#!/usr/bin/env python

import os
import re
import gzip
import json
import pysam
import subprocess
import numpy as np
import configparser
import pandas as pd
import pandas_plink
from collections import Counter

from flask import (
  url_for,
  redirect)

from __main__ import (
  lock,
  task_status)

from app.dirs import (
  root_dir,
  data_dir)

from modules.wpgma import wpgma

from modules.load_input import (
  load_gdata,
  load_config)

from modules.assoc_plot import (
  mplot,
  qqplot)

from app.task_worker import progress

tool_dir = os.path.join(
  root_dir,
  "modules",
  "tools")

# Bin size for making marker density histograms
mdbin = 1000000

###############################################################################
### Open a command-line process
###############################################################################

def popen(task_id, cmd):
    pr = subprocess.Popen(
      cmd,
      stdout=subprocess.PIPE,
      stderr=subprocess.STDOUT,
      text=True)
    for line in pr.stdout:
        with lock:
            task_status[task_id]["progress"] += line
    pr.wait()
    with lock:
        task_status[task_id]["progress"] += "\n"
    if pr.returncode != 0:
        raise subprocess.CalledProcessError(1, "")

###############################################################################
### Calculate marker set summary statistics
###############################################################################

def run_qc(genome, task_id, markers=None):
    config = load_config(genome)
    gdata = load_gdata(genome)

    plink_path = os.path.join(
      tool_dir,
      "plink2")

    if not markers:
        vcf_path = os.path.join(
          data_dir,
          genome,
          config["vcf"])

        stats_dir = os.path.join(
          data_dir,
          genome,
          "stats")

        plink_out = os.path.join(
          stats_dir,
          "gt")

    else:
        stats_dir = os.path.join(
          data_dir,
          genome,
          "assoc",
          markers,
          "stats")

        plink_out = os.path.join(
          data_dir,
          genome,
          "assoc",
          markers,
          "gt",
          "gt")

    if not os.path.exists(stats_dir):
        os.makedirs(stats_dir)

    assoc_dir = os.path.join(
      data_dir,
      genome,
      "assoc")
    if not os.path.exists(stats_dir):
        os.makedirs(stats_dir)

    if any((not os.path.exists(
      f"{plink_out}.{ext}")
      for ext in (
        "bed",
        "bim",
        "fam"))):
        # Create binary PLINK files
        progress(task_id, "Processing input genotype")
        popen(
          task_id,
          [plink_path,
           "--vcf", vcf_path,
           "--make-bed",
           "--set-all-var-ids", "@:#",
           "--allow-extra-chr",
           "--max-alleles", "2",
           "--polyploid-mode", "missing",
           "--vcf-half-call", "missing",           
           "--out", plink_out])

    if any((not os.path.exists(
      f"{plink_out}.{ext}")
      for ext in (
        "smiss",
        "vmiss",
        "het",
        "gcount",
        "afreq",
        "hardy"))):
        progress(task_id, "Calculating genotype statistics")
        popen(
          task_id,
          [plink_path,
           "--bfile", plink_out,
           "--missing",
           "--geno-counts",
           "--het",
           "--freq",
           "--hardy",
           "--allow-extra-chr",
           "--out", plink_out])

    metrics = (
      "smiss",
      "vmiss",
      "shet",
      "vhet",
      "maf",
      "hwe")
    exts = (
      "smiss",
      "vmiss",
      "het",
      "gcount",
      "afreq",
      "hardy")
    id_cols = (
      "IID",
      "ID",
      "IID",
      "ID",
      "ID",
      "ID")
    usecols = {
      "smiss": ["MISSING_CT", "OBS_CT"],
      "vmiss": ["MISSING_CT", "OBS_CT"],
      "shet": ["IID", "OBS_CT", "O(HOM)"],
      "vhet": ["ID", "HET_REF_ALT_CTS", "HOM_REF_CT", "TWO_ALT_GENO_CTS"],
      "maf": ["ALT_FREQS"],
      "hwe": ["P"]}

    progress(task_id, "Processing genotype statistics")
    for metric, ext, id_col in zip(metrics, exts, id_cols):
        qc_in = f"{plink_out}.{ext}"
        cnt = Counter()
        df = pd.read_csv(
          qc_in,
          sep="\t",
          chunksize=100000,
          usecols=usecols[metric])
        h5_out = os.path.join(
          stats_dir,
          f"{metric}.h5")

        for chk in df:
            if metric in ("smiss", "vmiss"):
                cnt.update(chk["MISSING_CT"] / chk["OBS_CT"])
            elif metric == "shet":
                chk["HET"] = (chk["OBS_CT"] - chk["O(HOM)"]) / chk["OBS_CT"]
                chk = chk.loc[:, ["IID", "HET"]]
                cnt.update(chk["HET"])
            elif metric == "vhet":
                chk["HET"] = (
                  chk["HET_REF_ALT_CTS"]
                  / (chk["HOM_REF_CT"]
                  + chk["HET_REF_ALT_CTS"]
                  + chk["TWO_ALT_GENO_CTS"]))
                chk = chk.loc[:, ["ID", "HET"]]
                cnt.update(chk["HET"])
            elif metric == "maf":
                chk["REF_FREQS"] = 1 - chk["ALT_FREQS"]
                cnt.update(chk.min(axis=1))
            elif metric == "hwe":
                cnt.update(-np.log10(chk["P"]))
            if metric in ("shet", "vhet"):
                chk.to_hdf(
                  h5_out,
                  key="het",
                  mode="a",
                  format="table",
                  append=True)

        cnt_vals = list(cnt.keys())
        weights = [cnt[v] for v in cnt_vals]
        hgts, bins = np.histogram(
          cnt_vals,
          weights=weights)
        # Get a larger number of bins
        cnt = {(-1 if np.isnan(v) else v): cnt[v] for v in cnt.keys()}
        cnt_vals = list(cnt.keys())
        weights = [cnt[v] for v in cnt_vals]
        hgts, bins = np.histogram(
          cnt_vals,
          bins=len(bins)*5,
          weights=weights)
        hgts = hgts.tolist()
        bins = bins.tolist()
        with open(os.path.join(stats_dir, metric), "w") as hist_out:
            json.dump([hgts, bins], hist_out)

        if metric == "smiss":
            samp_num = sum(weights)
        elif metric == "vmiss":
            var_num = sum(weights)

    # Keep full data for sample missingness and heterozygosity
    for metric, ext in (("smiss", "smiss"), ("shet", "het")):
        qc_in = f"{plink_out}.{ext}"
        df = pd.read_csv(
          qc_in,
          sep="\t",
          index_col=[0],
          usecols=["IID"] + usecols[metric])
        out_path = os.path.join(
          stats_dir,
          f"stats.{metric}")
        if metric == "smiss":
            ((df["MISSING_CT"] / df["OBS_CT"])
              .rename("SMISS")
              .to_csv(
                out_path,
                sep="\t",
                compression="gzip"))
        elif metric == "shet":
            (((df["OBS_CT"] - df["O(HOM)"]) / df["OBS_CT"])
              .rename("SHET")
              .to_csv(
                out_path,
                sep="\t",
                compression="gzip"))

    # Histogram of marker density
    progress(task_id, "Calculating marker density")
    concat = json.loads(config["concat"])
    ch_names = gdata["chr_names"]
    ch_sizes = dict(zip(ch_names, gdata["chr_sizes"]))
    ch_offsets = dict(zip(ch_names, gdata["chr_offsets"]))
    bin_offsets = dict()
    ch_names = [ch for ch in ch_names if ch not in concat]
    other_size = sum([ch_sizes.pop(ch) for ch in concat])
    for idx, ch in enumerate(ch_names):
        bin_offsets[ch] = sum([
          (ch_sizes[c] // mdbin * mdbin + mdbin)
          for c in ch_names[:idx]])
    if concat:
        bin_offsets["Other"] = sum([
          (ch_sizes[c] // mdbin * mdbin + mdbin)
          for c in ch_names])
        ch_sizes["Other"] = other_size
        ch_names.append("Other")

    hgts, bins, coords = list(), list(), list()
    bins.append(0)
    cnts = Counter()

    df = pd.read_csv(
      f"{plink_out}.bim",
      sep="\t",
      chunksize=100000,
      header=None,
      index_col=None,
      usecols=[0, 3])
    for chk in df:
        for ch, pos in chk.groupby(0):
            ch = str(ch)
            if ch in concat:
                offset = (
                  bin_offsets["Other"] +
                  sum(
                  [gdata["chr_sizes"][gdata["chr_names"].index(ch)]
                   for ch in concat[:concat.index(ch)]]))
                ch = "Other"
            else:
                offset = bin_offsets[ch]
            cnts.update(
              (pos[3] + offset)
              // mdbin
              * mdbin
              + mdbin)

    for ch in ch_names:
        size = ch_sizes[ch]
        for block in range(
          mdbin + bin_offsets[ch],
          size + bin_offsets[ch] + mdbin,
          mdbin):
            hgts.append(cnts.get(block, 0))
            bins.append(block)
            ch_block = block - bin_offsets[ch]
            coords.append(f"{ch}:{ch_block:,}-{ch_block + mdbin - 1:,}")

    with open(os.path.join(stats_dir, "dens"), "w") as dens_out:
        json.dump([hgts, bins, coords], dens_out)

    with open(os.path.join(stats_dir, "cnts"), "w") as cnts_out:
        json.dump([samp_num, var_num], cnts_out)

    # LD matrices per chromosome
    progress(task_id, "Calculating linkage disequilibrium")
    concat = json.loads(config["concat"])
    for ch in gdata["chr_names"]:
        if ch in concat:
            continue
        popen(
          task_id,
          [plink_path,
           "--bfile", plink_out,
           "--chr", str(ch),
           "--thin-count", "1000",
           "--r2-unphased", "square0",
           "--out", plink_out])

        r2_vals = pd.read_csv(
          f"{plink_out}.unphased.vcor2",
          sep="\t",
          header=None,
          index_col=None).round(3).fillna("")
        r2_idx = pd.read_csv(
          f"{plink_out}.unphased.vcor2.vars",
          sep="\t",
          header=None,
          index_col=None)
        r2_path = os.path.join(
          stats_dir,
          f"{ch}.r2")
        with gzip.open(r2_path, "wt") as r2_out:
            r2_out.write(
              json.dumps(
                [r2_vals.values.tolist(),
                 r2_idx[0].tolist()]))
    os.remove(f"{plink_out}.unphased.vcor2")
    os.remove(f"{plink_out}.unphased.vcor2.vars")

    # Clean up
    for fn in os.listdir(stats_dir):
        if (fn.startswith("gt")
          and not fn.endswith(
            (".gz",
             ".bed",
             ".bim",
             ".fam"))):
            os.remove(
              os.path.join(
                stats_dir,
                fn))

###############################################################################
### Apply filters to marker set
###############################################################################

def filt_qc(
  genome,
  data,
  task_id):
    gdata = load_gdata(genome)
    config = load_config(genome)
    vcf_file = config["vcf"]
    markers = data["name"]

    plink_path = os.path.join(
      tool_dir,
      "plink2")

    stats_dir = os.path.join(
      data_dir,
      genome,
      "stats")

    assoc_dir = os.path.join(
      data_dir,
      genome,
      "assoc",
      markers)

    filt_dir = os.path.join(
      assoc_dir,
      "gt")

    plink_in = os.path.join(
      stats_dir,
      "gt")

    plink_out = os.path.join(
      filt_dir,
      "gt")

    if not os.path.exists(filt_dir):
        os.makedirs(filt_dir)
    if not os.path.exists(stats_dir):
        os.makedirs(stats_dir)

    # Prepare heterozygosity filters
    if data["shet"]:
        shet_exc = os.path.join(
          filt_dir,
          "shet_exclude.txt")
        df = pd.read_hdf(
          os.path.join(stats_dir, "shet.h5"))
        df[df["HET"] > float(data["shet"])]["IID"].to_csv(
          shet_exc,
          index=False,
          header=False)
    if data["vhet"]:
        vhet_exc = os.path.join(
          filt_dir,
          "vhet_exclude.txt")
        df = pd.read_hdf(
          os.path.join(stats_dir, "vhet.h5"))
        df[df["HET"] > float(data["vhet"])]["ID"].to_csv(
          vhet_exc,
          index=False,
          header=False)

    # Prune by LD
    if data["prune"]:
        progress(task_id, "Pruning markers by LD")
        data["window"] = (
          str(data["window"]) + "kb"
          if str(data["window"]) else "200kb")

        popen(task_id,
          [plink_path,
          "--bfile", plink_in,
          "--indep-pairwise", data["window"], str(data["prune"]),
          "--allow-extra-chr",
          "--out", plink_out])

        prune_out_fn = plink_out + ".prune.out"
        with open(prune_out_fn, "r") as prune_out_in:
            prune_markers = prune_out_in.read().strip().split("\n")
            if vhet_exc:
                prune_rmv = len(
                  set(prune_markers).difference(
                    set(vhet_exc)))
            else:
                prune_rmv = len(prune_markers)

    else:
        prune_rmv = ""

    # Filter based on thresholds
    cmd = [plink_path, "--bfile", plink_in]

    cmd = cmd + ["--mind", str(data["smiss"])] if data["smiss"] else cmd
    cmd = cmd + ["--geno", str(data["vmiss"])] if data["vmiss"] else cmd

    cmd = cmd + ["--remove", shet_exc] if data["shet"] else cmd

    if vhet_exc or data["prune"]:
        cmd += ["--exclude"]
        cmd = cmd + [vhet_exc] if data["vhet"] else cmd
        cmd = cmd + [prune_out_fn] if data["prune"] else cmd

    cmd = cmd + ["--maf", str(data["maf"])] if data["maf"] else cmd
    cmd = cmd + ["--hwe", str(10**(-data["hwe"]))] if data["hwe"] else cmd

    cmd = cmd + ["--allow-extra-chr",
                 "--export", "vcf", "bgz", "id-paste=iid",
                 "--out", plink_out]

    progress(task_id, "Applying filters")
    popen(task_id, cmd)

    # Impute missing genotypes
    beagle_path = os.path.join(
      tool_dir,
      "beagle.22Jul22.46e.jar")

    vcf_in = plink_out + ".vcf.gz"
    vcf_root = os.path.join(
      filt_dir,
      "imp.gt")
    vcf_out = f"{vcf_root}.vcf.gz"

    if data["imp"] or 1:
        progress(task_id, "Imputing missing genotypes")
        popen(
          task_id,
          ["java",
           "-jar",
           beagle_path,
           f"gt={vcf_in}",
           f"out={vcf_root}"])

        os.remove(vcf_in)
        os.rename(vcf_out, vcf_in)

    vcf_out = vcf_in
    final_root = plink_out

    if data["hap"]:
        progress(task_id, "Constructing haplotype blocks")
        blocktool_path = os.path.join(
          tool_dir,
          "blocktool.py")
        # Make necessary genotype format files
        popen(
          task_id,
          [plink_path,
           "--vcf", vcf_out,
           "--allow-extra-chr",
           "--make-bed",
           "--export", "tped",
           "--out", final_root])

        # Make haplotype blocks
        popen(
          task_id,
          [blocktool_path,
           "block",
           "--type", "gab",
           "--gab-max-kb", "1000",
           "--gt-dir", filt_dir])

        # Genotype gased on blocks
        popen(
          task_id,
          [blocktool_path,
           "call",
           "--format", "per-hap",
           "--block-dir", filt_dir,
           "--gt-dir", filt_dir,
           "--transpose", "per-hap",
           "--maf", "0.05",
           "--min-hap", "1",
           "--tab"])

        os.remove(f"{final_root}.tped")
        os.remove(f"{final_root}.tfam")

    else:
        progress(task_id, "Processing filtered data")
        # Make necessary genotype format files
        popen(
          task_id,
          [plink_path,
           "--vcf", vcf_out,
           "--make-bed",
           "--allow-extra-chr",
           "--out", final_root])

    # Calculate K and Q
    progress(task_id, "Determining kinship and population structure")
    popen(
      task_id,
      [plink_path,
       "--bfile", final_root,
       "--make-rel", "square",
       "--pca", "15",
       "--allow-extra-chr",
       "--out", final_root])

    k_matr = pd.read_csv(
      f"{final_root}.rel",
      sep="\t",
      header=None,
      index_col=None)
#    k_idx = pd.read_csv(
#      f"{final_root}.rel.id",
#      sep="\t",
#      header=[0],
#      index_col=None)
#    k_matr.set_index(k_idx["IID"]).to_csv(
#      f"{final_root}.rel",
#      sep="\t",
#      header=False)

    gt_stats_dir = os.path.join(
      data_dir,
      genome,
      "assoc",
      markers,
      "stats")

    if not os.path.exists(gt_stats_dir):
        os.makedirs(gt_stats_dir)

    # Convert PCA eigenvectors to PVE
    sum_diag = sum(np.diag(k_matr))
    eigenval = pd.read_csv(
      f"{final_root}.eigenval",
      header    = None,
      index_col = None)
    eigenval = eigenval / sum_diag
    eigenval[0].to_json(
      os.path.join(
        gt_stats_dir,
        "pca"),
      orient="values")

    # Cluster kinship matrix
    plink1_path = os.path.join(
      tool_dir,
      "plink1.9")
    hamm_out = os.path.join(
      gt_stats_dir,
      "kin")

    popen(
      task_id,
      [plink1_path,
       "--bfile", final_root,
       "--distance", "square",
       "--allow-extra-chr",
       "--out", hamm_out])

    hamm_matrix = np.loadtxt(
      f"{hamm_out}.dist",
      delimiter="\t")
    hamm_labels = pd.read_csv(
      f"{hamm_out}.dist.id",
      sep="\t",
      header=None)[1].tolist()
    hamm_df = pd.DataFrame(hamm_matrix)
    hamm_df = hamm_df.set_axis(hamm_labels, axis=0)
    hamm_df = hamm_df.set_axis(hamm_labels, axis=1)

    tree, lab_order = wpgma(
      hamm_matrix,
      hamm_labels)

    hamm_df = hamm_df.loc[lab_order]
    hamm_df = hamm_df.loc[:, lab_order]
    with open(f"{hamm_out}", "w") as matr_out:
        json.dump(
          [hamm_df.values.tolist(),
           lab_order],
          matr_out)
    # Clean up
    for fn in (
      "kin.dist",
      "kin.dist.id",
      "kin.log",
      "kin.nosex"):
        os.remove(
          os.path.join(
            gt_stats_dir,
            fn))

    progress(task_id, "Saving summary")
    # Create file describing the filtering parameters
    params_fn = os.path.join(
      gt_stats_dir,
      "params")

    # Get original number of markers and samples
    (bim, fam, bed) = pandas_plink.read_plink(plink_in)
    var_all, samp_all = bed.shape
    (bim, fam, bed) = pandas_plink.read_plink(final_root)
    var_filt, samp_filt = bed.shape
    # Save parameters and counts in a file
    with open(params_fn, "w") as params_out:
        params = [
          ["",
           "Sample heterozygosity",
           "Sample missingness",
           "Variant heterozygosity",
           "LD window",
           "Variant LD",
           "Variant missingness",
           "Variant -log10(HWE)",
           "Variant MAF",
           "Markers retained/original",
           "Samples retained/original"],
          ["Threshold",
           float(data["shet"]) if data["shet"] else "-",
           float(data["smiss"]) if data["smiss"] else "-",
           float(data["vhet"]) if data["vhet"] else "-",
           data["window"] if data["window"] else "-",
           float(data["prune"]) if data["prune"] else "-",
           float(data["vmiss"]) if data["vmiss"] else "-",
           float(data["hwe"]) if data["hwe"] else "-",
           float(data["maf"]) if data["maf"] else "-",
           f"{var_filt}/{var_all}",
           f"{samp_filt}/{samp_all}"]]
        if data["hap"]:
            params[0].append("Haplotype count")
            params[1].append("-")
        json.dump(params, params_out)

    # Clean up
    for fn in os.listdir(filt_dir):
        if (fn.startswith("gt")
          and not fn.endswith(
            (".gz",
             ".bed",
             ".bim",
             ".fam",
             ".txt",
             ".eigenvec",
             ".rel"))):
            os.remove(os.path.join(filt_dir, fn))

    run_qc(
      genome,
      task_id,
      markers = data["name"])

    # Change chromosomes to numeric (required by BLINK)
    chr_rename = {
      ch: idx for idx, ch in enumerate(gdata["chr_names"], start=1)}
    (pd
      .read_csv(
        f"{final_root}.bim",
        sep="\t",
        index_col=None,
        header=None)
      .set_axis(
        ["CHR", "ID", "FAM",
         "POS", "REF", "ALT"],
        axis=1)
      .assign(CHR=lambda x: x["CHR"].astype(str).map(chr_rename))
      .to_csv(
        f"{final_root}.bim",
        sep="\t",
        header=False,
        index=False))

    init_path = os.path.join(
      assoc_dir,
      "init")
    open(init_path, "w").close()


###############################################################################
### Run association analysis and process output
###############################################################################

def run_assoc(genome, params, task_id):
    config  = load_config(genome)
    markers = params["markers_sgl"]
    model   = params["model_sgl"]

    progress(task_id, "Preparing input files")
    assoc_dir = os.path.join(
      data_dir,
      genome,
      "assoc",
      markers)

    out_dir = os.path.join(
      assoc_dir,
      "res")

    in_root = os.path.join(
      assoc_dir,
      "gt",
      "gt")

    phe_path = os.path.join(
      data_dir,
      genome,
      config["phe"])

    kq_root = os.path.join(
      data_dir,
      genome,
      "assoc",
      params["kq_from"] or markers,
      "gt",
      "gt")

    if model == "blink":
        # BLINK analyses either one or all traits
        # phenotype file rewritten each time to keep only selected traits
        out_fam = f"{in_root}.fam"
        fam_df = (pd
          .read_csv(
            out_fam,
            sep="\t",
            header=None,
            index_col=None)
          .set_axis(
            ["FID", "IID", "M", "P", "S", "PHE"],
            axis=1))
        phe_path = os.path.join(
          data_dir,
          genome,
          config["phe"])
        phe_df = pd.read_csv(
          phe_path,
          sep="\t",
          header=[0],
          index_col=[0])
        (phe_df
          .fillna("NaN")
          .loc[
            fam_df["IID"],
            params["traits"]]
          .to_csv(
            f"{in_root}.txt",
            sep="\t"))

        # Covariates file created to only keep top n selected ones
        q_matr = pd.read_csv(
            f"{kq_root}.eigenvec",
            sep="\t",
            header=[0],
            index_col=[0,1])
        (q_matr
#          .assign(I=1)
          .loc[:, q_matr.columns[:int(params["q_num"])]]
          .droplevel(0)
          .to_csv(
            f"{in_root}.cov",
            sep="\t"))

        blink_path = os.path.join(
          tool_dir,
          "BLINK",
          "blink")

        progress(task_id, "Running analysis")
        command = [
          blink_path,
          "--gwas",
          "--file", in_root,
          "--plink",
          "--out", os.path.join(out_dir, "BLINK"),
          "--trait", "0"]
        popen(
          task_id,
          command)

        for trait in params["traits"]:
            assoc_process(
              genome,
              f"BLINK_{trait}_GWAS_result.txt",
              f"{model}.{trait}",
              out_dir,
              ch_col="chr",
              pos_col="pos",
              pval_col="p_value",
              eff_col=None)

            # Save parameters
            run_params = [
              ["Markers", "Samples", "PCs", "Kinship", "K/Q source"],
              [sum(1 for l in open(f"{in_root}.bim")),
               phe_df.loc[fam_df["IID"]][trait].dropna().shape[0],
               params["q_num"],
               "No",
               params["kq_from"] or markers]]
            with open(
              os.path.join(
                out_dir,
                f"{model}.{trait}.params"),
              "w") as params_out:
                json.dump(run_params, params_out)

        os.remove(os.path.join(out_dir, "BLINK.log"))

    elif model == "rf":
        pass

    else:
        # GEMMA cannot accept row or column names in phenotype file
        # phenotype file rewritten each time to keep only selected traits
        out_fam = f"{in_root}.fam"
        fam_df = (pd
          .read_csv(
            out_fam,
            sep="\t",
            header=None,
            index_col=None)
          .set_axis(
            ["FID", "IID", "M", "P", "S", "PHE"],
            axis=1))
        phe_path = os.path.join(
          data_dir,
          genome,
          config["phe"])
        phe_df = pd.read_csv(
          phe_path,
          sep="\t",
          header=[0],
          index_col=[0])
        (phe_df
          .fillna("NA")
          # GEMMA treats -9 as a missing value                                    #################################
          .replace(-9, -9.000000000001)
          .replace("-9", "-9.000000000001")
          .loc[fam_df["IID"], ]
          .to_csv(
            f"{in_root}.txt",
            sep="\t",
            header=False,
            index=False))

        # Covariates file created to only keep top n selected ones
        q_matr = pd.read_csv(
            f"{kq_root}.eigenvec",
            sep="\t",
            header=[0],
            index_col=[0,1])
        (q_matr
          .assign(I=1)
          .loc[:, ["I"] + q_matr.columns.tolist()[:int(params["q_num"])]]
          .to_csv(
            f"{in_root}.cov",
            sep="\t",
            header=False,
            index=False))

        progress(task_id, "Running analysis")
        if model == "mtmlm":
            traits = params["trait_groups"]
        else:
            traits = zip([params["traits"]], params["traits"])
        for trait, trait_name in traits:
            n = [str(phe_df.columns.get_loc(t) + 1) for t in trait]
            command = [
              "gemma",
              "-bfile", in_root,
              "-p", f"{in_root}.txt",
              "-k", f"{kq_root}.rel",
              "-c", f"{in_root}.cov",
              "-n", *n,
              "-o", f"{model}.{trait_name}",
              "-outdir", out_dir]
            if model == "glm":
                command += ["-lm", "1"]
            else:
                command += ["-lmm", "1"]
            popen(
              task_id,
              command)

            assoc_process(
              genome,
              f"{model}.{trait_name}.assoc.txt",
              f"{model}.{trait_name}",
              out_dir,
              ch_col="chr",
              pos_col="ps",
              pval_col="p_wald",
              eff_col=None)

            # Save parameters
            run_params = [
              ["Markers", "Samples", "PCs", "Kinship", "K/Q source"],
              [sum(1 for l in open(f"{in_root}.bim")),
               phe_df.loc[fam_df["IID"]][trait].dropna().shape[0],
               params["q_num"],
               "No" if model == "glm" else "Yes",
               params["kq_from"] or markers]]
            if model == "mtmlm":
                run_params[0].append("Traits")
                run_params[1].append("<br>".join(trait))
            with open(
              os.path.join(
                out_dir,
                f"{model}.{trait_name}.params"),
              "w") as params_out:
                json.dump(run_params, params_out)

            os.remove(os.path.join(out_dir, f"{model}.{trait_name}.log.txt"))

def assoc_process(
  genome,
  assoc_fn,
  out_root,
  out_dir,
  ch_col,
  pos_col,
  pval_col,
  eff_col):
    gdata = load_gdata(genome)
    config = load_config(genome)
    assoc_out = os.path.join(
      out_dir,
      assoc_fn)
    processed_out = os.path.join(out_dir, out_root)
    columns = [c for c in (ch_col, pos_col, pval_col, eff_col) if c]

    # Output files
    mplot_fn    = f"{processed_out}.mplot.png"
    qqplot_fn   = f"{processed_out}.qqplot.png"
    out_feather = f"{processed_out}.feather"
    out_tsv     = f"{processed_out}.tsv"
    out_gz      = f"{processed_out}.tsv.gz"

    # Save raw results, rename numeric chromosomes back to original
    chr_rename = dict(enumerate(gdata["chr_names"], start=1))
    res = (
      pd
        .read_csv(
          assoc_out,
          sep="\t",
          header=[0],
          index_col=None)
        .loc[:, columns]
        .assign(**{ch_col: lambda x: x[ch_col].map(chr_rename)}))
    res.to_csv(
      out_tsv,
      sep="\t",
      index=False)
    pysam.tabix_compress(out_tsv, out_gz, force=True)
    pysam.tabix_index(
      out_gz,
      seq_col=0,
      start_col=1,
      end_col=1,
      line_skip=1, force=True)

    # Save results for MP visualisation
    res[pval_col] = -np.log10(res[pval_col])
    (res
      .sort_values(
        by=pval_col,
        ascending=False)
      .head(5000)
      .to_feather(out_feather))

    mplot_fig = mplot(
      genome,
      res,
      json.loads(config["concat"]),
      ch_col=ch_col,
      pos_col=pos_col,
      pval_col=pval_col)
    mplot_fig.savefig(
      mplot_fn,
      dpi=300,
      format="png")

    qqplot_fig = qqplot(
      res,
      pval_col=pval_col)
    qqplot_fig.savefig(
      qqplot_fn,
      dpi=300,
      format="png")

    os.remove(out_tsv)
    os.remove(assoc_out)