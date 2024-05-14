#!/usr/bin/env python

import os
import re
import gzip
import json
import subprocess
import numpy as np
import pandas as pd
from collections import Counter

from flask import (
  url_for,
  redirect)

from modules.wpgma import wpgma

from modules.load_input import (
  load_gdata,
  load_config)

from modules.assoc_plot import (
  mplot,
  qqplot)

from app.dirs import (
  root_dir,
  data_dir)

tool_dir = os.path.join(
  root_dir,
  "modules",
  "tools")
R_dir = os.path.join(
  root_dir,
  "modules",
  "R")

###############################################################################
### Calculate marker set summary statistics
###############################################################################

def run_qc(genome, markers=None):
    #yield "Summarizing genotype data"
    gconfig = load_config(genome)
    gdata = load_gdata(genome)

    plink_path = os.path.join(
      tool_dir,
      "plink2")

    if not markers:
        vcf_path = os.path.join(
          data_dir,
          genome,
          gconfig["vcf"])

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
        #yield "Processing input genotypes"
        # Create binary PLINK files
        pr = subprocess.Popen(
          [plink_path,
           "--vcf", vcf_path,
           "--make-bed",
           "--set-all-var-ids", "@:#",
           "--allow-extra-chr",
           "--out", plink_out],
          stdout=subprocess.PIPE,
          stderr=subprocess.PIPE)
        out, err = pr.communicate()
        if pr.returncode != 0:
            raise OSError(err.decode())

    if any((not os.path.exists(
      f"{plink_out}.{ext}")
      for ext in (
        "smiss",
        "vmiss",
        "het",
        "gcount",
        "afreq",
        "hardy"))):
        #yield "Calculating genotype statistics"
        pr = subprocess.Popen(
          [plink_path,
           "--bfile", plink_out,
           "--missing",
           "--geno-counts",
           "--het",
           "--freq",
           "--hardy",
           "--allow-extra-chr",
           "--out", plink_out],
          stdout=subprocess.PIPE,
          stderr=subprocess.PIPE)
        out, err = pr.communicate()
        if pr.returncode != 0:
            raise OSError(err.decode())

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

    #yield "Processing genotype statistics"
    for metric, ext, id_col in zip(metrics, exts, id_cols):
        qc_in = f"{plink_out}.{ext}"
        qc_df = pd.read_csv(
          qc_in,
          sep="\t")

        if metric in ("smiss", "vmiss"):
            qc_df["MIS_F"] = (qc_df["MISSING_CT"]
                              / qc_df["OBS_CT"])
            data = qc_df[[id_col, "MIS_F"]].sort_values(
              by="MIS_F",
              ascending=False)
        elif metric == "shet":
            qc_df["HET_F"] = ((qc_df["OBS_CT"]
                               - qc_df["O(HOM)"])
                              / qc_df["OBS_CT"])
            data = qc_df[[id_col, "HET_F"]].sort_values(
              by="HET_F",
              ascending=False)
        elif metric == "vhet":
            qc_df["HET_F"] = (qc_df["HET_REF_ALT_CTS"]
                              / (qc_df["HOM_REF_CT"]
                                 + qc_df["HET_REF_ALT_CTS"]
                                 + qc_df["TWO_ALT_GENO_CTS"]))
            data = qc_df[[id_col, "HET_F"]].sort_values(
            by="HET_F",
            ascending=False)
        elif metric == "maf":
            qc_df["REF_FREQS"] = 1 - qc_df["ALT_FREQS"]
            qc_df["MAF"] = qc_df[["REF_FREQS", "ALT_FREQS"]].min(axis=1)
            data = qc_df[[id_col, "MAF"]].sort_values(
              by="MAF")
        elif metric == "hwe":
            qc_df["log10P"] = -np.log10(qc_df["P"])
            data = qc_df[[id_col, "log10P"]].sort_values(
              by="log10P",
              ascending=False)

        data.iloc[:, 1] = data.iloc[:, 1].round(5)
        data.to_csv(
          f"{plink_out}.{metric}.csv.gz",
          index  = False,
          header = False,
          compression="gzip")

        hgts, bins = np.histogram(data.iloc[:, 1])
        with open(os.path.join(stats_dir, metric), "w") as hist_out:
            json.dump([hgts.tolist(), bins.tolist()], hist_out)

        if metric == "smiss":
            samp_num = data.shape[0]
        elif metric == "vmiss":
            var_num = data.shape[0]

    #yield "Calculating marker density"
    # Histogram of marker density
    chr_vis = gconfig["chr_vis"]
    if chr_vis:
        chr_vis = [
          ch.strip()
          for ch in chr_vis.split(",")]

    hgts, bins, coords = list(), list(), list()
    bins.append(0)
    ids = data["ID"].str.split(":", expand=True)
    for ch, pos in ids.groupby(0):
        if chr_vis and ch not in chr_vis:
            continue
        ch_idx = gdata["chr_names"].index(ch)
        ch_size = gdata["chr_sizes"][ch_idx]
        ch_offset = gdata["chr_offsets"][ch_idx]
        if ch_offset != 0:
            ch_offset = (
              ch_offset
              // 1000000
              * 1000000
              + 1000000)
        pos = (
          pos[1].astype(int)
          // 1000000
          * 1000000
          + 1000000)
        cnts = Counter(pos)
        max_block = (
          ch_size
          // 1000000
          * 1000000
          + 2000000)
        for block in range(
          1000000,
          max_block,
          1000000):
            hgts.append(cnts.get(block, 0))
            bins.append(block + ch_offset)
            coords.append(f"{ch}:{block - 1000000}-{block - 1}")
    if chr_vis:
        pos = ids[~ids[0].isin(chr_vis)]
        hgts.append(pos.shape[0])
        bins.append(
          max(
            ((gdata["chr_sizes"][-1]
              + gdata["chr_offsets"][-1])
             // 1000000
             * 1000000
             + 1000000),
            bins[-1] + 1000000))
        coords.append("Other")

    with open(os.path.join(stats_dir, "dens"), "w") as dens_out:
        json.dump([hgts, bins, coords], dens_out)

    with open(os.path.join(stats_dir, "cnts"), "w") as cnts_out:
        json.dump([samp_num, var_num], cnts_out)

    #yield "Calculating linkage disequilibrium"
    # LD matrices per chromosome
    for ch in gconfig["chr_vis"] or gdata["chr_names"]:
        pr = subprocess.Popen(
          [plink_path,
           "--bfile", plink_out,
           "--chr", str(ch),
           "--thin-count", "1000",
           "--r2-unphased", "square0",
           "--out", plink_out],
          stdout=subprocess.PIPE,
          stderr=subprocess.PIPE)
        out, err = pr.communicate()
        if pr.returncode != 0:
            raise OSError(err.decode())

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
            os.remove(os.path.join(stats_dir, fn))

###############################################################################
### Apply filters to marker set
###############################################################################

def filt_qc(
  genome,
  data,
  shet_exc,
  vhet_exc):
    gdata = load_gdata(genome)
    gconfig = load_config(genome)
    vcf_file = gconfig["vcf"]
    markers  = data["name"]

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

    # Get samples and variants to filter out by heterozygosity
    if shet_exc:
        shet_out_fn = os.path.join(
          filt_dir,
          "shet_exclude.txt")
        with open(shet_out_fn, "w") as shet_out:
            for s in shet_exc:
                _ = shet_out.write(s + "\n")

    if vhet_exc:
        vhet_out_fn = os.path.join(
          filt_dir,
          "vhet_exclude.txt")
        with open(vhet_out_fn, "w") as vhet_out:
            for v in vhet_exc:
                _ = vhet_out.write(v + "\n")

    # Prune by LD
    if data["prune"]:
        #yield "Pruning input markers"
        data["window"] = (
          str(data["window"]) + "kb"
          if str(data["window"]) else "200kb")

        pr = subprocess.Popen(
          [plink_path,
          "--bfile", plink_in,
          "--indep-pairwise", data["window"], str(data["prune"]),
          "--allow-extra-chr",
          "--out", plink_out],
         stdout=subprocess.PIPE,
         stderr=subprocess.PIPE)
        out, err = pr.communicate()
        if pr.returncode != 0:
            raise OSError(err.decode())

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

    cmd = cmd + ["--remove", shet_out_fn] if shet_exc else cmd

    if vhet_exc or data["prune"]:
        cmd += ["--exclude"]
        cmd = cmd + [vhet_out_fn] if vhet_exc else cmd
        cmd = cmd + [prune_out_fn] if data["prune"] else cmd

    cmd = cmd + ["--maf", str(data["maf"])] if data["maf"] else cmd
    cmd = cmd + ["--hwe", str(10**(-data["hwe"]))] if data["hwe"] else cmd

    cmd = cmd + ["--allow-extra-chr",
                 "--export", "vcf", "bgz", "id-paste=iid",
                 "--out", plink_out]

    #yield "Applying filters"
    pr = subprocess.Popen(
      cmd,
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE)
    out, err = pr.communicate()
    if pr.returncode != 0:
        raise OSError(err.decode())

    # Impute missing genotypes
    beagle_path = os.path.join(
      tool_dir,
      "beagle.22Jul22.46e.jar")

    vcf_in = plink_out + ".vcf.gz"
    vcf_root = os.path.join(
      filt_dir,
      "imp.gt")
    vcf_out = f"{vcf_root}.vcf.gz"

    if data["imp"]:
        #yield "Imputing missing genotypes"
        pr = subprocess.Popen(
          ["java",
           "-jar",
           beagle_path,
           f"gt={vcf_in}",
           f"out={vcf_root}"],
          stdout=subprocess.PIPE,
          stderr=subprocess.PIPE)
        out, err = pr.communicate()
        if pr.returncode != 0:
            raise OSError(err.decode())

        os.remove(vcf_in)
        os.rename(vcf_out, vcf_in)

    vcf_out = vcf_in
    final_root = plink_out

    if data["hap"]:
        #yield "Constructing haplotype blocks"
        blocktool_path = os.path.join(
          tool_dir,
          "blocktool.py")
        # Make necessary genotype format files
        pr = subprocess.Popen(
          [plink_path,
           "--vcf", vcf_out,
           "--allow-extra-chr",
           "--make-bed",
           "--export", "tped",
           "--out", final_root],
          stdout=subprocess.PIPE,
          stderr=subprocess.PIPE)
        out, err = pr.communicate()
        if pr.returncode != 0:
            raise OSError(err.decode())

        # Make haplotype blocks
        pr = subprocess.Popen(
          [blocktool_path,
           "block",
           "--type", "gab",
           "--gab-max-kb", "1000",
           "--gt-dir", filt_dir],
          stdout=subprocess.PIPE,
          stderr=subprocess.PIPE)
        out, err = pr.communicate()
        if pr.returncode != 0:
            raise OSError(err.decode())

        # Genotype gased on blocks
        pr = subprocess.Popen(
          [blocktool_path,
           "call",
           "--format", "per-hap",
           "--block-dir", filt_dir,
           "--gt-dir", filt_dir,
           "--transpose", "per-hap",
           "--maf", "0.05",
           "--min-hap", "1",
           "--tab"],
          stdout=subprocess.PIPE,
          stderr=subprocess.PIPE)
        out, err = pr.communicate()
        if pr.returncode != 0:
            raise OSError(err.decode())

        os.remove(f"{final_root}.tped")
        os.remove(f"{final_root}.tfam")
    
    else:
        #yield "Processing filtered data"
        # Make necessary genotype format files
        pr = subprocess.Popen(
          [plink_path,
           "--vcf", vcf_out,
           "--make-bed",
           "--allow-extra-chr",
           "--export", "A", "ped",
           "--out", final_root],
          stdout=subprocess.PIPE,
          stderr=subprocess.PIPE)
        out, err = pr.communicate()
        if pr.returncode != 0:
            raise OSError(err.decode())

        os.rename(
          f"{final_root}.raw",
          f"{final_root}.varA.raw")
        # Zip genotype files
        with open(f"{final_root}.varA.raw", "r") as raw_in:
            with gzip.open(f"{final_root}.raw.gz", "wt") as raw_out:
                header = next(raw_in)
                header = [
                  var.rsplit("_", 1)[0]
                  for var in header.split("\t")]
                _ = raw_out.write("\t".join(header) + "\n")
                for line in raw_in:
                    _ = raw_out.write(line)

        with open(f"{final_root}.map", "r") as map_in:
            with gzip.open(f"{final_root}.map.gz", "wt") as map_out:
                _ = map_out.writelines(map_in)

        # Clean up
        os.remove(f"{final_root}.map")
        os.remove(f"{final_root}.ped")
        os.remove(f"{final_root}.varA.raw")
        os.remove(vcf_out)

    #yield "Calculating kinship and population structure"
    # Calculate K and Q
    pr = subprocess.Popen(
      [plink_path,
       "--bfile", final_root,
       "--make-rel", "square",
       "--pca", "15",
       "--allow-extra-chr",
       "--out", final_root],
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE)
    out, err = pr.communicate()
    if pr.returncode != 0:
        raise OSError(err.decode())

    k_matr = pd.read_csv(
      f"{final_root}.rel",
      sep="\t",
      header=None,
      index_col=None)
    k_idx = pd.read_csv(
      f"{final_root}.rel.id",
      sep="\t",
      header=[0],
      index_col=None)
    k_matr.set_index(k_idx["IID"]).to_csv(
      f"{final_root}.rel",
      sep="\t",
      header=False)

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

    pr = subprocess.Popen(
      [plink1_path,
       "--bfile", final_root,
       "--distance", "square",
       "--allow-extra-chr",
       "--out", hamm_out],
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE)
    out, err = pr.communicate()
    if pr.returncode != 0:
        raise OSError(err.decode())

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

    #yield "Saving filtering summary"
    # Create file describing the filtering parameters
    params_fn = os.path.join(
      gt_stats_dir,
      "params")
    with open(params_fn, "w") as params_out:
        with open(f"{plink_out}.log", "r") as log_in:
            log_lines = log_in.read().split(".\n")

            # Exclude header lines, keep only lines describing filters
            _ = log_lines.pop(0)

            samp_ln = [
              line for line in log_lines 
              if "samples" in line and "loaded" in line]
            samp_all = (int(
                          samp_ln[0].split(" ", 1)[0])
                          if len(samp_ln) == 1
                          else "")

            vars_ln = [
              line for line in log_lines 
              if "variants" in line and "loaded" in line]
            vars_all = (int(
                          vars_ln[0].split(" ", 1)[0])
                          if len(vars_ln) == 1
                          else "")

            samp_rem = samp_all
            vars_rem = vars_all
            if samp_all and vars_all:
                shet_rem = samp_rem - len(shet_exc) if shet_exc else samp_rem
                vhet_rem = vars_rem - len(vhet_exc) if vhet_exc else vars_rem
                samp_rem, vars_rem  = shet_rem, vhet_rem

                prune_rem = vars_rem - prune_rmv if prune_rmv else vars_rem
                vars_rem  = prune_rem

                smiss_line = [
                  line for line in log_lines 
                  if "--mind" in line]
                smiss_rem = (
                  samp_rem
                  - int(re.search("\d+", smiss_line[0])[0])
                  if len(smiss_line) == 1 else "")
                samp_rem = smiss_rem or samp_rem

                vmiss_line = [
                  line for line in log_lines 
                  if "--geno" in line]
                vmiss_rem = (
                  vars_rem
                  - int(re.search("\d+", vmiss_line[0])[0])
                  if len(vmiss_line) == 1 else "")
                vars_rem = vmiss_rem or vars_rem

                hwe_line = [
                  line for line in log_lines 
                  if "--hwe" in line]
                hwe_rem = (
                  vars_rem
                  - int(re.search("\d+", hwe_line[0])[0])
                  if len(hwe_line) == 1 else "")
                vars_rem = hwe_rem or vars_rem

                maf_line = [
                  line for line in log_lines 
                  if "--maf" in line]
                maf_rem = (
                  vars_rem
                  - int(re.search("\d+", maf_line[0])[0])
                  if len(maf_line) == 1 else "")
                vars_rem = maf_rem or vars_rem

                samp_ln = [
                  line for line in log_lines 
                  if "samples" in line
                  and "remaining" in line
                  and "filters" in line]
                samp_rem = (int(
                              samp_ln[0].split(" ", 1)[0])
                              if len(samp_ln) == 1 else "")

                vars_ln = [
                  line for line in log_lines 
                  if "variants" in line
                  and "remaining" in line
                  and "filters" in line]
                vars_rem = (int(
                              vars_ln[0].split(" ", 1)[0])
                              if len(vars_ln) == 1 else "")

                params = list()
                params.append(
                  ["",
                   "Initial sample count",
                   "Sample heterozygosity",
                   "Sample missingness",
                   "Initial marker count",
                   "Variant heterozygosity",
                   "LD window",
                   "Variant LD",
                   "Variant missingness",
                   "Variant -log10(HWE)",
                   "Variant MAF",
                   "Final sample count",
                   "Final marker count"])

                params.append(
                  ["Threshold",
                   "-",
                   float(data["shet"]) if data["shet"] else "-",
                   float(data["smiss"]) if data["smiss"] else "-",
                   "-",
                   float(data["vhet"]) if data["vhet"] else "-",
                   data["window"] if data["window"] else "-",
                   float(data["prune"]) if data["prune"] else "-",
                   float(data["vmiss"]) if data["vmiss"] else "-",
                   float(data["hwe"]) if data["hwe"] else "-",
                   float(data["maf"]) if data["maf"] else "-",
                   "-",
                   "-"])

                params.append(
                  ["Remaining",
                   samp_all,
                   shet_rem,
                   smiss_rem,
                   vars_all,
                   vhet_rem,
                   "-",
                   prune_rem,
                   vmiss_rem, 
                   hwe_rem,
                   maf_rem,
                   samp_rem,
                   vars_rem])

                if data["hap"]:
                    params[0].append("Haplotype count")
                    params[1].append("-")
                    with open(f"{final_root}.blockinfo", "r") as blocks_in:
                        block_cnt = sum([1 for l in blocks_in]) - 1
                    params[2].append(block_cnt)

                json.dump(params, params_out)

            else:
                json.dump([], params_out)

    # Clean up
    for fn in os.listdir(filt_dir):
        if (fn.startswith("gt")
          and not fn.endswith(
            (".gz",
             ".bed",
             ".bim",
             ".fam",
             ".eigenvec",
             ".rel"))):
            os.remove(os.path.join(filt_dir, fn))

    #yield "Finished applying filters"
    run_qc(
      genome,
      markers = data["name"])
    return "Done."

###############################################################################
### Run association analysis and process output
###############################################################################

def run_assoc(genome, params):
    gconfig = load_config(genome)
    chr_vis = gconfig["chr_vis"]
    if chr_vis:
        chr_vis = [
          ch.strip() for ch in chr_vis.split(",")]

    markers = params["markers_sgl"]
    model   = params["model_sgl"]
    traits  = params["traits"][:]

    if not markers or not model or not traits:
        return

    gapit_path = os.path.join(
      R_dir,
      "GAPIT.R")

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
      gconfig["phe"])

    kq_root = os.path.join(
      data_dir,
      genome,
      "assoc",
      params["kq_from"] or markers,
      "gt",
      "gt")

    command = [
      "Rscript", gapit_path,
      "-g", f"{in_root}.raw.gz",
      "-m", f"{in_root}.map.gz",
      "-p", phe_path,
      "--phenotype_tsv",
      "-q", f"{kq_root}.eigenvec",
      "--q_tsv",
      "--top_q", str(params["q_num"]),
      "-k", f"{kq_root}.rel",
      "--k_tsv",
      "-o", out_dir,
      "-l", model,
      "-t", ", ".join(traits)]
    if params["chrom"]:
        command += ["-c", params["chrom"]]
    if params["start"]:
        command += ["-s", params["start"]]
    if params["end"]:
        command += ["-e", params["end"]]

    pr = subprocess.Popen(
      command,
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE)
    out, err = pr.communicate()

    raw_dir = os.path.join(out_dir, "raw")
    traits = traits if isinstance(traits, list) else [traits]
    for trait in traits:
        out_root = f"{model}.{trait}"

        in_fn = os.path.join(
          raw_dir,
          f"GAPIT.Association.GWAS_Results.{model}.{trait}.csv")
        in_se = os.path.join(
          raw_dir,
          f"GAPIT.Association.GWAS_StdErr.{model}.{trait}.csv")
        if not (os.path.exists(in_fn) and os.path.exists(in_se)):
            os.remove(
              os.path.join(
                out_dir,
                f"{model}.{trait}.feather"))
            continue

        mplot_fn = os.path.join(
          out_dir,
          f"{model}.{trait}.mplot.png")
        qqplot_fn = os.path.join(
          out_dir,
          f"{model}.{trait}.qqplot.png")

        out_feather = os.path.join(
          out_dir,
          f"{out_root}.feather")
        out_tsv = os.path.join(
          out_dir,
          f"{out_root}.tsv")
        out_gz = os.path.join(
          out_dir,
          f"{out_root}.tsv.gz")

        res = pd.read_csv(in_fn)
        res = res[
          ["Chr",
           "Pos",
           "P.value",
           "Effect"]]

        # Get std error if available
        se = pd.read_csv(in_se)["std Error"]
        if not se.isna().all():
            res["std Error"] = se

        res.to_csv(
          out_tsv,
          sep   = "\t",
          index = False)

        subprocess.run(
          ["bgzip",
           "-f",
           "-c",
           out_tsv],
          stdout = open(out_gz, "w"),
          check  = True)

        subprocess.run(
          ["tabix",
           "-s", "1",
           "-b", "2",
           "-e", "2",
           "--skip-lines", "1",
           "-f",
           out_gz],
          check=True)

        res["P.value"] = -np.log10(res["P.value"])
        res.sort_values(
          by="P.value",
          ascending=False).head(5000).to_feather(out_feather)

        mplot_fig = mplot(
          genome,
          res,
          chr_vis)
        mplot_fig.savefig(
          mplot_fn,
          dpi    = 300,
          format = "png")

        qqplot_fig = qqplot(res)
        qqplot_fig.savefig(
          qqplot_fn,
          dpi    = 300,
          format = "png")            
        # Clean up
        os.remove(out_tsv)
        for fn in (
          "GAPIT.Association.PVE",
          "GAPIT.Association.GWAS_StdErr",
          "GAPIT.Association.GWAS_Results",
          "GAPIT.Association.Vairance_markers",
          "GAPIT.Association.Prediction_results"):
            fn_path = os.path.join(
              raw_dir,
              f"{fn}.{model}.{trait}.csv")
            if os.path.exists(fn_path):
                os.remove(fn_path)
        for fn in (
          "GAPIT.Association.Optimum",
          "GAPIT.Association.Manhattan_Chro",
          "GAPIT.Association.Manhattan_Geno",
          "GAPIT.Association.Significant_SNPs",
          "GAPIT.Phenotype.Distribution_Significantmarkers"):
            fn_path = os.path.join(
              raw_dir,
              f"{fn}.{model}.{trait}.pdf")
            if os.path.exists(fn_path):
                os.remove(fn_path)

    if pr.returncode != 0:
        raise OSError(err.decode())
    return "Done."