#!/usr/bin/env python

import io
import os
import re
import gzip
import json
import pysam
import random
import shutil
import pgenlib
import bed_reader
import subprocess
import numpy as np
import pandas as pd
import statsmodels.api as sm
from collections import Counter
from glimix_core.lmm import LMM
from collections import defaultdict
from numpy_sugar.linalg import economic_qs_linear

from flask_app import (
  get_lock,
  running_pr,
  master_lock,
  task_status)

from app.dirs import (
  tmp_dir,
  data_dir,
  root_dir,
  tool_dir,
  plink1_path,
  plink2_path)

from modules.load_input import (
  load_gdata,
  load_config)

from modules.wpgma import wpgma
from modules.assoc_plot import mqplot
from modules.block import (
  win_block,
  block_genotype)

from app.task_worker import progress

from modules.forms import models as models_dict
models_dict = dict(models_dict)

# Bin size for making marker density histograms
mdbin = 1000000

###############################################################################
### Calculate marker set summary statistics
###############################################################################

def run_qc(genome, task_id, markers=None):
    config = load_config(genome)
    gdata = load_gdata(genome)

    if not markers:
        in_path = os.path.join(
          data_dir,
          genome,
          config["vcf"])
        stats_dir = os.path.join(
          data_dir,
          genome,
          "stats")
        gt_out = os.path.join(
          stats_dir,
          "gt")
        stats_out = gt_out
    else:
        stats_dir = os.path.join(
          data_dir,
          genome,
          "assoc",
          markers,
          "stats")
        in_path = os.path.join(
          data_dir,
          genome,
          "assoc",
          markers,
          "gt",
          "gt")
        gt_out = in_path
        stats_out = os.path.join(
          data_dir,
          genome,
          "assoc",
          markers,
          "stats",
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
      f"{stats_out}.{ext}")
      for ext in (
        "smiss",
        "vmiss",
        "het",
        "gcount",
        "afreq",
        "hardy"))):
        progress(task_id, "Calculating genotype statistics")
        if markers:
            popen(
              task_id,
              [plink2_path,
               "--pfile", in_path, "zst",
               "--missing", "scols=nmiss,nobs", "vcols=nmiss,nobs",
               "--geno-counts", "cols=refalt1,nobs",
               "--het", "cols=het,nobs",
               "--freq", "cols=reffreq",
               "--hardy", "cols=p",
               "--allow-extra-chr",
               "--out", stats_out])
        else:
            popen(
              task_id,
              [plink2_path,
               "--vcf", in_path,
               "--missing", "scols=nmiss,nobs", "vcols=nmiss,nobs",
               "--geno-counts", "cols=refalt1,nobs",
               "--het", "cols=het,nobs",
               "--bad-freqs",
               "--freq", "cols=reffreq",
               "--hardy", "cols=p",
               "--make-pgen", "vzs",
               "pvar-cols=", "psam-cols=",
               #"--export", "vcf", "bgz",
               "--allow-extra-chr",
               "--set-all-var-ids", "@:#",
               "--max-alleles", "2",
               "--polyploid-mode", "missing",
               "--vcf-half-call", "missing",
               "--out", stats_out])

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
      "smiss": [1,2],
      "vmiss": [1,2],
      "shet": [0,1,3],
      "vhet": [0,1,2],
      "maf": [1],
      "hwe": [2]}

    progress(task_id, "Processing genotype statistics")
    for metric, ext, id_col in zip(metrics, exts, id_cols):
        qc_in = f"{stats_out}.{ext}"
        cnt = Counter()
        df = pd.read_csv(
          qc_in,
          sep="\t",
          skiprows=1,
          header=None,
          index_col=None,
          chunksize=100000,
          usecols=usecols[metric])
        het_out = os.path.join(
          stats_dir,
          f"{metric}.txt.gz")
        if os.path.exists(het_out):
            os.remove(het_out)

        for chk in df:
            chk.columns = pd.RangeIndex(chk.columns.size)
            if metric in ("smiss", "vmiss"):
                cnt.update(chk[0] / chk[1])
            elif metric == "shet":
                chk[3] = chk[1] / chk[2]
                chk = chk.loc[:, [0,3]]
                cnt.update(chk[3])
            elif metric == "vhet":
                chk[3] = chk[1] / chk[2]
                chk = chk.loc[:, [0,3]]
                cnt.update(chk[3])
            elif metric == "maf":
                chk[1] = 1 - chk[0]
                cnt.update(chk.min(axis=1))
            elif metric == "hwe":
                cnt.update(-np.log10(chk[0]))
            if metric in ("shet", "vhet"):
                chk.to_csv(
                  het_out,
                  mode="a",
                  sep="\t",
                  header=False,
                  index=False,
                  compression="gzip")
        cnt = {(-1 if np.isnan(v) else v): cnt[v] for v in cnt.keys()}
        cnt_vals = list(cnt.keys())
        if metric == "hwe":
            bins = list(range(
              int(min(min(cnt), 0)),
              int(max(cnt_vals))+2))
        else:
            bins = [v/100 for v in range(
              int(min(min(cnt), 0)),
              int(max(cnt_vals)*100)+2)]
        weights = [cnt[v] for v in cnt_vals]
        hgts, bins = np.histogram(
          cnt_vals,
          weights=weights,
          bins=bins)
        hgts = hgts.tolist()
        bins = bins.tolist()
        with open(os.path.join(stats_dir, metric), "w") as hist_out:
            json.dump([hgts, bins], hist_out)

        if metric == "smiss":
            samp_num = sum(weights)
        elif metric == "vmiss":
            var_num = sum(weights)

    # Keep full data for sample missingness and heterozygosity
    if not markers:
        samp_path = os.path.join(
          data_dir,
          genome,
          config["samp"])
        with master_lock:
            lock = get_lock(samp_path)
        with lock:
            samp_df = pd.read_csv(
              samp_path,
              sep="\t",
              header=0,
              index_col=0)
            for metric, ext in (("smiss", "smiss"), ("shet", "het")):
                qc_in = f"{stats_out}.{ext}"
                df = (pd
                  .read_csv(
                    qc_in,
                    sep="\t",
                    index_col=0,
                    usecols=[0] + usecols[metric])
                  .rename_axis("Acc"))
                if metric == "smiss":
                    samp_df = pd.concat([
                      samp_df,
                      ((df["MISSING_CT"] / df["OBS_CT"])
                        .rename("Missingness")
                        .round(5))],
                      axis=1)
                elif metric == "shet":
                    samp_df = pd.concat([
                      samp_df,
                      ((df["O(HET)"] / df["OBS_CT"])
                        .rename("Heterozygosity")
                        .round(5))],
                      axis=1)
            samp_df.to_csv(
              samp_path,
              sep="\t")

    # LD decay scatter plot
    progress(task_id, "Calculating LD per distance")
    max_dist = str(max(gdata["chr_sizes"]) // 2000)
    popen(
      task_id,
      [plink2_path,
       "--pfile", gt_out, "vzs",
       "--thin-count", "50000",
       "--r2-unphased", "zs", "cols=pos",
       "--ld-window-kb", max_dist,
       "--ld-window-r2", "0",
       "--allow-extra-chr",
       "--out", stats_out])
    ldd_fn = os.path.join(
      stats_dir,
      "ldd")
    bin_sums = defaultdict(float)
    bin_counts = defaultdict(int)

    for chunk in pd.read_csv(
      f"{stats_out}.vcor.zst",
      sep="\t",
      chunksize=1000000):
        chunk["dist"] = abs(chunk["#POS_A"] - chunk["POS_B"])
        chunk["bin"] = chunk["dist"] // 10000 * 10000 + 10000
        for bin_, group in chunk.groupby("bin"):
            bin_sums[bin_] += group["UNPHASED_R2"].sum()
            bin_counts[bin_] += group["UNPHASED_R2"].count()

    bins = sorted(bin_sums.keys())
    mean_r2 = [round(bin_sums[b] / bin_counts[b], 3) for b in bins]
    (pd
      .DataFrame({"bin": bins, "R2": mean_r2})
      .to_csv(
        ldd_fn,
        sep="\t"))

    # Histogram of marker density
    progress(task_id, "Calculating marker density")
    concat = json.loads(config["concat"])
    ch_names = gdata["chr_names"]
    ch_sizes = dict(zip(ch_names, gdata["chr_sizes"]))
    ch_offsets = dict(zip(ch_names, gdata["chr_offsets"]))
    other_offset = sum(ch_sizes[c] for c in ch_names if c not in concat)

    bin_offsets = dict()
    ch_names = [ch for ch in ch_names if ch not in concat]
    other_size = sum(ch_sizes.pop(ch) for ch in concat)

    for idx, ch in enumerate(ch_names):
        bin_offsets[ch] = sum(
          (ch_sizes[c] // mdbin * mdbin + mdbin)
          for c in ch_names[:idx])
    if concat:
        bin_offsets["Other"] = sum(
          (ch_sizes[c] // mdbin * mdbin + mdbin)
          for c in ch_names if c not in concat)
        ch_sizes["Other"] = other_size
        ch_names.append("Other")

    hgts, bins, coords = list(), list(), list()
    cnts = Counter()

    df = pd.read_csv(
      f"{gt_out}.pvar.zst",
      sep="\t",
      chunksize=100000,
      header=None,
      skiprows=1,
      index_col=None,
      usecols=[0, 1])

    for chk in df:
        for ch, pos in chk.groupby(0):
            ch = str(ch)
            if ch in concat:
                pos[1] += ch_offsets[ch] - other_offset
                ch = "Other"
            cnts.update(
              (pos[1] + bin_offsets[ch])
              // mdbin
              * mdbin
              + mdbin)

    for ch in ch_names:
        size = ch_sizes[ch] // mdbin * mdbin + mdbin
        for block in range(
          bin_offsets[ch] + mdbin,
          size + bin_offsets[ch] + mdbin,
          mdbin):
            hgts.append(cnts.get(block, 0))
            bins.append(block - mdbin // 2)
            ch_block = block - bin_offsets[ch]
            coords.append(f"{ch}:{ch_block - mdbin:,}-{ch_block - 1:,}")

    with open(os.path.join(stats_dir, "dens"), "w") as dens_out:
        json.dump([hgts, bins, coords], dens_out)

    with open(os.path.join(stats_dir, "cnts"), "w") as cnts_out:
        json.dump([samp_num, var_num], cnts_out)

    # LD matrices per chromosome
    progress(task_id, "Calculating linkage disequilibrium")
    calc_chr_ld(genome, task_id, markers=markers)

    if not markers:
        var_only_root = os.path.join(
          stats_dir,
          "var_only")
        # Create a vars-only vcf file
        if not os.path.exists(f"{var_only_root}.vcf.gz"):
            popen(
              task_id,
              [plink2_path,
               "--pfile", gt_out, "vzs",
               "--make-just-pvar",
               "cols=vcfheader,qual,filter,info",
               "--allow-extra-chr",
               "--out", var_only_root])
            os.rename(f"{var_only_root}.pvar", f"{var_only_root}.vcf")
            pysam.tabix_compress(
              f"{var_only_root}.vcf",
              f"{var_only_root}.vcf.gz")
            pysam.tabix_index(
              f"{var_only_root}.vcf.gz",
              preset="vcf")

            os.remove(f"{var_only_root}.vcf")
            os.remove(f"{var_only_root}.log")

        cnts_fn = os.path.join(
          stats_dir,
          "cnts")
        # Save number of markers and samples
        pgen = pgenlib.PgenReader(f"{gt_out}.pgen".encode("utf-8"))
        var_all, samp_all = pgen.get_variant_ct(), pgen.get_raw_sample_ct()
        with open(cnts_fn, "w") as cnts_out:
            json.dump(
              [var_all, samp_all],
              cnts_out)

    # Clean up
    for fn in os.listdir(stats_dir):
        if (fn.startswith("gt")
          and not fn.endswith(
            (".gz",
             ".pgen",
             ".pvar.zst",
             ".psam"
             ))):
            os.remove(
              os.path.join(
                stats_dir,
                fn))

def calc_chr_ld(
  genome,
  task_id,
  markers=None):
    config = load_config(genome)
    gdata = load_gdata(genome)

    if not markers:
        stats_dir = os.path.join(
          data_dir,
          genome,
          "stats")
        in_path = os.path.join(
          stats_dir,
          "gt")
        gt_out = stats_out = in_path
    else:
        stats_dir = os.path.join(
          data_dir,
          genome,
          "assoc",
          markers,
          "stats")
        in_path = os.path.join(
          data_dir,
          genome,
          "assoc",
          markers,
          "gt",
          "gt")
        gt_out = in_path
        stats_out = os.path.join(
          data_dir,
          genome,
          "assoc",
          markers,
          "stats",
          "gt")

    concat = json.loads(config["concat"])
    for ch, span in zip(gdata["chr_names"], gdata["chr_sizes"]):
        if ch in concat:
            continue
        popen(
          task_id,
          [plink2_path,
           "--pfile", gt_out, "vzs",
           "--chr", str(ch),
           "--bp-space", str(span // 2000),
           "--r2-unphased", "square0",
           "--out", stats_out])

        r2_vals = pd.read_csv(
          f"{stats_out}.unphased.vcor2",
          sep="\t",
          header=None,
          index_col=None).round(3).fillna("")
        r2_idx = pd.read_csv(
          f"{stats_out}.unphased.vcor2.vars",
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
    os.remove(f"{stats_out}.unphased.vcor2")
    os.remove(f"{stats_out}.unphased.vcor2.vars")

###############################################################################
### Apply filters to marker set
###############################################################################

def filt_qc(
  genome,
  data,
  task_id):
    gdata = load_gdata(genome)
    config = load_config(genome)
    markers = data["name"]

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

    if data["filt_from"]:
        plink_in = os.path.join(
          data_dir,
          genome,
          "assoc",
          data["filt_from"],
          "gt",
          "gt")
        stats_dir = os.path.join(
          data_dir,
          genome,
          "assoc",
          data["filt_from"],
          "stats")

    if data["sample_group"]:
        samp_inc = os.path.join(
          filt_dir,
          "samp_include.txt")
        with open(samp_inc, "w") as samp_out:
            for group in data["sample_group"]:
                _ = samp_out.write("\n".join(json.loads(group)) + "\n")

    if data["coord"]:
        coord = []
        for group in data["coord"]:
            coord += json.loads(group)

        coord_inc = f"{plink_out}.bed1"
        with open(coord_inc, "w") as coord_out:
            for locus in coord:
                chrom, start = locus.split(":")
                start = start.replace(",", "")
                if "-" in start:
                    start, end = start.split("-")
                else:
                    end = start
                _ = coord_out.write("\t".join((chrom, start, end)) + "\n")

    # Prepare heterozygosity filters
    if data["shet"]:
        shet_exc = os.path.join(
          filt_dir,
          "shet_exclude.txt")
        if os.path.exists(shet_exc):
            os.remove(shet_exc)
        df = pd.read_csv(
          os.path.join(stats_dir, "shet.txt.gz"),
          sep="\t",
          header=None,
          index_col=None,
          compression="gzip",
          chunksize=100000)
        for chk in df:
            chk.columns = pd.RangeIndex(chk.columns.size)
            chk[chk[1] > float(data["shet"])][0].to_csv(
              shet_exc,
              mode="a",
              index=False,
              header=False)

    if data["vhet"]:
        vhet_exc = os.path.join(
          filt_dir,
          "vhet_exclude.txt")
        if os.path.exists(vhet_exc):
            os.remove(vhet_exc)
        df = pd.read_csv(
          os.path.join(stats_dir, "vhet.txt.gz"),
          sep="\t",
          header=None,
          index_col=None,
          compression="gzip",
          chunksize=100000)
        for chk in df:
            chk.columns = pd.RangeIndex(chk.columns.size)
            chk[chk[1] > float(data["vhet"])][0].to_csv(
              vhet_exc,
              mode="a",
              index=False,
              header=False)

    filt_out = plink_out + ".filt"
    # Apply filters
    if any(data[f] for f in (
      "sample_group", "coord", "smiss",
      "vmiss", "shet", "vhet", "maf", "hwe")):
        cmd = [plink2_path, "--pfile", plink_in, "vzs"]
        cmd = cmd + ["--remove", shet_exc] if data["shet"] else cmd
        cmd = cmd + ["--keep", samp_inc] if data["sample_group"] else cmd
        cmd = cmd + ["--exclude", vhet_exc] if data["vhet"] else cmd
        cmd = cmd + ["--extract", "bed1", coord_inc] if data["coord"] else cmd
        cmd = cmd + ["--mind", str(data["smiss"])] if data["smiss"] else cmd
        cmd = cmd + ["--geno", str(data["vmiss"])] if data["vmiss"] else cmd
        cmd = cmd + ["--maf", str(data["maf"])] if data["maf"] else cmd
        cmd = cmd + ["--hwe", str(10**(-data["hwe"]))] if data["hwe"] else cmd
        cmd = cmd + (["--make-just-pvar", "zs", "cols="] if data["imp"] else ["--make-pgen", "vzs", "pvar-cols=", "psam-cols="])
        cmd = cmd + ["--export", "vcf", "bgz", "id-paste=iid"] if data["imp"] else cmd
        cmd = cmd + ["--allow-extra-chr", "--out", filt_out]
        progress(task_id, "Applying filters")
        popen(task_id, cmd)

    else:
        filt_out = plink_in

    # Impute missing genotypes
    if data["imp"]:
        imp_out = os.path.join(
          filt_dir,
          "gt.imp")
        imp_vcf = f"{imp_out}.vcf.gz"
        vcf_filt = filt_out + ".vcf.gz"

        # Determine one-marker chromosomes (removed, cannot be imputed)
        chr_d = dict()
        chr_data = pd.read_csv(
          filt_out + ".pvar.zst",
          sep="\t",
          usecols=[0,1],
          skiprows=1,
          index_col=None,
          chunksize=1000000)
        for chk in chr_data:
            for ch, vs in chk.groupby(0):
                if vs.shape[0] == 1:
                    chr_d.setdefault(ch, []).append(vs.reset_index().loc[0,1])
        if any (len(v) == 1 for v in chr_d.values()):
            chrs = []
            excl_path = os.path.join(filt_dir, "exclude.txt")
            with open(excl_path, "w") as excl:
                for v in chr_d.values():
                    if len(v) == 1:
                        chrs.append(v[0].split(":")[0])
                        _ = excl.write(v[0] + "\n")
            progress(
              task_id, 
              "WARNING: these chromosomes have only a single marker and were removed:\n" +\
              ", ".join(chrs))

        beagle_path = os.path.join(
          tool_dir,
          "beagle.jar")

        progress(task_id, "Imputing missing genotypes")
        cmd = [
          "java",
          "-jar",
          beagle_path,
          f"gt={vcf_filt}",
          f"out={imp_out}"]
        if any (len(v) == 1 for v in chr_d.values()):
            cmd.append(f"excludemarkers={excl_path}")
        popen(
          task_id,
          cmd)

        popen(
          task_id,
          [plink2_path,
           "--vcf", imp_vcf,
           "--make-pgen", "vzs", "pvar-cols=", "psam-cols=",
           "--allow-extra-chr",
           "--out", imp_out])
    else:
        imp_out = filt_out

    # Prune by LD
    if data["prune"]:
        progress(task_id, "Calculating LD for pruning")
        data["window"] = (
          str(data["window"]) + "kb"
          if str(data["window"]) else "200kb")
        popen(task_id,
          [plink2_path,
          "--pfile", imp_out, "vzs",
          "--indep-pairwise", data["window"], str(data["prune"]),
          "--allow-extra-chr",
          "--out", plink_out])
        prune_in = plink_out + ".prune.in"

        cmd = [
          plink2_path,
          "--pfile", imp_out, "vzs",
          "--extract", prune_in,
          "--allow-extra-chr",
          "--out", plink_out]
        cmd = cmd + (["--export", "vcf", "bgz", "id-paste=iid"] if data["hap"] and data["btype"] == "ld" else ["--make-pgen", "vzs", "pvar-cols=", "psam-cols="])
        progress(task_id, "Pruning markers by LD")
        popen(task_id, cmd)


    gt_stats_dir = os.path.join(
      data_dir,
      genome,
      "assoc",
      markers,
      "stats")
    if not os.path.exists(gt_stats_dir):
        os.makedirs(gt_stats_dir)

    if data["hap"]:
        progress(task_id, "Constructing haplotype blocks")
        # Make haplotype blocks
        if data["btype"] == "ld":
            popen(
              task_id,
              [plink1_path,
               "--vcf", f"{plink_out}.vcf.gz",
               "--allow-extra-chr",
               "--blocks", "no-pheno-req",
               "--blocks-max-kb", str(data["bsize"]),
               "--out", plink_out])
        else:
            win_block(
              plink_out,
              data["btype"],
              data["bspan"],
              data["bstep"],
              task_id)

        #os.remove(full_vcf)
        #os.remove(plink_out + ".vcf.gz")
        # Genotype based on blocks
        maf = 0.05
        min_haps = 0
        block_genotype(
          plink_out,
          gt_stats_dir,
          maf,
          min_haps,
          json.loads(config["concat"]),
          task_id)

        kin_in =  plink_out + ".allhaps"
        # Remove sequence length information from .bim file
        os.rename(plink_out + ".bim", plink_out + ".lengths.bim")
        (pd
          .read_csv(
            plink_out + ".lengths.bim",
            sep="\t",
            header=None,
            index_col=None,
            usecols=[0,1,2,3])
          .assign(A1 = "A")
          .assign(A2 = "T")
          .to_csv(
            plink_out + ".bim",
            sep="\t",
            header=False,
            index=False))

    else:
        progress(task_id, "Processing filtered data")
        # Make necessary genotype format files
        popen(
          task_id,
          [plink2_path,
           "--pfile", plink_out, "vzs",
           "--make-bed",
           "--allow-extra-chr",
           "--out", plink_out])
        kin_in = plink_out

    # Calculate K and Q
    progress(task_id, "Determining kinship and population structure")
    popen(
      task_id,
      [plink2_path,
       "--bfile", kin_in,
       "--make-rel", "square",
       "--pca", "15",
       "--allow-extra-chr",
       "--out", plink_out])
    rel_matrix = pd.read_csv(
      f"{plink_out}.rel",
      sep="\t",
      header=None,
      index_col=None)

    # Convert PCA eigenvectors to PVE
    sum_diag = sum(np.diag(rel_matrix))
    eigenval = pd.read_csv(
      f"{plink_out}.eigenval",
      header=None,
      index_col=None)
    eigenval = (eigenval / sum_diag * 100)[0].tolist()
    eigenvec = pd.read_csv(
      f"{plink_out}.eigenvec",
      sep="\t",
      header=0,
      index_col=0,
      usecols=[1,2,3,4])
    samples = eigenvec.index.tolist()
    eigenvec = [eigenvec[c].values.tolist() for c in eigenvec.columns]
    pca_path = os.path.join(
      gt_stats_dir,
      "pca")
    with open(pca_path, "w") as pca_out:
        json.dump([samples, eigenvec, eigenval], pca_out)

    # Cluster kinship matrix
    hamm_out = os.path.join(
      gt_stats_dir,
      "kin")

    popen(
      task_id,
      [plink1_path,
       "--bfile", kin_in,
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
    kin_matrix = (pd
      .DataFrame(hamm_matrix)
      .set_axis(hamm_labels, axis=0)
      .set_axis(hamm_labels, axis=1))

    tree, label_order = wpgma(
      hamm_matrix,
      hamm_labels)

    kin_matrix = kin_matrix.loc[label_order, label_order]
    with open(f"{hamm_out}", "w") as matr_out:
        json.dump(
          [kin_matrix.values.tolist(),
           label_order],
          matr_out)

    progress(task_id, "Saving summary")
    # Create file describing the filtering parameters
    params_fn = os.path.join(
      gt_stats_dir,
      "params")

    # Get original number of markers and samples
    with open(os.path.join(stats_dir, "cnts"), "r") as cnts_in:
        samp_all, var_all = json.load(cnts_in)
    # Get filtered number of markers and samples
    bed = bed_reader.open_bed(f"{plink_out}.bed")
    samp_filt, var_filt = bed.iid_count, bed.sid_count

    with open(params_fn, "w") as params_out:
        params = [
          ["Sample missingness",
           "Variant missingness",
           "Sample heterozygosity",
           "Variant heterozygosity",
           "Variant MAF",
           "Variant -log10(HWE)",
           "Variant LD",
           "LD window",
           "Samples retained/original",
           "Markers retained/original",
           "Is a haplotype set"],
          [float(data["smiss"]) if data["smiss"] else "—",
           float(data["vmiss"]) if data["vmiss"] else "—",
           float(data["shet"]) if data["shet"] else "—",
           float(data["vhet"]) if data["vhet"] else "—",
           float(data["maf"]) if data["maf"] else "—",
           float(data["hwe"]) if data["hwe"] else "—",
           float(data["prune"]) if data["prune"] else "—",
           data["window"] if data["window"] else "—",
           f"{samp_filt}/{samp_all}",
           f"{var_filt}/{var_all}",
           "Yes" if data["hap"] else "No"]]
        if data["hap"]:
            params[0].append("Block type")
            if data["btype"] == "ld":
                params[0].append("Max size (kb)")
                params[1].append("LD-based")
                params[1].append(data["bsize"])
            else:
                params[0].append("Block span")
                params[0].append("Block step")
                params[1].append("SNP number window" if data["btype"] == "snpwin" else "Base-pair window")
                params[1].append(data["bspan"])
                params[1].append(data["bstep"])
            
        json.dump(params, params_out)

    # Clean up
    for fn in (
      "kin.dist",
      "kin.dist.id",
      "kin.log",
      "kin.nosex"):
        fn_path = os.path.join(
          gt_stats_dir,
          fn)
        if os.path.exists(fn_path):
            os.remove(fn_path)

    # Clean up
    for fn in (
      "gt.blocks",
      "gt.blocks.det",
      "gt.eigenval",
      "gt.filt.log",
      "gt.filt.vcf.gz",
      #"gt.full.vcf.gz",
      "gt.lengths.bim",
      "gt.log",
      "gt.miss.log",
      "gt.miss.mindrem.id",
      "gt.miss.vcf.gz",
      "gt.nosex",
      "gt.pgen",
      "gt.prune.in",
      "gt.prune.out",
      "shet_exclude.txt",
      "vhet_exclude.txt"):
        fn_path = os.path.join(
          filt_dir,
          fn)
        if os.path.exists(fn_path):
            os.remove(fn_path)

#    if data["hap"]:
#        for ext in (
#          ".bed",
#          ".bim",
#          ".fam"):
#            fn_path = os.path.join(
#              filt_dir,
#              "gt.full" + ext)
#            if os.path.exists(fn_path):
#                os.remove(fn_path)
#            shutil.copy(
#              os.path.join(filt_dir, "gt" + ext),
#              fn_path)

        calc_chr_ld(
          genome,
          task_id,
          markers = data["name"])
    else:
        run_qc(
          genome,
          task_id,
          markers = data["name"])

    # Change chromosomes to numeric (required by BLINK)
    chr_rename = {
      ch: idx for idx, ch in enumerate(gdata["chr_names"], start=1)}
    (pd
      .read_csv(
        f"{plink_out}.bim",
        sep="\t",
        index_col=None,
        header=None)
      .set_axis(
        ["CHR", "ID", "FAM",
         "POS", "REF", "ALT"],
        axis=1)
      .assign(CHR=lambda x: x["CHR"].astype(str).map(chr_rename))
      .to_csv(
        f"{plink_out}.bim",
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

# Extract QTLs into Excel
def extract_qtls(genome, params, task_id):
    runs   = json.loads(params["qtl_runs"])
    window = int(params["window_ld"])
    thresh = float(params["plot_thr"])
    r2     = float(params["r2_ld"])
    _2r2   = (r2 + (1 - r2) / 2)

    qtls = [[
     "Marker set", "Model", "Trait", "Correction", "Peak marker",
     "Peak significance", f"Locus at LD ≥ {_2r2:.2f}", f"Locus at LD ≥ {r2}",
     "Minor allele frequency", "Effect size", "Percentage of variance explained"]]

    for run in runs:
        markers, model, trait = run

        progress(task_id, ", ".join((markers, models_dict[model], trait)))

        params_fn = os.path.join(
          data_dir,
          genome,
          "assoc",
          markers,
          "stats",
          "params")
        with open(params_fn) as params_in:
            is_hap = json.load(params_in)[1][10] == "Yes"
        res_fn = os.path.join(
          data_dir,
          genome,
          "assoc",
          markers,
          "res",
          f"{model}.{trait}.feather")
        res = pd.read_feather(res_fn)
        sig_col, ch_col, pos_col = res.columns[:3]
        if is_hap:
            end_col, idx_col = res.columns[3:5]
        else:
            end_col, idx_col = None, None
        res = res[res[sig_col] >= thresh]
        if is_hap:
            res["id"] = (
              res[ch_col].astype(str) + ":" +
              res[pos_col].astype(str) + "-" +
              res[end_col].astype(str) + "-" +
              res[idx_col].astype(str))
        else:
            res["id"] = (
              res[ch_col].astype(str) + ":" +
              res[pos_col].astype(str))

        for ch, data in res.groupby(ch_col):
            while not data.empty:
                peak = data.iloc[0]
                if is_hap:
                    peak_pos = f"{ch}:{peak[pos_col]:,}-{peak[end_col]:,}_{peak[idx_col]:,}"
                else:
                    peak_pos = f"{ch}-{peak[pos_col]:,}"
                data = data.iloc[1:]

                                                                                                        ##################### CORRECT IF NEEDED
                in_root = os.path.join(
                  data_dir,
                  genome,
                  "assoc",
                  markers,
                  "gt",
                  "gt.orig")
                out_root = os.path.join(
                  tmp_dir,
                  task_id)

                if is_hap:
                    start, end = peak[pos_col], peak[end_col]
                else:
                    start = end = peak[pos_col]
                popen(
                  task_id,
                  [plink2_path,
                   "--bfile", in_root,
                   "--chr", ch,
                   "--from-bp", str(start - window*1000),
                   "--to-bp", str(end + window*1000),
                   "--r2-unphased", "cols=id,pos",
                   "--ld-snp", peak["id"],
                   "--ld-window-kb", str(window),
                   "--ld-window-r2", str(r2),
                   #"--ld-window-r2", "0",
                   #"--silent",                                                                      ############## UNSILENT
                   "--out", out_root])

                ld = pd.read_csv(
                  f"{out_root}.vcor",
                  sep="\t",
                  header=0,
                  index_col=None,
                  usecols=[2,3,4])

                # QTL span calculation
                if is_hap:
                    if not ld.empty:
                        span_r2 = (
                          ch + ":" + 
                          f"{ld.iloc[0]['POS_B']:,}" + "-" +
                          f"{ld.iloc[-1]['ID_B'].str.split('-')[1].split('_')[0]:,}")
                    else:
                        span_r2 = peak_pos

                    ld_2r2 = ld[ld["UNPHASED_R2"] >= _2r2]
                    if not ld_2r2.empty:
                        span_2r2 = (
                          ch + ":" + 
                          f"{ld_2r2.iloc[0]['POS_B']:,}" + "-" +
                          f"{ld_2r2.iloc[-1]['ID_B'].str.split('-')[1].split('_')[0]:,}")
                    else:
                        span_2r2 = peak_pos
                    
                else:
                    if not ld.empty:
                        span_r2 = (
                          ch + ":" + 
                          f"{ld.iloc[0]['POS_B']:,}" + "-" +
                          f"{ld.iloc[-1]['POS_B']:,}")
                    else:
                        span_r2 = peak_pos

                    ld_2r2 = ld[ld["UNPHASED_R2"] >= _2r2]
                    if not ld_2r2.empty:
                        span_2r2 = (
                          ch + ":" + 
                          f"{ld_2r2.iloc[0]['POS_B']:,}" + "-" +
                          f"{ld_2r2.iloc[-1]['POS_B']:,}")
                    else:
                        span_2r2 = peak_pos

                    if model == "mtmlm":
                        effect = "-"
                        pve = "-"
                    else:
                        # PVE and effect calculation
                        bed = bed_reader.open_bed(f"{in_root}.bed")
                        config = load_config(genome)
                        phe_path = os.path.join(
                          data_dir,
                          genome,
                          config["phe"])
                        with master_lock:
                            lock = get_lock(phe_path)
                        with lock:
                            phe = (pd
                              .read_csv(
                                phe_path,
                                sep="\t",
                                header=0,
                                index_col=0)
                              .loc[bed.iid, trait]
                              .dropna())

                        marker_idx = [
                          idx for idx, marker in enumerate(open(f"{in_root}.bim"))
                          if marker and marker.split("\t", 2)[1] == peak["id"]]

                        samples = phe.index
                        sample_idxs = [
                          idx for idx, sample in enumerate(open(f"{in_root}.fam"))
                          if sample and sample.split("\t", 2)[1] in samples]
                        if len(marker_idx) != 1 or len(sample_idxs) != samples.size:
                            pve = "-"
                        else:
                            gts = (bed
                              .read(
                                dtype="int8",
                                index=marker_idx)
                              .flatten())

                            filt_gts = gts[~np.isnan(gts)]
                            af = filt_gts.sum() / (filt_gts.size * 2)
                            maf = round(min(af, 1 - af), 2)
                            gts = gts[sample_idxs]

                            params_fn = os.path.join(
                              data_dir,
                              genome,
                              "assoc",
                              markers,
                              "res",
                              f"{model}.{trait}.params")
                            with open(params_fn) as params_in:
                                q_num, k, kq_from, cov = json.load(params_in)[1][2:6]
                                k = k == "Yes"
                            corr_string = ""
                            if int(q_num):
                                corr_string += f"{q_num} PCs, "
                            if k:
                                corr_string += "K, "
                            if k or int(q_num):
                                corr_string += f"from {kq_from}"
                            if cov:
                                if k or int(k_num):
                                    corr_string += ", + "
                                corr_string += "covariates"

                            phe_corr = kq_correct(
                              genome,
                              phe,
                              samples,
                              k,
                              int(q_num),
                              markers,
                              cov)
                            phe_fit = sm.OLS(phe_corr, sm.add_constant(gts)).fit()
                            effect = round(phe_fit.params["x1"].astype(float), 2)
                            pve = round(phe_fit.rsquared * 100, 2)

                qtls.append([
                  markers,
                  models_dict[model],
                  trait,
                  corr_string,
                  peak_pos,
                  round(peak[sig_col], 2),
                  span_2r2,
                  span_r2,
                  maf,
                  effect,
                  pve,
                  ])

                data = data[~data["id"].isin(ld["ID_B"])]

    #print(*qtls, sep="\n")
    qtl_df = pd.DataFrame(qtls[1:], columns=qtls[0])
    #excel_out = io.BytesIO()
    excel_out = os.path.join(tmp_dir, "QTL_data.xlsx")
    with pd.ExcelWriter(excel_out, engine="xlsxwriter") as writer:
        qtl_df.to_excel(writer, index=False, sheet_name="Significant QTLs")
    #excel_out.seek(0)

