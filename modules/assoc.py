#!/usr/bin/env python

import io
import os
import re
import pty
import gzip
import json
import pysam
import shutil
import bed_reader
import subprocess
import numpy as np
import configparser
import pandas as pd
import statsmodels.api as sm
from collections import Counter
from glimix_core.lmm import LMM
import matplotlib.pyplot as plt
from numpy_sugar.linalg import economic_qs_linear

from flask import (
  url_for,
  redirect)

from flask_app import (
  get_lock,
  master_lock,
  task_status)

from app.dirs import (
  tmp_dir,
  data_dir,
  root_dir,
  tool_dir,
  plink1_path,
  plink2_path)

from modules.rf import (
  run_rf,
  optim_rf)

from modules.load_input import (
  load_gdata,
  load_config)

from modules.wpgma import wpgma
from modules.assoc_plot import mqplot
from modules.block import (
  win_block,
  block_genotype)

from app.task_worker import progress

# Bin size for making marker density histograms
mdbin = 1000000

###############################################################################
### Open a command-line process
###############################################################################

def popen(task_id, cmd, show_stdout=True):
    master_fd, slave_fd = pty.openpty()
    pr = subprocess.Popen(
      cmd,
      stdout=slave_fd,
      stderr=subprocess.STDOUT,
      close_fds=True)
    os.close(slave_fd)

    if show_stdout:
        with master_lock:
            lock = get_lock("__task_status__")
        buffer = b""
        with os.fdopen(master_fd, "rb") as master:
            while True:
                try:
                    char = master.read(1)
                except OSError as e:
                    if e.errno == 5:
                        break
                    else:
                        raise
                else:
                    if not char:
                        break
                    if char in (b"\r", b"\n", b"\x08"):
                        if char == b"\r":
                            next_char = master.read(1)
                            if next_char and next_char == b"\n":
                                char = next_char
                            else:
                                char += next_char
                        with lock:
                            task_status[task_id]["progress"] += buffer.decode(
                              "utf-8",
                              errors="replace")
                        buffer = b""
                    buffer += char
        if buffer:
            with lock:
                task_status[task_id]["progress"] += buffer.decode("utf-8", errors="replace")

    pr.wait()
    if show_stdout:
        with lock:
            task_status[task_id]["progress"] += "\n"

    if (
      (cmd[0] != os.path.join(tool_dir, "pimass", "pimass-lin")
       and pr.returncode != 0) or
      (cmd[0] == os.path.join(tool_dir, "pimass", "pimass-lin")
       and pr.returncode == 0)):
        raise subprocess.CalledProcessError(1, "")
    return pr.returncode

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
               "--bfile", in_path,
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
               "--freq", "cols=reffreq",
               "--hardy", "cols=p",
               "--make-bed",
               "--export", "vcf", "bgz",
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
      f"{gt_out}.bim",
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
    calc_chr_ld(genome, task_id, markers=markers)

    if not markers:
        var_only_root = os.path.join(
          data_dir,
          genome,
          "var_only")
        # Create a vars-only vcf file
        if not os.path.exists(f"{var_only_root}.vcf.gz"):
            popen(
              task_id,
              [plink2_path,
               "--bfile", gt_out,
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
        bed = bed_reader.open_bed(f"{gt_out}.bed")
        var_all, samp_all = bed.iid_count, bed.sid_count
        with open(cnts_fn, "w") as cnts_out:
            json.dump(
              [var_all, samp_all],
              cnts_out)

    # Clean up
    for fn in os.listdir(stats_dir):
        if (fn.startswith("gt")
          and not fn.endswith(
            (".gz",
             ".bed",
             ".bim",
             ".fam"
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

    concat = json.loads(config["concat"])
    for ch, span in zip(gdata["chr_names"], gdata["chr_sizes"]):
        if ch in concat:
            continue
        popen(
          task_id,
          [plink2_path,
           "--bfile", gt_out,
           "--chr", str(ch),
           #"--thin-count", "1000",
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
    vcf_file = config["vcf"]
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

    if data["filt_ext"] and data["filt_from"]:
        plink_in = os.path.join(
          data_dir,
          genome,
          "assoc",
          data["filt_from"],
          "gt",
          "gt.full")
        stats_dir = os.path.join(
          data_dir,
          genome,
          "assoc",
          data["filt_from"],
          "stats")

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

    if data["samples_group"]:
        samp_inc = os.path.join(
          filt_dir,
          "samp_include.txt")
        with open(samp_inc, "w") as samp_out:
            for group in data["samples_group"]:
                _ = samp_out.write("\n".join(json.loads(group)) + "\n")

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

    plink_filt = plink_out + ".filt"
    # Apply filters
    if any(data[f] for f in (
      "samples_group", "smiss", "vmiss",
      "shet", "vhet", "maf", "hwe")):
        cmd = [plink2_path, "--vcf", f"{plink_in}.vcf.gz"]
        cmd = cmd + ["--remove", shet_exc] if data["shet"] else cmd
        cmd = cmd + ["--keep", samp_inc] if data["samples_group"] else cmd
        cmd = cmd + ["--exclude", vhet_exc] if data["vhet"] else cmd
        cmd = cmd + ["--mind", str(data["smiss"])] if data["smiss"] else cmd
        cmd = cmd + ["--geno", str(data["vmiss"])] if data["vmiss"] else cmd
        cmd = cmd + ["--maf", str(data["maf"])] if data["maf"] else cmd
        cmd = cmd + ["--hwe", str(10**(-data["hwe"]))] if data["hwe"] else cmd
        cmd = cmd + ["--make-bed"] if not data["imp"] else cmd
        cmd = cmd + ["--make-just-bim"] if data["imp"] else cmd
        cmd = cmd + [
          "--allow-extra-chr",
          "--export", "vcf", "bgz", "id-paste=iid",
          "--out", plink_filt]
        progress(task_id, "Applying filters")
        popen(task_id, cmd)

        if not data["imp"]:
            for ext in (".bed", ".bim", ".fam", ".vcf.gz"):
                os.rename(
                  plink_filt + ext,
                  plink_out + ".full" + ext)
            #plink_filt = plink_out + ".full"
    else:
        if data["imp"]:
            for ext in (".bed", ".bim", ".fam", ".vcf.gz"):
                if not os.path.exists(plink_filt + ext):
                    os.link(
                      plink_in + ext,
                      plink_filt + ext)
        else:
            plink_filt = plink_in
            for ext in (".bed", ".bim", ".fam", ".vcf.gz"):
                if not os.path.exists(plink_out + ".full" + ext):
                    os.link(
                      plink_filt + ext,
                      plink_out + ".full" + ext)

    full_root = os.path.join(
      filt_dir,
      "gt.full")
    full_vcf = f"{full_root}.vcf.gz"
    vcf_in = plink_filt + ".vcf.gz"

    # Impute missing genotypes
    if data["imp"]:
        # Determine one-marker chromosomes (removed, cannot be imputed)
        chr_d = dict()
        chr_data = pd.read_csv(
          plink_filt + ".bim",
          sep="\t",
          usecols=[0,1],
          header=None,
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
              "WARNING: these chromosomes have only a single marker:\n" +\
              ", ".join(chrs) + "\n" +\
              "and were removed.")

        beagle_path = os.path.join(
          tool_dir,
          "beagle.jar")

        progress(task_id, "Imputing missing genotypes")
        cmd = [
          "java",
          "-jar",
          beagle_path,
          f"gt={vcf_in}",
          f"out={full_root}"]
        if any (len(v) == 1 for v in chr_d.values()):
            cmd.append(f"excludemarkers={excl_path}")
        popen(
          task_id,
          cmd)

        popen(
          task_id,
          [plink2_path,
           "--vcf", full_vcf,
           "--make-bed",
           "--allow-extra-chr",
           "--out", full_root])

    # Prune by LD
    if data["prune"]:
        progress(task_id, "Calculating LD for pruning")
        data["window"] = (
          str(data["window"]) + "kb"
          if str(data["window"]) else "200kb")
        popen(task_id,
          [plink2_path,
          "--vcf", full_vcf,
          "--indep-pairwise", data["window"], str(data["prune"]),
          "--allow-extra-chr",
          "--out", plink_out])
        prune_in = plink_out + ".prune.in"

        cmd = [
          plink2_path,
          "--vcf", full_vcf,
          "--extract", prune_in,
          "--allow-extra-chr",
          "--export", "vcf", "bgz", "id-paste=iid",
          "--out", plink_out]
        progress(task_id, "Pruning markers by LD")
        popen(task_id, cmd)
    else:
        if not os.path.exists(plink_out + ".vcf.gz"):
            os.link(
              full_vcf,
              plink_out + ".vcf.gz")

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
        # Make necessary genotype format files
        popen(
          task_id,
          [plink2_path,
           "--vcf", f"{plink_out}.vcf.gz",
           "--allow-extra-chr",
           "--make-bpgen",
           "--out", plink_out])
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

        os.remove(full_vcf)
        os.remove(plink_out + ".vcf.gz")
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

        #pysam.tabix_index(
        #  output_file,
        #  preset="vcf")

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
           "--vcf", f"{plink_out}.vcf.gz",
           "--make-bed",
           "--allow-extra-chr",
           "--out", plink_out])
        kin_in = plink_out

    pysam.tabix_index(
      f"{plink_out}.vcf.gz",
      preset="vcf",
      force=True)

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

#    kin_matrix = (
#      rel_matrix
#        .set_axis(hamm_labels, axis=0)
#        .set_axis(hamm_labels, axis=1)
#        .loc[label_order, label_order])
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
      "gt.full.log",
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

    if data["hap"]:
        for ext in (
          ".bed",
          ".bim",
          ".fam",
          ".vcf.gz"):
            fn_path = os.path.join(
              filt_dir,
              "gt.full" + ext)
            if os.path.exists(fn_path):
                os.remove(fn_path)
            shutil.copy(
              os.path.join(filt_dir, "gt" + ext),
              fn_path)

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

def run_assoc(genome, params, task_id):
    config  = load_config(genome)
    gdata = load_gdata(genome)
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

    # Load samples and phenotype data
    fam_file = f"{in_root}.fam"
    fam_df = pd.read_csv(
      fam_file,
      sep="\t",
      header=None,
      index_col=None)
    samples = fam_df[1]

    phe_path = os.path.join(
      data_dir,
      genome,
      config["phe"])
    with master_lock:
        lock = get_lock(phe_path)
    with lock:
        phe_df = (pd
          .read_csv(
            phe_path,
            sep="\t",
            header=[0],
            index_col=[0])
         .loc[samples])

    if model == "blink":
        # BLINK analyses either one or all traits
        # phenotype file rewritten each time to keep only selected traits
        (phe_df
          .fillna("NaN")
          .loc[:, params["traits"]]
          .rename_axis("Acc")
          .to_csv(
            f"{in_root}.txt",
            sep="\t"))

        params["k"] = True
        q_matr, k_matr = kq_process(
          genome,
          markers,
          params["kq_from"] or markers,
          samples,
          False,
          int(params["q_num"]) if params["q"] else 0,
          params["cov"])

        if not q_matr.empty:
          q_matr.to_csv(
            f"{in_root}.cov",
            sep="\t")
        elif os.path.exists(f"{in_root}.cov"):
            os.remove(f"{in_root}.cov")

        blink_path = os.path.join(
          tool_dir,
          "BLINK",
          "blink_linux")
        progress(task_id, "Running analysis")
        command = [
          blink_path,
          "--gwas",
          "--file", in_root,
          "--plink",
          #"--bin_size", "5", "10", "2", "0.3",
          "--out", os.path.join(out_dir, "BLINK"),
          "--trait", "0"]
        popen(
          task_id,
          command)

        for trait in params["traits"]:
            assoc_process(
              genome,
              assoc_dir,
              f"BLINK_{trait}_GWAS_result.txt",
              f"{model}.{trait}",
              out_dir,
              ch_col="chr",
              pos_col="pos",
              sig_col="p_value")

            # Save parameters
            run_params = [
              ["Markers", "Samples", "PCs", "Kinship", "K/Q source"],
              [sum(1 for l in open(f"{in_root}.bim")),
               phe_df.loc[samples, trait].dropna().shape[0],
               params["q_num"] if params["q"] else 0,
               "Yes",
               params["kq_from"] or markers]]
            with open(
              os.path.join(
                out_dir,
                f"{model}.{trait}.params"),
              "w") as params_out:
                json.dump(run_params, params_out)

        os.remove(os.path.join(out_dir, "BLINK.log"))

    elif model == "rf":
        # Correct for kinship and population structure
        phe_df = phe_df.loc[samples, params["traits"]]

        if params["k"] and ((not params["q"] or params["q_num"] == 0) and not params["cov"]):
            raise OSError(
              "<b>Kiship correction must be used " +
              "together with covariate or PCA correction.<b>")

        q_matr, k_matr = kq_process(
          genome,
          markers,
          params["kq_from"] or markers,
          samples,
          params["k"],
          int(params["q_num"]) if params["q"] else 0,
          params["cov"])
        gt = bed_reader.open_bed(f"{in_root}.bed").read()

        if params["optimize"]:
            if len(params["traits"]) != 1 or len(params["traits_group"]) != 0:
                raise OSError("Single trait must be selected.")
            progress(task_id, "Running optimization")
        else:
            progress(task_id, "Running analysis")

        for trait in params["traits"]:
            progress(task_id, trait, bold=False)
            phe = (phe_df
              .loc[samples, trait]
              .dropna())
            phe_samples = phe.index
            phe_gt = gt[np.isin(samples, phe_samples), :]

            if not q_matr.empty:
                phe_samples = phe_samples.intersection(q_matr.index)
                phe = phe.loc[phe_samples]
                phe_q = q_matr.loc[phe_samples]
                phe_q = sm.add_constant(phe_q)
                if params["k"]:
                    phe_k = k_matr.loc[phe_samples, phe_samples]
                    QS = economic_qs_linear(phe_k)
                    lmm = LMM(
                      phe,
                      phe_q,
                      QS=QS)
                    lmm.fit()
                    phe_corr = phe - lmm.mean()
                else:
                    phe_corr = sm.OLS(phe, phe_q).fit().resid
            else:
                phe_corr = phe

            if params["optimize"]:
                mtryvals = [
                  int(v.strip())
                  for v in params["mtryvals"].split(", ")]
                ntreevals = [
                  int(v.strip())
                  for v in params["ntreevals"].split(", ")]
                res = optim_rf(
                  task_id,
                  phe_gt,
                  phe_corr,
                  mtryvals,
                  ntreevals)
                optim_path = os.path.join(
                  tmp_dir,
                  f"{task_id}.rfopt")
                with open(optim_path, "w") as optim_out:
                    json.dump([mtryvals, ntreevals, res], optim_out)
                return

            res = run_rf(
              task_id,
              phe_gt,
              phe_corr,
              params["mtry"],
              params["ntree"],
              params["niter"],
              params["nperm"])
            #res = pd.read_csv(f"data/{trait}_corr.txt", sep="\t", index_col=0)

            # Output files
            processed_out = os.path.join(out_dir, f"{model}.{trait}")
            mqplot_fn     = f"{processed_out}.mqplot.png"
            out_feather   = f"{processed_out}.feather"
            out_tsv       = f"{processed_out}.tsv"
            out_gz        = f"{processed_out}.tsv.gz"

            # Save raw results, rename numeric chromosomes back to original
            marker_map = (pd
              .read_csv(
                f"{in_root}.bim",
                sep="\t",
                header=None,
                index_col=None,
                usecols=[0,3])
              .set_axis(
                ["Chr", "Pos"],
                axis=1))
            chr_rename = dict(enumerate(gdata["chr_names"], start=1))
            res["Chr"]  = marker_map["Chr"].map(chr_rename)
            res["Pos"] = marker_map["Pos"]
            res = res.loc[:, ["Chr", "Pos", "Imp"]]
            res.to_csv(
              out_tsv,
              sep="\t",
              index=False)
            pysam.tabix_compress(
              out_tsv,
              out_gz,
              force=True)
            #pysam.tabix_index(
            #  out_gz,
            #  seq_col=0,
            #  start_col=1,
            #  end_col=1,
            #  line_skip=1, force=True)

            with open(os.path.join(assoc_dir, "stats", "params")) as params_in:
                is_hap = json.load(params_in)[1][10] == "Yes"
            if is_hap:
                blockmap = (pd
                  .read_csv(
                    os.path.join(assoc_dir, "gt", "gt.blockmap"),
                    sep="\t",
                    header=0,
                    index_col=None))
                res = (pd
                  .merge(
                    res,
                    blockmap,
                    left_on=["Chr", "Pos"],
                    right_on=["chr", "start"],
                    how="left")
                  .drop(
                    ["start"],
                    axis=1))
                res["idx"] = res.groupby(["Chr", "Pos"]).cumcount() + 1
                c = res.columns.tolist()
                (res[c[2:3] + c[:2] + c[-2:]]
                  .sort_values(
                    by="Imp",
                    ascending=False)
                  .head(5000)
                  .reset_index(drop=True)
                  .to_feather(out_feather))
            else:
                c = res.columns.tolist()
                (res[c[2:3] + c[:2]]
                  .sort_values(
                    by="Imp",
                    ascending=False)
                  .head(5000)
                  .reset_index(drop=True)
                  .to_feather(out_feather))

            mqplot_fig = mqplot(
              genome,
              res,
              mqplot_fn,
              json.loads(config["concat"]),
              ch_col="Chr",
              pos_col="Pos",
              sig_col="Imp",
              thresh=False,
              qqplot=False,
              end_col="end" if is_hap else None)

            os.remove(out_tsv)
            # Save parameters
            run_params = [
              ["Markers", "Samples", "PCs", "Kinship", "K/Q source",
               "Test markers", "Trees", "Iterations", "Permutations"],
              [sum(1 for l in open(f"{in_root}.bim")),
               phe_samples.shape[0],
               params["q_num"] if params["q"] else 0,
               "Yes" if params["k"] else "No",
               params["kq_from"] or markers,
               params["mtry"],
               params["ntree"],
               params["niter"],
               params["nperm"]]]
            with open(
              os.path.join(
                out_dir,
                f"{model}.{trait}.params"),
              "w") as params_out:
                json.dump(run_params, params_out)

    elif model in ("B", "C"):
        # GCTB requires space-separated .fam and cofactor/kinship files
        for ext in ("bed", "bim"):
            if not os.path.exists(f"{in_root}.gctb.{ext}"):
                os.link(f"{in_root}.{ext}", f"{in_root}.gctb.{ext}")

        (fam_df
          .assign(IID = samples.tolist())
          .loc[:, ["IID"] + fam_df.columns.tolist()[1:]]
          .to_csv(
            f"{in_root}.gctb.fam",
            sep=" ",
            header=False,
            index=False))
        (phe_df
          .assign(IID = phe_df.index.tolist())
          .loc[:, ["IID"] + params["traits"]]
          .to_csv(
            f"{in_root}.txt",
            header=False,
            sep=" "))

        q_matr, k_matr = kq_process(
          genome,
          markers,
          params["kq_from"] or markers,
          samples,
          params["k"],
          int(params["q_num"]) if params["q"] else 0,
          params["cov"])
        if not q_matr.empty:
            (q_matr
              .assign(IID = samples.tolist())
              .loc[
                :,
                ["IID"] + [c for c in q_matr.columns.tolist()]]
              .to_csv(
                f"{in_root}.cov",
                header=False,
                sep=" "))
        elif os.path.exists(f"{in_root}.cov"):
            os.remove(f"{in_root}.cov")

        if params["k"]:
            (k_matr
              .assign(IID = samples.tolist())
              .loc[:, ["IID"] + k_matr.columns.tolist()]
              .to_csv(
                f"{in_root}.kin",
                sep=" ",
                header=False,
                index=True))

        gctb_path = os.path.join(
          tool_dir,
          "gctb_2.5.2_Linux",
          "gctb")
        progress(task_id, "Running analysis")
        for idx, trait in enumerate(params["traits"], start=1):
            cmd = [
              gctb_path,
              "--bfile", f"{in_root}.gctb",
              "--pheno", f"{in_root}.txt",
              "--mpheno", str(idx),
              "--bayes", model,
              "--maf", "0",
              "--p-value", "1",
              "--out", os.path.join(out_dir, f"{model}.{trait}")]
            cmd = cmd + ["--burn-in", str(params["b_burn"])] if params["b_burn"] else cmd
            cmd = cmd + ["--chain-length", str(params["b_chain"])] if params["b_chain"] else cmd
            cmd = cmd + ["--thin", str(params["b_thin"])] if params["b_thin"] else cmd
            cmd = cmd + ["--covar", f"{in_root}.cov"] if not q_matr.empty else cmd
            cmd = cmd + ["--random-covar", f"{in_root}.kin"] if params["k"] else cmd
            popen(
              task_id,
              cmd)

            assoc_process(
              genome,
              assoc_dir,
              f"{model}.{trait}.snpRes",
              f"{model}.{trait}",
              out_dir,
              ch_col="Chrom",
              pos_col="Position",
              sig_col="PIP",
              eff_col=["A1Effect"],
              sep="\s+")

            # Save parameters
            run_params = [
              ["Markers", "Samples", "PCs", "Kinship", "K/Q source",
               "Burn-in iterations", "Total iterations", "Thinning"],
              [sum(1 for l in open(f"{in_root}.bim")),
               phe_df.loc[samples, trait].dropna().shape[0],
               params["q_num"] if params["q"] else 0,
               "No",
               params["kq_from"] or markers,
               params["b_burn"] if params["b_burn"] else 1000,
               params["b_chain"] if params["b_chain"] else 3000,
               params["b_thin"] if params["b_thin"] else 10]]
            with open(
              os.path.join(
                out_dir,
                f"{model}.{trait}.params"),
              "w") as params_out:
                json.dump(run_params, params_out)

        for ext in ("covRes", "parRes", "randCovRes"):
            fn = os.path.join(out_dir, f"{model}.{trait}.{ext}")
            if os.path.exists(fn):
                os.remove(fn)

    else:
        # GEMMA cannot accept row or column names in phenotype file
        # phenotype file rewritten each time to keep only selected traits
        (phe_df
          .fillna("NA")
          # GEMMA treats -9 as a missing value                                    #################################
          .replace(-9, -9.000000000001)
          .replace("-9", "-9.000000000001")
          .to_csv(
            f"{in_root}.txt",
            sep="\t",
            header=False,
            index=False))

        q_matr, k_matr = kq_process(
          genome,
          markers,
          params["kq_from"] or markers,
          samples,
          model != "glm",
          int(params["q_num"]) if params["q"] else 0,
          params["cov"])
        if not q_matr.empty:
            (q_matr
              .assign(I=1)
              .loc[:, ["I"] + q_matr.columns.tolist()]
              .to_csv(
                f"{in_root}.cov",
                sep="\t",
                header=False,
                index=False))
        elif os.path.exists(f"{in_root}.cov"):
            os.remove(f"{in_root}.cov")

        if not k_matr.empty:
            (k_matr
              .to_csv(
                f"{in_root}.kin",
                sep="\t",
                header=False,
                index=False))

#        k_matr = k_matr.loc[k_lab.isin(samples), k_lab.isin(samples)]
#        eigenvalues, eigenvectors = np.linalg.eigh(k_matr)
#        y = phe_df.loc[samples, params["traits"][0]]
#        sigma_g2 = np.var(eigenvalues)
#        sigma_e2 = np.var(y.dropna()) - sigma_g2
#        print(y)
#        print(sigma_g2 / (sigma_g2 + sigma_e2))

        progress(task_id, "Running analysis")
        if model == "mtmlm":
            traits = params["traits_group"]
        else:
            traits = [([t], t) for t in params["traits"]]
        for trait, trait_name in traits:
            n = [str(phe_df.columns.get_loc(t) + 1) for t in trait]
            command = [
              "gemma",
              "-bfile", in_root,
              #"-g", f"{in_root}.bimbam",
              #"-a", f"{in_root}.ann",
              "-p", f"{in_root}.txt",
              "-n", *n,
              "-o", f"{model}.{trait_name}",
              "-outdir", out_dir,
              "-miss", "1",
              "-maf", "0",
              "-r2", "1",
              "-hwe", "0"]
            if model == "glm":
                command += ["-lm", "1"]
            else:
                command += ["-lmm", "1"]
                command += ["-k", f"{in_root}.kin"]
            if not q_matr.empty:
                command += ["-c", f"{in_root}.cov"]
            progress(task_id, trait_name, bold=False)
            popen(
              task_id,
              command)

            assoc_process(
              genome,
              assoc_dir,
              f"{model}.{trait_name}.assoc.txt",
              f"{model}.{trait_name}",
              out_dir,
              ch_col="chr",
              pos_col="ps",
              sig_col="p_wald",
              eff_col=(
                ["beta"] if model in ("mlm", "glm")
                else ["beta_" + str(i) for i in range(1, len(n)+1)] if model == "mtmlm"
                else None),
              se_col="se" if model == "mlm" else None)

            # Save parameters
            run_params = [
              ["Markers", "Samples", "PCs", "Kinship", "K/Q source"],
              [sum(1 for l in open(f"{in_root}.bim")),
               phe_df.loc[samples, trait].dropna().shape[0],
               params["q_num"] if params["q"] else 0,
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

            log_file = os.path.join(out_dir, f"{model}.{trait_name}.log.txt")
            if os.path.exists(log_file):
                os.remove(log_file)

def kq_process(
  genome,
  markers,
  kq_from,
  samples,
  k,
  q_num,
  cov):
    kq_root = os.path.join(
      data_dir,
      genome,
      "assoc",
      kq_from,
      "gt",
      "gt")

    q_matr = pd.DataFrame()
    k_matr = pd.DataFrame()
    # Construct fixed effect covariates
    q_list = list()
    if q_num > 0:
        # Covariates file filtered to keep top n selected ones
        # and necessary samples
        q_matr = (pd
          .read_csv(
            f"{kq_root}.eigenvec",
            sep="\t",
            header=[0],
            index_col=[0,1])
          .droplevel(0))
        q_matr = q_matr.loc[:, q_matr.columns.tolist()[:q_num]]
        if not samples.isin(q_matr.index).all():
            raise OSError(
              "<b>Dataset selected for K/Q correction does not contain " +
              "all samples from the genotype dataset.<b>")
        q_list.append(q_matr)

    if cov:
        cov = pd.read_csv(
          io.StringIO(json.loads(cov)),
          sep="\t",
          header=None,
          index_col=0)
        if not samples.isin(cov.index).all():
            raise OSError(
              "<b>Additional covariates file does not contain " +
              "all samples from the genotype dataset.<b>")
        q_list.append(cov)
    if len(q_list) > 0:
        q_matr = (pd
          .concat(
            q_list,
            axis=1)
          .loc[samples])

    # Construct kinship matrix
    if k:
        k_matr = pd.read_csv(
          f"{kq_root}.rel",
          sep="\t",
          index_col=None,
          header=None)
        k_lab = pd.read_csv(
          f"{kq_root}.rel.id",
          sep="\t",
          usecols=[1],
          skiprows=[0],
          header=None,
          index_col=None)[1]
        if not k_lab.isin(samples).all():
            raise OSError(
              "<b>Dataset selected for K/Q correction does not contain " +
              "all samples from the genotype dataset.<b>")
        k_matr = (k_matr
          .set_axis(k_lab, axis=0)
          .set_axis(k_lab, axis=1)
          .loc[
            samples,
            samples])

    return q_matr, k_matr

def assoc_process(
  genome,
  assoc_dir,
  assoc_fn,
  out_root,
  out_dir,
  ch_col,
  pos_col,
  sig_col,
  sep="\t",
  eff_col=[],
  se_col=None):
    gdata = load_gdata(genome)
    config = load_config(genome)
    assoc_out = os.path.join(
      out_dir,
      assoc_fn)
    processed_out = os.path.join(out_dir, out_root)
    columns = [c for c in (ch_col, pos_col, sig_col, *eff_col) if c] #, se_col

    # Output files
    mqplot_fn   = f"{processed_out}.mqplot.png"
    out_feather = f"{processed_out}.feather"
    out_tsv     = f"{processed_out}.tsv"
    out_gz      = f"{processed_out}.tsv.gz"

    # Save raw results, rename numeric chromosomes back to original
    chr_rename = dict(enumerate(gdata["chr_names"], start=1))
    res = (
      pd
        .read_csv(
          assoc_out,
          sep=sep,
          header=[0],
          index_col=None)
        .loc[:, columns]
        .assign(**{ch_col: lambda x: x[ch_col].map(chr_rename)}))

    #res = res[~pd.isna(res[ch_col])]
    if eff_col:
        res[eff_col] = res[eff_col].abs()
    res.to_csv(
      out_tsv,
      sep="\t",
      index=False)
    pysam.tabix_compress(
      out_tsv,
      out_gz,
      force=True)
    #pysam.tabix_index(
    #  out_gz,
    #  seq_col=0,
    #  start_col=1,
    #  end_col=1,
    #  line_skip=1, force=True)

    if sig_col != "PIP":
        res[sig_col] = -np.log10(res[sig_col])
    # Save results for MP visualisation
    with open(os.path.join(assoc_dir, "stats", "params")) as params_in:
        is_hap = json.load(params_in)[1][10] == "Yes"
    if is_hap:
        blockmap = (pd
          .read_csv(
            os.path.join(assoc_dir, "gt", "gt.blockmap"),
            sep="\t",
            header=0,
            index_col=None))
        blockmap["chr"] = blockmap["chr"].astype(str)
        res = (pd
          .merge(
            res,
            blockmap,
            left_on=[ch_col, pos_col],
            right_on=["chr", "start"],
            how="left")
          .drop(
            ["start"],
            axis=1))
        res["idx"] = res.groupby([ch_col, pos_col]).cumcount() + 1
        c = res.columns.tolist()
        (res[c[2:3] + c[:2] + c[-2:] + eff_col]
          .sort_values(
            by=sig_col,
            ascending=False)
          .head(5000)
          .reset_index(drop=True)
          .to_feather(out_feather))
    else:
        c = res.columns.tolist()
        (res[c[2:3] + c[:2] + eff_col]
          .sort_values(
            by=sig_col,
            ascending=False)
          .head(5000)
          .reset_index(drop=True)
          .to_feather(out_feather))
    mqplot(
      genome,
      res,
      mqplot_fn,
      json.loads(config["concat"]),
      ch_col  = ch_col,
      pos_col = pos_col,
      sig_col = sig_col,
      qqplot  = sig_col != "PIP",
      thresh  = sig_col != "PIP",
      end_col = "end" if is_hap else None)

    os.remove(out_tsv)
    os.remove(assoc_out)

def assoc_finemap(
  genome,
  params,
  task_id):
    config  = load_config(genome)
    markers = params["mrk"]
    model   = params["mod"]
    trait   = params["trt"]
    coord   = params["coord"]
    n_loci  = len(coord.split(", "))

    out_root = os.path.join(
      tmp_dir,
      task_id)
    pimass_path = os.path.join(
      tool_dir,
      "pimass",
      "pimass-lin")
    gctb_path = os.path.join(
      tool_dir,
      "gctb_2.5.2_Linux",
      "gctb")
    markers_dir = os.path.join(
      data_dir,
      genome,
      "assoc",
      markers)
    gt_root = os.path.join(
      markers_dir,
      "gt",
      "gt.full")

    # Process phenotype file
    fam_file = f"{gt_root}.fam"
    samples = pd.read_csv(
        fam_file,
        sep="\s+",
        header=None,
        index_col=None,
        usecols=[1])[1]
    phe_path = os.path.join(
      data_dir,
      genome,
      config["phe"])
    with master_lock:
        lock = get_lock(phe_path)
    with lock:
        phe_df = pd.read_csv(
          phe_path,
          sep="\t",
          header=[0],
          index_col=[0])

    if params["fm_method"] != "pm":
        q_matr, k_matr = kq_process(
          genome,
          markers,
          params["kq_from"] or markers,
          samples,
          params["k"],
          int(params["q_num"]) if params["q"] else 0,
          params["cov"])

    # Extract genotype data at QTLs
    progress(task_id, "Extracting genotypes at QTLs")
    with open(f"{out_root}.bed1", "w") as coord_out:
        for locus in coord.split(", "):
            chrom, start = locus.split(":")
            start = start.replace(",", "")
            if "-" in start:
                start, end = start.split("-")
            else:
                end = start
            _ = coord_out.write("\t".join((chrom, start, end)) + "\n")

    if params["fm_method"] in ("B", "C", "glm", "mlm"):
        popen(
          task_id,
          [plink2_path,
           "--bfile", gt_root,
           "--extract", "bed1", f"{out_root}.bed1",
           "--make-bed",
           "--allow-extra-chr",
           "--out", out_root])

        num_samples = (phe_df
          .loc[samples, trait]
          .dropna()
          .size)
        num_vars = sum(1 for _ in open(f"{out_root}.bim"))

        if params["fm_method"] in ("B", "C"):
            fam_path = f"{out_root}.fam"
            fam_df = pd.read_csv(
              fam_path,
              sep="\t",
              header=None,
              index_col=None)
            (fam_df
              .assign(IID = samples.tolist())
              .loc[:, ["IID"] + fam_df.columns.tolist()[1:]]
              .to_csv(
                fam_path,
                sep=" ",
                header=False,
                index=False))

            (phe_df
              .assign(IID = phe_df.index.tolist())
              .loc[:,
                ["IID", trait]]
              .to_csv(
                f"{out_root}.txt",
                header=False,
                sep=" "))

            if not q_matr.empty:
                (q_matr
                  .assign(IID = samples.tolist())
                  .loc[
                    :,
                    ["IID"] + [c for c in q_matr.columns.tolist()]]
                  .to_csv(
                    f"{out_root}.cov",
                    header=False,
                    sep=" "))
            if params["k"]:
                (k_matr
                  .assign(IID = samples.tolist())
                  .loc[:, ["IID"] + k_matr.columns.tolist()]
                  .to_csv(
                    f"{out_root}.kin",
                    sep=" ",
                    header=False,
                    index=True))

            progress(task_id, "Fine-mapping with Bayes" + params["fm_method"])
            cmd = [
              gctb_path,
              "--bfile", f"{out_root}",
              "--pheno", f"{out_root}.txt",
              "--bayes", params["fm_method"],
              "--maf", "0",
              "--p-value", "1",
              "--out", out_root]
            cmd = cmd + ["--burn-in", str(params["b_burn"])] if params["b_burn"] else cmd
            cmd = cmd + ["--chain-length", str(params["b_chain"])] if params["b_chain"] else cmd
            cmd = cmd + ["--thin", str(params["b_thin"])] if params["b_thin"] else cmd
            cmd = cmd + ["--covar", f"{out_root}.cov"] if not q_matr.empty else cmd
            cmd = cmd + ["--random-covar", f"{out_root}.kin"] if params["k"] else cmd
            popen(
              task_id,
              cmd)

            pip_in = f"{out_root}.snpRes"
            pip_out = os.path.join(
              tmp_dir,
              f"{task_id}.fm")

            (pd
              .read_csv(
                pip_in,
                sep="\s+",
                header=0,
                index_col=None,
                usecols=[2,3,10])
              .set_axis(
                ["chr", "pos", "val"],
                axis=1)
              .sort_values(
                by="val",
                ascending=False)
              .head(5000)
              .reset_index(drop=True)
              .to_feather(
                pip_out))

            # Save parameters
            run_params = [
              ["Marker set", "Method", "Trait", "Loci", "Markers", "Samples"],
              [markers, "Bayes" + params["fm_method"], trait,
               "<br>".join(coord.split(", ")), num_vars, num_samples]]
            with open(
              os.path.join(
                tmp_dir,
                f"{task_id}.params"),
              "w") as params_out:
                json.dump(run_params, params_out)

        else:
            (phe_df
              .fillna("NA")
              # GEMMA treats -9 as a missing value                                    #################################
              .replace(-9, -9.000000000001)
              .replace("-9", "-9.000000000001")
              .loc[samples, trait]
              .to_csv(
                f"{out_root}.txt",
                sep="\t",
                header=False,
                index=False))

            if not q_matr.empty:
                (q_matr
                  .assign(I=1)
                  .loc[:, ["I"] + q_matr.columns.tolist()]
                  .to_csv(
                    f"{out_root}.cov",
                    sep="\t",
                    header=False,
                    index=False))
            if not k_matr.empty:
                (k_matr
                  .to_csv(
                    f"{out_root}.kin",
                    sep="\t",
                    header=False,
                    index=False))

            progress(task_id, f"Fine-mapping with {params['fm_method'].upper()}")
            command = [
              "gemma",
              "-bfile", out_root,
              "-p", f"{out_root}.txt",
              "-o", task_id,
              "-outdir", "tmp",
              "-miss", "1",
              "-maf", "0",
              "-r2", "1",
              "-hwe", "0"]
            if params["fm_method"] == "glm":
                command += ["-lm", "1"]
            else:
                command += ["-lmm", "1"]
                command += ["-k", f"{out_root}.kin"]
            if not q_matr.empty:
                command += ["-c", f"{out_root}.cov"]
            popen(
              task_id,
              command)

            res_in = f"{out_root}.assoc.txt"
            res_out = os.path.join(
              tmp_dir,
              f"{task_id}.fm")

            (pd
              .read_csv(
                res_in,
                sep="\t",
                header=0,
                index_col=None,
                usecols=[0,2,11])
              .set_axis(
                ["chr", "pos", "val"],
                axis=1)
              .assign(val=lambda d: -np.log10(d["val"]))
              .sort_values(
                by="val",
                ascending=False)
              .head(5000)
              .reset_index(drop=True)
              .to_feather(
                res_out))

            # Save parameters
            run_params = [
              ["Marker set", "Method", "Trait", "Loci", "Markers", "Samples"],
              [markers, params["fm_method"].upper(), trait, "<br>".join(coord.split(", ")), num_vars, num_samples]]
            with open(
              os.path.join(
                tmp_dir,
                f"{task_id}.params"),
              "w") as params_out:
                json.dump(run_params, params_out)

    else:
        popen(
          task_id,
          [plink2_path,
           "--bfile", gt_root,
           "--extract", "bed1", f"{out_root}.bed1",
           "--export", "Av",
           "--allow-extra-chr",
           "--out", out_root])

        if params["fm_method"] == "rf":
            if params["optimize"]:
                progress(task_id, "Running optimization")
            else:
                progress(task_id, "Fine-mapping with Random forest")
            phe = (phe_df
              .loc[samples, trait]
              .dropna())
            phe_samples = phe.index
            gt_df = (pd
              .read_csv(
                f"{out_root}.traw",
                sep="\t",
                header=0,
                index_col=None))

            if not q_matr.empty:
                phe_samples = phe_samples.intersection(q_matr.index)
                phe = phe.loc[phe_samples]
                phe_q = q_matr.loc[phe_samples]
                phe_q = sm.add_constant(phe_q)
                if params["k"]:
                    phe_k = k_matr.loc[phe_samples, phe_samples]
                    QS = economic_qs_linear(phe_k)
                    lmm = LMM(
                      phe,
                      phe_q,
                      QS=QS)
                    lmm.fit()
                    phe_corr = phe - lmm.mean()
                else:
                    phe_corr = sm.OLS(phe, phe_q).fit().resid
            else:
                phe_corr = phe

            gt = (gt_df
              .drop(
                ["CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT"],
                axis=1)
              .set_axis(
                samples,
                axis=1)
              .loc[:, phe_samples]
              .transpose())

            if params["optimize"]:
                mtryvals = [
                  int(v.strip())
                  for v in params["mtryvals"].split(", ")]
                ntreevals = [
                  int(v.strip())
                  for v in params["ntreevals"].split(", ")]
                res = optim_rf(
                  task_id,
                  gt,
                  phe_df,
                  mtryvals,
                  ntreevals)
                optim_path = os.path.join(
                  tmp_dir,
                  f"{task_id}.rfopt")
                with open(optim_path, "w") as optim_out:
                    json.dump([mtryvals, ntreevals, res], optim_out)
                return

            res = run_rf(
              task_id,
              gt,
              phe_corr,
              params["mtry"],
              params["ntree"],
              params["niter"],
              params["nperm"])
            imp_out = os.path.join(
              tmp_dir,
              f"{task_id}.fm")
            gt_map = gt_df[["CHR", "POS"]]
            res["Chr"]  = gt_map["CHR"]
            res["Pos"] = gt_map["POS"]
            #res["Imp"] = res["Imp"] / res["Imp"].max()
            (res
              .loc[:, ["Chr", "Pos", "Imp"]]
              .set_axis(
                ["chr", "pos", "val"],
                axis=1)
              .sort_values(
                by="val",
                ascending=False)
              .head(5000)
              .to_feather(
                imp_out))

            # Save parameters
            run_params = [
              ["Marker set", "Method", "Trait", "Loci",
               "Markers", "Samples", "Test markers",
               "Trees per forest", "Iterations", "Permutations"],
              [markers, "Random forest", trait, 
               "<br>".join(coord.split(", ")), len(gt.columns), phe_df.size,
               params["mtry"], params["ntree"], params["niter"], params["nperm"]]]
            with open(
              os.path.join(
                tmp_dir,
                f"{task_id}.params"),
              "w") as params_out:
                json.dump(run_params, params_out)

        elif params["fm_method"] == "pm":
            (phe_df
              .loc[samples, trait]
              .fillna("NA")
              .to_csv(
                f"{out_root}.phe",
                sep="\t",
                header=False,
                index=False))
            num_samples = (phe_df
              .loc[samples, trait]
              .dropna()
              .size)

            # Create a BIMBAM-format genotype file
            gt = (pd
              .read_csv(
                f"{out_root}.traw",
                sep="\t",
                header=[0],
                index_col=None)
              .drop(
                ["CHR", "(C)M", "POS"],
                axis=1))
            samp_cols = gt.columns.tolist()[3:]
            allele_cnt = len(samp_cols) * 2
            marker_cnt = gt.shape[0]
            gt["AF"] = gt.iloc[:, 3:].sum(axis=1) / allele_cnt
            switch_mask = gt["AF"] > 0.5
            gt.loc[switch_mask, samp_cols] = (gt.loc[switch_mask, samp_cols] - 2).abs()
            gt["COUNTED"] = "T"
            gt["ALT"] = "A"
            (gt
              .drop(["AF"], axis=1)
              .to_csv(
                f"{out_root}.bimbam",
                sep="\t",
                header=False,
                index=False))

            progress(task_id, "Calculating marker PIP")
            popen(
              task_id,
              [pimass_path,
               "-g", f"{out_root}.bimbam",
               "-p", f"{out_root}.phe",
               "-o", task_id,
               #"-r", "0",
               # MCMC sampling parameters
               "-w", str(params["pm_w"] or marker_cnt),
               "-s", str(params["pm_s"] or marker_cnt * 10),
               "-num", str(params["pm_num"] or 10),
               # Min / max expected heritability
               "-hmin", str(params["pm_hmin"] or 0.3),
               "-hmax", str(params["pm_hmax"] or 0.8),
               # Min / max expected causal SNPs
               "-pmin", str(params["pm_pmin"] or n_loci // 3),
               "-pmax", str(params["pm_pmax"] or n_loci * 3),
               # Number of causal SNPs to start with
               "-nstart", str(params["pm_nstart"] or n_loci),
               # Min / max number of SNPs to include in model
               "-smin", str(params["pm_smin"] or marker_cnt // 10),
               "-smax", str(params["pm_smax"] or marker_cnt // 2),
               ])
            pip_in = os.path.join(
              "output",
              f"{task_id}.mcmc.txt")
            pip_out = os.path.join(
              tmp_dir,
              f"{task_id}.fm")

            (pd
              .read_csv(
                pip_in,
                sep="\s+",
                header=0,
                index_col=0,
                usecols=[0,4])
              .assign(
                chr=lambda df: df.index.str.split(":").str[0],
                pos=lambda df: df.index.str.split(":").str[1].astype(int))
              [["chr", "pos", "postrb"]]
              .set_axis(
                ["chr", "pos", "val"],
                axis=1)
              .sort_values(
                by="val",
                ascending=False)
              .head(5000)
              .to_feather(
                pip_out))

            # Save parameters
            run_params = [
              ["Marker set", "Method", "Trait", "Loci", "Markers", "Samples", 
               "Burn-in steps", "Sampling steps", "Thinning", "Min heritability",
               "Max heritability", "Min causal SNPs", "Max causal SNPs", 
               "Init causal SNPs", "Model min SNPs", "Model max SNPs"],
              [markers, "PiMASS", trait,
               "<br>".join(coord.split(", ")), gt.shape[0], num_samples,
               params["pm_w"] or marker_cnt,
               params["pm_s"] or marker_cnt * 10,
               params["pm_num"] or 10,
               float(params["pm_hmin"] or 0.3),
               float(params["pm_hmax"] or 0.8),
               params["pm_pmin"] or n_loci // 3,
               params["pm_pmax"] or n_loci * 3,
               params["pm_nstart"] or n_loci,
               params["pm_smin"] or marker_cnt // 10,
               params["pm_smax"] or marker_cnt // 2,
               ]]
            with open(
              os.path.join(
                tmp_dir,
                f"{task_id}.params"),
              "w") as params_out:
                json.dump(run_params, params_out)
