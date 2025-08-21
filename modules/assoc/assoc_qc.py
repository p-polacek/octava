#!/usr/bin/env python

import os
import sys
import gzip
import json
import pysam
import tarfile

# import pgenlib
import subprocess
import portalocker
import numpy as np
import pandas as pd
from collections import Counter
from collections import defaultdict

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from app.dirs import data_dir, tmp_dir, plink2_path

from modules.mem import get_mem

from modules.load_input import load_gdata, load_config

# Bin size for making marker density histograms
mdbin = 1000000



def run_qc(genome, task_id, markers=None, data_ext=False):
    config = load_config(genome)
    gdata = load_gdata(genome)

    if not markers and data_ext:
        global data_dir
        ext_data_dir = data_dir
        data_dir = tmp_dir
        os.link(
            os.path.join(ext_data_dir, genome, config["vcf"]),
            os.path.join(data_dir, genome, config["vcf"])
        )

    if not markers:
        in_path = os.path.join(data_dir, genome, config["vcf"])
        stats_dir = os.path.join(data_dir, genome, "stats")
        gt_out = os.path.join(stats_dir, "gt")
        stats_out = gt_out
    else:
        stats_dir = os.path.join(data_dir, genome, "assoc", markers, "stats")
        in_path = os.path.join(data_dir, genome, "assoc", markers, "gt", "gt")
        gt_out = in_path
        stats_out = os.path.join(data_dir, genome, "assoc", markers, "stats", "gt")

    if not os.path.exists(stats_dir):
        os.makedirs(stats_dir)

    if not os.path.exists(stats_dir):
        os.makedirs(stats_dir)

    if any(
        (
            not os.path.exists(f"{stats_out}.{ext}")
            for ext in ("smiss", "vmiss", "het", "gcount", "afreq", "hardy")
        )
    ):
        print("<b>Calculating genotype statistics</b>")
        if markers:
            subprocess.run(
                [
                    plink2_path,
                    "--memory",
                    f"{get_mem()}",
                    "--pfile",
                    in_path,
                    "vzs",
                    "--missing",
                    "scols=nmiss,nobs",
                    "vcols=nmiss,nobs",
                    "--geno-counts",
                    "cols=refalt1,nobs",
                    "--het",
                    "cols=het,nobs",
                    "--freq",
                    "cols=reffreq",
                    "--hardy",
                    "cols=p",
                    "--allow-extra-chr",
                    "--out",
                    stats_out,
                ],
                check=True,
            )
        else:
            subprocess.run(
                [
                    plink2_path,
                    "--memory",
                    f"{get_mem()}",
                    "--vcf",
                    in_path,
                    "--missing",
                    "scols=nmiss,nobs",
                    "vcols=nmiss,nobs",
                    "--geno-counts",
                    "cols=refalt1,nobs",
                    "--het",
                    "cols=het,nobs",
                    "--bad-freqs",
                    "--freq",
                    "cols=reffreq",
                    "--hardy",
                    "cols=p",
                    "--make-pgen",
                    "vzs",
                    "pvar-cols=",
                    "psam-cols=",
                    # "--export", "vcf", "bgz",
                    "--allow-extra-chr",
                    "--set-all-var-ids",
                    "@:#",
                    "--max-alleles",
                    "2",
                    "--import-max-alleles",
                    "255",
                    "--polyploid-mode",
                    "missing",
                    "--vcf-half-call",
                    "missing",
                    "--out",
                    stats_out,
                ],
                check=True,
            )

    metrics = ("smiss", "vmiss", "shet", "vhet", "maf", "hwe")
    exts = ("smiss", "vmiss", "het", "gcount", "afreq", "hardy")
    id_cols = ("IID", "ID", "IID", "ID", "ID", "ID")
    usecols = {
        "smiss": [1, 2],
        "vmiss": [1, 2],
        "shet": [0, 1, 3],
        "vhet": [0, 1, 2],
        "maf": [1],
        "hwe": [2],
    }

    print("<b>Processing genotype statistics</b>")
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
            usecols=usecols[metric],
        )
        het_out = os.path.join(stats_dir, f"{metric}.txt.gz")
        if os.path.exists(het_out):
            os.remove(het_out)

        for chk in df:
            chk.columns = pd.RangeIndex(chk.columns.size)
            if metric in ("smiss", "vmiss"):
                cnt.update(chk[0] / chk[1])
            elif metric == "shet":
                chk[3] = chk[1] / chk[2]
                chk = chk.loc[:, [0, 3]]
                cnt.update(chk[3])
            elif metric == "vhet":
                chk[3] = chk[1] / chk[2]
                chk = chk.loc[:, [0, 3]]
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
                    compression="gzip",
                )
        cnt = {(-1 if np.isnan(v) else v): cnt[v] for v in cnt.keys()}
        cnt_vals = list(cnt.keys())
        if metric == "hwe":
            bins = list(range(int(min(min(cnt), 0)), int(max(cnt_vals)) + 2))
        else:
            bins = [
                v / 100
                for v in range(int(min(min(cnt), 0)), int(max(cnt_vals) * 100) + 2)
            ]
        weights = [cnt[v] for v in cnt_vals]
        hgts, bins = np.histogram(cnt_vals, weights=weights, bins=bins)
        hgts = hgts.tolist()
        bins = bins.tolist()
        with open(os.path.join(stats_dir, metric), "w") as hist_out:
            json.dump([hgts, bins], hist_out)

        if metric == "smiss":
            samp_num = sum(weights)
        elif metric == "vmiss":
            var_num = sum(weights)

    with open(os.path.join(stats_dir, "cnts"), "w") as cnts_out:
        json.dump([samp_num, var_num], cnts_out)

    # Keep full data for sample missingness and heterozygosity
    if not markers:
        smiss_in = f"{stats_out}.smiss"
        smiss_df = pd.read_csv(
            smiss_in, sep="\t", index_col=0, usecols=[0] + usecols["smiss"]
        ).rename_axis("Acc")
        smiss_df = (
            (smiss_df["MISSING_CT"] / smiss_df["OBS_CT"])
            .rename("Missingness")
            .round(5)
        )

        shet_in = f"{stats_out}.het"
        shet_df = pd.read_csv(
            shet_in, sep="\t", index_col=0, usecols=[0] + usecols["shet"]
        ).rename_axis("Acc")
        shet_df = (
            (shet_df["O(HET)"] / shet_df["OBS_CT"])
            .rename("Heterozygosity")
            .round(5)
        )

        if data_ext:
            samp_path = os.path.join(data_dir, genome, "stats", "samp_ext.txt")
            pd.concat([smiss_df, shet_df], axis=1).to_csv(samp_path, sep="\t")
        else:
            samp_path = os.path.join(data_dir, genome, config["samp"])
            with portalocker.Lock(samp_path, "r+", flags=portalocker.LOCK_EX) as locked_f:
                samp_df = pd.read_csv(locked_f, sep="\t", header=0, index_col=0)
                samp_df = pd.concat(
                    [samp_df, smiss_df, shet_df],
                    axis=1)
                locked_f.seek(0)
                samp_df.to_csv(locked_f, sep="\t")
                locked_f.truncate()

    # LD decay scatter plot
    print("<b>Calculating LD per distance</b>")
    max_dist = str(max(gdata["chr_sizes"]) // 2000)
    subprocess.run(
        [
            plink2_path,
            "--memory",
            f"{get_mem()}",
            "--pfile",
            gt_out,
            "vzs",
            "--thin-count",
            "10000",
            "--r2-unphased",
            "zs",
            "cols=pos",
            "--ld-window-kb",
            max_dist,
            "--ld-window-r2",
            "0",
            "--allow-extra-chr",
            "--out",
            stats_out,
        ],
        check=True,
    )
    ldd_fn = os.path.join(stats_dir, "ldd")
    bin_sums = defaultdict(float)
    bin_counts = defaultdict(int)

    for chunk in pd.read_csv(f"{stats_out}.vcor.zst", sep="\t", chunksize=1000000):
        chunk["dist"] = abs(chunk["#POS_A"] - chunk["POS_B"])
        chunk["bin"] = chunk["dist"] // 10000 * 10000 + 10000
        for bin_, group in chunk.groupby("bin"):
            bin_sums[bin_] += group["UNPHASED_R2"].sum()
            bin_counts[bin_] += group["UNPHASED_R2"].count()

    bins = sorted(bin_sums.keys())
    mean_r2 = [round(bin_sums[b] / bin_counts[b], 3) for b in bins]
    (pd.DataFrame({"bin": bins, "R2": mean_r2}).to_csv(ldd_fn, sep="\t"))

    # Histogram of marker density
    print("<b>Calculating marker density</b>")
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
            (ch_sizes[c] // mdbin * mdbin + mdbin) for c in ch_names[:idx]
        )
    if concat:
        bin_offsets["Other"] = sum(
            (ch_sizes[c] // mdbin * mdbin + mdbin) for c in ch_names if c not in concat
        )
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
        usecols=[0, 1],
    )

    for chk in df:
        for ch, pos in chk.groupby(0):
            ch = str(ch)
            if ch in concat:
                pos[1] += ch_offsets[ch] - other_offset
                ch = "Other"
            cnts.update((pos[1] + bin_offsets[ch]) // mdbin * mdbin + mdbin)

    for ch in ch_names:
        size = ch_sizes[ch] // mdbin * mdbin + mdbin
        for block in range(
            bin_offsets[ch] + mdbin, size + bin_offsets[ch] + mdbin, mdbin
        ):
            hgts.append(cnts.get(block, 0))
            bins.append(block - mdbin // 2)
            ch_block = block - bin_offsets[ch]
            coords.append(f"{ch}:{ch_block - mdbin:,}-{ch_block - 1:,}")

    with open(os.path.join(stats_dir, "dens"), "w") as dens_out:
        json.dump([hgts, bins, coords], dens_out)

    # LD matrices per chromosome
    print("<b>Calculating linkage disequilibrium</b>")
    calc_chr_ld(genome, markers=markers)

    if not markers:
        var_only_root = os.path.join(stats_dir, "var_only")
        # Create a vars-only vcf file
        if not os.path.exists(f"{var_only_root}.vcf.gz"):
            subprocess.run(
                [
                    plink2_path,
                    "--memory",
                    f"{get_mem()}",
                    "--pfile",
                    gt_out,
                    "vzs",
                    "--make-just-pvar",
                    "cols=vcfheader,qual,filter,info",
                    "--allow-extra-chr",
                    "--out",
                    var_only_root,
                ],
                check=True,
            )
            os.rename(f"{var_only_root}.pvar", f"{var_only_root}.vcf")
            pysam.tabix_compress(f"{var_only_root}.vcf", f"{var_only_root}.vcf.gz")
            pysam.tabix_index(f"{var_only_root}.vcf.gz", preset="vcf")

            os.remove(f"{var_only_root}.vcf")
            os.remove(f"{var_only_root}.log")

    # Clean up
    for fn in os.listdir(stats_dir):
        if fn.startswith("gt") and not fn.endswith(
            (".gz", ".pgen", ".pvar.zst", ".psam")
        ):
            os.remove(os.path.join(stats_dir, fn))

    if not markers and data_ext:
        init_path = os.path.join(stats_dir, "init")
        open(init_path, "w").close()

        tar_path = os.path.join(ext_data_dir, f"{task_id}.out.tar.gz")
        with tarfile.open(tar_path, "w:gz") as tar:
            arcpath = os.path.join(ext_data_dir, genome, "stats")
            tar.add(
                stats_dir,
                arcname=arcpath
            )
        print(f"Your files are ready at {tar_path}")

def calc_chr_ld(genome, markers=None):
    config = load_config(genome)
    gdata = load_gdata(genome)

    if not markers:
        stats_dir = os.path.join(data_dir, genome, "stats")
        in_path = os.path.join(stats_dir, "gt")
        gt_out = stats_out = in_path
    else:
        stats_dir = os.path.join(data_dir, genome, "assoc", markers, "stats")
        in_path = os.path.join(data_dir, genome, "assoc", markers, "gt", "gt")
        gt_out = in_path
        stats_out = os.path.join(data_dir, genome, "assoc", markers, "stats", "gt")

    # Determine empty chromosomes (skipped for LD calculation)
    chr_s = set()
    chr_data = pd.read_csv(
        gt_out + ".pvar.zst",
        sep="\t",
        usecols=[0],
        skiprows=1,
        header=None,
        index_col=None,
        chunksize=1000000,
    )
    for chk in chr_data:
        chr_s.update(chk[0].unique().astype(str))

    concat = json.loads(config["concat"])
    for ch, span in zip(gdata["chr_names"], gdata["chr_sizes"]):

        r2_path = os.path.join(stats_dir, f"{ch}.r2.gz")
        if ch in concat:
            continue
        if ch not in chr_s:
            with gzip.open(r2_path, "wt") as r2_out:
                r2_out.write(json.dumps([[], []]))
            continue

        subprocess.run(
            [
                plink2_path,
                "--memory",
                f"{get_mem()}",
                "--pfile",
                gt_out,
                "vzs",
                "--chr",
                str(ch),
                "--bp-space",
                str(span // 2000),
                "--r2-unphased",
                "square0",
                "--out",
                stats_out,
            ],
            check=True,
        )

        r2_vals = (
            pd.read_csv(
                f"{stats_out}.unphased.vcor2", sep="\t", header=None, index_col=None
            )
            .round(3)
            .fillna("")
        )
        r2_idx = pd.read_csv(
            f"{stats_out}.unphased.vcor2.vars", sep="\t", header=None, index_col=None
        )
        with gzip.open(r2_path, "wt") as r2_out:
            r2_out.write(json.dumps([r2_vals.values.tolist(), r2_idx[0].tolist()]))

        os.remove(f"{stats_out}.unphased.vcor2")
        os.remove(f"{stats_out}.unphased.vcor2.vars")

if __name__ == "__main__":
    genome, task_id = sys.argv[1:3]
    run_qc(genome, task_id, markers=None, data_ext=True)
