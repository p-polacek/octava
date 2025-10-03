#!/usr/bin/env python
import os
import gzip
import json
import pysam
import pgenlib
import bed_reader
import numpy as np
import pandas as pd


def win_block(gt_root, btype, span, step, task_id):
    bim = pd.read_csv(
        f"{gt_root}.bim", sep="\t", header=None, index_col=None, usecols=[0, 1, 3]
    ).set_axis(["Chr", "ID", "Pos"], axis=1)

    if btype == "bpwin":
        total = bim.groupby("Chr")["Pos"].max().sum()
    elif btype == "snpwin":
        total = bim.shape[0]
    done_cnt = 0

    pct_done = 0
    block = None
    last_block = None
    out = open(f"{gt_root}.blocks", "w")

    for ch, data in bim.groupby("Chr", sort=False):
        if btype == "bpwin":
            last_marker = data.iloc[-1, 2]
            data = data.set_index("Pos")
        elif btype == "snpwin":
            data.reset_index(inplace=True, drop=True)
            last_marker = data.shape[0] - 1

        for start in range(0, last_marker + 1, step):
            pct = int((start + done_cnt) / total * 100)
            if pct > pct_done:
                pct_done = pct
                print(f"\rProgress: {pct}%", end="")
            end = start + span
            block = data.loc[start:end, "ID"].tolist()
            if len(block) < 2:
                continue
            if last_block and block == last_block:
                continue
            _ = out.write(" ".join(["*"] + block) + "\n")
        done_cnt += last_marker + 1

    print("\rProgress: 100%\n")
    out.close()


def block_genotype(
    gt_root,
    stats_dir,
    #  out_format,
    #  transpose,
    maf,
    min_haps,
    concat,
    chr_names,
    #  group_minor,
    #  include_singles
    task_id,
):
    min_haps = max(1, min_haps)
    samples = pd.read_csv(
        f"{gt_root}.fam", sep="\t", header=None, index_col=None, usecols=[1]
    )[1].to_numpy()
    sample_cnt = samples.size
    allele_cnt = sample_cnt * 2
    allele_thresh = allele_cnt * maf
    # Summary statistics and marker info
    smiss, shet = {s: 0 for s in samples}, {s: 0 for s in samples}
    vmiss, vhet = list(), list()
    spans, mcnts, hcnts, fhcnts, afs = list(), list(), list(), list(), dict()
    blocks, bchrs, chrs, ps, ids, starts, ends = (
        list(),
        list(),
        list(),
        list(),
        list(),
        list(),
        list(),
    )
    # Load files
    in_b = open(f"{gt_root}.blocks")
    total = sum(1 for _ in open(f"{gt_root}.blocks"))
    pct_done = 0
    gt = pgenlib.PgenReader(f"{gt_root}.pgen".encode("utf8"))
    out_gt = list()
    out_bed = f"{gt_root}.bed"
    out_gt_allhaps = list()
    out_bed_allhaps = f"{gt_root}.allhaps.bed"
    chrs_all, ps_all, ids_all = list(), list(), list()
    # Marker to index map
    marker_map = {
        marker: idx
        for (idx, marker) in pd.read_csv(
            f"{gt_root}.bim", sep="\t", header=None, usecols=[1]
        ).to_records()
    }

    vcf_out = pysam.BGZFile(gt_root + ".vcf.gz", "w")
    _ = vcf_out.write("##fileformat=VCFv4.2\n".encode())
    columns = [
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
    ] + samples.tolist()
    _ = vcf_out.write(("\t".join(columns) + "\n").encode())

    print("<b>Genotyping haplotype blocks</b>")
    # Skip leading "*", split on space
    for block_idx, block in enumerate(in_b):
        pct = int(block_idx / total * 100)
        if pct > pct_done:
            pct_done = pct
            print(f"\rProgress: {pct}%", end="")
        m1, mn = block.split(None, 2)[1], block.rsplit(None, 1)[1].strip()
        start, end = marker_map[m1], marker_map[mn] + 1
        mcnt = end - start

        block_gt = np.empty([mcnt, allele_cnt], np.int32)
        phased = np.empty([mcnt, sample_cnt], np.int8)
        gt.read_alleles_and_phasepresent_range(start, end, block_gt, phased)
        block_gt = block_gt.transpose()
        haps, cnts = np.unique(block_gt, axis=0, return_counts=True)
        filt_haps = haps[cnts > allele_thresh]
        filt_cnts = cnts[cnts > allele_thresh]

        # Skip block if not enough haplotypes pass MAF threshold
        if len(filt_haps) < min_haps:
            continue

        # Save summary statistics
        ch, first_pos = m1.split(":")
        last_pos = mn.split(":")[1]
        span = int(last_pos) - (int(first_pos) - 1)
        blocks.append([ch, first_pos, last_pos, span])
        spans.append(span)
        mcnts.append(mcnt)
        hcnts.append(haps.shape[0])
        fhcnts.append(filt_haps.shape[0])
        for idx, af in enumerate(np.sort(filt_cnts / allele_cnt)[::-1], start=1):
            afs.setdefault(idx, list()).append(af)

        filt_haps = filt_haps[filt_cnts.argsort()[::-1]]
        for idx in range(filt_haps.shape[0]):
            hap_chromatids = np.all(block_gt == filt_haps[idx], axis=1)
            hap_gt = np.sum(hap_chromatids.reshape(-1, 2), axis=1)
            out_gt.append(hap_gt)
            # Save marker info
            chrs.append(ch)
            ps.append(first_pos)
            _id = f"{m1}-{last_pos}_{idx+1}"
            ids.append(_id)
        bchrs.append(ch)
        starts.append(first_pos)
        ends.append(last_pos)

        # Get haplotype block VCF entry
        matches = (block_gt[:, None, :] == filt_haps[None, :, :]).all(axis=2)
        idx_array = matches.argmax(axis=1) + 1
        idx_array[~matches.any(axis=1)] = -1
        idx_array = idx_array.reshape(-1, 2)
        idx_array = np.where(idx_array == -1, ".", idx_array.astype(str))

        # Determine misingness and heterozygosity
        block_gts = np.stack(out_gt[-filt_haps.shape[0] :]).transpose()
        block_vmiss = np.all(block_gts == 0, axis=1)
        vmiss.append(sum(block_vmiss) / sample_cnt)
        for s, m in zip(samples, block_vmiss):
            smiss[s] += m
        block_vhet = np.any(block_gts[~block_vmiss] == 1, axis=1)
        vhet.append(sum(block_vhet) / block_vhet.size)
        for s, m in zip(samples[~block_vmiss], block_vhet):
            shet[s] += m

        # Keep all haplotypes for kinship and PC calculation
        for idx in range(haps.shape[0]):
            hap_chromatids = np.all(block_gt == haps[idx], axis=1)
            hap_gt = np.sum(hap_chromatids.reshape(-1, 2), axis=1)
            out_gt_allhaps.append(hap_gt)

            chrs_all.append(ch)
            ps_all.append(first_pos)
            _id = f"{m1}-{last_pos}_{idx+1}"
            ids_all.append(_id)

        # write VCF entry
        _ = vcf_out.write(
            (
                "\t".join(
                    (
                        ch,
                        first_pos,
                        _id,
                        "A" * span,
                        ",".join(["T" * span] * len(filt_haps)),
                        ".",
                        "PASS",
                        ".",
                        "GT",
                        *[f"{a1}|{a2}" for a1, a2 in idx_array],
                    )
                )
                + "\n"
            ).encode()
        )

    print("\rProgress: 100%\n")
    print("<b>Saving genotyped blocks</b>")

    in_b.close()
    # Filtered data
    (
        pd.DataFrame({"chr": bchrs, "start": starts, "end": ends}).to_csv(
            gt_root + ".blockmap", sep="\t", index=False
        )
    )
    if not out_gt:
        raise OSError("<b>No haplotype blocks were constructed.</b>")
    bed_reader.to_bed(
        out_bed,
        np.stack(out_gt, axis=0).transpose(),
        properties={
            "iid": samples,
            "pheno": [-9] * sample_cnt,
            "chromosome": chrs,
            "sid": ids,
            "bp_position": ps,
        },
    )

    bim = pd.read_csv(f"{gt_root}.bim", sep="\t", header=None, index_col=None)
    (
        bim.set_axis(["CHR", "ID", "_", "POS", "A1", "A2"], axis=1)
        .assign(A1="A")
        .assign(A2="T")
        .to_csv(f"{gt_root}.bim", sep="\t", header=False, index=False)
    )
    # Full data for kinship and PC calculation
    bed_reader.to_bed(
        out_bed_allhaps,
        np.stack(out_gt_allhaps, axis=0).transpose(),
        properties={
            "iid": samples,
            "chromosome": chrs_all,
            "sid": ids_all,
            "bp_position": ps_all,
            "pheno": [-9]*sample_cnt
        },
    )
    pd.read_csv(f"{gt_root}.fam", sep=" ", header=None, index_col=None).to_csv(
        f"{gt_root}.fam", sep="\t", header=False, index=False
    )

    block_map = pd.DataFrame(data=blocks)
    for chrom, data in block_map.groupby(0):
        if chrom not in concat:
            data.loc[:, data.columns[1:]].to_csv(
                os.path.join(stats_dir, f"{chrom}.bdist.gz"),
                sep="\t",
                index=False,
                header=False,
                compression="gzip",
            )
    block_num = block_map.shape[0]

    block_chrs = block_map[0].unique()
    for chrom in chr_names:
        if chrom not in block_chrs:
            (gzip.open(os.path.join(stats_dir, f"{chrom}.bdist.gz"), "wb").close())

    smiss = {k: v / block_num for k, v in smiss.items()}
    shet = {k: v / block_num for k, v in shet.items()}
    # Save histogram data
    for stat, ext in zip(
        (list(smiss.values()), list(shet.values()), vmiss, vhet),
        ("smiss", "shet", "vmiss", "vhet"),
    ):
        bins = [
            v / 100 for v in range(int(min(min(stat), 0)), int(max(stat) * 100) + 2)
        ]
        hgts, bins = np.histogram(stat, bins=bins)
        with open(os.path.join(stats_dir, ext), "w") as hist_out:
            json.dump([hgts.tolist(), bins.tolist()], hist_out)
    for stat, ext in zip((smiss, shet), ("smiss", "shet")):
        (
            pd.DataFrame(list(stat.items()))
            .set_index(0)
            .loc[samples, :]
            .to_csv(
                os.path.join(stats_dir, f"stats.{ext}.gz"), sep="\t", compression="gzip"
            )
        )

    for stat, ext in zip((spans, mcnts, hcnts, fhcnts), ("bsp", "bmc", "bhc", "bfhc")):
        step = 1000 if stat == spans else 1
        hgts, bins = np.histogram(stat, bins=range(0, max(stat) + step, step))
        with open(os.path.join(stats_dir, ext), "w") as hist_out:
            json.dump([hgts.tolist(), bins.tolist()], hist_out)
    afs = [
        np.histogram(afs[k], bins=[i / 100 for i in range(0, 101, 1)])
        for k in sorted(afs.keys())
    ]

    with open(os.path.join(stats_dir, "cnts"), "w") as cnts_out:
        json.dump([sample_cnt, len(ids)], cnts_out)

    bins = afs[0][1].tolist()
    afs = [i[0].tolist() for i in afs]
    with open(os.path.join(stats_dir, "baf"), "w") as hist_out:
        json.dump([afs, bins], hist_out)
