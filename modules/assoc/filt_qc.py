#!/usr/bin/env python

import os
import sys
import json
import pysam
import shutil
import bed_reader
import subprocess
import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from app.dirs import data_dir, tool_dir, plink1_path, plink2_path

from modules.block import win_block, block_genotype

from modules.mem import get_mem

from modules.wpgma import wpgma

from modules.load_input import load_gdata, load_config

from modules.assoc.assoc_qc import run_qc, calc_chr_ld


def filt_qc(genome, data, task_id):

    gdata = load_gdata(genome)
    config = load_config(genome)
    markers = data["name"]

    stats_dir = os.path.join(data_dir, genome, "stats")
    assoc_dir = os.path.join(data_dir, genome, "assoc", markers)
    filt_dir = os.path.join(assoc_dir, "gt")
    plink_in = os.path.join(stats_dir, "gt")
    plink_out = os.path.join(filt_dir, "gt")

    if not os.path.exists(filt_dir):
        os.makedirs(filt_dir)
    if not os.path.exists(stats_dir):
        os.makedirs(stats_dir)

    gt_stats_dir = os.path.join(data_dir, genome, "assoc", markers, "stats")
    if not os.path.exists(gt_stats_dir):
        os.makedirs(gt_stats_dir)

    print("<b>Preparing filters</b>")
    if data["filt_from"]:
        plink_in = os.path.join(
            data_dir, genome, "assoc", data["filt_from"], "gt", "gt"
        )
        stats_dir = os.path.join(data_dir, genome, "assoc", data["filt_from"], "stats")

    if data["sample_group"]:
        samp_inc = os.path.join(filt_dir, "samp_include.txt")
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
        shet_exc = os.path.join(filt_dir, "shet_exclude.txt")
        if os.path.exists(shet_exc):
            os.remove(shet_exc)
        df = pd.read_csv(
            os.path.join(stats_dir, "shet.txt.gz"),
            sep="\t",
            header=None,
            index_col=None,
            chunksize=100000,
        )
        for chk in df:
            chk.columns = pd.RangeIndex(chk.columns.size)
            chk[chk[1] > float(data["shet"])][0].to_csv(
                shet_exc, mode="a", index=False, header=False
            )

    vhet_exc = None
    indels_exc = None

    if data["vhet"]:
        vhet_exc = os.path.join(filt_dir, "vhet_exclude.txt")
        if os.path.exists(vhet_exc):
            os.remove(vhet_exc)
        df = pd.read_csv(
            os.path.join(stats_dir, "vhet.txt.gz"),
            sep="\t",
            header=None,
            index_col=None,
            chunksize=100000,
        )
        for chk in df:
            chk.columns = pd.RangeIndex(chk.columns.size)
            chk[chk[1] > float(data["vhet"])][0].to_csv(
                vhet_exc, mode="a", index=False, header=False
            )

    if (
        (data["min_size"] and data["min_size"] > 1) or
        (data["max_size"] and data["max_size"] > 1)):

        min_size = data["min_size"] or -np.inf
        max_size = data["max_size"] or np.inf            

        indels_exc = os.path.join(filt_dir, "indels_exclude.txt")
        df = pd.read_csv(
            plink_in + ".pvar.zst",
            sep="\t",
            header=None,
            index_col=None,
            skiprows=1,
            usecols=[2,3,4],
            chunksize=100000,
        )
        for chk in df:
            var_size = np.maximum(chk[3].str.len(), chk[4].str.len())
            (chk[
                (var_size < min_size) |
                (var_size > max_size)][2]
             .to_csv(
                 indels_exc, mode="a", index=False, header=False
            ))

    else:
        snps_only = data["max_size"] == 1

    exc = []
    if vhet_exc:
        exc.append(vhet_exc)
    if indels_exc:
        exc.append(indels_exc)

    filt_out = plink_out + ".filt"
    # Apply filters
    if any(
        data[f]
        for f in (
            "sample_group",
            "min_size",
            "max_size",
            "coord",
            "smiss",
            "vmiss",
            "shet",
            "vhet",
            "maf",
            "hwe",
        )
    ):
        cmd = [plink2_path, "--pfile", plink_in, "vzs"]
        cmd = cmd + ["--memory", f"{get_mem()}"]
        cmd = cmd + ["--remove", shet_exc] if data["shet"] else cmd
        cmd = cmd + ["--keep", samp_inc] if data["sample_group"] else cmd
        cmd = cmd + ["--exclude", *exc] if exc else cmd
        cmd = cmd + ["--snps-only"] if not indels_exc and snps_only else cmd
        cmd = cmd + ["--extract", "bed1", coord_inc] if data["coord"] else cmd
        cmd = cmd + ["--mind", str(data["smiss"])] if data["smiss"] else cmd
        cmd = cmd + ["--geno", str(data["vmiss"])] if data["vmiss"] else cmd
        cmd = cmd + ["--maf", str(data["maf"])] if data["maf"] else cmd
        cmd = cmd + ["--hwe", str(10 ** (-data["hwe"]))] if data["hwe"] else cmd
        cmd = cmd + (
            ["--make-just-pvar", "zs", "cols="]
            if data["imp"]
            else ["--make-pgen", "vzs", "pvar-cols=", "psam-cols="]
        )
        cmd = cmd + ["--export", "vcf", "bgz", "id-paste=iid"] if data["imp"] else cmd
        cmd = cmd + ["--allow-extra-chr", "--out", filt_out]
        print("<b>Applying filters</b>")
        subprocess.run(cmd, check=True)

    else:
        filt_out = plink_in

    # Impute missing genotypes
    if data["imp"]:
        imp_out = plink_out + ".imp"
        imp_vcf = f"{imp_out}.vcf.gz"

        vcf_filt = filt_out + ".vcf.gz"
        if not os.path.exists(vcf_filt):
            subprocess.run(
                [
                    plink2_path,
                    "--memory",
                    f"{get_mem()}",
                    "--pfile",
                    filt_out,
                    "vzs",
                    "--export",
                    "vcf",
                    "bgz",
                    "id-paste=iid",
                    "--allow-extra-chr",
                    "--out",
                    filt_out,
                ],
                check=True,
            )

        # Determine one-marker chromosomes (removed, cannot be imputed)
        chr_d = dict()
        chr_data = pd.read_csv(
            filt_out + ".pvar.zst",
            sep="\t",
            usecols=[0, 1],
            skiprows=1,
            header=None,
            index_col=None,
            chunksize=1000000,
        )
        for chk in chr_data:
            for ch, vs in chk.groupby(0):
                if vs.shape[0] == 1:
                    chr_d.setdefault(ch, []).append(vs.reset_index().loc[0, 1])
        if any(len(v) == 1 for v in chr_d.values()):
            chrs = []
            excl_path = os.path.join(filt_dir, "exclude.txt")
            with open(excl_path, "w") as excl:
                for v in chr_d.values():
                    if len(v) == 1:
                        chrs.append(v[0].split(":")[0])
                        _ = excl.write(v[0] + "\n")
            print(
                "<b>WARNING: these chromosomes have only a single marker and were removed:</b>\n"
                + ", ".join(chrs)
            )

        beagle_path = os.path.join(tool_dir, "beagle.jar")

        print("<b>Imputing missing genotypes</b>")
        cmd = [
            "java",
            "-jar",
            f"-Xms{get_mem()}m",
            f"-Xmx{get_mem()}m",
            beagle_path,
            f"gt={vcf_filt}",
            f"out={imp_out}",
        ]
        if any(len(v) == 1 for v in chr_d.values()):
            cmd.append(f"excludemarkers={excl_path}")
        subprocess.run(cmd, check=True)

        subprocess.run(
            [
                plink2_path,
                "--memory",
                f"{get_mem()}",
                "--vcf",
                imp_vcf,
                "--make-pgen",
                "vzs",
                "pvar-cols=",
                "psam-cols=",
                "--allow-extra-chr",
                "--out",
                imp_out,
            ],
            check=True,
        )
    else:
        imp_out = filt_out

    # Prune by LD
    if data["prune"]:
        print("<b>Calculating LD for pruning</b>")
        data["window"] = str(data["window"]) + "kb" if str(data["window"]) else "200kb"
        cmd = [
            plink2_path,
            "--memory",
            f"{get_mem()}",
            "--pfile",
            imp_out,
            "vzs",
            "--indep-pairwise",
            data["window"],
            str(data["prune"]),
            "--allow-extra-chr",
            "--out",
            plink_out,
        ]
        subprocess.run(cmd, check=True)
        prune_in = plink_out + ".prune.in"

        cmd = [
            plink2_path,
            "--memory",
            f"{get_mem()}",
            "--pfile",
            imp_out,
            "vzs",
            "--extract",
            prune_in,
            "--allow-extra-chr",
            "--out",
            plink_out,
            "--make-pgen",
            "vzs",
            "pvar-cols=",
            "psam-cols=",
        ]
        print("<b>Pruning markers by LD</b>")
        subprocess.run(cmd, check=True)
    else:
        if data["hap"]:
            for ext in (".pgen", ".psam", ".pvar.zst"):
                if not os.path.exists(plink_out + ext):
                    shutil.copy(imp_out + ext, plink_out + ext)
        else:
            for ext in (".pgen", ".psam", ".pvar.zst"):
                if not os.path.exists(plink_out + ext):
                    os.link(imp_out + ext, plink_out + ext)

    print("<b>Processing filtered data</b>")
    # Make necessary genotype format files
    cmd = [
        plink2_path,
        "--memory",
        f"{get_mem()}",
        "--pfile",
        plink_out,
        "vzs",
        "--make-bed",
        "--allow-extra-chr",
        "--out",
        plink_out,
    ]
    subprocess.run(cmd, check=True)

    if data["hap"]:
        print("<b>Constructing haplotype blocks</b>")
        # Make haplotype blocks
        if data["btype"] == "ld":
            subprocess.run(
                [
                    plink1_path,
                    "--memory",
                    f"{get_mem()}",
                    "--bfile",
                    plink_out,
                    "--allow-extra-chr",
                    "--blocks",
                    "no-pheno-req",
                    "no-small-max-span",
                    "--blocks-max-kb",
                    str(data["bsize"]),
                    "--blocks-min-maf",
                    "0",
                    "--out",
                    plink_out,
                ],
                check=True,
            )
        else:
            win_block(plink_out, data["btype"], data["bspan"], data["bstep"], task_id)

        # os.remove(full_vcf)
        # os.remove(plink_out + ".vcf.gz")
        # Genotype based on blocks
        maf = 0.05
        min_haps = 0
        block_genotype(
            plink_out,
            gt_stats_dir,
            maf,
            min_haps,
            json.loads(config["concat"]),
            gdata["chr_names"],
            task_id,
        )
        kin_in = plink_out + ".allhaps"

        subprocess.run(
            [
                plink2_path,
                "--memory",
                f"{get_mem()}",
                "--bfile",
                plink_out,
                "--allow-extra-chr",
                "--make-pgen",
                "vzs",
                "pvar-cols=",
                "psam-cols=",
                "--out",
                plink_out,
            ],
            check=True,
        )
    else:
        kin_in = plink_out

    # Calculate K and Q
    print("<b>Determining kinship and population structure</b>")
    cmd = [
        plink2_path,
        "--memory",
        f"{get_mem()}",
        "--bfile",
        kin_in,
        "--make-rel",
        "square",
        "--pca",
        "15",
        "--allow-extra-chr",
        "--out",
        plink_out,
    ]
    if not data["hap"]:
        cmd += ["--export", "vcf", "bgz"]
    subprocess.run(cmd, check=True)

    pysam.tabix_index(plink_out + ".vcf.gz", preset="vcf")

    rel_matrix = pd.read_csv(f"{plink_out}.rel", sep="\t", header=None, index_col=None)
    # Convert PCA eigenvectors to PVE
    sum_diag = sum(np.diag(rel_matrix))
    eigenval = pd.read_csv(f"{plink_out}.eigenval", header=None, index_col=None)
    eigenval = (eigenval / sum_diag * 100)[0].tolist()
    eigenvec = pd.read_csv(
        f"{plink_out}.eigenvec",
        sep="\t",
        header=0,
        index_col=0,
        usecols=list(range(1, 16)),
    )
    samples = eigenvec.index.tolist()
    eigenvec = [eigenvec[c].values.tolist() for c in eigenvec.columns]
    pca_path = os.path.join(gt_stats_dir, "pca")
    with open(pca_path, "w") as pca_out:
        json.dump([samples, eigenvec, eigenval], pca_out)

    # Cluster kinship matrix
    hamm_out = os.path.join(gt_stats_dir, "kin")

    subprocess.run(
        [
            plink1_path,
            "--memory",
            f"{get_mem()}",
            "--bfile",
            kin_in,
            "--distance",
            "square",
            "--allow-extra-chr",
            "--out",
            hamm_out,
        ],
        check=True,
    )

    hamm_matrix = np.loadtxt(f"{hamm_out}.dist", delimiter="\t")
    hamm_labels = pd.read_csv(f"{hamm_out}.dist.id", sep="\t", header=None)[1].tolist()
    kin_matrix = (
        pd.DataFrame(hamm_matrix)
        .set_axis(hamm_labels, axis=0)
        .set_axis(hamm_labels, axis=1)
    )

    tree, label_order = wpgma(hamm_matrix, hamm_labels)

    kin_matrix = kin_matrix.loc[label_order, label_order]
    with open(f"{hamm_out}", "w") as matr_out:
        json.dump([kin_matrix.values.tolist(), label_order], matr_out)

    print("<b>Saving summary</b>")
    # Create file describing the filtering parameters
    params_fn = os.path.join(gt_stats_dir, "params")

    # Get original number of markers and samples
    with open(os.path.join(stats_dir, "cnts"), "r") as cnts_in:
        samp_all, var_all = json.load(cnts_in)
    # Get filtered number of markers and samples
    bed = bed_reader.open_bed(f"{plink_out}.bed")
    samp_filt, var_filt = bed.iid_count, bed.sid_count

    with open(params_fn, "w") as params_out:
        params = [
            [
                "Sample missingness",
                "Variant missingness",
                "Sample heterozygosity",
                "Variant heterozygosity",
                "Variant MAF",
                "Variant -log10(HWE)",
                "Variant LD",
                "LD window",
                "Samples retained/original",
                "Markers retained/original",
                "Is a haplotype set",
                "From markers",
            ],
            [
                float(data["smiss"]) if data["smiss"] else "—",
                float(data["vmiss"]) if data["vmiss"] else "—",
                float(data["shet"]) if data["shet"] else "—",
                float(data["vhet"]) if data["vhet"] else "—",
                float(data["maf"]) if data["maf"] else "—",
                float(data["hwe"]) if data["hwe"] else "—",
                float(data["prune"]) if data["prune"] else "—",
                data["window"] if data["window"] else "—",
                f"{samp_filt}/{samp_all}",
                f"{var_filt}/{var_all}",
                "Yes" if data["hap"] else "No",
                data["filt_from"],
            ],
        ]
        if data["hap"]:
            params[0].append("Block type")
            if data["btype"] == "ld":
                params[0].append("Max size (kb)")
                params[1].append("LD-based")
                params[1].append(data["bsize"])
            else:
                params[0].append("Block span")
                params[0].append("Block step")
                params[1].append(
                    "SNP number window"
                    if data["btype"] == "snpwin"
                    else "Base-pair window"
                )
                params[1].append(data["bspan"])
                params[1].append(data["bstep"])

        json.dump(params, params_out)

    if data["hap"]:
        hap_fn = os.path.join(gt_stats_dir, "is_hap")
        open(hap_fn, "w").close()

    # Clean up
    for fn in ("kin.dist", "kin.dist.id", "kin.log", "kin.nosex"):
        fn_path = os.path.join(gt_stats_dir, fn)
        if os.path.exists(fn_path):
            os.remove(fn_path)

    # Clean up
    for fn in (
        "gt.blocks",
        "gt.blocks.det",
        "gt.eigenval",
        # "gt.vcf.gz",
        "gt.filt.log",
        "gt.filt.vcf.gz",
        "gt.filt.pvar.zst",
        "gt.imp.log",
        "gt.imp.pgen",
        "gt.imp.psam",
        "gt.imp.pvar.zst",
        "gt.imp.vcf.gz",
        "gt.log",
        "gt.allhaps.bed",
        "gt.allhaps.bim",
        "gt.allhaps.fam",
        "gt.miss.log",
        "gt.miss.mindrem.id",
        "gt.miss.vcf.gz",
        "gt.nosex",
        "gt.prune.in",
        "gt.prune.out",
        "shet_exclude.txt",
        "vhet_exclude.txt",
    ):
        fn_path = os.path.join(filt_dir, fn)
        if os.path.exists(fn_path):
            os.remove(fn_path)

    if data["hap"]:
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

        calc_chr_ld(genome, markers=data["name"])
    else:
        run_qc(genome, markers=data["name"])

    # Change chromosomes to numeric (required by BLINK)
    chr_rename = {ch: idx for idx, ch in enumerate(gdata["chr_names"], start=1)}
    (
        pd.read_csv(f"{plink_out}.bim", sep="\t", index_col=None, header=None)
        .set_axis(["CHR", "ID", "FAM", "POS", "REF", "ALT"], axis=1)
        .assign(CHR=lambda x: x["CHR"].astype(str).map(chr_rename))
        .to_csv(f"{plink_out}.bim", sep="\t", header=False, index=False)
    )

    init_path = os.path.join(gt_stats_dir, "init")
    open(init_path, "w").close()


if __name__ == "__main__":
    genome, params, task_id = sys.argv[1:4]
    params = json.loads(params)
    filt_qc(genome, params, task_id)
