import os
import re
import json
import subprocess
import numpy as np
import pandas as pd
#import plotly.express as px
from itertools import product

from modules.getQueryOptions import getFile
from modules.pl_manhattan import manhattanplot, qqplot

# Define working directories
#root_dir = os.getcwd()
root_dir = "./"
data_dir = os.path.join(root_dir, "data")
tool_dir = os.path.join(root_dir,
                        "modules",
                        "tools")
R_dir = os.path.join(root_dir,
                     "modules",
                     "R")


def run_QC(genome):
    vcf_file = getFile(genome, "vcf")

    plink_path = os.path.join(tool_dir,
                              "plink2")

    stats_dir = os.path.join(data_dir,
                             genome,
                             "stats")

    if not os.path.exists(stats_dir):
        os.makedirs(stats_dir)

    plink_out = os.path.join(stats_dir,
                             "plink_out")

    if any((not os.path.exists(f"{plink_out}.{ext}")
            for ext in ("bed", "bim", "fam"))):
        # Create binary PLINK files
        subprocess.run([plink_path,
                        "--vcf", vcf_file,
                        "--make-bed",
                        "--set-all-var-ids", "@:#",
                        "--allow-extra-chr",
                        "--out", plink_out],
                       check=True)

    if any((not os.path.exists(f"{plink_out}.{ext}")
            for ext in ("smiss", "vmiss", "het", "gcount", "afreq", "hardy"))):
        # Compute summary statistics:
        # missingness per marker and sample,
        # heterozygosity per marker and sample,
        # minor allele frequency,
        # Hardy-Weinberg equilibrium
        subprocess.run([plink_path,
                        "--bfile", plink_out,
                        "--missing",
                        "--geno-counts",
                        "--het",
                        "--freq",
                        "--hardy",
                        "--allow-extra-chr",
                        "--out", plink_out],
                       check=True)

    metrics = ("smiss", "vmiss", "shet", "vhet", "maf", "hwe")
    exts = ("smiss", "vmiss", "het", "gcount", "afreq", "hardy")
    id_cols = ("IID", "ID", "IID", "ID", "ID", "ID")

    for metric, ext, id_col in zip(metrics, exts, id_cols):
        qc_in = f"{plink_out}.{ext}"
        qc_df = pd.read_csv(qc_in, sep="\t")

        if metric in ("smiss", "vmiss"):
            qc_df["MIS_F"] = qc_df["MISSING_CT"] / qc_df["OBS_CT"]
            data = qc_df[[id_col, "MIS_F"]].sort_values(by="MIS_F", ascending=False)
        elif metric == "shet":
            qc_df["HET_F"] = (qc_df["OBS_CT"] - qc_df["O(HOM)"]) / qc_df["OBS_CT"]
            data = qc_df[[id_col, "HET_F"]].sort_values(by="HET_F", ascending=False)
        elif metric == "vhet":
            qc_df["HET_F"] = qc_df["HET_REF_ALT_CTS"] / (qc_df["HOM_REF_CT"] +
                                                         qc_df["HET_REF_ALT_CTS"] +
                                                         qc_df["TWO_ALT_GENO_CTS"])
            data = qc_df[[id_col, "HET_F"]].sort_values(by="HET_F", ascending=False)
        elif metric == "maf":
            qc_df["REF_FREQS"] = 1 - qc_df["ALT_FREQS"]
            qc_df["MAF"] = qc_df[["REF_FREQS", "ALT_FREQS"]].min(axis=1)
            data = qc_df[[id_col, "MAF"]].sort_values(by="MAF")
        elif metric == "hwe":
            qc_df["log10P"] = -np.log10(qc_df["P"])
            data = qc_df[[id_col, "log10P"]].sort_values(by="log10P", ascending=False)

        data.iloc[:, 1] = data.iloc[:, 1].round(5)
        data.to_csv(f"{plink_out}.{metric}.csv", index=False, header=False)

        heights, bins = np.histogram(data.iloc[:, 1])
        with open(f"{plink_out}.{metric}.hist.json", "w") as fig_out:
            json.dump([heights.tolist(), bins.tolist()], fig_out)

    return None


def filter_QC(genome, data, shet_exclude, vhet_exclude):
    vcf_file = getFile(genome, "vcf")

    plink_path = os.path.join(tool_dir,
                              "plink2")

    stats_dir  = os.path.join(data_dir,
                              genome,
                              "stats")

    GWAS_dir   = os.path.join(data_dir,
                              genome,
                              "GWAS",
                              data["name"])

    filt_dir   = os.path.join(GWAS_dir,
                              "imp_filt")

    plink_in =   os.path.join(stats_dir,
                              "plink_out")

    plink_out =  os.path.join(filt_dir,
                              "plink_out")

    if not os.path.exists(filt_dir):
        os.makedirs(filt_dir)
    if not os.path.exists(stats_dir):
        os.makedirs(stats_dir)

    # Create lists of samples and variants to filter out based on heterozygosity
    if shet_exclude:
        shet_out_fn = os.path.join(filt_dir, "shet_exclude.txt")
        with open(shet_out_fn, "w") as shet_out:
            for s in shet_exclude:
                _ = shet_out.write(s + "\n")
    if vhet_exclude:
        vhet_out_fn = os.path.join(filt_dir, "vhet_exclude.txt")
        with open(vhet_out_fn, "w") as vhet_out:
            for v in vhet_exclude:
                _ = vhet_out.write(v + "\n")

    # Prune by LD
    if data["prune"]:
        data["window"] = str(data["window"]) + "kb" if data["window"] else "200kb"

        subprocess.run([plink_path,
                        "--bfile", plink_in,
                        "--indep-pairwise", data["window"], str(data["prune"]),
                        "--allow-extra-chr",
                        "--out", plink_out],
                       check=True)

        prune_out_fn = plink_out + ".prune.out"
        with open(prune_out_fn, "r") as prune_out_in:
            prune_markers = prune_out_in.read().strip().split("\n")
            if vhet_exclude:
                prune_removed = len(set(prune_markers).difference(set(vhet_exclude)))
            else:
                prune_removed = len(prune_markers)

    else:
        prune_removed = ""

    # Filter based on thresholds:
    # missingness per marker and sample,
    # heterozygosity per marker and sample,
    # minor allele frequency,
    # Hardy-Weinberg equilibrium

    cmd = [plink_path,
           "--bfile", plink_in]

    cmd = cmd + ["--mind", str(data["smiss"])] if data["smiss"] else cmd
    cmd = cmd + ["--geno", str(data["vmiss"])] if data["vmiss"] else cmd

    cmd = cmd + ["--remove", shet_out_fn] if shet_exclude else cmd

    if vhet_exclude or data["prune"]:
        cmd += ["--exclude"]
        cmd = cmd + [vhet_out_fn] if vhet_exclude else cmd
        cmd = cmd + [prune_out_fn] if data["prune"] else cmd

    cmd = cmd + ["--maf", str(data["maf"])] if data["maf"] else cmd
    cmd = cmd + ["--hwe", str(10**(-data["hwe"]))] if data["hwe"] else cmd

    cmd = cmd + ["--allow-extra-chr",
                 "--export", "vcf", "bgz", "id-paste=iid",
                 "--out", plink_out]

    subprocess.run(cmd,
                   check=True)
    
    # Impute missing genotypes
    beagle_path = os.path.join(tool_dir,
                               "beagle.22Jul22.46e.jar")

    beagle_in = plink_out + ".vcf.gz"
    beagle_root = os.path.join(filt_dir,
                               "imp.filt")

    subprocess.run(["java", "-jar", beagle_path,
                    f"gt={beagle_in}",
                    f"out={beagle_root}"],
                   check=True)

    # Convert to tab-separated genotype file
    final_root = beagle_root
    subprocess.run([plink_path,
                    "--vcf", f"{beagle_root}.vcf.gz",
                    "--allow-extra-chr",
                    "--export", "A", "ped",
                    "--out", final_root],
                   check=True)

    os.rename(f"{final_root}.raw", f"{final_root}.varA.raw")

    with open(f"{final_root}.varA.raw", "r") as raw_in:
        with open(f"{final_root}.raw", "w") as raw_out:
            header = next(raw_in)
            header = [var.rsplit("_", 1)[0] for var in header.split("\t")]
            _ = raw_out.write("\t".join(header) + "\n")
            for line in raw_in:
                _ = raw_out.write(line)

    os.remove(f"{final_root}.ped")
    os.remove(f"{final_root}.varA.raw")

    # Create file describing the filtering parameters
    with open(f"{final_root}.params.json", "w") as params_out:
        with open(f"{plink_out}.log", "r") as log_in:
            log_lines = log_in.read().split(".\n")

            # Exclude header lines, keep only lines describing filters
            _ = log_lines.pop(0)

            samples_line = [line for line in log_lines 
                            if "samples" in line and "loaded" in line]
            samples_all = (int(samples_line[0].split(" ", 1)[0])
                           if len(samples_line) == 1 else "")

            vars_line = [line for line in log_lines 
                        if "variants" in line and "loaded" in line]
            vars_all = (int(vars_line[0].split(" ", 1)[0])
                       if len(vars_line) == 1 else "")

            samples_remaining = samples_all
            vars_remaining = vars_all
            if samples_all and vars_all:
                shet_remaining = samples_remaining - len(shet_exclude) if shet_exclude else samples_remaining
                vhet_remaining = vars_remaining - len(vhet_exclude) if vhet_exclude else vars_remaining
                samples_remaining, vars_remaining  = shet_remaining, vhet_remaining

                prune_remaining = vars_remaining - prune_removed if prune_removed != None else vars_remaining
                vars_remaining  = prune_remaining

                smiss_line = [line for line in log_lines 
                              if "--mind" in line]
                smiss_remaining = (samples_remaining - int(re.search("\d+", smiss_line[0])[0])
                                   if len(smiss_line) == 1 else "")
                samples_remaining = smiss_remaining or samples_remaining

                vmiss_line = [line for line in log_lines 
                              if "--geno" in line]
                vmiss_remaining = (vars_remaining - int(re.search("\d+", vmiss_line[0])[0])
                                   if len(vmiss_line) == 1 else "")
                vars_remaining = vmiss_remaining or vars_remaining

                hwe_line = [line for line in log_lines 
                            if "--hwe" in line]
                hwe_remaining = (vars_remaining - int(re.search("\d+", hwe_line[0])[0])
                               if len(hwe_line) == 1 else "")
                vars_remaining = hwe_remaining or vars_remaining

                maf_line = [line for line in log_lines 
                            if "--maf" in line]
                maf_remaining = (vars_remaining - int(re.search("\d+", maf_line[0])[0])
                               if len(maf_line) == 1 else "")
                vars_remaining = maf_remaining or vars_remaining

                samples_line = [line for line in log_lines 
                                if "samples" in line and "remaining" in line and "filters" in line]
                samples_remaining = (int(samples_line[0].split(" ", 1)[0])
                                     if len(samples_line) == 1 else "")

                vars_line = [line for line in log_lines 
                            if "variants" in line and "remaining" in line and "filters" in line]
                vars_remaining = (int(vars_line[0].split(" ", 1)[0])
                                  if len(vars_line) == 1 else "")

                params = list()
                params.append(["", "Initial sample count",
                               "Sample heterozygosity", "Sample missingness",
                                "Initial marker count", "Variant heterozygosity", "LD window", "Variant LD",
                               "Variant missingness", "Variant -log10(HWE)",
                               "Variant MAF", "Final sample count", "Final marker count"])

                params.append(["Threshold", "-",
                               float(data["shet"]) if data["shet"] else "-",
                               float(data["smiss"]) if data["smiss"] else "-",
                               "-",
                               float(data["vhet"]) if data["vhet"] else "-",
                               data["window"] if data["window"] else "-",
                               float(data["prune"]) if data["prune"] else "-",
                               float(data["vmiss"]) if data["vmiss"] else "-",
                               float(data["hwe"]) if data["hwe"] else "-",
                               float(data["maf"]) if data["maf"] else "-", "-", "-"])
                params.append(["Remaining", samples_all, shet_remaining, smiss_remaining,
                               vars_all, vhet_remaining, "-", prune_remaining, vmiss_remaining, 
                               hwe_remaining, maf_remaining, samples_remaining, vars_remaining])

                json.dump(params, params_out)

            else:
                json.dump([], params_out)


def run_GWAS(genome, params):
    markers = params["markers"]
    models = params["models"]
    traits = params["traits"][:]
    for trait_group in params["trait_groups"]:
        for trait in trait_group.split(","):
            if trait not in params["traits"]:
                traits.append(trait)

    gapit_path = os.path.join(R_dir,
                              "GAPIT.R")

    for marker_set in markers:
        GWAS_dir = os.path.join(data_dir,
                                genome,
                                "GWAS",
                                marker_set)

        out_dir = os.path.join(GWAS_dir,
                               "assoc")

        in_root = os.path.join(GWAS_dir,
                               "imp_filt",
                               "imp.filt")

        phe_file = getFile(genome, "phe")
        q_file = os.path.join(data_dir,
                              genome,
                              "plink2.top3.eigenvec.txt")

        """
        command = ["Rscript", gapit_path,
                   "-g", f"{in_root}.raw",
                   "-m", f"{in_root}.map",
                   "-p", phe_file,
                   "-q", q_file,
                   "--q_tsv",
                   "--phenotype_tsv",
                   "-o", out_dir,
                   "-l", ",".join(models),
                   "-t", ",".join(traits)]
        if params["chrom"]:
            command += ["-c", params["chrom"]]
        if params["start"]:
            command += ["-s", params["start"]]
        if params["end"]:
            command += ["-e", params["end"]]
                   
        subprocess.run(command, check=True)
        """
        raw_dir = os.path.join(out_dir, "raw")

        models = models if isinstance(models, list) else [models]
        traits = traits if isinstance(traits, list) else [traits]

        for model, trait in product(models, traits):
            out_root = f"GAPIT.processed.{model}.{trait}"
            if all([os.path.exists(os.path.join(out_dir,
                                                f"{out_root}.{ext}")
                                   )
                    for ext in ("tsv",
                                "tsv.gz",
                                "tsv.gz.tbi")]):
                continue

            in_fn = os.path.join(raw_dir,
                                 f"GAPIT.Association.GWAS_Results.{model}.{trait}.csv")
            out_tsv = os.path.join(out_dir,
                                   f"{out_root}.tsv")
            out_gz = os.path.join(out_dir,
                                  f"{out_root}.tsv.gz")

            res = pd.read_csv(in_fn)
            cols = ["Chr", "Pos", "P.value", "Effect"]
            res[cols].to_csv(out_tsv, sep="\t", index=False)

            subprocess.run(["bgzip", "-f", "-c", out_tsv],
                           stdout=open(out_gz, "w"),
                           check=True)

            subprocess.run(["tabix",
                            "-s", "1",
                            "-b", "2",
                            "-e", "2",
                            "--skip-lines", "1",
                            "-f",
                            out_gz],
                           check=True)

            mplot_fn = os.path.join(out_dir, f"{model}.{trait}.mplot.png")
            qqplot_fn = os.path.join(out_dir, f"{model}.{trait}.qqplot.png")

            if not os.path.exists(mplot_fn) or not os.path.exists(qqplot_fn):

                spans_fn = os.path.join("data", genome, "spans.txt")
                mplot_fig = manhattanplot(out_tsv, spans_fn)
                mplot_fig.savefig(mplot_fn, format="png")

                qqplot_fig = qqplot(out_tsv)
                qqplot_fig.savefig(qqplot_fn, format="png")

                df = pd.read_csv(out_tsv, sep="\t")
                df["P.value"] = -np.log10(df["P.value"]).round(5)
                df = df.sort_values(by="P.value", ascending=False)
                df["Effect"] = df["Effect"].round(5)
                df.to_csv(out_tsv, sep="\t", index=False)
