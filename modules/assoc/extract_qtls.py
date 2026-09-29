#!/usr/bin/env python

import os
import sys
import json
import bed_reader
import subprocess
import portalocker
import numpy as np
import pandas as pd
import statsmodels.api as sm
from glimix_core.lmm import LMM
from numpy_sugar.linalg import economic_qs_linear

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from app.dirs import tmp_dir, data_dir, plink2_path

from modules.load_input import load_config

from modules.assoc.kq import kq_correct, kq_process

from modules.models import models as models_dict

models_dict = dict(models_dict)


def extract_qtls(genome, params, task_id):
    runs = json.loads(params["qtl_runs"])

    # Constants
    window = int(params["window_ld"])
    thresh = float(params["plot_thr"])
    r2 = float(params["r2_ld"])
    r2_markers = params["markers_sgl"]
    curr_locus = params["curr_locus"]

    # Header
    qtls = [
        [
            "Marker set",
            "Model",
            "Trait",
            "Correction",
            "Narrow sense heritability",
            "Peak marker",
            "Peak significance",
            "QTL span at r^2=0.8",
            "QTL span at r^2=0.4",
            f"QTL span at r^2={r2}",
            "Minor allele frequency",
            "Effect size",
            "Percentage of variance explained",
        ]
    ]
    pd.DataFrame(columns=qtls[0]).to_csv("QTL_data.txt", sep="\t", index=False)

    # Fixed input
    config = load_config(genome)
    phe_path = os.path.join(data_dir, genome, config["phe"])
    with portalocker.Lock(phe_path, "r", flags=portalocker.LOCK_SH) as locked_f:
        phe_df = pd.read_csv(phe_path, sep="\t", header=0, index_col=0)

    prev_markers = None
    pvar = None
    psam = None

    num_runs = len(runs)
    for run_idx, run in enumerate(runs, start=1):
        model_markers, model, trait = run
        print(
            ", ".join((model_markers, models_dict[model], trait)) + f" ({run_idx}/{num_runs})"
        )

        model_root = os.path.join(data_dir, genome, "assoc", model_markers, "gt", "gt")
        if prev_markers != model_markers:
            psam = pd.read_csv(f"{model_root}.psam", index_col=0).index
            pvar = pd.read_csv(f"{model_root}.pvar.zst", sep="\t", index_col=0, usecols=[2]).index
        bed = bed_reader.open_bed(f"{model_root}.bed")

        prev_markers = model_markers

        res_root = os.path.join(data_dir, genome, "assoc", model_markers, "res", f"{model}.{trait}")
        r2_root = os.path.join(data_dir, genome, "assoc", r2_markers or model_markers)

        model_hap = os.path.join(data_dir, genome, "assoc", model_markers, "stats", "is_hap")
        r2_hap = os.path.join(r2_root, "stats", "is_hap")
        is_hap = os.path.exists(model_hap)
        if is_hap != os.path.exists(r2_hap):
            print("Invalid marker set for QTL analysis.")
            continue

        res = pd.read_feather(f"{res_root}.feather")
        # Extract only current QTLs
        if params["fmt"] == "coord" and curr_locus:
            loci = curr_locus.split()
            mask = pd.Series(False, index=res.index)

            for locus in loci:
                if ":" not in locus:
                    continue
                chrom, span = locus.replace(",", "").split(":")
                if "-" in span:
                    start, end = map(int, span.split("-"))
                else:
                    start = end = int(span)

                mask |= (
                  (res["Chr"] == chrom) &
                  (res["Pos"] >= start) &
                  (res["Pos"] <= end)
                )
            res = res[mask]

        sig_col, ch_col, pos_col = res.columns[:3]
        if is_hap:
            end_col, idx_col = res.columns[3:5]
        else:
            end_col, idx_col = None, None
        res = res[res[sig_col] >= thresh]
        if is_hap:
            res["id"] = (
                res[ch_col].astype(str)
                + ":"
                + res[pos_col].astype(str)
                + "-"
                + res[end_col].astype(str)
                + "_"
                + res[idx_col].astype(str)
            )
        else:
            res["id"] = res[ch_col].astype(str) + ":" + res[pos_col].astype(str)

        # Get correction parameters
        with open(f"{res_root}.params") as params_in:
            q_num, k, kq_from, cov = json.load(params_in)[1][2:6]
        k = (k == "Yes")
        corr_string = ""
        if int(q_num):
            corr_string += f"{q_num} PCs, "
        if k:
            corr_string += "K, "
        if k or int(q_num):
            corr_string += f"from {kq_from}"
        if cov:
            if k or int(q_num):
                corr_string += ", + "
            corr_string += "covariates"

        # Calculate trait heritability
        phe = phe_df.loc[psam, trait].dropna()
        samples_mask = np.isin(psam, phe.index)

        q_matr, k_matr = kq_process(genome, model_markers, phe.index, True, int(q_num), cov, data_dir)

        null_cov = sm.add_constant(q_matr)
        QS = economic_qs_linear(k_matr)
        lmm0 = LMM(phe, null_cov["const"], QS=QS)
        lmm0.fit(verbose=False)

        sigma_g2_null = lmm0.v0   # genetic variance
        sigma_e2_null = lmm0.v1   # residual/noise variance

        hrt = round(sigma_g2_null / (sigma_g2_null + sigma_e2_null), 3)

        if res.empty:
            qtls.append(
                [
                    model_markers,
                    models_dict[model],
                    trait,
                    corr_string,
                    hrt,
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                ]
            )
            (
                pd.DataFrame([qtls[-1]], columns=qtls[0]).to_csv(
                    "QTL_data.txt", sep="\t", mode="a", header=False, index=False
                )
            )

        for ch, data in res.groupby(ch_col):
            while not data.empty:
                peak = data.iloc[0]
                if is_hap:
                    peak_pos = (
                        f"{ch}:{peak[pos_col]:,}-{peak[end_col]:,}_{peak[idx_col]:,}"
                    )
                else:
                    peak_pos = f"{ch}:{peak[pos_col]:,}"
                data = data.iloc[1:]

                in_root = os.path.join(r2_root, "gt", "gt")
                out_root = os.path.join(tmp_dir, task_id)

                if is_hap:
                    start, end = peak[pos_col], peak[end_col]
                else:
                    start = end = peak[pos_col]

                subprocess.run(
                    [
                        plink2_path,
                        "--pfile",
                        in_root,
                        "vzs",
                        "--chr",
                        ch,
                        "--from-bp",
                        str(start - window * 1000),
                        "--to-bp",
                        str(end + window * 1000),
                        "--r2-unphased",
                        "cols=id,pos",
                        "--ld-snp",
                        peak["id"],
                        "--ld-window-kb",
                        str(window),
                        "--ld-window-r2",
                        str(min(r2, 0.4)),
                        # "--ld-window-r2", "0",
                        "--silent",
                        "--out",
                        out_root,
                    ],
                    check=True,
                )

                ld = pd.read_csv(
                    f"{out_root}.vcor",
                    sep="\t",
                    header=0,
                    index_col=None,
                    usecols=[2, 3, 4],
                )

                ld.loc[-1] = [peak[pos_col], peak["id"], 1]
                #ld.index = ld.index + 1
                ld = ld.sort_values(by="POS_B").reset_index()

                # QTL span calculation
                ld_r2 = ld[ld["UNPHASED_R2"] >= r2]
                if ld_r2.empty:
                    span_r2 = peak_pos
                else:
                    if is_hap:
                        span_r2 = (
                            ch
                            + ":"
                            + f"{ld_r2.iloc[0]['POS_B']:,}"
                            + "-"
                            + f"{int(ld_r2.iloc[-1]['ID_B'].split('-')[1].split('_')[0]):,}"
                        )
                    else:
                        span_r2 = (
                            ch
                            + ":"
                            + f"{ld_r2.iloc[0]['POS_B']:,}"
                            + "-"
                            + f"{ld_r2.iloc[-1]['POS_B']:,}"
                        )

                ld_r4 = ld[ld["UNPHASED_R2"] >= 0.4]
                if ld_r4.empty:
                    span_r4 = peak_pos
                else:
                    if is_hap:
                        span_r4 = (
                            ch
                            + ":"
                            + f"{ld_r4.iloc[0]['POS_B']:,}"
                            + "-"
                            + f"{int(ld_r4.iloc[-1]['ID_B'].split('-')[1].split('_')[0]):,}"
                        )
                    else:
                        span_r4 = (
                            ch
                            + ":"
                            + f"{ld_r4.iloc[0]['POS_B']:,}"
                            + "-"
                            + f"{ld_r4.iloc[-1]['POS_B']:,}"
                        )

                ld_r8 = ld[ld["UNPHASED_R2"] >= 0.8]
                if ld_r8.empty:
                    span_r8 = peak_pos
                else:
                    if is_hap:
                        span_r8 = (
                            ch
                            + ":"
                            + f"{ld_r8.iloc[0]['POS_B']:,}"
                            + "-"
                            + f"{int(ld_r8.iloc[-1]['ID_B'].split('-')[1].split('_')[0]):,}"
                        )
                    else:
                        span_r8 = (
                            ch
                            + ":"
                            + f"{ld_r8.iloc[0]['POS_B']:,}"
                            + "-"
                            + f"{ld_r8.iloc[-1]['POS_B']:,}"
                        )

                # MAF calculation
                marker_idx = pvar.get_loc(peak["id"])
                gts = bed.read(index=marker_idx).flatten()[samples_mask]
                gts = gts[~np.isnan(gts)]
                af = gts.sum() / (gts.size * 2)
                maf = round(min(af, 1 - af), 2)
                gts_matr = pd.DataFrame(gts, columns=[trait], index=phe.index)

                if model == "mtmlm":
                    effect = "-"
                    pve = "-"
                else:
                    # PVE and effect calculation
                    full_cov = sm.add_constant(pd.concat([gts_matr, q_matr], axis=1))

                    snp_idx = list(full_cov.columns).index(trait)
                    gt = full_cov.iloc[:, snp_idx]
                    var_y = float(phe.var())

                    if k:
                        lmm = LMM(phe, full_cov, QS=QS)
                        lmm.fit(verbose=False)

                        snp_beta = lmm.beta[snp_idx]

                        sigma_g2 = lmm.v0   # genetic variance
                        sigma_e2 = lmm.v1   # residual/noise variance

                        var_snp = float(snp_beta**2 * gt.var())

                        effect = round(snp_beta, 2)
                        pve = round((var_snp / var_y) * 100, 2)

                    else:
                        full_fit = sm.OLS(phe, full_cov, missing="drop").fit()
                        #null_fit = sm.OLS(phe, null_cov, missing="drop").fit()

                        effect = round(full_fit.params[trait].astype(float), 2)
                        #pve = round((full_fit.rsquared - null_fit.rsquared) * 100, 2)

                        snp_beta = full_fit.params[trait]
                        var_snp = float(snp_beta**2 * gt.var())
                        pve = round((var_snp / var_y) * 100, 2)

                qtls.append(
                    [
                        model_markers,
                        models_dict[model],
                        trait,
                        corr_string,
                        hrt,
                        peak_pos,
                        round(peak[sig_col], 2),
                        span_r8,
                        span_r4,
                        span_r2,
                        maf,
                        effect,
                        pve,
                    ]
                )

                data = data[~data["id"].isin(ld_r2["ID_B"])]
                pd.DataFrame([qtls[-1]], columns=qtls[0]).to_csv(
                    "QTL_data.txt", sep="\t", mode="a", header=False, index=False
                )

    qtl_df = pd.DataFrame(qtls[1:], columns=qtls[0])
    excel_out = os.path.join(tmp_dir, "QTL_data.xlsx")
    with pd.ExcelWriter(excel_out, engine="xlsxwriter") as writer:
        qtl_df.to_excel(writer, index=False, sheet_name="Significant QTLs")


if __name__ == "__main__":
    genome, params, task_id = sys.argv[1:4]
    params = json.loads(params)
    extract_qtls(genome, params, task_id)
