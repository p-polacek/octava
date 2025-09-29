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

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from app.dirs import tmp_dir, data_dir, plink2_path

from modules.load_input import load_config

from modules.assoc.kq import kq_correct

from modules.models import models as models_dict

models_dict = dict(models_dict)


def extract_qtls(genome, params, task_id):
    runs = json.loads(params["qtl_runs"])
    window = int(params["window_ld"])
    thresh = float(params["plot_thr"])
    r2 = float(params["r2_ld"])

    qtls = [
        [
            "Marker set",
            "Model",
            "Trait",
            "Correction",
            #"Trait heritability",
            "Peak marker",
            "Peak significance",
            "QTL span at r^2=0.8",
            "QTL span at r^2=0.4",
            f"QTL span at r^2={r2}",
            "Minor allele frequency",
            "Corrected effect size",
            "Corrected percentage of variance explained",
            "Uncorrected effect size",
            "Uncorrected percentage of variance explained",
        ]
    ]

    num_runs = len(runs)
    pd.DataFrame(columns=qtls[0]).to_csv("QTL_data.txt", sep="\t", index=False)
    for run_idx, run in enumerate(runs, start=1):
        markers, model, trait = run
        print(
            ", ".join((markers, models_dict[model], trait)) + f" ({run_idx}/{num_runs})"
        )

        markers_root = os.path.join(data_dir, genome, "assoc", markers)
        hap_fn = os.path.join(markers_root, "stats", "is_hap")
        is_hap = os.path.exists(hap_fn)

        res_fn = os.path.join(markers_root, "res", f"{model}.{trait}.feather")
        res = pd.read_feather(res_fn)
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
        params_fn = os.path.join(markers_root, "res", f"{model}.{trait}.params")
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
            if k or int(q_num):
                corr_string += ", + "
            corr_string += "covariates"

        # Calculate trait heritability
        # heritability = ""

        if res.empty:
            qtls.append(
                [
                    markers,
                    models_dict[model],
                    trait,
                    corr_string,
                    # heritability,
                    "",
                    "",
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

                in_root = os.path.join(markers_root, "gt", "gt")
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

                # ld.loc[-1] = [peak["ps"], peak["id"], 1]
                # ld.index = ld.index + 1
                # ld = ld.sort_index()

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
                bed = bed_reader.open_bed(f"{in_root}.bed")
                config = load_config(genome)
                phe_path = os.path.join(data_dir, genome, config["phe"])
                with portalocker.Lock(
                    phe_path, "r", flags=portalocker.LOCK_SH
                ) as locked_f:
                    phe = (
                        pd.read_csv(phe_path, sep="\t", header=0, index_col=0)
                        .loc[bed.iid, trait]
                        .dropna()
                    )

                samples_mask = np.isin(bed.iid, phe.index)
                marker_idx = np.where(bed.sid == peak["id"])[0]
                if marker_idx.size != 1 or samples_mask.sum() != len(phe.index):
                    pve = "-"
                else:
                    gts = bed.read(index=marker_idx).flatten()

                    filt_gts = gts[samples_mask]
                    filt_gts = filt_gts[~np.isnan(filt_gts)]
                    af = filt_gts.sum() / (filt_gts.size * 2)
                    maf = round(min(af, 1 - af), 2)
                    gts = gts[samples_mask]

                if model == "mtmlm":
                    effect = "-"
                    pve = "-"
                else:
                    # PVE and effect calculation
                    phe_corr = kq_correct(
                        genome, phe, phe.index, False, int(q_num), markers, cov, data_dir  # k,
                    )
                    phe_fit = sm.OLS(
                        phe_corr, sm.add_constant(gts), missing="drop"
                    ).fit()
                    effect = round(phe_fit.params["x1"].astype(float), 2)
                    pve = round(phe_fit.rsquared * 100, 2)

                    phe_fit = sm.OLS(phe, sm.add_constant(gts), missing="drop").fit()
                    effect_u = round(phe_fit.params["x1"].astype(float), 2)
                    pve_u = round(phe_fit.rsquared * 100, 2)

                qtls.append(
                    [
                        markers,
                        models_dict[model],
                        trait,
                        corr_string,
                        # heritability,
                        peak_pos,
                        round(peak[sig_col], 2),
                        span_r8,
                        span_r4,
                        span_r2,
                        maf,
                        effect,
                        pve,
                        effect_u,
                        pve_u,
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
