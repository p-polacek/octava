#!/usr/bin/env python

import os
import json
import numpy as np
import pandas as pd

from modules.load_input import load_gdata, load_config

from modules.assoc_plot import mqplot

from modules.models import models as models_dict

models_dict = dict(models_dict)


def assoc_process(
    genome,
    assoc_dir,
    assoc_fn,
    mqplot_title,
    out_root,
    out_dir,
    ch_col,
    pos_col,
    sig_col,
    sep="\t",
    eff_col=None,
    sig_log=True,
    se_col=None,
    delete=True,
):
    gdata = load_gdata(genome)
    config = load_config(genome)
    assoc_out = os.path.join(out_dir, assoc_fn)
    processed_out = os.path.join(out_dir, out_root)
    columns = [c for c in (ch_col, pos_col, sig_col, eff_col) if c]  # , se_col

    # Output files
    mqplot_fn = f"{processed_out}.mqplot.png"
    out_feather = f"{processed_out}.feather"
    out_eff_feather = f"{processed_out}.eff.feather"
    out_gz = f"{processed_out}.tsv.gz"

    # Save raw results, rename numeric chromosomes back to original
    chr_rename = dict(enumerate(gdata["chr_names"], start=1))
    res = (
        pd.read_csv(assoc_out, sep=sep, header=0, index_col=None)
        .loc[:, columns]
        .assign(**{ch_col: lambda x: x[ch_col].map(chr_rename)})
    )

    # res = res[~pd.isna(res[ch_col])]
    res.to_csv(out_gz, sep="\t", index=False)
    if eff_col:
        res[eff_col] = res[eff_col].abs()

    if sig_log:
        res[sig_col] = -np.log10(res[sig_col])
    # Save results for MP visualisation
    hap_fn = os.path.join(assoc_dir, "stats", "is_hap")
    is_hap = os.path.exists(hap_fn)
    if is_hap:
        blockmap = pd.read_csv(
            os.path.join(assoc_dir, "gt", "gt.blockmap"),
            sep="\t",
            header=0,
            index_col=None,
        )
        blockmap["chr"] = blockmap["chr"].astype(str)
        res = pd.merge(
            res,
            blockmap,
            left_on=[ch_col, pos_col],
            right_on=["chr", "start"],
            how="left",
        ).drop(["start"], axis=1)
        res["idx"] = res.groupby([ch_col, pos_col]).cumcount() + 1
        c = res.columns.tolist()
        (
            res[c[2:3] + c[:2] + c[-2:]]
            .sort_values(by=sig_col, ascending=False)
            .head(5000)
            .reset_index(drop=True)
            .to_feather(out_feather)
        )
        c = res.columns.tolist()
        if eff_col:
            (
                res[[eff_col] + c[:2] + c[-2:]]
                .sort_values(by=eff_col, ascending=False)
                .head(5000)
                .reset_index(drop=True)
                .to_feather(out_eff_feather)
            )
    else:
        c = res.columns.tolist()
        (
            res[c[2:3] + c[:2]]
            .sort_values(by=sig_col, ascending=False)
            .head(5000)
            .reset_index(drop=True)
            .to_feather(out_feather)
        )
        if eff_col:
            (
                res[[eff_col] + c[:2]]
                .sort_values(by=eff_col, ascending=False)
                .head(5000)
                .reset_index(drop=True)
                .to_feather(out_eff_feather)
            )
    mqplot(
        genome,
        res,
        mqplot_fn,
        mqplot_title,
        json.loads(config["concat"]),
        ch_col=ch_col,
        pos_col=pos_col,
        sig_col=sig_col,
        qqplot=sig_col != "PIP",
        thresh=sig_col != "PIP",
        end_col="end" if is_hap else None,
    )

    if delete:
        os.remove(assoc_out)
