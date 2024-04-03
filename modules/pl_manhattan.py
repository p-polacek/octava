# based on https://github.com/PriyaLakr/Manhattan-plot

import os
import sys
import numpy as np
import pandas as pd
from itertools import cycle
import matplotlib.pyplot as plt

def read_file(in_fn):
    if in_fn.endswith(".csv"):
    	in_df = pd.read_csv(in_fn, sep=",")
    elif in_fn.endswith((".tsv", ".txt")):
        in_df = pd.read_csv(in_fn, sep="\t")
    else:
        return None
    return in_df


def manhattanplot(in_fn, spans_fn, ch_col="Chr", pos_col="Pos", pval_col="P.value"):
    in_df = read_file(in_fn)

    with open(spans_fn, "r") as spans_in:
        last_pos = 0
        spans = [line.split("\t")
                 for line in spans_in.read().strip().split("\n")]
        spans = [(ch, int(span)) for ch, span in spans]
        offsets = dict()
        for ch, span in spans:
            offsets[ch] = last_pos
            last_pos += span

    in_df[ch_col] = in_df[ch_col].astype(str)
    in_df[pos_col] = in_df[pos_col] + in_df[ch_col].map(offsets)
    in_df[pval_col] = -np.log10(in_df[pval_col])

    fig = plt.figure(figsize=(6, 4))
    colors = cycle(["#014e60", "#47b023"]) # Wageningen colors #3a8f1d

    plt.grid(visible=True, zorder=-1, axis="y", color="#dbdbdb")
    for ch, data in in_df.groupby(ch_col, sort=False):
        color = next(colors)
        plt.scatter(x=data[pos_col], y=data[pval_col], color=color, s=4, label=ch, zorder=10)

    plt.axhline(y=-np.log10(0.05/in_df.shape[0]), color="red", lw=1, zorder=0)

    plt.xticks(ticks=[(span + offsets[ch]) - span//2 for ch, span in spans],
               labels=list(range(len(spans))))
    plt.xlabel("Chromosomes")
    plt.ylabel("$-log_{10}$(p) values")

    plt.tight_layout()

    return fig


def qqplot(in_fn, pval_col="P.value"):
    in_df = read_file(in_fn)
    in_df[pval_col] = -np.log10(in_df[pval_col])
    in_df = in_df.sort_values(by=pval_col)

    num_tests = in_df.shape[0]
    theor_quant = -np.log10(np.linspace(1, 1/num_tests, num=num_tests))

    fig = plt.figure(figsize=(6, 4))
    plt.scatter(x=theor_quant, y=in_df[pval_col], s=4, color="#014e60")
    plt.plot([0, theor_quant[-1]], [0, theor_quant[-1]], color="#47b023", zorder=0)

    plt.xlabel("Expected $-log_{10}$(p)")
    plt.ylabel("Observed $-log_{10}$(p)")
    plt.tight_layout()    

    return fig


if __name__ == "__main__":
    pass