# based on https://github.com/PriyaLakr/Manhattan-plot

import os
import sys
import numpy as np
import pandas as pd
from itertools import cycle
import matplotlib.pyplot as plt
from modules.getQueryOptions import getGenomeData

def manhattanplot(genome, in_df, ch_col="Chr", pos_col="Pos", pval_col="P.value"):

    genome_data = getGenomeData(genome)
    chroms  = genome_data["Chromosomes"]
    spans   = {k: v["size"] for k, v in genome_data["Chromosome_data"].items()}
    offsets = {k: v["offset"] for k, v in genome_data["Chromosome_data"].items()}

    in_df[ch_col] = in_df[ch_col].astype(str)
    in_df[pos_col] = in_df[pos_col] + in_df[ch_col].map(offsets)

    fig = plt.figure(figsize=(6, 4))
    colors = cycle(["#014e60", "#47b023"]) # Wageningen colors #3a8f1d

    plt.grid(visible=True, zorder=-1, axis="y", color="#dbdbdb")
    for ch, data in in_df.groupby(ch_col, sort=False):
        color = next(colors)
        plt.scatter(x=data[pos_col], y=data[pval_col], color=color, s=4, label=ch, zorder=10)

    plt.axhline(y=-np.log10(0.05/in_df.shape[0]), color="red", lw=1, zorder=0)

    plt.xticks(ticks=list(range(len(chroms))),
               labels=[(spans[ch] + offsets[ch]) - spans[ch]//2 for ch in chroms])
    plt.xlabel("Chromosomes")
    plt.ylabel("$-log_{10}$(p) values")

    plt.tight_layout()

    return fig


def qqplot(in_df, pval_col="P.value"):

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