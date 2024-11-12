#!/usr/bin/env python

# based on https://github.com/PriyaLakr/Manhattan-plot

import os
import sys
import numpy as np
import pandas as pd
from itertools import cycle
import matplotlib.pyplot as plt

from modules.load_input import load_gdata

def mqplot(
  genome,
  in_df,
  out_fn,
  concat,
  ch_col="Chr",
  pos_col="Pos",
  pval_col="P.value",
  thresh=True,
  qqplot=True):
    gdata   = load_gdata(genome)
    chroms  = gdata["chr_names"]
    sizes   = {ch: size for ch, size in zip(chroms, gdata["chr_sizes"])}
    offsets = {ch: off for ch, off in zip(chroms, gdata["chr_offsets"])}

    in_df[ch_col] = in_df[ch_col].astype(str)
    in_df[pos_col] = in_df[pos_col] + in_df[ch_col].map(offsets)
    if concat:
        in_df[ch_col] = in_df[ch_col].apply(
          lambda ch: (
            ch if ch not in concat else "Other"))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    colors = cycle(["#014e60", "#47b023"])

    # Manhattan plot
    ax1.grid(
      visible=True,
      zorder=-1,
      axis="y",
      color="#dbdbdb")
    for ch, data in in_df.groupby(
      ch_col,
      sort=False):
        color = next(colors)
        ax1.scatter(
          x      = data[pos_col],
          y      = data[pval_col],
          color  = color,
          label  = ch,
          s      = 4,
          zorder = 10)

    if thresh:
        ax1.axhline(
          y      = -np.log10(0.05/in_df.shape[0]),
          color  = "red",
          lw     = 1,
          zorder = 0)

    ticks = [
      sizes[ch]//2 + offsets[ch]
      for ch in chroms if ch not in concat]
    labels = [ch for ch in chroms if ch not in concat]

    if concat:
        ticks.append((
          sum([sizes[ch] for ch in chroms if ch in concat]) // 2 +
          sum([sizes[ch] for ch in chroms if ch not in concat])
          ))
        labels.append("Other")

    ax1.set_xticks(ticks=ticks, labels=labels, rotation=45, ha="right")
    ax1.set_xlabel("Chromosome")
    ax1.set_ylabel("$-log_{10}$(p) values")

    # QQ plot
    if qqplot:
        in_df = in_df.sort_values(
          by=pval_col)
        num_tests = in_df.shape[0]
        theor_quant = -np.log10(
                        np.linspace(
                          1,
                          1/num_tests,
                          num = num_tests))
        ax2.scatter(
          x     = theor_quant,
          y     = in_df[pval_col],
          color = "#014e60",
          s     = 4)
        ax2.plot(
          [0, theor_quant[-1]],
          [0, theor_quant[-1]],
          color  = "#47b023",
          zorder = 0)
        ax2.set_xlabel("Expected $-log_{10}$(p)")
        ax2.set_ylabel("Observed $-log_{10}$(p)")

    plt.tight_layout()

    fig.savefig(
      out_fn,
      dpi=300,
      format="png")
