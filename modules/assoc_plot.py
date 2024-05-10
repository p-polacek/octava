#!/usr/bin/env python

# based on https://github.com/PriyaLakr/Manhattan-plot

import os
import sys
import numpy as np
import pandas as pd
from itertools import cycle
import matplotlib.pyplot as plt

from modules.load_input import load_gdata

def mplot(
  genome,
  in_df,
  chr_vis,
  ch_col="Chr",
  pos_col="Pos",
  pval_col="P.value"):
    gdata  = load_gdata(genome)
    chroms = gdata["chr_names"]
    sizes  = gdata["chr_sizes"]
    sizes  = {
      ch: siz 
      for ch, siz in zip(chroms, sizes)}
    offsets = gdata["chr_offsets"]
    offsets = {
      ch: off 
      for ch, off in zip(chroms, offsets)}

    in_df[ch_col] = in_df[ch_col].astype(str)
    in_df[pos_col] = in_df[pos_col] + in_df[ch_col].map(offsets)
    if chr_vis:
        in_df[ch_col] = in_df[ch_col].apply(
          lambda ch: (
            ch if ch in chr_vis else "Other"))

    fig = plt.figure(figsize=(6, 4))
    colors = cycle(["#014e60", "#47b023"])

    plt.grid(
      visible=True,
      zorder=-1,
      axis="y",
      color="#dbdbdb")
    for ch, data in in_df.groupby(
      ch_col,
      sort=False):
        color = next(colors)
        plt.scatter(
          x      = data[pos_col],
          y      = data[pval_col],
          color  = color,
          label  = ch,
          s      = 4,
          zorder = 10)

    plt.axhline(
      y      = -np.log10(0.05/in_df.shape[0]),
      color  = "red",
      lw     = 1,
      zorder = 0)

    if chr_vis:
        plt.xticks(
          ticks = list(range(len(chr_vis) + 1)),
          labels = chr_vis + ["Other"])        
    else:
        plt.xticks(
          ticks = list(range(len(chroms))),
          labels = chroms)
        #(sizes[ch] + offsets[ch]) - sizes[ch]//2
        #for ch in chroms])
    plt.xlabel("Chromosomes")
    plt.ylabel("$-log_{10}$(p) values")

    #plt.tight_layout()

    return fig

def qqplot(
  in_df,
  pval_col="P.value"):
    in_df = in_df.sort_values(
      by=pval_col)

    num_tests = in_df.shape[0]
    theor_quant = -np.log10(
                    np.linspace(
                      1,
                      1/num_tests,
                      num = num_tests))

    fig = plt.figure(figsize=(6, 4))
    plt.scatter(
      x     = theor_quant,
      y     = in_df[pval_col],
      color = "#014e60",
      s     = 4)
    plt.plot(
      [0, theor_quant[-1]],
      [0, theor_quant[-1]],
      color  = "#47b023",
      zorder = 0)

    plt.xlabel("Expected $-log_{10}$(p)")
    plt.ylabel("Observed $-log_{10}$(p)")
    plt.tight_layout()    

    return fig
