#!/usr/bin/env python

# based on https://github.com/PriyaLakr/Manhattan-plot

import numpy as np
import scipy.stats as st
from itertools import cycle
import matplotlib.pyplot as plt

from modules.load_input import load_gdata


def mqplot(
    genome,
    in_df,
    out_fn,
    plot_title,
    concat,
    ch_col="Chr",
    pos_col="Pos",
    sig_col="P.value",
    thresh=True,
    qqplot=True,
    end_col=None,
):
    gdata = load_gdata(genome)
    chroms = gdata["chr_names"]
    sizes = {ch: size for ch, size in zip(chroms, gdata["chr_sizes"])}
    offsets = {ch: off for ch, off in zip(chroms, gdata["chr_offsets"])}

    in_df[ch_col] = in_df[ch_col].astype(str)
    in_df[pos_col] = in_df[pos_col] + in_df[ch_col].map(offsets)

    if end_col:
        start_col = "__start__"
        in_df[start_col] = in_df.loc[:, pos_col]
        in_df[end_col] = in_df[end_col] + in_df[ch_col].map(offsets)
        in_df[pos_col] = in_df[start_col] + (in_df[end_col] - in_df[start_col]) / 2
        in_df["span"] = in_df[end_col] - in_df[start_col]

    if concat:
        in_df[ch_col] = in_df[ch_col].apply(
            lambda ch: (ch if ch not in concat else "Other")
        )

    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(10, 4),
        gridspec_kw={"width_ratios": [0.6, 0.4]}
    )

    #colors = ["#014e60", "#47b023"]
    colors = ["#e27267", "#60cdcd"]
    color_cycle = cycle(colors)

    # Manhattan plot
    ax1.grid(visible=True, zorder=-1, axis="y", color="#dbdbdb")

    in_df = in_df.sort_values(by=pos_col)
    for ch, data in in_df.groupby(ch_col, sort=False):
        color = next(color_cycle)
        ax1.scatter(
            x=data[pos_col], y=data[sig_col], color=color, label=ch, s=4, zorder=10
        )
        if end_col:
            ax1.hlines(
                y=data[sig_col],
                xmin=data[start_col],
                xmax=data[end_col],
                color=color,
                linewidth=0.5,
            )

    if thresh:
        if thresh is True:
            ax1.axhline(y=-np.log10(0.05 / in_df.shape[0]), color="red", lw=1, zorder=0)
        else:
            ax1.axhline(y=thresh, color="red", lw=1, zorder=0)

    ticks = [sizes[ch] // 2 + offsets[ch] for ch in chroms if ch not in concat]
    labels = [ch for ch in chroms if ch not in concat]

    if concat:
        ticks.append(
            (
                sum([sizes[ch] for ch in chroms if ch in concat]) // 2
                + sum([sizes[ch] for ch in chroms if ch not in concat])
            )
        )
        labels.append("Other")

    ax1.set_xticks(ticks=ticks, labels=labels, rotation=45, ha="right")
    ax1.set_xlabel("Chromosome")

    max_sig = in_df[sig_col].max()
    min_sig = in_df[sig_col].min()
    min_pos = in_df[pos_col].min()
    max_pos = in_df[pos_col].max()
    
    ylabel = (
        "Importance" if sig_col == "Imp" else
        "PIP" if sig_col == "PIP" else
        r"$-log_{10}$(p-value)"
    )
    ax1.set_ylabel(ylabel)
    ax1.set_xlim([-max_pos * 0.02 + min_pos, max_pos * 1.02])
    if not (np.isnan(max_sig) or max_sig == np.inf):
        ax1.set_ylim([-max_sig * 0.04 + min_sig, max_sig * 1.04])

    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # QQ plot
    if qqplot:
        in_df = in_df.sort_values(by=sig_col).dropna()
        num_tests = in_df.shape[0]
        theor_quant = -np.log10(np.linspace(1, 1 / num_tests, num=num_tests))

        i = np.arange(1, num_tests+1)
        a = num_tests - i + 1   # note the flip vs. Beta(i, n-i+1)

        eps = np.finfo(float).tiny  # guard against 0/1
        p_lo = np.clip(st.beta.ppf(0.025, a, i), eps, 1 - eps)
        p_hi = np.clip(st.beta.ppf(0.975, a, i), eps, 1 - eps)

        # After -log10, the ordering flips: smaller p -> larger y
        y_band_hi = -np.log10(p_lo)  # this is the *upper* envelope in y
        y_band_lo = -np.log10(p_hi)  # this is the *lower* envelope in y
        max_y = max(max(y_band_hi), max_sig)

        ax2.scatter(x=theor_quant, y=in_df[sig_col], color=colors[0], s=4)
        ax2.plot([0, theor_quant[-1]], [0, theor_quant[-1]], color=colors[1], zorder=0)

        ax2.fill_between(
            theor_quant,
            y_band_lo,
            y_band_hi,
            color="#51b4c226",
            edgecolor="none",
            zorder=0
        )

        ax2.set_xlabel("Expected $-log_{10}$(p-value)")
        ax2.set_ylabel("Observed $-log_{10}$(p-value)")

        ax2.set_xlim([-theor_quant[-1] * 0.04, theor_quant[-1] * 1.04])
        if not (np.isnan(max_sig) or max_sig == np.inf):
            ax2.set_ylim([-max_y * 0.04 + min_sig, max_y * 1.04])

        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)

    else:
        ax2.set_visible(False)

    plt.suptitle(plot_title)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.2)

    fig.savefig(out_fn, dpi=300, format="png")
    plt.close(fig)
