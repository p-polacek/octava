#!/usr/bin/env python

# based on https://github.com/PriyaLakr/Manhattan-plot

import numpy as np
import scipy.stats as st
from itertools import cycle
import plotly.graph_objects as go
from plotly.subplots import make_subplots

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

    fig = make_subplots(
        rows=1, cols=2,
        column_widths=[0.6, 0.4],
        horizontal_spacing=0.1
    )
    #colors = ["#014e60", "#47b023"]
    colors = ["#e27267", "#60cdcd"]
    color_cycle = cycle(colors)

    # Manhattan plot
    in_df = in_df.sort_values(by=pos_col)
    for ch, data in in_df.groupby(ch_col, sort=False):
        color = next(color_cycle)
        fig.add_trace(
            go.Scatter(
                x=data[pos_col],
                y=data[sig_col],
                mode="markers",
                marker=dict(color=color, size=4),
                showlegend=False  # disable legend for each chromosome
            ),
            row=1, col=1
        )
        if end_col:
            for y, x0, x1 in zip(data[sig_col], data[start_col], data[end_col]):
                fig.add_trace(
                    go.Scatter(
                        x=[x0, x1],
                        y=[y, y],
                        mode="lines",
                        line=dict(color=color, width=0.5),
                        showlegend=False
                    ),
                    row=1, col=1
                )

    if thresh:
        fig.add_hline(
            y=-np.log10(0.05 / in_df.shape[0]),
            line=dict(color="red", width=1),
            row=1, col=1
        )

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

    max_sig = in_df[sig_col].max()
    max_pos = in_df.loc[in_df.shape[0] - 1, pos_col]
    fig.update_xaxes(
        title_text="Chromosome",
        tickmode="array",
        tickvals=ticks,
        ticktext=labels,
        tickangle=-20,
        ticks="outside",
        ticklen=5,
        tickwidth=1,
        tickcolor="black",
        showline=True,
        linewidth=1,
        linecolor="black",
        mirror=False,
        range=[-max_pos * 0.02, max_pos * 1.02],
        row=1,
        col=1,
    )

    ylabel = (
        "Importance"
        if not qqplot
        else "PIP" if sig_col == "PIP" else r"-log<sub>10</sub>(p-value)"
    )
    fig.update_yaxes(
        title_text=ylabel,
        ticks="outside",
        ticklen=5,
        tickwidth=1,
        tickcolor="black",
        showgrid=True,
        gridcolor="#dbdbdb",
        gridwidth=1,
        showline=True,
        linewidth=1,
        linecolor="black",
        mirror=False,
        range=[-max_sig * 0.04, max_sig * 1.04],
        row=1,
        col=1,
    )

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

        fig.add_trace(
            go.Scatter(
                x=theor_quant,
                y=in_df[sig_col],
                mode="markers",
                marker=dict(color=colors[0], size=4),
                showlegend=False
            ),
            row=1, col=2
        )
        fig.add_trace(
            go.Scatter(
                x=[0, theor_quant[-1]],
                y=[0, theor_quant[-1]],
                mode="lines",
                line=dict(color="black", width=1),
                showlegend=False
            ),
            row=1, col=2
        )

        fig.add_trace(
            go.Scatter(x=theor_quant, y=y_band_hi, mode="lines",
                       line=dict(width=0), showlegend=False),
            row=1, col=2
        )
        fig.add_trace(
            go.Scatter(x=theor_quant, y=y_band_lo, mode="lines",
                       fill="tonexty", line=dict(width=0),
                       #fillcolor="#ededed", showlegend=False),
                       fillcolor="rgba(81,180,194,0.15)", showlegend=False),
            row=1, col=2
        )
        
        fig.update_xaxes(
            title_text=r"Expected -log<sub>10</sub>(p-value)",
            showline=True,
            linewidth=1,
            linecolor="black",
            mirror=False,
            ticks="outside",
            ticklen=5,
            tickwidth=1,
            tickcolor="black",
            range=[-theor_quant[-1] * 0.04, theor_quant[-1] * 1.04],
            row=1,
            col=2
        )
        fig.update_yaxes(
            title_text=r"Observed -log<sub>10</sub>(p-value)",
            showline=True,
            linewidth=1,
            linecolor="black",
            mirror=False,
            ticks="outside",
            ticklen=5,
            tickwidth=1,
            tickcolor="black",
            range=[-max_y * 0.04, max_y * 1.04],
            row=1,
            col=2
        )
    else:
        fig.update_xaxes(visible=False, row=1, col=2)
        fig.update_yaxes(visible=False, row=1, col=2)

    fig.update_layout(
        title=dict(
            text=plot_title,
            x=0.5,
            xanchor="center",
            yanchor="top"
        ),
        plot_bgcolor="white",
        height=400, width=1000,
        margin=dict(l=60, r=20, t=80, b=50)
    )

    fig.write_image(out_fn, format="png", scale=3)
