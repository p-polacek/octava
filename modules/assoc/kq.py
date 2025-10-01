#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
from glimix_core.lmm import LMM
from numpy_sugar.linalg import economic_qs_linear

#from app.dirs import data_dir


def kq_process(genome, kq_from, samples, k, q_num, cov, fn_data_dir):
    kq_root = os.path.join(fn_data_dir, genome, "assoc", kq_from, "gt", "gt")

    q_matr = pd.DataFrame()
    k_matr = pd.DataFrame()
    # Construct fixed effect covariates
    q_list = list()
    if q_num > 0:
        # Covariates file filtered to keep top n selected ones
        # and necessary samples
        q_matr = pd.read_csv(
            f"{kq_root}.eigenvec", sep="\t", header=0, index_col=[0, 1]
        ).droplevel(0)
        q_matr = q_matr.loc[:, q_matr.columns.tolist()[:q_num]]
        if samples is not None and not samples.isin(q_matr.index).all():
            raise OSError(
                "<b>Dataset selected for K/Q correction does not contain "
                + "all samples from the genotype dataset.<b>"
            )
        q_list.append(q_matr)

    if cov:
        for cov_fn in cov:
            cov = pd.read_csv(
                os.path.join(fn_data_dir, genome, "covariates", cov_fn),
                sep="\t",
                header=0,
                index_col=0,
            )
            if samples is not None and not samples.isin(cov.index).all():
                raise OSError(
                    "<b>Additional covariates file does not contain "
                    + "all samples from the genotype dataset.<b>"
                )
            q_list.append(cov)

    if len(q_list) > 0:
        q_matr = pd.concat(q_list, axis=1)

    # Construct kinship matrix
    if k:
        k_matr = pd.read_csv(f"{kq_root}.rel", sep="\t", index_col=None, header=None)
        k_lab = pd.read_csv(
            f"{kq_root}.rel.id",
            sep="\t",
            usecols=[1],
            skiprows=[0],
            header=None,
            index_col=None,
        )[1]
        if samples is not None and not samples.isin(k_lab).all():
            raise OSError(
                "<b>Dataset selected for K/Q correction does not contain "
                + "all samples from the genotype dataset.<b>"
            )
        k_matr = k_matr.set_axis(k_lab, axis=0).set_axis(k_lab, axis=1)

    if samples is not None:
        if not q_matr.empty:
            q_matr = q_matr.loc[samples]
        if not k_matr.empty:
            k_matr = k_matr.loc[samples, samples]

    return q_matr, k_matr


def kq_correct(genome, raw_phe, samples, k, q_num, kq_from, cov, fn_data_dir):
    q_matr, k_matr = kq_process(genome, kq_from, samples, k, q_num, cov, fn_data_dir)

    corr_samples = raw_phe.dropna().index
    if not q_matr.empty:
        corr_samples = corr_samples.intersection(q_matr.index)
    if not k_matr.empty:
        corr_samples = corr_samples.intersection(k_matr.index)

    if not q_matr.empty:
        phe_q = sm.add_constant(q_matr.loc[corr_samples])
    else:
        phe_q = pd.DataFrame(1, index=corr_samples, columns=["I"])

    if not k_matr.empty:
        phe_k = k_matr.loc[corr_samples, corr_samples]

    raw_phe = raw_phe.loc[corr_samples]
    if k:
        QS = economic_qs_linear(phe_k)
        lmm = LMM(raw_phe, phe_q, QS=QS)
        lmm.fit(verbose=False)

        beta = lmm.beta
        delta = lmm.delta
        sigma_g2 = lmm.v0
        #sigma_e2 = lmm.v1

        # Compute V = sigma_g^2 * (K + delta * I)
        V = sigma_g2 * (phe_k + delta * np.eye(len(corr_samples)))

        # Solve V_inv (y - Xβ)
        resid_fixed = raw_phe - phe_q @ beta
        V_inv_resid = np.linalg.solve(V, resid_fixed)

        # Genetic component (BLUP): g = σ_g^2 * K * V⁻¹ (y - Xβ)
        g_hat = sigma_g2 * phe_k @ V_inv_resid

        # Fitted values: Xβ + g_hat
        y_fitted = phe_q @ beta + g_hat

        # Residuals: y - fitted
        phe_corr = raw_phe - y_fitted

    else:
        phe_corr = sm.OLS(raw_phe, phe_q).fit().resid

    #    phe_corr = raw_phe - lmm.mean()
    return phe_corr
