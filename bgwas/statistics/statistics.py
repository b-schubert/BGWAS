"""
Helper functions for statistical
analysis

:author: Benjamin Schubert
"""
from logging import warning


import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy import linalg

from statsmodels.sandbox.stats.multicomp import multipletests

import bgwas.contrib.fastlmm as fastlmm
import bgwas.contrib.pylmm as pylmm


def apc(df, i="i", j="j", colname="fn", outcol="cn"):
    """
    Average product correction of a co-evolutionary matrix


    :param df: pandas dataframe holding tuples of positions
    :param i: the column name of position i
    :param j: the coulmn name of position j
    :param colname: the column that should be corrected
    :param outcol: the column to which the corrected values are written
    :return:
    """
    col_means = df.groupby(by=[i])[colname].mean()
    row_means = df.groupby(by=[j])[colname].mean()
    matrix_mean = df[colname].mean()

    i = df[i]
    j = df[j]

    apc = (col_means.loc[i].values*row_means.loc[j].values) / matrix_mean
    df[outcol] = df[colname] - apc
    return df


def test_epistasis_fastlmm(pairs, X, y, mapping=None, alpha=0.05, n_grid_h2=100):
    """
    Calculates association statistics for a pair of SNPs encoded in the vector X of size n.
    We calculate the full rank model y ~ x1+ x2 + x1*x2
    and test with a likelihood ration test for epistatic association

    ASSERTION:
    (1) pairs-index can directly query X
    (2) X is standardized (necessary for fastLMM to work properly)

    :param pairs: list of pairs
    :param X: np.matrix of encoded genotypes (entry must match index pairs)
    :param y: np.array (transformed) phenotype
    :param mapping: a dictionary that maps from external genomic positions to internal positions to query X
    :param alpha: significance level for FDR correction (default: 0.05)
    :param n_grid_h2: number of grid parameter for heritability estimate
    :return: pandas dataframe with columns i,j,D,p_val,r,q_val
    """
    def infer(X, UX, UUX):
        lmm.X = X
        lmm.UX = UX
        lmm.UUX = UUX
        res = lmm.nLLeval(delta=internal_delta, REML=False)
        return -res["nLL"]

    res_dic = {"i": [], "j": [], "D": [], "p_val": [], "phi": []}
    if mapping is None:
        mapping = {}

    M = X.shape[1]
    # delete genotypic data with missing phenotypic data
    v = np.isnan(y)
    keep = True ^ v
    if v.sum():
        warning("Cleaning the phenotype vector by removing %d individuals...\n" % (v.sum()))
        y = y[keep]
        X = X[keep, :]

    # train heritage parameter
    lmm = fastlmm.LMM()
    print("Calculating K")
    lmm.setG(G0=X)

    X0 = np.ones(len(y)).reshape(len(y), 1)
    print("Setting parameters")
    lmm.sety(y)
    lmm.setX(X0)

    # sid_count == number of snps ??
    print("optimizing h")

    log_delta = lmm.find_log_delta(REML=False, sid_count=M, nGrid=n_grid_h2)['log_delta']
    internal_delta = np.exp(log_delta) * M

    for ii, jj in pairs:
        i = mapping.get(ii, ii)
        j = mapping.get(jj, jj)

        # update exclude field to prevent
        # proximal contamination (see Supplement Note 2: An Efficient Algorithm for Avoiding Proximal Contamination)
        # available at: http://www.nature.com/nmeth/journal/v9/n6/extref/nmeth.2037-S1.pdf
        # exclude SNPs from the RRM in the likelihood evaluation
        lmm.exclude_idx = [i, j]

        # rotate data
        X_e = np.column_stack((X0[:, 0], X[:, i], X[:, j], X[:, i] * X[:, j]))
        UX = lmm.U.T.dot(X_e)
        k = lmm.S.shape[0]
        N = X.shape[0]
        if k < N:
            UUX = X_e - lmm.U.dot(UX)
            UUX_null = UUX[:, :-1]
        else:
            UUX = None
            UUX_null = None

        # first calculate null model with no interaction term
        ll_null = infer(X_e[:, :-1], UX[:, :-1], UUX_null)

        # then calculate the the alternative
        ll_alt = infer(X_e, UX, UUX)

        # calculate chi-sqrt test with df=1
        D = -2. * (ll_null - ll_alt)
        ps = stats.chi2.sf(D, df=1)
        r = np.sqrt(D/X.shape[0])

        res_dic["i"].append(ii)
        res_dic["j"].append(jj)
        res_dic["D"].append(D)
        res_dic["p_val"].append(ps)
        res_dic["phi"].append(r)

    rejected, q_val, _, _ = multipletests(res_dic["p_val"], alpha=alpha, method='fdr_bh')
    res_dic["q_val"] = q_val
    return pd.DataFrame.from_dict(res_dic)[["i", "j", "D", "p_val", "phi", "q_val"]]


def test_epistasis_pylmm(pairs, X, y, mapping=None, alpha=0.05, n_grid_h2=100):
    """
    Calculates association statistics for a pair of SNPs encoded in the vector X of size n.
    We calculate the full rank model y ~ x1+ x2 + x1*x2
    and test with a likelihood ration test for epistatic association

    ASSERTION:
    (1) pairs-index can directly query X
    (2) X is standardized (necessary for fastLMM to work properly)

    :param pairs: list of pairs
    :param X: np.matrix of encoded genotypes (entry must match index pairs)
    :param y: np.array (transformed) phenotype
    :param mapping: a dictionary that maps from external genomic positions to internal positions to query X
    :param alpha: significance level for FDR correction (default: 0.05)
    :param n_grid_h2: number of grid parameter for heritability estimate
    :return: pandas dataframe with columns i,j,D,p_val,r,q_val
    """
    # fit the null model
    mapping = {} if mapping is None else mapping
    res_dic = {"i": [], "j": [], "D": [], "p_val": [], "phi": [], "beta":[]}
    K = pylmm.calculateKinship(X)
    v = np.isnan(y)
    keep = True - v
    if v.sum():
        print("Cleaning the phenotype vector by removing %d individuals...\n" % (v.sum()))
        y = y[keep]
        K = K[keep, :][:, keep]
        X = X[keep, :]

    Kva, Kve = linalg.eigh(K)

    L = pylmm.LMM(y, K, Kva, Kve, verbose=True)
    L.fit(ngrids=n_grid_h2, REML=False)

    for ii, jj in pairs:
        i = mapping.get(ii, ii)
        j = mapping.get(jj, jj)
        x1 = X[:, i]
        x2 = X[:, j]
        xs = np.column_stack((x1, x2, x1 * x2))
        D, ps, r = L.epistasis(xs)
        res_dic["i"].append(ii)
        res_dic["j"].append(jj)
        res_dic["D"].append(D)
        res_dic["p_val"].append(ps)
        res_dic["phi"].append(r)

    rejected, q_val, _, _ = multipletests(res_dic["p_val"], alpha=alpha, method='fdr_bh')
    res_dic["q_val"] = q_val
    return pd.DataFrame.from_dict(res_dic)[["i", "j", "D", "beta", "p_val", "phi", "q_val"]]


def test_single_pylmm(sites, X, y, mapping=None, alpha=0.05, n_grid_h2=100, reml=True, refit=False):
    """
    Calculates association statistics for a single SNP encoded in the vector X of size n using a linear mixed model
    and calculate a Wald-test on the coefficients.

    :param sites: indices to be tested
    :param X: Genome matrix n x m
    :param y: phenotype n x 1
    :param mapping: dict mapping from external (sites) indices to internal (mattrix index)
    :param alpha: significant threshold
    :param n_grid_h2: number of grid parameter for heritability estimate
    :param reml: Boolean indicating whether the ML or REML should be estimated
    :param refit: Boolean indicating whether variance terms are refitted for each site
    :return: pandas dataframe with columns i,D,p_val,r,q_val
    """
    # fit the null model
    mapping = {} if mapping is None else mapping
    res_dic = {"i": [],  "T": [], "p_val": [], "beta": []}
    K = pylmm.calculateKinship(X)
    v = np.isnan(y)
    keep = True - v
    if v.sum():
        print("Cleaning the phenotype vector by removing %d individuals...\n" % (v.sum()))
        y = y[keep]
        K = K[keep, :][:, keep]
        X = X[keep, :]

    Kva, Kve = linalg.eigh(K)

    L = pylmm.LMM(y, K, Kva, Kve, verbose=True)
    L.fit(ngrids=n_grid_h2, REML=reml)
    n = X.shape[0]

    for ii in sites:
        i = mapping.get(ii, ii)
        xs = X[:, i].reshape((n, 1))
        if refit:
            L.fit(X=xs, REML=reml)

        ts, ps, ds, beta, beta_var = L.association(xs, REML=reml, returnBeta=True)
        res_dic["i"].append(ii)
        res_dic["T"].append(ts)
        res_dic["p_val"].append(ps)
        res_dic["beta"].append(beta)

    rejected, q_val, _, _ = multipletests(res_dic["p_val"], alpha=alpha, method='fdr_bh')
    res_dic["q_val"] = q_val
    return pd.DataFrame.from_dict(res_dic)[["i", "T", "p_val", "beta", "q_val"]]