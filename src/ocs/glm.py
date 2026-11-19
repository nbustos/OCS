"""Mass-univariate GLMs on seed maps (grayordinates or ROI means): design, OLS t, FDR, permutation FWE.

One row per subject (between-subject maps at a timepoint, or post - pre change maps), so rows are
exchangeable and plain permutation / sign-flipping is valid.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def design(df: pd.DataFrame, test: str, covariates: list[str] = ()) -> tuple[np.ndarray, list[str], int, dict]:
    """Design matrix [N x p] with intercept, the tested term and covariates.

    Numeric columns are mean-centred (the intercept is then the mean at average covariate values).
    Categorical columns (object/category/bool) become dummies (first sorted level = reference).
    ``test='intercept'`` tests the mean (e.g. mean post - pre change). A categorical test must have
    exactly 2 levels; its column codes level[1] vs level[0].
    Returns X, column names, index of the tested column, info (levels / reference).
    """
    cols, names, info = [np.ones(len(df))], ["intercept"], {}
    test_idx = 0
    for term in ([test] if test != "intercept" else []) + [c for c in covariates if c != test]:
        v = df[term]
        if pd.api.types.is_bool_dtype(v) or not pd.api.types.is_numeric_dtype(v):
            levels = sorted(v.astype(str).unique())
            if term == test and len(levels) != 2:
                raise ValueError(f"test term {term!r} needs exactly 2 levels, has {levels}")
            for lev in levels[1:]:
                cols.append((v.astype(str) == lev).astype(float).to_numpy()); names.append(f"{term}[{lev}]")
            if term == test:
                test_idx = len(names) - 1
                info = {"reference": levels[0], "level": levels[1]}
        else:
            x = v.astype(float).to_numpy()
            cols.append(x - x.mean()); names.append(term)
            if term == test:
                test_idx = len(names) - 1
    X = np.column_stack(cols)
    if np.linalg.matrix_rank(X) < X.shape[1]:
        raise ValueError(f"design is rank deficient: {names}")
    return X, names, test_idx, info


def ols_t(Y: np.ndarray, X: np.ndarray, j: int) -> tuple[np.ndarray, np.ndarray, int]:
    """OLS of every column of Y [N x G] on X [N x p]; returns (beta_j, t_j, df). Columns with NaN -> NaN."""
    Y = np.asarray(Y, dtype=np.float64)
    n, p = X.shape
    dof = n - p
    if dof < 1:
        raise ValueError(f"not enough rows: N={n}, parameters={p}")
    ok = np.isfinite(Y).all(axis=0)
    beta = np.full(Y.shape[1], np.nan); t = np.full(Y.shape[1], np.nan)
    XtX_inv = np.linalg.inv(X.T @ X)
    B = XtX_inv @ X.T @ Y[:, ok]
    res = Y[:, ok] - X @ B
    sigma2 = (res ** 2).sum(axis=0) / dof
    se = np.sqrt(sigma2 * XtX_inv[j, j])
    beta[ok] = B[j]
    with np.errstate(divide="ignore", invalid="ignore"):
        t[ok] = B[j] / se
    return beta, t, dof


def p_two_sided(t: np.ndarray, dof: int) -> np.ndarray:
    return 2 * stats.t.sf(np.abs(t), dof)


def fdr_bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values (NaN kept)."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(pv); out[order] = np.minimum(ranked, 1)
    q[ok] = out
    return q


def perm_fwe(Y: np.ndarray, X: np.ndarray, j: int, t_obs: np.ndarray, n_perm: int,
             seed: int = 0) -> tuple[np.ndarray, int]:
    """Max-|t| family-wise p-values by Freedman-Lane permutation of reduced-model residuals.

    Testing the intercept (one-sample, e.g. post - pre change) uses sign-flipping instead of row
    permutation. Returns (p_fwe [G], number of permutations actually used). Exhaustive when the
    number of distinct permutations is <= n_perm.
    """
    Y = np.asarray(Y, dtype=np.float64)
    ok = np.isfinite(Y).all(axis=0) & np.isfinite(t_obs)
    Yo = Y[:, ok]
    n = len(Y)
    Z = np.delete(X, j, axis=1)
    if Z.shape[1]:
        Bz = np.linalg.lstsq(Z, Yo, rcond=None)[0]
        fit_z, res_z = Z @ Bz, Yo - Z @ Bz
    else:
        fit_z, res_z = np.zeros_like(Yo), Yo
    signflip = bool(np.allclose(X[:, j], 1))
    rng = np.random.default_rng(seed)
    if signflip:
        n_all = 2 ** n
        draws = (np.array([[1 if (k >> i) & 1 else -1 for i in range(n)] for k in range(n_all)])
                 if n_all <= n_perm else rng.choice([-1, 1], size=(n_perm, n)))
    else:
        from math import factorial
        from itertools import permutations
        n_all = factorial(n)
        draws = (np.array(list(permutations(range(n)))) if n_all <= n_perm
                 else np.array([rng.permutation(n) for _ in range(n_perm)]))
    XtX_inv = np.linalg.inv(X.T @ X); dof = n - X.shape[1]
    max_t = np.empty(len(draws))
    for k, d in enumerate(draws):
        Ys = fit_z + (res_z * d[:, None] if signflip else res_z[d])
        B = XtX_inv @ X.T @ Ys
        res = Ys - X @ B
        se = np.sqrt((res ** 2).sum(axis=0) / dof * XtX_inv[j, j])
        with np.errstate(divide="ignore", invalid="ignore"):
            max_t[k] = np.nanmax(np.abs(B[j] / se))
    p = np.full(len(t_obs), np.nan)
    p[ok] = (1 + (max_t[None, :] >= np.abs(t_obs[ok])[:, None] - 1e-12).sum(axis=1)) / (1 + len(draws))
    return p, len(draws)


def fit(Y: np.ndarray, df: pd.DataFrame, test: str, covariates: list[str] = (),
        n_perm: int = 0, seed: int = 0) -> dict:
    """Fit one model to Y [N x G] (rows aligned with df). Returns maps and metadata."""
    X, names, j, info = design(df, test, list(covariates))
    beta, t, dof = ols_t(Y, X, j)
    p = p_two_sided(t, dof)
    out = dict(beta=beta, t=t, p=p, q=fdr_bh(p), dof=dof, n=len(df), columns=names, tested=names[j], **info)
    if n_perm:
        out["p_fwe"], out["n_perm"] = perm_fwe(Y, X, j, t, n_perm, seed)
    return out
