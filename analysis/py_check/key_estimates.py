#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
key_estimates.py -- independent second-programmer check (Python)
=================================================================

Purpose
-------
Protocol v1.0 (Section 5.6) requires that the key numbers of the NHANES III
analysis are programmed a second time, independently of the primary R code.
This script was written ONLY from the written specification and the protocol
(`analysis/prereg/protocol_v1.0_en.md`); the R code and any results were not seen.

What is computed (NHANES III)
-----------------------------
1. Part 1 ITT primary estimate: effect of randomized session Z (1 = afternoon/
   evening, 0 = morning) on ln(GLR), weighted by the MEC weight.
     - SE by Taylor linearization (stratified-cluster sandwich, as survey::svyglm)
     - SE by Fay BRR (52 replicate weights, rho = 0.3)
2. CACE (Wald estimator) = beta_ITT / (P(D=1|Z=1) - P(D=1|Z=0)), SE by BRR.
3. Part 2: weighted Cox proportional-hazards models (Breslow ties, written from
   scratch because lifelines is not available) for all-cause and cancer death,
   with the observed exposure log2(GLR) and the session-standardized exposure
   log2(GLR) - (CACE/ln 2) * D.  Delta = (beta_std - beta_obs) / beta_obs, with a
   paired BRR SE (the CACE is re-estimated inside every replicate).  The analytic
   attenuation prediction lambda is reported alongside.

Safety (blinding)
-----------------
* Default input is the SCRAMBLED file.  Any input whose file name does not
  contain "SCRAMBLED" (e.g. *_REAL*) is refused BEFORE it is opened, unless
  `analysis/prereg/OSF_REGISTRATION.txt` exists and contains "osf.io".
* All console output and the JSON are labelled "SCRAMBLED – NOT REAL" when the
  input is scrambled.

Usage
-----
    python3 key_estimates.py                    # scrambled data (default)
    python3 key_estimates.py --input FILE.csv.gz
    python3 key_estimates.py --skip-validation  # skip the statsmodels checks

Output: analysis/py_check/key_estimates_SCRAMBLED.json (or ..._REAL.json)
"""

import argparse
import datetime as _dt
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import linalg, stats
from scipy.interpolate import BSpline

# ---------------------------------------------------------------------------
# Paths and fixed constants
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent            # .../analysis/py_check
PROJECT_ROOT = SCRIPT_DIR.parent.parent                  # project root
DEFAULT_INPUT = PROJECT_ROOT / "analysis" / "derived" / "nhanes3_SCRAMBLED.csv.gz"
OSF_FILE = PROJECT_ROOT / "analysis" / "prereg" / "OSF_REGISTRATION.txt"

N_REP = 52                                   # Fay BRR replicate weights WTPXRP1-52
FAY_RHO = 0.3                                # Fay coefficient
BRR_DENOM = N_REP * (1.0 - FAY_RHO) ** 2     # 52 * 0.49
DF_DESIGN = 49                               # 98 PSUs - 49 strata
T975 = float(stats.t.ppf(0.975, DF_DESIGN))  # t quantile used for all CIs
COX_TOL = 1e-9                               # Newton-Raphson: stop when max|step| < 1e-9
LN2 = float(np.log(2.0))
BANNER = "SCRAMBLED – NOT REAL"
REP_COLS = [f"wtpxrp{r}" for r in range(1, N_REP + 1)]
SEED = 20260926                              # only used to draw integer test weights

# Factor levels used in Part 2.  The first level present in the data is the
# reference; the choice of reference does not affect the exposure coefficient.
FACTOR_LEVELS = {
    "sex": ["1", "2"],
    "race_eth": ["1", "2", "3", "4"],
    "phase": ["1", "2"],
    "smoke": ["never", "former", "current", "missing"],
    "bmi_missing": ["0", "1"],
    "educ_cat": ["<9", "9-11", "12", ">12", "missing"],
    "pir_cat": ["<1", "1-<2", "2-<4", ">=4", "missing"],
}


# ---------------------------------------------------------------------------
# 0. Input guard (runs before any data file is opened)
# ---------------------------------------------------------------------------
def osf_registered() -> bool:
    """True only if the OSF registration file exists and contains 'osf.io'."""
    try:
        return OSF_FILE.is_file() and ("osf.io" in OSF_FILE.read_text(encoding="utf-8", errors="ignore"))
    except OSError:
        return False


def guard_input(path: Path) -> bool:
    """Decide from the FILE NAME ONLY whether the input is scrambled.

    Returns True for a scrambled file.  For anything else (a *_REAL* file or any
    unrecognised name) the script stops unless the OSF registration exists.
    """
    name = path.name
    if not (name.endswith(".csv") or name.endswith(".csv.gz")):
        sys.exit(f"REFUSED: input must be a .csv or .csv.gz file, got '{name}'.")
    upper = name.upper()
    is_scrambled = ("SCRAMBLED" in upper) and ("REAL" not in upper)
    if is_scrambled:
        return True
    if not osf_registered():
        sys.exit(
            f"REFUSED: '{name}' is not a SCRAMBLED file. Real data may be analysed only after "
            f"OSF registration: '{OSF_FILE}' must exist and contain 'osf.io'. Nothing was read."
        )
    return False


def as_bool(s: pd.Series) -> pd.Series:
    """Parse a True/False column that may arrive as bool or as text."""
    if s.dtype == bool:
        return s
    m = s.astype(str).str.strip().str.lower().map(
        {"true": True, "false": False, "1": True, "0": False, "1.0": True, "0.0": False})
    if m.isna().any():
        raise ValueError("unparseable boolean values in in_itt")
    return m.astype(bool)


# ---------------------------------------------------------------------------
# 1. Survey helpers
# ---------------------------------------------------------------------------
def wls_taylor(y, X, w, strata, psu):
    """Weighted least squares + stratified-cluster sandwich variance.

    A = X'WX;  u_hi = sum over PSU i of stratum h of w_j x_j e_j (e = WLS residual);
    B = sum_h n_h/(n_h-1) sum_i (u_hi - ubar_h)(u_hi - ubar_h)';  V = A^-1 B A^-1.
    This is the linearization used by R survey::svyglm (no fpc).
    """
    Xw = X * w[:, None]
    A = X.T @ Xw
    beta = linalg.solve(A, Xw.T @ y)
    e = y - X @ beta
    u = Xw * e[:, None]
    tot = pd.DataFrame(u).groupby([np.asarray(strata), np.asarray(psu)]).sum()
    B = np.zeros((X.shape[1], X.shape[1]))
    n_psu = []
    for _, grp in tot.groupby(level=0):
        U = grp.to_numpy()
        nh = U.shape[0]
        if nh < 2:
            raise ValueError("a stratum has a single PSU (lonely PSU) - not expected in NHANES III")
        n_psu.append(nh)
        dev = U - U.mean(axis=0)
        B += nh / (nh - 1.0) * (dev.T @ dev)
    Ainv = linalg.inv(A)
    return beta, Ainv @ B @ Ainv, n_psu


def wmean_by_group(y, g, W):
    """Weighted means of y within g==1 and g==0 for one weight vector (n,) or
    a matrix of weights (n, R).  Returns (mean_g1, mean_g0), each scalar or (R,)."""
    W = np.asarray(W, float)
    g1 = (g == 1).astype(float)
    g0 = 1.0 - g1
    m1 = (g1 * y) @ W / (g1 @ W)
    m0 = (g0 * y) @ W / (g0 @ W)
    return m1, m0


def brr_var(theta_r, theta):
    """Fay BRR variance, centred at the full-sample estimate (as specified)."""
    theta_r = np.asarray(theta_r, float)
    return float(np.sum((theta_r - theta) ** 2) / BRR_DENOM)


def brr_var_mean_centred(theta_r):
    """Alternative (R survey default mse=FALSE): centred at the replicate mean."""
    theta_r = np.asarray(theta_r, float)
    return float(np.sum((theta_r - theta_r.mean()) ** 2) / BRR_DENOM)


def pct(b):
    return (np.exp(b) - 1.0) * 100.0


# ---------------------------------------------------------------------------
# 2. Natural cubic spline bases
# ---------------------------------------------------------------------------
def ns_esl(x, knots):
    """Natural cubic spline basis of Hastie, Tibshirani & Friedman (ESL, 2nd ed.,
    eqs. 5.4-5.5) for K knots xi_1 < ... < xi_K (boundary knots included):
        N_1 = 1, N_2 = x, N_{k+2} = d_k(x) - d_{K-1}(x), k = 1..K-2,
        d_k(x) = [(x - xi_k)_+^3 - (x - xi_K)_+^3] / (xi_K - xi_k).
    Returns N_2..N_K (K-1 columns); the constant N_1 is dropped because the Cox
    baseline hazard absorbs it.  Together with a constant this spans the same
    space as R splines::ns(x, knots = interior, Boundary.knots = boundary)."""
    x = np.asarray(x, float)
    xi = np.asarray(knots, float)
    K = xi.size
    if np.any(np.diff(xi) <= 0):
        raise ValueError(f"knots must be strictly increasing: {xi}")

    def d(k):  # 0-based k = 0..K-2
        return (np.clip(x - xi[k], 0.0, None) ** 3 - np.clip(x - xi[K - 1], 0.0, None) ** 3) / (xi[K - 1] - xi[k])

    d_last = d(K - 2)
    cols = [x] + [d(k) - d_last for k in range(K - 2)]
    return np.column_stack(cols)


def ns_rstyle(x, interior, boundary):
    """Second, independent construction used only as a check: the algorithm of
    R splines::ns (cubic B-splines on the augmented knot vector, first B-spline
    dropped (intercept = FALSE), then restricted to the null space of the
    second-derivative constraints at the two boundary knots)."""
    x = np.asarray(x, float)
    a, b = float(boundary[0]), float(boundary[1])
    t = np.concatenate([[a] * 4, np.asarray(interior, float), [b] * 4])
    m = t.size - 4
    spl = BSpline(t, np.eye(m), 3)
    B = spl(x)                                        # n x m B-spline basis
    const = spl.derivative(2)(np.array([a, b]))       # 2 x m second derivatives
    B, const = B[:, 1:], const[:, 1:]
    N = linalg.null_space(const)
    return B @ N


def span_residual(A, B):
    """Max relative residual of regressing each column of B on [1, A]."""
    X = np.column_stack([np.ones(A.shape[0]), A])
    coef, *_ = linalg.lstsq(X, B)
    res = B - X @ coef
    return float(np.max(np.abs(res)) / np.max(np.abs(B)))


# ---------------------------------------------------------------------------
# 3. Weighted Cox proportional hazards model, Breslow ties, Newton-Raphson
# ---------------------------------------------------------------------------
def cox_breslow(time_, event, X, w, beta0=None, tol=COX_TOL, max_iter=100):
    """Maximise the weighted Breslow partial log-likelihood

        l(b) = sum_i w_i d_i [x_i'b - log sum_{j: t_j >= t_i} w_j exp(x_j'b)]

    by Newton-Raphson with step halving.  Tied event times share the full risk
    set (Breslow).  Converged when the max absolute Newton step is < tol.
    Scaling all weights by a constant does not change the estimate, so weights
    are normalised to mean 1 for numerical stability.
    """
    time_ = np.asarray(time_, float)
    event = np.asarray(event, float)
    X = np.asarray(X, float)
    w = np.asarray(w, float)
    if np.any(~np.isfinite(w)) or np.any(w <= 0):
        raise ValueError("weights must be finite and positive")
    n, p = X.shape

    order = np.argsort(time_, kind="mergesort")
    t = time_[order]
    d = event[order]
    Xs = X[order]
    ws = w[order] / w.mean()
    uniq, first_idx, inv = np.unique(t, return_index=True, return_inverse=True)
    inv = inv.ravel()
    wd = ws * d
    dw_all = np.bincount(inv, weights=wd, minlength=uniq.size)   # weighted deaths per time
    is_ev = dw_all > 0
    fi = first_idx[is_ev]           # risk set of an event time = sorted rows fi .. n-1
    dw = dw_all[is_ev]
    sum_wdx = Xs.T @ wd

    def evaluate(beta):
        eta = Xs @ beta
        c = eta.max()                                   # overflow guard (cancels out)
        v = ws * np.exp(eta - c)
        S0 = np.cumsum(v[::-1])[::-1][fi]               # sum over risk set of w exp(eta)
        loglik = wd @ eta - dw @ (np.log(S0) + c)
        haz = np.zeros(uniq.size)
        haz[is_ev] = dw / S0
        vH = v * np.cumsum(haz)[inv]                    # v_j * sum_{event t <= t_j} dw/S0
        U = sum_wdx - Xs.T @ vH                         # score
        S1 = np.cumsum((v[:, None] * Xs)[::-1], axis=0)[::-1][fi]
        xbar = S1 / S0[:, None]
        # information = sum_t dw_t [S2_t/S0_t - xbar_t xbar_t'], with
        # sum_t dw_t S2_t/S0_t = sum_j v_j H_j x_j x_j'
        info = (Xs * vH[:, None]).T @ Xs - (xbar * dw[:, None]).T @ xbar
        return loglik, U, info

    beta = np.zeros(p) if beta0 is None else np.array(beta0, float).copy()
    ll, U, info = evaluate(beta)
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        try:
            step = linalg.solve(info, U, assume_a="pos")
        except linalg.LinAlgError:
            step = linalg.solve(info, U)
        max_step = float(np.max(np.abs(step)))
        new_beta = beta + step
        new_ll, new_U, new_info = evaluate(new_beta)
        n_half = 0
        while (not np.isfinite(new_ll) or new_ll < ll - 1e-10 * max(1.0, abs(ll))) and n_half < 40:
            step = step * 0.5
            new_beta = beta + step
            new_ll, new_U, new_info = evaluate(new_beta)
            n_half += 1
        beta, ll, U, info = new_beta, new_ll, new_U, new_info
        if max_step < tol:
            converged = True
            break
    if not converged:
        raise RuntimeError("Cox Newton-Raphson did not converge")
    return {"beta": beta, "loglik": float(ll), "info": info, "n_iter": it,
            "max_abs_score": float(np.max(np.abs(U)))}


# ---------------------------------------------------------------------------
# 4. Part 2 covariate design
# ---------------------------------------------------------------------------
def educ_category(e):
    out = np.full(e.shape, "missing", dtype=object)
    ok = ~np.isnan(e)
    out[ok & (e < 9)] = "<9"
    out[ok & (e >= 9) & (e < 12)] = "9-11"     # educ_years is integer valued: 9, 10, 11
    out[ok & (e == 12)] = "12"
    out[ok & (e > 12)] = ">12"
    return out


def pir_category(x):
    out = np.full(x.shape, "missing", dtype=object)
    ok = ~np.isnan(x)
    out[ok & (x < 1)] = "<1"
    out[ok & (x >= 1) & (x < 2)] = "1-<2"
    out[ok & (x >= 2) & (x < 4)] = "2-<4"
    out[ok & (x >= 4)] = ">=4"
    return out


def build_covariates(M, spline="esl"):
    """Core covariate set of the protocol, all knots and the BMI imputation value
    computed once in the full mortality sample M (unweighted)."""
    age = M["age"].to_numpy(float)
    age_knots = np.concatenate([[age.min()], np.percentile(age, [25, 50, 75]), [age.max()]])

    bmi = M["bmi"].to_numpy(float)
    miss = np.isnan(bmi)
    bmi_median = float(np.median(bmi[~miss]))
    bmi_filled = np.where(miss, bmi_median, bmi)
    bmi_knots = np.concatenate([[bmi_filled.min()],
                                np.percentile(bmi_filled, [100.0 / 3.0, 200.0 / 3.0]),
                                [bmi_filled.max()]])

    if spline == "esl":
        age_b = ns_esl(age, age_knots)
        bmi_b = ns_esl(bmi_filled, bmi_knots)
    else:
        age_b = ns_rstyle(age, age_knots[1:-1], age_knots[[0, -1]])
        bmi_b = ns_rstyle(bmi_filled, bmi_knots[1:-1], bmi_knots[[0, -1]])

    educ = M["educ_years"].to_numpy(float)
    if np.any(educ[~np.isnan(educ)] != np.round(educ[~np.isnan(educ)])):
        raise ValueError("educ_years has non-integer values; the 9-11 category would be ambiguous")

    factors = {
        "sex": M["sex"].astype(int).astype(str).to_numpy(),
        "race_eth": M["race_eth"].astype(int).astype(str).to_numpy(),
        "phase": M["phase"].astype(int).astype(str).to_numpy(),
        "smoke": M["smoke"].astype(str).str.strip().to_numpy(),
        "bmi_missing": np.where(miss, "1", "0"),
        "educ_cat": educ_category(educ),
        "pir_cat": pir_category(M["pir"].to_numpy(float)),
    }
    cols = [age_b, bmi_b]
    names = [f"ns_age_{k + 1}" for k in range(age_b.shape[1])] + \
            [f"ns_bmi_{k + 1}" for k in range(bmi_b.shape[1])]
    for f, levels in FACTOR_LEVELS.items():
        vals = factors[f]
        unexpected = set(np.unique(vals)) - set(levels)
        if unexpected:
            raise ValueError(f"unexpected levels in {f}: {unexpected}")
        present = [lv for lv in levels if np.any(vals == lv)]
        for lv in present[1:]:                           # first present level = reference
            cols.append((vals == lv).astype(float)[:, None])
            names.append(f"{f}[{lv}]")
    C = np.hstack(cols)
    info = {"age_knots": age_knots.tolist(), "bmi_knots": bmi_knots.tolist(),
            "bmi_median_imputed": bmi_median, "n_bmi_missing": int(miss.sum()),
            "covariate_columns": names,
            "level_counts": {f: {lv: int(np.sum(v == lv)) for lv in FACTOR_LEVELS[f]}
                             for f, v in factors.items()}}
    return C, names, factors, info


def zero_event_levels(factors, event):
    """A factor level whose members have no events has an MLE of -infinity (monotone
    likelihood).  The exact limit of all other coefficients equals the fit
    without those rows, so they are removed.  (R coxph would stop at a large
    negative coefficient with a warning; the exposure coefficient is the same to
    far more than 4 decimals.)  Positive weights in every BRR replicate make this
    condition weight-independent."""
    drop = np.zeros(event.size, bool)
    notes = []
    for f, vals in factors.items():
        for lv in np.unique(vals):
            m = vals == lv
            if event[m].sum() == 0:
                drop |= m
                notes.append(f"{f}={lv}: n={int(m.sum())}, 0 events -> rows removed (limit of MLE)")
    return drop, notes


def orthonormal_block(C):
    """Centre the covariates and replace them by an orthonormal basis of the same
    column space (pivoted QR, rank-revealing).  This is a linear re-parametrisation
    of the nuisance covariates only, so the exposure coefficient is unchanged,
    while the Newton-Raphson problem becomes well conditioned."""
    Cc = C - C.mean(axis=0)
    Q, R, piv = linalg.qr(Cc, mode="economic", pivoting=True)
    dg = np.abs(np.diag(R))
    rank = int(np.sum(dg > 1e-10 * dg[0]))
    return Q[:, :rank] * np.sqrt(C.shape[0]), rank


def design(x, Zcov):
    """Exposure first (centred; centring does not change its coefficient)."""
    return np.column_stack([x - x.mean(), Zcov])


# ---------------------------------------------------------------------------
# 5. Validation of the Cox routine against statsmodels PHReg
# ---------------------------------------------------------------------------
def phreg_fit(time_, event, X):
    """Fit statsmodels PHReg (Breslow) to its own optimum.

    statsmodels uses an undamped Newton method; from beta = 0 it can overshoot
    and overflow (it does so on the cancer model, where one covariate level has
    a single member).  In that case the optimum of statsmodels' OWN likelihood
    is found with BFGS from beta = 0 and then polished with statsmodels' Newton
    (tol 1e-12).  Nothing from the hand-written routine is used as a start."""
    import warnings
    from statsmodels.duration.hazard_regression import PHReg

    mod = PHReg(time_, X, status=event, ties="breslow")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = mod.fit(method="newton", tol=1e-12, maxiter=200)
        params = np.asarray(res.params, float)
        route = "newton from 0"
        if not (np.all(np.isfinite(params)) and np.isfinite(mod.loglike(params))):
            rb = mod.fit(method="bfgs", maxiter=5000, gtol=1e-8)
            res = mod.fit(method="newton", start_params=np.asarray(rb.params, float), tol=1e-12, maxiter=200)
            params = np.asarray(res.params, float)
            route = "bfgs from 0, then newton polish (undamped newton from 0 overflowed)"
    return params, float(mod.loglike(params)), route


def validate_cox(time_, event, X):
    """(a) weights = 1: compare with statsmodels PHReg(ties='breslow').
    (b) integer frequency weights: compare with PHReg on the row-expanded data
        (the weighted Breslow likelihood equals the expanded-data likelihood)."""
    ones = np.ones(time_.size)
    mine = cox_breslow(time_, event, X, ones)
    sm_par, sm_ll, route = phreg_fit(time_, event, X)
    diff_unw = float(np.max(np.abs(mine["beta"] - sm_par)))       # NaN propagates -> FAIL

    rng = np.random.default_rng(SEED)
    k = rng.integers(1, 4, size=time_.size)               # weights 1, 2 or 3
    mine_w = cox_breslow(time_, event, X, k.astype(float))
    idx = np.repeat(np.arange(time_.size), k)
    smw_par, _, route_w = phreg_fit(time_[idx], event[idx], X[idx])
    diff_w = float(np.max(np.abs(mine_w["beta"] - smw_par)))
    return {"max_abs_diff_unweighted": diff_unw,
            "exposure_coef_mine_unweighted": float(mine["beta"][0]),
            "exposure_coef_statsmodels_unweighted": float(sm_par[0]),
            "loglik_mine_unweighted": mine["loglik"], "loglik_statsmodels_unweighted": sm_ll,
            "statsmodels_route_unweighted": route,
            "max_abs_diff_integer_freq_weights_vs_expanded": diff_w,
            "statsmodels_route_freq_weights": route_w,
            "n_params": int(X.shape[1])}


# ---------------------------------------------------------------------------
# 6. Main
# ---------------------------------------------------------------------------
def to_py(o):
    """Convert numpy types for json (floats keep full precision via repr)."""
    if isinstance(o, dict):
        return {str(k): to_py(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [to_py(v) for v in o]
    if isinstance(o, np.ndarray):
        return [to_py(v) for v in o.tolist()]
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else None     # strict JSON: NaN/inf -> null
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", default=str(DEFAULT_INPUT), help="input csv(.gz); default = SCRAMBLED file")
    ap.add_argument("--skip-validation", action="store_true", help="skip the statsmodels PHReg checks")
    args = ap.parse_args()
    t_start = time.time()

    in_path = Path(args.input).expanduser().resolve()
    scrambled = guard_input(in_path)                       # may exit before reading anything
    label = BANNER if scrambled else "REAL DATA (post-registration run)"
    print("=" * 72)
    print(f"  {label}")
    print(f"  input: {in_path}")
    print("=" * 72)

    df = pd.read_csv(in_path, low_memory=False)
    if scrambled:
        if "SCRAMBLED" not in df.columns or not (df["SCRAMBLED"] == 1).all():
            if not osf_registered():
                sys.exit("REFUSED: file is named SCRAMBLED but its SCRAMBLED flag is not 1 for all rows; "
                         "treating it as real data, which is not allowed before OSF registration.")
            scrambled, label = False, "REAL DATA (post-registration run)"

    # ---------------- Part 1: ITT ----------------
    itt = df[as_bool(df["in_itt"])].copy()
    lnG = np.log(itt["glr"].to_numpy(float))
    if np.any(~np.isfinite(lnG)):
        raise ValueError("glr must be > 0 in the ITT sample")
    Z = itt["arm_pmeve"].to_numpy(float)
    D = itt["session_pmeve"].to_numpy(float)
    w = itt["w_mec"].to_numpy(float)
    Wr = itt[REP_COLS].to_numpy(float)
    strata = itt["strata"].to_numpy()
    psu = itt["psu"].to_numpy()
    assert set(np.unique(Z)) <= {0.0, 1.0} and set(np.unique(D)) <= {0.0, 1.0}
    assert np.all(((itt["arm"] == "PMEVE") == (Z == 1)).to_numpy()), "arm and arm_pmeve disagree"
    assert np.all(w > 0) and np.all(Wr > 0)

    X1 = np.column_stack([np.ones(Z.size), Z])
    b_wls, V_tay, n_psu = wls_taylor(lnG, X1, w, strata, psu)
    beta_itt = float(b_wls[1])
    m1, m0 = wmean_by_group(lnG, Z, w)
    assert abs((m1 - m0) - beta_itt) < 1e-12, "WLS slope must equal the weighted mean difference"
    se_taylor = float(np.sqrt(V_tay[1, 1]))

    m1r, m0r = wmean_by_group(lnG, Z, Wr)
    beta_r = m1r - m0r
    se_brr = float(np.sqrt(brr_var(beta_r, beta_itt)))

    ci_taylor = [pct(beta_itt - T975 * se_taylor), pct(beta_itt + T975 * se_taylor)]
    ci_brr = [pct(beta_itt - T975 * se_brr), pct(beta_itt + T975 * se_brr)]

    # ---------------- CACE ----------------
    p1, p0 = wmean_by_group(D, Z, w)
    cace = beta_itt / (p1 - p0)
    p1r, p0r = wmean_by_group(D, Z, Wr)
    cace_r = beta_r / (p1r - p0r)
    cace_se = float(np.sqrt(brr_var(cace_r, cace)))

    # ---------------- Part 2: mortality ----------------
    M = itt[(itt["eligstat"] == 1) & itt["fu_months_exam"].notna()].copy()
    tM = M["fu_months_exam"].to_numpy(float)
    ev = {"all": M["death_all"].to_numpy(float), "cancer": M["death_cancer"].to_numpy(float)}
    for k_, e_ in ev.items():
        if np.any(np.isnan(e_)) or not set(np.unique(e_)) <= {0.0, 1.0}:
            raise ValueError(f"event indicator {k_} must be 0/1 without missing values in M")
    wM = M["w_mec"].to_numpy(float)
    WrM = M[REP_COLS].to_numpy(float)
    DM = M["session_pmeve"].to_numpy(float)
    lnG_M = np.log(M["glr"].to_numpy(float))
    x_obs = lnG_M / LN2                                   # log2(GLR)
    x_std = x_obs - (cace / LN2) * DM                     # morning-equivalent log2(GLR)

    C, cov_names, factors, cov_info = build_covariates(M, spline="esl")

    # analytic attenuation prediction (weighted, in M)
    p_M = float(np.sum(wM * DM) / np.sum(wM))
    mu_M = float(np.sum(wM * lnG_M) / np.sum(wM))
    var_M = float(np.sum(wM * (lnG_M - mu_M) ** 2) / np.sum(wM))
    lam = 1.0 - cace ** 2 * p_M * (1.0 - p_M) / var_M
    predicted_delta = 1.0 / lam - 1.0
    n_M = lnG_M.size
    var_M_svy = var_M * n_M / (n_M - 1.0)                 # survey::svyvar convention
    predicted_delta_svyvar = 1.0 / (1.0 - cace ** 2 * p_M * (1.0 - p_M) / var_M_svy) - 1.0

    part2 = {}
    supp_part2 = {}
    val_design = {}
    for outcome, e_ in ev.items():
        drop, notes = zero_event_levels(factors, e_)
        keep = ~drop
        Zcov, rank = orthonormal_block(C[keep])
        if rank != C.shape[1]:
            notes.append(f"covariate block rank {rank} < {C.shape[1]} columns (dependent columns absorbed)")
        tk, ek, wk, WRk = tM[keep], e_[keep], wM[keep], WrM[keep]
        xo, xs = x_obs[keep], x_std[keep]

        fo = cox_breslow(tk, ek, design(xo, Zcov), wk)
        fs = cox_breslow(tk, ek, design(xs, Zcov), wk, beta0=fo["beta"])
        b_obs, b_std = float(fo["beta"][0]), float(fs["beta"][0])
        delta = (b_std - b_obs) / b_obs

        # paired BRR: re-estimate the CACE in each replicate, re-standardise, refit both
        b_obs_r = np.empty(N_REP)
        b_std_r = np.empty(N_REP)
        iters = []
        for r in range(N_REP):
            xs_r = x_obs - (cace_r[r] / LN2) * DM
            fo_r = cox_breslow(tk, ek, design(xo, Zcov), WRk[:, r], beta0=fo["beta"])
            fs_r = cox_breslow(tk, ek, design(xs_r[keep], Zcov), WRk[:, r], beta0=fs["beta"])
            b_obs_r[r], b_std_r[r] = fo_r["beta"][0], fs_r["beta"][0]
            iters += [fo_r["n_iter"], fs_r["n_iter"]]
        delta_r = (b_std_r - b_obs_r) / b_obs_r
        delta_se = float(np.sqrt(brr_var(delta_r, delta)))

        part2[outcome] = {
            "beta_obs": b_obs, "beta_std": b_std,
            "hr_obs": float(np.exp(b_obs)), "hr_std": float(np.exp(b_std)),
            "delta": delta, "delta_se_brr": delta_se,
            "delta_ci": [delta - T975 * delta_se, delta + T975 * delta_se],
            "predicted_delta": predicted_delta,
        }
        supp_part2[outcome] = {
            "n_used_in_cox": int(keep.sum()), "events_used": int(ek.sum()),
            "monotone_likelihood_notes": notes,
            "se_brr_beta_obs": float(np.sqrt(brr_var(b_obs_r, b_obs))),
            "se_brr_beta_std": float(np.sqrt(brr_var(b_std_r, b_std))),
            "delta_se_brr_centred_at_replicate_mean": float(np.sqrt(brr_var_mean_centred(delta_r))),
            "full_sample_newton_iterations": [fo["n_iter"], fs["n_iter"]],
            "full_sample_max_abs_score": [fo["max_abs_score"], fs["max_abs_score"]],
            "replicate_newton_iterations_range": [int(min(iters)), int(max(iters))],
            "replicate_delta": delta_r.tolist(),
            "replicate_cace": cace_r.tolist(),
        }
        val_design[outcome] = (tk, ek, design(xo, Zcov), wk, keep)

    # ---------------- checks of the implementation ----------------
    checks = {}
    # (i) the ESL basis spans the same space (with a constant) as an R-style ns basis
    C_r, _, _, _ = build_covariates(M, spline="rstyle")
    age_e, bmi_e = C[:, :4], C[:, 4:7]
    age_r, bmi_r = C_r[:, :4], C_r[:, 4:7]
    checks["spline_span_rel_residual_age"] = max(span_residual(age_e, age_r), span_residual(age_r, age_e))
    checks["spline_span_rel_residual_bmi"] = max(span_residual(bmi_e, bmi_r), span_residual(bmi_r, bmi_e))
    # (ii) exposure coefficient is invariant to the spline basis (ESL vs R-style ns)
    tk, ek, _, wk, keep = val_design["all"]
    Zr, _ = orthonormal_block(C_r[keep])
    b_alt = float(cox_breslow(tk, ek, design(x_obs[keep], Zr), wk)["beta"][0])
    checks["beta_obs_all_basis_invariance_abs_diff"] = abs(b_alt - part2["all"]["beta_obs"])
    # (iii) exposure coefficient is invariant to using the raw (non-orthonormalised) covariates
    Xraw = np.column_stack([x_obs[keep] - x_obs[keep].mean(), C[keep] - C[keep].mean(axis=0)])
    try:
        b_raw = float(cox_breslow(tk, ek, Xraw, wk)["beta"][0])
        checks["beta_obs_all_raw_design_abs_diff"] = abs(b_raw - part2["all"]["beta_obs"])
    except Exception as exc:  # pragma: no cover
        checks["beta_obs_all_raw_design_abs_diff"] = f"failed: {exc}"

    cox_val = None
    if not args.skip_validation:
        vals = {}
        for outcome in ("all", "cancer"):
            tk, ek, Xd, _, _ = val_design[outcome]
            vals[outcome] = validate_cox(tk, ek, Xd)
        diffs = np.array([v["max_abs_diff_unweighted"] for v in vals.values()]
                         + [v["max_abs_diff_integer_freq_weights_vs_expanded"] for v in vals.values()])
        # NaN-safe: any non-finite difference is a failure (np.max propagates NaN)
        cox_val = float(np.max(diffs[:2])) if np.all(np.isfinite(diffs[:2])) else float("nan")
        checks["cox_validation_detail"] = vals
        checks["cox_validation_pass_6_decimals"] = bool(np.all(np.isfinite(diffs)) and np.all(diffs < 5e-7))

    # ---------------- assemble output ----------------
    out = {
        "_label": label,
        "_input": str(in_path.relative_to(PROJECT_ROOT)) if in_path.is_relative_to(PROJECT_ROOT) else str(in_path),
        "n_itt": int(len(itt)), "n_z1": int(np.sum(Z == 1)), "n_z0": int(np.sum(Z == 0)),
        "beta_itt": beta_itt, "se_taylor": se_taylor, "se_brr": se_brr,
        "pct": float(pct(beta_itt)), "ci_taylor": ci_taylor, "ci_brr": ci_brr,
        "p1": float(p1), "p0": float(p0), "cace": float(cace), "cace_se_brr": cace_se,
        "cace_pct": float(pct(cace)),
        "n_mort": int(len(M)), "events_all": int(ev["all"].sum()), "events_cancer": int(ev["cancer"].sum()),
        "all": part2["all"], "cancer": part2["cancer"],
        "cox_validation_max_abs_diff": cox_val,
        "_supplementary": {
            "units": "pct, ci_taylor, ci_brr, cace_pct in percent; delta, delta_ci, predicted_delta as proportions",
            "t_quantile_0975_df49": T975,
            "brr_variance": "sum_r (theta_r - theta)^2 / (52 * (1-0.3)^2), centred at full-sample estimate",
            "se_brr_centred_at_replicate_mean": float(np.sqrt(brr_var_mean_centred(beta_r))),
            "cace_se_brr_centred_at_replicate_mean": float(np.sqrt(brr_var_mean_centred(cace_r))),
            "cace_ci_brr_pct": [float(pct(cace - T975 * cace_se)), float(pct(cace + T975 * cace_se))],
            "psu_per_stratum_itt": sorted(set(n_psu)), "n_strata_itt": len(n_psu),
            "weighted_mean_lnG_z1": float(m1), "weighted_mean_lnG_z0": float(m0),
            "lambda": lam, "p_D_in_M": p_M, "var_w_lnG_in_M": var_M,
            "predicted_delta_with_svyvar_n_over_n_minus_1": predicted_delta_svyvar,
            "covariates": cov_info,
            "part2": supp_part2,
            "checks": checks,
        },
        "_software": {"python": platform.python_version(), "numpy": np.__version__,
                      "pandas": pd.__version__, "scipy": scipy.__version__},
        "_run": {"timestamp": _dt.datetime.now().isoformat(timespec="seconds"),
                 "seconds": round(time.time() - t_start, 1)},
    }
    try:
        import statsmodels
        out["_software"]["statsmodels"] = statsmodels.__version__
    except ImportError:  # pragma: no cover
        pass

    out_name = "key_estimates_SCRAMBLED.json" if scrambled else "key_estimates_REAL.json"
    out_path = SCRIPT_DIR / out_name
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(to_py(out), fh, indent=2, ensure_ascii=False)

    # ---------------- console summary ----------------
    print(f"ITT n={out['n_itt']} (Z=1 {out['n_z1']}, Z=0 {out['n_z0']})")
    print(f"beta_ITT={beta_itt:.6f}  SE_Taylor={se_taylor:.6f}  SE_BRR={se_brr:.6f}  "
          f"pct={out['pct']:.4f}%  CI_T=[{ci_taylor[0]:.4f}, {ci_taylor[1]:.4f}]  "
          f"CI_BRR=[{ci_brr[0]:.4f}, {ci_brr[1]:.4f}]")
    print(f"p1={p1:.6f} p0={p0:.6f}  CACE={cace:.6f}  SE_BRR={cace_se:.6f}  pct={out['cace_pct']:.4f}%")
    print(f"Part 2: n={out['n_mort']}, events all={out['events_all']}, cancer={out['events_cancer']}")
    for oc in ("all", "cancer"):
        r_ = part2[oc]
        print(f"  {oc:6s} beta_obs={r_['beta_obs']:.6f} beta_std={r_['beta_std']:.6f} "
              f"HR_obs={r_['hr_obs']:.6f} HR_std={r_['hr_std']:.6f} Delta={r_['delta']:.6f} "
              f"SE={r_['delta_se_brr']:.6f} CI=[{r_['delta_ci'][0]:.6f}, {r_['delta_ci'][1]:.6f}] "
              f"pred={r_['predicted_delta']:.6f}")
    if cox_val is not None:
        ok = checks["cox_validation_pass_6_decimals"]
        print(f"Cox validation vs statsmodels PHReg (weights=1): max|diff|={cox_val:.3e}; "
              f"integer-weight vs expanded-data check included -> {'PASS' if ok else 'FAIL'} (>= 6 decimals)")
    print(f"written: {out_path}")
    print("=" * 72)
    print(f"  {label}")
    print("=" * 72)


if __name__ == "__main__":
    main()
