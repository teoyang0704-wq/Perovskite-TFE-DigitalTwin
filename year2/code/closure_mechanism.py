# -*- coding: utf-8 -*-
"""Elimination versus constriction: the diagnostic, and what data it needs.

f = N pi r^2.  Permeation gives f(d); decorated defect counts give N(d).  Both
observables carry a nucleation term on top of a floor,

    y(d) = A exp(-d / L) + C ,

so the decay length must be fitted together with the floor, not read off a
straight line -- otherwise points in the floor-dominated regime bias L upwards.
The ratio of the two fitted decay lengths names the mechanism:

    L_count = L_flux   defects are eliminated, radii unchanged
    L_count > L_flux   defects also constrict, radii shrink with thickness
"""
import numpy as np
from scipy.optimize import least_squares


def fit_decay(d, y, sigma_dex=0.15):
    """Fit y = A exp(-d/L) + C in log space. Returns (L, sigma_L)."""
    d, y = np.asarray(d, float), np.asarray(y, float)
    if len(d) < 4:
        return np.nan, np.nan

    def resid(p):
        A, L, C = 10 ** p[0], np.exp(p[1]), 10 ** p[2]
        return (np.log10(A * np.exp(-d / L) + C) - np.log10(y)) / sigma_dex

    p0 = [np.log10(y[0]), np.log(max(d[1] - d[0], 1.0)), np.log10(y[-1] * 0.9)]
    try:
        out = least_squares(resid, p0, method="lm", max_nfev=4000)
    except Exception:
        return np.nan, np.nan
    L = float(np.exp(out.x[1]))
    J = out.jac
    try:
        cov = np.linalg.inv(J.T @ J)
        sL = L * float(np.sqrt(max(cov[1, 1], 0.0)))
    except np.linalg.LinAlgError:
        sL = np.nan
    return L, sL


def counts_are_flat(d, N, sigma_dex=0.15):
    """Chi-square test of the hypothesis that the defect count does not change."""
    lN = np.log10(np.asarray(N, float))
    chi2 = np.sum(((lN - lN.mean()) / sigma_dex) ** 2)
    return chi2 < len(lN) + 2 * np.sqrt(2 * len(lN))     # ~2 sigma of chi2_(n-1)


def mechanism(Lc, sLc, Lf, sLf, k=2.0, flat_counts=False):
    if flat_counts and np.isfinite(Lf) and Lf > 0:
        # counts unchanged while flux decays: constriction, no fit needed
        return np.inf, np.nan, "CONSTRICTION"
    if not np.all(np.isfinite([Lc, sLc, Lf, sLf])) or Lf <= 0:
        return np.nan, np.nan, "indeterminate"
    ratio = Lc / Lf
    err = ratio * np.hypot(sLc / Lc, sLf / Lf)
    if ratio - k * err > 1.0:
        return ratio, err, "CONSTRICTION"
    if ratio + k * err < 1.0:
        return ratio, err, "radii grow -- model falsified"
    return ratio, err, "ELIMINATION (or undecided)"


def simulate(mode, ds, n_rep=250, sigma_dex=0.15, L_true=1.2, seed=0,
             floor_ratio=3e-6):
    """Fraction of synthetic datasets on which the diagnostic answers correctly.

    floor_ratio sets how far below the d=0 value the particulate floor sits;
    3e-6 reproduces the Wu and Groner films.
    """
    rng = np.random.default_rng(seed)
    L_N, L_f = (L_true, L_true) if mode == "elimination" else (2 * L_true, L_true)
    ds = np.asarray(ds, float)
    right = 0
    for _ in range(n_rep):
        N = 1e4 * (np.exp(-ds / L_N) + floor_ratio) * 10 ** (sigma_dex * rng.standard_normal(len(ds)))
        F = 1e-2 * (np.exp(-ds / L_f) + floor_ratio) * 10 ** (sigma_dex * rng.standard_normal(len(ds)))
        lc, slc = fit_decay(ds, N, sigma_dex)
        lf, slf = fit_decay(ds, F, sigma_dex)
        _, _, v = mechanism(lc, slc, lf, slf, flat_counts=counts_are_flat(ds, N, sigma_dex))
        right += (v == "CONSTRICTION") == (mode == "constriction")
    return right / n_rep


if __name__ == "__main__":
    print("Power of the diagnostic, with the particulate floor included")
    print("  scatter 0.15 dex, true closure length 1.2 nm, 250 realisations\n")
    print(f"  {'thickness points':36s} {'elimination':>12s} {'constriction':>13s}")
    cases = (("15, 20, 30, 50 nm", [15, 20, 30, 50]),
             ("Klumbies 15-100 nm", [15, 20, 30, 50, 75, 100]),
             ("5, 10, 15, 20 nm", [5, 10, 15, 20]),
             ("Groner 2.5, 5, 10, 26 nm", [2.5, 5, 10, 26]),
             ("2.5, 5, 10, 15, 20, 26 nm", [2.5, 5, 10, 15, 20, 26]))
    for label, ds in cases:
        e = simulate("elimination", ds, seed=3)
        c = simulate("constriction", ds, seed=11)
        print(f"  {label:36s} {100*e:11.0f}% {100*c:12.0f}%")
    print("\n  Reading: the second column is how often a truly eliminating film is")
    print("  not falsely called constricting, the third how often constriction is")
    print("  detected. A series that never goes below 15 nm can still be read, but")
    print("  its power to *detect* constriction is what decides whether one paper")
    print("  is enough or whether a sub-10 nm series has to be found as well.")
