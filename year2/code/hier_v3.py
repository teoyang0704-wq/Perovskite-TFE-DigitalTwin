# -*- coding: utf-8 -*-
"""Year-2, module 2c: the closure length is a property of the film, not of the chemistry.

Module 2b assumed one closure length for every ALD alumina film.  Adding
Buelow's plasma-ALD data broke that assumption visibly: the model missed their
25 nm point by 1.4 decades, and the strain propagated into the cracking onset,
which fell to 29 nm and left the design window empty in one draw out of six.

The physical reading is straightforward.  Closure length measures how quickly
nucleation islands merge into a continuous film, which depends on substrate
preparation, precursor chemistry and whether the process is thermal or
plasma-enhanced.  It is not a constant of Al2O3.  Here it is given the same
hierarchical treatment as the particulate floor:

    log10 d_close,i ~ Normal(log10 mu_dc, sigma_dc)

Only the cracking onset, which is mechanical, stays shared.
"""
import numpy as np
import hier_v2 as V2
from hier_v2 import (SINGLE, LADDER, FILMS, NF, IDX, MATERIAL, arr, f_shape,
                     ladder_wvtr, stretch_sample, P_PAR, CLOSURE_EPS, q)

# 0 lf0 | 1 lmu_dc | 2 lC | 3 d_crit | 4 lmu_f | 5 sigma_f | 6..10 lf_res
# 11 lK | 12 lm_pV3D3 | 13 lLpath | 14 Ea | 15 lm_PP | 16 sigma_dc | 17..21 ldc_i
NDIM = 17 + NF


def log_posterior(p):
    lf0, lmu_dc, lC, d_crit, lmu_f, sig_f = p[:6]
    lfres = p[6:6 + NF]
    lK, lm_pv, lL, Ea, lm_pp, sig_dc = p[6+NF], p[7+NF], p[8+NF], p[9+NF], p[10+NF], p[11+NF]
    ldc = p[12+NF:12+2*NF]
    if not (0 < sig_f < 3 and 0 < sig_dc < 1.5 and 10 < d_crit < 120 and 5 < Ea < 120):
        return -np.inf
    f0, C, K = 10 ** lf0, 10 ** lC, 10 ** lK
    fres, dc = 10 ** lfres, 10 ** ldc
    if not (1e-11 < fres.min() and fres.max() < 1e-3 and 0.05 < dc.min() and dc.max() < 20):
        return -np.inf
    mult = {"pV3D3": 10 ** lm_pv, "PP": 10 ** lm_pp}

    lp = 0.0
    lp += -0.5 * ((lf0 - np.log10(0.372)) / 0.6) ** 2
    lp += -0.5 * ((lmu_dc - np.log10(1.3)) / 0.4) ** 2
    lp += -0.5 * ((lC - np.log10(40.0)) / 0.8) ** 2
    lp += -0.5 * ((d_crit - 46.0) / 12.0) ** 2
    lp += -0.5 * ((lmu_f + 7.4) / 1.2) ** 2
    lp += -0.5 * (sig_f / 0.6) ** 2
    lp += -0.5 * (sig_dc / 0.35) ** 2
    lp += -0.5 * ((lm_pv - np.log10(30.0)) / 0.5) ** 2
    lp += -0.5 * ((lm_pp - np.log10(3.0)) / 1.0) ** 2
    lp += -0.5 * ((lL - np.log10(5e-3)) / 0.3) ** 2
    lp += -0.5 * ((lK - np.log10(2e4)) / 2.0) ** 2
    lp += -0.5 * ((Ea - 45.0) / 20.0) ** 2
    lp += np.sum(-0.5 * ((lfres - lmu_f) / sig_f) ** 2 - np.log(sig_f))
    lp += np.sum(-0.5 * ((ldc - lmu_dc) / sig_dc) ** 2 - np.log(sig_dc))

    for film, d, val, s_, T, cens in SINGLE:
        i = IDX[film]
        mod = K * arr(T, Ea) * f_shape(d, f0, dc[i], C, d_crit, fres[i])
        r = (np.log10(mod) - np.log10(val)) / s_
        lp += -0.5 * r ** 2 if not cens else (0.0 if r > 0 else -0.5 * r ** 2)
    for film, n, d_in, d_org, d_first, val, s_, da, T in LADDER:
        i = IDX[film]
        P = mult.get(MATERIAL.get(film, ""), 1.0) * P_PAR * arr(T, Ea)
        f = f_shape(d_in, f0, dc[i], C, d_crit, fres[i])
        mod = ladder_wvtr(n, d_in, d_org, d_first, f, P, da, film == "Lee", 10 ** lL)
        lp += -0.5 * ((np.log10(mod) - np.log10(val)) / s_) ** 2
    return lp if np.isfinite(lp) else -np.inf


if __name__ == "__main__":
    rng = np.random.default_rng(31)
    nw, nsteps, burn = 72, 7000, 2500
    start = np.array([np.log10(0.372), np.log10(1.3), np.log10(40.0), 46.0,
                      -7.4, 0.4, *([-7.4] * NF), np.log10(2e4), np.log10(25.0),
                      np.log10(5e-3), 45.0, np.log10(3.0), 0.25,
                      *([np.log10(1.3)] * NF)])
    p0 = start + 0.02 * rng.standard_normal((nw, NDIM))
    p0[:, 5], p0[:, 11 + NF] = np.abs(p0[:, 5]), np.abs(p0[:, 11 + NF])
    print(f"sampling ({NF} films, 2 temperatures, {NDIM} parameters) ...")
    chain, ar = stretch_sample(log_posterior, p0, nsteps, rng)
    flat = chain[burn:].reshape(-1, NDIM)
    print(f"acceptance {ar:.2f}, {len(flat)} samples")
    for k, nm in ((4, "log mu_f"), (5, "sigma_f"), (3, "d_crit"), (9 + NF, "Ea")):
        h1, h2 = chain[burn:nsteps//2], chain[nsteps//2:]
        print(f"  split-half shift {nm:10s} "
              f"{abs(h1[...,k].mean()-h2[...,k].mean())/chain[burn:][...,k].std():.2f} sd")

    print("\nposterior (median [90 % CI])")
    for k, nm, sc, un in ((9+NF, "activation energy", None, "kJ/mol"),
                          (3, "cracking onset d_crit", None, "nm"),
                          (5, "spread in floor", None, "decades"),
                          (11+NF, "spread in closure length", None, "decades")):
        lo, md, hi = q(flat[:, k]); print(f"  {nm:26s} {md:7.2f}  [{lo:.2f}, {hi:.2f}] {un}")
    print("\n  per-film closure length and defect floor")
    for i, f in enumerate(FILMS):
        a = q(flat[:, 12 + NF + i], lambda x: 10 ** x)
        b = q(flat[:, 6 + i], lambda x: 10 ** x)
        print(f"    {f:11s} d_close {a[1]:5.2f} [{a[0]:.2f}, {a[2]:.2f}] nm    "
              f"f_res {b[1]:.2e} [{b[0]:.1e}, {b[2]:.1e}]")

    rn = rng.standard_normal(len(flat))
    lf_new = flat[:, 4] + flat[:, 5] * rn
    ldc_new = flat[:, 1] + flat[:, 11 + NF] * rng.standard_normal(len(flat))
    d_close = 10 ** ldc_new * np.log(10 ** flat[:, 0] / (CLOSURE_EPS * 10 ** lf_new))
    w = flat[:, 3] - d_close
    print("\n  design window for a film nobody has measured")
    for lab, v in (("lower edge", d_close), ("upper edge", flat[:, 3]), ("width", w)):
        lo, md, hi = q(v); print(f"    {lab:12s} {md:6.1f}  [{lo:.1f}, {hi:.1f}] nm")
    print(f"    empty window in {100*np.mean(w <= 0):.1f} % of draws")

    print("\n  worst residuals (log10) at the posterior median")
    pm = np.median(flat, axis=0)
    f0, C, d_crit, K, Ea = 10**pm[0], 10**pm[2], pm[3], 10**pm[6+NF], pm[9+NF]
    fres, dcv = 10**pm[6:6+NF], 10**pm[12+NF:12+2*NF]
    mult = {"pV3D3": 10**pm[7+NF], "PP": 10**pm[10+NF]}
    res = []
    for film, d, val, s_, T, cens in SINGLE:
        i = IDX[film]
        m = K*arr(T, Ea)*f_shape(d, f0, dcv[i], C, d_crit, fres[i])
        r = np.log10(m/val)
        if cens and r > 0: r = 0.0
        res.append((abs(r), f"{film} single {d:.0f}nm @{T:.0f}C", r))
    for film, n, d_in, d_org, d_first, val, s_, da, T in LADDER:
        i = IDX[film]
        P = mult.get(MATERIAL.get(film, ""), 1.0)*P_PAR*arr(T, Ea)
        m = ladder_wvtr(n, d_in, d_org, d_first,
                        f_shape(d_in, f0, dcv[i], C, d_crit, fres[i]), P, da,
                        film == "Lee", 10**pm[8+NF])
        res.append((abs(np.log10(m/val)), f"{film} ladder n={n} @{T:.0f}C", np.log10(m/val)))
    for _, lab, r in sorted(res, reverse=True)[:5]:
        print(f"    {lab:28s} {r:+.2f}")
    rms = np.sqrt(np.mean([r**2 for _, _, r in res]))
    print(f"    rms over all {len(res)} observations: {rms:.3f} decades")
    np.save("hier_v3_posterior.npy", flat)
    print("\nsaved hier_v3_posterior.npy")
