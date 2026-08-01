# -*- coding: utf-8 -*-
"""Year-2, module 2f: a sixth film, and a cap the model was missing.

Module 2e left the closure length as the binding uncertainty, estimated from
essentially two informative films.  Groner, George, McLean and Carcia (48th SVC
Technical Conference Proceedings, 2005, open access) published a four-point
thickness series that runs straight through the steep part of the curve:

    2.5 nm  no improvement over the uncoated polymer (~1 g m^-2 day^-1)
    5   nm  about one order of magnitude below uncoated  (~0.1)
    10  nm  ~2e-3
    26  nm  ~1e-3

Adding it exposed a gap in the single-layer model.  `WVTR = K f(d)` has no
upper bound, so for a barely-nucleated film it predicts fluxes larger than the
bare substrate can pass, which is unphysical.  Nothing in the earlier data sat
near that limit, so the omission was invisible until now.  The fix is the
substrate in series, exactly as in the multilayer model:

    WVTR = [ 1/(K f(d)) + 1/W_bare ]^-1

Groner's test is a tritiated-water cell at 100 % RH with the film clamped
against a Viton o-ring; the temperature is not stated and is taken as ambient,
which is carried as an assumption rather than a measurement.
"""
import numpy as np
import hier_v2 as V2
from hier_v2 import arr, f_shape, ladder_wvtr, stretch_sample, P_PAR, CLOSURE_EPS, q

FILMS = ["Wu high-q", "Wu poor", "Lee", "Carcia", "Buelow", "Groner"]
NF = len(FILMS)
IDX = {f: i for i, f in enumerate(FILMS)}
MATERIAL = {"Lee": "pV3D3", "Buelow": "PP"}
# bare-substrate WVTR that caps a barely-nucleated film [g m^-2 day^-1]
W_BARE = {"Wu high-q": 10.0, "Wu poor": 10.0, "Lee": 10.0, "Carcia": 0.5,
          "Buelow": 0.5, "Groner": 1.0}

BF = V2.BUELOW_FLOOR
SINGLE = [
    ("Wu high-q", 15.0, 6.7e-3, 0.10, 38.0, False),
    ("Wu high-q", 20.0, 0.7e-3, 0.10, 38.0, False),
    ("Wu high-q", 30.0, 0.8e-3, 0.10, 38.0, False),
    ("Wu high-q", 50.0, 1.3e-3, 0.10, 38.0, False),
    ("Wu high-q", 60.0, 4.7e-3, 0.10, 38.0, False),
    ("Carcia",    25.0, 1.7e-5, 0.15, 38.0, False),
    ("Carcia",    25.0, 6.5e-5, 0.15, 60.0, False),
    ("Buelow",    50.0, 5.3e-3 - BF, 0.20, 60.0, False),
    ("Buelow",   100.0, 6.4e-3 - BF, 0.15, 60.0, False),
    ("Buelow",    25.0, 5.3e-2 - BF, 0.30, 60.0, True),
    ("Groner",     2.5, 1.0e0,  0.30, 25.0, False),
    ("Groner",     5.0, 1.0e-1, 0.30, 25.0, False),
    ("Groner",    10.0, 2.0e-3, 0.20, 25.0, False),
    ("Groner",    26.0, 1.0e-3, 0.20, 25.0, False),
]
LADDER = list(V2.LADDER)

# 0 lf0 | 1 lmu_dc | 2 lC | 3 d_crit | 4 lmu_f | 5 sigma_f | 6..11 lf_res
# 12 lK | 13 lm_pv | 14 lLpath | 15 Ea | 16 lm_pp | 17 sigma_dc | 18..23 ldc
NDIM = 18 + NF


def single_wvtr(d, K, Ea, T, f0, dc, C, d_crit, f_res, w_bare):
    g = K * arr(T, Ea) * f_shape(d, f0, dc, C, d_crit, f_res)
    return 1.0 / (1.0 / g + 1.0 / w_bare)          # substrate in series


def log_posterior(p):
    lf0, lmu_dc, lC, d_crit, lmu_f, sig_f = p[:6]
    lfres = p[6:6 + NF]
    lK, lm_pv, lL, Ea, lm_pp, sig_dc = (p[6+NF], p[7+NF], p[8+NF], p[9+NF],
                                        p[10+NF], p[11+NF])
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
        mod = single_wvtr(d, K, Ea, T, f0, dc[i], C, d_crit, fres[i], W_BARE[film])
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
    rng = np.random.default_rng(41)
    nw, nsteps, burn = 80, 8000, 3000
    start = np.array([np.log10(0.372), np.log10(1.3), np.log10(40.0), 46.0,
                      -7.4, 0.4, *([-7.4] * NF), np.log10(2e4), np.log10(25.0),
                      np.log10(5e-3), 45.0, np.log10(3.0), 0.25,
                      *([np.log10(1.3)] * NF)])
    p0 = start + 0.02 * rng.standard_normal((nw, NDIM))
    p0[:, 5], p0[:, 11 + NF] = np.abs(p0[:, 5]), np.abs(p0[:, 11 + NF])
    print(f"sampling ({NF} films, 3 temperatures, {NDIM} parameters) ...")
    chain, ar = stretch_sample(log_posterior, p0, nsteps, rng)
    flat = chain[burn:].reshape(-1, NDIM)
    print(f"acceptance {ar:.2f}, {len(flat)} samples")
    h1, h2 = chain[burn:nsteps//2], chain[nsteps//2:]
    for k, nm in ((1, "mu_dc"), (11+NF, "sigma_dc"), (3, "d_crit"), (9+NF, "Ea")):
        print(f"  split-half {nm:9s} {abs(h1[...,k].mean()-h2[...,k].mean())/chain[burn:][...,k].std():.2f} sd")

    print("\nposterior (median [90 % CI])")
    for k, nm, un in ((9+NF, "activation energy", "kJ/mol"),
                      (3, "cracking onset", "nm"),
                      (5, "spread in floor", "dec"),
                      (11+NF, "spread in closure length", "dec")):
        lo, md, hi = q(flat[:, k]); print(f"  {nm:26s} {md:7.2f}  [{lo:.2f}, {hi:.2f}] {un}")
    print("\n  per-film closure length (nm) and defect floor")
    for i, f in enumerate(FILMS):
        a = q(flat[:, 12 + NF + i], lambda x: 10 ** x)
        b = q(flat[:, 6 + i], lambda x: 10 ** x)
        print(f"    {f:11s} {a[1]:5.2f} [{a[0]:.2f}, {a[2]:5.2f}]    {b[1]:.2e} [{b[0]:.1e}, {b[2]:.1e}]")

    rn1, rn2 = rng.standard_normal(len(flat)), rng.standard_normal(len(flat))
    lf_new = flat[:, 4] + flat[:, 5] * rn1
    dc_new = 10 ** (flat[:, 1] + flat[:, 11 + NF] * rn2)
    N = np.log(10 ** flat[:, 0] / (CLOSURE_EPS * 10 ** lf_new))
    thr = flat[:, 3] / N
    lo, md, hi = q(thr)
    print(f"\n  closure lengths needed for closure: N = {q(N)[1]:.0f} [{q(N)[0]:.0f}, {q(N)[2]:.0f}]")
    print(f"  critical closure length     {md:5.2f}  [{lo:.2f}, {hi:.2f}] nm"
          f"   (was 2.01 [1.41, 2.71] with five films)")
    print("\n  P(a usable window exists)")
    for i, f in enumerate(FILMS):
        dc_i = 10 ** flat[:, 12 + NF + i]
        Ni = np.log(10 ** flat[:, 0] / (CLOSURE_EPS * 10 ** flat[:, 6 + i]))
        print(f"    {f:11s} {100*np.mean(dc_i < flat[:, 3]/Ni):5.1f} %")
    print(f"    {'unmeasured':11s} {100*np.mean(dc_new < thr):5.1f} %")
    lo, md, hi = q(dc_new * N)
    print(f"\n  window lower edge, unmeasured film   {md:5.1f} [{lo:.1f}, {hi:.1f}] nm")

    pm = np.median(flat, axis=0)
    f0, C, d_crit, K, Ea = 10**pm[0], 10**pm[2], pm[3], 10**pm[6+NF], pm[9+NF]
    fres, dcv = 10**pm[6:6+NF], 10**pm[12+NF:12+2*NF]
    mult = {"pV3D3": 10**pm[7+NF], "PP": 10**pm[10+NF]}
    res = []
    for film, d, val, s_, T, cens in SINGLE:
        i = IDX[film]
        m = single_wvtr(d, K, Ea, T, f0, dcv[i], C, d_crit, fres[i], W_BARE[film])
        r = np.log10(m/val)
        if cens and r > 0: r = 0.0
        res.append((abs(r), f"{film} single {d:.0f}nm", r))
    for film, n, d_in, d_org, d_first, val, s_, da, T in LADDER:
        i = IDX[film]
        P = mult.get(MATERIAL.get(film, ""), 1.0)*P_PAR*arr(T, Ea)
        m = ladder_wvtr(n, d_in, d_org, d_first,
                        f_shape(d_in, f0, dcv[i], C, d_crit, fres[i]), P, da,
                        film == "Lee", 10**pm[8+NF])
        res.append((abs(np.log10(m/val)), f"{film} ladder n={n}", np.log10(m/val)))
    print("\n  worst residuals (log10)")
    for _, lab, r in sorted(res, reverse=True)[:4]:
        print(f"    {lab:26s} {r:+.2f}")
    print(f"    rms over {len(res)} observations: "
          f"{np.sqrt(np.mean([r**2 for _,_,r in res])):.3f} decades")
    np.save("hier_v4_posterior.npy", flat)
    print("\nsaved hier_v4_posterior.npy")
