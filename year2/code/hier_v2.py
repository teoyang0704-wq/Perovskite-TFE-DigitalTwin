# -*- coding: utf-8 -*-
"""Year-2, module 2b: five films, two temperatures, one population.

Module 2 pooled four films all measured near 38 C, so the activation energy of
the organic transport could stay a prior.  Adding Buelow's plasma-ALD films,
measured at 60 C, forces the issue: a temperature has to be carried explicitly,
and the activation energy becomes a parameter the data can speak to.

Two things make that possible without inventing information:

  * Carcia measured *the same film* at 38 C and 60 C (1.7e-5 and 6.5e-5), which
    constrains the activation energy directly rather than through a prior;
  * Buelow measured their fixture floor with an impermeable glass lid
    (6e-4 g m^-2 day^-1), so their readings can be corrected before use --
    the module-1 finding feeding back into how the data are consumed.

Everything else follows module 2: one shared defect shape for ALD alumina, one
particulate floor per film drawn from a population, and a permeability
multiplier for each organic material (parylene-C, pV3D3, plasma polymer).
"""
import numpy as np

DAY, KG_TO_G, RGAS = 86400.0, 1e3, 8.314
P_PAR, T_REF = 7.5e-13, 311.15         # parylene-C at 38 C
R_SUB, R_PIN = 2.88e7, 50e-9
BUELOW_FLOOR = 6.0e-4                  # measured with a glass lid, same fixture

FILMS = ["Wu high-q", "Wu poor", "Lee", "Carcia", "Buelow"]
NF = len(FILMS)
MATERIAL = {"Lee": "pV3D3", "Buelow": "PP"}      # others: parylene-C

# --- bare single layers: (film, d_nm, WVTR, sigma_dex, T_C, censored_low) ---
SINGLE = [
    ("Wu high-q", 15.0, 6.7e-3, 0.10, 38.0, False),
    ("Wu high-q", 20.0, 0.7e-3, 0.10, 38.0, False),
    ("Wu high-q", 30.0, 0.8e-3, 0.10, 38.0, False),
    ("Wu high-q", 50.0, 1.3e-3, 0.10, 38.0, False),
    ("Wu high-q", 60.0, 4.7e-3, 0.10, 38.0, False),
    ("Carcia",    25.0, 1.7e-5, 0.15, 38.0, False),
    ("Carcia",    25.0, 6.5e-5, 0.15, 60.0, False),   # same film, second temperature
    ("Buelow",    50.0, 5.3e-3 - BUELOW_FLOOR, 0.20, 60.0, False),
    ("Buelow",   100.0, 6.4e-3 - BUELOW_FLOOR, 0.15, 60.0, False),
    ("Buelow",    25.0, 5.3e-2 - BUELOW_FLOOR, 0.30, 60.0, True),   # "> 1 order worse"
]
# --- ladders: (film, n_inorg, d_in, d_org, d_first, WVTR, sigma, da, T_C) ---
LADDER = [
    ("Wu high-q", 1, 50.0, 500.0,   0.0, 1.70e-4, 0.10, 1.00, 38.0),
    ("Wu poor",   1, 50.0, 500.0,   0.0, 1.6,     0.30, 1.00, 38.0),
    ("Wu poor",   2, 50.0, 500.0,   0.0, 1.6e-1,  0.30, 1.00, 38.0),
    ("Wu poor",   3, 50.0, 500.0,   0.0, 1.3e-1,  0.30, 1.00, 38.0),
    ("Lee",       1, 21.5, 100.0, 800.0, 3.0e-3,  0.10, 0.90, 38.0),
    ("Lee",       2, 21.5, 100.0, 800.0, 6.6e-4,  0.10, 0.90, 38.0),
    ("Lee",       3, 21.5, 100.0, 800.0, 5.4e-4,  0.10, 0.90, 38.0),
    ("Lee",       4, 21.5, 100.0, 800.0, 5.3e-4,  0.10, 0.90, 38.0),
    # 1.5 and 3.5 dyads = 2 and 4 oxide layers; fixture floor already removed
    ("Buelow",    2, 25.0, 100.0,   0.0, 3.6e-3 - BUELOW_FLOOR, 0.15, 0.90, 60.0),
    ("Buelow",    4, 25.0, 100.0,   0.0, 1.2e-3 - BUELOW_FLOOR, 0.20, 0.90, 60.0),
]
PAD, PERIM = 3e-3, 4 * 3e-3
CLOSURE_EPS = 0.01

# 0 lf0 | 1 ldc | 2 lC | 3 d_crit | 4 lmu | 5 sigma | 6..10 lf_res | 11 lK
# 12 l_m_pV3D3 | 13 l_Lpath | 14 Ea_kJ | 15 l_m_PP
NDIM = 16
IDX = {f: i for i, f in enumerate(FILMS)}


def arr(T_C, Ea_kJ):
    return np.exp(-Ea_kJ * 1e3 / RGAS * (1.0 / (T_C + 273.15) - 1.0 / T_REF))


def f_shape(d, f0, dc, C, d_crit, f_res):
    return f0 * np.exp(-d / dc) + f_res * (1.0 + C * max(0.0, (d - d_crit) / d_crit) ** 2)


def ladder_wvtr(n, d_in, d_org, d_first, f, P_org, da, edge, L_path):
    d_in, d_org, d_first = d_in * 1e-9, d_org * 1e-9, d_first * 1e-9
    s = R_PIN * np.sqrt(np.pi / f)
    tau2 = 1.0 + s ** 2 * np.log(max(s / R_PIN, 1.001)) / (2 * np.pi * d_org ** 2)
    R = R_SUB + n * d_in / (f * P_org) + (n - 1) * d_org * tau2 / P_org
    if d_first > 0:
        R += d_first / P_org
    w = da / R * KG_TO_G * DAY
    if edge:
        t_org = d_first + (n - 1) * d_org
        w += P_org * da * (PERIM * t_org) / (L_path * PAD ** 2) * KG_TO_G * DAY
    return w


def log_posterior(p):
    lf0, ldc, lC, d_crit, lmu, sig = p[:6]
    lfres, lK, lm_pv, lL, Ea, lm_pp = p[6:6+NF], p[6+NF], p[7+NF], p[8+NF], p[9+NF], p[10+NF]
    if not (0.0 < sig < 3.0 and 10.0 < d_crit < 120.0 and 5.0 < Ea < 120.0):
        return -np.inf
    f0, dc, C, K = 10 ** lf0, 10 ** ldc, 10 ** lC, 10 ** lK
    fres, L_path = 10 ** lfres, 10 ** lL
    mult = {"pV3D3": 10 ** lm_pv, "PP": 10 ** lm_pp}
    if not (0.05 < dc < 10.0 and 1e-11 < fres.min() and fres.max() < 1e-3):
        return -np.inf

    lp = 0.0
    lp += -0.5 * ((lf0 - np.log10(0.372)) / 0.6) ** 2
    lp += -0.5 * ((ldc - np.log10(1.07)) / 0.35) ** 2
    lp += -0.5 * ((lC - np.log10(40.0)) / 0.8) ** 2
    lp += -0.5 * ((d_crit - 46.0) / 12.0) ** 2
    lp += -0.5 * ((lmu + 7.4) / 1.2) ** 2
    lp += -0.5 * (sig / 0.6) ** 2
    lp += -0.5 * ((lm_pv - np.log10(30.0)) / 0.5) ** 2
    lp += -0.5 * ((lm_pp - np.log10(3.0)) / 1.0) ** 2       # plasma polymer, unknown
    lp += -0.5 * ((lL - np.log10(5e-3)) / 0.3) ** 2
    lp += -0.5 * ((lK - np.log10(2e4)) / 2.0) ** 2
    lp += -0.5 * ((Ea - 45.0) / 20.0) ** 2                  # broad, data-driven
    lp += np.sum(-0.5 * ((lfres - lmu) / sig) ** 2 - np.log(sig))

    for film, d, val, s_, T, cens in SINGLE:
        mod = K * arr(T, Ea) * f_shape(d, f0, dc, C, d_crit, fres[IDX[film]])
        r = (np.log10(mod) - np.log10(val)) / s_
        lp += -0.5 * r ** 2 if not cens else (0.0 if r > 0 else -0.5 * r ** 2)
    for film, n, d_in, d_org, d_first, val, s_, da, T in LADDER:
        P = mult.get(MATERIAL.get(film, ""), 1.0) * P_PAR * arr(T, Ea)
        f = f_shape(d_in, f0, dc, C, d_crit, fres[IDX[film]])
        mod = ladder_wvtr(n, d_in, d_org, d_first, f, P, da, film == "Lee", L_path)
        lp += -0.5 * ((np.log10(mod) - np.log10(val)) / s_) ** 2
    return lp if np.isfinite(lp) else -np.inf


def stretch_sample(logpost, p0, nsteps, rng, a=2.0):
    nw, nd = p0.shape
    p, lp = p0.copy(), np.array([logpost(x) for x in p0])
    chain, acc = np.empty((nsteps, nw, nd)), 0
    for t in range(nsteps):
        for half in (0, 1):
            me = np.where(np.arange(nw) % 2 == half)[0]
            you = np.where(np.arange(nw) % 2 != half)[0]
            for i in me:
                j = you[rng.integers(len(you))]
                z = ((a - 1.0) * rng.random() + 1.0) ** 2 / a
                q = p[j] + z * (p[i] - p[j])
                lq = logpost(q)
                if np.log(rng.random()) < (nd - 1) * np.log(z) + lq - lp[i]:
                    p[i], lp[i], acc = q, lq, acc + 1
        chain[t] = p
    return chain, acc / (nsteps * nw)


def q(v, s=lambda x: x):
    return np.percentile(s(v), [5, 50, 95])


if __name__ == "__main__":
    rng = np.random.default_rng(23)
    nw, nsteps, burn = 56, 6000, 2000
    start = np.array([np.log10(0.372), np.log10(1.07), np.log10(40.0), 46.0,
                      -7.4, 0.4, *([-7.4] * NF), np.log10(2e4),
                      np.log10(25.0), np.log10(5e-3), 45.0, np.log10(3.0)])
    p0 = start + 0.02 * rng.standard_normal((nw, NDIM))
    p0[:, 5] = np.abs(p0[:, 5])
    print("sampling (5 films, 2 temperatures, 16 parameters) ...")
    chain, ar = stretch_sample(log_posterior, p0, nsteps, rng)
    flat = chain[burn:].reshape(-1, NDIM)
    print(f"acceptance {ar:.2f}, {len(flat)} samples")
    h1, h2 = chain[burn:4000], chain[4000:]
    for k in (4, 5, 14):
        print(f"  split-half shift in param {k}: "
              f"{abs(h1[...,k].mean()-h2[...,k].mean())/chain[burn:][...,k].std():.2f} sd")

    print("\nposterior (median [90 % CI])")
    lo, md, hi = q(flat[:, 14]); print(f"  activation energy       {md:6.1f}  [{lo:.0f}, {hi:.0f}] kJ/mol")
    lo, md, hi = q(flat[:, 5]);  print(f"  between-film spread     {md:6.2f}  [{lo:.2f}, {hi:.2f}] decades")
    lo, md, hi = q(flat[:, 3]);  print(f"  cracking onset d_crit   {md:6.1f}  [{lo:.1f}, {hi:.1f}] nm")
    lo, md, hi = q(flat[:, 1], lambda x: 10 ** x)
    print(f"  closure length d_close  {md:6.2f}  [{lo:.2f}, {hi:.2f}] nm")
    for i, f in enumerate(FILMS):
        lo, md, hi = q(flat[:, 6 + i], lambda x: 10 ** x)
        print(f"  f_res {f:11s}      {md:.2e}  [{lo:.1e}, {hi:.1e}]")
    for name, k in (("pV3D3", 12), ("plasma polymer", 15)):
        lo, md, hi = q(flat[:, k], lambda x: 10 ** x)
        print(f"  {name:15s} perm.  {md:6.1f}  [{lo:.1f}, {hi:.1f}] x parylene")

    rn = rng.standard_normal(len(flat))
    lf_new = flat[:, 4] + flat[:, 5] * rn
    lo, md, hi = q(lf_new, lambda x: 10 ** x)
    print(f"\n  new unmeasured film     {md:.2e}  [{lo:.1e}, {hi:.1e}]")
    d_close = 10 ** flat[:, 1] * np.log(10 ** flat[:, 0] / (CLOSURE_EPS * 10 ** lf_new))
    lo, md, hi = q(d_close); print(f"  window lower edge       {md:6.1f}  [{lo:.1f}, {hi:.1f}] nm")
    lo, md, hi = q(flat[:, 3]); print(f"  window upper edge       {md:6.1f}  [{lo:.1f}, {hi:.1f}] nm")
    w = flat[:, 3] - d_close
    lo, md, hi = q(w); print(f"  window width            {md:6.1f}  [{lo:.1f}, {hi:.1f}] nm"
                             f"   (empty in {100*np.mean(w<=0):.1f} % of draws)")

    # residuals, to see whether one shared shape really fits five films
    print("\n  worst residuals (log10 units) at the posterior median")
    pm = np.median(flat, axis=0)
    f0, dc, C, d_crit = 10**pm[0], 10**pm[1], 10**pm[2], pm[3]
    K, Ea = 10**pm[11], pm[14]
    fres = 10**pm[6:6+NF]
    mult = {"pV3D3": 10**pm[12], "PP": 10**pm[15]}
    res = []
    for film, d, val, s_, T, cens in SINGLE:
        m = K*arr(T, Ea)*f_shape(d, f0, dc, C, d_crit, fres[IDX[film]])
        res.append((abs(np.log10(m/val)), f"{film} single {d:.0f} nm @{T:.0f}C", np.log10(m/val)))
    for film, n, d_in, d_org, d_first, val, s_, da, T in LADDER:
        P = mult.get(MATERIAL.get(film, ""), 1.0)*P_PAR*arr(T, Ea)
        f = f_shape(d_in, f0, dc, C, d_crit, fres[IDX[film]])
        m = ladder_wvtr(n, d_in, d_org, d_first, f, P, da, film == "Lee", 10**pm[13])
        res.append((abs(np.log10(m/val)), f"{film} ladder n={n} @{T:.0f}C", np.log10(m/val)))
    for a_, lab, r in sorted(res, reverse=True)[:5]:
        print(f"    {lab:28s} {r:+.2f}")
    np.save("hier_v2_posterior.npy", flat)
    print("\nsaved hier_v2_posterior.npy")
