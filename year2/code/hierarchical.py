# -*- coding: utf-8 -*-
"""Year-2, module 2: hierarchical calibration across films and laboratories.

Year-1 calibrated on one film in two point-estimate stages: a shape from a
single-layer thickness series, a scale from one multilayer anchor.  Everything
downstream then inherited that one film's cleanliness.  Module 1 showed the
between-film spread is far smaller than we thought (1.1x, not 15x, once the
sample edge is accounted for), which makes a proper hierarchical treatment
worthwhile: the defect floor of each film is drawn from a population, the
population is inferred from all films at once, and a laboratory that has never
measured anything can be given a predictive distribution rather than a number.

Model
-----
Defect area fraction of an inorganic layer of thickness d (nm):

    f(d) = f0 exp(-d/dc)  +  f_res  +  C f_res [(d - d_crit)/d_crit]_+^2
           ^ nucleation      ^ floor   ^ cracking

f0, dc, C, d_crit are properties of the ALD chemistry and are shared by every
film.  f_res is the particulate floor of one deposition campaign and is the
quantity that varies between films:

    log10 f_res,i  ~  Normal(log10 mu_f, sigma_f)          i = 1 .. n_films

sigma_f is the between-film spread in decades -- the number that decides how
far any calibrated rule can travel.

Observables
-----------
bare single layers      WVTR = K f(d)                (defect-limited, K shared)
multilayer ladders      series resistance network, plus the edge channel of
                        module 1 where the test geometry has an open perimeter

Sampling
--------
No MCMC package is available in this environment, so the affine-invariant
ensemble sampler (stretch move) is implemented here in about forty lines.  It
is the algorithm behind emcee; writing it out keeps every step auditable.

Author: Teo Yang.
"""
import numpy as np

DAY, KG_TO_G = 86400.0, 1e3
P_PAR = 7.5e-13          # parylene-C permeability [kg m^-1 s^-1 per activity]
R_SUB = 2.88e7           # PET substrate resistance [m^2 s kg^-1]
R_PIN = 50e-9            # pinhole radius [m]
FILMS = ["Wu high-q", "Wu poor", "Lee", "Carcia"]
NF = len(FILMS)

# ------------------------------------------------------------------ data
# (kind, film, geometry..., value, sigma_dex, driving activity)
SINGLE = [   # bare ALD on polymer: WVTR = K * f(d)
    ("Wu high-q", 15.0, 6.7e-3, 0.10), ("Wu high-q", 20.0, 0.7e-3, 0.10),
    ("Wu high-q", 30.0, 0.8e-3, 0.10), ("Wu high-q", 50.0, 1.3e-3, 0.10),
    ("Wu high-q", 60.0, 4.7e-3, 0.10),
    ("Carcia",    25.0, 1.7e-5, 0.15),
]
LADDER = [   # (film, n, d_in nm, d_org nm, d_first nm, value, sigma, da, edge?)
    ("Wu high-q", 1, 50.0, 500.0, 0.0, 1.70e-4, 0.10, 1.00, False),
    ("Wu poor",   1, 50.0, 500.0, 0.0, 1.6,     0.30, 1.00, False),
    ("Wu poor",   2, 50.0, 500.0, 0.0, 1.6e-1,  0.30, 1.00, False),
    ("Wu poor",   3, 50.0, 500.0, 0.0, 1.3e-1,  0.30, 1.00, False),
    ("Lee",       1, 21.5, 100.0, 800.0, 3.0e-3, 0.10, 0.90, True),
    ("Lee",       2, 21.5, 100.0, 800.0, 6.6e-4, 0.10, 0.90, True),
    ("Lee",       3, 21.5, 100.0, 800.0, 5.4e-4, 0.10, 0.90, True),
    ("Lee",       4, 21.5, 100.0, 800.0, 5.3e-4, 0.10, 0.90, True),
]
PAD, PERIM = 3e-3, 4 * 3e-3          # Lee calcium pad
CLOSURE_EPS = 0.01                   # closure complete when nucleation = 1 % of floor

# parameter vector
# 0 log10 f0 | 1 log10 dc | 2 log10 C | 3 d_crit | 4 log10 mu_f | 5 sigma_f
# 6..9 log10 f_res per film | 10 log10 K | 11 log10 m(pV3D3) | 12 log10 L_path
NDIM = 13
NAMES = ["log10 f0", "log10 dc", "log10 C", "d_crit", "log10 mu_f", "sigma_f",
         *[f"log10 f_res[{f}]" for f in FILMS], "log10 K", "log10 m", "log10 Lpath"]


def f_shape(d_nm, f0, dc, C, d_crit, f_res):
    crack = np.maximum(0.0, (d_nm - d_crit) / d_crit) ** 2
    return f0 * np.exp(-d_nm / dc) + f_res * (1.0 + C * crack)


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
    lf0, ldc, lC, d_crit, lmu, sig = p[0], p[1], p[2], p[3], p[4], p[5]
    lfres = p[6:6 + NF]
    lK, lm, lL = p[6 + NF], p[7 + NF], p[8 + NF]
    if not (0.0 < sig < 3.0 and 10.0 < d_crit < 90.0):
        return -np.inf
    f0, dc, C = 10 ** lf0, 10 ** ldc, 10 ** lC
    fres = 10 ** lfres
    K, m, L_path = 10 ** lK, 10 ** lm, 10 ** lL
    if not (0.05 < dc < 10.0 and 1e-10 < fres.min() and fres.max() < 1e-4):
        return -np.inf

    lp = 0.0
    # weakly informative priors on shared physics (Year-1 values, generous width)
    lp += -0.5 * ((lf0 - np.log10(0.372)) / 0.6) ** 2
    lp += -0.5 * ((ldc - np.log10(1.07)) / 0.35) ** 2
    lp += -0.5 * ((lC - np.log10(40.0)) / 0.8) ** 2
    lp += -0.5 * ((d_crit - 44.0) / 8.0) ** 2
    lp += -0.5 * ((lmu + 7.4) / 1.2) ** 2
    lp += -0.5 * (sig / 0.6) ** 2                      # half-normal on the spread
    lp += -0.5 * ((lm - np.log10(30.0)) / 0.5) ** 2
    lp += -0.5 * ((lL - np.log10(5e-3)) / 0.3) ** 2
    lp += -0.5 * ((lK - np.log10(2e4)) / 2.0) ** 2
    # hierarchical layer
    lp += np.sum(-0.5 * ((lfres - lmu) / sig) ** 2 - np.log(sig))

    idx = {f: i for i, f in enumerate(FILMS)}
    for film, d, val, s_ in SINGLE:
        mod = K * f_shape(d, f0, dc, C, d_crit, fres[idx[film]])
        lp += -0.5 * ((np.log10(mod) - np.log10(val)) / s_) ** 2
    for film, n, d_in, d_org, d_first, val, s_, da, edge in LADDER:
        P_org = m * P_PAR if film == "Lee" else P_PAR
        f = f_shape(d_in, f0, dc, C, d_crit, fres[idx[film]])
        mod = ladder_wvtr(n, d_in, d_org, d_first, f, P_org, da, edge, L_path)
        lp += -0.5 * ((np.log10(mod) - np.log10(val)) / s_) ** 2
    return lp if np.isfinite(lp) else -np.inf


# ------------------------------------------------------- ensemble sampler
def stretch_sample(logpost, p0, nsteps, rng, a=2.0):
    """Affine-invariant ensemble sampler (Goodman & Weare stretch move)."""
    nw, nd = p0.shape
    p = p0.copy()
    lp = np.array([logpost(x) for x in p])
    chain = np.empty((nsteps, nw, nd))
    acc = 0
    for t in range(nsteps):
        for half in (0, 1):
            me = np.where(np.arange(nw) % 2 == half)[0]
            you = np.where(np.arange(nw) % 2 != half)[0]
            for i in me:
                j = you[rng.integers(len(you))]
                u = rng.random()
                z = ((a - 1.0) * u + 1.0) ** 2 / a           # g(z) ~ 1/sqrt(z)
                q = p[j] + z * (p[i] - p[j])
                lq = logpost(q)
                if np.log(rng.random()) < (nd - 1) * np.log(z) + lq - lp[i]:
                    p[i], lp[i] = q, lq
                    acc += 1
        chain[t] = p
    return chain, acc / (nsteps * nw)


def summarise(flat, name, scale=lambda x: x, unit=""):
    v = scale(flat)
    lo, md, hi = np.percentile(v, [5, 50, 95])
    print(f"  {name:26s} {md:10.3g}  [{lo:.3g}, {hi:.3g}] {unit}")
    return md, lo, hi


if __name__ == "__main__":
    rng = np.random.default_rng(11)
    nw, nsteps, burn = 52, 6000, 2000
    start = np.array([np.log10(0.372), np.log10(1.07), np.log10(40.0), 44.0,
                      -7.4, 0.4, *([-7.4] * NF), np.log10(2e4),
                      np.log10(25.0), np.log10(5e-3)])
    p0 = start + 0.02 * rng.standard_normal((nw, NDIM))
    p0[:, 5] = np.abs(p0[:, 5])

    print("sampling ...")
    chain, arate = stretch_sample(log_posterior, p0, nsteps, rng)
    flat = chain[burn:].reshape(-1, NDIM)
    print(f"acceptance rate {arate:.2f}, {flat.shape[0]} posterior samples")

    # crude convergence check: split-half agreement of the population parameters
    h1, h2 = chain[burn:nsteps//2 + burn//2], chain[nsteps//2 + burn//2:]
    for k in (4, 5):
        m1, m2 = h1[..., k].mean(), h2[..., k].mean()
        sd = chain[burn:][..., k].std()
        print(f"  split-half shift in {NAMES[k]:12s} = {abs(m1-m2)/sd:.2f} sd")

    print("\nposterior (median [90 % credible interval])")
    summarise(flat[:, 0], "f0", lambda x: 10 ** x)
    summarise(flat[:, 1], "d_close", lambda x: 10 ** x, "nm")
    summarise(flat[:, 3], "d_crit", unit="nm")
    summarise(flat[:, 4], "population median f_res", lambda x: 10 ** x)
    sig_md, sig_lo, sig_hi = summarise(flat[:, 5], "between-film spread", unit="decades")
    for i, f in enumerate(FILMS):
        summarise(flat[:, 6 + i], f"f_res  {f}", lambda x: 10 ** x)
    summarise(flat[:, 6 + NF + 1], "pV3D3 permeability", lambda x: 10 ** x, "x parylene")

    # ---- posterior predictive for a film nobody has measured ----
    rn = rng.standard_normal(len(flat))
    lf_new = flat[:, 4] + flat[:, 5] * rn
    print("\nposterior predictive for a new, unmeasured film")
    summarise(lf_new, "f_res of a new film", lambda x: 10 ** x)
    dens = 10 ** lf_new / (np.pi * R_PIN ** 2) / 1e6
    summarise(np.log10(dens), "defect density", lambda x: 10 ** x, "per mm2")

    # ---- what the spread does to the design window ----
    d_close_new = 10 ** flat[:, 1] * np.log(10 ** flat[:, 0] /
                                            (CLOSURE_EPS * 10 ** lf_new))
    print("\ndesign window under between-film uncertainty")
    summarise(d_close_new, "lower edge (closure)", unit="nm")
    summarise(flat[:, 3], "upper edge (cracking)", unit="nm")
    width = flat[:, 3] - d_close_new
    summarise(width, "window width", unit="nm")
    print(f"  probability the window is empty for a new film: "
          f"{100*np.mean(width <= 0):.1f} %")
    print("\n  note: a cleaner film has a *lower* floor, so its nucleation term must"
          "\n  decay further before closure is complete -- cleaner films need thicker"
          "\n  inorganic layers, and their design window is narrower.")

    np.save("hier_posterior.npy", flat)
    print("\nsaved hier_posterior.npy")
