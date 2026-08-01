# -*- coding: utf-8 -*-
"""
step6_optimize.py — 6-objective Pareto exploration of TFE geometry
(design vars: d_org, d_inorg, n_pairs) on the CALIBRATED system
(ALD Al2O3 / parylene C, Step-5 parameters read from engine_parameters.csv).

Method
------
* Fast evaluator: analytic resistance ladder (machine-precision-equal to the
  PDE engine at steady state, verified in Steps 2-5) + exact discrete
  Frisch time-lag  t_lag = sum_j C_j*R_L,j*R_R,j / R_tot  (reduces to
  L^2/6D for a uniform slab), + peak/night radiative-convective energy
  balance, + incoherent optics, + deposition-time cost proxy (measured
  rates), + Basquin fatigue with thickness-dependent channel-cracking
  strength sigma_c(d) = sigma_c0*sqrt(d_ref/d).
* Optimizer: self-implemented NSGA-II (no external deps).
* Rigor check: exhaustive grid enumeration (feasible in 3-D) -> true
  Pareto front; NSGA-II convergence quantified by IGD.
Environment for lifetime: 38 C / 90 %RH (calibration-adjacent; Arrhenius
rescaling to 85/85 shifts absolute T80, not geometry ranking).
"""
import numpy as np, csv, os

RG, TREF = 8.314, 298.15
OUT = "/home/claude"

# ---------------- calibrated parameters from the DB --------------------
def load_params(path="/mnt/user-data/outputs/step1_db/engine_parameters.csv"):
    P = {}
    for r in csv.DictReader(open(path)):
        P[(r["engine_symbol"], r["material_id"])] = float(r["chosen_value"])
    return P

PB = load_params()
f0      = PB[("f0", "AL2O3_ALD")]
d_close = PB[("d_close_nm", "AL2O3_ALD")]
f_res   = PB[("f_res", "AL2O3_ALD")]
d_crit  = PB[("d_crit_nm", "AL2O3_ALD")]
f_c0    = 1.51e-6   # = C*f_res, corrected anchor
r_pin   = PB[("r_pinhole", "AL2O3_ALD")]
D_par, S_par = PB[("D_ref", "PARYLENE_C_CVD")], PB[("S", "PARYLENE_C_CVD")]
Ea_par  = 40e3
D_lat, Ea_lat, S_in = 1e-21, 60e3, 0.10
M_CRIT  = PB[("M_crit", "PSC_DEVICE")]          # g/m2 (placeholder, H)
RATE_ALD, RATE_PAR = 0.098, 6.0                  # nm/min (measured, WU2018)

T_ENV, RH_ENV = 38 + 273.15, 0.90               # lifetime environment
arr = lambda x, Ea, T: x * np.exp(-Ea / RG * (1 / T - 1 / TREF))

def f_pin(d_nm):
    return (f0 * np.exp(-d_nm / d_close) + f_res
            + f_c0 * max(0.0, (d_nm - d_crit) / d_crit) ** 2)

# ---------------- fast 6-objective evaluator ----------------------------
RHO_IN, RHO_ORG = 3000.0, 1289.0
N_AIR, N_ORG, N_IN, N_DEV = 1.0, 1.639, 1.65, 2.4
AL_ORG, AL_IN = 5e2, 1e3                        # 1/m absorption
E_IN, NU_IN, CTE_IN, CTE_SUB = 150e9, 0.24, 5e-6, 30e-6
SIG_C0, D_REF_SIG, M_FAT = 350e6, 30.0, 8.0
H_TOP, H_BOT, EPS_TOP, ETA_H, G_PEAK = 12.0, 6.0, 0.88, 0.75, 1000.0
SB = 5.670e-8

def transmittance(d_org, d_in, n):
    seq = [N_AIR] + [N_ORG, N_IN] * n + [N_DEV]
    T = 1.0
    for a, b in zip(seq[:-1], seq[1:]):
        T *= 1 - ((a - b) / (a + b)) ** 2
    return T * np.exp(-(AL_ORG * d_org + AL_IN * d_in) * 1e-9 * n)

def solve_T(Ta, G, Topt):
    T = Ta + 5.0
    Tsky = Ta - 15.0
    for _ in range(60):
        F = G * Topt * ETA_H + EPS_TOP * SB * (Tsky**4 - T**4) - (H_TOP + H_BOT) * (T - Ta)
        dF = -4 * EPS_TOP * SB * T**3 - (H_TOP + H_BOT)
        T -= F / dF
    return T

def moisture(d_org, d_in, n, T=T_ENV, da=RH_ENV):
    """Resistance ladder + discrete Frisch lag on the TFE cell chain
    (device sink | [inorg/org] x n | ambient)."""
    P_org = arr(D_par, Ea_par, T) * S_par
    P_in = arr(D_lat, Ea_lat, T) * S_in + f_pin(d_in) * P_org
    fv = f_pin(d_in)
    s = r_pin * np.sqrt(np.pi / fv)
    t2 = 1 + s**2 * np.log(max(s / r_pin, np.e)) / (2 * np.pi * (d_org * 1e-9) ** 2)
    r_cells, c_cells = [], []
    for i in range(n):
        r_cells.append(d_in * 1e-9 / P_in);  c_cells.append(S_in * d_in * 1e-9)
        tt = 1.0 if i == n - 1 else t2
        r_cells.append(d_org * 1e-9 * tt / P_org); c_cells.append(S_par * d_org * 1e-9)
    r = np.array(r_cells); c = np.array(c_cells)
    R = r.sum()
    RL = np.cumsum(r) - r / 2                    # inlet(device) -> cell centre
    RR = R - RL                                  # cell centre -> ambient
    t_lag = float((c * RL * RR).sum() / R)       # s
    J = da / R                                   # kg/m2/s
    wvtr = J * 1e3 * 86400                       # g/m2/day
    t80_h = t_lag / 3600 + M_CRIT / wvtr * 24
    return wvtr, t80_h, t_lag / 3600

def durability_years(d_in, dT):
    sig_a = E_IN / (1 - NU_IN) * abs(CTE_IN - CTE_SUB) * dT / 2
    sig_c = SIG_C0 * np.sqrt(D_REF_SIG / d_in)
    yrs = (sig_c / sig_a) ** M_FAT / 365.0
    if d_in > d_crit:
        yrs *= 0.1                               # pre-crack risk (strained regime)
    return min(yrs, 100.0)

def evaluate(x):
    d_org, d_in, n = float(x[0]), float(x[1]), int(round(x[2]))
    Topt = transmittance(d_org, d_in, n)
    wvtr, t80, tlag = moisture(d_org, d_in, n)
    Tpk = solve_T(42 + 273.15, G_PEAK, Topt)
    Tng = solve_T(15 + 273.15, 0.0, Topt)
    dT = Tpk - Tng
    dur = durability_years(d_in, dT)
    cost = n * (d_in / RATE_ALD + d_org / RATE_PAR)          # minutes
    wt = n * (RHO_IN * d_in + RHO_ORG * d_org) * 1e-9 * 1e3  # g/m2
    kpi = dict(d_org=d_org, d_in=d_in, n=n, T80_h=t80, wvtr=wvtr, t_lag_h=tlag,
               Tmax_C=Tpk - 273.15, Topt=Topt, cost_min=cost, weight=wt, dur_yr=dur)
    F = np.array([-t80, Tpk - 273.15, -Topt, cost, wt, -dur])  # all minimized
    return F, kpi

BOUNDS = np.array([[100.0, 1000.0], [15.0, 120.0], [1.0, 6.0]])  # d_org, d_in, n

# ---------------- exhaustive grid -> true Pareto front ------------------
def pareto_mask(F):
    N = len(F); mask = np.ones(N, bool)
    for i in range(N):
        if not mask[i]: continue
        d = np.all(F <= F[i], axis=1) & np.any(F < F[i], axis=1)
        if d.any(): mask[i] = False
    return mask

def run_grid():
    xs = [(o, i, n) for o in np.arange(100, 1001, 50)
          for i in np.arange(15, 121, 3) for n in range(1, 7)]
    F = np.empty((len(xs), 6)); K = []
    for j, x in enumerate(xs):
        F[j], k = evaluate(x); K.append(k)
    m = pareto_mask(F)
    return np.array(xs), F, K, m

# ---------------- self-implemented NSGA-II ------------------------------
rng = np.random.default_rng(42)

def nds(F):
    N = len(F); S = [[] for _ in range(N)]; nn = np.zeros(N, int)
    fronts = [[]]; rank = np.zeros(N, int)
    for p in range(N):
        for q in range(N):
            if p == q: continue
            if np.all(F[p] <= F[q]) and np.any(F[p] < F[q]): S[p].append(q)
            elif np.all(F[q] <= F[p]) and np.any(F[q] < F[p]): nn[p] += 1
        if nn[p] == 0: rank[p] = 0; fronts[0].append(p)
    i = 0
    while fronts[i]:
        nxt = []
        for p in fronts[i]:
            for q in S[p]:
                nn[q] -= 1
                if nn[q] == 0: rank[q] = i + 1; nxt.append(q)
        i += 1; fronts.append(nxt)
    return fronts[:-1], rank

def crowding(F, idx):
    n = len(idx); d = np.zeros(n)
    for m in range(F.shape[1]):
        o = np.argsort(F[idx, m]); fm = F[idx, m][o]
        rng_m = fm[-1] - fm[0]
        d[o[0]] = d[o[-1]] = np.inf
        if rng_m > 0:
            d[o[1:-1]] += (fm[2:] - fm[:-2]) / rng_m
    return d

def sbx(a, b, eta=15.0):
    u = rng.random(len(a)); beta = np.where(u <= .5, (2*u)**(1/(eta+1)), (1/(2-2*u))**(1/(eta+1)))
    return .5*((1+beta)*a+(1-beta)*b), .5*((1-beta)*a+(1+beta)*b)

def pmut(x, eta=20.0, p=1/3):
    for k in range(len(x)):
        if rng.random() < p:
            u = rng.random()
            dlt = (2*u)**(1/(eta+1))-1 if u < .5 else 1-(2-2*u)**(1/(eta+1))
            x[k] = np.clip(x[k]+dlt, 0, 1)
    return x

def denorm(u): return BOUNDS[:, 0] + u * (BOUNDS[:, 1] - BOUNDS[:, 0])

def nsga2(pop=80, gen=80):
    X = rng.random((pop, 3))
    F = np.array([evaluate(denorm(x))[0] for x in X])
    for g in range(gen):
        fronts, rank = nds(F)
        cd = np.zeros(pop)
        for fr in fronts: cd[fr] = crowding(F, fr)
        def tour():
            i, j = rng.integers(0, pop, 2)
            return i if (rank[i], -cd[i]) < (rank[j], -cd[j]) else j
        kids = []
        while len(kids) < pop:
            a, b = X[tour()].copy(), X[tour()].copy()
            if rng.random() < .9: a, b = sbx(a, b)
            kids += [pmut(a), pmut(b)]
        Xk = np.clip(np.array(kids[:pop]), 0, 1)
        Fk = np.array([evaluate(denorm(x))[0] for x in Xk])
        XA, FA = np.vstack([X, Xk]), np.vstack([F, Fk])
        fronts, _ = nds(FA)
        sel = []
        for fr in fronts:
            if len(sel) + len(fr) <= pop: sel += fr
            else:
                cdf = crowding(FA, fr)
                sel += [fr[i] for i in np.argsort(-cdf)[:pop - len(sel)]]; break
        X, F = XA[sel], FA[sel]
    fronts, _ = nds(F)
    return X[fronts[0]], F[fronts[0]]

# ---------------- run ----------------------------------------------------
if __name__ == "__main__":
    import time
    t0 = time.time()
    Xg, Fg, Kg, mg = run_grid()
    print("grid: %d designs, true Pareto front size %d  (%.1f s)"
          % (len(Xg), mg.sum(), time.time() - t0))
    t0 = time.time()
    Xn, Fn = nsga2()
    print("NSGA-II: front size %d  (%.1f s)" % (len(Fn), time.time() - t0))

    # IGD of NSGA front vs true grid front (normalized)
    Ft = Fg[mg]
    lo, hi = Ft.min(0), Ft.max(0); sc = np.where(hi > lo, hi - lo, 1)
    d = [np.min(np.linalg.norm((Fn - t) / sc, axis=1)) for t in Ft]
    print("convergence IGD (norm.): mean %.3f  max %.3f" % (np.mean(d), np.max(d)))

    # objective spans on the true front (effective dimensionality)
    names = ["-T80_h", "Tmax_C", "-Topt", "cost_min", "weight", "-dur_yr"]
    print("\ntrue-front objective spans:")
    for k, nm in enumerate(names):
        print("  %-9s %.4g .. %.4g" % (nm, Ft[:, k].min(), Ft[:, k].max()))

    # representative designs
    Kt = [Kg[i] for i in np.where(mg)[0]]
    t80s = np.array([k["T80_h"] for k in Kt])
    Fn_norm = (Ft - lo) / sc
    knee = int(np.argmin(np.linalg.norm(Fn_norm, axis=1)))
    best = int(np.argmax(t80s))
    cheap = int(np.argmin([k["cost_min"] if k["T80_h"] > 0.3 * t80s.max() else 1e18 for k in Kt]))
    print("\nrepresentative designs (true front):")
    for tag, i in (("MAX-T80", best), ("KNEE", knee), ("ECONOMY(T80>30%max)", cheap)):
        k = Kt[i]
        print("  %-20s d_org=%4.0f d_in=%3.0f n=%d | T80=%6.0f h  Topt=%.3f  "
              "cost=%5.0f min  wt=%4.2f g/m2  dur=%5.1f yr  Tmax=%.1f C"
              % (tag, k["d_org"], k["d_in"], k["n"], k["T80_h"], k["Topt"],
                 k["cost_min"], k["weight"], k["dur_yr"], k["Tmax_C"]))

    # save datasets (Step-7 SHAP input)
    with open(os.path.join(OUT, "grid_all.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(Kg[0].keys()) + ["pareto"])
        w.writeheader()
        for k, m in zip(Kg, mg): w.writerow({**k, "pareto": int(m)})
    with open(os.path.join(OUT, "pareto_front.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(Kt[0].keys())); w.writeheader()
        for k in Kt: w.writerow(k)

    # ---------------- figures ------------------------------------------
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Kg_arr = {k: np.array([d[k] for d in Kg]) for k in Kg[0]}
    Kt_arr = {k: np.array([d[k] for d in Kt]) for k in Kt[0]}
    fig, ax = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    fig.suptitle("Step 6 — 6-objective Pareto exploration (calibrated Al$_2$O$_3$/parylene C system)",
                 fontsize=13)

    a = ax[0, 0]
    a.scatter(Kg_arr["cost_min"], Kg_arr["T80_h"], s=6, c="lightgrey", label="all designs")
    sc1 = a.scatter(Kt_arr["cost_min"], Kt_arr["T80_h"], c=Kt_arr["n"], cmap="viridis",
                    s=40, label="true Pareto (grid)")
    nk = [evaluate(denorm(x))[1] for x in Xn]
    a.scatter([k["cost_min"] for k in nk], [k["T80_h"] for k in nk],
              marker="x", c="red", s=30, label="NSGA-II front")
    plt.colorbar(sc1, ax=a, label="n_pairs")
    a.set(xlabel="cost proxy [deposition min]", ylabel="T80 [h] @38C/90%RH",
          title="(a) Lifetime vs cost — n_pairs ladder")
    a.legend(fontsize=8); a.grid(alpha=.3)

    b = ax[0, 1]
    sc2 = b.scatter(Kt_arr["weight"], Kt_arr["T80_h"], c=Kt_arr["d_in"], cmap="plasma", s=40)
    plt.colorbar(sc2, ax=b, label="d_inorg [nm]")
    b.set(xlabel="added weight [g/m$^2$]", ylabel="T80 [h]",
          title="(b) Lifetime vs weight")
    b.grid(alpha=.3)

    c = ax[1, 0]
    sc3 = c.scatter(Kt_arr["dur_yr"], Kt_arr["T80_h"], c=Kt_arr["d_in"], cmap="plasma", s=40)
    plt.colorbar(sc3, ax=c, label="d_inorg [nm]")
    c.set(xlabel="fatigue durability [yr]", ylabel="T80 [h]",
          title="(c) Moisture lifetime vs thermo-mech durability")
    c.grid(alpha=.3)

    d_ = ax[1, 1]
    cols = ["T80_h", "Topt", "cost_min", "weight", "dur_yr"]
    M = np.array([[k[c2] for c2 in cols] for k in Kt], float)
    Mn = (M - M.min(0)) / np.where(M.max(0) > M.min(0), M.max(0) - M.min(0), 1)
    for row, nn_ in zip(Mn, Kt_arr["n"]):
        d_.plot(range(len(cols)), row, alpha=.5,
                color=plt.cm.viridis((nn_ - 1) / 5))
    d_.set_xticks(range(len(cols)))
    d_.set_xticklabels(["T80", "T_opt", "cost", "weight", "dur"], fontsize=9)
    d_.set(title="(d) Pareto set, parallel coordinates (color = n_pairs)",
           ylabel="normalized")
    d_.grid(alpha=.3)
    fig.savefig(os.path.join(OUT, "fig_step6_pareto.png"), dpi=150)
    print("\nsaved: grid_all.csv, pareto_front.csv, fig_step6_pareto.png")
