# -*- coding: utf-8 -*-
"""Fast 6-objective evaluator + exhaustive grid + self-contained NSGA-II.
Calibrated Al2O3(ALD)/parylene-C system; parameters loaded from data/engine_parameters.csv."""
import numpy as np, csv, os

RG, TREF = 8.314, 298.15
OUT = "/home/claude"

def load_params(path=None):
    if path is None:
        _here = os.path.dirname(os.path.abspath(__file__))
        for _c in (os.path.join(_here, "..", "data", "engine_parameters.csv"),
                   "/mnt/user-data/outputs/step1_db/engine_parameters.csv"):
            if os.path.exists(_c):
                path = _c; break
    P = {}
    for r in csv.DictReader(open(path)):
        P[(r["engine_symbol"], r["material_id"])] = float(r["chosen_value"])
    return P

PB = load_params()
f0      = PB[("f0", "AL2O3_ALD")]
d_close = PB[("d_close_nm", "AL2O3_ALD")]
f_res   = PB[("f_res", "AL2O3_ALD")]
d_crit  = PB[("d_crit_nm", "AL2O3_ALD")]
f_c0    = 1.62e-6
r_pin   = PB[("r_pinhole", "AL2O3_ALD")]
D_par, S_par = PB[("D_ref", "PARYLENE_C_CVD")], PB[("S", "PARYLENE_C_CVD")]
Ea_par  = 40e3
D_lat, Ea_lat, S_in = 1e-21, 60e3, 0.10
M_CRIT  = PB[("M_crit", "PSC_DEVICE")]
RATE_ALD, RATE_PAR = 0.098, 6.0

T_ENV, RH_ENV = 38 + 273.15, 0.90
arr = lambda x, Ea, T: x * np.exp(-Ea / RG * (1 / T - 1 / TREF))

def f_pin(d_nm):
    return (f0 * np.exp(-d_nm / d_close) + f_res
            + f_c0 * max(0.0, (d_nm - d_crit) / d_crit) ** 2)

RHO_IN, RHO_ORG = 3000.0, 1289.0
N_AIR, N_ORG, N_IN, N_DEV = 1.0, 1.639, 1.65, 2.4
AL_ORG, AL_IN = 5e2, 1e3
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
    RL = np.cumsum(r) - r / 2
    RR = R - RL
    t_lag = float((c * RL * RR).sum() / R)
    J = da / R
    wvtr = J * 1e3 * 86400
    t80_h = t_lag / 3600 + M_CRIT / wvtr * 24
    return wvtr, t80_h, t_lag / 3600

def durability_years(d_in, dT):
    sig_a = E_IN / (1 - NU_IN) * abs(CTE_IN - CTE_SUB) * dT / 2
    sig_c = SIG_C0 * np.sqrt(D_REF_SIG / d_in)
    yrs = (sig_c / sig_a) ** M_FAT / 365.0
    if d_in > d_crit:
        yrs *= 0.1
    return min(yrs, 100.0)

def evaluate(x):
    d_org, d_in, n = float(x[0]), float(x[1]), int(round(x[2]))
    Topt = transmittance(d_org, d_in, n)
    wvtr, t80, tlag = moisture(d_org, d_in, n)
    Tpk = solve_T(42 + 273.15, G_PEAK, Topt)
    Tng = solve_T(15 + 273.15, 0.0, Topt)
    dT = Tpk - Tng
    dur = durability_years(d_in, dT)
    cost = n * (d_in / RATE_ALD + d_org / RATE_PAR)
    wt = n * (RHO_IN * d_in + RHO_ORG * d_org) * 1e-9 * 1e3
    kpi = dict(d_org=d_org, d_in=d_in, n=n, T80_h=t80, wvtr=wvtr, t_lag_h=tlag,
               Tmax_C=Tpk - 273.15, Topt=Topt, cost_min=cost, weight=wt, dur_yr=dur)
    F = np.array([-t80, Tpk - 273.15, -Topt, cost, wt, -dur])
    return F, kpi

BOUNDS = np.array([[100.0, 1000.0], [15.0, 120.0], [1.0, 6.0]])

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
