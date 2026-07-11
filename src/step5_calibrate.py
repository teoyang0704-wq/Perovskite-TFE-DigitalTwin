# -*- coding: utf-8 -*-
"""
step5_calibrate.py — Calibration & zero-refit validation against Wu et al.,
RSC Adv. 2018, 8, 5721 (Al2O3(ALD,60C)/parylene C dyads on PET, 38C/100%RH).

Stage A (shape) : fit f(d) = f_res*[A*exp(-d/dc) + 1 + C*((d-d0)/d0)_+^2]
                  to Table 1 single-layer WVTR(d) up to a multiplicative
                  constant  (assumption: bare-layer WVTR ∝ f(d); per-pinhole
                  conductance is access-limited, thickness-independent).
Stage B (scale) : set absolute f_res so the 1-dyad point matches exactly.
Validation      : predict 2- and 3-dyad WVTR with ZERO refit —
                  (i) steady-state resistance ladder,
                  (ii) full transient PDE with finite test duration
                       (MOCON runs may sit in the lag phase for high barriers).
"""
import numpy as np
from scipy.optimize import least_squares
import sys
sys.path.insert(0, '/home/claude')
from tfe_physics_engine import ImplicitDiffusion1D

R_GAS, TREF = 8.314, 298.15
T = 38 + 273.15           # MOCON Aquatran II condition
DA = 1.0                  # activity difference at 100%RH vs dry sensor
G2 = lambda J: J * 1e3 * 86400.0   # kg/m2/s -> g/m2/day
K2 = lambda W: W / 1e3 / 86400.0

def arr(x_ref, Ea, T):
    return x_ref * np.exp(-Ea / R_GAS * (1.0 / T - 1.0 / TREF))

# ---------------- materials (parylene C system of Wu 2018) ------------
# Parylene C: datasheet-typical permeability ~0.08 g*mm/(m2*day) -> ~9e-13
# kg/(m*s); split D/S placeholders (quality C, flagged in DB).
D_par, Ea_par, S_par = 5.0e-13, 40e3, 1.5
D_lat, Ea_lat, S_in = 1.0e-21, 60e3, 0.10
r_pin = 50e-9
P_par = arr(D_par, Ea_par, T) * S_par
P_lat = arr(D_lat, Ea_lat, T) * S_in
# PET from measured bare WVTR (BM015, quality A): 3.0 g/m2/day @38C
J_PET = K2(3.0)
R_PET = DA / J_PET
L_PET, S_PET = 125e-6, 5.0
P_PET = L_PET / R_PET

# ---------------- data (from DB) --------------------------------------
d_tab = np.array([15., 20., 30., 50., 60.])                 # Table 1 [nm]
W_tab = np.array([6.7e-3, 7.0e-4, 8.0e-4, 1.3e-3, 4.7e-3])  # g/m2/day
dyads_meas = {1: 1.7e-4, 2: 3.6e-5, 3: 7.7e-6}              # Fig.3 digitized
DIN, DORG = 50.0, 500.0  # CORRECTED 2026-07-10: Wu full text - Fig.3 dyad series uses 50 nm Al2O3 (primal thickness), parylene 500 nm

# ---------------- Stage A: shape fit -----------------------------------
def fshape(d, A, dc, C, d0):
    return A * np.exp(-d / dc) + 1.0 + C * np.maximum(0.0, (d - d0) / d0) ** 2

def resid(p):
    lK, lA, dc, lC, d0 = p
    return np.log10(10**lK * fshape(d_tab, 10**lA, dc, 10**lC, d0)) - np.log10(W_tab)

p0 = [np.log10(7e-4), np.log10(300), 4.0, np.log10(3), 35.0]
fit = least_squares(resid, p0, bounds=([-6, 0, 1.0, -2, 20], [0, 7, 15.0, 4, 60]))
lK, lA, dc, lC, d0 = fit.x
A, C = 10**lA, 10**lC
res_dec = resid(fit.x)
print("Stage A shape fit: A=%.3g dc=%.2f nm C=%.3g d_crit=%.1f nm" % (A, dc, C, d0))
print("  residuals [decades]:", np.round(res_dec, 3), " max|.|=%.3f" % np.max(np.abs(res_dec)))

# ---------------- Stage B: absolute scale from 1-dyad ------------------
R_top = DORG * 1e-9 / P_par
R1_needed = DA / K2(dyads_meas[1])
R_in_needed = R1_needed - R_PET - R_top
f30 = (DIN * 1e-9 / R_in_needed - P_lat) / P_par
f_res = f30 / fshape(DIN, A, dc, C, d0)
f = lambda d: f_res * fshape(d, A, dc, C, d0)
n_mm2 = f30 / (np.pi * r_pin**2) / 1e6
s30 = r_pin * np.sqrt(np.pi / f30)
print("Stage B: f(30nm)=%.3g  f_res=%.3g  -> pinhole density %.2f /mm2, spacing %.0f um"
      % (f30, f_res, n_mm2, s30 * 1e6))

# consistency: implied per-area bare-pinhole conductance vs open-vapor limit
g_bare30 = K2(np.interp(30, d_tab, W_tab)) / DA
P_vap = 2.6e-5 * 0.046
print("  bare-hole conductance / open-vapor limit = %.2g (should be <1, spreading-limited)"
      % (g_bare30 * (DIN * 1e-9) / (f30 * P_vap)))

# ---------------- forward model ---------------------------------------
def tau2(d_org_nm, fval):
    s = r_pin * np.sqrt(np.pi / fval)
    return 1.0 + s**2 * np.log(max(s / r_pin, np.e)) / (2 * np.pi * (d_org_nm * 1e-9)**2)

def R_ladder(n, d_in=DIN, d_org=DORG):
    fv = f(d_in)
    R_in = d_in * 1e-9 / (P_lat + fv * P_par)
    R_sand = d_org * 1e-9 * tau2(d_org, fv) / P_par
    return R_PET + n * R_in + (n - 1) * R_sand + R_top

print("\nZero-refit steady-state prediction (Fig.3 high-quality):")
print("  tau^2(sandwiched parylene) = %.3g" % tau2(DORG, f30))
pred_ss = {}
for n in (1, 2, 3):
    pred_ss[n] = G2(DA / R_ladder(n))
    err = np.log10(pred_ss[n] / dyads_meas[n])
    print("  %d dyad: model %.2e vs meas %.2e  (%+.2f decades)"
          % (n, pred_ss[n], dyads_meas[n], err))

# ---------------- transient PDE: finite-test apparent WVTR -------------
def build_cells(n, d_in=DIN, d_org=DORG):
    fv = f(d_in)
    dx, K, Cc = [], [], []
    npet = 24
    dx += [L_PET / npet] * npet; K += [P_PET] * npet; Cc += [S_PET] * npet
    for i in range(n):
        dx += [d_in * 1e-9 / 8] * 8
        K += [P_lat + fv * P_par] * 8; Cc += [S_in] * 8
        top = (i == n - 1)
        t2 = 1.0 if top else tau2(d_org, fv)
        dx += [d_org * 1e-9 / 8] * 8
        K += [P_par / t2] * 8; Cc += [S_par] * 8
    return map(np.array, (dx, K, Cc))

def transient(n, t_end_d=30.0):
    dx, K, Cc = build_cells(n)
    sol = ImplicitDiffusion1D(dx)
    a = np.zeros(len(dx)); t, dt, M = 0.0, 1.0, 0.0
    ts, Js = [], []
    while t < t_end_d * 86400:
        dtn = min(dt, t_end_d * 86400 - t)
        a, Jl, _ = sol.step(a, dtn, K, Cc, ("dirichlet", 0.0), ("dirichlet", DA))
        t += dtn; M += Jl * dtn
        ts.append(t); Js.append(Jl)
        dt = min(dt * 1.2, 900.0)
    ts, Js = np.array(ts), np.array(Js)
    Jss = DA / R_ladder(n)
    tlag = ts[-1] - M / Jss
    return ts, G2(Js), G2(Jss), tlag / 3600.0

print("\nTransient (finite MOCON test) — apparent WVTR:")
rows = []
for n in (1, 2, 3):
    ts, Jt, Jss, tlag = transient(n)
    ap = {h: float(np.interp(h * 3600, ts, Jt)) for h in (24, 72, 168)}
    rows.append((n, ts, Jt, Jss, tlag, ap))
    print("  %d dyad: t_lag=%.1f h | J(24h)=%.2e J(72h)=%.2e J(168h)=%.2e | J_ss=%.2e | meas=%.2e"
          % (n, tlag, ap[24], ap[72], ap[168], Jss, dyads_meas[n]))

# ---------------- figure ------------------------------------------------
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)
fig.suptitle("Step 5 — calibration (Table 1) and zero-refit validation (Fig. 3), Wu et al. RSC Adv. 2018", fontsize=12)

dd = np.linspace(12, 65, 200)
ax[0].semilogy(d_tab, W_tab, "ks", ms=8, label="Table 1 (meas)")
ax[0].semilogy(dd, 10**lK * fshape(dd, A, dc, C, d0), "b-", label="fitted f(d) shape")
ax[0].set(xlabel="Al$_2$O$_3$ thickness [nm]", ylabel="single-layer WVTR [g m$^{-2}$ day$^{-1}$]",
          title="(a) Stage A: closure + cracking shape")
ax[0].legend(); ax[0].grid(alpha=0.3, which="both")

ns = np.array([1, 2, 3])
ax[1].semilogy(ns, [dyads_meas[n] for n in ns], "ko", ms=9, label="measured (digitized)")
ax[1].semilogy(ns, [pred_ss[n] for n in ns], "b^-", label="model steady-state")
ax[1].semilogy(ns, [r[5][72] for r in rows], "rv--", label="model apparent @72 h")
ax[1].axhline(5e-5, color="grey", ls=":", label="stated instrument limit")
ax[1].plot([1], [dyads_meas[1]], "y*", ms=16, label="scale anchor (fitted)")
ax[1].set(xticks=[1, 2, 3], xlabel="number of dyads",
          ylabel="WVTR [g m$^{-2}$ day$^{-1}$]", title="(b) Zero-refit dyad prediction")
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3, which="both")

for (n, ts, Jt, Jss, tlag, ap), c in zip(rows, ("tab:blue", "tab:green", "tab:red")):
    ax[2].loglog(ts / 3600, Jt, color=c, label=f"{n} dyad (t_lag={tlag:.0f} h)")
    ax[2].axhline(Jss, color=c, ls=":", alpha=0.6)
    ax[2].plot(72, dyads_meas[n], "o", color=c, mfc="none", ms=9)
ax[2].axhline(5e-5, color="grey", ls=":")
ax[2].set(xlabel="test time [h]", ylabel="apparent WVTR [g m$^{-2}$ day$^{-1}$]",
          title="(c) Lag-phase effect (circles: measured @assumed 72 h)")
ax[2].legend(fontsize=8); ax[2].grid(alpha=0.3, which="both")
fig.savefig("/home/claude/fig_step5_calibration.png", dpi=150)
print("\nfigure saved; calibrated params: f0=%.3g d_close=%.2f nm f_res=%.3g f_crack0=%.3g d_crit=%.1f nm"
      % (A * f_res, dc, f_res, C * f_res, d0))
