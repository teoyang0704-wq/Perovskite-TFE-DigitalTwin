# -*- coding: utf-8 -*-
"""si_compute.py — computations backing Supporting Information (Figs. S1–S3, Tables)."""
import numpy as np, time, sys, os, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

sys.path.insert(0, "/home/claude")
OUT = "/mnt/user-data/outputs/SI"
os.makedirs(OUT, exist_ok=True)
C = dict(blue="#0072B2", green="#009E73", orange="#E69F00", red="#D55E00", grey="#666")
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "savefig.bbox": "tight"})
def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=300); fig.savefig(f"{OUT}/{name}.pdf"); plt.close(fig)

env = {"python": sys.version.split()[0]}
for m in ("numpy", "scipy", "pandas", "matplotlib"):
    env[m] = __import__(m).__version__
runt = {}

# ---------- Fig. S1: identifiability ridge (profile over d_c) -----------
d_tab = np.array([15., 20., 30., 50., 60.])
W_tab = np.array([6.7e-3, 7.0e-4, 8.0e-4, 1.3e-3, 4.7e-3])
def fshape(d, A, dc, Cc, d0):
    return A*np.exp(-d/dc) + 1.0 + Cc*np.maximum(0.0, (d-d0)/d0)**2
def profile(dc):
    def r(p):
        lK, lA, lC, d0 = p
        return np.log10(10**lK*fshape(d_tab, 10**lA, dc, 10**lC, d0)) - np.log10(W_tab)
    fit = least_squares(r, [np.log10(7e-4), 4, np.log10(3), 40],
                        bounds=([-6, 0, -2, 20], [0, 12, 4, 60]))
    lK, lA, lC, d0 = fit.x
    f15_over_floor = fshape(15, 10**lA, dc, 10**lC, d0)  # /f_res-normalized
    return 2*fit.cost, 10**lA, f15_over_floor, d0
dcs = np.linspace(1.0, 15.0, 29)
prof = np.array([profile(dc) for dc in dcs])
print("S1 profile: SSR range %.3e–%.3e (decades^2); f(15)/f_floor range %.2f–%.2f; d_crit range %.1f–%.1f"
      % (prof[:,0].min(), prof[:,0].max(), prof[:,2].min(), prof[:,2].max(),
         prof[:,3].min(), prof[:,3].max()))
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.6))
ax[0].plot(dcs, prof[:, 0], color=C["blue"])
ax[0].set(xlabel="fixed $d_c$ [nm]", ylabel="profile SSR [decade$^2$]",
          title="(a) Flat likelihood valley: ($f_0$,$d_c$) ridge")
ax0b = ax[0].twinx(); ax0b.semilogy(dcs, prof[:, 1], color=C["orange"], ls="--")
ax0b.set_ylabel("compensating $f_0/f_{res}$", color=C["orange"])
ax0b.spines["right"].set_visible(True)
ax[1].plot(dcs, prof[:, 2], color=C["green"])
ax[1].set(xlabel="fixed $d_c$ [nm]", ylabel="$f(15\\,\\mathrm{nm})/f_{res}$",
          title="(b) What IS identified: $f$(15 nm)")
save(fig, "FigS1_identifiability_ridge")

# ---------- Fig. S3: test-duration sensitivity ---------------------------
t0 = time.time(); import step5_calibrate as S5; runt["step5_import_s"] = round(time.time()-t0, 1)
fig, ax = plt.subplots(figsize=(4.8, 3.2))
band = {}
tgrid = np.linspace(24, 240, 200) * 3600
errs = {}
for (n, ts, Jt, Jss, tlag, ap), col in zip(S5.rows[1:], (C["green"], C["red"])):
    e = np.log10(np.interp(tgrid, ts, Jt) / S5.dyads_meas[n])
    errs[n] = e
    ax.plot(tgrid/86400, e, color=col, label=f"{n} dyad")
emax = np.maximum(np.abs(errs[2]), np.abs(errs[3]))
ok = tgrid[emax <= 0.4] / 86400
band = (float(ok.min()), float(ok.max())) if len(ok) else (np.nan, np.nan)
print("S3 duration window with max|err|<=0.4 dec:", f"{band[0]:.1f}–{band[1]:.1f} days")
for d in (1,2,3,5,7):
    i2 = np.argmin(np.abs(tgrid-d*86400))
    print(f"   t={d} d: err2={errs[2][i2]:+.2f}, err3={errs[3][i2]:+.2f}")
ax.axhspan(-0.4, 0.4, color="grey", alpha=.15)
ax.axvspan(band[0], band[1], color=C["blue"], alpha=.12)
ax.axvline(3, ls="--", c="k", lw=.8)
ax.set(xlabel="assumed test duration [days]", ylabel="prediction error [decades]",
       title=f"Duration sensitivity: |err|≤0.4 dec for {band[0]:.1f}–{band[1]:.1f} d")
ax.legend(); save(fig, "FigS3_duration_sensitivity")

# ---------- Fig. S2: MC convergence --------------------------------------
t0 = time.time(); import step8_uq as S8; runt["step8_MC500_s"] = round(time.time()-t0, 1)
Ns = np.arange(25, 501, 25)
keys = ("closure", "crack", "a_rel")
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4))
conv = {}
for a, k in zip(ax, keys):
    vals = np.array([r[k] for r in S8.res])
    med = [np.nanmedian(vals[:n]) for n in Ns]
    lo = [np.nanpercentile(vals[:n], 2.5) for n in Ns]
    hi = [np.nanpercentile(vals[:n], 97.5) for n in Ns]
    a.plot(Ns, med, color=C["blue"]); a.fill_between(Ns, lo, hi, alpha=.2, color=C["blue"])
    a.set(xlabel="N samples", title=k)
    conv[k] = abs(med[-1] - med[Ns.tolist().index(250)])
print("S2 convergence |median(500)-median(250)|:", {k: round(v, 3) for k, v in conv.items()})
save(fig, "FigS2_mc_convergence")

# ---------- Table S1: full exact Sobol ------------------------------------
import pandas as pd, step7_xai as S7
df = pd.read_csv("/mnt/user-data/outputs/step6_optimization/grid_all.csv")
rows = []
for o in S7.OBJS:
    Sb, ST, V = S7.exact_sobol(df, o)
    rows.append(dict(objective=o, V=V, **{k: round(v, 3) for k, v in Sb.items()},
                     **{f"ST_{k}": round(v, 3) for k, v in ST.items()}))
pd.DataFrame(rows).to_csv(f"{OUT}/TableS1_sobol_full.csv", index=False)
print("\nTable S1 (Sobol, first-order & interactions):")
print(pd.DataFrame(rows).to_string(index=False))

json.dump(dict(env=env, runtimes=runt, seeds=dict(nsga2=42, uq=7),
               duration_band_days=band), open(f"{OUT}/repro_meta.json", "w"), indent=1)
print("\nenv:", env, "\nruntimes:", runt)
