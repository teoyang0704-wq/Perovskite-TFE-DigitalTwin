# -*- coding: utf-8 -*-
"""Interactive TFE Digital Twin - v2 (environment control + honesty band +
constraint optimizer). Main file for Streamlit Community Cloud."""
import os, sys, math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "optimization"))
from step6_optimize import (f_pin, arr, D_par, S_par, Ea_par, D_lat, Ea_lat,
                            S_in, r_pin, M_CRIT, d_crit, evaluate)

REPO = "https://github.com/teoyang0704-wq/Perovskite-TFE-DigitalTwin"
FLOOR = 97.0
G1, G2 = 22.5, 44.0
T_CAL = 311.15

def lifetime_env(d_org, d_in, n, T_C, RH, mfac):
    T, da = T_C + 273.15, RH / 100.0
    P_org = arr(D_par, Ea_par, T) * S_par
    P_lat = arr(D_lat, Ea_lat, T) * S_in
    fv = f_pin(d_in)
    P_in = P_lat + fv * P_org
    s = r_pin * math.sqrt(math.pi / fv)
    t2 = 1 + s**2 * math.log(max(s / r_pin, math.e)) / (2 * math.pi * (d_org * 1e-9) ** 2)
    r_l, c_l = [], []
    for i in range(int(n)):
        r_l += [d_in * 1e-9 / P_in]; c_l += [S_in * d_in * 1e-9]
        tt = 1.0 if i == n - 1 else t2
        r_l += [d_org * 1e-9 * tt / P_org]; c_l += [S_par * d_org * 1e-9]
    R = sum(r_l); RL, tl = 0.0, 0.0
    for r_i, c_i in zip(r_l, c_l):
        a, b = RL, R - RL
        tl += c_i / r_i * (a * b * r_i + (b - a) * r_i**2 / 2 - r_i**3 / 3)
        RL += r_i
    t_lag_h = tl / R / 3600.0
    wvtr = da / R * 1e3 * 86400.0
    t80 = t_lag_h + (M_CRIT * mfac) / wvtr * 24.0
    return t80, wvtr, t_lag_h

def ea_band(T_C):
    return math.exp(10e3 / 8.314 * abs(1 / (T_C + 273.15) - 1 / T_CAL))

st.set_page_config(page_title="TFE Digital Twin", page_icon="🛡️", layout="wide")
st.title("🛡️ Perovskite-TFE Digital Twin - design & environment explorer")
st.caption(f"Literature-calibrated, defect-mediated 1-D twin | zero-refit validated "
           f"(±0.25 decade) | [paper & code]({REPO}) | research prototype")

DEFAULTS = dict(d_org=100, d_in=30, n=3, T_C=38, RH=90)
for k, v in DEFAULTS.items():
    st.session_state.setdefault(k, v)

with st.sidebar:
    st.header("Design")
    st.slider("Organic interlayer d_org [nm]", 100, 1000, step=25, key="d_org")
    st.slider("Inorganic barrier d_in [nm]", 15, 120, step=1, key="d_in")
    st.slider("Dyads n", 1, 6, key="n")
    st.header("Environment")
    st.slider("Temperature [degC]", 15, 85, key="T_C")
    st.slider("Relative humidity [%]", 30, 100, key="RH")
    pc = st.columns(3)
    def preset(t, rh): st.session_state.update(T_C=t, RH=rh)
    pc[0].button("38/90", on_click=preset, args=(38, 90), use_container_width=True)
    pc[1].button("85/85", on_click=preset, args=(85, 85), use_container_width=True)
    pc[2].button("25/50", on_click=preset, args=(25, 50), use_container_width=True)
    apply_floor = st.toggle("Device-anchored M_crit floor (x97)", value=True)
    st.markdown(f"---\n**Design window (G1-G2):** {G1}-{G2} nm")

d_org, d_in, n = st.session_state.d_org, st.session_state.d_in, st.session_state.n
T_C, RH = st.session_state.T_C, st.session_state.RH
mfac = FLOOR if apply_floor else 1.0
t80, wvtr, tlag = lifetime_env(d_org, d_in, n, T_C, RH, mfac)
_, k38 = evaluate((float(d_org), float(d_in), int(n)))

if d_in < G1:
    st.error(f"G1 violated: d_in < {G1} nm - pinhole closure incomplete.")
elif d_in > G2:
    st.error(f"G2 violated: d_in > {G2} nm - cracking regime; lifetime AND durability collapse.")
else:
    st.success(f"Inside the calibrated window ({G1}-{G2} nm).")

b = ea_band(T_C)
c = st.columns(4)
c[0].metric(f"Lifetime @ {T_C} degC / {RH} %RH", f"{t80:,.0f} h", f"= {t80/8760:,.1f} yr")
c[1].metric("Extrapolation uncertainty", f"x/ {b:.2f}",
            "calibration-adjacent" if b < 1.15 else "Arrhenius-extrapolated",
            delta_color="off")
t85, _, _ = lifetime_env(d_org, d_in, n, 85, 85, mfac)
c[2].metric("Same design @ 85/85", f"{t85:,.0f} h", f"= {t85/8760:,.1f} yr", delta_color="off")
c[3].metric("Steady WVTR here", f"{wvtr:.2e} g/m2/day", f"lag {tlag:.1f} h", delta_color="off")
c2 = st.columns(4)
c2[0].metric("Fatigue durability (diurnal proxy)", f"{k38['dur_yr']:.1f} yr",
             "pre-crack penalty!" if d_in > d_crit else "chamber-independent", delta_color="off")
c2[1].metric("Cost proxy", f"{k38['cost_min']:,.0f} min")
c2[2].metric("Weight", f"{k38['weight']:.2f} g/m2")
c2[3].metric("T_opt / T_max", f"{k38['Topt']*100:.1f} % | {k38['Tmax_C']:.1f} degC",
             "geometry-flat", delta_color="off")
if RH > 95:
    st.info("RH > 95%: Henry-linear sorption may underestimate uptake - paper Sec. 3.5(vii).")

@st.cache_data
def load_grid():
    return pd.read_csv(os.path.join(ROOT, "optimization", "grid_all.csv"))

df = load_grid()
env_fac = (t80 - tlag) / max((k38["T80_h"] - k38["t_lag_h"]) * mfac, 1e-9)
df_t80 = df["t_lag_h"] + (df["T80_h"] - df["t_lag_h"]) * mfac * env_fac

left, right = st.columns([1.5, 1])
with left:
    st.subheader(f"4,104-design map @ {T_C} degC / {RH} %RH")
    fig, ax = plt.subplots(figsize=(7.6, 4.0))
    ax.scatter(df["cost_min"], df_t80, s=4, c="#dddddd")
    m4 = df["pareto4"] == 1
    sc = ax.scatter(df.loc[m4, "cost_min"], df_t80[m4], c=df.loc[m4, "n"],
                    cmap="viridis", s=16)
    ax.scatter([k38["cost_min"]], [t80], marker="*", s=340, c="#D55E00",
               edgecolors="k", zorder=5)
    ax.set(xlabel="cost proxy [deposition min]", ylabel="lifetime [h]", yscale="log")
    fig.colorbar(sc, ax=ax, label="n dyads")
    st.pyplot(fig, clear_figure=True)
with right:
    st.subheader("Optimize under my constraints")
    cmax = st.number_input("max cost [min]", 100, 9000, 1500, 100)
    wmax = st.number_input("max weight [g/m2]", 0.2, 10.0, 2.0, 0.1)
    dmin = st.number_input("min durability [yr]", 0.0, 100.0, 10.0, 1.0)
    inwin = st.checkbox("respect G1-G2 window", True)
    m = (df.cost_min <= cmax) & (df.weight <= wmax) & (df.dur_yr >= dmin)
    if inwin:
        m &= df.d_in.between(G1, G2)
    if m.any():
        i = df_t80[m].idxmax(); r = df.loc[i]
        st.success(f"**Best: {int(r.n)} x ({int(r.d_in)}/{int(r.d_org)}) nm** - "
                   f"{df_t80[i]:,.0f} h (={df_t80[i]/8760:.1f} yr), "
                   f"{r.cost_min:,.0f} min, {r.weight:.2f} g/m2, {r.dur_yr:.1f} yr dur.")
        def _apply(dorg=int(r.d_org), din=int(r.d_in), nn=int(r.n)):
            st.session_state.update(d_org=dorg, d_in=din, n=nn)
        st.button("Apply this design", on_click=_apply)
    else:
        st.warning("No design satisfies these constraints - relax one of them.")

with st.expander("Why trust this? (validation & limits)"):
    p3 = os.path.join(ROOT, "figures", "Fig3_validation.png")
    if os.path.exists(p3):
        st.image(p3, caption="Zero-refit validation vs Wu 2018.")
    st.markdown(
        "- Solver verified to **0.45%** vs analytic; every mechanism survives an "
        "**ablation study** (SI S7).\n"
        "- Thresholds carry 95% CIs from **N=500 Monte-Carlo through calibration**: "
        "22.5 [21.0-26.0] / 44.0 [37.7-48.0] nm.\n"
        "- **Honest limits**: calibrated at 38 degC on Al2O3/parylene-C; environment "
        "sliders extrapolate by Arrhenius (band above; rankings invariant); durability "
        "is a Basquin proxy; organic planarization unmodeled (100 nm bound).\n"
        f"- Paper, SI, provenance, source: [repository]({REPO}). Built with AI "
        "assistance (Anthropic Claude); physics verified by the author.")
st.caption("(c) 2026 Teo Yang | MIT | v0.9.2 | https://tfe-twin.streamlit.app")
