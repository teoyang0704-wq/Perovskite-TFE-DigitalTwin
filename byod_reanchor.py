# -*- coding: utf-8 -*-
"""BYOD re-anchoring: calibrate the twin on *your* film, and audit your own
measurements.

The published calibration describes one laboratory's deposition campaign. Year-2
established that the design window is a property of a campaign rather than of
alumina, so applying our numbers to your film is the wrong thing to do. This
page does the right thing instead: you supply a thickness series, it re-fits the
two quantities that vary between campaigns, and it returns the window for your
film with credible intervals.

It also runs three diagnostics on your measurements that no other tool performs:
the series bound, the fixture floor, and the steady-state check.
"""
import io
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st

RG, TREF = 8.314, 298.15
D_PAR, S_PAR, EA_PAR = 5.0e-13, 1.5, 40e3
D_LAT, S_IN, EA_LAT = 1e-21, 0.10, 60e3
R_SUB_DEFAULT = 2.88e7
CLOSURE_EPS = 0.01
REPO = "https://github.com/teoyang0704-wq/Perovskite-TFE-DigitalTwin"

arr = lambda x, Ea, T: x * np.exp(-Ea / RG * (1 / T - 1 / TREF))


# --------------------------------------------------------------- physics
def f_shape(d_nm, f0, d_close, f_res, C, d_crit):
    crack = np.maximum(0.0, (d_nm - d_crit) / d_crit) ** 2
    return f0 * np.exp(-d_nm / d_close) + f_res * (1.0 + C * crack)


def wvtr_single(d_nm, p, T_C, RH, R_sub):
    """Steady WVTR of a bare barrier on a substrate, g m^-2 day^-1."""
    T, da = T_C + 273.15, RH / 100.0
    P_org = arr(D_PAR, EA_PAR, T) * S_PAR
    P_lat = arr(D_LAT, EA_LAT, T) * S_IN
    f = f_shape(d_nm, p["f0"], p["d_close"], p["f_res"], p["C"], p["d_crit"])
    P_in = P_lat + f * P_org
    R = R_sub + d_nm * 1e-9 / P_in
    return da / R * 1e3 * 86400.0


def fit_film(d_nm, wvtr, T_C, RH, R_sub, n_boot=400, seed=0):
    """Re-fit the two campaign-specific quantities on a user's thickness series.

    Shared physics (nucleation amplitude, cracking coefficient and onset) is held
    at the published values; the closure length and particulate floor are free,
    because those are what Year-2 found to vary between campaigns.
    """
    d_nm, wvtr = np.asarray(d_nm, float), np.asarray(wvtr, float)
    base = dict(f0=0.372, C=40.0, d_crit=44.0)
    grid_dc = np.logspace(np.log10(0.15), np.log10(12.0), 140)
    grid_fr = np.logspace(-10, -4.0, 140)
    best, surf = None, np.zeros((len(grid_dc), len(grid_fr)))
    for i, dc in enumerate(grid_dc):
        for j, fr in enumerate(grid_fr):
            p = dict(base, d_close=dc, f_res=fr)
            pred = np.array([wvtr_single(d, p, T_C, RH, R_sub) for d in d_nm])
            err = float(np.sum((np.log10(pred) - np.log10(wvtr)) ** 2))
            surf[i, j] = err
            if best is None or err < best[0]:
                best = (err, dc, fr)
    err, dc, fr = best
    rms = math.sqrt(err / len(d_nm))

    rng = np.random.default_rng(seed)
    boot = []
    for _ in range(n_boot):
        k = rng.integers(0, len(d_nm), len(d_nm))
        bb = None
        for dcx in grid_dc[::3]:
            for frx in grid_fr[::3]:
                p = dict(base, d_close=dcx, f_res=frx)
                pred = np.array([wvtr_single(d, p, T_C, RH, R_sub) for d in d_nm[k]])
                e = float(np.sum((np.log10(pred) - np.log10(wvtr[k])) ** 2))
                if bb is None or e < bb[0]:
                    bb = (e, dcx, frx)
        boot.append((bb[1], bb[2]))
    boot = np.array(boot)
    return dict(base, d_close=dc, f_res=fr), rms, boot


def window(p, eps=CLOSURE_EPS):
    """Lower edge (closure complete) and upper edge (cracking onset), nm."""
    lo = p["d_close"] * math.log(p["f0"] / (eps * p["f_res"]))
    return lo, p["d_crit"]


# --------------------------------------------------------------- diagnostics
def series_bound(pairs):
    """pairs: list of (n, wvtr). Returns rows and a verdict."""
    rows, worst = [], None
    ns = sorted(pairs)
    for i, a in enumerate(ns):
        for b in ns[i + 1:]:
            ratio, bound = pairs[b] / pairs[a], a / b
            rows.append(dict(n1=a, n2=b, measured=ratio, bound=bound,
                             verdict="VIOLATED" if ratio > bound else "ok",
                             approach=ratio / bound))
            if worst is None or ratio / bound > worst:
                worst = ratio / bound
    return rows, worst


def fixture_floor(sensor_mm, seal, t_seal_um, L_path_mm, T_C, RH):
    """Floor imposed by the path joining sample to instrument.

    Two regimes. With a lid, the path runs through the adhesive bead, and the
    cross-section is the bondline height. With the barrier deposited directly on
    the sensor there is no adhesive: the lateral path is the stack's own organic
    interlayers, whose total thickness is of order a micrometre, and the seal
    thickness entered by the user is ignored in favour of that.
    """
    # Effective permeabilities calibrated against the four published fixture
    # floors (6e-4, 1e-4, 5e-5 with a lid; 2e-6 with the barrier on the sensor).
    # A geometry-only calculation from tabulated sealant permeabilities
    # overpredicts those floors by 6-50x, because cells contain a getter and
    # because bondline dimensions are almost never reported; these effective
    # values absorb that, and are order-of-magnitude estimates only.
    P = {"none (barrier on sensor)": 2.9e-13,
         "polyisobutylene": 1.6e-15,
         "epoxy": 1.6e-13,
         "unknown / not reported": 1.6e-13}[seal]
    a = sensor_mm * 1e-3
    t = 1.1e-6 if seal.startswith("none") else t_seal_um * 1e-6
    J = P * (RH / 100.0) * (4 * a * t) / (L_path_mm * 1e-3)
    return J / (a * a) * 1e3 * 86400.0


def steady_state_fraction(t_hours, lag_hours):
    tau = max(lag_hours, 1e-9) * 6 / math.pi ** 2
    x = max(t_hours / tau, 1e-9)
    s = 1.0 + 2 * sum((-1) ** m * math.exp(-m * m * x) for m in range(1, 120))
    return max(min(s, 1.0), 0.0)


# --------------------------------------------------------------- UI
def render():
    st.header("Calibrate on your own film")
    st.caption(
        "The published calibration describes one laboratory's campaign. Year-2 of "
        "this project showed the design window is a property of a deposition "
        "campaign, not of alumina — so the right move is to re-fit on your data, "
        "not to borrow ours. Five single-layer points are enough: they narrow the "
        "window's lower edge by roughly eighty-fold."
    )

    tab1, tab2 = st.tabs(["1 · Re-anchor on your film", "2 · Audit your measurements"])

    # ---------------- tab 1
    with tab1:
        c1, c2 = st.columns([2, 1])
        with c2:
            T_C = st.number_input("Test temperature (°C)", 20.0, 90.0, 38.0, 1.0)
            RH = st.number_input("Relative humidity (%)", 10.0, 100.0, 90.0, 5.0)
            R_sub = st.number_input("Substrate resistance (m² s kg⁻¹)",
                                    value=float(R_SUB_DEFAULT), format="%.3e",
                                    help="Bare-substrate WVTR converted to a resistance; "
                                         "leave at the default for PET unless you measured it.")
        with c1:
            st.markdown("**Your single-layer thickness series**")
            st.caption("At least three points; at least two below 20 nm, or the closure "
                       "length is not identifiable (Year-2, module 2f).")
            default = pd.DataFrame({"thickness_nm": [15.0, 20.0, 30.0, 50.0, 60.0],
                                    "wvtr_g_m2_day": [6.7e-3, 0.7e-3, 0.8e-3, 1.3e-3, 4.7e-3]})
            df = st.data_editor(default, num_rows="dynamic", use_container_width=True,
                                key="series")

        up = st.file_uploader("…or upload a CSV with columns thickness_nm, wvtr_g_m2_day",
                              type=["csv"])
        if up is not None:
            df = pd.read_csv(up)

        if st.button("Re-anchor on this film", type="primary"):
            d = df["thickness_nm"].astype(float).values
            w = df["wvtr_g_m2_day"].astype(float).values
            if len(d) < 3:
                st.error("At least three thicknesses are needed.")
                return
            thin = int((d < 20).sum())
            with st.spinner("Fitting your campaign…"):
                p, rms, boot = fit_film(d, w, T_C, RH, R_sub)
            lo, hi = window(p)
            los = [window(dict(p, d_close=b[0], f_res=b[1]))[0] for b in boot]
            lo_ci = np.percentile(los, [5, 95])

            st.success(f"Fitted your film: closure length {p['d_close']:.2f} nm, "
                       f"particulate floor {p['f_res']:.2e} — rms {rms:.3f} decades")
            m1, m2, m3 = st.columns(3)
            m1.metric("Window lower edge", f"{lo:.1f} nm",
                      f"90 % CI {lo_ci[0]:.1f}–{lo_ci[1]:.1f}")
            m2.metric("Window upper edge", f"{hi:.1f} nm", "cracking onset (shared)")
            m3.metric("Usable width", f"{hi - lo:.1f} nm",
                      "negative means no thickness works")

            if thin < 2:
                st.warning(
                    f"Only {thin} of your points lie below 20 nm. Closure is measurable "
                    "only where the barrier is still poor; with points this thick the "
                    "closure length is set by the prior, not by your data, and the "
                    "lower edge above should be treated as indicative."
                )
            if p["d_close"] > 2.0:
                st.error(
                    f"Closure length {p['d_close']:.2f} nm exceeds the critical value of "
                    "about 2 nm. For this campaign the film cracks before its defects "
                    "close, and **no inorganic thickness gives a usable barrier**. "
                    "The lever is deposition conditions, not geometry."
                )

            fig, ax = plt.subplots(figsize=(6.2, 3.6), dpi=160)
            grid = np.linspace(max(2, d.min() * 0.6), d.max() * 1.25, 250)
            ax.plot(grid, [wvtr_single(x, p, T_C, RH, R_sub) for x in grid],
                    color="#1B3A5C", lw=1.8, label="your re-anchored film")
            ax.plot(d, w, "o", color="#D55E00", ms=7, label="your measurements")
            if hi > lo:
                ax.axvspan(lo, hi, color="#4C9F70", alpha=.15, label="your design window")
            ax.set(xlabel="inorganic thickness (nm)",
                   ylabel="WVTR (g m$^{-2}$ day$^{-1}$)", yscale="log")
            ax.legend(fontsize=7, frameon=False)
            fig.tight_layout()
            st.pyplot(fig)

            out = pd.DataFrame([{"closure_length_nm": p["d_close"],
                                 "particulate_floor": p["f_res"],
                                 "window_lower_nm": lo, "window_lower_lo": lo_ci[0],
                                 "window_lower_hi": lo_ci[1], "window_upper_nm": hi,
                                 "fit_rms_decades": rms, "T_C": T_C, "RH": RH}])
            st.download_button("Download your calibration (CSV)",
                               out.to_csv(index=False).encode(),
                               "my_film_calibration.csv", "text/csv")

    # ---------------- tab 2
    with tab2:
        st.markdown("**Three checks on your own measurements.** None of these needs "
                    "our model; they need only numbers you already have.")

        st.divider()
        st.markdown("### A · Series bound")
        st.caption("Layers in series cannot improve transmission more slowly than 1/n. "
                   "A violation means something conducts in parallel with your stack, "
                   "or your layers are not equivalent.")
        lad = st.data_editor(pd.DataFrame({"n_dyads": [1, 2, 3], "wvtr": [3.0e-3, 6.6e-4, 5.4e-4]}),
                             num_rows="dynamic", key="ladder", use_container_width=True)
        if st.button("Run the bound"):
            pairs = {float(r.n_dyads): float(r.wvtr) for r in lad.itertuples()}
            rows, worst = series_bound(pairs)
            st.dataframe(pd.DataFrame(rows), use_container_width=True)
            if worst and worst > 1:
                st.error(f"Bound violated (worst case {worst:.2f}× the bound). The most "
                         "common cause is lateral ingress through the seal joining your "
                         "sample to the instrument — check B below.")
            elif worst and worst > 0.9:
                st.warning(f"Within {100*(1-worst):.0f} % of the bound. Not a violation, "
                           "but close enough that measurement scatter could cross it.")
            else:
                st.success("Bound satisfied with margin.")

        st.divider()
        st.markdown("### B · What your fixture can resolve")
        c1, c2, c3 = st.columns(3)
        sensor = c1.number_input("Sensor / pad side (mm)", 1.0, 200.0, 10.0, 1.0)
        seal = c2.selectbox("How is the sample sealed?",
                            ["epoxy", "polyisobutylene", "none (barrier on sensor)",
                             "unknown / not reported"])
        t_seal = c3.number_input("Seal bondline height (µm)", 1.0, 2000.0, 25.0, 5.0,
                                 help="Ignored when the barrier is deposited directly on the "
                                      "sensor: the lateral path is then the stack's own organic "
                                      "layers, taken as ~1.1 µm.")
        c4, c5 = st.columns(2)
        L_path = c4.number_input("Bondline width / in-plane path (mm)", 0.1, 50.0, 2.0, 0.5)
        target = c5.number_input("Value you want to report (g m⁻² day⁻¹)",
                                 value=1e-5, format="%.2e")
        floor = fixture_floor(sensor, seal, t_seal, L_path, 38.0, 90.0)
        st.metric("Estimated floor of your fixture", f"{floor:.2e} g m⁻² day⁻¹")
        if target < floor:
            st.error(f"Your target is **{floor/target:.0f}× below** your fixture's floor. "
                     "A value there would be your seal, not your film. Enlarge the pad, "
                     "change the sealant, or deposit the barrier directly on the sensor.")
        else:
            st.success(f"Your target sits {target/floor:.1f}× above the floor — resolvable, "
                       "though a margin under about 2× is uncomfortable.")
        st.caption("Order-of-magnitude only: bondline dimensions and getter capacity are "
                   "rarely reported, and our attempt to predict four published floors from "
                   "geometry alone was 6–50× off. Run an impermeable control to measure it.")

        st.divider()
        st.markdown("### C · Did your run reach steady state?")
        c1, c2 = st.columns(2)
        t_run = c1.number_input("Run duration (h)", 1.0, 5000.0, 72.0, 6.0)
        lag = c2.number_input("Lag time of your stack (h)", 0.1, 5000.0, 60.0, 5.0,
                              help="From the intercept of the transient, or estimated.")
        frac = steady_state_fraction(t_run, lag)
        st.metric("Fraction of the eventual flux you are seeing", f"{100*frac:.0f} %")
        if frac < 0.9:
            st.error(f"You are reading {100*frac:.0f} % of the steady-state flux. Your "
                     f"barrier looks about {1/max(frac,1e-6):.1f}× better than it is, and "
                     "because lag grows with dyad number the error grows along a ladder — "
                     "which bends it and inflates any exponent fitted to it.")
        else:
            st.success("Close enough to steady state to report.")

    st.divider()
    st.caption(f"Method and code: {REPO} · the diagnostics implement results from the "
               "Year-2 and Year-3 analyses in that repository.")


if __name__ == "__main__":
    st.set_page_config(page_title="BYOD re-anchoring", layout="wide")
    render()
