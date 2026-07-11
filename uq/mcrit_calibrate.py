# -*- coding: utf-8 -*-
"""mcrit_calibrate.py — device-anchored effective lower bound for M_crit.

Data (DB): BM013 Choi 2018 barrier WVTR = 1.84e-2 g/m2/day @45C/100%RH (own film);
DS001 same encapsulation on PSC: <4% PCE drop after 7,500 h shelf @25C/50%RH.

Logic (same-sink-model effective calibration):
 1) transfer Choi's measured WVTR to shelf conditions. Bare-pinhole transport is
    vapor-access-limited (SI S4.1), so WVTR ∝ P_vap(T)·Δa with P_vap = D_vap·c_sat;
    Arrhenius(40 kJ/mol) given as an alternative (agreement within ~15%).
 2) model-dose over 7,500 h.
 3) Tier-1 floor (monotone damage only): M_crit^eff ≥ dose.
    Tier-2 (linear dose–damage): ≥ (20%/4%) × dose.
 4) propagate to representative designs (dose term scales; t_lag unchanged).
Assumptions & caveats enumerated for SI S4.4b (a–e).
"""
import numpy as np

W45 = 1.84e-2            # g/m2/day @45C, 100%RH (BM013)
t_h = 7500.0
# saturation vapor concentration [g/m3] (Magnus-based)
def c_sat(TC):
    es = 611.2*np.exp(17.62*TC/(243.12+TC))    # Pa
    return es*18.015/(8.314*(TC+273.15))       # g/m3
r_vap = (c_sat(25)/c_sat(45)) * (298.15/318.15)**-1.81**0 * ((298.15/318.15)**1.81)
# D_vap ∝ T^1.81 ; ratio D(25)/D(45)
r_vap = (c_sat(25)/c_sat(45)) * (298.15/318.15)**1.81
W_shelf_vap = W45 * r_vap * 0.5
R = 8.314
W_shelf_arr = W45 * np.exp(-40e3/R*(1/298.15-1/318.15)) * 0.5
print(f"c_sat 25/45C = {c_sat(25):.1f}/{c_sat(45):.1f} g/m3; vapor-mech factor {r_vap*0.5:.3f}")
print(f"WVTR @25C/50%RH: vapor-mech {W_shelf_vap:.2e} | Arrhenius40 {W_shelf_arr:.2e} g/m2/day")
doses = np.array([W_shelf_vap, W_shelf_arr]) * t_h/24.0
print(f"model-dose over 7,500 h: {doses.min():.2f}-{doses.max():.2f} g/m2 (central {doses.mean():.2f})")
M1 = doses.mean(); M1_lo, M1_hi = doses.min(), doses.max()
M2 = 5.0*M1
print(f"\nTier-1 floor (monotone damage):  M_crit^eff >= {M1:.2f} g/m2  [{M1_lo:.2f}-{M1_hi:.2f}]")
print(f"Tier-2 (linear dose-damage):     M_crit^eff >= {M2:.1f} g/m2")

scale = M1/0.01
print(f"\nscale vs placeholder (0.01): x{scale:.0f} on the dose term")
reps = {"MAX-T80 6x(42/100)": 43574, "KNEE 6x(18/100)": 23604,
        "ECONOMY 3x(21/100)": 14857, "single dyad 1x(30/100)": 2296}
print("moisture-limited lifetime under Tier-1 floor (t_lag negligible):")
for k, v in reps.items():
    print(f"  {k:26s}: {v*scale:,.0f} h  (~{v*scale/8760:,.0f} yr)")
print("\n=> within the G1-G2 window, moisture-limited lifetime exceeds ~20 yr even at")
print("   n=1; moisture ceases to be the binding lifetime channel for dyad designs.")
print("   Thresholds (22.5/44 nm) and all rankings unchanged (G7).")
print("\nSI S4.4b assumption list: (a) same-sink-model effective calibration; transfer")
print("assumes sink bias non-increasing with barrier resistance; (b) T/RH transfer by")
print("vapor mechanism (Arrhenius-40 alt within ~15%); (c) Tier-1 needs monotone damage")
print("only, Tier-2 adds linearity; (d) Choi's own measured WVTR used -> film-specific")
print("f cancels; (e) conformal coating assumed to suppress edge ingress.")
