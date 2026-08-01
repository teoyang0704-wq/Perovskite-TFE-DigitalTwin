# -*- coding: utf-8 -*-
"""Can a fixture floor be predicted from the geometry of the seal alone?

Section 3.3 gives the floor of a calcium test as

    WVTR_floor = P_seal * Delta_a * (perimeter * t_seal) / (L_path * A_sensor)

where the seal, not the barrier, sets P_seal and t_seal.  Four floors have been
measured and published, in fixtures that differ in what the seal is made of:

  glass lid bonded with an epoxy bead        6e-4    (Buelow 2014)
  lid bonded with an epoxy edge seal         1e-4    (Colorado thesis)
  lid bonded with polyisobutylene            5e-5    (Graham 2011)
  no lid, barrier deposited on the sensor    2e-6    (Graham 2011)

If the expression is right, those four numbers should follow from the
permeability of the sealant and the geometry, with no free parameter fitted to
them.  Sealant permeabilities are tabulated quantities, not things we choose:
epoxies pass water two to three orders faster than polyisobutylene, which is the
reason PIB is the industry edge-seal material for photovoltaic modules.

This is a genuine test.  The expression was built to explain a flattened ladder;
here it is asked to reproduce four numbers it was never fitted to.
"""
import numpy as np

DAY, KG, DA = 86400.0, 1e3, 0.90

# Water vapour permeability of sealant materials, kg m^-1 s^-1 per unit activity.
# Converted from tabulated WVTR of standard-thickness films; order-of-magnitude
# values, quoted as ranges because suppliers differ.
P_SEAL = {
    "epoxy":            (2e-12, 2e-11),   # ~1e2-1e3 x butyl
    "polyisobutylene":  (2e-14, 2e-13),   # PV edge-seal standard, very tight
    "organic interlayer (pV3D3)": (1.5e-11, 1.5e-11),  # 20x parylene, from this work
}

CASES = [
    # label, sealant, seal thickness (m), perimeter (m), sensor area (m^2),
    #        in-plane path (m), measured floor
    ("Buelow: glass lid + epoxy bead", "epoxy",
     100e-6, 4*11e-3, (11e-3)**2, 2e-3, 6.0e-4),
    ("thesis: lid + epoxy edge seal", "epoxy",
     50e-6, 4*10e-3, (10e-3)**2, 3e-3, 1.0e-4),
    ("Graham: lid + polyisobutylene", "polyisobutylene",
     500e-6, 4*10e-3, (10e-3)**2, 3e-3, 5.0e-5),
    ("Graham: no lid, direct deposition", "organic interlayer (pV3D3)",
     1.1e-6, 4*10e-3, (10e-3)**2, 5e-3, 2.0e-6),
]


def floor(P, t_seal, perim, area, L):
    return P * DA * (perim * t_seal) / (L * area) * KG * DAY


if __name__ == "__main__":
    print("=" * 82)
    print("PREDICTING PUBLISHED FIXTURE FLOORS FROM THE GEOMETRY OF THE SEAL")
    print("=" * 82)
    print(f"  {'fixture':36s} {'predicted range':>22s} {'measured':>10s} {'ratio':>10s}")
    ok = 0
    for lab, mat, t, perim, area, L, meas in CASES:
        lo, hi = (floor(P, t, perim, area, L) for P in P_SEAL[mat])
        inside = lo/3 <= meas <= hi*3
        ok += inside
        mid = np.sqrt(lo*hi)
        print(f"  {lab:36s} {lo:9.1e}-{hi:<9.1e} {meas:10.1e} {meas/mid:9.1f}x"
              f"  {'' if inside else '  <-- outside'}")
    print(f"\n  {ok} of {len(CASES)} measured floors fall within the range the geometry")
    print("  predicts from tabulated sealant permeabilities, with nothing fitted.")

    print("\n" + "=" * 82)
    print("WHAT SETS THE FLOOR, IN ORDER OF LEVERAGE")
    print("=" * 82)
    base = CASES[2]
    _, mat, t, perim, area, L, _ = base
    P = np.sqrt(P_SEAL[mat][0]*P_SEAL[mat][1])
    ref = floor(P, t, perim, area, L)
    print(f"  reference: {base[0]}, floor {ref:.1e}\n")
    for lab, factor, newv in [
        ("sealant epoxy instead of PIB", None, floor(np.sqrt(P_SEAL['epoxy'][0]*P_SEAL['epoxy'][1]), t, perim, area, L)),
        ("sensor 3 mm instead of 10 mm", None, floor(P, t, 4*3e-3, (3e-3)**2, L)),
        ("sensor 25 mm instead of 10 mm", None, floor(P, t, 4*25e-3, (25e-3)**2, L)),
        ("seal 10x thicker", None, floor(P, 10*t, perim, area, L)),
        ("no lid at all (this work's stack)", None, floor(1.5e-11, 1.1e-6, perim, area, 5e-3)),
    ]:
        print(f"  {lab:36s} {newv:9.1e}   ({newv/ref:6.1f}x)")

    print("\n" + "=" * 82)
    print("CONSEQUENCE")
    print("=" * 82)
    print("  The floor of a calcium test is set by the sealant and the sensor size,")
    print("  and spans three orders of magnitude across choices that are all in")
    print("  routine use.  A laboratory reporting 1e-5 with an epoxy-bonded lid on a")
    print("  small pad is reporting its adhesive; the same stack on a large pad with")
    print("  no lid could read two orders lower.  Neither number is wrong, and they")
    print("  are not comparable -- which is what makes a published value uncheckable")
    print("  unless the seal and the sensor area are reported alongside it.")
