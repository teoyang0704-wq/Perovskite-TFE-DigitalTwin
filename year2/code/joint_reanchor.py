# -*- coding: utf-8 -*-
"""Year-2 experiment 3: joint re-anchoring, and design rules for finite devices.

Experiment 2 fitted the Lee ladder with two independent scales -- one for the
series path, one for the edge -- and noted the degeneracy honestly: the organic
permeability P_org multiplies both, while the defect fraction f multiplies only
the series term.  Forcing a single P_org therefore *breaks* the degeneracy
rather than worsening it, because the edge term carries no f:

    series  ~  f * P_org          edge  ~  P_org

Fitting (m = P_org / P_parylene, f) jointly to the four-point ladder gives both.
Three independent observables then have to agree on m: the lag times (Year-1,
SI S4.5), the flattening of the ladder, and the absolute level of the ladder.

The second half of the script asks what the edge channel means for design.
Year-1 derived its rules for an unbounded area.  Real devices have a perimeter,
so the benefit of adding a dyad saturates once the edge term dominates -- and
that happens at a device size, not at a stack thickness.
"""
import numpy as np

DAY, KG_TO_G = 86400.0, 1e3
P_PAR = 7.5e-13            # parylene-C permeability D*S [kg m^-1 s^-1 per activity]
DA = 0.90                  # driving activity, 38 C / 90 %RH

# Lee et al. stack geometry (AEM 8, 1701928)
D_IN, D_ORG, D_FIRST, R_PIN = 21.5e-9, 100e-9, 800e-9, 50e-9
LEE = {1: 3.0e-3, 2: 6.6e-4, 3: 5.4e-4, 4: 5.3e-4}
PAD, L_PATH = 3e-3, 5e-3   # calcium pad side, in-plane path length (assumed)

# Year-1 anchors for comparison
F_WU = 3.72e-8             # calibrated defect fraction, Wu high-quality film
F_MODEB = 5.6e-7           # Lee defect fraction inferred without an edge term
M_LAG = 38.0               # organic permeability multiplier from the lag times


def series_wvtr(n, f, m, area):
    """Through-the-stack WVTR of an n-dyad Lee stack [g m^-2 day^-1]."""
    P = m * P_PAR
    s = R_PIN * np.sqrt(np.pi / f)                       # mean pinhole spacing
    tau2 = 1.0 + s ** 2 * np.log(s / R_PIN) / (2 * np.pi * D_ORG ** 2)
    R = n * D_IN / (f * P * area)                        # pinhole-limited inorganics
    R += (n - 1) * D_ORG * tau2 / (P * area)             # interlayer detours
    R += D_FIRST / (P * area)                            # 800 nm passivation, no detour
    return DA / R / area * KG_TO_G * DAY


def edge_wvtr(n, m, area, perimeter, L_path=L_PATH):
    """Lateral ingress from the sample edge [g m^-2 day^-1]."""
    t_org = D_FIRST + (n - 1) * D_ORG
    J = m * P_PAR * DA * (perimeter * t_org) / L_path
    return J / area * KG_TO_G * DAY


def joint_fit():
    area, perim = PAD ** 2, 4 * PAD
    ns = np.array(sorted(LEE), dtype=float)
    y = np.array([LEE[int(n)] for n in ns])
    best = None
    for lf in np.linspace(-9.0, -5.5, 141):              # f from 1e-9 to 3e-6
        for lm in np.linspace(-0.5, 2.6, 141):           # m from 0.3 to 400
            f, m = 10 ** lf, 10 ** lm
            pred = np.array([series_wvtr(n, f, m, area)
                             + edge_wvtr(n, m, area, perim) for n in ns])
            err = np.sum((np.log10(pred) - np.log10(y)) ** 2)
            if best is None or err < best[0]:
                best = (err, f, m, pred)
    err, f, m, pred = best
    print("=" * 72)
    print("JOINT RE-ANCHORING of the Lee ladder (one organic permeability)")
    print("=" * 72)
    print(f"  fitted defect fraction f      = {f:.2e}")
    print(f"  fitted permeability  m        = {m:.1f} x parylene-C")
    print(f"  rms residual                  = {np.sqrt(err/len(y)):.3f} decades")
    print("\n   n   measured    model      series      edge    edge share")
    for n, ym, yp in zip(ns, y, pred):
        s_ = series_wvtr(n, f, m, area)
        e_ = edge_wvtr(n, m, area, perim)
        print(f"  {int(n)}   {ym:.2e}   {yp:.2e}   {s_:.2e}  {e_:.2e}  {100*e_/yp:5.1f} %")

    print("\n  cross-checks")
    print(f"    m from lag times (independent)      : {M_LAG:.0f}"
          f"   -> agreement within {max(m,M_LAG)/min(m,M_LAG):.1f}x")
    print(f"    f without an edge term (mode B)     : {F_MODEB:.1e}"
          f"   -> revised downward by {F_MODEB/f:.0f}x")
    print(f"    f of the Year-1 calibration film    : {F_WU:.1e}"
          f"   -> Lee film is {f/F_WU:.1f}x that density")
    s = R_PIN * np.sqrt(np.pi / f)
    print(f"    implied pinhole spacing             : {s*1e6:.0f} um "
          f"({f/(np.pi*R_PIN**2)/1e6:.1f} defects per mm2)")

    big = 1.0                                            # 1 m x 1 m coupon
    print("\n  what the same stack would read on a large sample (no edge limit)")
    for n in (2, 3, 4):
        print(f"    n={n}: series only = {series_wvtr(n, f, m, big*big):.2e} "
              f"g m^-2 day^-1   (reported {LEE[n]:.1e})")
    return f, m


def dyads_worth_adding(f, m, sizes=(3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0), n_max=8,
                       gain_threshold=0.10):
    """For a square device of side a, how many dyads still help?

    A dyad is 'worth adding' while it lowers the total WVTR by more than
    `gain_threshold` (10 % by default).  Once the edge term dominates, further
    dyads change nothing -- and they make the edge slightly worse.
    """
    print("\n" + "=" * 72)
    print("DESIGN RULE: how many dyads a finite device can actually use")
    print("=" * 72)
    print("  device side   useful dyads   floor WVTR      limited by")
    for a in sizes:
        area, perim = a * a, 4 * a
        tot = [series_wvtr(n, f, m, area) + edge_wvtr(n, m, area, perim)
               for n in range(1, n_max + 1)]
        useful = 1
        for k in range(1, n_max):
            if (tot[k - 1] - tot[k]) / tot[k - 1] > gain_threshold:
                useful = k + 1
            else:
                break
        n_star = useful
        e = edge_wvtr(n_star, m, area, perim)
        s_ = series_wvtr(n_star, f, m, area)
        lim = "edge" if e > s_ else "stack"
        print(f"  {a*1e3:7.0f} mm   {useful:^12d}   {tot[n_star-1]:.2e}    {lim}")
    print("\n  Below roughly a centimetre the stack stops being the limit: the\n"
          "  perimeter is. Adding dyads there buys nothing, and thicker organic\n"
          "  interlayers make it worse, since they widen the lateral cross-section.")


def figure(f, m, path="fig_finite_device.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ns = np.arange(1, 9)
    fig, ax = plt.subplots(figsize=(5.4, 4.0), dpi=200)
    for a, c in zip((3e-3, 1e-2, 1e-1, 1.0),
                    ("#D55E00", "#8A5AA8", "#4C9F70", "#1B3A5C")):
        area, perim = a * a, 4 * a
        tot = [series_wvtr(n, f, m, area) + edge_wvtr(n, m, area, perim) for n in ns]
        ax.plot(ns, tot, "o-", ms=4, lw=1.6, color=c,
                label=f"device {a*1e3:.0f} mm")
    ax.plot(ns, [series_wvtr(n, f, m, 1.0) for n in ns], "k--", lw=1.4,
            label="unbounded area (Year-1)")
    ax.set(xlabel="number of dyads $n$", ylabel="WVTR  (g m$^{-2}$ day$^{-1}$)",
           yscale="log", xticks=list(ns))
    ax.legend(fontsize=7.5, frameon=False)
    ax.set_title("Below a centimetre, the perimeter sets the barrier", fontsize=9)
    fig.tight_layout(); fig.savefig(path)
    print(f"\nsaved {path}")


if __name__ == "__main__":
    f, m = joint_fit()
    dyads_worth_adding(f, m)
    figure(f, m)
