# -*- coding: utf-8 -*-
"""Year-2 experiment 2: the dyad-independent channel is the sample edge.

Experiment 1 showed that published dyad ladders can be flatter than the series
bound G(n2)/G(n1) <= n1/n2, which no arrangement of defects can produce.  Two
explanations survived: a parallel ingress channel that does not scale with dyad
number, or barrier layers that degrade with deposition order.

This script tests the first.  Water can reach the sensor by diffusing laterally
through the organic interlayers from the cut edge of the sample, bypassing every
inorganic layer.  That channel is set by the *test geometry*, not by the stack:

    J_edge = P_org * da * (perimeter * t_org) / L_path                [kg/s]
    WVTR_edge = J_edge / A_sensor                                     [g/m2/day]

It therefore predicts a sharp, checkable signature:

  * Ca tests with a small sensor pad (large perimeter-to-area ratio, exposed
    edge) should show the violation;
  * coulometric tests such as MOCON, where the film is clamped and the measured
    area is edge-sealed, should not;
  * the channel grows with total organic thickness, hence weakly *with* dyad
    number, making the ladder flatter still.

The database contains one Ca-test ladder (Lee) and two MOCON ladders (Wu, high
and poor quality), which is exactly the comparison needed.
"""
import numpy as np
from percolation3d import Stack

DAY = 86400.0
KG_TO_G = 1e3

# ---------------------------------------------------------------- ladders
LADDERS = {
    "Wu high-quality (MOCON)":  {1: 1.70e-4, 2: 3.60e-5, 3: 7.70e-6},
    "Wu poor-quality (MOCON)":  {1: 1.6,     2: 1.6e-1,  3: 1.3e-1},
    "Lee (Ca test, 3 mm pad)":  {1: 3.0e-3,  2: 6.6e-4,  3: 5.4e-4, 4: 5.3e-4},
}
QUALITY = {
    "Wu high-quality (MOCON)": "explicit values in text",
    "Wu poor-quality (MOCON)": "digitised from a log-scale figure (+-30 %)",
    "Lee (Ca test, 3 mm pad)": "explicit values in figure legend",
}


def survey_bound(ladders=LADDERS):
    """Check every pair (n1 < n2) against G(n2)/G(n1) <= n1/n2."""
    print("=" * 74)
    print("LADDER SURVEY: violations of the series bound  G(n2)/G(n1) <= n1/n2")
    print("=" * 74)
    for name, d in ladders.items():
        ns = sorted(d)
        flags = []
        for i, n1 in enumerate(ns):
            for n2 in ns[i + 1:]:
                ratio, bound = d[n2] / d[n1], n1 / n2
                flags.append((n1, n2, ratio, bound, ratio > bound))
        nv = sum(f[-1] for f in flags)
        print(f"\n{name}   [{QUALITY[name]}]")
        for n1, n2, r, b, v in flags:
            mark = "VIOLATED" if v else "ok"
            print(f"   {n1}->{n2}:  measured {r:6.3f}   bound {b:6.3f}   {mark}")
        print(f"   => {nv} of {len(flags)} pairs violate the bound")


# ---------------------------------------------------------------- edge model
def t_org_total(n, d_first=800e-9, d_rest=100e-9):
    """Total organic thickness of an n-dyad stack (first layer is thicker)."""
    return d_first + max(n - 1, 0) * d_rest


def wvtr_edge(P_org, t_org, perimeter, area, L_path, da=0.90):
    """Apparent WVTR contributed by lateral ingress from the sample edge."""
    J = P_org * da * (perimeter * t_org) / L_path          # kg/s
    return J / area * KG_TO_G * DAY


def fit_lee(P_par=7.5e-13, pad=3e-3, L_path=5e-3):
    """Split the Lee ladder into a series part and an edge part.

    G(n) = alpha * G_1D(n) + WVTR_edge(n);  free parameters: alpha and the
    organic permeability multiplier m = P_org / P_parylene.
    """
    d = LADDERS["Lee (Ca test, 3 mm pad)"]
    ns = np.array(sorted(d), dtype=float)
    y = np.array([d[int(n)] for n in ns])

    area, perim = pad ** 2, 4 * pad
    base = Stack(n_inorg=1, d_in=21.5e-9, d_org=100e-9, f=5.6e-7,
                 r_pin=50e-9, P_org=P_par)

    def series(n):
        st = Stack(n_inorg=int(n), d_in=base.d_in, d_org=base.d_org,
                   f=base.f, r_pin=base.r_pin, P_org=base.P_org)
        return st.conductance_1d(area) * 0.90 / area * KG_TO_G * DAY

    S = np.array([series(n) for n in ns])
    E1 = np.array([wvtr_edge(P_par, t_org_total(int(n)), perim, area, L_path)
                   for n in ns])

    # least squares in log space over (alpha, m), both positive
    best = None
    for la in np.linspace(-3, 1, 121):
        for lm in np.linspace(-1, 2.6, 121):
            pred = 10 ** la * S + 10 ** lm * E1
            err = np.sum((np.log10(pred) - np.log10(y)) ** 2)
            if best is None or err < best[0]:
                best = (err, 10 ** la, 10 ** lm, pred)
    err, alpha, m, pred = best
    print("\n" + "=" * 74)
    print(f"EDGE-CHANNEL FIT to the Lee ladder  (pad {pad*1e3:.0f} mm, "
          f"lateral path {L_path*1e3:.0f} mm)")
    print("=" * 74)
    print(f"  series scale alpha = {alpha:.3g}   organic permeability "
          f"m = P_org/P_parylene = {m:.1f}")
    print(f"  rms error in log10 = {np.sqrt(err/len(y)):.3f} decades")
    print("   n   measured    predicted    series      edge     edge share")
    for n, ym, yp, s_, e_ in zip(ns, y, pred, alpha * S, m * E1):
        print(f"  {int(n)}   {ym:.2e}   {yp:.2e}   {s_:.2e}  {e_:.2e}   "
              f"{100*e_/yp:5.1f} %")
    print(f"\n  independent estimate from the lag analysis (Year-1, SI S4.5):"
          f" P_org/P_parylene ~ 38")
    print(f"  edge-channel estimate here: {m:.0f}  ->  agreement within a "
          f"factor {max(m,38)/min(m,38):.1f}")
    return alpha, m


def resolution_limit(m, P_par=7.5e-13, L_path=5e-3, da=0.90):
    """Smallest WVTR a Ca test of a given pad size can resolve before its own
    edge channel dominates."""
    print("\n" + "=" * 74)
    print("MEASUREMENT DESIGN RULE: edge-limited resolution of a Ca test")
    print("=" * 74)
    print("  pad side   organic 1.1 um        organic 0.5 um")
    for pad in (1e-3, 3e-3, 10e-3, 25e-3, 50e-3):
        vals = [wvtr_edge(m * P_par, t, 4 * pad, pad ** 2, L_path, da)
                for t in (1.1e-6, 0.5e-6)]
        print(f"  {pad*1e3:5.0f} mm   {vals[0]:.2e} g/m2/day   "
              f"{vals[1]:.2e} g/m2/day")
    print("\n  A stack whose true WVTR lies below the value for its pad size is\n"
          "  measuring its own edge, not its barrier.")


if __name__ == "__main__":
    survey_bound()
    alpha, m = fit_lee()
    resolution_limit(m)


def figure(alpha, m, P_par=7.5e-13, pad=3e-3, L_path=5e-3, path="fig_edge_channel.png"):
    """Decomposition of the Lee ladder into series and edge contributions."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = LADDERS["Lee (Ca test, 3 mm pad)"]
    ns = np.array(sorted(d), dtype=float)
    y = np.array([d[int(n)] for n in ns])
    area, perim = pad ** 2, 4 * pad

    def series(n):
        st = Stack(n_inorg=int(n), d_in=21.5e-9, d_org=100e-9, f=5.6e-7,
                   r_pin=50e-9, P_org=P_par)
        return st.conductance_1d(area) * 0.90 / area * KG_TO_G * DAY

    S = alpha * np.array([series(n) for n in ns])
    E = m * np.array([wvtr_edge(P_par, t_org_total(int(n)), perim, area, L_path)
                      for n in ns])
    fig, ax = plt.subplots(figsize=(5.4, 4.0), dpi=200)
    ax.plot(ns, y, "s", color="#D55E00", ms=9, label="Lee et al. (measured)")
    ax.plot(ns, S + E, "-", color="#1B3A5C", lw=2, label="model: series + edge")
    ax.plot(ns, S, "--", color="#4C9F70", lw=1.6, label="series (through the stack)")
    ax.plot(ns, E, ":", color="#8A5AA8", lw=1.8, label="edge (lateral, through organics)")
    ax.set(xlabel="number of dyads $n$", ylabel="WVTR  (g m$^{-2}$ day$^{-1}$)",
           yscale="log", xticks=[1, 2, 3, 4])
    ax.legend(fontsize=7.5, frameon=False)
    ax.set_title("A 3 mm calcium pad measures its own edge beyond two dyads",
                 fontsize=9)
    fig.tight_layout(); fig.savefig(path)
    print(f"\nsaved {path}")
