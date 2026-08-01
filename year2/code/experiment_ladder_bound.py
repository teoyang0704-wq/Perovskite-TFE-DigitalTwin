# -*- coding: utf-8 -*-
"""Year-2 experiment 1: how flat can a dyad ladder be?

Every mechanism in the Year-1 twin, and every defect arrangement in the 3-D
network, places the n inorganic layers *in series*.  Adding a dyad can
therefore never help by less than the series limit.  The flattest possible
ladder is the one in which lateral detours contribute nothing at all, i.e. all
transport runs through particle-punched columnar channels, giving

        G(n) = G_1 / n        ->        G(n2)/G(n1) = n1 / n2.

This script measures the ladder produced by the 3-D network as the columnar
fraction phi is swept from 0 to 1, confirms that it saturates at the 1/n bound,
and compares both with published ladders.  A measured ladder flatter than 1/n
cannot be explained by any series arrangement of defects and requires a
parallel, dyad-independent ingress channel.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from percolation3d import Stack, place_defects, effective_conductance

NX, L_OVER_S, N_REAL = 160, 8.0, 4
PHIS = (0.0, 0.03, 0.1, 0.3, 1.0)
N_MAX = 4

# published ladders, WVTR in g m^-2 day^-1 at 38 C / 90 %RH
LEE = {1: 3.0e-3, 2: 6.6e-4, 3: 5.4e-4, 4: 5.3e-4}          # AEM 8, 1701928, Fig. 5a
WU_HQ = {1: 1.70e-4, 2: 3.60e-5, 3: 7.70e-6}                 # RSC Adv. 8, 5721, Fig. 3


def model_ladder(base: Stack, phi: float):
    L = L_OVER_S * base.spacing
    out = []
    for n in range(2, N_MAX + 1):
        st = Stack(n_inorg=n, d_in=base.d_in, d_org=base.d_org, f=base.f,
                   r_pin=base.r_pin, P_org=base.P_org)
        g = []
        for k in range(N_REAL):
            rng = np.random.default_rng(2000 + 91 * k + n)
            lay, col = place_defects(st, L, rng, phi, "poisson")
            g.append(effective_conductance(st, L, NX, lay, col))
        out.append(float(np.mean(g)))
    return np.array(out)


if __name__ == "__main__":
    lee = Stack(n_inorg=4, d_in=21.5e-9, d_org=100e-9, f=5.6e-7,
                r_pin=50e-9, P_org=7.5e-13)
    ns = np.arange(2, N_MAX + 1)

    print("model ladders, normalised to n = 2")
    curves = {}
    for phi in PHIS:
        g = model_ladder(lee, phi)
        curves[phi] = g / g[0]
        print(f"  phi={phi:<5.2f} " + "  ".join(f"n={n}:{v:.3f}"
                                                for n, v in zip(ns, curves[phi])))
    bound = 2.0 / ns
    print(f"  series bound (1/n) " + "  ".join(f"n={n}:{v:.3f}"
                                               for n, v in zip(ns, bound)))

    lee_rel = np.array([LEE[n] / LEE[2] for n in ns])
    print("\nmeasured, normalised to n = 2")
    print("  Lee et al.        " + "  ".join(f"n={n}:{v:.3f}"
                                             for n, v in zip(ns, lee_rel)))
    wu_rel = np.array([WU_HQ[n] / WU_HQ[2] for n in (2, 3)])
    print("  Wu et al. (high q) " + "  ".join(f"n={n}:{v:.3f}"
                                              for n, v in zip((2, 3), wu_rel)))
    viol = lee_rel > bound
    print("\nviolation of the series bound (measured flatter than 1/n):")
    for n, m, b, v in zip(ns, lee_rel, bound, viol):
        print(f"  n={n}: measured {m:.3f} vs bound {b:.3f}  ->  "
              f"{'VIOLATED' if v else 'consistent'}")

    # residual dyad-independent channel implied by the violation
    # G(n) = A/n + C  fitted to the Lee ladder (A series part, C parallel leak)
    G_meas = np.array([LEE[n] for n in ns])
    M = np.c_[1.0 / ns, np.ones_like(ns, dtype=float)]
    A, C = np.linalg.lstsq(M, G_meas, rcond=None)[0]
    print(f"\nfit  G(n) = A/n + C  ->  A = {A:.2e},  C = {C:.2e} g m^-2 day^-1")
    print(f"  dyad-independent channel accounts for {100*C/G_meas[-1]:.0f}% "
          f"of the 4-dyad reading")

    fig, ax = plt.subplots(figsize=(5.4, 4.0), dpi=200)
    for phi in PHIS:
        ax.plot(ns, curves[phi], "o-", ms=4, lw=1.2, label=f"3-D model, $\\phi$={phi}")
    ax.plot(ns, bound, "k--", lw=1.6, label="series bound $1/n$")
    ax.plot(ns, lee_rel, "s-", color="#D55E00", ms=7, lw=2, label="Lee et al. (measured)")
    ax.set(xlabel="number of dyads $n$", ylabel="WVTR relative to $n=2$",
           yscale="log", xticks=list(ns))
    ax.legend(fontsize=7, frameon=False)
    ax.set_title("Measured ladders are flatter than any series model allows", fontsize=9)
    fig.tight_layout()
    fig.savefig("fig_ladder_bound.png")
    print("\nsaved fig_ladder_bound.png")
