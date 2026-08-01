# -*- coding: utf-8 -*-
"""Year-2, module 2d: what it is worth measuring before you design.

Module 2c let the closure length vary between films, which fixed the fit but
destroyed the predictive design window: for a film nobody has measured, the
lower edge of the window is uncertain by two orders of magnitude and the window
is empty in four draws out of ten.

That is not a failure of the inference; it is the answer to a question Year-1
never asked. Year-1 reported a window as though it were a property of alumina.
It is not: it is a property of a deposition campaign, and it can be pinned down
only by data from that campaign. The useful question is therefore how much data
is enough -- and the posterior already contains the answer, because two of the
five films come with a single-layer thickness series and three do not.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from hier_v2 import FILMS, NF, CLOSURE_EPS

FLAT = np.load("hier_v3_posterior.npy")
LF0, LMU_DC, D_CRIT, LMU_F, SIG_F = 0, 1, 3, 4, 5
SIG_DC = 11 + NF
HAS_SERIES = {"Wu high-q": "5 thicknesses", "Buelow": "3 thicknesses",
              "Wu poor": "none", "Lee": "none", "Carcia": "one thickness, two temperatures"}


def window_for(dc, lf_res, flat):
    """Lower edge of the design window given closure length and defect floor."""
    return dc * np.log(10 ** flat[:, LF0] / (CLOSURE_EPS * 10 ** lf_res))


def q(v):
    return np.percentile(v, [5, 50, 95])


if __name__ == "__main__":
    print("=" * 74)
    print("DESIGN WINDOW, FILM BY FILM")
    print("=" * 74)
    print("  film          data on that film            lower edge (nm)      width")
    rows = {}
    for i, f in enumerate(FILMS):
        dc = 10 ** FLAT[:, 12 + NF + i]
        lo_edge = window_for(dc, FLAT[:, 6 + i], FLAT)
        w = FLAT[:, D_CRIT] - lo_edge
        a, b = q(lo_edge), q(w)
        rows[f] = (lo_edge, w)
        print(f"  {f:12s} {HAS_SERIES[f]:26s} {a[1]:5.1f} [{a[0]:4.1f}, {a[2]:5.1f}]"
              f"   {b[1]:5.1f} [{b[0]:.1f}, {b[2]:.1f}]")

    rng = np.random.default_rng(5)
    lf_new = FLAT[:, LMU_F] + FLAT[:, SIG_F] * rng.standard_normal(len(FLAT))
    dc_new = 10 ** (FLAT[:, LMU_DC] + FLAT[:, SIG_DC] * rng.standard_normal(len(FLAT)))
    lo_new = window_for(dc_new, lf_new, FLAT)
    w_new = FLAT[:, D_CRIT] - lo_new
    a, b = q(lo_new), q(w_new)
    print(f"  {'unmeasured':12s} {'nothing':26s} {a[1]:5.1f} [{a[0]:4.1f}, {a[2]:5.1f}]"
          f"   {b[1]:5.1f} [{b[0]:.1f}, {b[2]:.1f}]")

    print("\n" + "=" * 74)
    print("THE VALUE OF MEASURING YOUR OWN THICKNESS SERIES")
    print("=" * 74)
    ref = rows["Wu high-q"][0]
    f_meas = (q(ref)[2] - q(ref)[0])
    f_un = (a[2] - a[0])
    print(f"  90 % interval on the window's lower edge")
    print(f"    with a five-point thickness series : {f_meas:6.1f} nm")
    print(f"    with nothing                       : {f_un:6.1f} nm")
    print(f"    ratio                              : {f_un/f_meas:6.0f}x narrower")
    print(f"  probability of an empty window")
    print(f"    measured film (Wu high quality)    : {100*np.mean(rows['Wu high-q'][1] <= 0):.1f} %")
    print(f"    unmeasured film                    : {100*np.mean(w_new <= 0):.1f} %")
    print("\n  Five single-layer depositions and five calcium or coulometric readings")
    print("  are enough to turn an unusable prior into a usable window.  That is the")
    print("  experiment the twin should be asked to justify, and it is cheap.")

    fig, ax = plt.subplots(figsize=(6.4, 4.0), dpi=200)
    bins = np.linspace(0, 80, 120)
    ax.hist(np.clip(lo_new, 0, 80), bins=bins, density=True, alpha=.45,
            color="#8A5AA8", label="film nobody has measured")
    ax.hist(np.clip(rows["Wu high-q"][0], 0, 80), bins=bins, density=True, alpha=.6,
            color="#1B3A5C", label="film with a 5-point thickness series")
    ax.axvspan(*np.percentile(FLAT[:, D_CRIT], [5, 95]), color="#D55E00", alpha=.15)
    ax.axvline(np.median(FLAT[:, D_CRIT]), color="#D55E00", lw=1.6)
    ax.text(np.median(FLAT[:, D_CRIT]) + 1, ax.get_ylim()[1] * .85,
            "cracking onset\n(shared, mechanical)", fontsize=7, color="#D55E00")
    ax.set(xlabel="lower edge of the design window (nm)", ylabel="density", xlim=(0, 80))
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    ax.set_title("The window is a property of your deposition, not of alumina", fontsize=9.5)
    fig.tight_layout(); fig.savefig("fig_window_value.png")
    print("\nsaved fig_window_value.png")
