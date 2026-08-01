# -*- coding: utf-8 -*-
"""Year-2, module 2e: what survived as universal.

Module 2d showed that the lower edge of the design window is a property of a
deposition campaign, not of alumina, and that two of five published films have
no window at all.  That reads like the loss of a universal rule.  It is not:
the rule moved up one level of abstraction.

Closure is complete when the nucleation term has fallen to a fixed fraction of
the particulate floor,

    d_closure = d_close * ln( f0 / (eps * f_res) ) = N * d_close,

and N = ln(f0 / (eps f_res)) is a logarithm, so it barely moves: a film a
hundred times cleaner needs only ln(100) = 4.6 more closure lengths.  The
window therefore exists if and only if

    d_crit / d_close  >  N,

a dimensionless criterion with no thickness in it.  Equivalently there is a
critical closure length, d_close < d_crit / N, above which no thickness works:
the film cracks before its defects have closed.

This script computes N and the critical closure length from the posterior, and
checks the criterion film by film.
"""
import numpy as np
from hier_v2 import FILMS, NF, CLOSURE_EPS

FLAT = np.load("hier_v3_posterior.npy")
LF0, D_CRIT, LMU_F, SIG_F, SIG_DC, LMU_DC = 0, 3, 4, 5, 11 + NF, 1


def q(v):
    return np.percentile(v, [5, 50, 95])


if __name__ == "__main__":
    rng = np.random.default_rng(17)
    lf_new = FLAT[:, LMU_F] + FLAT[:, SIG_F] * rng.standard_normal(len(FLAT))

    N_films = np.log(10 ** FLAT[:, LF0] / (CLOSURE_EPS * 10 ** FLAT[:, 6:6 + NF].T)).T
    N_new = np.log(10 ** FLAT[:, LF0] / (CLOSURE_EPS * 10 ** lf_new))

    print("=" * 72)
    print("HOW MANY CLOSURE LENGTHS DOES CLOSURE TAKE?")
    print("=" * 72)
    lo, md, hi = q(N_new)
    print(f"  N = ln(f0 / 0.01 f_res)   {md:5.1f}  [{lo:.1f}, {hi:.1f}]   (dimensionless)")
    print("  per film:")
    for i, f in enumerate(FILMS):
        a = q(N_films[:, i])
        print(f"    {f:12s} {a[1]:5.1f}  [{a[0]:.1f}, {a[2]:.1f}]")
    print("\n  N is a logarithm, so a film 100x cleaner needs only 4.6 more closure")
    print("  lengths.  Almost all of the between-film variation in the window comes")
    print("  from the closure length itself, not from the defect floor.")

    dcrit_over_N = FLAT[:, D_CRIT] / N_new
    lo, md, hi = q(dcrit_over_N)
    print("\n" + "=" * 72)
    print("THE UNIVERSAL CRITERION")
    print("=" * 72)
    print(f"  a window exists iff  d_close < d_crit / N")
    print(f"  critical closure length   {md:5.2f}  [{lo:.2f}, {hi:.2f}] nm")
    print("\n  film          closure length (nm)     P(window exists)")
    for i, f in enumerate(FILMS):
        dc = 10 ** FLAT[:, 12 + NF + i]
        thr = FLAT[:, D_CRIT] / N_films[:, i]
        a = q(dc)
        print(f"  {f:12s}  {a[1]:5.2f} [{a[0]:.2f}, {a[2]:.2f}]        "
              f"{100*np.mean(dc < thr):5.1f} %")
    dc_new = 10 ** (FLAT[:, LMU_DC] + FLAT[:, SIG_DC] * rng.standard_normal(len(FLAT)))
    print(f"  {'unmeasured':12s}  {q(dc_new)[1]:5.2f} [{q(dc_new)[0]:.2f}, {q(dc_new)[2]:.2f}]"
          f"        {100*np.mean(dc_new < dcrit_over_N):5.1f} %")

    print("\n" + "=" * 72)
    print("WHAT IS STILL UNIVERSAL")
    print("=" * 72)
    lo, md, hi = q(FLAT[:, D_CRIT])
    print(f"  upper edge (cracking onset)  {md:5.1f} [{lo:.1f}, {hi:.1f}] nm  "
          f"- mechanical, shared by assumption")
    print(f"  closure takes N = {q(N_new)[1]:.0f} closure lengths, whatever the film")
    print(f"  a window exists iff d_close < {q(dcrit_over_N)[1]:.1f} nm")
    print("  and, from module 1, independent of film: the series bound on ladders,")
    print("  the linear lifetime gain per dyad, the preference for thin organic")
    print("  interlayers, and the device size below which the perimeter dominates.")
    print("\n  What is *not* universal is a number in nanometres for the lower edge.")
