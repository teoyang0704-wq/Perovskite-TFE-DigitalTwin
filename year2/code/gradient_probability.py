# -*- coding: utf-8 -*-
"""How likely is the observed ladder if the layers merely differ at random?

Section 3.2 notes that the ladder can be reproduced by four independent layer
resistances -- a quality gradient -- and that this fit is unfalsifiable because
it has one parameter per data point.  That is an argument from parsimony.  This
script replaces it with a probability.

The question is not whether *some* set of four resistances reproduces the
ladder; one always exists.  It is whether the particular pattern the ladder
requires is likely to arise from layers that differ only by ordinary process
scatter.  Inverting the measured ladder gives the resistance each added layer
contributed:

    layer 1 (+ substrate)  333       layer 3   337
    layer 2               1182       layer 4    35

so the fourth layer contributes 34 times less than the second, and the sequence
is not merely scattered but monotonically collapsing after layer 2.

Two things are needed to turn that into a probability: a distribution for
layer-to-layer scatter, and a definition of the pattern.  For the first we do
not invent a number.  The published thickness series of a single laboratory,
measured on nominally identical films, gives the scatter of a real process
directly; we use the residual scatter of that series about its fitted trend as
the layer-to-layer sigma, and repeat the whole calculation across a wide range
of sigma so that the conclusion does not depend on it.

For the second we require two properties that the observed ladder has and that
a parallel channel produces automatically:

  P1  the ratio of the largest to smallest layer contribution is at least 34;
  P2  the contributions after the second layer are monotonically decreasing.

P2 matters because random scatter produces large ratios often enough, but
rarely in a monotone sequence -- and a fixture channel produces exactly a
monotone collapse, because each added layer contributes less against a fixed
parallel path.
"""
import numpy as np

# contributions implied by the Lee ladder, in units of the first
OBS = np.array([333.3, 1181.8, 336.7, 34.9])
RATIO_OBS = OBS.max() / OBS.min()          # 34
N_TRIAL = 200_000


def process_sigma_from_wu():
    """Layer-to-layer scatter of a real process, from a published series.

    Wu et al. report five nominally identical depositions at different
    thicknesses; the scatter of those values about the fitted thickness trend
    bounds how much two layers from one process differ.
    """
    d = np.array([15.0, 20.0, 30.0, 50.0, 60.0])
    w = np.array([6.7e-3, 0.7e-3, 0.8e-3, 1.3e-3, 4.7e-3])
    # fit the non-monotonic trend the paper itself describes (closure + cracking)
    A = np.c_[np.ones_like(d), np.exp(-d / 1.2), np.maximum(0, (d - 44) / 44) ** 2]
    coef, *_ = np.linalg.lstsq(A, np.log10(w), rcond=None)
    resid = np.log10(w) - A @ coef
    return float(np.std(resid, ddof=1))


def trial_probability(sigma_dex, n_trial=N_TRIAL, seed=0, monotone=True):
    """P(a random four-layer stack shows the observed pattern)."""
    rng = np.random.default_rng(seed)
    R = 10 ** (sigma_dex * rng.standard_normal((n_trial, 4)))
    ratio = R.max(axis=1) / R.min(axis=1)
    hit = ratio >= RATIO_OBS
    if monotone:
        mono = (R[:, 1] > R[:, 2]) & (R[:, 2] > R[:, 3])
        hit &= mono
    return hit.mean()


if __name__ == "__main__":
    s_wu = process_sigma_from_wu()
    print("=" * 76)
    print("IS THE LADDER PLAUSIBLE AS ORDINARY PROCESS SCATTER?")
    print("=" * 76)
    print(f"  layer contributions implied by the measured ladder: "
          f"{np.round(OBS,0)}")
    print(f"  largest / smallest = {RATIO_OBS:.0f}")
    print(f"\n  layer-to-layer scatter of a real published process "
          f"(residual of a five-point series): sigma = {s_wu:.2f} decades")

    print("\n  probability that four layers differing only by scatter reproduce")
    print("  a spread of >=34x, and do so monotonically after layer 2:\n")
    print(f"  {'sigma (dec)':>12s} {'P(spread >= 34x)':>18s} {'P(spread and monotone)':>24s}")
    for sig in (s_wu, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5):
        p_any = trial_probability(sig, monotone=False, seed=1)
        p_mon = trial_probability(sig, monotone=True, seed=2)
        tag = "  <- measured process" if abs(sig - s_wu) < 1e-9 else ""
        print(f"  {sig:12.2f} {p_any:18.4f} {p_mon:24.5f}{tag}")

    p = trial_probability(s_wu, monotone=True, seed=3)
    print(f"\n  At the scatter of a real process, the pattern arises by chance with")
    print(f"  probability {p:.4f}" + (f" (about 1 in {1/p:.0f})" if p > 0 else " (none in %d trials)" % N_TRIAL))
    print("\n  The same pattern is the generic signature of a parallel channel: with a")
    print("  fixed path beside the stack, each added layer necessarily contributes")
    print("  less than the last, so P2 holds by construction rather than by chance.")

    print("\n" + "=" * 76)
    print("WHAT THIS DOES AND DOES NOT SETTLE")
    print("=" * 76)
    print("  It does not exclude a systematic gradient -- a process that genuinely")
    print("  degrades with each layer would produce the same pattern deliberately,")
    print("  and no probability argument can rule that out.  What it excludes is the")
    print("  weaker and more common claim that ordinary run-to-run variation could")
    print("  have produced the ladder.  Distinguishing a systematic gradient from a")
    print("  fixture channel still requires inverting the deposition order.")
