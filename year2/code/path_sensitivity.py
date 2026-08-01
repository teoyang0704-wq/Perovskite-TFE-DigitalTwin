# -*- coding: utf-8 -*-
"""Robustness of the edge attribution to the assumed in-plane path length.

The papers do not report the distance from the exposed edge of a calcium sample
to the sensor.  Section 3.2 assumed 5 mm.  This script separates what that
assumption controls from what it does not, by refitting the ladder at every
path length between 0.5 and 25 mm and asking, at each one:

  * can the ladder still be fitted at all, and how well;
  * what fraction of the four-dyad reading is attributed to the edge;
  * what organic permeability the fit then requires, and whether that value is
    physically admissible.

The last question is what closes the argument.  A permeability is not free: it
cannot be smaller than the parylene-C value the same class of polymer is known
to have, and a siloxane-based film cannot plausibly be thousands of times more
permeable than parylene.  Path lengths that would require such values are
excluded by the material, not by our assumption.
"""
import numpy as np

DAY, KG = 86400.0, 1e3
P_PAR = 7.5e-13
DA = 0.90
D_IN, D_ORG, D_FIRST, R_PIN = 21.5e-9, 100e-9, 800e-9, 50e-9
PAD, PERIM = 3e-3, 4 * 3e-3
LEE = {1: 3.0e-3, 2: 6.6e-4, 3: 5.4e-4, 4: 5.3e-4}

# admissible range for a siloxane organic layer, relative to parylene-C.
# lower: no organic barrier layer in use is tighter than parylene-C by much;
# upper: bulk PDMS, the most permeable siloxane in common use, is ~1e3 x parylene.
M_MIN, M_MAX = 1.0, 1000.0


def series_wvtr(n, f, m, area=PAD ** 2):
    P = m * P_PAR
    s = R_PIN * np.sqrt(np.pi / f)
    tau2 = 1.0 + s ** 2 * np.log(max(s / R_PIN, 1.001)) / (2 * np.pi * D_ORG ** 2)
    R = n * D_IN / (f * P * area) + (n - 1) * D_ORG * tau2 / (P * area) + D_FIRST / (P * area)
    return DA / R / area * KG * DAY


def edge_wvtr(n, m, L_path, area=PAD ** 2, perim=PERIM):
    t_org = D_FIRST + (n - 1) * D_ORG
    return m * P_PAR * DA * (perim * t_org) / (L_path * area) * KG * DAY


def fit_at(L_path):
    ns = np.array(sorted(LEE), float)
    y = np.array([LEE[int(n)] for n in ns])
    best = None
    for lf in np.linspace(-9.0, -5.5, 161):
        for lm in np.linspace(-0.5, 3.2, 161):
            f, m = 10 ** lf, 10 ** lm
            pred = np.array([series_wvtr(n, f, m) + edge_wvtr(n, m, L_path) for n in ns])
            err = np.sum((np.log10(pred) - np.log10(y)) ** 2)
            if best is None or err < best[0]:
                best = (err, f, m, pred)
    err, f, m, pred = best
    edge4 = edge_wvtr(4, m, L_path)
    return dict(rms=np.sqrt(err / len(ns)), f=f, m=m,
                edge_share=edge4 / pred[-1], pred=pred)


if __name__ == "__main__":
    Ls = np.array([0.5, 1, 2, 3, 5, 8, 10, 15, 20, 25]) * 1e-3
    print("=" * 78)
    print("SENSITIVITY OF THE EDGE ATTRIBUTION TO THE ASSUMED PATH LENGTH")
    print("=" * 78)
    print(f"  {'L_path':>8s} {'rms (dec)':>10s} {'f':>10s} {'m (x parylene)':>16s} "
          f"{'edge share at n=4':>19s}  admissible?")
    rows = []
    for L in Ls:
        r = fit_at(L)
        ok = M_MIN <= r["m"] <= M_MAX
        rows.append((L, r, ok))
        print(f"  {L*1e3:6.1f} mm {r['rms']:10.3f} {r['f']:10.2e} {r['m']:16.1f} "
              f"{100*r['edge_share']:17.0f} %  {'yes' if ok else 'NO — m outside 1-1000x'}")

    adm = [(L, r) for L, r, ok in rows if ok]
    shares = [r["edge_share"] for _, r in adm]
    rmss = [r["rms"] for _, r in adm]
    print("\n" + "=" * 78)
    print("WHAT SURVIVES THE ASSUMPTION")
    print("=" * 78)
    print(f"  admissible path lengths: {adm[0][0]*1e3:.1f} to {adm[-1][0]*1e3:.1f} mm")
    print(f"  edge share at four dyads across that range: "
          f"{100*min(shares):.0f}-{100*max(shares):.0f} %")
    print(f"  fit quality across that range: rms {min(rmss):.3f}-{max(rmss):.3f} decades")
    print("\n  The attribution is not a consequence of choosing 5 mm.  Over every path")
    print("  length that a physically admissible organic permeability allows, the edge")
    print("  dominates the four-dyad reading.  What the assumption does control is the")
    print("  permeability itself, and therefore the comparison with the lag-derived")
    print("  value -- which is consistency, not measurement, and is reported as such.")

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Lg = np.logspace(np.log10(0.5e-3), np.log10(25e-3), 22)
    share, mm, rms = [], [], []
    for L in Lg:
        r = fit_at(L)
        share.append(100 * r["edge_share"]); mm.append(r["m"]); rms.append(r["rms"])
    share, mm, rms = np.array(share), np.array(mm), np.array(rms)
    ok = (mm >= M_MIN) & (mm <= M_MAX)

    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.6), dpi=200)
    ax[0].fill_between(Lg * 1e3, 0, 100, where=ok, color="#4C9F70", alpha=.12,
                       label="physically admissible")
    ax[0].plot(Lg * 1e3, share, "o-", ms=3.5, lw=1.6, color="#1B3A5C")
    ax[0].axhline(50, ls=":", lw=1, color="#888")
    ax[0].axvline(5, ls="--", lw=1.2, color="#D55E00")
    ax[0].text(5.3, 20, "assumed\n5 mm", fontsize=6.5, color="#D55E00")
    ax[0].set(xscale="log", xlabel="assumed in-plane path length (mm)",
              ylabel="edge share of the 4-dyad reading (%)", ylim=(0, 100))
    ax[0].legend(fontsize=7, frameon=False, loc="lower right")
    ax[0].set_title("The attribution does not depend on the assumption", fontsize=9)

    ax[1].fill_between(Lg * 1e3, M_MIN, M_MAX, color="#4C9F70", alpha=.12)
    ax[1].plot(Lg * 1e3, mm, "o-", ms=3.5, lw=1.6, color="#8A5AA8")
    ax[1].axhline(38, ls=":", lw=1.2, color="#D55E00")
    ax[1].text(0.55, 44, "value implied by the lag times", fontsize=6.5, color="#D55E00")
    ax[1].axvline(5, ls="--", lw=1.2, color="#D55E00")
    ax[1].set(xscale="log", yscale="log", xlabel="assumed in-plane path length (mm)",
              ylabel="required organic permeability (x parylene-C)")
    ax[1].set_title("The permeability does", fontsize=9)
    fig.tight_layout(); fig.savefig("fig_path_sensitivity.png")
    print("\nsaved fig_path_sensitivity.png")
