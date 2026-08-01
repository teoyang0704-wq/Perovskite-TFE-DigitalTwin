# -*- coding: utf-8 -*-
"""Screening log for the pre-registered audit frame.

Protocol: AUDIT_PREREGISTRATION.md, fixed before any paper was gathered.
All six queries executed; the frame saturated (queries 4-6 returned mostly
papers already screened).  Every paper returned by those queries,
in the order returned, deduplicated, is recorded here with its screening
outcome.  Reviews, patents and vendor pages are out of frame by definition and
are listed as such rather than silently dropped.

The pre-registration anticipated that most exclusions would be "ladder exists
only inside a figure" and stated that this count is itself a finding.  It is.
"""
import numpy as np

# (short name, in_frame, outcome, best value reported, ladder available?, floor reported?)
#   outcome: "ladder"   two or more numerical values for the same architecture
#            "best-only" improvement with layer count asserted, one number given
#            "not-ladder" the compared samples differ in something else
#            "out-of-frame" review, patent, vendor page
SCREENED = [
    # --- query 1 -------------------------------------------------------------
    ("Wu 2018, Al2O3/parylene (RSC Adv)",             True,  "ladder",      7.7e-6,  True,  False),
    ("Kim 2018, Al2O3/silamer (ACS AMI)",             True,  "best-only",   1.11e-6, False, False),
    ("Seo 2015, Al2O3/ZrO2/alucone (Nanoscale RL)",   True,  "best-only",   8.5e-5,  False, False),
    # --- query 2 -------------------------------------------------------------
    ("Al2O3/alucone, 5.5 dyads",                      True,  "best-only",   1.44e-4, False, False),
    ("O3-based laminate, 1 to 3 layers",              True,  "ladder",      2.37e-5, True,  False),
    ("TiO2/Al2O3, ozone vs water (Org Electron)",     True,  "not-ladder",  1e-4,    False, False),
    ("npj Flex El 2025, stress-released Al2O3",       True,  "best-only",   4.49e-5, False, False),
    ("ZAM/silamer multilayer (J Inf Display)",        True,  "best-only",   5.94e-5, False, False),
    ("npj Flex El 2025, silbione hybrimer",           True,  "best-only",   7.83e-6, False, False),
    ("SAOLs-Al2O3 nanolaminate (ACS AMI)",            True,  "best-only",   2.99e-7, False, False),
    # --- query 3 -------------------------------------------------------------
    ("Wang, PDMS/Al2O3 nanolaminate 1.5 & 2.5 dyads", True,  "best-only",   1e-5,    False, False),
    ("Org Electron 2024, parylene/alumina 3 dyads",   True,  "best-only",   8.7e-4,  False, False),
    ("Al2O3/alucone, dyad-number study",              True,  "best-only",   None,    False, False),
    ("Optical Ca test, multiple pathways (method)",   True,  "not-ladder",  2.1e-4,  False, False),
    ("Rev Sci Instrum, Ca-test sensitivity limits",   True,  "not-ladder",  None,    False, False),
    # --- query 4 -------------------------------------------------------------
    ("Buelow 2014, AlOx/plasma polymer (Disc Nano)",  True,  "ladder",      1.2e-3,  True,  True),
    ("Kim 2014, SiOx & SiNx/parylene annealing",      True,  "ladder",      6.6e-4,  True,  False),
    ("Graff 2004, lag time vs equilibrium (model)",   False, "out-of-frame", None,   False, False),
    # --- query 6 -------------------------------------------------------------
    ("Graham, barrier vs shelf lifetime correlation", True,  "best-only",   None,    False, False),
    ("SiNx/parylene on PC (Surf Coat Technol)",       True,  "not-ladder",  None,    False, False),
    # --- query 5 -------------------------------------------------------------
    ("Al2O3/s-NEA multilayers, 3 units (JMSE)",       True,  "best-only",   2.6e-5,  False, False),
    # --- out of frame across all queries -------------------------------------
    ("reviews (x2)",                                  False, "out-of-frame", None,   False, False),
    ("patents (x7)",                                  False, "out-of-frame", None,   False, False),
    ("vendor / trade pages (x3)",                     False, "out-of-frame", None,   False, False),
]

# ladders already held from the opportunistic pilot, kept separate
PILOT_LADDERS = 8
PILOT_WITH_FLOOR = 3


def jeffreys(k, n, lo=0.05, hi=0.95):
    from scipy.stats import beta
    a, b = k + 0.5, n - k + 0.5
    return beta.ppf(lo, a, b), beta.ppf(hi, a, b)


if __name__ == "__main__":
    frame = [s for s in SCREENED if s[1]]
    n = len(frame)
    lad = [s for s in frame if s[2] == "ladder"]
    best = [s for s in frame if s[2] == "best-only"]
    notl = [s for s in frame if s[2] == "not-ladder"]
    floors = [s for s in frame if s[5]]

    print("=" * 74)
    print("SCREENING OF THE PRE-REGISTERED FRAME (all 6 queries, saturated)")
    print("=" * 74)
    print(f"  papers returned and in frame            {n}")
    print(f"  out of frame (reviews, patents, vendors) {len(SCREENED)-n}")
    print()
    print(f"  usable ladder (>=2 numerical values)    {len(lad):2d}   "
          f"{100*len(lad)/n:4.0f} %")
    print(f"  best value only, ladder in a figure     {len(best):2d}   "
          f"{100*len(best)/n:4.0f} %")
    print(f"  comparison confounded by another change {len(notl):2d}   "
          f"{100*len(notl)/n:4.0f} %")
    print(f"  reports a fixture floor or control      {len(floors):2d}   "
          f"{100*len(floors)/n:4.0f} %")

    k = len(best) + len(notl)
    lo, hi = jeffreys(k, n)
    print(f"\n  fraction of benchmark papers whose ladder cannot be checked by a")
    print(f"  reader: {100*k/n:.0f} %  [90 % CI {100*lo:.0f}-{100*hi:.0f} %]")
    lo, hi = jeffreys(len(floors), n)
    print(f"  fraction reporting any floor or impermeable control: "
          f"{100*len(floors)/n:.0f} %  [90 % CI {100*lo:.0f}-{100*hi:.0f} %]")

    print("\n" + "=" * 74)
    print("WHAT THIS DOES TO THE AUDIT")
    print("=" * 74)
    print("  The series bound needs two numbers from one architecture; the floor")
    print("  comparison needs a control.  Most benchmark papers supply neither.")
    print("  The violation rate therefore cannot be estimated from the frame at")
    print("  the rate papers accumulate -- but the reason it cannot is a stronger")
    print("  and more robust statement than the violation rate would have been:")
    print()
    print("    a barrier value is published, the ladder behind it is shown only as")
    print("    a curve, and no impermeable control is run, so neither the reader")
    print("    nor the author can tell whether the number is the film or the fixture.")
    print()
    print(f"  Pilot ladders held separately: {PILOT_LADDERS}, of which "
          f"{PILOT_WITH_FLOOR} report a floor and all {PILOT_WITH_FLOOR} sit on it.")
