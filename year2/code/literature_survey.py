# -*- coding: utf-8 -*-
"""Year-2 experiment 4: the edge floor has already been measured, twice.

Module 1 inferred a dyad-independent channel from the shape of a ladder.
A literature search then turned up two groups who measured that channel
directly, reported the number, and did not draw the conclusion:

  * Buelow et al., Nanoscale Res. Lett. 9, 223 (2014), open access.  A glass
    lid -- an essentially impermeable barrier -- was run through the same
    calcium test as their multilayers and read about 6e-4 g m^-2 day^-1 at
    60 C / 90 %RH.  Their best 3.5-dyad multilayer read 1.2e-3, which the
    authors describe as "only by a factor of 2 higher" than the glass lid.
    Read the other way: half of that measurement is the test fixture.

  * Colorado thesis on calcium testing of ALD barriers: an ALD film and a
    glass lid both returned about 1e-4 g m^-2 day^-1, "assigned to H2O
    permeation through the epoxy edge seal", with the explicit warning that
    the calcium test cannot evaluate barriers at or below 5e-5 and that the
    epoxy edge seal is the largest limiting factor.

Both of those fixtures seal a lid with epoxy, so their edge channel runs
through the adhesive.  Lee et al. deposit the encapsulation directly onto the
calcium with no lid and no epoxy, so nothing in their geometry can leak except
the encapsulation's own organic interlayers -- which is the channel module 1
proposed and quantified.  The phenomenon is therefore established; what is new
here is that it also occurs in lid-free, directly deposited stacks, that it can
be read off the shape of a dyad ladder, and that the series bound tells you
when it is happening without a control sample.
"""
import numpy as np

# ladders: name -> {n: WVTR}, plus measured floor where the authors provided one
LADDERS = {
    "Wu high-quality (MOCON, sealed)": (
        {1: 1.70e-4, 2: 3.60e-5, 3: 7.70e-6}, None, "38 C / 100 %RH"),
    "Wu poor-quality (MOCON, sealed)": (
        {1: 1.6, 2: 1.6e-1, 3: 1.3e-1}, None, "38 C / 100 %RH"),
    "Lee (Ca, 3 mm pad, no lid)": (
        {1: 3.0e-3, 2: 6.6e-4, 3: 5.4e-4, 4: 5.3e-4}, None, "38 C / 90 %RH"),
    "Buelow (Ca, glass lid + epoxy)": (
        {1.5: 3.6e-3, 3.5: 1.2e-3}, 6.0e-4, "60 C / 90 %RH"),
    # Groner 2005: 26 nm Al2O3 on Kapton, one side vs both sides = two layers
    # in series, measured in a tritiated-water cell with a clamped Viton o-ring
    "Groner (tritium cell, o-ring)": (
        {1: 1.0e-3, 2: 7.0e-4}, None, "100 %RH, ambient"),
    # added 2026-08-01: five-point ladder, RSC Adv. 9 (2019) 20884
    "Al2O3/alucone (Ca, 25 C/60 %RH)": (
        {1.5: 1.74e-2, 2.5: 2.47e-3, 3.5: 6.41e-4, 4.5: 2.23e-4, 5.5: 1.44e-4},
        None, "25 C/60 %RH"),
    # two-point ladder; the second value sits at 99 % of its bound
    "SiNx/Al2O3 (Ca)": (
        {1.5: 2.6e-4, 2.5: 1.55e-4}, None, "Ca test"),
}
# independently measured fixture floors (impermeable control through the same test)
FLOORS = {
    "Buelow, glass lid + epoxy bead": (
        6.0e-4, "glass lid run through the same calcium fixture; 11x11 mm aperture"),
    "Bertrand thesis, lid + epoxy edge seal": (
        1.0e-4, "an ALD film and a glass lid both read this; authors assign it to the epoxy seal"),
    "Graham, lid + polyisobutylene sealant": (
        5.0e-5, "separate test of the edge sealant alone"),
    "Graham, no lid, barrier on the sensor": (
        2.0e-6, "same laboratory and method, no adhesive seal"),
    "Lee (inferred, not measured)": (
        3.9e-4, "no lid: inferred here from the ladder shape; lateral path is the stack's own organics"),
}


def survey():
    print("=" * 76)
    print("EXTENDED SURVEY: series bound  G(n2)/G(n1) <= n1/n2")
    print("=" * 76)
    for name, (d, floor, cond) in LADDERS.items():
        ns = sorted(d)
        print(f"\n{name}   [{cond}]")
        nv = 0
        for i, n1 in enumerate(ns):
            for n2 in ns[i + 1:]:
                r, b = d[n2] / d[n1], n1 / n2
                nv += r > b
                print(f"   {n1}->{n2}:  measured {r:6.3f}   bound {b:6.3f}   "
                      f"{'VIOLATED' if r > b else 'ok'}")
        if floor:
            top = d[max(ns)]
            print(f"   measured fixture floor {floor:.1e} = {100*floor/top:.0f} % "
                  f"of the {max(ns)}-dyad reading")
            # what the ladder looks like once the measured floor is removed
            corr = {n: d[n] - floor for n in ns}
            if min(corr.values()) > 0:
                r = corr[max(ns)] / corr[min(ns)]
                b = min(ns) / max(ns)
                print(f"   after subtracting the floor: {r:.3f} vs bound {b:.3f}"
                      f"   -> {'still violated' if r > b else 'consistent'}")
        print(f"   => {nv} violation(s)")


def floors():
    print("\n" + "=" * 76)
    print("DIRECTLY MEASURED OR INFERRED FIXTURE FLOORS")
    print("=" * 76)
    for k, (v, note) in FLOORS.items():
        print(f"  {k:28s} {v:.1e} g m^-2 day^-1   ({note})")
    print("\n  Three independent instances spanning 1e-4 to 6e-4 g m^-2 day^-1.")
    print("  Two were measured with an impermeable control; the third was inferred")
    print("  from the ladder shape alone and lands in the same range.")


def implication():
    print("\n" + "=" * 76)
    print("WHAT THIS MEANS FOR PUBLISHED VALUES")
    print("=" * 76)
    print("  A calcium test reports the sum of the stack and its fixture.  Any")
    print("  published multilayer WVTR within a factor of a few of 1e-4 that was")
    print("  measured on a small pad, with or without an epoxy-sealed lid, is a")
    print("  candidate for being fixture-limited.  Two diagnostics need no extra")
    print("  hardware:")
    print("    1. run the series bound on the dyad ladder; a violation is proof")
    print("       that something parallel to the stack is contributing;")
    print("    2. compare the best value with an impermeable control measured in")
    print("       the same fixture -- the control is the floor, not zero.")
    print("\n  And one that does: invert the deposition order.  A fixture channel")
    print("  is unchanged by it; degrading layer quality reverses.")


if __name__ == "__main__":
    survey()
    floors()
    implication()
