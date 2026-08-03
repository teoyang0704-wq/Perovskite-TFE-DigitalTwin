# -*- coding: utf-8 -*-
"""Year-2, module 3: an audit of published multilayer barrier ladders.

Two questions, both answerable from numbers that papers already print.

  Q1  How often does a published ladder violate the series bound
      G(n2)/G(n1) <= n1/n2, which no arrangement of defects can escape?
      This needs only the ladder, not the fixture geometry.

  Q2  How often does a published value sit at or below the floor of its own
      test fixture -- the value that fixture returns for an impermeable
      sample?  This needs a floor, which is available whenever the authors
      measured a control or stated their resolution.

Every entry below is a number stated in the text or abstract of the source, or
read from a legend; nothing is digitised from a curve except where marked.  The
registry is deliberately small and verifiable rather than large and soft.

Author: Teo Yang.
"""
import numpy as np

# ---------------------------------------------------------------- registry
# key: (ladder {n: WVTR}, fixture floor or None, seal type, conditions, note)
LADDERS = {
    "Wu 2018 high quality": (
        {1: 1.70e-4, 2: 3.60e-5, 3: 7.70e-6}, None, "MOCON, clamped",
        "38 C/100 %RH", "values in text"),
    "Wu 2018 poor quality": (
        {1: 1.6, 2: 1.6e-1, 3: 1.3e-1}, None, "MOCON, clamped",
        "38 C/100 %RH", "digitised from a log figure, +-30 %"),
    "Lee 2018": (
        {1: 3.0e-3, 2: 6.6e-4, 3: 5.4e-4, 4: 5.3e-4}, None, "Ca, no lid",
        "38 C/90 %RH", "values in figure legend"),
    "Buelow 2014": (
        {1.5: 3.6e-3, 3.5: 1.2e-3}, 6.0e-4, "Ca, glass lid + epoxy",
        "60 C/90 %RH", "floor measured with a glass lid"),
    "Groner 2005 (1 vs 2 sides)": (
        {1: 1.0e-3, 2: 7.0e-4}, None, "tritium cell, Viton o-ring",
        "100 %RH, ambient", "double-sided coating = two layers in series"),
    "Graham 2011 SiNx/parylene": (
        {1: 4.0e-3, 4: 4.0e-5}, 5.0e-5, "Ca, lid + polyisobutylene sealant",
        "20 C/50 %RH", "1-dyad inferred from 'two orders between 1 and 4 dyads'"),
    "Graham 2011 ZnO/Al2O3 nanolaminate": (
        {1: 7.0e-4, 10: 5.0e-5}, 5.0e-5, "Ca, lid + polyisobutylene sealant",
        "20 C/50 %RH", "10-dyad value is censored at the floor (<5e-5)"),
    "Kim 2014 SiOx/parylene": (
        {1: 3.1e-2, 3: 6.6e-4}, None, "Ca", "20 C/50 %RH", "values in text"),
    "Al2O3/alucone 2019 (five points)": (
        {1.5: 1.74e-2, 2.5: 2.47e-3, 3.5: 6.41e-4, 4.5: 2.23e-4, 5.5: 1.44e-4},
        None, "Ca", "25 C/60 %RH", "values in text"),
    "SiNx/Al2O3 hygroscopic interlayer": (
        {1.5: 2.6e-4, 2.5: 1.55e-4}, None, "Ca", "n/s", "values in text"),
}

# fixture floors reported by the authors themselves, from an impermeable
# control or from a separate test of the seal
FLOORS = {
    "Buelow 2014": (6.0e-4, "glass lid run through the same calcium fixture"),
    "Colorado thesis": (1.0e-4, "ALD film and glass lid both read this; "
                                "attributed to the epoxy edge seal"),
    "Graham 2011, lid + PIB sealant": (5.0e-5, "separate test of the edge sealant"),
    "Graham 2011, direct deposition": (2.0e-6, "same laboratory, no lid, no sealant"),
    "Lee 2018 (inferred)": (3.9e-4, "not measured; inferred here from the ladder shape"),
}


def bound_test(ladder):
    """Return (n_pairs, n_violations, worst_case)."""
    ns = sorted(ladder)
    npair = nviol = 0
    worst = None
    for i, n1 in enumerate(ns):
        for n2 in ns[i + 1:]:
            r, b = ladder[n2] / ladder[n1], n1 / n2
            npair += 1
            if r > b:
                nviol += 1
                if worst is None or (r - b) > worst[2]:
                    worst = (n1, n2, r - b, r, b)
    return npair, nviol, worst


if __name__ == "__main__":
    print("=" * 78)
    print("Q1  SERIES BOUND  G(n2)/G(n1) <= n1/n2")
    print("=" * 78)
    tot_p = tot_v = 0
    viol_papers = []
    for name, (lad, floor, seal, cond, note) in LADDERS.items():
        p, v, w = bound_test(lad)
        tot_p += p; tot_v += v
        tag = "VIOLATED" if v else "ok"
        if v:
            viol_papers.append(name)
        extra = ""
        if w:
            extra = f"  worst {w[0]}->{w[1]}: {w[3]:.3f} vs {w[4]:.3f}"
        print(f"  {name:36s} {seal:32s} {tag:9s}{extra}")
    print(f"\n  {len(viol_papers)} of {len(LADDERS)} ladders violate the bound "
          f"({100*len(viol_papers)/len(LADDERS):.0f} %); "
          f"{tot_v} of {tot_p} pairs ({100*tot_v/tot_p:.0f} %)")
    sealed = [n for n, (l, f, s, c, no) in LADDERS.items() if "clamped" in s]
    print(f"  violating: {', '.join(viol_papers)}")
    print(f"  of the {len(sealed)} clamped, edge-sealed-by-design measurements "
          f"({', '.join(sealed)}),")
    print(f"  {sum(n in viol_papers for n in sealed)} violate.")

    print("\n" + "=" * 78)
    print("Q2  VALUES AT OR BELOW THEIR OWN FIXTURE FLOOR")
    print("=" * 78)
    n_at = n_check = 0
    for name, (lad, floor, seal, cond, note) in LADDERS.items():
        if floor is None:
            continue
        best_n = max(lad)
        best = lad[best_n]
        n_check += 1
        ratio = best / floor
        at = ratio <= 2.0
        n_at += at
        print(f"  {name:36s} best {best:.1e}  floor {floor:.1e}  "
              f"= {ratio:4.1f}x floor   {'AT THE FLOOR' if at else 'above floor'}")
    print(f"\n  {n_at} of {n_check} ladders with a known floor report their best "
          f"value within 2x of it.")
    print("  The remaining ladders have no reported floor, which is itself the")
    print("  finding: their best values cannot be checked by anyone, including")
    print("  their authors.")

    print("\n" + "=" * 78)
    print("FIXTURE FLOORS THE AUTHORS THEMSELVES REPORTED")
    print("=" * 78)
    for k, (v, how) in FLOORS.items():
        print(f"  {v:.1e}  {k:34s} {how}")
    vals = [v for k, (v, _) in FLOORS.items() if "inferred" not in k]
    print(f"\n  measured floors span {min(vals):.0e} to {max(vals):.0e} "
          f"g m^-2 day^-1, a factor {max(vals)/min(vals):.0f}, set entirely by how")
    print("  the sample is sealed to the instrument.  Graham's laboratory measured")
    print("  both ends of that range with one method: 5e-5 through a "
          "polyisobutylene\n  edge seal, 2e-6 with the barrier deposited straight "
          "onto the sensor.")
