# -*- coding: utf-8 -*-
"""Closing the disjunction with published data.

A ladder that violates the series bound admits two explanations: a transport
path parallel to the stack, or barrier layers whose quality falls with
deposition order.  The manuscript has treated this as open, pending an
experiment that inverts the deposition order.

That experiment is the cleanest test, but it is not the only one.  The two
hypotheses differ in every observable that depends on something other than the
stack, and three of those observables are already in the literature.

  T3  Lag time.  A gradient changes how much resistance each layer adds; it
      cannot change the fact that the lag is the stack's own transient, of order
      10 h.  A lateral path has a lag set by its own length, of order 10^3 h.

  T5  Sealing arrangement.  A gradient is a property of the deposited stack and
      is unaffected by how the sample is attached to the instrument.  A parallel
      path through the seal changes when the seal changes.

  T6  Single-layer consistency.  A gradient acts only from the second layer
      onward -- the first layer is by definition the reference.  A parallel
      channel contaminates every reading including n = 1, and by a predictable
      amount.

Each is evaluated below against published numbers.
"""
import numpy as np

DAY, KG, DA = 86400.0, 1e3, 0.90
P_PAR = 7.5e-13
D_IN, D_ORG, D_FIRST, R_PIN, S_PAR = 21.5e-9, 100e-9, 800e-9, 50e-9, 1.5
PAD, PERIM, L_PATH = 3e-3, 4 * 3e-3, 5e-3
LEE = {1: 3.0e-3, 2: 6.6e-4, 3: 5.4e-4, 4: 5.3e-4}
M = 22.9                      # organic permeability multiplier from the flux fit
F_FIT = 4.2e-8                # defect fraction from the joint fit


# ----------------------------------------------------------------- T3
def t3():
    print("=" * 74)
    print("T3  LAG TIME  --  can a gradient produce a 400 h lag?")
    print("=" * 74)
    P = M * P_PAR
    s = R_PIN * np.sqrt(np.pi / F_FIT)
    tau2 = 1 + s ** 2 * np.log(s / R_PIN) / (2 * np.pi * D_ORG ** 2)

    def stack_lag(n, degrade=1.0):
        """Frisch lag with layer k resistance scaled by degrade**k."""
        r, c = [], []
        for i in range(n):
            r.append(D_IN / (F_FIT * P) * degrade ** i); c.append(0.10 * D_IN)
            if i < n - 1:
                r.append(D_ORG * tau2 / P); c.append(S_PAR * D_ORG)
        r.append(D_FIRST / P); c.append(S_PAR * D_FIRST)
        r, c = np.array(r), np.array(c)
        R = r.sum(); RL = np.cumsum(r) - r / 2
        return float((c * RL * (R - RL)).sum() / R) / 3600.0

    print("  A gradient is modelled as layer k having resistance scaled by g**k.")
    print("  The value of g needed to reproduce the flux ladder is about 1/34 per")
    print("  layer; larger g means a milder gradient.\n")
    print(f"  {'gradient g':>12} {'lag at n=4 (h)':>16}")
    for g in (1.0, 0.5, 0.1, 1/34, 0.01):
        print(f"  {g:12.3f} {stack_lag(4, g):16.1f}")
    lag_edge = L_PATH ** 2 / (6.0 * (M * P_PAR / S_PAR)) / 3600.0
    print(f"\n  lateral path                {lag_edge:8.0f}")
    print(f"  measured (Lee et al.)       {'>400':>8}")
    print("\n  No gradient reproduces the measured lag: making later layers worse")
    print("  *reduces* the stack's transient, because the resistance-capacitance")
    print("  product that sets the lag shrinks. The gradient hypothesis moves the")
    print("  prediction away from the measurement, not toward it.")


# ----------------------------------------------------------------- T5
def t5():
    print("\n" + "=" * 74)
    print("T5  SEALING ARRANGEMENT  --  one laboratory, one method, two seals")
    print("=" * 74)
    print("  Graham (54th SVC Proc., 2011) reports, from the same laboratory and")
    print("  the same calcium method:")
    print("     barrier deposited directly on the sensor, no lid   2e-6 g/m2/day")
    print("     lid attached with a polyisobutylene edge seal      5e-5 g/m2/day")
    print(f"     ratio                                              {5e-5/2e-6:.0f}x")
    print("\n  The stacks are made by the same group with the same processes. A")
    print("  gradient in layer quality is a property of the deposited stack and")
    print("  cannot change by a factor of 25 when the sample is attached to the")
    print("  instrument differently. A parallel path through the seal does exactly")
    print("  that, and its magnitude tracks the sealant: polyisobutylene, the")
    print("  material chosen industrially for photovoltaic edge seals precisely")
    print("  because it is tight, still admits 25 times more than no seal at all.")
    print("\n  This observation is independent of any model in this work.")


# ----------------------------------------------------------------- T6
def t6():
    print("\n" + "=" * 74)
    print("T6  SINGLE-LAYER CONSISTENCY  --  is the n = 1 reading contaminated too?")
    print("=" * 74)
    print("  A gradient defines the first layer as the reference: by construction")
    print("  it predicts nothing about the n = 1 value. A parallel channel predicts")
    print("  the n = 1 value is contaminated by the same edge term as the others.\n")
    t_org1 = D_FIRST
    edge1 = M * P_PAR * DA * (PERIM * t_org1) / (L_PATH * PAD ** 2) * KG * DAY
    print(f"  edge contribution at n = 1 (no free parameters left): {edge1:.2e}")
    print(f"  measured n = 1:                                       {LEE[1]:.2e}")
    print(f"  edge share at n = 1:                                  {100*edge1/LEE[1]:.0f} %")
    print("\n  The prediction is that the single-dyad value is only ~10 % edge, so")
    print("  the stack dominates there and the ladder's first step is nearly clean.")
    print("  That is what makes the first step steep (a factor 4.5 from n=1 to n=2)")
    print("  while later steps flatten. A gradient has no mechanism that makes the")
    print("  first step steep and the later ones flat; it must be tuned to do so.")

    print("\n  Quantitatively: the observed ladder falls 4.5x from n=1 to n=2, then")
    print("  1.2x and 1.02x. The edge model reproduces the whole shape with two")
    print("  parameters (rms 0.011 decades). A gradient needs one parameter per")
    print("  step and reproduces it exactly by construction, predicting nothing.")


# ----------------------------------------------------------------- verdict
def verdict():
    print("\n" + "=" * 74)
    print("WHAT THE THREE TESTS TOGETHER SUPPORT")
    print("=" * 74)
    print("  T3  A gradient cannot produce the measured lag; it predicts a shorter")
    print("      transient than a uniform stack, while the measurement is two orders")
    print("      longer. The lateral path predicts the right order with no freedom.")
    print("  T5  The floor moves by 25x when only the seal changes, within one")
    print("      laboratory. A stack property cannot do that.")
    print("  T6  The parallel channel predicts a clean first step and flat later")
    print("      ones, which is the observed shape; a gradient must be fitted to it.")
    print("\n  These are independent of one another and of the flux fit. Taken")
    print("  together they exclude a pure quality gradient as the explanation of the")
    print("  observed ladders. What they do not exclude is a gradient existing")
    print("  *alongside* the edge channel and contributing part of the flattening;")
    print("  the deposition-order experiment remains the way to bound that share.")


if __name__ == "__main__":
    t3(); t5(); t6(); verdict()
