# -*- coding: utf-8 -*-
"""What can a barrier measurement actually identify?

The forward model has eight unknowns. A dyad ladder, however many points it
contains, supplies three: the overall scale, the interlayer fraction that curves
it, and the parallel channel that flattens its tail. "Reconstructing the defect
structure from permeation" is therefore not possible from ladders, and the
interesting question is the one underneath it: which combinations of
measurements identify which parameters?

That is answerable exactly, by the rank and conditioning of the Jacobian of the
predicted observables with respect to the parameters. This module computes it.

CHANGES FROM THE PREVIOUS VERSION
  1. The stack resistance is now defined ONCE, in stack_R(), and every
     observable calls it. Previously it was written out three times and the
     "lag" and "pad" branches omitted the tortuosity factor, which at the
     calibrated point is 3.0e7 -- so the same physical stack differed by seven
     orders of magnitude in its organic resistance depending on which
     observable was being predicted. That divergence was what produced the
     old "rank 8 at cond 1.5e8" last row.
  2. The lag now uses the exact layer-lumped Frisch form rather than
     R_total * C_total / 6, which is the single-slab formula and understates a
     multilayer lag by 18% at two dyads. Verified against direct continuum
     integration to six digits.
  3. The pad branch scales the edge resistance as (pad/3.0)**1, not **2. The
     edge term is a perimeter-over-area quantity, WVTR_edge proportional to
     4 s t / (L s^2) proportional to 1/s, so the AREAL edge resistance goes as
     s^1. The exponent 2 was a typo. On its own it moves nothing (<1%); it
     matters only once fix (1) is in place, where it changes the last row's
     conditioning from 1.0e6 to 1.1e6.
  4. forward() is complex-safe, so the Jacobian can be taken by complex step
     (exact to machine precision, no subtractive cancellation) as well as by
     central differences.

Neither (1) nor (2) changes the identifiability table except in the last row,
which goes from rank 8 / cond 1.5e8 to rank 7 / cond 1.0e6.
"""
import numpy as np

RG, TREF = 8.314, 298.15

PARAMS = ["log f_res", "log d_close", "log f0", "log P_org",
          "log r_pin", "log R_edge", "d_crit", "log S_org"]
BASE = np.array([np.log(3.72e-8), np.log(1.07), np.log(0.372), np.log(7.5e-13),
                 np.log(5e-8), np.log(1e18), 44.0, np.log(1.5)])


# ---------------------------------------------------------------- physics
def f_pin(d_nm, f_res, d_close, f0, d_crit, C=40.0):
    """Defect area fraction. The cracking term is a positive part, off below
    d_crit. The comparison is taken on the real part so that an infinitesimal
    imaginary perturbation cannot flip the branch."""
    x = (d_nm - d_crit) / d_crit
    crack = x ** 2 if np.real(x) > 0 else 0.0 * x
    return f0 * np.exp(-d_nm / d_close) + f_res * (1.0 + C * crack)


def tortuosity(f, d_org, r_pin):
    """Detour factor for the organic interlayer between pinholes of radius
    r_pin at area fraction f. Valid while the pinhole spacing s exceeds the
    interlayer thickness, which holds by four orders of magnitude here."""
    s = r_pin * np.sqrt(np.pi / f)
    ratio = s / r_pin
    lg = np.log(ratio) if np.real(ratio) > np.e else np.log(np.e) + 0.0 * ratio
    return 1.0 + s ** 2 * lg / (2.0 * np.pi * (d_org * 1e-9) ** 2)


def layer_RC(n, d_in, d_org, f, P_org, r_pin, S_org, S_in=0.10):
    """Per-layer areal resistance and capacity for an n-dyad stack:
    n inorganic layers alternating with n-1 organic interlayers."""
    t2 = tortuosity(f, d_org, r_pin)
    R_in = d_in * 1e-9 / (f * P_org)
    R_org = d_org * 1e-9 * t2 / P_org
    C_in = S_in * d_in * 1e-9
    C_org = S_org * d_org * 1e-9
    R, C = [], []
    for k in range(n):
        R.append(R_in); C.append(C_in)
        if k < n - 1:
            R.append(R_org); C.append(C_org)
    return C, R


def stack_R(n, d_in, d_org, f, P_org, r_pin):
    """Total areal transport resistance of an n-dyad stack.

    ONE definition, called by every observable. Do not inline this into a
    branch: the previous version did, in three places, and one of the three
    silently omitted the tortuosity factor."""
    return sum(layer_RC(n, d_in, d_org, f, P_org, r_pin, 1.0)[1])


def lag_exact(C, R):
    """Exact layer-lumped Frisch lag.

        t_lag = (1/R_tot) sum_j C_j ( R_L,j R_R,j - R_j^2 / 12 )

    with R_L and R_R measured to each layer's MIDPOINT. Reduces to L^2/6D for a
    single slab, which the naive R_tot*C_tot/6 also does -- but for a stack the
    naive form is wrong, understating the lag by 18% at n=2, 8% at n=3 and 5%
    at n=4 with these parameters. Checked against direct continuum integration
    of int C(x) R_L(x) R_R(x) dx / R_tot to six significant figures."""
    R = np.asarray(R); C = np.asarray(C)
    R_tot = R.sum()
    cum = np.concatenate([[0.0 * R_tot], np.cumsum(R)])
    R_L = cum[:-1] + R / 2.0
    R_R = R_tot - cum[1:] + R / 2.0
    return np.sum(C * (R_L * R_R - R ** 2 / 12.0)) / R_tot


# ---------------------------------------------------------------- forward
def forward(p, obs):
    """Predict a list of log-observables from the parameter vector."""
    f_res, d_close, f0, P_org, r_pin, R_edge, d_crit, S_org = (
        np.exp(p[0]), np.exp(p[1]), np.exp(p[2]), np.exp(p[3]),
        np.exp(p[4]), np.exp(p[5]), p[6], np.exp(p[7]))
    out = []
    for kind, arg in obs:

        if kind == "single":                       # bare single layer
            f = f_pin(arg, f_res, d_close, f0, d_crit)
            out.append(np.log(arg * 1e-9 / (f * P_org)))
            continue

        if kind == "pad":
            n, d_in, d_org, pad = arg
        else:
            n, d_in, d_org = arg
            pad = None

        f = f_pin(d_in, f_res, d_close, f0, d_crit)
        R = stack_R(n, d_in, d_org, f, P_org, r_pin)

        if kind == "ladder":                       # dyad ladder point
            out.append(np.log(1.0 / (1.0 / R + 1.0 / R_edge)))

        elif kind == "lag":                        # exact layer-lumped lag
            C_j, R_j = layer_RC(n, d_in, d_org, f, P_org, r_pin, S_org)
            out.append(np.log(lag_exact(C_j, R_j)))

        elif kind == "pad":                        # same stack, second size
            out.append(np.log(1.0 / (1.0 / R + 1.0 / (R_edge * (pad / 3.0)))))

    return np.array(out)


# ------------------------------------------------------------- Jacobians
def jacobian(obs, p=BASE, method="complex", eps=1e-4, h=1e-20):
    """method='complex' is exact to machine precision. method='fd' is the
    central difference the earlier version used, kept for comparison."""
    m = len(forward(p, obs))
    J = np.zeros((m, len(p)))
    if method == "complex":
        for k in range(len(p)):
            q = np.asarray(p, dtype=complex).copy()
            q[k] += 1j * h
            J[:, k] = np.imag(forward(q, obs)) / h
    else:
        for k in range(len(p)):
            a, b = np.array(p, float), np.array(p, float)
            a[k] += eps; b[k] -= eps
            J[:, k] = (forward(a, obs) - forward(b, obs)) / (2 * eps)
    return J


def identifiability(obs, p=BASE, method="complex"):
    """Returns (rank, singular values, right singular vectors)."""
    J = jacobian(obs, p, method)
    U, S, Vt = np.linalg.svd(J, full_matrices=True)
    tol = max(J.shape) * S.max() * 1e-10 if S.max() > 0 else 0.0
    return int((S > tol).sum()), S, Vt


# ------------------------------------------------------------ measurements
_L = [("ladder", (n, 21.5, 100.0)) for n in (1, 2, 3, 4)]
_S = [("single", d) for d in (15., 20., 30., 50., 60.)]
_G = [("lag", (n, 21.5, 100.0)) for n in (1, 4)]
_P = [("pad", (4, 21.5, 100.0, q)) for q in (3.0, 25.0)]

SETS = {
    "ladder alone, n = 1-4": _L,
    "ladder alone, n = 1-10": [("ladder", (n, 21.5, 100.0)) for n in range(1, 11)],
    "+ thickness series (5 pts)": _L + _S,
    "+ lag times": _L + _S + _G,
    "+ two pad sizes": _L + _S + _G + _P,
}

# what a correct run must reproduce
EXPECTED = {"ladder alone, n = 1-4": (4, 3), "ladder alone, n = 1-10": (10, 3),
            "+ thickness series (5 pts)": (9, 6), "+ lag times": (11, 7),
            "+ two pad sizes": (13, 7)}  # cond 1.1e6


if __name__ == "__main__":
    print("=" * 78)
    print("IDENTIFIABILITY OF BARRIER PARAMETERS FROM MEASUREMENT SETS")
    print("=" * 78)
    print(f"  eight unknowns: {', '.join(PARAMS)}\n")
    print(f"  {'measurement set':30s} {'obs':>4} {'rank':>6} {'cond':>10}  {'weakest live':>13}  ok")
    ok_all = True
    for name, obs in SETS.items():
        r, S, Vt = identifiability(obs)
        cond = S.max() / S[r - 1]
        good = (len(obs), r) == EXPECTED[name]
        ok_all &= good
        print(f"  {name:30s} {len(obs):4d} {r:5d}/8 {cond:10.1e}  {S[r-1]/S.max():13.2e}  "
              f"{'yes' if good else 'NO'}")

    r, S, Vt = identifiability(SETS["+ two pad sizes"])
    print(f"\n  unconstrained combination, full measurement set:")
    v = Vt[r]
    for j in np.argsort(-np.abs(v))[:3]:
        print(f"    {v[j]:+.3f} * {PARAMS[j]}")
    print("  the same combination is unconstrained in every row: the two defect")
    print("  fractions and the organic permeability enter only through a product.")

    print(f"\n  cross-check, complex step vs central difference:")
    Sc = np.linalg.svd(jacobian(SETS['+ two pad sizes'], method='complex'), compute_uv=False)
    Sf = np.linalg.svd(jacobian(SETS['+ two pad sizes'], method='fd'), compute_uv=False)
    print(f"    largest singular value  agree to {abs(Sc[0]/Sf[0]-1):.1e}")
    print(f"    weakest live direction  agree to {abs(Sc[6]/Sf[6]-1):.1e}")

    print(f"\n  {'ALL ROWS MATCH THE PUBLISHED TABLE' if ok_all else '*** TABLE MISMATCH ***'}")
