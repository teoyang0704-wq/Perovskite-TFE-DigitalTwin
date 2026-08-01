# -*- coding: utf-8 -*-
"""Validation and first experiments for the 3-D defect-network module."""
import numpy as np
from percolation3d import Stack, place_defects, effective_conductance

rng0 = np.random.default_rng(7)


# ---------------------------------------------------------------- test 1
def test_two_electrode(nx_list=(64, 128, 256)):
    """One pinhole per layer, separation rho: compare with the analytic
    two-electrode sheet result  R = 2 R_through + ln(rho/r)/(pi P d)."""
    st = Stack(n_inorg=2, d_in=50e-9, d_org=500e-9, r_pin=50e-9, P_org=7.5e-13)
    L, rho = 400e-6, 50e-6
    R_through = st.d_in / (st.P_org * np.pi * st.r_pin ** 2)
    R_lat = np.log(rho / st.r_pin) / (np.pi * st.P_org * st.d_org)
    R_ana = 2 * R_through + R_lat
    print(f"  analytic: R_through={R_through:.2e}  R_lat={R_lat:.2e}  total={R_ana:.2e}")
    for nx in nx_list:
        layers = [np.array([[L / 2, L / 2]]), np.array([[L / 2 + rho, L / 2]])]
        G = effective_conductance(st, L, nx, layers)
        print(f"  nx={nx:4d}  R_network={1/G:.3e}   ratio to analytic = {R_ana*G:.3f}")


# ---------------------------------------------------------------- test 2
def lattice_vs_tau2(st, n_realise=1, nx=192, L_over_s=8.0):
    """Regular array, phi=0: does the network reproduce the Year-1 tau^2 form?"""
    L = L_over_s * st.spacing
    layers, col = place_defects(st, L, rng0, phi=0.0, arrangement="lattice")
    G = effective_conductance(st, L, nx, layers, col)
    G1 = st.conductance_1d(L * L)
    print(f"  N per layer={len(layers[0]):4d}  s={st.spacing*1e6:6.1f} um  "
          f"tau^2={st.tau2():.3e}")
    print(f"  G_network/G_1D = {G/G1:.3f}   (1.0 = analytic formula exact)")
    return G / G1


# ---------------------------------------------------------------- test 3
def randomness_bias(st, n_realise=8, nx=192, L_over_s=8.0):
    """Poisson vs regular placement at the same density, phi = 0."""
    L = L_over_s * st.spacing
    lay0, col0 = place_defects(st, L, rng0, 0.0, "lattice")
    lat = effective_conductance(st, L, nx, lay0, col0)
    vals = []
    for k in range(n_realise):
        rng = np.random.default_rng(100 + k)
        lay, col = place_defects(st, L, rng, 0.0, "poisson")
        vals.append(effective_conductance(st, L, nx, lay, col))
    v = np.array(vals) / lat
    print(f"  Poisson/lattice conductance ratio = {v.mean():.3f} +- {v.std():.3f}"
          f"  (n={n_realise})")
    return v.mean(), v.std()


# ---------------------------------------------------------------- test 4
def correlation_ladder(st_base, phis=(0.0, 0.003, 0.01, 0.03, 0.1),
                       n_max=4, n_realise=4, nx=160, L_over_s=8.0):
    """WVTR ladder versus dyad number for several columnar-defect fractions."""
    L = L_over_s * st_base.spacing
    print(f"  domain {L*1e3:.2f} mm, {int(st_base.density*L*L)} defects/layer")
    out = {}
    for phi in phis:
        row = []
        for n in range(2, n_max + 1):
            st = Stack(n_inorg=n, d_in=st_base.d_in, d_org=st_base.d_org,
                       f=st_base.f, r_pin=st_base.r_pin, P_org=st_base.P_org)
            g = []
            for k in range(n_realise):
                rng = np.random.default_rng(1000 + 37 * k + n)
                lay, col = place_defects(st, L, rng, phi, "poisson")
                g.append(effective_conductance(st, L, nx, lay, col))
            row.append(np.mean(g))
        out[phi] = np.array(row)
        rel = out[phi] / out[phi][0]
        print(f"  phi={phi:<6.3f}  G(n)/G(2) = " + "  ".join(f"{x:.3f}" for x in rel)
              + f"   [G(2)={out[phi][0]:.2e}]")
    return out


if __name__ == "__main__":
    print("\n[1] solver check against the two-electrode analytic result")
    test_two_electrode()

    wu = Stack(n_inorg=3, d_in=50e-9, d_org=500e-9, f=3.72e-8,
               r_pin=50e-9, P_org=7.5e-13)
    print("\n[2] regular lattice versus the Year-1 tau^2 formula (Wu geometry)")
    lattice_vs_tau2(wu)

    print("\n[3] effect of random (Poisson) placement at equal density")
    randomness_bias(wu)

    lee = Stack(n_inorg=4, d_in=21.5e-9, d_org=100e-9, f=5.6e-7,
                r_pin=50e-9, P_org=7.5e-13)
    print("\n[4] dyad ladder versus columnar-defect fraction (Lee geometry)")
    correlation_ladder(lee)
