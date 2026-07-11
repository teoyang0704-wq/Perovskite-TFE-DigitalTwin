# -*- coding: utf-8 -*-
"""
=====================================================================
 tfe_physics_engine.py
 TFE Digital Twin - Physics Engine v1.0   (Project Step 1)
---------------------------------------------------------------------
 1-D coupled moisture-diffusion / heat-conduction PDE solver for
 multilayer thin-film encapsulation (TFE) on perovskite solar cells.

 [Moisture]  Fick's 2nd law in water-ACTIVITY form
        S(x) * da/dt = d/dx( P(x,T) * da/dx ),   P = D * S
   - a = C/S (water activity, 0..1) is continuous across layer
     interfaces  ->  physically correct multilayer interface
     condition (continuity of water chemical potential). At the
     ambient boundary a simply equals the relative humidity RH(t).
   - Arrhenius temperature coupling:
        D(T) = D_ref * exp[ -Ea/R * (1/T - 1/T_ref) ]
   - DEFECT-MEDIATED EFFECTIVE MEDIA  (this is what makes pure
     GEOMETRY matter in a 1-D model -- in a perfect-film 1-D model
     the WVTR would depend only on total thickness, not n_pairs):
       inorganic layer (parallel paths, lattice + pinholes):
           P_in_eff(T) = D_lat(T)*S_in + f_pin * D_org(T)*S_org
           f_pin(d) = f0*exp(-d/d_close)      (nucleation closure)
                      + f_res                 (particle defect floor)
                      + f_c0*max(0,(d-d_crit)/d_crit)^2
                                              (channel-cracking penalty)
       organic decoupling layer (series tortuous path between
       laterally staggered pinholes; Graff et al., J. Appl. Phys.
       96, 1840 (2004) lag-time picture, folded into 1-D as a
       tortuosity correction):
           P_org_eff(T) = D_org(T)*S_org / tau^2
           tau^2 = 1 + s^2 * ln(s/r_pin) / (2*pi*d_org^2)
           s = mean pinhole spacing = r_pin * sqrt(pi / f_pin)

 [Heat]      rho*cp * dT/dt = d/dx( k dT/dx ) + q(x,t)
   - domain: substrate | device | TFE (sun-facing TFE)
   - solar heat source in the device layer, parasitic absorption
     in the TFE, convection at both faces + sky radiation on top.

 [Numerics]  cell-centred finite volume, harmonic-mean interface
   conductances, fully-implicit backward Euler (tridiagonal solve
   via scipy.linalg.solve_banded)
   -> unconditionally stable for the >1e5 permeability contrast
      between organic and inorganic layers.

 Units: SI internally. User-facing: nm, hours, g/m2/day, degC.
 Dependencies: numpy, scipy, matplotlib  (no FiPy / COMSOL needed).
=====================================================================
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Callable, List, Tuple, Dict, Optional
from scipy.linalg import solve_banded

__version__ = "1.0"

R_GAS = 8.314        # J/(mol K)
SIGMA = 5.670e-8     # W/(m2 K4)
T_REF = 298.15       # K
DAY   = 86400.0      # s


# ------------------------------------------------------------------
# 1. Materials
# ------------------------------------------------------------------
@dataclass
class Material:
    name: str
    D_ref: float        # water diffusivity at T_REF          [m2/s]
    Ea_D: float         # diffusion activation energy         [J/mol]
    S: float            # water solubility at unit activity   [kg/m3]
    k: float            # thermal conductivity                [W/mK]
    rho: float          # density                             [kg/m3]
    cp: float           # specific heat                       [J/kgK]
    n_opt: float        # refractive index (visible mean)     [-]
    alpha_opt: float    # optical absorption coefficient      [1/m]
    cte: float          # thermal expansion coefficient       [1/K]
    E_mod: float        # Young's modulus                     [Pa]
    nu: float           # Poisson ratio                       [-]
    cost_per_nm: float  # cost coefficient                    [arb/nm]
    emissivity: float = 0.85

    def D(self, T: float) -> float:
        """Arrhenius diffusivity."""
        return self.D_ref * np.exp(-self.Ea_D / R_GAS * (1.0 / T - 1.0 / T_REF))


# NOTE: D_ref of the organic equals the user's 1e-8 cm2/s = 1e-12 m2/s.
#       For the inorganic, the user's 1e-13 cm2/s was an EFFECTIVE value
#       (lattice + defects). Here we separate the two: dense ALD Al2O3
#       lattice diffusion is immeasurably small (< 1e-20 m2/s), so ALL
#       practical transport goes through pinholes/defects (DefectModel).
#       This separation is precisely what gives the model its geometry
#       (n_pairs / thickness) sensitivity.
ORGANIC = Material("organic (acrylate/pV3D3)",
                   D_ref=1.0e-12, Ea_D=30e3, S=10.0,
                   k=0.20, rho=1200, cp=1500,
                   n_opt=1.52, alpha_opt=2e4,
                   cte=60e-6, E_mod=3e9, nu=0.35,
                   cost_per_nm=0.10, emissivity=0.88)

INORGANIC = Material("inorganic (ALD Al2O3)",
                     D_ref=1.0e-21, Ea_D=60e3, S=0.10,
                     k=1.80, rho=3000, cp=880,
                     n_opt=1.65, alpha_opt=1e3,
                     cte=5e-6, E_mod=150e9, nu=0.24,
                     cost_per_nm=0.50, emissivity=0.65)

PET = Material("PET substrate", 1e-13, 40e3, 5.0, 0.20, 1380, 1200,
               1.57, 5e2, 30e-6, 4e9, 0.40, 0.0, 0.90)

DEVICE = Material("perovskite device", 1e-14, 50e3, 1.0, 0.50, 4000, 300,
                  2.40, 1e7, 40e-6, 20e9, 0.30, 0.0, 0.90)


# ------------------------------------------------------------------
# 2. Defect / tortuosity model  (geometry <-> transport link)
# ------------------------------------------------------------------
@dataclass
class DefectModel:
    """Pinhole statistics of the inorganic barrier layers.

    All four parameters are CALIBRATION handles against experimental
    Ca-test / WVTR data. Defaults reproduce literature orders of
    magnitude for ALD Al2O3 on polymer.
    """
    f0: float = 1e-3          # pinhole area fraction of an ultrathin film
    d_close_nm: float = 5.0   # nucleation-closure decay length [nm]
    f_res: float = 1e-9       # residual (particle) defect floor
    #   f_res = 1e-9 with r = 50 nm corresponds to ~0.1 pinhole/mm2;
    #   f(50 nm) ~ 5e-8 corresponds to ~2 pinholes/mm2 and ~0.4 mm
    #   spacing -- the range reported from Ca-test decoration of good
    #   ALD Al2O3. Together the defaults put a 3x(50/200 nm) dyad stack
    #   in the 1e-3..1e-4 g/m2/day WVTR class (85C / 25C respectively).
    f_crack0: float = 3e-5    # cracking penalty prefactor
    d_crit_nm: float = 120.0  # channel-cracking critical thickness [nm]
    r_pinhole: float = 50e-9  # effective pinhole radius [m]

    def pinhole_fraction(self, d_inorg_nm: float) -> float:
        d = float(d_inorg_nm)
        f_crack = self.f_crack0 * max(0.0, (d - self.d_crit_nm) / self.d_crit_nm) ** 2
        return self.f0 * np.exp(-d / self.d_close_nm) + self.f_res + f_crack

    def lateral_spacing(self, d_inorg_nm: float) -> float:
        """Mean centre-to-centre pinhole spacing s = r*sqrt(pi/f)."""
        f = self.pinhole_fraction(d_inorg_nm)
        return self.r_pinhole * np.sqrt(np.pi / f)


# ------------------------------------------------------------------
# 3. Generic 1-D implicit finite-volume diffusion solver
#    C(x) du/dt = d/dx( K(x) du/dx ) + src(x)
# ------------------------------------------------------------------
class ImplicitDiffusion1D:
    """Backward-Euler FV solver on a fixed (possibly non-uniform) grid.

    Used with (u=a, K=P=D*S, C=S)          for moisture, and
              (u=T, K=k,     C=rho*cp)     for heat.
    """

    def __init__(self, dx: np.ndarray):
        self.dx = np.asarray(dx, dtype=float)
        self.N = len(self.dx)

    @staticmethod
    def _bc(dx_b: float, K_b: float, bc: tuple) -> Tuple[float, float]:
        """Return (boundary conductance G_b, reservoir value u_b)."""
        kind = bc[0]
        if kind == "dirichlet":                    # ("dirichlet", value)
            return 2.0 * K_b / dx_b, bc[1]
        if kind == "robin":                        # ("robin", h, u_inf)
            h = bc[1]
            return 1.0 / (1.0 / h + dx_b / (2.0 * K_b)), bc[2]
        if kind == "adiabatic":                    # ("adiabatic",)
            return 0.0, 0.0
        raise ValueError(f"unknown BC {kind}")

    def step(self, u, dt, K, C, bc_left, bc_right, src=None):
        """One implicit time step. Returns (u_new, J_left, J_right),
        where J is the flux LEAVING the domain through that face
        (units of [K]*[u]/m, e.g. kg/m2/s for moisture, W/m2 for heat)."""
        dx, N = self.dx, self.N
        # harmonic-mean interface conductances (exact flux continuity)
        G = 1.0 / (dx[:-1] / (2.0 * K[:-1]) + dx[1:] / (2.0 * K[1:]))

        diag = C * dx / dt
        rhs = C * dx / dt * u
        if src is not None:
            rhs = rhs + src * dx

        diag[:-1] += G
        diag[1:] += G

        Gl, ul = self._bc(dx[0], K[0], bc_left)
        Gr, ur = self._bc(dx[-1], K[-1], bc_right)
        diag[0] += Gl
        rhs[0] += Gl * ul
        diag[-1] += Gr
        rhs[-1] += Gr * ur

        ab = np.zeros((3, N))
        ab[0, 1:] = -G          # super-diagonal
        ab[1, :] = diag         # diagonal
        ab[2, :-1] = -G         # sub-diagonal
        u_new = solve_banded((1, 1), ab, rhs)

        J_left = Gl * (u_new[0] - ul)
        J_right = Gr * (u_new[-1] - ur)
        return u_new, J_left, J_right


# ------------------------------------------------------------------
# 4. TFE stack geometry -> cell arrays
# ------------------------------------------------------------------
class TFEStack:
    """Cell arrays of the TFE. Index 0 touches the device,
    index -1 touches the ambient. Layer order (device -> ambient):
        [ inorganic | organic ] * n_pairs
    The topmost organic faces the ambient directly, hence no lateral
    pinhole detour (tau = 1) is applied to it."""

    def __init__(self, d_org_nm, d_inorg_nm, n_pairs,
                 defects: Optional[DefectModel] = None, cells_per_layer=8):
        self.d_org_nm = float(d_org_nm)
        self.d_inorg_nm = float(d_inorg_nm)
        self.n_pairs = int(round(n_pairs))
        self.defects = defects or DefectModel()
        self.cpl = int(cells_per_layer)

        self.f_pin = self.defects.pinhole_fraction(self.d_inorg_nm)
        self.s_lat = self.defects.lateral_spacing(self.d_inorg_nm)

        d_org_m = self.d_org_nm * 1e-9
        d_in_m = self.d_inorg_nm * 1e-9
        ratio = max(self.s_lat / self.defects.r_pinhole, np.e)
        self.tau2_sand = 1.0 + self.s_lat ** 2 * np.log(ratio) / (2.0 * np.pi * d_org_m ** 2)

        dx, is_org, tau2, edges = [], [], [], [0.0]
        for p in range(self.n_pairs):
            dx += [d_in_m / self.cpl] * self.cpl        # inorganic
            is_org += [False] * self.cpl
            tau2 += [1.0] * self.cpl
            edges.append(edges[-1] + d_in_m)
            top = (p == self.n_pairs - 1)               # organic
            dx += [d_org_m / self.cpl] * self.cpl
            is_org += [True] * self.cpl
            tau2 += [1.0 if top else self.tau2_sand] * self.cpl
            edges.append(edges[-1] + d_org_m)

        self.dx = np.array(dx)
        self.is_org = np.array(is_org)
        self.tau2 = np.array(tau2)
        self.layer_edges = np.array(edges)
        self.x_edges = np.concatenate([[0.0], np.cumsum(self.dx)])
        self.x_cent = 0.5 * (self.x_edges[:-1] + self.x_edges[1:])
        self.L = self.x_edges[-1]

        self.S = np.where(self.is_org, ORGANIC.S, INORGANIC.S)
        self.k_th = np.where(self.is_org, ORGANIC.k, INORGANIC.k)
        self.C_th = np.where(self.is_org, ORGANIC.rho * ORGANIC.cp,
                             INORGANIC.rho * INORGANIC.cp)

    # effective permeability field P(x,T) = D_eff * S  [kg/(m s)]
    def P_cells(self, T: float) -> np.ndarray:
        P_org = ORGANIC.D(T) * ORGANIC.S
        P_in = INORGANIC.D(T) * INORGANIC.S
        return np.where(self.is_org,
                        P_org / self.tau2,                # tortuous path
                        P_in + self.f_pin * P_org)        # lattice + pinholes

    def inorg_spans_nm(self) -> List[Tuple[float, float]]:
        e = self.layer_edges * 1e9
        return [(e[2 * p], e[2 * p + 1]) for p in range(self.n_pairs)]


# ------------------------------------------------------------------
# 5. Environment scenarios (dynamic boundary forcing)
# ------------------------------------------------------------------
@dataclass
class Environment:
    name: str
    T_amb: Callable[[float], float]     # [K]
    RH: Callable[[float], float]        # activity 0..1
    G_sun: Callable[[float], float]     # [W/m2]


def damp_heat_env(T_C=85.0, RH=0.85) -> Environment:
    """IEC 61215-style constant damp-heat chamber."""
    return Environment(f"damp heat {T_C:.0f}C/{RH*100:.0f}%RH",
                       lambda t: T_C + 273.15, lambda t: RH, lambda t: 0.0)


def diurnal_cycle(T_min_C=15.0, T_max_C=42.0, RH_min=0.30, RH_max=0.90,
                  G_max=1000.0, t_peak_h=14.0) -> Environment:
    """Sinusoidal day/night: T peaks at t_peak_h, RH anti-correlated,
    solar flux is a half-sine between 06:00 and 18:00."""
    Tm = 0.5 * (T_min_C + T_max_C) + 273.15
    Ta = 0.5 * (T_max_C - T_min_C)
    Rm = 0.5 * (RH_min + RH_max)
    Ra = 0.5 * (RH_max - RH_min)

    def T_amb(t):
        return Tm + Ta * np.cos(2 * np.pi * (t / DAY - t_peak_h / 24.0))

    def RH(t):
        r = Rm - Ra * np.cos(2 * np.pi * (t / DAY - t_peak_h / 24.0))
        return float(np.clip(r, 0.05, 0.98))

    def G(t):
        h = (t / 3600.0) % 24.0
        return G_max * max(0.0, np.sin(np.pi * (h - 6.0) / 12.0))

    return Environment("diurnal cycle", T_amb, RH, G)


# ------------------------------------------------------------------
# 6. Digital twin: coupled runs + 6-objective KPI extraction
# ------------------------------------------------------------------
class TFEDigitalTwinPDE:
    """PDE-based replacement of the analytic ``simulate_tfe``.

    Design variables:  d_org_nm, d_inorg_nm, n_pairs   (pure geometry)
    ``evaluate_design()`` returns all 6 objectives -> ready to be
    wrapped by pymoo/NSGA-II in project Step 2.
    """

    def __init__(self, d_org_nm, d_inorg_nm, n_pairs,
                 defects: Optional[DefectModel] = None,
                 substrate: Material = PET, sub_thk=125e-6, dev_thk=500e-9,
                 cells_per_layer=8, n_sub_cells=12, n_dev_cells=2,
                 h_top=12.0, h_bot=6.0, eta_heat=0.75):
        self.stack = TFEStack(d_org_nm, d_inorg_nm, n_pairs,
                              defects, cells_per_layer)
        self.substrate, self.sub_thk, self.dev_thk = substrate, sub_thk, dev_thk
        self.h_top, self.h_bot, self.eta_heat = h_top, h_bot, eta_heat

        # thermal grid: [substrate | device | TFE], x=0 = module back
        dx = [sub_thk / n_sub_cells] * n_sub_cells + [dev_thk / n_dev_cells] * n_dev_cells
        kk = [substrate.k] * n_sub_cells + [DEVICE.k] * n_dev_cells
        CC = [substrate.rho * substrate.cp] * n_sub_cells + \
             [DEVICE.rho * DEVICE.cp] * n_dev_cells
        self.th_dx = np.concatenate([dx, self.stack.dx])
        self.th_k = np.concatenate([kk, self.stack.k_th])
        self.th_C = np.concatenate([CC, self.stack.C_th])
        self.i_dev = slice(n_sub_cells, n_sub_cells + n_dev_cells)
        self.i_tfe = slice(n_sub_cells + n_dev_cells, None)
        self.eps_top = ORGANIC.emissivity if self.stack.is_org[-1] \
            else INORGANIC.emissivity

        self._moist = ImplicitDiffusion1D(self.stack.dx)
        self._heat = ImplicitDiffusion1D(self.th_dx)

    # ---------------- moisture-only run (isothermal chamber) --------
    def run_damp_heat(self, T_C=85.0, RH=0.85, t_end_h=600.0,
                      dt0=0.5, dt_max=600.0, ramp=1.15,
                      snap_frac=(0.002, 0.02, 0.1, 0.4, 1.0)) -> Dict:
        """Fick's-2nd-law PDE under constant T/RH.
        Time step ramps geometrically (implicit scheme = stable)."""
        st = self.stack
        T = T_C + 273.15
        K = st.P_cells(T)                       # constant (isothermal)
        a = np.zeros(len(st.dx))
        t_end = t_end_h * 3600.0
        snap_t = [f * t_end for f in snap_frac]
        snaps: List[Tuple[float, np.ndarray]] = []

        t, dt, M = 0.0, dt0, 0.0
        ts, Js, Ms = [0.0], [0.0], [0.0]
        while t < t_end - 1e-9:
            dtn = min(dt, t_end - t)
            a, Jl, _ = self._moist.step(a, dtn, K, st.S,
                                        ("dirichlet", 0.0),      # device = perfect sink (Ca-test)
                                        ("dirichlet", RH))       # ambient activity = RH
            t += dtn
            M += Jl * dtn
            ts.append(t); Js.append(Jl); Ms.append(M)
            while snap_t and t >= snap_t[0] - 1e-9:
                snaps.append((snap_t.pop(0), a.copy()))
            dt = min(dt * ramp, dt_max)

        ts = np.array(ts); Js = np.array(Js); Ms = np.array(Ms)
        J_ss = float(np.mean(Js[int(0.9 * len(Js)):]))           # steady flux
        flat = float(Js[-1] / J_ss) if J_ss > 0 else np.nan
        t_lag = t_end - Ms[-1] / max(J_ss, 1e-300)               # time-lag method
        return dict(t=ts, wvtr=Js * 1e3 * 86400.0,               # g/m2/day
                    M=Ms * 1e3,                                  # g/m2
                    wvtr_ss=J_ss * 1e3 * 86400.0,
                    steady_flatness=flat, t_lag_h=t_lag / 3600.0,
                    snapshots=snaps, x_nm=st.x_cent * 1e9,
                    T_C=T_C, RH=RH)

    @staticmethod
    def t80_hours(run: Dict, M_crit_g_m2=0.01) -> float:
        """T80 = time for cumulative ingress to reach the critical
        water dose M_crit (calibration parameter of the perovskite)."""
        M, t = run["M"], run["t"]
        if M[-1] >= M_crit_g_m2:
            return float(np.interp(M_crit_g_m2, M, t)) / 3600.0
        extra_days = (M_crit_g_m2 - M[-1]) / max(run["wvtr_ss"], 1e-300)
        return t[-1] / 3600.0 + extra_days * 24.0

    def steady_activity_profile(self, RH: float, T: float) -> np.ndarray:
        """Analytic steady-state activity profile a(x) from the series
        resistance ladder (used to pre-condition long-horizon runs)."""
        st = self.stack
        r_half = st.dx / (2.0 * st.P_cells(T))       # half-cell resistances
        Rc = np.empty(len(st.dx))
        Rc[0] = r_half[0]
        for i in range(1, len(st.dx)):
            Rc[i] = Rc[i - 1] + r_half[i - 1] + r_half[i]
        R_tot = Rc[-1] + r_half[-1]
        return RH * Rc / R_tot

    # ---------------- coupled diurnal run (heat <-> moisture) -------
    def run_diurnal(self, env: Environment, days=3.0, dt=60.0,
                    init="steady") -> Dict:
        """Operator-split coupling per step:
        (1) implicit heat step with solar source + radiation,
        (2) update Arrhenius D(T) with mean TFE temperature,
        (3) implicit moisture step with dynamic RH(t) boundary.
        init='steady' pre-conditions the moisture field so the run
        shows the periodic 'breathing' regime rather than the initial
        lag transient (which run_damp_heat characterises)."""
        st = self.stack
        Nth = len(self.th_dx)
        T = np.full(Nth, env.T_amb(0.0))
        if init == "steady":
            RH0 = 0.5 * (env.RH(0.0) + env.RH(DAY / 2))
            a = self.steady_activity_profile(RH0, float(np.mean(
                [env.T_amb(0.0), env.T_amb(DAY / 2)])))
        else:
            a = np.zeros(len(st.dx))
        T_opt = self.transmittance()
        M, t, t_end = 0.0, 0.0, days * DAY
        rec = {k: [] for k in ("t", "T_dev", "T_tfe", "T_amb", "G", "RH",
                               "wvtr", "wvtr_uptake", "M")}
        while t < t_end - 1e-6:
            t += dt
            Ta, G, RH = env.T_amb(t), env.G_sun(t), env.RH(t)
            # -- heat --
            q = np.zeros(Nth)
            q[self.i_dev] += G * T_opt * self.eta_heat / self.dev_thk
            q[self.i_tfe] += G * (1.0 - T_opt) / st.L        # parasitic abs.
            Tsky = Ta - 15.0                                  # clear-sky approx
            q[-1] += self.eps_top * SIGMA * (Tsky ** 4 - T[-1] ** 4) / self.th_dx[-1]
            T, _, _ = self._heat.step(T, dt, self.th_k, self.th_C,
                                      ("robin", self.h_bot, Ta),
                                      ("robin", self.h_top, Ta), src=q)
            # -- moisture with D(T) --
            T_tfe = float(np.mean(T[self.i_tfe]))
            Km = st.P_cells(T_tfe)
            a, Jl, Jr = self._moist.step(a, dt, Km, st.S,
                                         ("dirichlet", 0.0), ("dirichlet", RH))
            M += Jl * dt
            rec["t"].append(t); rec["T_dev"].append(float(np.mean(T[self.i_dev])))
            rec["T_tfe"].append(T_tfe); rec["T_amb"].append(Ta)
            rec["G"].append(G); rec["RH"].append(RH)
            rec["wvtr"].append(Jl * 1e3 * 86400.0)          # into device
            rec["wvtr_uptake"].append(-Jr * 1e3 * 86400.0)  # from ambient
            rec["M"].append(M * 1e3)
        return {k: np.array(v) for k, v in rec.items()}

    # ---------------- analytic KPIs ---------------------------------
    def transmittance(self) -> float:
        """Incoherent estimate: Fresnel interface losses x Beer-Lambert.
        (Coherent transfer-matrix upgrade planned for a later step.)"""
        st = self.stack
        n_seq = [1.0]                                # air
        for _ in range(st.n_pairs):                  # from ambient downwards
            n_seq += [ORGANIC.n_opt, INORGANIC.n_opt]
        n_seq.append(DEVICE.n_opt)
        Topt = 1.0
        for na, nb in zip(n_seq[:-1], n_seq[1:]):
            Topt *= 1.0 - ((na - nb) / (na + nb)) ** 2
        ab = (ORGANIC.alpha_opt * st.d_org_nm +
              INORGANIC.alpha_opt * st.d_inorg_nm) * 1e-9 * st.n_pairs
        return float(Topt * np.exp(-ab))

    def cost(self) -> float:                         # user's cost model kept
        st = self.stack
        return (st.d_org_nm * ORGANIC.cost_per_nm +
                st.d_inorg_nm * INORGANIC.cost_per_nm) * st.n_pairs

    def added_weight_g_m2(self) -> float:
        st = self.stack
        return (ORGANIC.rho * st.d_org_nm +
                INORGANIC.rho * st.d_inorg_nm) * 1e-9 * st.n_pairs * 1e3

    def durability(self, dT_cycle: float,
                   sigma_crit_inorg=350e6, sigma_crit_org=60e6,
                   m_fatigue=8.0) -> Dict:
        """Thermo-mechanical fatigue from the diurnal Delta-T:
        biaxial film stress amplitude sigma_a = E/(1-nu)*|d_alpha|*dT/2,
        Basquin life N_f = (sigma_crit / sigma_a)^m  [cycles = days]."""
        out, Nf_min = {}, np.inf
        for mat, sc in ((INORGANIC, sigma_crit_inorg),
                        (ORGANIC, sigma_crit_org)):
            sig_a = mat.E_mod / (1 - mat.nu) * abs(mat.cte - self.substrate.cte) \
                    * dT_cycle / 2.0
            Nf = (sc / max(sig_a, 1.0)) ** m_fatigue
            out[mat.name] = dict(sigma_a_MPa=sig_a / 1e6, N_f_cycles=Nf)
            Nf_min = min(Nf_min, Nf)
        out["N_f_cycles"] = Nf_min
        out["years"] = Nf_min / 365.0
        out["cracked_as_deposited"] = \
            self.stack.d_inorg_nm > self.stack.defects.d_crit_nm
        return out

    # ---------------- one-call 6-objective evaluation ----------------
    def evaluate_design(self, damp_hours=600.0, diurnal_days=3.0,
                        env: Optional[Environment] = None,
                        M_crit_g_m2=0.01) -> Dict:
        run_dh = self.run_damp_heat(t_end_h=damp_hours)
        t80 = self.t80_hours(run_dh, M_crit_g_m2)
        run_di = self.run_diurnal(env or diurnal_cycle(), days=diurnal_days)
        last = run_di["t"] > (run_di["t"][-1] - DAY)             # final day
        T_max = float(run_di["T_dev"][last].max())
        dT = float(run_di["T_dev"][last].max() - run_di["T_dev"][last].min())
        dur = self.durability(dT)
        return dict(
            lifetime_T80_h=t80,                    # obj 1 (max)
            max_temp_C=T_max - 273.15,             # obj 2 (min)
            transmittance=self.transmittance(),    # obj 3 (max)
            cost=self.cost(),                      # obj 4 (min)
            weight_g_m2=self.added_weight_g_m2(),  # obj 5 (min)
            durability_years=dur["years"],         # obj 6 (max)
            wvtr_g_m2_day=run_dh["wvtr_ss"],
            dT_cycle_K=dT, durability=dur,
            runs=dict(damp=run_dh, diurnal=run_di))


# ------------------------------------------------------------------
# 7. Drop-in replacement for the user's original analytic API
# ------------------------------------------------------------------
def simulate_tfe_compat(d_org, d_inorg, n_pairs,
                        damp_hours=400.0, diurnal_days=2.0) -> Dict:
    twin = TFEDigitalTwinPDE(d_org, d_inorg, n_pairs)
    k = twin.evaluate_design(damp_hours, diurnal_days)
    return {"lifetime": k["lifetime_T80_h"], "max_temp": k["max_temp_C"],
            "transmittance": k["transmittance"], "cost": k["cost"],
            "weight": k["weight_g_m2"], "durability": k["durability_years"]}


# ------------------------------------------------------------------
# 8. Verification against the analytic Fickian slab solution
# ------------------------------------------------------------------
def validate_single_slab(L=1e-6, D=1e-14, S=10.0, a1=0.85,
                         N=60, t_end_fac=4.0, n_terms=300):
    """Numerical breakthrough flux vs the exact Fourier-series solution
        J(t)/J_ss = 1 + 2*sum_n (-1)^n exp(-n^2 pi^2 D t / L^2)."""
    dx = np.full(N, L / N)
    solver = ImplicitDiffusion1D(dx)
    K, C = np.full(N, D * S), np.full(N, S)
    t_char = L * L / D
    dt = t_char / 800.0
    a = np.zeros(N)
    ts, Jn = [], []
    t = 0.0
    while t < t_end_fac * t_char:
        a, Jl, _ = solver.step(a, dt, K, C,
                               ("dirichlet", 0.0), ("dirichlet", a1))
        t += dt
        ts.append(t); Jn.append(Jl)
    ts, Jn = np.array(ts), np.array(Jn)
    J_ss = D * S * a1 / L
    n = np.arange(1, n_terms + 1)[:, None]
    Ja = J_ss * (1 + 2 * np.sum((-1.0) ** n *
                                np.exp(-D * n ** 2 * np.pi ** 2 * ts[None, :] / L ** 2),
                                axis=0))
    mask = ts > 0.05 * t_char
    err = float(np.max(np.abs(Jn[mask] - Ja[mask]) / J_ss))
    return ts / t_char, Jn / J_ss, Ja / J_ss, err


# ------------------------------------------------------------------
# 9. Demo / figures
# ------------------------------------------------------------------
def _print_kpis(name: str, k: Dict):
    print(f"\n[{name}]  6-objective KPI summary")
    print(f"  1. T80 lifetime @85C/85%RH : {k['lifetime_T80_h']:10.0f} h"
          f"   (steady WVTR = {k['wvtr_g_m2_day']:.3e} g/m2/day)")
    print(f"  2. Max device temperature  : {k['max_temp_C']:10.1f} degC"
          f"   (diurnal dT = {k['dT_cycle_K']:.1f} K)")
    print(f"  3. Optical transmittance   : {k['transmittance']:10.3f}")
    print(f"  4. Cost (user model)       : {k['cost']:10.1f} arb")
    print(f"  5. Added weight            : {k['weight_g_m2']:10.2f} g/m2")
    print(f"  6. Fatigue durability      : {k['durability_years']:10.1f} yr"
          f"   (cracked as-deposited: "
          f"{k['durability']['cracked_as_deposited']})")


def main(make_figures=True, out_prefix="."):
    import os
    np.set_printoptions(precision=3)
    print("=" * 68)
    print(" TFE Digital Twin - Physics Engine v1.0  (PDE prototype demo)")
    print("=" * 68)

    # --- (A) numerical verification -------------------------------
    tv, Jn, Ja, err = validate_single_slab()
    print(f"\n(A) Verification vs analytic Fickian slab: "
          f"max flux error = {err*100:.2f} %  "
          f"({'PASS' if err < 0.02 else 'CHECK'})")

    # --- (B) baseline design: full transient behaviour -------------
    base = TFEDigitalTwinPDE(d_org_nm=200, d_inorg_nm=50, n_pairs=3)
    st = base.stack
    print(f"\n(B) Baseline stack: 3 x (Al2O3 50 nm / organic 200 nm), "
          f"L = {st.L*1e9:.0f} nm, {len(st.dx)} cells")
    print(f"    defect model: f_pin = {st.f_pin:.2e}, "
          f"pinhole spacing s = {st.s_lat*1e6:.1f} um, "
          f"tau^2(sandwiched organic) = {st.tau2_sand:.0f}")
    kpi = base.evaluate_design(damp_hours=600, diurnal_days=3)
    _print_kpis("baseline 3x(50/200) nm", kpi)
    run_dh = kpi["runs"]["damp"]
    run_di = kpi["runs"]["diurnal"]
    print(f"    damp-heat lag time (time-lag method): "
          f"{run_dh['t_lag_h']*60:.1f} min, "
          f"steady flatness J_end/J_ss = {run_dh['steady_flatness']:.3f}")

    # --- (C) geometry study: same 150 nm inorganic budget ----------
    print("\n(C) Geometry effect at FIXED inorganic budget (150 nm total),"
          "\n    d_org = 200 nm  ->  pure n_pairs / thickness trade-off:")
    designs = [("1 x 150 nm", 1, 150.0),
               ("3 x  50 nm", 3, 50.0),
               ("5 x  30 nm", 5, 30.0)]
    comp = []
    for label, n, d_in in designs:
        tw = TFEDigitalTwinPDE(200.0, d_in, n)
        r = tw.run_damp_heat(t_end_h=600)
        t80 = tw.t80_hours(r)
        comp.append((label, r["wvtr_ss"], t80,
                     tw.stack.defects.pinhole_fraction(d_in)))
        print(f"    {label}:  WVTR = {r['wvtr_ss']:.3e} g/m2/day,  "
              f"T80 = {t80:7.0f} h,  f_pin = {comp[-1][3]:.2e}")
    print("    -> non-monotonic optimum: pinhole closure favours thick "
          "layers,\n       cracking + defect decoupling favour multilayers "
          "(geometry novelty).")

    if not make_figures:
        return kpi, comp

    # --- figures ----------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    fig.suptitle("TFE Digital Twin v1 - moisture PDE engine "
                 "(damp heat 85C/85%RH)", fontsize=13)

    a0 = ax[0, 0]
    a0.plot(tv, Ja, "k-", lw=2, label="analytic (Fourier series)")
    a0.plot(tv[::12], Jn[::12], "o", ms=4, mfc="none", label="FVM (this work)")
    a0.set(xlabel=r"$t\,/\,(L^2/D)$", ylabel=r"$J(t)\,/\,J_{ss}$",
           title=f"(a) Verification, max err = {err*100:.2f}%")
    a0.legend(); a0.grid(alpha=0.3)

    a1 = ax[0, 1]
    cmap = plt.cm.viridis(np.linspace(0.1, 0.9, len(run_dh["snapshots"])))
    for (tsnap, prof), c in zip(run_dh["snapshots"], cmap):
        a1.plot(run_dh["x_nm"], prof, color=c,
                label=f"t = {tsnap/3600:.2g} h")
    for x0, x1 in st.inorg_spans_nm():
        a1.axvspan(x0, x1, color="grey", alpha=0.25)
    a1.set(xlabel="depth from device [nm]", ylabel="water activity a = C/S",
           title="(b) Activity profiles (grey = inorganic)")
    a1.legend(fontsize=8); a1.grid(alpha=0.3)

    a2 = ax[1, 0]
    m = run_dh["t"] > 0
    a2.loglog(run_dh["t"][m] / 3600, run_dh["wvtr"][m], "b-")
    a2.axhline(run_dh["wvtr_ss"], color="r", ls="--",
               label=f"steady = {run_dh['wvtr_ss']:.2e}")
    a2.set(xlabel="time [h]", ylabel="WVTR [g m$^{-2}$ day$^{-1}$]",
           title="(c) Breakthrough curve (baseline)")
    a2.legend(); a2.grid(alpha=0.3, which="both")

    a3 = ax[1, 1]
    labels = [c[0] for c in comp]
    xpos = np.arange(len(comp))
    b1 = a3.bar(xpos - 0.18, [c[1] for c in comp], width=0.36,
                color="steelblue", label="WVTR")
    a3.set_yscale("log")
    a3.set_ylabel("WVTR [g m$^{-2}$ day$^{-1}$]", color="steelblue")
    a3b = a3.twinx()
    b2 = a3b.bar(xpos + 0.18, [c[2] for c in comp], width=0.36,
                 color="darkorange", label="T80")
    a3b.set_ylabel("T80 [h]", color="darkorange")
    a3.set_xticks(xpos); a3.set_xticklabels(labels)
    a3.set_title("(d) Fixed 150 nm inorganic budget")
    a3.grid(alpha=0.3, axis="y")
    f1 = os.path.join(out_prefix, "fig1_damp_heat_engine.png")
    fig.savefig(f1, dpi=150)

    fig2, bx = plt.subplots(3, 1, figsize=(10, 9), sharex=True,
                            constrained_layout=True)
    fig2.suptitle("Coupled heat <-> moisture under diurnal cycling "
                  "(baseline design)", fontsize=13)
    th = run_di["t"] / 3600.0
    bx[0].plot(th, run_di["T_amb"] - 273.15, "k--", label="T ambient")
    bx0b = bx[0].twinx()
    bx0b.fill_between(th, run_di["G"], color="gold", alpha=0.4)
    bx0b.set_ylabel("solar flux [W/m$^2$]", color="goldenrod")
    bx[0].plot(th, run_di["RH"] * 100 / 2, "c:",
               label="RH/2 [%] (right scale ~)")
    bx[0].set_ylabel("T [degC]"); bx[0].legend(loc="upper left")
    bx[0].set_title("(a) Environmental forcing")

    bx[1].plot(th, run_di["T_dev"] - 273.15, "r-", label="device")
    bx[1].plot(th, run_di["T_tfe"] - 273.15, "b:", label="TFE mean")
    bx[1].axhline(kpi["max_temp_C"], color="r", ls="--", alpha=0.5)
    bx[1].set_ylabel("T [degC]"); bx[1].legend()
    bx[1].set_title(f"(b) Device temperature "
                    f"(max = {kpi['max_temp_C']:.1f} degC, "
                    f"dT cycle = {kpi['dT_cycle_K']:.1f} K)")

    bx[2].plot(th, run_di["wvtr_uptake"], color="steelblue", alpha=0.8,
               label="uptake flux @ ambient face")
    bx[2].plot(th, run_di["wvtr"], "r-", lw=2,
               label="ingress flux @ device face")
    bx2b = bx[2].twinx()
    bx2b.plot(th, run_di["M"] * 1e3, "g--", label="cumulative")
    bx2b.set_ylabel("cumulative ingress [mg/m$^2$]", color="g")
    bx[2].set_ylabel("flux [g m$^{-2}$ day$^{-1}$]")
    bx[2].set_xlabel("time [h]")
    bx[2].legend(loc="upper left", fontsize=8)
    bx[2].set_title("(c) Arrhenius D(T) + RH(t) forcing: the stack acts as a "
                    "low-pass filter (device-side flux is damped)")
    for b in bx:
        b.grid(alpha=0.3)
    f2 = os.path.join(out_prefix, "fig2_diurnal_coupled.png")
    fig2.savefig(f2, dpi=150)
    print(f"\nFigures saved: {f1}, {f2}")
    return kpi, comp


if __name__ == "__main__":
    main()
