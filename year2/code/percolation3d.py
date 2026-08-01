# -*- coding: utf-8 -*-
"""
Year-2 module: three-dimensional defect-network model of multilayer barriers.

Purpose
-------
The 1-D twin (Year-1) represents the organic interlayer by an analytic detour
factor

    tau^2 = 1 + s^2 ln(s/r) / (2 pi d_org^2),      s = mean pinhole spacing,

which assumes (i) pinholes in successive inorganic layers are statistically
independent and (ii) each pinhole feeds its neighbour through a single lateral
path.  This module drops both assumptions and solves the actual transport
problem on a three-dimensional network, so that the analytic factor can be
tested rather than trusted.

Physical picture and scale separation
-------------------------------------
Three lengths are widely separated in every real stack:

    pinhole radius r (~50 nm)  <<  organic thickness d_org (~1e2 nm)
                               <<  pinhole spacing s (~1e5 nm)

so the organic interlayer transports laterally like a thin conducting sheet
(sheet conductance P_org * d_org per square), while each pinhole is a point
terminal connecting two adjacent sheets.  We therefore discretise every organic
interlayer as a periodic 2-D grid and connect the grids through pinholes, with
the sub-grid constriction resistance of each pinhole added analytically.  This
is exact in the limit r << Delta << s and avoids a full 3-D mesh, which would
need ~1e8 cells for a single realisation.

Network elements (conductances, SI: kg m^-2 s^-1 per unit activity)
-------------------------------------------------------------------
sheet cell-to-cell     G_sheet = P_org * d_org
through-pinhole        G_hole  = 1 / (R_through + n_c * R_constrict)
    R_through   = d_in / (P_org * pi r^2)          (hole filled with organic)
    R_constrict = ln(r_cell / r) / (2 pi P_org d_org),  r_cell = Delta / sqrt(pi)
    n_c         = number of sheet-side constrictions (1 at the boundaries, 2 inside)
bulk inorganic         G_bulk  = P_lat * Delta^2 / d_in     (per cell, optional)

Defect correlation
------------------
A fraction phi of defects is *columnar*: a particle that punctures every
inorganic layer at the same (x, y).  The remaining fraction is independent
between layers.  phi = 0 reproduces the Year-1 assumption; phi > 0 introduces
short circuits that no additional dyad can close, and is the candidate
explanation for the WVTR plateau seen in published dyad ladders.

Author: Teo Yang.  Released with the project repository.
"""

from dataclasses import dataclass, field
import numpy as np
from scipy.sparse import coo_matrix, identity
from scipy.sparse.linalg import factorized


# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------
@dataclass
class Stack:
    """Geometry and material parameters of one multilayer stack."""
    n_inorg: int = 3           # number of inorganic (pinhole-bearing) layers
    d_in: float = 50e-9        # inorganic thickness [m]
    d_org: float = 500e-9      # organic interlayer thickness [m]
    f: float = 3.72e-8         # pinhole area fraction (Year-1 calibration)
    r_pin: float = 50e-9       # pinhole radius [m]
    P_org: float = 7.5e-13     # organic permeability D*S [kg m^-1 s^-1 per activity]
    P_lat: float = 0.0         # intact-oxide permeability (0 = pinhole-dominated)

    @property
    def density(self) -> float:
        """Pinhole number density [m^-2]."""
        return self.f / (np.pi * self.r_pin ** 2)

    @property
    def spacing(self) -> float:
        """Mean pinhole spacing s = 1/sqrt(density) [m]."""
        return 1.0 / np.sqrt(self.density)

    def tau2(self) -> float:
        """Year-1 analytic detour factor of the organic interlayer."""
        s, r = self.spacing, self.r_pin
        return 1.0 + s ** 2 * np.log(s / r) / (2 * np.pi * self.d_org ** 2)

    def conductance_1d(self, area: float) -> float:
        """Year-1 prediction for the same stack and area [kg s^-1 per activity]."""
        P_in = self.P_lat + self.f * self.P_org
        R = self.n_inorg * self.d_in / (P_in * area)
        R += (self.n_inorg - 1) * self.d_org * self.tau2() / (self.P_org * area)
        return 1.0 / R


# --------------------------------------------------------------------------
# defect placement
# --------------------------------------------------------------------------
def place_defects(stack: Stack, L: float, rng, phi: float = 0.0,
                  arrangement: str = "poisson"):
    """Return a list of (N_k, 2) position arrays, one per inorganic layer.

    phi          fraction of defects that are columnar (identical in every layer)
    arrangement  'poisson' (random) or 'lattice' (regular square array)
    """
    n_total = int(round(stack.density * L * L))
    if n_total < 1:
        raise ValueError("domain too small: no defects. Increase L.")
    n_col = int(round(phi * n_total))
    n_ind = n_total - n_col

    def draw(n):
        if n <= 0:
            return np.empty((0, 2))
        if arrangement == "lattice":
            m = int(np.ceil(np.sqrt(n)))
            xs = (np.arange(m) + 0.5) * L / m
            X, Y = np.meshgrid(xs, xs)
            return np.c_[X.ravel(), Y.ravel()][:n]
        return rng.random((n, 2)) * L

    columnar = draw(n_col)                       # same positions in every layer
    layers = [draw(n_ind) for _ in range(stack.n_inorg)]
    return layers, columnar


# --------------------------------------------------------------------------
# network assembly and solution
# --------------------------------------------------------------------------
def effective_conductance(stack: Stack, L: float, nx: int, layers,
                          columnar=None, return_detail: bool = False):
    """Solve the sheet network and return the effective conductance [kg s^-1].

    Boundary conditions: activity 1 above the top inorganic layer, 0 below the
    bottom one; periodic in x and y.  Nodes are the cells of the (n_inorg - 1)
    interior organic sheets.
    """
    columnar = np.empty((0, 2)) if columnar is None else np.asarray(columnar)
    n_col = len(columnar)
    n_sheets = stack.n_inorg - 1
    if n_sheets < 1:
        raise ValueError("need at least two inorganic layers for a sheet network")
    delta = L / nx
    ncell = nx * nx
    N = n_sheets * ncell + n_col * n_sheets     # sheet cells + channel nodes

    G_sheet = stack.P_org * stack.d_org                      # per cell face
    r_cell = delta / np.sqrt(np.pi)
    R_through = stack.d_in / (stack.P_org * np.pi * stack.r_pin ** 2)
    R_constr = np.log(max(r_cell / stack.r_pin, 1.001)) / (2 * np.pi * stack.P_org * stack.d_org)
    G_bulk_cell = stack.P_lat * delta ** 2 / stack.d_in if stack.P_lat > 0 else 0.0

    rows, cols, vals = [], [], []
    diag = np.zeros(N)
    rhs = np.zeros(N)

    def cell_index(sheet, ix, iy):
        return sheet * ncell + (iy % nx) * nx + (ix % nx)

    # ---- in-plane sheet conduction (periodic) ----
    idx = np.arange(ncell)
    ix, iy = idx % nx, idx // nx
    for dx, dy in ((1, 0), (0, 1)):
        jx, jy = (ix + dx) % nx, (iy + dy) % nx
        j = jy * nx + jx
        for sheet in range(n_sheets):
            a = sheet * ncell + idx
            b = sheet * ncell + j
            rows.extend(a); cols.extend(b); vals.extend([-G_sheet] * ncell)
            rows.extend(b); cols.extend(a); vals.extend([-G_sheet] * ncell)
            np.add.at(diag, a, G_sheet)
            np.add.at(diag, b, G_sheet)

    # ---- vertical links through inorganic layers ----
    #  layer k connects sheet k-1 to sheet k  (sheet -1 = source, sheet n_sheets = sink)
    for k, pts in enumerate(layers):
        upper, lower = k - 1, k                     # sheet indices
        at_top, at_bot = (upper < 0), (lower >= n_sheets)
        n_c = 2 - int(at_top) - int(at_bot)         # constrictions on sheet sides
        G_hole = 1.0 / (R_through + n_c * R_constr)

        cx = np.clip((pts[:, 0] / delta).astype(int), 0, nx - 1)
        cy = np.clip((pts[:, 1] / delta).astype(int), 0, nx - 1)
        cid = cy * nx + cx
        # several defects may fall in one cell: conductances add
        uniq, counts = np.unique(cid, return_counts=True)
        g = G_hole * counts

        if at_top:                                  # source -> sheet 0
            a = 0 * ncell + uniq
            np.add.at(diag, a, g)
            np.add.at(rhs, a, g * 1.0)              # activity 1 at source
        elif at_bot:                                # sheet n-1 -> sink (activity 0)
            a = (n_sheets - 1) * ncell + uniq
            np.add.at(diag, a, g)
        else:                                       # sheet k-1 <-> sheet k
            a = upper * ncell + uniq
            b = lower * ncell + uniq
            rows.extend(a); cols.extend(b); vals.extend(-g)
            rows.extend(b); cols.extend(a); vals.extend(-g)
            np.add.at(diag, a, g)
            np.add.at(diag, b, g)

        if G_bulk_cell > 0:                         # intact-oxide leakage, all cells
            if at_top:
                a = np.arange(ncell)
                np.add.at(diag, a, G_bulk_cell); np.add.at(rhs, a, G_bulk_cell)
            elif at_bot:
                a = (n_sheets - 1) * ncell + np.arange(ncell)
                np.add.at(diag, a, G_bulk_cell)
            else:
                a = upper * ncell + np.arange(ncell)
                b = lower * ncell + np.arange(ncell)
                rows.extend(a); cols.extend(b); vals.extend([-G_bulk_cell] * ncell)
                rows.extend(b); cols.extend(a); vals.extend([-G_bulk_cell] * ncell)
                np.add.at(diag, a, G_bulk_cell); np.add.at(diag, b, G_bulk_cell)

    # ---- columnar (particle) channels: continuous tubes through every layer ----
    #  A particle punctures all inorganic layers at one (x, y), so water crosses
    #  the stack without ever spreading laterally.  Each channel gets its own
    #  node at every sheet level: vertical links carry R_through only, while a
    #  constriction link couples the channel to the surrounding sheet cell.
    if n_col:
        base = n_sheets * ncell
        cx = np.clip((columnar[:, 0] / delta).astype(int), 0, nx - 1)
        cy = np.clip((columnar[:, 1] / delta).astype(int), 0, nx - 1)
        cid = cy * nx + cx
        G_v = 1.0 / R_through           # hole-to-hole, no spreading
        G_c = 1.0 / R_constr            # channel <-> local sheet
        for lvl in range(n_sheets):
            node = base + lvl * n_col + np.arange(n_col)
            cell = lvl * ncell + cid
            rows.extend(node); cols.extend(cell); vals.extend([-G_c] * n_col)
            rows.extend(cell); cols.extend(node); vals.extend([-G_c] * n_col)
            np.add.at(diag, node, G_c); np.add.at(diag, cell, G_c)
            if lvl == 0:                                   # source -> channel top
                np.add.at(diag, node, G_v); np.add.at(rhs, node, G_v * 1.0)
            else:                                          # channel segment
                prev = base + (lvl - 1) * n_col + np.arange(n_col)
                rows.extend(node); cols.extend(prev); vals.extend([-G_v] * n_col)
                rows.extend(prev); cols.extend(node); vals.extend([-G_v] * n_col)
                np.add.at(diag, node, G_v); np.add.at(diag, prev, G_v)
            if lvl == n_sheets - 1:                        # channel -> sink
                np.add.at(diag, node, G_v)

    rows.extend(range(N)); cols.extend(range(N)); vals.extend(diag)
    A = coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsc()
    A = A + 1e-30 * identity(N, format="csc")       # guard against isolated cells
    solve = factorized(A)
    a_field = solve(rhs)

    # current into the sink = sum over bottom-layer defects of G * activity
    pts = layers[-1]
    cx = np.clip((pts[:, 0] / delta).astype(int), 0, nx - 1)
    cy = np.clip((pts[:, 1] / delta).astype(int), 0, nx - 1)
    uniq, counts = np.unique(cy * nx + cx, return_counts=True)
    G_hole_bot = 1.0 / (R_through + 1 * R_constr)
    I = float(np.sum(G_hole_bot * counts * a_field[(n_sheets - 1) * ncell + uniq]))
    if G_bulk_cell > 0:
        I += float(G_bulk_cell * np.sum(a_field[(n_sheets - 1) * ncell:n_sheets * ncell]))
    if n_col:
        tail = n_sheets * ncell + (n_sheets - 1) * n_col
        I += float((1.0 / R_through) * np.sum(a_field[tail:tail + n_col]))

    if return_detail:
        return I, a_field[:n_sheets * ncell].reshape(n_sheets, nx, nx)
    return I                                        # activity drop is unity


def wvtr(G: float, area: float, da: float = 0.90) -> float:
    """Convert conductance to WVTR in g m^-2 day^-1 for a given driving activity."""
    return G * da / area * 1e3 * 86400.0
