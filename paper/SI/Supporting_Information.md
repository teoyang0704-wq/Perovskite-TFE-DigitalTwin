# Supporting Information

**for** "Geometry-first design of organic/inorganic multilayer thin-film encapsulation via a literature-calibrated digital twin" (draft v2)

**Correction note (2026-07-10).** Full-text verification of Ref. [7] established that the Fig.-3 dyad series uses **50 nm** Al2O3 ('primal thickness'), not 30 nm; the Stage-B anchor and every downstream quantity were regenerated. Net effects: f_res −6.7% (design space +7% lifetimes), validation errors improve to +0.11/−0.01 decade, thresholds and all structural conclusions unchanged. The article states its detection limit but **no run duration**; all duration statements now use the insensitivity window (S4.3).

Contents: S1 Governing equations, numerics and time-lag derivations · S2 Steady-state ladder bound · S3 Defect model: derivation, validity, identifiability · S4 Calibration protocol and uncertainty quantification · S5 Extended tables · S6 Reproducibility. Figures S1–S3; Tables S1–S4.

---

## S1. Governing equations, numerics, and the layer-lumped time lag

### S1.1 Activity formulation
With Henry-law sorption C = S(x)·a (a: water activity), Fick's second law ∂C/∂t = ∂x(D ∂x C) becomes S ∂a/∂t = ∂x(P ∂x a), P ≡ D·S. Continuity of chemical potential across an interface implies continuity of *a* (not of C), so the activity field is smooth through the stack and equals RH at the ambient face; this is the standard multilayer interface condition and removes solubility discontinuities from the numerics. The Henry assumption is revisited in the main-text limitations (water clustering near saturation).

### S1.2 Finite-volume discretization and verification
Cell-centered volumes with interface conductances G_{i+1/2} = [Δx_i/(2P_i) + Δx_{i+1}/(2P_{i+1})]⁻¹ (exact flux matching for piecewise-constant P; the discrete series resistance Σ 1/G equals Σ Δx/P identically), fully implicit backward Euler, tridiagonal solve. Unconditional stability accommodates the >10⁵ permeability contrast between layers. Verification against the exact Fourier-series breakthrough flux of a uniform slab, J(t)/J_ss = 1 + 2Σₙ(−1)ⁿ exp(−n²π²Dt/L²): maximum relative flux error **0.45%** for t > 0.05·L²/D (main Fig. 2/3 pipeline; script `src/tfe_physics_engine.py`).

### S1.3 Exact layer-lumped Frisch time lag
For 1-D diffusion with position-dependent resistance density ρ(x) = 1/P(x) and capacity c(x) = S(x), the flux time lag at the outlet after a step at the inlet is (Frisch [22])
t_lag = (1/R_tot) ∫₀ᴸ c(x) R_L(x) R_R(x) dx, with R_L(x) = ∫₀ˣ ρ dx′, R_R = R_tot − R_L.
For a stack of uniform layers (resistance r_i, capacity c_i, entered at cumulative resistance a_i, with b_i = R_tot − a_i), the integral evaluates in closed form per layer:
Δt_i = (c_i/r_i)[a_i b_i r_i + (b_i − a_i) r_i²/2 − r_i³/3],  t_lag = Σᵢ Δt_i / R_tot.
Uniform-slab check: a=0, b=R=L/P, r=R, c=SL ⇒ t_lag = SL·(R²/2 − R²/3)/R·(1/…) = S L²/(6P) = **L²/6D** ✓. This closed form makes the lag exact without time stepping and is what the fast evaluator uses.

## S2. Steady-state ladder bound (why geometric per-dyad scaling is impossible in-model)

For n dyads the ladder is R_n ≈ R_1 + (n−1)(R_in + R_sand) (substrate and top-organic terms are negligible for the validation stack: R_PET = 2.9×10⁷, R_top = 3.4×10⁵ ≪ R_in = 5.1×10¹¹ s·m²·kg⁻¹). Hence
R₃/R₁ = 1 + 2(1 + ρ), ρ ≡ R_sand/R_in = (d_org/d_in)·[f + r² ln(s/r)/(2 d_org²)],
using R_in = d_in/(f P_org), R_sand = d_org τ²/P_org and τ² − 1 = s² ln(s/r)/(2π d_org²) with s² = π r²/f (so f cancels except in the residual f-term and the logarithm). For the validation geometry (d_org = 500 nm, d_in = 50 nm) and admissible parameters (r ≤ 100 nm from the ×/÷2 prior; ln(s/r) ≤ 12 for any s ≤ 1 cm; f ≤ 10⁻⁶ negligible):
ρ ≤ (500/50)·[10⁻⁶ + (100 nm)²·12/(2·(500 nm)²)] = 10 × 0.24 ≈ **2.4**, ⇒ R₃/R₁ ≤ 1 + 2(1+2.4) = **7.8**.
The measured geometric improvement requires R₃/R₁ ≈ 22 (÷4.7 per dyad). Therefore no admissible parameter set reproduces the dyad series at steady state; the discrepancy is structural, resolved only by the finite-duration (lag-phase) simulation of §3.1.

## S3. Defect model: derivation, validity, identifiability

### S3.1 Pinhole area fraction f(d)
f(d) = f₀e^(−d/d_c) + f_res + f_c0[(d−d_crit)/d_crit]₊². The exponential encodes nucleation-driven coalescence (closure); the floor f_res encodes particle-seeded defects that do not anneal with thickness—identified as the dominant defect source in Ref. [7] and consistent with defect-permeation correlations on coated polymers [17,18]; the quadratic onset encodes channel cracking, whose energy-release rate grows with film thickness (strained flexible-substrate regime of Ref. [7]).

### S3.2 Lateral (tortuous-path) reduction τ²
Between vertically offset pinholes of areal density 1/s² (s = r√(π/f)), permeant must traverse the organic interlayer laterally over ≈s through a conduit of thickness d_org. The per-unit-area lateral resistance of this stage is the disk-source spreading result R_lat ≈ s² ln(s/r)/(2π P_org d_org); demanding an equivalent 1-D layer of the same geometric thickness gives P_org^eff = P_org/τ² with **τ² = 1 + s² ln(s/r)/(2π d_org²)** — the 1-D reduction of the lag-time/tortuous-path picture of Graff et al. [8]. Validity requires s ≫ d_org and s ≫ r (lateral flow is thin-film-like): here s = 343–459 µm against d_org ≤ 1 µm and r ≈ 50 nm, satisfied by >2 orders of magnitude across the whole design space.

### S3.3 Identifiability (Fig. S1)
Profile analysis over fixed d_c (refitting all other parameters) shows the sum of squared residuals rising steeply away from the lower bound: SSR = 3.3×10⁻³ decade² at d_c = 1 nm versus 0.27 at 6 nm. Against the measurement-noise level Nσ² = 5·(0.08)² = 0.032, the admissible region is **d_c ≤ 1.8 nm**: the closure length is *upper-bounded, not identified*—closure completes below the first measured thickness (15 nm)—with f₀ compensating along the boundary (Fig. S1a). The identified invariant is the ratio **f(15 nm)/f_res = 9.1** (Fig. S1b), together with d_crit = 43.9 nm, both stable across the admissible region. Consequence adopted throughout: the design space is bounded at d_in ≥ 15 nm and no sub-15 nm statement is made.

## S4. Calibration protocol and uncertainty quantification

### S4.1 Two-stage calibration
Stage A fits {K, f₀/f_res, d_c, f_c0/f_res, d_crit} to the five-point single-layer series in log space (max residual 0.047 decade), under WVTR_bare(d) ∝ f(d). Consistency of that assumption: the implied per-area bare-pinhole conductance is **5.8×10⁻⁶** of the open-vapor limit 4P_vap/(πr), i.e., strongly access/spreading-limited [17,18], hence thickness-independent per hole. Stage B solves f_res from the 1-dyad point: f(50) = [d_in/(R₁ − R_PET − R_top) − P_lat]/P_org (d_in = 50 nm), giving f(50 nm) = 6.67×10⁻⁸ → 8.5 pinholes mm⁻², spacing 343 µm (anchor: the 50 nm dyad series).

### S4.2 Priors for Monte-Carlo (Table S4) and convergence
Each uncertain input, its distribution and rationale are listed in Table S4. The full Stage-A/B calibration is re-run per draw (N = 500, seed 7), preserving calibration-induced correlations. Convergence (Fig. S2): |median(N=500) − median(N=250)| = 0.0 nm (closure), 0.0 nm (crack onset), 0.003 (relative pair gain) — fully converged at N = 500.


### S4.4b Device-anchored effective floor for M_crit (added after v2 review)
Using the encapsulated-device shelf data of Ref. [4] (<4% PCE loss after 7,500 h at 25 degC/50% RH) together with the same film's measured barrier WVTR (1.84e-2 g m-2 day-1 at 45 degC/100% RH), transferred to shelf conditions by the vapor-access mechanism of S4.1 (factor 0.156; Arrhenius-40 kJ/mol alternative within 15%), the sink model accumulates a dose of 0.90-1.04 g m-2 over the test. Since the device remained within 4% of its initial PCE, the *same-model* effective threshold obeys **M_crit^eff >= 0.97 g m-2 (Tier-1: monotone dose-damage only)**; a linear dose-damage reading gives >= 4.9 g m-2 (Tier-2). Assumptions: (a) effective calibration within the identical sink model, transferable if the sink bias is non-increasing with barrier resistance; (b) T/RH transfer by the vapor mechanism; (c) Tier-1 requires only monotonicity; (d) the film's own measured WVTR is used, so film-specific f cancels; (e) conformal coating suppresses edge ingress. Consequence: the dose term of every lifetime scales by >=97x; within the G1-G2 window even a single dyad exceeds ~25 yr of moisture-limited life, so moisture ceases to be the binding lifetime channel for dyad designs, while thresholds and rankings are unchanged (main G7).

### S4.3 Test-duration sensitivity (Fig. S3)
Sweeping the assumed instrument run time: prediction errors (decades) for the {2,3}-dyad points are {−0.37, −1.34} at 1 d, {−0.03, −0.39} at 2 d, {+0.11, −0.01} at 3 d, {+0.22, +0.35} at 5 d, {+0.26, +0.51} at 7 d. Both points remain within ≤0.4 decade for assumed durations of **2.0–5.5 days**, spanning typical multi-day coulometric runs (the source states the 5×10⁻⁵ detection limit but not the duration); steady state (t→∞) fails per §S2 regardless.

## S5. Extended tables

**Table S1 — Exact Sobol decomposition (uniform grid measure), all six objectives** (file `TableS1_sobol_full.csv`):

| objective | Var(Y) | d_org | d_in | n | d_org×d_in | d_org×n | d_in×n | 3-way | ST_d_org | ST_d_in | ST_n |
|---|---|---|---|---|---|---|---|---|---|---|---|
| t₈₀ᵈᵒˢᵉ [h] | 3.48×10⁷ | 0.033 | 0.551 | 0.141 | 0.052 | 0.015 | 0.184 | 0.024 | 0.125 | 0.811 | 0.364 |
| durability [yr] | 3.43×10² | 0 | 1.000 | 0 | 0 | 0 | 0 | 0 | 0 | 1.000 | 0 |
| cost [min] | 3.34×10⁶ | 0.008 | 0.371 | 0.532 | 0.000 | 0.002 | 0.088 | ~0 | 0.009 | 0.459 | 0.622 |
| weight [g/m²] | 4.45 | 0.343 | 0.024 | 0.545 | ~0 | 0.082 | 0.006 | 0 | 0.425 | 0.030 | 0.633 |
| T_max [°C] | 4.8×10⁻⁴ | 0.334 | 0.017 | 0.565 | 0 | 0.079 | 0.004 | 0 | 0.414 | 0.021 | 0.649 |
| T_opt [–] | 5.7×10⁻⁷ | 0.334 | 0.017 | 0.565 | 0 | 0.079 | 0.004 | 0 | 0.414 | 0.021 | 0.649 |

(For T_max and T_opt the *total variance itself* is negligible—0.10 K and 0.35 %p spans—so their share rows describe a vanishing pie; see main §3.2.)

**Table S2 — Design rules with 95% CI and confidence grades** (English rendering of `guidelines_table.csv`):
G1 [HIGH] Inorganic thickness ≥ 22.5 nm (CI 21.0–26.0): pinhole-closure completion. · G2 [HIGH] < 44.0 nm (37.7–48.0): cracking onset; beyond it lifetime and durability collapse jointly (dominated region; strained flexible-substrate regime). · G3 [LOW-MED] Within the window, N_f ∝ d^(−m/2) (median exponent −3.6 [−4.9, −0.0]); a 10-yr fatigue target needs d_in ≤ 30 nm (16–49) and is unreachable at any thickness in 16% of draws under the fatigue prior (proxy weakness, primary-data target). · G4 [HIGH structure / MED absolute] Pair law linear in all 500 draws (min R² = 0.9998); relative gain 0.5 (0.4–0.5) of the 3-dyad lifetime per added dyad; absolute gain 7,592 h (1,321–45,276) is M_crit-prior-dominated. · G5 [MED, transport-only, conditional] Do not exceed the planarization-required organic thickness (R_sand ∝ 1/d_org for n ≥ 2; planarization function unmodeled [21]). · G6 [HIGH] Thermal/optical objectives are not geometric degrees of freedom in this material pair (0.10 K / 0.35 %p). · G7 [HIGH] Thresholds and rankings invariant under M_crit ×/÷3.

**Table S3 — Representative Pareto designs (baseline parameters)**
MAX-t₈₀: 6×(42 nm/100 nm), 46,857 h, 1.7 yr, 2,671 min, 1.53 g/m². · KNEE: 6×(18/100), 25,399 h, 50.2 yr, 1,202 min, 1.10 g/m². · ECONOMY: 3×(21/100), 15,983 h, 27.1 yr, 693 min, 0.58 g/m². (Absolute lifetimes carry the M_crit credible span, ×37 at KNEE: 24,106 h [3,803–139,482] under UQ.)

**Table S4 — UQ priors and rationale (all sampled per-draw before re-calibration)**
| input | prior | rationale |
|---|---|---|
| Table-1 WVTR (5 pts) | lognormal σ = 0.08 dec | 2-sig-fig reporting + coulometric repeatability |
| 1-dyad anchor | lognormal σ = 0.10 dec | digitization ±0.15 dec (~1.5σ) |
| D_parylene | ×/÷3 log-uniform | quality-C datasheet split; validated-insensitive (±0.03 dec at ×3) |
| S_parylene | 0.5–2× log-uniform | datasheet uptake <0.1 wt% ∧ joint dyad-fit constraint (§3.1) |
| Ea (organic) | 30–50 kJ/mol uniform | typical polymer range; near-inert at the calibration temperature |
| pinhole radius r | ×/÷2 log-uniform | Ca-test morphology spread |
| M_crit | ×/÷3 log-uniform | unmeasured perovskite dose threshold (declared transfer item) |
| σ_c0 | ×/÷1.5 log-uniform | thin-film channel-cracking strength spread |
| fatigue exponent m | 6–10 uniform | Basquin range for brittle films (proxy) |

**NSGA-II settings**: population 80, generations 80, SBX crossover (η = 15, p_c = 0.9), polynomial mutation (η = 20, p_m = 1/3), binary tournament on (rank, crowding), real-coded genes on [0,1] with n rounded at evaluation; seed 42. Benchmarks vs exhaustive front: GD mean 0.060 / max 0.246; IGD 0.168 (80-point archive vs 3,064-point six-objective front).

**Exact Shapley**: with features F = {d_org, d_in, n}, φ_i = (1/3)[v(i)−v(∅)] + (1/6)Σ_{j≠i}[v(ij)−v(j)] + (1/3)[v(F)−v(F\i)], where every coalition value v(S) is the factorial-grid group mean conditioned on x_S (interventional semantics, uniform measure). Efficiency identity Σφ = y − ȳ satisfied to ≤10⁻¹¹ of the response for every objective.

## S6. Reproducibility

Live interactive deployment of the calibrated twin: https://tfe-twin.streamlit.app

Environment: Python 3.12.3, numpy 2.4.4, scipy 1.17.1, pandas 3.0.2, matplotlib 3.10.8. Seeds: NSGA-II 42; UQ 7. Runtimes (single CPU): exhaustive grid (4,104 designs) 2.0 s; NSGA-II 20 s; calibration + transient validation suite ≈2 s; UQ N = 500 (with per-draw re-calibration) 14.1 s; full figure regeneration ≈2 min. Commands (repository root): `python src/tfe_physics_engine.py`; `python src/step5_calibrate.py`; `python optimization/step6_optimize.py`; `python explainability/step7_xai.py`; `python uq/step8_uq.py`; `python figures/make_paper_figures.py`. Data provenance: every database row carries source ID, measurement conditions, censoring flag (=, <, ~), quality grade (A/B/C) and calibration/validation role; digitization overlays and the digitizer script are archived under `data/provenance/`. Bibliography entries flagged [VERIFY] must be resolved against primary sources before submission; no such entry is used for any quantitative claim except where explicitly noted.

## S7. Ablation study: what each mechanism is pinned by

Each physics term was removed and the **remaining** model fairly re-calibrated (same two-stage protocol), then scored on (i) the five-point single-layer fit and (ii) the zero-refit dyad validation at the confirmed ~72 h duration (Fig. S4a; `ablation/ablation_table.csv`).

| ablation | Table-1 max residual [dec] | dyad-2 err | dyad-3 err |
|---|---|---|---|
| baseline | 0.047 | +0.11 | −0.01 |
| no closure (f₀=0) | **0.654** | +0.11 | −0.01 |
| no particle floor (f_res=0) | **0.433** | +0.11 | −0.01 |
| no cracking (f_c0=0) | **0.538** | +0.11 | −0.01 |
| no tortuous path (τ²=1) | 0.047 | **+0.24** | **+0.33** |
| no lag (steady only) | 0.047 | +0.29 | **+0.75** |
| no calibration (generic 1 mm⁻², no fit) | — | **−0.64** | −0.18* |

The structure is diagnostic: the three **shape terms are each pinned by a distinct single-layer feature** (15 nm point → closure; 20–30 nm plateau → particle floor; 50–60 nm rise → cracking) and, when removed, break the calibration by ≥0.43 decade while leaving the d=50 nm dyad predictions untouched — i.e., they carry the *design-space geometry sensitivity* (without them the 22.5/44 nm thresholds do not exist). Conversely the **transport terms are pinned by the multilayer transient**: removing τ² degrades the joint 72 h agreement (max|err| 0.33 vs 0.11) and worsens the steady-state limit to +0.87 decade at 3 dyads; removing the lag reproduces the structural +0.75-decade failure of §S2. (*) The uncalibrated model's −0.18 at 3 dyads remains a partial error cancellation — it misses the 1-dyad anchor by −0.92 decade and scales incorrectly with n — and must not be read as calibration being dispensable.

## S8. Additional sensitivity views (response surface, interactions)

Fig. S4b maps the full response surface T80(d_in, n) at d_org = 100 nm from the exhaustive grid, visualizing the G1–G2 window as a ridge; Fig. S4c overlays T80(d_in) for n = 1…6, showing the multiplicative amplification that exact ANOVA quantifies as the d_in×n interaction (S = 0.184, Table S1). These views are projections of the same exhaustive dataset (no additional sampling machinery required).
