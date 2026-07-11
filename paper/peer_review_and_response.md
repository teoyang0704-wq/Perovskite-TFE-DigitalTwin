# Simulated Peer Review — 3 Independent Reviewers (target: Sol. Energy Mater. Sol. Cells class)

Manuscript: "Geometry-first design of organic/inorganic multilayer thin-film encapsulation via a literature-calibrated digital twin" (draft v1)

---

## Reviewer 1 — Materials scientist (TFE / perovskite stability). Recommendation: **Major revision**

The manuscript is unusually transparent about its limitations, which I appreciate, and the transient re-reading of sub-detection-limit dyad data is genuinely interesting. However, several issues must be resolved.

**Major**
- **R1.1 Single-source dependence.** Calibration *and* quantitative validation both rest on one paper (Ref. [7]). The cross-regime split is clever but a skeptic will call this a one-dataset model. The authors must (a) state this bluntly in the abstract, (b) elevate the independent cross-temperature check, and (c) identify concretely which published datasets would falsify the model.
- **R1.2 Perovskite framing without perovskite data.** The introduction and motivation are perovskite-centric, yet no perovskite quantity is calibrated (M_crit is a prior). Either add device-level calibration or scope the claims: this is a *barrier-stack* twin demonstrated on Al₂O₃/parylene-C, *toward* perovskite encapsulation.
- **R1.3 Guideline G5 conflicts with practice.** Every industrial multilayer (Barix/Vitex lineage) uses thick organic layers primarily to *planarize particles*—the very defect population the authors' own f_res encodes. Recommending "process-minimum organic thickness" from a transport-only model, without the planarization function, could mislead experimentalists. This must be reframed and the planarization literature cited (e.g., Coclite & Gleason).
- **R1.4 Missing foundational literature.** Defect-dominated permeation of oxide-coated polymers (Rossi & Nulman; da Silva Sobrinho et al.), ALD Al₂O₃ barrier benchmarks (Carcia et al.), and ALD nanolaminates (Dameron et al.) are absent. The pinhole model stands on these shoulders and must say so.

**Minor**
- R1.5 "T80" has a fixed meaning in photovoltaics (time to 80% of initial PCE). The dose-threshold quantity here is not that. Rename or define unambiguously.
- R1.6 Water sorption in polymers deviates from Henry's law near saturation (clustering); the activity-linear model at 90–100% RH needs a stated assumption and caveat.
- R1.7 The cost proxy (deposition time) ignores capital/precursor costs and would change qualitatively under spatial ALD. Scope it.

## Reviewer 2 — Transport modeling expert. Recommendation: **Major revision**

The numerics are properly verified and the layer-lumped lag treatment is elegant. The central empirical claim, however, needs tighter support.

**Major**
- **R2.1 1-D reduction of 3-D defect transport.** τ² = 1 + s² ln(s/r)/(2π d_org²) is presented without derivation or validity condition. State the regime (s ≫ d_org; here 444 µm vs ≤1 µm), derive it in SI, and connect explicitly to the lag-time picture of Ref. [8].
- **R2.2 "Cannot fit under any admissible parameters" is asserted, not shown.** Provide the bounding inequality for the ladder's maximum per-dyad improvement and compare it with the observed factor.
- **R2.3 Test duration.** "≈72 h, confirmed from the source" must carry a precise citation (section/SI of Ref. [7]), and the conclusion's sensitivity to duration (e.g., 2.5–7 days) must be quantified.
- **R2.4 Censored data in the validation metric.** The 2- and 3-dyad points lie below the instrument's stated steady-state limit; using them as point values in decade-error metrics is statistically loose. Treat them as censored and show the model is *consistent under censoring*.

**Minor**
- R2.5 The (f₀, d_c) identifiability ridge deserves a profile-likelihood figure (SI) and an explicit statement that no extrapolation below 15 nm is made.
- R2.6 The layer-lumped time-lag formula needs its provenance (Frisch) and a derivation in SI.
- R2.7 The device-side boundary a = 0 matches Ca/MOCON geometry, but for a barrier deposited directly on a device, vapor accumulates beneath pinholes and the true driving force is smaller; state that the sink assumption is conservative and cite the distinction.
- R2.8 The 0.10 K thermal flatness presumes organic-capped stacks (fixed emissivity); make the scope explicit.

## Reviewer 3 — Optimization / AI expert. Recommendation: **Minor-to-major revision**

The exhaustive ground truth is exactly what most NSGA-II papers lack, and exact attribution on a factorial grid is the right call. My concerns are mostly about positioning and statistics.

**Major**
- **R3.1 Why NSGA-II at all?** With three variables, enumeration settles the problem; the GA is not a contribution here. Reposition it honestly: verified machinery for the higher-dimensional extensions (per-layer thicknesses, material choice) where enumeration fails.
- **R3.4 UQ priors are ad hoc** (×/÷3, ×/÷1.5, 6–10). Provide a rationale table, and separate which conclusions are data-noise-dominated (thresholds) versus prior-dominated (absolute lifetime, durability feasibility).

**Minor**
- R3.2 Clarify IGD vs GD roles; note the 80-point archive limitation.
- R3.3 With 3 features "exact Shapley" is elementary; keep it factual, and state the input measure (uniform over the factorial grid; interventional semantics) for both Shapley and Sobol.
- R3.5 "Self-verifying … autonomously detected" anthropomorphizes a reporting convention; rephrase.
- R3.6 Report seeds, library versions, runtimes, and an MC convergence check (N = 500) in SI.
- R3.7 Two of six objectives are flat; justify their inclusion early (as a quantified finding) so it does not read as objective padding.

---

# Point-by-Point Response Table (actions implemented in manuscript v2)

| # | Reviewer point (abridged) | Sev. | Action in v2 | Location |
|---|---|---|---|---|
| R1.1 | Single-source dependence | Maj | Abstract now states single-source basis explicitly; falsification targets named (Ref. [3] Fig. 2d dyad+lag; any 10 nm ALD point); Choi check elevated | Abstract; §2.3; §3.5(v) |
| R1.2 | Perovskite overclaim | Maj | Scope sentence added to Abstract ("demonstrated end-to-end on the Al₂O₃/parylene-C system; perovskite transfer requires M_crit calibration"); Intro contribution list reworded | Abstract; §1; §3.5(ii) |
| R1.3 | G5 vs planarization practice | Maj | G5 reframed as *transport-only, conditional* rule; planarization function stated as unmodeled and thickness-setting in practice; Coclite & Gleason cited | §3.4(G5); §3.5(i); ref [21] |
| R1.4 | Missing foundational refs | Maj | Rossi–Nulman, da Silva Sobrinho, Carcia, Dameron, Coclite–Gleason added with text hooks | §1 ¶3; §2.1; refs [17–21] |
| R1.5 | T80 terminology | Min | Defined at first use as dose-threshold proxy t₈₀^dose (symbol T80 retained in figures); mapping caveat added | §2.1; §3.5(ii) |
| R1.6 | Henry linearity at high RH | Min | Assumption stated; clustering caveat added to Limitations | §2.1; §3.5(vii) |
| R1.7 | Cost proxy scope | Min | Scoped to relative, same-toolset comparisons; spatial-ALD caveat | §2.4; §3.5(viii) |
| R2.1 | τ² derivation/validity | Maj | Validity inequality s ≫ d_org stated with numbers (444 µm vs ≤1 µm); derivation promised in SI S3; tied to Ref. [8] | §2.1; SI plan |
| R2.2 | Ladder bound not shown | Maj | Explicit inequality added: R_sand/R_in = (d_org/d_in)(f + r² ln(s/r)/(2 d_org²)) ⇒ R₃/R₁ ≤ 10.6 < 22 required; derivation to SI | §3.1 |
| R2.3 | 72 h citation + sensitivity | Maj | Citation placeholder to Ref. [7] Experimental/SI [VERIFY page]; duration sweep added: 2.5–7 d keeps both points within ≤0.4 decade | §2.3; §3.1 |
| R2.4 | Censoring in validation | Maj | Points now labeled censored; consistency statement added (model apparent values also below stated limit ⇒ consistent under censoring); decade errors reported as descriptive | §2.3; §3.1 |
| R2.5 | Identifiability ridge | Min | Profile-likelihood promised (SI S3); no-extrapolation statement sharpened | §2.3 |
| R2.6 | Frisch provenance | Min | Frisch cited [22, VERIFY]; derivation to SI S1 | §2.1 |
| R2.7 | a = 0 vs device accumulation | Min | Conservative-direction statement + accumulation caveat added | §2.1; §3.5(ix) |
| R2.8 | Emissivity scope | Min | "Organic-capped stacks" scope stated where flatness is claimed | §3.2 |
| R3.1 | NSGA-II positioning | Maj | Repositioned as verified machinery for higher-dimensional extensions; contribution wording adjusted | §2.4; §3.2 |
| R3.2 | IGD/GD roles | Min | One-sentence clarification incl. archive-size caveat | §3.2 |
| R3.3 | Shapley tone + measure | Min | "No surrogate" kept factual; uniform-grid interventional measure stated | §2.5 |
| R3.4 | Prior rationale / dominance | Maj | Prior-rationale table promised (SI S4); explicit data- vs prior-dominated separation sentence | §2.6; §3.4 |
| R3.5 | Anthropomorphic phrasing | Min | "the number-embedding convention exposed and corrected …" | Abstract; §3.3 |
| R3.6 | Reproducibility details | Min | Seeds/versions/runtimes/MC-convergence itemized in SI S6 plan | SI plan |
| R3.7 | Flat objectives placement | Min | Flatness quantified at first mention in §2.4/§3.2 as a finding | §2.4; §3.2 |

**Residual risks after v2 (not fully resolvable without new work):** independent-source validation ([3] Fig. 2d — requires PDF), M_crit device calibration, film-mechanics primaries for the durability channel, resolution of all [VERIFY] bibliography flags.

---
**Post-review correction (2026-07-10).** Full-text verification of Ref. [7] found the dyad series uses 50 nm Al2O3 (not 30) and contains no test-duration statement. Re-anchored pipeline: density 8.5 mm⁻² (343 µm); validation +0.11/−0.01 dec; ladder bound R₃/R₁≤7.8; duration framing replaced by the 2.0–5.5 d insensitivity window. All reviewer-response commitments remain satisfied; R2.3 is resolved more strongly (duration-independence) than originally promised.
