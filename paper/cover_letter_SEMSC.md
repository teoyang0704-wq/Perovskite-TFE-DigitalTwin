# Cover Letter (draft) — Solar Energy Materials and Solar Cells

[Date] · Teo Yang · Independent Researcher · teoyang0704@icloud.com

Dear Editor,

Please consider our manuscript, **"Geometry-first design of organic/inorganic multilayer thin-film encapsulation via a literature-calibrated digital twin: defect-mediated moisture physics, zero-refit validation, and confidence-qualified design rules,"** for publication as a full-length article in *Solar Energy Materials and Solar Cells*.

**What the paper reports.** We treat multilayer thin-film encapsulation (TFE) as a *geometry* design problem and build a lightweight digital twin whose sensitivity to layer thicknesses and dyad number derives from explicit pinhole-closure, particle-floor, cracking and lateral tortuous-path physics. Calibrated only on a published single-layer ALD-Al₂O₃ thickness series plus one multilayer point, the model predicts the remaining multilayer barrier data with zero refitting to +0.11/−0.01 decade—but only when a finite multi-day run is simulated (the source states its detection limit, not the duration): we show analytically that steady-state theory cannot reproduce the observed per-dyad improvement under any admissible parameters, identifying widely reported sub-detection-limit multilayer WVTR values as lag-phase-limited measurements. Exhaustive enumeration of the full design space provides a ground-truth six-objective Pareto front; exact Shapley/Sobol attribution is translated into causal chains that carry model-internal quantities (a convention that exposed and corrected one of our own over-interpretations); and Monte-Carlo uncertainty propagated through the *entire calibration pipeline* yields design rules with explicit confidence—chief among them an inorganic-thickness window of 22.5 [21.0–26.0] to 44.0 [37.5–48.0] nm and a saturation-free linear dyad law robust in all 500 draws.

**Why this is new.** Each ingredient exists—defect-dominated permeation theory, laminate models, NSGA-II, SHAP—but no prior work couples a validated transport PDE to encapsulation *geometry* search with reliability-centred objectives, validates it zero-refit against published multilayer data (censoring-consistently), or qualifies its design rules by uncertainty propagated through calibration itself. We deliberately claim no new material, law or algorithm: the contribution is the validated integration and the design map it produces, with named falsification targets for the community.

**Why this journal.** The work sits squarely in this journal's core territory—encapsulation and barrier performance, module-level moisture ingress, and accelerated-stability interpretation—and builds directly on studies published here (ALD encapsulation of perovskite devices; module moisture-ingress modelling). Its outputs are addressed to experimentalists: fabrication-ready thickness windows, a robust design map, and an explanation of a measurement artefact that affects how high-barrier WVTR data are read.

**Why now.** Perovskite photovoltaics are at the encapsulation-limited stage of commercialization, and the community has just articulated the need for multi-property TFE evaluation (2025 figure-of-merit proposal) while machine-learning studies proliferate without mechanism grounding. A transparent, physics-anchored, uncertainty-honest exemplar is timely.

All data, code, database provenance (including digitization overlays) and reproduction commands are released at https://github.com/teoyang0704-wq/Perovskite-TFE-DigitalTwin, with an interactive deployment of the calibrated twin at https://tfe-twin.streamlit.app. Portions of the code, analysis and text were developed with a large-language-model assistant (Anthropic Claude); all methods, results and claims were verified by the author, who takes full responsibility. The manuscript is original, not under consideration elsewhere, and reports no conflicts of interest.

Suggested reviewers: [3–5 researchers spanning (i) TFE/ALD barriers, (ii) permeation modelling, (iii) PV reliability or ML-for-materials; select from the fields of the cited literature].

Word count ≈ [n]; 7 figures; 3 tables; Supporting Information S1–S6.

Sincerely,
Teo Yang
