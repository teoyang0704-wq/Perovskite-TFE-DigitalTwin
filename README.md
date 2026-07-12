# Perovskite-TFE-DigitalTwin
Physics-informed digital-twin framework for thin-film encapsulation (TFE) geometry design:
literature-calibrated defect-mediated 1-D moisture/heat PDE engine, 6-objective Pareto
exploration (self-implemented NSGA-II + exhaustive ground truth), exact Shapley/Sobol
physics-mapped interpretation, and calibration-through uncertainty quantification.

**Live demo:** https://tfe-twin.streamlit.app  (interactive twin: design sliders, Pareto map, M_crit tiers)

## Reproduce
pip install -r requirements.txt
1) src/tfe_physics_engine.py      - PDE engine demo + verification (0.45% vs analytic)
2) src/step5_calibrate.py         - calibration & zero-refit validation (Wu, RSC Adv. 2018)
3) optimization/step6_optimize.py - Pareto front (grid ground truth + NSGA-II)
4) explainability/step7_xai.py    - exact Shapley/Sobol + causal-chain auto-report
5) uq/step8_uq.py                 - Monte-Carlo through calibration + robust guidelines

## Data provenance
data/*.csv: literature-sourced with per-row source IDs, censoring flags, quality grades,
calibration/validation role split; digitization overlays in data/provenance/.
Known TODO before submission: verify flagged (needs_check) bibliography; confirm MOCON
test duration in Wu 2018 SI; primary sources for parylene D/S and Al2O3 film mechanics.

## License
MIT (see LICENSE; fill in the copyright holder). Cite via CITATION.cff.
