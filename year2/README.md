# Year 2 — a series bound for multilayer permeation barriers

This directory contains everything behind the manuscript *"A series bound for
multilayer permeation barriers, and what it reveals about published
measurements"* (Teo Yang, under consideration). No new experimental data were
collected; every input is a value published elsewhere, recorded here with its
source.

The manuscript text is withheld while it is under review. Code, data and the
audit protocol are not.

## Contents

```
year2/
  AUDIT_PREREGISTRATION.md   protocol fixed before any paper was gathered
  RESEARCH_LOG.md            full development record, including failed attempts
  DISCOVERY_PROGRAM.md       a line of enquiry that was tested and abandoned
  code/                      22 scripts; every number in the paper comes from one
  data/                      machine-readable registries (CSV)
  figures/                   the four figures, as produced by the scripts
```

## Reproducing every number in the paper

Run from `year2/code` with numpy, scipy and matplotlib installed.

| claim in the paper | value | script |
|---|---|---|
| bound attained exactly for nonlinear diffusion | 0.500000 / 0.250000 | `universality.py` |
| bound satisfied by random walk and percolation | 0.487 / 0.352 | `universality.py` |
| network solver vs analytic two-electrode result | 3.6 % | `validate_3d.py` |
| grid and domain convergence | 3 digits / 2 % | `validate_3d.py` |
| random vs regular defect placement | 0.801 ± 0.007 | `validate_3d.py` |
| bound never crossed, 12 geometries | max excess 0.00000 | `experiment_ladder_bound.py` |
| joint re-anchoring of a published ladder | rms 0.011 decades | `joint_reanchor.py` |
| edge share at four dyads | 71–76 % | `path_sensitivity.py` |
| permeability required vs path length | 2.2–117 × parylene | `path_sensitivity.py` |
| process scatter from a published series | 0.06 decades | `gradient_probability.py` |
| gradient pattern in 200,000 draws | 0 occurrences | `gradient_probability.py` |
| lag: stack vs lateral path vs measured | 10 h / 101 h / >400 h | `referee_answers.py` |
| gradient shortens the lag as it strengthens | 13.7 → 11.6 h | `close_disjunction.py` |
| ladder survey, violations | 3 of 8 | `audit_registry.py` |
| audit, values a reader cannot check | 80 % [63–91 %] | `audit_screening.py` |
| audit, papers reporting a control | 1 of 20 | `audit_screening.py` |
| measured fixture floors, range | 2e-6 to 6e-4 (300×) | `literature_survey.py` |
| floors are not predictable from published information | 6–50× discrepancy | `floor_prediction.py` |

## Data files

`ladder_registry.csv` — every published ladder used, with the transmission rate
at each repeat count, the test conditions, and where the number appears in the
source.

`measured_floors.csv` — fixture floors, separated into those measured with an
impermeable control and the one inferred here.

`audit_screening_log.csv` — the complete screening log of the pre-registered
audit: every paper returned by the fixed queries, whether it fell in frame, its
screening outcome, and whether it reported a ladder or a control. Papers
excluded are listed with the reason, as the protocol requires.

## Scripts not used in the paper

`hierarchical.py`, `hier_v2.py`, `hier_v3.py`, `hier_v4.py`,
`hier_populations.py`, `window_value.py`, `universal_criterion.py` implement a
hierarchical calibration across six published films. That work is reported
separately. `closure_mechanism.py` and `DISCOVERY_PROGRAM.md` record a question
about defect closure that the available data cannot answer, together with the
measurement that would answer it; they are kept because the negative result is
worth having.

## Licence

MIT, as for the rest of this repository.
