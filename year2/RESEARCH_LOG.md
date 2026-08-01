# Year-2, module 1 — three-dimensional defect network

Status: solver implemented and validated; first experiment complete.
Author: Teo Yang. Started 2026-08-01.

## 1. Why this module exists

The Year-1 twin replaces the whole lateral problem inside an organic interlayer
with one analytic factor,

    tau^2 = 1 + s^2 ln(s/r) / (2 pi d_org^2),

derived by assuming that pinholes in successive inorganic layers are
statistically independent and that each one feeds its neighbour along a single
radial path. Those assumptions were never tested — they were adopted, calibrated
through, and declared as a limitation. This module tests them by solving the
transport problem on an explicit three-dimensional network of defects, so that
the analytic factor becomes a result rather than an assumption.

## 2. Method

Three lengths are widely separated in every real stack: pinhole radius
(r ~ 50 nm) << organic thickness (d_org ~ 1e2 nm) << pinhole spacing
(s ~ 1e5 nm). Exploiting that separation, each organic interlayer is treated as
a thin conducting sheet discretised on a periodic 2-D grid (sheet conductance
P_org * d_org per square), and each pinhole becomes a vertical terminal linking
two adjacent sheets, carrying the through-hole resistance plus the analytic
sub-grid constriction resistance. A full 3-D mesh resolving 50 nm holes over a
millimetre domain would need ~1e8 cells; this formulation needs ~1e5 and is
exact in the limit r << Delta << s.

Defect correlation enters through one parameter. A fraction phi of defects is
*columnar*: a particle that punctures every inorganic layer at the same (x, y),
represented as a continuous tube with its own node at each level, so that water
crosses the whole stack without ever spreading laterally. The remaining
defects are independent between layers. phi = 0 is the Year-1 assumption.

Files: `percolation3d.py` (solver), `validate_3d.py` (tests 1-3),
`experiment_ladder_bound.py` (experiment 1).

## 3. Validation

| test | result |
|---|---|
| two-electrode sheet problem, analytic answer known | network within **3.6 %**, converged in grid resolution (nx = 64/128/256 identical to 3 digits) |
| regular defect array vs the Year-1 tau^2 formula | **G_network / G_1D = 0.960** — the analytic factor is accurate to 4 % for regular arrays |
| pure columnar limit (phi = 1) vs analytic series bound 1/n | reproduces **0.667 / 0.500** exactly |

The Year-1 detour factor therefore survives its first direct test. That is a
real outcome: the calibrated model was not resting on a broken approximation.

## 4. First results

**(a) Random placement raises the barrier by ~20 %.** At equal defect density,
Poisson-distributed pinholes give conductance 0.801 +/- 0.007 times that of a
regular array (8 realisations). Randomness creates defect-poor regions whose
long detours cost more than close pairs save. The Year-1 formula, which is
implicitly a regular-array result, therefore over-predicts transmission by about
a quarter — a bias absorbed into the calibrated defect density, but one that
matters whenever the model is re-anchored on a new laboratory.

**(b) A structural bound on how flat a dyad ladder can be.** Every arrangement
of defects puts the n inorganic layers in series, so the flattest achievable
ladder is the one where lateral detours contribute nothing and all transport
runs through columnar channels:

    G(n) = G(1) / n     =>     G(n2) / G(n1) = n1 / n2.

The network confirms this: sweeping phi from 0 to 1 moves G(4)/G(2) from 0.353
to exactly 0.500 and no further.

**(c) Published ladders violate the bound.** Lee et al. measure
G(3)/G(2) = 0.818 and G(4)/G(2) = 0.803, against bounds of 0.667 and 0.500.
Wu et al.'s high-quality series (0.214 at n = 3) sits comfortably inside it.
Decomposing the Lee ladder as G(n) = A/n + C gives C = 3.8e-4 g m^-2 day^-1:
**a dyad-independent channel that accounts for 71 % of the four-dyad reading.**

## 5. What (c) means, stated honestly

The bound assumes n *identical* layers in series. The measurement is therefore
incompatible with that picture, and exactly two classes of explanation remain:

1. **a parallel ingress channel that does not scale with dyad number** — most
   plausibly lateral permeation from the sample perimeter in a Ca test, or
   substrate-side ingress; or
2. **systematically degrading layer quality** — later inorganic layers carrying
   more defects than earlier ones, so each addition helps less. Lee's stack is
   not layer-uniform (first organic layer 800 nm, the rest 100 nm), which makes
   this physically plausible.

Neither is in the Year-1 model, and the two are distinguishable by experiment:
a perimeter channel scales with the perimeter-to-area ratio of the test sample
and is unchanged by stack order, while degrading quality is independent of
sample size and inverts when the deposition order is inverted.

This is the Year-2 research question in one sentence: **the improvement of
multilayer barriers saturates faster than any series model permits, and the
reason is a channel nobody is currently measuring.**

## 6. Next steps

1. Sensitivity of the bound violation to r, d_org and density; is C universal
   across published ladders, or sample-geometry-specific?
2. Add an explicit perimeter channel to the network (finite domain with an open
   edge) and fit C from the geometry of the published Ca tests.
3. Re-analyse every ladder in the Year-1 database against the bound; count how
   many violate it. If most do, the field has been mis-reading dyad scaling.
4. Feed the 20 % randomness bias back into the Year-1 calibration and check
   whether the design window moves (it should not; the thresholds are ratios).
5. Only then: hierarchical Bayesian multi-laboratory calibration (module 2).

---

## 7. Experiment 2 — the dyad-independent channel is the sample edge

Script: `edge_ingress.py`. Figure: `fig_edge_channel.png`.

### 7.1 The survey splits by measurement method, not by stack

Every ladder in the Year-1 database was checked against the bound for all pairs
(n1 < n2):

| ladder | test method | pairs violating |
|---|---|---|
| Wu et al., high quality | MOCON coulometric, clamped and edge-sealed | **0 of 3** |
| Wu et al., poor quality | MOCON, values digitised from a log figure (+-30 %) | 1 of 3 (marginal: 0.812 vs 0.667) |
| Lee et al. | calcium test, 3 mm x 3 mm pad, exact values | **3 of 6** (worst 0.981 vs 0.750) |

This is the discriminating observation. Barrier layers that degrade with
deposition order would flatten a ladder regardless of how it is measured; a
channel that bypasses the stack laterally can only appear where the sample edge
is open to the atmosphere. The violation appears in the calcium test and not in
the clamped coulometric test.

### 7.2 The edge channel reproduces the Lee ladder quantitatively

Water reaching the sensor by diffusing in from the cut edge contributes

    WVTR_edge = P_org * da * (perimeter * t_org) / (L_path * A_sensor),

which depends on the test geometry and on the total organic thickness, and so
grows *slightly with dyad number* -- flatter than even the 1/n bound. Fitting
the Lee ladder as series + edge (two free scales) gives an rms residual of
**0.019 decades** over four points:

| n | measured | model | series | edge | edge share |
|---|---|---|---|---|---|
| 1 | 3.0e-3 | 2.9e-3 | 2.6e-3 | 2.7e-4 | 9.5 % |
| 2 | 6.6e-4 | 7.1e-4 | 4.0e-4 | 3.1e-4 | 43 % |
| 3 | 5.4e-4 | 5.6e-4 | 2.2e-4 | 3.4e-4 | 61 % |
| 4 | 5.3e-4 | 5.2e-4 | 1.5e-4 | 3.7e-4 | **72 %** |

The fit requires an organic permeability 22 times that of parylene-C. Year-1
extracted, from an entirely different observable -- the >400 h lag times of the
same stack -- an effective value of about 38 times parylene. Two independent
routes to the same unmeasured material property agree within a factor 1.7, and
they would coincide exactly for an in-plane path length of 8.7 mm rather than
the 5 mm assumed here (the papers do not state it).

### 7.3 Consequence: a measurement design rule

Edge ingress scales as 1/(pad side), so small sensor pads impose a floor:

| calcium pad | floor, 1.1 um organic | floor, 0.5 um organic |
|---|---|---|
| 1 mm | 1.1e-3 | 5.1e-4 |
| 3 mm | 3.7e-4 | 1.7e-4 |
| 10 mm | 1.1e-4 | 5.1e-5 |
| 25 mm | 4.5e-5 | 2.0e-5 |
| 50 mm | 2.3e-5 | 1.0e-5 |

A stack whose true WVTR lies below the entry for its pad size is measuring its
own edge. Lee's four-dyad reading of 5.3e-4 sits a factor 1.4 above the floor
of its own 3 mm pad -- the measurement had reached its geometric limit, and the
barrier itself is better than the number implies.

### 7.4 What is established, and what is not

Established: published ladders violate a bound that no series model can escape;
the violation tracks the measurement method; an edge channel of the right order
reproduces the flattening and independently recovers a material parameter
obtained from lag times.

Not established: the in-plane path length L_path is assumed, and the fitted
permeability scales inversely with it. The series and edge scales are partly
degenerate, since the unknown defect density also multiplies the series term;
a joint re-anchoring in which one organic permeability serves both is the next
refinement. Degrading layer quality is disfavoured by the method split but not
excluded -- inverting the deposition order would settle it experimentally.

### 7.5 Next

1. Joint re-anchoring: single P_org, free defect density, both ladders at once.
2. Add the edge channel to the Year-1 twin as an explicit term and re-derive the
   design rules for finite modules; the organic interlayer is now doubly
   penalised for thickness (lateral sheet conduction and edge cross-section),
   which sharpens rule G5 rather than softening it.
3. Extend the survey: every multilayer WVTR in the literature that was measured
   on a small calcium pad is a candidate for edge limitation.

---

## 8. Experiment 3 — joint re-anchoring, and what a finite device can use

Script: `joint_reanchor.py`. Figure: `fig_finite_device.png`.

### 8.1 The degeneracy breaks, and the fit improves

Forcing one organic permeability to serve both paths removes the ambiguity
flagged in §7.4, because the edge term carries no defect fraction:

    series ~ f * P_org        edge ~ P_org

Fitting (f, m = P_org / P_parylene) jointly to the four-point ladder gives
**f = 4.2e-8, m = 22.9, rms residual 0.011 decades** -- a better fit than the
two-free-scale version, with one parameter fewer.

| n | measured | model | series | edge | edge share |
|---|---|---|---|---|---|
| 1 | 3.0e-3 | 2.9e-3 | 2.6e-3 | 2.9e-4 | 9.8 % |
| 2 | 6.6e-4 | 6.8e-4 | 3.6e-4 | 3.2e-4 | 47 % |
| 3 | 5.4e-4 | 5.5e-4 | 1.9e-4 | 3.6e-4 | 65 % |
| 4 | 5.3e-4 | 5.2e-4 | 1.3e-4 | 3.9e-4 | **75 %** |

### 8.2 A correction to the Year-1 external validation

Year-1 concluded, from the same ladder, that the second laboratory's film
carried **15 times** the defect density of the calibration film, and read that
as the low-temperature deposition penalty the source itself documents. The
joint analysis says otherwise. Once the organic layer is allowed its own
permeability and the edge is accounted for, the inferred defect fraction falls
by a factor 13, to 4.2e-8 -- **1.1 times the calibration film**, a pinhole
spacing of 432 um, 5.4 defects per mm2.

The two laboratories' films are essentially equally clean. The apparent
fifteen-fold penalty was an artefact of two assumptions made explicit but not
tested in Year-1: that pV3D3 transports like parylene-C, and that all measured
flux crosses the stack.

The same correction says the reported barrier understates the stack. On a large
sample, free of edge limitation, the four-dyad stack would read
**1.3e-4 instead of 5.3e-4** g m^-2 day^-1 -- four times better than published.

Third consistency check on m: lag times gave ~38, the ladder shape gives 22.9,
and the two coincide for an in-plane path of 8.7 mm rather than the 5 mm
assumed. Three observables, one material parameter, agreement within a factor
1.7.

### 8.3 Design rule: dyads stop paying below a device size

Year-1 derived its rules for unbounded area, so the lifetime gain per dyad was
linear without end. With a perimeter in the model, the gain saturates -- and the
saturation is set by device size, not by stack design:

| device side | dyads still worth adding (>10 % gain) | floor WVTR | limited by |
|---|---|---|---|
| 3 mm | 3 | 5.5e-4 | **edge** |
| 10 mm | 4 | 2.5e-4 | stack |
| 30 mm | 6 | 1.3e-4 | stack |
| 100 mm | 8+ | 7.4e-5 | stack |
| 1000 mm | 8+ | 6.0e-5 | stack |

Below roughly a centimetre the perimeter, not the stack, is the barrier. This
is the first rule in the project that depends on the size of the object being
protected, and it sharpens Year-1's preference for thin organic interlayers:
they are penalised twice over, once for conducting laterally between pinholes
and once for widening the edge cross-section.

### 8.4 Consequence for the Year-1 manuscript

The submitted manuscript reports that the mode-B residuals against this ladder
are all one-sided, and attributes them to a decorrelated-defect model bounding
transmission from below. That attribution is now superseded: the residuals are
the sample edge, and they can be reproduced quantitatively. If the manuscript
returns for revision, replacing a declared but unexplained systematic with a
measured mechanism is a straightforward strengthening; no published claim needs
retracting, since every affected statement was made conditional on the parylene
proxy and on all flux crossing the stack.

### 8.5 Next

1. Literature survey: how many published multilayer WVTR values were measured
   on small calcium pads, and which of them sit below their own geometric floor?
2. Deposition-order inversion: the one experiment that would separate the edge
   channel from degrading layer quality outright.
3. Module 2: hierarchical Bayesian calibration across laboratories, now that the
   between-laboratory spread has shrunk from 15x to 1.1x.

---

## 9. Module 2 — hierarchical calibration across films

Scripts: `hierarchical.py` (model and sampler), `hier_populations.py` (two
populations, figure). Figure: `fig_hierarchical.png`.

Year-1 calibrated on one film with two point-estimate stages, so every rule
inherited that film's cleanliness. Here the defect floor of each film is a draw
from a population, all films are fitted at once, and the shared physics
(nucleation amplitude, closure length, cracking onset) is common to all:

    log10 f_res,i ~ Normal(log10 mu_f, sigma_f)

`sigma_f`, the between-film spread in decades, is the number that decides how
far any calibrated rule can travel. No sampler was available in this
environment, so the affine-invariant ensemble sampler is implemented in the
module (about forty lines); acceptance 0.32-0.37, split-half drift in the
population parameters below 0.05 sd.

### 9.1 Two populations, because one film is not like the others

Wu's poor-quality film was deposited deliberately at high base pressure to show
that particulates ruin a barrier. Pooling it with films meant for service
inflates the spread, so both are reported.

| | all published films | service-quality only |
|---|---|---|
| between-film spread sigma_f | **1.32** [0.96, 1.83] dec | **0.76** [0.45, 1.29] dec |
| new film, predictive f_res | 9.9e-8 [3.5e-10, 2.6e-5] | 1.8e-8 [5.4e-10, 7.0e-7] |
| window lower edge | 21.9 [14.2, 29.0] nm | **24.2** [19.9, 29.2] nm |
| window upper edge | 39.9 [32.8, 44.3] nm | **43.3** [36.8, 48.6] nm |
| window width | 17.6 [8.0, 27.0] nm | 18.8 [11.2, 25.9] nm |
| chance of an empty window | 0.3 % | 0.0 % |

Per-film floors (service population): Carcia 1.8e-9, Lee 3.2e-8, Wu high-quality
7.3e-8 -- a factor 40 between the cleanest and the least clean film that anyone
intended to use, and a factor 1600 up to the deliberately bad one.

### 9.2 The Year-1 window survives, with one correction

Year-1 reported 22.5 [21.0-26.0] to 44.0 [37.5-48.0] nm from a single film. The
hierarchical analysis over three independent films gives 24.2 [19.9, 29.2] to
43.3 [36.8, 48.6] nm. The upper edge is essentially unchanged, since cracking is
mechanical and film-independent; the lower edge moves up by about 2 nm, and its
interval widens.

The direction of that shift is the interesting part. Closure is complete when the
nucleation term has fallen to a fixed fraction of the particulate floor, so a
*cleaner* film -- a lower floor -- requires the nucleation term to decay further:

    d_closure = d_close * ln( f0 / (0.01 f_res) )

**A cleaner film needs a thicker inorganic layer before its defects are closed,
and its usable window is narrower.** Year-1's window came from a mid-range film
and was, by luck, close to the population median; a laboratory with unusually
clean deposition should design to about 29 nm rather than 22 nm.

### 9.3 Third independent estimate of the pV3D3 permeability

The hierarchical fit, which sees the whole dataset rather than one ladder,
returns 23.5 [7.9, 65.8] times parylene-C. Module 1 gave 22.9 from the ladder
shape; the Year-1 lag analysis gave about 38. Three routes, one material
parameter, all inside a factor two.

### 9.4 Honest limits

Four films is a thin population, and sigma_f is correspondingly uncertain -- its
90 % interval spans a factor of two in decades. The predictive interval for a
new film is therefore wide by construction, and that width is the result, not a
failure of it: with the literature as it stands, a laboratory that has not
measured its own film cannot know its defect floor to better than about two
decades. Narrowing that is a matter of adding films, not of better inference.

### 9.5 Next

1. Add every single-layer thickness series in the literature to the population;
   each new film tightens sigma_f directly.
2. Propagate the predictive floor through the Year-1 lifetime map: design rules
   conditional on "your film is as clean as the median" versus "worst case".
3. The literature survey of module 1 (which published values sit below their own
   edge-limited floor) is still open and needs web access.

---

## 10. Experiment 4 — the edge floor has already been measured, twice

Script: `literature_survey.py`.

A literature search for calcium-test multilayer data turned up two groups who
measured the channel of module 1 directly, published the number, and did not
draw the conclusion.

**Bülow et al., Nanoscale Res. Lett. 9, 223 (2014)** (open access) ran a glass
lid -- an essentially impermeable barrier -- through the same calcium fixture as
their multilayers. It read about 6e-4 g m^-2 day^-1 at 60 C / 90 %RH. Their
best 3.5-dyad stack read 1.2e-3, which the paper describes as "only by a factor
of 2 higher" than the glass lid. Read the other way: **half of that measurement
is the fixture.** Their ladder (3.6e-3 at 1.5 dyads, 1.2e-3 at 3.5) satisfies
the series bound both before and after the measured floor is removed, which is
what a ladder should look like when the floor is known and modest.

**A doctoral thesis on calcium testing of ALD barriers** reports an ALD film and
a glass lid both returning about 1e-4, "assigned to H2O permeation through the
epoxy edge seal", and states plainly that the calcium test cannot evaluate
barriers at or below 5e-5 because the epoxy edge seal is the limiting factor.

Both fixtures seal a lid with epoxy, so their channel runs through adhesive.
Lee et al. deposit the encapsulation straight onto the calcium with no lid and
no epoxy: nothing in that geometry can leak except the stack's own organic
interlayers, which is precisely the channel module 1 proposed and quantified.

| fixture floor | value | how obtained |
|---|---|---|
| Bülow, glass lid + epoxy | 6e-4 | measured, impermeable control |
| thesis, ALD or glass lid | 1e-4 | measured, attributed to the epoxy seal |
| Lee, no lid | 3.9e-4 | inferred here from the ladder shape alone |

Three independent instances inside one order of magnitude, two of them measured
with a control, the third inferred from shape alone and landing in the same
range.

### 10.1 What is new, stated precisely

The phenomenon of a fixture-limited calcium test is known for lids sealed with
epoxy, and at least two groups have said so in print. This work adds three
things:

1. it occurs in **lid-free, directly deposited** stacks, where the leak path is
   the encapsulation's own organic interlayers and no adhesive is involved;
2. it can be detected **from the shape of a dyad ladder alone**, with no control
   sample, because a fixture channel forces a violation of the series bound;
3. its size follows from the geometry, which turns it into a **design rule for
   the measurement** -- a sensor of a given size cannot resolve below a floor
   that can be calculated in advance.

### 10.2 Two diagnostics that need no new hardware

Run the series bound on any published ladder: a violation proves something is
parallel to the stack. And compare the best reading with an impermeable control
measured in the same fixture -- the control is the floor, not zero. The one
experiment that would settle the remaining ambiguity does need hardware:
inverting the deposition order leaves a fixture channel unchanged and reverses
degrading layer quality.

### 10.3 Data added to the population

Bülow's ladder and single-layer thicknesses (TALD 100 nm 6.4e-3, PEALD 25 nm
4.4e-3, at 60 C) are candidates for the hierarchical population of module 2,
but they are measured at 60 C rather than 38 C, so folding them in requires the
activation energy to be carried as a fitted parameter rather than a prior. That
is the next increment.

---

## 11. Module 2b-d — five films, two temperatures, and a correction to Year-1

Scripts: `hier_v2.py` (activation energy free, Buelow added), `hier_v3.py`
(closure length hierarchical), `window_value.py` (per-film windows).
Figure: `fig_window_value.png`.

### 11.1 Adding a fifth film broke the shared-shape assumption, visibly

Buelow's films are plasma-ALD on PEN, measured at 60 C, so the activation energy
had to become a parameter. Carcia's two temperatures on one film supplied the
constraint directly, giving **42 [28, 55] kJ/mol** -- consistent with the 52.5
obtained from that pair alone in module 1, and with the 40 +/- 10 prior Year-1
had assumed.

With one closure length shared by all films, the fit missed Buelow's 25 nm point
by 1.4 decades and the strain propagated: the cracking onset was dragged down to
29 nm and the design window came out empty in 16 % of draws. Letting the closure
length vary per film -- it measures how fast nucleation islands merge, which
depends on substrate, precursor and whether the plasma is on, and is not a
constant of alumina -- removed the tension. Fit quality over all twenty
observations: **rms 0.198 decades**. The cracking onset returned to
**41.5 [31.5, 47.6] nm**, close to Year-1's 44.0.

Between-film spreads: **0.84 [0.46, 1.83] decades in the defect floor**, and
**0.53 [0.31, 0.91] decades in the closure length**.

### 11.2 The design window is a property of a deposition campaign

| film | data available on it | window lower edge (nm) | window width (nm) |
|---|---|---|---|
| Wu high quality | five thicknesses | **22.7 [20.9, 24.7]** | 18.8 [8.6, 25.0] |
| Buelow | three thicknesses | 73.7 [4.7, 86.0] | **negative** |
| Wu poor quality | none | 168 [42, 268] | **negative** |
| Lee | none | 11.0 [1.8, 26.5] | 29.6 [12.5, 41.4] |
| Carcia | one thickness, two temperatures | 15.4 [2.6, 29.3] | 25.2 [8.7, 39.9] |
| a film nobody has measured | nothing | 30.5 [2.4, 291] | 9.8 [-250, 40] |

Two readings of this table matter.

First, **Year-1's window is exactly right for the film Year-1 measured**: 22.7
[20.9, 24.7] nm here against 22.5 [21.0, 26.0] nm there. What was wrong was
presenting it as a property of the material rather than of that deposition run.

Second, **two of the five films have no usable window at all**. Wu's poor film
never closes before it cracks -- which is the point that paper was making --
and Buelow's plasma-ALD film is marginal by the same measure. A window is not
guaranteed to exist.

### 11.3 How much measurement is enough

| | 90 % interval on the lower edge | chance of no window |
|---|---|---|
| with a five-point thickness series | **3.8 nm** | 0.5 % |
| with nothing | 288 nm | 41 % |

**A five-point single-layer thickness series narrows the design window by a
factor of about eighty.** Five depositions and five barrier readings -- the
experiment Wu already did -- convert an unusable prior into a design rule. That
is the concrete prescription this project can now offer a laboratory, and it is
a better deliverable than a universal number that does not exist.

### 11.4 Honest limits

Only two of the five films come with a thickness series, so the spread in
closure length is estimated from essentially two informative draws and its
interval is correspondingly soft. Split-half convergence diagnostics for the
population parameters sit at 0.3-0.4 sd, adequate for the numbers quoted but not
for tighter claims; a longer chain is warranted before publication. The shared
cracking onset is an assumption that has not yet been stressed the way the
closure length was, and it will fail in the same way if a film with very
different mechanics is added.

### 11.5 Next

1. Longer chains and a proper convergence statistic across independent runs.
2. More thickness series from the literature -- each one directly tightens the
   closure-length population, which is now the binding uncertainty.
3. Propagate the per-film window into the Year-1 lifetime map, so that lifetime
   predictions carry the campaign-specific uncertainty rather than a single
   film's.

---

## 12. Module 2e — the universal rule did not disappear, it changed units

Script: `universal_criterion.py`.

Module 2d looked like the loss of a universal design rule. It was not. Closure
is complete when the nucleation term reaches a fixed fraction of the particulate
floor, so the closure thickness is a *number of closure lengths*:

    d_closure = N * d_close,        N = ln( f0 / (0.01 f_res) )

and N is a logarithm, so it barely moves. Across all five films N sits between
19 and 23; a film a hundred times cleaner needs only 4.6 more closure lengths.
Essentially all of the between-film variation in the window comes from the
closure length itself, not from the defect floor.

The window therefore exists if and only if

    **d_crit / d_close > N   ,  i.e.   d_close < 2.0 [1.4, 2.7] nm**

a criterion with no thickness in it. Above that closure length the film cracks
before its defects have closed, and no inorganic thickness works at all.

| film | closure length (nm) | P(a window exists) |
|---|---|---|
| Lee | 0.55 [0.09, 1.30] | 99.9 % |
| Wu high quality | 1.17 [0.98, 1.42] | 99.5 % |
| Carcia | 0.67 [0.11, 1.29] | 99.3 % |
| Buelow | 3.29 [0.24, 4.23] | **25.6 %** |
| Wu poor quality | 8.37 [3.34, 14.16] | **3.7 %** |
| a film nobody has measured | 1.50 [0.12, 14.3] | 58.6 % |

The criterion reproduces, from one number per film, exactly which films had a
usable window in the table of §11.2.

### 12.1 The revised statement of what this project knows

Universal, film-independent:

* the series bound on dyad ladders, and its use as a fixture-leak diagnostic;
* the linear lifetime gain per added dyad;
* the preference for thin organic interlayers (lateral conduction and edge
  cross-section both penalise thickness);
* the device size below which the perimeter, not the stack, is the barrier;
* closure takes N ~ 21 closure lengths;
* a window exists iff the closure length is below about 2 nm;
* the cracking onset, 41.5 [31.5, 47.6] nm -- shared by assumption, not yet
  stressed the way the closure length was.

Campaign-specific, and knowable only by measuring:

* the closure length, hence the lower edge of the window in nanometres;
* the particulate floor, hence the absolute lifetime.

Year-1 reported the second kind as though it were the first. The correction is
not that the number was wrong -- for Wu's film it was right to within a
nanometre -- but that it was a measurement of one deposition run.

---

## 13. Module 2f — a sixth film, and a cap the model was missing

Script: `hier_v4.py`.

Groner, George, McLean and Carcia (48th SVC Proceedings, 2005, open access)
published the thickness series the population needed, running straight through
the steep region: 2.5 nm gives no improvement over the bare polymer, 5 nm about
one decade, 10 nm 2e-3, 26 nm 1e-3 g m^-2 day^-1.

Adding it exposed a gap. `WVTR = K f(d)` has no upper bound, so for a
barely-nucleated film the model predicted fluxes far larger than the bare
substrate can pass. No earlier data point sat near that limit, so the omission
was invisible; Groner's 2.5 nm point sits exactly there. Putting the substrate
in series -- as the multilayer model already did -- fixes it.

With six films and 24 observations the fit reaches **rms 0.167 decades**.

| quantity | five films | six films |
|---|---|---|
| spread in closure length | 0.53 [0.31, 0.91] dec | **0.50 [0.32, 0.81]** dec |
| critical closure length | 2.01 [1.41, 2.71] nm | **2.11 [1.51, 2.77]** nm |
| P(window) for an unmeasured film | 58.6 % | **65.5 %** |
| cracking onset | 41.5 [31.5, 47.6] nm | 42.8 [33.7, 48.3] nm |
| activation energy | 42 [28, 55] kJ/mol | 45 [24, 62] kJ/mol |

Per-film closure lengths now: Groner 0.61 [0.44, 0.88], Carcia 0.59, Lee 0.68,
Wu high quality 1.20, Buelow 3.46, Wu poor 8.79 nm. Four films of six have a
usable window with near-certainty; Buelow's plasma film is at 21 % and Wu's
deliberately poor film at 0 %.

One added film moves a population by very little, which is the honest lesson:
the criterion tightened from 2.01 to 2.11 nm with a slightly narrower interval,
and the predictive probability of having any window at all rose by seven points.
Pinning this down is a matter of ten films, not two.

### 13.1 A third measurement technique, the same signature

Groner also coated Kapton on both sides -- two barrier layers in series -- and
the WVTR fell only from 1e-3 to 7e-4, a ratio of 0.70 against a series bound of
0.50. **Violated.** Their lag time meanwhile rose from about one day to four or
five, and the paper cites Graff's lag-versus-equilibrium argument for exactly
this reason.

That makes three techniques -- coulometric MOCON, calcium, and a tritiated-water
cell -- of which the latter two show ladders flatter than any series arrangement
allows, and in both cases the sample is sealed against a fixture: an epoxy bead
in one, a clamped Viton o-ring in the other. The MOCON ladders, where the
measured area is edge-sealed by design, do not violate the bound.

---

## 14. Module 3 — first pass of the literature audit

Script: `audit_registry.py`. Eight ladders entered, every value taken from text,
abstract or figure legend rather than digitised from a curve (one exception,
marked).

### 14.1 Q1 — the series bound

| ladder | seal | verdict |
|---|---|---|
| Wu 2018 high quality | MOCON, clamped | ok |
| Wu 2018 poor quality | MOCON, clamped | violated (2->3: 0.812 vs 0.667; digitised) |
| Lee 2018 | calcium, **no lid** | violated (2->4: 0.803 vs 0.500) |
| Buelow 2014 | calcium, glass lid + epoxy | ok |
| Groner 2005, one vs two sides | tritium cell, Viton o-ring | violated (0.700 vs 0.500) |
| Graham 2011, SiNx/parylene | calcium, lid + polyisobutylene | ok |
| Graham 2011, ZnO/Al2O3 | calcium, lid + polyisobutylene | ok |
| Kim 2014, SiOx/parylene | calcium | ok |

**Three of eight ladders violate a bound that no arrangement of defects can
escape** -- 29 % of all comparable pairs.

### 14.2 Q2 — values sitting on their own floor

Only three of the eight papers report a fixture floor, obtained by running an
impermeable control or by testing the seal alone. In **all three**, the best
value in the ladder lands within a factor of two of that floor:

| ladder | best value | own floor | ratio |
|---|---|---|---|
| Buelow 2014 | 1.2e-3 | 6.0e-4 | 2.0 |
| Graham 2011, SiNx/parylene | 4.0e-5 | 5.0e-5 | 0.8 |
| Graham 2011, ZnO/Al2O3 | 5.0e-5 (censored) | 5.0e-5 | 1.0 |

The other five ladders report no floor at all, which means their best values
cannot be checked by anyone, including their authors.

### 14.3 The floors, in the authors' own words

| floor | source | how obtained |
|---|---|---|
| 6e-4 | Buelow 2014 | glass lid through the same calcium fixture |
| 1e-4 | Colorado thesis | ALD film and glass lid both read this; attributed to the epoxy edge seal |
| 5e-5 | Graham 2011 | separate test of the polyisobutylene edge sealant |
| 2e-6 | Graham 2011 | same laboratory, barrier deposited straight onto the sensor |
| 3.9e-4 | Lee 2018 | not measured; inferred here from the ladder shape |

Graham's group measured both ends of a **factor-300 range with one method**,
differing only in whether the sample was sealed to a lid or the barrier was
deposited directly on the sensor. The floor of a calcium test is a property of
its seal, and it spans more than two decades across the literature.

### 14.4 Why these numbers are not yet a result

The sample is biased, and the bias runs the wrong way. Every one of these eight
papers was found by searching for edge effects, fixture limits and dyad ladders
-- that is, by hunting where the phenomenon lives. Papers that measure their
fixture floor are exactly the papers whose authors suspected they were near it,
which is why three out of three sit on it. Quoting 38 % as the violation rate of
the field would be indefensible.

**The audit only becomes a result with an unbiased sample.** The protocol has to
be fixed before the data are gathered:

1. define the frame -- for example every experimental paper reporting a
   multilayer barrier WVTR in a stated set of journals and years, or the full
   reference list of a recent review;
2. define inclusion -- at least two stack thicknesses or dyad counts with
   numerical WVTR values, any measurement method;
3. record for each: the ladder, the method, the seal, the sensor area, whether a
   floor or control is reported, and whether values are censored;
4. run the bound and the floor comparison on everything included, with no
   discretion after the fact.

Everything needed to execute that already exists in this module; what is missing
is the sampling, which is deliberate work rather than opportunistic search. A
frame of twenty to forty papers is enough to state a rate with a credible
interval, and the answer decides the scale of the claim: a violation rate near
10 % is a cautionary note, a rate near half is a statement about the field.

---

## 15. Module 3, second pass — the pre-registered frame, and what it found instead

Files: `AUDIT_PREREGISTRATION.md` (fixed before any paper was gathered),
`audit_screening.py`.

The pilot of §14 could not be quoted, because its eight papers were found by
searching for edge effects. The protocol was therefore written first: a fixed
set of six queries, none of which mentions edges, seals, fixtures, leakage or
detection limits; every returned paper screened in the order returned; inclusion
and extraction rules fixed in advance; thresholds for what each outcome would
mean fixed in advance.

All six queries have been executed and the frame saturated: queries 4-6 returned
almost only papers already screened. Twenty in-frame papers, four out of frame
(two reviews, seven patents, three vendor pages).

| screening outcome | count | share |
|---|---|---|
| usable ladder, two or more numerical values for one architecture | 4 | 20 % |
| **best value only; the ladder appears as a curve with no numbers** | 12 | **60 %** |
| comparison confounded by another change, or a method paper | 4 | 20 % |
| **reports a fixture floor or an impermeable control** | **1** | **5 %** |

**Eighty per cent of benchmark papers publish a number whose ladder a reader
cannot check** (90 % interval 63-91 %), and one of twenty ran an impermeable
control (90 % interval 1-18 %). The estimates are unchanged from the first ten
papers, with intervals now half as wide.

### 15.1 The audit found a different result from the one it was designed to find

The plan was to estimate a violation rate. That is not going to be estimable at
the rate papers accumulate, because the series bound needs two numbers from one
architecture and most papers give one. But the reason it is not estimable is a
stronger statement than the rate would have been, and it needed no digitising,
no modelling and no assumptions:

> A barrier value is published. The ladder behind it is shown as a curve with no
> numbers. No impermeable control is run. Neither the reader nor the author can
> tell whether the number describes the film or the fixture.

Both diagnostics this project produced -- the series bound and the floor
comparison -- are cheap, need no new hardware, and are unusable on most of the
literature purely because of what is left out of it. The pilot shows what
happens when the information is there: of eight ladders, three violate the
bound, and of the three papers that measured their own floor, all three sit on
it.

This was pre-registered as a possible outcome: the protocol states that papers
excluded because their ladder exists only inside a figure are to be counted,
"because the count of papers whose ladders exist only inside a figure is itself
a reportable finding."

### 15.2 Status and what remains

The frame is complete and saturated at twenty papers, which is the size the
accessible benchmark literature on multilayer ladders turns out to be through
this route. A larger frame would need full-text access to paywalled figures
rather than more searching. The result points at a recommendation rather than an
indictment:

1. publish the ladder numerically, not only as a curve;
2. run one impermeable control in the same fixture and publish its value;
3. state the sensor area and how the sample is sealed.

None of the three costs an experiment worth mentioning, and together they make
every published barrier value checkable by its reader. That is the paper this
audit is turning into, and it is a more useful one than a violation rate.

---

## 16. The "stagnation" angle: tested, and what it turned into

Scripts: `floor_prediction.py`. Attempted 2026-08-01.

The idea was to convert the paper from a negative claim ("some published values
are contaminated") into a positive one ("barrier performance stopped improving
years ago and nobody noticed"). It requires being able to say what a given
fixture's floor is, so that reported bests can be compared against it. We tried
to predict the four measured floors from the geometry of the seal and tabulated
sealant permeabilities, with nothing fitted.

**It failed, and the failure is informative.** Predictions exceed measurements
by factors of 6 to 50, all in the same direction, with a scatter of 0.34 decades
about a geometric mean of 20x. A literature check found the reason, or most of
it: real calcium cells contain a **getter** -- a BaO or CaO desiccant bonded to
the lid to absorb residual and ingressing moisture -- which no permeation model
without it can be expected to match. Seal geometry is also specified in the
trade by bondline *height* and *width*, quantities that papers almost never
report and that we had to guess.

So the floor of a fixture cannot currently be predicted from published
information. That is not a failure of the expression; it is a consequence of
three quantities going unreported: bondline dimensions, getter presence and
capacity, and sensor area. The result strengthens the paper's recommendation
rather than adding a new claim, and it adds a fourth item to it.

**The stagnation claim is therefore not available.** Reported bests span 1e-6 to
1e-3 with a year-to-performance correlation of -0.22 across 2004-2025, which is
consistent with stagnation -- and equally consistent with a mixture of
incomparable fixtures, which is exactly what the audit says the literature is.
Distinguishing the two requires knowing each paper's floor, which requires
information the papers do not contain. We record the attempt and stop, rather
than presenting a correlation that cannot be interpreted.
