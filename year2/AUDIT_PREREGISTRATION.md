# Pre-registration: audit of published multilayer barrier ladders

Written **before** any paper was gathered under this protocol, and not modified
afterwards. The eight ladders already in `audit_registry.py` were found
opportunistically while chasing the edge channel and are **excluded** from the
frame; they are reported separately as the pilot that motivated the audit.

Date: 2026-08-01. Author: Teo Yang.

---

## 1. Why a frame is needed

The pilot found that three of eight published ladders violate the series bound
and that three of three papers reporting a fixture floor sit on it. Those
numbers cannot be quoted, because every one of those papers was found by
searching for edge effects and fixture limits. Selection was informative with
respect to the outcome. The rate has to be re-measured on a sample chosen
without reference to the phenomenon.

## 2. Frame

Papers returned by the fixed query set in section 3, taken **in the order the
search engine returns them**, deduplicated by DOI or title, until the frame is
exhausted or forty papers have been screened. Search ranking favours
well-cited work, so the frame is a sample of the results the field treats as
benchmarks rather than of all published barrier measurements. That is stated as
a limitation and is arguably the population of interest: these are the numbers
other people design against.

No query in the set mentions edges, seals, fixtures, leakage, artefacts or
detection limits. This is the debiasing step and it is the reason the query list
is fixed in advance.

## 3. Query set (fixed)

1. `multilayer thin film encapsulation WVTR dyads Al2O3 organic layer values`
2. `nanolaminate barrier water vapor transmission rate number of layers OLED encapsulation`
3. `flexible encapsulation barrier WVTR versus number of dyads calcium test values`
4. `organic inorganic multilayer permeation barrier WVTR improvement per dyad`
5. `thin film encapsulation quantum dot OLED multilayer WVTR g/m2/day dyad`
6. `SiNx parylene multilayer barrier water vapor transmission rate dyads`

## 4. Inclusion criteria

A paper is included if all hold:

* it is a primary experimental report (not a review, not a patent);
* it reports water vapour transmission rates for **two or more** different dyad
  counts, layer counts, or repetitions of the same multilayer architecture;
* at least two of those values are given **numerically** in the abstract, body
  text, a table or a figure legend.

## 5. Exclusion criteria

* values available only as points on a plot;
* the two values differ in something other than the number of repeat units
  (for example a change of material or of deposition temperature);
* both values censored (reported only as "<" a limit).

Papers excluded for each reason are counted, because the count of papers whose
ladders exist only inside a figure is itself a reportable finding.

## 6. Recorded for every included paper

ladder (n : WVTR), measurement method, how the sample is sealed to the
instrument if stated, sensor area if stated, whether an impermeable control or
instrument floor is reported, whether any value is censored, temperature and
humidity, and the location of the numbers in the paper.

## 7. Analysis, fixed in advance

* **Q1** For every pair (n1 < n2) in every ladder, test G(n2)/G(n1) <= n1/n2.
  Primary statistic: the fraction of *papers* with at least one violating pair,
  with a Jeffreys binomial interval. Secondary: the fraction of pairs.
* **Q2** Among papers reporting a floor, the fraction whose best value lies
  within a factor of two of it.
* **Q3** The fraction of included papers that report a floor or control at all.
* Digitised values, where unavoidable, are flagged and the analysis is repeated
  with them excluded.

## 8. What each outcome would mean, fixed in advance

* violation rate below 15 % -- a cautionary note about specific measurements;
* violation rate 15-35 % -- a substantial minority of benchmark ladders carry a
  parallel channel; worth a specialist paper on its own;
* violation rate above 35 % -- a statement about the field, and the diagnostic
  becomes something referees should apply routinely.

Stating this now removes the temptation to decide afterwards which threshold
was the interesting one.
