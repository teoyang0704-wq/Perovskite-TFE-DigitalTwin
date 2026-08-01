# Year-2 discovery program: what closure actually does

Written 2026-08-01, after attempting the nucleation-closure correlation and
failing for a reason worth recording.

---

## 1. The attempt, and why it failed

The plan was to test whether the closure length inferred from permeation is
predictable from ALD growth parameters -- nucleation delay, growth per cycle,
deposition temperature -- which would let barrier quality be forecast from
growth data alone.

It cannot be tested with what we have. Of six films in the hierarchical
population, only **two** have closure lengths that the data actually identify:

| film | 90 % width of log10 d_close | thickness points below 20 nm | verdict |
|---|---|---|---|
| Wu high quality | 0.16 dex | 2 (15, 20 nm) | identified, 1.20 nm |
| Groner | 0.30 dex | 3 (2.5, 5, 10 nm) | identified, 0.61 nm |
| Buelow | 1.58 dex | 0 | prior-dominated |
| Carcia | 1.10 dex | 0 | prior-dominated |
| Lee | 1.05 dex | 0 | prior-dominated |
| Wu poor | 0.39 dex | 0 | degenerate with its floor, not meaningful |

The lesson is sharper than the failure. **A thickness series only constrains
closure if it has points below about 20 nm.** Buelow measured three thicknesses
and none of them informs closure, because 25, 50 and 100 nm all sit after the
nucleation term has died. Most published thickness series are in that regime,
which is why this quantity has never been pinned down: the field measures where
the barrier is good, and closure happens where the barrier is bad.

Two points do hint at a temperature trend -- Groner at 120 C gives 0.61 nm,
Wu at 60 C gives 1.20 nm, in the direction expected if warmer deposition drives
more trimethylaluminium into the polymer subsurface and seeds a denser
nucleation population -- but two points with different substrates and different
test methods is an anecdote, not evidence.

**Data needed to revive this:** five to eight films with at least two thickness
points below 20 nm, quantitative permeation, and a stated deposition
temperature.

---

## 2. The better question the attempt uncovered

Our model writes the defect area fraction as

    f(d) = f0 exp(-d / d_close) + f_res

and treats it as one number. But the area fraction is a product of two
independent things:

    f = N * pi * r^2         N = defect number density,  r = defect radius

**Nobody has established which of the two the exponential decay belongs to.**
Two mechanisms are physically distinct and predict identical permeation:

* **elimination** -- nucleation islands merge and defects vanish one by one.
  N decays exponentially, r stays roughly constant.
* **constriction** -- defects survive but narrow as material closes in from
  their rims. N stays roughly constant, r shrinks.

Both give f ~ exp(-d/d_close), so permeation alone cannot separate them. That is
why the question has stayed open: the standard observable is degenerate.

### The diagnostic

It stops being degenerate the moment somebody counts defects on the same films.
Writing the decay lengths of the two observables separately,

    WVTR(d) ~ exp(-d / L_flux)          (permeation)
    N(d)    ~ exp(-d / L_count)         (decorated defect counts)

then, since f = N pi r^2,

    r(d)^2 ~ exp[ -d ( 1/L_flux - 1/L_count ) ]

so that

    L_flux = L_count      pure elimination, radii unchanged
    L_flux < L_count      constriction as well; radii shrink with thickness
    L_flux > L_count      would mean radii *grow*; would falsify the picture

**A single number -- the ratio L_count / L_flux -- decides between two
mechanisms of thin-film closure.** This is measurable from data that already
exists.

### Why it matters beyond bookkeeping

If closure works by constriction, the surviving defects are getting narrower,
and a molecule stops fitting through them at a thickness that depends on the
molecule. The critical thickness of a barrier would then be a property of the
**permeant**, not of the film:

    d_critical(sigma) decreases as the kinetic diameter sigma increases

so oxygen (3.46 A) would appear to close earlier than water (2.65 A) on the very
same film, and helium (2.6 A) later than both. If closure works by elimination,
all permeants close together and the critical thickness is a film property.

There is a suggestive hint already. Groner et al. report, for one set of films,
oxygen below the coulometric detection limit at 5 nm while water is still at
0.1 g m^-2 day^-1 at the same thickness and does not reach its floor until about
26 nm. That is the direction constriction predicts. It is not evidence, because
the oxygen number is censored -- below a detection limit is not a measurement --
but it is exactly the comparison that would settle it.

### Why this is a discovery-shaped question and the earlier ones were not

The edge channel was a correction to how the field reads its instruments. The
series bound is a theorem about circuits. The hierarchical calibration is
better statistics. All three are useful and none is a new natural phenomenon.

"Does a growing film close its defects by removing them or by narrowing them,
and does that make the barrier's critical thickness a property of the gas rather
than the film?" is a question about a physical process, it has two mutually
exclusive answers, it is decidable with existing published data, and its answer
changes what a barrier specification means. If constriction wins, then every
critical thickness in the literature is quoted for one permeant and is not
transferable to another -- and multi-permeant permeation becomes a way to
measure the defect size distribution of a film without imaging it.

---

## 3. What data would decide it -- checked before asking for any

The diagnostic was written first and its statistical power measured on synthetic
data, with the particulate floor included and the decay length fitted jointly
with the floor. That check overturned the obvious plan.

| thickness series | pure constriction found | partial constriction found | pure elimination **misread as constriction** |
|---|---|---|---|
| 15-100 nm (six points) | 97 % | **0 %** | **78 %** |
| 2.5-26 nm (four points) | 98 % | 49 % | 1 % |
| 2-30 nm (six points) | 97 % | **84 %** | 1 % |

Two things follow, and the first is a warning.

**A series that starts at 15 nm cannot answer this question, and worse, it
answers it wrongly.** Once every measured thickness sits in the floor-dominated
regime, the defect counts look flat whether or not they are, and the test calls
constriction on a film that is purely eliminating in about four cases out of
five. Klumbies et al. measured defect counts and permeation on the same films --
the only such dataset found -- but from 15 to 100 nm. Fetching it and running the
diagnostic would have produced a confident, publishable, wrong answer.

**The decisive measurement is well defined and appears not to exist**: six
thicknesses spanning roughly 2 to 30 nm, on one deposition campaign, with defect
density counted by electrodeposition decoration or equivalent *and* permeation
measured on the same samples. Groner's series has the right thickness range but
no defect counts; Klumbies has the counts but the wrong range. Nobody has put
the two together, because a 2 nm alumina film is not a barrier and nobody has had
a reason to characterise one carefully.

That is a modest experiment -- six ALD runs and two characterisations -- and it
is the single most valuable thing a laboratory collaboration could contribute to
this project. It is also the right thing to ask for, in preference to anything
that merely extends the model.

### Still worth obtaining, for other reasons

1. **Klumbies et al., Organic Electronics 15 (2014) 3242.** Not for the
   mechanism test, but because their AFM result -- 5 nm films leaking at flat,
   particle-free sites -- is direct microscopic support for splitting the defect
   fraction into a nucleation term and a particulate floor, which our model
   assumed without evidence.
2. **Any paper with WVTR and OTR on the same films at three or more
   thicknesses**, neither value censored. This tests the permeant-size
   corollary independently of the counting route.
3. **Thickness series below 20 nm with a stated deposition temperature**, to
   revive section 1 -- currently blocked at two identified films.

## 4. Honest assessment

The elimination-versus-constriction question is the best discovery candidate
this project has produced, and it is genuinely open. It is also modest in scope:
answering it is a strong specialist paper, not a Nature paper, because it
resolves a mechanism rather than opening a field. The permeant-size corollary is
the part with real reach -- if critical thickness turns out to be a property of
the gas, every barrier specification in flexible electronics is quoted against
an unstated variable.

The risk has moved. It is no longer that the answer might be negative -- a clean
"elimination" would close a real question and is worth publishing. It is that
the decisive data does not exist and cannot be manufactured from literature, so
this line stalls at a well-posed question with a designed experiment attached
unless a laboratory runs it.

What the attempt has already produced, without any new data: a diagnostic whose
power is known in advance, a demonstration that the obvious dataset would give
the wrong answer four times in five, and a specification of the six-point
experiment that would settle it. Knowing which measurement to ask for, and
knowing which one would have misled us, is the useful output of a failed
search.
