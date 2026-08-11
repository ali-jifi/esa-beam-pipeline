# candidate labeling standard

Codified 2026-07-31 after two disagreement-review rounds, so the standard
survives the labeler. Applies to all `candidates/*.jsonl` records.

## class axes

Beam-ness (morphology) and population (context) are independent axes.
A morphologically perfect beam in the wrong population is a negative
(exhibit: `e_2013-08-24T0900_123020_17`, PSBL beamlet, R=32, organized
dispersion, still negative). Positive means cold ionospheric outflow,
not "looks like a beam."

## positive requires ALL of

1. **Discrete line** (clear spectral peak, Eb/dE well above ~1.5) OR an
   **organized chain terminating in a confirmed parent** (falling/rising
   energy ramp that connects, step by step, to a confirmed beam).
2. **Cold stack**: ne/ni >= ~1.2 and cold electrons (Te well under ~200 eV).
   Amplitude override (2026-07-31, formulation fixed 2026-08-02): strong
   morphology (flux_z >= ~10 or prom >= ~1) with cold Te overrides the
   ne/ni bar — sub-unity ne/ni on a spin snapshot co-occurs with
   unambiguous beams, so the bar keeps its force only in the
   weak-amplitude regime. The arms are evaluated as EPISODE-MAX, not
   single-sample: confirmed positives fluctuate a median 50% in flux_z
   (38% in prom) between adjacent samples, so hard single-sample constants
   are false precision. Isolated candidates have no episode to take a max
   over and cannot claim the override (case law: 214308 negative on
   exactly this). 005223 and 224646 are the logged test cases of the
   weak-regime bar for the 2015 scale-out — if the bar systematically
   kills real episode members there, revisit the bar with fresh evidence,
   not case-by-case mercy. Watchlist of episode-anchored sub-unity
   positives kept under the override: 090859, 082144, 095728, 130516,
   103800, 150937, 235219.

   **Regime conditioning (2026-08-03, closes the instrumental check)**:
   ne/ni is not noise (median 16% adjacent-sample variation vs 50% for
   flux_z; dips ANTIcorrelate with counting-floor flags) — it is a smooth
   regime variable measuring the hidden cold-ion fraction, riding beta
   (deep lobe median 1.87, plasma-sheet-adjacent 1.06). The bar is
   therefore regime-conditioned:
   - beta < 3e-3: bar in full force (ne/ni class auc .83 below the cut;
     positives run ~1.65 mid-band).
   - beta >= 3e-3: bar INAPPLICABLE (auc collapses .73 -> .57 between
     2e-3 and 3e-3, chance by 6e-3; sweep run with the adjudicated cases
     excluded). Population calls there rest on Te + episode + line.
   - Boundary frozen at beta = 3e-3, the first sampled cut at chance.
     Tie-break: beta within ~35% of the boundary (q75 of adjacent-sample
     beta variation) = boundary case, evaluated under BOTH band logics;
     if the outcomes differ, register. Precedent: 235219 (beta 4.0e-3,
     assigned transition band, survives mid-band treatment also).
   - ne/ni < ~0.5 at ANY beta = moment QC flag, not evidence (055405's
     0.09, 090859's 0.45) — extreme collapse is a different phenomenon
     from a marginal dip.
   - Evidence hierarchy, not mercy: mid-band sub-1.2 positives
     (145247/145609/025248) are held because episode membership + visual
     confirmation outrank one context feature — the codified ordering,
     distinct from the case-by-case mercy ruled out for 005223/224646.
   - FALSIFIABLE PREDICTION for the 2015 scale-out: mid-band sub-1.2
     positives should be RARE (positives there run ~1.65). If the
     scale-out produces many, either the band statistics or the
     2017-01-31 episode labels have a problem. CONFIRMED on the first
     peir batch (2026-08-03): 0/31 episode positives below 1.2, lowest
     exactly 1.20; watchlist gains no members.
   - Watchlist through the regime lens: the clause explains only
     103800 (transition) and 235219 (boundary). 090859 (0.45) falls to
     the QC threshold. 082144 (0.55), 095728, 130516, 150937 are
     mid/deep-band anomalies that remain the amplitude override's
     independent caseload — the override is NOT absorbed by the regime
     clause. Stage-2 instrumental check (moments-level, spacecraft
     potential) targeted at deep/mid-band extreme dips during scale-out.
   - Circularity exposure: same class as the flux_z floor — derived once
     (deferred cases excluded from the derivation), applied forward.
3. Adjacency alone NEVER suffices. Being minutes from a confirmed episode
   at a similar energy is supporting context, not evidence. No chain to
   the parent, no promotion (see borderlines 065940, 084444).

Episode membership (inside a confirmed episode window, matching energy,
confirmed neighbors on both sides) plus the cold stack can promote a faint
member (105909, 115827). Leading edges promote only with an organized
chain into the confirmed parent (115002, 012524 — low confidence).

**Isolation clause (2026-08-02)**: isolation is a confidence modifier,
not a kill rule. Isolated single-sample candidates that pass the full
standard label positive at LOW confidence, never medium+, and carry an
`isolated: true` record field so they can be excluded as a sensitivity
check in model A validation. The flag DROPS when a companion is later
confirmed (115410/115056 pair, de-isolated 2026-08-02 when the companion
labeled positive). Codified because the rule set contained no isolation
requirement while isolation was doing real intuitive work in the queue
reviews. Current flagged: 065940, 092349, 140154, 081324.

## negative kill rules

- **psbl_boundary_beam**: broad keV beam (E_b >~ 2 keV), earthward
  (dir_x_bx=+1), warm electrons (Te >~ 200-3000 eV), often off-lobe
  (beta > 0.1). Any two of direction/energy/Te suffice; ne/ni >= 1.2 does
  NOT rescue direction+energy (122343). A tailward hot broad beam is still
  negative on Te + Eb/dE (063627).
- **artifact / ratio-collapse**: single-spin (duration 1), ne/ni < ~1.2,
  floored cone or floored perp, no episode anchor on the day. The x5
  amplitude of an anti bar does not rescue ne/ni=0.67 (013535).
- **no discrete line**: Eb/dE ~ 1, prominence < ~0.3 without chain
  support. Cold stack does not rescue linelessness (071059, 102026).
- **warm electrons**: Te >~ 200 eV kills even tailward in-lobe candidates
  unless the record is chain-connected to a cold confirmed parent.
- **flux_z floor** (derived 2026-07-31 from the labeled distributions, not
  from the cases it adjudicates): positive and artifact flux_z interleave
  through 2-5 with no clean gap, so flux_z alone cannot adjudicate. Rule:
  flux_z < 3.5 (p5 of labeled positives) -> negative unless chain-anchored
  or episode-sandwiched (confirmed positives on both sides at episode
  energy — stronger evidence than a chain; codified 2026-08-02, keeps
  112254 and 105909). Applies FORWARD only (register sweeps, future
  batches); resolved labels are frozen
  (candidates/labels_frozen_20260731.csv) and never relabeled by rules
  derived from them.

## excluded third class: psbl_embedded_cold_line

FROZEN DEFINITION (written before any targeted search, 2026-08-02, so
the search cannot quietly redraw the boundary): narrow <=~250 eV para
line, R >~ 20, Te >~ 1 keV, off-lobe or beta > 0.05. Morphologically
perfect, population-ambiguous — plausibly cold outflow threading the
PSBL flank. Neither positive nor psbl_boundary_beam: label null,
`class: psbl_embedded_cold_line` record field, excluded from training
(user decision 2026-07-31). The class field is a CURATION field, never a
model feature — leaking it into model A would rebuild the model B
circularity.

Members: 123020 (reclassified from negative, prior label preserved in
the freeze), 125854, 130057, 132729, 124235, 130259 (all e_2013-08-24).

Case law (2026-08-02, first targeted encounter with the boundary):
- 124235 ADMITTED with R=13.5, below the ~20 in the definition — core
  criteria (cold narrow line in hot off-lobe plasma) carry; the R
  relaxation is explicit precedent.
- 105424 EXCLUDED on three edges (336 eV > 250, Te=772 < 1 keV,
  beta=0.043 < 0.05, in-lobe) despite a strong tailward line — the
  frozen boundary is not stretched; held on register for
  class-definition review.
- Candidate member to adjudicate: 085502 (tailward, ne/ni=1.9,
  Te=1.5 keV, fits neither class).

## catalog post-filter

The deployed catalog (scripts/score_candidates.py) hard-drops
dir_x_bx=+1 (earthward) from both tiers regardless of score — the shared
blind spot of both classifiers. 055405, the one labeled positive it cost,
was relabeled negative 2026-08-02 (pre-standard label; fails the current
standard on energy and Te; its ne/ni=0.09 is a moment QC flag, not
evidence). None (no B data) passes.

The catalog also carries a `class` column (2026-08-02): class-tagged
records (psbl_embedded_cold_line) are excluded from tiers 1/2 and listed
in their own block — preserving the distinction between "ruled out"
(earthward hard drop, a known different thing) and "not yet ruled" (an
unresolved population worth a targeted pass).

Physics note (advisor, 2026-08-02): earthward cold beams CAN be genuine
ionospheric outflow — tailward ions may mirror near the opposite
hemisphere and bounce, giving simultaneous earthward+tailward beams. The
filter is a scope restriction (clean tailward-only sample, avoids
earthward ambiguity with PSBL return flow), not a claim that earthward
cold beams are unphysical. Frame it that way in the paper.

## borderline convention

Unresolvable records stay `label: null` with a `label_note` explaining the
tension. They are excluded from training and from future review queues.

Register as of 2026-08-02 post-cutout-review: **105424 alone, by
design** — the primary test case for the 2015 coldline boundary review.
085502 resolved negative (keV energy dispositive per 122028/055405
despite ne/ni=1.9 and dur=10). Standing entries resolved-in-place:
071247 (permanent), 012337 (below floor but chain-anchored, exemption
holds).

Consistency audit 2026-08-02 (uniform standard over all pre-standard
positives): 12 flagged, resolved as — 3 relabeled negative on cutout
review (142819 contamination-floor + Te=554; 114715 two-population
structure conflated into one record; 101132 broad step below floor);
2 kept via the episode-sandwich floor exemption (112254, 105909);
7 weak-regime ne/ni cases (103800, 235219, 103116, 080959, 145247,
145609, 025248) resolved 2026-08-03 by the regime clause — all keep
their positive labels with regime annotations (4 transition-band incl
the 235219 boundary precedent, 3 mid-band on evidence hierarchy). The
instrumental check is CLOSED at stage 1; stage 2 (moments-level, for
deep/mid-band extreme dips) rides along with the 2015 scale-out.

Sweep history: 2026-08-02 register sweep resolved 7 negative (045426,
102026, 145024, 145708 by floor, 103455, 134429, 214308 by episode-max
isolation) and promoted 2 positive (065940 under the isolation clause;
084444 on chain evidence completed by the 085917 promotion). Earlier:
130957 (2.0) and 095150 (2.4) swept negative by the flux_z-floor rule;
145553 resolved negative (earthward per JSONL); 123020 moved to the
psbl_embedded_cold_line class; 145024 note — tailward confirmed by
audit, did not flip, died on ne/ni.

## record fields

- `label`: positive | negative | null
- `label_source`: auto_rule | manual_unanimous | manual_promoted_review |
  manual_disagreement_review | manual_borderline_resolution |
  flux_z_floor_rule | manual_catalog_review | manual_class_adjudication |
  manual_episode_review; suffix `+chain_evidence` flags labels whose
  evidence included chain features (circularity audit for model A).
  manual_episode_review (peir era, 2026-08-03): verdicts are given at
  EPISODE level and propagate to every member sample; `group=<episode_id>`
  in the note ties members together.
- `class`: curation-only population tag (psbl_embedded_cold_line);
  excluded from tiers and NEVER a model feature.
- `isolated`: true on positives promoted under the isolation clause;
  excludable as a sensitivity check in model A validation.
- `label_confidence`: high | medium | low
- `label_note`: free text; `group=<event-id>` tags episode membership for
  positives promoted via episode/chain evidence. CV grouping is probe+date
  (automatic from candidate_id), which already isolates parents with their
  extensions and episode members.

## tail2015 / peir episode era (2026-08-03)

First peir batch: 48 episodes reviewed — 31 positive, 11 negative, 6
borderline (tail2015 register: a_120457, d_125737, a_120851, e_125829,
a_125552, a_121009; r2 added d_2015-12-21_124615, a_2015-12-20_114541,
d_2015-12-30_133538, d_2015-12-21_125056 floor-protected,
e_2015-12-24_123747; r3 added a_2015-12-30_132342, the best
cold-on-negative-day case). Constraints codified from that review, binding for
the episode-model phase:

- **CV splits by DAY, not probe or episode** — multi-probe same-day
  blocks are one physical event (2015-12-31 a/d/e; 11-03 four correlated
  negatives from one crossing). Effective sample size of the 48 is ~10
  events. Implemented in train_model.load_labeled.
- **Episode feature tables aggregate over members, never rep samples** —
  five confirmed positives have rep rows with strict-gate FAILs; fitting
  on reps would teach the model to reject valid beams.
- **NaN flux_z policy is explicit** — episodes with no defined flux_z
  (slow-survey/starved baselines) carry a missing indicator; never
  silently imputed for episode fits.
- **This batch is a sanity set, not a generalization test** — Te
  separates pos/neg perfectly (<=110 vs >=120, no overlap) so model B is
  a ceiling artifact by construction; flux_z_max alone is also linearly
  separating (>=21.1 vs <=16.2). Expect near-perfect AUCs that will not
  survive harder episodes.
- New negative note vocabulary: plasma_sheet_anisotropy (warm broad
  para excess in ps-adjacent plasma, the peir-era analog of
  psbl_boundary_beam at lower energy).

Round 2 (2026-08-03, 48 adversarial episodes: 27 pos / 16 neg / 5
borderline) — structural outcomes:
- **THE DAY CONFOUND**: only 2 of 18 labeled days contain both classes;
  the label is nearly a function of the date. Any feature with more
  between-day than within-day variance looks discriminating, and
  leave-one-day-out CV cannot protect (16/18 folds are single-class).
  Episode-model AUCs are optimistic until round 3 (within-day
  stratification: negatives hunted on positive days and vice versa)
  either breaks the confound or proves the model is a day classifier.
- **Duration features poisoned in this label set**: r2 stratification
  paired long hot negatives against short cold positives, so n_samples /
  duration coefficients read sampling, not physics (fitted -1.80 on
  log n). Dropped from episode model A.
- flux_z and prominence got honest (AUC 1.000/0.990 -> 0.859/0.837
  pooled) and are the realistic ceiling; Te separation (0.982) is real
  physics but stays out of model A; model B stays an unreported
  artifact.
- Sub-unity ne/ni: two batches, zero sub-unity positives (6 more in r2,
  all neg/borderline). Sub-unity is established negative-side
  transition-band signal; watchlist stays empty.
- Amplitude override: first three confirmed fires, all beta < 4e-3.
- Te is breakable but barely: stratum A (Te 100-250) produced exactly
  one hot positive (125729, Te 207 + beta 8e-4 + discrete line).
- Constant-protection precedent: 125056 left borderline BECAUSE it sits
  on the floor constant (fz 3.9) — labeling cases that sit on a
  pending-re-derivation constant biases the re-derivation.

Round 3 (2026-08-04, 48 within-day episodes: 0 pos / 47 neg / 1
borderline) — THE CONFOUND RESOLVED:
- **It was beta, not the date.** Zero of 58 positives sit above beta
  0.0272; beta < 0.03 captures 58/58 positives admitting 16/74
  negatives (AUC 0.991 alone). The negative days were days the
  spacecraft never entered the lobe. The genuine null: 24 positive-hunts
  on negative days returned zero (one borderline, 132342, on register).
  The episode model is a LOBE DETECTOR plus a weak morphology term —
  frame it that way everywhere.
- **The corridor population**: 10 narrow tailward strong-amplitude
  spikes in warm/hot ambient (fz up to 178.6, prom to 1.83), labeled
  negative because a deployed detector fires on them. Cost accepted as
  the true measurement: flux_z is DEAD as a discriminator (negative
  class now out-fluxes every positive; refit coefficient went negative);
  prominence is the only detection-only feature left carrying weight.
  Honest episode-model number: ~0.80-0.82 day-blocked.
- Duration features: permanently dropped, never audited (r3 negatives
  run to 123 samples).
- Next-batch rule: stratify on BETA WITHIN THE LOBE (beta < 0.03),
  never on days (deck E).

Final sitting, decks C/D/E (2026-08-11) — CLOSEOUT:
- Deck C (32 disagreements: 31 neg, 1 borderline -> register): the two
  strata were two defects, not disagreements. Model-high failures are
  direction/energy blindness — **direction and the E_b window stay HARD
  GATES upstream, never features** (unlearnable on 72 positives).
  Rules-high failures confirm the project thesis: the strict AND-chain
  fires on warm high-beta noise a fitted model rejects. New named
  excluded population: **earthward keV-beam family** (cold electrons,
  low-mod beta, very high flux_z, E_b 4-7 keV, earthward — reconnection
  outflow / PSBL earthward beam; one line in the paper).
- Deck D (16): **coldline recurrence CONFIRMED** — 8 new members on 7
  distinct non-12-31 days, probes a/d/e. One-event objection closed.
  44% boundary rate -> TWO MORE CONSTANTS DECLARED: Eb/dE >= ~1.0 and
  prom >= ~0.7, plus R tightened to >= 25. Six-constant definition
  (E_b<=300, R>=25, Te>=1keV, beta>0.05/off-lobe, Eb/dE>=1.0,
  prom>=0.7) yields the FINAL CENSUS: 41 members over 13 days
  (coldline_final.csv). 2015-12-08 probe-d excluded at DAY level
  (both candidates ne/ni 0.40/0.69 — moment failure, not physics).
- Deck E (35: 14 pos, 21 neg): morphology separates inside the lobe
  (strong 53% vs weak 25% positive) but beta dominates (78% -> 0%
  across two decades). The cleanly separating quantity is SIGNED
  asym_dev: every positive in the deck has asym_dev > 0; five strong
  negatives are anti-parallel at asym_dev -0.80..-0.86. Added to the
  episode feature vector; it became the TOP coefficient (+1.80) and
  lifted the lobe-scoped model to AUC 0.956 / f1 .927 (109 lobe
  episodes, 20 day-folds). Report lobe-scoped numbers ONLY — a global
  AUC is a lobe detector's score. Context-morphology covary below
  beta 1e-3 and decouple above 1e-2; d_2015-12-31_131202 (ne/ni 4.63,
  no beam) is the canonical counterexample to context-only reasoning.
  Four TE breakers now exist (Te 207/236/274/276, all deep-lobe with
  visible lines) — the Te kill keeps its chain/episode exemption arm.
- Floors re-derived for the peir era (forward-only, episode level):
  flux_z floor 6.5 (p5 of 72 positives; old 3.5 separated nothing —
  positives bottom at 4.8, floor-region negatives at 2.7/3.1);
  prominence floor 0.55 (observed positive minimum; tight — negatives
  reach 0.43 nearby, apply with episode-anchor exemption only).
- 2015 catalog assembled (beam_catalog_2015.csv): hard gates upstream
  (tailward, E_b<600, beta<0.03) then model x prominence two-tier:
  tier 1 = 150 episodes (70 lpos / 6 lneg / 74 unlabeled), tier 2 =
  115. Spot-check queue: 24 unlabeled tier-1 over 14 days — the final
  review of the arc.

Spot-check + arc closure (2026-08-11, 24 episodes: 8 pos / 14 neg / 2
borderline):
- **Measured contamination**: 33% at threshold 0.5 (4/12 decided, CI
  14-61%), 0% at 0.85+ (0/8, CI 0-32%). Threshold moved to 0.82 —
  costs nothing in this sample (lowest positive p=0.881, highest
  negative p=0.764). Catalog recut with it plus a Te < 500 eV scope arm
  at admission (TE breakers max at 276, corridor leaks start at 596 —
  clean daylight): **final tier 1 = 89 episodes (53 lpos / 0 lneg / 36
  unlabeled), tier 2 = 142.**
- **BLIND CAVEAT (binding)**: spot-check verdicts were NOT blind to
  model probabilities (prob column visible during review). The perfect
  verdict/prob separation is NOT validation and must never be quoted.
  One re-run deck with probs withheld is REQUIRED before any separation
  number goes near the paper. Contamination rates are less affected
  (driven by ne/ni and Te) but inherit the caveat.
- **The prediction test — the methods-section sentence**: 15/15
  context-only pre-verdicts correct; morphology changed exactly 1
  verdict in 24 (105418, the fifth TE breaker, carried by R 49.5 +
  cross-probe energy match). The detector is a REGIME CLASSIFIER WITH A
  THIN MORPHOLOGY CORRECTION: regime does the work, morphology
  resolves the residual.
- Scoping leak traced and closed: beta < 0.03 alone admits warm
  corridor plasma because beta and Te decouple in the transition region
  (deck E's own result); 5 tier-1 episodes had te_med > 500, all fall
  out at 0.82, and the Te scope arm closes the class at admission.
- Intra-rater consistency: d_2015-12-21_124615 re-presented blind to
  its round-2 verdict, same call both times, same reasoning. n=1 —
  state as a limitation, not a strength.
- Multi-probe events for statistics: 12-31 08:33 triple conjunction
  (d 171 eV / e 295 / a 454 within 64 s — energy ordering worth a
  dispersion analysis) and the 12-20 11:25-11:40 five-episode
  three-probe cluster each count as ONE event.
- Isolation clause calibrated: single-spin 3/13 positive vs
  multi-sample 5/11, all singles at low/medium confidence.
- Register adds: e_2015-11-13_141929; 124615 re-confirmed.

Coldline 2015 census adjudication (2026-08-03): 10 confirmed + 5
boundary + 1 not_coldline of the top-16 strict; 5 strict hits pulled on
dir_x_bx=+1 (scope rule applies to class membership). CORRECTIONS OF
RECORD: the class core is ONE EVENT (the 2015-12-31 d/e conjunction —
R>=50 & prom>=1 leaves 7 candidates, all in one 25-min window); the
beta 0.05 / Te 800 / R floors are FILTER CONSTANTS not observed edges
(0/911 hits below beta 0.051 by construction); the corridor beta
0.027-0.051 between lobe positives and coldline hits contains only
plasma_sheet_anisotropy negatives. 108 hits = 108 filter passes, not a
population, until the negative control runs (corridor candidates the
filter REJECTS at beta 0.05-0.2, Te>1000, checked for visual
distinctness). 105424 stays HELD — it is the type specimen by
construction and cannot also evidence the class. Paper status: the
12-31 conjunction is a defensible case study; no "recurring population"
claim without the control and a multi-day core.

Control outcome (2026-08-04): the deck was not a fair control (it
varied E_b and direction alongside breadth/R — 8/16 earthward, 14/16
above 300 eV). What it proved anyway: the class SURVIVES against broad
low-R humps at comparable energy (R gate earns its keep at the bottom)
and DOES NOT survive against narrow high-R lines at 581-1010 eV, which
are morphologically indistinguishable from the core and excluded by an
energy gate that was never declared. **FOURTH CONSTANT DECLARED:
E_b <= ~300 eV** (strict pool tops at 175.3, wide at 303.6), justified
physically: ionospheric ions accelerated by polar-cap potentials arrive
at tens to a few hundred eV; narrow keV parallel lines in the plasma
sheet are more plausibly reconnection outflow or bounce-resonant
populations. The 581-1010 eV narrow lines are a morphologically
identical but DISTINCT population, reported separately. The honest
claim: morphology does not establish ionospheric origin — energy does.
Direct evidence the gates cut a continuous structure: the 12-31 14:37
three-probe conjunction is split by the beta gate (d at beta 0.0265 /
Te 956 falls below both gates while e and a sit above); all three
labeled corridor negatives.

## mechanics

Labels live inside the jsonls. Before any batch regen: back up to
`candidates/pre_regen/`, regen, then `scripts/migrate_labels.py` (joins by
candidate_id, falls back probe+t_center rounded 0.1 s, carries notes).
Verify orphan count and that gold auto-rule counts reproduce before
deleting backups.
