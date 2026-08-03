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
   103800, 150937, 235219. Open check: do ne/ni dips anticorrelate with
   anything instrumental.
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
145609, 025248) deferred pending the ne/ni instrumental check — the
only open analysis gating them.

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
  flux_z_floor_rule | manual_catalog_review | manual_class_adjudication;
  suffix `+chain_evidence` flags labels whose evidence included chain
  features (circularity audit for model A).
- `class`: curation-only population tag (psbl_embedded_cold_line);
  excluded from tiers and NEVER a model feature.
- `isolated`: true on positives promoted under the isolation clause;
  excludable as a sensitivity check in model A validation.
- `label_confidence`: high | medium | low
- `label_note`: free text; `group=<event-id>` tags episode membership for
  positives promoted via episode/chain evidence. CV grouping is probe+date
  (automatic from candidate_id), which already isolates parents with their
  extensions and episode members.

## mechanics

Labels live inside the jsonls. Before any batch regen: back up to
`candidates/pre_regen/`, regen, then `scripts/migrate_labels.py` (joins by
candidate_id, falls back probe+t_center rounded 0.1 s, carries notes).
Verify orphan count and that gold auto-rule counts reproduce before
deleting backups.
