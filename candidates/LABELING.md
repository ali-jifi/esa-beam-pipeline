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
   Amplitude override (2026-07-31): strong morphology (flux_z >= ~10 or
   prom >= ~1) with cold Te overrides the ne/ni bar — sub-unity ne/ni on a
   3 s spin co-occurs with unambiguous beams, so the bar keeps its force
   only in the weak-amplitude regime. Watchlist of episode-anchored
   sub-unity positives kept under the override: 090859, 082144, 095728,
   130516, 103800, 150937, 235219. Open check: do ne/ni dips anticorrelate
   with anything instrumental.
3. Adjacency alone NEVER suffices. Being minutes from a confirmed episode
   at a similar energy is supporting context, not evidence. No chain to
   the parent, no promotion (see borderlines 065940, 084444).

Episode membership (inside a confirmed episode window, matching energy,
confirmed neighbors on both sides) plus the cold stack can promote a faint
member (105909, 115827). Leading edges promote only with an organized
chain into the confirmed parent (115002, 012524 — low confidence).

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
  flux_z < 3.5 (p5 of labeled positives) -> negative unless chain-anchored.
  Applies FORWARD only (register sweeps, future batches); resolved labels
  are frozen (candidates/labels_frozen_20260731.csv) and never relabeled
  by rules derived from them.

## excluded third class: psbl_embedded_cold_line

Intense narrow 100-250 eV para lines (R in the tens, flux_z up to ~50)
embedded in Te ~1.5 keV off-lobe plasma. Morphologically perfect,
population-ambiguous — plausibly cold outflow threading the PSBL flank.
Neither positive nor psbl_boundary_beam: label null, class tag in the
note, excluded from training (user decision 2026-07-31). Members: 123020
(reclassified from negative, prior label preserved in the freeze),
125854, 130057, 132729 (all e_2013-08-24). Candidate member to
adjudicate: 085502 (tailward, ne/ni=1.9, Te=1.5 keV, fits neither class).

## catalog post-filter

The deployed catalog (scripts/score_candidates.py) hard-drops
dir_x_bx=+1 (earthward) from both tiers regardless of score — the shared
blind spot of both classifiers; costs one labeled positive (055405, a
pre-standard label that contradicts the codified rules on four axes,
flagged for relabel). None (no B data) passes.

## borderline convention

Unresolvable records stay `label: null` with a `label_note` explaining the
tension. They are excluded from training and from future review queues.
Register as of 2026-07-31 post-tier-1-review (12): 071247 (permanent),
085502 (candidate psbl_embedded_cold_line member), 012337 (below floor
but chain-anchored, exemption holds), 045426, 084444, 065940, 102026,
145024 (tailward confirmed, does not flip), 145708, and from the tier-1
review: 134429, 214308, 103455. Most resolve mechanically under the
fixed floor (3.5) and the ne/ni amendment — pending one sweep pass.
Swept negative by the flux_z-floor rule: 130957 (2.0), 095150 (2.4).
Resolved off register earlier: 145553 (negative, earthward dir_x_bx=+1
per JSONL). 123020 moved to the psbl_embedded_cold_line class.

## record fields

- `label`: positive | negative | null
- `label_source`: auto_rule | manual_unanimous | manual_promoted_review |
  manual_disagreement_review | manual_borderline_resolution |
  flux_z_floor_rule | manual_catalog_review; suffix `+chain_evidence`
  flags labels whose evidence included chain features (circularity audit
  for model A).
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
