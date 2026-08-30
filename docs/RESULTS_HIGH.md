# High-layer results — validation, methods findings, and discoveries

The interpretive results built on the catalog described in RESULTS_SURFACE.md. Quotability rules at the end are binding.

## 1. Blind validation (the credibility result)

A 24-episode deck was adjudicated with model probabilities and tiers withheld (key sealed). After probe + time-overlap dedup — episode_id dedup is insufficient because episode boundaries are pipeline artifacts — the deck contained 20 independent structures.

| quantity | value |
|---|---|
| blind verdict-vs-probability AUC | **0.969** |
| pooled contamination at thr 0.82 | **1/14 = 7% (95% CI 0–34%)** |
| day structure | reproduced blind — positives landed only on 12-20 / 12-31 of 14 days sampled |

The threshold itself was set earlier by a spot-check (contamination 33% at thr 0.5 → 0% at 0.82); that check saw model scores, so only the blind numbers are quotable for separation.

## 2. The central methods finding

**The detector is a plasma-regime classifier with a thin morphology correction.** Context-only information predicted 15/15 review verdicts; morphology changed exactly 1 verdict in 24. The honest counterexample is d_131202 (ne/ni 4.63, context screams beam, no beam) — context is dominant, not infallible.

Supporting structure:

- **Hard-gates doctrine**: direction (tailward) and E_b window are upstream gates, never model features — 11/16 of one review round's model false positives were pure direction failures the gates now remove.
- **Lobe-scoped model: AUC 0.956 / f1 0.927** (109 lobe episodes, 20 day-folds). Signed asymmetry deviation (asym_dev) is the top feature (+1.80); anti-parallel negatives sit at −0.8. Pooled (non-lobe-scoped) numbers are regime-inflated and are not reported.
- **flux_z demoted to feature-only**: the blind deck exposed baseline absorption — the running baseline eats persistent beams (factor ~50 on one continuous structure: fz 0.64 vs 30.1, same spectral line, 90 s apart). Explains its AUC decay 1.0 → 0.79. Prominence + R carry all morphology gating.
- **Named failure mode — event-boundary episodes**: morphology outlives the plasma regime at structure exits. The type case, e_2015-12-31_083709 (p = 0.837), is the 08:33 conjunction's exit and is the single labeled false positive in tier 1.

## 3. Multi-spacecraft events (the physics showcase)

Nine physical structures (see RESULTS_SURFACE.md), six multi-probe, three simultaneous on all three probes:

- **12-31 08:33 triple conjunction**: 5 members (4 tier-1 + 1 tier-2) across a/d/e, core three within 64 s, total duration bounded at ~4–6 min. Catalog beam energies are 171 (d) / 295 (e) / 454 (a) eV — but those are three *different times* spanning that 64 s, and each probe's energy moves within the event: probe a also detects 175 eV at 08:33:01 (sub-threshold, not in catalog) and probe e runs 295 → 57 → 441 eV. The spread is confounded between probe position and time and is **not** a per-probe ordering. **Cross-probe time-of-flight is not resolvable at survey cadence** (corrected 2026-08-17): measured member separations are ΔX 0.34, ΔY 0.06, ΔZ 0.97 RE, so the field-aligned baseline (ΔX, since B ≈ ±X in the lobe) implies only ≈10 s of arrival spread against a ~100 s peif sampling interval. Dispersion analysis requires peir (3.2 s). Its exit episode is the catalog's known false positive — the event-boundary failure mode illustrated inside the showcase event.
- **12-20 11:21 structure**: 62 min, 39 episodes, all three probes — the Artemyev+ 2020 Figure 1 day, recovered independently at 30× time resolution.
- **12-31 11:59 structure**: 71 min, 25 episodes, all three probes — the persistent structure that exposed the flux_z baseline-absorption failure.

Occurrence statistics count each conjunction/cluster as one event.

## 4. The coldline population (secondary discovery)

A distinct recurring population: narrow cold field-aligned ion lines embedded in hot (Te ≥ 1 keV) boundary-layer plasma. **41 members over 13 days** (coldline_final.csv) under a six-constant frozen definition: E_b ≤ 300 eV, R ≥ 25, Te ≥ 1 keV, beta > 0.05, E_b/dE ≥ 1.0, prom ≥ 0.7.

| quantity | value |
|---|---|
| census | 41 members / 13 days (Nov 13 – Dec 31) |
| concentration | 12-31 holds 20 members (a/d/e, incl. the 25-min d/e conjunction = class core); 12-07 has a 3-probe day |
| properties | E_b median 171 eV (33–296), Te median 1.4 keV, R median 34 |
| controls passed | negative control (corridor rejects look different) + recurrence (8 members on 7 non-12-31 days) |

Discipline: the type specimen (e_2013-08-24 123020) is held out and never counted as class evidence; the morphologically identical 581–1010 eV narrow-line population is reported separately — morphology does not establish ionospheric origin, energy does. Origin unresolved (plausibly cold outflow threading the PSBL flank); no mechanism claim.

## 5. Excluded populations (scope boundary, one line each)

1. **psbl_boundary_beam** — broad hot keV beams (E_b ≳ 2 keV, Te 200–500 eV, beta 0.01–0.3): boundary-layer, not outflow.
2. **earthward keV-beam family** — E_b 4–7 keV, cold electrons, very high flux_z, earthward: reconnection outflow / PSBL return.
3. **corridor population** — narrow tailward strong-amplitude spikes in warm ambient (beta 0.027–0.051): plasma-sheet transition region.
4. **plasma_sheet_anisotropy** — warm broad parallel enhancements: anisotropy, not beams.

The earthward exclusion is a scope restriction, not physics — bounce returns mean earthward cold beams can be genuine outflow (tailward-only gives a clean sample at the cost of the bounce population).

## 6. Quotability rules (binding)

- Separation/contamination numbers: blind deck only (AUC 0.969, 7%). Spot-check numbers appear only as the threshold-setting narrative.
- Model performance: lobe-scoped only (AUC 0.956). Never pooled.
- Occurrence: physical events (9), not episodes (89), for any rate-of-nature claim; episodes are pipeline segments.
- Coldline: census + case-study framing; no origin claim; type specimen excluded from evidence.
- Limitations to state: single labeler, intra-rater consistency n = 1, tailward-only scope, 2015 Oct–Dec coverage (verified complete for 2015).
