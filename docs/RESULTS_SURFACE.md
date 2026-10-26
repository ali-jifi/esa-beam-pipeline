# Surface-layer results — 2015 catalog statistics

Descriptive statistics of the 2015 beam catalog (candidates/tail2015/analysis/beam_catalog_2015.csv). All numbers computed directly from the catalog files 2026-09-19, after the catalog was regenerated with the direction-invariant episode model (dominant-cone features, refit on the unchanged 2015 labels; see LABELING.md refit addendum). The 2026-08-11 catalog (89 / 142) is kept as beam_catalog_2015_frozen20260811.csv. Interpretation, validation, and physics-level results are in RESULTS_HIGH.md.

## Input coverage

| quantity | value |
|---|---|
| survey scope | all of 2015; tail region x_GSM < −8 Re, \|y\| < \|x\| |
| qualifying coverage | mid-Oct through Dec only (365-day scan verified; Jan–Sep has zero qualifying intervals) |
| intervals / probe-hours | 246 intervals, 2,908 probe-hours (Oct 63 / Nov 90 / Dec 93 intervals); 2,684 angle-resolved hours |
| probes | THEMIS a, d, e (b/c are ARTEMIS by 2015) |
| data | 3.2 s reduced-mode (peir), 50 angle × 24 energy bins |

## Detection funnel

868k candidate snapshots → 40k beam-flagged → 11,871 grouped episodes → hard gates (tailward, E_b ≤ 600 eV, beta < 0.03, Te < 500 eV) + model ≥ 0.82 + prominence ≥ 0.55 → **tier 1 = 86 episodes**, tier 2 (review) = 145.

## Tier 1 — 86 episodes

Label state: 66 positive, 0 negative, 20 unlabeled. The known event-boundary false positive (e_2015-12-31_083709, p = 0.837 under the 2026-08-11 model) no longer clears the 0.82 cut and sits in tier 2. Relative to the 2026-08-11 catalog: 81 kept, 8 dropped (4 labeled positive, 1 negative, 3 unlabeled), 5 added (4 labeled positive, 1 unlabeled); tiers 1 and 2 together are the same 231 episodes.

### By day

| day | episodes | probes |
|---|---|---|
| 2015-11-13 | 11 | d, e |
| 2015-12-06 | 1 | d |
| 2015-12-20 | 40 | a, d, e |
| 2015-12-21 | 7 | d, e |
| 2015-12-31 | 27 | a, d, e |

5 days total. By probe: d 40 / e 30 / a 16.

### Episode properties

| quantity | min | p25 | median | mean | p75 | max |
|---|---|---|---|---|---|---|
| duration (s) | 0* | 14.8 | 47.4 | 99.3 | 128 | 659.4 |
| n_samples | 1 | 4 | 8 | 15.9 | 21 | 131 |
| E_b median (eV) | 33 | 57 | 99 | 133 | 171 | 442 |
| model prob | 0.834 | 0.958 | 0.996 | 0.972 | 1.000 | 1.000 |

*zero duration = single-sample episode (one 3.2 s snapshot).

All energies in this file are measured (instrument frame) unless marked otherwise. True energies are higher by qΦ; see "Beam energies in the true frame" below.

Beam energy bands: <100 eV: 52 episodes, 100–200 eV: 20, 200–300 eV: 10, 300–500 eV: 4, >500 eV: 0. The population is dominantly sub-200 eV — squarely ionospheric.

### Plasma context (episode medians)

| quantity | min | median | max |
|---|---|---|---|
| Te (eV) | 18 | 54 | 391 |
| beta | 0.0003 | 0.0009 | 0.0272 |
| ne/ni | 1.06 | 2.22 | 5.17 |

Deep-lobe context throughout: beta median 9e-4, cold electrons, heavy-ion composition signature (ne/ni median 2.2; pure proton plasma reads ~1).

## Occurrence

Total tier-1 beam time: **2.37 h of 2,908 scanned probe-hours ≈ 0.08%** of tail-region residence (0.09% of the 2,684 angle-resolved hours). Longest single episode ~11 min.

## Physical events (cross-probe time merge)

Merging tier-1 episodes across probes by time proximity gives **8 distinct physical structures** at 30–60 min gap tolerance (10 at 10 min, where the two hour-scale structures each split once):

| event | span | episodes | probes | E_b range (eV, measured) |
|---|---|---|---|---|
| 11-13 14:08 | 21 min | 11 | d, e | 33–99 |
| 12-06 11:03 | single | 1 | d | 442 |
| 12-20 04:57 | single | 1 | d | 99 |
| 12-20 11:21 | 62 min | 39 | a, d, e | 34–304 |
| 12-21 10:53 | 9 min | 3 | d, e | 234–296 |
| 12-21 12:57 | 15 min | 4 | d, e | 296–441 |
| 12-31 08:32 | <1 min | 2 | d, e | 171–295 |
| 12-31 11:54 | 75 min | 25 | a, d, e | 57–441 |

Two structures are simultaneous three-probe detections; five of eight are multi-probe. The two longest (12-20 11:21, 12-31 11:54) are hour-scale persistent structures containing most of the catalog's episodes — episode count and physical event count are very different measures of occurrence. Against the 2026-08-11 catalog, the 12-21 16:49 single-episode structure and probe a's member of the 12-31 08:32 conjunction fell to tier 2.

## Beam energies in the true frame (added 2026-10-17)

Ions lose qΦ climbing the positive spacecraft potential, so the true energy is E_b + qΦ (singly charged; method and caveats in the README section "Spacecraft potential correction"). The numbers below use the current pooled catalog (`candidates/all/beam_catalog_all.csv`, tail2015). It has **89 tier-1 episodes** under the 2026-09-24 model and tier-1 rules, not the 86 of the 2026-09-19 tables above.

| quantity | min | p25 | median | mean | p75 | max |
|---|---|---|---|---|---|---|
| Φ (V) | 14.1 | 22.9 | 26.4 | 27.4 | 30.1 | 45.8 |
| E_b measured (eV) | 46 | 58 | 99 | 142 | 171 | 583 |
| E_b + qΦ (eV) | 71 | 86 | 126 | 169 | 204 | 624 |
| per-beam 1σ error (eV) | 9.3 | 11.1 | 17.5 | 22.8 | 29.5 | 46.9 |

| band | <100 eV | 100–200 | 200–300 | 300–500 | >500 |
|---|---|---|---|---|---|
| measured | 48 | 25 | 12 | 3 | 1 |
| true frame | 25 | 37 | 12 | 14 | 1 |

Median E_b/Te goes from 2.87 measured to 3.25 corrected (Te is already potential-corrected). In the true frame 62 of 89 beams (70%) are still below 200 eV, so the population stays dominantly sub-200 eV. The 40 eV tier-1 floor stays on measured energy and corresponds to about 66 eV true at the 2015 median potential. 86 of 89 beams have Φ ≤ 0.5 E_b (`fit_ok`), none has Φ > 0.8 E_b, and none is in Earth's shadow.

## Tier 2 — 145 episodes (review tier)

22 positive / 57 negative / 66 unlabeled, spread over 24 days, probes d 73 / e 43 / a 29. Energies broader and hotter (E_b median 134 eV, max 598) — the tier-2 mix is where the excluded populations live.
