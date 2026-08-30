# Surface-layer results — 2015 catalog statistics

Descriptive statistics of the 2015 beam catalog (candidates/tail2015/analysis/beam_catalog_2015.csv). All numbers computed directly from the catalog files 2026-08-17. Interpretation, validation, and physics-level results are in RESULTS_HIGH.md.

## Input coverage

| quantity | value |
|---|---|
| survey scope | all of 2015; tail region x_GSM < −8 Re, \|y\| < \|x\| |
| qualifying coverage | mid-Oct through Dec only (365-day scan verified; Jan–Sep has zero qualifying intervals) |
| intervals / probe-hours | 246 intervals, 2,908 probe-hours (Oct 63 / Nov 90 / Dec 93 intervals) |
| probes | THEMIS a, d, e (b/c are ARTEMIS by 2015) |
| data | 3.2 s reduced-mode (peir), 50 angle × 24 energy bins |

## Detection funnel

868k candidate snapshots → 40k beam-flagged → 11,871 grouped episodes → hard gates (tailward, E_b ≤ 600 eV, beta < 0.03, Te < 500 eV) + model ≥ 0.82 + prominence ≥ 0.55 → **tier 1 = 89 episodes**, tier 2 (review) = 142.

## Tier 1 — 89 episodes

Label state: 66 positive, 1 negative, 22 unlabeled. The single negative is the known event-boundary blind false positive (e_2015-12-31_083709, p = 0.837 — also the tier-1 probability minimum).

### By day

| day | episodes | probes |
|---|---|---|
| 2015-11-13 | 12 | d, e |
| 2015-12-06 | 1 | d |
| 2015-12-20 | 40 | a, d, e |
| 2015-12-21 | 7 | d, e |
| 2015-12-31 | 29 | a, d, e |

5 days total. By probe: d 38 / e 33 / a 18.

### Episode properties

| quantity | min | p25 | median | mean | p75 | max |
|---|---|---|---|---|---|---|
| duration (s) | 0* | 8.7 | 34.7 | 87.8 | 115.1 | 659.4 |
| n_samples | 1 | 3 | 6 | 14.2 | 20 | 131 |
| E_b median (eV) | 33 | 57 | 99 | 135 | 171 | 454 |
| model prob | 0.837 | 0.960 | 0.994 | 0.968 | 0.997 | 1.000 |

*zero duration = single-sample episode (one 3.2 s snapshot).

Beam energy bands: <100 eV: 53 episodes, 100–200 eV: 22, 200–300 eV: 7, 300–500 eV: 7, >500 eV: 0. The population is dominantly sub-200 eV — squarely ionospheric.

### Plasma context (episode medians)

| quantity | min | median | max |
|---|---|---|---|
| Te (eV) | 18 | 54 | 391 |
| beta | 0.0003 | 0.0009 | 0.0272 |
| ne/ni | 1.06 | 2.14 | 5.17 |

Deep-lobe context throughout: beta median 9e-4, cold electrons, heavy-ion composition signature (ne/ni median 2.1; pure proton plasma reads ~1).

## Occurrence

Total tier-1 beam time: **2.17 h of 2,908 probe-hours ≈ 0.075%** of tail-region residence. Longest single episode ~11 min.

## Physical events (cross-probe time merge)

Merging tier-1 episodes across probes by time proximity gives **9 distinct physical structures** — stable at every gap tolerance from 10 to 60 min:

| event | span | episodes | probes | E_b range (eV) |
|---|---|---|---|---|
| 11-13 14:08 | 24 min | 12 | d, e | 33–99 |
| 12-06 11:03 | 2 min | 1 | d | 442 |
| 12-20 04:57 | 2 min | 1 | d | 99 |
| 12-20 11:21 | 62 min | 39 | a, d, e | 34–304 |
| 12-21 10:54 | 9 min | 2 | d, e | 234–295 |
| 12-21 12:58 | 15 min | 4 | d, e | 296–441 |
| 12-21 16:49 | 2 min | 1 | e | 170 |
| 12-31 08:32 | 6 min | 4 | a, d, e | 171–454 |
| 12-31 11:59 | 71 min | 25 | a, d, e | 57–175 |

Three structures are simultaneous three-probe detections; six of nine are multi-probe. The two longest (12-20 11:21, 12-31 11:59) are hour-scale persistent structures containing most of the catalog's episodes — episode count and physical event count are very different measures of occurrence.

## Tier 2 — 142 episodes (review tier)

22 positive / 56 negative / 64 unlabeled, spread over 24 days, probes d 75 / e 40 / a 27. Energies broader and hotter (E_b median 134 eV, max 598) — the tier-2 mix is where the excluded populations live.
