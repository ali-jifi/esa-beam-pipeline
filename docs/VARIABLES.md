# Variable reference

Every variable in the pipeline, candidate records, models, episode layer, and output catalogs. Flux units are always energy flux, eV/(cm² s sr eV). "dex" = log10 units (0.3 dex ≈ factor 2).

## 1. Pitch-angle spectra (`PitchAngleSpectra`, beam_pipeline.py)

Built by `compute_pa_spectra()` from the 3D distribution + FGM B in DSL. All flux arrays are (ntime, nenergy).

| variable | meaning |
|---|---|
| times | timestep centers, unix seconds |
| energy | energy bin centers, eV |
| omni | solid-angle-weighted mean flux over all angles |
| para / anti / perp | cone-averaged flux: 0-30° / 150-180° / 75-105° pitch angle |
| para_sig / anti_sig / perp_sig | Poisson sigma of each cone average, from var(flux) = onecount × flux per bin |
| perp_floor | one-count noise floor of the perp cone average; R's denominator is clamped here |
| para_floor / anti_floor | one-count floors for the field-aligned cones; a sampled-but-empty cone reads as this upper limit, not NaN |
| pa_coverage_para / _anti / _perp | fraction of solid angle each cone actually sampled that spin |

## 2. Per-timestep detector features (`FeatureTable`, beam_pipeline.py)

Computed by `extract_features()`, one value per timestep.

| variable | unit | meaning |
|---|---|---|
| e_peak | eV | energy of max omni flux |
| e_beam | eV | energy of max \|asymmetry\| |
| width | – | normalized spectral width of the omni spectrum (narrow = beam-like) |
| asymmetry | – | (para − anti) / (para + anti) at e_beam, signed; + = parallel beam |
| para_to_omni | – | para / omni flux ratio at e_beam |
| energy_ratio | – | bulk flow energy / thermal energy from moments (dead in scoring, w = 0) |
| peak_prom | dex | prominence of the narrow spectral line inside the coherent run (0.3 = 2× above local baseline) |
| peak_width | bins | FWHM of that line |
| e_line | eV | energy of the line = **E_b**, the beam energy |
| de_line | eV | line FWHM = **ΔE** |
| eb_over_de | – | E_b / ΔE, monochromaticity; logged only, never gated |
| r_beam | – | **R**: flux-weighted mean of dominant cone / perp cone over the run bins |
| pa_max_ratio | – | max(para) / max(anti) over run bins, reported dominant/sub (≥ 1) |
| sig_margin | σ | weakest \|para − anti\| / sigma over the run bins |
| asym_baseline | – | rolling-median asymmetry over the run band, 3 h window, lobe-masked |
| asym_dev | – | asymmetry minus asym_baseline |
| flux_z | MAD-z | dominant-cone flux vs rolling median/MAD baseline, 6 h window |
| flux_z_perp | MAD-z | same z for the perp cone; separates beam from bulk compression |
| duration | steps | chain length of adjacent same-band candidates (1 step ≈ 100 s peif / 3.2 s peir) |
| chain_e_slope | dex/step | d log10(E) along the chain, TOF dispersion |
| chain_e_scatter | dex | rms log10(E) residual of the chain fit, drift coherence |
| coherent_ok | bool | a coherent directional run was found (asym + R + sigma gates over ≥ min_bins adjacent bins) |
| hyst_promoted | bool | run accepted at the lo Poisson bar (1.5σ) via a same-band significant neighbor |
| perp_depleted | bool | R denominator was clamped at the perp one-count floor |
| cone_floored | bool | run used an empty para/anti cone read at its one-count floor |
| pa_ok_both / pa_ok_para | bool | cone coverage sufficient for asym / for p2o |

## 3. Candidate JSONL record (`generate_candidates.py`)

One record per coherent-run candidate. Files: `candidates/*.jsonl` (literature batch, peif) and `candidates/tail2015/*.jsonl` (2015 batch, peir).

Top level:

| field | meaning |
|---|---|
| candidate_id | `{probe}_{trange-tag}_{HHMMSS}_{energy-bin-index}` |
| probe | THEMIS probe a-e |
| t_center / t_center_ut | timestep center, unix s / ISO UT |
| trange | interval the record came from |
| interval_tag | batch tag: beam_period / plasma_sheet / unknown |
| profile | strict (production params) or survey (0.5× sigma and R bars) |
| t_idx | timestep index into the interval's spectra |
| is_beam | heuristic classifier verdict at the active profile |
| direction | +1 parallel, −1 anti-parallel, 0 unknown |
| label / label_source / label_confidence / label_note | hand label and provenance; null until labeled |

`features{}` — snapshot of the FeatureTable at t (names as in section 2, except): `R` = r_beam, `E_b` = e_line, `dE` = de_line, `Eb_over_dE` = eb_over_de, `duration_steps` = duration, plus `beam_score` (the heuristic weighted score).

`gates{}` — strict-threshold pass/fail over the candidate band ±1 bin, all at production values regardless of profile:

| gate | passes if |
|---|---|
| sig | ≥ min_bins adjacent bins clear the lo Poisson bar on both \|para−anti\| and dominant−perp |
| asym | ≥ min_bins adjacent bins with \|asym\| ≥ 0.2 |
| dir | ≥ min_bins adjacent bins with R ≥ 1.2 |
| min_bins | ≥ min_bins adjacent bins clearing all three at once |

`context{}` — plasma environment at t:

| field | meaning |
|---|---|
| Bx_sign | sign of GSM Bx (hemisphere/lobe sense) |
| ne_ni_ratio | electron/ion density ratio; ≥ 2 flags O+ (composition proxy, ESA sees only protons cleanly) |
| dir_x_bx | direction × sign(Bx): **−1 = tailward = outflow-consistent**, +1 = earthward |
| beta | plasma beta = 0.4027 n[cm⁻³] T[eV] / B[nT]²; lobe ≲ 0.03 |
| in_lobe | beta < 0.1 |
| te_ev | electron temperature, eV; lobe indicator is cold (< 500 eV) |
| perp_depleted / cone_floored | copied bool flags from the FeatureTable |

## 4. Snapshot model features (`train_model.py`)

Transforms applied to record fields before the logistic fit. NaN → column median + a companion `{name}_missing` 0/1 indicator column.

| model feature | built from | transform |
|---|---|---|
| R, sig_margin, peak_prom, width, para_to_omni, asym_dev, flux_z_perp | features | as-is |
| abs_asym | asymmetry | absolute value |
| eb_over_de | Eb_over_dE | as-is |
| log_E_b / log_dE | E_b / dE | log10 |
| flux_z | flux_z | hinged: min(z, 5) |
| flux_z_excess | flux_z | hinged: max(z − 5, 0); moderate excursion reads beam, extreme reads population change |
| duration | duration_steps | × cadence → seconds, cadence-invariant across peif/peir |
| chain_e_slope | chain_e_slope | ÷ cadence → dex/second |
| chain_e_scatter | chain_e_scatter | as-is |
| hyst_promoted | hyst_promoted | 0/1 |
| anchor_fwd / anchor_bwd | candidate pool | log10(1 + seconds) to nearest later / earlier is_beam neighbor within 0.3 dex in energy, same probe+day |
| anchor_score | candidate pool | best heuristic beam_score among energy-matched neighbors within ±600 s |
| n_near | candidate pool | count of energy-matched neighbors within ±600 s, self excluded |
| log_te | context te_ev | log10; **model A only** (promoted 2026-07-31, provenance caveat: Te is a labeling feature) |
| ne_ni_ratio, log_beta, dir_x_bx | context | **model B only** — labels leaned on these, B is a ceiling not a result |

Feature sets: **A0** = morphology audit (no Te), **A** = A0 + log_te (deployed, used by `score_candidates.py`), **B** = A + context. Fitted coefficients: `candidates/analysis/model_report.md`.

## 5. Episode record (`group_episodes.py` → `episodes.jsonl`)

Beam-flagged snapshots grouped by: gap ≤ 60 s, same direction, adjacent members within 0.3 dex in E_b.

| field | meaning |
|---|---|
| episode_id | `{probe}_{date}_{HHMMSS}` of the first member |
| probe / file | probe and source jsonl |
| t_start_ut / t_end_ut / duration_s | first/last member time and span |
| n_samples | member count |
| direction / dir_x_bx | from the first member |
| e_b_med / e_b_min / e_b_max | member E_b stats, eV |
| te_med / ne_ni_med / beta_med | member context medians |
| flux_z_max / prom_max | member maxima |
| rep_id | member with the highest beam_score (used for cutouts) |
| member_ids | all member candidate_ids |

## 6. Episode model and catalog (`build_episode_catalog.py`)

Hard gates, applied upstream of the model, never features: dir_x_bx = −1 (tailward), e_b_med < 600 eV, beta_med < 0.03, te_med < 500 eV.

Episode model features (aggregated over members; NaN → median + missing indicator):

| feature | meaning |
|---|---|
| log_n | log10 n_samples |
| prom_max / flux_z_max / r_max / sig_max | member maxima of peak_prom, flux_z, R, sig_margin |
| p2o_med / fzp_med / asym_dev_med | member medians of para_to_omni, flux_z_perp, asym_dev |
| log_eb_med | log10 e_b_med |
| eb_span_dex | log10(e_b_max / e_b_min), energy wander of the episode |

Catalog cut (frozen 2026-08-11): model fires at prob ≥ 0.82, prominence fires at prom_max ≥ 0.55. **Tier 1 = both fire, tier 2 = exactly one.**

## 7. Output CSV columns

`candidates/analysis/model_scores.csv`: candidate_id, prob (model A probability), is_beam (heuristic), label.

`candidates/analysis/beam_catalog.csv` (snapshot-level): tier, class, candidate_id, t_center_ut, prob, is_beam, direction, E_b, R, peak_prom, flux_z, duration_steps, label, label_confidence.

`candidates/tail2015/analysis/episode_features.csv`: episode_id, day, label, n_samples, log_duration, then the episode model features (section 6) plus log_te, ne_ni, log_beta context columns.

`candidates/tail2015/analysis/beam_catalog_2015.csv` (the 2015 catalog): tier, episode_id, t_start_ut, duration_s, n_samples, prob, prom_max, e_b_med, te_med, ne_ni_med, beta_med, label. prob is the lobe-scoped episode model probability; the blind-validated operating point is 0.82.
