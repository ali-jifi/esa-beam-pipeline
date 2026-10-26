# Data Analysis Project #1

Detect and characterize field-aligned ion beams in THEMIS ESA particle data by reducing 3D distributions into pitch-angle-resolved energy spectra and classifying directional, narrow-line beam signatures against isotropic plasma-sheet background. On top of the per-timestep detector sits a candidate/labeling layer, a trained logistic model, and an episode-level catalog builder used to produce a validated 2015 tail-season beam catalog.

Built with care, by Ali Jifi-Bahlool

# esa_plotting

Python package for plotting and analyzing THEMIS satellite Electrostatic Analyzer (ESA) data, with an ion beam detection pipeline built on top of pyspedas.

## Core package - `src/esa_plotting/`

- **config.py** - sets up where satellite data lives on disk. `set_data_dir()` checks for a path passed in, then falls back to the `THM_DATA_DIR` env var, then defaults to `./data`. Also defines plot defaults (color range, energy range, color map).
- **probes.py** - THEMIS has 5 probes labeled a-e. Defines that tuple and a helper `eflux_var()` that builds the tplot variable name for a given probe/species combo, such as `tha_peif_en_eflux`.
- **loaders.py** - thin wrappers around pyspedas to load three types of THEMIS data:
  - `load_esa()` - energy flux spectra (the main data product)
  - `load_esd()` - 3D particle distributions
  - `load_fgm()` - magnetic field data

  All return tplot variable names that pyspedas stores internally.
- **plotting.py** - two plotting helpers:
  - `configure_eflux_panel()` - sets up a single energy flux spectrogram panel (log scales, color range, labels)
  - `stack_plot()` - takes multiple tplot variables, filters out ones that didn't load, and renders them stacked into a PNG
- **beam_pipeline.py** - the beam detection pipeline (below).
- **__init__.py** - re-exports the public api from all of the above.

## Beam detection pipeline - `beam_pipeline.py`

Detects ion beams in THEMIS data through 6 phases.

### Phase 0 - Data acquisition

Loads the data products for a given probe/time range:

- `load_esd_distribution()` - reads the CDF directly (not through pyspedas) to get the full 3D ion distribution including the angle/energy lookup tables that pyspedas doesn't expose. Supports both `peif` (full-distribution survey, ~97-122 s cadence = 32 spins, 32 energies x 88 angles) and `peir` (reduced mode, 3.2 s cadence, 24 energies x 50 angles) via `datatype`. Those are the populated modes; the CDF arrays themselves are 32 x 176 for `peif` and 32 x 88 for `peir`, with the unused slots masked off. Parses time ranges, clips to the requested interval, and extracts the dominant energy/angle mode. Also pulls the calibration arrays (`geom_factor`, `gf`, `eff`, `integ_t`) and builds the per-bin one-count eflux level `1/(integ_t * geom_factor * gf * eff)`, the inverse of `thm_convert_esa_units`, used downstream as the cone noise floors.
- `load_bfield_dsl()` - loads magnetic field in DSL (Despun Spacecraft L) coordinates, trying survey/low/high cadence in order.
- `load_bfield_gsm()` / `load_state_gsm()` - GSM field and spacecraft position, used for plasma beta and tail-region context.
- `load_moments()` - loads bulk plasma parameters (density, velocity, temperature).

### Phase 1 - Spectral reduction

Takes the full 3D ion distribution and collapses it into 1D energy spectra based on pitch angle relative to the magnetic field.

Pitch angle is the angle between a particle's velocity and the local B-field direction. 0° = moving along B (field-aligned/parallel), 180° = moving opposite to B (anti-parallel), 90° = perpendicular.

The four spectra:

- **Omnidirectional** - average flux at each energy across all angles, weighted by solid angle (domega). This is what a standard spectrogram shows.
- **Parallel (0-30°)** - only bins where the particle velocity is roughly field-aligned. A beam streaming along B shows up here.
- **Anti-parallel (150-180°)** - only bins moving opposite to B. A beam coming from the other direction shows up here.
- **Perpendicular (75-105°)** - bins moving roughly across B. A field-aligned beam depletes this cone, so it serves as a clean directional background, used as the denominator of R below.

A beam is a narrow, directional population, it lights up in one PA gate but not the other. Plasma sheet ions are roughly isotropic, so all curves overlap. Comparing parallel and anti-parallel spectra surfaces the directional asymmetry that defines a beam.

`compute_pa_spectra()`: for each timestep, interpolates the magnetic field to the distribution time, computes pitch angles for every angle bin, then sorts flux into omni / para / anti / perp energy spectra. It also propagates the per-bin one-count level through the same solid-angle weighting to get per-cone noise floors, and propagates Poisson counting uncertainty (`var(flux) = onecount * flux` per bin) into per-cone sigmas (`para_sig` / `anti_sig` / `perp_sig`) used by the significance gates below. **Cone one-count floors**: a sampled-but-empty para/anti cone reads as its one-count upper limit instead of NaN - without this the pipeline is blind to the strongest lobe beams, since a pure anti-parallel beam *should* have an empty parallel cone. For the geometry, the CDF look-direction angles are the particle travel direction (verified against L2 velocity moments), so they dot directly with the B-field unit vector. The work is vectorized over chunks of timesteps, with results bit-identical to a per-timestep loop (see Optimization).

### Phase 2 - Feature extraction (`extract_features()`)

For each timestep, computes the discriminating features.

Spectral + moment features:

- `e_peak` - energy of peak flux
- `width` - how broad the spectrum is (narrow = beam-like)
- `asymmetry` - `(para - anti) / (para + anti)` near the peak, picked from the bin with max `|asym|` across a window. Positive = field-aligned beam.
- `para_to_omni` - how much stronger the parallel flux is vs omnidirectional
- `energy_ratio` - bulk flow energy / thermal energy from moments (high = directed flow)

Coherent-run + spectral-line features:

- `coherent_ok` - bool, a real coherent directional run was found (both cones sampled, enough adjacent bins where `|asym|` and **R** clear their per-bin thresholds AND both are statistically significant). Significance uses a **hysteresis pair of Poisson bars**: bins clear at `n_sigma_lo=1.5` if they sit beside a same-energy significant neighbor, isolated bins must clear `n_sigma_hi=2.5` alone (benchmarked +9.5 recall / -0.3 precision vs a single bar; promoted bins tracked via `hyst_promoted`). **R = dominant cone / perpendicular cone**, not dominant/omni: omni includes the beam cone so it dilutes the enhancement, while the perp cone is a clean background a field-aligned beam depletes. The perp denominator is clamped at its one-count noise floor so a depleted perp can't blow R up. Runs above `beam_e_max=7000` eV are rejected.
- `perp_depleted` / `cone_floored` - bools, the run leaned on a perp clamped at its floor, or on an empty para/anti cone read at its one-count floor. Tracked for calibration, not gated.
- `peak_prom` - log10 prominence of the narrow spectral line found inside the coherent run (0.3 = 2x above local baseline)
- `peak_width` - width of that line in bins at half its log10 prominence (scipy `find_peaks`, so not a linear-flux FWHM)
- `e_line` - energy of the line in eV (this is **E_beam**)
- `de_line` - **ΔE**, the line width in eV at half its log10 prominence, interpolated between bins. It is not a true linear-flux FWHM and sits near the bin width below about 100 eV
- `eb_over_de` - `e_line / ΔE` = beam monochromaticity; logged only, never gated (its median-vs-R trend was shown to be pure selection)
- `r_beam` - flux-weighted mean R over the coherent-run bins
- `pa_max_ratio` - `max(para) / max(anti)` over the run bins, reported dominant/sub (>= 1)
- `sig_margin` - weakest `|para - anti| / sigma` over the run bins

Temporal-context features (logged only, consumed by the model layer, never gated in the heuristic):

- `asym_baseline` / `asym_dev` - rolling-median asymmetry over the run band (3 h window) and the deviation from it
- `flux_z` / `flux_z_perp` - dominant-cone and perp-cone flux vs rolling median/MAD baseline (6 h window); separates beam enhancement from bulk compression
- `duration` - chain length of adjacent same-band candidates in steps
- `chain_e_slope` / `chain_e_scatter` - d log10(E) per step along the chain and the rms residual of the fit (time-of-flight dispersion / drift coherence)

`lobe_baseline_mask()` restricts baseline statistics to lobe-like intervals (beta cut) so the rolling baselines aren't poisoned by plasma-sheet crossings.

Spectral-line detection is local to the coherent run, not global. It scans only the dominant cone, compresses to finite/positive bins, takes `log10(flux)`, and runs `scipy.signal.find_peaks` with a bounded prominence window (`peak_wlen=5`) and a width cap (`peak_width_max=4.0` bins). A peak only counts if it sits inside the directional run band (±1 bin slop). The idea: a beam = the directional region is also a narrow line; a prominent line elsewhere (e.g. the anti-parallel plasma-sheet peak) is rejected.

### Phase 3 - Classification (`classify_beams()`)

Single-path score classifier, wrapped in an AND-gate.

**Score** - weighted sum of normalized feature scores. Weights: `w_r_beam=0.50`, `w_peak_prom=0.30`, `w_asymmetry=0.10`, `w_width=0.05`, `w_para_to_omni=0.05` (`w_energy_ratio=0.0`, dead code). Each component ramps to 0.5 at its own threshold and caps at 1; the R ramp anchors on `r_beam_min=8.0`, which splits the labeled classes. Beam candidate if score clears `score_threshold=0.6`.

**AND-gate:** final `is_beam = score_ok and gate`, where `gate = coherent_ok and peak_prom >= peak_prom_min`. Both detectors, the directional coherent run AND a narrow spectral line at the same energy, must agree, which kills noise that fires only one signal alone.

The **hard-rule fallback is gone** (removed after label calibration): the bypass predecided 39/74 labeled cases and froze the score. Asymmetry still guards upstream through the per-bin gate and the hysteresis pass, so `asymmetry_min` / `width_max` / `para_to_omni_min` now act only as score-ramp anchors.

Line-detection params on `ClassifierParams`: `peak_prom_min=0.3` (log10, so 0.3 = 2x above local baseline), `peak_width_max=4.0` (FWHM cap in bins), `peak_wlen=5` (local prominence window).

Beam direction is tagged (+1 parallel, -1 anti-parallel) from the asymmetry sign.

### Phase 4 - Temporal smoothing

- `smooth_labels()` - requires N consecutive beam-flagged timesteps to keep a beam interval. Default `min_consecutive` is 1 (keep isolated beams) because the AND-gate already enforces precision; raise it to suppress more.
- `threshold_sensitivity()` - sweeps each threshold parameter and reports how beam count changes, so you can see how stable the classification is.

### Phase 5 - Plotting

- `plot_feature_timeseries()` - 8-panel overview: omni spectrogram, E_peak, width, asymmetry, peak prominence, beam score, a classification color bar (red = parallel beam, blue = anti-parallel, orange = unknown direction, gray = no beam), and a bottom spectrogram with beam detections overlaid at `e_peak`. Takes `ClassifierParams` so every threshold guide-line is driven by the actual params.
- `plot_curated_snapshots()` - picks representative timesteps (confirmed beams, plasma sheet, borderline cases) and plots the three-curve energy spectra at each.
- `diagnose_window()` - dumps per-timestep spectra and features, including the per-bin omni/para/anti/perp flux and R, plus the line, R, and context features above.

Per-beam outputs (one record per flagged timestep, written when plotting is on):

- `write_beam_table()` - `<prefix>_beams.csv`: UT, direction, `e_beam`, `delta_e`, `eb_over_de`, `r_beam`, `pa_max_ratio`, `sig_margin`, asymmetry, `e_peak`, beam score.
- `plot_beam_histograms()` - `<prefix>_histograms.png`, 4 panels: R, E_beam, ΔE, and E_beam/ΔE over the flagged timesteps.
- `plot_threshold_comparison()` - `<prefix>_threshold_compare.png`, small multiples: one omni spectrogram per `coherent_dir_min` (R) value with beam dots overlaid. Re-runs `extract_features` per value since the R gate lives there.

`run_pipeline()` ties it all together, returns everything in a `PipelineResult` dataclass.

## Candidate, model, and catalog layers - `scripts/`

The per-timestep pipeline is deliberately permissive: it *generates candidates*, it doesn't judge them. The layers above it turn candidates into a labeled dataset, a trained model, and an episode-level catalog.

### Candidate generation and labeling

- `generate_candidates.py` - wraps the pipeline without changing it; emits one JSONL record per coherent-run candidate plus a PNG spectrum cutout each. Two profiles: **strict** (production `ClassifierParams`, single source of truth) and **survey** (sigma bars and R gate at 0.5x, cast wide and let the strict gate bitmask sort it). Each record carries features, per-gate pass/fail at strict thresholds, and plasma context (beta from moments + GSM B, ne/ni ratio, Te, `dir_x_bx` = beam direction x sign(Bx), where -1 = tailward = outflow-consistent). Batch mode via `--events` CSV; `--datatype peir` for reduced-mode 3.2 s data.
- `contact_sheet.py` - tiles candidate cutouts into contact sheets for fast eyeball labeling.
- `rerender_cutouts.py` - redraws cutout PNGs from existing JSONLs, records never touched. `--existing` limits to cutouts already on disk, `--catalog` puts the model prob in the cutout title (episode catalogs map through the `episodes.jsonl` beside them), `--cut-dir` / `--datatype peir` for the 2015 layout.
- `migrate_labels.py` - carries labels from a backup of candidate JSONLs into regenerated ones (join on candidate_id, fallback probe + t_center).
- `analyze_candidates.py` - pools candidate JSONLs (deduped on physical timestep) and plots feature distributions for auto-label calibration.
- `candidates/LABELING.md` - the labeling standard: kill rules, derived constants, case law, frozen-labels discipline. Labels live inline in the JSONLs; `labels_frozen_*.csv` snapshots them.

### Snapshot model

- `train_model.py` - logistic regression on the labeled candidate dataset, event-level `GroupKFold` splits, benchmarked against the heuristic classifier. Feature set `FEATS_A` is detection + context morphology (R, sig_margin, peak_prom, flux_z family, duration, chain features, asym_dev, anchor/neighbor episode context, log_te); `FEATS_B` adds the labeling-adjacent context features (ne/ni, beta, direction) and is kept separate as an audit trail.
- `score_candidates.py` - two-tier snapshot catalog: heuristic generates, model judges. Tier 1 = both fire (high confidence), tier 2 = exactly one fires (review queue). Threshold self-calibrates off the out-of-fold PR curve unless given.
- `analyze_r_monochromaticity.py` - diagnostic for why median E_beam/ΔE falls as the R gate tightens (answer: selection, which is why `eb_over_de` is never gated).

### Episode layer and the 2015 catalog

- `group_episodes.py` - groups beam-flagged snapshots into episodes: gap <= 60 s, same direction, adjacent members within 0.3 dex in energy. Writes `episodes.jsonl` with per-episode stats (duration, E_b range, context medians, representative member).
  - **Why 0.3 dex (a factor of 2).** Adjacent reduced-mode energy bins near 60 to 100 eV differ by a factor of about 1.7 (0.24 dex), and the beam peak wobbles by one bin between spins from counting noise. 0.3 dex is the smallest tolerance that absorbs that wobble. A sensitivity test on tail2015 (2026-09-27, frozen model, tier-1 rules applied) shows factor 2 sits on a plateau: loosening to 2.5 (0.4 dex) leaves tier 1 unchanged, while tightening to 1.6 (0.2 dex), below one bin step, splits beams into fragments.

    | energy-link factor            |     1.6    | 2.0 (used) | 2.5  |
    |-------------------------------|------------|---------|---------|
    | episodes                      | 16,297     | 11,871  | 11,343  |
    | tier 1                        | 214        | 89      | 89      |
    | tier-1 beam time, h           | 1.70       | 2.48    | 2.48    |
    | tier-1 structures, same probe | 24         | 25      | 25      |
    | tier-1 events, all probes     | 9          | 8       | 8       |
    | days with tier 1              | 6          | 6       | 6       |
    | labeled tier-1 negatives      | 3 of 204   | 1 of 85 | 1 of 85 |

    Events, active days and contamination are insensitive to the choice, so they are the robust units for comparing seasons. Episode counts, and rates built on them, depend on it. The contamination row is in-sample, since most 2015 labels trained the model. Only 2015 was tested.
- `build_episode_catalog.py` - episode features, lobe-scoped logistic model, and the two-tier episode catalog. **Hard gates are upstream of the model, never features**: tailward only (`dir_x_bx = -1`), E_b < 600 eV (outflow energy window), beta < 0.03 (lobe scope), Te < 500 eV (closes the corridor leak, and a missing `te_med` is admitted rather than rejected). All read the episode median except direction, which comes from the first member. Catalog cut at model prob 0.82 (blind-validated: contamination 1/14 at this threshold) plus prominence >= 0.55. Model features are defined on the beam's dominant cone (`d2o_med`, `asym_dev_dir`), not the parallel cone: the 2015 lobe set is 269/282 parallel beams, and the original para-cone features flipped sign or collapsed for anti-parallel beams, which scored the northern-lobe seasons near zero. Refit 2026-09-19 on the 2015 labels (day-blocked AUC 0.938). The default run applies the frozen model in `candidates/tail2015/analysis/episode_model.joblib` to whichever `--season` is given; `--fit` refits on that season's labels and overwrites the model file.
- `render_episode_sheets.py` - renders episode representative cutouts and tiles contact sheets from a review-queue CSV.
- `render_beam_members.py` - per-member cutouts for every catalog beam, one folder per episode under `analysis/beams/tier{n}/<episode_id>/`, refreshes the rep cutouts in `analysis/cutouts/` on the way.

Outputs land in `candidates/tail2015/analysis/`: `beam_catalog_2015.csv` (tier 1 = 86 episodes / tier 2 = 145 under the 2026-09-19 refit; the 2026-08-11 catalog with 89 / 142 is kept beside it as `beam_catalog_2015_frozen20260811.csv`), `episodes.jsonl`, the blind-validation deck, review-round queues and sheets, and the coldline census (`coldline_final.csv`, 41 members under a six-constant frozen definition).

### Multi-season runs

- `scan_tail_intervals.py` - builds `events_tail<year>.csv` from L1 state: x_GSM < -8 Re and |y| < |x| on probes a/d/e, merged into hourly intervals per probe-day (reproduces `events_tail2015.csv` exactly). Pad the scan window a month either side of the expected season and let the geometry decide.
- `group_episodes.py`, `build_episode_catalog.py`, `add_episode_positions.py`, `add_episode_scpot.py`, `link_events.py`, `render_beam_members.py`, `render_episode_sheets.py` all take `--season tail<year>` (default `tail2015`) and read/write `candidates/tail<year>/`. Per season: `generate_candidates.py --events events_tail<year>.csv --datatype peir --prefetch`, then `generate_candidates.py --events events_tail<year>.csv --profile survey --datatype peir --no-cutouts --no-update --workers 12 --out candidates/tail<year>` (see Optimization), then the chain above in that order. Generation checkpoints per interval, so rerunning the same command resumes.
- `pool_seasons.py` - concatenates the per-season catalogs and episodes into `candidates/all/` and writes `coverage.csv`. It joins each season's potential sidecar (`analysis/episode_scpot.jsonl`) beside the measured energies. Episodes join on `episode_id` and `rep_id`. In the catalog, which has no `rep_id`, an `episode_id` shared by two episodes is left blank rather than joined to the wrong row. Coverage is counted in angle-resolved hours (reduced-mode samples with the 50-angle packing inside each processed interval, cached per season in `analysis/coverage_intervals.csv`): in the early seasons fast survey covers only about half of each tail interval, so scanned hours overstate exposure. Tier 1 per 100 angle-resolved hours is the cross-season number.
- Reduced-mode intervals outside the fast-survey schedule fail in the loader with `no angle-resolved peir samples`, are logged as `[fail]`, and count as unprocessed. Expect a block of them in the padding months.

### Current model and catalog status (2026-09-24)

- The frozen model has 16 inputs: the ten shape features above, log Te and log beta, and four features computed from the snapshot files (`log_iso`, beam-flagged snapshots within 600 s of the episode on the same interval; `duty`, members over snapshots in the episode span; `e_stab`, share of members within 0.12 dex of the median beam energy; `eb_drift`, beam-energy slope in dex per minute). It is fit on every user label from 2011 to 2016. Held-out-season AUC is 0.974 pooled and 0.905 for 2013, the weakest season. Earlier catalogs are kept beside each season's catalog as `_frozen20260919` (direction-invariant model), `_tebeta20260923` (plus Te and beta) and `_shape20260924` (plus the four shape features, before the tier-1 rules below).
- Tier 1 also requires at least 2 samples and a median beam energy above 40 eV. These are part of the catalog definition, not tuned cuts: a single spin, or a peak in the lowest energy bin (about 34 eV), cannot be a measured beam. Episodes that fail either move to tier 2. Across the 15 seasons this moved 232 of 1,409 tier-1 episodes, 157 of them single spins, leaving 1,177.
- Tier-1 contamination from the blind decks is roughly 5 to 13% in the quiet seasons, 2011 to 2019, and roughly 25 to 35% in the solar-maximum seasons, 2022 to 2024. These are working estimates, not yet rescored under the tier-1 rules above or under the pending density-ratio and 300 to 600 eV rulings. Every catalog row carries the model probability, so a stricter cut can be applied downstream.
- Reduced-mode detection has a flux threshold, because a 3-second spin collects too few counts for very faint cold beams. On published events this costs about 4% of in-scope beams, typically the faint edges of events that are otherwise detected.

### Spacecraft potential correction (2026-10-17)

A sunlit spacecraft in thin lobe plasma charges positive. Incoming ions climb the potential barrier and lose qΦ, so the measured beam energy sits below the true energy: E_true = E_meas + qΦ for singly charged ions (McFadden et al. 2008, first-results paper, Fig. 3 and p. 486). In the 2011 to 2016 tier-1 beams, Φ has a median of 30 V (16 to 47 V at the 5th to 95th percentile). That is a median of 0.3 of the measured E_b, and 0.5 in the 40 to 60 eV bin, so the correction matters most for the coldest beams.

- **Script.** `add_episode_scpot.py --season tail<year>` runs after `add_episode_positions.py`. It takes the EFI-derived `peir_sc_pot` from the ESA L2 file (falling back to `peif_sc_pot`) at each member spin and adds it to that spin's peak energy before taking episode statistics, because the potential moves by a median of 8.5 V within one beam. It writes `analysis/episode_scpot.jsonl`, keyed by `episode_id` and `rep_id` (two peaks in one spin can share an `episode_id`), and `pool_seasons.py` joins it into the pooled catalog and episodes. The raw `e_b_med` stays the measured value and the corrected values sit beside it. Columns are listed in `docs/VARIABLES.md`.
- **What stays in the measured frame.** All gates (the 40 eV tier-1 floor and the 600 eV window) and the model feature `log_eb_med` stay on measured energy, because detection happens in the instrument frame. Rewriting the floor in true energy would make it a potential-dependent cut and change catalog membership. At the median potential the 40 eV floor corresponds to about 70 eV true (65 to 75 eV interquartile). Each beam's own bound is in `eb_floor_true` (40 eV + Φ), the truncation bound for any fit. The correction does not recover beams below the barrier.
- **What uses the true frame.** Corrected E_b, E_b/Te and all comparisons across seasons or probes. Te is the L2 reduced-mode electron moment, which the ground processing corrects for the potential (L2 quality flag 1 means "no spacecraft potential available"), so E_b/Te is consistent only with corrected E_b. The tier-1 median E_b/Te goes from 1.95 measured to 2.8 corrected.
- **Per-beam error.** `e_b_corr_err` combines two terms in quadrature: σ_E² = (f·E_meas)² + σ_Φ².
  - The first term is rounding to the reduced-mode energy bins. A beam pinned to one bin has an error of w/√12, where w is the mean gap to the two neighbouring bins. The bins are 0.238 dex apart up to the ~300 eV bin, 0.174 dex from there to the ~450 eV bin and 0.12 dex above. So f is about 17% below 290 eV, 15% at the ~300 eV bin, 10% at the ~450 eV bin and 8% at the ~580 eV bin (band edges 290, 366 and 514 eV). 30 of the 303 tier-1 beams are above 290 eV. Taking the episode median doesn't shrink this term when the true energy stays in one bin.
  - The second term, σ_Φ = 5 V, is the EFI-to-plasma calibration in thin plasma. Eq. 1 of the McFadden instrument paper uses a scale factor of 1.15 and a 1.0 V offset fitted on one probe. Probes differ by about 5%, with a ratio that isn't constant, and the offset grows to about 2 V at low density, so a few volts at 30 V is plausible.
  - At 58 eV the rounding term is about 10 eV and dominates.
- **Flags.** These beams stay in the catalog but should be left out of quantitative fits, or dropped as a sensitivity test.
  - `sheath_flag`, Φ > 0.8 E_b (8 of 303 tier-1): the sheath bends trajectories as well as energies, so arrival direction, pitch angle and dir×Bx are unreliable.
  - `fit_ok`, Φ ≤ 0.5 E_b (241 of 303): the beams to use in quantitative fits.
  - `eclipse_flag`, episode midpoint in Earth's shadow from ephemeris: x < 0 and √(y²+z²) < 1.1 RE, which also covers penumbra (GSE and GSM share x, so the stored GSM position works). In shadow the EFI potential is invalid even when it reads positive: values there pile up at 1.0 V or near 98 V. The fly-wheel spin phase also drifts, which corrupts pitch angle and dir×Bx. Flag from geometry, never from the sign of Φ. No 2011 to 2016 tier-1 beam is in shadow. 20 tier-2 beams are (16 in 2013, 4 in 2014).
- **Velocity.** E_b and E_b + qΦ are energies per charge, and the ESA does not separate species. Any velocity is reported for both H+ and O+: v = √(2qE/m), which is 13.8·√E km/s for H+ and 3.47·√E km/s for O+, with E in eV. Species is not assigned from ne/ni alone, because cold ions hidden below the barrier produce the same density mismatch as O+ (McFadden et al. 2008, §2.2).

Checks and open caveats (2011 to 2016, 2026-10-15 to 10-17):

- **Potential against density.** In lobe episodes outside shadow, Φ ≈ 4.5 − 37.4·log₁₀(ne) on probe-day medians (R² 0.55).
  - The higher potential on probe A in the tier-1 sample is explained by density: A sits in thinner plasma, and at the same density it isn't higher than D or E.
  - Of the 11 V raw gap between 2011 to 2013 and 2014 to 2016, about 6 V remains at the same density. F10.7 is a small positive effect (0.030 ± 0.015 V per sfu, about 1.5 V across the range) and doesn't change the era term.
  - The leftover is not a clean step. 2013 D and E sit 8 and 13 V high for the whole season, independent of distance from the shadow, and 2016 sits 4 to 7 V low. It may be instrumental. The EFI operations history for those years hasn't been checked, and the documented bias changes cover 2007 to 2008 only. Treat it as a per-era systematic of about σ_Φ.
  - Lobe ne is not fully independent of Φ, because cold electrons and photoelectrons both arrive near eΦ.
- **Season comparison.** In the raw corrected medians the gap between 2011 to 2013 and 2014 to 2016 widens (log E_b difference −0.069 dex measured, −0.092 corrected). Most of that is truncation: 2011 to 2013 has a higher true-frame floor (median 74 eV against 67 eV). Above a common 90 eV true floor the difference is not significant (p = 0.37). The truncated regression on `eb_floor_true` that would settle it hasn't been built.
- **Line width.** Below about 100 eV measured, `de_line` is limited by the instrument and bin width. A model with an instrumental and a thermal term fits better than either alone, but the fitted kT_∥ ranges from 0.9 to 6 eV depending on model choices, because `de_line` is a width at half the log prominence, not a true linear FWHM. Only kT_∥ ≲ 6 eV is robust, well below the 10 to 20 eV oxygen moment temperature in Artemyev et al. 2020.
- **Te bias.** 28% of tier-1 beams have Te < Φ. Their electron population sits near the photoelectron cutoff, where wide bins and photoelectrons bias the moment, so lean on their E_b/Te least. A 10 V error in Φ at fixed density moves Te by about −7%.
- **L2 quality flags.** The meanings in the CDF TEXT attribute are out of date. Per themis.ssl.berkeley.edu/esa_flag.shtml they are: 1 no spacecraft potential, 2 counter saturation, 4 solar wind outside solar-wind mode, 8 slow-survey single angle bin, 16 Earth shadow, 32 lunar shadow, 64 maneuver, 128 ne and ni differ by more than 2×.
  - Flag 128 is set on 86% of tier-1 beams. Their median ne/ni is 2.1, against 1.2 for the rest.
  - Flag 2 is set on 62% of lobe episodes, together with high Φ and low density. One unverified reading is photoelectrons saturating the low-energy electron counters.
  - No tier-1 beam carries flag 16. It agrees only partly with the ephemeris flag.
- **Artemyev et al. 2020 frame.** The paper never says whether its ion energies are potential-corrected. Its event criteria (|Bx| > 0.9 B_lobe, tailward flow faster than 10 km/s, n_O > n_p, parallel flux above 5× average, all held for at least 5 min) involve no energy, so recall against Table S1 is frame-free. Its W_O is a bulk flow energy from oxygen moments, not a spectral peak, so it isn't comparable to E_b in either frame.
  - Anchor check on 2015-12-20 THEMIS-D: tier-1 corrected E_b has a median of about 109 eV at Φ ≈ 25 V. That is half the 210 eV implied by O+ at the 50 km/s of their Fig. 1h.
  - That is within the factor-of-2 tolerance a moment-to-peak comparison allows. The comparison is ambiguous either way, because Fig. 1h may be the proton-assumed moment.
- **Predicting E_b from location and activity (2026-10-19), null result.**
  - Data: tier-1 beams with `fit_ok` set and outside shadow, from 2011 to 2019 plus 2022 to 2024 at prob ≥ 0.90, 718 beams on 107 days.
  - Method: day medians, 1-min OMNI averaged over the 30 minutes before each beam, and season-out cross-validation. Peak energy is interpolated below bin resolution with a 3-point log parabola on the beam cone (62% of spins fit, otherwise the bin centre).
  - Result: no predictor beats the training mean out of sample (season-out R² ≤ 0.03). This covers Y, Newell coupling, AE, SYM-H and Φ.
  - Promising hints did not replicate. In 2011 to 2016 there was Newell coupling over 0 to 30 min (ρ +0.41) and AE 3 to 6 h before the beam (ρ +0.43), but neither held on 2017 to 2019 or 2022 to 2024.
  - δE carries no physics. It is a median 1.6 bin widths, no beam is wider than 3, and it shows no trend with any covariate.
  - Energy statistics against Φ need a cut that every beam passes whatever its potential. `fit_ok` (E_b ≥ 2Φ) removes the ~57 eV bin once Φ > 29 V, which held 44% of low-Φ beams and 0% of high-Φ ones. Adding Φ to a bin-pinned energy also correlates corrected E_b with Φ by construction. Together these produced a spurious E_b–Φ correlation (ρ 0.7, season-out R² 0.5) that disappears when only measured E_b ≥ 99 eV and Φ ≤ 45 V are kept.

### Plotting and utility scripts

- `sanity_check.py` - smoke test that loads one day of probe A data and renders a test PNG
- `plot_single_probe.py` - CLI to plot one probe's energy flux for a given date
- `plot_multi_probe.py` - CLI to plot all probes stacked for a given date
- `run_beam_pipeline.py` - CLI to run the full per-timestep pipeline with configurable thresholds. Each run writes the per-beam CSV and the histogram / threshold-comparison PNGs unless `--no-plots` is set.
- `plot_eb_maps.py` - E_b, δE and E_b/δE maps (Y vs X and |Bx| vs R) plus per-probe histograms with Te, one marker per episode (`--seasons`, `--tier`, `--rebuild`). E_b and E_b/Te are plotted in the true frame (`e_b_corr`, with E_b/Te taken per member then median), and the CSV cache `candidates/all/eb_map_<seasons>_t<tier>.csv` keeps the measured `e_b` beside them. Figures go to `figures/eb_maps_<seasons>_t<tier>.png`. The cache rebuilds when any source catalog, episode, position or potential file is newer.
- `app.py` - tkinter front end for the pipeline, spectrogram, and candidate CLIs (`python scripts/app.py`, no extra deps). Pick a probe and interval by calendar + UT times or from the `events_*.csv` lists, set params from the `strict` / `survey` presets or save your own to `configs/ui_presets.json`, toggle the diagnose dump and threshold-compare sweep, then run. The script log streams into the window and the output figures (wheel to zoom, drag to pan), beams table, and diagnose text load in-app, with a browser for earlier runs in `figures/`. Defaults come from `ClassifierParams` so they can't drift from the code.

Command flags:
        Flag	Type	Default	Description
        --probe	choice a-e	a	THEMIS probe
        --trange	2 args	2019-05-01 2019-05-02	start/end times, hour syntax ok (2019-05-01/06:00)
        --hours	2 ints	none	hour window on start date, e.g. 6 12, end rolls to next day if <= start
        --energy-cutoff	float	30.0	low-energy cutoff (eV)
        --min-consecutive	int	1	min consecutive beam steps to keep
        --asym-threshold	float	0.2	asymmetry score-ramp anchor
        --width-threshold	float	0.8	width score-ramp anchor
        --p2o-threshold	float	1.3	para-to-omni score-ramp anchor
        --score-threshold	float	0.6	beam score threshold
        --min-coverage	float	0.01	min PA cone solid-angle coverage
        --n-sigma-lo	float	1.5	lo poisson bar, runs form here, needs a neighbor
        --n-sigma-hi	float	2.5	hi poisson bar, isolated runs stand alone
        --beam-e-max	float	7000.0	band-energy ceiling for any accepted run (eV)
        --coherent-asym-min	float	0.2	per-bin |asym| threshold for coherent run
        --coherent-dir-min	float	1.2	per-bin dominant-cone/perp (R) threshold
        --coherent-min-bins	int	2	min adjacent bins for coherent beam
        --peak-prom-min	float	0.3	log10 prominence for spectral line score
        --peak-width-max	float	4.0	max line FWHM in bins
        --threshold-compare-values	floats	1.0 1.2 1.5 2.0	R values for the threshold-comparison plot
        --no-plots	flag	off	skip plotting
        --diagnose	0 or 2 args	none	dump per-bin spectra/features, UT window e.g. 06:00 07:00, no args = use --hours/trange window

## Optimization

Candidate generation (`generate_candidates.py`) is the slow step. Serial timings are measured on a fixed 12-interval sample of `events_tail2016.csv` (every tenth row, 129 of the season's 1,095 probe-hours, probes a/d/e) with the data already on disk, and the serial season column scales the sample by probe-hours. The 12-worker column is measured on full-season runs. Every change was checked by a byte-for-byte diff of the jsonl output: sample runs against a run of the unchanged code, and season runs against the existing `candidates/tail2016` files, which match that run. Measured 2026-10-25 and 2026-10-26 on the 16-core, 32-thread, 63 GB machine.

Each row includes the changes above it:

| Change | Sample, serial | Per interval | 2016 season, serial (est.) | 2016 season, 12 workers |
|---|---|---|---|---|
| Baseline | 1194 s | 99.5 s | 2.8 h | n/a |
| `--no-update` | 1122 s | 93.5 s | 2.6 h | n/a |
| `OMP_NUM_THREADS=1` | 1110 s | 92.5 s | 2.6 h | 22.9 min |
| FGM reuse | 958 s | 79.8 s | 2.3 h | not run |
| Vectorized `compute_pa_spectra` | 368 s | 30.7 s | 52 min | 6.3 min |
| ESA `varnames` | 346 s | 28.8 s | 49 min | **5.2 min** |

With every change in, a full 2016 season is about 32× faster than the original serial baseline, and 12 workers give 9.4× over serial with the same code. Every 12-worker season run was byte-identical (114/114). On the 12-interval sample, 12 workers took 205 s against 1110 s serial (5.4×), less than the season's 6.9×, because 12 intervals on 12 workers is a single wave, bounded by the slowest interval. A second-season check, 12 intervals spread across 2013 run with every change in, is byte-identical to `candidates/tail2013` (12/12).

Running a season:

```
python scripts/generate_candidates.py --events events_tail<year>.csv --datatype peir --prefetch
python scripts/generate_candidates.py --events events_tail<year>.csv --profile survey --datatype peir --no-cutouts --no-update --workers 12 --out candidates/tail<year>
```

Where an interval's time went before these changes (cProfile, one 13 h interval): about 60% in `compute_pa_spectra`, about 25% in pyspedas loading and time clipping, about 8% in `extract_features`, and 4% in server checks. Unprofiled, the old `compute_pa_spectra` alone took 89 to 96 s on that interval, and the vectorized one takes 10 s.

How each change works:

- **Interval checkpoint.** `run_interval` skips any interval whose `.jsonl` already exists. It writes `.jsonl.tmp` and renames it at the end, so a killed run never leaves a partial file that would be skipped on resume.
- **Day filter on the ESD glob.** `load_esd_distribution` keeps only the CDFs whose date falls in the trange, because the glob spans every year on disk and each file's epoch conversion costs about 0.3 s.
- **Download first, then compute (`--prefetch`, `--no-update`).** Without `--no-update`, every interval makes four server round trips (ESD, FGM, ESA, then FGM again for GSM). Each fetches the remote directory index and confirms that the local file is current, and pyspedas sleeps 1 s on each call. `--prefetch` walks the event list once per probe-day span and makes the same three `themis` calls with `downloadonly=True`, then exits. `--no-update` threads `no_update=True` through the four loaders, so they read the local cache and never contact the server. Measured saving is 6.1 s per interval (6%), and outputs are identical. The main purpose is to make parallel workers safe, since two workers fetching the same daily CDF at once can corrupt it. A file the prefetch missed fails as `FileNotFoundError` or "no data" instead of downloading, so check the `[fail]` list after each season. Prefetch on a full cache costs about 6 s per span (75 s for the 12-interval sample).
- **One BLAS thread per process.** `generate_candidates.py` sets `OMP_NUM_THREADS=1` (with `setdefault`, so it can be overridden) before numpy is imported. By default numpy's OpenBLAS starts 24 threads per process, and 12 workers would carry 288 threads on 32 hardware threads, which busy-wait after each call. The pipeline's only BLAS calls are `np.linalg.norm` on 3-component vectors, so serially the change is within noise (1%). It exists to keep the parallel run from oversubscribing. Spawned workers inherit the environment.
- **Parallel intervals (`--workers N`).** Intervals share no state and write to their own files (`{probe}_{start}_{profile}`), so they can run concurrently. With `--workers` above 1, `run_pool` submits the event rows to a `ProcessPoolExecutor`. Each worker returns only a record count, and the parent prints the same `[fail]` lines and tally as the serial loop. Rows are submitted longest first (intervals run 2 to 13 h), so a long interval doesn't start last and hold up the end of the run. `--workers` above 1 refuses to run without `--no-update`, so workers never download. The default of 1 keeps the serial loop unchanged.
- **Processes, not threads.** The loaders call `themis.fgm` and `themis.esa` and then read the result back by name with `get_data`, from pytplot's global variable store. Two threads on the same probe would overwrite each other's variables between the load and the read, and one would silently compute on the wrong interval's field. Each process has its own store, and its own `THM_DATA_DIR`.
- **Memory.** Tail intervals run 2 to 13 h (median 11 h). A 13 h reduced-mode interval is about 16,000 spins. Its `eflux` array, (time, 32, 88) float32, is about 180 MB, and `onecount`, the same shape in float64, is about 360 MB. With every change in, measured peak is about 2 GB for one process and 13.0 GB for the whole 12-worker tree.
- **Why not 12×.** Per-interval time grows under load: the slowest sample interval took about 1.7× its serial time with 12 workers running. That points at shared memory bandwidth and lower all-core clocks, not at the scheduling.
- **spacepy rc race.** pyspedas imports spacepy, and spacepy 0.7.0 drops the `[spacepy]` header from `~/.spacepy/spacepy.rc` whenever it adds missing defaults (a slice assignment with `nextsec = None` replaces the whole file). After that, every import treats the file as bad, truncates it and rewrites it. Serially that's harmless. With 12 workers importing at once, one reads the file mid-truncation, crashes with `UnboundLocalError`, and the whole pool fails with `BrokenProcessPool`. With the header in place and every default key present, spacepy writes nothing, so the fix is a `[spacepy]` line at the top of the rc. `--workers` above 1 checks for the header and refuses to start without it. A spacepy upgrade that adds a new default key will strip it again.
- **FGM reuse.** `load_bfield_dsl` already loads the L2 FGM file into tplot, GSM variables included. `run_interval` now calls `load_bfield_gsm(..., reload=False)`, which reads the GSM variable from tplot instead of reading the CDF and time clipping again. The ESA moments load in between touches no FGM variable, and a check on two intervals showed the reused arrays are identical to a fresh reload. The default stays `reload=True`, so `run_beam_pipeline.py --diagnose` is unchanged. Measured saving is 12.7 s per interval (14%), more than the 5 s a standalone reload takes. The extra probably comes from reloading into a tplot store that already holds the interval's ESD, FGM and ESA variables, but that hasn't been isolated.
- **Vectorized `compute_pa_spectra`.** The old version looped in Python over every (timestep, energy) pair, about 520k iterations of about 30 small numpy calls on a 13 h interval, so call overhead dominated. The new one works in 512-step time chunks on (time × energy) rows of angle bins: pitch angles and every cone mask come from one pass, and all sums go through `_masked_sums`. Bit-identical output was a requirement, and a plain zero-filled masked sum breaks it, because numpy's pairwise summation then groups the adds differently and shifts the last bit. `_masked_sums` therefore packs each row's selected elements into a contiguous prefix and reduces equal-count rows together, which reproduces `np.sum` on the compressed slice exactly. That was tested on float32 and float64 for every row length from 1 to 176. Every expression also keeps the old operand dtypes (`eflux`, `theta`, `phi` and `domega` are float32, `onecount` is float64), and |B| stays a per-row `np.linalg.norm`. Checked against the old loop on four reduced-mode intervals (2011 and 2016, probes a/d/e) and one fast-survey (176-bin) interval: every output array is bit-identical, and the function is 9 to 10× faster. Temporaries are about 300 MB per chunk. At that stage the 12-worker peak was 14.6 GB, against 14.1 GB before. All callers of `compute_pa_spectra` get the speedup, including the cutout and deck renderers.
- **Only the needed ESA variables.** Without `varnames`, `themis.esa(level="l2")` loads and time clips all 104 variables in the ESA L2 file. `run_interval` passes `varnames` from the `ESA_VARS` constant, the 7 variables it actually reads: `peif_density`, `peif_velocity_gsm` and its fallback `peif_velocity_dsl`, `peif_avgtemp`, `peif_vthermal`, and `peer_density` and `peer_avgtemp` for ne/ni and Te. `load_moments` defaults to `varnames=None` (everything), so its other callers are unchanged. pyspedas does no ESA-specific post-processing, so `varnames` changes how much is loaded, not the values. On five intervals (2011, 2015, 2016), restricted and full loads give bit-identical arrays for all 7 variables and an identical moments dict. The load drops from 5.5 s to 0.6 s. In the pipeline the measured saving is smaller: 1.8 s per interval serially (6%), and 17% on the 12-worker season, where loads compete for I/O. A variable that isn't in `ESA_VARS` comes back as `None`, so add it there before `run_interval` reads it.

Caveats for the parallel run:

- Cutouts are rendered after the jsonl rename, so a worker killed mid-cutouts leaves an interval that resumes as skipped with an incomplete PNG set. Season runs use `--no-cutouts`, so they aren't affected.
- Ctrl-C on Windows can leave orphan workers. The checkpoint makes a restart safe.
- Killed workers leave `.jsonl.tmp` files behind. They are harmless and get overwritten on the rerun.

## Data and event lists

- `events_artemyev2020.csv` - literature event list (Artemyev+ 2020 Table S1) used as the labeling batch source
- `events_tail2015.csv` - the 2015 tail-season scan: 246 tail-region intervals (x_GSM < -8 Re, |y| < |x|), 2,908 probe-hours, probes a/d/e
- `events_tail<year>.csv` - the other seasons, same columns, written by `scan_tail_intervals.py`
- `events_2015.csv`, `events_smoketest.csv` - supporting batches
- `candidates/` - candidate JSONLs + cutouts for the literature batch, `LABELING.md`, frozen label snapshots, `analysis/` outputs (model report, PR curve, catalog, sheets)
- `candidates/tail2015/` - the 2015 batch and its `analysis/` (episode catalog, blind deck, coldline census)
- `candidates/tail<year>/` - the other seasons, same layout; `candidates/all/` - pooled catalog, episodes, coverage, and the `eb_map_*.csv` caches
- ESA L2 files (`data/th?/l2/esa/`) also supply the EFI-derived spacecraft potential (`peir_sc_pot`) and the electron moments (Te, ne) used by the potential correction
- `paper/skeleton.md` - working paper outline (JGR Space Physics target)
- `docs/VARIABLES.md` - reference sheet for every variable: detector features, candidate record fields, model transforms, episode features, catalog columns

## Tests

- `test_smoke.py` - basic import and unit tests for the probe/variable helpers

## Credits

Built on [pyspedas](https://github.com/spedas/pyspedas). THEMIS ESA data courtesy of the THEMIS mission (NASA) and the instrument teams.

Ad Astra per Aspera

## License

MIT (see [LICENSE](LICENSE))
