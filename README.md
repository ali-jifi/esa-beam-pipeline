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

`compute_pa_spectra()`: for each timestep, interpolates the magnetic field to the distribution time, computes pitch angles for every angle bin, then sorts flux into omni / para / anti / perp energy spectra. It also propagates the per-bin one-count level through the same solid-angle weighting to get per-cone noise floors, and propagates Poisson counting uncertainty (`var(flux) = onecount * flux` per bin) into per-cone sigmas (`para_sig` / `anti_sig` / `perp_sig`) used by the significance gates below. **Cone one-count floors**: a sampled-but-empty para/anti cone reads as its one-count upper limit instead of NaN - without this the pipeline is blind to the strongest lobe beams, since a pure anti-parallel beam *should* have an empty parallel cone. `_compute_pitch_angles()` does the geometry; the CDF look-direction angles are particle travel direction (verified against L2 velocity moments), so they dot directly with the B-field unit vector.

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
- `peak_width` - FWHM of that line in bins
- `e_line` - energy of the line in eV (this is **E_beam**)
- `de_line` - **ΔE**, the line FWHM in eV
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
- `build_episode_catalog.py` - episode features, lobe-scoped logistic model, and the two-tier episode catalog. **Hard gates are upstream of the model, never features**: tailward only (`dir_x_bx = -1`), E_b < 600 eV (outflow energy window), beta < 0.03 (lobe scope), Te < 500 eV (closes the corridor leak, and a missing `te_med` is admitted rather than rejected). All read the episode median except direction, which comes from the first member. Catalog cut at model prob 0.82 (blind-validated: contamination 1/14 at this threshold) plus prominence >= 0.55. Constants are the frozen 2026-08-11 closeout values; rerunning refits on current labels.
- `render_episode_sheets.py` - renders episode representative cutouts and tiles contact sheets from a review-queue CSV.
- `render_beam_members.py` - per-member cutouts for every catalog beam, one folder per episode under `analysis/beams/tier{n}/<episode_id>/`, refreshes the rep cutouts in `analysis/cutouts/` on the way.

Outputs land in `candidates/tail2015/analysis/`: `beam_catalog_2015.csv` (tier 1 = 89 episodes / tier 2 = 142), `episodes.jsonl`, the blind-validation deck, review-round queues and sheets, and the coldline census (`coldline_final.csv`, 41 members under a six-constant frozen definition).

### Plotting and utility scripts

- `sanity_check.py` - smoke test that loads one day of probe A data and renders a test PNG
- `plot_single_probe.py` - CLI to plot one probe's energy flux for a given date
- `plot_multi_probe.py` - CLI to plot all probes stacked for a given date
- `run_beam_pipeline.py` - CLI to run the full per-timestep pipeline with configurable thresholds. Each run writes the per-beam CSV and the histogram / threshold-comparison PNGs unless `--no-plots` is set.

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

## Data and event lists

- `events_artemyev2020.csv` - literature event list (Artemyev+ 2020 Table S1) used as the labeling batch source
- `events_tail2015.csv` - the 2015 tail-season scan: 246 tail-region intervals (x_GSM < -8 Re, |y| < |x|), 2,908 probe-hours, probes a/d/e
- `events_2015.csv`, `events_smoketest.csv` - supporting batches
- `candidates/` - candidate JSONLs + cutouts for the literature batch, `LABELING.md`, frozen label snapshots, `analysis/` outputs (model report, PR curve, catalog, sheets)
- `candidates/tail2015/` - the 2015 batch and its `analysis/` (episode catalog, blind deck, coldline census)
- `paper/skeleton.md` - working paper outline (JGR Space Physics target)
- `docs/VARIABLES.md` - reference sheet for every variable: detector features, candidate record fields, model transforms, episode features, catalog columns

## Tests

- `test_smoke.py` - basic import and unit tests for the probe/variable helpers

## Credits

Built on [pyspedas](https://github.com/spedas/pyspedas). THEMIS ESA data courtesy of the THEMIS mission (NASA) and the instrument teams.

Ad Astra per Aspera

## License

MIT (see [LICENSE](LICENSE))
