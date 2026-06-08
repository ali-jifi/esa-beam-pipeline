# Data Analysis Project #1

Detect and characterize field-aligned ion beams in THEMIS ESA particle data by reducing 3D distributions into pitch-angle-resolved energy spectra and classifying directional, narrow-line beam signatures against isotropic plasma-sheet background.

Built by Ali Jifi

# esa_plotting

Python package for plotting and analyzing THEMIS satellite Electrostatic Analyzer (ESA) data, with an ion beam detection pipeline built on top of pyspedas.

## Core package - `src/esa_plotting/`

- **Config.py** - sets up where satellite data lives on disk. `set_data_dir()` checks for a path passed in, then falls back to the `THM_DATA_DIR` env var, then defaults to `./data`. Also defines plot defaults (color range, energy range, color map).
- **Probes.py** - THEMIS has 5 probes labeled a-e. Defines that tuple and a helper `eflux_var()` that builds the tplot variable name for a given probe/species combo, such as `tha_peif_en_eflux`.
- **Loaders.py** - thin wrappers around pyspedas to load three types of THEMIS data:
  - `load_esa()` - energy flux spectra (the main data product)
  - `load_esd()` - 3D particle distributions
  - `load_fgm()` - magnetic field data
  
  All return tplot variable names that pyspedas stores internally.
- **Plotting.py** - two plotting helpers:
  - `configure_eflux_panel()` - sets up a single energy flux spectrogram panel (log scales, color range, labels)
  - `stack_plot()` - takes multiple tplot variables, filters out ones that didn't load, and renders them stacked into a PNG
- **__init__.py** - re-exports the public api from all of the above.

## Beam detection pipeline - `beam_pipeline.py`

Detects ion beams in THEMIS data through 6 phases.

### Phase 0 - Data acquisition

Loads three data products for a given probe/time range:

- `load_esd_distribution()` - reads the CDF directly (not through pyspedas) to get the full 3D ion distribution including the angle/energy lookup tables that pyspedas doesn't expose. Parses time ranges, clips to the requested interval, and extracts the dominant energy/angle mode.
- `load_bfield_dsl()` - loads magnetic field in DSL (Despun Spacecraft L) coordinates, trying survey/low/high cadence in order.
- `load_moments()` - loads bulk plasma parameters (density, velocity, temperature).

### Phase 1 - Spectral reduction

Takes the full 3D ion distribution (flux at 32 energies x 176 angle bins) and collapses it into three 1D energy spectra based on pitch angle relative to the magnetic field.

Pitch angle is the angle between a particle's velocity and the local B-field direction. 0° = moving along B (field-aligned/parallel), 180° = moving opposite to B (anti-parallel), 90° = perpendicular.

The three spectra:

- **Omnidirectional** - average flux at each energy across all angles, weighted by solid angle (domega). This is what a standard spectrogram shows.
- **Parallel (0-30°)** - only bins where the particle velocity is roughly field-aligned. A beam streaming along B shows up here.
- **Anti-parallel (150-180°)** - only bins moving opposite to B. A beam coming from the other direction shows up here.

A beam is a narrow, directional population, it lights up in one PA gate but not the other. Plasma sheet ions are roughly isotropic, so all three curves overlap. Comparing parallel and anti-parallel spectra surfaces the directional asymmetry that defines a beam.

`compute_pa_spectra()`: for each timestep, interpolates the magnetic field to the distribution time, computes pitch angles for every angle bin, then sorts flux into omni / para / anti energy spectra. `_compute_pitch_angles()` does the geometry by converting instrument look directions to particle velocities (opposite direction), then dotting with the B-field unit vector.

### Phase 2 - Feature extraction (`extract_features()`)

For each timestep, computes the discriminating features.

Spectral + moment features:

- `e_peak` - energy of peak flux
- `width` - how broad the spectrum is (narrow = beam-like)
- `asymmetry` - `(para - anti) / (para + anti)` near the peak, picked from the bin with max `|asym|` across a window. Positive = field-aligned beam.
- `para_to_omni` - how much stronger the parallel flux is vs omnidirectional
- `energy_ratio` - bulk flow energy / thermal energy from moments (high = directed flow)

Coherent-run + spectral-line features:

- `coherent_ok` - bool, a real coherent directional run was found (both cones sampled, enough adjacent bins where `|asym|` and dominant-cone/omni clear their per-bin thresholds)
- `peak_prom` - log10 prominence of the narrow spectral line found inside the coherent run (0.3 = 2x above local baseline)
- `peak_width` - FWHM of that line in bins
- `e_line` - energy of the line in eV

Spectral-line detection is local to the coherent run, not global. It scans only the dominant cone (para if `asymmetry >= 0`, else anti), compresses to finite/positive bins, takes `log10(flux)`, and runs `scipy.signal.find_peaks` with a bounded prominence window (`peak_wlen=5`) and a width cap (`peak_width_max=4.0` bins). A peak only counts if it sits inside the directional run band (±1 bin slop). The idea: a beam = the directional region is also a narrow line; a prominent line elsewhere (e.g. the anti-parallel plasma-sheet peak) is rejected.

### Phase 3 - Classification (`classify_beams()`)

Two-path heuristic classifier, wrapped in an AND-gate:

1. **Score-based** - weighted sum of normalized feature scores. Weights: `w_asymmetry=0.35`, `w_width=0.25`, `w_para_to_omni=0.25`, `w_peak_prom=0.15` (peak prominence took the slot energy_ratio vacated). `w_energy_ratio=0.0`. The moments energy_ratio score is dead code, kept for now. Beam candidate if score clears the threshold.
2. **Hard rule fallback** - beam candidate if asymmetry exceeds threshold AND either width is narrow enough or para_to_omni is high enough. Catches strong beams that might miss the score threshold.

**AND-gate:** final `is_beam = (score_ok or hard_ok) and gate`, where `gate = coherent_ok and peak_prom >= peak_prom_min`. Both detectors, the directional coherent run AND a narrow spectral line at the same energy, must agree, which kills noise that fires only one signal alone.

Line-detection params on `ClassifierParams`: `peak_prom_min=0.3` (log10, so 0.3 = 2x above local baseline), `peak_width_max=4.0` (FWHM cap in bins), `peak_wlen=5` (local prominence window).

Beam direction is tagged (+1 parallel, -1 anti-parallel) from the asymmetry sign.

### Phase 4 - Temporal smoothing

- `smooth_labels()` - requires N consecutive beam-flagged timesteps to keep a beam interval. Default `min_consecutive` is 1 (keep isolated beams) because the AND-gate already enforces precision; raise it to suppress more.
- `threshold_sensitivity()` - sweeps each threshold parameter and reports how beam count changes, so you can see how stable the classification is.

### Phase 5 - Plotting

- `plot_feature_timeseries()` - 7-panel overview with spectrogram, features, beam score, and a color bar showing classification (red = parallel beam, blue = anti-parallel, orange = unknown direction, gray = no beam). Takes `ClassifierParams` so every threshold guide-line is driven by the actual params instead of hardcoded values. Panel 4 plots **Peak Prominence** (the AND-gate line) with its threshold, replacing the old E_flow/E_th panel; the score panel has a threshold line.
- `plot_curated_snapshots()` - picks representative timesteps (confirmed beams, plasma sheet, borderline cases) and plots the three-curve energy spectra at each.
- `diagnose_window()` - dumps per-timestep spectra and features, including `peak_prom`, `peak_width`, `e_line`, and `coherent_ok`.

`run_pipeline()` ties it all together via loading data, runs the phases, optionally saves plots, returns everything in a `PipelineResult` dataclass. Default `min_consecutive` is 1, threads `peak_width_max` / `peak_wlen` through to `extract_features`, and passes `params` into the plotter.

## Scripts - `scripts/`

- `sanity_check.py` - smoke test that loads one day of probe A data and renders a test PNG
- `plot_single_probe.py` - CLI to plot one probe's energy flux for a given date
- `plot_multi_probe.py` - CLI to plot all probes stacked for a given date
- `run_beam_pipeline.py` - CLI to run the full beam detection pipeline with configurable thresholds. Flags include `--peak-prom-min` (default 0.3), `--peak-width-max` (default 4.0), and `--min-consecutive` (default 1, where 1 = keep isolated).
Command flags:
        Flag	Type	Default	Description
        --probe	choice a-e	a	THEMIS probe
        --trange	2 args	2019-05-01 2019-05-02	start/end times
        --energy-cutoff	float	30.0	low-energy cutoff (eV)
        --min-consecutive	int	1	min consecutive beam steps to keep
        --asym-threshold	float	0.2	asymmetry threshold
        --width-threshold	float	0.8	width threshold
        --p2o-threshold	float	1.3	para-to-omni ratio threshold
        --score-threshold	float	0.4	beam score threshold
        --min-coverage	float	0.01	min PA cone solid-angle coverage
        --beam-flux-floor	float	0.1	min omni flux frac of peak for asym scan
        --coherent-asym-min	float	0.2	per-bin |asym| threshold for coherent run
        --coherent-dir-min	float	1.2	per-bin dominant-cone/omni threshold
        --coherent-min-bins	int	2	min adjacent bins for coherent beam
        --peak-prom-min	float	0.3	log10 prominence for spectral line score
        --peak-width-max	float	4.0	max line FWHM in bins
        --no-plots	flag	off	skip plotting
        --diagnose	2 args	none	dump per-bin spectra/features for UT window, e.g. 06:00 07:00

## Tests

- `test_smoke.py` - basic import and unit tests for the probe/variable helpers

## Credits

Built on [pyspedas](https://github.com/spedas/pyspedas). THEMIS ESA data courtesy of the THEMIS mission (NASA) and the instrument teams.

Ad astra per aspera

## License

MIT (see [LICENSE](LICENSE))
