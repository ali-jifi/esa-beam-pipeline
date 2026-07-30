# run ion beam detection pipeline
import argparse
from pathlib import Path

from esa_plotting.config import set_data_dir
from esa_plotting.beam_pipeline import (run_pipeline, ClassifierParams, diagnose_window, apply_hours, load_state_gsm, load_bfield_gsm)

FIGURES = Path(__file__).resolve().parents[1] / "figures"

def main() -> None:
    p = argparse.ArgumentParser(description="THEMIS ion beam detection pipeline")
    p.add_argument("--probe", default="a", choices=list("abcde"))
    p.add_argument("--trange", nargs=2, default=["2019-05-01", "2019-05-02"],
                   help="Start and end times, e.g. 2019-05-01 2019-05-02, "
                        "hour syntax also works, e.g. 2019-05-01/06:00")
    p.add_argument("--hours", nargs=2, type=int, metavar=("H_START", "H_END"),
                   help="Hour window on the start date, e.g. 6 12; "
                        "end rolls to next day if <= start")
    p.add_argument("--energy-cutoff", type=float, default=30.0,
                   help="Low-energy cutoff in eV (default: 30)")
    p.add_argument("--min-consecutive", type=int, default=1,
                   help="Min consecutive beam steps to keep, 1 = keep isolated (default: 1)")
    p.add_argument("--asym-threshold", type=float, default=0.2,
                   help="Asymmetry threshold (default: 0.2)")
    p.add_argument("--width-threshold", type=float, default=0.8,
                   help="Width threshold (default: 0.8)")
    p.add_argument("--p2o-threshold", type=float, default=1.3,
                   help="Para-to-omni ratio threshold (default: 1.3)")
    p.add_argument("--score-threshold", type=float, default=0.6,
                   help="Beam score threshold (default: 0.6)")
    p.add_argument("--min-coverage", type=float, default=0.01,
                   help="Min PA cone solid-angle coverage (default: 0.01)")
    p.add_argument("--n-sigma-lo", type=float, default=1.5,
                   help="Lo poisson bar, runs form here, needs a neighbor (default: 1.5)")
    p.add_argument("--n-sigma-hi", type=float, default=2.5,
                   help="Hi poisson bar, isolated runs stand alone here (default: 2.5)")
    p.add_argument("--beam-e-max", type=float, default=7000.0,
                   help="Band-energy ceiling for any accepted run (default: 7000)")
    p.add_argument("--coherent-asym-min", type=float, default=0.2,
                   help="Per-bin |asym| threshold for coherent run (default: 0.2)")
    p.add_argument("--coherent-dir-min", type=float, default=1.2,
                   help="Per-bin dominant-cone/omni threshold (default: 1.2)")
    p.add_argument("--coherent-min-bins", type=int, default=2,
                   help="Min adjacent bins to count as a coherent beam (default: 2)")
    p.add_argument("--peak-prom-min", type=float, default=0.3,
                   help="Log10 prominence threshold for spectral line score (default: 0.3)")
    p.add_argument("--peak-width-max", type=float, default=4.0,
                   help="Max line FWHM in bins (default: 4)")
    p.add_argument("--threshold-compare-values", nargs="+", type=float,
                   default=[1.0, 1.2, 1.5, 2.0],
                   help="R (coherent_dir_min) values for the threshold-comparison plot")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--diagnose", nargs="*", metavar="UT",
                   help="Dump per-bin spectra and features. Give a UT window, "
                        "e.g. 06:00 07:00, or no args to use the --hours/trange window")
    args = p.parse_args()

    trange = apply_hours(args.trange, args.hours)
    data_dir = set_data_dir()

    params = ClassifierParams(
        asymmetry_min=args.asym_threshold,
        width_max=args.width_threshold,
        para_to_omni_min=args.p2o_threshold,
        score_threshold=args.score_threshold,
        min_coverage=args.min_coverage,
        n_sigma_lo=args.n_sigma_lo,
        n_sigma_hi=args.n_sigma_hi,
        beam_e_max=args.beam_e_max,
        coherent_asym_min=args.coherent_asym_min,
        coherent_dir_min=args.coherent_dir_min,
        coherent_min_bins=args.coherent_min_bins,
        peak_prom_min=args.peak_prom_min,
        peak_width_max=args.peak_width_max,
    )

    result = run_pipeline(
        probe=args.probe,
        trange=trange,
        data_dir=data_dir,
        params=params,
        min_consecutive=args.min_consecutive,
        energy_cutoff_low=args.energy_cutoff,
        figures_dir=str(FIGURES) if not args.no_plots else None,
        threshold_compare_values=tuple(args.threshold_compare_values),
    )

    n = len(result.features.times)
    n_beam = result.classification_smoothed.is_beam.sum()
    print(f"\n=== Summary ===")
    print(f"Total timesteps: {n}")
    print(f"Beam timesteps (smoothed): {n_beam} ({100*n_beam/max(n,1):.1f}%)")

    if args.diagnose is not None:
        if len(args.diagnose) == 2:
            ut0, ut1 = args.diagnose
        elif len(args.diagnose) == 0:
            # no args, inherit the hour window, full day if none
            ut0 = trange[0].split("/")[1] if "/" in trange[0] else "00:00"
            ut1 = trange[1].split("/")[1] if "/" in trange[1] else "23:59:59"
        else:
            p.error("--diagnose takes 0 or 2 args")
        # dated diagnose file in the same per-day folder as the figures
        date_str = trange[0].split("/")[0]
        diag_dir = FIGURES / date_str
        diag_dir.mkdir(parents=True, exist_ok=True)
        # hour tag keeps sub-day windows from clobbering each other
        tag = "-".join(s.split("/")[1].replace(":", "") for s in trange if "/" in s)
        tag = f"_{tag}" if tag else ""
        diag_path = diag_dir / f"th{args.probe}_beam_{date_str}{tag}_diagnose.txt"
        # gsm position and b-field for context, skip if load fails
        pos_gsm = b_gsm = None
        try:
            pos_gsm = load_state_gsm(args.probe, trange, data_dir)
            b_gsm = load_bfield_gsm(args.probe, trange, data_dir)
        except Exception as e:
            print(f"[diag] gsm context unavailable: {e}")
        diagnose_window(result.spectra, result.features, result.classification,
                        params, ut0, ut1,
                        out_path=str(diag_path),
                        pos_gsm=pos_gsm, b_gsm=b_gsm)

    print(f"\n=== Threshold Sensitivity ===")
    for param, info in result.sensitivity.items():
        print(f"  {param}: {list(zip(info['values'], info['beam_counts']))}")

if __name__ == "__main__":
    main()