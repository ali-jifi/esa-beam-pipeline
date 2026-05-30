# run ion beam detection pipeline
import argparse
from pathlib import Path

from esa_plotting.config import set_data_dir
from esa_plotting.beam_pipeline import run_pipeline, ClassifierParams, diagnose_window

FIGURES = Path(__file__).resolve().parents[1] / "figures"


def main() -> None:
    p = argparse.ArgumentParser(description="THEMIS ion beam detection pipeline")
    p.add_argument("--probe", default="a", choices=list("abcde"))
    p.add_argument("--trange", nargs=2, default=["2019-05-01", "2019-05-02"],
                   help="Start and end times, e.g. 2019-05-01 2019-05-02")
    p.add_argument("--energy-cutoff", type=float, default=30.0,
                   help="Low-energy cutoff in eV (default: 30)")
    p.add_argument("--min-consecutive", type=int, default=2,
                   help="Min consecutive beam steps for smoothing (default: 2)")
    p.add_argument("--asym-threshold", type=float, default=0.2,
                   help="Asymmetry threshold (default: 0.2)")
    p.add_argument("--width-threshold", type=float, default=0.8,
                   help="Width threshold (default: 0.8)")
    p.add_argument("--p2o-threshold", type=float, default=1.3,
                   help="Para-to-omni ratio threshold (default: 1.3)")
    p.add_argument("--score-threshold", type=float, default=0.4,
                   help="Beam score threshold (default: 0.4)")
    p.add_argument("--min-coverage", type=float, default=0.01,
                   help="Min PA cone solid-angle coverage (default: 0.01)")
    p.add_argument("--beam-flux-floor", type=float, default=0.1,
                   help="Min omni flux as frac of peak for asym scan (default: 0.1)")
    p.add_argument("--coherent-asym-min", type=float, default=0.2,
                   help="Per-bin |asym| threshold for coherent run (default: 0.2)")
    p.add_argument("--coherent-dir-min", type=float, default=1.2,
                   help="Per-bin dominant-cone/omni threshold (default: 1.2)")
    p.add_argument("--coherent-min-bins", type=int, default=2,
                   help="Min adjacent bins to count as a coherent beam (default: 2)")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--diagnose", nargs=2, metavar=("UT_START", "UT_END"),
                   help="Dump per-bin spectra and features for UT window, e.g. 06:00 07:00")
    args = p.parse_args()

    data_dir = set_data_dir()

    params = ClassifierParams(
        asymmetry_min=args.asym_threshold,
        width_max=args.width_threshold,
        para_to_omni_min=args.p2o_threshold,
        score_threshold=args.score_threshold,
        min_coverage=args.min_coverage,
        beam_flux_floor=args.beam_flux_floor,
        coherent_asym_min=args.coherent_asym_min,
        coherent_dir_min=args.coherent_dir_min,
        coherent_min_bins=args.coherent_min_bins,
    )

    result = run_pipeline(
        probe=args.probe,
        trange=args.trange,
        data_dir=data_dir,
        params=params,
        min_consecutive=args.min_consecutive,
        energy_cutoff_low=args.energy_cutoff,
        figures_dir=str(FIGURES) if not args.no_plots else None,
    )

    n = len(result.features.times)
    n_beam = result.classification_smoothed.is_beam.sum()
    print(f"\n=== Summary ===")
    print(f"Total timesteps: {n}")
    print(f"Beam timesteps (smoothed): {n_beam} ({100*n_beam/max(n,1):.1f}%)")

    if args.diagnose:
        diag_path = Path(__file__).resolve().parents[1] / "figures" / "diagnose.txt"
        diagnose_window(result.spectra, result.features, result.classification,
                        params, args.diagnose[0], args.diagnose[1],
                        out_path=str(diag_path))

    print(f"\n=== Threshold Sensitivity ===")
    for param, info in result.sensitivity.items():
        print(f"  {param}: {list(zip(info['values'], info['beam_counts']))}")


if __name__ == "__main__":
    main()
