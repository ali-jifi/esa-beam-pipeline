"""Multi-probe ESA ion eflux stack plot. Example: python scripts/plot_multi_probe.py --date 2009-02-27"""
import argparse
from datetime import datetime, timedelta
from pathlib import Path

from esa_plotting import PROBES, configure_eflux_panel, eflux_var, load_esa, set_data_dir
from esa_plotting.plotting import stack_plot

FIGURES = Path(__file__).resolve().parents[1] / "figures"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--date", default="2009-02-27", help="YYYY-MM-DD")
    p.add_argument("--probes", default="abcde", help="probe letters concatenated, e.g. abcde")
    args = p.parse_args()

    set_data_dir()
    probes = list(args.probes)
    start = args.date
    end = (datetime.strptime(start, "%Y-%m-%d") + timedelta(days=1)).strftime("%Y-%m-%d")
    trange = [start, end]

    load_esa(probes, trange)

    vars_to_plot = [eflux_var(pr) for pr in probes]
    for v in vars_to_plot:
        configure_eflux_panel(v)

    FIGURES.mkdir(exist_ok=True)
    out = FIGURES / f"multi_probe_eflux_{start.replace('-', '')}"
    plotted = stack_plot(vars_to_plot, title=f"THEMIS ESA ion eflux {start}", out_png=str(out))
    print(f"[OK] wrote {out}.png with {len(plotted)} panels")


if __name__ == "__main__":
    main()
