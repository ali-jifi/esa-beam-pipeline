"""Single-probe ESA ion eflux plot. Example: python scripts/plot_single_probe.py --probe a --date 2019-05-01"""
import argparse
from datetime import datetime, timedelta
from pathlib import Path

from esa_plotting import configure_eflux_panel, eflux_var, load_esa, set_data_dir
from pyspedas import tplot, tplot_options

FIGURES = Path(__file__).resolve().parents[1] / "figures"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--probe", default="a", choices=list("abcde"))
    p.add_argument("--date", default="2019-05-01", help="YYYY-MM-DD")
    args = p.parse_args()

    set_data_dir()
    start = args.date
    end = (datetime.strptime(start, "%Y-%m-%d") + timedelta(days=1)).strftime("%Y-%m-%d")
    trange = [start, end]

    load_esa(args.probe, trange)
    v = eflux_var(args.probe)

    configure_eflux_panel(v, ytitle=f"ESA i+ th{args.probe}\nE [eV]")
    tplot_options("title", f"THEMIS-{args.probe.upper()} ESA ion eflux {start}")

    FIGURES.mkdir(exist_ok=True)
    out = FIGURES / f"th{args.probe}_peif_eflux_{start.replace('-', '')}"
    tplot(v, save_png=str(out), xsize=12, ysize=3, display=False)
    print(f"[OK] wrote {out}.png")


if __name__ == "__main__":
    main()
