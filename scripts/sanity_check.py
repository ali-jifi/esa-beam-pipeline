"""Smoke test: imports, env, download, render one PNG. Run from repo root."""
import sys
from pathlib import Path

from esa_plotting import configure_eflux_panel, eflux_var, load_esa, set_data_dir
from esa_plotting.loaders import pyspedas_version
from pyspedas import get_data, tplot, tplot_names, tplot_options

FIGURES = Path(__file__).resolve().parents[1] / "figures"


def check(label, condition, detail=""):
    mark = "OK " if condition else "FAIL"
    print(f"[{mark}] {label}" + (f" ({detail})" if detail else ""))
    if not condition:
        sys.exit(1)


def main() -> None:
    check("imports", True, f"pyspedas {pyspedas_version()}")

    data_dir = set_data_dir()
    check("data dir exists", Path(data_dir).exists(), data_dir)

    probe = "a"
    trange = ["2019-05-01", "2019-05-02"]
    print(f"\nLoading THEMIS-{probe.upper()} ESA L2 for {trange[0]}...")
    loaded = load_esa(probe, trange)
    check("esa() returned variables", bool(loaded), f"{len(loaded)} vars")

    target = eflux_var(probe)
    check(f"{target} in tplot store", target in tplot_names(quiet=True))

    d = get_data(target)
    check(
        f"{target} has data",
        d is not None and len(d.times) > 0,
        f"{len(d.times)} time samples" if d is not None else "no data",
    )

    FIGURES.mkdir(exist_ok=True)
    out = FIGURES / "test_plot_output.png"
    if out.exists():
        out.unlink()

    configure_eflux_panel(target)
    tplot_options("title", f"TEST  THEMIS-{probe.upper()}  {target}  {trange[0]}")
    tplot(target, save_png=str(out.with_suffix("")), xsize=12, ysize=3, display=False)

    check("PNG written", out.exists(), str(out.resolve()))
    check("PNG nonempty", out.stat().st_size > 1000, f"{out.stat().st_size} bytes")
    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
