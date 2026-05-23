from typing import Iterable

from pyspedas import data_exists, options, tplot, tplot_options

from esa_plotting.config import DEFAULT_COLORMAP, DEFAULT_YRANGE, DEFAULT_ZRANGE


def configure_eflux_panel(
    var: str,
    zrange: tuple[float, float] = DEFAULT_ZRANGE,
    yrange: tuple[float, float] = DEFAULT_YRANGE,
    colormap: str = DEFAULT_COLORMAP,
    ytitle: str | None = None,
    ztitle: str = "eflux\n[eV/cm^2-s-sr-eV]",
) -> None:
    options(var, "spec", 1)
    options(var, "ylog", 1)
    options(var, "zlog", 1)
    options(var, "yrange", list(yrange))
    options(var, "zrange", list(zrange))
    options(var, "Colormap", colormap)
    options(var, "ztitle", ztitle)
    if ytitle is None:
        probe = var.split("_")[0].upper()
        ytitle = f"{probe}\nEnergy [eV]"
    options(var, "ytitle", ytitle)


def stack_plot(
    variables: Iterable[str],
    title: str,
    out_png: str,
    xsize: int = 12,
    ysize: int = 3,
    vertical_spacing: float = 0.05,
    display: bool = False,
) -> list[str]:
    """Filter to vars that loaded, set common options, render to PNG. Returns plotted vars."""
    plotted = [v for v in variables if data_exists(v)]
    missing = [v for v in variables if v not in plotted]
    if missing:
        print(f"[WARN] no data for: {missing}")

    tplot_options("title", title)
    tplot_options("vertical_spacing", vertical_spacing)
    tplot(plotted, save_png=out_png, xsize=xsize, ysize=ysize, display=display)
    return plotted
