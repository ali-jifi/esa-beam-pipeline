from typing import Iterable

import pyspedas
from pyspedas.projects import themis


def load_esa(
    probes: str | Iterable[str],
    trange: tuple[str, str] | list[str],
    level: str = "l2",
    time_clip: bool = True,
) -> list[str]:
    """Load THEMIS ESA data and return the tplot variables created."""
    return themis.esa(
        probe=list(probes) if not isinstance(probes, str) else probes,
        trange=list(trange),
        level=level,
        notplot=False,
        time_clip=time_clip,
    )


def pyspedas_version() -> str:
    return pyspedas.version()
