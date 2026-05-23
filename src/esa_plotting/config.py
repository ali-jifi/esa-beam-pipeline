import os
from pathlib import Path

DEFAULT_ZRANGE = (1e3, 1e8)
DEFAULT_YRANGE = (5, 30000)
DEFAULT_COLORMAP = "spedas"


def set_data_dir(path: str | os.PathLike | None = None) -> str:
    """Resolve and export THM_DATA_DIR. Falls back to env var, then ./data."""
    if path is None:
        path = os.environ.get("THM_DATA_DIR")
    if path is None:
        path = Path(__file__).resolve().parents[2] / "data"
    path = str(Path(path).resolve())
    os.environ["THM_DATA_DIR"] = path
    return path
