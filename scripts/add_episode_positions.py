# joins gsm spacecraft position onto episodes as a sidecar
# separate file on purpose so frozen candidate records and labels stay untouched
# reads analysis/episodes.jsonl, writes analysis/episode_positions.jsonl keyed by episode_id
import json
import logging
import os
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.WARNING)

ROOT = Path(__file__).resolve().parents[1]
ANA = ROOT / "candidates" / "tail2015" / "analysis"
os.environ["THM_DATA_DIR"] = str(ROOT / "data")

from pyspedas import get_data  # noqa: E402
from pyspedas.projects import themis  # noqa: E402

RE_KM = 6371.2
PAD_S = 1800.0
# known 12-31 triple conjunction, sanity anchor for the join
CHECK_DAY = "2015-12-31"
CHECK_WIN = ("08:25", "08:45")


def load_pos(probe, t0, t1):
    # state cadence is ~60 s so pad wide and let interp cover the episode span
    a = datetime.fromtimestamp(t0 - PAD_S, tz=timezone.utc)
    b = datetime.fromtimestamp(t1 + PAD_S, tz=timezone.utc)
    trange = [a.strftime("%Y-%m-%d/%H:%M"), b.strftime("%Y-%m-%d/%H:%M")]
    var = f"th{probe}_pos_gsm"
    themis.state(probe=probe, trange=trange, time_clip=True,
                 varnames=[var], no_update=True)
    d = get_data(var)
    return (d.times, d.y / RE_KM) if d is not None else None


def main():
    eps = [json.loads(line) for line in open(ANA / "episodes.jsonl")]
    for e in eps:
        e["_t0"] = datetime.fromisoformat(e["t_start_ut"]).timestamp()
        e["_t1"] = datetime.fromisoformat(e["t_end_ut"]).timestamp()

    groups = {}
    for e in eps:
        groups.setdefault((e["probe"], e["t_start_ut"][:10]), []).append(e)

    out, missing = [], []
    for (probe, day), grp in sorted(groups.items()):
        t0 = min(e["_t0"] for e in grp)
        t1 = max(e["_t1"] for e in grp)
        try:
            pos = load_pos(probe, t0, t1)
        except Exception as ex:
            pos = None
            print(f"[warn] state load failed {probe} {day}: {ex}")
        if pos is None:
            missing.extend(e["episode_id"] for e in grp)
            continue
        pt, py = pos
        for e in grp:
            tm = 0.5 * (e["_t0"] + e["_t1"])
            # np.interp clamps silently so reject out-of-range explicitly
            if tm < pt[0] or tm > pt[-1]:
                missing.append(e["episode_id"])
                continue
            x, y, z = (float(np.interp(tm, pt, py[:, k])) for k in range(3))
            out.append({
                "episode_id": e["episode_id"],
                "probe": probe,
                "t_mid_ut": datetime.fromtimestamp(tm, tz=timezone.utc).isoformat(),
                "t_mid": round(tm, 3),
                "x_gsm": round(x, 3),
                "y_gsm": round(y, 3),
                "z_gsm": round(z, 3),
                "r_gsm": round(float(np.sqrt(x * x + y * y + z * z)), 3),
            })

    with open(ANA / "episode_positions.jsonl", "w") as fh:
        for r in out:
            fh.write(json.dumps(r) + "\n")

    print(f"positions joined: {len(out)}/{len(eps)} episodes "
          f"({len(groups)} probe-days)")
    if missing:
        print(f"missing: {len(missing)} -> {missing[:5]}")

    chk = sorted((r for r in out
                  if r["t_mid_ut"][:10] == CHECK_DAY
                  and CHECK_WIN[0] <= r["t_mid_ut"][11:16] <= CHECK_WIN[1]),
                 key=lambda r: r["t_mid_ut"])
    print(f"\nanchor check, {CHECK_DAY} {CHECK_WIN[0]}-{CHECK_WIN[1]}:")
    for r in chk:
        print(f"  {r['episode_id']}  th{r['probe']}  {r['t_mid_ut'][11:19]}  "
              f"x={r['x_gsm']:>7.2f} y={r['y_gsm']:>7.2f} z={r['z_gsm']:>6.2f} "
              f"r={r['r_gsm']:.2f}")
    if len(chk) > 1:
        xs = [r["x_gsm"] for r in chk]
        ys = [r["y_gsm"] for r in chk]
        print(f"  spread: dX={max(xs)-min(xs):.2f} RE  dY={max(ys)-min(ys):.2f} RE")
    print(f"[ok] wrote {ANA / 'episode_positions.jsonl'}")


if __name__ == "__main__":
    main()
