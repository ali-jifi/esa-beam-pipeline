# tail-region interval scan, x_gsm < -8 re and |y| < |x|, hourly intervals per probe-day
# reconstructs the scan behind events_tail2015.csv, that one was never committed
import argparse
import csv
import os
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
# pyspedas reads this at import so it has to land first
os.environ["THM_DATA_DIR"] = str(ROOT / "data")

from pyspedas import get_data  # noqa: E402
from pyspedas.projects import themis  # noqa: E402

RE_KM = 6371.2
FMT = "%Y-%m-%d/%H:%M"


def day_mask(probe, day):
    t0 = day.strftime("%Y-%m-%d")
    t1 = (day + timedelta(days=1)).strftime("%Y-%m-%d")
    var = f"th{probe}_pos_gsm"
    themis.state(probe=probe, trange=[t0, t1], varnames=[var], time_clip=True)
    d = get_data(var)
    if d is None:
        return None, None
    x, y = d.y[:, 0] / RE_KM, d.y[:, 1] / RE_KM
    return d.times, (x < -8) & (np.abs(y) < np.abs(x))


def runs_to_intervals(times, mask, min_h=2):
    # contiguous true runs, floor both ends to the hour, drop short ones
    # floor on the end is what reproduces events_tail2015.csv, ceil runs an hour long
    out = []
    i = 0
    while i < len(mask):
        if not mask[i]:
            i += 1
            continue
        j = i
        while j < len(mask) and mask[j]:
            j += 1
        a = datetime.fromtimestamp(times[i], tz=timezone.utc).replace(
            minute=0, second=0, microsecond=0)
        b = datetime.fromtimestamp(times[j - 1], tz=timezone.utc).replace(
            minute=0, second=0, microsecond=0)
        if (b - a) >= timedelta(hours=min_h):
            out.append((a, b))
        i = j
    return out


def merge(ivs):
    # join runs on the same day when the gap is under an hour
    out = []
    for a, b in ivs:
        if out and a - out[-1][1] < timedelta(hours=1):
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def main():
    p = argparse.ArgumentParser(description="tail-region interval scan")
    p.add_argument("--season", default="tail2015", help="interval_tag and output name")
    p.add_argument("--start", required=True, help="first day, yyyy-mm-dd")
    p.add_argument("--end", required=True, help="stop day, exclusive, yyyy-mm-dd")
    p.add_argument("--probes", nargs="+", default=["a", "d", "e"])
    p.add_argument("--min-hours", type=int, default=2)
    p.add_argument("--out", default=None)
    args = p.parse_args()
    args.out = args.out or str(ROOT / f"events_{args.season}.csv")

    start = datetime.strptime(args.start, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    end = datetime.strptime(args.end, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    rows = []
    for probe in args.probes:
        day = start
        while day < end:
            times, mask = day_mask(probe, day)
            if times is None:
                print(f"[warn] no state th{probe} {day:%Y-%m-%d}")
            else:
                last = day.replace(hour=23, minute=59)
                for a, b in merge(runs_to_intervals(times, mask, args.min_hours)):
                    b = min(b, last)
                    if b <= a:
                        continue
                    rows.append({"probe": probe, "trange_start": a.strftime(FMT),
                                 "trange_end": b.strftime(FMT),
                                 "interval_tag": args.season})
            day += timedelta(days=1)
    rows.sort(key=lambda r: (r["trange_start"], r["probe"]))

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["probe", "trange_start", "trange_end", "interval_tag"])
        w.writeheader()
        w.writerows(rows)

    for probe in args.probes:
        sub = [r for r in rows if r["probe"] == probe]
        hrs = sum((datetime.strptime(r["trange_end"], FMT)
                   - datetime.strptime(r["trange_start"], FMT)).total_seconds() / 3600
                  for r in sub)
        print(f"th{probe}: {len(sub)} intervals, {hrs:.1f} h")
    print(f"[ok] wrote {args.out} ({len(rows)} intervals)")


if __name__ == "__main__":
    main()
