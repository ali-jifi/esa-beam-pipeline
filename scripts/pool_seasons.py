# pool per-season catalogs, episodes and coverage into candidates/all
# tier 1 per 100 angle-resolved probe-hours is the cross-season number, never the raw
# count, scanned hours overstate coverage where fast survey only covers part of an interval
import csv
import glob
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CAND = ROOT / "candidates"
OUT = CAND / "all"
FMT = "%Y-%m-%d/%H:%M"
COV_COLS = ["probe", "trange_start", "trange_end", "hours", "processed", "hours_angle"]


def read(path):
    return list(csv.DictReader(open(path, newline="")))


def write(path, rows, cols=None):
    if not rows:
        print(f"[skip] {path.name}, nothing to write")
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols or list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"[ok] wrote {path} ({len(rows)} rows)")


def angle_hours(probe, t0, t1):
    # hours of >=50 bin reduced-mode samples inside the window, from the peir cdf on disk
    # this is what the loader keeps, the 1-angle omni packing is dropped
    from cdflib import CDF
    from cdflib.epochs import CDFepoch
    fs = sorted(glob.glob(str(ROOT / "data" / f"th{probe}" / "l2" / "esd" / t0.strftime("%Y")
                              / f"th{probe}_l2_esa_peir_{t0:%Y%m%d}_*.cdf")))
    if not fs:
        return None
    c = CDF(fs[-1])
    t = np.asarray(CDFepoch.unixtime(c.varget("epoch")))
    nb = np.asarray(c.varget("nbins"))[c.varget("an_ind")]
    u0, u1 = (x.replace(tzinfo=timezone.utc).timestamp() for x in (t0, t1))
    m = (t >= u0) & (t <= u1)
    if m.sum() < 2:
        return 0.0
    dt = float(np.median(np.diff(t[m])))
    return float((nb[m] >= 50).sum() * dt / 3600)


def season_coverage(season):
    # per-interval cache in analysis/, an interval's angle hours never change once its cdf is down
    ev = ROOT / f"events_{season}.csv"
    if not ev.exists():
        return []
    ana = CAND / season / "analysis"
    ana.mkdir(parents=True, exist_ok=True)
    cache_p = ana / "coverage_intervals.csv"
    cache = {(r["probe"], r["trange_start"]): r for r in read(cache_p)} if cache_p.exists() else {}
    done = {p.name[:17] for p in (CAND / season).glob("*_survey.jsonl")}
    rows, new = [], 0
    for r in read(ev):
        a, b = datetime.strptime(r["trange_start"], FMT), datetime.strptime(r["trange_end"], FMT)
        key = (r["probe"], r["trange_start"])
        stem = f"{r['probe']}_{r['trange_start'].replace('/', 'T').replace(':', '')}"
        processed = stem in done
        c = cache.get(key)
        if c is None or c["hours_angle"] == "" or c["processed"] != str(processed):
            try:
                ha = angle_hours(r["probe"], a, b)
            except Exception as e:
                print(f"[warn] {season} {stem}: {type(e).__name__}: {e}")
                ha = None
            c = {"probe": r["probe"], "trange_start": r["trange_start"], "trange_end": r["trange_end"],
                 "hours": round((b - a).total_seconds() / 3600, 2), "processed": str(processed),
                 "hours_angle": "" if ha is None else round(ha, 3)}
            new += 1
        rows.append(c)
    if new:
        write(cache_p, rows, COV_COLS)
    return rows


def main():
    OUT.mkdir(exist_ok=True)
    cats, eps, cataloged = [], [], set()
    for ana in sorted(CAND.glob("tail*/analysis")):
        season = ana.parent.name
        cat = ana / f"beam_catalog_{season[4:]}.csv"
        if cat.exists():
            cataloged.add(season)
            cats += [{"season": season, **r} for r in read(cat)]
        else:
            print(f"[warn] {season}: no {cat.name}")
        ep = ana / "episodes.csv"
        if ep.exists():
            eps += [{"season": season, **r} for r in read(ep)]
    write(OUT / "beam_catalog_all.csv", cats)
    write(OUT / "episodes_all.csv", eps)

    cov = defaultdict(lambda: {"n_intervals": 0, "hours": 0.0, "n_processed": 0,
                               "hours_processed": 0.0, "hours_angle": 0.0, "n_unknown": 0})
    for f in sorted(ROOT.glob("events_tail*.csv")):
        season = f.stem[len("events_"):]
        for r in season_coverage(season):
            c = cov[(season, r["probe"])]
            h = float(r["hours"])
            c["n_intervals"] += 1
            c["hours"] += h
            if r["processed"] == "True":
                c["n_processed"] += 1
                c["hours_processed"] += h
                if r["hours_angle"] != "":
                    c["hours_angle"] += float(r["hours_angle"])
                else:
                    c["n_unknown"] += 1
    write(OUT / "coverage.csv",
          [{"season": s, "probe": pr, **{k: (round(v, 2) if isinstance(v, float) else v) for k, v in c.items()}}
           for (s, pr), c in sorted(cov.items())])

    tot = defaultdict(lambda: defaultdict(float))
    for (s, _), c in cov.items():
        for k in ("hours", "hours_processed", "hours_angle", "n_unknown"):
            tot[s][k] += c[k]
    t1 = defaultdict(int)
    for r in cats:
        if r["tier"] == "1":
            t1[r["season"]] += 1
    print(f"\n{'season':10} {'tier1':>6} {'scanned':>8} {'processed':>10} {'angle_h':>8} {'t1/100h':>8}")
    for s in sorted(set(tot) | cataloged):
        h = tot[s]["hours_angle"]
        flag = "*" if tot[s]["n_unknown"] else " "
        rate = f"{100 * t1[s] / h:8.2f}" if h and s in cataloged else f"{'-':>8}"
        print(f"{s:10} {t1[s]:6d} {tot[s]['hours']:8.1f} {tot[s]['hours_processed']:10.1f} {h:7.1f}{flag} {rate}")
    print("* = some processed intervals had no cdf on disk, angle hours incomplete")


if __name__ == "__main__":
    main()
