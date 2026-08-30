# links episodes into physical events, single-linkage union-find over a
# time + flow-direction predicate with dY as a reject-only filter
# cross-probe TOF is NOT used, it is sub-cadence in peif survey, see docs/RESULTS_HIGH.md
# needs scripts/add_episode_positions.py to have run first
import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ANA = ROOT / "candidates" / "tail2015" / "analysis"

# bridges the ~400 s slow-survey gap, below that the 12-31 11:59 structure splits in two
# anchors are identical from 480 through 900 s so this sits at the start of the stable plateau
TAU = 480.0
# cross-tail separation beyond which probes sit on different field lines
DY_MAX = 2.0
# anchors that any acceptable linking must reproduce, see docs/RESULTS_HIGH.md
ANCHORS = {
    "12-31 08:33 triple": ["d_2015-12-31_083258", "e_2015-12-31_083302",
                           "a_2015-12-31_083402", "e_2015-12-31_083437",
                           "e_2015-12-31_083709"],
    "12-20 11:21 structure": ["d_2015-12-20_112131", "e_2015-12-20_112345",
                              "a_2015-12-20_112847"],
    "12-31 11:59 structure": ["d_2015-12-31_115454", "e_2015-12-31_120218",
                              "a_2015-12-31_120457"],
}


class DSU:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def load(catalog_only):
    eps = {}
    for line in open(ANA / "episodes.jsonl"):
        e = json.loads(line)
        eps[e["episode_id"]] = e
    pos = {}
    for line in open(ANA / "episode_positions.jsonl"):
        p = json.loads(line)
        pos[p["episode_id"]] = p
    cat = {}
    if catalog_only:
        with open(ANA / "beam_catalog_2015.csv") as fh:
            cat = {r["episode_id"]: r for r in csv.DictReader(fh)}
        ids = [i for i in cat if i in eps and i in pos]
    else:
        ids = [i for i in eps if i in pos]

    items = []
    for i in ids:
        e, p = eps[i], pos[i]
        items.append({
            "episode_id": i,
            "probe": e["probe"],
            "t0": datetime.fromisoformat(e["t_start_ut"]).timestamp(),
            "t1": datetime.fromisoformat(e["t_end_ut"]).timestamp(),
            "dir": e["dir_x_bx"],
            "e_b": e.get("e_b_med"),
            "n_samples": e["n_samples"],
            "x": p["x_gsm"], "y": p["y_gsm"], "z": p["z_gsm"],
            "tier": int(cat[i]["tier"]) if i in cat else 0,
            "prob": float(cat[i]["prob"]) if i in cat else None,
            "label": cat[i]["label"] if i in cat else "",
        })
    items.sort(key=lambda r: r["t0"])
    return items


def _gap(a, b):
    # 0 when the spans overlap, else the dead time between them
    return max(0.0, max(a["t0"], b["t0"]) - min(a["t1"], b["t1"]))


def link(items, tau, dy_max):
    d = DSU(len(items))
    for i, a in enumerate(items):
        for j in range(i + 1, len(items)):
            b = items[j]
            # sorted by t0 so once b starts too late nothing after it can match
            if b["t0"] - a["t1"] > tau:
                break
            if a["dir"] != b["dir"] or _gap(a, b) > tau:
                continue
            if a["probe"] != b["probe"] and abs(a["y"] - b["y"]) > dy_max:
                continue
            d.union(i, j)
    return d


def groups_of(items, dsu):
    g = {}
    for idx, it in enumerate(items):
        g.setdefault(dsu.find(idx), []).append(it)
    return list(g.values())


def summarize(members):
    members = sorted(members, key=lambda m: m["t0"])
    probes = sorted({m["probe"] for m in members})
    xs = [m["x"] for m in members]
    ys = [m["y"] for m in members]
    zs = [m["z"] for m in members]
    ebs = [m["e_b"] for m in members if m["e_b"] is not None]
    t0, t1 = members[0]["t0"], max(m["t1"] for m in members)
    ut = datetime.fromtimestamp(t0, tz=timezone.utc)
    return {
        "event_id": f"{''.join(probes)}_{ut.strftime('%Y-%m-%d_%H%M%S')}",
        "probes": "".join(probes),
        "n_probes": len(probes),
        "n_episodes": len(members),
        "t_start_ut": ut.isoformat(),
        "duration_s": round(t1 - t0, 1),
        "dir_x_bx": members[0]["dir"],
        "n_tier1": sum(1 for m in members if m["tier"] == 1),
        "n_tier2": sum(1 for m in members if m["tier"] == 2),
        "prob_max": max([m["prob"] for m in members if m["prob"] is not None],
                        default=None),
        "e_b_min": round(min(ebs), 1) if ebs else None,
        "e_b_max": round(max(ebs), 1) if ebs else None,
        "dx_span": round(max(xs) - min(xs), 3),
        "dy_span": round(max(ys) - min(ys), 3),
        "dz_span": round(max(zs) - min(zs), 3),
        "n_pos": sum(1 for m in members if m["label"] == "positive"),
        "n_neg": sum(1 for m in members if m["label"] == "negative"),
        "member_ids": [m["episode_id"] for m in members],
    }


def null_test(items, tau, dy_max, n, seed, span_h=6.0):
    # shift each probe-day block by a random offset, kills cross-probe alignment
    # while preserving within-probe clustering
    # caveat: positions travel with the episode so the dY filter is approximate here
    rng = np.random.default_rng(seed)
    blocks = {}
    for it in items:
        day = datetime.fromtimestamp(it["t0"], tz=timezone.utc).strftime("%Y-%m-%d")
        blocks.setdefault((it["probe"], day), []).append(it)
    out = []
    for _ in range(n):
        shifted = []
        for grp in blocks.values():
            off = float(rng.uniform(-span_h * 3600, span_h * 3600))
            for it in grp:
                s = dict(it)
                s["t0"] += off
                s["t1"] += off
                shifted.append(s)
        shifted.sort(key=lambda r: r["t0"])
        gs = groups_of(shifted, link(shifted, tau, dy_max))
        out.append(sum(1 for g in gs if len({m["probe"] for m in g}) > 1))
    return np.array(out)


def main():
    p = argparse.ArgumentParser(description="link episodes into multi-probe events")
    p.add_argument("--tau", type=float, default=TAU,
                   help="max dead time between linked episodes, s")
    p.add_argument("--dy-max", type=float, default=DY_MAX,
                   help="max cross-tail probe separation, Re")
    p.add_argument("--all", action="store_true",
                   help="link every episode instead of catalog members only")
    p.add_argument("--null", type=int, default=0,
                   help="run N time-scrambled trials for the coincidence floor")
    p.add_argument("--seed", type=int, default=0)
    # not events_2015.csv, that name is taken by the root interval scan
    p.add_argument("--out", default=str(ANA / "linked_events_2015.csv"))
    args = p.parse_args()

    items = load(catalog_only=not args.all)
    evs = [summarize(g) for g in groups_of(items, link(items, args.tau, args.dy_max))]
    evs.sort(key=lambda e: e["t_start_ut"])

    multi = [e for e in evs if e["n_probes"] > 1]
    print(f"{len(items)} episodes -> {len(evs)} events "
          f"(tau={args.tau:.0f}s dy_max={args.dy_max} Re, "
          f"{'all' if args.all else 'catalog'})")
    print(f"  multi-probe: {len(multi)}  "
          f"(3-probe: {sum(1 for e in multi if e['n_probes'] == 3)})")
    print(f"  singletons:  {sum(1 for e in evs if e['n_episodes'] == 1)}")

    keys = [k for k in evs[0] if k != "member_ids"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(evs)
    with open(Path(args.out).with_suffix(".jsonl"), "w") as fh:
        for e in evs:
            fh.write(json.dumps(e) + "\n")

    print("\nanchor recovery:")
    for name, want in ANCHORS.items():
        hit = [e for e in evs if any(w in e["member_ids"] for w in want)]
        if not hit:
            print(f"  {name:24} NOT FOUND")
            continue
        for e in hit:
            got = sum(1 for w in want if w in e["member_ids"])
            print(f"  {name:24} {e['event_id']:26} {got}/{len(want)} anchors, "
                  f"{e['n_episodes']} eps, {e['probes']}, {e['duration_s']:.0f}s")

    print("\nlargest multi-probe events:")
    for e in sorted(multi, key=lambda e: -e["n_episodes"])[:8]:
        print(f"  {e['event_id']:26} {e['probes']:3} n={e['n_episodes']:>3} "
              f"dur={e['duration_s']:>7.0f}s E_b {e['e_b_min']}-{e['e_b_max']} "
              f"dY={e['dy_span']:.2f} t1={e['n_tier1']}")

    if args.null:
        c = null_test(items, args.tau, args.dy_max, args.null, args.seed)
        obs = len(multi)
        # z is only meaningful when the null actually varies
        z = (obs - c.mean()) / c.std() if c.std() > 0 else float("inf")
        print(f"\nnull test ({args.null} time-scrambled trials):")
        print(f"  observed multi-probe events: {obs}")
        print(f"  scrambled: mean {c.mean():.1f} +/- {c.std():.1f} "
              f"(min {c.min()}, max {c.max()})")
        print(f"  excess {obs - c.mean():+.1f}, z = {z:.1f}, "
              f"trials >= observed: {(c >= obs).sum()}/{args.null}")

    print(f"\n[ok] wrote {args.out}")


if __name__ == "__main__":
    main()
