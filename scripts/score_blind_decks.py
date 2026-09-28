# score labeled blind decks against their sealed keys
# run only after every label in a deck is final, opening a key ends that decks blindness
import argparse
import csv
import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.stats import beta as beta_dist
from sklearn.metrics import roc_auc_score

CAND = Path(__file__).resolve().parents[1] / "candidates"
THR = 0.82
# same-probe episodes closer than this count as one structure, same tau link_events uses
TAU = 480.0


def ci(k, n, a=0.05):
    # clopper pearson
    if n == 0:
        return (float("nan"), float("nan"))
    lo = 0.0 if k == 0 else beta_dist.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta_dist.ppf(1 - a / 2, k + 1, n - k)
    return lo, hi


def load(deck_csv):
    ana = deck_csv.parent
    name = deck_csv.stem[len("blind_deck_"):]
    key = {r["episode_id"]: r for r in csv.DictReader(open(ana / f"blind_key_{name}_SEALED.csv"))}
    eps = {}
    for line in open(ana / "episodes.jsonl"):
        e = json.loads(line)
        eps[e["episode_id"]] = e
    rows = []
    for r in csv.DictReader(open(deck_csv)):
        e, k = eps[r["episode_id"]], key[r["episode_id"]]
        rows.append({"deck": name, "episode_id": r["episode_id"], "probe": e["probe"],
                     "t0": datetime.fromisoformat(e["t_start_ut"]).timestamp(),
                     "t1": datetime.fromisoformat(e["t_end_ut"]).timestamp(),
                     "label": (r.get("label") or "").strip().lower(),
                     "tier": int(k["tier"]), "prob": float(k["prob"]),
                     "te": float(r["te_med"]) if r["te_med"] not in ("", "None") else None,
                     "beta": float(r["beta_med"])})
    return rows


def structures(rows):
    # same probe and within tau is one structure, positive if any member is
    out = []
    by = defaultdict(list)
    for r in rows:
        by[(r["deck"], r["probe"])].append(r)
    for grp in by.values():
        grp.sort(key=lambda r: r["t0"])
        cur = [grp[0]]
        for r in grp[1:]:
            if r["t0"] - cur[-1]["t1"] <= TAU:
                cur.append(r)
            else:
                out.append(cur)
                cur = [r]
        out.append(cur)
    return out


def report(name, rows, scope=None):
    if scope:
        rows = [r for r in rows if (r["te"] is None or r["te"] < scope[0]) and r["beta"] < scope[1]]
    dec = [r for r in rows if r["label"] in ("positive", "negative")]
    t1 = [r for r in dec if r["prob"] >= THR and r["tier"] == 1]
    neg = sum(1 for r in t1 if r["label"] == "negative")
    lo, hi = ci(neg, len(t1))
    y = [1 if r["label"] == "positive" else 0 for r in dec]
    auc = roc_auc_score(y, [r["prob"] for r in dec]) if len(set(y)) == 2 else float("nan")
    st = [s for s in structures(t1)]
    sneg = sum(1 for s in st if not any(m["label"] == "positive" for m in s))
    slo, shi = ci(sneg, len(st))
    npos = sum(y)
    # only meaningful for census decks, sampled decks fix the tier mix by design
    rec1 = sum(1 for r in dec if r["label"] == "positive" and r["tier"] == 1)
    rec12 = sum(1 for r in dec if r["label"] == "positive" and r["tier"] in (1, 2))
    print(f"{name:14} n={len(rows):3d} decided={len(dec):3d} pos={npos:3d} borderline={sum(1 for r in rows if r['label']=='borderline'):2d}"
          f" | auc {auc:.3f} | tier1 contamination {neg}/{len(t1)} = {100*neg/max(len(t1),1):3.0f}% (CI {100*lo:.0f}-{100*hi:.0f})"
          f" | by structure {sneg}/{len(st)} = {100*sneg/max(len(st),1):3.0f}% (CI {100*slo:.0f}-{100*shi:.0f})"
          f" | pos in tier1 {rec1}/{npos}, tier1+2 {rec12}/{npos}")


def main():
    p = argparse.ArgumentParser(description="score blind decks")
    p.add_argument("--decks", nargs="*", help="deck names like 2024 2012b, default all with a key")
    p.add_argument("--partial", action="store_true", help="allow decks with unlabeled rows")
    p.add_argument("--scope", nargs=2, type=float, metavar=("TE", "BETA"),
                   help="also score under a tighter scope, te_med < TE and beta_med < BETA")
    args = p.parse_args()

    files = sorted(CAND.glob("tail*/analysis/blind_deck_*.csv"))
    files = [f for f in files if (f.parent / f"blind_key_{f.stem[11:]}_SEALED.csv").exists()]
    if args.decks:
        files = [f for f in files if f.stem[11:] in args.decks]
    allrows = []
    for f in files:
        rows = load(f)
        missing = sum(1 for r in rows if not r["label"])
        if missing and not args.partial:
            print(f"[skip] {f.stem[11:]}: {missing} unlabeled rows, key stays sealed")
            continue
        allrows += [r for r in rows if r["label"]]
    if not allrows:
        return
    for name in sorted({r["deck"] for r in allrows}):
        report(name, [r for r in allrows if r["deck"] == name])
    report("pooled", allrows)
    if args.scope:
        print(f"\nunder scope te < {args.scope[0]:g}, beta < {args.scope[1]:g}")
        for name in sorted({r["deck"] for r in allrows}):
            report(name, [r for r in allrows if r["deck"] == name], args.scope)
        report("pooled", allrows, args.scope)
        pos = [r for r in allrows if r["label"] == "positive"]
        kept = [r for r in pos if (r["te"] is None or r["te"] < args.scope[0]) and r["beta"] < args.scope[1]]
        print(f"labeled positives kept by the scope: {len(kept)}/{len(pos)}")


if __name__ == "__main__":
    main()
