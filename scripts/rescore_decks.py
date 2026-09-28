# rescore the blind-deck labels against the current catalogs, labels stay as they are, keys unused
import csv
from collections import defaultdict
from datetime import datetime
from pathlib import Path

from score_blind_decks import ci, TAU

CAND = Path(__file__).resolve().parents[1] / "candidates"
# deck, season, who labeled, role for the current model
DECKS = [("2011b", "tail2011", "user", "training"), ("2012b", "tail2012", "user", "training"),
         ("2013", "tail2013", "user", "training"), ("2014census", "tail2014", "user", "training"),
         ("2015held", "tail2015", "user", "training"), ("2016", "tail2016", "user", "training"),
         ("2017", "tail2017", "assistant", "design"), ("2019v2", "tail2019", "assistant", "design"),
         ("2023", "tail2023", "assistant", "design"), ("2022", "tail2022", "assistant", "fresh"),
         ("2024v2", "tail2024", "assistant", "fresh")]
QUIET = {"tail2011", "tail2012", "tail2013", "tail2014", "tail2015", "tail2016", "tail2017", "tail2019"}


def load():
    rows = []
    for name, season, who, role in DECKS:
        ana = CAND / season / "analysis"
        cat = {r["episode_id"]: r for r in csv.DictReader(open(ana / f"beam_catalog_{season[4:]}.csv"))}
        eps = {r["episode_id"]: r for r in csv.DictReader(open(ana / "episodes.csv"))}
        for r in csv.DictReader(open(ana / f"blind_deck_{name}.csv")):
            lab = r["label"].strip().lower()
            if lab not in ("positive", "negative"):
                continue
            e, c = eps[r["episode_id"]], cat.get(r["episode_id"])
            rows.append({"deck": name, "season": season, "who": who, "role": role, "label": lab,
                         "tier": int(c["tier"]) if c else 0, "probe": e["probe"],
                         "t0": datetime.fromisoformat(e["t_start_ut"]).timestamp(),
                         "t1": datetime.fromisoformat(e["t_end_ut"]).timestamp()})
    return rows


def structures(t1):
    by = defaultdict(list)
    for r in t1:
        by[(r["deck"], r["probe"])].append(r)
    out = []
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


def line(name, rows):
    t1 = [r for r in rows if r["tier"] == 1]
    neg = sum(r["label"] == "negative" for r in t1)
    lo, hi = ci(neg, len(t1))
    st = structures(t1) if t1 else []
    sneg = sum(not any(m["label"] == "positive" for m in s) for s in st)
    pos = [r for r in rows if r["label"] == "positive"]
    p1, p12 = sum(r["tier"] == 1 for r in pos), sum(r["tier"] in (1, 2) for r in pos)
    print(f"{name:26} decided {len(rows):4d} | tier1 contamination {neg:3d}/{len(t1):<4d} {100*neg/max(len(t1),1):4.0f}% "
          f"(CI {100*lo:3.0f}-{100*hi:3.0f}) | by structure {sneg}/{len(st)} | pos in tier1 {p1}/{len(pos)}, tier1+2 {p12}/{len(pos)}")


def main():
    rows = load()
    for name, *_ in DECKS:
        line(name, [r for r in rows if r["deck"] == name])
    print()
    line("quiet 2011-2019", [r for r in rows if r["season"] in QUIET])
    line("solar max 2022-2024", [r for r in rows if r["season"] not in QUIET])
    for role in ("training", "design", "fresh"):
        line(f"role: {role}", [r for r in rows if r["role"] == role])
    for who in ("user", "assistant"):
        line(f"labeled by {who}", [r for r in rows if r["who"] == who])
    line("all", rows)


if __name__ == "__main__":
    main()
