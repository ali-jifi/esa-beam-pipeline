# blind validation deck for one season: mixed tiers, spread over days, probs withheld,
# key written to a separate sealed file, cutouts kept apart from the prob-titled ones
import argparse
import csv
import json
import random
from collections import defaultdict
from pathlib import Path

from esa_plotting.config import set_data_dir
from esa_plotting.beam_pipeline import (load_esd_distribution, load_bfield_dsl,
                                        compute_pa_spectra)
from generate_candidates import render_cutout

COLS = ["episode_id", "t_start_ut", "duration_s", "n_samples", "prom_max", "e_b_med",
        "te_med", "ne_ni_med", "beta_med", "rep_id", "label", "label_confidence", "label_note"]


def spread(rows, n, rng):
    # round robin over days so the big days cant fill the deck
    by_day = defaultdict(list)
    for r in rows:
        by_day[r["t_start_ut"][:10]].append(r)
    for v in by_day.values():
        rng.shuffle(v)
    days = list(by_day)
    rng.shuffle(days)
    out = []
    while len(out) < n and any(by_day.values()):
        for d in days:
            if by_day[d] and len(out) < n:
                out.append(by_day[d].pop())
    return out


def main():
    p = argparse.ArgumentParser(description="blind validation deck")
    p.add_argument("--season", required=True)
    p.add_argument("--n1", type=int, default=18, help="tier 1 draws")
    p.add_argument("--n2", type=int, default=12, help="tier 2 draws")
    p.add_argument("--n0", type=int, default=0, help="untiered draws, needs a census catalog")
    p.add_argument("--catalog", default=None, help="catalog csv, default is the season one")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--tag", default="", help="suffix after the year in output names")
    p.add_argument("--exclude", nargs="*", default=[],
                   help="csvs whose episode_ids were already shown to the labeler")
    p.add_argument("--datatype", default="peir", choices=["peif", "peir"])
    args = p.parse_args()

    year = args.season[4:]
    name = f"{year}{args.tag}"
    tail = Path(__file__).resolve().parents[1] / "candidates" / args.season
    ana = tail / "analysis"
    cut = ana / "blind_cutouts"
    cut.mkdir(exist_ok=True)

    seen = set()
    for f in args.exclude:
        seen |= {r["episode_id"] for r in csv.DictReader(open(f))}
    cat = [r for r in csv.DictReader(open(args.catalog or ana / f"beam_catalog_{year}.csv"))
           if r["episode_id"] not in seen]
    eps = {}
    for line in open(ana / "episodes.jsonl"):
        e = json.loads(line)
        eps[e["episode_id"]] = e
    rng = random.Random(args.seed)
    deck = (spread([r for r in cat if r["tier"] == "1"], args.n1, rng)
            + spread([r for r in cat if r["tier"] == "2"], args.n2, rng)
            + spread([r for r in cat if r["tier"] == "0"], args.n0, rng))
    rng.shuffle(deck)

    with open(ana / f"blind_deck_{name}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        for r in deck:
            w.writerow({**{k: r.get(k, "") for k in COLS},
                        "rep_id": eps[r["episode_id"]]["rep_id"],
                        "label": "", "label_confidence": "", "label_note": ""})
    with open(ana / f"blind_key_{name}_SEALED.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["episode_id", "tier", "prob"])
        w.writeheader()
        w.writerows({"episode_id": r["episode_id"], "tier": r["tier"], "prob": r["prob"]}
                    for r in deck)
    n1 = sum(1 for r in deck if r["tier"] == "1")
    print(f"[ok] deck {name}: {len(deck)} episodes ({n1} tier 1) over "
          f"{len({r['t_start_ut'][:10] for r in deck})} days, {len(seen)} excluded")

    data_dir = set_data_dir()
    by_file = defaultdict(set)
    for r in deck:
        e = eps[r["episode_id"]]
        by_file[e["file"]].add(e["rep_id"])
    for fname, want in by_file.items():
        recs = {}
        with open(tail / fname) as fh:
            for line in fh:
                rec = json.loads(line)
                if rec["candidate_id"] in want:
                    recs[rec["candidate_id"]] = rec
        a = next(iter(recs.values()))
        dist = load_esd_distribution(a["probe"], a["trange"], data_dir, datatype=args.datatype)
        bt, bd = load_bfield_dsl(a["probe"], a["trange"], data_dir)
        spectra = compute_pa_spectra(dist, bt, bd)
        for rid, rec in recs.items():
            # prob stays none on purpose, thats the blinding
            render_cutout(spectra, rec, cut / f"{rid}.png", prob=None)
        print(f"[ok] {fname}: {len(recs)} cutouts", flush=True)

    from PIL import Image
    ids = [eps[r["episode_id"]]["rep_id"] for r in deck]
    per = 16
    for pi in range(0, len(ids), per):
        tiles = []
        for cid in ids[pi:pi + per]:
            im = Image.open(cut / f"{cid}.png")
            tiles.append(im.resize((im.width // 2, im.height // 2)))
        w_, h_ = tiles[0].width, tiles[0].height
        sheet = Image.new("RGB", (4 * w_, ((len(tiles) + 3) // 4) * h_), "white")
        for k, t in enumerate(tiles):
            sheet.paste(t, ((k % 4) * w_, (k // 4) * h_))
        out = ana / f"blind{name}_sheet{pi // per + 1}.png"
        sheet.save(out)
        print(f"[ok] {out.name}: {len(tiles)}")


if __name__ == "__main__":
    main()
