# redraw cutout pngs from existing candidate jsonls, records never touched
import argparse
import csv
import json
from pathlib import Path

from esa_plotting.config import set_data_dir
from esa_plotting.beam_pipeline import (
    load_esd_distribution, load_bfield_dsl, compute_pa_spectra,
)
from generate_candidates import render_cutout, CANDIDATES


def load_probs(catalog):
    # candidate_id -> model prob, episode-keyed catalogs fan out over member_ids
    # via the episodes.jsonl sitting next to the catalog
    rows = list(csv.DictReader(open(catalog)))
    if not rows or "candidate_id" in rows[0]:
        return {r["candidate_id"]: float(r["prob"]) for r in rows}
    ep = {r["episode_id"]: float(r["prob"]) for r in rows}
    out = {}
    for line in open(Path(catalog).parent / "episodes.jsonl"):
        e = json.loads(line)
        if e["episode_id"] in ep:
            out.update({m: ep[e["episode_id"]] for m in e["member_ids"]})
    return out


def main():
    p = argparse.ArgumentParser(description="re-render candidate cutouts")
    p.add_argument("--out", default=str(CANDIDATES))
    p.add_argument("--cut-dir", default=None, help="default <out>/cutouts")
    p.add_argument("--datatype", default="peif", choices=["peif", "peir"])
    p.add_argument("--existing", action="store_true",
                   help="only redraw cutouts already on disk")
    p.add_argument("--catalog", default=None,
                   help="csv with prob column, shown in the cutout title")
    args = p.parse_args()

    data_dir = set_data_dir()
    cand_dir = Path(args.out)
    cut_dir = Path(args.cut_dir) if args.cut_dir else cand_dir / "cutouts"
    cut_dir.mkdir(exist_ok=True)
    probs = load_probs(args.catalog) if args.catalog else None
    have = {f.stem for f in cut_dir.glob("*.png")} if args.existing else None

    for jl in sorted(cand_dir.glob("*.jsonl")):
        recs = [json.loads(l) for l in open(jl)]
        if have is not None:
            recs = [r for r in recs if r["candidate_id"] in have]
        if not recs:
            continue
        probe, trange = recs[0]["probe"], recs[0]["trange"]
        dist = load_esd_distribution(probe, trange, data_dir, datatype=args.datatype)
        b_times, b_dsl = load_bfield_dsl(probe, trange, data_dir)
        spectra = compute_pa_spectra(dist, b_times, b_dsl)
        for rec in recs:
            cid = rec["candidate_id"]
            prob = probs.get(cid, float("nan")) if probs is not None else None
            render_cutout(spectra, rec, cut_dir / f"{cid}.png", prob=prob)
        print(f"[ok] {jl.name}: {len(recs)} cutouts re-rendered")


if __name__ == "__main__":
    main()
