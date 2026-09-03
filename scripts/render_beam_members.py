# per-member cutouts for every catalog beam, one folder per episode under
# analysis/beams/tier{n}/<episode_id>/, rep cutouts in analysis/cutouts/ refreshed on the way
import argparse
import csv
import json
import shutil
from pathlib import Path

from esa_plotting.config import set_data_dir
from esa_plotting.beam_pipeline import (load_esd_distribution, load_bfield_dsl,
                                        compute_pa_spectra)
from generate_candidates import render_cutout
from rerender_cutouts import load_probs

TAIL = Path(__file__).resolve().parents[1] / "candidates" / "tail2015"
ANA = TAIL / "analysis"


def main():
    p = argparse.ArgumentParser(description="per-member cutouts per catalog beam")
    p.add_argument("--catalog", default=str(ANA / "beam_catalog_2015.csv"))
    p.add_argument("--out", default=str(ANA / "beams"))
    p.add_argument("--datatype", default="peir", choices=["peif", "peir"])
    p.add_argument("--file", default=None, help="only jsonls whose name contains this")
    args = p.parse_args()

    data_dir = set_data_dir()
    cut = ANA / "cutouts"
    reps = {f.stem for f in cut.glob("*.png")}
    probs = load_probs(args.catalog)
    tier = {r["episode_id"]: r["tier"] for r in csv.DictReader(open(args.catalog))}
    dest = {}
    for line in open(ANA / "episodes.jsonl"):
        e = json.loads(line)
        if e["episode_id"] in tier:
            d = Path(args.out) / f"tier{tier[e['episode_id']]}" / e["episode_id"]
            dest.update({m: d for m in e["member_ids"]})
    want = set(dest) | reps

    for jl in sorted(TAIL.glob("*.jsonl")):
        if args.file and args.file not in jl.name:
            continue
        # id sits first on every line, cheap filter before the json parse
        recs = [json.loads(l) for l in open(jl) if l.split('"')[3] in want]
        if not recs:
            continue
        a = recs[0]
        dist = load_esd_distribution(a["probe"], a["trange"], data_dir,
                                     datatype=args.datatype)
        bt, bd = load_bfield_dsl(a["probe"], a["trange"], data_dir)
        spectra = compute_pa_spectra(dist, bt, bd)
        for r in recs:
            cid = r["candidate_id"]
            prob = probs.get(cid, float("nan"))
            if cid in dest:
                dest[cid].mkdir(parents=True, exist_ok=True)
                out = dest[cid] / f"{cid}.png"
                render_cutout(spectra, r, out, prob=prob)
                if cid in reps:
                    shutil.copy(out, cut / f"{cid}.png")
            else:
                render_cutout(spectra, r, cut / f"{cid}.png", prob=prob)
        print(f"[ok] {jl.name}: {len(recs)} cutouts", flush=True)


if __name__ == "__main__":
    main()
