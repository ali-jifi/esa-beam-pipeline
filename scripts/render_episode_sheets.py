# render episode rep cutouts and tile contact sheets from a queue csv
# csv needs episode_id or rep_id + file columns, missing cutouts get rendered
import argparse
import csv
import json
from pathlib import Path

from esa_plotting.config import set_data_dir
from esa_plotting.beam_pipeline import (load_esd_distribution, load_bfield_dsl,
                                        compute_pa_spectra)
from generate_candidates import render_cutout

TAIL = Path(__file__).resolve().parents[1] / "candidates" / "tail2015"
ANA = TAIL / "analysis"


def main():
    p = argparse.ArgumentParser(description="episode contact sheets")
    p.add_argument("--csv", required=True)
    p.add_argument("--prefix", required=True, help="sheet filename prefix")
    p.add_argument("--datatype", default="peir", choices=["peif", "peir"])
    p.add_argument("--per-page", type=int, default=16)
    args = p.parse_args()

    data_dir = set_data_dir()
    cut = ANA / "cutouts"
    cut.mkdir(exist_ok=True)
    rows = list(csv.DictReader(open(args.csv)))
    eps = {json.loads(l)["episode_id"]: json.loads(l)
           for l in open(ANA / "episodes.jsonl")}
    jobs = []
    for r in rows:
        rid = r.get("rep_id") or eps[r["episode_id"]]["rep_id"]
        fname = r.get("file") or eps[r["episode_id"]]["file"]
        jobs.append((rid, fname))

    by_file = {}
    for rid, fname in jobs:
        if not (cut / f"{rid}.png").exists():
            by_file.setdefault(fname, []).append(rid)
    for fname, want in by_file.items():
        want = set(want)
        recs = {}
        for line in open(TAIL / fname):
            r = json.loads(line)
            if r["candidate_id"] in want:
                recs[r["candidate_id"]] = r
        if not recs:
            print(f"[warn] no recs in {fname}")
            continue
        a = next(iter(recs.values()))
        dist = load_esd_distribution(a["probe"], a["trange"], data_dir,
                                     datatype=args.datatype)
        bt, bd = load_bfield_dsl(a["probe"], a["trange"], data_dir)
        spectra = compute_pa_spectra(dist, bt, bd)
        for rid, r in recs.items():
            render_cutout(spectra, r, cut / f"{rid}.png")
        print(f"[ok] {fname}: {len(recs)} cutouts")

    from PIL import Image
    ids = [rid for rid, _ in jobs]
    pages = [ids[i:i + args.per_page] for i in range(0, len(ids), args.per_page)]
    for pi, page in enumerate(pages, 1):
        tiles = []
        for cid in page:
            fp = cut / f"{cid}.png"
            if not fp.exists():
                print(f"[warn] missing {cid}")
                continue
            im = Image.open(fp)
            tiles.append(im.resize((im.width // 2, im.height // 2)))
        if not tiles:
            continue
        w, h = tiles[0].width, tiles[0].height
        nrows = (len(tiles) + 3) // 4
        sheet = Image.new("RGB", (4 * w, nrows * h), "white")
        for k, t in enumerate(tiles):
            sheet.paste(t, ((k % 4) * w, (k // 4) * h))
        out = ANA / f"{args.prefix}_sheet{pi}.png"
        sheet.save(out)
        print(f"[ok] {out.name}: {len(tiles)}")


if __name__ == "__main__":
    main()
