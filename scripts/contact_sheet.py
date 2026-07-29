# tile candidate cutouts into contact sheets for fast eyeball labeling
import argparse
import csv
from pathlib import Path

from PIL import Image

CANDIDATES = Path(__file__).resolve().parents[1] / "candidates"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--csv", default=str(CANDIDATES / "analysis" / "auto_positive_candidates.csv"))
    p.add_argument("--out", default=str(CANDIDATES / "analysis"))
    p.add_argument("--cols", type=int, default=4)
    p.add_argument("--per-page", type=int, default=16)
    p.add_argument("--scale", type=float, default=0.5)
    args = p.parse_args()

    with open(args.csv, newline="") as fh:
        ids = [row["candidate_id"] for row in csv.DictReader(fh)]
    cut_dir = CANDIDATES / "cutouts"
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = Path(args.csv).stem

    pages = [ids[i:i + args.per_page] for i in range(0, len(ids), args.per_page)]
    for pi, page in enumerate(pages, 1):
        tiles = []
        for cid in page:
            fp = cut_dir / f"{cid}.png"
            if not fp.exists():
                print(f"[warn] missing cutout {cid}")
                continue
            im = Image.open(fp)
            tiles.append(im.resize((int(im.width * args.scale),
                                    int(im.height * args.scale))))
        if not tiles:
            continue
        w, h = tiles[0].width, tiles[0].height
        rows = (len(tiles) + args.cols - 1) // args.cols
        sheet = Image.new("RGB", (args.cols * w, rows * h), "white")
        for k, tile in enumerate(tiles):
            sheet.paste(tile, ((k % args.cols) * w, (k // args.cols) * h))
        out = out_dir / f"{stem}_sheet{pi}.png"
        sheet.save(out)
        print(f"[ok] {out.name}: {len(tiles)} cutouts")


if __name__ == "__main__":
    main()
