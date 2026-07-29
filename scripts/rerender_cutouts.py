# redraw cutout pngs from existing candidate jsonls, records never touched
import argparse
import json
from pathlib import Path

from esa_plotting.config import set_data_dir
from esa_plotting.beam_pipeline import (
    load_esd_distribution, load_bfield_dsl, compute_pa_spectra,
)
from generate_candidates import render_cutout, CANDIDATES


def main():
    p = argparse.ArgumentParser(description="re-render candidate cutouts")
    p.add_argument("--out", default=str(CANDIDATES))
    args = p.parse_args()

    data_dir = set_data_dir()
    cand_dir = Path(args.out)
    cut_dir = cand_dir / "cutouts"
    cut_dir.mkdir(exist_ok=True)

    for jl in sorted(cand_dir.glob("*.jsonl")):
        recs = [json.loads(l) for l in open(jl)]
        if not recs:
            continue
        probe, trange = recs[0]["probe"], recs[0]["trange"]
        dist = load_esd_distribution(probe, trange, data_dir)
        b_times, b_dsl = load_bfield_dsl(probe, trange, data_dir)
        spectra = compute_pa_spectra(dist, b_times, b_dsl)
        for rec in recs:
            render_cutout(spectra, rec, cut_dir / f"{rec['candidate_id']}.png")
        print(f"[ok] {jl.name}: {len(recs)} cutouts re-rendered")


if __name__ == "__main__":
    main()
