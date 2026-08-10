# group is_beam samples into episodes, gap<=60s, same direction,
# adjacent members within 0.3 dex in energy
# episode_ids here are referenced by group= tags in labeled records
import json
import glob
from pathlib import Path

import numpy as np

TAIL = Path(__file__).resolve().parents[1] / "candidates" / "tail2015"
OUT = TAIL / "analysis"
OUT.mkdir(exist_ok=True)
GAP = 60.0
DEX = 0.3

episodes = []
for f in sorted(TAIL.glob("*.jsonl")):
    recs = [json.loads(l) for l in open(f)]
    beams = sorted((r for r in recs if r["is_beam"]),
                   key=lambda r: r["t_center"])
    cur = []
    def flush():
        if not cur:
            return
        eb = [r["features"].get("E_b") for r in cur if r["features"].get("E_b")]
        te = [r["context"].get("te_ev") for r in cur
              if r["context"].get("te_ev") and np.isfinite(r["context"]["te_ev"])]
        nn = [r["context"].get("ne_ni_ratio") for r in cur
              if r["context"].get("ne_ni_ratio") and np.isfinite(r["context"]["ne_ni_ratio"])]
        fz = [r["features"].get("flux_z") for r in cur
              if r["features"].get("flux_z") is not None]
        pr = [r["features"].get("peak_prom") for r in cur
              if r["features"].get("peak_prom") is not None]
        bt = [r["context"].get("beta") for r in cur
              if r["context"].get("beta") and np.isfinite(r["context"]["beta"])]
        rep = max(cur, key=lambda r: r["features"].get("beam_score") or 0)
        episodes.append({
            "episode_id": f"{cur[0]['probe']}_{cur[0]['t_center_ut'][:10]}_{cur[0]['t_center_ut'][11:19].replace(':','')}",
            "probe": cur[0]["probe"],
            "file": f.name,
            "t_start_ut": cur[0]["t_center_ut"],
            "t_end_ut": cur[-1]["t_center_ut"],
            "duration_s": round(cur[-1]["t_center"] - cur[0]["t_center"], 1),
            "n_samples": len(cur),
            "direction": cur[0]["direction"],
            "dir_x_bx": cur[0]["context"].get("dir_x_bx"),
            "e_b_med": round(float(np.median(eb)), 1) if eb else None,
            "e_b_min": round(min(eb), 1) if eb else None,
            "e_b_max": round(max(eb), 1) if eb else None,
            "te_med": round(float(np.median(te)), 1) if te else None,
            "ne_ni_med": round(float(np.median(nn)), 2) if nn else None,
            "beta_med": float(np.median(bt)) if bt else None,
            "flux_z_max": round(max(fz), 1) if fz else None,
            "prom_max": round(max(pr), 2) if pr else None,
            "rep_id": rep["candidate_id"],
            "member_ids": [r["candidate_id"] for r in cur],
        })
    for r in beams:
        if not cur:
            cur = [r]
            continue
        prev = cur[-1]
        eb_a, eb_b = prev["features"].get("E_b"), r["features"].get("E_b")
        e_ok = (eb_a and eb_b and abs(np.log10(eb_b / eb_a)) <= DEX)
        if (r["t_center"] - prev["t_center"] <= GAP
                and r["direction"] == prev["direction"] and e_ok):
            cur.append(r)
        else:
            flush()
            cur = [r]
    flush()

with open(OUT / "episodes.jsonl", "w") as fh:
    for e in episodes:
        fh.write(json.dumps(e) + "\n")

# cold tailward selection for the validation sheets
cold = [e for e in episodes
        if e["dir_x_bx"] == -1 and e["te_med"] is not None and e["te_med"] < 200
        and e["e_b_med"] is not None and e["e_b_med"] < 600]
cold.sort(key=lambda e: -e["n_samples"])
import csv
keys = [k for k in episodes[0] if k != "member_ids"]
with open(OUT / "episodes.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
    w.writeheader()
    w.writerows(episodes)
with open(OUT / "episodes_cold_tailward.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
    w.writeheader()
    w.writerows(cold)

n1 = sum(1 for e in episodes if e["n_samples"] == 1)
print(f"episodes: {len(episodes)} total ({n1} single-sample)")
print(f"cold tailward (<600 eV, Te<200, dir_x_bx=-1): {len(cold)}")
multi = [e for e in cold if e["n_samples"] >= 3]
print(f"  of which >=3 samples: {len(multi)}")
print("top 10 cold tailward by samples:")
for e in cold[:10]:
    print(f"  {e['episode_id']}: n={e['n_samples']} dur={e['duration_s']}s "
          f"E_b {e['e_b_min']}-{e['e_b_max']} Te={e['te_med']} ne/ni={e['ne_ni_med']} fz_max={e['flux_z_max']}")
