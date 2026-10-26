# joins s/c potential onto episodes as a sidecar, per spin, so e_b gets corrected before the episode median
# ions lose q*phi climbing the positive sheath, so e_true = e_meas + phi (singly charged)
# gates and model features stay on measured e_b, detection happens in the instrument frame
# reads analysis/episodes.jsonl + member snapshots, writes analysis/episode_scpot.jsonl keyed by episode_id
import argparse
import glob
import json
from pathlib import Path

import numpy as np
from cdflib import CDF

from build_episode_catalog import EB_FLOOR

p = argparse.ArgumentParser()
p.add_argument("--season", default="tail2015")
# efi-to-plasma calibration in thin lobe plasma is good to several volts, carried as 1 sigma
p.add_argument("--sig", type=float, default=5.0)
a = p.parse_args()
ROOT = Path(__file__).resolve().parents[1]
CAND = ROOT / "candidates" / a.season
ANA = CAND / "analysis"
# max gap to the nearest potential sample, peir is spin cadence, peif ~100 s
MAX_DT = {"peir": 10.0, "peif": 300.0}
# peir steps 0.238 dex up to the ~300 eV bin, 0.174 to ~450, then 0.12, a peak pinned to one bin carries
# w/sqrt(12) rounding, w = mean gap to both neighbours, the episode median doesnt shrink it for a steady beam
# bins ~300 (295-304 by probe), ~450, ~580, w 0.207, 0.147, 0.12, edges are geometric midpoints between bins
QUANT = ((290.0, 0.17), (366.0, 0.15), (514.0, 0.10), (np.inf, 0.08))
SHEATH = 0.8   # phi/e_b above this the sheath bends trajectories, direction gate unreliable
FIT_MAX = 0.5  # phi/e_b above this, keep in catalog but out of quantitative fits
# earth shadow cylinder from ephemeris, gse and gsm share x so sqrt(y2+z2) is frame free
# 1.1 re covers penumbra and atmosphere at tail distances, efi potential and spin phase are bad in shadow
SHADOW_R = 1.1


def load_pot(probe, day):
    fs = sorted(glob.glob(str(ROOT / "data" / f"th{probe}" / "l2" / "esa" / day[:4]
                              / f"th{probe}_l2_esa_{day}_v*.cdf")))
    if not fs:
        return None
    c = CDF(fs[-1])
    out = {}
    for dt in MAX_DT:
        try:
            t, v = np.asarray(c.varget(f"th{probe}_{dt}_time")), np.asarray(c.varget(f"th{probe}_{dt}_sc_pot"))
        except Exception:
            continue
        ok = np.isfinite(t) & np.isfinite(v)
        o = np.argsort(t[ok])
        out[dt] = (t[ok][o], v[ok][o])
    return out


def phi_at(pot, t):
    for dt, (tt, vv) in pot.items():
        if not len(tt):
            continue
        i = np.clip(np.searchsorted(tt, t), 1, len(tt) - 1)
        if min(abs(tt[i] - t), abs(tt[i - 1] - t)) <= MAX_DT[dt]:
            return float(np.interp(t, tt, vv))
    return None


def quant(eb):
    return next(q for lim, q in QUANT if eb < lim)


def med(v):
    v = [x for x in v if x is not None]
    return round(float(np.median(v)), 2) if v else None


def main():
    eps = [json.loads(line) for line in open(ANA / "episodes.jsonl")]
    pos = {r["episode_id"]: r for r in map(json.loads, open(ANA / "episode_positions.jsonl"))}
    need = {m for e in eps for m in e["member_ids"]}
    recs = {}
    for f in {e["file"] for e in eps}:
        for line in open(CAND / f):
            r = json.loads(line)
            if r["candidate_id"] in need:
                recs[r["candidate_id"]] = (r["t_center"], r["features"].get("E_b"), r["context"].get("te_ev"))

    cache, out, nopot = {}, [], 0
    for e in eps:
        key = (e["probe"], e["t_start_ut"][:10].replace("-", ""))
        if key not in cache:
            cache[key] = load_pot(*key)
            if cache[key] is None:
                print(f"[warn] no esa l2 cdf for th{key[0]} {key[1]}")
        pot = cache[key]
        phis, eb_c, ebte_m, ebte_c = [], [], [], []
        for m in e["member_ids"]:
            t, eb, te = recs.get(m, (None, None, None))
            ph = phi_at(pot, t) if pot and t is not None else None
            phis.append(None if ph is None else round(ph, 2))
            # same truthy E_b filter as group_episodes so medians line up with e_b_med
            if not eb:
                continue
            if ph is not None:
                eb_c.append(eb + ph)
            if te:
                ebte_m.append(eb / te)
                if ph is not None:
                    ebte_c.append((eb + ph) / te)
        ph_ok = [x for x in phis if x is not None]
        if not ph_ok:
            nopot += 1
        pm = med(ph_ok)
        frac = round(pm / e["e_b_med"], 3) if pm is not None and e.get("e_b_med") else None
        ps = pos.get(e["episode_id"])
        r_yz = float(np.hypot(ps["y_gsm"], ps["z_gsm"])) if ps else None
        out.append({
            "episode_id": e["episode_id"],
            "rep_id": e["rep_id"],
            "scpot_med": pm,
            "scpot_min": round(min(ph_ok), 2) if ph_ok else None,
            "scpot_max": round(max(ph_ok), 2) if ph_ok else None,
            "scpot_frac": frac,
            "e_b_corr_med": med(eb_c),
            "e_b_corr_min": round(min(eb_c), 2) if eb_c else None,
            "e_b_corr_max": round(max(eb_c), 2) if eb_c else None,
            "e_b_corr_err": round(float(np.hypot(quant(e["e_b_med"]) * e["e_b_med"], a.sig)), 2) if eb_c else None,
            # true-frame floor the tier-1 gate implies for this beam, truncation bound for fits
            "eb_floor_true": round(EB_FLOOR + pm, 2) if pm is not None else None,
            "eb_over_te_meas": med(ebte_m),
            "eb_over_te_corr": med(ebte_c),
            "sheath_flag": frac is not None and frac > SHEATH,
            "fit_ok": frac is not None and frac <= FIT_MAX,
            "r_yz": round(r_yz, 3) if r_yz is not None else None,
            "eclipse_flag": r_yz is not None and r_yz < SHADOW_R and pos[e["episode_id"]]["x_gsm"] < 0,
            "member_scpot": phis,
        })

    with open(ANA / "episode_scpot.jsonl", "w") as fh:
        for r in out:
            fh.write(json.dumps(r) + "\n")
    print(f"scpot joined: {len(out) - nopot}/{len(out)} episodes, {len(cache)} probe-days, sig {a.sig} V")
    print(f"eclipse: {sum(r['eclipse_flag'] for r in out)} in shadow, {sum(r['r_yz'] is None for r in out)} without position")
    print(f"[ok] wrote {ANA / 'episode_scpot.jsonl'}")


if __name__ == "__main__":
    main()
