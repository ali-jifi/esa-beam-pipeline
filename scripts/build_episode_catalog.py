# episode features, lobe-scoped model, and two-tier episode catalog
# constants below are the frozen 2026-08-11 closeout values, see LABELING.md
# rerunning REFITS on current labels and overwrites the catalog, use --out to test
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

TAIL = Path(__file__).resolve().parents[1] / "candidates" / "tail2015"
ANA = TAIL / "analysis"

# hard gates, upstream of the model, never features (deck c doctrine)
DIR_GATE = -1          # tailward only
EB_MAX = 600.0         # outflow energy window, ev
BETA_MAX = 0.03        # lobe scope
TE_MAX = 500.0         # closes the corridor leak, te breakers max at 276
# catalog cut
PROB_THR = 0.82        # blind-validated, contamination 1/14 at this cut
PROM_BAR = 0.55        # positive prominence minimum

FEATS = ["prom_max", "flux_z_max", "r_max", "sig_max", "p2o_med",
         "log_eb_med", "eb_span_dex", "fzp_med", "asym_dev_med", "log_n"]


def build_features():
    eps = [json.loads(l) for l in open(ANA / "episodes.jsonl")]
    recs = {}
    for f in TAIL.glob("*.jsonl"):
        for line in open(f):
            r = json.loads(line)
            recs[r["candidate_id"]] = r
    rows = []
    for e in eps:
        ms = [recs[m] for m in e["member_ids"] if m in recs]
        def mx(key):
            vals = [m["features"].get(key) for m in ms]
            vals = [v for v in vals if v is not None and np.isfinite(v)]
            return max(vals) if vals else None
        def med(key):
            vals = [m["features"].get(key) for m in ms]
            vals = [v for v in vals if v is not None and np.isfinite(v)]
            return float(np.median(vals)) if vals else None
        lab = None
        for m in ms:
            if m.get("label") in ("positive", "negative"):
                lab = m["label"]
                break
        rows.append({
            "episode_id": e["episode_id"],
            "day": e["episode_id"].split("_")[1],
            "label": lab, "meta": e,
            "log_n": np.log10(e["n_samples"]),
            "prom_max": mx("peak_prom"), "flux_z_max": mx("flux_z"),
            "r_max": mx("R"), "sig_max": mx("sig_margin"),
            "p2o_med": med("para_to_omni"),
            "log_eb_med": np.log10(e["e_b_med"]) if e.get("e_b_med") else None,
            "eb_span_dex": (np.log10(e["e_b_max"] / e["e_b_min"])
                            if e.get("e_b_max") and e.get("e_b_min") else None),
            "fzp_med": med("flux_z_perp"),
            "asym_dev_med": med("asym_dev"),
        })
    return rows


def gates_ok(r):
    e = r["meta"]
    te = e.get("te_med")
    return (e.get("dir_x_bx") == DIR_GATE
            and e.get("e_b_med") is not None and e["e_b_med"] < EB_MAX
            and e.get("beta_med") is not None and e["beta_med"] < BETA_MAX
            and (te is None or te < TE_MAX))


def mat(rws):
    raw = np.array([[r[n] if r[n] is not None else np.nan for n in FEATS]
                    for r in rws], dtype=float)
    med_ = np.nanmedian(raw, axis=0)
    miss = np.isnan(raw)
    return np.hstack([np.where(miss, med_, raw), miss.astype(float)])


def main():
    p = argparse.ArgumentParser(description="build the episode catalog")
    p.add_argument("--out", default=str(ANA / "beam_catalog_2015.csv"))
    args = p.parse_args()

    rows = build_features()
    lobe = [r for r in rows if gates_ok(r)]
    lab = [r for r in lobe if r["label"]]
    y = np.array([1 if r["label"] == "positive" else 0 for r in lab])
    days = np.array([r["day"] for r in lab])
    X = mat(lab)

    # day-blocked oof report, splits by DAY never probe or episode
    from sklearn.metrics import roc_auc_score
    probs = np.full(len(y), np.nan)
    gkf = GroupKFold(n_splits=min(5, len(set(days))))
    for tr, te in gkf.split(X, y, groups=days):
        m = make_pipeline(StandardScaler(),
                          LogisticRegression(class_weight="balanced", max_iter=1000))
        m.fit(X[tr], y[tr])
        probs[te] = m.predict_proba(X[te])[:, 1]
    print(f"lobe-scoped oof auc {roc_auc_score(y, probs):.3f} "
          f"({len(y)} labeled, {y.sum()} pos, {len(set(days))} day-folds)")

    full = make_pipeline(StandardScaler(),
                         LogisticRegression(class_weight="balanced", max_iter=1000))
    full.fit(X, y)
    pall = full.predict_proba(mat(lobe))[:, 1]

    cat = []
    for r, pr in zip(lobe, pall):
        mfire, pfire = pr >= PROB_THR, (r["prom_max"] or 0) >= PROM_BAR
        tier = 1 if (mfire and pfire) else (2 if (mfire or pfire) else 0)
        if tier:
            e = r["meta"]
            cat.append({"tier": tier, "episode_id": r["episode_id"],
                        "t_start_ut": e["t_start_ut"], "duration_s": e["duration_s"],
                        "n_samples": e["n_samples"], "prob": round(float(pr), 3),
                        "prom_max": r["prom_max"], "e_b_med": e.get("e_b_med"),
                        "te_med": e.get("te_med"), "ne_ni_med": e.get("ne_ni_med"),
                        "beta_med": e.get("beta_med"), "label": r["label"] or ""})
    cat.sort(key=lambda c: (c["tier"], c["t_start_ut"]))
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cat[0].keys()))
        w.writeheader()
        w.writerows(cat)
    for t in (1, 2):
        sub = [c for c in cat if c["tier"] == t]
        lp = sum(1 for c in sub if c["label"] == "positive")
        ln = sum(1 for c in sub if c["label"] == "negative")
        print(f"tier {t}: {len(sub)} ({lp} lpos, {ln} lneg, {len(sub)-lp-ln} unlabeled)")
    print(f"[ok] wrote {args.out}")


if __name__ == "__main__":
    main()
