# two-tier beam catalog, heuristic generates candidates, model a judges
# tier 1 = both fire (high confidence), tier 2 = exactly one fires (review)
import argparse
import csv
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_recall_curve
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from train_model import CANDIDATES, FEATS_A, load_all, matrix, out_of_fold_probs


def make_model():
    return make_pipeline(StandardScaler(),
                         LogisticRegression(class_weight="balanced",
                                            max_iter=1000))


def main():
    p = argparse.ArgumentParser(description="score candidates into a catalog")
    p.add_argument("--threshold", type=float, default=None,
                   help="model threshold, default = oof best-f1")
    p.add_argument("--out", default=str(CANDIDATES / "analysis" / "beam_catalog.csv"))
    args = p.parse_args()

    recs = load_all()
    X, _ = matrix(recs, FEATS_A)
    lab = np.array([{"positive": 1, "negative": 0}.get(r.get("label"), -1)
                    for r in recs])
    tr = lab >= 0
    print(f"{len(recs)} candidates, {tr.sum()} labeled")

    thr = args.threshold
    if thr is None:
        # self-calibrate off the oof pr curve, 0.5 is arbitrary at this n
        events = np.array(["_".join(r["candidate_id"].split("_")[:2]).split("T")[0]
                           for r in recs])
        oof = out_of_fold_probs(X[tr], lab[tr], events[tr])
        prec, rec, thrs = precision_recall_curve(lab[tr], oof)
        f1s = 2 * prec * rec / np.clip(prec + rec, 1e-9, None)
        j = int(np.nanargmax(f1s[:-1]))
        thr = float(thrs[j])
        print(f"oof best-f1 threshold {thr:.2f} "
              f"(P {prec[j]:.3f} R {rec[j]:.3f} f1 {f1s[j]:.3f})")

    model = make_model()
    model.fit(X[tr], lab[tr])
    probs = model.predict_proba(X)[:, 1]
    heur = np.array([bool(r["is_beam"]) for r in recs])
    mfire = probs >= thr

    def tier(i):
        if mfire[i] and heur[i]:
            return 1
        if mfire[i] or heur[i]:
            return 2
        return 0

    tiers = np.array([tier(i) for i in range(len(recs))])
    # hard earthward post-filter, user decision 2026-07-31: dir_x_bx=+1 is
    # psbl territory whatever the scores say, None (no b data) passes
    dropped = sum(1 for i, r in enumerate(recs)
                  if tiers[i] and r["context"].get("dir_x_bx") == 1)
    for i, r in enumerate(recs):
        if tiers[i] and r["context"].get("dir_x_bx") == 1:
            tiers[i] = 0
    print(f"earthward post-filter dropped {dropped}")
    rows = sorted((i for i in range(len(recs)) if tiers[i]),
                  key=lambda i: (tiers[i], recs[i]["t_center"]))
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tier", "candidate_id", "t_center_ut", "prob", "is_beam",
                    "direction", "E_b", "R", "peak_prom", "flux_z",
                    "duration_steps", "label", "label_confidence"])
        for i in rows:
            r, f = recs[i], recs[i]["features"]
            def g(k):
                v = f.get(k)
                return f"{v:.3g}" if v is not None and np.isfinite(v) else ""
            w.writerow([tiers[i], r["candidate_id"], r["t_center_ut"],
                        f"{probs[i]:.3f}", int(heur[i]), r["direction"],
                        g("E_b"), g("R"), g("peak_prom"), g("flux_z"),
                        g("duration_steps"), r.get("label") or "",
                        r.get("label_confidence") or ""])

    for t in (1, 2):
        m = tiers == t
        lp = int((lab[m] == 1).sum())
        ln = int((lab[m] == 0).sum())
        print(f"tier {t}: {m.sum()} candidates "
              f"({lp} labeled pos, {ln} labeled neg, {m.sum() - lp - ln} unlabeled)")
    print(f"[ok] wrote {args.out}")


if __name__ == "__main__":
    main()
