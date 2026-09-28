# episode features, lobe-scoped model, and two-tier episode catalog
# constants below are the frozen 2026-08-11 closeout values, see LABELING.md
# --fit refits on current labels and saves the model, default applies the frozen one, use --out to test
import argparse
import csv
import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

TAIL = Path(__file__).resolve().parents[1] / "candidates" / "tail2015"
ANA = TAIL / "analysis"
MODEL = ANA / "episode_model.joblib"

# hard gates, upstream of the model, never features (deck c doctrine)
DIR_GATE = -1          # tailward only
EB_MAX = 600.0         # outflow energy window, ev
BETA_MAX = 0.03        # lobe scope
TE_MAX = 500.0         # closes the corridor leak, te breakers max at 276
# catalog cut
PROB_THR = 0.82        # blind-validated, contamination 1/14 at this cut
PROM_BAR = 0.55        # positive prominence minimum
ISO_WIN = 600.0        # s, window for counting nearby beam samples
# tier 1 definition, not tuned, a single spin or a peak in the lowest bin (~34 ev) cant be a measured beam
MIN_SAMPLES = 2
EB_FLOOR = 40.0        # ev

# d2o and asym_dev_dir are dominant-cone versions of para_to_omni and asym_dev,
# the 2015 fit saw parallel beams almost only and both flip or collapse in the
# other lobe where tailward beams are anti-parallel
FEATS = ["prom_max", "flux_z_max", "r_max", "sig_max", "d2o_med",
         "log_eb_med", "eb_span_dex", "fzp_med", "asym_dev_dir", 
         "log_n", "log_te", "log_beta", "log_iso", "duty", "e_stab", "eb_drift"]

def build_features():
    eps = [json.loads(l) for l in open(ANA / "episodes.jsonl")]
    recs = {}
    for f in TAIL.glob("*.jsonl"):
        for line in open(f):
            r = json.loads(line)
            recs[r["candidate_id"]] = r
    # snapshot times per interval file, for isolation and duty cycle
    t_all, t_beam = {}, {}
    for r in recs.values():
        k = r["candidate_id"].rsplit("_", 2)[0]
        t_all.setdefault(k, []).append(r["t_center"])
        if r["is_beam"]:
            t_beam.setdefault(k, []).append(r["t_center"])
    t_all = {k: np.sort(v) for k, v in t_all.items()}
    t_beam = {k: np.sort(v) for k, v in t_beam.items()}
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
        d2o, adev = [], []
        for m in ms:
            f, d = m["features"], m.get("direction") or 0
            p2o, asym, dev = f.get("para_to_omni"), f.get("asymmetry"), f.get("asym_dev")
            # anti/omni = para/omni * anti/para, and anti/para = (1-asym)/(1+asym)
            if d == -1 and p2o is not None and asym is not None and asym > -1:
                p2o = p2o * (1 - asym) / (1 + asym)
            d2o.append(p2o)
            adev.append(dev * d if dev is not None and d else None)
        def medv(vals):
            vals = [v for v in vals if v is not None and np.isfinite(v)]
            return float(np.median(vals)) if vals else None
        k = e["file"].split("_survey")[0]
        mt = np.array([m["t_center"] for m in ms])
        t0, t1 = mt.min(), mt.max()
        tb, ta = t_beam.get(k, np.array([])), t_all[k]
        near = np.searchsorted(tb, t1 + ISO_WIN, "right") - np.searchsorted(tb, t0 - ISO_WIN)
        span = np.searchsorted(ta, t1, "right") - np.searchsorted(ta, t0)
        eb = np.array([m["features"].get("E_b") or np.nan for m in ms], float)
        ok = np.isfinite(eb) & (eb > 0)
        le = np.log10(eb[ok])
        e_stab = float(np.mean(np.abs(le - np.median(le)) <= 0.12)) if ok.sum() else None
        # dex per minute, needs 3 points to mean anything
        drift = (float(np.polyfit((mt[ok] - t0) / 60.0, le, 1)[0])
                 if ok.sum() >= 3 and np.ptp(mt[ok]) > 0 else None)

        rows.append({
            "episode_id": e["episode_id"],
            "day": e["episode_id"].split("_")[1],
            "label": lab, "meta": e,
            "log_n": np.log10(e["n_samples"]),
            "prom_max": mx("peak_prom"), "flux_z_max": mx("flux_z"),
            "r_max": mx("R"), "sig_max": mx("sig_margin"),
            "d2o_med": medv(d2o),
            "log_eb_med": np.log10(e["e_b_med"]) if e.get("e_b_med") else None,
            "eb_span_dex": (np.log10(e["e_b_max"] / e["e_b_min"])
                            if e.get("e_b_max") and e.get("e_b_min") else None),
            "fzp_med": med("flux_z_perp"),
            "asym_dev_dir": medv(adev),
            "log_te": np.log10(e["te_med"]) if e.get("te_med") else None,
            "log_beta": np.log10(max(e["beta_med"], 1e-5)) if e.get("beta_med") is not None else None,
            "log_iso": np.log10(1 + max(near - len(ms), 0)), "duty": len(ms) / max(span, 1),
            "e_stab": e_stab, "eb_drift": drift,

        })
    return rows

def gates_ok(r):
    e = r["meta"]
    te = e.get("te_med")
    return (e.get("dir_x_bx") == DIR_GATE
            and e.get("e_b_med") is not None and e["e_b_med"] < EB_MAX
            and e.get("beta_med") is not None and e["beta_med"] < BETA_MAX
            and (te is None or te < TE_MAX))

def mat(rws, med_=None):
    raw = np.array([[r[n] if r[n] is not None else np.nan for n in FEATS]
                    for r in rws], dtype=float)
    if med_ is None:
        med_ = np.nanmedian(raw, axis=0)
    miss = np.isnan(raw)
    return np.hstack([np.where(miss, med_, raw), miss.astype(float)]), med_

def deck_labels():
    # user labels from blind decks with a key, plus the 2011 and 2012 tuning decks
    files = [f for f in ANA.glob("blind_deck_*.csv")
             if (ANA / f"blind_key_{f.stem[len('blind_deck_'):]}_SEALED.csv").exists()]
    files += list(ANA.glob("tier1_deck_*_labeled.csv"))
    labs = {}
    for f in files:
        for r in csv.DictReader(open(f)):
            if r.get("label") in ("positive", "negative"):
                labs[r["episode_id"]] = r["label"]
    return labs

def main():
    p = argparse.ArgumentParser(description="build the episode catalog")
    p.add_argument("--out", default=None)
    p.add_argument("--season", default="tail2015")
    p.add_argument("--fit", action="store_true",
                   help="refit on current labels and save the model")
    p.add_argument("--model", default=str(MODEL),
                   help="frozen model path, default is the 2015 one")
    p.add_argument("--fit-seasons", nargs="+",
                   default=["tail2011", "tail2012", "tail2013", "tail2014", "tail2015", "tail2016"])

    args = p.parse_args()
    global TAIL, ANA
    TAIL = Path(__file__).resolve().parents[1] / "candidates" / args.season
    ANA = TAIL / "analysis"
    args.out = args.out or str(ANA / f"beam_catalog_{args.season[4:]}.csv")

    rows = build_features()
    lobe = [r for r in rows if gates_ok(r)]
    if args.fit:
        root = TAIL.parent
        pool, lab, grp = [], [], []
        for s in args.fit_seasons:
            TAIL, ANA = root / s, root / s / "analysis"
            dl = deck_labels()
            for r in build_features():
                if not gates_ok(r):
                    continue
                pool.append(r)
                r["label"] = dl.get(r["episode_id"]) or r["label"]
                if r["label"]:
                    lab.append(r)
                    grp.append(s)
        TAIL, ANA = root / args.season, root / args.season / "analysis"
        _, lobe_med = mat(pool)
        X, _ = mat(lab, med_=lobe_med)
        y = np.array([1 if r["label"] == "positive" else 0 for r in lab])
        g = np.array(grp)
        prom = np.array([r["prom_max"] or 0 for r in lab])

        # leave one season out, the honest number for a new season
        from sklearn.metrics import roc_auc_score
        probs = np.full(len(y), np.nan)
        for s in args.fit_seasons:
            te = g == s
            m = make_pipeline(StandardScaler(),
                            LogisticRegression(class_weight="balanced", max_iter=1000))
            m.fit(X[~te], y[~te])
            probs[te] = m.predict_proba(X[te])[:, 1]
            auc = roc_auc_score(y[te], probs[te]) if len(set(y[te])) == 2 else float("nan")
            print(f"held out {s}: n={te.sum()} pos={y[te].sum()} auc {auc:.3f}")
        print(f"pooled loso auc {roc_auc_score(y, probs):.3f}")
        for thr in (0.7, 0.75, 0.82, 0.9):
            t1 = (probs >= thr) & (prom >= PROM_BAR)
            print(f"thr {thr}: tier1 neg {(t1 & (y == 0)).sum()}/{t1.sum()}, pos in tier1 {(t1 & (y == 1)).sum()}/{y.sum()}")

        full = make_pipeline(StandardScaler(),
                            LogisticRegression(class_weight="balanced", max_iter=1000))
        full.fit(X, y)
        joblib.dump({"model": full, "medians": lobe_med, "feats": FEATS}, args.model)
        print(f"[ok] saved model {args.model}")
    else:
        saved = joblib.load(args.model)
        full, lobe_med = saved["model"], saved["medians"]
    Xl, _ = mat(lobe, med_=lobe_med)
    pall = full.predict_proba(Xl)[:, 1]

    cat = []
    for r, pr in zip(lobe, pall):
        mfire, pfire = pr >= PROB_THR, (r["prom_max"] or 0) >= PROM_BAR
        tier = 1 if (mfire and pfire) else (2 if (mfire or pfire) else 0)
        if tier == 1 and (r["meta"]["n_samples"] < MIN_SAMPLES or (r["meta"].get("e_b_med") or 0) <= EB_FLOOR):
            tier = 2
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
