# train logistic regression on the labeled candidate dataset
# event-level splits, benchmarked against the heuristic classifier

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (precision_recall_curve, auc, precision_score, recall_score, f1_score)
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

CANDIDATES = Path(__file__).resolve().parents[1] / "candidates" 
OUT = CANDIDATES / "analysis"

# model A: detection features only
# labeling leaned on ne_ni/te/beta/dir_x_bx so those are model B only
# flux_z hinged at 5, moderate excursion reads beam, extreme reads
# population change, a single linear coef cannot express both
# anchor_dt/anchor_score/n_near are label-free episode context, built from
# the candidate pool itself (heuristic is_beam/score, never labels)
FEATS_A0 = ["R", "sig_margin", "peak_prom", "abs_asym", "eb_over_de", "log_E_b",
            "log_dE", "flux_z", "flux_z_excess", "flux_z_perp", "duration",
            "chain_e_slope", "chain_e_scatter", "width", "para_to_omni",
            "asym_dev", "hyst_promoted", "anchor_fwd", "anchor_bwd",
            "anchor_score", "n_near"]
# log_te promoted 2026-07-31, user decision, provenance caveat documented:
# te is a labeling feature, validation of te-dependent calls rests on the
# anchors, a0 stays te-free as the morphology audit trail
FEATS_A = FEATS_A0 + ["log_te"]
FEATS_B = FEATS_A + ["ne_ni_ratio", "log_beta", "dir_x_bx"]

NBR_WIN = 600.0   # episode window, s
NBR_DEX = 0.3     # energy match, factor 2


def add_neighbor_feats(recs):
    # per probe+day, self excluded, energy-matched neighbors only
    by_day = {}
    for r in recs:
        day = (r["probe"], r["candidate_id"].split("_")[1].split("T")[0])
        by_day.setdefault(day, []).append(r)
    for group in by_day.values():
        for r in group:
            e, t = r["features"].get("E_b"), r["t_center"]
            if not e or e <= 0:
                r["_nbr"] = {}
                continue
            dt_f, dt_b, sc, n = None, None, None, 0
            for o in group:
                if o is r:
                    continue
                oe = o["features"].get("E_b")
                if not oe or oe <= 0 or abs(np.log10(oe / e)) > NBR_DEX:
                    continue
                d = o["t_center"] - t
                if o["is_beam"]:
                    # split by side, a leading edge has anchors ahead only
                    if d >= 0 and (dt_f is None or d < dt_f):
                        dt_f = d
                    if d <= 0 and (dt_b is None or -d < dt_b):
                        dt_b = -d
                if abs(d) <= NBR_WIN:
                    n += 1
                    s = o["features"].get("beam_score")
                    if s is not None and (sc is None or s > sc):
                        sc = s
            def _log(v):
                return np.log10(1 + v) if v is not None else None
            dt_a = min(x for x in (dt_f, dt_b, np.inf) if x is not None)
            r["_nbr"] = {
                "anchor_dt": _log(dt_a) if np.isfinite(dt_a) else None,
                "anchor_fwd": _log(dt_f),
                "anchor_bwd": _log(dt_b),
                "anchor_score": sc,
                "n_near": float(n),
            }


def load_all():
    # survey files first so their records win the dedupe,
    # neighbor feats need the full pool so they are computed here
    files = sorted(CANDIDATES.glob("*.jsonl"),
                   key=lambda p: (0 if "survey" in p.name else 1, p.name))
    seen = set()
    recs = []
    for jl in files:
        for line in open(jl):
            r = json.loads(line)
            key = (r["probe"], round(r["t_center"], 1))
            if key in seen:
                continue
            seen.add(key)
            recs.append(r)
    add_neighbor_feats(recs)
    return recs


def load_labeled():
    recs, events = [], []
    for r in load_all():
        if r.get("label") not in ("positive", "negative"):
            continue
        # probe + date, hour tags merge, works with or without the T part
        event = "_".join(r["candidate_id"].split("_")[:2]).split("T")[0]
        recs.append(r)
        events.append(event)
    return recs, np.array(events)

def log10(v):
    return np.log10(v) if v is not None and v > 0 else None

def raw_features(r):
    f, c = r["features"], r["context"]
    nbr = r.get("_nbr", {})
    a = f.get("asymmetry")
    z = f.get("flux_z")
    return {
        "R": f.get("R"),
        "sig_margin": f.get("sig_margin"),
        "peak_prom": f.get("peak_prom"),
        "abs_asym": abs(a) if a is not None else None,
        "eb_over_de": f.get("Eb_over_dE"),
        "log_E_b": log10(f.get("E_b")),
        "log_dE": log10(f.get("dE")),
        "flux_z": min(z, 5.0) if z is not None else None,
        "flux_z_excess": max(z - 5.0, 0.0) if z is not None else None,
        "flux_z_perp": f.get("flux_z_perp"),
        "duration": f.get("duration_steps"),
        "chain_e_slope": f.get("chain_e_slope"),
        "chain_e_scatter": f.get("chain_e_scatter"),
        "width": f.get("width"),
        "para_to_omni": f.get("para_to_omni"),
        "asym_dev": f.get("asym_dev"),
        "hyst_promoted": float(f.get("hyst_promoted") or 0),
        "anchor_dt": nbr.get("anchor_dt"),
        "anchor_fwd": nbr.get("anchor_fwd"),
        "anchor_bwd": nbr.get("anchor_bwd"),
        "anchor_score": nbr.get("anchor_score"),
        "n_near": nbr.get("n_near"),
        "ne_ni_ratio": c.get("ne_ni_ratio"),
        "log_te": log10(c.get("te_ev")),
        "log_beta": log10(c.get("beta")),
        "dir_x_bx": c.get("dir_x_bx"),
    }

def matrix(recs, names):
    # median impute + was-missing indicator
    raw = np.array([[raw_features(r)[n] if raw_features(r)[n] is not None
                     else np.nan for n in names] for r in recs], dtype=float)
    bad = [n for j, n in enumerate(names) if np.all(np.isnan(raw[:, j]))]
    if bad:
        raise ValueError(f"all-nan feature columns: {bad}")
    med = np.nanmedian(raw, axis=0)
    miss = np.isnan(raw)
    filled = np.where(miss, med, raw)
    X = np.hstack([filled, miss.astype(float)])
    cols = names + [f"{n}_missing" for n in names]
    return X, cols

def out_of_fold_probs(X, y, events, n_splits=5):
    probs = np.full(len(y), np.nan)
    gkf = GroupKFold(n_splits=n_splits)
    for tr, te in gkf.split(X, y, groups=events):
        model = make_pipeline(StandardScaler(),
                              LogisticRegression(class_weight="balanced",
                                                 max_iter=1000))
        model.fit(X[tr], y[tr])
        probs[te] = model.predict_proba(X[te])[:, 1]
    return probs

def scores(y, pred):
    return (precision_score(y, pred, zero_division=0),
            recall_score(y, pred), f1_score(y, pred))

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    recs, events = load_labeled()
    y = np.array([1 if r["label"] == "positive" else 0 for r in recs])
    print(f"{len(recs)} labeled, {y.sum()} pos / {(1 - y).sum()} neg, "
          f"{len(set(events))} events")

    lines = ["# beam classifier model report", "",
             f"n = {len(recs)} labeled candidates "
             f"({y.sum()} positive, {(1 - y).sum()} negative), "
             f"{len(set(events))} events, GroupKFold(5) by event", ""]

    fig, ax = plt.subplots(figsize=(7, 5))
    results = {}
    for tag, names in (("A0_morphology_audit", FEATS_A0),
                       ("A_deployed", FEATS_A), ("B_all_features", FEATS_B)):
        X, cols = matrix(recs, names)
        probs = out_of_fold_probs(X, y, events)
        pred = (probs >= 0.5).astype(int)
        p, r_, f1 = scores(y, pred)
        prec, rec, thr = precision_recall_curve(y, probs)
        pr_auc = auc(rec, prec)
        # best-f1 operating point off the oof curve, 0.5 is arbitrary at this n
        f1s = 2 * prec * rec / np.clip(prec + rec, 1e-9, None)
        j = int(np.nanargmax(f1s[:-1]))
        results[tag] = (probs, pred)
        lines += [f"## model {tag}",
                  f"out-of-fold: precision {p:.3f}, recall {r_:.3f}, "
                  f"f1 {f1:.3f}, pr-auc {pr_auc:.3f}",
                  f"best-f1 threshold {thr[j]:.2f}: precision {prec[j]:.3f}, "
                  f"recall {rec[j]:.3f}, f1 {f1s[j]:.3f}", ""]
        ax.plot(rec, prec, label=f"{tag} (auc {pr_auc:.2f})")

        # coefficients from a full fit, standardized to rank importance
        model = make_pipeline(StandardScaler(),
                              LogisticRegression(class_weight="balanced",
                                                 max_iter=1000))
        model.fit(X, y)
        coefs = model.named_steps["logisticregression"].coef_[0]
        order = np.argsort(-np.abs(coefs))
        lines += ["| feature | coef |", "|---|---|"]
        lines += [f"| {cols[i]} | {coefs[i]:+.2f} |" for i in order
                  if abs(coefs[i]) > 0.05]
        lines.append("")

    # circularity audit, chain feats were label evidence for these records,
    # rerun model a without them and confirm the result holds
    chain_ev = np.array(["chain_evidence" in (r.get("label_source") or "")
                         for r in recs])
    if chain_ev.any():
        sub = ~chain_ev
        Xs, _ = matrix([r for r, k in zip(recs, sub) if k], FEATS_A)
        ps = out_of_fold_probs(Xs, y[sub], events[sub])
        p, r_, f1 = scores(y[sub], (ps >= 0.5).astype(int))
        lines += ["## circularity audit",
                  f"model A minus {chain_ev.sum()} chain_evidence-labeled "
                  f"record(s): precision {p:.3f}, recall {r_:.3f}, "
                  f"f1 {f1:.3f}", ""]

    # heuristic baselines on same labeled set
    for name, pred in (
            ("heuristic is_beam", np.array([r["is_beam"] for r in recs], int)),
            ("all strict gates", np.array([all(r["gates"].values())
                                           for r in recs], int))):
        p, r_, f1 = scores(y, pred)
        lines += [f"## baseline: {name}",
                  f"precision {p:.3f}, recall {r_:.3f}, f1 {f1:.3f}", ""]
        ax.plot(r_, p, "o", label=name)

    # every error is 1 cutout away from physical explanation
    probs, pred = results["A_deployed"]
    lines += ["## model A errors (out-of-fold, threshold 0.5)", "",
              "| candidate | label | prob |", "|---|---|---|"]
    for i in np.where(pred != y)[0]:
        lines.append(f"| {recs[i]['candidate_id']} | {recs[i]['label']} "
                     f"| {probs[i]:.2f} |")
    lines += ["", "## caveats",
              "- small n, weak-supervision labels, single-catalog provenance",
              "- model B includes the features the labels were built from, "
              "treat it as a ceiling not a result"]

    ax.set_xlabel("recall")
    ax.set_ylabel("precision")
    ax.legend(fontsize=8)
    ax.set_title("out-of-fold precision-recall, event-level splits")
    fig.tight_layout()
    fig.savefig(OUT / "model_pr_curve.png", dpi=110)

    (OUT / "model_report.md").write_text("\n".join(lines))
    print(f"[ok] wrote {OUT / 'model_report.md'} and model_pr_curve.png")

if __name__ == "__main__":
    main()