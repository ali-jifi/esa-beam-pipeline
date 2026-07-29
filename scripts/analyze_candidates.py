# pool candidate jsonls and plot feature distributions for auto-label calibration
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CANDIDATES = Path(__file__).resolve().parents[1] / "candidates"

CONFIRMED_POS = ["d_2015-12-20T1115_survey", ]  # interval stem of the anchor beams
ANCHOR_BEAM_UTS = ("11:48:30", "11:53:33", "11:56:55", "12:08:41")


def load(cand_dir):
    # dedupe on physical timestep, overlapping intervals emit the same candidate
    # under different stems, survey files first so their records win
    files = sorted(cand_dir.glob("*.jsonl"),
                   key=lambda p: (0 if "survey" in p.name else 1, p.name))
    seen, recs = set(), []
    for jl in files:
        for line in open(jl):
            r = json.loads(line)
            key = (r["probe"], round(r["t_center"], 1))
            if key in seen:
                continue
            seen.add(key)
            recs.append(r)
    return recs


def col(recs, group, key):
    out = []
    for r in recs:
        v = r[group].get(key) if group else r.get(key)
        out.append(np.nan if v is None else float(v))
    return np.array(out)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default=str(CANDIDATES / "analysis"))
    args = p.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    recs = load(CANDIDATES)
    strict = np.array([all(r["gates"].values()) for r in recs])
    print(f"{len(recs)} unique candidates, {strict.sum()} pass all strict gates")

    panels = [
        ("R", col(recs, "features", "R"), "log"),
        ("E_b [eV]", col(recs, "features", "E_b"), "log"),
        ("sig_margin", col(recs, "features", "sig_margin"), "linear"),
        ("peak_prom", col(recs, "features", "peak_prom"), "linear"),
        ("|asymmetry|", np.abs(col(recs, "features", "asymmetry")), "linear"),
        ("Eb_over_dE", col(recs, "features", "Eb_over_dE"), "linear"),
        ("ne_ni_ratio", col(recs, "context", "ne_ni_ratio"), "linear"),
        ("Te [eV]", col(recs, "context", "te_ev"), "log"),
        ("beta", col(recs, "context", "beta"), "log"),
    ]

    fig, axes = plt.subplots(3, 3, figsize=(14, 10))
    for ax, (name, v, scale) in zip(axes.flat, panels):
        fin = np.isfinite(v)
        if scale == "log":
            fin &= v > 0
            bins = np.logspace(np.log10(np.nanpercentile(v[fin], 0.5)),
                               np.log10(np.nanpercentile(v[fin], 99.5)), 40)
            ax.set_xscale("log")
        else:
            bins = np.linspace(np.nanpercentile(v[fin], 0.5),
                               np.nanpercentile(v[fin], 99.5), 40)
        ax.hist(v[fin & ~strict], bins=bins, color="0.7", label="gate fail")
        ax.hist(v[fin & strict], bins=bins, color="crimson", alpha=0.6,
                label="all strict gates")
        ax.set_xlabel(name)
        ax.set_yscale("log")
    axes.flat[0].legend(fontsize=8)
    fig.suptitle(f"candidate feature distributions, n={len(recs)}")
    fig.tight_layout()
    png = out_dir / "feature_distributions.png"
    fig.savefig(png, dpi=110)
    plt.close(fig)
    print(f"[ok] wrote {png}")

    # dir_x_bx composition among strict-pass
    dxb = col(recs, "context", "dir_x_bx")
    for grp, m in (("strict-pass", strict), ("gate-fail", ~strict)):
        fin = np.isfinite(dxb) & m
        n_tail = int(np.sum(dxb[fin] == -1))
        n_earth = int(np.sum(dxb[fin] == 1))
        print(f"{grp}: tailward {n_tail}, earthward {n_earth}, "
              f"unknown {int(m.sum() - fin.sum())}")

    # split quantiles for the discriminators
    print(f"\n{'feature':>12} {'split':>12} {'p10':>8} {'p50':>8} {'p90':>8}")
    for name, v, _ in panels:
        for grp, m in (("strict-pass", strict), ("gate-fail", ~strict)):
            fin = np.isfinite(v) & m
            if fin.sum() < 5:
                continue
            q = np.nanpercentile(v[fin], [10, 50, 90])
            print(f"{name:>12} {grp:>12} {q[0]:>8.2f} {q[1]:>8.2f} {q[2]:>8.2f}")

    # candidate auto-positive rule, thresholds from the confirmed-beam signature
    ne_ni = col(recs, "context", "ne_ni_ratio")
    te = col(recs, "context", "te_ev")
    tag = np.array([r["interval_tag"] for r in recs])
    rule = (strict & (dxb == -1) & (ne_ni >= 2.0) & (te <= 200)
            & (tag == "beam_period"))
    print(f"\nproposed auto-positive rule: strict AND tailward AND "
          f"ne_ni>=2 AND Te<=200 AND beam_period -> {rule.sum()} candidates")

    # anchors must land on the right side
    for r, is_pos in zip(recs, rule):
        ut = r["t_center_ut"][11:19]
        if r["candidate_id"].startswith("d_2015-12-20T1115") and ut in ANCHOR_BEAM_UTS:
            print(f"  confirmed beam {ut}: rule={'POS' if is_pos else 'NEG'}")
        if r.get("label") == "negative":
            print(f"  labeled negative {r['candidate_id']}: "
                  f"rule={'POS (BAD)' if is_pos else 'NEG (correct)'}")

    # silver tier, one notch looser, minus gold, for hand labeling
    silver = (strict & (dxb == -1) & (ne_ni >= 1.5) & (te <= 500)
              & (tag == "beam_period") & ~rule)
    print(f"silver tier (ne_ni>=1.5, Te<=500, minus gold): {silver.sum()} candidates")

    def write_csv(name, mask):
        rows = [r for r, m in zip(recs, mask) if m]
        out_csv = out_dir / name
        with open(out_csv, "w") as fh:
            fh.write("candidate_id,ut,probe,E_b,R,sig_margin,ne_ni,te_ev,beta\n")
            for r in rows:
                f, c = r["features"], r["context"]
                fh.write(f"{r['candidate_id']},{r['t_center_ut']},{r['probe']},"
                         f"{f['E_b']},{f['R']},{f['sig_margin']},"
                         f"{c['ne_ni_ratio']},{c['te_ev']},{c['beta']}\n")
        print(f"[ok] wrote {out_csv} ({len(rows)} rows)")

    write_csv("auto_positive_candidates.csv", rule)
    write_csv("silver_candidates.csv", silver)


if __name__ == "__main__":
    main()
