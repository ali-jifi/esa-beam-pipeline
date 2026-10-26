# maps of e_b, de, eb over de on y vs x and bx vs r, one marker per episode by probe, plus per probe hists with te
import argparse
import csv
import json
import logging
import os
import warnings
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import (FixedLocator, FuncFormatter, LogLocator, NullFormatter, NullLocator,
                               ScalarFormatter)

warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.WARNING)

ROOT = Path(__file__).resolve().parents[1]
CAND = ROOT / "candidates"
os.environ["THM_DATA_DIR"] = str(ROOT / "data")

from pyspedas import get_data  # noqa: E402
from pyspedas.projects import themis  # noqa: E402

MARKERS = {"a": "^", "d": "o", "e": "s"}
COLORS = {"a": "#2a78d6", "d": "#eb6834", "e": "#1baf7a"}
# peir energy table spacing, e_b medians sit on these steps or between two
ESA_DEX = 0.238
PAD_S = 600.0
# e_b is measured, e_b_corr adds the per-spin s/c potential, eb_over_te uses e_b_corr since l2 te is potential-corrected
COLS = ["season", "episode_id", "probe", "t_start_ut", "n_samples",
        "e_b", "e_b_corr", "scpot", "scpot_min", "scpot_max", "de", "eb_over_de", "te", "eb_over_te",
        "x_gsm", "y_gsm", "z_gsm", "r_xy", "bx_gsm"]


def load_bx(probe, t0, t1):
    a = datetime.fromtimestamp(t0 - PAD_S, tz=timezone.utc)
    b = datetime.fromtimestamp(t1 + PAD_S, tz=timezone.utc)
    trange = [a.strftime("%Y-%m-%d/%H:%M:%S"), b.strftime("%Y-%m-%d/%H:%M:%S")]
    themis.fgm(probe=probe, trange=trange, level="l2", time_clip=True, no_update=True)
    # time_clip leaves a var unclipped when it has no data in range, so check coverage
    for cad in ("fgs", "fgl"):
        d = get_data(f"th{probe}_{cad}_gsm")
        if d is not None and d.times[0] <= t0 and d.times[-1] >= t1:
            return d.times, d.y[:, 0]
    return None


def build(year, tier):
    ana = CAND / f"tail{year}" / "analysis"
    ids = {r["episode_id"] for r in csv.DictReader(open(ana / f"beam_catalog_{year}.csv"))
           if r["tier"] == str(tier)}
    eps = [e for e in map(json.loads, open(ana / "episodes.jsonl")) if e["episode_id"] in ids]
    pos = {p["episode_id"]: p for p in map(json.loads, open(ana / "episode_positions.jsonl"))}
    scp = {p["rep_id"]: p for p in map(json.loads, open(ana / "episode_scpot.jsonl"))}

    recs = {}
    for f in {e["file"] for e in eps}:
        for line in open(CAND / f"tail{year}" / f):
            r = json.loads(line)
            recs[r["candidate_id"]] = r

    groups = {}
    for e in eps:
        groups.setdefault((e["probe"], e["t_start_ut"][:10]), []).append(e)

    rows = []
    for (probe, day), grp in sorted(groups.items()):
        ts = [(datetime.fromisoformat(e["t_start_ut"]).timestamp(),
               datetime.fromisoformat(e["t_end_ut"]).timestamp()) for e in grp]
        try:
            bx = load_bx(probe, min(t[0] for t in ts), max(t[1] for t in ts))
        except Exception as ex:
            bx = None
            print(f"[warn] fgm load failed {probe} {day}: {ex}")
        for e, (t0, t1) in zip(grp, ts):
            p = pos.get(e["episode_id"])
            if p is None:
                print(f"[warn] no position {e['episode_id']}")
                continue
            miss = [m for m in e["member_ids"] if m not in recs]
            if miss:
                print(f"[warn] {len(miss)} members missing from candidates {e['episode_id']}")
                continue
            ms =[recs[m]["features"] for m in e["member_ids"]]
            te = [recs[m]["context"].get("te_ev") for m in e["member_ids"]]
            sp = scp[e["rep_id"]]
            for m, t, ph in zip(ms, te, sp["member_scpot"]):
                m["te"] = t
                m["phi"] = ph
                m["eb_c"] = m["E_b"] + ph if m.get("E_b") and ph is not None else None
                m["eb_te"] = m["eb_c"] / t if m["eb_c"] and t else None

            def med(k):
                v = [m[k] for m in ms if m.get(k) is not None and np.isfinite(m[k])]
                return float(np.median(v)) if v else None

            b = None
            # stale tplot var from an earlier day survives a failed load, reject it
            if bx is not None and not (bx[0][0] <= t0 and bx[0][-1] >= t1):
                print(f"[warn] fgm does not cover {e['episode_id']}")
            elif bx is not None:
                sel = (bx[0] >= t0) & (bx[0] <= t1)
                # single sample episodes fall between fgs points, interp at mid
                b = (float(np.nanmedian(bx[1][sel])) if sel.any()
                     else float(np.interp(0.5 * (t0 + t1), bx[0], bx[1])))
            rows.append({
                "season": f"tail{year}", "episode_id": e["episode_id"], "probe": probe,
                "t_start_ut": e["t_start_ut"], "n_samples": e["n_samples"],
                "e_b": med("E_b"), "e_b_corr": med("eb_c"), "scpot": sp["scpot_med"],
                "scpot_min": sp["scpot_min"], "scpot_max": sp["scpot_max"],
                "de": med("dE"), "eb_over_de": med("Eb_over_dE"),
                "te": med("te"), "eb_over_te": med("eb_te"),
                "x_gsm": p["x_gsm"], "y_gsm": p["y_gsm"], "z_gsm": p["z_gsm"],
                "r_xy": round(float(np.hypot(p["x_gsm"], p["y_gsm"])), 3),
                "bx_gsm": None if b is None else round(b, 3),
            })
    return rows


def fcol(rows, k):
    return np.array([float(r[k]) if r[k] not in (None, "") else np.nan for r in rows])


def plot(rows, out, title):
    x, y = fcol(rows, "x_gsm"), fcol(rows, "y_gsm")
    r, bx = fcol(rows, "r_xy"), fcol(rows, "bx_gsm")
    probe = np.array([row["probe"] for row in rows])
    qs = [("e_b_corr", r"$E_b + q\Phi_{sc}$ (eV)", True), ("de", r"$\delta E_b$ FWHM (eV)", True),
          ("eb_over_de", r"$E_b/\delta E_b$", None)]

    fig = plt.figure(figsize=(15, 16), constrained_layout=True)
    top, bot = fig.subfigures(2, 1, height_ratios=[9, 7])
    axes = top.subplots(2, 3)
    for j, (k, lab, log) in enumerate(qs):
        c = fcol(rows, k)
        ok = np.isfinite(c) & (c > 0)
        axes[0, j].set_title(lab)
        axes[0, j].set_xlabel(r"$X_{GSM}$ ($R_E$)")
        axes[0, j].set_ylabel(r"$Y_{GSM}$ ($R_E$)")
        axes[0, j].invert_xaxis()
        axes[1, j].set_xlabel(r"$R=\sqrt{X^2+Y^2}$ ($R_E$)")
        axes[1, j].set_ylabel(r"$|B_{x,GSM}|$ (nT)")
        if not ok.any():
            for ax in axes[:, j]:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", va="center")
            continue
        lo, hi = np.nanpercentile(c[ok], [1, 99])
        # log only when the quantity spans more than a decade
        if log is None:
            log = hi / lo > 10
        norm = LogNorm(lo, hi) if log else Normalize(lo, hi)
        for i, (xv, yv) in enumerate([(x, y), (r, np.abs(bx))]):
            ax = axes[i, j]
            for pr, mk in MARKERS.items():
                m = ok & (probe == pr) & np.isfinite(xv) & np.isfinite(yv)
                sc = ax.scatter(xv[m], yv[m], c=c[m], norm=norm, cmap="viridis", marker=mk,
                                s=36, edgecolors="k", linewidths=0.4)
            ax.grid(alpha=0.25, lw=0.5)
        cb = top.colorbar(sc, ax=axes[:, j], location="bottom", shrink=0.8, label=lab,
                          extend="both")
        if log:
            dec = range(int(np.floor(np.log10(lo))), int(np.ceil(np.log10(hi))) + 1)
            ticks = [t for t in (m * 10.0 ** p for p in dec for m in (1, 2, 3, 5)) if lo <= t <= hi]
            cb.ax.xaxis.set_major_locator(FixedLocator(ticks))
            cb.ax.xaxis.set_major_formatter(ScalarFormatter())
            cb.ax.xaxis.set_minor_locator(NullLocator())
    handles = [Line2D([], [], ls="", marker=mk, mfc="0.7", mec="k", ms=7,
                      label=f"th{pr} (n={(probe == pr).sum()})") for pr, mk in MARKERS.items()]
    axes[0, 0].legend(handles=handles, loc="best", fontsize=8, framealpha=0.8)
    hists(bot, rows, probe)
    fig.suptitle(title)
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


def step_bins(v):
    # one bin per esa step centered on the median, finer bins split the step clusters
    lv = np.log10(v[np.isfinite(v) & (v > 0)])
    k = np.median(lv)
    n = np.arange(np.round((lv.min() - k) / ESA_DEX), np.round((lv.max() - k) / ESA_DEX) + 2)
    return 10 ** (k + ESA_DEX * (n - 0.5))


def hists(sf, rows, probe):
    qs = [("e_b_corr", r"$E_b + q\Phi_{sc}$ (eV)", "log"),
          ("de", r"$\delta E_b$ FWHM (eV)", "step"),
          ("eb_over_de", r"$E_b/\delta E_b$", "lin"),
          ("te", r"$T_e$ (eV)", "log"),
          ("eb_over_te", r"$(E_b + q\Phi_{sc})/T_e$", "log")]
    axes = sf.subplots(2, 3).ravel()
    for ax, (key, lab, bins) in zip(axes, qs):
        v = fcol(rows, key)
        ok = np.isfinite(v) & (v > 0)
        ax.set_xlabel(lab)
        if not ok.any():
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", va="center")
            continue
        if bins == "step":
            bins = step_bins(v)
        else:
            lo, hi = v[ok].min(), v[ok].max()
            bins = (np.logspace(np.log10(lo), np.log10(hi), 16) if bins == "log"
                    else np.linspace(lo, hi, 16))
        for pr in MARKERS:
            m = ok & (probe == pr)
            ax.hist(v[m], bins=bins, weights=np.full(m.sum(), 1 / max(m.sum(), 1)),
                    histtype="step", lw=2, color=COLORS[pr])
        md = np.median(v[ok])
        ax.axvline(md, color="0.3", ls="--", lw=1)
        ax.text(md, 1.0, f" median {md:.3g}", transform=ax.get_xaxis_transform(),
                va="top", ha="left", fontsize=8, color="0.3")
        if key != "eb_over_de":
            ax.set_xscale("log")
            ax.xaxis.set_major_locator(LogLocator(subs=(1, 2, 5)))
            ax.xaxis.set_major_formatter(FuncFormatter(lambda t, _: f"{t:g}"))
            ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_ylabel("fraction of episodes")
        ax.grid(alpha=0.25, lw=0.5)
    handles = [Line2D([], [], color=COLORS[pr], lw=2, label=f"th{pr} (n={(probe == pr).sum()})")
               for pr in MARKERS]
    axes[-1].axis("off")
    axes[-1].legend(handles=handles, loc="center", title="per-probe fractions, dashed = all-probe median")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seasons", nargs="+", default=["2011", "2012", "2013", "2014", "2015", "2016"])
    ap.add_argument("--tier", type=int, default=1)
    ap.add_argument("--rebuild", action="store_true")
    a = ap.parse_args()
    tag = f"{'-'.join(a.seasons)}_t{a.tier}"
    cache = CAND / "all" / f"eb_map_{tag}.csv"
    srcs = [CAND / f"tail{s}" / "analysis" / n for s in a.seasons
            for n in (f"beam_catalog_{s}.csv", "episodes.jsonl", "episode_positions.jsonl",
                      "episode_scpot.jsonl")]
    # rebuild when any catalog was refit after the cache was written
    stale = cache.exists() and any(p.stat().st_mtime > cache.stat().st_mtime for p in srcs if p.exists())
    if cache.exists() and not a.rebuild and not stale:
        rows = list(csv.DictReader(open(cache)))
    else:
        rows = [row for s in a.seasons for row in build(s, a.tier)]
        with open(cache, "w", newline="") as f:
            w = csv.DictWriter(f, COLS)
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {cache} ({len(rows)} episodes)")
    plot(rows, ROOT / "figures" / f"eb_maps_{tag}.png",
         f"tier-{a.tier} beam episodes, {', '.join(a.seasons)}, marker = probe")


if __name__ == "__main__":
    main()
