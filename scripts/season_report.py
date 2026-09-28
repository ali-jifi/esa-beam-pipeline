# aggregate season outputs into one interactive html report
# python scripts/season_report.py --seasons 2014 2015 2016
import argparse
import csv
import html
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

import build_episode_catalog as bc
import rescore_decks as rd

ROOT = Path(__file__).resolve().parents[1]
CAND = ROOT / "candidates"
# claude design export with placeholders, see render()
TEMPLATE = Path(__file__).parent / "season_report_design.html"
BANDS = [(40, 60), (60, 120), (120, 250), (250, 400), (400, 600)]


def fnum(v):
    try:
        x = float(v)
        return x if np.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def gated(e):
    te, eb, beta = fnum(e["te_med"]), fnum(e["e_b_med"]), fnum(e["beta_med"])
    return (e["dir_x_bx"] == str(bc.DIR_GATE) and eb is not None and eb < bc.EB_MAX
            and beta is not None and beta < bc.BETA_MAX and (te is None or te < bc.TE_MAX))


def med(v):
    v = [x for x in v if x is not None]
    return round(float(np.median(v)), 4) if v else None


def season(year, cov, decks):
    s = f"tail{year}"
    ana = CAND / s / "analysis"
    eps = {r["episode_id"]: r for r in csv.DictReader(open(ana / "episodes.csv"))}
    cat = list(csv.DictReader(open(ana / f"beam_catalog_{year}.csv")))
    pos = {}
    for line in open(ana / "episode_positions.jsonl"):
        p = json.loads(line)
        pos[p["episode_id"]] = p
    links = list(csv.DictReader(open(ana / f"linked_events_{year}.csv")))
    t1 = [r for r in cat if r["tier"] == "1"]
    t2 = [r for r in cat if r["tier"] == "2"]
    c = [r for r in cov if r["season"] == s]
    hours = sum(float(r["hours_angle"]) for r in c)
    ev1 = [l for l in links if int(l["n_tier1"]) > 0]
    multi = [l for l in ev1 if int(l["n_probes"]) >= 2]
    t1_in_multi = sum(int(l["n_tier1"]) for l in multi)
    day1 = Counter(r["t_start_ut"][:10] for r in t1)
    day2 = Counter(r["t_start_ut"][:10] for r in t2)
    dur = [fnum(r["duration_s"]) or 0 for r in t1]
    by_day_time = defaultdict(float)
    for r in t1:
        by_day_time[r["t_start_ut"][:10]] += fnum(r["duration_s"]) or 0
    top3 = sum(sorted(by_day_time.values(), reverse=True)[:3]) / max(sum(dur), 1)
    ebs = [fnum(r["e_b_med"]) for r in t1]
    bands = [sum(1 for e in ebs if e is not None and lo <= e < hi) for lo, hi in BANDS]
    dirs = Counter(eps[r["episode_id"]]["direction"] for r in t1)
    probes = Counter(eps[r["episode_id"]]["probe"] for r in t1)
    points = []
    for r in t1:
        p = pos.get(r["episode_id"])
        if p:
            points.append({"id": r["episode_id"], "t": r["t_start_ut"][:19].replace("T", " "),
                           "x": p["x_gsm"], "y": p["y_gsm"], "z": p["z_gsm"],
                           "eb": fnum(r["e_b_med"]), "te": fnum(r["te_med"]), "beta": fnum(r["beta_med"]),
                           "dur": fnum(r["duration_s"]), "prob": fnum(r["prob"]), "n": int(r["n_samples"])})
    days = sorted(set(day1) | set(day2) | {e["t_start_ut"][:10] for e in eps.values()})
    events = [{"id": l["event_id"], "probes": l["probes"], "t": l["t_start_ut"][:16].replace("T", " "),
               "dur": round(float(l["duration_s"]) / 60, 1), "n1": int(l["n_tier1"]), "n2": int(l["n_tier2"]),
               "eb": f'{float(l["e_b_min"]):.0f}-{float(l["e_b_max"]):.0f}', "prob": fnum(l["prob_max"])}
              for l in ev1]
    dk = [r for r in decks if r["season"] == s]
    d1 = [r for r in dk if r["tier"] == 1]
    dpos = [r for r in dk if r["label"] == "positive"]
    lo, hi = rd.ci(sum(r["label"] == "negative" for r in d1), len(d1))
    return {
        "season": year, "intervals": sum(int(r["n_intervals"]) for r in c),
        "processed": sum(int(r["n_processed"]) for r in c), "hours": round(hours, 1),
        "episodes": len(eps), "gated": sum(gated(e) for e in eps.values()),
        "t1": len(t1), "t2": len(t2), "rate": round(100 * len(t1) / max(hours, 1), 2),
        "beam_h": round(sum(dur) / 3600, 2), "days": len(day1), "top3_share": round(top3, 3),
        "med_dur": med(dur), "med_eb": med(ebs), "med_te": med([fnum(r["te_med"]) for r in t1]),
        "med_beta": med([fnum(r["beta_med"]) for r in t1]),
        "med_prob": med([fnum(r["prob"]) for r in t1]),
        "p99": sum((fnum(r["prob"]) or 0) >= 0.99 for r in t1) / max(len(t1), 1),
        "par": dirs.get("1", 0), "anti": dirs.get("-1", 0), "probes": dict(sorted(probes.items())),
        "events": len(ev1), "multi": len(multi), "multi_frac": round(t1_in_multi / max(len(t1), 1), 3),
        "bands": bands, "points": points,
        "daily": [{"d": d, "t1": day1.get(d, 0), "t2": day2.get(d, 0)} for d in days],
        "top_events": sorted(events, key=lambda e: (-e["n1"], -e["dur"]))[:40],
        "val": {"decks": sorted({r["deck"] for r in dk}), "who": sorted({r["who"] for r in dk}),
                "decided": len(dk), "t1": len(d1), "t1_neg": sum(r["label"] == "negative" for r in d1),
                "ci": [round(lo, 3), round(hi, 3)], "pos": len(dpos),
                "pos_t1": sum(r["tier"] == 1 for r in dpos), "pos_t12": sum(r["tier"] in (1, 2) for r in dpos),
                "roles": sorted({r["role"] for r in dk})},
    }


def pct(x):
    # half up like the page, python rounds half to even
    return int(100 * x + 0.5)


def joined(parts):
    return parts[0] if len(parts) == 1 else ", ".join(parts[:-1]) + " and " + parts[-1]


def insights(S):
    # short findings in the design's voice, one row each
    out = []
    top, low = max(S, key=lambda s: s["rate"]), min(S, key=lambda s: s["rate"])
    if len(S) > 1:
        out.append(f"{top['season']} has the highest tier-1 rate: {top['rate']} per 100 h, against {low['rate']} in {low['season']}.")
    else:
        out.append(f"{top['season']} tier-1 rate: {top['rate']} per 100 h.")
    conc = max(S, key=lambda s: s["top3_share"])
    out.append(f"Beam time clusters in a few days. In {conc['season']}, the three busiest days hold "
               f"{pct(conc['top3_share'])}% of it.")
    out.append("Beams appeared on only " + joined(
        [f"{s['days']} of {len(s['daily'])}{' scanned days' if i == 0 else ''} in {s['season']}" for i, s in enumerate(S)]) + ".")
    durs = [s["med_dur"] for s in S if s["med_dur"] is not None]
    head = "Median beam duration is about a minute" if durs and all(30 <= d <= 120 for d in durs) else "Median beam duration"
    out.append(head + ": " + joined([f"{s['med_dur']:.0f} s in {s['season']}" for s in S if s["med_dur"] is not None]) + ".")
    low = sorted(S, key=lambda s: s["bands"][0] / max(s["t1"], 1), reverse=True)
    if len(low) > 1 and low[0]["t1"]:
        a, b = low[0], low[1]
        out.append(f"{a['season']} beams are the coldest and slowest: {pct(a['bands'][0] / a['t1'])}% sit at "
                   f"{BANDS[0][0]} to {BANDS[0][1]} eV, against {pct(b['bands'][0] / max(b['t1'], 1))}% in {b['season']}.")
    mb = max(S, key=lambda s: s["multi_frac"])
    out.append(f"{pct(mb['multi_frac'])}% of {mb['season']} tier-1 episodes are multi-probe, "
               f"but the probes fly close together.")
    hp = sorted([s for s in S if s["t1"]], key=lambda s: s["p99"], reverse=True)
    if len(hp) > 1:
        out.append(f"{pct(hp[0]['p99'])}% of {hp[0]['season']} tier-1 episodes score a model probability of 0.99 or more, "
                   f"against {pct(hp[-1]['p99'])}% in {hp[-1]['season']}.")
    lab = [s for s in S if s["val"]["t1"]]
    if lab:
        roles = {r for s in lab for r in s["val"]["roles"]}
        caveat = ", but every deck is a training deck, so this is optimistic" if roles == {"training"} else ""
        if all(s["val"]["t1_neg"] == 0 for s in lab):
            out.append(f"No labeled tier-1 negatives in any season{caveat}.")
        else:
            parts = ", ".join(f"{s['season']} {s['val']['t1_neg']} of {s['val']['t1']}" for s in lab)
            out.append(f"Labeled tier-1 negatives: {parts}{caveat}.")
    return out


def findings_html(items):
    rows = []
    for i, t in enumerate(items):
        top = "" if i == 0 else "; border-top: 1px solid var(--line)"
        rows.append(f'<div class="row" style="display: grid; grid-template-columns: 60px minmax(0, 1fr){top}">'
                    f'<div style="border-right: 1px solid var(--line); padding: 12px 0; text-align: center; '
                    f'font-family: var(--mono); font-size: 12px">{i + 1:02d}</div>'
                    f'<div style="padding: 11px 16px">{html.escape(t)}</div></div>')
    return "\n".join(rows)


def render(data, seasons):
    # fill the claude design bundle, the page lives as a json string inside it
    b = TEMPLATE.read_text(encoding="utf-8")
    m = re.search(r'(<script type="__bundler/template">\s*)(.*?)(\s*</script>)', b, re.S)
    page = json.loads(m.group(2))
    ys = seasons[0] if len(seasons) == 1 else f"{seasons[0]} to {seasons[-1]}"
    md = data["model"]
    footer = (f"python scripts/season_report.py --seasons {' '.join(seasons)} · model {md['feats']} inputs "
              f"· cut {md['thr']} · prominence {md['prom']} · tier 1: {md['min_samples']}+ samples, "
              f"E_b > {md['eb_floor']:g} eV")
    for k, v in (("__DATA__", json.dumps(data)), ("__SUBTITLE__", html.escape(f"Tail seasons {ys}")),
                 ("__FINDINGS__", findings_html(data["insights"])), ("__FOOTER__", html.escape(footer))):
        assert page.count(k) == 1, k
        page = page.replace(k, v)
    # a literal </ would close the bundle script tag early
    enc = json.dumps(page, ensure_ascii=False).replace("</", "<\\u002F")
    return b[:m.start(2)] + enc + b[m.end(2):]


def main():
    p = argparse.ArgumentParser(description="interactive season report")
    p.add_argument("--seasons", nargs="+", default=["2014", "2015", "2016"])
    p.add_argument("--out", default=None)
    args = p.parse_args()
    cov = list(csv.DictReader(open(CAND / "all" / "coverage.csv")))
    decks = rd.load()
    S = [season(y, cov, decks) for y in args.seasons]
    data = {"seasons": S, "bands": [f"{lo}-{hi}" for lo, hi in BANDS], "insights": insights(S),
            "model": {"thr": bc.PROB_THR, "prom": bc.PROM_BAR, "min_samples": bc.MIN_SAMPLES,
                      "eb_floor": bc.EB_FLOOR, "feats": len(bc.FEATS)}}
    out = Path(args.out or CAND / "all" / f"report_{'_'.join(args.seasons)}.html")
    out.write_text(render(data, args.seasons), encoding="utf-8", newline="")
    for s in S:
        print(f"{s['season']}: t1 {s['t1']} t2 {s['t2']} rate {s['rate']} beam_h {s['beam_h']} multi {s['multi_frac']}")
    print(f"[ok] wrote {out}")


if __name__ == "__main__":
    main()
