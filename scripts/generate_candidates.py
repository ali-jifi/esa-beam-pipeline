# labeled candidate generation, wraps the pipeline without changing it
# emits one jsonl record per coherent-run candidate plus a png cutout each
import argparse
import csv
import json
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from esa_plotting.config import set_data_dir
from esa_plotting.beam_pipeline import (
    load_esd_distribution, load_bfield_dsl, load_bfield_gsm, load_moments,
    compute_pa_spectra, extract_features, classify_beams, ClassifierParams,
    plot_spectra_snapshot, apply_hours, lobe_baseline_mask,
)

CANDIDATES = Path(__file__).resolve().parents[1] / "candidates"


def _fstag(s: str) -> str:
    # windows-safe filename tag from a trange string
    return s.replace("/", "T").replace(":", "")

# strict = production values, single source of truth
STRICT = ClassifierParams()


def make_profile(name: str) -> ClassifierParams:
    if name == "strict":
        return STRICT
    if name == "survey":
        # beam_flux_floor is gone, the sigma bars are its replacement so they take
        # the 0.5x relax, dir_min at 0.5x drops below 1 so that gate is effectively
        # open, thats intended for survey, cast wide and let the strict bitmask sort it
        return replace(STRICT,
                       n_sigma_lo=0.5 * STRICT.n_sigma_lo,
                       n_sigma_hi=0.5 * STRICT.n_sigma_hi,
                       coherent_dir_min=0.5 * STRICT.coherent_dir_min,
                       coherent_min_bins=2)
    raise ValueError(f"unknown profile {name}")


def _bin_arrays(spectra, t, energy_cutoff):
    # per-bin asym, R, and sigma margins, same formulas as extract_features
    valid_e = spectra.energy >= energy_cutoff
    e_valid = spectra.energy[valid_e]
    para = spectra.para[t, valid_e]
    anti = spectra.anti[t, valid_e]
    perp = spectra.perp[t, valid_e]
    pfloor = spectra.perp_floor[t, valid_e]
    psig = spectra.para_sig[t, valid_e]
    asig = spectra.anti_sig[t, valid_e]
    ppsig = spectra.perp_sig[t, valid_e]

    # empty-cone floors, mirrors extract_features so gates see the same values
    if spectra.para_floor is not None:
        pf = spectra.para_floor[t, valid_e]
        emp = ~np.isfinite(para) & np.isfinite(pf)
        para = np.where(emp, pf, para)
        psig = np.where(emp, pf, psig)
    if spectra.anti_floor is not None:
        af = spectra.anti_floor[t, valid_e]
        emp = ~np.isfinite(anti) & np.isfinite(af)
        anti = np.where(emp, af, anti)
        asig = np.where(emp, af, asig)
    emp = ~np.isfinite(perp) & np.isfinite(pfloor)
    perp = np.where(emp, pfloor, perp)
    ppsig = np.where(emp, pfloor, ppsig)

    with np.errstate(invalid="ignore", divide="ignore"):
        denom = para + anti
        asym = np.where(denom > 0, (para - anti) / denom, np.nan)
        perp_eff = np.where(np.isfinite(perp), np.maximum(perp, pfloor), np.nan)
        dom = np.where(np.isfinite(asym) & (asym >= 0), para,
                       np.where(np.isfinite(asym), anti, np.nan))
        dom_sig = np.where(np.isfinite(asym) & (asym >= 0), psig,
                           np.where(np.isfinite(asym), asig, np.nan))
        r = np.where(perp_eff > 0, dom / perp_eff, np.nan)
        sig_d = np.sqrt(psig ** 2 + asig ** 2)
        asym_sg = np.where(sig_d > 0, np.abs(para - anti) / sig_d, np.nan)
        rn = np.sqrt(dom_sig ** 2 + ppsig ** 2)
        r_sg = np.where(rn > 0, (dom - perp_eff) / rn, np.nan)
    signs = np.where(asym > 0, 1, np.where(asym < 0, -1, 0)).astype(int)
    return e_valid, asym, r, asym_sg, r_sg, signs


def _best_run(qual, signs):
    # longest same-sign run with 1-bin gap tolerance, mirrors extract_features
    best, cur, cur_sign, gap = [], [], 0, 0
    for k in range(len(qual)):
        if qual[k] and signs[k] != 0:
            if cur_sign == 0 or signs[k] == cur_sign:
                cur.append(k)
                cur_sign = signs[k]
                gap = 0
            else:
                if len(cur) > len(best):
                    best = cur
                cur, cur_sign, gap = [k], signs[k], 0
        else:
            if cur:
                gap += 1
                if gap > 1:
                    if len(cur) > len(best):
                        best = cur
                    cur, cur_sign, gap = [], 0, 0
    if len(cur) > len(best):
        best = cur
    return best


def _has_adjacent(mask, signs, min_bins):
    # min_bins adjacent true bins with consistent asym sign
    count, sign = 0, 0
    for k in range(len(mask)):
        if mask[k] and signs[k] != 0 and (sign == 0 or signs[k] == sign):
            count += 1
            sign = signs[k]
            if count >= min_bins:
                return True
        else:
            count = 1 if (mask[k] and signs[k] != 0) else 0
            sign = signs[k] if count else 0
    return False


def strict_gates(spectra, t, energy_cutoff, band_lo, band_hi):
    # per-gate pass/fail at strict thresholds over the candidate band +/-1 bin
    e_valid, asym, r, asym_sg, r_sg, signs = _bin_arrays(spectra, t, energy_cutoff)
    lo = max(0, band_lo - 1)
    hi = min(len(asym), band_hi + 2)
    sl = slice(lo, hi)
    with np.errstate(invalid="ignore"):
        asym_q = np.isfinite(asym) & (np.abs(asym) >= STRICT.coherent_asym_min)
        dir_q = np.isfinite(r) & (r >= STRICT.coherent_dir_min)
        # lo bar only, hysteresis needs neighbor timesteps and isnt mirrored here
        sig_q = (np.isfinite(asym_sg) & (asym_sg >= STRICT.n_sigma_lo) &
                 np.isfinite(r_sg) & (r_sg >= STRICT.n_sigma_lo))
    mb = STRICT.coherent_min_bins
    s = signs[sl]
    return {
        "sig": _has_adjacent(sig_q[sl], s, mb),
        "asym": _has_adjacent(asym_q[sl], s, mb),
        "dir": _has_adjacent(dir_q[sl], s, mb),
        "min_bins": _has_adjacent((asym_q & dir_q & sig_q)[sl], s, mb),
    }


def beta_at(tc, moments, b_gsm):
    # plasma beta from ion moments and |B|, beta = 0.403 n[cc] T[eV] / B[nT]^2
    if b_gsm is None or "density" not in moments or "temperature" not in moments:
        return None
    n = np.interp(tc, moments["density_times"], moments["density"])
    T = np.interp(tc, moments["temp_times"], moments["temperature"])
    bv = [np.interp(tc, b_gsm[0], b_gsm[1][:, i]) for i in range(3)]
    b2 = bv[0] ** 2 + bv[1] ** 2 + bv[2] ** 2
    if not (np.isfinite(b2) and b2 > 0):
        return None
    beta = 0.4027 * n * T / b2
    return float(beta) if np.isfinite(beta) else None


LOBE_BETA_MAX = 0.1


def build_records(probe, trange, tag, profile_name, spectra, features, cls,
                  energy_cutoff, b_gsm, ne_ni, moments, te=None):
    params = make_profile(profile_name)
    records = []
    for t in range(len(spectra.times)):
        if not features.coherent_ok[t]:
            continue
        # active-profile band via the same qual + run logic as extract_features
        e_valid, a, r, asg, rsg, sg = _bin_arrays(spectra, t, energy_cutoff)
        with np.errstate(invalid="ignore"):
            qual = (np.isfinite(a) & np.isfinite(r) &
                    (np.abs(a) >= params.coherent_asym_min) &
                    (r >= params.coherent_dir_min) &
                    np.isfinite(asg) & (asg >= params.n_sigma_lo) &
                    np.isfinite(rsg) & (rsg >= params.n_sigma_lo))
        run = _best_run(qual, sg)
        if run:
            band_lo, band_hi = min(run), max(run)
        else:
            # shouldnt happen, fall back to the bin nearest e_beam
            band_lo = band_hi = int(np.argmin(np.abs(e_valid - features.e_beam[t])))

        el = features.e_line[t]
        e_idx = (int(np.argmin(np.abs(spectra.energy - el)))
                 if np.isfinite(el)
                 else int(np.argmin(np.abs(spectra.energy - features.e_beam[t]))))
        tc = float(spectra.times[t])
        ut = datetime.fromtimestamp(tc, tz=timezone.utc)

        bx_sign = None
        if b_gsm is not None:
            bx = float(np.interp(tc, b_gsm[0], b_gsm[1][:, 0]))
            bx_sign = int(np.sign(bx))
        ratio = None
        if ne_ni is not None:
            ratio = float(np.interp(tc, ne_ni[0], ne_ni[1]))
        beta = beta_at(tc, moments, b_gsm)
        te_ev = float(np.interp(tc, te[0], te[1])) if te is not None else None

        def f(x):
            return float(x) if np.isfinite(x) else None

        records.append({
            "candidate_id": f"{probe}_{_fstag(trange[0])}_"
                            f"{ut.strftime('%H%M%S')}_{e_idx}",
            "probe": probe,
            "t_center": tc,
            "t_center_ut": ut.isoformat(),
            "trange": list(trange),
            "interval_tag": tag,
            "profile": profile_name,
            "t_idx": t,
            "features": {
                "R": f(features.r_beam[t]),
                "E_b": f(features.e_line[t]),
                "dE": f(features.de_line[t]),
                "Eb_over_dE": f(features.eb_over_de[t]),
                "peak_prom": f(features.peak_prom[t]),
                "sig_margin": f(features.sig_margin[t]),
                "asymmetry": f(features.asymmetry[t]),
                "asym_baseline": f(features.asym_baseline[t]),
                "asym_dev": f(features.asym_dev[t]),
                "flux_z": f(features.flux_z[t]),
                "flux_z_perp": f(features.flux_z_perp[t]),
                "duration_steps": f(features.duration[t]),
                "chain_e_slope": f(features.chain_e_slope[t]),
                "chain_e_scatter": f(features.chain_e_scatter[t]),
                "width": f(features.width[t]),
                "para_to_omni": f(features.para_to_omni[t]),
                "hyst_promoted": bool(features.hyst_promoted[t]),
                "beam_score": f(cls.beam_score[t]),
            },
            "gates": strict_gates(spectra, t, energy_cutoff, band_lo, band_hi),
            "is_beam": bool(cls.is_beam[t]),
            "direction": int(cls.beam_direction[t]),
            "context": {
                "Bx_sign": bx_sign,
                "ne_ni_ratio": ratio,
                "perp_depleted": bool(features.perp_depleted[t]),
                "cone_floored": bool(features.cone_floored[t]),
                # dir*bx flow sense, -1 tailward, +1 earthward
                "dir_x_bx": (int(cls.beam_direction[t] * bx_sign)
                             if bx_sign is not None and cls.beam_direction[t] != 0
                             else None),
                "beta": beta,
                "in_lobe": (beta < LOBE_BETA_MAX) if beta is not None else None,
                "te_ev": te_ev,
            },
            "label": None,
            "label_source": None,
            "label_confidence": None,
        })
    return records


def render_cutout(spectra, rec, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    t_idx = rec["t_idx"]
    tc = rec["t_center"]
    fig = plt.figure(figsize=(10, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.2], width_ratios=[2, 1])

    # spectrogram +/-10 min, fixed color scale and axes across all cutouts
    ax0 = fig.add_subplot(gs[0, :])
    win = (spectra.times >= tc - 600) & (spectra.times <= tc + 600)
    idxs = np.where(win)[0]
    e_ok = spectra.energy > 0
    e_plot = spectra.energy[e_ok]
    times_dt = [datetime.fromtimestamp(x, tz=timezone.utc)
                for x in spectra.times[idxs]]
    z = spectra.omni[np.ix_(idxs, e_ok)].T
    z = np.where(z > 0, z, np.nan)
    pcm = ax0.pcolormesh(times_dt, e_plot, z,
                         norm=plt.matplotlib.colors.LogNorm(vmin=1e3, vmax=1e8),
                         cmap="jet", shading="auto")
    ax0.set_yscale("log")
    ax0.set_ylim(5, 30000)
    ax0.set_ylabel("Energy [eV]")
    ax0.axvline(datetime.fromtimestamp(tc, tz=timezone.utc),
                color="w", ls="--", lw=1.5)
    if np.isfinite(rec["features"]["E_b"] or np.nan):
        ax0.plot(datetime.fromtimestamp(tc, tz=timezone.utc),
                 rec["features"]["E_b"], "wo", mec="k", ms=8)
    ax0.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    fig.colorbar(pcm, ax=ax0, label="eflux", pad=0.01)
    ax0.set_title(rec["candidate_id"])

    # 1d spectra at t_center via the existing certified snapshot plotter
    ax1 = fig.add_subplot(gs[1, 0])
    plot_spectra_snapshot(spectra, t_idx, label=rec["t_center_ut"][11:19], ax=ax1)
    ax1.set_xlim(5, 30000)
    ax1.set_ylim(1e3, 1e8)

    # feature vector + gate bitmask text box
    ax2 = fig.add_subplot(gs[1, 1])
    ax2.axis("off")
    ft, gt, cx = rec["features"], rec["gates"], rec["context"]
    # one key column and right-justified values across all sections so the
    # numbers land inline, monospace makes the right edges a clean rail
    kw = max(len(k) for k in [*ft, *gt, *cx]) + 1
    def fmt(v):
        s = ("nan" if v is None
             else f"{v:.3g}" if isinstance(v, float) else str(v))
        return f"{s:>9}"
    txt = ["features"]
    txt += [f"  {k:<{kw}}= {fmt(v)}" for k, v in ft.items()]
    txt += ["gates (strict)"]
    txt += [f"  {k:<{kw}}= {fmt('PASS' if v else 'FAIL')}" for k, v in gt.items()]
    txt += ["context"]
    txt += [f"  {k:<{kw}}= {fmt(v)}" for k, v in cx.items()]
    txt += [f"is_beam({rec['profile']}) = {rec['is_beam']}  dir={rec['direction']}"]
    ax2.text(0.0, 1.0, "\n".join(txt), va="top", family="monospace", fontsize=9)

    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


def _load_ne_ni(probe, moments):
    # electron/ion density ratio, peer vars come in with the same esa load
    from pyspedas import get_data
    try:
        d = get_data(f"th{probe}_peer_density")
        if d is None or "density" not in moments:
            return None
        ni = np.interp(d.times, moments["density_times"], moments["density"])
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = np.where(ni > 0, d.y / ni, np.nan)
        return d.times, ratio
    except Exception:
        return None


def _load_te(probe):
    # electron temp, papers lobe indicator is cold (<500 ev) electrons
    from pyspedas import get_data
    try:
        d = get_data(f"th{probe}_peer_avgtemp")
        return (d.times, d.y) if d is not None else None
    except Exception:
        return None


def run_interval(probe, trange, tag, profile_name, data_dir, out_dir,
                 energy_cutoff=30.0, cutouts=True, datatype="peif"):
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{probe}_{_fstag(trange[0])}_{profile_name}"
    jsonl = out_dir / f"{stem}.jsonl"
    # interval-level checkpoint, downloads are slow so the job must restart clean
    if jsonl.exists():
        print(f"[skip] {jsonl.name} exists")
        return None

    params = make_profile(profile_name)
    dist = load_esd_distribution(probe, trange, data_dir, datatype=datatype)
    b_times, b_dsl = load_bfield_dsl(probe, trange, data_dir)
    moments = load_moments(probe, trange, data_dir)
    spectra = compute_pa_spectra(dist, b_times, b_dsl)
    b_mask = lobe_baseline_mask(spectra.times, moments, b_times, b_dsl)
    features = extract_features(spectra, moments,
                                baseline_mask=b_mask,
                                energy_cutoff_low=energy_cutoff,
                                pa_coverage_threshold=params.min_coverage,
                                coherent_asym_min=params.coherent_asym_min,
                                coherent_dir_min=params.coherent_dir_min,
                                coherent_min_bins=params.coherent_min_bins,
                                n_sigma_lo=params.n_sigma_lo,
                                n_sigma_hi=params.n_sigma_hi,
                                beam_e_max=params.beam_e_max,
                                peak_width_max=params.peak_width_max,
                                peak_wlen=params.peak_wlen)
    cls = classify_beams(features, params)

    b_gsm = None
    try:
        b_gsm = load_bfield_gsm(probe, trange, data_dir)
    except Exception as e:
        print(f"[warn] b_gsm unavailable: {e}")
    ne_ni = _load_ne_ni(probe, moments)
    te = _load_te(probe)

    records = build_records(probe, trange, tag, profile_name, spectra,
                            features, cls, energy_cutoff, b_gsm, ne_ni, moments,
                            te=te)

    tmp = jsonl.with_suffix(".jsonl.tmp")
    with open(tmp, "w") as fh:
        for rec in records:
            fh.write(json.dumps(rec) + "\n")
    tmp.rename(jsonl)
    print(f"[ok] {jsonl.name}: {len(records)} candidates "
          f"({sum(r['is_beam'] for r in records)} is_beam under {profile_name})")

    if cutouts:
        cut_dir = out_dir / "cutouts"
        cut_dir.mkdir(exist_ok=True)
        for rec in records:
            render_cutout(spectra, rec, cut_dir / f"{rec['candidate_id']}.png")
        print(f"[ok] {len(records)} cutouts -> {cut_dir}")
    return records


def main():
    p = argparse.ArgumentParser(description="labeled beam candidate generation")
    p.add_argument("--events", help="csv with probe,trange_start,trange_end,interval_tag")
    p.add_argument("--probe", default="a")
    p.add_argument("--trange", nargs=2, default=["2019-05-01", "2019-05-02"])
    p.add_argument("--hours", nargs=2, type=int)
    p.add_argument("--tag", default="unknown",
                   choices=["beam_period", "plasma_sheet", "unknown"])
    p.add_argument("--profile", default="strict", choices=["strict", "survey"])
    p.add_argument("--energy-cutoff", type=float, default=30.0)
    p.add_argument("--out", default=str(CANDIDATES))
    p.add_argument("--no-cutouts", action="store_true")
    p.add_argument("--datatype", default="peif", choices=["peif", "peir"])
    args = p.parse_args()

    data_dir = set_data_dir()
    out_dir = Path(args.out)

    if args.events:
        with open(args.events, newline="") as fh:
            rows = list(csv.DictReader(fh))
        fails = 0
        for row in rows:
            trange = [row["trange_start"], row["trange_end"]]
            try:
                run_interval(row["probe"], trange, row["interval_tag"],
                             args.profile, data_dir, out_dir,
                             energy_cutoff=args.energy_cutoff,
                             cutouts=not args.no_cutouts,
                             datatype=args.datatype)
            except Exception as e:
                # one bad interval must not kill a multi-day batch
                fails += 1
                print(f"[fail] {row['probe']} {trange[0]}: {type(e).__name__}: {e}")
        if fails:
            print(f"[warn] {fails}/{len(rows)} intervals failed")
    else:
        trange = apply_hours(args.trange, args.hours)
        run_interval(args.probe, trange, args.tag, args.profile,
                     data_dir, out_dir,
                     energy_cutoff=args.energy_cutoff,
                     cutouts=not args.no_cutouts,
                     datatype=args.datatype)


if __name__ == "__main__":
    main()
