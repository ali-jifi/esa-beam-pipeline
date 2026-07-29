# probe why median E_beam/dE falls as the R gate tightens
import argparse
import numpy as np
from scipy.stats import spearmanr

from esa_plotting.config import set_data_dir
from esa_plotting.beam_pipeline import (
    load_esd_distribution, load_bfield_dsl, load_moments,
    compute_pa_spectra, extract_features, classify_beams, ClassifierParams,
    apply_hours,
)


def _med(x):
    x = x[np.isfinite(x)]
    return float(np.median(x)) if x.size else float("nan")


def _runs(mask):
    # contiguous flagged stretches, one beam interval each
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return []
    return np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)


def _spearman_block(r, pairs, indent="  "):
    for name, y in pairs:
        rho, pv = spearmanr(r, y)
        print(f"{indent}spearman(r_beam, {name:>5}) = {rho:+.3f}  (p={pv:.3g})")


def _extract(spectra, moments, params, thr, cutoff):
    feat = extract_features(spectra, moments,
                            energy_cutoff_low=cutoff,
                            pa_coverage_threshold=params.min_coverage,
                            coherent_asym_min=params.coherent_asym_min,
                            coherent_dir_min=thr,
                            coherent_min_bins=params.coherent_min_bins,
                            n_sigma_lo=params.n_sigma_lo,
                            n_sigma_hi=params.n_sigma_hi,
                            pair_e_max=params.pair_e_max,
                            peak_width_max=params.peak_width_max,
                            peak_wlen=params.peak_wlen)
    return feat, classify_beams(feat, params)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--probe", default="a")
    # default event, differs from pipeline default on purpose
    p.add_argument("--trange", nargs=2, default=["2009-02-27", "2009-02-28"])
    p.add_argument("--thresholds", nargs="+", type=float, default=[1.0, 1.2, 1.5, 2.0])
    p.add_argument("--energy-cutoff", type=float, default=30.0)
    p.add_argument("--hours", nargs=2, type=int, metavar=("H_START", "H_END"),
                   help="hour window on the start date, e.g. 6 12")
    args = p.parse_args()

    trange = apply_hours(args.trange, args.hours)
    data_dir = set_data_dir()
    params = ClassifierParams()

    dist = load_esd_distribution(args.probe, trange, data_dir)
    b_times, b_dsl = load_bfield_dsl(args.probe, trange, data_dir)
    moments = load_moments(args.probe, trange, data_dir)
    spectra = compute_pa_spectra(dist, b_times, b_dsl)

    # one extract per threshold, reused by every section below
    results = {thr: _extract(spectra, moments, params, thr, args.energy_cutoff)
               for thr in args.thresholds}

    print(f"\n=== THEMIS-{args.probe.upper()} {trange[0]} -> {trange[1]} ===")
    print(f"{'R>=':>6} {'n':>4} {'med Eb/dE':>10} {'med dE':>9} {'med Eb':>9} "
          f"{'med r_beam':>11} {'perp_depl%':>11}")
    for thr in args.thresholds:
        feat, cls = results[thr]
        m = cls.is_beam
        n = int(m.sum())
        eb, de, el = feat.eb_over_de[m], feat.de_line[m], feat.e_line[m]
        r, pdep = feat.r_beam[m], feat.perp_depleted[m]
        pct_depl = 100.0 * np.mean(pdep) if n else float("nan")
        print(f"{thr:>6.2f} {n:>4d} {_med(eb):>10.3f} {_med(de):>9.1f} "
              f"{_med(el):>9.1f} {_med(r):>11.3f} {pct_depl:>10.1f}%")

    # split the median trend into selection vs re-measurement
    thr_lo, thr_hi = args.thresholds[0], args.thresholds[-1]
    if thr_hi > thr_lo:
        feat_lo, cls_lo = results[thr_lo]
        feat_hi, cls_hi = results[thr_hi]
        both = cls_lo.is_beam & cls_hi.is_beam
        print(f"\nselection vs re-measurement, R>={thr_lo} vs R>={thr_hi}")
        print(f"  {int(cls_lo.is_beam.sum())} beams at loose gate, "
              f"{int(both.sum())} survive tight gate")
        if both.any():
            eb_all = feat_lo.eb_over_de[cls_lo.is_beam]
            eb_lo = feat_lo.eb_over_de[both]
            eb_hi = feat_hi.eb_over_de[both]
            print(f"  selection: med Eb/dE {_med(eb_all):.3f} all -> "
                  f"{_med(eb_lo):.3f} survivors, both at loose gate")
            print(f"  re-measurement: survivors {_med(eb_lo):.3f} loose -> "
                  f"{_med(eb_hi):.3f} tight")
            ok = np.isfinite(eb_lo) & np.isfinite(eb_hi)
            if ok.any():
                changed = eb_lo[ok] != eb_hi[ok]
                if changed.any():
                    dmed = float(np.median(eb_hi[ok][changed] - eb_lo[ok][changed]))
                    print(f"  {100 * np.mean(changed):.0f}% of survivor lines "
                          f"re-measured, med paired delta {dmed:+.3f}")
                else:
                    print("  no survivor lines re-measured, trend is pure selection")

    # correlations across the loose-gate beam population
    feat, cls = results[thr_lo]
    m = cls.is_beam
    r = feat.r_beam[m]
    eb, de, el = feat.eb_over_de[m], feat.de_line[m], feat.e_line[m]
    ok = np.isfinite(r) & np.isfinite(eb) & np.isfinite(de) & np.isfinite(el)
    r, eb, de, el = r[ok], eb[ok], de[ok], el[ok]
    print(f"\nbeam population at R>={thr_lo}, n={r.size}")
    if r.size >= 5:
        # per-timestep, descriptive only, consecutive steps arent independent
        _spearman_block(r, (("Eb/dE", eb), ("dE", de), ("Eb", el)))
        pdep = feat.perp_depleted[m][ok]
        if pdep.any() and (~pdep).any():
            print(f"  median r_beam: perp_depleted={_med(r[pdep]):.3f}  "
                  f"not_depleted={_med(r[~pdep]):.3f}  "
                  f"(depleted frac {100 * np.mean(pdep):.0f}%)")
            # clamped r_beam is a lower bound not a measurement, redo without it
            if (~pdep).sum() >= 5:
                print(f"  excluding depleted, n={int((~pdep).sum())}")
                _spearman_block(r[~pdep], (("Eb/dE", eb[~pdep]),))
        else:
            print(f"  perp_depleted uniform ({100 * np.mean(pdep):.0f}%), no split")
    else:
        print("  too few beams for correlation")

    # per-interval medians, one sample per beam so p-values arent inflated
    runs = _runs(m)
    print(f"\nper-interval correlation, {len(runs)} beam intervals")
    if len(runs) >= 5:
        r_run = np.array([_med(feat.r_beam[ix]) for ix in runs])
        eb_run = np.array([_med(feat.eb_over_de[ix]) for ix in runs])
        ok_run = np.isfinite(r_run) & np.isfinite(eb_run)
        if ok_run.sum() >= 5:
            _spearman_block(r_run[ok_run], (("Eb/dE", eb_run[ok_run]),))
        else:
            print("  too few finite intervals")
    else:
        print("  too few intervals for inference, pool more dates via the beam csvs")


if __name__ == "__main__":
    main()
