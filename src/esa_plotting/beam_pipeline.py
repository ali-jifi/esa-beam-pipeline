# ion beam detection pipeline for themis esa data

from __future__ import annotations

import glob
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from cdflib import CDF
from pyspedas import get_data
from pyspedas.projects import themis
from scipy.interpolate import interp1d
from scipy.signal import find_peaks


@dataclass
class ESDDistribution:
    # 3d ion dist from themis esa esd files
    times: np.ndarray          # (ntime,) unix timestamps
    eflux: np.ndarray          # (ntime, 32, 176) energy flux
    energy: np.ndarray         # (32,) energy vals in ev for active mode
    theta: np.ndarray          # (32, 176) polar angle in instrument coords, deg
    phi: np.ndarray            # (32, 176) azimuthal angle in instrument coords, deg
    domega: np.ndarray         # (32, 176) solid angle element
    bins_mask: np.ndarray      # (ntime, 32, 176) valid bin mask
    phi_offset: np.ndarray     # (ntime,) spin-phase correction
    en_ind: np.ndarray         # (ntime,) energy mode idx per timestep
    an_ind: np.ndarray         # (ntime,) angle mode idx per timestep
    onecount: np.ndarray       # (ntime, 32, 176) one-count eflux level per bin


def load_esd_distribution(probe: str, trange: list[str], data_dir: str) -> ESDDistribution:
    # reads cdf directly cuz pyspedas doesnt expose angle/energy lookup tables
    import os
    os.environ["THM_DATA_DIR"] = data_dir

    themis.esd(probe=probe, trange=trange, datatype="peif",
               time_clip=True, downloadonly=True)

    pattern = str(
        Path(data_dir)
        / f"th{probe}" / "l2" / "esd" / "*"
        / f"th{probe}_l2_esa_peif_*.cdf"
    )
    cdf_files = sorted(glob.glob(pattern))
    if not cdf_files:
        raise FileNotFoundError(f"No ESD CDF files found matching {pattern}")

    all_times, all_eflux, all_bins = [], [], []
    all_en_ind, all_an_ind, all_phi_offset = [], [], []
    all_eff, all_integ_t = [], []
    energy_table = phi_table = theta_table = domega_table = None
    gf_table = geom_factor = None

    for cdf_path in cdf_files:
        cdf = CDF(cdf_path)

        epoch = cdf.varget("epoch")
        from cdflib.epochs import CDFepoch
        unix_times = CDFepoch.unixtime(epoch)

        t0_str, t1_str = trange
        from datetime import datetime, timezone
        t0_unix = datetime.strptime(t0_str.split("/")[0], "%Y-%m-%d").replace(
            tzinfo=timezone.utc
        ).timestamp()
        if "/" in t0_str:
            parts = t0_str.split("/")[1].split(":")
            t0_unix += int(parts[0]) * 3600
            if len(parts) > 1:
                t0_unix += int(parts[1]) * 60
            if len(parts) > 2:
                t0_unix += float(parts[2])

        t1_unix = datetime.strptime(t1_str.split("/")[0], "%Y-%m-%d").replace(
            tzinfo=timezone.utc
        ).timestamp()
        if "/" in t1_str:
            parts = t1_str.split("/")[1].split(":")
            t1_unix += int(parts[0]) * 3600
            if len(parts) > 1:
                t1_unix += int(parts[1]) * 60
            if len(parts) > 2:
                t1_unix += float(parts[2])

        mask = (unix_times >= t0_unix) & (unix_times <= t1_unix)
        if not np.any(mask):
            continue

        unix_times = np.asarray(unix_times)[mask]
        eflux = cdf.varget("eflux")[mask]
        bins_data = cdf.varget("bins")[mask]
        en_ind = cdf.varget("en_ind")[mask]
        an_ind = cdf.varget("an_ind")[mask]
        phi_offset_data = cdf.varget("phi_offset")[mask]
        eff_data = cdf.varget("eff")[mask]
        integ_t_data = cdf.varget("integ_t")[mask]

        all_times.append(unix_times)
        all_eflux.append(eflux)
        all_bins.append(bins_data)
        all_en_ind.append(en_ind)
        all_an_ind.append(an_ind)
        all_phi_offset.append(phi_offset_data)
        all_eff.append(eff_data)
        all_integ_t.append(integ_t_data)

        if energy_table is None:
            energy_table = cdf.varget("energy")     # (32, 5)
            theta_table = cdf.varget("theta")        # (32, 176, 3)
            phi_table = cdf.varget("phi")             # (32, 176, 3)
            domega_table = cdf.varget("domega")       # (32, 176, 3)
            gf_table = cdf.varget("gf")               # (32, 176, 3)
            geom_factor = float(cdf.varget("geom_factor"))

        del cdf

    if not all_times:
        raise ValueError(f"No ESD data in time range {trange}")

    times = np.concatenate(all_times)
    eflux = np.concatenate(all_eflux)
    bins_mask = np.concatenate(all_bins)
    en_ind = np.concatenate(all_en_ind)
    an_ind = np.concatenate(all_an_ind)
    phi_offset = np.concatenate(all_phi_offset)
    eff = np.concatenate(all_eff)
    integ_t = np.concatenate(all_integ_t)

    mode_en = int(np.median(en_ind))
    mode_an = int(np.median(an_ind))
    energy = energy_table[:, mode_en]
    theta = theta_table[:, :, mode_an]
    phi = phi_table[:, :, mode_an]
    domega = domega_table[:, :, mode_an]
    gf_mode = gf_table[:, :, mode_an]

    # one-count eflux per bin, inverts the thm_convert_esa_units eflux scale
    # eflux = counts / (integ_t * geom_factor * gf * eff), so one count is 1/that
    # this is the perp channels noise floor downstream, decoupled from signal amplitude
    gf_eff = geom_factor * gf_mode[None, :, :] * eff
    with np.errstate(divide="ignore", invalid="ignore"):
        onecount = np.where(gf_eff > 0, 1.0 / (integ_t[:, None, None] * gf_eff), np.nan)

    return ESDDistribution(
        times=times, eflux=eflux, energy=energy,
        theta=theta, phi=phi, domega=domega,
        bins_mask=bins_mask, phi_offset=phi_offset,
        en_ind=en_ind, an_ind=an_ind, onecount=onecount,
    )


def load_bfield_dsl(probe: str, trange: list[str], data_dir: str) -> tuple[np.ndarray, np.ndarray]:
    # returns (times, b_dsl) where b_dsl is (n,3)
    import os
    os.environ["THM_DATA_DIR"] = data_dir
    themis.fgm(probe=probe, trange=trange, level="l2", time_clip=True)
    varname = f"th{probe}_fgs_dsl"
    d = get_data(varname)
    if d is None:
        varname = f"th{probe}_fgl_dsl"
        d = get_data(varname)
    if d is None:
        varname = f"th{probe}_fgh_dsl"
        d = get_data(varname)
    if d is None:
        raise ValueError(f"No FGM DSL data found for probe {probe}")
    return d.times, d.y


def load_moments(probe: str, trange: list[str], data_dir: str) -> dict:
    import os
    os.environ["THM_DATA_DIR"] = data_dir
    themis.esa(probe=probe, trange=trange, level="l2", time_clip=True)

    prefix = f"th{probe}_peif"
    result = {}

    d = get_data(f"{prefix}_density")
    if d is not None:
        result["density_times"] = d.times
        result["density"] = d.y

    d = get_data(f"{prefix}_velocity_gsm")
    if d is None:
        d = get_data(f"{prefix}_velocity_dsl")
    if d is not None:
        result["velocity_times"] = d.times
        result["velocity"] = d.y

    d = get_data(f"{prefix}_avgtemp")
    if d is not None:
        result["temp_times"] = d.times
        result["temperature"] = d.y

    d = get_data(f"{prefix}_vthermal")
    if d is not None:
        result["vthermal"] = d.y

    return result


@dataclass
class PitchAngleSpectra:
    # pa-gated energy spectra over a time interval
    times: np.ndarray          # (ntime,)
    energy: np.ndarray         # (nenergy,) energy bin centers in ev
    omni: np.ndarray           # (ntime, nenergy) omnidirectional flux
    para: np.ndarray           # (ntime, nenergy) field-aligned 0-30 deg
    anti: np.ndarray           # (ntime, nenergy) anti-field-aligned 150-180 deg
    perp: np.ndarray           # (ntime, nenergy) perpendicular 75-105 deg
    perp_floor: np.ndarray     # (ntime, nenergy) one-count noise floor of perp cone avg
    pa_coverage_para: np.ndarray  # (ntime,) frac of solid angle in para cone
    pa_coverage_anti: np.ndarray  # (ntime,) frac of solid angle in anti cone
    pa_coverage_perp: np.ndarray  # (ntime,) frac of solid angle in perp cone


def _compute_pitch_angles(theta_inst: np.ndarray, phi_inst: np.ndarray,
                          b_dsl: np.ndarray, phi_offset: float) -> np.ndarray:
    # pa btwn particle vel and b, vel is opposite to look dir
    theta_rad = np.deg2rad(theta_inst)
    # phi_offset corrects spin phase from inst frame to dsl
    phi_rad = np.deg2rad(phi_inst + phi_offset)

    # sphere to cartesian, theta is elevation so cos goes on x/y
    # need vector form to dot with b-field below
    look_x = np.cos(theta_rad) * np.cos(phi_rad)
    look_y = np.cos(theta_rad) * np.sin(phi_rad)
    look_z = np.sin(theta_rad)

    # detector sees particles coming in, so v = -look
    vx, vy, vz = -look_x, -look_y, -look_z

    b_mag = np.linalg.norm(b_dsl)
    # skip when b is too small, pa is undefined
    if b_mag < 0.1:
        return np.full_like(theta_inst, np.nan)

    bhat = b_dsl / b_mag
    # cos(pa) = vhat . bhat, both are unit vectors
    # dot product of unit vectors gives cos of angle between them
    cos_pa = vx * bhat[0] + vy * bhat[1] + vz * bhat[2]
    # clip floating point overshoot before arccos
    cos_pa = np.clip(cos_pa, -1.0, 1.0)
    return np.rad2deg(np.arccos(cos_pa))


def compute_pa_spectra(
    dist: ESDDistribution,
    b_times: np.ndarray,
    b_dsl: np.ndarray,
    para_range: tuple[float, float] = (0.0, 30.0),
    anti_range: tuple[float, float] = (150.0, 180.0),
    perp_range: tuple[float, float] = (75.0, 105.0),
) -> PitchAngleSpectra:
    # reduces 3d dist to pa-gated energy spectra per timestep
    b_interp = interp1d(b_times, b_dsl, axis=0, kind="linear",
                        bounds_error=False, fill_value="extrapolate")

    ntime = len(dist.times)
    nenergy = len(dist.energy)
    valid_energy_mask = dist.energy > 0

    omni = np.full((ntime, nenergy), np.nan)
    para = np.full((ntime, nenergy), np.nan)
    anti = np.full((ntime, nenergy), np.nan)
    perp = np.full((ntime, nenergy), np.nan)
    perp_floor = np.full((ntime, nenergy), np.nan)
    cov_para = np.zeros(ntime)
    cov_anti = np.zeros(ntime)
    cov_perp = np.zeros(ntime)

    for t in range(ntime):
        b_vec = b_interp(dist.times[t])
        phi_off = dist.phi_offset[t]
        bins_valid = dist.bins_mask[t]
        flux = dist.eflux[t]

        for e in range(nenergy):
            if not valid_energy_mask[e]:
                continue

            valid = bins_valid[e].astype(bool)
            if not np.any(valid):
                continue

            f = flux[e, valid]
            dw = dist.domega[e, valid]
            oc = dist.onecount[t, e, valid]

            th_bins = dist.theta[e, valid]
            ph_bins = dist.phi[e, valid]

            pa = _compute_pitch_angles(th_bins, ph_bins, b_vec, phi_off)

            nan_mask = np.isnan(pa)
            if np.all(nan_mask):
                continue

            good = ~nan_mask & np.isfinite(f) & (f > 0)
            if not np.any(good):
                continue

            # solid-angle weighted avg, bigger bins count more
            # equal-weight avg would bias toward small high-latitude bins
            total_weight = dw[good].sum()
            if total_weight > 0:
                omni[t, e] = np.sum(f[good] * dw[good]) / total_weight

            # gate by pa range then weighted avg over the cone only
            # isolates the directional pop we care about
            in_para = good & (pa >= para_range[0]) & (pa <= para_range[1])
            if np.any(in_para):
                w_para = dw[in_para].sum()
                para[t, e] = np.sum(f[in_para] * dw[in_para]) / w_para
                # coverage = frac of solid angle inside the cone
                # low coverage means the avg is unreliable, drop it downstream
                cov_para[t] = max(cov_para[t], w_para / total_weight)

            in_anti = good & (pa >= anti_range[0]) & (pa <= anti_range[1])
            if np.any(in_anti):
                w_anti = dw[in_anti].sum()
                anti[t, e] = np.sum(f[in_anti] * dw[in_anti]) / w_anti
                cov_anti[t] = max(cov_anti[t], w_anti / total_weight)

            # perp cone is the clean background for R, beam depletes it not enhances
            in_perp = good & (pa >= perp_range[0]) & (pa <= perp_range[1])
            if np.any(in_perp):
                w_perp = dw[in_perp].sum()
                perp[t, e] = np.sum(f[in_perp] * dw[in_perp]) / w_perp
                # one-count noise floor through the same solid-angle weighting
                perp_floor[t, e] = np.sum(oc[in_perp] * dw[in_perp]) / w_perp
                cov_perp[t] = max(cov_perp[t], w_perp / total_weight)

    return PitchAngleSpectra(
        times=dist.times,
        energy=dist.energy,
        omni=omni,
        para=para,
        anti=anti,
        perp=perp,
        perp_floor=perp_floor,
        pa_coverage_para=cov_para,
        pa_coverage_anti=cov_anti,
        pa_coverage_perp=cov_perp,
    )


@dataclass
class FeatureTable:
    # per-timestep features for beam classification
    times: np.ndarray
    e_peak: np.ndarray           # energy of max flux, ev
    e_beam: np.ndarray           # energy where max |asym| was found, ev
    width: np.ndarray            # normalized spectral width
    asymmetry: np.ndarray        # signed asym at e_beam, picked from max |asym| across bins
    para_to_omni: np.ndarray     # parallel cone / omni ratio at e_beam
    energy_ratio: np.ndarray     # e_flow / e_th from moments
    peak_prom: np.ndarray        # log10 prominence of narrow line inside the coherent run
    peak_width: np.ndarray       # fwhm of that line in bins
    e_line: np.ndarray           # energy of the line, ev (E_beam)
    de_line: np.ndarray          # line fwhm in ev (delta energy)
    eb_over_de: np.ndarray       # e_line / line fwhm, beam monochromaticity
    r_beam: np.ndarray           # flux-weighted mean R (dominant cone / perp) over the run
    pa_max_ratio: np.ndarray     # max para / max anti over run bins, dominant/sub, flux transfer
    coherent_ok: np.ndarray      # bool, a real coherent directional run was found
    perp_depleted: np.ndarray    # bool, run had sub-floor perp flux so R denom was clamped
    pa_ok_both: np.ndarray       # bool, both cones sampled, asym is trustworthy
    pa_ok_para: np.ndarray       # bool, para cone sampled, p2o is trustworthy


PROTON_MASS_KG = 1.6726219e-27
EV_PER_JOULE = 6.242e18
BOLTZMANN_EV = 8.617e-5


def extract_features(
    spectra: PitchAngleSpectra,
    moments: dict,
    energy_cutoff_low: float = 30.0,
    pa_coverage_threshold: float = 0.01,
    beam_flux_floor: float = 0.1,
    coherent_asym_min: float = 0.2,
    coherent_dir_min: float = 1.2,
    coherent_min_bins: int = 2,
    peak_width_max: float = 4.0,
    peak_wlen: int = 5,
) -> FeatureTable:
    ntime = len(spectra.times)
    energy = spectra.energy

    valid_e = energy >= energy_cutoff_low
    e_valid = energy[valid_e]

    e_peak = np.full(ntime, np.nan)
    e_beam = np.full(ntime, np.nan)
    width = np.full(ntime, np.nan)
    asymmetry = np.full(ntime, np.nan)
    para_to_omni = np.full(ntime, np.nan)
    energy_ratio = np.full(ntime, np.nan)
    peak_prom = np.full(ntime, np.nan)
    peak_width = np.full(ntime, np.nan)
    e_line = np.full(ntime, np.nan)
    de_line = np.full(ntime, np.nan)
    eb_over_de = np.full(ntime, np.nan)
    r_beam = np.full(ntime, np.nan)
    pa_max_ratio = np.full(ntime, np.nan)
    coherent_ok = np.zeros(ntime, dtype=bool)
    perp_depleted = np.zeros(ntime, dtype=bool)
    pa_ok_both = np.zeros(ntime, dtype=bool)
    pa_ok_para = np.zeros(ntime, dtype=bool)

    if "velocity" in moments and "temperature" in moments:
        vel_interp = interp1d(moments["velocity_times"], moments["velocity"],
                              axis=0, kind="nearest",
                              bounds_error=False, fill_value=np.nan)
        temp_interp = interp1d(moments["temp_times"], moments["temperature"],
                               axis=0, kind="nearest",
                               bounds_error=False, fill_value=np.nan)
    else:
        vel_interp = temp_interp = None

    for t in range(ntime):
        omni_t = spectra.omni[t, valid_e]
        para_t = spectra.para[t, valid_e]
        anti_t = spectra.anti[t, valid_e]
        perp_t = spectra.perp[t, valid_e]
        perp_floor_t = spectra.perp_floor[t, valid_e]

        # per-feature coverage: asym needs both cones, p2o only needs para
        para_ok = spectra.pa_coverage_para[t] >= pa_coverage_threshold
        anti_ok = spectra.pa_coverage_anti[t] >= pa_coverage_threshold
        pa_ok_para[t] = para_ok
        pa_ok_both[t] = para_ok and anti_ok

        if np.all(np.isnan(omni_t)):
            continue
        omni_finite = np.where(np.isfinite(omni_t), omni_t, 0.0)
        # peak energy = argmax of omni spectrum
        # characterizes where the bulk population sits in energy
        idx_peak = np.argmax(omni_finite)
        e_peak[t] = e_valid[idx_peak]

        # width = std/mean of flux distribution over energy, narrow beam -> small
        total_flux = np.nansum(omni_finite)
        if total_flux > 0 and e_peak[t] > 0:
            # flux-weighted mean energy, like center of mass
            # high-flux bins dominate so it tracks the populated part of the spectrum
            e_mean = np.nansum(omni_finite * e_valid) / total_flux
            # flux-weighted variance
            # 2nd moment of the flux dist, measures energy spread
            e_var = np.nansum(omni_finite * (e_valid - e_mean) ** 2) / total_flux
            # normalize by e_peak so width is dimensionless
            width[t] = np.sqrt(e_var) / e_peak[t]

        # scan all bins for max |asym|, catches beams that ride on a plasma sheet
        # peak-window approach misses beams when e_peak is the plasma sheet not the beam
        peak_omni = omni_finite[idx_peak]
        # flux floor mask, low-flux bins are noisy and produce spurious asym near +/-1
        flux_mask = omni_finite >= beam_flux_floor * peak_omni

        # per-bin signed asym
        denom_asym = para_t + anti_t
        with np.errstate(invalid="ignore", divide="ignore"):
            asym_bins = np.where(denom_asym > 0,
                                 (para_t - anti_t) / denom_asym, np.nan)
        # per-bin p2o, feeds the para_to_omni feature only
        with np.errstate(invalid="ignore", divide="ignore"):
            p2o_bins = np.where(omni_finite > 0,
                                para_t / omni_finite, np.nan)

        # R = dominant cone / perp cone, the coherent_dir gate ratio
        # omni includes the beam cone so it partially cancels the signal, perp is a
        # clean background since a field-aligned beam depletes the perpendicular cone
        # clamp perp at its own one-count noise floor, not the omni peak, so a depleted
        # perp cant blow R up and perp_depleted stays comparable across intervals
        perp_subfloor = np.isfinite(perp_t) & (perp_t < perp_floor_t)
        with np.errstate(invalid="ignore", divide="ignore"):
            perp_eff = np.where(np.isfinite(perp_t),
                                np.maximum(perp_t, perp_floor_t), np.nan)
            r_para = np.where(perp_eff > 0, para_t / perp_eff, np.nan)
            r_anti = np.where(perp_eff > 0, anti_t / perp_eff, np.nan)

        # dominant cone / perp, direction-agnostic enhancement
        # high value means the dominant cone is brighter than the perp background
        with np.errstate(invalid="ignore"):
            dir_enhanced = np.where(
                np.isfinite(asym_bins) & (asym_bins >= 0), r_para,
                np.where(np.isfinite(asym_bins) & (asym_bins < 0), r_anti, np.nan)
            )

        # a bin qualifies if flux-floor cleared, asym magnitude clears gate,
        # and the dominant cone is enhanced over the perp background
        with np.errstate(invalid="ignore"):
            qual = (flux_mask &
                    np.isfinite(asym_bins) &
                    np.isfinite(dir_enhanced) &
                    (np.abs(asym_bins) >= coherent_asym_min) &
                    (dir_enhanced >= coherent_dir_min))

        # find longest run of same-sign qualifying bins, allow gap of 1 non-qual bin
        # gap tolerance catches beams whose qualifying bins are interrupted by
        # a single nan or a bin that just barely misses the dir/asym gate
        n = len(asym_bins)
        signs = np.where(asym_bins > 0, 1,
                         np.where(asym_bins < 0, -1, 0)).astype(int)
        best_idxs = []   # qualifying bin indices in the best run
        cur_idxs = []
        cur_sign = 0
        gap = 0
        max_gap = 1
        for k in range(n):
            if qual[k] and signs[k] != 0:
                if cur_sign == 0 or signs[k] == cur_sign:
                    cur_idxs.append(k)
                    cur_sign = signs[k]
                    gap = 0
                else:
                    # opposite sign, close current run and start fresh
                    if len(cur_idxs) > len(best_idxs):
                        best_idxs = cur_idxs
                    cur_idxs = [k]
                    cur_sign = signs[k]
                    gap = 0
            else:
                # non-qualifying bin, extend gap if a run is active
                if cur_idxs:
                    gap += 1
                    if gap > max_gap:
                        if len(cur_idxs) > len(best_idxs):
                            best_idxs = cur_idxs
                        cur_idxs = []
                        cur_sign = 0
                        gap = 0
        if len(cur_idxs) > len(best_idxs):
            best_idxs = cur_idxs

        if pa_ok_both[t] and len(best_idxs) >= coherent_min_bins:
            coherent_ok[t] = True
            # flux-weighted avg over only the qualifying bins in the run
            # gaps are excluded so their non-beam contribution doesnt pollute the avg
            idx_arr = np.array(best_idxs)
            # flag if any run bin leaned on a clamped (depleted) perp denominator
            perp_depleted[t] = bool(np.any(perp_subfloor[idx_arr]))
            w = omni_finite[idx_arr]
            w_sum = w.sum()
            asymmetry[t] = np.sum(asym_bins[idx_arr] * w) / w_sum
            e_beam[t] = np.sum(e_valid[idx_arr] * w) / w_sum
            # flux-weighted mean R over the run, per-beam directional enhancement
            r_beam[t] = np.sum(dir_enhanced[idx_arr] * w) / w_sum
            # max para vs max anti within the beam band, dominant/sub = flux transfer
            run_para = para_t[idx_arr]
            run_anti = anti_t[idx_arr]
            mp = np.nanmax(run_para) if np.any(np.isfinite(run_para)) else np.nan
            ma = np.nanmax(run_anti) if np.any(np.isfinite(run_anti)) else np.nan
            if np.isfinite(mp) and np.isfinite(ma) and mp > 0 and ma > 0:
                pa_max_ratio[t] = max(mp, ma) / min(mp, ma)
            if pa_ok_para[t]:
                p2o_vals = p2o_bins[idx_arr]
                p2o_ok = np.isfinite(p2o_vals)
                if np.any(p2o_ok):
                    w_p = w[p2o_ok]
                    para_to_omni[t] = np.sum(p2o_vals[p2o_ok] * w_p) / w_p.sum()
        else:
            # fallback to peak +/-2 window when no coherent run found
            lo = max(0, idx_peak - 2)
            hi = min(len(e_valid), idx_peak + 3)
            f_para = np.nanmean(para_t[lo:hi])
            f_anti = np.nanmean(anti_t[lo:hi])
            f_omni = np.nanmean(omni_finite[lo:hi])
            if pa_ok_both[t] and np.isfinite(f_para) and np.isfinite(f_anti) and (f_para + f_anti) > 0:
                asymmetry[t] = (f_para - f_anti) / (f_para + f_anti)
                e_beam[t] = e_peak[t]
            if pa_ok_para[t] and np.isfinite(f_para) and f_omni > 0:
                para_to_omni[t] = f_para / f_omni

        # spectral line detection, local to the coherent run not global
        # a beam = the directional region is ALSO a narrow line, so we only
        # look in the dominant cone within the run band, a prominent line
        # elsewhere (eg the anti plasma sheet peak) is not the beam
        if coherent_ok[t]:
            # dominant cone is the one the run points at
            cone = para_t if asymmetry[t] >= 0 else anti_t
            fin = np.isfinite(cone) & (cone > 0)
            if fin.sum() >= 3:
                # compress to finite bins, log flux so prominence = brightness ratio
                idx_map = np.where(fin)[0]
                logf = np.log10(cone[fin])
                # wlen bounds the prominence window so a narrow line on a broad
                # pedestal gets a local prominence and width, not inflated by far valleys
                pks, props = find_peaks(logf, prominence=0.05, wlen=peak_wlen,
                                        width=(None, peak_width_max))
                # run band in valid_e index space, allow 1 bin of slop
                lo_b = min(best_idxs) - 1
                hi_b = max(best_idxs) + 1
                # energies on the compressed axis, for converting fwhm bins to ev
                e_comp = e_valid[idx_map]
                comp_x = np.arange(len(e_comp))
                for j, pk in enumerate(pks):
                    b = idx_map[pk]
                    # line must sit inside the directional run to count
                    if b < lo_b or b > hi_b:
                        continue
                    prom = props["prominences"][j]
                    if not np.isfinite(peak_prom[t]) or prom > peak_prom[t]:
                        peak_prom[t] = prom
                        peak_width[t] = props["widths"][j]
                        e_line[t] = e_valid[b]
                        # de = line fwhm in ev from the half-max crossings,
                        # eb/de is how monoenergetic the beam line is
                        e_left = np.interp(props["left_ips"][j], comp_x, e_comp)
                        e_right = np.interp(props["right_ips"][j], comp_x, e_comp)
                        de = abs(e_right - e_left)
                        de_line[t] = de
                        eb_over_de[t] = e_valid[b] / de if de > 0 else np.nan

        # e_flow/e_th, ratio of bulk kinetic to thermal energy, beams have higher ratio
        if vel_interp is not None and temp_interp is not None:
            v = vel_interp(spectra.times[t])
            T = temp_interp(spectra.times[t])
            if np.all(np.isfinite(v)) and np.isfinite(T) and T > 0:
                v_mag = np.linalg.norm(v)
                # 1/2 m v^2 in joules, then convert to ev, v is km/s so *1e3
                # puts e_flow in same units as T so the ratio is meaningful
                e_flow = 0.5 * PROTON_MASS_KG * (v_mag * 1e3) ** 2 * EV_PER_JOULE
                energy_ratio[t] = e_flow / T

    return FeatureTable(
        times=spectra.times,
        e_peak=e_peak,
        e_beam=e_beam,
        width=width,
        asymmetry=asymmetry,
        para_to_omni=para_to_omni,
        energy_ratio=energy_ratio,
        peak_prom=peak_prom,
        peak_width=peak_width,
        e_line=e_line,
        de_line=de_line,
        eb_over_de=eb_over_de,
        r_beam=r_beam,
        pa_max_ratio=pa_max_ratio,
        coherent_ok=coherent_ok,
        perp_depleted=perp_depleted,
        pa_ok_both=pa_ok_both,
        pa_ok_para=pa_ok_para,
    )


@dataclass
class ClassifierParams:
    # spectral features are primary, moments are soft bc they mix beam + plasma sheet
    asymmetry_min: float = 0.2
    width_max: float = 0.8
    para_to_omni_min: float = 1.3
    energy_ratio_min: float = 0.5
    score_threshold: float = 0.4
    min_coverage: float = 0.01
    beam_flux_floor: float = 0.1
    # per-bin gates for coherent-region detection
    coherent_asym_min: float = 0.2   # min |asym| per bin to qualify
    coherent_dir_min: float = 1.2    # min dominant cone / omni per bin
    coherent_min_bins: int = 2       # min adjacent bins for a run
    # spectral line gates, prominence in log10 so 0.3 = 2x above local baseline
    peak_prom_min: float = 0.3
    peak_width_max: float = 4.0      # fwhm cap in bins, rejects broad bumps
    peak_wlen: int = 5               # local window for prominence, bounds pedestal inflation
    # weights, spectral > moments
    w_asymmetry: float = 0.35
    w_width: float = 0.25
    w_para_to_omni: float = 0.25
    # peak prom takes the slot energy ratio vacated, er stays dead
    w_peak_prom: float = 0.15
    w_energy_ratio: float = 0.0


@dataclass
class ClassificationResult:
    times: np.ndarray
    is_beam: np.ndarray          # bool
    beam_score: np.ndarray       # 0-1 continuous score
    beam_direction: np.ndarray   # +1 parallel, -1 anti-parallel, 0 unknown


def classify_beams(features: FeatureTable,
                   params: ClassifierParams | None = None) -> ClassificationResult:
    if params is None:
        params = ClassifierParams()

    ntime = len(features.times)
    is_beam = np.zeros(ntime, dtype=bool)
    beam_score = np.zeros(ntime)
    beam_direction = np.zeros(ntime, dtype=int)

    for t in range(ntime):
        er = features.energy_ratio[t]
        asym = features.asymmetry[t]
        w = features.width[t]
        p2o = features.para_to_omni[t]

        # no full-timestep gate, nan-features score 0 and naturally drop out

        # component scores in [0,1], 1 means feature fully satisfies its threshold
        # ramp: 0 at asym=0, 0.5 at threshold, capped at 1 when asym = 2*threshold
        # smooth ramp lets borderline features still contribute to total score
        s_asym = 0.0
        if np.isfinite(asym):
            s_asym = np.clip(abs(asym) / params.asymmetry_min, 0, 2) / 2

        # ramp: 1 at width=0, 0 at width_max, 0 beyond
        # inverted cuz beams have narrow spectra, lower width = stronger signal
        s_width = 0.0
        if np.isfinite(w):
            s_width = np.clip((params.width_max - w) / params.width_max, 0, 1)

        # ramp shifted so p2o=1 (no enhancement) gives 0, threshold gives 0.5
        # baseline is 1 not 0 cuz p2o<1 just means para cone is dimmer than avg, not anti-beam
        s_p2o = 0.0
        if np.isfinite(p2o):
            s_p2o = np.clip((p2o - 1.0) / (params.para_to_omni_min - 1.0), 0, 2) / 2

        # dead code, will remove/edit later
        s_er = 0.0
        if np.isfinite(er):
            s_er = np.clip(er / params.energy_ratio_min, 0, 2) / 2

        # ramp on log prominence, 0.5 at threshold, width already gated in find_peaks
        s_peak = 0.0
        prom = features.peak_prom[t]
        if np.isfinite(prom):
            s_peak = np.clip(prom / params.peak_prom_min, 0, 2) / 2

        # weighted sum, weights total to 1 so score stays in [0,1]
        beam_score[t] = (params.w_asymmetry * s_asym +
                         params.w_width * s_width +
                         params.w_para_to_omni * s_p2o +
                         params.w_peak_prom * s_peak +
                         params.w_energy_ratio * s_er)

        # score-based, catches weak beams spread across features
        score_ok = beam_score[t] >= params.score_threshold

        # hard rule fallback, catches strong beams even if score is borderline
        asym_ok = np.isfinite(asym) and abs(asym) >= params.asymmetry_min
        width_ok = np.isfinite(w) and w <= params.width_max
        p2o_ok = np.isfinite(p2o) and p2o >= params.para_to_omni_min
        hard_ok = asym_ok and (width_ok or p2o_ok)

        # and-gate, directional run must also be a narrow line at same energy
        # both detectors agree or its not a beam, kills noise that fires one alone
        gate = (features.coherent_ok[t] and np.isfinite(prom) and
                prom >= params.peak_prom_min)

        is_beam[t] = (score_ok or hard_ok) and gate

        if np.isfinite(asym):
            if asym > 0:
                beam_direction[t] = 1
            elif asym < 0:
                beam_direction[t] = -1

    return ClassificationResult(
        times=features.times,
        is_beam=is_beam,
        beam_score=beam_score,
        beam_direction=beam_direction,
    )


def smooth_labels(result: ClassificationResult,
                  min_consecutive: int = 1) -> ClassificationResult:
    # requires min_consecutive flagged steps to keep a beam interval
    # default 1 = keep isolated beams, and-gate already handles precision
    smoothed = np.zeros_like(result.is_beam)
    n = len(smoothed)

    run_start = None
    for i in range(n + 1):
        if i < n and result.is_beam[i]:
            if run_start is None:
                run_start = i
        else:
            if run_start is not None:
                run_len = i - run_start
                if run_len >= min_consecutive:
                    smoothed[run_start:i] = True
                run_start = None

    return ClassificationResult(
        times=result.times,
        is_beam=smoothed,
        beam_score=result.beam_score,
        beam_direction=result.beam_direction,
    )


def diagnose_window(
    spectra: PitchAngleSpectra,
    features: FeatureTable,
    classification: ClassificationResult,
    params: ClassifierParams,
    ut_start: str,
    ut_end: str,
    out_path: str | None = None,
) -> None:
    # dump per-bin spectra and features for timesteps in [ut_start, ut_end]
    # writes to file if out_path given, prints summary to stdout either way
    from datetime import datetime, timezone

    def parse_ut(s: str, ref_unix: float) -> float:
        ref_dt = datetime.fromtimestamp(ref_unix, tz=timezone.utc)
        parts = s.split(":")
        h = int(parts[0])
        m = int(parts[1]) if len(parts) > 1 else 0
        sec = float(parts[2]) if len(parts) > 2 else 0.0
        dt = ref_dt.replace(hour=h, minute=m, second=int(sec), microsecond=0)
        return dt.timestamp()

    t0 = parse_ut(ut_start, spectra.times[0])
    t1 = parse_ut(ut_end, spectra.times[0])
    mask = (spectra.times >= t0) & (spectra.times <= t1)
    idxs = np.where(mask)[0]

    if len(idxs) == 0:
        print(f"[diag] no timesteps in {ut_start}-{ut_end}")
        return

    lines = []
    lines.append(f"=== Diagnostic: {ut_start}-{ut_end} UT, {len(idxs)} timesteps ===")
    lines.append(f"params: asym_min={params.asymmetry_min} width_max={params.width_max} "
                 f"p2o_min={params.para_to_omni_min} score_thr={params.score_threshold} "
                 f"flux_floor={params.beam_flux_floor}")

    for t in idxs:
        dt = datetime.fromtimestamp(spectra.times[t], tz=timezone.utc)
        lines.append(f"\n--- {dt.strftime('%H:%M:%S')} (idx={t}) ---")
        lines.append(f"  e_peak={features.e_peak[t]:.1f}eV e_beam={features.e_beam[t]:.1f}eV "
                     f"width={features.width[t]:.3f}")
        lines.append(f"  asym={features.asymmetry[t]:.3f} p2o={features.para_to_omni[t]:.3f} "
                     f"e_ratio={features.energy_ratio[t]:.3f}")
        lines.append(f"  peak_prom={features.peak_prom[t]:.3f} "
                     f"peak_width={features.peak_width[t]:.2f}bins "
                     f"e_line={features.e_line[t]:.1f}eV "
                     f"coherent_ok={features.coherent_ok[t]}")
        lines.append(f"  eb_over_de={features.eb_over_de[t]:.2f} "
                     f"de_line={features.de_line[t]:.1f}eV "
                     f"r_beam={features.r_beam[t]:.3f} "
                     f"pa_max_ratio={features.pa_max_ratio[t]:.3f} "
                     f"perp_depleted={features.perp_depleted[t]}")
        lines.append(f"  pa_cov: para={spectra.pa_coverage_para[t]:.3f} "
                     f"anti={spectra.pa_coverage_anti[t]:.3f} "
                     f"perp={spectra.pa_coverage_perp[t]:.3f} "
                     f"(both_ok={features.pa_ok_both[t]} para_ok={features.pa_ok_para[t]})")
        lines.append(f"  score={classification.beam_score[t]:.3f} "
                     f"is_beam={classification.is_beam[t]} dir={classification.beam_direction[t]}")

        e = spectra.energy
        omni = spectra.omni[t]
        para = spectra.para[t]
        anti = spectra.anti[t]
        omni_finite = np.where(np.isfinite(omni), omni, 0.0)
        peak = np.max(omni_finite) if np.any(omni_finite > 0) else 0.0
        floor = params.beam_flux_floor * peak

        perp = spectra.perp[t]
        pfloor = spectra.perp_floor[t]
        lines.append(f"  per-bin (omni >= {floor:.2e}):")
        lines.append(f"    {'E[eV]':>9} {'omni':>10} {'para':>10} {'anti':>10} "
                     f"{'perp':>10} {'asym':>7} {'p2o':>6} {'R':>6}")
        for b in range(len(e)):
            if not np.isfinite(omni[b]) or omni[b] < floor:
                continue
            denom = (para[b] if np.isfinite(para[b]) else 0) + \
                    (anti[b] if np.isfinite(anti[b]) else 0)
            ab = (para[b] - anti[b]) / denom if denom > 0 else np.nan
            pb = para[b] / omni[b] if omni[b] > 0 and np.isfinite(para[b]) else np.nan
            # R = dominant cone / clamped perp, mirrors the coherent_dir gate
            pe = perp[b]
            pf = pfloor[b]
            if np.isfinite(pe):
                perp_eff = max(pe, pf) if np.isfinite(pf) else pe
                dom = para[b] if (np.isfinite(ab) and ab >= 0) else anti[b]
                rb = dom / perp_eff if perp_eff > 0 and np.isfinite(dom) else np.nan
            else:
                rb = np.nan
            lines.append(f"    {e[b]:>9.1f} {omni[b]:>10.2e} "
                         f"{para[b] if np.isfinite(para[b]) else float('nan'):>10.2e} "
                         f"{anti[b] if np.isfinite(anti[b]) else float('nan'):>10.2e} "
                         f"{pe if np.isfinite(pe) else float('nan'):>10.2e} "
                         f"{ab:>7.3f} {pb:>6.3f} {rb:>6.3f}")

    text = "\n".join(lines)
    if out_path:
        with open(out_path, "w") as f:
            f.write(text + "\n")
        # short stdout summary
        n_flagged = sum(int(classification.is_beam[t]) for t in idxs)
        max_score = max((classification.beam_score[t] for t in idxs), default=0.0)
        max_asym = max((abs(features.asymmetry[t]) for t in idxs
                        if np.isfinite(features.asymmetry[t])), default=0.0)
        print(f"[diag] wrote {out_path}")
        print(f"[diag] window summary: {len(idxs)} steps, {n_flagged} flagged, "
              f"max_score={max_score:.3f} max|asym|={max_asym:.3f}")
    else:
        print(text)


def threshold_sensitivity(features: FeatureTable,
                          param_ranges: dict | None = None) -> dict:
    # varies thresholds, reports label stability
    if param_ranges is None:
        param_ranges = {
            "asymmetry_min": [0.1, 0.15, 0.2, 0.3, 0.4],
            "width_max": [0.5, 0.6, 0.8, 1.0, 1.2],
            "para_to_omni_min": [1.1, 1.2, 1.3, 1.5, 1.8],
        }

    results = {}
    base = ClassifierParams()
    base_result = classify_beams(features, base)
    base_count = base_result.is_beam.sum()

    for param_name, values in param_ranges.items():
        counts = []
        for val in values:
            kwargs = {
                "asymmetry_min": base.asymmetry_min,
                "width_max": base.width_max,
                "para_to_omni_min": base.para_to_omni_min,
                "energy_ratio_min": base.energy_ratio_min,
            }
            kwargs[param_name] = val
            p = ClassifierParams(**kwargs)
            r = classify_beams(features, p)
            counts.append(int(r.is_beam.sum()))
        results[param_name] = {"values": values, "beam_counts": counts,
                               "base_count": int(base_count)}

    return results


def plot_spectra_snapshot(
    spectra: PitchAngleSpectra,
    time_idx: int,
    label: str = "",
    ax=None,
):
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 5))

    e = spectra.energy
    valid = e > 0

    ax.loglog(e[valid], spectra.omni[time_idx, valid], "k-", lw=2, label="Omni")
    ax.loglog(e[valid], spectra.para[time_idx, valid], "r-", lw=1.5, label="0°–30° (para)")
    ax.loglog(e[valid], spectra.anti[time_idx, valid], "b-", lw=1.5, label="150°–180° (anti)")

    ax.set_xlabel("Energy [eV]")
    ax.set_ylabel("Energy Flux [eV/cm²-s-sr-eV]")
    ax.legend(fontsize=9)
    if label:
        ax.set_title(label)
    ax.grid(True, alpha=0.3)
    return ax


def _score_to_size(score, thr):
    # marker size grows with score above thr, borderline small, strong large
    frac = np.clip((score - thr) / max(1.0 - thr, 1e-6), 0.0, 1.0)
    return 12.0 + 70.0 * frac


def plot_feature_timeseries(
    features: FeatureTable,
    classification: ClassificationResult,
    spectra: PitchAngleSpectra,
    out_png: str,
    params: ClassifierParams,
    title: str = "",
):
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from datetime import datetime, timezone

    times_dt = [datetime.fromtimestamp(t, tz=timezone.utc) for t in features.times]

    fig, axes = plt.subplots(8, 1, figsize=(14, 21), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1, 1, 1, 1, 1, 0.5, 3]})

    # omni spectrogram
    ax = axes[0]
    e_valid = spectra.energy > 0
    e_plot = spectra.energy[e_valid]
    z = spectra.omni[:, e_valid].T
    z = np.where(z > 0, z, np.nan)
    pcm = ax.pcolormesh(times_dt, e_plot, z,
                        norm=plt.matplotlib.colors.LogNorm(vmin=1e3, vmax=1e8),
                        cmap="jet", shading="auto")
    ax.set_yscale("log")
    ax.set_ylabel("Energy [eV]")
    ax.set_ylim(5, 30000)
    # colorbar via divider so it doesnt shrink the spectrogram data area
    # keeps x-axis aligned with the panels below
    div = make_axes_locatable(ax)
    cax = div.append_axes("right", size="1.5%", pad=0.05)
    fig.colorbar(pcm, cax=cax, label="eflux")
    if title:
        ax.set_title(title)

    # bottom panel, same spectrogram with beam detections overlaid at e_peak
    axb = axes[7]
    pcmb = axb.pcolormesh(times_dt, e_plot, z,
                          norm=plt.matplotlib.colors.LogNorm(vmin=1e3, vmax=1e8),
                          cmap="jet", shading="auto")
    axb.set_yscale("log")
    axb.set_ylabel("Energy [eV]")
    axb.set_ylim(5, 30000)
    divb = make_axes_locatable(axb)
    caxb = divb.append_axes("right", size="1.5%", pad=0.05)
    fig.colorbar(pcmb, cax=caxb, label="eflux")
    # markers at e_peak so they sit on the flux peak band, red para blue anti
    is_b = classification.is_beam
    t_arr = np.array(times_dt, dtype=object)
    for mask, col, lab in ((is_b & (classification.beam_direction > 0), "red", "para beam"),
                           (is_b & (classification.beam_direction < 0), "blue", "anti beam")):
        if np.any(mask):
            axb.scatter(t_arr[mask], features.e_peak[mask],
                        s=_score_to_size(classification.beam_score[mask], params.score_threshold),
                        c=col, edgecolors="white", linewidths=0.5, label=lab, zorder=3)
    if np.any(is_b):
        axb.legend(loc="upper right", fontsize=8, title="size ∝ beam score")

    # invisible spacer on the middle panels so they share the spectrogram width
    for other_ax in axes[1:7]:
        d = make_axes_locatable(other_ax)
        spacer = d.append_axes("right", size="1.5%", pad=0.05)
        spacer.axis("off")

    axes[1].semilogy(times_dt, features.e_peak, "k.", ms=2)
    axes[1].set_ylabel("E_peak [eV]")
    axes[1].set_ylim(10, 30000)

    axes[2].plot(times_dt, features.width, "k.", ms=2)
    axes[2].set_ylabel("Width")
    axes[2].axhline(params.width_max, color="r", ls="--", alpha=0.5, label="threshold")
    axes[2].legend(fontsize=8)

    axes[3].plot(times_dt, features.asymmetry, "k.", ms=2)
    axes[3].set_ylabel("Asymmetry")
    axes[3].axhline(params.asymmetry_min, color="r", ls="--", alpha=0.5)
    axes[3].axhline(-params.asymmetry_min, color="r", ls="--", alpha=0.5)
    axes[3].set_ylim(-1.1, 1.1)

    # peak prominence, the and-gate line, key precision filter
    axes[4].plot(times_dt, features.peak_prom, "k.", ms=2)
    axes[4].set_ylabel("Peak Prom")
    axes[4].axhline(params.peak_prom_min, color="r", ls="--", alpha=0.5, label="threshold")
    axes[4].legend(fontsize=8)

    axes[5].plot(times_dt, classification.beam_score, "k.", ms=2)
    axes[5].set_ylabel("Beam Score")
    axes[5].axhline(params.score_threshold, color="r", ls="--", alpha=0.5, label="threshold")
    axes[5].legend(fontsize=8)
    axes[5].set_ylim(0, 1.1)

    # classification color bar
    ax = axes[6]
    colors = []
    for i in range(len(features.times)):
        if classification.is_beam[i]:
            if classification.beam_direction[i] > 0:
                colors.append("red")
            elif classification.beam_direction[i] < 0:
                colors.append("blue")
            else:
                colors.append("orange")
        else:
            colors.append("lightgray")

    for i, (t, c) in enumerate(zip(times_dt, colors)):
        ax.axvspan(t, times_dt[min(i + 1, len(times_dt) - 1)],
                   color=c, alpha=0.8)
    ax.set_ylabel("Beam")
    ax.set_yticks([])

    # hourly major ticks, set on all panels via sharex
    axes[-1].xaxis.set_major_locator(mdates.HourLocator(interval=1))
    axes[-1].xaxis.set_minor_locator(mdates.MinuteLocator(byminute=[15, 30, 45]))
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axes[-1].set_xlabel("UT")
    # bolder hourly ticks across all panels
    for ax in axes:
        ax.tick_params(axis="x", which="major", width=1.8, length=7, color="black")
        ax.tick_params(axis="x", which="minor", width=0.8, length=3, color="gray")
    for lbl in axes[-1].get_xticklabels(which="major"):
        lbl.set_fontweight("bold")
    plt.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] wrote {out_png}")


def plot_curated_snapshots(
    spectra: PitchAngleSpectra,
    classification: ClassificationResult,
    features: FeatureTable,
    out_png: str,
    n_beam: int = 3,
    n_ps: int = 3,
    include_borderline: bool = True,
):
    import matplotlib.pyplot as plt
    from datetime import datetime, timezone

    beam_idx = np.where(classification.is_beam)[0]
    ps_idx = np.where(~classification.is_beam & features.pa_ok_both)[0]

    selected = []
    labels = []

    if len(beam_idx) > 0:
        step = max(1, len(beam_idx) // n_beam)
        for i in beam_idx[::step][:n_beam]:
            selected.append(i)
            t = datetime.fromtimestamp(spectra.times[i], tz=timezone.utc)
            labels.append(f"BEAM {t.strftime('%H:%M:%S')}")

    if len(ps_idx) > 0:
        step = max(1, len(ps_idx) // n_ps)
        for i in ps_idx[::step][:n_ps]:
            selected.append(i)
            t = datetime.fromtimestamp(spectra.times[i], tz=timezone.utc)
            labels.append(f"PS {t.strftime('%H:%M:%S')}")

    if include_borderline:
        borderline = np.where(
            (classification.beam_score > 0.3) &
            (classification.beam_score < 0.7) &
            features.pa_ok_both
        )[0]
        if len(borderline) > 0:
            step = max(1, len(borderline) // 2)
            for i in borderline[::step][:2]:
                if i not in selected:
                    selected.append(i)
                    t = datetime.fromtimestamp(spectra.times[i], tz=timezone.utc)
                    labels.append(f"BORDER {t.strftime('%H:%M:%S')}")

    if not selected:
        print("[WARN] No timesteps selected for snapshot plot")
        return

    ncols = min(3, len(selected))
    nrows = (len(selected) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
    if len(selected) == 1:
        axes = np.array([axes])
    axes = np.atleast_2d(axes)

    for idx, (sel, lbl) in enumerate(zip(selected, labels)):
        r, c = divmod(idx, ncols)
        plot_spectra_snapshot(spectra, sel, label=lbl, ax=axes[r, c])

    for idx in range(len(selected), nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r, c].set_visible(False)

    plt.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] wrote {out_png}")


def write_beam_table(
    features: FeatureTable,
    classification: ClassificationResult,
    out_csv: str,
) -> int:
    # one row per flagged timestep with the per-beam derived quantities
    import csv
    from datetime import datetime, timezone

    idxs = np.where(classification.is_beam)[0]
    cols = ["ut", "unix", "direction", "e_beam_eV", "delta_e_eV", "eb_over_de",
            "r_beam", "pa_max_ratio", "asymmetry", "e_peak_eV", "beam_score"]
    with open(out_csv, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(cols)
        for t in idxs:
            ut = datetime.fromtimestamp(features.times[t], tz=timezone.utc)
            wr.writerow([
                ut.strftime("%Y-%m-%d %H:%M:%S"),
                f"{features.times[t]:.0f}",
                int(classification.beam_direction[t]),
                f"{features.e_line[t]:.1f}",
                f"{features.de_line[t]:.1f}",
                f"{features.eb_over_de[t]:.3f}",
                f"{features.r_beam[t]:.3f}",
                f"{features.pa_max_ratio[t]:.3f}",
                f"{features.asymmetry[t]:.3f}",
                f"{features.e_peak[t]:.1f}",
                f"{classification.beam_score[t]:.3f}",
            ])
    print(f"[OK] wrote {out_csv} ({len(idxs)} beams)")
    return len(idxs)


def plot_beam_histograms(
    features: FeatureTable,
    classification: ClassificationResult,
    out_png: str,
    title: str = "",
):
    # 4-panel hist of the per-beam values over flagged timesteps
    import matplotlib.pyplot as plt

    is_b = classification.is_beam

    def _hist(ax, data, label, logx=False):
        d = data[is_b]
        d = d[np.isfinite(d)]
        if d.size == 0:
            ax.text(0.5, 0.5, f"no data\n{label}", ha="center", va="center")
            ax.set_title(label)
            return
        if logx and np.all(d > 0) and d.min() < d.max():
            bins = np.logspace(np.log10(d.min()), np.log10(d.max()), 25)
            ax.hist(d, bins=bins)
            ax.set_xscale("log")
        else:
            ax.hist(d, bins=25)
        ax.set_xlabel(label)
        ax.set_ylabel("count")
        ax.set_title(f"{label}  (n={d.size}, median={np.median(d):.2f})")

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    _hist(axes[0, 0], features.r_beam, "R (flux-weighted, beam)")
    _hist(axes[0, 1], features.e_line, "E_beam = e_line [eV]", logx=True)
    _hist(axes[1, 0], features.de_line, "delta E = line FWHM [eV]", logx=True)
    _hist(axes[1, 1], features.eb_over_de, "E_beam / delta E")
    if title:
        fig.suptitle(title)
    plt.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] wrote {out_png}")


def plot_threshold_comparison(
    spectra: PitchAngleSpectra,
    moments: dict,
    params: ClassifierParams,
    out_png: str,
    thresholds: tuple = (1.0, 1.2, 1.5, 2.0),
    energy_cutoff_low: float = 30.0,
    title: str = "",
):
    # small multiples, left col omni spectrogram per coherent_dir_min (R) value
    # beam detections overlaid as dots at the omni flux peak (e_peak)
    # right col E_beam/delta_E (eb_over_de) of those beams, dots colored by value
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from datetime import datetime, timezone

    times_dt = [datetime.fromtimestamp(t, tz=timezone.utc) for t in spectra.times]
    e_valid = spectra.energy > 0
    e_plot = spectra.energy[e_valid]
    z = spectra.omni[:, e_valid].T
    z = np.where(z > 0, z, np.nan)
    t_arr = np.array(times_dt, dtype=object)

    n = len(thresholds)
    fig, axes = plt.subplots(n, 2, figsize=(22, 3.2 * n), sharex=True)
    axes = np.atleast_2d(axes)

    # R threshold lives in extract_features, so re-extract per value
    results = []
    for thr in thresholds:
        feat = extract_features(spectra, moments,
                                energy_cutoff_low=energy_cutoff_low,
                                pa_coverage_threshold=params.min_coverage,
                                beam_flux_floor=params.beam_flux_floor,
                                coherent_asym_min=params.coherent_asym_min,
                                coherent_dir_min=thr,
                                coherent_min_bins=params.coherent_min_bins,
                                peak_width_max=params.peak_width_max,
                                peak_wlen=params.peak_wlen)
        results.append((feat, classify_beams(feat, params)))

    # shared eb_over_de y-range so rows stay comparable, median is per-row
    eb_all = [f.eb_over_de[c.is_beam & np.isfinite(f.eb_over_de)] for f, c in results]
    eb_all = np.concatenate(eb_all) if eb_all else np.array([])
    if eb_all.size and eb_all.max() > eb_all.min():
        pad = 0.05 * (eb_all.max() - eb_all.min())
        eb_lo, eb_hi = float(eb_all.min() - pad), float(eb_all.max() + pad)
    else:
        eb_lo, eb_hi = 0.0, 1.0

    for i, (thr, (feat, cls)) in enumerate(zip(thresholds, results)):
        ax = axes[i, 0]
        pcm = ax.pcolormesh(times_dt, e_plot, z,
                            norm=plt.matplotlib.colors.LogNorm(vmin=1e3, vmax=1e8),
                            cmap="jet", shading="auto")
        ax.set_yscale("log")
        ax.set_ylabel("Energy [eV]")
        ax.set_ylim(5, 30000)
        div = make_axes_locatable(ax)
        cax = div.append_axes("right", size="1.5%", pad=0.05)
        fig.colorbar(pcm, cax=cax, label="eflux")
        is_b = cls.is_beam
        for mask, col in ((is_b & (cls.beam_direction > 0), "red"),
                          (is_b & (cls.beam_direction < 0), "blue")):
            if np.any(mask):
                ax.scatter(t_arr[mask], feat.e_peak[mask],
                           s=_score_to_size(cls.beam_score[mask], params.score_threshold),
                           c=col, edgecolors="white", linewidths=0.5, zorder=3)
        ax.set_title(f"coherent_dir_min (R) >= {thr}   ({int(is_b.sum())} beams)")

        # right col, eb_over_de of detected beams colored by direction like left
        rax = axes[i, 1]
        rdiv = make_axes_locatable(rax)
        rcax = rdiv.append_axes("right", size="1.5%", pad=0.05)
        rcax.set_visible(False)  # spacer, keeps right panel width matching left
        ebok = np.isfinite(feat.eb_over_de)
        for mask, col in ((is_b & ebok & (cls.beam_direction > 0), "red"),
                          (is_b & ebok & (cls.beam_direction < 0), "blue")):
            if np.any(mask):
                rax.scatter(t_arr[mask], feat.eb_over_de[mask], s=20, c=col,
                            edgecolors="white", linewidths=0.5, zorder=3)
        ebmask = is_b & ebok
        if eb_all.size:
            rax.set_ylim(eb_lo, eb_hi)
        if np.any(ebmask):
            # per-row median of this threshold's detected beams
            row_med = float(np.median(feat.eb_over_de[ebmask]))
            rax.axhline(row_med, color="gray", ls="--", lw=0.8, alpha=0.7, zorder=2)
            rax.text(0.99, row_med, f"median {row_med:.2f}", color="gray",
                     fontsize=8, ha="right", va="bottom",
                     transform=rax.get_yaxis_transform())
        else:
            rax.text(0.5, 0.5, "no beams", ha="center", va="center",
                     transform=rax.transAxes)
        rax.set_ylabel("E_beam / delta E")
        rax.set_title(f"R >= {thr}: E_beam/delta_E")

    # time ticks under every row, smaller font since they now repeat
    for row in range(n):
        for j in (0, 1):
            a = axes[row, j]
            a.xaxis.set_major_locator(mdates.HourLocator(interval=1))
            a.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
            a.tick_params(axis="x", labelbottom=True, labelsize=7)
    for j in (0, 1):
        axes[-1, j].set_xlabel("UT")
    if title:
        fig.suptitle(title)
    # direction key centered below the suptitle, sign convention from classify_beams
    fig.text(0.49, 0.965, "red = parallel", color="red", fontsize=9,
             ha="right", va="top")
    fig.text(0.51, 0.965, "blue = antiparallel", color="blue", fontsize=9,
             ha="left", va="top")
    fig.text(0.5, 0.945, "dot size in left col grows with beam score", color="gray",
             fontsize=7, ha="center", va="top")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] wrote {out_png}")


@dataclass
class PipelineResult:
    spectra: PitchAngleSpectra
    features: FeatureTable
    classification: ClassificationResult
    classification_smoothed: ClassificationResult
    sensitivity: dict = field(default_factory=dict)


def run_pipeline(
    probe: str,
    trange: list[str],
    data_dir: str,
    params: ClassifierParams | None = None,
    min_consecutive: int = 1,
    energy_cutoff_low: float = 30.0,
    figures_dir: str | None = None,
    threshold_compare_values: tuple = (1.0, 1.2, 1.5, 2.0),
) -> PipelineResult:
    if params is None:
        params = ClassifierParams()
    print(f"=== Phase 0: Loading data for THEMIS-{probe.upper()} {trange} ===")
    dist = load_esd_distribution(probe, trange, data_dir)
    print(f"  ESD: {len(dist.times)} timesteps, {len(dist.energy)} energy bins")

    b_times, b_dsl = load_bfield_dsl(probe, trange, data_dir)
    print(f"  FGM: {len(b_times)} samples")

    moments = load_moments(probe, trange, data_dir)
    print(f"  Moments: density={len(moments.get('density', []))}, "
          f"velocity={len(moments.get('velocity', []))}")

    print(f"\n=== Phase 1: Pitch-angle spectral reduction ===")
    spectra = compute_pa_spectra(dist, b_times, b_dsl)
    n_valid_para = np.sum(spectra.pa_coverage_para > 0)
    n_valid_anti = np.sum(spectra.pa_coverage_anti > 0)
    print(f"  Para coverage: {n_valid_para}/{len(spectra.times)} timesteps")
    print(f"  Anti coverage: {n_valid_anti}/{len(spectra.times)} timesteps")

    print(f"\n=== Phase 2: Feature extraction ===")
    features = extract_features(spectra, moments,
                                energy_cutoff_low=energy_cutoff_low,
                                pa_coverage_threshold=params.min_coverage,
                                beam_flux_floor=params.beam_flux_floor,
                                coherent_asym_min=params.coherent_asym_min,
                                coherent_dir_min=params.coherent_dir_min,
                                coherent_min_bins=params.coherent_min_bins,
                                peak_width_max=params.peak_width_max,
                                peak_wlen=params.peak_wlen)
    n_finite = np.sum(np.isfinite(features.asymmetry))
    print(f"  Features computed: {n_finite}/{len(features.times)} with valid asymmetry")

    print(f"\n=== Phase 3: Classification ===")
    classification = classify_beams(features, params)
    n_beam = classification.is_beam.sum()
    print(f"  Raw beams: {n_beam}/{len(features.times)} timesteps")

    print(f"\n=== Phase 4: Temporal smoothing ===")
    smoothed = smooth_labels(classification, min_consecutive=min_consecutive)
    n_smooth = smoothed.is_beam.sum()
    print(f"  Smoothed beams: {n_smooth}/{len(features.times)} timesteps")

    sensitivity = threshold_sensitivity(features)
    print(f"  Sensitivity analysis complete")

    if figures_dir is not None:
        print(f"\n=== Phase 5: Plotting ===")
        date_str = trange[0].split("/")[0]
        # per-day subfolder keeps the top-level figures dir clean
        fig_dir = Path(figures_dir) / date_str
        fig_dir.mkdir(parents=True, exist_ok=True)

        prefix = f"th{probe}_beam_{date_str}"

        plot_feature_timeseries(
            features, smoothed, spectra,
            str(fig_dir / f"{prefix}_overview.png"),
            params,
            title=f"THEMIS-{probe.upper()} Beam Detection {trange[0]}",
        )

        plot_curated_snapshots(
            spectra, smoothed, features,
            str(fig_dir / f"{prefix}_snapshots.png"),
        )

        write_beam_table(
            features, smoothed,
            str(fig_dir / f"{prefix}_beams.csv"),
        )

        plot_beam_histograms(
            features, smoothed,
            str(fig_dir / f"{prefix}_histograms.png"),
            title=f"THEMIS-{probe.upper()} Beam Distributions {trange[0]}",
        )

        plot_threshold_comparison(
            spectra, moments, params,
            str(fig_dir / f"{prefix}_threshold_compare.png"),
            thresholds=tuple(threshold_compare_values),
            energy_cutoff_low=energy_cutoff_low,
            title=f"THEMIS-{probe.upper()} R-threshold Comparison {trange[0]}",
        )

    return PipelineResult(
        spectra=spectra,
        features=features,
        classification=classification,
        classification_smoothed=smoothed,
        sensitivity=sensitivity,
    )