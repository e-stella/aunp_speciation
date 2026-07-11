"""Load experimental UV-Vis spectra for fitting.

Supported layouts (CSV/TSV, with header):
  1) Wide temperature series: first column = wavelength (nm); each remaining
     column = extinction at one temperature, header carrying the temperature
     e.g. "5C", "20 C", "T=35", "50degC". Temperatures parsed from headers.
  2) Single spectrum: first column wavelength, second column extinction.

Returns wavelength (nm), spectra (nT, nwl), and temperatures in Kelvin (or None
if not a temperature series). Baseline/units are left as-is; the fitters work on
normalized shapes, so absolute scaling and constant offsets are tolerated (though
a proper baseline subtraction upstream is recommended).
"""

from __future__ import annotations
import csv
import re
import numpy as np


def _parse_temperature(label):
    """Extract a temperature in Celsius from a column header, else None."""
    m = re.search(r"(-?\d+(?:\.\d+)?)\s*(?:deg)?\s*C", label, re.IGNORECASE)
    if m:
        return float(m.group(1))
    m = re.search(r"T\s*=?\s*(-?\d+(?:\.\d+)?)", label, re.IGNORECASE)
    if m:
        return float(m.group(1))
    m = re.search(r"(-?\d+(?:\.\d+)?)", label)
    return float(m.group(1)) if m else None


def water_density(temperature_C):
    """Density of air-free water (kg/m^3), Kell 1975 polynomial (0-100 C)."""
    t = np.asarray(temperature_C, dtype=float)
    num = (999.83952 + 16.945176 * t - 7.9870401e-3 * t**2
           - 46.170461e-6 * t**3 + 105.56302e-9 * t**4
           - 280.54253e-12 * t**5)
    return num / (1.0 + 16.879850e-3 * t)


def antoine_psat(temperature_C):
    """Saturation vapour pressure of water (mmHg), Antoine equation:
    log10 p = 8.07131 - 1730.63/(233.426 + T[degC])  (valid ~1-100 C)."""
    t = np.asarray(temperature_C, dtype=float)
    return 10.0 ** (8.07131 - 1730.63 / (233.426 + t))


def _read_wide(path, delimiter=None):
    """Read a wide CSV -> (wavelength, Y (ncols, nwl), column labels)."""
    with open(path, newline="") as f:
        sample = f.read(2048)
        f.seek(0)
        if delimiter is None:
            delimiter = "\t" if sample.count("\t") > sample.count(",") else ","
        reader = csv.reader(f, delimiter=delimiter)
        header = next(reader)
        rows = [r for r in reader if r and r[0].strip() != ""]
    arr = np.array([[float(x) for x in r] for r in rows], dtype=float)
    return arr[:, 0], arr[:, 1:].T, header[1:]


def subtract_blanks(wl, Y, temps_C, blanks_path):
    """Subtract per-temperature water-blank spectra (scripts/export_blanks.py).

    The 2011 Cary blanks are NON-ZERO (~0.03 in the far red — exactly where
    the pedestal lives) and the instrument's baseline correction did NOT
    remove them. Blank columns are matched to sample columns by temperature
    and interpolated onto the sample wavelength grid. Raises if any sample
    temperature has no blank.
    """
    bwl, B, blabels = _read_wide(blanks_path)
    btemps = [_parse_temperature(l) for l in blabels]
    out = np.array(Y, dtype=float, copy=True)
    for i, t in enumerate(temps_C):
        js = [j for j, bt in enumerate(btemps)
              if bt is not None and abs(bt - t) < 0.01]
        if not js:
            raise ValueError(f"no blank column for T={t} C in {blanks_path}")
        out[i] = Y[i] - np.interp(wl, bwl, B[js[0]])
    return out


def _band_ratio(wl_h, yh, wl_c, yc, band):
    """Mean of cooling/heating over a wavelength band (cooling interpolated)."""
    m = (wl_h >= band[0]) & (wl_h <= band[1])
    return float(np.mean(np.interp(wl_h[m], wl_c, yc) / yh[m]))


def check_concentration_drift(heating, cooling, blue=(400.0, 470.0),
                              red=(700.0, 800.0), warn_pct=0.5, verbose=True):
    """Matched-temperature branch-ratio diagnostic. RUN THIS whenever a
    dataset has both a heating and a cooling branch — it would have caught
    the C500 evaporation in 2011.

    heating/cooling: (wl, Y, temps_C) tuples (ideally blank-subtracted).
    For each temperature present in BOTH branches, computes the ratio
    cooling/heating and decomposes it into:
      - a wavelength-FLAT percent offset (blue-band mean) => pure
        multiplicative CONCENTRATION drift (evaporation/dilution);
      - a wavelength-STRUCTURED residual (red-band minus blue-band, pp)
        => irreversible chemistry / pedestal residue.
    The two can COEXIST (they do in C500: flat ~+5.6% @15C plus ~+1.1 pp
    extra in the red); both components are reported, no single verdict is
    forced. If the flat component exceeds warn_pct percent anywhere, a
    WARNING is printed: normalize='density' is blind to this drift.

    Returns dict(temps, flat_pct, structure_pp, drift_detected,
    structured_detected).
    """
    wl_h, Yh, th = heating
    wl_c, Yc, tc = cooling
    temps, flat, struct = [], [], []
    for i, t in enumerate(th):
        js = [j for j, u in enumerate(tc) if abs(u - t) < 0.01]
        if not js:
            continue
        j = js[0]
        rb = _band_ratio(wl_h, Yh[i], wl_c, Yc[j], blue)
        rr = _band_ratio(wl_h, Yh[i], wl_c, Yc[j], red)
        temps.append(t)
        flat.append(100.0 * (rb - 1.0))
        struct.append(100.0 * (rr - rb))
    temps, flat, struct = map(np.array, (temps, flat, struct))
    drift = bool(np.any(np.abs(flat) > warn_pct))
    structured = bool(np.any(np.abs(struct) > warn_pct))
    if verbose:
        print("concentration-drift diagnostic (cooling/heating at matched T):")
        print("    T(C)   flat offset %   structure (red-blue, pp)")
        for t, f, s in zip(temps, flat, struct):
            print(f"   {t:5.0f}      {f:+6.2f}            {s:+5.2f}")
        if drift:
            print("  *** WARNING: wavelength-FLAT branch offset detected => the "
                  "sample CONCENTRATION drifted during the run (evaporation or "
                  "dilution). normalize='density' CANNOT see this and will be "
                  "biased — use normalize='evaporation' (needs both branches). "
                  "NB sealing alone is NOT sufficient (the sealed 2011 C500 "
                  "cell still lost ~5.6%): minimise headspace and weigh the "
                  "cell before/after the run. ***")
        if structured:
            print("  note: a wavelength-STRUCTURED component is also present "
                  "(irreversible chemistry / pedestal residue) — the flat and "
                  "structured parts coexist; treat them separately.")
    return dict(temps=temps, flat_pct=flat, structure_pp=struct,
                drift_detected=drift, structured_detected=structured)


def _exposures(th, tc, scan_order=None):
    """Cumulative Antoine vapour-pressure-weighted exposure E per column.

    Default (scan_order=None): simple reconstruction — one exposure unit per
    loaded column, heating columns in file order then cooling columns in file
    order. scan_order: path to a scan_order.csv (index,kind,branch,
    temperature_C — see scripts/export_blanks.py) giving the TRUE chronology
    incl. replicates and blanks; E for a column is then the cumulative
    exposure at that (branch, T)'s LAST sample scan. The convention only
    rescales the fitted alpha; c/c0 is convention-invariant to <0.2% (checked
    on C500).
    Returns (E_heat dict, E_cool dict) keyed by temperature.
    """
    if scan_order is None:
        # chronological reconstruction is independent of file COLUMN order:
        # a monotonic ramp means heating ran ascending, cooling descending
        E, Eh, Ec = 0.0, {}, {}
        for t in sorted(th):
            E += float(antoine_psat(t))
            Eh[float(t)] = E
        for t in sorted(tc, reverse=True):
            E += float(antoine_psat(t))
            Ec[float(t)] = E
        return Eh, Ec
    Eh, Ec, E = {}, {}, 0.0
    with open(scan_order, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            t = float(row["temperature_C"])
            E += float(antoine_psat(t))
            if row["kind"] == "sample":
                (Eh if row["branch"] == "heating" else Ec)[t] = E
    return Eh, Ec


def fit_evaporation_alpha(heating, cooling, scan_order=None,
                          band=(400.0, 470.0), verbose=True):
    """Calibrate the 1-parameter evaporation model MODEL-FREE (no optics, no
    basis) from the matched-temperature branch offsets:

        cooling(lam,T)/heating(lam,T) ~= 1 + alpha * (E_cool(T) - E_heat(T))

    with E the cumulative Antoine exposure (see _exposures). alpha is fitted
    by least squares over the matched temperatures, on the BLUE band
    (400-470 nm, centred on Haiss's 450 nm concentration reference) where the
    branch offset is purest concentration — NOT the red, which carries ~1 pp
    of extra (pedestal/irreversible) offset that biases alpha low by ~25%.
    Per-region alphas are returned as a robustness check: blue-vs-red
    divergence measures the pedestal's leakage into the calibration.

    Returns dict(alpha, scatter, residuals, temps, dE, per_region_alpha,
    E_heat, E_cool).
    """
    wl_h, Yh, th = heating
    wl_c, Yc, tc = cooling
    Eh, Ec = _exposures(th, tc, scan_order)

    regions = {"blue 400-470": (400.0, 470.0), "450 alone": (447.0, 453.0),
               "peak 500-560": (500.0, 560.0), "red 700-800": (700.0, 800.0),
               "full 400-800": (400.0, 800.0)}

    def alpha_for(bnd):
        r, dE = [], []
        for i, t in enumerate(th):
            js = [j for j, u in enumerate(tc) if abs(u - t) < 0.01]
            if not js:
                continue
            r.append(_band_ratio(wl_h, Yh[i], wl_c, Yc[js[0]], bnd) - 1.0)
            dE.append(Ec[float(tc[js[0]])] - Eh[float(t)])
        r, dE = np.array(r), np.array(dE)
        a = float(np.sum(r * dE) / np.sum(dE**2))
        per_pair = r / dE
        scatter = float(np.std(per_pair) / abs(a)) if a != 0 else np.inf
        return a, scatter, r - a * dE, dE

    per_region = {k: alpha_for(b)[:2] for k, b in regions.items()}
    alpha, scatter, resid, dE = alpha_for(band)
    temps = [t for t in th if any(abs(u - t) < 0.01 for u in tc)]
    if verbose:
        print(f"evaporation alpha (band {band[0]:.0f}-{band[1]:.0f} nm): "
              f"{alpha:.3g}  scatter {100*scatter:.0f}%  "
              f"residuals +/-{100*np.max(np.abs(resid)):.2f}%")
        print("  per-region robustness check:")
        for k, (a, s) in per_region.items():
            print(f"    {k:14s}: alpha={a:.3g}  scatter {100*s:.0f}%")
        ab = per_region["blue 400-470"][0]
        ar = per_region["red 700-800"][0]
        if abs(ab - ar) > 0.15 * abs(ab):
            print(f"  note: blue/red alphas diverge {100*abs(ab-ar)/abs(ab):.0f}% "
                  "-> the red carries non-concentration (pedestal/irreversible) "
                  "offset; this is the pedestal's leakage measure, and why the "
                  "calibration uses the blue.")
    return dict(alpha=alpha, scatter=scatter, residuals=resid, temps=temps,
                dE=dE, per_region_alpha=per_region, E_heat=Eh, E_cool=Ec)


def load_series(path, delimiter=None, wavelength_range=(420, 800),
                normalize=None, anchor_nm=400.0, companion=None,
                blanks=None, companion_blanks=None, scan_order=None,
                branch=None):
    """Load a spectra file. Returns (wavelength_nm, spectra, temps_K_or_None).

    wavelength_range=(min_nm, max_nm) restricts the returned arrays to that
    (inclusive) window; pass None to keep the full range. The default
    (420, 800) nm trims the deep-UV interband region the Etchegoin dielectric
    does not model and the 350 nm lamp-changeover artifact seen on Cary
    instruments (see "Known limitations" #6 in CLAUDE.md).

    normalize (default None = no-op, keeps synthetic/pre-cleaned data as-is):
      "density" — RECOMMENDED for temperature series. Multiplicative dilution
      correction from water thermal expansion only:
      A_norm(l, T) = A(l, T) * rho(T_ref) / rho(T), Kell (1975) density
      polynomial, T_ref the FIRST ext_ column's temperature. No anchor
      wavelength, no assumption about any spectral point being invariant.
      Requires parseable temperatures in every column header.
      "evaporation" — SALVAGE MODE for runs that lost solvent (limitation
      #14; NB the sealed 2011 C500 cell still lost ~5.6% — a seal alone does
      not guarantee anything). c(T_i)/c(T_ref) = [rho(T_i)/rho(T_ref)] *
      [1 + alpha*E_i], with
      E_i the cumulative Antoine vapour-pressure-weighted scan exposure and
      alpha ONE parameter FITTED per dataset (never hardcoded) from the
      matched-temperature heating-vs-cooling branch offsets (model-free —
      see fit_evaporation_alpha). REQUIRES `companion=` (the other branch's
      RAW file); with one branch only it falls back to "density" with a LOUD
      warning (evaporation is uncorrectable from a single branch). Assumes
      the entire flat branch offset is concentration — irreversible
      aggregation would be misread as evaporation. "density" is the right
      mode only when the concentration is KNOWN stable — verify by weighing
      the cell before/after the run (and minimise headspace), don't assume
      it from a seal.
      "mult_400nm" — DEPRECATED / BIASED. Anchors every spectrum to equal
      extinction at anchor_nm. Its premise (A(400) tracks concentration only)
      is measurably violated on the C500 series: A(400) RISES +3.18% over
      15->75 C while Kell expansion predicts -2.43%; forcing A(400)=const is a
      flat rescale that inflates the apparent plasmon-peak change ~8.5x and
      manufactures/destroys isosbestic structure (limitation #12). Kept for
      reproducing old results only.

    blanks / companion_blanks: paths to per-temperature blank CSVs
    (scripts/export_blanks.py) subtracted from this file / the companion
    BEFORE anything else — the 2011 blanks carry a real ~0.03 far-red offset.
    scan_order: optional path to scan_order.csv for the true chronology
    (default: simple reconstruction — heating ascending then cooling
    descending; changes c/c0 by <0.2% on C500).

    ORDER IS FIXED: blank subtraction, then normalization (the anchor must
    still be in the array), THEN the wavelength_range clip. anchor_nm and
    wavelength_range are independent — the 400 nm anchor deliberately sits
    outside the 420-800 nm fit window, so feed RAW full-range files (e.g.
    390-900 nm), not pre-clipped.
    """
    wl, Y, col_labels = _read_wide(path, delimiter)
    temps_C = [_parse_temperature(lbl) for lbl in col_labels]

    def need_temps(mode):
        if any(t is None for t in temps_C):
            raise ValueError(f"normalize={mode!r} needs a parseable temperature "
                             f"in every column header; got {col_labels}")

    if blanks is not None:
        need_temps("blanks subtraction")
        Y = subtract_blanks(wl, Y, temps_C, blanks)

    if normalize == "evaporation" and companion is None:
        print("*** WARNING: normalize='evaporation' needs BOTH branches "
              "(companion=...). Evaporation is UNCORRECTABLE from a single "
              "branch — falling back to 'density', which corrects thermal "
              "expansion only and is BIASED if the sample lost solvent. ***")
        normalize = "density"

    if normalize == "density":
        need_temps("density")
        rho = water_density(np.array(temps_C))
        Y = Y * (rho[0] / rho)[:, None]          # T_ref = first ext_ column
    elif normalize == "evaporation":
        need_temps("evaporation")
        wl2, Y2, labels2 = _read_wide(companion, delimiter)
        temps2 = [_parse_temperature(lbl) for lbl in labels2]
        if any(t is None for t in temps2):
            raise ValueError("companion file needs parseable temperatures")
        if companion_blanks is not None:
            Y2 = subtract_blanks(wl2, Y2, temps2, companion_blanks)
        elif blanks is not None:
            print("  warning: blanks subtracted from the loaded file but no "
                  "companion_blanks given — the alpha calibration ratio mixes "
                  "blank-subtracted and raw branches.")
        # Which branch is THIS file? Prefer an explicit branch=; the fallback
        # heuristic (ascending columns = heating) FAILS when a cooling file was
        # saved with ascending column order (true of the C500 cooling CSV).
        if branch is None:
            this_is_heating = temps_C[0] < temps_C[-1]
            print(f"  branch not given — assuming this file is "
                  f"{'HEATING' if this_is_heating else 'COOLING'} from its "
                  "column order (pass branch='heating'/'cooling' to be sure; "
                  "column order is NOT always chronological).")
        else:
            if branch not in ("heating", "cooling"):
                raise ValueError("branch must be 'heating' or 'cooling'")
            this_is_heating = branch == "heating"
        heat = (wl, Y, temps_C) if this_is_heating else (wl2, Y2, temps2)
        cool = (wl2, Y2, temps2) if this_is_heating else (wl, Y, temps_C)
        print("*** normalize='evaporation' is a SALVAGE mode for runs that "
              "lost solvent: it assumes the ENTIRE flat branch offset is "
              "concentration; irreversible aggregation would be misread as "
              "evaporation. Next time minimise headspace and weigh the cell "
              "before/after (the sealed 2011 cell still lost ~5.6%). ***")
        check_concentration_drift(heat, cool)
        cal = fit_evaporation_alpha(heat, cool, scan_order=scan_order)
        a, Eh, Ec = cal["alpha"], cal["E_heat"], cal["E_cool"]
        th, tc_ = heat[2], cool[2]
        rho_ref = water_density(th[0])
        E_ref = Eh[float(th[0])]
        Emap = Eh if this_is_heating else Ec
        c_over_c0 = np.array([
            (water_density(t) / rho_ref) * (1 + a * Emap[float(t)]) / (1 + a * E_ref)
            for t in temps_C])
        print("  applied c/c0 (this file):",
              " ".join(f"{t:g}C:{c:.4f}" for t, c in zip(temps_C, c_over_c0)))
        Y = Y / c_over_c0[:, None]
    elif normalize is not None:
        if normalize != "mult_400nm":
            raise ValueError(f"unknown normalize mode: {normalize!r} (expected "
                             "None, 'density', 'evaporation' or 'mult_400nm')")
        if not (wl.min() <= anchor_nm <= wl.max()):
            raise ValueError(
                f"anchor_nm={anchor_nm} nm is outside the loaded data "
                f"({wl.min():.1f}-{wl.max():.1f} nm). Normalization needs the "
                "anchor in the RAW file (it runs before the wavelength_range "
                "clip) — feed the full-range RAW CSV, not a pre-clipped one.")
        ia = int(np.argmin(np.abs(wl - anchor_nm)))
        Y = Y * (Y[0, ia] / Y[:, ia])[:, None]   # T_ref = first ext_ column
    if wavelength_range is not None:
        lo, hi = wavelength_range
        mask = (wl >= lo) & (wl <= hi)
        wl = wl[mask]
        Y = Y[:, mask]
    if Y.shape[0] >= 2 and all(t is not None for t in temps_C):
        temps_K = np.array(temps_C) + 273.15
        order = np.argsort(temps_K)
        return wl, Y[order], temps_K[order]
    return wl, Y, None


def save_series(path, wl, spectra, temps_C, header_prefix="ext"):
    """Write a wide temperature-series CSV (used to make example data)."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["wavelength_nm"] + [f"{header_prefix}_{t:g}C" for t in temps_C])
        for i, lam in enumerate(wl):
            w.writerow([f"{lam:.1f}"] + [f"{spectra[j, i]:.6f}"
                                         for j in range(len(temps_C))])
