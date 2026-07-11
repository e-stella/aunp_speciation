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


def load_series(path, delimiter=None, wavelength_range=(420, 800),
                normalize=None, anchor_nm=400.0):
    """Load a spectra file. Returns (wavelength_nm, spectra, temps_K_or_None).

    wavelength_range=(min_nm, max_nm) restricts the returned arrays to that
    (inclusive) window; pass None to keep the full range. The default
    (420, 800) nm trims the deep-UV interband region the Etchegoin dielectric
    does not model and the 350 nm lamp-changeover artifact seen on Cary
    instruments (see "Known limitations" #6 in CLAUDE.md).

    normalize (default None = no-op, keeps synthetic/pre-cleaned data as-is):
      "mult_400nm" — multiplicative concentration normalization anchored at
      anchor_nm: A_norm(l, T) = A(l, T) * A(anchor, T_ref) / A(anchor, T),
      with T_ref the FIRST ext_ column in the file. At ~400 nm gold extinction
      is interband-dominated and ~speciation-invariant, so it tracks total
      concentration; a multiplicative scale is the Beer-Lambert-correct
      drift correction (see README_experimental_data.md).

    ORDER IS FIXED: normalization runs FIRST (the anchor must still be in the
    array), THEN the wavelength_range clip. anchor_nm and wavelength_range are
    independent — the 400 nm anchor deliberately sits outside the 420-800 nm
    fit window, so feed RAW full-range files (e.g. 390-900 nm), not pre-clipped.
    """
    with open(path, newline="") as f:
        sample = f.read(2048)
        f.seek(0)
        if delimiter is None:
            delimiter = "\t" if sample.count("\t") > sample.count(",") else ","
        reader = csv.reader(f, delimiter=delimiter)
        header = next(reader)
        rows = [r for r in reader if r and r[0].strip() != ""]
    arr = np.array([[float(x) for x in r] for r in rows], dtype=float)
    wl = arr[:, 0]
    Y = arr[:, 1:].T  # (n_columns, n_wl)
    if normalize is not None:
        if normalize != "mult_400nm":
            raise ValueError(f"unknown normalize mode: {normalize!r} "
                             "(expected None or 'mult_400nm')")
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
    col_labels = header[1:]
    temps_C = [_parse_temperature(lbl) for lbl in col_labels]
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
