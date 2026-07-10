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


def load_series(path, delimiter=None, wavelength_range=(420, 800)):
    """Load a spectra file. Returns (wavelength_nm, spectra, temps_K_or_None).

    wavelength_range=(min_nm, max_nm) restricts the returned arrays to that
    (inclusive) window; pass None to keep the full range. The default
    (420, 800) nm trims the deep-UV interband region the Etchegoin dielectric
    does not model and the 350 nm lamp-changeover artifact seen on Cary
    instruments (see "Known limitations" #7 in CLAUDE.md).
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
