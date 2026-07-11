"""Export the per-temperature water-blank scans (and the true scan order)
from the 2011 Cary run's final.xls `raw` sheet.

Usage:
    python scripts/export_blanks.py /path/to/final.xls

Writes, next to the RAW sample CSVs (experimental/ESK_2011/):
  - blanks_heating_390-900.csv  (wavelength_nm + one blank_<T>C column per T)
  - blanks_cooling_390-900.csv
  - scan_order.csv              (chronological: scan index, kind, branch, T)

Why: the blanks are NON-ZERO (~0.03-0.05 in the far red — exactly where the
pedestal lives) and the instrument's baseline correction did NOT remove them.
Subtract them per temperature before fitting (io_data.load_series(blanks=...)).
Their smallness of T-drift (~-0.002 over 15->75 C) is itself the evidence that
the pedestal's T-growth is not instrumental.

Needs `xlrd` (pip install xlrd). The .xls lives outside the repo (2011 archive).
"""
import sys, os
import numpy as np

HERE = os.path.dirname(__file__)
OUTDIR = os.path.join(HERE, "..", "experimental", "ESK_2011")
XLS = sys.argv[1] if len(sys.argv) > 1 else \
    "/Users/eskoh/Documents/Projects/aunp-speciation_old_exp/final.xls"

import xlrd
wb = xlrd.open_workbook(XLS)
sh = wb.sheet_by_name("raw")

GRID = np.arange(390.0, 901.0, 1.0)


def read_scan(c):
    """Scan at column pair (c, c+1) -> (label, ext on the 1-nm GRID)."""
    label = str(sh.cell_value(0, c)).strip()
    wl, ab = [], []
    for r in range(2, sh.nrows):
        v, a = sh.cell_value(r, c), sh.cell_value(r, c + 1)
        if v == "" or a == "":
            break
        wl.append(float(v)); ab.append(float(a))
    wl, ab = np.asarray(wl), np.asarray(ab)
    order = np.argsort(wl)
    return label, np.interp(GRID, wl[order], ab[order])


# ---- parse all scans in chronological (column) order ----
scans = [read_scan(c) for c in range(0, sh.ncols, 2)]
print(f"parsed {len(scans)} scans from {os.path.basename(XLS)}")

# classify: 'Cell1 Blank' -> blank; '<T> deg' -> heating sample; '<T> deg-2'
# (or trailing '-2') -> cooling sample. Blanks belong to the NEXT sample's T.
seq = []   # (kind, branch, T_C, spectrum)
for label, spec in scans:
    if "blank" in label.lower():
        seq.append(["blank", None, None, spec])
    else:
        t = label.lower().replace("deg", "").strip()
        cooling = t.endswith("-2")
        t = float(t[:-2] if cooling else t)
        seq.append(["sample", "cooling" if cooling else "heating", t, spec])
# assign blanks forward to the next sample scan's (branch, T)
for i, s in enumerate(seq):
    if s[0] == "blank":
        for j in range(i + 1, len(seq)):
            if seq[j][0] == "sample":
                s[1], s[2] = seq[j][1], seq[j][2]
                break

# ---- average blanks per (branch, T), preserve first-appearance order ----
def write_branch(branch):
    temps, blanks = [], {}
    for kind, br, t, spec in seq:
        if kind == "blank" and br == branch:
            blanks.setdefault(t, []).append(spec)
            if t not in temps:
                temps.append(t)
    path = os.path.join(OUTDIR, f"blanks_{branch}_390-900.csv")
    with open(path, "w") as f:
        f.write("wavelength_nm," + ",".join(f"blank_{t:g}C" for t in temps) + "\n")
        M = np.column_stack([np.mean(blanks[t], axis=0) for t in temps])
        for i, lam in enumerate(GRID):
            f.write(f"{lam:.1f}," + ",".join(f"{x:.10f}" for x in M[i]) + "\n")
    n = {t: len(blanks[t]) for t in temps}
    print(f"wrote {path}  ({len(temps)} temps; scans per T: {n})")
    return temps, M


write_branch("heating")
write_branch("cooling")

# ---- true chronological scan order (for the evaporation E_i, optional) ----
path = os.path.join(OUTDIR, "scan_order.csv")
with open(path, "w") as f:
    f.write("index,kind,branch,temperature_C\n")
    for i, (kind, br, t, _) in enumerate(seq):
        f.write(f"{i},{kind},{br},{t:g}\n")
print(f"wrote {path}  ({len(seq)} scans)")

# quick diagnostics the docs cite
b15 = np.mean([s for k, b, t, s in seq if k == "blank" and b == "heating" and t == 15], axis=0)
b75 = np.mean([s for k, b, t, s in seq if k == "blank" and b == "heating" and t == 75], axis=0)
i7, i8 = np.argmin(abs(GRID - 700)), np.argmin(abs(GRID - 790))
print(f"blank far-red level @700/790 nm (15C): {b15[i7]:.4f} / {b15[i8]:.4f}")
print(f"blank T-drift 15->75C @700/790 nm: {b75[i7]-b15[i7]:+.4f} / {b75[i8]-b15[i8]:+.4f}")
