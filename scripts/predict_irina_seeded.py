"""Zero-free-parameter spectrum PREDICTION for the seeded-growth GNP series
(Oct + Nov 2024, citrate/H2O2 seeded growth, unpurified; data in
experimental/Irina/).

For each sample the model input is ONLY the TEM mean diameter and TEM
polydispersity (from the accompanying *_TEM_measurements files); the optics is
the production stack (Mie, jc gold, gamma_S A_surf=0.25, water). Nothing is
fitted — both data and model are normalized to their own plasmon peak, so the
comparison is pure shape prediction.

Also serves as the CSV-vs-pptx consistency check: panels (a)/(c) replot the
CSVs in the same style as the UV-Vis slides in GNP_Oct2024.pptx /
GNP_Nov2024.pptx for visual comparison (peak positions, ordering, tails).

Known sample caveats (from the slides / TEM notes):
  - samples are UNPURIFIED: H2O2 saturates the detector below ~300 nm and the
    residual reaction matrix may shift the medium index slightly above water;
  - the Nov "55 nm" sample is flagged by the experimenter as aggregated
    ("need to redo at different concentration") — the prediction is EXPECTED
    to fail there; it is the negative control.

Run:  python scripts/predict_irina_seeded.py
Writes outputs/fig17_irina_seeded.png and prints a per-sample table.
"""
import sys, os, csv, io
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np

from aunp_speciation import dielectric
dielectric.use_gold_model("jc")
from aunp_speciation.spectra import monomer_polydisperse

ROOT = os.path.join(os.path.dirname(__file__), "..")
IRINA = os.path.join(ROOT, "experimental", "Irina")
OUT = os.path.join(ROOT, "outputs")
A_SURF = 0.25
WLM = np.arange(400.0, 1001.0, 1.0)      # model grid

# TEM inputs (mean nm, sd nm) transcribed from *_TEM_measurements files.
# 'names' maps the CSV column headers (synthesis TARGET names) to the
# TEM-anchored display names used on the slides' size-analysis pages.
SETS = {
    "Oct 2024": dict(
        csv="Scan_10_22_2024_5_43_15_PM.csv",
        tem={"GNP12nm": (11.55, 0.91), "GNP20nm": (21.58, 2.18),
             "GNP40nm": (41.59, 3.34), "GNP60nm": (61.27, 4.50)},
        names={}),
    "Nov 2024": dict(
        csv="Scan_11_16_2024_4_37_22_PM.csv",
        tem={"Seed 1(Abs)": (10.85, 0.83), "30nm 1(Abs)": (22.90, 1.46),
             "40nm 1(Abs)": (38.52, 3.37), "60nm 1(Abs)": (54.61, 4.89)},
        names={"Seed 1(Abs)": "GNP11nm (seed)", "30nm 1(Abs)": "GNP23nm",
               "40nm 1(Abs)": "GNP39nm", "60nm 1(Abs)": "GNP55nm"}),
}


def load_utf16_tsv(path):
    with open(path, "r", encoding="utf-16") as f:
        rows = [r for r in csv.reader(f, delimiter="\t") if r and r[0].strip()]
    hdr = [c.strip() for c in rows[0] if c.strip()]
    data = []
    for r in rows[1:]:
        try:
            vals = [float(c) for c in r[:len(hdr)] if c.strip() != ""]
        except ValueError:
            continue
        if len(vals) == len(hdr):
            data.append(vals)
    arr = np.array(data)
    o = np.argsort(arr[:, 0])
    return hdr[1:], arr[o, 0], arr[o, 1:]


def peak_of(wl, y, lo=480.0, hi=650.0):
    m = (wl >= lo) & (wl <= hi)
    k = np.convolve(y[m], np.ones(5) / 5, mode="same")   # 5-pt smooth vs noise
    return wl[m][np.argmax(k)]


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25, "axes.axisbelow": True})
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
COLS = ("#D62728", "#2CA02C", "#17BECF", "#9467BD")   # match the slides

print(f"{'set':9s} {'sample':14s} {'TEM D':>7s} {'poly%':>6s} "
      f"{'peak data':>10s} {'peak model':>11s} {'dpk':>5s} "
      f"{'A700/pk data':>13s} {'model':>6s}")
for row, (label, S) in enumerate(SETS.items()):
    hdr, wl, Y = load_utf16_tsv(os.path.join(IRINA, S["csv"]))
    axd, axm = axes[row]
    for i, (name, col) in enumerate(zip(hdr, COLS)):
        D, sd = S["tem"][name]
        disp = S["names"].get(name, name.replace("(Abs)", "").strip())
        poly = 100 * sd / D
        m = (wl >= 300) & (wl <= 1000)
        y = Y[m, i]
        pk_idx = (wl[m] >= 400) & (wl[m] <= 700)
        y = y / y[pk_idx].max()
        axd.plot(wl[m], y, color=col, lw=1.2,
                 label=f"{disp} (TEM {D:.1f}±{sd:.1f})")
        # zero-free-parameter prediction
        mod = monomer_polydisperse(WLM, D, poly, "water", None, True, A_SURF)
        mod = mod / mod.max()
        mfit = (wl >= 400) & (wl <= 1000)
        yd = Y[mfit, i] / Y[mfit, i].max()
        axm.plot(wl[mfit], yd, "-", lw=1.0, color=col, alpha=0.45)
        axm.plot(WLM, mod, "--", color=col, lw=1.5)
        pk_d = peak_of(wl[mfit], yd)
        pk_m = peak_of(WLM, mod)
        i700d = np.argmin(np.abs(wl[mfit] - 700.0))
        r700d = yd[i700d]
        r700m = mod[np.argmin(np.abs(WLM - 700.0))]
        print(f"{label:9s} {disp:14s} "
              f"{D:7.2f} {poly:6.1f} {pk_d:10.0f} {pk_m:11.0f} "
              f"{pk_m-pk_d:+5.0f} {r700d:13.3f} {r700m:6.3f}")
    axd.set_xlabel("wavelength (nm)")
    axd.set_ylabel("extinction / peak")
    axd.set_title(f"({'ac'[row]}) {label} data replotted from CSV "
                  "(compare to the pptx UV-Vis slide)", fontsize=9)
    axd.legend(frameon=False, fontsize=7)
    axd.set_xlim(300, 1000); axd.set_ylim(0, 1.25)
    axm.set_xlabel("wavelength (nm)")
    axm.set_ylabel("extinction / peak")
    axm.set_title(f"({'bd'[row]}) {label}: TEM-pinned prediction, "
                  "ZERO fitted parameters", fontsize=9)
    axm.set_xlim(400, 1000); axm.set_ylim(0, 1.25)
    from matplotlib.lines import Line2D
    axm.legend(handles=[
        Line2D([], [], color="0.4", lw=1.0, alpha=0.6, label="solid = measured"),
        Line2D([], [], color="0.2", lw=1.5, ls="--", label="dashed = model")],
        frameon=False, fontsize=7.5)
fig.suptitle("Seeded-growth GNP series: CSV replot vs pptx + zero-free-"
             "parameter Mie prediction (jc, γ_S A=0.25, TEM-pinned D & poly)",
             fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig17_irina_seeded.png"))
print("\nwrote outputs/fig17_irina_seeded.png")
