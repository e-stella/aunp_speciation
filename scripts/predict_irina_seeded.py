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
from aunp_speciation.mie import mie_ab
from aunp_speciation.dielectric import gold_epsilon, medium_index

ROOT = os.path.join(os.path.dirname(__file__), "..")
IRINA = os.path.join(ROOT, "experimental", "Irina")
OUT = os.path.join(ROOT, "outputs")
A_SURF = 0.25
WLM = np.arange(400.0, 1001.0, 1.0)      # model grid

# CALIBRATED damping (frozen 2026-07-12, scripts/calibrate_damping.py):
# gamma_eff = S_CAL*gamma_bulk + A_CAL*hbar*v_F/R, one global pair fitted
# jointly on the five CTAC monomer samples (mean RMS 4.6% -> 2.5%), plus the
# documented -2.7 nm offset (CLAUDE.md #9). The seeded-growth series below is
# OUT-OF-SAMPLE for this calibration (different chemistry: citrate vs CTAC).
S_CAL, A_CAL, SHIFT_CAL = 0.05, 0.6, 2.7
_E_EV = 1239.842 / WLM
_EPS_JC = gold_epsilon(WLM)
_N_MED = medium_index("water", None).real
_OM = getattr(dielectric, "_OMEGA_P_EV", 8.45)
_G0 = getattr(dielectric, "_GAMMA_BULK_EV", 0.044)


def _xsec_cal(D):
    g = S_CAL * _G0 + A_CAL * 0.9215 / (D / 2.0)
    eps = _EPS_JC + _OM**2/(_E_EV**2 + 1j*_G0*_E_EV) \
        - _OM**2/(_E_EV**2 + 1j*g*_E_EV)
    m = np.sqrt(eps) / _N_MED
    C = np.zeros_like(WLM)
    for i, l0 in enumerate(WLM):
        k = 2.0 * np.pi * _N_MED / l0
        a, b = mie_ab(m[i], k * D / 2.0)
        n = np.arange(1, len(a) + 1)
        C[i] = (2.0 * np.pi / k**2) * np.sum((2*n+1) * np.real(a + b))
    return C


def calibrated_model(D0, poly_pct):
    sd = D0 * poly_pct / 100.0
    ds = np.linspace(D0 - 3*sd, D0 + 3*sd, 15)
    w = np.exp(-0.5*((ds - D0)/sd)**2); w /= w.sum()
    y = np.sum([wi * _xsec_cal(d) for d, wi in zip(ds, w)], axis=0)
    return np.interp(WLM, WLM - SHIFT_CAL, y)    # -2.7 nm offset

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
      f"{'pk data':>10s} {'dpk_bas':>6s} {'dpk_cal':>6s} "
      f"{'A700 d':>8s} {'bas':>6s} {'cal':>6s} {'RMSb':>6s} {'RMSc':>6s}")
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
        # zero-free-parameter predictions: baseline stack + CTAC-calibrated
        mod = monomer_polydisperse(WLM, D, poly, "water", None, True, A_SURF)
        mod = mod / mod.max()
        mcal = calibrated_model(D, poly)
        mcal = mcal / mcal.max()
        mfit = (wl >= 400) & (wl <= 1000)
        yd = Y[mfit, i] / Y[mfit, i].max()
        axm.plot(wl[mfit], yd, "-", lw=1.0, color=col, alpha=0.45)
        axm.plot(WLM, mod, ":", color=col, lw=1.0, alpha=0.8)
        axm.plot(WLM, mcal, "--", color=col, lw=1.5)
        pk_d = peak_of(wl[mfit], yd)
        pk_m = peak_of(WLM, mod)
        pk_c = peak_of(WLM, mcal)
        i700d = np.argmin(np.abs(wl[mfit] - 700.0))
        i700m = np.argmin(np.abs(WLM - 700.0))
        # RMS on peak-normalized curves, max(lo,450)-700 nm
        lo = max(S_ANCH := 450.0, 480.0 if label.startswith("Nov") or label.startswith("Oct") else 450.0)
        band_d = (wl[mfit] >= lo) & (wl[mfit] <= 700)
        band_m = (WLM >= lo) & (WLM <= 700)
        ydb = np.interp(WLM[band_m], wl[mfit][band_d], yd[band_d])
        rms_b = 100*np.sqrt(np.mean((mod[band_m]/mod[band_m].max() - ydb/ydb.max())**2))
        rms_c = 100*np.sqrt(np.mean((mcal[band_m]/mcal[band_m].max() - ydb/ydb.max())**2))
        print(f"{label:9s} {disp:14s} "
              f"{D:7.2f} {poly:6.1f} {pk_d:10.0f} {pk_m-pk_d:+6.0f} {pk_c-pk_d:+6.0f} "
              f"{yd[i700d]:8.3f} {mod[i700m]:6.3f} {mcal[i700m]:6.3f} "
              f"{rms_b:6.2f} {rms_c:6.2f}")
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
        Line2D([], [], color="0.2", lw=1.0, ls=":", label="dotted = baseline model"),
        Line2D([], [], color="0.2", lw=1.5, ls="--",
               label="dashed = CTAC-calibrated (out-of-sample)")],
        frameon=False, fontsize=7)
fig.suptitle("Seeded-growth GNP series: TEM-pinned zero-free-parameter predictions — "
             "baseline vs CTAC-calibrated damping (s=0.05, A=0.6, −2.7 nm; out-of-sample)",
             fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig17_irina_seeded.png"))
print("\nwrote outputs/fig17_irina_seeded.png")
