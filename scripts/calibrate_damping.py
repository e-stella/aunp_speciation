"""Joint damping calibration on the CTAC monomer set (all five samples).

Model: gamma_eff(D) = s * gamma_bulk + A * hbar*v_F / (D/2), with ONE global
(s, A) pair for the whole size series (A is a property of the interface, s of
the crystallinity — neither may depend on size; the 1/R law carries the size
dependence). Grid-scans (s, A), minimizing the summed peak-normalized RMS
(450-700 nm) over the five TEM-pinned per-particle histograms, with the
documented -2.7 nm offset (CLAUDE.md #9) applied.

The result is meant to be FROZEN and validated out-of-sample on the
seeded-growth citrate series (scripts/predict_irina_seeded.py) — the CTAC set
is the calibration set, the seeded set the test set.

Run:  python scripts/calibrate_damping.py     (~2 min)
"""
import sys, os, csv
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import pandas as pd

from aunp_speciation import dielectric
dielectric.use_gold_model("jc")
from aunp_speciation.mie import mie_ab
from aunp_speciation.dielectric import gold_epsilon, medium_index

ROOT = os.path.join(os.path.dirname(__file__), "..")
WL = np.arange(420.0, 801.0, 1.0)
E_EV = 1239.842 / WL
EPS_JC = gold_epsilon(WL)
N_MED = medium_index("water", None).real
OM = getattr(dielectric, "_OMEGA_P_EV", 8.45)
G0 = getattr(dielectric, "_GAMMA_BULK_EV", 0.044)
HBAR_VF = 0.9215
SHIFT_NM = 2.7

S_GRID = (0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0)
A_GRID = (0.25, 0.4, 0.5, 0.6, 0.75, 1.0)


def drude(g):
    return OM ** 2 / (E_EV ** 2 + 1j * g * E_EV)


def xsec(D, s, A):
    eps = EPS_JC + drude(G0) - drude(s * G0 + A * HBAR_VF / (D / 2.0))
    m = np.sqrt(eps) / N_MED
    C = np.zeros_like(WL)
    for i, l0 in enumerate(WL):
        k = 2.0 * np.pi * N_MED / l0
        a, b = mie_ab(m[i], k * D / 2.0)
        n = np.arange(1, len(a) + 1)
        C[i] = (2.0 * np.pi / k ** 2) * np.sum((2 * n + 1) * np.real(a + b))
    return C


tem = pd.read_excel(os.path.join(ROOT, "experimental/GNP_CTAC_TEM_UV.xlsx"),
                    sheet_name="Sheet1")
with open(os.path.join(ROOT, "experimental/ctac/ctac_uvvis.csv")) as f:
    rows = list(csv.reader(f))
uv = np.array([[float(x) for x in r] for r in rows[1:]])
o = np.argsort(uv[:, 0]); uv = uv[o]

SAMPLES = []
for i, col in enumerate(tem.columns):
    d = tem[col].dropna().to_numpy(float)
    edges = np.arange(d.min() - 0.25, d.max() + 0.75, 0.5)
    cnt, _ = np.histogram(d, edges)
    ctr = 0.5 * (edges[:-1] + edges[1:])
    keep = cnt > 0
    SAMPLES.append(dict(lab=rows[0][i + 1].replace("Abs_GNP_", "").replace(" ", ""),
                        bins=ctr[keep], w=cnt[keep] / cnt.sum(),
                        y=np.interp(WL, uv[:, 0], uv[:, i + 1])))

MSK = (WL >= 450) & (WL <= 700)


def rms(ym, y):
    a = ym[MSK] / ym[MSK].max()
    b = y[MSK] / y[MSK].max()
    return 100 * np.sqrt(np.mean((a - b) ** 2))


best = (1e9, None, None, None)
for s in S_GRID:
    for A in A_GRID:
        per = []
        for S in SAMPLES:
            ym = np.sum([wi * xsec(d, s, A) for d, wi in zip(S["bins"], S["w"])],
                        axis=0)
            ym = np.interp(WL, WL - SHIFT_NM, ym)      # documented offset
            per.append(rms(ym, S["y"]))
        tot = float(np.mean(per))
        if tot < best[0]:
            best = (tot, s, A, per)
print(f"grid: s in {S_GRID}, A in {A_GRID}, shift -{SHIFT_NM} nm applied")
print(f"BEST: s = {best[1]}, A_surf = {best[2]}  (mean RMS {best[0]:.2f}%)")
for S, r in zip(SAMPLES, best[3]):
    print(f"  {S['lab']:6s}: RMS {r:.2f}%")
print("\nBaseline for comparison (s=1, A=0.25, no shift):")
for S in SAMPLES:
    ym = np.sum([wi * xsec(d, 1.0, 0.25) for d, wi in zip(S["bins"], S["w"])],
                axis=0)
    print(f"  {S['lab']:6s}: RMS {rms(ym, S['y']):.2f}%")
