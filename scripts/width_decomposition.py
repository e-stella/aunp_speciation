"""Width decomposition: WHY is the modeled plasmon band broader than the
measured one on verified-clean monomer samples — and how much of the gap does
each ingredient close?

The measured CTAC spectra are NARROWER than the TEM-pinned model (sizing fits
rail polydispersity at the grid floor). Candidate contributions, tested here
cumulatively on the CTAC samples (the set with true per-particle TEM lists):

  A  baseline: Gaussian(TEM mean, RAW sd) + jc + gamma_S(A=0.25)   [status quo]
  B  per-particle histogram instead of the Gaussian (uses every measured
     particle at its own diameter; outliers contribute their true extinction
     instead of inflating a symmetric width)
  C  B + bulk-damping scale s: gamma_eff = s*gamma_bulk + gamma_S. s<1 models
     single-crystal colloidal gold being LESS damped than the evaporated
     polycrystalline films behind the J&C table (Reddy 2016 measured a ~6x
     single-vs-poly damping difference). s is calibrated per sample here as a
     DIAGNOSTIC; adopting one global value is a follow-up decision.
  D  C + the documented -2.7 nm peak-position offset (CLAUDE.md #9)

Metrics on peak-normalized curves, 450-700 nm: full width at 0.75*max (FW75),
red-side half width (peak -> red 0.5*max), and RMS vs data.

Run:  python scripts/width_decomposition.py   (system env, Mie only, ~2 min)
Writes outputs/fig19_width_decomposition.png.
"""
import sys, os, csv
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np

from aunp_speciation import dielectric
dielectric.use_gold_model("jc")
from aunp_speciation.mie import mie_ab
from aunp_speciation.dielectric import gold_epsilon, medium_index

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "outputs")
WL = np.arange(420.0, 801.0, 1.0)
A_SURF = 0.25
OMEGA_P = getattr(dielectric, "_OMEGA_P_EV", 8.45)
GAMMA0 = getattr(dielectric, "_GAMMA_BULK_EV", 0.044)
HBAR_VF = 0.9215        # eV*nm; gamma_S = A * HBAR_VF / R  (0.143 eV @ D=12.9, A=1)
S_SCAN = np.arange(0.05, 1.01, 0.05)

E_EV = 1239.842 / WL
EPS_JC = gold_epsilon(WL)                       # base J&C table (spline)
N_MED = medium_index("water", None).real


def drude(gam):
    return OMEGA_P ** 2 / (E_EV ** 2 + 1j * gam * E_EV)


def cross_section(D, s):
    """Mie extinction with gamma_eff = s*gamma_bulk + gamma_S(D)."""
    g_new = s * GAMMA0 + A_SURF * HBAR_VF / (D / 2.0)
    eps = EPS_JC + drude(GAMMA0) - drude(g_new)
    m = np.sqrt(eps) / N_MED
    a_rad = D / 2.0
    C = np.zeros_like(WL)
    for i, l0 in enumerate(WL):
        k = 2.0 * np.pi * N_MED / l0
        a, b = mie_ab(m[i], k * a_rad)
        n = np.arange(1, len(a) + 1)
        C[i] = (2.0 * np.pi / k ** 2) * np.sum((2 * n + 1) * np.real(a + b))
    return C


# ------------------------------------------------------- per-particle TEM
import pandas as pd
tem = pd.read_excel(os.path.join(ROOT, "experimental/GNP_CTAC_TEM_UV.xlsx"),
                    sheet_name="Sheet1")
with open(os.path.join(ROOT, "experimental/ctac/ctac_uvvis.csv")) as f:
    rows = list(csv.reader(f))
uv_hdr = rows[0]
uv = np.array([[float(x) for x in r] for r in rows[1:]])
o = np.argsort(uv[:, 0]); uv = uv[o]
WLD = uv[:, 0]

SAMPLES = []
for i, col in enumerate(tem.columns):
    d = tem[col].dropna().to_numpy(float)
    lab = uv_hdr[i + 1].replace("Abs_GNP_", "").replace(" ", "")
    # 0.5 nm bins of the actual particle list
    edges = np.arange(d.min() - 0.25, d.max() + 0.75, 0.5)
    cnt, _ = np.histogram(d, edges)
    ctr = 0.5 * (edges[:-1] + edges[1:])
    keep = cnt > 0
    SAMPLES.append(dict(lab=lab, D=d.mean(), sd=d.std(ddof=1),
                        bins=ctr[keep], w=cnt[keep] / cnt.sum(),
                        y=np.interp(WL, WLD, uv[:, i + 1])))


def gauss_bins(D0, sd):
    ds = np.linspace(D0 - 3 * sd, D0 + 3 * sd, 15)
    w = np.exp(-0.5 * ((ds - D0) / sd) ** 2)
    return ds, w / w.sum()


def ensemble(bins, w, s):
    return np.sum([wi * cross_section(d, s) for d, wi in zip(bins, w)], axis=0)


def norm(y):
    m = (WL >= 450) & (WL <= 700)
    return y / y[m].max()


def widths(y):
    """FW at 0.75*max and red-side HWHM, on the peak-normalized curve."""
    m = (WL >= 450) & (WL <= 720)
    w, yy = WL[m], y[m] / y[m].max()
    ipk = np.argmax(yy)
    def x_at(level, side):
        seg = slice(ipk, None) if side > 0 else slice(None, ipk + 1)
        ws, ys = w[seg], yy[seg]
        if side < 0:
            ws, ys = ws[::-1], ys[::-1]
        below = np.where(ys <= level)[0]
        if len(below) == 0:
            return np.nan
        j = below[0]
        return np.interp(level, [ys[j], ys[j - 1]], [ws[j], ws[j - 1]])
    fw75 = x_at(0.75, +1) - x_at(0.75, -1)
    rhw = x_at(0.50, +1) - w[ipk]
    return fw75, rhw


def rms_vs(y_model, y_data):
    m = (WL >= 450) & (WL <= 700)
    a, b = norm(y_model)[m], norm(y_data)[m]
    return 100 * np.sqrt(np.mean((a - b) ** 2))


print(f"{'sample':7s} {'':>10s} {'FW75':>6s} {'redHW':>6s} {'RMS%':>6s}   "
      "(FW75/redHW in nm; RMS on peak-normalized 450-700)")
RESULTS = []
for S in SAMPLES:
    fw_d, rh_d = widths(S["y"])
    row = dict(lab=S["lab"], data=(fw_d, rh_d))
    print(f"{S['lab']:7s} {'data':>10s} {fw_d:6.1f} {rh_d:6.1f} {'—':>6s}")
    gb, gw = gauss_bins(S["D"], S["sd"])
    variants = {
        "A": ("gauss(raw sd)", ensemble(gb, gw, 1.0)),
        "B": ("per-particle", ensemble(S["bins"], S["w"], 1.0)),
    }
    # C: scan damping scale on the per-particle histogram
    best = (1e9, None, None)
    for s in S_SCAN:
        yC = ensemble(S["bins"], S["w"], s)
        r = rms_vs(yC, S["y"])
        if r < best[0]:
            best = (r, s, yC)
    variants["C"] = (f"+damping s={best[1]:.2f}", best[2])
    yD = np.interp(WL, WL - 2.7, best[2])          # shift model 2.7 nm blue
    variants["D"] = ("+2.7 nm shift", yD)
    row["variants"] = {}
    for key, (lab2, ym) in variants.items():
        fw, rh = widths(ym)
        r = rms_vs(ym, S["y"])
        row["variants"][key] = (fw, rh, r, ym, lab2)
        print(f"{'':7s} {key:>10s} {fw:6.1f} {rh:6.1f} {r:6.2f}   {lab2}")
    RESULTS.append(row)

# ------------------------------------------------------------------- figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25, "axes.axisbelow": True})
fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))

ax = axes[0]
S = SAMPLES[2]; R = RESULTS[2]                     # 31 nm example
ax.plot(WL, norm(S["y"]), "k.", ms=2.5, label="data (CTAC 31 nm)")
for (key, (fw, rh, r, ym, lab2)), c in zip(R["variants"].items(),
                                      ("#0072B2", "#009E73", "#D55E00", "#CC79A7")):
    ax.plot(WL, norm(ym), color=c, lw=1.4, label=f"{key} {lab2}  (RMS {r:.2f}%)")
ax.set_xlim(440, 720); ax.set_xlabel("wavelength (nm)")
ax.set_ylabel("extinction / peak")
ax.set_title("(a) cumulative ingredients, 31 nm sample", fontsize=9)
ax.legend(frameon=False, fontsize=7)

ax = axes[1]
labs = [r["lab"] for r in RESULTS]
x = np.arange(len(RESULTS))
wd = 0.2
for k, (vkey, c) in enumerate(zip("ABCD",
                                  ("#0072B2", "#009E73", "#D55E00", "#CC79A7"))):
    ax.bar(x + (k - 1.5) * wd, [r["variants"][vkey][2] for r in RESULTS],
           wd, color=c, label=vkey)
ax.set_xticks(x); ax.set_xticklabels(labs)
ax.set_ylabel("RMS vs data (% of peak, 450–700 nm)")
ax.set_title("(b) residual closure A → D per sample", fontsize=9)
ax.legend(frameon=False, fontsize=7.5)

ax = axes[2]
for r in RESULTS:
    fw_d = r["data"][0]
    vals = [r["variants"][k][0] for k in "ABCD"]
    ax.plot(range(4), np.array(vals) - fw_d, "o-", ms=4, lw=1,
            label=f"{r['lab']} (data {fw_d:.0f} nm)")
ax.axhline(0, color="0.3", lw=0.8)
ax.set_xticks(range(4)); ax.set_xticklabels(["A", "B", "C", "D"])
ax.set_ylabel("model FW75 − data FW75 (nm)")
ax.set_title("(c) width surplus by variant", fontsize=9)
ax.legend(frameon=False, fontsize=7)

fig.suptitle("Plasmon-width decomposition on CTAC monomers: Gaussian-sd "
             "inflation vs per-particle histogram vs bulk-damping scale",
             fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig19_width_decomposition.png"))
print("\nwrote outputs/fig19_width_decomposition.png")
