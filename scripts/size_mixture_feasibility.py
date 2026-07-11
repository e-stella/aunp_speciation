"""Feasibility screen for the size-mixing validation experiment (29 + 42 nm).

Question: if two nominal sizes from the CTAC series are mixed at known ratios,
can the pinned exact-optics basis retrieve the mixing fraction from UV-Vis —
and can a dimer of the small size be told apart from a monomer of the large
size? Everything is computed per unit GOLD VOLUME (so NNLS weights convert
directly to gold fractions), on the exact T-matrix basis cached by
scripts/fit_ctac_validation.py (jc, gamma_S A_surf=0.25, gap 3.5 nm,
polydispersity pinned to TEM).

Outputs:
  - pairwise cosine-similarity matrix (monomer/dimer/trimer x 29.25/42.0 nm)
  - mixing-fraction recovery vs noise level (NNLS, monomer-pair basis)
  - masquerade test: pure 29 nm dimer fitted with a monomer-only basis
  - outputs/fig15_mixture_feasibility.png

Run:  python scripts/size_mixture_feasibility.py   (system env; needs the
experimental/ctac/ctac_exact_basis.npz cache)
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import nnls

from aunp_speciation import dielectric
dielectric.use_gold_model("jc")
from aunp_speciation.spectra import monomer_polydisperse

ROOT = os.path.join(os.path.dirname(__file__), "..")
CACHE = os.path.join(ROOT, "experimental", "ctac", "ctac_exact_basis.npz")
OUT = os.path.join(ROOT, "outputs")
A_SURF = 0.25
WLF = np.arange(420.0, 801.0, 1.0)
N_NOISE = 400
NOISE_LEVELS = (0.001, 0.003, 0.01)     # fraction of peak: 0.1% / 0.3% / 1%

z = np.load(CACHE)
WL_TM, D_TEM, POLY = z["wl"], z["D_tem"], z["poly"]
IDX = {"29": 1, "42": 4}                # D_tem = [7.78, 29.25, 31.45, 37.89, 42.0]


def per_gold(name):
    """Extinction per unit gold volume on the 1 nm grid."""
    sp, tag = name.split("_")
    i = IDX[tag]
    vol = (np.pi / 6.0) * D_TEM[i] ** 3
    if sp == "m":
        return monomer_polydisperse(WLF, D_TEM[i], POLY[i], "water",
                                    None, True, A_SURF) / vol
    k, key = {"d": (2, "dimer"), "t": (3, "trimer_linear")}[sp]
    return CubicSpline(WL_TM, z[f"cube__{key}"][i])(WLF) / (k * vol)


NAMES = ["m_29", "m_42", "d_29", "d_42", "t_29", "t_42"]
B = {n: per_gold(n) for n in NAMES}


def cos(a, b):
    return float(a @ b / np.sqrt((a @ a) * (b @ b)))


print(f"pinned inputs: D = {D_TEM[IDX['29']]:.2f} / {D_TEM[IDX['42']]:.2f} nm, "
      f"poly = {POLY[IDX['29']]:.1f} / {POLY[IDX['42']]:.1f} %  "
      f"(gap 3.5 nm, A_surf={A_SURF}, jc)")
print("\ncosine-similarity matrix (per gold volume, 420-800 nm):")
print("        " + "".join(f"{n:>8s}" for n in NAMES))
CM = np.zeros((len(NAMES), len(NAMES)))
for i, a in enumerate(NAMES):
    row = ""
    for j, b in enumerate(NAMES):
        CM[i, j] = cos(B[a], B[b])
        row += f"{CM[i, j]:8.4f}"
    print(f"{a:8s}" + row)

print("\nkey pairs:")
for a, b, why in [("m_29", "m_42", "the size-mixture experiment"),
                  ("m_42", "d_29", "does a 29-dimer mimic a 42-monomer?"),
                  ("m_42", "d_42", "42 nm speciation (reference)"),
                  ("d_29", "d_42", "cross-size dimers")]:
    c = cos(B[a], B[b])
    print(f"  cos({a},{b}) = {c:.4f}  orth {100*np.sqrt(1-c*c):.1f}%   [{why}]")

# ---------------- mixing-fraction recovery (monomer-pair basis) ------------
rng = np.random.default_rng(7)
b29, b42 = B["m_29"], B["m_42"]
Bmono = np.column_stack([b29, b42])
fracs = np.array([0.0, 0.02, 0.05, 0.10, 0.20, 0.50, 0.80, 1.0])
print("\nmixing-fraction recovery: f = gold fraction of 42 nm in a 29+42 mix")
print(f"{'f_true':>7s}" + "".join(f"  rec@{100*s:.1f}%noise" for s in NOISE_LEVELS))
rec_tbl = np.zeros((len(fracs), len(NOISE_LEVELS), 2))
for i, f in enumerate(fracs):
    y0 = f * b42 + (1 - f) * b29
    scale = y0.max()
    for j, s in enumerate(NOISE_LEVELS):
        rec = np.empty(N_NOISE)
        for k in range(N_NOISE):
            y = y0 + rng.normal(0.0, s * scale, y0.shape)
            w, _ = nnls(Bmono, y)
            rec[k] = w[1] / max(w.sum(), 1e-30)
        rec_tbl[i, j] = rec.mean(), rec.std()
    print(f"{f:7.2f}" + "".join(
        f"   {rec_tbl[i, j, 0]:.3f}±{rec_tbl[i, j, 1]:.3f}"
        for j in range(len(NOISE_LEVELS))))

# detection limit: smallest f whose recovered distribution clears 2 sd of f=0
print("\n2σ detection limits for the 42 nm component:")
for j, s in enumerate(NOISE_LEVELS):
    floor = rec_tbl[0, j, 0] + 2 * rec_tbl[0, j, 1]
    above = [f for i, f in enumerate(fracs[1:], 1)
             if rec_tbl[i, j, 0] - 2 * rec_tbl[i, j, 1] > floor]
    lim = f"~{min(above):.2f}" if above else ">0.8"
    print(f"  noise {100*s:.1f}% of peak: f_42 detectable above {lim} "
          f"(f=0 baseline reads {rec_tbl[0, j, 0]:.3f}±{rec_tbl[0, j, 1]:.3f})")

# ---------------- masquerade test: pure 29-dimer vs monomer basis ----------
print("\nmasquerade test: pure 29 nm DIMER population fitted with monomers only")
y0 = B["d_29"]
scale = y0.max()
for s in (0.003,):
    w, _ = nnls(Bmono, y0)
    resid = Bmono @ w - y0
    rms = np.sqrt(np.mean(resid ** 2)) / scale
    print(f"  best monomer-only fit: {100*w[0]/w.sum():.1f}% '29-mono' + "
          f"{100*w[1]/w.sum():.1f}% '42-mono', RMS {100*rms:.2f}% of peak "
          f"(vs noise {100*s:.1f}%) -> "
          f"{'DISTINGUISHABLE (structured residual)' if rms > 3*s else 'NOT distinguishable at this noise'}")
# and with the dimer column available, does NNLS find it?
B3 = np.column_stack([b29, b42, B["d_29"]])
w, _ = nnls(B3, y0 + rng.normal(0, 0.003 * scale, y0.shape))
print(f"  with d_29 in the basis: weights (m29/m42/d29 gold) = "
      f"{w[0]/w.sum():.3f}/{w[1]/w.sum():.3f}/{w[2]/w.sum():.3f}")

# ------------------------------------------------------------------- figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25, "axes.axisbelow": True})
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
ax = axes[0]
for n, c in zip(NAMES, ("#0072B2", "#D55E00", "#56B4E9", "#E69F00",
                        "#009E73", "#CC79A7")):
    ax.plot(WLF, B[n] / B[n].max(), color=c, lw=1.5,
            ls="-" if n.startswith("m") else "--" if n.startswith("d") else ":",
            label=n.replace("_", " ") + " nm")
ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("extinction / peak (per gold)")
ax.set_title("(a) per-gold basis shapes (exact optics, pinned)", fontsize=9)
ax.legend(frameon=False, fontsize=7.5, ncol=2)

ax = axes[1]
im = ax.imshow(CM, vmin=0.9, vmax=1.0, cmap="viridis")
ax.set_xticks(range(len(NAMES))); ax.set_yticks(range(len(NAMES)))
ax.set_xticklabels([n.replace("_", "") for n in NAMES], fontsize=8)
ax.set_yticklabels([n.replace("_", "") for n in NAMES], fontsize=8)
for i in range(len(NAMES)):
    for j in range(len(NAMES)):
        ax.text(j, i, f"{CM[i, j]:.3f}", ha="center", va="center", fontsize=6.5,
                color="w" if CM[i, j] < 0.97 else "k")
ax.set_title("(b) cosine similarity (1 = indistinguishable)", fontsize=9)
ax.grid(False)

ax = axes[2]
cols = ("#0072B2", "#D55E00", "#009E73")
for j, (s, c) in enumerate(zip(NOISE_LEVELS, cols)):
    ax.errorbar(fracs, rec_tbl[:, j, 0], yerr=2 * rec_tbl[:, j, 1],
                fmt="o-", ms=3.5, lw=1.2, capsize=2.5, color=c,
                label=f"noise {100*s:.1f}% of peak")
ax.plot([0, 1], [0, 1], "k--", lw=0.8)
ax.set_xlabel("true gold fraction of 42 nm"); ax.set_ylabel("recovered fraction")
ax.set_title("(c) 29+42 nm mixture retrieval (±2σ)", fontsize=9)
ax.legend(frameon=False, fontsize=7.5)
fig.suptitle("Size-mixture feasibility: 29.25 + 42.0 nm CTAC, exact basis, "
             "TEM-pinned, per-gold NNLS", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig15_mixture_feasibility.png"))
print("\nwrote outputs/fig15_mixture_feasibility.png")
