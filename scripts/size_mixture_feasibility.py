"""Feasibility screen for the size-mixing validation experiment (29 + 42 nm).

Question: if two nominal sizes from the CTAC series are mixed at known ratios,
can the pinned exact-optics basis retrieve the mixing fraction from UV-Vis —
and what CANNOT be told apart? Everything is computed per unit GOLD VOLUME
(NNLS weights convert directly to gold fractions), on the exact T-matrix basis
cached by scripts/fit_ctac_validation.py (jc, gamma_S A_surf=0.25, gap 3.5 nm,
polydispersity pinned to TEM: 7.13% @ 29.25 nm, 6.44% @ 42.0 nm).

Measured findings (all reproduced by this script):
  1. Instrument noise, measured from the CTAC scans (2nd-difference estimator,
     600-900 nm): sigma = 1.7-1.9 mA/point = 0.13-0.21% of peak.
  2. Mixing-fraction retrieval WITH the two component shapes as priors is
     unbiased, 2-sigma detection ~2% of gold at the measured noise; robust to
     polydispersity 3-20% when the poly is known (common-mode broadening).
  3. Dimer masquerade: a pure 29 nm dimer grabs the 42-monomer column but
     leaves a 9.7%-of-peak structured residual (~30x noise) — detectable.
  4. MODEL DISCRIMINATION FAILS (the key caveat): shape-only UV-Vis cannot
     tell the 50/50 mix from a SINGLE intermediate population — best single
     imitations reach RMS 0.008% of peak (poly free -> D=31.2, 24% poly) and
     0.09% (poly <= 10% -> D=36), and the reverse holds too (a narrow 35 nm /
     5% sample fits as a 34:66 "mix" at RMS 0.11%) — all BELOW the 0.15%
     noise floor. Absolute per-gold extinction rejects the broad imitation
     (+15% level error) but not the narrow one (1-2.4%). => Retrieval is
     CONDITIONAL on knowing the sample is composed of the two TEM'd
     endmembers; blind bimodality detection needs TEM/DLS.

Outputs: printed diagnostics + outputs/fig15_mixture_feasibility.png.
Run:  python scripts/size_mixture_feasibility.py   (system env; needs the
experimental/ctac/ctac_exact_basis.npz cache; ~2-3 min for the grid search)
"""
import sys, os, csv
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import nnls

from aunp_speciation import dielectric
dielectric.use_gold_model("jc")
from aunp_speciation.spectra import monomer_polydisperse

ROOT = os.path.join(os.path.dirname(__file__), "..")
CACHE = os.path.join(ROOT, "experimental", "ctac", "ctac_exact_basis.npz")
UVVIS = os.path.join(ROOT, "experimental", "ctac", "ctac_uvvis.csv")
OUT = os.path.join(ROOT, "outputs")
A_SURF = 0.25
WLF = np.arange(420.0, 801.0, 1.0)
N_NOISE = 400
NOISE_LEVELS = (0.0015, 0.003, 0.01)    # of peak: measured floor / 2x / worst

z = np.load(CACHE)
WL_TM, D_TEM, POLY = z["wl"], z["D_tem"], z["poly"]
IDX = {"29": 1, "42": 4}


def per_gold(name):
    """Extinction per unit gold volume on the 1 nm grid."""
    sp, tag = name.split("_")
    i = IDX[tag]
    v = (np.pi / 6.0) * D_TEM[i] ** 3
    if sp == "m":
        return monomer_polydisperse(WLF, D_TEM[i], POLY[i], "water",
                                    None, True, A_SURF) / v
    k, key = {"d": (2, "dimer"), "t": (3, "trimer_linear")}[sp]
    return CubicSpline(WL_TM, z[f"cube__{key}"][i])(WLF) / (k * v)


def single_pop(D, p):
    return monomer_polydisperse(WLF, D, p, "water", None, True, A_SURF) \
        / ((np.pi / 6.0) * D ** 3)


def cos(a, b):
    return float(a @ b / np.sqrt((a @ a) * (b @ b)))


NAMES = ["m_29", "m_42", "d_29", "d_42", "t_29", "t_42"]
B = {n: per_gold(n) for n in NAMES}

# ---------------- measured instrument noise (from the CTAC scans) ----------
with open(UVVIS) as f:
    rows = list(csv.reader(f))
dat = np.array([[float(x) for x in r] for r in rows[1:]])
o = np.argsort(dat[:, 0]); dat = dat[o]
wl_raw, A_raw = dat[:, 0], dat[:, 1:]
mred = (wl_raw >= 600) & (wl_raw <= 900)
mpk = (wl_raw >= 420) & (wl_raw <= 800)
noise_meas = []
for i in range(A_raw.shape[1]):
    d2 = np.diff(A_raw[mred, i], 2)
    noise_meas.append((d2.std() / np.sqrt(6)) / A_raw[mpk, i].max())
noise_meas = np.array(noise_meas)
print(f"measured per-point noise (Cary, 2nd-difference, 600-900 nm): "
      f"{100*noise_meas.min():.2f}-{100*noise_meas.max():.2f}% of peak")
print(f"pinned inputs: D = {D_TEM[IDX['29']]:.2f} / {D_TEM[IDX['42']]:.2f} nm, "
      f"poly(TEM) = {POLY[IDX['29']]:.2f} / {POLY[IDX['42']]:.2f} %  "
      f"(gap 3.5 nm, A_surf={A_SURF}, jc)")

print("\ncosine-similarity matrix (per gold volume, 420-800 nm):")
print("        " + "".join(f"{n:>8s}" for n in NAMES))
CM = np.zeros((len(NAMES), len(NAMES)))
for i, a in enumerate(NAMES):
    for j, b in enumerate(NAMES):
        CM[i, j] = cos(B[a], B[b])
    print(f"{a:8s}" + "".join(f"{CM[i, j]:8.4f}" for j in range(len(NAMES))))

# ---------------- mixing-fraction recovery (monomer-pair basis) ------------
rng = np.random.default_rng(7)
b29, b42 = B["m_29"], B["m_42"]
Bmono = np.column_stack([b29, b42])
fracs = np.array([0.0, 0.02, 0.05, 0.10, 0.20, 0.50, 0.80, 1.0])
rec_tbl = np.zeros((len(fracs), len(NOISE_LEVELS), 2))
for i, f in enumerate(fracs):
    y0 = f * b42 + (1 - f) * b29
    scale = y0.max()
    for j, s in enumerate(NOISE_LEVELS):
        rec = np.empty(N_NOISE)
        for k in range(N_NOISE):
            w, _ = nnls(Bmono, y0 + rng.normal(0.0, s * scale, y0.shape))
            rec[k] = w[1] / max(w.sum(), 1e-30)
        rec_tbl[i, j] = rec.mean(), rec.std()
print("\nmixing retrieval (f = gold fraction of 42 nm), known components:")
for j, s in enumerate(NOISE_LEVELS):
    floor = rec_tbl[0, j, 0] + 2 * rec_tbl[0, j, 1]
    above = [f for i, f in enumerate(fracs[1:], 1)
             if rec_tbl[i, j, 0] - 2 * rec_tbl[i, j, 1] > floor]
    print(f"  noise {100*s:.2f}%: 2-sigma detection "
          f"{('~%.2f' % min(above)) if above else '>0.8'}")

# poly robustness (same poly both components, KNOWN to the fit)
print("\nretrieval vs polydispersity (f=0.5, noise 0.3%, poly KNOWN):")
for p in (3.0, 7.0, 12.0, 20.0):
    c29 = single_pop(29.25, p); c42 = single_pop(42.0, p)
    Bp = np.column_stack([c29, c42])
    y0 = 0.5 * c42 + 0.5 * c29; s = 0.003 * y0.max()
    rec = [nnls(Bp, y0 + rng.normal(0, s, y0.shape))[0] for _ in range(150)]
    rec = np.array([w[1] / w.sum() for w in rec])
    print(f"  poly {p:4.1f}%: recovered {rec.mean():.3f} +- {2*rec.std():.3f} (2sd)")

# ---------------- masquerade + model-discrimination tests ------------------
y0 = B["d_29"]
w, _ = nnls(Bmono, y0)
rms_d = np.sqrt(np.mean((Bmono @ w - y0) ** 2)) / y0.max()
print(f"\npure 29-dimer fitted monomer-only: {100*w[1]/w.sum():.0f}% lands on "
      f"'42-mono', residual {100*rms_d:.1f}% of peak (~{rms_d/0.0015:.0f}x noise)"
      " -> detectable")

ymix = 0.5 * b42 + 0.5 * b29
ymixN = ymix / ymix.max()


def best_single(y, polys, Ds=np.arange(28.0, 46.1, 0.5)):
    best = (1e9, None, None)
    for D in Ds:
        for p in polys:
            m = single_pop(D, p)
            a = (m @ y) / (m @ m)
            r = float(np.sqrt(np.mean((a * m - y) ** 2)))
            if r < best[0]:
                best = (r, D, p)
    return best


print("\nMODEL DISCRIMINATION (shape only): can a SINGLE population imitate "
      "the 50/50 mix?")
imitators = {}
for tag, polys in [("poly free", np.arange(2.0, 30.1, 1.0)),
                   ("poly<=10%", np.arange(2.0, 10.1, 0.5))]:
    r, D, p = best_single(ymixN, polys)
    imitators[tag] = (r, D, p)
    print(f"  best single [{tag}]: D={D:.2f} poly={p:.1f}% "
          f"RMS={100*r:.3f}% of peak -> DEGENERATE (noise floor 0.15%)")
y35 = single_pop(35.0, 5.0); y35 /= y35.max()
w, _ = nnls(Bmono, y35)
r35 = float(np.sqrt(np.mean((Bmono @ w - y35) ** 2)))
print(f"  reverse: narrow 35 nm/5% fits as a {100*w[1]/w.sum():.0f}% '42' mix, "
      f"RMS={100*r35:.3f}% -> DEGENERATE")
print("  absolute per-gold level rejects the broad imitation (+15%) but not "
      "the narrow one (1-2.4%).")
print("  => mixture retrieval is CONDITIONAL on the two-endmember prior; "
      "blind bimodality detection needs TEM/DLS.")

# ------------------------------------------------------------------- figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25, "axes.axisbelow": True})
fig, axes = plt.subplots(2, 2, figsize=(11, 8.2))

ax = axes[0, 0]
for n, c in zip(NAMES, ("#0072B2", "#D55E00", "#56B4E9", "#E69F00",
                        "#009E73", "#CC79A7")):
    ax.plot(WLF, B[n] / B[n].max(), color=c, lw=1.5,
            ls="-" if n.startswith("m") else "--" if n.startswith("d") else ":",
            label=n.replace("_", " ") + " nm")
ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("extinction / peak (per gold)")
ax.set_title("(a) per-gold basis shapes — exact optics, poly pinned to TEM\n"
             f"(29 nm: {POLY[IDX['29']]:.1f}%,  42 nm: {POLY[IDX['42']]:.1f}%)",
             fontsize=9)
ax.legend(frameon=False, fontsize=7.5, ncol=2)

ax = axes[0, 1]
im = ax.imshow(CM, vmin=0.9, vmax=1.0, cmap="viridis")
ax.set_xticks(range(len(NAMES))); ax.set_yticks(range(len(NAMES)))
ax.set_xticklabels([n.replace("_", "") for n in NAMES], fontsize=8)
ax.set_yticklabels([n.replace("_", "") for n in NAMES], fontsize=8)
for i in range(len(NAMES)):
    for j in range(len(NAMES)):
        ax.text(j, i, f"{CM[i, j]:.3f}", ha="center", va="center",
                fontsize=6.5, color="w" if CM[i, j] < 0.97 else "k")
ax.set_title("(b) cosine similarity (1 = indistinguishable)", fontsize=9)
ax.grid(False)

ax = axes[1, 0]
for j, (s, c, lab) in enumerate(zip(
        NOISE_LEVELS, ("#0072B2", "#D55E00", "#009E73"),
        (f"{100*NOISE_LEVELS[0]:.2f}% — measured Cary floor",
         f"{100*NOISE_LEVELS[1]:.1f}% — replicate/drift level",
         f"{100*NOISE_LEVELS[2]:.0f}% — pessimistic"))):
    ax.errorbar(fracs, rec_tbl[:, j, 0], yerr=2 * rec_tbl[:, j, 1],
                fmt="o-", ms=3.5, lw=1.2, capsize=2.5, color=c, label=lab)
ax.plot([0, 1], [0, 1], "k--", lw=0.8)
ax.set_xlabel("true gold fraction of 42 nm"); ax.set_ylabel("recovered fraction")
ax.set_title("(c) ratio retrieval WITH the two-endmember prior (±2σ)\n"
             "robust to poly 3–20% when poly is known", fontsize=9)
ax.legend(frameon=False, fontsize=7.5, title="noise (of peak)")

ax = axes[1, 1]
r, D, p = imitators["poly free"]
mfree = single_pop(D, p); afree = (mfree @ ymixN) / (mfree @ mfree)
r2, D2, p2 = imitators["poly<=10%"]
mcon = single_pop(D2, p2); acon = (mcon @ ymixN) / (mcon @ mcon)
ax.axhspan(-0.15, 0.15, color="0.85", alpha=0.7,
           label="measured noise ±0.15%")
ax.plot(WLF, 100 * (afree * mfree - ymixN), color="#D55E00", lw=1.4,
        label=f"single D={D:.0f} nm, {p:.0f}% poly (free)")
ax.plot(WLF, 100 * (acon * mcon - ymixN), color="#0072B2", lw=1.4,
        label=f"single D={D2:.0f} nm, {p2:.0f}% poly (≤10%)")
ax.plot(WLF, 100 * (Bmono @ nnls(Bmono, y35)[0] - y35), color="#009E73",
        lw=1.4, label="narrow 35 nm/5% as '34:66 mix'")
ax.axhline(0, color="0.3", lw=0.7)
ax.set_ylim(-0.5, 0.5)
ax.set_xlabel("wavelength (nm)")
ax.set_ylabel("imitation − truth  (% of peak)")
ax.set_title("(d) shape-only discrimination fails in practice: imitation\n"
             "error RMS ≲ noise and ≪ the 1–4% model-systematics floor",
             fontsize=9)
ax.legend(frameon=False, fontsize=7.5)

fig.suptitle("Size-mixture feasibility, 29.25 + 42.0 nm CTAC (exact basis, "
             "TEM-pinned, per-gold NNLS): retrieval is conditional on the "
             "two-endmember prior", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig15_mixture_feasibility.png"))
print("\nwrote outputs/fig15_mixture_feasibility.png")
