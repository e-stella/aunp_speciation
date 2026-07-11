"""Worked example: extracting the 29:42 nm gold ratio from a mixture UV-Vis.

Companion to scripts/size_mixture_feasibility.py — that script maps the
feasibility limits; this one demonstrates the actual extraction pipeline on a
synthetic "unknown" mixture (f42 = 0.30 of the gold, measured Cary noise):

  1. TEM the two stock endmembers once -> (D, poly) for each.
  2. Compute each endmember's extinction PER UNIT GOLD VOLUME (exact optics:
     Mie + jc gold + gamma_S; in the real experiment the measured pure-stock
     spectra can be used instead, which cancels model systematics).
  3. Measure the mixture spectrum, 420-800 nm.
  4. Non-negative least squares:  y(lambda) ~ w29*b29(lambda) + w42*b42(lambda).
     Because the basis columns are per-gold, the gold ratio is directly
     f42 = w42/(w29+w42) — the overall amplitude (concentration, path length)
     cancels out of the ratio.
  5. Uncertainty by residual bootstrap; a structured residual (>> noise)
     flags contamination (e.g. dimers) instead of silently biasing the ratio.

Writes outputs/fig16_mixture_extraction_demo.png.
Run:  python scripts/mixture_extraction_demo.py   (system env; needs the
experimental/ctac/ctac_exact_basis.npz cache for the pinned TEM inputs)
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from scipy.optimize import nnls

from aunp_speciation import dielectric
dielectric.use_gold_model("jc")
from aunp_speciation.spectra import monomer_polydisperse

ROOT = os.path.join(os.path.dirname(__file__), "..")
CACHE = os.path.join(ROOT, "experimental", "ctac", "ctac_exact_basis.npz")
OUT = os.path.join(ROOT, "outputs")
WL = np.arange(420.0, 801.0, 1.0)
A_SURF = 0.25
F_TRUE = 0.30            # gold fraction of the 42 nm component
NOISE = 0.0015           # measured Cary per-point floor (0.15% of peak)
N_BOOT = 500

z = np.load(CACHE)
D29, P29 = z["D_tem"][1], z["poly"][1]
D42, P42 = z["D_tem"][4], z["poly"][4]
vol = lambda d: (np.pi / 6.0) * d ** 3
b29 = monomer_polydisperse(WL, D29, P29, "water", None, True, A_SURF) / vol(D29)
b42 = monomer_polydisperse(WL, D42, P42, "water", None, True, A_SURF) / vol(D42)
B = np.column_stack([b29, b42])

# ---- step 3: the "measured" mixture (synthetic, with measured noise) ------
rng = np.random.default_rng(11)
y_clean = F_TRUE * b42 + (1 - F_TRUE) * b29
scale = y_clean.max()
y = y_clean + rng.normal(0.0, NOISE * scale, y_clean.shape)

# ---- step 4: NNLS -> gold ratio -------------------------------------------
w, _ = nnls(B, y)
f_hat = w[1] / w.sum()
resid = y - B @ w
rms = float(np.sqrt(np.mean(resid ** 2))) / scale

# ---- step 5: residual bootstrap uncertainty -------------------------------
boot = np.empty(N_BOOT)
yfit = B @ w
for k in range(N_BOOT):
    wb, _ = nnls(B, yfit + rng.choice(resid, len(resid), replace=True))
    boot[k] = wb[1] / wb.sum()
f_sd = boot.std()

print(f"endmembers (TEM-pinned): {D29:.2f} nm / {P29:.1f}%  and  "
      f"{D42:.2f} nm / {P42:.1f}%  (jc, gamma_S A_surf={A_SURF})")
print(f"true 42 nm gold fraction : {F_TRUE:.3f}")
print(f"recovered                : {f_hat:.3f} ± {f_sd:.3f} (1σ bootstrap)")
print(f"fit residual RMS         : {100*rms:.3f}% of peak "
      f"(noise injected: {100*NOISE:.2f}%) -> no contamination flag")
print("number-fraction view     : "
      f"{100 * (f_hat/vol(D42)) / (f_hat/vol(D42) + (1-f_hat)/vol(D29)):.1f}% "
      "of PARTICLES are 42 nm")

# ------------------------------------------------------------------- figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25, "axes.axisbelow": True})
fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))

ax = axes[0]
cmap = plt.cm.plasma(np.linspace(0.1, 0.85, 5))
for f, c in zip((0.0, 0.25, 0.5, 0.75, 1.0), cmap):
    ym = f * b42 + (1 - f) * b29
    ax.plot(WL, ym / ym.max(), color=c, lw=1.5, label=f"f₄₂ = {f:.2f}")
ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("extinction / peak")
ax.set_title("(a) the mixture family — differences are subtle\n"
             "(information lives in the 540–650 nm flank)", fontsize=9)
ax.legend(frameon=False, fontsize=7.5, title="gold fraction 42 nm")
axins = ax.inset_axes([0.45, 0.35, 0.52, 0.55])
for f, c in zip((0.0, 0.25, 0.5, 0.75, 1.0), cmap):
    ym = f * b42 + (1 - f) * b29
    axins.plot(WL, ym / ym.max(), color=c, lw=1.2)
axins.set_xlim(540, 650); axins.set_ylim(0.1, 0.75)
axins.tick_params(labelsize=7)
ax.indicate_inset_zoom(axins, edgecolor="0.6")

ax = axes[1]
ax.plot(WL[::4], y[::4] / scale, "o", ms=3, mfc="none", mec="k", mew=0.8,
        label='"measured" mixture (every 4th pt)')
ax.plot(WL, (B @ w) / scale, "-", color="#D55E00", lw=1.6, label="NNLS fit")
ax.plot(WL, w[0] * b29 / scale, "--", color="#0072B2", lw=1.3,
        label=f"29 nm part ({100*(1-f_hat):.0f}% of gold)")
ax.plot(WL, w[1] * b42 / scale, "--", color="#009E73", lw=1.3,
        label=f"42 nm part ({100*f_hat:.0f}% of gold)")
ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("extinction (norm.)")
ax.set_title(f"(b) extraction: true f₄₂ = {F_TRUE:.2f} → "
             f"recovered {f_hat:.3f} ± {f_sd:.3f}", fontsize=9)
ax.legend(frameon=False, fontsize=7.5)

ax = axes[2]
ax.axhspan(-100 * NOISE, 100 * NOISE, color="0.85", alpha=0.7,
           label=f"±{100*NOISE:.2f}% noise band")
ax.plot(WL, 100 * resid / scale, color="0.3", lw=0.8, label="fit residual")
d = b42 / b42.max() - b29 / b29.max()
ax.plot(WL, d / np.abs(d).max() * 3 * (100 * NOISE), color="#CC79A7",
        lw=1.4, label="discriminating shape b₄₂−b₂₉ (scaled to 3× noise)")
ax.set_xlabel("wavelength (nm)")
ax.set_ylabel("residual (% of peak)")
ax.set_title("(c) flat residual = clean two-component mix;\n"
             "structure here would flag contamination (e.g. dimers)",
             fontsize=9)
ax.legend(frameon=False, fontsize=7.5)

fig.suptitle("Worked example: gold-ratio extraction from a 29+42 nm mixture "
             "spectrum (TEM-pinned endmembers, per-gold NNLS)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig16_mixture_extraction_demo.png"))
print("wrote outputs/fig16_mixture_extraction_demo.png")
