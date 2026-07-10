"""Self-test of the Layer 3 inverse fit.

Generate a synthetic 'experimental' spectrum from known parameters (+ noise),
fit it, and check we recover the mean size, polydispersity, and multimer
fractions. Writes outputs/fig5_inverse_fit.png.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib.pyplot as plt

from aunp_speciation.spectra import species_basis, mix
from aunp_speciation.fitting import fit_spectrum

OUT = os.path.join(os.path.dirname(__file__), "..", "outputs")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(0)

# ---- ground truth ----
SPECIES = ("monomer", "dimer", "trimer_linear")
D_true, P_true, GAP = 12.0, 4.0, 1.0
w_true = {"monomer": 0.70, "dimer": 0.15, "trimer_linear": 0.075}
K = {"monomer": 1, "dimer": 2, "trimer_linear": 3}

wl = np.arange(400, 760, 1.0)
basis_true = species_basis(wl, D_true, P_true, GAP, "water", species=SPECIES)
truth = mix(basis_true, w_true)
truth /= truth.max()
noise = 0.005
data = truth + rng.normal(0, noise, truth.shape)

# true gold fractions
kw = np.array([K[s] * w_true[s] for s in SPECIES])
gold_true = {s: kw[i] / kw.sum() for i, s in enumerate(SPECIES)}

# ---- fit ----
print("fitting...")
r = fit_spectrum(wl, data, species=SPECIES, gap_nm=GAP, fit_stride=2)

agg_true = gold_true["dimer"] + gold_true["trimer_linear"]
print(f"\n{'param':<26}{'truth':>10}{'fit':>16}{'verdict':>10}")
print(f"{'diameter (nm)':<26}{D_true:>10.2f}{r.diameter_nm:>10.2f} ± {r.diameter_sd:>4.1f}"
      f"{r.identifiability['diameter']:>10}")
print(f"{'polydispersity (%)':<26}{P_true:>10.2f}{r.pct_poly:>10.2f} ± {r.pct_poly_sd:>4.0f}"
      f"{r.identifiability['pct_poly']:>10}")
print(f"{'AGGREGATED gold frac':<26}{agg_true:>10.3f}"
      f"{r.aggregated_gold:>10.3f} ± {r.aggregated_gold_sd:>4.2f}"
      f"{r.identifiability['aggregated_gold']:>10}")
print(f"{'  monomer gold frac':<26}{gold_true['monomer']:>10.3f}"
      f"{r.gold_fractions['monomer']:>10.3f} ± {r.gold_fractions_sd['monomer']:>4.2f}")
print(f"\n  dimer/trimer split: {r.identifiability['dimer_vs_trimer_split']} "
      f"(collinear red tails — needs the temperature/isosbestic series to resolve)")
print(f"\nresidual RMS = {r.residual_rms:.4f}  (noise = {noise})   success={r.success}")

# ---- figure ----
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 6), height_ratios=[3, 1], sharex=True)
ax1.plot(wl, data, ".", ms=3, color="0.55", label="synthetic data (+noise)")
ax1.plot(r.wavelength, r.model, "-", color="#D55E00", lw=2, label="best fit")
txt = (f"D = {r.diameter_nm:.1f}±{r.diameter_sd:.1f} nm\n"
       f"poly = {r.pct_poly:.1f}% (weakly constrained)\n"
       f"aggregated gold = {r.aggregated_gold*100:.0f}±{r.aggregated_gold_sd*100:.0f}%\n"
       f"(truth {agg_true*100:.0f}%)")
ax1.text(0.97, 0.95, txt, transform=ax1.transAxes, ha="right", va="top",
         fontsize=9, bbox=dict(boxstyle="round", fc="white", ec="0.7"))
ax1.set_ylabel("extinction (norm.)"); ax1.legend(frameon=False)
ax1.set_title("Layer 3 fit: aggregated fraction recovered; size/poly under-determined")
resid = np.interp(r.wavelength, wl, data) - r.model
ax2.axhline(0, color="0.7", lw=1)
ax2.plot(r.wavelength, resid, color="#0072B2", lw=1)
ax2.set_ylabel("residual"); ax2.set_xlabel("wavelength (nm)")
ax2.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig5_inverse_fit.png"), dpi=130)
print("wrote outputs/fig5_inverse_fit.png")
