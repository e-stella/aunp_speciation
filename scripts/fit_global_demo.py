"""Self-test of the global multi-temperature fit.

Generate a synthetic temperature series (same sample: shared size/poly/gap,
conserved gold, van 't Hoff populations), add noise, fit jointly, and check we
recover size, polydispersity, thermodynamics, and the per-temperature gold
fractions — the degeneracies that a single spectrum could not resolve.

Writes outputs/fig7_global_fit.png
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib.pyplot as plt

from aunp_speciation.spectra import species_basis
from aunp_speciation.equilibrium import association_constants, gold_fractions
from aunp_speciation.fit_global import fit_temperature_series

OUT = os.path.join(os.path.dirname(__file__), "..", "outputs")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(1)

SPECIES = ("monomer", "dimer", "trimer_linear")
K = {"monomer": 1, "dimer": 2, "trimer_linear": 3}
D_true, P_true, GAP = 12.0, 4.0, 1.0
C_TOT = 0.01   # dimensionless total-gold scale (mostly-monomer regime)
# ground-truth thermodynamics (exothermic association; heating -> more monomer)
TH = dict(dH2=-25.0, dS2=-0.06, dH3=-45.0, dS3=-0.14)
temps_C = np.array([5.0, 20.0, 35.0, 50.0, 65.0])
temps_K = temps_C + 273.15
wl = np.arange(460, 700, 4.0)
noise = 0.006

# ---- synthesize the T-series (CDA backend for speed) ----
basis = species_basis(wl, D_true, P_true, GAP, "water", species=SPECIES,
                      backend="cda", n_sizes=7)
B = np.column_stack([basis[s] for s in SPECIES])
kvec = np.array([K[s] for s in SPECIES], dtype=float)


def canon(s):
    return "trimer" if s.startswith("trimer") else s


true_frac = {s: [] for s in SPECIES}
data = []
A_true = 1.0
for T in temps_K:
    K2, K3 = association_constants(T, **TH)
    gf = gold_fractions(C_TOT, K2, K3)
    f = np.array([gf[canon(s)] for s in SPECIES])
    for j, s in enumerate(SPECIES):
        true_frac[s].append(f[j])
    data.append(A_true * (B @ (f / kvec)))
data = np.array(data)
data = data / data.max() + rng.normal(0, noise, data.shape)
true_frac = {s: np.array(v) for s, v in true_frac.items()}
agg_true = true_frac["dimer"] + true_frac["trimer_linear"]

# ---- global fit ----
print("global fitting (this runs the CDA basis build each iteration)...")
r = fit_temperature_series(temps_K, wl, data, species=SPECIES, gap_nm=GAP,
                           backend="cda", n_sizes=5, fit_stride=2, max_nfev=120, C_tot=C_TOT)

print(f"\n{'param':<16}{'truth':>10}{'fit':>16}")
print(f"{'diameter (nm)':<16}{D_true:>10.2f}{r.diameter_nm:>10.2f} ± {r.param_sd['diameter']:.2f}")
print(f"{'poly (%)':<16}{P_true:>10.2f}{r.pct_poly:>10.2f} ± {r.param_sd['pct_poly']:.2f}")
print(f"{'dH2 (kJ/mol)':<16}{TH['dH2']:>10.1f}{r.dH2:>10.1f} ± {r.param_sd['dH2']:.1f}")
print(f"{'dS2 (kJ/mol/K)':<16}{TH['dS2']:>10.3f}{r.dS2:>10.3f} ± {r.param_sd['dS2']:.3f}")
agg_fit = r.gold_fractions["dimer"] + r.gold_fractions["trimer_linear"]
print("\naggregated gold fraction vs T:")
for i, Tc in enumerate(temps_C):
    print(f"  {Tc:4.0f} C : truth {agg_true[i]:.3f}   fit {agg_fit[i]:.3f}")
print(f"\nresidual RMS = {r.residual_rms:.4f}  (noise {noise})  success={r.success}")

# ---- figure ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
cmap = plt.cm.coolwarm(np.linspace(0, 1, len(temps_K)))
for i, Tc in enumerate(temps_C):
    ax1.plot(wl, data[i], ".", ms=2.5, color=cmap[i])
    ax1.plot(r.wavelength, r.models[i], "-", color=cmap[i], lw=1.6,
             label=f"{Tc:.0f} °C")
ax1.set_xlabel("wavelength (nm)"); ax1.set_ylabel("extinction (norm.)")
ax1.set_title("Global fit to the temperature series (points=data, lines=fit)")
ax1.legend(frameon=False, fontsize=8, ncol=2)

ax2.plot(temps_C, agg_true, "o-", color="0.3", label="truth (aggregated)")
ax2.plot(temps_C, agg_fit, "s--", color="#D55E00", label="fit (aggregated)")
for s, c in [("monomer", "#0072B2"), ("dimer", "#D55E00"), ("trimer_linear", "#009E73")]:
    ax2.plot(temps_C, true_frac[s], "-", color=c, alpha=0.4)
    ax2.plot(temps_C, r.gold_fractions[s], "x", color=c, label=f"fit {s}")
ax2.set_xlabel("temperature (°C)"); ax2.set_ylabel("gold fraction")
ax2.set_title(f"Recovered speciation vs T  (D={r.diameter_nm:.1f} nm, poly={r.pct_poly:.1f}%)")
ax2.legend(frameon=False, fontsize=8, ncol=2)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig7_global_fit.png"), dpi=130)
print("wrote outputs/fig7_global_fit.png")
