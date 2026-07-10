"""Generate the three key demonstration figures.

  fig1_species_spectra.png  - per-particle spectra of monomer/dimer/trimers
  fig2_broadening.png       - monomer-only (4% poly) vs realistic speciation mix
  fig3_isosbestic.png       - temperature series showing an isosbestic point
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib.pyplot as plt

from aunp_speciation.spectra import monomer_polydisperse, species_basis, mix
from aunp_speciation.clusters import species_spectrum
from aunp_speciation.equilibrium import association_constants, solve_populations

OUT = os.path.join(os.path.dirname(__file__), "..", "outputs")
os.makedirs(OUT, exist_ok=True)

D = 12.0          # nm nominal diameter
POLY = 4.0        # % polydispersity
GAP = 1.0         # nm ligand/citrate surface gap
MED = "water"
wl = np.arange(400, 780, 1.0)

# colorblind-friendly palette (Okabe-Ito)
C = {"mono": "#0072B2", "dimer": "#D55E00", "triL": "#009E73",
     "triT": "#CC79A7", "mix": "#000000", "poly": "#0072B2"}
plt.rcParams.update({"figure.dpi": 130, "font.size": 11, "axes.grid": True,
                     "grid.alpha": 0.25, "axes.axisbelow": True})


def fwhm(x, y):
    ymax = y.max(); half = ymax / 2
    idx = np.where(y >= half)[0]
    return x[idx[-1]] - x[idx[0]]


def red_tail_index(x, y):
    """Red-side intensity (560-700 nm) relative to peak height — the signature
    that discriminates speciation from polydispersity (peak barely moves)."""
    band = (x >= 560) & (x <= 700)
    return np.trapezoid(y[band], x[band]) / y.max()


# ---------------- Figure 1: species spectra (per particle) ----------------
mono = species_spectrum("monomer", D, wl, n_medium=MED)
dim = species_spectrum("dimer", D, wl, GAP, MED) / 2
triL = species_spectrum("trimer_linear", D, wl, GAP, MED) / 3
triT = species_spectrum("trimer_triangular", D, wl, GAP, MED) / 3

fig, ax = plt.subplots(figsize=(7, 4.5))
for y, lab, c in [(mono, "monomer", C["mono"]), (dim, "dimer", C["dimer"]),
                  (triL, "trimer (linear)", C["triL"]),
                  (triT, "trimer (triangular)", C["triT"])]:
    ax.plot(wl, y / mono.max(), label=lab, color=c, lw=2)
ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("extinction per particle (norm.)")
ax.set_title(f"{D:.0f} nm Au, gap {GAP:.0f} nm — coupling adds red-side intensity")
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig1_species_spectra.png")); plt.close(fig)

# ---------------- Figure 2: broadening from speciation ----------------
basis = species_basis(wl, D, POLY, GAP, MED,
                      species=("monomer", "dimer", "trimer_linear"))
mono_only = basis["monomer"]
# realistic mixture: number concentrations (most particles still monomers)
weights = {"monomer": 0.70, "dimer": 0.15, "trimer_linear": 0.075}
mixture = mix(basis, weights)

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(wl, mono_only / mono_only.max(), color=C["poly"], lw=2,
        label=f"monomer only, {POLY:.0f}% poly  (red-tail idx {red_tail_index(wl, mono_only):.1f} nm)")
ax.plot(wl, mixture / mixture.max(), color=C["mix"], lw=2.2,
        label=f"+ 15% dimer, 7.5% trimer  (red-tail idx {red_tail_index(wl, mixture):.1f} nm)")
ax.fill_between(wl, mono_only / mono_only.max(), mixture / mixture.max(),
                where=(mixture / mixture.max() >= mono_only / mono_only.max()),
                color=C["dimer"], alpha=0.15, label="excess (red tail + broadening)")
ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("extinction (norm.)")
ax.set_title("Speciation broadens & red-tails the spectrum beyond polydispersity")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig2_broadening.png")); plt.close(fig)

# ---------------- Figure 3: temperature series -> isosbestic point ----------------
C_tot = 1.0
temps_C = np.array([5, 15, 25, 35, 45, 55, 65])
cmap = plt.cm.plasma(np.linspace(0.1, 0.85, len(temps_C)))

# per-cluster, polydispersity-averaged basis (per gold normalization below)
b_mono = basis["monomer"]                    # already per-particle, poly-averaged
b_dim = basis["dimer"]
b_tri = basis["trimer_linear"]

fig, ax = plt.subplots(figsize=(7, 4.5))
curves = []
for T_C, col in zip(temps_C, cmap):
    T = T_C + 273.15
    K2, K3 = association_constants(T)
    c1, c2, c3 = solve_populations(C_tot, K2, K3)
    # extinction per unit total gold (gold is conserved -> isosbestic)
    E = (c1 * b_mono + c2 * b_dim + c3 * b_tri) / C_tot
    curves.append(E)
    ax.plot(wl, E, color=col, lw=1.8, label=f"{T_C} °C")
curves = np.array(curves)
# find isosbestic: wavelength of minimum spread across temperatures
spread = curves.std(axis=0)
iso = wl[np.argmin(spread[(wl > 500) & (wl < 650)]) + np.where(wl > 500)[0][0]]
ax.axvline(iso, color="0.4", ls="--", lw=1)
ax.annotate(f"isosbestic ≈ {iso:.0f} nm", xy=(iso, curves[:, np.argmin(np.abs(wl-iso))].mean()),
            xytext=(iso + 25, ax.get_ylim()[1]*0.7), fontsize=9,
            arrowprops=dict(arrowstyle="->", color="0.4"))
ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("extinction per total gold (a.u.)")
ax.set_title("Temperature shifts the equilibrium → isosbestic point")
ax.legend(frameon=False, ncol=2, fontsize=9, title="temperature")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig3_isosbestic.png")); plt.close(fig)

# ---------------- Figure 4: gap sensitivity (the dominant lever) ----------------
fig, ax = plt.subplots(figsize=(7, 4.5))
gaps = [5.0, 2.0, 1.0, 0.5]
gcolors = plt.cm.viridis(np.linspace(0.15, 0.8, len(gaps)))
ax.plot(wl, mono / mono.max(), color="0.4", lw=2, ls="--", label="monomer")
for g, col in zip(gaps, gcolors):
    d = species_spectrum("dimer", D, wl, g, MED) / 2
    ax.plot(wl, d / mono.max(), color=col, lw=2, label=f"dimer, gap {g:.1f} nm")
ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("extinction per particle (norm.)")
ax.set_title("Interparticle gap is the dominant lever on the red tail")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig4_gap_sensitivity.png")); plt.close(fig)

print("wrote figures to outputs/")
print(f"  monomer-only red-tail index  = {red_tail_index(wl, mono_only):.2f} nm")
print(f"  speciation-mix red-tail index = {red_tail_index(wl, mixture):.2f} nm  "
      f"(x{red_tail_index(wl, mixture)/red_tail_index(wl, mono_only):.2f})")
print(f"  isosbestic ~ {iso:.0f} nm")
