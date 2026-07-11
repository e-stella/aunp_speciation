"""Generate the four key demonstration figures.

  fig1_species_spectra.png  - per-particle spectra of monomer/dimer/trimers
  fig2_broadening.png       - monomer-only (4% poly) vs realistic speciation mix
  fig3_isosbestic.png       - temperature series -> isosbestic point, computed
                              with the PRODUCTION stack (jc + gamma_S + cached
                              exact T-matrix optics + eps(T)); two panels:
                              speciation-only (stationary crossing, red of the
                              peak) vs +eps(T) (crossing shifts red but stays
                              nearly stationary — the measured 24 nm drift
                              needs evolving cluster geometry, CLAUDE.md #12)
  fig4_gap_sensitivity.png  - dimer red tail vs interparticle gap

Figs 1/2/4 keep the legacy prototype configuration (etchegoin, CDA) — they
illustrate mechanisms, not magnitudes. Fig 3 was upgraded because its output
(the crossing wavelength) is compared against experiment: the legacy CDA/
etchegoin version put the crossing BLUE of the peak, an artifact of CDA's
weak shift + amplitude enhancement (see CLAUDE.md limitation #1).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib.pyplot as plt

from aunp_speciation.spectra import monomer_polydisperse, species_basis, mix
from aunp_speciation.clusters import species_spectrum
from aunp_speciation.equilibrium import (association_constants,
                                         solve_populations, gold_fractions)

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
# PRODUCTION STACK: jc dielectric (12.9 nm peak ~523 nm), gamma_S surface
# damping, cached exact T-matrix optics when available (CDA fallback), and
# endothermic K(T) so heating grows the clusters (the lab-hypothesis
# direction). Panel (a) freezes eps at the 20 C reference — speciation alone
# then gives a TRUE stationary isosbestic, RED of the monomer peak (~532 nm;
# the legacy CDA/etchegoin version put it BLUE of the peak — an artifact).
# Panel (b) adds gold+water eps(T): the crossing shifts red ~7 nm but stays
# NEARLY STATIONARY (drift <1 nm). Measured finding, not a bug: a fixed-
# geometry equilibrium + eps(T) cannot reproduce the measured 544->568 nm
# drift — the drift requires evolving cluster geometry (CLAUDE.md #12).
from aunp_speciation import dielectric
from aunp_speciation.basis_cache import load_cache

C_tot = 1.0
D3, POLY3, GAP3 = 12.9, 7.0, 1.0       # the C500 sample (TEM 12.9 nm +-7%)
wl3 = np.arange(420.0, 801.0, 1.0)     # stay inside the cache grid (420-800)
temps_C = np.array([15, 25, 35, 45, 55, 65, 75], float)
cmap = plt.cm.plasma(np.linspace(0.1, 0.85, len(temps_C)))
# endothermic demo thermodynamics: mostly monomer at 15 C, ~40% clustered gold
# at 75 C (demonstration values, not fitted constants)
THERMO = dict(dH2=30.0, dS2=0.081, dH3=50.0, dS3=0.13)
SPECIES3 = ("monomer", "dimer", "trimer_linear")

_prev_gold = dielectric.current_gold_model()
dielectric.use_gold_model("jc")
backend3, size_corr3, optics_lab = "cda", True, "CDA fallback (no cache)"
try:
    _cache = load_cache(os.path.join(OUT, "tmatrix_basis.npz"))
    if _cache.gold_model == "jc" and _cache.temps is not None:
        backend3, size_corr3 = _cache.species_fn, _cache.size_correction
        optics_lab = "exact T-matrix (cached)"
except FileNotFoundError:
    pass


def mix_per_gold(b, T_K):
    """Population-weighted extinction per unit total gold (gold conserved)."""
    K2, K3 = association_constants(T_K, **THERMO)
    c1, c2, c3 = solve_populations(C_tot, K2, K3)
    return (c1 * b["monomer"] + c2 * b["dimer"]
            + c3 * b["trimer_linear"]) / C_tot


def rising_crossings(wl_, curves_):
    """Bleach->growth zero of E(T_i)-E(T_0) for each T_i>T_0, 500-700 nm."""
    ref, out = curves_[0], []
    m = (wl_ >= 500) & (wl_ <= 700)
    w = wl_[m]
    for E in curves_[1:]:
        d = (E - ref)[m]
        idx = np.where((d[:-1] < 0) & (d[1:] >= 0))[0]
        if len(idx) == 0:
            out.append(np.nan)
            continue
        i = idx[0]
        out.append(w[i] - d[i] * (w[i + 1] - w[i]) / (d[i + 1] - d[i]))
    return np.array(out)


# panel (a): eps frozen at the 20 C reference — populations are the ONLY
# T-dependence, so the crossing must be stationary
basis_fixed = species_basis(wl3, D3, POLY3, GAP3, MED, species=SPECIES3,
                            backend=backend3, n_sizes=7, temperature_C=None,
                            size_correction=size_corr3)
curves_a = np.array([mix_per_gold(basis_fixed, T + 273.15) for T in temps_C])
# panel (b): full production physics — gold eps(T) + water n(T) per T
curves_b = np.array([
    mix_per_gold(species_basis(wl3, D3, POLY3, GAP3, MED, species=SPECIES3,
                               backend=backend3, n_sizes=7, temperature_C=T,
                               size_correction=size_corr3), T + 273.15)
    for T in temps_C])

iso_a = rising_crossings(wl3, curves_a)
iso_b = rising_crossings(wl3, curves_b)
peak3 = wl3[np.argmax(curves_a[0])]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
for ax, curves, iso, title in [
        (axes[0], curves_a, iso_a,
         "(a) speciation only (ε at 20 °C): stationary isosbestic"),
        (axes[1], curves_b, iso_b,
         "(b) + gold & water ε(T): crossing shifts red, drift <1 nm")]:
    for E, T_C, col in zip(curves, temps_C, cmap):
        ax.plot(wl3, E, color=col, lw=1.6, label=f"{T_C:.0f} °C")
    ax.axvline(peak3, color="0.6", ls=":", lw=1)
    ax.set_xlabel("wavelength (nm)")
    ax.set_title(title, fontsize=10)
    # the crossing is invisible at full scale — zoom on it
    axins = ax.inset_axes([0.52, 0.45, 0.45, 0.5])
    for E, col in zip(curves, cmap):
        axins.plot(wl3, E, color=col, lw=1.2)
    axins.scatter(iso, np.interp(iso, wl3, curves[0]), s=22, marker="x",
                  color="k", zorder=5)
    band = (wl3 >= 518) & (wl3 <= 556)
    axins.set_xlim(518, 556)
    axins.set_ylim(curves[:, band].min() * 0.99, curves[:, band].max() * 1.01)
    if np.all(np.isnan(iso)):
        # CDA fallback: the under-coupled dimer is brighter per gold across
        # the whole window, so there is no bleach->growth crossing (#1)
        axins.set_title("no crossing 500–700 nm (CDA artifact, #1)",
                        fontsize=8)
    else:
        axins.set_title(
            f"crossing {np.nanmin(iso):.1f}–{np.nanmax(iso):.1f} nm",
            fontsize=8)
    axins.tick_params(labelsize=7)
    ax.indicate_inset_zoom(axins, edgecolor="0.6")
axes[0].annotate(f"peak {peak3:.0f} nm", xy=(peak3, curves_a[0].max()),
                 xytext=(428, curves_a[0].max() * 0.99), fontsize=9,
                 color="0.4")
axes[1].text(0.02, 0.03,
             "measured C500 crossing drifts 544→568 nm —\n"
             "fixed-geometry equilibrium + ε(T) can't reproduce that:\n"
             "the drift needs evolving cluster geometry (#12)",
             transform=axes[1].transAxes, fontsize=8, color="0.25")
axes[0].set_ylabel("extinction per total gold (a.u.)")
axes[0].legend(frameon=False, ncol=2, fontsize=8, loc="lower left",
               title="temperature")
fig.suptitle(f"Endothermic monomer⇌dimer⇌trimer, D={D3} nm, "
             f"{optics_lab}, jc + γ_S", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig3_isosbestic.png")); plt.close(fig)
dielectric.use_gold_model(_prev_gold)   # figs 1/2/4 keep the legacy config

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
print(f"  fig3 [{optics_lab}]: monomer peak {peak3:.1f} nm")
if np.all(np.isnan(iso_a)):
    print("    speciation-only: NO crossing in 500-700 nm — CDA artifact "
          "(#1); build outputs/tmatrix_basis.npz for the exact result")
else:
    print(f"    speciation-only isosbestic ~ {np.nanmean(iso_a):.1f} nm "
          f"(spread {np.nanstd(iso_a):.2f} nm) — should sit RED of the peak")
print(f"    with eps(T): crossing drifts {iso_b[0]:.1f} -> {iso_b[-1]:.1f} nm")
for T_C in (temps_C[0], temps_C[-1]):
    K2, K3 = association_constants(T_C + 273.15, **THERMO)
    f = gold_fractions(C_tot, K2, K3)
    print(f"    gold fractions @{T_C:.0f} C: mono {f['monomer']:.2f} / "
          f"dimer {f['dimer']:.2f} / trimer {f['trimer']:.2f}")
