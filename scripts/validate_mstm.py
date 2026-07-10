"""Validate the exact T-matrix engine against Mie (monomer) and CDA (clusters).

Run with the pinned venv:  mstm-env/bin/python scripts/validate_mstm.py
Writes outputs/fig6_mstm_vs_cda.png
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from aunp_speciation.mie import monomer_cross_sections
from aunp_speciation.clusters import species_spectrum as cda_spec
from aunp_speciation.clusters_tmatrix import (
    species_spectrum_tmatrix as mstm_spec, available)

assert available(), "treams not available — use mstm-env/bin/python"
OUT = os.path.join(os.path.dirname(__file__), "..", "outputs")
os.makedirs(OUT, exist_ok=True)

D, GAP, MED = 12.0, 1.0, "water"
wl = np.arange(460, 700, 8.0)
LMAX = 5

# 1) monomer: Mie vs treams  (validation of both)
mie = monomer_cross_sections(D, wl, MED)["ext"]
tm_mono = mstm_spec("monomer", D, wl, GAP, MED)
rel = np.max(np.abs(tm_mono - mie) / mie.max())
print(f"[monomer] Mie vs treams max rel diff: {rel:.2e}  "
      f"(peaks: Mie {wl[mie.argmax()]:.0f} nm, treams {wl[tm_mono.argmax()]:.0f} nm)")

# 2) lmax convergence for a contact dimer
print("[dimer lmax convergence @ gap 1 nm]")
for L in (2, 3, 5):
    d = mstm_spec("dimer", D, np.array([600.0]), GAP, MED, lmax=L)[0]
    print(f"    lmax={L}: Cext(600nm)={d:.1f} nm^2")

# 3) dimer & trimer: CDA vs exact
fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True)
for ax, sp, title in [(axes[0], "dimer", "dimer"),
                      (axes[1], "trimer_linear", "linear trimer")]:
    cda = cda_spec(sp, D, wl, GAP, MED)
    tm = mstm_spec(sp, D, wl, GAP, MED, lmax=LMAX)
    npart = 2 if sp == "dimer" else 3
    ax.plot(wl, mie / mie.max(), "k--", lw=1.3, label="monomer (Mie)")
    ax.plot(wl, cda / npart / mie.max(), color="#0072B2", lw=2,
            label="CDA (1 dipole/sphere)")
    ax.plot(wl, tm / npart / mie.max(), color="#D55E00", lw=2,
            label="exact T-matrix (treams)")
    ax.set_title(f"{title}: CDA under-couples vs exact")
    ax.set_xlabel("wavelength (nm)"); ax.legend(frameon=False, fontsize=9)
    # quantify red-band (600-700) per-particle intensity
    band = (wl >= 600) & (wl <= 700)
    r_cda = np.trapz(cda[band] / npart, wl[band])
    r_tm = np.trapz(tm[band] / npart, wl[band])
    print(f"[{sp}] red-band(600-700) per-particle: CDA={r_cda:.0f}  "
          f"exact={r_tm:.0f}  exact/CDA = x{r_tm/max(r_cda,1e-9):.2f}")
axes[0].set_ylabel("extinction per particle (norm.)")
fig.suptitle(f"{D:.0f} nm Au, gap {GAP:.0f} nm — exact multipole vs point-dipole")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig6_mstm_vs_cda.png"), dpi=130)
print("wrote outputs/fig6_mstm_vs_cda.png")
