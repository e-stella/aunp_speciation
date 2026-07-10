"""Quick physics sanity checks before plotting."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from aunp_speciation.mie import mie_ab, monomer_cross_sections, dipole_polarizability
from aunp_speciation.clusters import species_spectrum

# 1) Rayleigh limit check: a_1 ~ -i (2/3) x^3 (m^2-1)/(m^2+2) for small x
m = 1.5 + 0.0j
x = 0.01
a1 = mie_ab(m, x, nmax=4)[0][0]
a1_ray = -1j * (2.0 / 3.0) * x**3 * (m**2 - 1) / (m**2 + 2)
print(f"[1] Rayleigh a1 check: mie={a1:.3e}  analytic={a1_ray:.3e}  "
      f"ratio={abs(a1/a1_ray):.4f}")

# 2) Monomer LSPR peak location for 12 nm Au in water (expect ~515-525 nm)
wl = np.arange(400, 800, 1.0)
ext = monomer_cross_sections(12.0, wl, "water")["ext"]
peak = wl[np.argmax(ext)]
print(f"[2] 12 nm Au monomer LSPR peak in water: {peak:.0f} nm  (expect ~515-525)")

# 3) CDA monomer must equal Mie monomer (consistency)
cda_mono = species_spectrum("monomer", 12.0, wl, n_medium="water")
rel = np.max(np.abs(cda_mono - ext) / np.max(ext))
print(f"[3] CDA-monomer vs Mie-monomer max rel diff: {rel:.2e}  (expect ~0)")

# 4) Dimer should red-shift / develop red-side intensity vs monomer
dimer = species_spectrum("dimer", 12.0, wl, gap_nm=1.0, n_medium="water")
peak_d = wl[np.argmax(dimer)]
# red-side intensity ratio (600-700 nm integral) per particle
band = (wl >= 600) & (wl <= 700)
red_mono = ext[band].sum()
red_dim = (dimer[band] / 2).sum()  # per-particle
print(f"[4] dimer peak {peak_d:.0f} nm; red-band(600-700) intensity per particle: "
      f"monomer={red_mono:.3g}  dimer={red_dim:.3g}  "
      f"enhancement x{red_dim/max(red_mono,1e-30):.1f}")
print("OK")
