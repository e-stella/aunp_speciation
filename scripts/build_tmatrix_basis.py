"""Precompute the exact T-matrix species basis on a (diameter, gap) grid.

RUN WITH THE VENV:   mstm-env/bin/python scripts/build_tmatrix_basis.py
Writes outputs/tmatrix_basis.npz, loadable from the system interpreter for fast,
exact-optics fits (see basis_cache.py).

Keep the grid modest; each grid point is a treams cluster solve. Widen once you
have time/compute.
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from aunp_speciation import dielectric
from aunp_speciation.basis_cache import build_grid, save_cache
from aunp_speciation.clusters_tmatrix import available

assert available(), "treams not available — use mstm-env/bin/python"
OUT = os.path.join(os.path.dirname(__file__), "..", "outputs", "tmatrix_basis.npz")

# Measured J&C constants — same model fit_real.py uses (cache records this;
# fit_real refuses a mismatched cache).
dielectric.use_gold_model("jc")

# wl MUST cover the full fit window (420-800 nm): the cache interpolator
# linearly EXTRAPOLATES outside its grid (fill_value=None), which silently
# corrupted fits when the old 470-675 nm grid met 420-800 nm data.
diameters = np.array([11., 12., 13., 14., 15.])   # nm (TEM 12.9 +- 7%, margin)
gaps = np.array([0.5, 1.0, 2.0, 3.0])             # nm surface-surface
wl = np.arange(420, 801, 10.0)                    # nm
species = ("dimer", "trimer_linear")
LMAX = 6   # cache-build accuracy/speed tradeoff (bump to 8 for production)

t0 = time.time()
print(f"building {len(diameters)}x{len(gaps)}x{len(wl)} grid for {species} "
      f"@ lmax={LMAX}, gold_model={dielectric.current_gold_model()} ...")
grid = build_grid(diameters, gaps, wl, species=species, lmax=LMAX)
save_cache(OUT, grid)
print(f"wrote {OUT}  ({time.time()-t0:.0f} s)")
