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
from aunp_speciation.basis_cache import build_grid, save_cache
from aunp_speciation.clusters_tmatrix import available

assert available(), "treams not available — use mstm-env/bin/python"
OUT = os.path.join(os.path.dirname(__file__), "..", "outputs", "tmatrix_basis.npz")

diameters = np.array([10., 11., 12., 13., 14.])   # nm
gaps = np.array([0.5, 1.0, 2.0])                  # nm surface-surface
wl = np.arange(470, 675, 10.0)                    # nm
species = ("dimer", "trimer_linear")
LMAX = 6   # cache-build accuracy/speed tradeoff (bump to 8 for production)

t0 = time.time()
print(f"building {len(diameters)}x{len(gaps)}x{len(wl)} grid for {species} "
      f"@ lmax={LMAX} ...")
grid = build_grid(diameters, gaps, wl, species=species, lmax=LMAX)
save_cache(OUT, grid)
print(f"wrote {OUT}  ({time.time()-t0:.0f} s)")
