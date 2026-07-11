"""Precomputed exact-optics (T-matrix) basis, with interpolation.

The exact T-matrix backend is accurate but too slow to call inside a fit loop.
This module precomputes per-cluster extinction spectra on a grid of primary
diameters x interparticle gaps (built once with the mstm-env venv), saves them
to .npz, and provides an interpolator with the SAME signature as
`clusters.species_spectrum` — so the exact optics can be used inside fits from
the *system* interpreter (no treams needed at fit time).

Build (in the venv):
    mstm-env/bin/python scripts/build_tmatrix_basis.py
Use (system python):
    from aunp_speciation.basis_cache import load_cache
    cache = load_cache("outputs/tmatrix_basis.npz")
    basis = species_basis(wl, D, poly, gap, backend=cache.species_fn)
"""

from __future__ import annotations
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .mie import monomer_cross_sections


def build_grid(diameters_nm, gaps_nm, wavelength_nm,
               species=("dimer", "trimer_linear", "trimer_triangular"),
               n_medium="water", lmax=8):
    """Compute exact per-cluster spectra over (species, D, gap, lambda).

    Returns a dict ready for np.savez. Requires treams (run in mstm-env).
    Monomer is not cached (Mie is exact & instant at fit time).
    """
    from .clusters_tmatrix import species_spectrum_tmatrix, available
    if not available():
        raise RuntimeError("build_grid needs treams — run with mstm-env/bin/python")
    D = np.asarray(diameters_nm, float)
    G = np.asarray(gaps_nm, float)
    wl = np.asarray(wavelength_nm, float)
    data = {}
    for sp in species:
        cube = np.zeros((len(D), len(G), len(wl)))
        for i, d in enumerate(D):
            for j, g in enumerate(G):
                cube[i, j] = species_spectrum_tmatrix(sp, d, wl, g, n_medium, lmax)
        data[sp] = cube
    from .dielectric import current_gold_model
    return dict(diameters=D, gaps=G, wavelength=wl,
                species=np.array(species), medium=np.array(str(n_medium)),
                gold_model=np.array(current_gold_model()),
                lmax=lmax, **{f"cube__{k}": v for k, v in data.items()})


def save_cache(path, grid):
    np.savez_compressed(path, **grid)


class CachedBasis:
    """Interpolated exact-optics species spectra."""

    def __init__(self, npz):
        self.D = npz["diameters"]
        self.G = npz["gaps"]
        self.wl = npz["wavelength"]
        self.medium = str(npz["medium"])
        # dielectric the cluster cubes were built with (older caches: unknown).
        # The fit-time monomer (Mie) uses the CURRENT module default — keep the
        # two consistent or the basis mixes dielectrics.
        self.gold_model = str(npz["gold_model"]) if "gold_model" in npz else "unknown"
        self.species_list = [str(s) for s in npz["species"]]
        self._interp = {}
        for sp in self.species_list:
            cube = npz[f"cube__{sp}"]
            self._interp[sp] = RegularGridInterpolator(
                (self.D, self.G, self.wl), cube,
                bounds_error=False, fill_value=None)  # linear extrapolation

    def species_fn(self, species, diameter_nm, wavelength_nm, gap_nm=1.0,
                   n_medium="water"):
        """Signature-compatible with clusters.species_spectrum."""
        wl = np.atleast_1d(np.asarray(wavelength_nm, float))
        if species == "monomer":
            return monomer_cross_sections(diameter_nm, wl, n_medium)["ext"]
        interp = self._interp[species]
        pts = np.column_stack([np.full(wl.shape, diameter_nm),
                               np.full(wl.shape, gap_nm), wl])
        return interp(pts)


def load_cache(path):
    with np.load(path, allow_pickle=True) as npz:
        return CachedBasis(npz)
