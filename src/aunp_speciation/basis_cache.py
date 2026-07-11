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
               n_medium="water", lmax=8, temperatures_C=None,
               size_correction=False, A_surf=1.0):
    """Compute exact per-cluster spectra over (species[, T], D, gap, lambda).

    Returns a dict ready for np.savez. Requires treams (run in mstm-env).
    Monomer is not cached (Mie is exact & instant at fit time).

    temperatures_C: iterable of temperatures for the eps(T) axis (limitation
    #11) — cubes become 4-D (nT, nD, nG, nwl) and are interpolated in T at
    fit time. None = legacy single-reference-T 3-D cube.
    size_correction/A_surf: bake the (T-independent) surface damping gamma_S
    into each sphere (limitation #2). Recorded; the fit-time monomer uses the
    same flags.
    """
    from .clusters_tmatrix import species_spectrum_tmatrix, available
    if not available():
        raise RuntimeError("build_grid needs treams — run with mstm-env/bin/python")
    D = np.asarray(diameters_nm, float)
    G = np.asarray(gaps_nm, float)
    wl = np.asarray(wavelength_nm, float)
    temps = None if temperatures_C is None else np.asarray(temperatures_C, float)
    data = {}
    for sp in species:
        if temps is None:
            cube = np.zeros((len(D), len(G), len(wl)))
            for i, d in enumerate(D):
                for j, g in enumerate(G):
                    cube[i, j] = species_spectrum_tmatrix(
                        sp, d, wl, g, n_medium, lmax,
                        size_correction=size_correction, A_surf=A_surf)
        else:
            cube = np.zeros((len(temps), len(D), len(G), len(wl)))
            for t, Tc in enumerate(temps):
                for i, d in enumerate(D):
                    for j, g in enumerate(G):
                        cube[t, i, j] = species_spectrum_tmatrix(
                            sp, d, wl, g, n_medium, lmax, temperature_C=Tc,
                            size_correction=size_correction, A_surf=A_surf)
        data[sp] = cube
    from .dielectric import current_gold_model, EPS_T_MODEL
    out = dict(diameters=D, gaps=G, wavelength=wl,
               species=np.array(species), medium=np.array(str(n_medium)),
               gold_model=np.array(current_gold_model()),
               size_correction=np.array(bool(size_correction)),
               A_surf=np.array(float(A_surf)),
               lmax=lmax, **{f"cube__{k}": v for k, v in data.items()})
    if temps is not None:
        out["temperatures_C"] = temps
        out["eps_t_model"] = np.array(EPS_T_MODEL)
    return out


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
        # eps(T) axis (limitation #11): None for legacy fixed-eps caches
        self.temps = npz["temperatures_C"] if "temperatures_C" in npz else None
        self.eps_t_model = str(npz["eps_t_model"]) if "eps_t_model" in npz else None
        # surface damping baked into the cubes (limitation #2)
        self.size_correction = bool(npz["size_correction"]) \
            if "size_correction" in npz else False
        self.A_surf = float(npz["A_surf"]) if "A_surf" in npz else 1.0
        self.species_list = [str(s) for s in npz["species"]]
        self._interp = {}
        for sp in self.species_list:
            cube = npz[f"cube__{sp}"]
            axes = (self.D, self.G, self.wl) if self.temps is None else \
                   (self.temps, self.D, self.G, self.wl)
            self._interp[sp] = RegularGridInterpolator(
                axes, cube,
                bounds_error=False, fill_value=None)  # linear extrapolation

    def species_fn(self, species, diameter_nm, wavelength_nm, gap_nm=1.0,
                   n_medium="water", temperature_C=None, size_correction=None,
                   A_surf=None):
        """Signature-compatible with clusters.species_spectrum.

        The cubes' baked-in size_correction/A_surf are authoritative; passing
        a conflicting value raises (rebuild the cache instead). The fit-time
        monomer is computed live with the SAME flags + temperature.
        """
        if size_correction is not None and bool(size_correction) != self.size_correction:
            raise ValueError(
                f"cache was built with size_correction={self.size_correction}; "
                "rebuild it rather than mixing basis physics.")
        if temperature_C is not None and temperature_C != 20.0 and self.temps is None:
            raise ValueError(
                "cache has no temperature axis — rebuild with "
                "build_grid(..., temperatures_C=[...]) for eps(T) fits.")
        wl = np.atleast_1d(np.asarray(wavelength_nm, float))
        if species == "monomer":
            return monomer_cross_sections(
                diameter_nm, wl, n_medium, size_correction=self.size_correction,
                A_surf=self.A_surf, temperature_C=temperature_C)["ext"]
        interp = self._interp[species]
        if self.temps is None:
            pts = np.column_stack([np.full(wl.shape, diameter_nm),
                                   np.full(wl.shape, gap_nm), wl])
        else:
            tc = 20.0 if temperature_C is None else float(temperature_C)
            pts = np.column_stack([np.full(wl.shape, tc),
                                   np.full(wl.shape, diameter_nm),
                                   np.full(wl.shape, gap_nm), wl])
        return interp(pts)


def load_cache(path):
    with np.load(path, allow_pickle=True) as npz:
        return CachedBasis(npz)
