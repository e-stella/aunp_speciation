"""Assemble ensemble spectra from species + populations, with polydispersity.

Total extinction (per unit volume, arbitrary units):
    sigma(lambda) = sum_species  n_species * <C_ext,species(lambda)>

where n_species is a *number* concentration of that cluster type. Polydispersity
is applied to the primary-particle diameter by integrating the monomer (and,
approximately, the cluster) spectra over a Gaussian size distribution.
"""

from __future__ import annotations
import numpy as np

from .mie import monomer_cross_sections
from .clusters import species_spectrum as _species_spectrum_cda


def _get_species_fn(backend):
    """Return the per-cluster spectrum function for the chosen optics backend.

    backend='cda'      -> fast point-dipole coupled-dipole (clusters.py)
    backend='tmatrix'  -> exact multipole T-matrix (clusters_tmatrix.py, needs
                          the mstm-env venv). Quantitative but slower.
    backend=callable   -> a custom species_fn(species, D, wl, gap, medium), e.g.
                          the interpolator from basis_cache (exact optics, fast).
    """
    if callable(backend):
        return backend
    if backend == "cda":
        return _species_spectrum_cda
    if backend == "tmatrix":
        from .clusters_tmatrix import species_spectrum_tmatrix, available
        if not available():
            raise RuntimeError("tmatrix backend requires treams (mstm-env venv)")
        return species_spectrum_tmatrix
    raise ValueError(f"unknown backend {backend!r}")


def gaussian_sizes(mean_nm, pct_polydispersity, n=15):
    """Return (diameters, weights) sampling a Gaussian size distribution."""
    sigma = mean_nm * pct_polydispersity / 100.0
    if sigma <= 0:
        return np.array([mean_nm]), np.array([1.0])
    ds = np.linspace(mean_nm - 3 * sigma, mean_nm + 3 * sigma, n)
    w = np.exp(-0.5 * ((ds - mean_nm) / sigma) ** 2)
    w /= w.sum()
    return ds, w


def monomer_polydisperse(wavelength_nm, mean_nm, pct_poly, n_medium="water"):
    """Monomer extinction averaged over a Gaussian size distribution."""
    ds, w = gaussian_sizes(mean_nm, pct_poly)
    total = np.zeros_like(np.asarray(wavelength_nm, dtype=float))
    for d, wi in zip(ds, w):
        total += wi * monomer_cross_sections(d, wavelength_nm, n_medium)["ext"]
    return total


def species_basis(wavelength_nm, mean_nm, pct_poly, gap_nm=1.0, n_medium="water",
                  species=("monomer", "dimer", "trimer_linear", "trimer_triangular"),
                  backend="cda", n_sizes=7):
    """Per-cluster extinction spectra for each species (polydispersity applied).

    backend selects the optics engine ('cda' fast, 'tmatrix' exact/slower).
    Polydispersity is applied exactly to the monomer and approximately to the
    clusters by averaging over the size distribution.
    """
    species_fn = _get_species_fn(backend)
    basis = {}
    for sp in species:
        if sp == "monomer":
            basis[sp] = monomer_polydisperse(wavelength_nm, mean_nm, pct_poly, n_medium)
        else:
            ds, w = gaussian_sizes(mean_nm, pct_poly, n=n_sizes)
            acc = np.zeros_like(np.asarray(wavelength_nm, dtype=float))
            for d, wi in zip(ds, w):
                acc += wi * species_fn(sp, d, wavelength_nm, gap_nm, n_medium)
            basis[sp] = acc
    return basis


def mix(basis, weights):
    """Weighted sum of species spectra. weights: dict species->number conc."""
    wl_len = len(next(iter(basis.values())))
    total = np.zeros(wl_len)
    for sp, spec in basis.items():
        total += weights.get(sp, 0.0) * spec
    return total
