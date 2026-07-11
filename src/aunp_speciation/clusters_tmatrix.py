"""Exact multi-sphere T-matrix backend (MSTM/GMM-equivalent) via `treams`.

This is the quantitative optical engine: unlike the single-dipole coupled-dipole
model in `clusters.py`, it includes all multipole orders, so it captures the
strong near-contact coupling that CDA underestimates.

`treams` needs an older numpy/scipy than the main environment (scipy < 1.13).
Install into a dedicated venv:

    python -m venv mstm-env
    mstm-env/bin/pip install "numpy==1.26.4" "scipy==1.11.4" treams

Then run T-matrix scripts with `mstm-env/bin/python`. The API here mirrors
`clusters.species_spectrum(...)` so it is a drop-in for building the Layer-2
species basis; `spectra.species_basis(..., backend=...)` selects the engine.

Orientation averaging is exact (treams `xs_ext_avg`), not a 3-direction
approximation.
"""

from __future__ import annotations
import numpy as np

from .dielectric import gold_epsilon, medium_index
from .clusters import GEOMETRIES

try:
    import treams
    _HAVE_TREAMS = True
except Exception:  # pragma: no cover
    _HAVE_TREAMS = False


def available():
    return _HAVE_TREAMS


def cluster_cross_section_tmatrix(positions, diameter_nm, wavelength_nm,
                                  n_medium="water", lmax=8, temperature_C=None,
                                  size_correction=False, A_surf=1.0):
    """Exact orientation-averaged extinction (nm^2) per cluster via treams.

    temperature_C: gold eps(T) + water n(T); None = 20 C reference (no-op).
    size_correction: give each sphere the surface-damped (gamma_S) gold eps."""
    if not _HAVE_TREAMS:
        raise RuntimeError(
            "treams not importable in this interpreter; use the mstm-env venv."
        )
    nmed = medium_index(n_medium, temperature_C).real
    lam = np.atleast_1d(np.asarray(wavelength_nm, dtype=float))
    radius = diameter_nm / 2.0
    positions = np.asarray(positions, dtype=float)
    positions = positions - positions.mean(axis=0)  # center for global expansion
    mat_med = treams.Material(nmed**2)

    from .dielectric import gold_epsilon_sized
    Cext = np.zeros(lam.shape)
    for i, l0 in enumerate(lam):
        k0 = 2.0 * np.pi / l0
        if size_correction:
            eps_au = gold_epsilon_sized(l0, diameter_nm, A_surf, temperature_C)
        else:
            eps_au = gold_epsilon(l0, temperature_C=temperature_C)
        mat_au = treams.Material(complex(eps_au))
        sphere = treams.TMatrix.sphere(lmax, k0, radius, [mat_au, mat_med])
        if len(positions) == 1:
            Cext[i] = sphere.xs_ext_avg
        else:
            cl = treams.TMatrix.cluster([sphere] * len(positions), positions)
            cl = cl.interaction.solve()
            # expand local cluster T-matrix to a global (origin-centered) one
            # so the exact orientation average (xs_ext_avg) is defined
            g = cl.expand(treams.SphericalWaveBasis.default(lmax))
            Cext[i] = g.xs_ext_avg
    return Cext


def species_spectrum_tmatrix(species, diameter_nm, wavelength_nm,
                             gap_nm=1.0, n_medium="water", lmax=8,
                             temperature_C=None, size_correction=False,
                             A_surf=1.0):
    """Per-cluster orientation-averaged extinction for a named species (exact)."""
    pos = GEOMETRIES[species](diameter_nm, gap_nm)
    return cluster_cross_section_tmatrix(pos, diameter_nm, wavelength_nm,
                                         n_medium, lmax, temperature_C,
                                         size_correction, A_surf)
