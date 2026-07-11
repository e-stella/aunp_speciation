"""Coupled-dipole approximation (CDA) for small clusters of spheres.

At 12 nm each particle is well described by a single dipole (its Mie a_1
polarizability), so a dimer/trimer is a 2-3 dipole coupled linear system.
This is fast, analytic, and accurate in the quasi-static-ish regime relevant
here. MSTM (exact multisphere T-matrix) is the recommended cross-check; see
CLAUDE.md.

Orientation averaging: cluster extent (~24-36 nm) << wavelength, so we neglect
the incident phase variation across the cluster and average the response over
three orthogonal polarizations. Dipole-dipole interaction retains full retarded
Green tensor between sites.
"""

from __future__ import annotations
import numpy as np

from .mie import dipole_polarizability
from .dielectric import medium_index


# ----- geometry generators: center positions (nm), gap = surface-surface -----

def dimer(diameter_nm, gap_nm):
    s = diameter_nm + gap_nm
    return np.array([[0.0, 0.0, 0.0], [s, 0.0, 0.0]])


def trimer_linear(diameter_nm, gap_nm):
    s = diameter_nm + gap_nm
    return np.array([[0.0, 0.0, 0.0], [s, 0.0, 0.0], [2 * s, 0.0, 0.0]])


def trimer_triangular(diameter_nm, gap_nm):
    s = diameter_nm + gap_nm
    return np.array([[0.0, 0.0, 0.0], [s, 0.0, 0.0], [s / 2, s * np.sqrt(3) / 2, 0.0]])


GEOMETRIES = {
    "monomer": lambda D, g: np.array([[0.0, 0.0, 0.0]]),
    "dimer": dimer,
    "trimer_linear": trimer_linear,
    "trimer_triangular": trimer_triangular,
}


def _green_tensor(r_vec, k):
    """Retarded dipole-dipole interaction tensor A_ij (3x3), Draine convention.

    A_ij = exp(ikr)/r [ k^2 (rr - I) + (1 - ikr)/r^2 (3 rr - I) ]
    """
    r = np.linalg.norm(r_vec)
    rhat = r_vec / r
    rr = np.outer(rhat, rhat)
    I = np.eye(3)
    term1 = k**2 * (rr - I)
    term2 = (1.0 - 1j * k * r) / r**2 * (3.0 * rr - I)
    return np.exp(1j * k * r) / r * (term1 + term2)


def cluster_cross_section(positions, diameter_nm, wavelength_nm, n_medium="water",
                          temperature_C=None, size_correction=False, A_surf=1.0):
    """Orientation-averaged extinction cross section (nm^2) per cluster.

    positions: (N,3) sphere centers in nm. Vectorized over wavelength.
    temperature_C: gold eps(T) + water n(T); None = 20 C reference (no-op).
    size_correction: apply the (T-independent) surface-damping gamma_S to each
    sphere's polarizability.
    """
    nm = medium_index(n_medium, temperature_C).real
    lam = np.atleast_1d(np.asarray(wavelength_nm, dtype=float))
    N = len(positions)
    alpha = dipole_polarizability(diameter_nm, lam, n_medium, temperature_C,
                                  size_correction, A_surf)  # (nwl,)

    Cext = np.zeros(lam.shape)
    for iw, l0 in enumerate(lam):
        k = 2.0 * np.pi * nm / l0
        a = alpha[iw]
        # Build 3N x 3N interaction matrix M
        M = np.zeros((3 * N, 3 * N), dtype=complex)
        for i in range(N):
            M[3 * i:3 * i + 3, 3 * i:3 * i + 3] = np.eye(3) / a
            for j in range(N):
                if i == j:
                    continue
                M[3 * i:3 * i + 3, 3 * j:3 * j + 3] = _green_tensor(
                    positions[i] - positions[j], k
                )
        Minv = np.linalg.inv(M)
        # Average over 3 orthogonal incident polarizations (uniform field)
        cext_pol = 0.0
        for d in range(3):
            Einc = np.zeros(3 * N, dtype=complex)
            Einc[d::3] = 1.0
            p = Minv @ Einc
            cext_pol += 4.0 * np.pi * k * np.imag(np.vdot(Einc, p)).real
        Cext[iw] = cext_pol / 3.0
    return Cext


def species_spectrum(species, diameter_nm, wavelength_nm, gap_nm=1.0,
                     n_medium="water", temperature_C=None,
                     size_correction=False, A_surf=1.0):
    """Per-cluster orientation-averaged extinction spectrum for a named species."""
    pos = GEOMETRIES[species](diameter_nm, gap_nm)
    return cluster_cross_section(pos, diameter_nm, wavelength_nm, n_medium,
                                 temperature_C, size_correction, A_surf)
