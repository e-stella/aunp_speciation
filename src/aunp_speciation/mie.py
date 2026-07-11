"""Mie theory for a single sphere.

Computes Mie coefficients a_n, b_n (Bohren & Huffman convention) using
spherical Bessel functions, valid for complex refractive index. Provides:
  - monomer extinction/scattering/absorption cross sections
  - the dipolar polarizability alpha (from a_1) used by the coupled-dipole
    cluster model.

Everything is done in the surrounding medium: size parameter and wavenumber
use the medium refractive index; the particle enters as the *relative* index
m = n_particle / n_medium.
"""

from __future__ import annotations
import numpy as np
from scipy.special import jv, yv

from .dielectric import gold_index, medium_index


def _sph_jn_complex(n, z):
    """Spherical Bessel j_n(z) for complex z (array of orders n)."""
    z = complex(z)
    return np.sqrt(np.pi / (2.0 * z)) * jv(n + 0.5, z)


def _sph_jn_yn_real(n, x):
    """Spherical j_n(x), y_n(x) for real x (array of orders n)."""
    x = float(x)
    j = np.sqrt(np.pi / (2.0 * x)) * jv(n + 0.5, x)
    y = np.sqrt(np.pi / (2.0 * x)) * yv(n + 0.5, x)
    return j, y


def mie_ab(m, x, nmax=None):
    """Mie coefficients a_n, b_n for relative index m, size parameter x.

    Returns arrays a[0..nmax-1], b[...] corresponding to orders n=1..nmax.
    """
    if nmax is None:
        nmax = int(np.ceil(x + 4.0 * x ** (1.0 / 3.0) + 2.0)) + 2
    n = np.arange(1, nmax + 1)

    mx = m * x
    # psi_n(y) = y j_n(y); xi_n(x) = x h_n^(1)(x) = x (j_n + i y_n)
    jn_mx = _sph_jn_complex(n, mx)
    jn_mx_m1 = _sph_jn_complex(n - 1, mx)
    jx, yx = _sph_jn_yn_real(n, x)
    jx_m1, yx_m1 = _sph_jn_yn_real(n - 1, x)

    psi_mx = mx * jn_mx
    psi_x = x * jx
    hx = jx + 1j * yx
    hx_m1 = jx_m1 + 1j * yx_m1
    xi_x = x * hx

    # derivatives: d/dz[z f_n(z)] = z f_{n-1}(z) - n f_n(z)
    psi_mx_d = mx * jn_mx_m1 - n * jn_mx
    psi_x_d = x * jx_m1 - n * jx
    xi_x_d = x * hx_m1 - n * hx

    a = (m * psi_mx * psi_x_d - psi_x * psi_mx_d) / (
        m * psi_mx * xi_x_d - xi_x * psi_mx_d
    )
    b = (psi_mx * psi_x_d - m * psi_x * psi_mx_d) / (
        psi_mx * xi_x_d - m * xi_x * psi_mx_d
    )
    return a, b


def monomer_cross_sections(diameter_nm, wavelength_nm, n_medium="water",
                           size_correction=False, A_surf=1.0,
                           temperature_C=None):
    """Extinction, scattering, absorption cross sections (nm^2) of one sphere.

    Vectorized over wavelength. Returns dict of arrays. If size_correction,
    applies the small-particle surface-scattering damping (see dielectric.py).
    temperature_C (default None = 20 C reference) applies BOTH the gold bulk
    thermal damping (limitation #11) and the water n(T) (limitation #13).
    """
    nm = medium_index(n_medium, temperature_C).real   # propagation: real part
    lam = np.atleast_1d(np.asarray(wavelength_nm, dtype=float))
    a_rad = diameter_nm / 2.0

    Cext = np.zeros_like(lam)
    Csca = np.zeros_like(lam)
    d_corr = diameter_nm if size_correction else None
    m_arr = gold_index(lam, d_corr, A_surf, temperature_C=temperature_C) / nm
    for i, l0 in enumerate(lam):
        k = 2.0 * np.pi * nm / l0          # wavenumber in medium (1/nm)
        x = k * a_rad
        a, b = mie_ab(m_arr[i], x)
        n = np.arange(1, len(a) + 1)
        pref = 2.0 * np.pi / k**2
        Cext[i] = pref * np.sum((2 * n + 1) * np.real(a + b))
        Csca[i] = pref * np.sum((2 * n + 1) * (np.abs(a) ** 2 + np.abs(b) ** 2))
    return {"ext": Cext, "sca": Csca, "abs": Cext - Csca, "wavelength": lam}


def dipole_polarizability(diameter_nm, wavelength_nm, n_medium="water",
                          temperature_C=None, size_correction=False, A_surf=1.0):
    """Dynamic dipolar polarizability alpha(lambda) from Mie a_1.

    alpha = (3i / (2 k^3)) a_1, with k the medium wavenumber. This includes
    finite-size (dynamic depolarization + radiative) corrections and is
    consistent with the Mie dipole term. Units: nm^3. Vectorized over lambda.
    """
    nm = medium_index(n_medium, temperature_C).real
    lam = np.atleast_1d(np.asarray(wavelength_nm, dtype=float))
    a_rad = diameter_nm / 2.0
    alpha = np.zeros(lam.shape, dtype=complex)
    d_corr = diameter_nm if size_correction else None
    m_arr = gold_index(lam, d_corr, A_surf, temperature_C=temperature_C) / nm
    for i, l0 in enumerate(lam):
        k = 2.0 * np.pi * nm / l0
        x = k * a_rad
        a1 = mie_ab(m_arr[i], x, nmax=4)[0][0]
        alpha[i] = (3j / (2.0 * k**3)) * a1
    return alpha
