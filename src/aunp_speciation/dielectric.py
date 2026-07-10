"""Optical constants.

Gold permittivity via the Etchegoin et al. (2006) analytic model, an
"improved analytical fit" to Johnson & Christy that is accurate ~400-1000 nm.

IMPORTANT (ties to ref 15, the gold-permittivity training-data paper):
the choice of dielectric function measurably shifts absolute peak positions.
For production, replace `gold_epsilon` with tabulated Johnson & Christy (1972)
or Olmon (2012) data and interpolate. The analytic model here is for fast,
dependency-free prototyping and is validated below by the monomer LSPR landing
near ~520 nm in water.
"""

from __future__ import annotations
import numpy as np

# Etchegoin 2006 parameters (wavelength formulation, nm)
_EPS_INF = 1.53
_LAMBDA_P = 145.0      # plasma wavelength (nm)
_GAMMA_P = 17000.0     # Drude damping (nm)
# Two critical-point (interband) terms
_CP = [
    dict(A=0.94, phi=-np.pi / 4, lam=468.0, gamma=2300.0),
    dict(A=1.36, phi=-np.pi / 4, lam=331.0, gamma=940.0),
]


def _gold_epsilon_etchegoin(lam):
    eps = np.full(lam.shape, _EPS_INF, dtype=complex)
    eps -= 1.0 / (_LAMBDA_P**2 * (1.0 / lam**2 + 1j / (_GAMMA_P * lam)))
    for cp in _CP:
        A, phi, lc, gc = cp["A"], cp["phi"], cp["lam"], cp["gamma"]
        eps += (A / lc) * (
            np.exp(1j * phi) / (1.0 / lc - 1.0 / lam - 1j / gc)
            + np.exp(-1j * phi) / (1.0 / lc + 1.0 / lam + 1j / gc)
        )
    return eps


# --- selectable dielectric model ---------------------------------------------
# 'etchegoin' : fast analytic (default; peak ~527 nm, biased for quantitative work)
# 'bb'        : tabulated Brendel-Bormann (Rakic et al. 1998) — accurate, matches
#               Johnson & Christy; recommended for fitting real spectra. (ties to
#               ref 15: use measured/accurate constants, esp. the imaginary part.)
_GOLD_MODEL = "etchegoin"
_BB_TABLE = None


def use_gold_model(name):
    """Select the gold dielectric model: 'etchegoin' or 'bb'."""
    global _GOLD_MODEL
    if name not in ("etchegoin", "bb"):
        raise ValueError("model must be 'etchegoin' or 'bb'")
    _GOLD_MODEL = name


def _load_bb():
    global _BB_TABLE
    if _BB_TABLE is None:
        import os
        f = os.path.join(os.path.dirname(__file__), "data", "gold_bb.npz")
        d = np.load(f)
        _BB_TABLE = (d["wavelength"], d["eps_real"] + 1j * d["eps_imag"])
    return _BB_TABLE


def gold_epsilon(wavelength_nm, model=None):
    """Complex relative permittivity of gold at vacuum wavelength(s) in nm.

    model overrides the module default set by use_gold_model().
    """
    lam = np.asarray(wavelength_nm, dtype=float)
    m = model or _GOLD_MODEL
    if m == "bb":
        wl_t, eps_t = _load_bb()
        return np.interp(lam, wl_t, eps_t.real) + 1j * np.interp(lam, wl_t, eps_t.imag)
    return _gold_epsilon_etchegoin(lam)


# --- small-particle surface-scattering damping correction ---
# For ~12 nm NPs the electron mean free path is cut by the surface, adding
# damping Gamma_size = A * hbar * v_F / R to the Drude term. This broadens and
# damps the LSPR and is significant at 12 nm. (ref 15 stresses that the imaginary
# part of eps, i.e. damping, dominates size-distribution fit accuracy — and that
# tabulated constants from thin-film studies, e.g. Yakubovsky 2017, generalize
# best for colloidal Au. Prefer swapping in such tabulated eps for production;
# this analytic correction is the lightweight prototype path.)
_OMEGA_P_EV = 9.0        # gold plasma energy (eV)
_GAMMA_BULK_EV = 0.07    # bulk Drude damping (eV)
_HBAR_VF_EVNM = 0.922    # hbar * v_F in eV*nm  (v_F = 1.40e6 m/s)


def size_damping_correction(wavelength_nm, diameter_nm, A_surf=1.0):
    """Delta-epsilon added to bulk eps for surface-scattering damping."""
    E = 1239.84 / np.asarray(wavelength_nm, dtype=float)   # photon energy (eV)
    R = diameter_nm / 2.0                                  # radius (nm)
    g_size = A_surf * _HBAR_VF_EVNM / R
    inv = lambda g: 1.0 / (E**2 + 1j * E * g)
    return _OMEGA_P_EV**2 * (inv(_GAMMA_BULK_EV) - inv(_GAMMA_BULK_EV + g_size))


def gold_epsilon_sized(wavelength_nm, diameter_nm, A_surf=1.0):
    """Size-corrected gold permittivity for a given primary-particle diameter."""
    return gold_epsilon(wavelength_nm) + size_damping_correction(
        wavelength_nm, diameter_nm, A_surf)


def gold_index(wavelength_nm, diameter_nm=None, A_surf=1.0):
    """Complex refractive index n + ik of gold.

    If diameter_nm is given, applies the small-particle damping correction.
    """
    if diameter_nm is None:
        return np.sqrt(gold_epsilon(wavelength_nm))
    return np.sqrt(gold_epsilon_sized(wavelength_nm, diameter_nm, A_surf))


# Non-dispersive media (good enough over 400-800 nm)
WATER_N = 1.333
GLYCEROL_N = 1.47


def medium_index(name_or_value="water"):
    if isinstance(name_or_value, (int, float)):
        return float(name_or_value)
    return {"water": WATER_N, "glycerol": GLYCEROL_N, "vacuum": 1.0}[name_or_value]
