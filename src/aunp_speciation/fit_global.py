"""Layer 3 (global): joint fit of a temperature series of spectra.

A single spectrum under-determines size, polydispersity, and the dimer/trimer
split (see fitting.py). The temperature series breaks these degeneracies because
all spectra share ONE sample: the primary size, polydispersity, and gap are
common to every temperature, the total gold is conserved, and the populations
move with T according to van 't Hoff. The isosbestic point is the visible
consequence of that shared structure.

Model for the spectrum at temperature T_i (same sample, conserved gold):

    E_i(lambda) = A * sum_s  [ f_s(T_i) / k_s ] * C_s(lambda; D, poly, gap)

  - C_s : per-cluster extinction basis (Layer 2), shared across T (depends on
          D, polydispersity, gap only).
  - f_s(T_i) : gold fraction in species s at T_i, from the equilibrium model
          with van 't Hoff K(T) (equilibrium.py).
  - k_s : particles per cluster (1,2,3) -> converts gold fraction to number.
  - A : single global amplitude (path length x total gold), profiled out
        linearly across ALL temperatures at once (keeps cross-T amplitude
        information that encodes the speciation shift).

Fitted parameters: D, polydispersity, dH2, dS2, dH3, dS3  (gap fixed by default;
enable by widening its bound). C_tot is fixed to 1 (absorbed into A / K).
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.optimize import least_squares

from .spectra import species_basis
from .equilibrium import association_constants, gold_fractions

_K = {"monomer": 1, "dimer": 2, "trimer_linear": 3, "trimer_triangular": 3}


@dataclass
class GlobalFitResult:
    diameter_nm: float
    pct_poly: float
    dH2: float
    dS2: float
    dH3: float
    dS3: float
    amplitude: float
    species: tuple
    temps_K: np.ndarray
    gold_fractions: dict            # species -> array over temperatures
    models: np.ndarray              # (nT, nwl) best-fit spectra (fit grid)
    wavelength: np.ndarray
    residual_rms: float
    param_sd: dict
    success: bool
    message: str = ""


def _shapes_for(theta, wlf, species, gap, med, temps, C_tot, backend, n_sizes):
    """Build the (nT, nwl) model-shape matrix (before global amplitude)."""
    D, p, dH2, dS2, dH3, dS3 = theta
    basis = species_basis(wlf, D, p, gap, med, species=species,
                          backend=backend, n_sizes=n_sizes)
    B = np.column_stack([basis[s] for s in species])   # (nwl, nspec)
    k = np.array([_K[s] for s in species], dtype=float)
    shapes = np.zeros((len(temps), len(wlf)))
    fracs = {s: np.zeros(len(temps)) for s in species}
    for i, T in enumerate(temps):
        K2, K3 = association_constants(T, dH2, dS2, dH3, dS3)
        gf = gold_fractions(C_tot, K2, K3)
        f = np.array([gf[_canon(s)] for s in species])
        for j, s in enumerate(species):
            fracs[s][i] = f[j]
        shapes[i] = B @ (f / k)
    return shapes, fracs


def _canon(s):
    # map species name to gold_fractions keys (monomer/dimer/trimer)
    if s.startswith("trimer"):
        return "trimer"
    return s


def fit_temperature_series(
    temps_K,
    wavelength_nm,
    spectra,
    species=("monomer", "dimer", "trimer_linear"),
    gap_nm=1.0,
    n_medium="water",
    C_tot=1.0,
    backend="cda",
    n_sizes=5,
    x0=None,
    bounds=None,
    fit_stride=1,
    max_nfev=80,
):
    """Jointly fit spectra measured at several temperatures. Returns GlobalFitResult.

    spectra: (nT, nwl) array aligned with temps_K and wavelength_nm.
    """
    temps = np.asarray(temps_K, dtype=float)
    wl = np.asarray(wavelength_nm, dtype=float)
    Y = np.asarray(spectra, dtype=float)
    wlf = wl[::fit_stride]
    Yf = Y[:, ::fit_stride]
    scale = Yf.max()
    Yf = Yf / scale
    data_stack = Yf.reshape(-1)

    if x0 is None:
        x0 = [11.0, 5.0, -40.0, -0.10, -70.0, -0.20]
    if bounds is None:
        bounds = ([8.0, 0.5, -120.0, -0.5, -200.0, -0.8],
                  [20.0, 12.0, 0.0, 0.0, 0.0, 0.0])

    def residual(theta):
        shapes, _ = _shapes_for(theta, wlf, species, gap_nm, n_medium,
                                temps, C_tot, backend, n_sizes)
        s = shapes.reshape(-1)
        A = float(np.dot(s, data_stack) / max(np.dot(s, s), 1e-30))
        return A * s - data_stack

    res = least_squares(residual, x0=x0, bounds=bounds, method="trf",
                        x_scale=[10, 5, 40, 0.1, 70, 0.2],
                        diff_step=3e-3, max_nfev=max_nfev)
    theta = res.x
    shapes, fracs = _shapes_for(theta, wlf, species, gap_nm, n_medium,
                                temps, C_tot, backend, n_sizes)
    s = shapes.reshape(-1)
    A = float(np.dot(s, data_stack) / max(np.dot(s, s), 1e-30))
    models = (A * shapes)
    resid = models.reshape(-1) - data_stack
    m, n = len(data_stack), len(theta)
    s2 = 2.0 * res.cost / max(m - n, 1)
    # pseudo-inverse: identifiable params get real error bars; flat directions
    # (e.g. polydispersity, which UV-Vis barely constrains) get large ones.
    cov = np.linalg.pinv(res.jac.T @ res.jac) * s2
    sd = np.sqrt(np.abs(np.diag(cov)))
    names = ["diameter", "pct_poly", "dH2", "dS2", "dH3", "dS3"]

    return GlobalFitResult(
        diameter_nm=float(theta[0]), pct_poly=float(theta[1]),
        dH2=float(theta[2]), dS2=float(theta[3]),
        dH3=float(theta[4]), dS3=float(theta[5]),
        amplitude=A * scale, species=tuple(species), temps_K=temps,
        gold_fractions=fracs, models=models, wavelength=wlf,
        residual_rms=float(np.sqrt(np.mean(resid**2))),
        param_sd={names[i]: float(sd[i]) for i in range(n)},
        success=bool(res.success), message=res.message,
    )
