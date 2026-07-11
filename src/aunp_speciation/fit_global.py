"""Layer 3 (global): joint fit of a temperature series of spectra.

A single spectrum under-determines size, polydispersity, and the dimer/trimer
split (see fitting.py). The temperature series breaks these degeneracies because
all spectra share ONE sample: the primary size, polydispersity, and gap are
common to every temperature, the total gold is conserved, and the populations
move with T according to van 't Hoff. The isosbestic point is the visible
consequence of that shared structure.

Model for the spectrum at temperature T_i (same sample, conserved gold):

    E_i(lambda) = A * sum_s [ f_s(T_i) / k_s ] * C_s(lambda; D, poly, gap)
                  + A_sca_i * (lambda / 550)^(-n_sca)

  - C_s : per-cluster extinction basis (Layer 2), shared across T (depends on
          D, polydispersity, gap only).
  - f_s(T_i) : gold fraction in species s at T_i, from the equilibrium model
          with van 't Hoff K(T) (equilibrium.py).
  - k_s : particles per cluster (1,2,3) -> converts gold fraction to number.
  - A : single global amplitude (path length x total gold), profiled out
        linearly across ALL temperatures at once (keeps cross-T amplitude
        information that encodes the speciation shift).
  - A_sca_i >= 0, n_sca: broadband scattering pedestal from the UN-ENUMERATED
        minority population of large aggregates (limitation #8: the real red
        wing is a near-flat pedestal no monomer/dimer/trimer basis produces).
        One amplitude per temperature (it grows with T), ONE shared exponent
        n_sca fitted across the series, bounded [0, 6]. n_sca is the
        discriminator: ~4 => Rayleigh-like scattering from large aggregates;
        ~0 => flat instrumental offset. All A_sca_i and A are profiled jointly
        by NNLS (they are linear); only n_sca adds a nonlinear parameter.
        DESIGN NOTE for >=60 nm work: the species basis C_s (Mie/T-matrix
        extinction) already carries each enumerated species' own scattering —
        at 60 nm even monomers scatter 15-30% of extinction. A_sca must ONLY
        absorb what is not in the basis (the large-aggregate tail); keep the
        basis exact so the pedestal term does not steal the species'
        scattering, and expect partial degeneracy there.

Fitted parameters: D, polydispersity, dH2, dS2, dH3, dS3, n_sca (gap fixed by
default; enable by widening its bound). C_tot is fixed to 1 (absorbed into A/K).
"""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from scipy.optimize import least_squares, nnls

from .spectra import species_basis
from .equilibrium import association_constants, gold_fractions

_K = {"monomer": 1, "dimer": 2, "trimer_linear": 3, "trimer_triangular": 3}

_SCA_REF_NM = 550.0   # reference wavelength: A_sca is the pedestal height there


def _sca_shape(wl, n_sca):
    """Broadband-scattering pedestal shape (lambda/550)^(-n), unit at 550 nm."""
    return (np.asarray(wl, float) / _SCA_REF_NM) ** (-n_sca)


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
    # scattering pedestal: per-T amplitude at 550 nm, in the same normalized
    # units as `models` (i.e. fraction of the series maximum), + shared exponent
    a_sca: np.ndarray = field(default=None)
    n_sca: float = np.nan


def _profile_linear(shapes, sca_vec, data_stack):
    """Jointly profile the linear amplitudes by NNLS: one global species
    amplitude + one pedestal amplitude per temperature, all constrained >= 0.
    Returns (coeffs, model_stack); coeffs[0] = A, coeffs[1:] = A_sca per T."""
    nT, nwl = shapes.shape
    M = np.zeros((nT * nwl, 1 + nT))
    M[:, 0] = shapes.reshape(-1)
    for i in range(nT):
        M[i * nwl:(i + 1) * nwl, 1 + i] = sca_vec
    c, _ = nnls(M, data_stack)
    return c, M @ c


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
        # Two-sided: dH and dS may take either sign. One-sided (<=0) bounds
        # would forbid endothermic/entropy-driven association (dH>0, dS>0),
        # i.e. thermally-INDUCED aggregation — which real samples show.
        bounds = ([8.0, 0.5, -150.0, -0.6, -150.0, -0.6],
                  [20.0, 12.0, 150.0, 0.6, 150.0, 0.6])
    # n_sca appended as the 7th parameter; 6-long x0/bounds from older callers
    # are extended automatically. Start neutral (2.0) — do NOT seed at the
    # Rayleigh value, n must be earned from the data.
    x0 = list(x0)
    if len(x0) == 6:
        x0 = x0 + [2.0]
    if len(bounds[0]) == 6:
        bounds = (list(bounds[0]) + [0.0], list(bounds[1]) + [6.0])

    def residual(theta):
        shapes, _ = _shapes_for(theta[:6], wlf, species, gap_nm, n_medium,
                                temps, C_tot, backend, n_sizes)
        sca_vec = _sca_shape(wlf, theta[6])
        _, model_stack = _profile_linear(shapes, sca_vec, data_stack)
        return model_stack - data_stack

    res = least_squares(residual, x0=x0, bounds=bounds, method="trf",
                        x_scale=[10, 5, 40, 0.1, 70, 0.2, 2],
                        diff_step=3e-3, max_nfev=max_nfev)
    theta = res.x
    shapes, fracs = _shapes_for(theta[:6], wlf, species, gap_nm, n_medium,
                                temps, C_tot, backend, n_sizes)
    sca_vec = _sca_shape(wlf, theta[6])
    c, model_stack = _profile_linear(shapes, sca_vec, data_stack)
    A, a_sca = float(c[0]), c[1:]
    models = model_stack.reshape(shapes.shape)
    resid = model_stack - data_stack
    m, n = len(data_stack), len(theta)
    s2 = 2.0 * res.cost / max(m - n, 1)
    # pseudo-inverse: identifiable params get real error bars; flat directions
    # (e.g. polydispersity, which UV-Vis barely constrains; or n_sca when
    # A_sca ~ 0) get large ones.
    cov = np.linalg.pinv(res.jac.T @ res.jac) * s2
    sd = np.sqrt(np.abs(np.diag(cov)))
    names = ["diameter", "pct_poly", "dH2", "dS2", "dH3", "dS3", "n_sca"]

    return GlobalFitResult(
        diameter_nm=float(theta[0]), pct_poly=float(theta[1]),
        dH2=float(theta[2]), dS2=float(theta[3]),
        dH3=float(theta[4]), dS3=float(theta[5]),
        amplitude=A * scale, species=tuple(species), temps_K=temps,
        gold_fractions=fracs, models=models, wavelength=wlf,
        residual_rms=float(np.sqrt(np.mean(resid**2))),
        param_sd={names[i]: float(sd[i]) for i in range(n)},
        success=bool(res.success), message=res.message,
        a_sca=a_sca, n_sca=float(theta[6]),
    )
