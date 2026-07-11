"""Layer 3: inverse fit of an experimental extinction spectrum.

Given a measured extinction spectrum, recover:
  - nominal (mean) primary-particle diameter
  - % polydispersity of the primary particle
  - number fractions of monomer / dimer / trimer (-> gold fractions)
with uncertainties.

Strategy: separable / variable-projection least squares.
  * Nonlinear parameters theta = (mean_diameter, pct_polydispersity, n_sca) set
    the shape of each species' basis spectrum and of the broadband-scattering
    pedestal (lambda/550)^(-n_sca) (see fit_global.py for the pedestal's
    rationale and the >=60 nm degeneracy caution).
  * For each theta the species mixing weights AND the pedestal amplitude are
    LINEAR, so we profile them out with non-negative least squares (NNLS) at
    every evaluation.
  * An outer `least_squares` optimizes theta on the profiled residual.

Covariances: theta from the outer Jacobian; weights from a linear model at the
optimum; gold fractions by first-order propagation from the weights.

The overall amplitude is arbitrary (extinction in a.u.); fractions are
scale-invariant. The absolute-scale weight is absorbed into the linear weights.
"""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from scipy.optimize import least_squares, nnls

from .spectra import species_basis
from .fit_global import _sca_shape

# particles per cluster, for converting number fractions -> gold fractions
_K = {"monomer": 1, "dimer": 2, "trimer_linear": 3, "trimer_triangular": 3}


@dataclass
class FitResult:
    diameter_nm: float
    diameter_sd: float
    pct_poly: float
    pct_poly_sd: float
    species: tuple
    number_fractions: dict          # clusters, sums to 1
    gold_fractions: dict            # gold atoms, sums to 1
    gold_fractions_sd: dict
    aggregated_gold: float          # robust: 1 - monomer gold fraction
    aggregated_gold_sd: float
    identifiability: dict           # per-param "ok"/"weak" verdicts
    weights: np.ndarray
    model: np.ndarray               # best-fit spectrum on the fit grid
    wavelength: np.ndarray
    residual_rms: float
    success: bool
    message: str = ""
    # broadband-scattering pedestal A_sca*(lambda/550)^(-n): amplitude at
    # 550 nm (normalized units) + exponent. See fit_global.py's model docstring
    # (incl. the >=60 nm degeneracy caution). A single spectrum constrains
    # n_sca weakly — trust the temperature-series fit for it.
    a_sca: float = 0.0
    n_sca: float = np.nan


def _build_basis_matrix(wl, D, p, gap, med, species, backend, n_sizes):
    """(n_wl, n_species) matrix of per-cluster species spectra."""
    basis = species_basis(wl, D, p, gap, med, species=species,
                          backend=backend, n_sizes=n_sizes)
    return np.column_stack([basis[s] for s in species])


def fit_spectrum(
    wavelength_nm,
    ext,
    species=("monomer", "dimer", "trimer_linear"),
    gap_nm=1.0,
    n_medium="water",
    d_bounds=(8.0, 20.0),
    p_bounds=(0.5, 12.0),
    d0=11.0,
    p0=6.0,
    fit_stride=2,
    backend="cda",
    n_sizes=7,
):
    """Fit an extinction spectrum. Returns a FitResult.

    fit_stride subsamples the wavelength grid for speed during optimization
    (the CDA cluster solve is the cost); the returned model is on that grid.
    backend: optics engine, same options as species_basis ('cda', 'tmatrix',
    or a cached interpolator's species_fn) — MUST match whatever produced the
    data being fitted, else the mismatch is absorbed by D and the scattering
    pedestal (demonstrated: exact-optics data + CDA basis rails D and fakes
    a_sca even on clean synthetics).
    """
    wl = np.asarray(wavelength_nm, dtype=float)
    y = np.asarray(ext, dtype=float)
    wlf = wl[::fit_stride]
    yf = y[::fit_stride]
    yf = yf / yf.max()  # scale-free target

    def weights_for(theta):
        D, p, n_sca = theta
        B = _build_basis_matrix(wlf, D, p, gap_nm, n_medium, species,
                                backend, n_sizes)
        # extra non-negative column: the scattering pedestal (amplitude linear,
        # profiled with the species weights; only the exponent is nonlinear)
        B = np.column_stack([B, _sca_shape(wlf, n_sca)])
        w, _ = nnls(B, yf)
        return B, w

    def residual(theta):
        B, w = weights_for(theta)
        return B @ w - yf

    res = least_squares(
        residual,
        x0=[d0, p0, 2.0],
        bounds=([d_bounds[0], p_bounds[0], 0.0],
                [d_bounds[1], p_bounds[1], 6.0]),
        method="trf",
        x_scale=[10.0, 5.0, 2.0],
        diff_step=1e-3,
    )
    D_hat, p_hat, nsca_hat = res.x
    B, w = weights_for(res.x)
    model = B @ w
    a_sca = float(w[-1])
    B, w = B[:, :-1], w[:-1]   # species-only blocks for the fraction math
    m, n = len(yf), 3
    dof = max(m - n, 1)
    s2 = 2.0 * res.cost / dof
    # theta covariance from Gauss-Newton approx (pinv: n_sca is a flat
    # direction whenever the profiled pedestal amplitude is ~0)
    cov_theta = np.linalg.pinv(res.jac.T @ res.jac) * s2
    D_sd, p_sd, _ = np.sqrt(np.abs(np.diag(cov_theta)))

    # weight covariance from linear model at theta*
    resid = model - yf
    sw2 = np.sum(resid**2) / max(m - len(species), 1)
    try:
        cov_w = sw2 * np.linalg.inv(B.T @ B)
    except np.linalg.LinAlgError:
        cov_w = np.full((len(species), len(species)), np.nan)

    # number fractions
    wsum = w.sum() if w.sum() > 0 else 1.0
    num_frac = {s: w[i] / wsum for i, s in enumerate(species)}
    # gold fractions g_s = k_s w_s / sum(k w)
    k = np.array([_K[s] for s in species], dtype=float)
    kw = k * w
    G = kw.sum() if kw.sum() > 0 else 1.0
    gold = kw / G
    # propagate: d g_i / d w_j = (k_i delta_ij G - k_i w_i k_j) / G^2
    Jg = np.zeros((len(species), len(species)))
    for i in range(len(species)):
        for j in range(len(species)):
            Jg[i, j] = (k[i] * (i == j) * G - kw[i] * k[j]) / G**2
    cov_g = Jg @ cov_w @ Jg.T
    gold_sd = np.sqrt(np.abs(np.diag(cov_g)))

    # robust aggregated (non-monomer) gold fraction + propagated sd
    agg = np.array([0.0 if s == "monomer" else 1.0 for s in species])
    aggregated = float(agg @ gold)
    aggregated_sd = float(np.sqrt(np.abs(agg @ cov_g @ agg)))

    # identifiability verdicts (relative uncertainty thresholds)
    def verdict(val, sd, rel=0.5):
        if not np.isfinite(sd):
            return "weak"
        return "ok" if sd <= rel * max(abs(val), 1e-9) else "weak"

    ident = {
        "diameter": verdict(D_hat, D_sd, 0.25),
        "pct_poly": verdict(p_hat, p_sd, 0.5),
        "aggregated_gold": verdict(aggregated, aggregated_sd, 0.5),
        "dimer_vs_trimer_split": "weak",  # collinear red tails; needs T-series
    }

    return FitResult(
        diameter_nm=float(D_hat),
        diameter_sd=float(D_sd),
        pct_poly=float(p_hat),
        pct_poly_sd=float(p_sd),
        species=tuple(species),
        number_fractions={s: float(num_frac[s]) for s in species},
        gold_fractions={s: float(gold[i]) for i, s in enumerate(species)},
        gold_fractions_sd={s: float(gold_sd[i]) for i, s in enumerate(species)},
        aggregated_gold=aggregated,
        aggregated_gold_sd=aggregated_sd,
        identifiability=ident,
        weights=w,
        model=model,
        wavelength=wlf,
        residual_rms=float(np.sqrt(np.mean(resid**2))),
        success=bool(res.success),
        message=res.message,
        a_sca=a_sca,
        n_sca=float(nsca_hat),
    )
