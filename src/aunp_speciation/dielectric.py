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
# 'bb'        : tabulated Brendel-Bormann (Rakic et al. 1998) model fit. CAUTION:
#               in the 500-540 nm window its eps1 is ~0.4 less negative than
#               measured J&C, so the 13 nm monomer peak lands ~527 nm — same red
#               bias as 'etchegoin' (measured here; do not assume it "matches J&C").
# 'jc'        : tabulated Johnson & Christy (1972) measured n,k — 49 points,
#               187.9-1937 nm, embedded below (public-domain digitization from
#               refractiveindex.info). 13 nm monomer peak ~523 nm in bare water.
# 'yakubovsky25' / 'yakubovsky53' : tabulated thin-film n,k from Yakubovsky,
#               Arsenin, Stebunov, Fedyanin & Volkov, Opt. Express 25, 25574
#               (2017), for 25 nm and 53 nm e-beam-evaporated films (171 points,
#               300-2000 nm; CC0 digitization via refractiveindex.info,
#               main/Au/nk/Yakubovsky-{25nm,53nm}.yml; data file
#               data/gold_yakubovsky.npz). Ref 18 (Klinavicius 2025, JPCC 129,
#               17616) finds these, applied SIZE-MATCHED, are the best gold eps
#               for colloidal AuNP retrieval: the thin film's measured eps2
#               already embeds the surface-scattering/mean-free-path damping
#               that bulk J&C needs a correction formula for. Size rule (make it
#               explicit — use gold_model_for_diameter(), or pick a film):
#               25 nm film for particle RADIUS < 25 nm, 53 nm film otherwise.
_GOLD_MODEL = "etchegoin"
_BB_TABLE = None

# Johnson & Christy, Phys. Rev. B 6, 4370 (1972), gold: (vacuum lambda nm, n, k).
# CC0 digitization via refractiveindex.info (main/Au/nk/Johnson.yml).
_JC_NM_N_K = np.array([
    (187.9, 1.280, 1.188), (191.6, 1.320, 1.203), (195.3, 1.340, 1.226),
    (199.3, 1.330, 1.251), (203.3, 1.330, 1.277), (207.3, 1.300, 1.304),
    (211.9, 1.300, 1.350), (216.4, 1.300, 1.387), (221.4, 1.300, 1.427),
    (226.2, 1.310, 1.460), (231.3, 1.300, 1.497), (237.1, 1.320, 1.536),
    (242.6, 1.320, 1.577), (249.0, 1.330, 1.631), (255.1, 1.330, 1.688),
    (261.6, 1.350, 1.749), (268.9, 1.380, 1.803), (276.1, 1.430, 1.847),
    (284.4, 1.470, 1.869), (292.4, 1.490, 1.878), (300.9, 1.530, 1.889),
    (310.7, 1.530, 1.893), (320.4, 1.540, 1.898), (331.5, 1.480, 1.883),
    (342.5, 1.480, 1.871), (354.2, 1.500, 1.866), (367.9, 1.480, 1.895),
    (381.5, 1.460, 1.933), (397.4, 1.470, 1.952), (413.3, 1.460, 1.958),
    (430.5, 1.450, 1.948), (450.9, 1.380, 1.914), (471.4, 1.310, 1.849),
    (495.9, 1.040, 1.833), (520.9, 0.620, 2.081), (548.6, 0.430, 2.455),
    (582.1, 0.290, 2.863), (616.8, 0.210, 3.272), (659.5, 0.140, 3.697),
    (704.5, 0.130, 4.103), (756.0, 0.140, 4.542), (821.1, 0.160, 5.083),
    (892.0, 0.170, 5.663), (984.0, 0.220, 6.350), (1088.0, 0.270, 7.150),
    (1216.0, 0.350, 8.145), (1393.0, 0.430, 9.519), (1610.0, 0.560, 11.210),
    (1937.0, 0.920, 13.780),
])


_NK_SPLINES = {}


def _nk_spline(name, wl_t, n_t, k_t):
    """Cubic spline of (n, k) vs wavelength, cached by name. Cubic, not linear:
    the J&C grid is ~25 nm apart in the visible and linear interpolation puts
    kinks at the nodes that the LSPR maximum snaps to (measured: peak pinned to
    the 520.9 nm node). Same treatment for every tabulated dataset."""
    if name not in _NK_SPLINES:
        from scipy.interpolate import CubicSpline
        _NK_SPLINES[name] = CubicSpline(wl_t, np.stack([n_t, k_t], axis=1), axis=0)
    spl = _NK_SPLINES[name]
    return lambda lam: spl(lam).T


def _jc_spline():
    wl_t, n_t, k_t = _JC_NM_N_K.T
    return _nk_spline("jc", wl_t, n_t, k_t)


_YK_TABLE = None


def _yk_spline(film):
    """film: '25' or '53' (nm nominal film thickness)."""
    global _YK_TABLE
    if _YK_TABLE is None:
        import os
        f = os.path.join(os.path.dirname(__file__), "data", "gold_yakubovsky.npz")
        _YK_TABLE = np.load(f)
    t = _YK_TABLE
    return _nk_spline(f"yk{film}", t["wavelength"], t[f"n{film}"], t[f"k{film}"])


_MODELS = ("etchegoin", "bb", "jc", "yakubovsky25", "yakubovsky53",
           "reddy_p200")


def use_gold_model(name):
    """Select the gold dielectric model: 'etchegoin', 'bb', 'jc',
    'yakubovsky25', or 'yakubovsky53'.

    'yakubovsky' (no film) is deliberately NOT accepted: the film choice is a
    size-matching decision that must be visible in the calling code — use
    gold_model_for_diameter(D) to apply the ref-18 rule, or name a film.
    """
    global _GOLD_MODEL
    if name == "yakubovsky":
        raise ValueError(
            "'yakubovsky' is size-matched: call "
            "use_gold_model(gold_model_for_diameter(D)) to apply the rule "
            "(25 nm film for radius < 25 nm, else 53 nm film), or select "
            "'yakubovsky25' / 'yakubovsky53' explicitly.")
    if name not in _MODELS:
        raise ValueError(f"model must be one of {_MODELS}")
    _GOLD_MODEL = name


def gold_model_for_diameter(diameter_nm, family="yakubovsky"):
    """Ref-18 size-matching rule (Klinavicius 2025): Yakubovsky 25 nm-film
    eps for particle radius < 25 nm, the 53 nm-film eps otherwise. Returns the
    model NAME so the choice is explicit at the call site."""
    if family != "yakubovsky":
        raise ValueError("size matching is defined for family='yakubovsky'")
    return "yakubovsky25" if diameter_nm / 2.0 < 25.0 else "yakubovsky53"


def current_gold_model():
    """Name of the active gold dielectric model (for provenance/caches)."""
    return _GOLD_MODEL


def _load_bb():
    global _BB_TABLE
    if _BB_TABLE is None:
        import os
        f = os.path.join(os.path.dirname(__file__), "data", "gold_bb.npz")
        d = np.load(f)
        _BB_TABLE = (d["wavelength"], d["eps_real"] + 1j * d["eps_imag"])
    return _BB_TABLE


def gold_epsilon(wavelength_nm, model=None, temperature_C=None):
    """Complex relative permittivity of gold at vacuum wavelength(s) in nm.

    model overrides the module default set by use_gold_model().
    temperature_C (default None = T_REF_C = 20 C, a mathematical no-op) adds
    the bulk thermal Drude-damping correction (limitation #11): every base
    dataset is tabulated/fitted near room T; heating is applied as a delta on
    top. See thermal_damping_correction.
    """
    lam = np.asarray(wavelength_nm, dtype=float)
    m = model or _GOLD_MODEL
    if m == "bb":
        wl_t, eps_t = _load_bb()
        eps = np.interp(lam, wl_t, eps_t.real) + 1j * np.interp(lam, wl_t, eps_t.imag)
    elif m == "jc":
        nk = _jc_spline()(lam)
        eps = (nk[0] + 1j * nk[1])**2
    elif m in ("yakubovsky25", "yakubovsky53"):
        nk = _yk_spline(m[-2:])(lam)
        eps = (nk[0] + 1j * nk[1])**2
    elif m == "reddy_p200":
        # natively T-dependent (full DCP incl. interband) — do NOT add the
        # Drude-only thermal retune on top
        return gold_epsilon_reddy(lam, 23.0 if temperature_C is None
                                  else temperature_C)
    elif m == "yakubovsky":
        use_gold_model(m)  # raises with the size-matching guidance
    else:
        eps = _gold_epsilon_etchegoin(lam)
    if temperature_C is not None and temperature_C != T_REF_C:
        eps = eps + thermal_damping_correction(lam, temperature_C)
    return eps


# --- Drude damping corrections: size (T-independent) + temperature (bulk) ----
# The plasmon relaxation rate decomposes (Chetoui et al. 2026; lit-map §G) as
#     gamma(T) = gamma_bulk(T) + gamma_S,
# where gamma_bulk is the BULK electron-phonon damping (strongly T-dependent)
# and gamma_S = A*hbar*v_F/R is the size/surface (Landau) damping —
# EXPLICITLY temperature-independent. Keep the two separable: do NOT import a
# monolithic "nanoparticle eps(T)" (it would carry someone else's gamma_S and
# environment; gamma_S scales 1/R).
#
# Drude baseline (Olmon et al., PRB 86, 235147 (2012), single-crystal Au;
# AOT 2022 ellipsometry review): hbar*omega_p = 8.45 eV, gamma_bulk(RT) =
# 44 meV. (Previously 9.0 eV / 70 meV generic textbook values — updated;
# affects only the size-corrected path, which is off by default in fits.)
#
# gamma_bulk(T): Holstein electron-phonon form endorsed and used by Reddy,
# Guler, Kildishev, Boltasseva & Shalaev, Opt. Mater. Express 6, 2776 (2016)
# ["Temperature-dependent optical properties of gold thin films" — measured
# Drude+2CP fits 23-500 C; their single-crystal 200 nm film gives
# Gamma_D = 0.0534 eV (23 C) -> 0.0898 eV (500 C)]:
#     gamma_ep(T) = gamma_0 * [ 2/5 + 4 (T/Th)^5 * Int_0^{Th/T} z^4/(e^z-1) dz ],
# Th = Debye temperature of gold = 170 K. Near room T this slope
# (~1.2e-4 eV/K anchored at 44 meV) matches Reddy's measured first-cycle
# single-crystal slope (~1.1e-4 eV/K) to ~12%. We implement the full integral
# (not the high-T limit) and anchor at gamma_bulk(20 C) = 44 meV.
#
# The correction is applied as a RETUNE of the Drude damping on top of any
# base dataset: delta_eps = omega_p^2 * [inv(gamma_ref) - inv(gamma_new)],
# inv(g) = 1/(E^2 + iEg). To leading order in the visible (E >> gamma) this is
# ~ +i * omega_p^2 * (gamma_new - gamma_ref)/E^3 — i.e. the delta depends on
# the CHANGE in damping and only weakly on the assumed reference, so it
# composes cleanly with tabulated bases (jc/yakubovsky) whose internal gamma
# is not exactly 44 meV.
_OMEGA_P_EV = 8.45       # gold Drude plasma energy (eV) — Olmon 2012
_GAMMA_BULK_EV = 0.044   # bulk Drude damping at T_REF (eV) — Olmon 2012
_HBAR_VF_EVNM = 0.922    # hbar * v_F in eV*nm  (v_F = 1.40e6 m/s)
_THETA_D_K = 170.0       # Debye temperature of gold (K) — Reddy 2016
T_REF_C = 20.0           # reference: tabulated datasets are ~room-T
EPS_T_MODEL = "holstein170K-olmon44meV"   # recorded by basis caches


def _holstein_factor(T_K):
    """Holstein e-ph damping shape f(T) = 2/5 + 4 (T/Th)^5 Int_0^{Th/T} z^4/(e^z-1) dz."""
    T_K = float(T_K)
    x_D = _THETA_D_K / T_K
    z = np.linspace(1e-9, x_D, 400)
    trapz = getattr(np, "trapezoid", None) or np.trapz  # numpy<2 compat (mstm-env)
    integral = trapz(z**4 / np.expm1(z), z)
    return 0.4 + 4.0 * (T_K / _THETA_D_K) ** 5 * integral


_F_REF = None  # lazy _holstein_factor(T_REF)

# --- gamma_bulk(T) scaling: a FACTOR-SEVERAL SYSTEMATIC, bracketed -----------
# The three defensible scalings disagree strongly over 15-75 C (numbers are
# dGamma/dT near room T, and the relative rise over 15->75 C w.r.t. 44 meV):
#   'holstein'     : full Holstein integral, Th=170 K -> 1.47e-4 eV/K, +20.1%
#                    (pure electron-phonon theory; NB the high-T limit gives
#                    ~1.2e-4 — an earlier doc claim of "~12% agreement with
#                    Reddy" compared THAT number, not the implemented 1.47e-4,
#                    and was WRONG: Holstein is 36% steeper than Reddy-single)
#   'reddy-single' : Reddy 2016 Table 6, 200 nm single-crystal film, measured
#                    Gamma_D 0.0534 (23 C) -> 0.0725 eV (200 C)
#                    -> 1.08e-4 eV/K, +14.7% over 15->75 C (w.r.t. 44 meV)
#   'reddy-poly'   : Reddy 2016 Table 3, 200 nm POLYcrystalline film,
#                    0.0471 (23 C) -> 0.0489 eV (100 C)
#                    -> 0.23e-4 eV/K, +3.2% over 15->75 C
# Reconciliation: a MEASURED Gamma_D includes T-independent grain-boundary/
# defect scattering that damps the RELATIVE rise; Holstein is the pure e-ph
# term. Citrate-grown AuNPs are polycrystalline, so 'reddy-poly' is arguably
# the closest analogue — but the truth for a colloidal particle is unknown.
# DO NOT leave this as a fixed assumption: bracket fits with all three
# (use_eps_t_scaling) and report how much the conclusions move.
_EPS_T_SCALINGS = {
    "holstein": None,           # full integral (implemented below)
    "reddy-single": 1.08e-4,    # linear dGamma/dT (eV/K), anchored at T_REF
    "reddy-poly": 0.23e-4,
}
_EPS_T_SCALING = "holstein"


def use_eps_t_scaling(name):
    """Select the gamma_bulk(T) scaling: 'holstein' (default), 'reddy-single',
    or 'reddy-poly'. NB: T-matrix caches record the eps(T) model they were
    built with (EPS_T_MODEL); non-holstein scalings will not match a
    holstein-built cache — use the CDA backend for bracket fits, or rebuild."""
    global _EPS_T_SCALING
    if name not in _EPS_T_SCALINGS:
        raise ValueError(f"scaling must be one of {tuple(_EPS_T_SCALINGS)}")
    _EPS_T_SCALING = name


def current_eps_t_scaling():
    return _EPS_T_SCALING


def gamma_bulk_ev(temperature_C=None):
    """Bulk electron-phonon Drude damping (eV) at temperature_C.

    temperature_C=None means the reference (20 C) -> exactly _GAMMA_BULK_EV.
    Scaling selected by use_eps_t_scaling(); see the bracket comment above —
    the magnitude of dGamma/dT is a factor-several systematic between pure
    e-ph theory (Holstein) and Reddy 2016's measured films.
    """
    if temperature_C is None or temperature_C == T_REF_C:
        return _GAMMA_BULK_EV
    slope = _EPS_T_SCALINGS[_EPS_T_SCALING]
    if slope is not None:   # linear Reddy-measured scaling
        return _GAMMA_BULK_EV + slope * (temperature_C - T_REF_C)
    global _F_REF
    if _F_REF is None:
        _F_REF = _holstein_factor(T_REF_C + 273.15)
    return _GAMMA_BULK_EV * _holstein_factor(temperature_C + 273.15) / _F_REF


def _drude_damping_delta(wavelength_nm, gamma_new_ev):
    """Delta-eps that retunes the Drude damping from _GAMMA_BULK_EV to gamma_new_ev."""
    E = 1239.84 / np.asarray(wavelength_nm, dtype=float)   # photon energy (eV)
    inv = lambda g: 1.0 / (E**2 + 1j * E * g)
    return _OMEGA_P_EV**2 * (inv(_GAMMA_BULK_EV) - inv(gamma_new_ev))


def size_damping_correction(wavelength_nm, diameter_nm, A_surf=1.0):
    """Delta-epsilon added to bulk eps for surface-scattering damping
    (gamma_S = A*hbar*v_F/R — temperature-INDEPENDENT)."""
    g_size = A_surf * _HBAR_VF_EVNM / (diameter_nm / 2.0)
    return _drude_damping_delta(wavelength_nm, _GAMMA_BULK_EV + g_size)


def thermal_damping_correction(wavelength_nm, temperature_C):
    """Delta-epsilon for BULK heating: gamma_bulk(T_REF) -> gamma_bulk(T).

    Zero at T_REF_C by construction. This is limitation #11's fix: heating
    raises eps2 -> the plasmon broadens, drops and redshifts with NO
    aggregation. Implemented from electron-phonon theory (Holstein, Th=170 K)
    anchored at Olmon's 44 meV — NOT from a fitted eps(T) table; stated
    explicitly per the task spec. Validated against Reddy 2016's measured
    Drude broadening slope (~12% agreement near room T).
    """
    return _drude_damping_delta(wavelength_nm, gamma_bulk_ev(temperature_C))


# --- full Reddy DCP eps(lambda, T) — interband T-dependence included --------
# Reddy et al., OME 6, 2776 (2016), Table S1 (arXiv:1604.00064 SI p.33):
# 200-nm POLYCRYSTALLINE film, first cycle. Drude + 2 critical points,
# phases fixed at -pi/4:
#   eps = eps_inf - wp^2/(w^2 + i*Gd*w)
#         + sum_j C_j*E_j*( e^{-i pi/4}/(E_j - w - i g_j)
#                          + e^{+i pi/4}/(E_j + w + i g_j) ),  w in eV.
# Unlike the Drude-only thermal retune (thermal_damping_correction), this
# also moves the INTERBAND parameters with T (eps_inf 2.27->2.45 and the
# 2.62 eV CP damping g2 0.256->0.273 already over 23->100 C). Parameters are
# interpolated linearly in T (extrapolated slightly for 15 C < 23 C).
# NB: an evaporated-film eps — absolute peak position inherits film bias
# (like 'bb'/'yakubovsky'); use it to quantify the T-RESPONSE, not for
# absolute peak calibration.
_REDDY_P200 = dict(
    T=np.array([23.0, 100.0, 200.0, 300.0, 400.0, 500.0]),
    eps_inf=np.array([2.27, 2.45, 2.04, 2.06, 2.03, 1.90]),
    wp=np.array([8.856, 8.863, 9.113, 9.012, 8.978, 8.959]),
    gd=np.array([0.0471, 0.0489, 0.0496, 0.0678, 0.0695, 0.0815]),
    C1=np.array([2.31, 2.28, 2.64, 2.49, 2.56, 2.54]),
    g1=np.array([1.215, 1.198, 1.464, 1.357, 1.407, 1.369]),
    E1=np.array([3.082, 3.060, 3.270, 3.218, 3.208, 3.209]),
    C2=np.array([0.226, 0.224, 0.395, 0.371, 0.375, 0.394]),
    g2=np.array([0.256, 0.273, 0.359, 0.375, 0.398, 0.429]),
    E2=np.array([2.625, 2.620, 2.631, 2.607, 2.599, 2.578]),
)


def _reddy_param(key, temperature_C):
    """Linear interp in T; linear extrapolation below 23 C (first segment)."""
    T, v = _REDDY_P200["T"], _REDDY_P200[key]
    t = float(temperature_C)
    if t < T[0]:
        return float(v[0] + (v[1] - v[0]) / (T[1] - T[0]) * (t - T[0]))
    return float(np.interp(t, T, v))


def gold_epsilon_reddy(wavelength_nm, temperature_C=23.0, drude_only=False):
    """Full Reddy DCP eps(lambda, T) for 200-nm polycrystalline gold.

    drude_only=True freezes the interband (eps_inf, wp, both CPs) at the
    23 C row and lets ONLY Gamma_D move with T — isolates the interband
    contribution to the thermal response by difference.
    """
    lam = np.asarray(wavelength_nm, dtype=float)
    w = 1239.84 / lam
    tP = 23.0 if drude_only else temperature_C
    gd = _reddy_param("gd", temperature_C)
    eps = (_reddy_param("eps_inf", tP)
           - _reddy_param("wp", tP) ** 2 / (w**2 + 1j * gd * w))
    phi = -np.pi / 4
    for j in ("1", "2"):
        C, g, E0 = (_reddy_param(k + j, tP) for k in ("C", "g", "E"))
        eps = eps + C * E0 * (np.exp(1j * phi) / (E0 - w - 1j * g)
                              + np.exp(-1j * phi) / (E0 + w + 1j * g))
    return eps


def gold_epsilon_sized(wavelength_nm, diameter_nm, A_surf=1.0,
                       temperature_C=None):
    """Size- and temperature-corrected gold permittivity.

    Single retune with gamma_new = gamma_bulk(T) + gamma_S — the correct
    composition (the two deltas are NOT additive; inv() is nonlinear).
    """
    g_size = A_surf * _HBAR_VF_EVNM / (diameter_nm / 2.0)
    g_new = gamma_bulk_ev(temperature_C) + g_size
    return gold_epsilon(wavelength_nm) + _drude_damping_delta(wavelength_nm, g_new)


def gold_index(wavelength_nm, diameter_nm=None, A_surf=1.0,
               temperature_C=None):
    """Complex refractive index n + ik of gold at temperature_C.

    If diameter_nm is given, applies the (T-independent) small-particle
    damping correction on top of the thermal one.
    """
    if diameter_nm is None:
        return np.sqrt(gold_epsilon(wavelength_nm, temperature_C=temperature_C))
    return np.sqrt(gold_epsilon_sized(wavelength_nm, diameter_nm, A_surf,
                                      temperature_C))


# --- media ---------------------------------------------------------------
WATER_N = 1.333          # legacy fixed value (~22 C, 589 nm)
GLYCEROL_N = 1.47

# Water n(T) at 589 nm (CRC Handbook); dn/dT dispersion across 420-800 nm is
# small (<10% of dn/dT) and neglected. n falls ~1.3334 (15 C) -> ~1.3240 (75 C):
# a BLUE shift of the plasmon on heating — opposite in sign to the observed
# red-shift (limitation #13: ignoring it MASKS part of the real signal).
_WATER_N_T = np.array([
    (10, 1.33369), (15, 1.33339), (20, 1.33299), (25, 1.33251), (30, 1.33192),
    (40, 1.33051), (50, 1.32894), (60, 1.32718), (70, 1.32511), (80, 1.32287),
])

# Water absorption k(lambda), 400-900 nm. Sources: Pope & Fry, Appl. Opt. 36,
# 8710 (1997) for 400-700 nm; Hale & Querry, Appl. Opt. 12, 555 (1973) above
# 700 nm (the ~740-760 nm O-H overtone bump seen in the blank scans).
# k = a*lambda/(4*pi); values are ~1e-10..5e-7 — utterly negligible for the
# PARTICLE cross sections (and blank-referenced spectra cancel the bulk path
# absorption), but tabulated here so the blank/reference channel can be
# modeled. NOT wired into Mie: an absorbing host makes the extinction cross
# section ill-defined (active research); use for diagnostics only.
_WATER_K = np.array([
    (400, 2.11e-10), (450, 3.30e-10), (500, 8.12e-10), (550, 2.47e-9),
    (600, 1.06e-8), (650, 1.76e-8), (700, 3.48e-8), (725, 8.7e-8),
    (750, 1.56e-7), (775, 1.48e-7), (800, 1.25e-7), (850, 2.93e-7),
    (900, 4.86e-7),
])


def water_index_T(temperature_C):
    """Real refractive index of water vs temperature (CRC 589 nm table)."""
    t, n = _WATER_N_T.T
    return float(np.interp(temperature_C, t, n))


def water_k(wavelength_nm):
    """Imaginary index k of water (tabulated, 400-900 nm; see _WATER_K note)."""
    wl_t, k_t = _WATER_K.T
    return np.interp(np.asarray(wavelength_nm, dtype=float), wl_t, k_t)


def medium_index(name_or_value="water", temperature_C=None):
    """Refractive index of the medium, as a COMPLEX number.

    - numeric input: returned as complex(value) (back-compat pass-through).
    - "water": real part from the CRC n(T) table when temperature_C is given
      (default None -> legacy fixed 1.333); imaginary part 0 — water's k(λ)
      is available separately via water_k() (see _WATER_K for why it is not
      folded into the particle optics).
    Callers doing propagation math should take .real explicitly.
    """
    if isinstance(name_or_value, (int, float, complex)):
        return complex(name_or_value)
    if name_or_value == "water" and temperature_C is not None:
        return complex(water_index_T(temperature_C))
    return complex({"water": WATER_N, "glycerol": GLYCEROL_N,
                    "vacuum": 1.0}[name_or_value])
