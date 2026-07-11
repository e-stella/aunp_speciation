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


_MODELS = ("etchegoin", "bb", "jc", "yakubovsky25", "yakubovsky53")


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


def gold_epsilon(wavelength_nm, model=None):
    """Complex relative permittivity of gold at vacuum wavelength(s) in nm.

    model overrides the module default set by use_gold_model().
    """
    lam = np.asarray(wavelength_nm, dtype=float)
    m = model or _GOLD_MODEL
    if m == "bb":
        wl_t, eps_t = _load_bb()
        return np.interp(lam, wl_t, eps_t.real) + 1j * np.interp(lam, wl_t, eps_t.imag)
    if m == "jc":
        nk = _jc_spline()(lam)
        return (nk[0] + 1j * nk[1])**2
    if m in ("yakubovsky25", "yakubovsky53"):
        nk = _yk_spline(m[-2:])(lam)
        return (nk[0] + 1j * nk[1])**2
    if m == "yakubovsky":
        use_gold_model(m)  # raises with the size-matching guidance
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
