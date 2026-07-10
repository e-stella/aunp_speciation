"""Generate an example experimental-style temperature-series CSV to demo the
real-data driver. Writes data/example_series.csv.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from aunp_speciation.spectra import species_basis
from aunp_speciation.equilibrium import association_constants, gold_fractions
from aunp_speciation.io_data import save_series

rng = np.random.default_rng(7)
DATA = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATA, exist_ok=True)

SPECIES = ("monomer", "dimer", "trimer_linear")
K = {"monomer": 1, "dimer": 2, "trimer_linear": 3}
D, P, GAP, CTOT = 12.0, 4.0, 1.0, 0.01
TH = dict(dH2=-25.0, dS2=-0.06, dH3=-45.0, dS3=-0.14)
temps_C = np.array([5., 20., 35., 50., 65.])
# stay within the cached wavelength grid so data & fit use consistent optics
wl = np.arange(470, 672, 2.0)

# use the exact cached optics if available (so fit_real, which prefers the cache,
# is validated self-consistently); else CDA
CACHE = os.path.join(os.path.dirname(__file__), "..", "outputs", "tmatrix_basis.npz")
if os.path.exists(CACHE):
    from aunp_speciation.basis_cache import load_cache
    backend = load_cache(CACHE).species_fn
    print("example data uses EXACT cached optics")
else:
    backend = "cda"
    print("example data uses CDA optics (cache not found)")
basis = species_basis(wl, D, P, GAP, "water", species=SPECIES, backend=backend, n_sizes=7)
B = np.column_stack([basis[s] for s in SPECIES])
kvec = np.array([K[s] for s in SPECIES], float)
canon = lambda s: "trimer" if s.startswith("trimer") else s

spectra = []
for Tc in temps_C:
    K2, K3 = association_constants(Tc + 273.15, **TH)
    gf = gold_fractions(CTOT, K2, K3)
    f = np.array([gf[canon(s)] for s in SPECIES])
    spectra.append(B @ (f / kvec))
spectra = np.array(spectra)
spectra = spectra / spectra.max()
spectra += rng.normal(0, 0.004, spectra.shape)  # measurement noise

save_series(os.path.join(DATA, "example_series.csv"), wl, spectra, temps_C)
print("wrote data/example_series.csv  (truth: D=12 nm, poly=4%, dH2=-25, dS2=-0.06)")
