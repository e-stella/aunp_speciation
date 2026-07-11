# CLAUDE.md — AuNP Speciation / UV-Vis broadening model

Context for AI coding sessions (Claude Code etc.). Read this first.

## The scientific problem
12 nm gold nanoparticles at ~4% polydispersity show a UV-Vis extinction peak
that is **much broader and more red-tailed** than an FDTD/Mie model accounting
only for the 4% size spread. Hypothesis (Prof. Dragnea's lab): in liquid the
particles form a reversible equilibrium of **monomers ⇌ dimers ⇌ trimers ⇌
multimers**; the UV-Vis is an ensemble/time-integrated snapshot, and plasmonic
coupling in the clusters adds red-shifted intensity. Supporting evidence: a
**temperature-dependent isosbestic point** in the measured spectra.

**Goal:** Monte-Carlo the cluster populations → compute a simulated ensemble
spectrum → compare to experiment → invert into a tool that reports nominal
size + size distribution (%) AND multimer fractions, to replace/augment TEM.
This reframes UV-Vis width as *polydispersity + speciation*, not size alone.

## Architecture (three layers)
- **Layer 1 — `equilibrium.py`**: monomer/dimer/trimer association with
  van 't Hoff K(T); produces temperature-dependent populations → isosbestic.
- **Layer 2 — optical forward model** (two selectable backends):
  - `dielectric.py`: gold ε(λ), Etchegoin analytic model + optional small-particle
    surface-scattering damping (`gold_epsilon_sized`, `size_damping_correction`).
    Per ref 18, the imaginary part (damping) dominates fit accuracy and tabulated
    thin-film constants (Yakubovsky 2017) generalize best for colloidal Au — swap
    those in for production. NOTE: at 12 nm the size-damping correction alone
    broadens the monomer FWHM ~79→117 nm, so speciation broadening must be judged
    against a monomer model that already includes it (else over-attribution).
  - `mie.py`: single-sphere Mie cross sections + dipole polarizability α from a₁
    (`size_correction=True` applies the damping correction).
  - `clusters.py`: **CDA backend** — coupled point-dipole dimer/trimer, orientation
    averaged. Fast; UNDER-couples at contact (lower bound).
  - `clusters_tmatrix.py`: **exact T-matrix backend** (MSTM/GMM via `treams`) —
    all multipoles, exact orientation average. Validated: reproduces Mie for a
    monomer to 1e-15; gives 1.6× (dimer) / 2.8× (trimer) more red-tail than CDA
    at gap 1 nm; converged by lmax≈8 (sub-0.5 nm gaps need more). Needs the
    `mstm-env` venv (numpy 1.26 / scipy 1.11). Select via
    `spectra.species_basis(..., backend="tmatrix")`.
  - `spectra.py`: polydispersity integration + population-weighted mixing;
    `backend=` chooses the optics engine.
- **Layer 3 — `fitting.py` (skeleton built)**: separable/variable-projection
  least squares. Nonlinear params (mean diameter, % polydispersity) set the
  species basis shapes; species mixing weights are profiled out by NNLS at each
  step; outer `least_squares` optimizes the nonlinear params. Returns per-param
  uncertainties + identifiability verdicts. `scripts/fit_demo.py` is a synthetic
  recovery self-test.
  **Key finding from the self-test (real, not a bug):** a *single* 12 nm
  spectrum robustly constrains the **aggregated (non-monomer) gold fraction**
  (verdict "ok") but only weakly constrains mean size, barely constrains
  polydispersity (matches Bilén ref 19), and cannot separate dimer from trimer
  (collinear red tails). → Break these degeneracies by (a) using absolute
  extinction + the A_spr/A_450 ratio (Haiss, ref 14) to pin size/concentration,
  and (b) a **global fit across the temperature series** (shared size/gap,
  T-dependent K's) exploiting the isosbestic constraint.
- **Layer 3 — `fit_global.py` (global multi-temperature fit, BUILT)**: jointly
  fits a T-series that shares one sample — common (diameter, polydispersity,
  gap), conserved gold, one global amplitude profiled out linearly, and van 't
  Hoff K(T) driving populations. `scripts/fit_global_demo.py` self-test
  (mostly-monomer synthetic series) recovers mean diameter ≈12.7 nm (vs 12),
  thermodynamics ΔH₂≈−27 (vs −25 kJ/mol), and the aggregated-gold-fraction-vs-T
  curve to within noise. Polydispersity stays unconstrained (flat direction,
  huge SD — real, matches Bilén). Monomer-vs-aggregated split is robust; dimer-
  vs-trimer is still weak even globally → resolve with the exact-backend trimer
  *peak shape* and/or a concentration series.

## Known limitations (important — do not silently trust magnitudes)
1. **CDA under-couples at near-contact** (use the T-matrix backend for magnitudes).
   One dipole per sphere omits higher multipoles that dominate for gap/D ≲ 0.1,
   so CDA red-tails are a **lower bound** (measured 1.6×/2.8× low for dimer/
   trimer). The exact backend (`clusters_tmatrix.py`, treams) is wired in; use
   MSTM (exact multisphere
   T-matrix; Python `treams`/`smuthi`, or Mackowski MSTM) — verified to match
   CDA for a monomer and in the weak-coupling limit. This is the top priority
   upgrade for Layer 2.
2. **Analytic gold dielectric** shifts the absolute peak ~5–10 nm (currently
   527 nm vs ~520 nm expected). Swap in tabulated J&C / Yakubovsky (ref 18).
   A size-damping correction exists (`gold_epsilon_sized`) but is NOT yet wired
   into the size-integral of `spectra.species_basis` (mie has `size_correction=`
   flag, default off) — wire it so each size in the polydispersity integral uses
   its own ε. Big effect: FWHM 79→117 nm at 12 nm.
3. **Cluster polydispersity** is applied approximately (size-scaled). Refine
   with per-size cluster runs if needed.
4. Cluster geometry set is minimal (dimer, linear/triangular trimer). Extend
   with more multimer motifs + a gap distribution.
5. Global fit currently uses CDA basis in the loop (fast). For quantitative
   results precompute a T-matrix basis on a (D, gap) grid and interpolate, so
   the exact optics enter the fit without the per-iteration treams cost.
6. **Wavelength range + normalization live in the data layer (`io_data.py`).**
   The Etchegoin dielectric is only validated ~400–1000 nm and real Cary spectra
   carry deep-UV interband absorption plus a 350 nm lamp-changeover artifact, so
   `io_data.load_series` trims to `wavelength_range=(420, 800)` by default
   (`fit_real.py`: `--range MIN MAX`).
   Real UV-Vis T-series also need **concentration normalization**, applied IN
   CODE (not baked into CSVs) so the preprocessing is explicit and reproducible.
   IMPLEMENTED: `load_series(..., normalize=None, anchor_nm=400)`.
   - `normalize=None` (DEFAULT): no-op — preserves synthetic tests and
     `example_series.csv`, which are already clean and must NOT be re-normalized.
   - `normalize="mult_400nm"`: multiplicative anchor at `anchor_nm`,
     `A_norm(λ,T) = A(λ,T)·A(anchor,T_ref)/A(anchor,T)`, T_ref = first column.
     Corrects concentration/drift (Beer–Lambert multiplicative; 400 nm is the
     ~speciation-invariant dielectric point for Au). Removes a real ~3.3%
     inter-scan drift on the C500 heating series. See README_experimental_data.
   ORDER IS FIXED: normalize FIRST (anchor must still be in the array), THEN clip
   to `wavelength_range`; raises ValueError if the anchor is outside the loaded
   data. `anchor_nm` and `wavelength_range` are SEPARATE: the anchor lives in the
   raw data but deliberately outside the fit window, so the input CSV must span
   the anchor (feed the RAW 390–900 nm files, not pre-clipped).
   `fit_real.py`: `--normalize {none,mult_400nm}`, `--anchor NM`.
   Residual caveats: `fitting.py`/`fit_global.py` still weight the kept range
   uniformly — no per-wavelength noise model; and the 400 nm "speciation-
   invariant" assumption degrades if a broadband scattering pedestal (see #8)
   also lifts 400 nm.
7. **[FIXED] One-sided thermodynamic sign bounds in `fit_global.py`.** Default
   bounds used to force dH2,dS2,dH3,dS3 ≤ 0, which cannot represent thermally-
   induced aggregation (endothermic/entropy-driven, dH>0 & dS>0) — on the real
   C500 heating series the aggregated fraction fell with T (backwards). Now
   two-sided: dH∈[−150,+150] kJ/mol, dS∈[−0.6,+0.6] kJ/(mol·K); the
   `equilibrium.association_constants` docstring now covers both sign regimes.
   Synthetic self-test (exothermic truth) still recovers correctly.
8. **The real red wing is mostly a FLAT PEDESTAL no small-cluster basis can fit
   (supersedes the CDA-only framing).** Measured on C500 after mult_400nm
   normalization: 700 nm/peak ≈ 0.18→0.21 growing with T, but still 0.156→0.186
   at 790 nm — a nearly wavelength-flat pedestal ≈15–19% of peak. Best available
   basis red-tails (jc, D=12.9): CDA dimer ~0.010; EXACT T-matrix dimer 0.021
   (gap 1 nm), trimer_linear 0.068/0.023 (700/790, gap 0.5) — an order of
   magnitude short at 790 nm. The tail-minus-pedestal is ~T-constant, so the
   GROWTH of the red wing is the pedestal growing (reversibly — it tracks T on
   both heating and cooling branches): broadband scattering from large
   multimers/turbidity, NOT dimer/trimer coupling. Consequences: (a) fitting the
   raw shape rails D/poly (the fit fakes the pedestal with size); (b) subtracting
   a per-T baseline (mean 780–800 nm) before `fit_temperature_series` gives
   physical results: D = 12.5±1.5 (heating) / 12.6 (cooling) vs TEM 12.9, RMS
   0.125→0.037, aggregated gold ~0.33–0.42 roughly T-stable. NEXT STEP: add a
   scattering/baseline component to the fit model (e.g. a flat or λ^-n term, or
   a large-multimer species) instead of ad-hoc subtraction.
9. **[FIXED] Monomer peak ~4 nm red of experiment.** Root cause was the gold
   dielectric DATASET, not the medium: tabulated Johnson & Christy (model
   `'jc'`, embedded in `dielectric.py`, cubic-spline interpolated — linear
   interpolation pins the peak to the sparse J&C grid nodes) puts the 12.9 nm
   monomer at 522.8 nm in bare water (target 523.1). No ligand-shell index
   needed. NB: the Brendel-Bormann table (`'bb'`) is ~4–6 nm red-biased in the
   500–540 nm window just like Etchegoin — do not use either for calibration.
   `fit_real.py` sets `'jc'`; caches record their `gold_model` and `fit_real`
   refuses a mismatched cache.

## Conventions
- Units: lengths in nm, cross sections in nm². Wavelengths are vacuum λ₀.
- All optics done *in the medium*: size parameter/wavenumber use n_medium; the
  particle enters as relative index m = n_Au/n_medium.
- CDA and Mie are consistent for N=1 (regression-tested in `scripts/verify.py`).
- Keep everything dependency-light (numpy/scipy/matplotlib). Add `treams` when
  wiring MSTM.

## Entry points
Two interpreters: the **system** python for CDA/Mie/fits, and the **mstm-env**
venv for the exact T-matrix backend (treams needs numpy 1.26 / scipy 1.11).

- `python scripts/verify.py` — physics sanity checks.
- `python scripts/demo.py` — figs 1–4 (species, broadening, isosbestic, gap).
- `python scripts/fit_demo.py` — fig 5, single-spectrum inverse fit.
- `python scripts/fit_global_demo.py` — fig 7, global multi-temperature fit.
- `mstm-env/bin/python scripts/validate_mstm.py` — fig 6, exact-vs-CDA (needs venv).
  Recreate venv (MUST use python ≤3.12 — numpy 1.26.4 does not run on 3.13; a
  3.10 lives at /Library/Frameworks): `python3.10 -m venv mstm-env &&
  mstm-env/bin/pip install "numpy==1.26.4" "scipy==1.11.4" treams matplotlib`.
- `mstm-env/bin/python scripts/build_tmatrix_basis.py` — precompute the exact
  T-matrix basis over a (diameter, gap) grid -> outputs/tmatrix_basis.npz
  (~5.5 min at lmax=6, D 11–15, gaps 0.5–3, wl 420–800, gold_model='jc').
- `python scripts/make_example_data.py` — write data/example_series.csv (uses the
  cache if present, so the demo is optics-consistent).
- `python scripts/fit_real.py [file.csv] [--range MIN MAX] [--normalize
  {none,mult_400nm}] [--anchor NM]` — fig 8; load a real UV-Vis file and fit.
  Sets gold model 'jc'; uses the cached EXACT optics if present (fast, no treams
  at fit time; refuses a cache built with a different gold_model), else CDA.
  Validated: on the self-consistent example (jc basis) it recovers D=11.6 nm
  (vs 12, ±0.1 formal — mildly biased by n_sizes/stride differences between
  generation and fit), ΔH₂=-25.2 (vs -25), RMS at the noise floor.
  Real C500 series: `python scripts/fit_real.py
  experimental/ESK_2011/aunp_heating_RAW_390-900.csv --normalize mult_400nm` —
  but see limitation #8 (baseline pedestal) before trusting the raw-shape fit.

## The cached-exact-optics pattern (important)
The exact T-matrix is too slow inside a fit loop, so build it ONCE on a grid
(`basis_cache.build_grid`, run in mstm-env) and interpolate at fit time
(`basis_cache.load_cache(...).species_fn`, pure numpy — runs in the SYSTEM env).
Pass that callable as `backend=` to `species_basis` / the fitters. Keep the
optics backend used to *fit* the same as any used to *simulate* test data — a
mismatch (CDA data vs exact fit) biases size & aggregation (demonstrated).
Two gotchas learned the hard way:
- The interpolator EXTRAPOLATES linearly outside its grid (`fill_value=None`)
  with no warning — in wavelength AND diameter. An old 470–675 nm cache fed
  420–800 nm data silently corrupted fits (and D wandered past the grid edge).
  Build the wl grid to cover the full fit window and keep D bounds within the
  grid (pass explicit `bounds=` to the fitters).
- The cluster cubes embed the gold dielectric ACTIVE AT BUILD TIME, while the
  monomer is computed live at fit time with the current module default. Caches
  now record `gold_model`; `fit_real.py` and `make_example_data.py` enforce/
  adopt it. Rebuild the cache after changing the dielectric.

## Next steps (suggested order)
1. **Add a scattering/baseline component to the fit model** (limitation #8): a
   flat or λ^-n pedestal term (or a large-multimer species) profiled out with
   the species weights, replacing the ad-hoc per-T 780–800 nm subtraction that
   currently makes the C500 fits physical (D≈12.5 vs TEM 12.9, RMS 0.037).
   Then revisit the mult_400nm anchor (the pedestal lifts 400 nm too).
2. Wire the size-damping correction into the polydispersity integral (tabulated
   J&C 'jc' is DONE and calibrated to 522.8 nm; size-damping remains off there).
3. Add the Haiss A_spr/A_450 ratio + absolute extinction to pin size &
   concentration (breaks the size/poly degeneracy the T-series leaves open).
4. Real-data driver: DONE end-to-end (load → normalize → clip → global fit;
   limitations 6/7/9 fixed, 8 diagnosed). Remaining: promote the baseline-
   subtracted diagnostic into `fit_real.py` proper once step 1 lands.
5. Validate against TEM-characterized samples + the temperature/isosbestic
   series (first pass done: D 12.5–12.6 both branches vs TEM 12.9±7%; pedestal
   reversible with T on heating AND cooling).
6. Optional: NN surrogate trained on Layer-1+2 spectra for instant inference.
