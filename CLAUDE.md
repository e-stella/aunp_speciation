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
2. **Gold dielectric: `yakubovsky` IMPLEMENTED and measured — ε is a real but
   SMALL red-tail confound, and a LARGE speciation confound.** Models now:
   `etchegoin`, `bb`, `jc`, `yakubovsky25`, `yakubovsky53` (Yakubovsky 2017
   thin-film n,k, 300–2000 nm, `data/gold_yakubovsky.npz`, CC0 via
   refractiveindex.info; cubic-spline like `jc`). Ref-18 size-matching rule is
   explicit: `gold_model_for_diameter(D)` → 25 nm film for R<25 nm, else 53 nm
   film; bare `'yakubovsky'` deliberately raises (no silent film switching).
   Measured at D=12.9 in water (peak | 700/peak | 790/peak):
   etchegoin 526.9|.0090|.0045 · bb 527.0|.0188|.0091 · jc 522.8|.0124|.0066 ·
   yk25 515.6|.0166|.0079 · yk53 515.8|.0122|.0059. So: **yk25 raises the
   monomer red tail vs jc (+34% @700, +20% @790) — but that is only ~1–2% of
   the measured pedestal (#8): the red excess is NOT an ε artifact.** Fitted on
   C500, switching jc→yk25 moves A_sca by only ~−5% and n_sca by ~0.1. BUT
   yk25/53 put the peak 7.5 nm BLUE of experiment (evaporated-film ε₁ ≠
   colloidal Au), and in the global fit that bias masquerades as aggregation
   (agg jumps to a non-credible 0.7–0.8 with lower RMS — more-flexible-but-
   wrong). ⇒ **keep `jc` for the 13 nm C500 fits; use yakubovsky as the
   ε-systematics bound and for ≥50 nm work (where ref 18's validation lives).**
   Also: a size-damping correction exists (`gold_epsilon_sized`) but is NOT yet
   wired into the size-integral of `spectra.species_basis` (mie has
   `size_correction=` flag, default off) — wire it so each size in the
   polydispersity integral uses its own ε. Big effect: FWHM 79→117 nm at 12 nm.
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
8. **Scattering pedestal: MODELED — and it is NOT Rayleigh.** The C500 red wing
   is a near-flat pedestal (700/peak 0.18→0.21, 790/peak 0.156→0.186, growing
   with T) that no monomer/dimer/trimer basis reaches (exact T-matrix dimer
   ~0.021 @700; order of magnitude short at 790). IMPLEMENTED: the forward
   model in `fit_global.py` (and `fitting.py`) now includes
   A_sca,T·(λ/550)^(−n_sca): per-T amplitudes profiled by NNLS together with
   the global amplitude (linear, ≥0), ONE shared exponent fitted in [0,6].
   This replaces (deprecates) the ad-hoc per-T 780–800 nm subtraction.
   **Findings on C500 (exact jc basis):** n_sca = 0.73±0.01 (heating) /
   0.79±0.03 (cooling); ε-robust (yk25: 0.63–0.65) — **decisively NOT n≈4
   Rayleigh**: the pedestal decays too weakly for small aggregates and points
   to scatterers ≳λ (µm-scale flocs; large-particle Mie flattens toward λ⁰),
   and/or partly instrumental stray light. A_sca ≈ 0.21–0.22 of series max
   @550, grows ~+2–3% over 15→75 °C and is REVERSIBLE across branches
   (hysteresis ~1–2.5%, cooling slightly above heating) — supporting a
   physical scatterer over pure drift. Fit quality: RMS 0.125 (no pedestal) →
   0.037 (ad-hoc subtraction) → 0.0225 (fitted pedestal). Speciation: agg
   small-cluster gold ≈ 0.41–0.47, ~T-stable/mildly decreasing; D = 11.2±0.7
   (heating) / 12.2±0.6 (cooling). CAVEATS: (a) with n≈0.7 the pedestal
   extends under the plasmon peak, so it partially trades against D/poly — D
   sits ~2σ below TEM 12.9; pin size via absolute extinction/Haiss (#10), not
   the normalized shape. (b) Guardrails passed: on clean synthetics A_sca
   returns ≈0 and old parameters are recovered (matched-backend caveat:
   `fit_spectrum` now takes `backend=` — an optics mismatch between data and
   basis gets absorbed by D AND fake a_sca). (c) See the fit_global.py
   docstring for the ≥60 nm design note: the species basis keeps its own
   Mie/T-matrix scattering; A_sca only absorbs the un-enumerated population.
9. **[FIXED] Monomer peak ~4 nm red of experiment.** Root cause was the gold
   dielectric DATASET, not the medium: tabulated Johnson & Christy (model
   `'jc'`, embedded in `dielectric.py`, cubic-spline interpolated — linear
   interpolation pins the peak to the sparse J&C grid nodes) puts the 12.9 nm
   monomer at 522.8 nm in bare water (target 523.1). No ligand-shell index
   needed. NB: the Brendel-Bormann table (`'bb'`) is ~4–6 nm red-biased in the
   500–540 nm window just like Etchegoin — do not use either for calibration.
   `fit_real.py` sets `'jc'`; caches record their `gold_model` and `fit_real`
   refuses a mismatched cache.
10. **Single-spectrum size retrieval fails below ~35 nm — INDEPENDENTLY CONFIRMED,
    and the red tail biases even the Haiss rescue.** Validated on the CTAC size
    series (5 monomer samples w/ paired TEM: 7.8 / 29.2 / 31.5 / 37.9 / 42.0 nm):
    shape-only Mie fit returns 15.7 / 19.7 / 20.2 / 30.7 / 40.2 nm — i.e. only the
    42 nm sample is recovered. This is an INFORMATION limit, not a bug: normalized
    small-particle spectra are near size-degenerate (Haiss ref 14 states their own
    peak-position method is valid only ≥25 nm). The **Haiss A_spr/A_450 ratio only
    partly helps**: it separates 8 nm from ~30 nm but SATURATES across 29–42 nm
    (~2.39–2.43), and measured ratios sit ABOVE clean-monomer Mie because the
    excess red tail (#8) inflates A_spr — so the red tail contaminates the size
    readout too. ⇒ Size must be pinned by the TEMPERATURE SERIES (absolute
    extinction + van 't Hoff), not a single normalized spectrum.
    **Also validated (good news):** the Mie scattering channel is sound —
    σ_sca/σ_ext at peak rises 0.4% (8 nm) → 6.7% (42 nm), the expected
    (size)⁶-vs-(size)³ scaling. For planned ~60 nm work scattering is 15–30% of
    extinction and even monomers scatter, which will make the pedestal/species
    separation (#8) harder than it is at 13 nm — design the scattering term for
    that now.
    **Red-tail excess is reproducible in clean samples:** with D FIXED to TEM
    ground truth, data ride above Mie across ~560–660 nm in every CTAC sample,
    including the pristine 29/31 nm ones — independently reproducing the
    persistent-red-tail-vs-simulation effect reported in the lab.

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
0. **Add `yakubovsky` size-matched thin-film ε to `dielectric.py`** (limitation
   #2, ref 18) and re-check the monomer peak + red tail. Cheap, and it RULES OUT
   an ε artifact as the source of the static red-side excess before #1 attributes
   that excess to a scattering population. Do this BEFORE trusting #1's amplitude.
1. **Scattering term: DONE** (see #8; n_sca≈0.7, not Rayleigh). Follow-ups:
   (a) replace the phenomenological λ^-n with a real large-aggregate Mie
   species (a few 100 nm–1 µm effective sizes) to test the floc hypothesis
   against the fitted n; (b) revisit the mult_400nm anchor — the pedestal
   lifts 400 nm too, so the anchor slightly over-corrects concentration;
   (c) the pedestal–(D, poly) degeneracy leaves D ~2σ low vs TEM: pin size
   with absolute extinction/Haiss (#10).
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
