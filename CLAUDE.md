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
   UPDATE: the size-damping correction IS now wired through the whole basis —
   `species_basis(..., size_correction=True, A_surf=)` gives each size in the
   polydispersity integral its own γ_S ε, and both cluster backends accept it
   (the T-matrix cache bakes it in and records the flag). Required for ε(T)
   fits (#11). Remaining knob: A_surf is assumed 1.0 (matches Chetoui's
   γ_S = 0.143 eV at D=12.9) but trades against D/poly in fits — with it the
   fitted D falls to ~10.5–10.8 vs TEM 12.9; consider calibrating A_surf on a
   TEM-pinned sample.
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
   **EXCLUSIONS CLOSED (what the pedestal is NOT):**
   (i) NOT small-aggregate Rayleigh — n_sca = 0.73/0.79, not 4 (above).
   (ii) NOT a permittivity *dataset* artifact — yakubovsky raises the model tail
        but accounts for only ~1–2% of the measured pedestal (#2).
   (iii) NOT instrumental / NOT a thermal-optical artifact of the cell. The 2011
        `raw` sheet holds **2–4 water-blank scans at EVERY temperature**, interleaved
        before each sample scan. **Blanks are flat with T:** far-red level drifts only
        −0.0016 (700 nm) / −0.0020 (790 nm) over 15→75 °C — ~15× smaller than the
        sample pedestal's ≈ +0.03 growth, and in the OPPOSITE direction. Excludes,
        for the *T-dependent* component: stray light, convection/schlieren, and
        gas-bubble formation from dissolved gas (pure water carries the same gas load
        and shows nothing). NB a *static* ~0.03 instrumental offset IS present in the
        far red and contributes to the pedestal's absolute baseline — not its T-growth.
   (iv) NOT nanobubbles on the particles (lit-map §F): nanobubbles require
        **hydrophobic** surfaces (citrate/TEG AuNP are deliberately hydrophilic), high
        curvature **suppresses** them on gold, observed surface nanobubbles are
        120–145 nm (too large for a 13 nm sphere), and heating to 40–60 °C *reduces*
        their number. The *laser-driven* plasmonic-nanobubble literature is a different
        regime and must not be conflated with mild bulk heating.
   **(v) FIFTH EXCLUSION MEASURED — ε(T) (see #11) — BUT THE RESIDUAL-PEDESTAL
   CHARACTERIZATION IS ⚠️ PROVISIONAL (see #14: evaporation).** Under ε(T) +
   γ_S basis + Kell density norm: A_sca drops ~0.22 → 0.13–0.17 of series max
   (≈40% of the old pedestal was gold-physics + basis-width artifact — this
   part stands), and the surviving exponent collapses n_sca = 0.73 →
   0.00±0.06, a flat floor growing +31% over 15→75 °C. **CAUTION: a
   wavelength-flat, T-growing brightness term is exactly the signature of the
   Kell normalization's evaporation blindness (#14) — Kell over-brightens
   high-T spectra by ~+2.5%(45 °C) to ~+4.2%(75 °C). The n_sca=0.00 landing
   exactly on the bound, and its move 0.73→0.00 when the normalization
   changed, are red flags that the pedestal term is absorbing a normalization
   error. Do NOT treat the flat pedestal or its +31% growth as physical until
   an evaporation-aware concentration term is in.** Floc-test result (stands
   as a shape statement): 0.7–1 µm Maxwell-Garnett–Mie spheres reproduce a
   flat pedestal unimposed (equiv n = 0.00/−0.27, shape-RMS 0.088–0.12);
   ≤200 nm flocs give n ≥ +1.3. Next discriminator: evaporation-aware
   normalization first, then DLS/filtration if a pedestal still survives.
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

11. **[IMPLEMENTED — AND THE VERDICT IS LARGELY AGAINST THE OLD READING] Gold ε(λ,T).**
    `gold_epsilon(..., temperature_C=)` now applies a bulk thermal Drude-damping
    retune on top of any base dataset: γ_bulk(T) from the **Holstein electron–phonon
    form (Θ_D = 170 K), anchored at Olmon's 44 meV @ 20 °C** — implemented from
    theory, stated explicitly (the Reddy 2016 measured single-crystal slope
    ~1.1e-4 eV/K agrees to ~12%; their Drude+2CP tables were obtained from the
    open-access full text, and the author list is verified: Reddy, Guler,
    Kildishev, Boltasseva, Shalaev, *OME* 6, 2776). γ_S stays T-independent per
    Chetoui and composes with γ_bulk(T) in a SINGLE Drude retune
    (`gold_epsilon_sized(..., temperature_C=)`) — the deltas are not additive.
    Threaded through mie/clusters/clusters_tmatrix/spectra/fitters; the T-matrix
    cache gained a temperature axis ([15,45,75] °C, linear interp) and records
    `eps_t_model`; `fit_real.py` refuses stale caches; `--fixed-eps` reproduces
    legacy fits. Drude baseline updated to Olmon (ω_p=8.45 eV, γ=44 meV; was
    9.0/0.07 textbook values).
    **MEASURED (D=12.9, jc, 15→75 °C; observed: peak −2.7%, @700 +15.6%, @790
    +19.1%):** γ_bulk rises 20.1% (Holstein) → without γ_S the basis responds
    peak −1.7%, @700 +10.5%, @790 +10.9% (~57–67% of everything); WITH γ_S the
    physical response is peak −1.3%, @700 +3.7%, @790 +3.8% (~20–47%). The
    napkin's "~40%" was directionally right; the exact share depends on the γ_S
    dilution it itself identified. **γ_S must be in the basis for ε(T) fits**
    (`size_correction=True`, now wired through the whole basis incl. the
    polydispersity integral — closes #2's open item) or the fit over-attributes
    T-changes to gold ~3×. Water n(T) adds peak −2.85%/blue −0.4 nm by itself.
    **REFIT VERDICT — ⚠️ PROVISIONAL, pending evaporation-aware concentration
    correction (#14):** with ε(T)+n(T) in the model and the Kell density
    normalization, ΔH₂ collapses to ≈0 (+0.9±1.2 heating, +4.1 cooling), the
    aggregated fraction is ~T-flat (0.24 heating / 0.10 cooling), and the
    thermally-driven speciation signal essentially vanishes — BUT the Kell
    normalization is blind to the measured ~5.6% evaporation over the run
    (#14), which injects a spurious wavelength-flat, T-growing brightness
    error of exactly the kind the pedestal term absorbs. The claim "the
    thermally-driven speciation does not survive" is therefore NOT settled;
    re-decide after an evaporation-aware (time-monotonic) concentration term.
    What stands regardless: the ε(T)/γ_S/n(T) physics, its measured basis
    response, and the ~3× γ_S-dilution requirement. RMS improves
    0.0225→0.0166 (fit quality, not verdict). **Honest caveat: the
    ε(T)-ON and ε(T)-OFF fits reach the SAME RMS (0.0166 vs 0.0164)** — the
    data cannot statistically distinguish "gold heats" from "speciation
    shifts"; ε(T) belongs in the model because it is mandatory physics, not
    because the fit prefers it. Also: with γ_S(A_surf=1) the fitted D falls to
    10.5–10.8 (below TEM 12.9 and the cache grid) — a real A_surf/D/poly
    degeneracy; pin size via absolute extinction (#10) and consider
    calibrating A_surf.
12. **[IMPLEMENTED] Normalization re-decided: `normalize="density"` (Kell 1975).**
    `load_series` now supports a pure dilution correction
    A·ρ(T_ref)/ρ(T) (Kell polynomial, verified ρ(15)=999.10, ρ(75)=974.85 —
    NB the denominator coefficient is 16.87985e-3; a digit transposition was
    caught against these reference values). No anchor-wavelength assumption;
    needs parseable temps in every column header. **`mult_400nm` is DEPRECATED/
    BIASED** (premise measurably violated: A(400) rises +3.18% while Kell
    predicts −2.43%; flat rescale inflates the apparent peak change ~8.5×);
    kept only for reproducing old results. `fit_real.py --normalize density`.
    **ISOSBESTIC RE-TEST — ⚠️ PROVISIONAL, pending evaporation-aware correction
    (#14).** Under the Kell density normalization the speciation residual has
    no common zero-crossing (small-ΔT curves show dozens of noise-level
    crossings; ΔT≥20 °C curves NONE in 500–700 nm). BUT this was tested ONLY
    under Kell, which is blind to the measured ~5.6% evaporation; under a
    correct multiplicative concentration correction, crossings DO exist and
    drift 562 → 577 nm. So "the 575 nm isosbestic does not exist" is NOT
    established either way — what IS established is that its location/existence
    is normalization-dependent, so it cannot carry the two-state argument
    until the concentration history is pinned down.
13. **[IMPLEMENTED] Medium n(T) + water k(λ).** `medium_index(name, temperature_C)
    -> complex`; water n from the CRC 589 nm table (1.3334@15 → 1.3229@80 °C;
    dn/dT dispersion neglected), threaded through all optics. Water k(λ)
    tabulated 400–900 nm (`water_k()`; Pope & Fry 1997 + Hale & Querry 1973,
    incl. the 740–760 nm O–H bump seen in blanks) but NOT wired into Mie —
    k ≤ 5e-7 is negligible for particle cross sections, absorbing-host Mie is
    ill-defined, and blank-referenced spectra cancel the bulk path absorption.
    **Prediction verified:** freezing n(T) moves D/ΔH₂ slightly (it had been
    masking part of the T-signal — n(T) alone is peak −2.85%, blue −0.4 nm)
    and leaves A_sca essentially untouched (0.131→0.167 vs 0.132→0.172) — no
    loud flag needed; nothing is misunderstood there.
14. **[CRITICAL — INVALIDATES THE #8/#11/#12 VERDICTS AS STATED] `normalize=
    "density"` is blind to EVAPORATION; a time-monotonic concentration term is
    required.** Branch-offset evidence (raw data, heating vs cooling at MATCHED
    temperatures): the cooling branch is systematically BRIGHTER — at 15 °C
    +5.63% @400, +5.65% @523, +5.95% @700, +6.78% @790 nm; shrinking to +1.9%
    at 65 °C. Conclusive because (a) the offset is nearly wavelength-
    independent in percent ⇒ a pure multiplicative CONCENTRATION change, not
    chemistry; (b) it scales with the TIME GAP between scans (heating-15 °C is
    the first scan, cooling-15 °C the last ⇒ largest offset; the 65 °C scans
    are adjacent ⇒ smallest). The sample concentrated ~5.6% over the run by
    evaporation. Kell corrects only thermal expansion (concentration DOWN
    2.4%) while the true concentration went UP ⇒ Kell OVER-BRIGHTENS high-T
    spectra by ~+2.5% (45 °C) to ~+4.2% (75 °C), monotonically in T. A
    spurious wavelength-flat brightness error growing with T is EXACTLY what
    the fitted "flat pedestal (n_sca=0.00) growing +31%" looks like — and
    n_sca moving 0.73 (anchor norm) → exactly 0.00 (Kell) when only the
    normalization changed is itself a red flag that the pedestal term absorbs
    normalization error. FIX (not yet implemented): an evaporation-aware
    multiplicative concentration model, e.g. c(t) monotonic in SCAN TIME/ORDER
    (identifiable from the matched-T branch offsets), composed with Kell;
    then re-run the #8/#11/#12 verdicts.

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
  T-matrix basis over a (T, diameter, gap) grid -> outputs/tmatrix_basis.npz
  (~19 min at lmax=6: T [15,45,75] °C, D 11–15, gaps 0.5–3, wl 420–800,
  gold_model='jc', size_correction=True baked in; cache records eps_t_model
  and the flags — fit_real refuses mismatches).
- `python scripts/make_example_data.py` — write data/example_series.csv (uses the
  cache if present, so the demo is optics-consistent).
- `python scripts/fit_real.py [file.csv] [--range MIN MAX] [--normalize
  {none,density,mult_400nm}] [--anchor NM] [--fixed-eps]` — fig 8; load a real
  UV-Vis file and fit with gold ε(T) + water n(T) per temperature (--fixed-eps
  for legacy fixed-ε fits). Sets gold model 'jc'; uses the cached EXACT optics
  if present (refuses a cache with mismatched gold_model / eps_t_model / no T
  axis), else CDA with γ_S ON. Validated: on the self-consistent example
  (jc + ε(T) + γ_S basis) it recovers D=11.91 nm (vs 12), ΔH₂=-26.2 (vs -25),
  RMS at the noise floor, A_sca≈0 with n_sca flagged unidentified.
  Real C500 series: `python scripts/fit_real.py
  experimental/ESK_2011/aunp_heating_RAW_390-900.csv --normalize density` —
  see #11's verdict before interpreting the speciation numbers.

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

**A–D + F: code DONE; the scientific verdicts are ⚠️ PROVISIONAL pending #14
(evaporation).** Under Kell the fits said: thermally-driven speciation gone
(ΔH₂≈0), pedestal flat (n=0) growing +31%, no isosbestic — but Kell is blind
to the measured ~5.6% evaporation, whose signature (wavelength-flat,
T-monotonic brightness) is degenerate with all three claims. What stands
regardless: the ε(T)/γ_S/n(T) physics and infrastructure, the measured basis
T-response, the γ_S-dilution requirement, and the floc shape result.

**A0. [NEW — DO FIRST] Evaporation-aware concentration correction (#14):**
   model c(t) as monotonic in scan time/order, identified from the matched-T
   heating-vs-cooling offsets (+5.6% @15 °C → +1.9% @65 °C), composed with
   Kell; then RE-RUN the #8/#11/#12 verdicts (speciation survival, isosbestic
   562→577 nm drift, pedestal shape/growth).

**E. Break the A_surf–D–poly degeneracy (#2, #10, #11):** pin size with ABSOLUTE
   extinction + Haiss A_spr/A_450 (not the normalized shape), and calibrate
   A_surf on a TEM-pinned monomer sample. With γ_S(A_surf=1) in the basis the
   fitted D rails to 10.5–10.8 vs TEM 12.9 — the biggest remaining model knob.

**G. EXPERIMENTAL (lab, not code) — now THE decisive test:** DLS across the ramp
   and/or a filtered/centrifuged control. The fit cannot distinguish µm-flocs
   from stray light (both are λ-flat); only the T-reversible growth + flat
   blanks argue for a physical scatterer. No amount of fitting settles it.

**H. Model-selection honesty:** ε(T)-ON and ε(T)-OFF fits tie on RMS (0.0166 vs
   0.0164). If a referee asks "how do you know it's the gold and not the
   particles", the answer is "because bulk gold measurably does this (Reddy),
   at exactly the magnitude used" — an a-priori argument, not a fit-quality
   one. Keep it that way; do not tune ε(T) parameters against the C500 data.

I. Optional: NN surrogate trained on Layer-1+2 spectra for instant inference.
