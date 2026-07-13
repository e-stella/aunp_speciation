# CLAUDE.md — AuNP Speciation / UV-Vis broadening model

Context for AI coding sessions (Claude Code etc.). Read this first.

## The scientific problem
12 nm gold nanoparticles at ~4% polydispersity show a UV-Vis extinction peak
that is **much broader and more red-tailed** than an FDTD/Mie model accounting
only for the 4% size spread. Hypothesis (Dragnea group, Indiana University): in liquid the
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
   **WIDTH DECOMPOSITION MEASURED (2026-07-12,
   `scripts/width_decomposition.py`, fig19; CTAC per-particle TEM):** the
   TEM-pinned model band is BROADER than clean-monomer data at 29–42 nm and
   NARROWER at 8 nm, and the two ends calibrate the two damping terms
   separately:
   (a) 29–42 nm: swapping the Gaussian(raw TEM sd) for the true per-particle
   histogram changes the width by ~nothing (polydispersity is NOT the width
   driver at these sizes); scaling the BULK Drude damping down
   (γ_eff = s·γ_bulk + γ_S, s ≲ 0.2, degenerate below that) closes the
   +2.5 nm FW75 surplus and halves RMS 5→2–3%; adding the −2.7 nm offset
   (#9) reaches RMS 1.7–2.3% (from 3.7–5.1% baseline). Physical reading:
   J&C's evaporated-film ε₂ over-damps single-crystal-like colloids —
   equivalent to ε₂ ~15% (550 nm) to ~35% (650 nm) too high; do not read
   s literally as γ_bulk→2 meV (it proxies any ε₂ overestimate; Reddy's
   single-vs-poly factor ~6 is the anchor).
   (b) 8 nm CTAC: the opposite — data BROADER than the A=0.25 model by
   ~10 nm FW75; A_surf calibrates to ≈0.5–0.8 (best-RMS 0.5, best-width
   0.75) — chemical-interface damping above Berciaud's 0.25 for this
   chemistry/size. ⇒ γ_S calibration belongs at small D (where it
   dominates), bulk-ε₂ calibration at large D. JOINT CALIBRATION
   (`scripts/calibrate_damping.py`): ONE global pair (s=0.05, A_surf=0.6)
   fits all five CTAC sizes simultaneously (mean RMS 4.6→2.5%) — the 1/R law
   carries the size dependence; no per-size values needed.
   **OUT-OF-SAMPLE TRANSFER TESTED (frozen CTAC pair → seeded-growth citrate
   series, fig17): the −2.7 nm offset TRANSFERS (peak errors go from
   +1..+4 to −3..+1, centered on 0), the damping pair does NOT** — the
   calibrated (narrower) basis worsens the peak-region RMS on the cleanest
   seeded samples (GNP39nm 7.3→10.4%): seeded-growth H₂O₂/citrate particles
   are measurably MORE damped than CTAC ones, consistent with
   defect-richer/multiply-twinned overgrowth. ⇒ (s, A_surf) are
   PER-SYNTHESIS-ROUTE material parameters: calibrate once per route
   (ligand sets A, crystallinity sets s); the plasmon width doubles as a
   crystallinity probe. Production default remains uncalibrated
   (s=1, A per caller) until per-route adoption is decided.
   **SEEDED-GROWTH ROUTE CALIBRATED (2026-07-13,
   `scripts/calibrate_damping_seeded.py`, fig20; per HANDOFF.md): (s=1.75,
   A_surf=1.0) + the −2.7 nm offset**, jointly on the four verified-clean
   monomers (Oct GNP12 + Nov 11/23/39 nm; hand-picked per-particle TEM from
   ../tem-particle-metrics/outputs/dm3/gt/): clean-set mean RMS 7.74% (s=1,
   A=0.25 baseline) / 7.92% (CTAC pair) → 3.00%. OPPOSITE direction from
   CTAC: this route is MORE damped than jc bulk (s>1 proxies ε₂ ABOVE the
   evaporated film — defect-rich/multiply-twinned H₂O₂ overgrowth — where
   CTAC was single-crystal-like at s=0.05), confirming fig17's transfer
   failure and the per-route reading. The s/A split sits on a shallow
   diagonal ridge (γ_eff is the constrained combination: ≈0.23 eV at
   11.5 nm, ≈0.16 eV at 23 nm across variants; GNP12-only rails to s=3.0 —
   s is a flat direction at small D); no Oct-vs-Nov batch split (nov-only
   picks the same ridge, GNP12 near its optimum). Peak residuals after
   −2.7 nm: −1..+2 nm — no medium-index flag from the unpurified matrix.
   API: `dielectric.set_bulk_damping_scale(s)` (module-level like
   `use_gold_model`; default 1.0 = old behavior; applies inside
   `gold_epsilon_sized`/`size_damping_correction`; NB T-matrix caches bake
   damping at build time and do NOT record s — rebuild or use CDA).
   **SPECIATION VERDICT (`scripts/fit_speciation_seeded.py`, fig21;
   TEM-pinned NNLS, mirrors ../tem-particle-metrics fit_speciation.py):
   the Oct GNP12 13% "aggregation" was the line-width artifact (13→0%,
   likewise Nov11 6→0%), but Oct GNP20 does NOT collapse (41→45%; 41% at
   the GNP12-only ridge point) and GNP40/60 are stable (40/68%) — their
   residual misfit is red-asymmetric, a real aggregation signature, while
   the same recalibrated monomer nails every clean sample (Nov23 0%,
   Nov39 7%, RMS all improve; GNP55 stays 100% — negative controls pass).
   **EXACT T-MATRIX RE-FIT DONE (2026-07-13,
   `scripts/build_seeded_tmatrix_basis.py` + `fit_speciation_seeded_exact.py`,
   fig22): the Oct fractions HOLD OR RISE under exact optics.** Per-sample
   exact bases (dimer/trimer_linear/tetramer_linear — tetramer added to
   GEOMETRIES; calibrated damping baked in; caches
   outputs/seeded_tmatrix_basis_*.npz record s_bulk; probe-chosen lmax per
   gap): GNP20 **31%** agg gold (31–38% across gaps 1/2/3.5; RMS halves vs
   CDA 0.066→0.038), GNP40 **46%** (41–46%), GNP60 **81%** (75–81%), GNP55
   100% (RMS 0.108 — its smooth 580 plateau is BEYOND any single-gap
   dimer/trimer/tetramer basis: needs motif+gap distributions / bigger
   flocs; don't quote a fraction). ⚠ The intuited "CDA under-couples ⇒ its
   fractions are upper bounds" was FALSIFIED at large D: exact ≥ CDA for
   GNP40/60 — exact chain modes redistribute red intensity in shaped bands,
   not a uniform boost; quote the CDA↔exact spread as the backend
   systematic. GNP60's exact fit shows a ~620 nm chain-mode bump the smooth
   data lacks (single-motif discreteness, limitation #4). NUMERICAL LESSONS
   (in build_seeded_tmatrix_basis.py): (a) treams chain (trimer/tetramer)
   interaction solves go UNSTABLE (negative/spiky Cext) at lmax=12 for
   gap≤1 nm even where the dimer is fine — the dimer convergence proxy does
   NOT catch it; cap chains at lmax=10 near contact and check positivity;
   (b) gap-1 bonding-mode structure is REAL and sharp — build those columns
   on a 5 nm wl grid (per-gap `wl_gap<g>` axis in the cache).
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
   **(v) FIFTH EXCLUSION MEASURED — ε(T) (see #11) — AND THE PEDESTAL
   RE-MEASURED UNDER normalize="evaporation" + BLANKS (#14 fixed): the flat
   pedestal is REAL, and about half its T-growth was the Kell artifact.**
   Corrected numbers (exact γ_S basis, evaporation norm, blank-subtracted —
   i.e. with the static instrumental floor REMOVED from the data): the
   pedestal EXISTS (A_sca ≈ 0.13–0.19 of series max) and its T-growth is
   +9–16% across branches and ε(T) treatments — vs +31% under Kell:
   **roughly half the previously-fitted growth was evaporation-normalization
   error; the rest is a real, branch-consistent, reversible growth.** The
   concentration-drift diagnostic independently confirms a structured
   (non-concentration) branch component of ~1–1.8 pp in the red — the
   pedestal's own reversibility residue. **Its SHAPE, however, is
   ε(T)-model-dependent (see #11): n_sca = 0.00 under jc+Drude-retune vs
   0.85 under the full Reddy DCP — and RMS prefers the DCP (0.0138 vs
   0.0166). Treat the pedestal exponent as UNSETTLED in [0, ~0.9]; do not
   base the floc-vs-straylight call on it.** Floc-test (shape statement):
   0.7–1 µm Maxwell-Garnett–Mie spheres reproduce flat-to-weakly-decaying
   pedestals unimposed (equiv n = 0.00/−0.27 at fill 0.35; ~+1.3 at 200 nm);
   n≈0.85 sits between the large-floc and small-floc regimes (~300–500 nm
   effective sizes). Remaining discriminator (lab): DLS across the ramp /
   filtration — fitting cannot settle it; the reversible T-growth + flat
   blanks argue physical.
9. **[FIXED] Monomer peak ~4 nm red of experiment.** Root cause was the gold
   dielectric DATASET, not the medium: tabulated Johnson & Christy (model
   `'jc'`, embedded in `dielectric.py`, cubic-spline interpolated — linear
   interpolation pins the peak to the sparse J&C grid nodes) puts the 12.9 nm
   monomer at 522.8 nm in bare water (target 523.1). No ligand-shell index
   needed. NB: the Brendel-Bormann table (`'bb'`) is ~4–6 nm red-biased in the
   500–540 nm window just like Etchegoin — do not use either for calibration.
   `fit_real.py` sets `'jc'`; caches record their `gold_model` and `fit_real`
   refuses a mismatched cache.
   **RE-ANCHORED 2026-07-12 (`scripts/dielectric_reanchor.py`, fig18): jc
   carries a +2.7 ± 1.1 nm RED bias against verified-clean monomers.** Tested
   against NINE clean anchors (CTAC 7.8–42 nm + seeded-growth citrate
   Oct/Nov 2024; two chemistries): jc +2.7±1.1 / reddy_p200 +3.2±0.9 /
   etchegoin +7.4 / bb +7.9 / yakubovsky25 −4.3. jc REMAINS the best fixed
   table (do not switch), but the original C500 anchor match was likely
   flattered by C500's own aggregation red-shift (~15–30% aggregated gold) —
   treat model-peak expectations for clean monomers as jc − ~2.7 nm.
   Transferability (leave-one-out): the PEAK-SHIFT component of an empirical
   correction transfers (residual +1.0 nm mean); a FULL-SHAPE multiplicative
   c(λ) does NOT (LOO RMS worsens for 7/9 anchors; far-red ratios are
   ill-conditioned) — for shape-critical work use measured endmember spectra
   instead. Sizing consequence re-confirmed: shape-only D retrieval stays
   information-limited below ~35 nm with or without the correction and with
   poly pinned (a 22.9 nm sample fits at 7.8 nm); ≥38 nm is good to ~10%.
   Absolute extinction (#10, Next-step E) remains the required rescue.
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
    **[SUPERSEDED 2026-07-11] "Red-tail excess is reproducible in clean
    samples" — NOT reproduced under the calibrated model.** The earlier claim
    (data above Mie at ~560–660 nm in every CTAC sample) traced to the
    bb-era analysis. Re-measured with jc + TEM-pinned D AND poly
    (fig13_redtail_excess.png): the 29–42 nm samples sit 0.5–2.5% BELOW the
    model across 560–660 nm under both peak- and LSQ-normalization, robust to
    trimming the TEM low-side outliers (raw TEM poly 4.4–7.7% falls to
    2.8–4.9% trimmed — the TEM sd's are outlier-inflated, part of why the
    data look taller/narrower than the model). The 8 nm sample's +6.8% excess
    is inside the γ_S systematic (+3.1% at A_surf=0.25, −8.5% at A_surf=1).
    ⇒ The red-tail anomaly is NOT generic to gold colloids: it is (so far)
    specific to the citrate/TEG C500 sample — which strengthens the
    surface-chemistry-contrast hypothesis. The C500 far-red floor itself is
    damping-robust: measured A(790)/A(peak) = 0.146 blank-subtracted vs
    monomer model 0.0066 (no γ_S) / 0.0101 (A_surf=0.25) / 0.0228 (A_surf=1)
    — a ≥6× excess under even the most conservative damping, ~14–22× at the
    fit-preferred levels.
    **VALIDATION RUN (2026-07-11, `scripts/fit_ctac_validation.py`, mstm-env):
    exact-basis speciation fit on the CTAC series, D/poly PINNED to TEM,
    A_surf=0.25, gap=3.5 nm** (CTAC bilayer is ~3.2–4 nm, so smaller gaps are
    sterically impossible for this ligand; conveniently the T-matrix also
    converges there at lmax=8 for D≥29 nm, whereas gap 1 nm is unconverged
    even at lmax=10 AND unphysical for CTAC).
    - **Collinearity diagnostic** (cosine of per-gold-atom basis vectors,
      420–800 nm): CDA 12.9 nm/gap 1 = 0.9965 (reproduced); **EXACT 12.9/gap 1
      = 0.9593** — exact optics DO break the CDA collinearity at near-contact.
      BUT the gap dominates: EXACT 12.9/gap 3.5 = **0.9979** — collinear
      again. Per sample (gap 3.5): 8 nm **0.9997** (hopeless); 29/31/38/42 nm
      = 0.969/0.963/0.940/0.925 — identifiable, improving with size.
    - **Fits (negative control PASSED):** 29–42 nm all return 100% monomer,
      dimer gold <0.5% (profile bound at RMS+1%; <2.6% at RMS+5%) — with an
      identifiable basis the pipeline does NOT report spurious aggregates on clean
      monomer samples. The 8 nm sample nominally returns 96/0/4 (%m/d/t) but
      its profile bound is 26% dimer at RMS+1% — NOT identifiable; any
      speciation percentage at ~8 nm must carry that caveat.
    - **CONSEQUENCE FOR C500 (#8/#11): identifiability at ~13 nm is
      GAP-CONDITIONAL.** The aggregated-gold ≈0.15–0.30 numbers rest on the
      assumed 1 nm citrate/TEG gap (cos 0.96); if the true gap is ≳2–3 nm they
      become unidentifiable (cos→0.998). Constraining the gap (TEM of dried
      dimers / ligand-length argument) joins A_surf in Next-step E.
    - Residual caveat: the pinned-TEM fits at 29–42 nm misfit the PEAK (data
      taller & narrower than the model; RMS 0.035–0.047 in A units; pedestal
      takes zero weight) — same shape-systematics family as the red-tail
      excess above (jc interband, TEM-poly tails, CTAC medium index ≠ pure
      water). The %-fractions are robust to it (profile bounds above), but do
      not read the RMS as noise-limited.
    - **Size-mixture screen (fig15, `scripts/size_mixture_feasibility.py`,
      2026-07-11):** per-gold cosines on the same exact basis — m29/m42 =
      0.9976, yet NNLS retrieval of the 42 nm gold fraction in a synthetic
      29+42 mix is UNBIASED with 2σ detection ≈ 2% of gold at ≤0.3% iid noise
      (≈10% at 1%; structured peak-shape systematics are the practical floor —
      the mixing experiment measures it). m42/d29 = 0.9803, and a pure
      29-dimer fitted monomer-only leaves a 9.7%-of-peak structured residual
      (~30× noise) before grabbing the m42 column — dimer contamination is
      detectable, not silent. WARNING: d42 vs t29 = 0.9984 — cross-size
      AGGREGATE species are mutually degenerate; mixing experiments must keep
      aggregation suppressed for the two-monomer reading to stay clean.
      ADDENDUM (measured; CORRECTED 2026-07-11): shape-only BLIND analysis
      cannot detect the mixture as such — and the degeneracy is DEEPER than
      the first pass suggested. A single Gaussian population imitates the
      50/50 mix at RMS 0.01% of peak with poly free (D=32, 22% poly), and
      STILL at RMS 0.09% with poly constrained ≤10% (D=36) — i.e. the
      "implausibly-large-poly flag" is UNRELIABLE at this size gap. The
      reverse also holds: a narrow 35 nm/5% sample fits as a 34:66 "mix" at
      RMS 0.11%. All are below the measured per-point noise (0.13–0.21% of
      peak, 2nd-difference on the CTAC scans) and far below the 1–4%
      systematics floor. Absolute per-gold extinction rejects the BROAD
      imitation (+15% level error) but not the narrow one (1–2.4%). ⇒ Ratio
      retrieval REQUIRES the two-endmember prior (TEM'd stocks, as the
      mixing experiment provides); blind bimodality detection needs TEM/DLS,
      not UV-Vis. Retrieval is otherwise robust to polydispersity 3–20%
      when the widths are KNOWN (common-mode broadening; 2σ ±0.01).

11. **[IMPLEMENTED — AND THE VERDICT IS LARGELY AGAINST THE OLD READING] Gold ε(λ,T).**
    `gold_epsilon(..., temperature_C=)` now applies a bulk thermal Drude-damping
    retune on top of any base dataset: γ_bulk(T) from the **Holstein electron–phonon
    form (Θ_D = 170 K), anchored at Olmon's 44 meV @ 20 °C** — implemented from
    theory, stated explicitly. **⚠️ CORRECTION (2026-07-11): the claimed "~12%
    agreement" with Reddy is WRONG — recheck.** Reddy 2016 (author list verified:
    Reddy, Guler, Kildishev, Boltasseva, Shalaev, *OME* 6, 2776; Drude+2CP tables
    from the open-access full text) MEASURED Γ_D:
      200 nm SINGLE-crystal (Table 6): 0.0534 eV (23 °C) → 0.0725 (200 °C)
        ⇒ dΓ/dT = **1.08e-4 eV/K**  (≈ +12.3% over our 15→75 °C ramp)
      200 nm POLY-crystalline (Table 3): 0.0471 eV (23 °C) → 0.0489 (100 °C)
        ⇒ dΓ/dT = **0.023e-4 eV/K** (≈ +3.0% over our ramp)
    Holstein anchored at Olmon 44 meV @20 °C rises 20.1% over 60 K
        ⇒ dΓ/dT ≈ **1.47e-4 eV/K** — i.e. **36% steeper than Reddy single-crystal
    and ~6× steeper than Reddy poly**, not 12%.
    Reconciliation: *measured* Γ_D includes T-INDEPENDENT grain-boundary/defect
    scattering, which damps the RELATIVE rise; Holstein is the pure γ_ep term.
    **Citrate-grown AuNPs are POLYcrystalline**, so Table 3 is arguably the better
    analogue. **BRACKETED (2026-07-11): the Drude-slope systematic turns out NOT
    to matter, but the INTERBAND treatment does.** Fits with γ_ep(T) scaled to
    Holstein / Reddy-single / Reddy-poly (`use_eps_t_scaling`; CDA+γ_S basis,
    evaporation norm, heating) give near-IDENTICAL conclusions (D 8.16–8.24,
    ΔH₂ −10..−13, n_sca 0.00, A_sca +19%, RMS 0.0179 for all three — γ_S
    dominates γ_total, so even the 6× slope spread moves the basis by too
    little to matter). BUT the FULL Reddy DCP ε(λ,T) (`'reddy_p200'`, Table S1
    interpolated in T — interband parameters move too: ε∞ 2.27→2.45, γ₂
    0.256→0.273 already by 100 °C) is a genuine outlier: peak response −2.56%
    vs Drude-only −0.33% (the interband DOMINATES), its 12.9 nm peak lands at
    522.8 nm (same as jc — usable), and in the CDA fit it returns n_sca 0.70,
    halves the aggregated fraction and improves RMS (0.0142). **EXACT-basis
    refits under 'reddy_p200' (cache outputs/tmatrix_basis_reddy.npz,
    evaporation norm + blanks) confirm the interband sensitivity: heating
    D=10.66, ΔH₂=+0.5±0.9, n_sca=0.86±0.02, agg 0.289→0.295, RMS 0.0138;
    cooling n_sca=0.84, RMS 0.0137 (ΔH₂ rails/unidentified on 6 points).
    RMS prefers the full DCP on both branches (0.0138 vs 0.0166, −17%).**
    ⇒ ROBUST across ε(T) treatments: no thermal speciation signal; pedestal
    exists with reversible ~+10–16% growth; D rails low. ε(T)-MODEL-DEPENDENT
    (the dominant remaining systematic): the pedestal SHAPE (n_sca 0.00 under
    jc+Drude-retune vs 0.85 under full DCP) and the aggregated-fraction LEVEL
    (0.15–0.30). γ_S stays T-independent per
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
    back-of-envelope "~40%" estimate was directionally right; the exact share
    depends on the γ_S dilution it itself identified. **γ_S must be in the basis for ε(T) fits**
    (`size_correction=True`, now wired through the whole basis incl. the
    polydispersity integral — closes #2's open item) or the fit over-attributes
    T-changes to gold ~3×. Water n(T) adds peak −2.85%/blue −0.4 nm by itself.
    **REFIT VERDICT — now under normalize="evaporation" + blanks (#14 fixed;
    the Kell-era numbers are superseded):** the speciation conclusion SURVIVES
    the corrected normalization. With ε(T)+n(T)+γ_S and the evaporation-aware
    concentration correction: ΔH₂ = −15.2±16.3 (heating; consistent with 0) /
    +2.8±2.9 (cooling); aggregated gold ≈ 0.15–0.21, nearly T-flat (heating
    0.181→0.150, cooling 0.186→0.208 — opposite weak trends, no robust
    thermal direction); **branch consistency improves dramatically vs Kell
    (0.18 vs 0.19 at 15 °C; Kell said 0.24 vs 0.10)** — a strong sign the
    normalization is now right. **There is NO robust thermally-driven
    monomer⇌dimer⇌trimer signal in the C500 series.** What the temperature
    does move: gold+water ε(T) physics (the peak region) and the pedestal's
    T-growth (#8). D falls to 10.3–10.6 with γ_S(A_surf=1) — the A_surf/D/poly
    degeneracy stands (pin size via absolute extinction, #10; calibrate
    A_surf). RMS 0.0166/0.0172 — identical to Kell: fit quality CANNOT choose
    the normalization; the branch-offset physics does. **Caveats retained:
    (a) ε(T)-ON and ε(T)-OFF fits tie on RMS (0.0166 vs 0.0164) — the data
    cannot statistically distinguish "gold heats" from "particles aggregate";
    ε(T) is in the model because bulk gold measurably does this, an a-priori
    argument. (b) The γ(T) magnitude itself is a factor-6 systematic
    (Holstein vs Reddy-poly, above) — bracketed below.**
12. **[IMPLEMENTED] Normalization re-decided: `normalize="density"` (Kell 1975).**
    `load_series` now supports a pure dilution correction
    A·ρ(T_ref)/ρ(T) (Kell polynomial, verified ρ(15)=999.10, ρ(75)=974.85 —
    NB the denominator coefficient is 16.87985e-3; a digit transposition was
    caught against these reference values). No anchor-wavelength assumption;
    needs parseable temps in every column header. **`mult_400nm` is DEPRECATED/
    BIASED** (premise measurably violated: A(400) rises +3.18% while Kell
    predicts −2.43%; flat rescale inflates the apparent peak change ~8.5×);
    kept only for reproducing old results. `fit_real.py --normalize density`.
    **ISOSBESTIC — RESOLVED under normalize="evaporation" (+blanks): a
    crossing EXISTS but it DRIFTS, so it is NOT a two-state isosbestic.**
    Measured (data-level, ΔT-vs-15 °C curves, blank-subtracted): heating
    crossing walks ~544 nm (30–50 °C) → 550.5 (55) → 556 (65) → **568 nm
    (75 °C)**; cooling sits ~546–554 nm, essentially flat. (Without blank
    subtraction the same drift appears ~10 nm redder — 562→577 nm — matching
    the hand-analysis; the blank offset shifts the crossing, another reason
    blanks matter.) The speciation residual (pedestal removed) also has
    bipolar structure with drifting crossings (heating ~618→634 nm). VERDICTS:
    (a) the earlier Kell-based "no crossing anywhere" was a normalization
    artifact — under Kell the difference curves are positive everywhere and no
    crossing is possible; (b) a TRUE isosbestic is stationary — a ~24 nm
    drifting crossing instead indicates an EVOLVING cluster geometry
    (cf. Cardellini 2025, ref 11: λ_iso tracks interparticle spacing), or a
    residual mixture of >2 effective species; (c) the historical single
    "575 nm isosbestic" value is not recoverable — it was read off
    mult_400nm-normalized curves. The two-state argument cannot rest on a
    stationary isosbestic; the drift itself is the finding.
    **MODEL CROSS-CHECK (demo.py fig3, regenerated with the production stack —
    exact cached T-matrix, jc + γ_S, endothermic K(T), D=12.9):** a
    fixed-geometry monomer⇌dimer⇌trimer series gives a STATIONARY crossing at
    ~532 nm (spread 0.25 nm), RED of the 524 nm peak; adding full gold+water
    ε(T) shifts it to ~539 nm but it STAYS stationary (drift <1 nm over
    15→75 °C). Two consequences: (a) the legacy demo's blue-side crossing
    (~524 nm, CDA/etchegoin) was an artifact of CDA under-coupling (weak shift
    + amplitude enhancement pulls the crossing up the blue flank) — fixed;
    (b) the measured ~24 nm drift CANNOT be produced by a fixed-geometry
    equilibrium even with ε(T) in the model — independent support for the
    evolving-geometry reading above.
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
    normalization error.
    **FIX IMPLEMENTED — `normalize="evaporation"` (SALVAGE mode) + the
    reusable diagnostic.** Three pieces, in order of importance:
    (1) `io_data.check_concentration_drift(heating, cooling)` — run on EVERY
    dual-branch dataset; decomposes the matched-T branch ratio into a
    wavelength-FLAT part (concentration drift; WARNS that "density" is blind
    to it) and a wavelength-STRUCTURED part (irreversible chemistry/pedestal
    residue) — the two COEXIST on C500 (flat +6.1% @15 °C shrinking to +1.9%
    @65 °C, plus ~±1–1.8 pp structure). Would have caught this in 2011.
    (2) `normalize="density"` stays the RECOMMENDED mode when concentration is
    VERIFIED stable (weigh the cell before/after; minimise headspace). NB the
    2011 C500 cell WAS sealed and still lost ~5.6% — a seal alone proves
    nothing; always run the drift diagnostic on dual-branch data.
    (3) `normalize="evaporation"`: c(T_i)/c(T_ref) = [ρ(T_i)/ρ(T_ref)]·
    [1+α·E_i], E_i = cumulative Antoine-psat scan exposure (simple ramp
    reconstruction by default; optional true scan_order.csv — convention only
    rescales α, c/c0 invariant to <0.3%); α is ONE parameter FITTED model-free
    from the blue-band (400–470 nm, Haiss-450-centred) matched-T offsets —
    never hardcoded; needs BOTH branches (`companion=`, `branch=`), else falls
    back to density with a LOUD warning. Calibrated on C500: α = 8.6e-6
    (blank-subtracted; raw 8.4e-6 vs hand-analysis 8.0e-6), scatter 6%,
    residuals ±0.35%; per-region robustness confirms the blue is purest
    (blue 8.6 / peak 7.4 / red 7.8e-6 with 31% scatter — the blue/red
    divergence is the pedestal's leakage measure). c/c0 reproduces the
    hand-computed table to ≤0.65% (heating 75 °C: 1.023 vs 1.015; cooling
    15 °C: 1.058 vs 1.054) — decisively away from Kell's 0.976: the
    degeneracy IS broken. BLANKS: `scripts/export_blanks.py` pulls the
    per-temperature water blanks from final.xls (non-zero: 0.031/0.036 at
    700/790 nm; T-drift only −0.003) into companion CSVs;
    `load_series(blanks=, companion_blanks=)` subtracts them BEFORE
    normalization. Re-run verdicts: see #8/#11/#12.

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
  Fig 3 uses the PRODUCTION stack (jc + γ_S + ε(T); exact cached optics if
  outputs/tmatrix_basis.npz exists, else CDA fallback — the crossing then sits
  less red). Figs 1/2/4 stay in the legacy etchegoin/CDA prototype config.
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
- `mstm-env/bin/python scripts/fit_ctac_validation.py` — fig12; CTAC size-series
  validation (D/poly pinned to TEM, A_surf=0.25, gap 3.5 nm): collinearity
  diagnostic + pinned exact-basis speciation fits (see #10). ~20 min first run;
  exact basis cached in experimental/ctac/.
- `python scripts/size_mixture_feasibility.py` — fig15; 29+42 nm size-mixture
  identifiability screen (per-gold cosine matrix + NNLS retrieval vs noise;
  system env, reuses the ctac cache). Pre-registration for the mixing
  experiment; results recorded in #10.
- `python scripts/calibrate_damping_seeded.py` — fig20; seeded-growth-route
  (s, A_surf) grid calibration on the four verified-clean monomers (needs
  ../tem-particle-metrics for the per-particle TEM gt). ~3 min; see #2.
- `python scripts/fit_speciation_seeded.py` — fig21; the HANDOFF.md
  deliverable: TEM-pinned monomer/dimer/trimer NNLS on all 8 seeded samples,
  first-cut (s=1, A=0.25) vs recalibrated (1.75, 1.0, −2.7 nm). Verdict in #2.
- `mstm-env/bin/python scripts/build_seeded_tmatrix_basis.py --probe|--sample
  SID` — exact per-sample cluster bases for the 4 aggregated seeded samples
  (calibrated damping baked in; --probe first: lmax convergence + the gap-1
  chain-instability guard). ~30–75 min per sample; shard and run in parallel.
- `python scripts/fit_speciation_seeded_exact.py` — fig22; CDA-vs-exact
  speciation on those caches, with per-gap sensitivity. Results in #2.
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
Two failure modes found during development:
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

**A0–D + F: DONE, verdicts now MEASURED under normalize="evaporation" +
blanks** (see #8/#11/#12/#14). Where the science stands — robust across ε(T)
treatments: no thermally-driven monomer⇌dimer⇌trimer signal (ΔH₂≈0 or
unidentified; agg ~T-flat); the pedestal is REAL with half its former
T-growth an artifact (+9–16% survives, reversible); the isosbestic crossing
exists but DRIFTS 544→568 nm (not two-state; drift = evolving geometry,
cf. Cardellini ref 11). Remaining model systematics, in order: (1) the
INTERBAND ε(T) treatment — Drude-retune vs full Reddy DCP changes the
pedestal exponent (0.00 vs 0.85) and the agg level (0.15–0.30), and RMS
prefers the DCP (γ-slope choice within the Drude retune is immaterial —
bracketed, #11); (2) the A_surf–D–poly degeneracy (D rails to ~10 vs TEM
12.9). Decide the interband treatment (or report both), then E/G below.

**E. Break the A_surf–D–poly degeneracy (#2, #10, #11):** pin size with ABSOLUTE
   extinction + Haiss A_spr/A_450 (not the normalized shape), and calibrate
   A_surf on a TEM-pinned monomer sample. With γ_S(A_surf=1) in the basis the
   fitted D rails to 10.5–10.8 vs TEM 12.9 — the biggest remaining model knob.

**G. EXPERIMENTAL (lab, not code) — now THE decisive test:** DLS across the ramp
   and/or a filtered/centrifuged control. The fit cannot distinguish µm-flocs
   from stray light (both are λ-flat); only the T-reversible growth + flat
   blanks argue for a physical scatterer. No amount of fitting settles it.

**H. Model-selection caveat:** ε(T)-ON and ε(T)-OFF fits tie on RMS (0.0166 vs
   0.0164). To the question "is it the gold or the particles?", the answer is
   "bulk gold measurably does this (Reddy), at exactly the magnitude used" —
   an a-priori argument, not a fit-quality one. Keep it that way; do not tune ε(T) parameters against the C500 data.

I. Optional: NN surrogate trained on Layer-1+2 spectra for instant inference.
