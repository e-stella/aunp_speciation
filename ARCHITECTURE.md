# A guided tour of the code

This walks through every module in the order the physics flows — from "what
colour is gold" up to "fit a real spectrum." For each file you get: the idea, the
key functions, and the papers behind it. If you read top to bottom you will
understand the whole tool.

The big picture in one sentence: **we model a gold colloid as a temperature-
dependent mixture of monomers, dimers and trimers, compute each species'
UV-Vis spectrum, add them up, and fit that model to experimental spectra to
recover size, polydispersity and the amount of clustering.**

The code is organised in three conceptual layers:

- **Layer 1 — populations**: how much monomer / dimer / trimer is present, and
  how that shifts with temperature. → `equilibrium.py`
- **Layer 2 — optics**: the extinction spectrum of a single sphere and of small
  clusters. → `dielectric.py`, `mie.py`, `clusters.py`, `clusters_tmatrix.py`,
  `spectra.py`
- **Layer 3 — inference**: given a measured spectrum (or a temperature series),
  fit the model parameters. → `fitting.py`, `fit_global.py`
- plus **plumbing**: `basis_cache.py` (make the exact optics fast enough to fit
  with) and `io_data.py` (load real files).

---

## Layer 2a — the colour of gold: `dielectric.py`

**Idea.** Everything optical starts from the metal's complex permittivity
ε(λ) = ε₁ + iε₂ (equivalently the refractive index n + ik). ε₁ controls where the
plasmon resonance sits; ε₂ (the imaginary/absorptive part) controls how broad and
damped it is.

**Available models** (`use_gold_model(name)`; the cache records which one was used and
`fit_real.py` refuses a mismatch):
- `etchegoin` — analytic Drude + 2 interband critical points. Fast, no data files.
  **Peak ~527 nm: ~4–6 nm too red. Do not use for calibration.**
- `bb` — tabulated Brendel–Bormann (Rakić 1998). **Same ~4–6 nm red bias as etchegoin.**
- `jc` — tabulated Johnson & Christy (1972), 49 points, cubic spline. **PRODUCTION for
  the 13 nm system**: puts the 12.9 nm monomer at 522.8 nm (target 523.1).
- `yakubovsky25` / `yakubovsky53` — tabulated thin-film n,k (Yakubovsky 2017), size-matched
  (25 nm film for R < 25 nm, 53 nm film above), per lit-map ref 18. Bare `yakubovsky`
  deliberately raises — the film choice must be explicit. Peak sits ~7.5 nm BLUE of
  experiment at 12.9 nm, which the global fit converts into fake aggregation, so it is
  the **systematics bound and the ≥50 nm option**, not production at 13 nm.

**Damping decomposes** as γ = γ^(ep)(T) + γ_S, where γ^(ep) is the bulk electron–phonon
term and γ_S = 3ℏk_F/(4Rmc) is the size-dependent surface/Landau damping. For a 13 nm
particle γ_S ≈ 0.143 eV DOMINATES γ_bulk ≈ 0.044 eV. A small-particle (size-damping)
correction exists as `gold_epsilon_sized` / `A_surf`.

**ε(λ, T) IS IMPLEMENTED.** `gold_epsilon(..., temperature_C=)` applies a bulk
thermal Drude-damping retune on top of any base dataset: γ_bulk(T) from the
Holstein electron–phonon form (Debye Θ_D = 170 K), anchored at Olmon's 44 meV
@ 20 °C — implemented from theory (stated explicitly in the docstring), and
validated against Reddy 2016's measured single-crystal Drude slope (~12%
agreement near room T). γ_bulk rises 20% over 15→75 °C. γ_S stays
T-independent (Chetoui) and composes with γ_bulk(T) in a SINGLE retune
(`gold_epsilon_sized(λ, D, A_surf, temperature_C)`) — the two deltas are not
additive. Temperature threads through mie/clusters/clusters_tmatrix/spectra
and the fitters; the T-matrix cache carries a temperature axis and records
`eps_t_model`. Measured impact at D=12.9 (15→75 °C): with γ_S in the basis,
peak −1.3%, @700 +3.7% — about a fifth to half of the observed changes; without
γ_S the basis over-responds ~3× (γ_S dilutes the thermal term), which is why
`size_correction=True` is REQUIRED for quantitative ε(T) fits. Verdict on the
C500 series: with ε(T)+n(T) in the model the thermally-driven speciation
signal collapses (ΔH₂≈0) — see CLAUDE.md #11.

**Key functions.** `gold_epsilon(λ, model=, temperature_C=)`,
`size_damping_correction` and `gold_epsilon_sized(λ, D, A_surf, temperature_C)`,
`gold_index(...)` (returns n + ik), `gamma_bulk_ev(T)`,
`thermal_damping_correction(λ, T)`; `medium_index("water", temperature_C)` →
COMPLEX (CRC n(T): 1.3334@15 → 1.3229@80 °C — heating blue-shifts, masking part
of any red-shift signal); `water_k(λ)` (tabulated 400–900 nm water absorption,
incl. the 740–760 nm O–H bump; diagnostics only — not wired into Mie, where an
absorbing host is ill-defined and blank-referencing cancels the bulk path).

**Why the small-particle correction matters.** In a 12 nm particle the conduction
electrons hit the surface before completing a mean free path, which broadens the
resonance. At 12 nm this alone widens the monomer peak from ~79 to ~117 nm FWHM —
so it is a real confound you must include before blaming broadening on clustering.

**Papers.**
- Analytic gold model: Etchegoin, Le Ru & Meyer, *J. Chem. Phys.* 124, 164705 (2006).
- Measured optical constants (production choice): Johnson & Christy, *Phys. Rev. B*
  6, 4370 (1972); Olmon et al., *Phys. Rev. B* 86, 235147 (2012) — also the
  Drude baseline used here (ħω_p = 8.45 eV, γ_bulk(RT) = 44 meV); Yakubovsky et
  al., *Opt. Express* 25, 25574 (2017).
- Small-particle (surface-scattering) damping: Kreibig & Vollmer, *Optical
  Properties of Metal Clusters* (Springer, 1995); Coronado & Schatz, *J. Chem.
  Phys.* 119, 3926 (2003).
- **Temperature dependence of ε:** Reddy, Guler, Kildishev, Boltasseva &
  Shalaev, *Opt. Mater. Express* 6, 2776 (2016) — measured Drude+2CP ε(λ,T),
  23–500 °C (ε₂ nearly doubles by 500 °C; ε₁ marginal); Chetoui et al. 2026 —
  the γ(T) = γ^(ep)(T) + γ_S decomposition, γ_S explicitly T-independent
  (lit-map §G).
- **Which ε to actually use** (very relevant here): Klinavičius et al., *J. Phys.
  Chem. C* 129, 17616 (2025) — ref 18; they show the imaginary part dominates
  size-fit accuracy and thin-film Yakubovsky ε generalises best for colloids
  (measured here: at ~13 nm, peak-anchored fitting still favors `jc` — see
  CLAUDE.md #2).

---

## Layer 2b — one sphere: `mie.py`

**Idea.** A single sphere has an exact analytic solution (Mie theory). We compute
the Mie coefficients aₙ, bₙ and from them the extinction/scattering/absorption
cross sections. We also extract the **dipole polarizability** α from the first
coefficient a₁ — this is the bridge to the cluster model, where each sphere is
treated as one polarizable dipole.

**Key functions.** `mie_ab(m, x)` (the coefficients; m = relative index, x = size
parameter), `monomer_cross_sections(diameter, λ)` (the single-particle spectrum;
set `size_correction=True` to use the small-particle ε), `dipole_polarizability(...)`
(α from a₁, used by `clusters.py`).

**Sanity checks** live in `scripts/verify.py`: the small-size limit of a₁ matches
the Rayleigh formula exactly, and the 12 nm peak lands near 520–527 nm.

**Papers.**
- Mie theory / the aₙ, bₙ recursion: Bohren & Huffman, *Absorption and Scattering
  of Light by Small Particles* (Wiley, 1983) — the standard reference.
- Size (and concentration) from a single-sphere UV-Vis spectrum, the classical
  baseline our tool extends: Haiss, Thanh et al., *Anal. Chem.* 79, 4215 (2007);
  Amendola & Meneghetti, *J. Phys. Chem. C* 113, 4277 (2009).

---

## Layer 2c — small clusters, the fast way: `clusters.py`

**Idea (coupled-dipole approximation, CDA).** When two particles are close their
plasmons couple: the mode polarised *along* the pair red-shifts and brightens.
Because a 12 nm particle is far smaller than the wavelength, we can approximate
each sphere as a single dipole (its Mie α) and solve the small linear system of
dipoles interacting through the retarded dipole field. Cheap and analytic.

**Key functions.** Geometry builders `dimer`, `trimer_linear`, `trimer_triangular`
(positions given a surface-to-surface `gap`); `_green_tensor` (the dipole–dipole
interaction); `cluster_cross_section(...)` / `species_spectrum(species, ...)`
(orientation-averaged extinction per cluster).

**What it gets right and wrong.** It captures the qualitative red-shift/red-tail,
and it is exact for a single particle (regression-checked against `mie.py`). But
one dipole per sphere **under-estimates coupling at near-contact** (it omits
higher multipoles) — so its red-tail is a *lower bound*. That is why the exact
backend below exists.

**Papers.**
- CDA / discrete-dipole formalism: Draine, *Astrophys. J.* 333, 848 (1988);
  Draine & Flatau, *J. Opt. Soc. Am. A* 11, 1491 (1994) (the DDSCAT code).
- The distance-decay of plasmon coupling (the physics of the red-shift; the
  "plasmon ruler"): Jain, Huang & El-Sayed, *Nano Lett.* 7, 2080 (2007) — in the
  project library.
- Coupling/red-tail in real Au aggregates, and why point-dipole under-counts it:
  Kelesidis et al., *Anal. Chem.* 94, 5310 (2022) — in the project library.

---

## Layer 2d — small clusters, the exact way: `clusters_tmatrix.py`

**Idea (multi-sphere T-matrix, a.k.a. MSTM / GMM).** For clusters of spheres there
is an *exact* semi-analytic method that includes all multipole orders, not just
the dipole. We use the `treams` library. This is the quantitative engine: at a
1 nm gap it produces 1.6× (dimer) / 2.8× (trimer) more red-tail than CDA, and a
genuinely resolved longitudinal peak.

**Key functions.** `species_spectrum_tmatrix(species, ...)` — same signature as
the CDA version, so it is a drop-in. `available()` tells you whether `treams`
imported.

**Practical note.** `treams` needs older numpy/scipy, so it lives in a separate
`mstm-env` virtualenv (see README). It is also slower — hence `basis_cache.py`,
which precomputes it. Validation: it reproduces Mie for a monomer to 1e-15.

**Papers.**
- Multi-sphere T-matrix method: Mackowski & Mishchenko, *J. Opt. Soc. Am. A* 13,
  2266 (1996); generalized multiparticle Mie: Xu, *Appl. Opt.* 34, 4573 (1995).
- The `treams` implementation: Beutel, Fernandez-Corbaton & Rockstuhl, *Comput.
  Phys. Commun.* 297, 109076 (2024).
- Why exact multipoles matter for touching Au spheres: Kelesidis et al., *Anal.
  Chem.* 94, 5310 (2022) — they needed 10⁵ dipoles/particle in DDA to converge.

---

## Layer 1 — how much of each species: `equilibrium.py`

**Idea.** Treat clustering as a reversible chemical equilibrium
monomer ⇌ dimer ⇌ trimer, with association constants K₂, K₃. Temperature enters
through the van 't Hoff relation K(T) = exp(−(ΔH − TΔS)/RT). Given the total gold
and the K's, we solve mass balance for how many monomers/dimers/trimers exist.
When the spectrum is a population-weighted sum of species that share crossing
wavelengths, changing T pivots the family of spectra through a fixed
**isosbestic point** — the thermodynamic fingerprint the lab observed.

**Key functions.** `association_constants(T, dH2, dS2, dH3, dS3)` (van 't Hoff),
`solve_populations(C_tot, K2, K3)` (mass balance), `gold_fractions(...)` (fraction
of gold as monomer/dimer/trimer).

**Papers.**
- Thermally reversible, temperature-controlled Au aggregation with the
  corresponding UV-Vis and isosbestic behaviour: e.g. *Nano Lett.* 13, 5844 (2013)
  (H-bonded thermoreversible AuNP assembly); *Langmuir* 27, 14761 (2011)
  (thermosensitive polymer-coated AuNP).
- Statistical-mechanics of reversible self-assembly (the modelling analogue —
  distributions of assembly states in equilibrium): Hagan, *Phys. Rev. E* /
  related 2008; Elrad & Hagan, *Nano Lett.* 8, 3850 (2008) — both in the project
  library.
- Aggregation kinetics / cluster-size populations (the kinetic counterpart):
  Kim et al., *J. Colloid Interface Sci.* 318, 238 (2008); Jungblut & Eychmüller,
  "Modeling nanoparticle aggregation" (2019) — in the project library.

---

## Layer 2e — assembling the ensemble spectrum: `spectra.py`

**Idea.** Glue Layers 1 and 2 together. Take the per-species spectra (from CDA or
the exact backend), average the monomer over a Gaussian **size distribution**
(polydispersity), and add the species weighted by their populations. This is the
"forward model": parameters in → predicted spectrum out.

**Key functions.** `gaussian_sizes(mean, pct)` (the size distribution),
`monomer_polydisperse(...)`, `species_basis(..., backend=...)` (the per-species
spectra; `backend` is `'cda'`, `'tmatrix'`, or a callable such as the cached
interpolator), `mix(basis, weights)` (population-weighted sum).

**A key finding baked into the design.** Polydispersity has only a *weak* imprint
on the spectrum — so broadening beyond a few percent size spread is evidence of
clustering, not size. That is exactly the hypothesis the tool tests.

**Papers.**
- Polydispersity barely broadens a monomer spectrum, and coupling is the missing
  ingredient: Bilén et al., *J. Phys. Chem. C* 128, 13909 (2024) — ref 19;
  Klinavičius et al., *J. Phys. Chem. C* 128, 9662 (2024) — ref 17.

---

## Layer 3a — fit one spectrum: `fitting.py`

**Idea.** Invert the forward model for a single measured spectrum. We use
*separable* (variable-projection) least squares: the nonlinear parameters (mean
diameter, % polydispersity) set the basis shapes, and the species mixing weights —
which are linear — are solved at each step by non-negative least squares (NNLS).
An outer optimiser handles the two nonlinear parameters. We also report
**identifiability**: which parameters the data actually constrain.

**Key functions.** `fit_spectrum(wl, ext, ...)` → a `FitResult` with diameter,
polydispersity, per-species and **aggregated** gold fractions, uncertainties, and
verdicts.

**The honest result.** From one spectrum you robustly get the *aggregated* gold
fraction; size, polydispersity, and the dimer-vs-trimer split are under-determined
(they trade off). That motivates the global fit next.

**Papers.**
- The inverse "spectrum → size distribution" problem and the "replace TEM" goal:
  Haiss et al., *Anal. Chem.* 79, 4215 (2007); Klinavičius et al. 2024 (ref 17);
  Bilén et al. 2024 (ref 19); Glaubitz et al., *J. Phys. Chem. C* 128, 421 (2024)
  — ref 20.

---

## Layer 3b — fit a whole temperature series: `fit_global.py`

**Idea.** The degeneracies a single spectrum leaves open are broken by fitting the
*whole temperature series at once*. All spectra share one sample, so they share
diameter, polydispersity and gap; the total gold is conserved (one global
amplitude); and the populations move with T by van 't Hoff. The isosbestic
structure is now a constraint, not just a curiosity. This recovers size,
thermodynamics (ΔH, ΔS), and the aggregation-vs-temperature curve together.

**Key functions.** `fit_temperature_series(temps_K, wl, spectra, ...)` →
`GlobalFitResult` with the shared parameters, the fitted ΔH/ΔS, and per-temperature
gold fractions.

**Papers.**
- Van 't Hoff / global ("global analysis") fitting of spectral titrations is
  standard chemometrics; the isosbestic point as a two-state signature is textbook
  physical chemistry. The physical system (temperature-driven reversible AuNP
  equilibria with an isosbestic point) is the reversible-aggregation literature
  cited under `equilibrium.py`.

---

## Plumbing 1 — make the exact optics fast: `basis_cache.py`

**Idea.** The exact T-matrix is too slow to call inside a fit loop, so precompute
it *once* on a grid of (diameter, gap) — in the `mstm-env` venv — save to `.npz`,
and interpolate at fit time from the ordinary interpreter. The interpolator has
the same signature as `species_spectrum`, so you pass it straight into
`species_basis(..., backend=cache.species_fn)`. Result: fits run on exact optics
with no `treams` cost per iteration.

**Key functions.** `build_grid(...)` / `save_cache(...)` (build step, venv),
`load_cache(path)` → `CachedBasis` with `.species_fn` (fit-time interpolator).

**Gotcha (documented from experience):** the optics backend used to *fit* must
match any used to *simulate* test data, or size and aggregation come out biased.

---

## Plumbing 2 — read real files: `io_data.py`

**Idea.** Load experimental UV-Vis from CSV/TSV. It handles a single spectrum
(wavelength, extinction) or a wide **temperature series** (wavelength column +
one extinction column per temperature, with the temperature parsed from the
header, e.g. `ext_5C`). Baseline/units are left to you; the fitters work on
normalised shapes.

**Key functions.** `load_series(path)` → (wavelength, spectra, temps_K or None),
`save_series(...)` (used to write the bundled example).

---

## The `scripts/` folder — runnable entry points

Each script is a small driver that produces one figure or one check:

- `verify.py` — physics sanity checks (run after touching `mie.py` /
  `clusters.py` / `dielectric.py`).
- `demo.py` — figures 1–4: species spectra, speciation broadening, the isosbestic
  point, the gap lever.
- `fit_demo.py` — figure 5: single-spectrum fit self-test.
- `validate_mstm.py` — figure 6: exact-vs-CDA (run with `mstm-env/bin/python`).
- `fit_global_demo.py` — figure 7: global temperature-fit self-test.
- `build_tmatrix_basis.py` — precompute the exact-optics cache (venv).
- `make_example_data.py` — write `data/example_series.csv`.
- `fit_real.py` — load a real UV-Vis file and fit it (figure 8). **Start here to
  fit your own data:** `python scripts/fit_real.py your_file.csv`.

---

## Suggested reading order for a newcomer

1. `dielectric.py` → `mie.py` (one sphere), run `scripts/verify.py`.
2. `clusters.py` then `clusters_tmatrix.py` (two/three spheres), run
   `scripts/validate_mstm.py`.
3. `equilibrium.py` (populations), then `spectra.py` (put it together), run
   `scripts/demo.py`.
4. `fitting.py` then `fit_global.py` (inference), run `scripts/fit_demo.py` and
   `scripts/fit_global_demo.py`.
5. `basis_cache.py` + `io_data.py` + `scripts/fit_real.py` (fit real data).

`CLAUDE.md` has the developer-oriented version (limitations, next steps, the
two-interpreter setup); this file is the conceptual tour.
