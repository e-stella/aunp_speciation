# Stochastic Multimer Formation & UV–Vis Broadening in Gold Nanoparticle Colloids
### A Literature Map and Simulation-Tool Roadmap

**Context:** 12 nm AuNP, 4% polydispersity. The measured UV–Vis is broader / red-tailed than an FDTD extinction that accounts only for the 4% size spread. Working hypothesis: reversible monomer ⇌ dimer ⇌ trimer ⇌ multimer equilibria in liquid, so the UV–Vis is an ensemble / time-integrated snapshot of coexisting states. Temperature-dependent UV–Vis shows an **isosbestic point**, consistent with two (or few) interconverting optical species. Goal: Monte-Carlo the cluster population → compute a simulated spectrum → compare to experiment → invert into a tool that reports nominal size + size distribution (%), as an alternative to TEM.

---

## Why the hypothesis is physically sound (the short version)

A single 12 nm Au sphere has a narrow, symmetric dipolar LSPR (~520 nm) with negligible radiative broadening at this size. When two particles come within ~1 diameter, near-field coupling splits the mode: the longitudinal (bonding) mode red-shifts strongly and gains oscillator strength, while the transverse mode blue-shifts weakly. Because dimers/trimers form at **random orientations and gaps**, an ensemble of them contributes a smear of red-shifted intensity — exactly a red tail + apparent broadening on top of the monomer peak. A pure size-distribution model (FDTD with only the 4% spread) cannot reproduce this because it never introduces coupling. The isosbestic point is the thermodynamic fingerprint: it appears when population shifts between two species whose spectra cross at a fixed wavelength — the classic signature of a two-state equilibrium.

**Three independent ML papers corroborate the effect as everyone else's "error term."** Klinavičius (ref 17) had to *sonicate agglomerates away* to get good monomer fits and explicitly warns that aggregation red-shifts and broadens the spectrum; Bilén (ref 19) built polydisperse spectra by additive mixing "not considering… plasmon coupling" and attributes the sim-vs-experiment gap to exactly that; and Bilén further finds that **polydispersity barely imprints on a monomer-only spectrum** (a 1–8 nm spread produces only minor broadening). That last point is close to a direct measurement of the central hypothesis here: if realistic polydispersity barely broadens the monomer spectrum, the large excess width observed at only 4% spread cannot be polydispersity — it must be coupling / speciation.

---

## A. Plasmon coupling in dimers/trimers — the broadening / red-shift mechanism

1. **Jain, El-Sayed et al., "On the Universal Scaling Behavior of the Distance Decay of Plasmon Coupling in Metal Nanoparticle Pairs: A Plasmon Ruler Equation," *Nano Lett.* 2007.** — *(reviewed)* The backbone equation. Fractional peak shift Δλ/λ₀ ≈ a·exp[−(gap/D)/τ] with a near-universal decay length τ ≈ 0.2 in units of particle diameter. Fast analytic map from interparticle gap → red-shift; the natural cheap "physics engine" inside a Monte-Carlo / inverse tool. Coupling is only significant at near-contact (gap ≲ 0.2–0.3·D); beyond ~1 diameter a pair is optically a monomer. https://pubs.acs.org/doi/10.1021/nl071008a

2. **"Aggregation affects optical properties and photothermal heating of gold nanospheres," *Sci. Rep.* 2020.** — *(open access)* DDA (DDSCAT) of small clusters N ≤ 5 in idealized close-packed contact, plus fractal aggregates (N = k_f (R_g/a)^D_f). Finds a secondary red-shifted "longitudinal-like" peak emerging in 2–4 particle clusters and **broadening + red-shift of the primary peak for all sizes** — precisely the observed spectral change. Also ties aggregation to reduced photothermal heating. Best single template for the forward optical model. https://www.nature.com/articles/s41598-020-79393-w

3. **"Surface Plasmon Coupling in Dimers of Gold Nanoparticles: Experiment and Theory for Ideal (Spherical) and Nonideal (Faceted) Building Blocks," *ACS Photonics* 2018.** — Experiment + theory on how real (faceted) dimers couple; relevant because 12 nm particles are slightly faceted and facet/gap geometry changes the shift. https://pubs.acs.org/doi/10.1021/acsphotonics.8b01424

4. **"Manipulating the confinement of electromagnetic field in size-specific gold nanoparticle dimers and trimers," *RSC Adv.* 2019.** — Systematic dimer vs trimer field/spectral behavior; useful for building distinct trimer basis spectra (linear vs bent vs triangular). https://pubs.rsc.org/en/content/articlehtml/2019/ra/c9ra07346a

5. **"Plasmon Coupling in Nanorod Assemblies: Optical Absorption, DDA Simulation, and Exciton-Coupling Model," *J. Phys. Chem. B* 2006.** — The plasmon-hybridization / exciton-coupling picture that lets you *predict* coupled-mode positions analytically instead of running full DDA every time. https://pubs.acs.org/doi/10.1021/jp063879z

6. **"Size Dependence of the Plasmon Ruler Equation for Two-Dimensional Metal Nanosphere Arrays," *J. Phys. Chem. C* 2011.** — How the ruler constant drifts with particle size / cluster order; needed to keep the analytic model honest at 12 nm. https://pubs.acs.org/doi/10.1021/jp2055415

---

## B. Reversible aggregation, spectral signatures & equilibrium validation

### B1 — Reversible temperature/stimulus-driven aggregation (reversibility + redshift/broadening signature)

*These establish that a monomer↔cluster switch produces the λmax redshift + FWHM broadening the forward model predicts, and that it is reversible. NONE of the three demonstrates an isosbestic point — do not cite them for that. Two are kinetically FROZEN, not equilibrium populations (caution for the equilibrium solver).*

7. **Kelesidis et al., *Anal. Chem.* 2022.** — *(reviewed in full)* Ensemble extinction ↔ agglomerate interparticle distance; closest precedent for reading cluster structure from a bulk UV–Vis spectrum. DOI: 10.1021/acs.analchem.1c05145

8. **Heo, Miesch, Emrick & Hayward, *Nano Lett.* 2013, 13, 5297.** — *(reviewed)* Solid-state: 4.5 nm P(S-r-2VP)-AuNPs in P(S-r-4VPh); phenol–pyridine H-bonds disperse, heating (200 °C) aggregates, cooling (120 °C) redisperses. λmax 527→537 nm + FWHM broadening, both ~linear in aggregate size; FTIR confirms H-bond break/reform. Reversible at 4.1 vol%, only partial at 14.2 vol%. Role: reversibility + shift↔size anchor (solid limit). DOI: 10.1021/nl402813q

9. **Durand-Gasselin, Sanson & Lequeux, *Langmuir* 2011, 27, 12329.** — *(reviewed)* Aqueous: 6 nm poly(EO-st-PO)-AuNPs, LCST aggregation above Tagg; free Pluronic P123 arrests growth to limited-size clusters only when BOTH Tagg and cmt are exceeded. λSPR linear in log(P123/AuNP). Reversible on cooling, but aggregates KINETICALLY FROZEN at fixed T (no rearrangement over weeks). Names dimer/trimer control as unsolved future work (= this regime). DOI: 10.1021/la2023852

10. **Trantakis, Bolisetty, Mezzenga & Sturla, *Langmuir* 2013, 29, 10824.** — *(reviewed)* Cleanest two-state dispersed↔aggregated system: 20 nm AuNPs, linker strand aggregates, complementary target disassembles at room temp, driven by duplex ΔTm (61.9 vs 31.3 °C); single-base specificity. Readout via A530/A700 ratio + SPR shift. NOT thermal (temp only tunes rate), NOT isosbestic, clusters kinetically formed. Role: two-state reversible switch with a well-defined thermodynamic driver. DOI: 10.1021/la401211u

### B2 — Isosbestic points & thermal multimer equilibrium (equilibrium-model validation)

*Backs the van 't Hoff + mass-balance + isosbestic-pivot assumptions in the association-constant / spectra layer. Verified against full PDFs. NB: no single paper combines thermal + isosbestic + small-N — see Gap note.*

11. **Cardellini, Marchetti, et al., "Nanoplasmonic Isosbestics Uncover Mesoscale Assembly of Gold Nanoparticles on Soft Templates," *J. Am. Chem. Soc.* 2025, 147, 20008–20022.** — *(reviewed in full)* PRIMARY isosbestic-mechanism ref. Experiment + electrodynamics simulation prove nanoplasmonic isosbestics arise from a conserved two-species equilibrium (individual AuNP ↔ clustered AuNP), and that the isosbestic WAVELENGTH is a sigmoidal function of interparticle spacing (λiso = 537 nm @ ⟨sp⟩ = 0.3 nm → 527 nm @ 0.8 nm; ⟨sp⟩ ≈ 0.3 nm recovered for AuNP–DOPC). ~12 nm citrate AuNPs templated on liposomes of tuned stiffness. Directly justifies "isosbestic position → geometry" in the inverse model. Equilibrium is TEMPLATE/CONCENTRATION-driven, NOT thermal. DOI: 10.1021/jacs.5c05189

12. **Mezzasalma, Kruse, Iturrospe Ibarra, Arbe & Grzelczak, "Statistical Thermodynamics in Reversible Clustering of Gold Nanoparticles. A First Step Towards Nanocluster Heat Engines," *J. Colloid Interface Sci.* 2022, 628, 205–214.** — *(reviewed in full)* PRIMARY thermal-equilibrium ref. 13.4 nm BSPP AuNPs, reversible thermal clustering (UV–Vis–NIR recovery on T-reversal). Validates the van 't Hoff EXTRACTION METHODOLOGY (transient equilibrium mapping → ΔH⁰ = 188 kJ/mol, ΔS⁰ = 461 J/(mol·K) on the dispersion branch), NOT a plug-in two-state K. Caveats: INVERSE-T (clusters on cooling, disperses >40 °C) and genuine thermal HYSTERESIS (framed as a thermodynamic cycle). No isosbestic shown. DOI: 10.1016/j.jcis.2022.07.037

13. **Kruse, Rao, Sánchez-Iglesias, Montaño-Priede, Seifert, Grzelczak et al., "Temperature-Modulated Reversible Clustering of Gold Nanorods Driven by Small Surface Ligands," *Chem. Eur. J.* 2024, 30, e202302793.** — *(reviewed in full)* OPTIONAL / anisotropic extension. Same group as ref 12. Thermoresponsive AuNR clustering with small BSPP ligands (no polymer/DNA); real-time longitudinal LSPR intensity → interparticle distance, cluster size, melting/freezing temps, hysteresis. Notes the roster of thermoresponsive reversible-clustering systems is "surprisingly short" (supports novelty). Rod-specific readout, no isosbestic. DOI: 10.1002/chem.202302793

### Gap note

No reference here combines a thermally-driven isosbestic multimer equilibrium in the small-N (dimer/trimer) regime. Ref 11 = isosbestic mechanism (concentration-templated); refs 12/13 = thermal equilibrium (no isosbestic). The combined thermal + isosbestic + few-particle demonstration is an open niche — which reinforces the tool's novelty claim. Ref 11's λiso↔⟨sp⟩ sigmoidal is the closest existing anchor for the inverse map, but it is ensemble-averaged over a broad cluster population, not a per-N basis, so it calibrates direction/geometry, not a species-resolved basis spectrum.

---

## C. Size & concentration from UV–Vis — the classical (monomer-only) tool foundation

14. **Haiss, Thanh et al., "Determination of Size and Concentration of Gold Nanoparticles from UV−Vis Spectra," *Anal. Chem.* 2007.** — *(reviewed in full)* The canonical calibration (peak position & the ratio A_spr/A_450 → diameter; extinction → concentration). The baseline "monomer-only" method being extended here. **Verified details that matter for this project:**
    - **Permittivity:** bulk-gold n(λ) from **Johnson & Christy (1972)**, spline-fitted for continuous evaluation, then **corrected for the reduced mean free path** of conduction electrons in small particles. (Bulk measurement + analytic small-particle correction — contrast with ref 18's size-matched thin-film approach.)
    - **Peak-position method fails below ~25 nm:** theory matches experiment for **25–120 nm**, but below 25 nm the observed peak sits lower than predicted; attributed to the rising surface-to-bulk atom ratio below 20 nm. (Independently reproduced here: normalized single spectra of 8–31 nm AuNP are near-degenerate in size.)
    - **Haiss already knew the red side is aggregation-contaminated.** They state that long-wavelength absorbance ratios can be used to *assess hydrosol quality or monitor aggregation* but are **not suitable for size determination**, and that theory–experiment agreement improves for ratios taken **below 600 nm** — which is exactly why the reference is 450 nm. → The A_spr/A_450 design is a deliberate *avoidance* of the red-side signal this project aims to *model*. Strong support for the central thesis, and a direct citation for "red-side extinction is unmodeled in monomer-only tools."
    - **Stated limitation:** calculations assume **monodisperse, perfectly spherical** particles.
    https://pubs.acs.org/doi/10.1021/ac0702084

15. **"Determination of Size and Concentration of Gold Nanoparticles from Extinction Spectra," *Anal. Chem.* 2008.** — Full-spectrum (not just peak) inversion using Mie theory; the natural forward model to fit against. https://pubs.acs.org/doi/10.1021/ac800834n

16. **"Size Evaluation of Gold Nanoparticles by UV−vis Spectroscopy," *J. Phys. Chem. C* 2009.** — Complementary calibration and limits of the peak-based approach. https://pubs.acs.org/doi/10.1021/jp8082425

---

## D. Machine-learning / inverse recovery of size distribution — the "replace TEM" precedent

17. **Klinavičius et al., "Deep Learning Methods for Colloidal Silver Nanoparticle Concentration and Size Distribution Determination from UV–Vis Extinction Spectra," *J. Phys. Chem. C* 2024.** — *(reviewed in full)* Tandem DNN: "DipoleNet" isolates the dipolar LSPR component, "ColloidNet" maps it → log-normal size distribution + concentration; trained on effective-medium (Maxwell-Garnett-Mie) / Mie simulated spectra; μ recovered to ~1.2% (large NPs), 6.1% overall. **Two caveats that matter here:** (i) all colloids were **sonicated before UV–Vis and TEM to break up agglomerates** — i.e. good matches are in a forced-monomer state; (ii) the paper explicitly warns aggregation red-shifts, broadens, and lowers the spectrum, causing over-predicted size / under-predicted concentration. This is the *control experiment* for the hypothesis, not a counterexample. Also: σ (distribution width) was consistently **over**estimated. https://pubs.acs.org/doi/10.1021/acs.jpcc.4c02459

18. **Klinavičius, Tamulevičienė & Tamulevičius, "How to Select the Proper Gold Permittivity for Generating Training Data for Deep Neural Networks for Analysis of the Nanoparticle Colloid Size Distributions from Their UV–Vis–NIR Extinction Spectra?," *J. Phys. Chem. C* 2025, 129, 17616–17631.** — *(reviewed in full)* **The permittivity-choice paper — directly actionable for `dielectric.py`.** Key findings:
    - Gold permittivities measured by different authors **vary significantly — especially the imaginary part ε₂** — and this has a **profound effect** on DNN training and prediction accuracy. ε choice is not a detail.
    - **Best performers for colloidal AuNP were NOT bulk Johnson & Christy** but **Yakubovsky's thin-film measurements**, applied *size-matched*: the **25 nm-film ε for particles below 25 nm radius, the 53 nm-film ε for larger**. Lemarchand's 11 nm-film data also performed well.
    - Of all permittivities tested, those with **medium-to-high LSPR quality factor** performed best.
    - Philosophy: a thin film's *measured* ε already embeds some reduced-mean-free-path/surface-scattering physics that a nanoparticle experiences — obtained empirically rather than via a correction formula (contrast ref 14's bulk-J&C + analytic correction).
    → **Consequence for this tool:** `dielectric.py` currently offers `etchegoin` / `bb` / `jc`. A **`yakubovsky` (size-matched thin-film) option should be added** and tested; a wrong ε₂ directly changes damping and therefore the red tail, so ε is a **confound that must be ruled out before attributing static red-side excess to speciation**. (Note: the *temperature-dependent, reversible* part of the red-side growth cannot be a dielectric artifact.)
    https://pubs.acs.org/doi/10.1021/acs.jpcc.5c03973

19. **Bilén et al., "Machine Learning-Based Interpretation of Optical Properties of Colloidal Gold with Convolutional Neural Networks," *J. Phys. Chem. C* 2024.** — *(reviewed in full)* CNN takes a UV–Vis spectrum → predicts mean diameter + polydispersity (σ) of a Gaussian, for spheres and rods; trained on in-silico spectra. **Forward model = MEEP/FDTD single-particle spectra + additive mixing of 100,000 sampled sizes, explicitly "not considering… plasmon coupling."** Two key results: (i) mean diameter strongly imprints on the spectrum but **polydispersity barely does**, and becomes the hardest parameter under noise; (ii) the sim-vs-experiment gap is blamed on "not accounting for capping agents or plasmon coupling." The authors also hit FDTD numerical instabilities on some spectra. https://pubs.acs.org/doi/10.1021/acs.jpcc.4c02971

20. **Glaubitz et al., "Leveraging Machine Learning for Size and Shape Analysis of Nanoparticles: A Shortcut to Electron Microscopy," *J. Phys. Chem. C* 2024.** — *(reviewed)* Explicitly frames optical→size ML as a TEM replacement; reports accuracy vs EM. The framing paper. https://pubs.acs.org/doi/10.1021/acs.jpcc.3c05938

---

## E. Cluster-population statistics & simulation engines (for the Monte-Carlo core)

21. **Kim et al., "Kinetics of Gold Nanoparticle Aggregation: Experiments and Modeling," *J. Colloid Interface Sci.* 2008.** — Smoluchowski population-balance modeling of AuNP cluster-size distributions; the kinetic counterpart to the equilibrium MC. https://pubmed.ncbi.nlm.nih.gov/18022182/

22. **"Modeling Nanoparticle Aggregation," *Chemical Modelling* (RSC book chapter).** — Review of DLVO + Brownian / Monte-Carlo aggregation modeling; a menu of methods for generating cluster ensembles. https://books.rsc.org/books/edited-volume/1453/chapter/953996/Modeling-nanoparticle-aggregation

---

## Deep-read: Kelesidis et al. 2022 (ref 7) — the closest forward model, and where it stops

**What it does.** Two-stage physics pipeline. (1) Discrete element modeling (DEM) simulates Brownian coagulation of 1000 monodisperse Au primary particles (dp = 20–80 nm) in water at 300 K → realistic fractal agglomerates of 2–15 particles with quantified morphology (fractal dimension Df, prefactor, radius of gyration). (2) DDA (DDSCAT, **100,000 dipoles/particle**) computes extinction, sweeping interparticle distance s = 1–50 nm. Regresses longitudinal-peak λ_l vs s/dp into a new plasmon ruler (their eq. 4), more accurate below 10 nm than the older Jain-type dimer ruler (which underestimates the redshift by ~40% at small s).

**Findings that shape the forward model:**

- **Contact is everything.** At s = 50 nm a "dimer" is optically identical to a single sphere at all λ. Coupling switches on only at near-contact (consistent with Jain's 0.2·D decay length). For 12 nm citrate particles with a ~1 nm ligand shell, only genuine **contact** dimers/trimers contribute red-shifted intensity → the equilibrium model should be a *contact-association* model, not a general pair-distance integral.
- **λ_l encodes the gap, not the cluster order.** 7-mer → 15-mer shifts λ_l by only ~10%; increasing Df barely moves λ_l. What Df *does* control is **broadening**. Consequence: peak *position* is nearly degenerate in cluster order — so dimer/trimer *fractions* must be read from **relative amplitude of the coupled feature + lineshape breadth + monomer:coupled intensity ratio**, i.e. full-spectrum information (exactly what the monomer-only NNs discard).
- **Resolution trap.** Coarse DDA / T-matrix underestimate the redshift by up to 55%; the authors needed 100k dipoles/particle. → For spheres, prefer an exact multisphere method (MSTM) over DDA to avoid this entirely.

**Where it stops (all in this project's favor):** calibrated for 20–80 nm (not 12 nm — coupling ∝ volume ∝ D³, so at 12 nm expect a modest red **shoulder/tail near ~550 nm, not a resolved NIR peak**, which matches the observed broadening); reads s for a *single known agglomerate* rather than decomposing a bulk ensemble into monomer/dimer/trimer *fractions*; and has *no* thermodynamics (static coating-set gap, no equilibrium, no isosbestic point). This project's contribution is precisely these three extensions.

---

## Optical simulation layer — solver choice (MEEP/FDTD vs DDA vs MSTM vs CDA)

The key simplifier: **12 nm ≪ λ (~520 nm), so we are deep in the quasi-static regime.** That makes the cheap methods accurate and the heavy ones unnecessary as the workhorse. Ladder from cheapest to heaviest:

- **Mie theory (analytic) — monomers.** Exact closed form (`PyMieScatt` / `miepython` / `scattnlay`); instant; polydispersity is just an integral over the size distribution. No reason to use FDTD for single spheres.
- **Coupled-dipole approximation (CDA) — dimers/trimers, recommended first build.** At 12 nm each particle ≈ one dipole (Mie a₁ / Clausius-Mossotti polarizability); a cluster is a 2–3 dipole coupled linear system, solved analytically, with cheap / closed-form orientation averaging. Excellent at 12 nm (strains only at large sizes). Essentially free → can live inside the fit loop.
- **Multisphere T-matrix (MSTM / GMM) — accurate reference for sphere clusters.** Semi-analytic, essentially *exact* for spheres, with **analytic orientation averaging** (ideal for randomly-oriented clusters in solution). Mackowski's MSTM (Fortran) or Python `treams` / `smuthi`. Strictly better than DDA *for spheres* (no staircasing, far cheaper for equal accuracy).
- **DDA (DDSCAT / ADDA / pyGDM) — only if shape matters.** Its purpose is arbitrary shapes; for spheres it pays a staircasing penalty (Kelesidis needed 100k dipoles/particle). Use only for faceting / neck formation.
- **FDTD (MEEP) — most general, most expensive.** Time-domain grid solver: one pulse → full broadband spectrum, any geometry, any dispersive material. Costs for *this* problem: grid staircasing of curved surfaces (needs fine resolution + subpixel smoothing + convergence study); gold enters as a Drude–Lorentz fit that must be validated vs Johnson–Christy in 450–700 nm; metal FDTD is finicky (sharp surface fields, small time steps, careful PML — the numerical instabilities Bilén reported); **no analytic orientation averaging** (every gap × orientation is a separate 3D run); minutes–hours per run. Fine for building a training library *offline*; too slow inside an iterative fit.

**Recommendation (hybrid).** Primary engine = Mie (monomers) + CDA/MSTM (clusters): exact for spheres, orientation-averages analytically, fast enough to precompute a basis or run in the fit loop. Use **MEEP for the two jobs it's genuinely best at**: (1) independent cross-validation of the MSTM/CDA cluster spectra (frequency-domain semi-analytic vs time-domain grid agreement = strong check), and (2) the escape hatch if non-spherical effects (faceting, flattened contact facets, partial coalescence) turn out to matter. Since Bilén (ref 19) showed FDTD produces usable training libraries, MEEP is a reasonable offline library / validation tool — just keep FDTD out of the inner fit loop.

---

## Connections to related work

- **Jain et al. 2007 (plasmon ruler)** — item A1; the analytic heart of the scheme.
- **Hagan 2008 / Elrad & Hagan 2008 (viral capsid assembly)** — methodological analogs: statistical-mechanics / Monte-Carlo of reversible self-assembly with size control and polymorphism. The "distribution of assembly states in equilibrium" formalism transfers directly to monomer/dimer/trimer populations.
- **Photothermal detection literature (Selmke & Cichos; Adhikari; and related thesis work)** — the *Sci. Rep.* 2020 paper (A2) links aggregation state to photothermal heating, so the multimer model also predicts a photothermal-contrast signature to cross-check against photothermal detection methods.

---

## Proposed tool architecture (three layers)

**Layer 1 — Population model (the Monte Carlo).** (a) *Reaction / equilibrium MC*: monomer⇌dimer⇌trimer⇌… as an association network with equilibrium constants K₂, K₃… (or a single binding free energy ΔG + coordination penalty). Sample populations; temperature enters via K(T) = exp(−ΔG/kT), which *predicts the isosbestic point*. Because coupling is a contact effect (Kelesidis), model *contact* association specifically. (b) *Configurational MC*: particles in a box with DLVO + van-der-Waals + ligand-steric potential; Metropolis / Brownian sampling → realistic cluster geometries. More physical, more expensive. Start with (a) for a first version.

**Layer 2 — Optical forward model (spectrum of each species).** Precompute orientation-averaged extinction spectra: monomer (Mie), contact dimer, trimer (linear/bent/triangular), and a few multimer motifs — using **CDA** (fast, great at 12 nm) with **MSTM** as the exact reference and **MEEP/FDTD** as cross-check. Convolve each with the 4% size polydispersity. Total: σ(λ) = Σᵢ nᵢ · ⟨σᵢ(λ)⟩, weighted by Layer-1 populations. (Note: at 12 nm expect red *shoulders/tails*, not resolved NIR peaks.)

**Layer 3 — Inverse fit (the deliverable tool).** Fit an experimental spectrum for {nominal size, size-distribution width, multimer fractions / K's}:

- *Least-squares / Bayesian* over the Layer-2 basis (interpretable, error bars, small data) — best first version.
- *Neural network* trained on Layer-1+2 simulated spectra (refs 17–20 do this for monomers; extend the training set with cluster basis spectra + speciation labels). Instant inference once trained.

**Validation plan:** (i) reproduce the measured broad spectrum and its temperature series *including the isosbestic point*; (ii) recover known size distributions on samples with matching TEM; (iii) show the fitted monomer size distribution matches TEM while the residual broadening is absorbed by multimer fractions — the scientific punchline: **UV–Vis width ≠ polydispersity; it's polydispersity + speciation.**

---

## Empirical findings from lab data (analysis sessions, 2026-07)

Two in-house datasets were analyzed against the monomer forward model. Both show
**red-side extinction in excess of clean-monomer Mie prediction** — the effect the
literature above treats as an error term. A chain of controls has since *excluded*
most of the mundane explanations; what remains is a genuine, particle-derived,
reversible scattering population of **unexpectedly large** scatterers.

### Dataset 1 — C500 temperature series
Cary 100, 2011; citrate/TEG AuNP, TEM 12.9 nm ± 7%, 45× diluted into the linear
absorbance regime from the optical-force stock. Reproduces report Figs 16–17.

- **Normalization:** multiplicative anchor at **400 nm** (gold ≈ dielectric there →
  extinction tracks total cross-section, ~speciation-invariant; Beer–Lambert makes
  the correction multiplicative, not additive). Removes ~3% inter-scan drift.
- **Isosbestic crossing ≈ 575 nm**, recovered only after 400 nm normalization
  (un-normalized data gives a spurious ~542 nm crossing). Consistent with a
  predominantly two-state monomer ⇌ dimer equilibrium; slight red drift of the
  crossing at the highest T suggests minor higher-order aggregation.
- **The red wing is a near-flat "pedestal", not a coupled-plasmon shoulder.**
  Red-tail/peak (700/523 nm) ≈ 0.18→0.21 growing with T, still ~0.16–0.19 at 790 nm.
  Best basis red-tails fall an order of magnitude short (CDA dimer ~0.010;
  **exact T-matrix** dimer ~0.021, linear trimer ~0.023 at 790 nm). The
  tail-minus-pedestal is ~T-constant ⇒ **the growth is the pedestal**, and it is
  **reversible on both heating and cooling branches**.
- Fitting the raw shape without a pedestal term rails D/polydispersity (the fit
  fakes the pedestal with size), corrupting BOTH size AND the species fractions
  that are the end goal.

### The pedestal: what it is NOT (four independent exclusions)

**(a) NOT small-aggregate Rayleigh scattering.** A λ^(−n) scattering term was added
to the fit with **n fitted, not fixed** (bounded [0,6], started neutrally at 2) as an
explicit falsifiable test: n ≈ 4 ⇒ Rayleigh from small aggregates; n ≈ 0 ⇒ flat
offset. **Measured: n_sca ≈ 0.63–0.79** across both branches and both ε bases —
**decisively not Rayleigh.** The pedestal decays far too weakly for sub-λ scatterers.
n ≈ 0.7 points to scatterers **≳ λ** (i.e. hundreds of nm to µm), where Mie
scattering flattens toward λ⁰. *The small-multimer-Rayleigh hypothesis is rejected
by the data; it was not tuned to fit the story.*
A_sca is reversible (tracks T on both branches, ~+2–3% over 15→75 °C, hysteresis
~1–2.5%). Adding the term: RMS 0.037 → 0.0225, agg fraction 0.41–0.47.
⚠️ **Open degeneracy:** the near-flat pedestal trades against diameter — D falls to
11.2 ± 0.7 / 12.2 ± 0.6 nm (~2σ below TEM 12.9), vs 12.5–12.6 under the old ad-hoc
baseline subtraction. **D, A_sca, and the species fractions are therefore not yet
independently trustworthy.**

**(b) NOT a permittivity artifact.** Yakubovsky size-matched thin-film ε (ref 18)
*does* raise the model red tail vs J&C (+34% @700 nm, +20% @790 nm) — so ε error is
a real contributor — but it accounts for only **~1–2% of the measured pedestal**, and
re-estimating on the yk25 basis moved A_sca by just ~−5%. Trade-off found:
yakubovsky's peak sits **7.5 nm blue** of experiment at 12.9 nm, which the global fit
converts into *fake aggregation* (agg 0.7–0.8) ⇒ **`jc` stays production for the
13 nm system; yakubovsky is the systematics bound and the ≥50 nm option.**

**(c) NOT instrumental, and NOT thermal-optical artifact of the cell.** The 2011
`raw` sheet contains **2–4 water-blank scans at every temperature**, interleaved
before each sample scan — a proper blank temperature series. **The blanks are flat
with T:** far-red blank level drifts by only **−0.0016 (700 nm) / −0.0020 (790 nm)**
over 15→75 °C, i.e. **~15× smaller than the sample pedestal's ≈ +0.03 growth, and in
the OPPOSITE direction.** This experimentally excludes, for the *temperature-dependent*
component: stray light, convection/schlieren in the heated cuvette, and gas-bubble
formation from dissolved gas (pure water carries the same dissolved-gas load and
shows nothing). *A static instrumental offset (~0.03 abs) is present in the far red
and does contribute to the pedestal's absolute baseline — but not to its T-dependence,
which is where the science is.*

**(d) NOT nanobubbles on the particle surfaces.** Checked against the steady-state
(non-laser) surface-nanobubble literature — see section F. Every specific detail cuts
against it: nanobubbles require **hydrophobic** surfaces (citrate/TEG AuNP are
deliberately hydrophilic), **high curvature suppresses** their formation on gold,
observed surface nanobubbles are **120–145 nm** (far too large to sit on a 13 nm
sphere), and **heating to 40–60 °C *reduces* their number** rather than growing them —
the opposite of the pedestal's behavior. Note the *laser-driven* plasmonic-nanobubble
literature is a different regime entirely (particle superheated far above the bulk
liquid) and must not be conflated with mild bulk heating to 75 °C.

**⇒ What remains:** a genuine, nanoparticle-derived, reversible population of
scatterers **comparable to or larger than the wavelength** (≳ several hundred nm).
Not the dimers/trimers being speciated — a distinct, larger channel. Small clusters
shape 520–600 nm; this population makes the far-red floor.
**Decisive next test is experimental, not more fitting:** DLS across the temperature
ramp (does µm-scale material appear reversibly on heating?), and/or a
filtered/centrifuged control (does the pedestal survive removal of large material?).

### Dataset 2 — CTAC size series
Irina; 5 samples, TEM 7.8 / 29.2 / 31.5 / 37.9 / 42.0 nm, 4–9% polydispersity, single
room-T spectra; well-dispersed CTAC-capped monomers.

- **Confirms the red-tail excess in "good" samples.** With diameter **fixed to TEM
  ground truth** (amplitude the only free parameter), the data ride above the Mie
  curve across ~560–660 nm in **every** sample — including the pristine 29/31 nm ones.
  Independently reproduces the persistent-red-tail-vs-simulation effect reported in
  the lab, in a *different* surface chemistry and a *different* sample set.
- **Single-spectrum size retrieval fails below ~35 nm** (shape-only fit: TEM 7.8 →
  15.7 nm; 29.2 → 19.7; 31.5 → 20.2; 37.9 → 30.7; but 42.0 → 40.2 ✓). An *information*
  limit, not a bug — matches ref 14's own ≥25 nm validity statement.
- **The Haiss A_spr/A_450 ratio is only a partial rescue:** it separates 8 nm from
  ~30 nm but **saturates across 29–42 nm** (ratios 2.39–2.43), and the measured ratios
  sit **above** clean-monomer Mie — the excess red tail inflates A_spr, so **the red
  tail biases the Haiss size readout upward**, exactly the contamination ref 14 itself
  warned about.
- **Scattering channel validated across size** (Mie, jc): σ_sca/σ_ext at peak = 0.4%
  (8 nm) → 0.8% (29–31) → 3.0% (38) → **6.7% (42 nm)** — the expected (size)⁶-vs-(size)³
  scaling. **For the planned ~60 nm work** scattering is 15–30% of extinction and even
  *monomers* scatter appreciably, which will partially degenerate with the A_sca term
  and make the pedestal/species separation harder than it is at 13 nm.

### Known model gaps exposed by this analysis
- **Medium index is frozen and purely real.** `WATER_N = 1.333`, no T-dependence and
  **no imaginary part** — water is modeled as lossless. Water's index falls ~1.3334
  (15 °C) → ~1.3200 (75 °C), and water genuinely absorbs in the red/NIR (visible in the
  blank scans as a rise toward the red with a bump near 740–760 nm). Both should be
  added. **Falsifiable prediction: neither will explain the pedestal** — a ~1% index
  change shifts the *resonance* by ~1 nm; it cannot create a flat far-red floor. If
  adding n(T) absorbs a meaningful fraction of A_sca, something is badly misunderstood
  and it must be flagged loudly.
  (Note the medium n(T) effect is a *blue* shift on heating, opposite in sign to the
  observed peak red-shift — so it is currently *masking* part of the speciation signal.)

---

## F. Steady-state (non-laser) surface nanobubbles — the excluded hypothesis

Consulted to test whether the reversible T-dependent pedestal could be gas bubbles on
the AuNP surfaces. **Conclusion: no** (see exclusion (d) above). Retained because the
exclusion is load-bearing for the interpretation.

21. **"Surface nanobubbles: Theory, simulation, and experiment. A review," *Adv. Colloid Interface Sci.* 2019.** — Establishes that surface nanobubbles **do** form spontaneously by simple immersion, without lasers; gas concentration ~100–110% and liquid temperature ~25–45 °C favor formation; lifetimes reach hours-to-days, far exceeding the Epstein–Plesset prediction from the high Laplace pressure. https://www.sciencedirect.com/science/article/abs/pii/S0001868619300685

22. **"On the stability of bubbles trapped at a solid–liquid interface: A thermodynamical approach."** — The phase-space study on a **hydrophobised silicon** substrate with independently controlled temperature and gas concentration; NBs occupy a distinct region at ~100–110% gas concentration; micropancakes appear at higher T and gas concentration; supersaturation is **not** required for nucleation.

23. **"Solid-liquid interfacial nanobubble nucleation dynamics influenced by surface hydrophobicity and gas oversaturation," *ScienceDirect* 2024.** — **The load-bearing exclusion.** MD simulations show a *distinct preference* for interfacial nanobubbles to form on **hydrophobic** rather than hydrophilic surfaces, agreeing with experiment (gas enriches at the hydrophobic surface, disrupting the hydration layer). Citrate/TEG-capped AuNP are deliberately **hydrophilic** ⇒ wrong surface chemistry. https://www.sciencedirect.com/science/article/abs/pii/S0167732224018178

24. **"Thermal Conductance of the Gold–Water Interface…," PMC12814529.** — States directly that **high interfacial curvature enhances conductance, SUPPRESSES nanobubble formation**, and maintains efficient heat transfer even at high temperatures. Supports the Young–Laplace argument (2γ/r makes nm-radius bubbles hard to sustain). ⚠️ *Caveat: derived in a **laser-processing** context; the curvature result is general interfacial physics, but the transfer to mild bulk heating is supporting, not decisive.* https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12814529/

25. **"Experimental investigation on nucleation of surface nanobubbles with polycarbonate hydrophobic surface tested under heterogeneous boiling," *ScienceDirect* 2023.** — AFM before/after heating water to 40/50/60 °C: the **number of bubbles DECREASED by ~5%** after heating; average surface nanobubble size **120–145 nm**; size "not considerably changed" on heating. ⇒ Wrong size scale (a 13 nm sphere cannot host a pinned 120 nm spherical cap) and wrong temperature response (heating reduces, not grows). ⚠️ *Caveat: flat **polycarbonate** substrate, not AuNP.* https://www.sciencedirect.com/science/article/abs/pii/S2214785323030249

**Evidence-quality note:** refs 21 & 23 (existence + hydrophobicity requirement) are
solid and directly on point; refs 24 & 25 are corroborating but come from adjacent
systems (laser-heated gold; flat polycarbonate). **No direct study of bubble formation
on citrate-capped AuNP under mild bulk heating appears to exist** — that gap is itself
worth knowing. Lead with the hydrophobicity argument and the in-house blank control
(exclusion (c)), which is the strongest single piece of evidence.

---

*Access status (2026-07-10): refs 7, 14, 17, 18, 19, 20 reviewed in full; refs 21-25 (nanobubble section F) consulted via search snippets, not full text
(ref 14 Haiss and ref 18 gold-permittivity paper now obtained and read).
Remaining summaries rely on abstracts + open-access sources.*
