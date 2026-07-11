# `aunp_speciation` — Modelling Stack: Choices, Alternatives, and Rationale

*Status: 2026-07-11. Written to document why each layer of the forward model is what it is,
what else was available, and what was rejected. Open issues are flagged inline.*

---

## Purpose of the model

Given a UV–Vis extinction spectrum (or a temperature series) of a gold nanoparticle colloid,
recover the **speciation** — what fraction of the gold is present as isolated monomers versus
dimers, trimers, and larger aggregates. The scientific motivation is that measured AuNP spectra
are consistently **broader and more red-tailed than a monodisperse-monomer model predicts**, even
for samples with excellent TEM-verified size uniformity. The classical tools (Haiss et al. 2007)
and the modern ML successors (Klinavičius 2024; Bilén 2024; Glaubitz 2024) all *assume* a pure
monomer population and treat this excess as an error term — sonicating it away or absorbing it
into a residual. This tool aims to model it instead.

The forward model is a sum over species, each with its own orientation-averaged extinction
cross-section, weighted by an equilibrium population, plus a broadband scattering baseline:

    A(λ,T) = c(T)·ℓ · Σ_s f_s(T)·σ_s(λ, T)  +  A_sca(T)·(λ/550)^(−n)

Everything below is a defence of the individual terms.

---

## Layer 1 — Gold permittivity ε_Au(λ)

**Chosen: tabulated Johnson & Christy (`jc`), cubic-spline interpolated.**

| option | verdict |
|---|---|
| `etchegoin` — analytic Drude + 2 critical points | ❌ peak at 526.9 nm; **4–6 nm too red** |
| `bb` — Brendel–Bormann (Rakić 1998) tabulated fit | ❌ same 4–6 nm red bias |
| **`jc` — Johnson & Christy (1972) measured n,k** | ✅ **522.8 nm vs 523.1 nm measured** |
| `yakubovsky25/53` — thin-film n,k, size-matched | ⚠️ peak 7.5 nm **blue**; the fit converts the mismatch into *fake aggregation*. Retained as a systematics bound and the ≥50 nm option. |

The analytic models are convenient (evaluable at any λ, no data files) but both place the plasmon
several nm red of experiment, because their ε₁ is slightly off near 520 nm. Since the entire
scientific claim rests on excess red-side extinction, a systematic red bias in the *baseline* model
is disqualifying.

Klinavičius et al. (2025) is the dedicated study of this question and recommends **Yakubovsky
thin-film ε, size-matched**, over bulk J&C — on the grounds that a thin film's measured ε already
embeds surface-scattering physics that bulk data needs a correction formula for. We tested it and
found it lands 7.5 nm blue at 12.9 nm, which the global fit then "explains" by inventing
aggregation. **We therefore keep `jc` as production for the 13 nm system and hold `yakubovsky` as
a documented systematic.** This disagreement with the literature is deliberate and recorded.

---

## Layer 2 — Size / surface damping (the Kreibig term)

**Chosen: γ(R) = γ_bulk + A·v_F/R, with A = 0.25.**

Below ~50 nm a gold particle is smaller than the bulk electron mean free path (~42 nm), so
conduction electrons scatter off the *surface*, broadening the plasmon. **A** is the dimensionless
probability that this scattering is diffuse.

- **A = 1** is the fully-diffuse limit — a theoretical *maximum*, not a measurement. The code
  originally defaulted here.
- **A = 0.25** is the measured value: Berciaud, Cognet, Tamarat & Lounis (*Nano Lett.* 2005),
  using **Photothermal Heterodyne Imaging of individual nanoparticles** (5/10/20/33 nm, >30
  single-particle spectra per size). Single-particle measurement gives the **homogeneous**
  linewidth, free of the polydispersity *and* aggregation broadening that make any ensemble
  determination of A degenerate. ⚠️ *Secondary literature widely misquotes this as 0.33; it is 0.25.*

**Why it matters more than it looks.** At D = 12.9 nm:

| A | γ_S | γ_bulk | which dominates? |
|---|---|---|---|
| 1.00 | 0.143 eV | 0.044 eV | γ_S, 3:1 |
| **0.25** | **0.036 eV** | 0.044 eV | **γ_bulk — it inverts** |

This determines whether the plasmon linewidth is controlled by the *surface* (temperature-independent)
or by the *bulk* (strongly temperature-dependent) — i.e. whether the observed thermal changes belong
to gold or to aggregation. **Open issue: the code currently ships A = 1.0. This is scheduled to change.**

---

## Layer 3 — Temperature dependence of gold, ε_Au(λ, T)

**Chosen: a bracketing pair, with fitted quantities reported as ranges — the production
default is a bulk thermal Drude retune (Holstein form, Olmon-anchored) on the J&C
baseline; the full Reddy et al. (2016) Drude + two-critical-point (DCP) model,
interpolated in T, is the bracketing variant.** The Reddy DCP fits better (RMS −17%)
but is thin-film-derived, and fit quality is deliberately not used as the selector
(see Layer 1's Yakubovsky case). The DCP defines the ε(T) *magnitude* argument below.

Gold's permittivity is *not* temperature-independent: heating raises the electron–phonon scattering
rate, increasing ε₂, which **broadens the plasmon, lowers its peak, and red-shifts it — reversibly,
with no aggregation whatsoever**. This is precisely the signature we were attributing to dimerization,
and the model originally ignored it entirely.

The damping decomposes as **γ(T) = γ^(ep)(T) + γ_S** (Chetoui et al. 2026), where the electron–phonon
term is a **bulk** property (strongly T-dependent) and the surface term γ_S is **temperature-independent**.
The two are additive and separable — so the correct recipe is **bulk ε(T) + the existing size correction**,
*not* a monolithic "nanoparticle ε(T)" (which would import another sample's γ_S and environment).

**Drude-only was insufficient.** Reddy's tables show the *interband* critical points also move with
temperature (ε∞ 2.27→2.45; the 2.62 eV CP damping 0.256→0.273 over 23→100 °C). Measured effect on the
12.9 nm peak over the 15→75 °C ramp:

| treatment | peak change | vs observed (−2.7%) |
|---|---|---|
| Drude-only | −0.33% | ~12% |
| **Full Reddy DCP** | **−2.56%** | **~95%** |

**The interband terms dominate.** Nearly all of the "monomer loss at 523 nm" that motivated the
two-state interpretation is gold thermal physics.

⚠️ **Open systematic.** The measured dΓ_D/dT spans a factor of ~6 depending on film type
(Reddy poly-crystalline 0.23e-4 eV/K; Reddy single-crystal 1.08e-4; Holstein theory ~1.45e-4).
Citrate-grown AuNPs are polycrystalline, arguing for the low end. This was previously dismissed as
"immaterial because γ_S dominates" — a conclusion that fails once A = 0.25 (Layer 2). **To be re-bracketed.**

---

## Layer 4 — Cluster optics: **the central choice**

**Chosen: exact multi-sphere T-matrix (`treams`).**

| method | accuracy for sphere clusters | cost | verdict |
|---|---|---|---|
| **CDA** — 1 dipole per particle | **fails** | fastest | ❌ see below |
| **DDA** — N dipoles per particle | good, but discretisation error | slow | ❌ strictly worse than T-matrix for spheres |
| **FDTD** (MEEP etc.) | general, meshing error | slowest | ❌ same |
| **Maxwell–Garnett effective medium** | mean-field; **discards near-field coupling** | trivial | ❌ throws away the physics of interest |
| **T-matrix / MSTM** | **exact** (analytic multipole) | fast *after* caching | ✅ |

**Why not DDA or FDTD?** For clusters of *spheres*, the multi-sphere T-matrix is **exact** — it solves
the full multipole scattering problem analytically, with no discretisation. DDA approximates a sphere
with a cubic dipole lattice; FDTD with a Yee mesh. Both introduce shape error and are slower. DDA and
FDTD earn their keep for *arbitrary* geometries (rods, stars, faceted particles) — for spheres they are
a downgrade. MSTM additionally provides **analytic orientation averaging**, which is exactly what a
cuvette of tumbling clusters requires.

**Why CDA failed — quantified.** The coupled-dipole approximation is DDA in its crudest limit (one
dipole per particle). It cannot represent the higher multipoles that dominate when surfaces are ~1 nm
apart. Measuring the cosine similarity between the per-gold-atom basis vectors:

| D | cos(monomer, dimer) | cos(monomer, trimer) |
|---|---|---|
| **12.9 nm** | **0.9965** | 0.9857 |
| 29.2 nm | 0.9844 | 0.9410 |
| 42.0 nm | 0.9608 | 0.8736 |

At 13 nm the CDA dimer basis vector is **99.65% collinear with the monomer's**. NNLS consequently
assigns it zero weight and returns 100% monomer, leaving the red tail entirely unexplained.
**A speciation fit on a CDA basis is not approximate — it is meaningless.** This is why the exact
T-matrix cache is mandatory, and why any run that silently falls back to CDA must be discarded.

Note the collinearity *worsens* as particles shrink — the 13 nm system is the hardest case, and even
on the exact basis the **dimer:trimer partition will be far less identifiable than the aggregated
fraction**. This should be reported with every speciation number.

**Measured update (2026-07-11, `scripts/fit_ctac_validation.py`):** the exact basis at
12.9 nm / gap 1 nm brings cos(monomer, dimer) down to **0.9593** — identifiable — but the gain is
**gap-conditional**: at a ligand-realistic 3.5 nm gap the exact 13 nm dimer is back to **0.9979**.
Identifiability at ~13 nm therefore rests on the assumed near-contact geometry, and any reported
fraction is conditional on it. At 29–42 nm the exact basis is well-conditioned at any plausible
gap (0.925–0.969 at gap 3.5), which is what makes the CTAC series a usable validation set.

**Implementation.** `treams` runs offline over a grid (species × diameter × gap × wavelength ×
temperature) and is cached (`tmatrix_basis.npz`); the fit interpolates. A live solve inside the
optimiser loop would be intractable.

---

## Layer 5 — Equilibrium

Species populations follow van 't Hoff, K(T) = exp[−(ΔH − TΔS)/RT], with total gold conserved.
**The thermodynamic sign bounds must permit both signs** — the original defaults forced ΔH, ΔS ≤ 0,
which makes thermally-*induced* aggregation impossible by construction. This was a genuine bug.

---

## Layer 6 — The scattering pedestal (λ^−n)

**Chosen: a phenomenological λ^(−n) baseline, with n fitted.**

The measured red wing carries ~16–19% of peak extinction at 790 nm; the *exact* T-matrix dimer and
trimer reach only ~2%. An order of magnitude is unaccounted for. The pedestal term absorbs it.

**This term is a construction of convenience, not a literature model, and it is documented as such.**
Its anchors are rigorous — **n = 4** (Rayleigh, particles ≪ λ) and **n = 0** (geometric limit, ≫ λ) —
and the named analogue is the **Ångström exponent** (Ångström 1929), used exactly this way as a
particle-size proxy in aerosol optics. But:

- The **n → floc-size inversion is non-monotonic**, strongly dependent on an unconstrained packing
  fraction (n ≈ 0.85 maps to anywhere from 160–345 nm), and  was computed via
  Maxwell–Garnett + Mie, which is **mean-field and therefore not plasmon-coupling aware**. That is
  ironic: near-field coupling is the entire reason gold aggregates have red tails.
- **⇒ No floc size should be reported from n_sca.** It is a shape descriptor only.

**Planned replacement:** a real coupling-aware **fractal aggregate** species (Sorensen 2001) —
parameterised by (N_primaries, fractal dimension D_f) rather than an unphysical fill fraction, and
solved with the same exact MSTM machinery. This is tractable (~200 primaries at D_f ≈ 1.8 spans
150–250 nm; a 6,000×6,000 matrix, seconds per wavelength), whereas a *dense* 213 nm floc would need
~1,500 primaries and a 45,000² matrix — infeasible. **Deferred pending DLS**, which measures the floc
size directly and would render the whole exercise unnecessary if no such population exists.

---

## Layer 7 — Data preprocessing (concentration)

Three modes, in order of scientific preference:

1. **`density`** — Kell (1975) water thermal expansion. **The correct, general mode.** No free
   parameters, works on a single branch. **Recommended for any properly sealed cell.**
2. **`evaporation`** — a **salvage mode** for the 2011 dataset, whose cuvette — although
   sealed — still lost solvent from the probed volume (leak past the cap and/or headspace
   condensation).
   Kell expansion × a one-parameter Antoine-vapour-pressure evaporation term, with α calibrated
   *model-free* on the matched-temperature heating-vs-cooling branch offsets. Requires both branches.
3. **`mult_400nm`** — the original anchor. **Deprecated.** Its premise (that A(400 nm) responds only to
   concentration) is violated, and forcing A(400) constant inflates the apparent plasmon-peak change
   ~8.5×.

The genuinely reusable contribution here is the **diagnostic**: comparing branches at matched
temperatures decomposes any discrepancy into a wavelength-**flat** part (concentration drift) and a
wavelength-**structured** part (chemistry). Every temperature-series dataset should run it. It would
have caught the 2011 evaporation at the bench.

**Lab recommendation: a seal alone is not enough — minimise headspace and weigh the cell
before and after the run.** The 2011 cuvette was sealed and still concentrated ~5.6%, which is
indistinguishable from a growing broadband scattering signal and required a two-branch salvage
correction to remove. Weighing turns any residual loss into a one-line correction.

---

## References

**Gold permittivity**
- Johnson, P. B.; Christy, R. W. *Phys. Rev. B* **1972**, 6, 4370. (tabulated n,k — production)
- Reddy, H.; Guler, U.; Kildishev, A. V.; Boltasseva, A.; Shalaev, V. M. *Opt. Mater. Express* **2016**, 6(9), 2776. (**Drude + 2 critical-point ε(λ,T) — the T-dependence**)
- Olmon, R. L. et al. *Phys. Rev. B* **2012**, 86, 235147. (Drude baseline: ħω_p = 8.45 eV, Γ_D = 44 meV)
- Rakić, A. D. et al. *Appl. Opt.* **1998**, 37, 5271. (Brendel–Bormann)
- Etchegoin, P. G. et al. *J. Chem. Phys.* **2006**, 125, 164705. (analytic model)
- Yakubovsky, D. I. et al. *Opt. Express* **2017**, 25, 25574. (thin-film n,k)
- Klinavičius, T. et al. *J. Phys. Chem. C* **2025**, 129, 17616. (which ε to use for colloidal AuNP)
- Zollner, S. et al. *Adv. Opt. Technol.* **2022**. (ellipsometry 10–700 K)

**Surface / size damping**
- Kreibig, U.; Vollmer, M. *Optical Properties of Metal Clusters*, Springer, **1995**.
- **Berciaud, S.; Cognet, L.; Tamarat, P.; Lounis, B. *Nano Lett.* 2005, 5(3), 515–518.** (**A = 0.25**, single-particle PHI)
- Chetoui et al. **2026**. (γ(T) = γ^(ep)(T) + γ_S, with γ_S temperature-independent)

**Cluster optics**
- Beutel, D.; Fernandez-Corbaton, I.; Rockstuhl, C. *Comput. Phys. Commun.* **2024**, 297, 109076. (**treams**)
- Mackowski, D. W.; Mishchenko, M. I. *J. Opt. Soc. Am. A* **1996**, 13, 2266. (MSTM; analytic orientation averaging)
- Xu, Y. *Appl. Opt.* **1995**, 34, 4573. (generalised multiparticle Mie)
- Jain, P. K.; Huang, W.; El-Sayed, M. A. *Nano Lett.* **2007**, 7, 2080. (plasmon ruler)
- Kelesidis, G. A. et al. **2022**. (DDA on gold agglomerates — the validation target for the floc model)

**Aggregate / baseline**
- Sorensen, C. M. *Aerosol Sci. Technol.* **2001**, 35, 648–687. (fractal aggregate light scattering)
- Ångström, A. **1929**. (the λ^−α exponent as a size proxy)

**Classical / comparison tools**
- Haiss, W.; Thanh, N. T. K. et al. *Anal. Chem.* **2007**, 79, 4215. (A_450 → concentration; A_spr/A_450 → size; explicitly notes the red side is aggregation-contaminated)
- Kell, G. S. *J. Chem. Eng. Data* **1975**, 20, 97. (water density vs T)
- Klinavičius 2024; Bilén 2024; Glaubitz 2024. (ML size-from-spectrum — all assume monomers)
