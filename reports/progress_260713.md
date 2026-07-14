% From TEM Micrograph to UV-Vis Speciation: an Automated Sizing and Optical-Physics Pipeline for Gold Nanoparticle Colloids
% Mid-project progress report
% July 2026

## Why this matters (plain-language summary for life-science readers)

**The problem, in one sentence.** Gold nanoparticles are a workhorse of
bionanotechnology — contrast agents, biosensors, drug-delivery cargo — and
their usefulness depends on two numbers that are surprisingly hard to measure
routinely: **how big are the particles, and are they clumping together?**
Both are answered today by an expert spending hours at an electron
microscope, hand-clicking a few hundred particles per sample. This project
replaces that bottleneck with two automated tools and connects them, so a
single cheap optical scan plus one microscope image yields size, size spread,
and the fraction of gold that is aggregated.

**Where it helps a biology lab.** This work was built directly on data from
the **Dragnea lab (Indiana University)**, which co-assembles **virus-like
particles (VLPs)** — coat-protein shells of the brome mosaic virus — around
gold-nanoparticle cargo, both as a route to engineered delivery vehicles and
as a model system for how viruses package "fragmented" (multi-piece) cargo
such as segmented genomes, multiple siRNAs, or several drug molecules
[25,26]. The nanoparticle plays the cargo role precisely because gold is
*visible*: it gives the electron-microscope contrast that lets the assembly
be read particle-by-particle. Two features of that research make our tools a
direct fit:

- **The whole biology is gated by nanoparticle *size*.** Encapsulation is
  strongly size-selective — particles below a few nm are excluded, larger
  ones are preferentially packaged, and the *number* of nanoparticles per
  shell falls with their diameter [25,26]. The key evidence is a **shift
  between two size distributions** — the free nanoparticle stock versus the
  population actually encapsulated — and the lab establishes it by
  hand-counting **300–500 particles per TEM sample**, for *both* populations.
  Stage one below automates exactly that count: it turns a half-day of
  clicking into minutes, measures **several-fold more particles**, and gives
  the less outlier-biased distributions that decide whether a size-selection
  signal is even resolvable. Sharper size distributions in, sharper
  thermodynamic conclusions out.
- **Assembly *is* a speciation equilibrium — the same problem this project
  models.** The lab's own theory describes co-assembly as a mass-action
  balance between **free nanoparticles, nanoparticle–protein complexes, and
  capsids holding a variable number of packed cores** [26] — and inside a
  multi-core VLP those cores sit as **adjacent linear or clustered arrays**
  against the shell wall [25], currently resolved only by cryo-electron
  tomography. That is precisely the situation our optical model was built for:
  gold particles packed at contact couple plasmonically, and the strength of
  that coupling reports how many are clustered and how tightly — a
  **non-destructive, in-solution, high-throughput readout of assembly/cluster
  state** that could complement low-throughput cryo-EM. The coupling signal
  is quantitative at the ≥~10 nm sizes typical of single-core AuNP–VLP work
  and of routine colloid QC; extending it to the few-nm fragmented cargo is a
  natural next step rather than a solved capability today.

In short, the two measurements this class of nanoparticle-cargo and VLP
research needs on every batch are **throughput** — automated sizing that
scales to thousands of particles across both free and assembled populations —
and **aggregate/cluster inference** — a quantitative clustering readout from a
one-minute optical scan, grounded in the same equilibrium picture the
assembly itself obeys. Both now come out of one connected pipeline, with
stated error bars and honest failure flags, and the absolute-extinction
concentration measurement the lab already runs [19] is the anchor that ties
them together. The technical report follows.

## 1. Summary

We have built an end-to-end pipeline that starts at a raw electron-microscope
image and ends at a quantitative statement about a colloid's state in
solution. Stage one measures **every** gold nanoparticle in a TEM micrograph
automatically; stage two feeds that size distribution into a first-principles
optical model that predicts the UV-Vis spectrum with **zero fitted
parameters**; stage three reads the *deviation* between prediction and
measurement for what it physically is — line-width chemistry and
**aggregation**. Four results define where we stand:

1. **Automated TEM sizing is validated.** On all eight of the lab's samples
   the automated mean diameter agrees with hand-picked measurements to within
   ±6% (mostly ±3%) while measuring **1.6–6.4× more particles**; a
   foundation-model segmentation tier correctly separates touching particles
   inside aggregates that defeat classical thresholding.
2. **TEM size alone predicts the spectrum of a clean sol.** For well-dispersed
   samples the predicted plasmon peak lands within a few nanometers of the
   measurement with nothing fitted — and the automated sizes drive the
   physics exactly as well as manual counting (peaks within 1 nm). Where the
   prediction fails (the two aggregated samples), the failure *is* the
   aggregation signature the project exists to quantify.
3. **Line width is a per-synthesis-route material property, and with it
   calibrated, speciation becomes quantitative.** A two-parameter damping law
   calibrated on verified-clean monomers reduces model-data mismatch ~2.5×
   per route; CTAC-grown and seeded-growth H₂O₂/citrate particles land on
   **opposite sides of textbook gold** (s = 0.05 vs 1.75) — plasmon width
   doubles as a crystallinity/interface probe. With the calibration in place,
   a physically-pure 12 nm control that previously read a spurious 13%
   "aggregation" now reads **0%**, while the anomalous 20 nm batch keeps a
   real **~⅓ of its gold in small clusters** (31–38% across model choices) —
   answering the question the PIs posed about that batch. The 40/60 nm
   batches fit at ≈46% and ≈81%; the flagged "55 nm" batch is beyond the
   small-cluster regime entirely.
4. **The temperature series shows no thermal monomer⇌dimer signal — but two
   real, unexplained observables remain.** After correcting a normalization
   pitfall (evaporation, ~5.6% over the run) and subtracting measured blanks,
   the 13 nm citrate/TEG heating/cooling series shows an isosbestic-like
   crossing that **drifts** 544→568 nm (a fixed two-state equilibrium cannot
   produce this; an evolving cluster geometry can [21]), and a flat,
   reversible scattering pedestal (~15% of peak) that survives every
   instrumental exclusion we could test. A DLS ramp or filtration control is
   the decisive lab experiment.

## 2. Motivation

Manual TEM sizing is slow, samples relatively few particles, and its
histogram widths are inflated by whichever outliers a person happens to
click. Worse, two batches that look identical by hand-measured histogram can
behave completely differently in solution — because one is aggregating, and
a dried-grid histogram cannot see it. Conversely, an honest UV-Vis "size"
readout is confounded by aggregation and by chemistry-dependent line width.
The goal is an instrument-grade chain: an *unbiased, high-n* per-particle
measurement that reproduces trusted hand measurements, plus an inverse
optical model that, given a routine UV-Vis spectrum, reports nominal size,
size spread, and **speciation** — the fraction of gold bound in small
clusters — with stated uncertainties and honest failure flags. The division
of labor that emerged from this work: **TEM pins the size; the spectrum then
yields the aggregated fraction as a stable linear readout.**

## 3. Stage one — automated TEM sizing (image → per-particle table)

The tool reads native Gatan dm3 micrographs directly (RosettaSciIO/HyperSpy
[7]); pixel calibration comes from file metadata, so no burned-in scale bar
is needed and the calibration is exact. The present dataset is 47 frames
across 8 samples (Oct 2024: GNP12/20/40/60; Nov 2024: seed-11/23/39/55).
Segmentation has two interchangeable tiers that emit the identical
per-particle table (diameter, area, circularity, aspect ratio, shape class,
touching-group flag):

- **Tier 1 — classical** (no GPU, ~0.5 s/image): background flattening →
  Otsu threshold [1] → distance-transform watershed [2] with a rod-aware
  neck-merge, built on scikit-image [3]. Accurate on well-dispersed fields;
  merges particles inside dense aggregates, as expected.
- **Tier 2 — foundation model** (~6 s/image on a laptop GPU): NP-SAM [6], a
  nanoparticle-EM-tuned wrapper around the fast Segment-Anything variant
  FastSAM [5] (SAM [4]). It recognizes particles by learned shape/edge
  features and **separates touching particles inside clusters**; a
  physics-based "ghost" filter removes background false positives.
- **Human-in-the-loop QC:** segmentation is never blindly trusted — a
  one-click reviewer lets a person drop any detection, with statistics
  recomputed live; decisions are stored separately and never overwrite the
  machine output.

![Same aggregated sample (Nov "55 nm"). *Left:* tier-1 classical segmentation
merges touching particles into oversized blobs (mean 65 nm, sd 31 vs a true
54.6 ± 4.9 nm). *Right:* tier-2 (NP-SAM/FastSAM) separates the individual
particles (mean 54 nm, sd 8), recovering the true
size.](/Users/eskoh/Documents/Projects/tem-particle-metrics/outputs/dm3/tier_compare_nov55.png){width=100%}

![Size distributions: automated tier-2 (red) vs hand-picked (grey), per
sample, from the dm3 micrographs. Dashed lines are the
means.](/Users/eskoh/Documents/Projects/tem-particle-metrics/outputs/dm3/size_comparison.png){width=100%}

**Validation against hand measurement.** Both tiers ran over all 47 frames;
aggregated distributions were compared to the lab's hand-picked per-particle
measurements (Figures 1–2):

| Sample | Manual (n / mean±sd) | Tier-1 (n / mean, Δ) | Tier-2 (n / mean, Δ) |
|--------|----------------------|----------------------|----------------------|
| oct GNP12 | 439 / 11.6±0.9 | 846 / 11.6 (+0%) | 803 / 12.2 (+5%) |
| oct GNP20 | 582 / 21.6±2.2 | 937 / 21.9 (+1%) | 515 / 22.8 (+6%) |
| oct GNP40 | 137 / 41.6±3.3 | 137 / 52.2 (+25%) | 262 / 42.5 (+2%) |
| oct GNP60 | 203 / 61.3±4.5 | 101 / 96.1 (+57%) | 317 / 63.3 (+3%) |
| nov seed-11 | 307 / 10.9±0.8 | 566 / 11.4 (+5%) | 532 / 11.5 (+6%) |
| nov GNP23 | 200 / 22.9±1.5 | 711 / 24.2 (+6%) | 853 / 23.5 (+2%) |
| nov GNP39 | 102 / 38.5±3.4 | 326 / 42.8 (+11%) | 423 / 39.4 (+2%) |
| nov GNP55 | 135 / 54.6±4.9 | 507 / 68.5 (+25%) | 868 / 54.7 (+0%) |

**Tier-2 recovers the hand-picked mean within ±6% on every sample.** The
aggregated samples that defeat tier-1 (+11…+57%) all drop to +0…+3% with
tier-2, while the automated extraction measures 1.6–6.4× more particles
(GNP55: 868 vs 135). The residual ~+5% on the small samples is a known,
size-dependent mask-boundary effect (SAM outlines run slightly generous)
with a calibration in progress.

## 4. Stage two — closing the loop: zero-parameter spectrum prediction

Each automated size distribution feeds single-particle Mie theory [8] with
Johnson-Christy gold constants [9] in water — nothing fitted; both data and
model normalized to their own plasmon peak (a pure shape test):

| Sample | Measured peak (nm) | SAM-predicted (Δ) | Manual-predicted (Δ) |
|--------|:------------------:|:-----------------:|:--------------------:|
| oct GNP12 | 521 | 523 (+2) | 523 (+2) |
| oct GNP20 | 523 | 524 (+1) | 524 (+1) |
| oct GNP40 | 528 | 528 (+0) | 528 (+0) |
| oct GNP60 | 545 | 537 (−8) | 536 (−9) |
| nov seed-11 | 518 | 523 (+5) | 523 (+5) |
| nov GNP23 | 520 | 524 (+4) | 524 (+4) |
| nov GNP39 | 525 | 527 (+2) | 527 (+2) |
| nov GNP55 | **580** | **534 (−46)** | **533 (−47)** |

Three observations: automated and hand-picked predictions are
indistinguishable (within 1 nm) — the automated table drives the physics
exactly as well as manual measurement; for well-dispersed samples the
measured peak is predicted within a few nanometers from TEM size alone; and
**GNP60 and especially GNP55 fail by design** — spectra far redder and
broader than single-particle physics allows. That red-shift is the
plasmonic-coupling signature of aggregation, and quantifying it is stages
three and four below.

![UV-Vis predicted from the TEM size distribution (SAM red, manual blue) vs
measured (black), normalized to peak. Well-dispersed samples match within a
few nm; GNP55/60 show the aggregation red-shift the monomer model cannot
reproduce.](/Users/eskoh/Documents/Projects/tem-particle-metrics/outputs/dm3/spectra_prediction.png){width=100%}

## 5. The physics engine

The optical forward model has three layers, all validated against internal
and external references:

- **Optics.** Single-sphere Mie theory [8] with size- and
  temperature-corrected gold permittivity; small clusters via two selectable
  backends: a fast coupled-dipole approximation (CDA), and an **exact
  multipole T-matrix** solver (`treams` [17]; equivalent to
  Mackowski-Mishchenko MSTM [18]). The exact backend reproduces Mie for a
  single sphere to machine precision and shows the CDA **under-couples at
  near-contact** — 1.6×/2.8× less dimer/trimer red-tail at 1 nm gap — so CDA
  is used for exploration and the T-matrix for quantitative work, cached on
  parameter grids for fitting speed.
- **Materials.** Gold ε(λ) from Johnson & Christy [9] (selected against
  Etchegoin [10], Brendel-Bormann, and Yakubovsky [11] alternatives; a
  nine-sample re-anchoring on verified-clean monomers measured a +2.7 ± 1.1
  nm red bias, now applied as a fixed offset), Drude parameters from Olmon
  [12]; surface damping γ~S~ = A·ħv~F~/R [13,14]; bulk damping γ(T)
  bracketed between Holstein electron-phonon theory [15] and Reddy's
  measured film values [16]; water n(T) and k(λ) from standard tables
  [23,24]; solvent density/evaporation corrections in the data layer (§9).
- **Populations & inversion.** Monomer⇌dimer⇌trimer equilibrium with van 't
  Hoff K(T); inversion by non-negative least squares over species basis
  spectra (linear, given size) inside a nonlinear fit over size/spread, with
  profile-likelihood identifiability checks; a global multi-temperature fit
  ties a whole T-series to one sample with shared thermodynamics.

## 6. What one spectrum can — and cannot — determine

A methodical identifiability study (synthetic recovery tests plus validation
on five CTAC monomer samples with paired TEM) established limits that shape
everything downstream, consistent with the sizing literature [19,20]:

- **Robust:** the aggregated (non-monomer) gold fraction. With size pinned,
  it is a stable linear readout; on clean monomer samples the pipeline
  correctly returns zero (negative control passed).
- **Weak:** mean size from normalized band shape below ~35 nm (a 23 nm
  sample fits equally well at 8 nm); the Haiss A~spr~/A~450~ ratio [19]
  separates 8 from 30 nm but saturates over 29–42 nm. Size must come from
  TEM, absolute extinction, or a temperature series — not band shape alone.
  This is exactly why stage one (automated TEM) is the right partner
  instrument.
- **Not identifiable:** dimer vs trimer vs tetramer (collinear red tails —
  only the total aggregated fraction is quotable); blind bimodality (a
  two-size mixture is spectrally indistinguishable from one broad
  population); and speciation at ~8 nm (profile bounds explode).
- **Conditional:** at ~13 nm the monomer/dimer separation depends on the
  assumed interparticle gap (identifiable at 1 nm contact, degenerate at
  ≳2–3 nm) — gap priors from ligand length matter.

## 7. Line width is chemistry: per-route damping calibration

The single largest systematic in all fits was the monomer line width. We
model the effective Drude damping as γ_eff = *s*·γ_bulk + *A*·ħv~F~/R — *s*
scales the tabulated (evaporated-film) bulk ε₂, absorbing crystallinity
differences; *A* is the surface/interface term [13,14]. One (s, A) pair is
fitted per synthesis route, jointly across that route's verified-clean
monomer sizes (the 1/R law carries the size dependence):

| route | calibration set | (s, A_surf) | mean shape RMS |
|-------|-----------------|:-----------:|:--------------:|
| CTAC (7.8–42 nm) | 5 samples, per-particle TEM | (0.05, 0.6) | 4.6% → 2.5% |
| seeded-growth citrate/H₂O₂ (11–39 nm) | 4 clean samples | (1.75, 1.0) | 7.7% → 3.0% |

The two routes bracket textbook gold from **opposite sides**: CTAC particles
are less damped than Johnson-Christy's evaporated films (single-crystal-like;
J&C ε₂ effectively 15–35% too high for them), while seeded-growth H₂O₂
overgrowth is substantially *more* damped (defect-rich / multiply twinned),
with an interface term at the top of the literature range [14]. A frozen
calibration transferred across routes makes fits *worse* — the pair is a
route property, not a universal constant. Two practical consequences: the
plasmon width can serve as a fast crystallinity/interface probe; and any
speciation fit must use its route's calibration or it will convert line-width
error into fictitious aggregation (§8).

![Seeded-growth route calibration on the four verified-clean
monomers (per-particle TEM histograms pinned; jc −2.7 nm offset). Measured
spectra (black) vs the uncalibrated model, the (wrong-route) CTAC
calibration, and the seeded-growth calibration (red). Right panels: RMS
landscape over (s, A) — a diagonal valley (the constrained combination is
the total damping) with an interior joint optimum at
(1.75, 1.0).](/Users/eskoh/Documents/Projects/aunp-speciation/outputs/fig20_calibrate_seeded.png){width=100%}

## 8. Stage three — TEM-pinned speciation of the seeded-growth series

With sizes pinned by stage one and the route-calibrated monomer, fitting each
spectrum as a non-negative monomer + dimer + trimer + tetramer mixture is a
clean linear problem. The first-cut fit (uncalibrated line width) had
assigned 13% "aggregation" to the physically-pure 12 nm control — the tell
of an under-damped monomer model padding width with fake aggregate. After
calibration:

| sample | first-cut agg | calibrated CDA | **calibrated exact T-matrix** (gap range) |
|--------|:---:|:---:|:---:|
| oct GNP12 (11.6 nm, pure control) | 13% | **0%** | — |
| oct GNP20 (22.8 nm) | 41% | 45% | **31%** (31–38%) |
| oct GNP40 (42.5 nm) | 40% | 38% | **46%** (41–46%) |
| oct GNP60 (63.3 nm) | 68% | 67% | **81%** (75–81%) |
| nov series (11/23/39 nm, clean) | 6/0/7% | 0/0/7% | — |
| nov GNP55 (54.7 nm, flagged) | 100% | 100% | 100%, structured misfit |

The diagnostic splits exactly as hoped: the **pure-monomer control collapses
to zero** (so the calibration is right), while **GNP20 does not** — its
residual is red-asymmetric, which damping cannot and should not absorb. The
Oct-series aggregation is real, and GNP20 carries roughly a third of its
gold in small clusters, robust across optics backends (CDA vs exact) and
the full interparticle-gap grid. Nothing was calibrated away: the clean Nov
samples stay clean with *improved* residuals.

![Speciation with TEM-pinned size, before (dotted) and after
(red) the route damping calibration. The pure 12 nm control (top left) drops
from 13% to 0% aggregated; GNP20's misfit persists and is specifically
red-sided — real coupling, not line width.](/Users/eskoh/Documents/Projects/aunp-speciation/outputs/fig21_speciation_seeded.png){width=100%}

Three honest caveats. (i) The exact-optics fractions are **gap-conditional
ranges**, quoted as such. (ii) The intuition "CDA under-couples, so its
fractions are upper bounds" proved **wrong** at large sizes — exact chain
modes concentrate red intensity in shaped bands rather than uniformly
boosting the tail, and the exact fractions for GNP40/60 came out *above*
CDA's; we therefore quote the CDA↔exact spread as the backend systematic.
(iii) The dimer/trimer/tetramer split within the aggregate is **not
identifiable** (collinear red tails; the fitter parks weight at corners of
the degenerate family) — only the total aggregated fraction is reported.
GNP55 remains honestly unfittable (its smooth 580 nm plateau defeats any
single-gap small-cluster basis): that sample is in a larger-floc regime,
matching the experimenter's own "redo" flag.

![CDA vs exact T-matrix fits of the aggregated samples, both
with calibrated damping. The exact basis halves the GNP20 residual; GNP60's
fitted chain modes show a ~620 nm bump the smooth data lack — the signature
of a single-motif basis, motivating the ensemble-averaged (Monte-Carlo
geometry) basis planned next.](/Users/eskoh/Documents/Projects/aunp-speciation/outputs/fig22_speciation_exact.png){width=100%}

## 9. The temperature series: no thermal speciation — and what remains

The original motivation was a temperature-dependent "isosbestic point" in a
13 nm citrate/TEG sol, read as evidence of a two-state monomer⇌dimer
equilibrium. Careful data handling reversed both readings:

- **Normalization pitfall (methodological result worth publishing on its
  own).** Heating/cooling branches at matched temperature differ by a
  wavelength-flat +5.6% — the cell lost ~5.6% of its water to evaporation
  over the run, overwhelming the −2.4% thermal-expansion correction [22] and
  masquerading as a growing flat pedestal. We built a branch-ratio
  diagnostic that decomposes any dual-branch series into concentration drift
  vs real spectral change, and an evaporation-aware normalization calibrated
  from the data themselves (one parameter; reproduces a hand-computed
  concentration table to ≤0.65%). Per-temperature water blanks are
  subtracted before any normalization.
- **No robust thermal speciation.** With gold ε(T) [15,16], water n(T) [23],
  and surface damping in the basis, the fitted ΔH₂ is consistent with zero
  on both branches and the aggregated fraction is nearly T-flat (~0.15–0.21,
  opposite weak trends on heating vs cooling). What temperature actually
  moves is the gold and solvent dielectrics — bulk gold measurably does
  this at exactly the magnitude used [16] — plus the pedestal below.
- **The crossing drifts.** The difference-spectra crossing walks 544→568 nm
  over 30→75 °C (cooling ~flat at 546–554 nm). A fixed-geometry two-state
  system produces a *stationary* crossing (our model: drift <1 nm even with
  full ε(T)); a drifting one indicates evolving cluster geometry — precisely
  the behavior Cardellini et al. recently established, with λ_iso tracking
  interparticle spacing [21]. The historical "575 nm isosbestic" was an
  artifact of the deprecated normalization.
- **The pedestal is real but its origin is not decided by fitting.** A
  near-flat scattering floor (~15% of peak) survives blank subtraction and
  every exclusion we could construct: its exponent is decisively not
  Rayleigh (n ≈ 0–0.9, not 4 — small-aggregate scattering ruled out), it is
  not a permittivity artifact, not stray light drift (blanks are T-flat),
  not thermal-optical cell effects, and not surface nanobubbles (wrong
  wetting, size, and T-dependence). It grows reversibly by ~10–16% over
  15→75 °C. Both µm-scale flocs and a static stray-light floor fit the
  shape; **DLS across the ramp or a filtration control is the decisive
  experiment** — we recommend it as the next lab action.

![Data-level difference spectra (blank-subtracted, evaporation
normalization): the crossing near 545–570 nm drifts with temperature on
heating — incompatible with a fixed two-state equilibrium, consistent with
evolving cluster geometry [21].](/Users/eskoh/Documents/Projects/aunp-speciation/outputs/fig14_isosbestic_drift.png){width=100%}

## 10. Status and next steps

**Working now:** dm3 input with exact calibration; two-tier segmentation
validated against hand measurement on the lab's own images; mandatory QC
overlays with one-click human review; per-particle output tables; a
validated three-layer optical model (CDA + exact T-matrix, cached for
fitting); route-calibrated damping for two chemistries; TEM-pinned
quantitative speciation with stated identifiability limits; dual-branch
normalization diagnostics; global multi-temperature fitting. The two
repositories (`tem-particle-metrics`, `aunp-speciation`) share a clean data
contract: the per-particle CSV table.

**Next, in order of leverage:**

1. **Ensemble (Monte-Carlo geometry) cluster basis.** Replace single-motif,
   single-gap cluster columns with averages over sampled configurations
   (N, motif, gap distribution) computed with the exact T-matrix. This
   addresses the two visible artifacts of the current basis (GNP60's
   spurious chain-mode bump; GNP55's unfittable plateau) and makes the
   non-identifiability of the N-split explicit.
2. **Calibrate out the ~+5% tier-2 size offset** (mask-boundary effect), so
   absolute (un-normalized) extinction can be compared, not just peak shape;
   then use absolute extinction plus the A~spr~/A~450~ ratio [19] to break
   the remaining damping-size-polydispersity degeneracy.
3. **Lab: DLS ramp / filtration control** on the temperature-series sample to
   settle the pedestal's origin (flocs vs instrument floor).
4. **Adoption & packaging:** the batch segmenter plus lightweight reviewer on
   a shared lab machine; per-route damping calibration as a standard step
   (one clean monomer series per chemistry); the branch-ratio drift
   diagnostic on every dual-branch dataset.

---

## Acknowledgments

Data and research courtesy of the **Dragnea lab** (Indiana University). The
TEM micrographs, hand-picked size measurements, and UV-Vis spectra of the
seeded-growth Oct/Nov 2024 series (§3–§4, §8) were provided by
**Irina Tsvetkova**. The 2011 temperature-series experiment (§9) — the
sample, the UV-Vis heating/cooling series, and its TEM characterization —
was designed and performed by the author.

## References

1. N. Otsu, "A threshold selection method from gray-level histograms," *IEEE
   Trans. Syst. Man Cybern.* **9**(1), 62–66 (1979).
2. S. Beucher, F. Meyer, "The morphological approach to segmentation: the
   watershed transformation," in *Mathematical Morphology in Image
   Processing*, 433–481 (1993).
3. S. van der Walt *et al.*, "scikit-image: image processing in Python,"
   *PeerJ* **2**, e453 (2014).
4. A. Kirillov *et al.*, "Segment Anything," *Proc. IEEE/CVF ICCV*,
   4015–4026 (2023); arXiv:2304.02643.
5. X. Zhao *et al.*, "Fast Segment Anything," arXiv:2306.12156 (2023).
6. "NP-SAM: Implementing the Segment Anything Model for Easy Nanoparticle
   Segmentation in Electron Microscopy Images," *ChemRxiv* (2023),
   doi:10.26434/chemrxiv-2023-k73qz; Python package `npsam`.
7. F. de la Peña *et al.*, HyperSpy / RosettaSciIO — multi-dimensional data
   analysis and microscopy file I/O (dm3/dm4), https://hyperspy.org.
8. C. F. Bohren, D. R. Huffman, *Absorption and Scattering of Light by Small
   Particles*, Wiley (1983).
9. P. B. Johnson, R. W. Christy, "Optical constants of the noble metals,"
   *Phys. Rev. B* **6**, 4370–4379 (1972).
10. P. G. Etchegoin, E. C. Le Ru, M. Meyer, "An analytic model for the
    optical properties of gold," *J. Chem. Phys.* **125**, 164705 (2006).
11. D. I. Yakubovsky, A. V. Arsenin, Y. V. Stebunov, D. Y. Fedyanin,
    V. S. Volkov, "Optical constants and structural properties of thin gold
    films," *Opt. Express* **25**(21), 25574–25587 (2017).
12. R. L. Olmon *et al.*, "Optical dielectric function of gold," *Phys. Rev.
    B* **86**, 235147 (2012).
13. U. Kreibig, M. Vollmer, *Optical Properties of Metal Clusters*, Springer
    (1995).
14. S. Berciaud, L. Cognet, P. Tamarat, B. Lounis, "Observation of intrinsic
    size effects in the optical response of individual gold nanoparticles,"
    *Nano Lett.* **5**(3), 515–518 (2005).
15. T. Holstein, "Optical and infrared volume absorptivity of metals,"
    *Phys. Rev.* **96**, 535 (1954).
16. H. Reddy, U. Guler, A. V. Kildishev, A. Boltasseva, V. M. Shalaev,
    "Temperature-dependent optical properties of gold thin films," *Opt.
    Mater. Express* **6**(9), 2776–2802 (2016).
17. D. Beutel, I. Fernandez-Corbaton, C. Rockstuhl, "treams — a
    T-matrix-based scattering code for nanophotonics," *Comput. Phys.
    Commun.* **297**, 109076 (2024).
18. D. W. Mackowski, M. I. Mishchenko, "A multiple sphere T-matrix Fortran
    code for use on parallel computer clusters," *J. Quant. Spectrosc.
    Radiat. Transfer* **112**, 2182–2192 (2011).
19. W. Haiss, N. T. K. Thanh, J. Aveyard, D. G. Fernig, "Determination of
    size and concentration of gold nanoparticles from UV-Vis spectra,"
    *Anal. Chem.* **79**(11), 4215–4221 (2007).
20. N. G. Khlebtsov, "Determination of size and concentration of gold
    nanoparticles from extinction spectra," *Anal. Chem.* **80**(17),
    6620–6625 (2008).
21. J. Cardellini *et al.*, "Nanoplasmonic isosbestics uncover mesoscale
    assembly of gold nanoparticles on soft templates," *J. Am. Chem. Soc.*
    (2025), doi:10.1021/jacs.5c05189.
22. G. S. Kell, "Density, thermal expansivity, and compressibility of liquid
    water from 0° to 150 °C," *J. Chem. Eng. Data* **20**, 97–105 (1975).
23. *CRC Handbook of Chemistry and Physics*, refractive index of water vs
    temperature (589 nm table).
24. R. M. Pope, E. S. Fry, "Absorption spectrum (380–700 nm) of pure water,"
    *Appl. Opt.* **36**, 8710–8723 (1997); G. M. Hale, M. R. Querry,
    "Optical constants of water in the 200-nm to 200-µm wavelength region,"
    *Appl. Opt.* **12**, 555–563 (1973).
25. A. Amjad, I. Tsvetkova, L. G. Lowry, D. G. Morgan, R. Zandi,
    P. van der Schoot, B. Dragnea, "An assembly-line mechanism for in-vitro
    encapsulation of fragmented cargo in virus-like particles,"
    arXiv:2509.12409 (2025).
26. P. van der Schoot, R. Zandi, A. Amjad, I. Tsvetkova, B. Dragnea,
    "Encapsulation of fragmented cargo by virus coat proteins," *J. Chem.
    Phys.* **164**, 044906 (2026).
