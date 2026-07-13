# HANDOFF — recalibrate monomer damping for the seeded-growth series

_From the tem-particle-metrics session, 2026-07-13. Context for a Claude
session working in this (aunp-speciation) repo._

## Task
Recalibrate the monomer surface-damping (γ_S) for the **seeded-growth
citrate/H₂O₂ series** (Oct + Nov 2024), then re-run the speciation fit and see
whether the Oct "aggregation" is real or a line-width artifact.

## Why
The current damping (`S_CAL=0.05, A_CAL=0.6`, frozen 2026-07-12 in
`scripts/calibrate_damping.py`) was fit on the **CTAC** samples — out-of-sample
for this chemistry. Symptom, from a first-cut speciation fit (TEM-pinned size,
monomer+dimer+trimer NNLS, done in ../tem-particle-metrics):

- The Oct UV-Vis spectra are broader than the damped-monomer model **on both
  the blue and red sides of the peak** (roughly symmetric).
- Even Oct **GNP12 (~11.55 nm, physically pure monomer) gets assigned ~13 %
  "aggregation"** — unphysical for a 12 nm sol. That is the tell: the model
  line-width is too *narrow* for this chemistry, so the fit pads width with
  fake aggregation (GNP20 → ~41 %, GNP40 → ~40 %).
- Symmetric broadening argues for **damping/line-broadening, not aggregation**
  (aggregation adds a *red-shifted* tail specifically).

## Plan
1. Use **Oct GNP12 as the pure-monomer reference** (it must be ~100 % monomer):
   fit γ_S (adapt `scripts/calibrate_damping.py`) so the model matches the
   *measured* GNP12 spectrum shape.
2. Re-fit the Oct series (GNP20/40/60) with the recalibrated damping.
   **Diagnostic:**
   - GNP20 "aggregation" collapses toward 0 → it was a **model artifact**
     (whole Oct series just under-damped); the real broadening cause (H₂O₂
     matrix?) is a symmetric line-width effect, not aggregation.
   - GNP20 stays anomalous after GNP12 is fixed → a **real, sample-specific**
     effect worth chasing.

## Pointers (this repo)
- Measured UV-Vis: `experimental/Irina/Scan_10_22*.csv` (Oct),
  `Scan_11_16*.csv` (Nov). Column↔sample map: Oct = GNP12/20/40/60nm;
  Nov "Seed/30/40/60" ↔ TEM 11/23/39/55 nm.
- Forward model: `spectra.monomer_polydisperse`, `spectra.species_basis`,
  `mie.monomer_cross_sections`; existing `scripts/calibrate_damping.py`,
  `scripts/predict_irina_seeded.py`.
- Env: `mstm-env` (Python 3.10).

## From the TEM side (../tem-particle-metrics)
- Per-sample TEM size (mean±sd) is validated to **within ±6 %** of hand
  measurement; use the TEM mean as the fixed size when fitting.
- First-cut speciation + method: `docs/speciation_fit.md`,
  `scripts/fit_speciation.py`. Ground-truth diameters: `outputs/dm3/gt/*.csv`.
- A ~+5 % tier-2 size-oversizing calibration is being applied on the TEM side
  (mask-boundary artifact); once landed, the TEM means shift down slightly.

## Note
An email to the PIs (Bogdan, Irina) already stated this recalibration should
resolve the Oct GNP20 question — so the diagnostic above is the deliverable.

---

# RESULTS (2026-07-13, aunp-speciation session)

**Verdict: the GNP12 13% was the line-width artifact (13→0%), but the Oct
GNP20/40/60 aggregation is REAL — it does not collapse under the
recalibrated monomer (GNP20 41→45%).** The two arms of the diagnostic split:
the route was indeed under-damped, fixing it fully cleans the small clean
samples, and what remains on GNP20/40/60 is a red-asymmetric excess — the
aggregation signature, not symmetric line width.

## Calibration (`scripts/calibrate_damping_seeded.py`, fig20)
- Same methodology as the CTAC calibration; fitted jointly on the four
  verified-clean monomers (Oct GNP12 + Nov 11/23/39) using the hand-picked
  per-particle TEM (`outputs/dm3/gt/`). GNP55/60 and Oct GNP20/40 excluded.
- **Route pair: s = 1.75, A_surf = 1.0** (−2.7 nm offset kept — it
  transfers). Clean-set mean RMS: baseline 7.74% / CTAC pair 7.92% →
  **3.00%**.
- Direction is OPPOSITE the CTAC result (s=0.05): seeded-growth H₂O₂/citrate
  particles are more damped than J&C bulk — defect-rich/multiply-twinned
  overgrowth, exactly the fig17 transfer-failure reading. The s/A split is a
  shallow ridge (γ_eff is what's constrained); GNP12-only rails to s=3.0
  because s is a flat direction at 11.6 nm — the joint fit is the production
  number.
- No Oct-vs-Nov batch split (nov-only picks the same ridge; GNP12 sits near
  its optimum) → route-level property, H₂O₂-matrix-as-batch-effect not
  supported. Peak residuals after −2.7 nm: −1..+2 nm (no medium-index shift
  from the unpurified matrix).

## Speciation re-fit (`scripts/fit_speciation_seeded.py`, fig21)
Same fit as your `fit_speciation.py` (SAM tier-2 sizes, CDA basis, gap scan,
NNLS 450–740 nm); first-cut config reproduced your table exactly, then
re-fit with the route pair:

| sample | agg first-cut | agg recalibrated | fit RMS |
|--------|:---:|:---:|:---:|
| oct_GNP12 | 13% | **0%** | 0.060 → 0.022 |
| oct_GNP20 | 41% | **45%** (41% at the s=3.0 ridge point) | 0.105 → 0.070 |
| oct_GNP40 | 40% | 40% | 0.074 → 0.055 |
| oct_GNP60 | 68% | 68% | 0.067 → 0.060 |
| nov_GNP11 | 6% | 0% | 0.060 → 0.026 |
| nov_GNP23 | 0% | 0% | 0.035 → 0.026 |
| nov_GNP39 | 7% | 7% | 0.029 → 0.012 |
| nov_GNP55 | 100% | 100% | 0.207 → 0.193 |

Guardrails: Nov clean samples stay clean with improved RMS; GNP55/60 stay
aggregated (extra damping did not eat real signal).

## Exact T-matrix re-fit (same day; fig22 — supersedes the caveat below)
Per-sample exact bases (dimer + linear trimer + linear tetramer, calibrated
damping baked in, gap grids {1,2,3.5} nm at GNP20 / {2,3.5} at the larger
three, lmax from a convergence probe):

| sample | agg CDA | agg EXACT (best gap) | gap range | fit RMS CDA→exact |
|--------|:---:|:---:|:---:|:---:|
| oct_GNP20 | 38% | **31%** | 31–38% | 0.066 → 0.038 |
| oct_GNP40 | 38% | **46%** | 41–46% | 0.041 → 0.030 |
| oct_GNP60 | 67% | **81%** | 75–81% | 0.036 → 0.045 |
| nov_GNP55 | 100% | 100% | — | 0.134 → 0.108 |

- **GNP20: ~⅓ of the gold is aggregated, robust across backends and gaps**
  (CDA 38–45%, exact 31–38%) — the quantitative version of the verdict.
- ⚠ The earlier "CDA fractions are upper bounds" reasoning was WRONG at
  large D: exact ≥ CDA for GNP40/60. Exact chain modes place red intensity
  in shaped bands rather than uniformly boosting the tail; quote the
  CDA↔exact spread as the backend systematic.
- GNP55 has no quotable fraction: its smooth 580 nm plateau is beyond any
  single-gap dimer/trimer/tetramer basis (RMS 0.108, structured residual) —
  consistent with the experimenter's own "aggregated, redo" flag; needs
  motif/gap distributions or DLS, not a longer basis.
- Numerical notes for reuse: treams chain solves go unstable (negative
  Cext) at lmax=12 for gap≤1 nm even where the dimer converges — cap chains
  at lmax=10 near contact; gap-1 bonding structure is real and sharp, build
  those columns at 5 nm wavelength sampling.

## Caveat on the absolute Oct fractions [SUPERSEDED by the exact re-fit
above — kept for the record]
CDA under-couples at near-contact, so 40–68% are UPPER bounds (CDA needs
more aggregate per unit red intensity). The T-matrix re-fit (your next-step
2) is still the missing piece for quantitative numbers — but the
artifact-vs-real question needed no exact optics.

## To adopt in ../tem-particle-metrics/scripts/fit_speciation.py
```python
from aunp_speciation import dielectric
dielectric.set_bulk_damping_scale(1.75)   # new module-level knob, default 1.0
# in species_basis(...): A_surf=1.0  (was 0.25)
# blue-shift each basis column 2.7 nm: np.interp(WL, WL - 2.7, col)
```
