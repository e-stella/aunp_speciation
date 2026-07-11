# C500 experimental UV-Vis T-series — data for `aunp_speciation` fitting

Source: `final.xls` (Cary 100 UV-Vis, "Isosbestic pt exp", Dragnea Lab, 2011).
Sample: citrate/TEG-functionalized AuNP, TEM diameter **12.9 nm ± 7%**, ~45×
diluted from the optical-force stock into the spectrophotometer's linear range.
Characterized in `C500_Final_Report_Eun_Sohl_Koh.pdf` (Figs 16–17 = this data).

## CANONICAL INPUT — use these (raw, full range, normalization done IN CODE)
- `aunp_heating_RAW_390-900.csv` — 18 temperatures, 15→75 °C (heating branch)
- `aunp_cooling_RAW_390-900.csv` —  6 temperatures, 65→15 °C (cooling branch)

These are the honest instrument data: per-temperature 2-replicate averages from
the `orgnz` sheet, 390–900 nm, **no normalization applied**. 400 nm IS present so
the normalization can be derived in code. Design intent: the repo's data layer
(`io_data.load_series`) does the preprocessing, not a pre-baked CSV.

### Preprocessing the data layer must apply (in this order)
1. **Multiplicative normalization anchored at 400 nm** (opt-in, `normalize=`):
   `A_norm(λ,T) = A(λ,T) · A(400,T_ref) / A(400,T)`, with `T_ref` the first/
   coldest column (15 °C).
   Rationale: at 400 nm gold is ~dielectric (interband), so extinction there
   tracks total cross-section (concentration) and is ~speciation-invariant;
   scaling (not offset) is the Beer–Lambert-correct concentration correction.
   Removes ~2–3 % inter-scan drift (solvent thermal expansion etc.).
2. **THEN clip to the fit range** `wavelength_range=(420, 800)` (inside the
   Etchegoin dielectric's validated window; clear of the 350 nm lamp changeover
   and deep-UV interband edge).
   ORDER MATTERS: normalize first (needs 400 nm in the array), clip second.

Anchor (`anchor_nm=400`) and fit range (`wavelength_range`) are SEPARATE knobs:
the anchor lives in the raw data but outside the fit window on purpose.

Format: column 1 `wavelength_nm`, then one `ext_<T>C` column per temperature.
Run: `python scripts/fit_real.py aunp_heating_RAW_390-900.csv --normalize mult_400nm`
(exact flag name per the io_data API Claude Code builds).

## PRE-NORMALIZED CONVENIENCE COPIES (optional)
- `aunp_heating_series.csv`, `aunp_cooling_series.csv`
  Same data with the 400 nm multiplicative normalization ALREADY applied and
  clipped to 420–800 nm. Fit-ready as-is, but the 400 nm anchor is gone so the
  normalization can't be re-derived or audited from them. Prefer the RAW files +
  in-code normalization for reproducibility; these exist only for a quick run.

## DO NOT USE
- `aunp_*_uncorrected.csv` (earliest files): raw averages, already clipped
  (no 400 nm point), NO normalization. Superseded — missing anchor AND missing
  normalization. Kept only for provenance.

## Physics notes for the fit
- Experimental plasmon peak ≈ **523 nm** (nearly T-stationary, 523.1→523.8 nm).
- Isosbestic crossing ≈ **575 nm** (report value; reproduced from this data
  after the 400 nm normalization).
- Red-tail/peak ratio (700/523 nm) ≈ 0.18→0.21, **growing with T** — the
  aggregation signal. Any model whose dimer/trimer basis lacks a comparable
  red wing (e.g. CDA) cannot fit this. See CLAUDE.md "Known limitations" 7–9.
