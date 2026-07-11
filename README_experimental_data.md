# C500 experimental UV-Vis T-series — data for `aunp_speciation` fitting

> ⚠️ **STATUS (2026-07-11): ALL THREE HEADLINE VERDICTS BELOW ARE PROVISIONAL —
> pending an evaporation-aware concentration correction (CLAUDE.md #14).**
> Matched-temperature branch comparison shows the sample CONCENTRATED ~5.6% over
> the run (cooling brighter than heating by +5.6% @15 °C, shrinking to +1.9%
> @65 °C; wavelength-independent in % ⇒ concentration, and scaling with the scan
> time gap ⇒ evaporation). Kell corrects thermal expansion only and therefore
> OVER-BRIGHTENS high-T spectra by +2.5–4.2% — a wavelength-flat, T-growing error
> degenerate with the fitted pedestal.
>
> **1. Gold ε(λ,T) IS implemented** (bulk Holstein γ_bulk(T) + T-independent γ_S)
> plus water n(T) — that physics stands. Under Kell normalization the
> thermally-driven speciation signal collapsed (ΔH₂ ≈ 0) and the surviving
> pedestal came out exactly FLAT (n≈0) growing +31% — **both PROVISIONAL**: the
> flat growing pedestal is exactly what Kell's evaporation blindness injects,
> and n_sca moving 0.73→0.00 with the normalization change is a red flag.
> See CLAUDE.md #8/#11/#14.
>
> **2. `normalize="mult_400nm"` is DEPRECATED/BIASED** (premise violated; peak
> change inflated ~8.5×). `normalize="density"` (Kell) is better-motivated but
> **blind to evaporation** (#14). The isosbestic verdict is **PROVISIONAL**:
> under Kell there is no crossing; under a correct multiplicative concentration
> correction crossings DO exist and drift 562→577 nm. Its existence/location is
> normalization-dependent — do not build the two-state argument on it yet.
>
> **The RAW files below remain the correct canonical input either way.**

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
1. **Concentration normalization** (opt-in, `normalize=`). RECOMMENDED:
   `normalize="density"` — pure dilution correction from water thermal
   expansion, `A_norm(λ,T) = A(λ,T) · ρ(T_ref)/ρ(T)` (Kell 1975;
   ρ 999.10 → 974.85 kg/m³ over 15→75 °C). No anchor-wavelength assumption;
   needs a parseable temperature in every column header.
   DEPRECATED: `normalize="mult_400nm"` (anchor at 400 nm) — premise violated
   (A(400) rises +3.18% while expansion predicts −2.43%; the flat rescale
   inflates the apparent peak change ~8.5×, CLAUDE.md #12). Kept only for
   reproducing old results.
2. **THEN clip to the fit range** `wavelength_range=(420, 800)` (clear of the
   350 nm lamp changeover and deep-UV interband edge).
   ORDER MATTERS: normalize first, clip second.

Anchor (`anchor_nm=400`, mult_400nm only) and fit range (`wavelength_range`)
are SEPARATE knobs: the anchor lives in the raw data but outside the fit
window on purpose.

Format: column 1 `wavelength_nm`, then one `ext_<T>C` column per temperature.
Run: `python scripts/fit_real.py aunp_heating_RAW_390-900.csv --normalize density`
(fits with gold ε(T) + water n(T) per temperature; `--fixed-eps` for legacy fits).

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
- Isosbestic crossing ≈ 575 nm (report value) — **status PROVISIONAL /
  normalization-dependent** (CLAUDE.md #12, #14): absent under Kell density
  normalization, present and drifting 562→577 nm under a multiplicative
  concentration correction. Not settled until the evaporation-aware
  correction is in.
- Red-tail/peak ratio (700/523 nm) ≈ 0.18→0.21, growing with T. Post-ε(T)
  decomposition (PROVISIONAL, #14): gold+water T-physics covers the
  peak-region changes; the red-wing growth appeared as a flat pedestal under
  Kell — possibly evaporation-injected in part. See CLAUDE.md #8/#11/#14.
- **Evaporation (measured, #14):** cooling branch brighter than heating at
  matched T (+5.63% @400 / +5.65% @523 / +5.95% @700 / +6.78% @790 nm at
  15 °C; +1.9% at 65 °C) — sample concentrated ~5.6% over the run.
