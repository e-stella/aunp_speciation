"""Fit a real UV-Vis file (or the bundled example).

Usage:
    python scripts/fit_real.py [path/to/spectra.csv] [--range MIN_NM MAX_NM]
                               [--normalize {none,density,evaporation,mult_400nm}]
                               [--anchor NM] [--fixed-eps]
                               [--companion CSV] [--branch heating|cooling]
                               [--blanks CSV] [--companion-blanks CSV]
                               [--scan-order CSV]

- Auto-detects a temperature series (>=2 T columns) -> global fit; otherwise a
  single-spectrum fit. The global fit uses gold eps(T) + water n(T) per
  temperature (limitations #11/#13) unless --fixed-eps is passed.
- --range restricts the fitted window (default 420 800 nm); pass "--range 0 1e9"
  to keep the full spectrum. See "Known limitations" #6 in CLAUDE.md.
- --normalize density (for runs with VERIFIED-stable concentration — weigh the
  cell before/after; a seal alone is not proof) corrects concentration for
  water thermal expansion (Kell 1975) — no anchor-wavelength assumption.
  --normalize evaporation (SALVAGE for runs that lost solvent, limitation #14;
  the sealed 2011 C500 cell still lost ~5.6%)
  additionally fits a 1-parameter evaporation model from the matched-T
  heating-vs-cooling offsets; REQUIRES --companion (the other branch) and
  --branch; use --blanks/--companion-blanks (export via
  scripts/export_blanks.py) and optionally --scan-order for true chronology.
  mult_400nm is DEPRECATED/BIASED (limitation #12; inflates the apparent peak
  change ~8.5x); kept only to reproduce old results.
  C500 example:
    python scripts/fit_real.py experimental/ESK_2011/aunp_heating_RAW_390-900.csv \
      --normalize evaporation --branch heating \
      --companion experimental/ESK_2011/aunp_cooling_RAW_390-900.csv \
      --blanks experimental/ESK_2011/blanks_heating_390-900.csv \
      --companion-blanks experimental/ESK_2011/blanks_cooling_390-900.csv \
      --scan-order experimental/ESK_2011/scan_order.csv
- Uses the precomputed EXACT T-matrix basis (outputs/tmatrix_basis.npz) if it
  exists, else falls back to the fast CDA optics. (Build the cache once with
  mstm-env/bin/python scripts/build_tmatrix_basis.py — it now carries the
  eps(T) axis and bakes in gamma_S; fit_real refuses stale caches.)
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from aunp_speciation import dielectric
from aunp_speciation.io_data import load_series
from aunp_speciation.fitting import fit_spectrum
from aunp_speciation.fit_global import fit_temperature_series

# Real-data driver -> measured optical constants. Tabulated J&C puts the 12.9 nm
# monomer at ~523 nm in water (matches the C500 series); the analytic Etchegoin
# and BB models are both ~4-6 nm red-biased here. NB: a T-matrix basis cache
# embeds whatever dielectric it was BUILT with — rebuild it under 'jc'
# (scripts/build_tmatrix_basis.py) or the cache silently overrides this choice.
dielectric.use_gold_model("jc")

HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "..", "outputs")
CACHE = os.path.join(OUT, "tmatrix_basis.npz")
SPECIES = ("monomer", "dimer", "trimer_linear")


def parse_args(argv):
    """Return (data_path, wavelength_range, normalize, anchor_nm, eps_T).
    Flags (`--range MIN MAX`, `--normalize MODE`, `--anchor NM`, `--fixed-eps`)
    may appear anywhere after the program name."""
    args = list(argv[1:])
    wl_range = (420.0, 800.0)
    normalize = None
    anchor = 400.0
    eps_T = True
    if "--range" in args:
        i = args.index("--range")
        wl_range = (float(args[i + 1]), float(args[i + 2]))
        del args[i:i + 3]
    if "--normalize" in args:
        i = args.index("--normalize")
        mode = args[i + 1]
        if mode not in ("none", "density", "evaporation", "mult_400nm"):
            sys.exit("--normalize must be 'none', 'density', 'evaporation' or "
                     f"'mult_400nm', got {mode!r}")
        normalize = None if mode == "none" else mode
        del args[i:i + 2]
    if "--anchor" in args:
        i = args.index("--anchor")
        anchor = float(args[i + 1])
        del args[i:i + 2]
    if "--fixed-eps" in args:
        eps_T = False
        args.remove("--fixed-eps")
    extras = {}
    for flag, key in (("--companion", "companion"), ("--blanks", "blanks"),
                      ("--companion-blanks", "companion_blanks"),
                      ("--scan-order", "scan_order"), ("--branch", "branch")):
        if flag in args:
            i = args.index(flag)
            extras[key] = args[i + 1]
            del args[i:i + 2]
    data = args[0] if args else os.path.join(HERE, "..", "data",
                                             "example_series.csv")
    return data, wl_range, normalize, anchor, eps_T, extras


DATA, WL_RANGE, NORMALIZE, ANCHOR, EPS_T, EXTRAS = parse_args(sys.argv)


def choose_backend():
    """Return (backend, size_correction) — the cache's baked-in flags win."""
    if os.path.exists(CACHE):
        from aunp_speciation.basis_cache import load_cache
        cache = load_cache(CACHE)
        print(f"using EXACT cached T-matrix optics: {os.path.basename(CACHE)} "
              f"(gold_model={cache.gold_model}, "
              f"wl {cache.wl.min():.0f}-{cache.wl.max():.0f} nm, "
              f"D {cache.D.min():.0f}-{cache.D.max():.0f} nm, "
              f"T axis={list(cache.temps) if cache.temps is not None else 'NONE'}, "
              f"size_corr={cache.size_correction})")
        if cache.gold_model != dielectric.current_gold_model():
            sys.exit(f"cache dielectric ({cache.gold_model}) != active "
                     f"({dielectric.current_gold_model()}); rebuild the cache "
                     "with scripts/build_tmatrix_basis.py or delete it.")
        if EPS_T and cache.temps is None:
            sys.exit("eps(T) fitting needs a cache with a temperature axis — "
                     "rebuild with scripts/build_tmatrix_basis.py, or pass "
                     "--fixed-eps to reproduce legacy fixed-eps fits.")
        if EPS_T and cache.eps_t_model != dielectric.EPS_T_MODEL:
            sys.exit(f"cache eps(T) model ({cache.eps_t_model}) != code "
                     f"({dielectric.EPS_T_MODEL}); rebuild the cache.")
        return cache.species_fn, cache.size_correction
    print("cache not found -> using fast CDA optics (lower-bound coupling; "
          "gamma_S ON)")
    return "cda", True


def main():
    wl, spectra, temps_K = load_series(DATA, wavelength_range=WL_RANGE,
                                       normalize=NORMALIZE, anchor_nm=ANCHOR,
                                       **EXTRAS)
    if NORMALIZE == "mult_400nm":
        print(f"applied {NORMALIZE} normalization (anchor {ANCHOR:g} nm) "
              "-- DEPRECATED/BIASED, see limitation #12; prefer 'density'")
    elif NORMALIZE:
        print(f"applied {NORMALIZE} normalization")
    backend, size_corr = choose_backend()
    print(f"loaded {os.path.basename(DATA)}: {spectra.shape[0]} spectrum(a), "
          f"{len(wl)} wavelengths, {wl.min():.0f}-{wl.max():.0f} nm")

    if temps_K is not None and spectra.shape[0] >= 2:
        print(f"temperature series: {temps_K - 273.15} C -> GLOBAL fit "
              f"(eps(T) {'ON' if EPS_T else 'OFF'}, gamma_S "
              f"{'ON' if size_corr else 'OFF'})")
        r = fit_temperature_series(temps_K, wl, spectra, species=SPECIES,
                                   backend=backend, n_sizes=5, fit_stride=2,
                                   C_tot=0.01, max_nfev=120,
                                   eps_temperature=EPS_T,
                                   size_correction=size_corr)
        agg = r.gold_fractions["dimer"] + r.gold_fractions["trimer_linear"]
        print(f"\n  diameter   = {r.diameter_nm:.2f} nm (± {r.param_sd['diameter']:.1f})")
        print(f"  poly       = {r.pct_poly:.2f} %  (± {r.param_sd['pct_poly']:.0f})")
        print(f"  dH2        = {r.dH2:.1f} kJ/mol (± {r.param_sd['dH2']:.1f})")
        print(f"  dS2        = {r.dS2:.3f} kJ/mol/K (± {r.param_sd['dS2']:.3f})")
        print(f"  n_sca      = {r.n_sca:.2f} (± {r.param_sd['n_sca']:.2f})   "
              "[~4 Rayleigh large-aggregate scattering; ~0 flat offset]")
        print("  gold fractions + scattering pedestal vs T:")
        print("      T      mono   dimer  trimer   agg   A_sca@550")
        for i, Tc in enumerate(temps_K - 273.15):
            print(f"     {Tc:5.0f} C {r.gold_fractions['monomer'][i]:7.3f}"
                  f"{r.gold_fractions['dimer'][i]:8.3f}"
                  f"{r.gold_fractions['trimer_linear'][i]:8.3f}"
                  f"{agg[i]:7.3f}   {r.a_sca[i]:.4f}")
        print(f"  residual RMS = {r.residual_rms:.4f}")
        # show only a few temperatures: overlaying all of them makes data and
        # fit indistinguishable. The FIT still uses every temperature.
        TC = temps_K - 273.15
        want = np.linspace(TC.min(), TC.max(), 4)
        show = sorted({int(np.argmin(np.abs(TC - t))) for t in want})
        cmap = plt.cm.coolwarm(np.linspace(0, 1, len(temps_K)))
        norm = spectra.max()
        fig, (ax, axr) = plt.subplots(
            2, 1, figsize=(7.5, 5.8), sharex=True,
            gridspec_kw={"height_ratios": [3.2, 1]})
        for i in show:
            ax.plot(wl[::4], spectra[i, ::4] / norm, "o", ms=3.5,
                    mfc="none", mec=cmap[i], mew=0.9)
            ax.plot(r.wavelength, r.models[i], "-", color=cmap[i], lw=1.7,
                    label=f"{TC[i]:.0f} °C")
            resid = np.interp(r.wavelength, wl, spectra[i] / norm) - r.models[i]
            axr.plot(r.wavelength, resid, "-", color=cmap[i], lw=0.9)
        from matplotlib.lines import Line2D
        proxies = [Line2D([], [], marker="o", mfc="none", mec="k", ls="",
                          ms=4.5, label="circles = measured (every 4th pt)"),
                   Line2D([], [], color="k", lw=1.7, label="lines = global fit")]
        leg1 = ax.legend(handles=proxies, frameon=False, fontsize=8.5,
                         loc="upper right")
        ax.add_artist(leg1)
        ax.legend(frameon=False, fontsize=8, loc="center right",
                  title="temperature")
        ax.set_ylabel("extinction (norm.)")
        ax.set_title(f"Global fit to real T-series — {len(show)} of "
                     f"{len(temps_K)} temperatures shown (fit uses all)")
        axr.axhspan(-r.residual_rms, r.residual_rms, color="0.85", alpha=0.7)
        axr.axhline(0, color="0.4", lw=0.7)
        axr.set_ylabel("data − fit")
        axr.set_xlabel("wavelength (nm)")
    else:
        print("single spectrum -> single-spectrum fit")
        tc = float(temps_K[0] - 273.15) if (temps_K is not None and EPS_T) else None
        r = fit_spectrum(wl, spectra[0], species=SPECIES, gap_nm=1.0,
                         backend=backend, n_sizes=5, temperature_C=tc,
                         size_correction=size_corr)
        print(f"\n  diameter = {r.diameter_nm:.2f} ± {r.diameter_sd:.1f} nm  "
              f"[{r.identifiability['diameter']}]")
        print(f"  poly     = {r.pct_poly:.2f} %  [{r.identifiability['pct_poly']}]")
        print(f"  aggregated gold = {r.aggregated_gold:.3f} ± "
              f"{r.aggregated_gold_sd:.3f}  [{r.identifiability['aggregated_gold']}]")
        print(f"  pedestal a_sca = {r.a_sca:.4f} @550 nm, n_sca = {r.n_sca:.2f} "
              "(weakly identified from one spectrum)")
        fig, ax = plt.subplots(figsize=(7, 4.4))
        ax.plot(wl, spectra[0] / spectra[0].max(), ".", ms=3, color="0.5", label="data")
        ax.plot(r.wavelength, r.model, "-", color="#D55E00", lw=2, label="fit")
        ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("extinction (norm.)")
        ax.legend(frameon=False)
    out = os.path.join(OUT, "fig8_fit_real.png")
    fig.tight_layout(); fig.savefig(out, dpi=130)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
