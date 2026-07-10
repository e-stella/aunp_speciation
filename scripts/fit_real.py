"""Fit a real UV-Vis file (or the bundled example).

Usage:
    python scripts/fit_real.py [path/to/spectra.csv] [--range MIN_NM MAX_NM]

- Auto-detects a temperature series (>=2 T columns) -> global fit; otherwise a
  single-spectrum fit.
- --range restricts the fitted window (default 420 800 nm); pass "--range 0 1e9"
  to keep the full spectrum. See "Known limitations" #7 in CLAUDE.md.
- Uses the precomputed EXACT T-matrix basis (outputs/tmatrix_basis.npz) if it
  exists, else falls back to the fast CDA optics. (Build the cache once with
  mstm-env/bin/python scripts/build_tmatrix_basis.py)
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from aunp_speciation.io_data import load_series
from aunp_speciation.fitting import fit_spectrum
from aunp_speciation.fit_global import fit_temperature_series

HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "..", "outputs")
CACHE = os.path.join(OUT, "tmatrix_basis.npz")
SPECIES = ("monomer", "dimer", "trimer_linear")


def parse_args(argv):
    """Return (data_path, wavelength_range). Supports an optional
    `--range MIN MAX` flag anywhere after the program name."""
    args = list(argv[1:])
    wl_range = (420.0, 800.0)
    if "--range" in args:
        i = args.index("--range")
        wl_range = (float(args[i + 1]), float(args[i + 2]))
        del args[i:i + 3]
    data = args[0] if args else os.path.join(HERE, "..", "data",
                                             "example_series.csv")
    return data, wl_range


DATA, WL_RANGE = parse_args(sys.argv)


def choose_backend():
    if os.path.exists(CACHE):
        from aunp_speciation.basis_cache import load_cache
        print(f"using EXACT cached T-matrix optics: {os.path.basename(CACHE)}")
        return load_cache(CACHE).species_fn
    print("cache not found -> using fast CDA optics (lower-bound coupling)")
    return "cda"


def main():
    wl, spectra, temps_K = load_series(DATA, wavelength_range=WL_RANGE)
    backend = choose_backend()
    print(f"loaded {os.path.basename(DATA)}: {spectra.shape[0]} spectrum(a), "
          f"{len(wl)} wavelengths, {wl.min():.0f}-{wl.max():.0f} nm")

    if temps_K is not None and spectra.shape[0] >= 2:
        print(f"temperature series: {temps_K - 273.15} C -> GLOBAL fit")
        r = fit_temperature_series(temps_K, wl, spectra, species=SPECIES,
                                   backend=backend, n_sizes=5, fit_stride=2,
                                   C_tot=0.01, max_nfev=120)
        agg = r.gold_fractions["dimer"] + r.gold_fractions["trimer_linear"]
        print(f"\n  diameter   = {r.diameter_nm:.2f} nm (± {r.param_sd['diameter']:.1f})")
        print(f"  poly       = {r.pct_poly:.2f} %  (± {r.param_sd['pct_poly']:.0f})")
        print(f"  dH2        = {r.dH2:.1f} kJ/mol (± {r.param_sd['dH2']:.1f})")
        print(f"  dS2        = {r.dS2:.3f} kJ/mol/K (± {r.param_sd['dS2']:.3f})")
        print("  aggregated gold fraction vs T:")
        for Tc, a in zip(temps_K - 273.15, agg):
            print(f"     {Tc:5.0f} C : {a:.3f}")
        print(f"  residual RMS = {r.residual_rms:.4f}")
        fig, ax = plt.subplots(figsize=(7, 4.6))
        cmap = plt.cm.coolwarm(np.linspace(0, 1, len(temps_K)))
        for i, Tc in enumerate(temps_K - 273.15):
            ax.plot(wl, spectra[i] / spectra.max(), ".", ms=2.5, color=cmap[i])
            ax.plot(r.wavelength, r.models[i], "-", color=cmap[i], lw=1.5,
                    label=f"{Tc:.0f} C")
        ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("extinction (norm.)")
        ax.set_title("Global fit to real T-series (points=data, lines=fit)")
        ax.legend(frameon=False, fontsize=8, ncol=2)
    else:
        print("single spectrum -> single-spectrum fit")
        r = fit_spectrum(wl, spectra[0], species=SPECIES, gap_nm=1.0)
        print(f"\n  diameter = {r.diameter_nm:.2f} ± {r.diameter_sd:.1f} nm  "
              f"[{r.identifiability['diameter']}]")
        print(f"  poly     = {r.pct_poly:.2f} %  [{r.identifiability['pct_poly']}]")
        print(f"  aggregated gold = {r.aggregated_gold:.3f} ± "
              f"{r.aggregated_gold_sd:.3f}  [{r.identifiability['aggregated_gold']}]")
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
