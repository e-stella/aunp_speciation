"""Exact T-matrix speciation re-fit for the aggregated seeded-growth samples.

Follow-up to fit_speciation_seeded.py (fig21): CDA under-couples at
near-contact, so its aggregated fractions are upper bounds. This re-fits Oct
GNP20/40/60 + Nov GNP55 with the exact multipole basis from
scripts/build_seeded_tmatrix_basis.py (route-calibrated damping s=1.75,
A_surf=1.0, -2.7 nm; dimer + linear trimer + linear tetramer), and reports
CDA-vs-exact aggregated gold fractions side by side. The monomer column is
Mie (exact for a sphere) computed live with the same damping.

Run:  python scripts/fit_speciation_seeded_exact.py     (system env)
Writes outputs/fig22_speciation_exact.png.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from scipy.optimize import nnls

from aunp_speciation import dielectric
dielectric.use_gold_model("jc")
from aunp_speciation.spectra import monomer_polydisperse, species_basis

import fit_speciation_seeded as F        # data loaders + SAM stats + config

ROOT = os.path.join(os.path.dirname(__file__), "..")
CACHE_TPL = os.path.join(ROOT, "outputs", "seeded_tmatrix_basis_{sid}.npz")
S_CAL, A_CAL, SHIFT = F.SEED_S, F.SEED_A, F.SEED_SHIFT
WL = F.WL
SPECIES = ("monomer", "dimer", "trimer_linear", "tetramer_linear")
N_GOLD = {"monomer": 1, "dimer": 2, "trimer_linear": 3, "tetramer_linear": 4}
TARGETS = ("oct_GNP20", "oct_GNP40", "oct_GNP60", "nov_GNP55")


def shift(col):
    return np.interp(WL, WL - SHIFT, col)


def agg_of(w):
    gold = np.array([w[i] * N_GOLD[sp] for i, sp in enumerate(SPECIES)])
    return 0.0 if gold.sum() == 0 else float(gold[1:].sum() / gold.sum())


def fit_basis(columns_by_gap, meas):
    """NNLS over the gap grid; columns_by_gap: {gap: (n_wl, n_species)}."""
    best = None
    for gap, M in sorted(columns_by_gap.items()):
        w, rnorm = nnls(M, meas)
        if best is None or rnorm < best["rnorm"]:
            best = dict(gap=gap, w=w, rnorm=rnorm, M=M, fit=M @ w)
    best["agg"] = agg_of(best["w"])
    return best


def main():
    from scipy.interpolate import CubicSpline
    dielectric.set_bulk_damping_scale(S_CAL)
    try:
        rows = []
        panels = {}
        for sid in TARGETS:
            z = np.load(CACHE_TPL.format(sid=sid))
            assert str(z["meta_gold_model"]) == "jc" \
                and float(z["meta_s_bulk"]) == S_CAL \
                and float(z["meta_A_surf"]) == A_CAL, "cache/calibration mismatch"
            wlc = z["wl"]
            gaps = sorted(float(k.split(":")[0])
                          for k in str(z["meta_gaps"]).split(","))
            mean, poly, _ = F.sam_stats(sid)
            meas = F.read_measured(*F.SAMPLES[sid][:2])
            mono = shift(monomer_polydisperse(WL, mean, poly, "water", None,
                                              True, A_CAL))
            # exact basis from the cache (coarse grid -> cubic spline; a gap
            # may carry its own finer axis as wl_gap<g>, e.g. the near-contact
            # gap-1 columns whose bonding-mode structure needs 5 nm sampling)
            ex_cols = {}
            for gap in gaps:
                wlg = z[f"wl_gap{gap:g}"] if f"wl_gap{gap:g}" in z.files else wlc
                cols = [mono] + [
                    shift(CubicSpline(wlg, z[f"{sp}/gap{gap:g}"])(WL))
                    for sp in SPECIES[1:]]
                ex_cols[gap] = np.vstack(cols).T
            exact = fit_basis(ex_cols, meas)
            sens = {}
            for gap, M in sorted(ex_cols.items()):
                wg, _ = nnls(M, meas)
                sens[gap] = (agg_of(wg),
                             float(np.sqrt(np.mean((M @ wg - meas) ** 2)) / meas.max()))
            exact["sens"] = sens
            # CDA with the SAME four species/gaps (like-for-like comparison)
            cda_cols = {}
            for gap in gaps:
                b = species_basis(WL, mean, poly, gap_nm=gap, species=SPECIES,
                                  backend="cda", size_correction=True,
                                  A_surf=A_CAL)
                cda_cols[gap] = np.vstack([shift(b[sp]) for sp in SPECIES]).T
            cda = fit_basis(cda_cols, meas)
            rows.append((sid, mean, cda, exact,
                         float(np.sqrt(np.mean((cda["fit"] - meas) ** 2)) / meas.max()),
                         float(np.sqrt(np.mean((exact["fit"] - meas) ** 2)) / meas.max())))
            panels[sid] = (meas, cda, exact, mean)
    finally:
        dielectric.set_bulk_damping_scale(1.0)

    print(f"calibrated damping (s={S_CAL}, A={A_CAL}, -{SHIFT} nm); "
          f"species: {'+'.join(SPECIES)}; per-sample gap grids from the cache\n")
    print(f"{'sample':11} {'SAM_nm':>7} | {'agg_CDA':>8} {'gap':>4} {'rms':>6} | "
          f"{'agg_EXACT':>9} {'gap':>4} {'rms':>6}")
    print("-" * 66)
    for sid, mean, cda, exact, r_c, r_e in rows:
        print(f"{sid:11} {mean:7.1f} | {cda['agg']*100:7.0f}% {cda['gap']:4.1f} "
              f"{r_c:6.3f} | {exact['agg']*100:8.0f}% {exact['gap']:4.1f} {r_e:6.3f}")
    print("\ngap sensitivity of the EXACT aggregated fraction (agg% @ rms):")
    for sid, mean, cda, exact, _, _ in rows:
        s = "  ".join(f"gap{g:g}: {a*100:.0f}% @ {r:.3f}"
                      for g, (a, r) in sorted(exact["sens"].items()))
        print(f"  {sid}: {s}")
    print("\nper-species gold split (EXACT fit):")
    for sid, mean, cda, exact, _, _ in rows:
        gold = np.array([exact["w"][i] * N_GOLD[sp] for i, sp in enumerate(SPECIES)])
        tot = gold.sum() if gold.sum() > 0 else 1.0
        split = "  ".join(f"{sp.split('_')[0]}={g/tot*100:.0f}%"
                          for sp, g in zip(SPECIES, gold))
        print(f"  {sid}: {split}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.grid": True,
                         "grid.alpha": 0.25, "axes.axisbelow": True})
    fig, axes = plt.subplots(1, 4, figsize=(20, 4.6))
    for ax, sid in zip(axes, TARGETS):
        meas, cda, exact, mean = panels[sid]
        scale = 1.0 / meas.max()
        ax.plot(WL, meas * scale, "k", lw=2.2, label="measured")
        ax.plot(WL, cda["fit"] * scale, color="#999", lw=1.3, ls=":",
                label=f"CDA (agg {cda['agg']*100:.0f}%, gap {cda['gap']:g})")
        ax.plot(WL, exact["fit"] * scale, color="#d1495b", lw=1.6,
                label=f"exact (agg {exact['agg']*100:.0f}%, gap {exact['gap']:g})")
        ax.plot(WL, exact["w"][0] * exact["M"][:, 0] * scale, "#3b7dd8",
                lw=1.2, ls="--", label="monomer part")
        ax.set_title(f"{sid}  (SAM {mean:.0f} nm)", fontsize=10)
        ax.set_xlabel("wavelength (nm)")
        ax.legend(fontsize=7)
    fig.suptitle("Aggregated seeded samples: CDA vs exact T-matrix basis, both with "
                 f"route-calibrated damping (s={S_CAL}, A={A_CAL}, −{SHIFT} nm), "
                 "monomer+dimer+trimer+tetramer NNLS", fontsize=11)
    fig.tight_layout()
    out = os.path.join(ROOT, "outputs", "fig22_speciation_exact.png")
    fig.savefig(out, dpi=110)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
