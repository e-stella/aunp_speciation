"""Exact T-matrix cluster basis for the seeded-growth speciation re-fit.

Builds orientation-averaged dimer / trimer_linear / tetramer_linear extinction
columns for the four aggregated seeded-growth samples (Oct GNP20/40/60, Nov
GNP55) at their SAM tier-2 sizes, with the ROUTE-CALIBRATED damping baked in
(s=1.75, A_surf=1.0 — CLAUDE.md #2; the module-level s knob is honored live
by clusters_tmatrix via gold_epsilon_sized). Cache: outputs/
seeded_tmatrix_basis.npz (monomer columns are NOT cached — Mie is cheap and
exact, compute at fit time with the same damping).

--probe mode runs the lmax-convergence study FIRST (dimer, the standard
proxy): successive-lmax relative changes at 600/700 nm per (D, gap). Fill
GAP_LMAX below from its output before building. CLAUDE.md #10 already warns
gap 1 nm is unconverged at lmax=10 for D>=29 — expect the small-gap columns
to be usable only at small D or high lmax.

Run (mstm-env; one shard per sample so the four can run in parallel,
~30-75 min each):
  mstm-env/bin/python scripts/build_seeded_tmatrix_basis.py --probe   (~5 min)
  mstm-env/bin/python scripts/build_seeded_tmatrix_basis.py --sample oct_GNP20
"""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np

from aunp_speciation import dielectric
dielectric.use_gold_model("jc")
from aunp_speciation.clusters_tmatrix import species_spectrum_tmatrix, available

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT_TPL = os.path.join(ROOT, "outputs", "seeded_tmatrix_basis_{sid}.npz")

# route calibration (scripts/calibrate_damping_seeded.py, 2026-07-13)
S_CAL, A_CAL = 1.75, 1.0
GOLD_MODEL = "jc"

# SAM tier-2 (mean nm, poly %) — identical to fit_speciation_seeded.py inputs
SAMPLES = {
    "oct_GNP20": (22.8, 11.9),
    "oct_GNP40": (42.5, 11.8),
    "oct_GNP60": (63.3, 12.7),
    "nov_GNP55": (54.7, 19.5),
}
SPECIES = ("dimer", "trimer_linear", "tetramer_linear")

# Per-sample (gap_nm -> lmax) build plan, from the --probe run (2026-07-13):
#   gap 1.0 converges (<1%) only at D=22.8 (l12 +0.5%); at D>=42.5 it is
#     UNCONVERGED even at l14 (+2.7/+7.6/+2.4%) -> dropped there.
#   gap 2.0: l10 +0.3% (22.8), l12 +0.4% (42.5), l14 +0.8% (54.7/63.3);
#   gap 3.5: l8 +0.1% (22.8), l10 +0.3-0.4% (42.5/54.7), l12 +0.3% (63.3).
# The l12/l10 choices for the two largest samples accept a ~1-2% residual on
# the largest polydispersity bins (l14 would double an already-hour-long
# build). Direction of that residual: under-converged = red tail LOW =
# fitted aggregate fraction slightly HIGH -> the exact fractions stay upper
# bounds, same sign as the CDA bias and much smaller.
# ⚠ STABILITY (found 2026-07-13): at gap 1.0 the TRIMER/TETRAMER interaction
# solve goes numerically unstable at l12 (negative/spiky Cext beyond ~650 nm
# at D=22.8) even though the DIMER is fine there — the probe's dimer proxy
# does NOT catch it. l10 is smooth. Cap chain species at l10 for gap <= 1.
PLAN = {
    "oct_GNP20": {1.0: 10, 2.0: 10, 3.5: 8},
    "oct_GNP40": {2.0: 12, 3.5: 10},
    "oct_GNP60": {2.0: 12, 3.5: 10},
    "nov_GNP55": {2.0: 12, 3.5: 10},
}

WL = np.arange(450.0, 741.0, 10.0)   # cluster columns are smooth; fit splines
N_SIZES = 3                          # cluster polydispersity is approximate
                                     # anyway (CLAUDE.md limitation #3)

PROBE_D = (22.8, 42.5, 54.7, 63.3)
PROBE_GAPS = (1.0, 2.0, 3.5)
PROBE_LMAX = (6, 8, 10, 12, 14)
PROBE_WL = (600.0, 700.0)


def gaussian_sizes(mean_nm, pct, n):
    sd = mean_nm * pct / 100.0
    ds = np.linspace(mean_nm - 2 * sd, mean_nm + 2 * sd, n)
    w = np.exp(-0.5 * ((ds - mean_nm) / sd) ** 2)
    return ds, w / w.sum()


def probe():
    print("dimer lmax convergence, calibrated damping "
          f"(s={S_CAL}, A={A_CAL}); values = Cext(600), Cext(700); "
          "pct = change vs previous lmax\n")
    for D in PROBE_D:
        for gap in PROBE_GAPS:
            prev = None
            line = f"D={D:5.1f} gap={gap:3.1f}: "
            for lmax in PROBE_LMAX:
                t = time.time()
                c = species_spectrum_tmatrix("dimer", D, PROBE_WL, gap_nm=gap,
                                             lmax=lmax, size_correction=True,
                                             A_surf=A_CAL)
                dt = time.time() - t
                if prev is None:
                    line += f"l{lmax}[{c[0]:.0f},{c[1]:.0f}] "
                else:
                    pct = 100 * np.max(np.abs(c - prev) / prev)
                    line += f"l{lmax}[{pct:+.1f}%,{dt:.1f}s] "
                    if pct < 1.0:      # converged; skip higher lmax
                        prev = c
                        break
                prev = c
            print(line, flush=True)
    print("\nfill GAP_LMAX from the first lmax with <1% change, then rerun "
          "without --probe")


def build(sid):
    D0, poly = SAMPLES[sid]
    plan = PLAN[sid]
    dielectric.set_bulk_damping_scale(S_CAL)
    data = {"wl": WL}
    meta = dict(gold_model=GOLD_MODEL, s_bulk=S_CAL, A_surf=A_CAL,
                n_sizes=N_SIZES, species=",".join(SPECIES),
                gaps=",".join(f"{g}:{l}" for g, l in sorted(plan.items())))
    t0 = time.time()
    ds, w = gaussian_sizes(D0, poly, N_SIZES)
    for gap, lmax in sorted(plan.items()):
        for sp in SPECIES:
            t = time.time()
            col = np.zeros_like(WL)
            for d, wi in zip(ds, w):
                col += wi * species_spectrum_tmatrix(
                    sp, d, WL, gap_nm=gap, lmax=lmax,
                    size_correction=True, A_surf=A_CAL)
            data[f"{sp}/gap{gap:g}"] = col
            print(f"{sid} {sp} gap={gap:g} lmax={lmax}: "
                  f"{time.time()-t:.0f} s  (total {(time.time()-t0)/60:.1f} min)",
                  flush=True)
    out = OUT_TPL.format(sid=sid)
    np.savez(out, **data, **{f"meta_{k}": v for k, v in meta.items()})
    print(f"\nwrote {out}")


if __name__ == "__main__":
    if not available():
        sys.exit("treams not importable — run with mstm-env/bin/python")
    if "--probe" in sys.argv:
        probe()
    elif "--sample" in sys.argv:
        build(sys.argv[sys.argv.index("--sample") + 1])
    else:
        sys.exit("usage: --probe | --sample <sid>  (one shard per sample; "
                 f"sids: {', '.join(SAMPLES)})")
