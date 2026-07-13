"""Speciation re-fit for the seeded-growth series: old vs recalibrated damping.

The HANDOFF.md deliverable. Mirrors ../tem-particle-metrics/scripts/
fit_speciation.py exactly (SAM tier-2 sizes, monomer+dimer+trimer_linear CDA
basis, gap scan, NNLS on the raw absorbance, 450-740 nm) and runs it twice:

  first-cut : s=1, A_surf=0.25, no offset — reproduces the sister repo's
              table (oct_GNP12 13% "agg" etc.) as the baseline.
  seeded-cal: the route-calibrated (s, A_surf) from
              scripts/calibrate_damping_seeded.py, plus the -2.7 nm jc offset
              (part of the calibration; CLAUDE.md #9).

Diagnostic (verbatim from HANDOFF.md): if oct_GNP20's "aggregation" collapses
toward 0 under the recalibrated monomer line-width, the Oct aggregation was a
model artifact; if it stays, it is a real sample-specific effect.

Run:  python scripts/fit_speciation_seeded.py
Writes outputs/fig21_speciation_seeded.png.
"""
import csv, io, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from scipy.optimize import nnls

from aunp_speciation import dielectric
dielectric.use_gold_model("jc")
from aunp_speciation.spectra import species_basis

ROOT = os.path.join(os.path.dirname(__file__), "..")
TEM = os.path.join(ROOT, "..", "tem-particle-metrics")
IRINA = os.path.join(ROOT, "experimental", "Irina")
OUT = os.path.join(ROOT, "outputs")

# route calibration (scripts/calibrate_damping_seeded.py, 2026-07-13):
# joint fit on the four verified-clean monomers (oct12 + nov11/23/39),
# mean RMS 7.74% -> 3.00%. The s/A split sits on a shallow ridge (gamma_eff
# is the constrained combination); the gnp12-only corner of the ridge is the
# robustness alternative below.
SEED_S, SEED_A, SEED_SHIFT = 1.75, 1.0, 2.7
ALT_PAIR = (3.0, 0.75)  # gnp12-only pair, robustness check on the Oct verdict

WL = np.arange(450.0, 741.0, 2.0)
SPECIES = ("monomer", "dimer", "trimer_linear")
N_GOLD = {"monomer": 1, "dimer": 2, "trimer_linear": 3}
GAPS = [0.5, 1.0, 2.0, 3.0]

OCT, NOV = "Scan_10_22_2024_5_43_15_PM.csv", "Scan_11_16_2024_4_37_22_PM.csv"
SAMPLES = {
    "oct_GNP12": (OCT, "GNP12nm", "Oct GNP12"),
    "oct_GNP20": (OCT, "GNP20nm", "Oct GNP20"),
    "oct_GNP40": (OCT, "GNP40nm", "Oct GNP40"),
    "oct_GNP60": (OCT, "GNP60nm", "Oct GNP60"),
    "nov_GNP11": (NOV, "Seed 1(Abs)", "Nov Seed 11"),
    "nov_GNP23": (NOV, "30nm 1(Abs)", "Nov GNP23"),
    "nov_GNP39": (NOV, "40nm 1(Abs)", "Nov GNP39"),
    "nov_GNP55": (NOV, "60nm 1(Abs)", "Nov GNP55"),
}


def sam_stats(sid):
    """SAM tier-2 (mean nm, polydispersity %, TEM touching fraction) — same
    inputs as the sister repo's first-cut fit."""
    diam, touch = [], 0
    with open(os.path.join(TEM, "outputs", "dm3", "tier2",
                           f"{sid}_particles.csv"), newline="") as f:
        for row in csv.DictReader(f):
            diam.append(float(row["diameter_nm"]))
            if row.get("touching_group_id", "").strip():
                touch += 1
    d = np.array(diam)
    return float(d.mean()), float(d.std(ddof=1) / d.mean() * 100), touch / len(d)


def read_measured(csv_name, column):
    raw = open(os.path.join(IRINA, csv_name), encoding="utf-16").read()
    rows = list(csv.reader(io.StringIO(raw), delimiter="\t"))
    header = [h.strip() for h in rows[0]]
    ci = header.index(column)
    wl, ab = [], []
    for r in rows[1:]:
        if len(r) > ci and r[0].strip():
            try:
                wl.append(float(r[0])); ab.append(float(r[ci]))
            except ValueError:
                pass
    wl, ab = np.array(wl), np.array(ab)
    o = np.argsort(wl)
    return np.interp(WL, wl[o], ab[o])


def fit_one(mean_nm, poly, meas, s, A, shift_nm):
    """NNLS species mixture at the best-fitting gap, under damping (s, A)."""
    dielectric.set_bulk_damping_scale(s)
    try:
        best = None
        for gap in GAPS:
            basis = species_basis(WL, mean_nm, poly, gap_nm=gap,
                                  species=SPECIES, backend="cda",
                                  size_correction=True, A_surf=A)
            M = np.vstack([np.interp(WL, WL - shift_nm, basis[sp])
                           for sp in SPECIES]).T
            w, rnorm = nnls(M, meas)
            if best is None or rnorm < best["rnorm"]:
                best = dict(gap=gap, w=w, rnorm=rnorm, M=M, fit=M @ w)
    finally:
        dielectric.set_bulk_damping_scale(1.0)
    gold = np.array([best["w"][i] * N_GOLD[sp] for i, sp in enumerate(SPECIES)])
    best["agg"] = 0.0 if gold.sum() == 0 else float(gold[1:].sum() / gold.sum())
    return best


def main():
    if SEED_S is None:
        sys.exit("fill SEED_S/SEED_A from scripts/calibrate_damping_seeded.py first")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.grid": True,
                         "grid.alpha": 0.25, "axes.axisbelow": True})
    fig, axes = plt.subplots(2, 4, figsize=(20, 9))

    print(f"seeded pair: s={SEED_S}, A_surf={SEED_A}, shift -{SEED_SHIFT} nm\n")
    print(f"{'sample':11} {'SAM_nm':>7} | {'agg_old':>8} {'gap':>4} {'rms':>6} | "
          f"{'agg_new':>8} {'gap':>4} {'rms':>6} | {'TEM_touch':>9}")
    print("-" * 78)
    rows = []
    for ax, (sid, (csv_name, col, title)) in zip(axes.ravel(), SAMPLES.items()):
        mean, poly, touch = sam_stats(sid)
        meas = read_measured(csv_name, col)
        old = fit_one(mean, poly, meas, 1.0, 0.25, 0.0)      # first-cut config
        new = fit_one(mean, poly, meas, SEED_S, SEED_A, SEED_SHIFT)
        r_old = float(np.sqrt(np.mean((old["fit"] - meas) ** 2)) / meas.max())
        r_new = float(np.sqrt(np.mean((new["fit"] - meas) ** 2)) / meas.max())
        rows.append((sid, old["agg"], new["agg"]))
        print(f"{sid:11} {mean:7.1f} | {old['agg']*100:7.0f}% {old['gap']:4.1f} "
              f"{r_old:6.3f} | {new['agg']*100:7.0f}% {new['gap']:4.1f} "
              f"{r_new:6.3f} | {touch*100:8.0f}%")

        scale = 1.0 / meas.max()
        ax.plot(WL, meas * scale, "k", lw=2.2, label="measured")
        ax.plot(WL, old["fit"] * scale, color="#999", lw=1.3, ls=":",
                label=f"old damping (agg {old['agg']*100:.0f}%)")
        ax.plot(WL, new["fit"] * scale, color="#d1495b", lw=1.6,
                label=f"recalibrated (agg {new['agg']*100:.0f}%)")
        mono = new["w"][0] * new["M"][:, 0]
        ax.plot(WL, mono * scale, "#3b7dd8", lw=1.2, ls="--", label="monomer part")
        ax.set_title(f"{title}  (SAM {mean:.0f} nm)", fontsize=11)
        ax.set_xlabel("wavelength (nm)")
        ax.legend(fontsize=8)

    fig.suptitle("Speciation with TEM-pinned size: first-cut (s=1, A=0.25) vs "
                 f"seeded-route calibrated damping (s={SEED_S}, A={SEED_A}, "
                 f"−{SEED_SHIFT} nm)", fontsize=13)
    fig.tight_layout()
    out = os.path.join(OUT, "fig21_speciation_seeded.png")
    fig.savefig(out, dpi=110)
    print(f"\nwrote {out}")

    print("\nHANDOFF diagnostic (agg%, first-cut -> recalibrated"
          + (" -> gnp12-only pair" if ALT_PAIR else "") + "):")
    for sid, a_old, a_new in rows:
        if sid.startswith("oct_"):
            line = f"  {sid}: {a_old*100:.0f}% -> {a_new*100:.0f}%"
            if ALT_PAIR:
                mean, poly, _ = sam_stats(sid)
                meas = read_measured(*SAMPLES[sid][:2])
                alt = fit_one(mean, poly, meas, *ALT_PAIR, SEED_SHIFT)
                line += f" -> {alt['agg']*100:.0f}%"
            print(line)


if __name__ == "__main__":
    main()
