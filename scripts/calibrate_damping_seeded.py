"""Joint damping calibration for the SEEDED-GROWTH citrate/H2O2 route
(Oct + Nov 2024 series; HANDOFF.md).

Same methodology as the frozen CTAC calibration (calibrate_damping.py):
gamma_eff(D) = s*gamma_bulk + A*hbar*v_F/(D/2), ONE global (s, A) pair per
synthesis route, grid-scanned to minimize the mean peak-normalized RMS
(450-700 nm) over TEM-pinned per-particle histograms, with the documented
-2.7 nm jc offset applied (CLAUDE.md #9 — the offset transfers across
chemistries; the damping pair does not, CLAUDE.md #2).

Calibration set = the VERIFIED-CLEAN monomers of the route only:
  oct_GNP12 (11.6 nm)  — the handoff's pure-monomer reference
  nov_GNP11 (10.9 nm), nov_GNP23 (22.9 nm), nov_GNP39 (38.5 nm)
  (first-cut speciation: 13/6/0/7 % "agg"; the 13% on a 12 nm sol is the
  line-width artifact under test). GNP55/60 (really aggregated) and Oct
  GNP20/40 (the samples under test) are EXCLUDED.

Three variants are reported to localize any Oct-vs-Nov batch effect:
  gnp12-only  — the handoff's literal plan (constrains mostly A: gamma_S
                dominates at small D, s is nearly free there)
  nov-only    — Nov batch alone (11/23/39 nm spans the 1/R lever arm)
  joint       — all four clean monomers (production candidate)

Per-particle diameters: hand-picked ground truth from the sister repo,
../tem-particle-metrics/outputs/dm3/gt/*.csv (preferred over SAM tier-2
until its +5% mask-boundary correction lands).

Run:  python scripts/calibrate_damping_seeded.py     (~3 min)
Writes outputs/fig20_calibrate_seeded.png.
"""
import sys, os, csv
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np

from aunp_speciation import dielectric
dielectric.use_gold_model("jc")
from aunp_speciation.mie import mie_ab
from aunp_speciation.dielectric import gold_epsilon, medium_index

ROOT = os.path.join(os.path.dirname(__file__), "..")
TEM_GT = os.path.join(ROOT, "..", "tem-particle-metrics", "outputs", "dm3", "gt")
IRINA = os.path.join(ROOT, "experimental", "Irina")
OUT = os.path.join(ROOT, "outputs")

WL = np.arange(420.0, 801.0, 1.0)
E_EV = 1239.842 / WL
EPS_JC = gold_epsilon(WL)
N_MED = medium_index("water", None).real
OM = getattr(dielectric, "_OMEGA_P_EV", 8.45)
G0 = getattr(dielectric, "_GAMMA_BULK_EV", 0.044)
HBAR_VF = 0.9215
SHIFT_NM = 2.7                       # jc red bias, transfers (CLAUDE.md #9/#2)

# First pass (0.05-1.5 x 0.1-1.25) railed at the (1.5, 1.0) corner — the
# seeded route is MORE damped than jc bulk, opposite the CTAC result; grid
# extended upward accordingly.
S_GRID = (0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0)
A_GRID = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)
CTAC_PAIR = (0.05, 0.6)              # frozen 2026-07-12 (CTAC route)
BASE_PAIR = (1.0, 0.25)              # uncalibrated production default

# sample -> (scan csv, column header, TEM gt file)
OCT, NOV = "Scan_10_22_2024_5_43_15_PM.csv", "Scan_11_16_2024_4_37_22_PM.csv"
CAL_SAMPLES = {
    "oct_GNP12": (OCT, "GNP12nm", "oct_GNP12.csv"),
    "nov_GNP11": (NOV, "Seed 1(Abs)", "nov_GNP11.csv"),
    "nov_GNP23": (NOV, "30nm 1(Abs)", "nov_GNP23.csv"),
    "nov_GNP39": (NOV, "40nm 1(Abs)", "nov_GNP39.csv"),
}
VARIANTS = {
    "gnp12-only": ("oct_GNP12",),
    "nov-only": ("nov_GNP11", "nov_GNP23", "nov_GNP39"),
    "joint": tuple(CAL_SAMPLES),
}


def load_utf16_tsv(path):
    with open(path, "r", encoding="utf-16") as f:
        rows = [r for r in csv.reader(f, delimiter="\t") if r and r[0].strip()]
    hdr = [c.strip() for c in rows[0] if c.strip()]
    data = []
    for r in rows[1:]:
        try:
            vals = [float(c) for c in r[:len(hdr)] if c.strip() != ""]
        except ValueError:
            continue
        if len(vals) == len(hdr):
            data.append(vals)
    arr = np.array(data)
    o = np.argsort(arr[:, 0])
    return hdr, arr[o, 0], arr[o, 1:]


def drude(g):
    return OM ** 2 / (E_EV ** 2 + 1j * g * E_EV)


def xsec(D, s, A):
    eps = EPS_JC + drude(G0) - drude(s * G0 + A * HBAR_VF / (D / 2.0))
    m = np.sqrt(eps) / N_MED
    C = np.zeros_like(WL)
    for i, l0 in enumerate(WL):
        k = 2.0 * np.pi * N_MED / l0
        a, b = mie_ab(m[i], k * D / 2.0)
        n = np.arange(1, len(a) + 1)
        C[i] = (2.0 * np.pi / k ** 2) * np.sum((2 * n + 1) * np.real(a + b))
    return C


def model(S, s, A):
    ym = np.sum([wi * xsec(d, s, A) for d, wi in zip(S["bins"], S["w"])], axis=0)
    return np.interp(WL, WL - SHIFT_NM, ym)     # -2.7 nm offset


MSK = (WL >= 450) & (WL <= 700)


def rms(ym, y):
    a = ym[MSK] / ym[MSK].max()
    b = y[MSK] / y[MSK].max()
    return 100 * np.sqrt(np.mean((a - b) ** 2))


def peak_of(wl, y, lo=480.0, hi=650.0):
    m = (wl >= lo) & (wl <= hi)
    k = np.convolve(y[m], np.ones(5) / 5, mode="same")
    return wl[m][np.argmax(k)]


# --- load calibration samples ------------------------------------------------
SAMPLES = {}
scans = {}
for sid, (scan, col, gtf) in CAL_SAMPLES.items():
    gt_path = os.path.join(TEM_GT, gtf)
    if not os.path.exists(gt_path):
        sys.exit(f"missing TEM ground truth {gt_path} — clone/build "
                 "../tem-particle-metrics first")
    d = np.loadtxt(gt_path)
    edges = np.arange(d.min() - 0.25, d.max() + 0.75, 0.5)
    cnt, _ = np.histogram(d, edges)
    ctr = 0.5 * (edges[:-1] + edges[1:])
    keep = cnt > 0
    if scan not in scans:
        scans[scan] = load_utf16_tsv(os.path.join(IRINA, scan))
    hdr, wl_d, Y = scans[scan]
    ci = hdr.index(col) - 1                     # hdr[0] is the wl column
    SAMPLES[sid] = dict(bins=ctr[keep], w=cnt[keep] / cnt.sum(),
                        y=np.interp(WL, wl_d, Y[:, ci]),
                        D=float(d.mean()), sd=float(d.std(ddof=1)), n=len(d))
    print(f"{sid}: TEM D = {d.mean():.2f} +/- {d.std(ddof=1):.2f} nm "
          f"(n={len(d)}), data peak {peak_of(WL, SAMPLES[sid]['y']):.0f} nm")

# --- RMS grid over all (s, A) for every sample (variants share it) -----------
print(f"\nscanning s in {S_GRID}\n         A in {A_GRID} "
      f"({len(S_GRID) * len(A_GRID)} pairs, shift -{SHIFT_NM} nm) ...")
RMS = {sid: np.zeros((len(S_GRID), len(A_GRID))) for sid in SAMPLES}
for i, s in enumerate(S_GRID):
    for j, A in enumerate(A_GRID):
        for sid, S in SAMPLES.items():
            RMS[sid][i, j] = rms(model(S, s, A), S["y"])

results = {}
for name, sids in VARIANTS.items():
    tot = np.mean([RMS[sid] for sid in sids], axis=0)
    i, j = np.unravel_index(np.argmin(tot), tot.shape)
    results[name] = (S_GRID[i], A_GRID[j], tot[i, j])

print(f"\n{'variant':12s} {'s':>5s} {'A_surf':>7s} {'meanRMS%':>9s}   per-sample RMS%")
for name, (s, A, r) in results.items():
    si, ai = S_GRID.index(s), A_GRID.index(A)
    per = "  ".join(f"{sid.split('_')[1]}={RMS[sid][si, ai]:.2f}"
                    for sid in CAL_SAMPLES)
    print(f"{name:12s} {s:5.2f} {A:7.2f} {r:9.2f}   {per}")
for label, (s, A) in (("BASELINE", BASE_PAIR), ("CTAC pair", CTAC_PAIR)):
    per, tot = [], []
    for sid, S in SAMPLES.items():
        r = rms(model(S, s, A), S["y"])
        per.append(f"{sid.split('_')[1]}={r:.2f}"); tot.append(r)
    print(f"{label:12s} {s:5.2f} {A:7.2f} {np.mean(tot):9.2f}   " + "  ".join(per))

# --- peak-position residuals at the joint pair (medium-index diagnostic) -----
sj, Aj, _ = results["joint"]
print(f"\npeak residuals (model - data, nm) at joint pair (s={sj}, A={Aj}), "
      f"after the -{SHIFT_NM} nm offset:")
for sid, S in SAMPLES.items():
    pk_m = peak_of(WL, model(S, sj, Aj))
    pk_d = peak_of(WL, S["y"])
    print(f"  {sid}: model {pk_m:.0f}  data {pk_d:.0f}  delta {pk_m - pk_d:+.1f}")

# --- figure -------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25, "axes.axisbelow": True})
fig, axes = plt.subplots(2, 3, figsize=(15, 8))
for ax, (sid, S) in zip(axes.ravel()[:4], SAMPLES.items()):
    yd = S["y"][MSK] / S["y"][MSK].max()
    ax.plot(WL[MSK], yd, "k", lw=2.0, label="measured")
    for (s, A), style, lab in (
            (BASE_PAIR, dict(ls=":", color="#888"), "baseline (1, 0.25)"),
            (CTAC_PAIR, dict(ls="--", color="#2CA02C"), "CTAC (0.05, 0.6)"),
            ((sj, Aj), dict(ls="-", color="#D62728"),
             f"seeded ({sj}, {Aj})")):
        ym = model(S, s, A)
        ax.plot(WL[MSK], ym[MSK] / ym[MSK].max(), lw=1.4, **style,
                label=f"{lab}: RMS {rms(ym, S['y']):.1f}%")
    ax.set_title(f"{sid}  (TEM {S['D']:.1f}±{S['sd']:.1f} nm, n={S['n']})",
                 fontsize=9)
    ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("ext / peak")
    ax.legend(fontsize=7, frameon=False)

for ax, name in zip(axes.ravel()[4:], ("joint", "gnp12-only")):
    tot = np.mean([RMS[sid] for sid in VARIANTS[name]], axis=0)
    im = ax.imshow(tot, origin="lower", aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(A_GRID)), [f"{a:g}" for a in A_GRID])
    ax.set_yticks(range(len(S_GRID)), [f"{s:g}" for s in S_GRID])
    ax.set_xlabel("A_surf"); ax.set_ylabel("s (bulk scale)")
    s_b, A_b, r_b = results[name]
    ax.plot(A_GRID.index(A_b), S_GRID.index(s_b), "r*", ms=14)
    ax.set_title(f"mean RMS%, {name} (best {r_b:.2f} @ s={s_b}, A={A_b})",
                 fontsize=9)
    ax.grid(False)
    fig.colorbar(im, ax=ax, shrink=0.85)

fig.suptitle("Seeded-growth route damping calibration on verified-clean monomers "
             f"(per-particle TEM histograms, jc −{SHIFT_NM} nm)", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig20_calibrate_seeded.png"))
print("\nwrote outputs/fig20_calibrate_seeded.png")
