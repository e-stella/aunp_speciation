"""Re-anchor the gold-dielectric calibration on VERIFIED-CLEAN monomers, learn
an empirical basis correction, and test sizing fits on the clean samples.

Background (CLAUDE.md #9 + 2026-07 finding): the production dataset 'jc' was
selected because it put the 12.9 nm model peak on the C500 sample's measured
peak — but C500 is now known to hold 15-30% aggregated gold, which red-shifts
its ensemble peak. Every verified-clean monomer sample since (CTAC series,
2024 seeded-growth series) peaks 1-4 nm BLUE of the TEM-pinned jc model. This
script:

  STEP 1  re-runs the dataset-selection test against NINE verified-clean
          anchors (CTAC 7.8/29.2/31.5/37.9/42.0; seeded Oct GNP12nm, Nov
          seed/GNP23nm/GNP39nm; excluded: Oct 20/40/60 nm — tail excess — and
          Nov GNP55nm — aggregated). Candidates: jc, etchegoin, bb,
          yakubovsky25 (size-matched film for all these R<25 nm), reddy_p200.
          Selection among FIXED published tables only; nothing is tuned.
  STEP 2  learns a multiplicative correction c(lambda) = median over anchors
          of data/(scaled jc model), smoothed; validated LEAVE-ONE-OUT (the
          held-out sample never contributes to its own correction).
  STEP 3  sizing fits (D, poly free; amplitude profiled) on each clean sample
          with the uncorrected and the LOO-corrected basis, vs TEM truth.
          NB the measured Cary noise is wavelength-flat (~1.8 mA), so uniform
          weighting is already statistically correct for these single spectra.

Fit windows: 420-800 nm (CTAC, aqueous) / 480-800 nm (seeded-growth, whose
unpurified matrix absorbs in the blue). gamma_S with A_surf=0.25 throughout.

Run:  python scripts/dielectric_reanchor.py    (system env, Mie only, ~1 min)
Writes outputs/fig18_reanchor.png and prints the verdict tables.
"""
import sys, os, csv
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np

from aunp_speciation import dielectric
from aunp_speciation.spectra import monomer_polydisperse
from aunp_speciation.mie import monomer_cross_sections

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "outputs")
A_SURF = 0.25
WLM = np.arange(420.0, 801.0, 1.0)
CANDIDATES = ("jc", "etchegoin", "bb", "yakubovsky25", "reddy_p200")


# ------------------------------------------------------------- data loading
def load_ctac():
    with open(os.path.join(ROOT, "experimental/ctac/ctac_uvvis.csv")) as f:
        rows = list(csv.reader(f))
    hdr = rows[0]
    dat = np.array([[float(x) for x in r] for r in rows[1:]])
    o = np.argsort(dat[:, 0]); dat = dat[o]
    with open(os.path.join(ROOT, "experimental/ctac/ctac_tem_stats.csv")) as f:
        tem = list(csv.reader(f))[1:]
    out = []
    for i, r in enumerate(tem):
        lab = hdr[i + 1].replace("Abs_GNP_", "CTAC ").replace(" ", "")
        out.append(dict(name=lab, wl=dat[:, 0], y=dat[:, i + 1],
                        D=float(r[1]), poly=float(r[2]), lo=420.0))
    return out


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
    return hdr[1:], arr[o, 0], arr[o, 1:]


def load_seeded():
    E = os.path.join(ROOT, "experimental", "Irina")
    out = []
    hdr, wl, Y = load_utf16_tsv(os.path.join(E, "Scan_10_22_2024_5_43_15_PM.csv"))
    out.append(dict(name="Oct GNP12nm", wl=wl, y=Y[:, hdr.index("GNP12nm")],
                    D=11.55, poly=100 * 0.91 / 11.55, lo=480.0))
    hdr, wl, Y = load_utf16_tsv(os.path.join(E, "Scan_11_16_2024_4_37_22_PM.csv"))
    for col, name, D, sd in [("Seed 1(Abs)", "Nov GNP11nm", 10.85, 0.83),
                             ("30nm 1(Abs)", "Nov GNP23nm", 22.90, 1.46),
                             ("40nm 1(Abs)", "Nov GNP39nm", 38.52, 3.37)]:
        out.append(dict(name=name, wl=wl, y=Y[:, hdr.index(col)],
                        D=D, poly=100 * sd / D, lo=480.0))
    return out


ANCHORS = load_ctac() + load_seeded()


def data_peak(a):
    m = (a["wl"] >= 480) & (a["wl"] <= 650)
    k = np.convolve(a["y"][m], np.ones(5) / 5, "same")
    return a["wl"][m][np.argmax(k)]


# ------------------------------ STEP 1: dataset selection on clean anchors
print("=" * 78)
print("STEP 1 — dataset selection against verified-clean monomer anchors")
print("        (model peak − data peak, nm; gamma_S A_surf=0.25, water)")
print("=" * 78)
models = {}
for cand in CANDIDATES:
    dielectric.use_gold_model(cand)
    models[cand] = [monomer_polydisperse(WLM, a["D"], a["poly"], "water",
                                         None, True, A_SURF) for a in ANCHORS]
dielectric.use_gold_model("jc")

print(f"{'anchor':14s} {'TEM D':>6s} {'pk data':>8s}"
      + "".join(f"{c:>13s}" for c in CANDIDATES))
deltas = {c: [] for c in CANDIDATES}
for i, a in enumerate(ANCHORS):
    pkd = data_peak(a)
    row = f"{a['name']:14s} {a['D']:6.1f} {pkd:8.1f}"
    for c in CANDIDATES:
        pkm = WLM[np.argmax(models[c][i])]
        deltas[c].append(pkm - pkd)
        row += f"{pkm - pkd:+13.1f}"
    print(row)
print(f"{'MEAN ± SD':29s}" + "".join(
    f"{np.mean(deltas[c]):+8.1f}±{np.std(deltas[c]):4.1f}" for c in CANDIDATES))

# ---------------------- STEP 2: empirical multiplicative basis correction
print()
print("=" * 78)
print("STEP 2 — empirical correction c(lambda) from clean anchors (jc basis),")
print("         leave-one-out validated")
print("=" * 78)


def ratio_curve(a, model):
    """data / (LSQ-scaled model) on the anchor's window, on the WLM grid."""
    m = (a["wl"] >= a["lo"]) & (a["wl"] <= 800)
    y = np.interp(WLM, a["wl"][m], a["y"][m], left=np.nan, right=np.nan)
    keep = ~np.isnan(y) & (WLM >= a["lo"])
    s = (model[keep] @ y[keep]) / (model[keep] @ model[keep])
    r = np.full_like(WLM, np.nan)
    r[keep] = y[keep] / (s * model[keep])
    return r


def smooth(x, n=31):
    k = np.ones(n) / n
    pad = np.r_[np.repeat(x[0], n // 2), x, np.repeat(x[-1], n // 2)]
    return np.convolve(pad, k, "valid")


RATIOS = np.array([ratio_curve(a, models["jc"][i]) for i, a in enumerate(ANCHORS)])


def correction_from(idx):
    med = np.nanmedian(RATIOS[idx], axis=0)
    med[np.isnan(med)] = 1.0
    return np.clip(smooth(med), 0.8, 1.2)


print(f"{'held-out anchor':14s} {'pk err befor':>13s} {'after':>7s} "
      f"{'RMS befor':>10s} {'after':>7s}   (RMS % of peak, 480-700 nm)")
loo_corr = {}
for i, a in enumerate(ANCHORS):
    c = correction_from([j for j in range(len(ANCHORS)) if j != i])
    loo_corr[a["name"]] = c
    m0, m1 = models["jc"][i], models["jc"][i] * c
    pkd = data_peak(a)
    band = (WLM >= 480) & (WLM <= 700)
    md = (a["wl"] >= 480) & (a["wl"] <= 700)
    y = np.interp(WLM[band], a["wl"][md], a["y"][md])
    out = []
    for m in (m0, m1):
        s = (m[band] @ y) / (m[band] @ m[band])
        out.append((WLM[np.argmax(m)] - pkd,
                    100 * np.sqrt(np.mean((s * m[band] - y) ** 2)) / y.max()))
    print(f"{a['name']:14s} {out[0][0]:+13.1f} {out[1][0]:+7.1f} "
          f"{out[0][1]:10.2f} {out[1][1]:7.2f}")
CORR_ALL = correction_from(list(range(len(ANCHORS))))

# ------------------------- STEP 3: sizing fits (D, poly free) vs TEM truth
print()
print("=" * 78)
print("STEP 3 — sizing fits (D, poly free; amplitude profiled) vs TEM")
print("=" * 78)
# precompute monomer spectra on a fine D grid once (fast poly integration)
DG = np.arange(5.0, 70.01, 0.25)
dielectric.use_gold_model("jc")
MGRID = np.array([monomer_cross_sections(d, WLM, "water", size_correction=True,
                                         A_surf=A_SURF)["ext"] for d in DG])


def poly_spectrum(D0, p):
    sig = D0 * p / 100.0
    w = np.exp(-0.5 * ((DG - D0) / sig) ** 2)
    w /= w.sum()
    return w @ MGRID


def size_fit(a, corr=None):
    m = (a["wl"] >= a["lo"]) & (a["wl"] <= 800)
    y = np.interp(WLM, a["wl"][m], a["y"][m], left=np.nan, right=np.nan)
    keep = ~np.isnan(y) & (WLM >= a["lo"])
    yk = y[keep]
    best = (1e9, None, None)
    for D0 in np.arange(6.0, 65.01, 0.25):
        for p in np.arange(2.0, 15.01, 0.5):
            mod = poly_spectrum(D0, p)
            if corr is not None:
                mod = mod * corr
            mk = mod[keep]
            s = (mk @ yk) / (mk @ mk)
            r = np.sqrt(np.mean((s * mk - yk) ** 2))
            if r < best[0]:
                best = (r, D0, p)
    return best


print(f"{'sample':14s} {'TEM D':>6s} | {'D uncorr':>9s} {'poly':>5s} | "
      f"{'D corr(LOO)':>11s} {'poly':>5s}")
res_unc, res_cor = [], []
for a in ANCHORS:
    _, D_u, p_u = size_fit(a)
    _, D_c, p_c = size_fit(a, loo_corr[a["name"]])
    res_unc.append(D_u); res_cor.append(D_c)
    print(f"{a['name']:14s} {a['D']:6.1f} | {D_u:9.1f} {p_u:5.1f} | "
          f"{D_c:11.1f} {p_c:5.1f}")

# ------------------------------------------------------------------- figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25, "axes.axisbelow": True})
fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))

ax = axes[0]
xs = np.arange(len(ANCHORS))
mk = dict(jc="o", etchegoin="s", bb="v", yakubovsky25="^", reddy_p200="d")
for c in CANDIDATES:
    ax.plot(xs, deltas[c], mk[c] + "-", ms=4, lw=0.8, label=c)
ax.axhline(0, color="0.3", lw=0.8)
ax.set_xticks(xs)
ax.set_xticklabels([a["name"].replace(" ", "\n") for a in ANCHORS], fontsize=6)
ax.set_ylabel("model peak − data peak (nm)")
ax.set_title("(a) dataset selection on clean anchors", fontsize=9)
ax.legend(frameon=False, fontsize=7.5)

ax = axes[1]
for r in RATIOS:
    ax.plot(WLM, r, color="0.7", lw=0.7)
ax.plot(WLM, CORR_ALL, color="#D55E00", lw=2, label="smoothed median c(λ)")
ax.axhline(1, color="0.3", lw=0.8)
ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("data / scaled jc model")
ax.set_ylim(0.75, 1.6)
ax.set_title("(b) per-anchor mismatch ratios → correction", fontsize=9)
ax.legend(frameon=False, fontsize=7.5)

ax = axes[2]
tem = [a["D"] for a in ANCHORS]
ax.plot([5, 65], [5, 65], "k--", lw=0.8)
ax.plot(tem, res_unc, "o", ms=6, mfc="none", mec="#0072B2",
        label="uncorrected jc basis")
ax.plot(tem, res_cor, "o", ms=5, color="#D55E00",
        label="LOO-corrected basis")
ax.set_xlabel("TEM diameter (nm)"); ax.set_ylabel("fitted diameter (nm)")
ax.set_title("(c) shape-only sizing fit vs TEM", fontsize=9)
ax.legend(frameon=False, fontsize=7.5)

fig.suptitle("Dielectric re-anchoring on verified-clean monomers + empirical "
             "basis correction (jc, γ_S A=0.25; LOO-validated)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig18_reanchor.png"))
print("\nwrote outputs/fig18_reanchor.png")
