"""Speciation-fit VALIDATION on the CTAC size series (5 samples, paired TEM).

The cleanest data available: 5 monomer samples with paired per-particle TEM,
no temperature ramp, no concentration drift, no evaporation. D and
polydispersity are PINNED to TEM (nothing size-like is fitted), A_surf=0.25,
exact T-matrix optics — so the ONLY question the fit answers is the species
split. Two things are reported, in order:

1. IDENTIFIABILITY DIAGNOSTIC (report FIRST): cosine similarity between the
   per-gold-atom basis vectors (monomer vs dimer/2 vs trimer/3) over the fit
   window. On CDA at 12.9 nm, gap 1, A_surf 0.25 this is 0.9966 — collinear:
   NNLS returns ~100% monomer and the red tail goes unexplained. If the EXACT
   basis is not meaningfully below ~0.99, the dimer fraction is NOT
   identifiable at that size and every speciation percentage must carry that
   caveat.
2. %monomer / %dimer / %trimer (gold fractions) per sample, from NNLS with
   the pinned exact basis; with and without the scattering-pedestal term.
   Fraction SDs from residual bootstrap at the best pedestal exponent.

GAP CHOICE (matters): fits use gap = 3.5 nm — the CTAC bilayer is ~3.2-4 nm,
so smaller gaps are sterically impossible for this ligand; conveniently the
T-matrix is also converged there (lmax=8 changes <0.6% vs lmax=10 at D=47).
Gap 1 nm at D>=29 is NOT converged even at lmax=10 (gap/D~0.02) AND is
unphysical for CTAC. The 12.9 nm/gap-1 reference IS converged at lmax=8 and
is computed for apples-to-apples comparison with the CDA 0.9966.

Run:  mstm-env/bin/python scripts/fit_ctac_validation.py
(~20 min on first run; the exact basis is cached in
experimental/ctac/ctac_exact_basis.npz and reused afterwards)

Inputs (exported from experimental/GNP_CTAC_TEM_UV.xlsx):
  experimental/ctac/ctac_uvvis.csv      - 300-900 nm, 1 nm, 5 samples
  experimental/ctac/ctac_tem_stats.csv  - per-sample TEM mean D, % poly
"""
import sys, os, csv, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import nnls

from aunp_speciation import dielectric
dielectric.use_gold_model("jc")
from aunp_speciation.spectra import monomer_polydisperse, gaussian_sizes
from aunp_speciation.clusters import species_spectrum as species_spectrum_cda

ROOT = os.path.join(os.path.dirname(__file__), "..")
CTAC = os.path.join(ROOT, "experimental", "ctac")
OUT = os.path.join(ROOT, "outputs")
CACHE = os.path.join(CTAC, "ctac_exact_basis.npz")

GAP = 3.5            # nm — CTAC bilayer (see module docstring)
GAP_REF = 1.0        # nm — 12.9 nm reference diagnostic (CDA comparison)
LMAX = 8
A_SURF = 0.25
N_SIZES = 5          # polydispersity samples for the cluster integral
WINDOW = (420.0, 800.0)
WL_TM = np.arange(420.0, 801.0, 10.0)   # treams grid; cubic-interp to 1 nm
CLUSTERS = ("dimer", "trimer_linear")
K_S = {"monomer": 1, "dimer": 2, "trimer_linear": 3}
N_BOOT = 200


def load_csv(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    return rows[0], np.array([[float(x) for x in r] for r in rows[1:]])


def cosine(a, b):
    return float(a @ b / np.sqrt((a @ a) * (b @ b)))


def build_exact_basis(D_tem, poly_pct):
    """Exact per-cluster spectra on WL_TM for each sample (treams; slow)."""
    from aunp_speciation.clusters_tmatrix import species_spectrum_tmatrix, available
    if not available():
        raise RuntimeError("needs treams — run with mstm-env/bin/python")
    cubes = {sp: np.zeros((len(D_tem), len(WL_TM))) for sp in CLUSTERS}
    for i, (D, poly) in enumerate(zip(D_tem, poly_pct)):
        ds, w = gaussian_sizes(D, poly, n=N_SIZES)
        for sp in CLUSTERS:
            t0 = time.time()
            acc = np.zeros_like(WL_TM)
            for d, wi in zip(ds, w):
                acc += wi * species_spectrum_tmatrix(
                    sp, d, WL_TM, GAP, "water", LMAX,
                    size_correction=True, A_surf=A_SURF)
            cubes[sp][i] = acc
            print(f"  built {sp:14s} D={D:5.2f} poly={poly:.1f}% "
                  f"({time.time()-t0:5.1f} s)", flush=True)
    # 12.9 nm / gap 1 reference dimer (single size, converged at lmax=8)
    ref = species_spectrum_tmatrix("dimer", 12.9, WL_TM, GAP_REF, "water",
                                   LMAX, size_correction=True, A_surf=A_SURF)
    ref35 = species_spectrum_tmatrix("dimer", 12.9, WL_TM, GAP, "water",
                                     LMAX, size_correction=True, A_surf=A_SURF)
    return cubes, ref, ref35


def get_basis(D_tem, poly_pct):
    meta = dict(gap=GAP, gap_ref=GAP_REF, lmax=LMAX, A_surf=A_SURF,
                n_sizes=N_SIZES, gold_model="jc")
    if os.path.exists(CACHE):
        z = np.load(CACHE)
        ok = (np.allclose(z["D_tem"], D_tem) and np.allclose(z["poly"], poly_pct)
              and np.allclose(z["wl"], WL_TM)
              and all(float(z[k]) == float(v) for k, v in meta.items()
                      if k != "gold_model")
              and str(z["gold_model"]) == meta["gold_model"])
        if ok:
            print(f"using cached exact basis {CACHE}")
            return ({sp: z[f"cube__{sp}"] for sp in CLUSTERS},
                    z["ref_dimer_129_gap1"], z["ref_dimer_129_gap35"])
        print("cache stale (parameters changed) — rebuilding")
    cubes, ref, ref35 = build_exact_basis(D_tem, poly_pct)
    np.savez_compressed(
        CACHE, D_tem=D_tem, poly=poly_pct, wl=WL_TM,
        ref_dimer_129_gap1=ref, ref_dimer_129_gap35=ref35,
        gold_model=np.array("jc"),
        **{k: np.array(v) for k, v in meta.items() if k != "gold_model"},
        **{f"cube__{sp}": c for sp, c in cubes.items()})
    return cubes, ref, ref35


def fit_sample(y, B):
    """NNLS species-only and species+pedestal (exponent grid-scanned)."""
    w0, _ = nnls(B, y)
    rms0 = float(np.sqrt(np.mean((B @ w0 - y) ** 2)))
    best = None
    for n in np.arange(0.0, 6.01, 0.25):
        Bp = np.column_stack([B, (wl_fit / 550.0) ** (-n)])
        w, _ = nnls(Bp, y)
        rms = float(np.sqrt(np.mean((Bp @ w - y) ** 2)))
        if best is None or rms < best[0]:
            best = (rms, n, w, Bp)
    return (w0, rms0), best


def gold_fractions(w):
    kw = np.array([K_S[sp] * wi for sp, wi in zip(SPECIES, w[:len(SPECIES)])])
    return kw / kw.sum() if kw.sum() > 0 else kw


# ---------------------------------------------------------------- load data
hdr_uv, uv = load_csv(os.path.join(CTAC, "ctac_uvvis.csv"))
with open(os.path.join(CTAC, "ctac_tem_stats.csv")) as f:
    tem_rows = list(csv.reader(f))[1:]          # sample,D_tem,poly_pct,n
labels = [h.replace("Abs_GNP_", "").replace(" ", "") for h in hdr_uv[1:]]
D_tem = np.array([float(r[1]) for r in tem_rows])
poly_pct = np.array([float(r[2]) for r in tem_rows])
wl_raw, A_raw = uv[:, 0], uv[:, 1:]
order = np.argsort(wl_raw)
wl_raw, A_raw = wl_raw[order], A_raw[order]
m = (wl_raw >= WINDOW[0]) & (wl_raw <= WINDOW[1])
wl_fit, A_fit = wl_raw[m], A_raw[m]
SPECIES = ("monomer",) + CLUSTERS

print(f"samples: {labels}")
print(f"pinned:  D = {np.round(D_tem, 2)}  poly% = {np.round(poly_pct, 2)}")
print(f"config:  gap={GAP} nm (CTAC bilayer), lmax={LMAX}, A_surf={A_SURF}, "
      f"jc, window {WINDOW[0]:.0f}-{WINDOW[1]:.0f} nm\n")

cubes, ref_dimer_g1, ref_dimer_g35 = get_basis(D_tem, poly_pct)

# per-sample basis on the 1 nm data grid: exact monomer (live Mie) + splined
# exact clusters; per-cluster columns (fraction math divides by k later)
bases = []
for i in range(len(labels)):
    cols = [monomer_polydisperse(wl_fit, D_tem[i], poly_pct[i], "water",
                                 None, True, A_SURF)]
    for sp in CLUSTERS:
        cols.append(CubicSpline(WL_TM, cubes[sp][i])(wl_fit))
    bases.append(np.column_stack(cols))

# ------------------------------------------------- DIAGNOSTIC 1: collinearity
print("=" * 76)
print("DIAGNOSTIC 1 — basis collinearity: cosine of per-gold-atom vectors")
print("=" * 76)
mono_129 = monomer_polydisperse(wl_fit, 12.9, 0.0, "water", None, True, A_SURF)
cda_129 = species_spectrum_cda("dimer", 12.9, wl_fit, GAP_REF, "water",
                               size_correction=True, A_surf=A_SURF) / 2
exact_129_g1 = CubicSpline(WL_TM, ref_dimer_g1)(wl_fit) / 2
exact_129_g35 = CubicSpline(WL_TM, ref_dimer_g35)(wl_fit) / 2
print(f"reference 12.9 nm (monodisperse), mono vs dimer:")
print(f"  gap 1.0 nm:  CDA {cosine(mono_129, cda_129):.4f}   "
      f"EXACT {cosine(mono_129, exact_129_g1):.4f}")
print(f"  gap 3.5 nm (CTAC-realistic):    EXACT {cosine(mono_129, exact_129_g35):.4f}")
print(f"\nper sample (EXACT basis, gap {GAP} nm, pinned D/poly):")
print(f"  {'sample':8s} {'D_tem':>6s} {'poly%':>6s} {'cos(m,d)':>9s} "
      f"{'cos(m,t)':>9s} {'cos(d,t)':>9s} {'orth(m,d)%':>11s} {'cond':>7s}")
coss = []
for i, lab in enumerate(labels):
    b = bases[i] / np.array([1.0, 2.0, 3.0])          # per gold atom
    c_md = cosine(b[:, 0], b[:, 1])
    c_mt = cosine(b[:, 0], b[:, 2])
    c_dt = cosine(b[:, 1], b[:, 2])
    cond = np.linalg.cond(b / np.linalg.norm(b, axis=0))
    coss.append(c_md)
    print(f"  {lab:8s} {D_tem[i]:6.2f} {poly_pct[i]:6.2f} {c_md:9.4f} "
          f"{c_mt:9.4f} {c_dt:9.4f} {100*np.sqrt(1-c_md**2):11.2f} {cond:7.1f}")
if max(coss) > 0.99:
    print("\n  *** WARNING: cos(mono,dimer) > 0.99 for at least one sample —")
    print("  *** dimer fraction NOT identifiable there; treat the percentages")
    print("  *** below as upper-bound decompositions, not measurements.")

# ------------------------------------------------------------------- fits
print()
print("=" * 76)
print("SPECIATION FITS — gold fractions, D/poly pinned to TEM")
print("=" * 76)
rng = np.random.default_rng(0)
results = []
for i, lab in enumerate(labels):
    y, B = A_fit[:, i], bases[i]
    (w0, rms0), (rms, n_sca, w, Bp) = fit_sample(y, B)
    g0, g = gold_fractions(w0), gold_fractions(w)
    # residual bootstrap at fixed best pedestal exponent
    resid = y - Bp @ w
    boots = np.empty((N_BOOT, len(SPECIES)))
    for k in range(N_BOOT):
        wb, _ = nnls(Bp, Bp @ w + rng.choice(resid, len(resid), replace=True))
        boots[k] = gold_fractions(wb)
    gsd = boots.std(axis=0)
    ped_share = (w[-1] * (wl_fit / 550.0) ** (-n_sca)).max() / y.max()
    results.append(dict(lab=lab, g=g, gsd=gsd, g0=g0, n_sca=n_sca, rms=rms,
                        rms0=rms0, w=w, Bp=Bp, ped=ped_share))
    print(f"{lab:8s} mono {100*g[0]:5.1f}±{100*gsd[0]:4.1f}%  "
          f"dimer {100*g[1]:5.1f}±{100*gsd[1]:4.1f}%  "
          f"trimer {100*g[2]:5.1f}±{100*gsd[2]:4.1f}%  "
          f"| pedestal n={n_sca:.2f} share {100*ped_share:4.1f}%  "
          f"RMS {rms:.5f}")
    print(f"{'':8s} species-only: mono/dimer/trimer "
          f"{100*g0[0]:.1f}/{100*g0[1]:.1f}/{100*g0[2]:.1f}%  RMS {rms0:.5f}")

# ------------------------------------------------------------------ figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25})
fig, axes = plt.subplots(2, 3, figsize=(13, 7))
for i, r in enumerate(results):
    ax = axes.flat[i]
    y = A_fit[:, i]
    ax.plot(wl_fit, y, "k.", ms=2, label="data")
    ax.plot(wl_fit, r["Bp"] @ r["w"], "-", color="#D55E00", lw=1.5, label="fit")
    for j, (sp, c) in enumerate(zip(SPECIES, ("#0072B2", "#009E73", "#CC79A7"))):
        if r["w"][j] > 0:
            ax.plot(wl_fit, r["w"][j] * r["Bp"][:, j], "--", color=c, lw=1,
                    label=sp)
    if r["w"][-1] > 0:
        ax.plot(wl_fit, r["w"][-1] * r["Bp"][:, -1], ":", color="0.5", lw=1,
                label=f"pedestal n={r['n_sca']:.1f}")
    ax.set_title(f"{r['lab']}: D={D_tem[i]:.1f} poly={poly_pct[i]:.1f}% | "
                 f"m/d/t {100*r['g'][0]:.0f}/{100*r['g'][1]:.0f}/"
                 f"{100*r['g'][2]:.0f}% | cos {coss[i]:.4f}", fontsize=8)
    ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("A")
    ax.legend(frameon=False, fontsize=7)
ax = axes.flat[5]
ax.bar(range(len(labels)), coss, color="#0072B2")
ax.axhline(0.9966, color="#D55E00", ls="--", lw=1, label="CDA 12.9 nm (0.9966)")
ax.axhline(0.99, color="0.3", ls=":", lw=1, label="identifiability ~0.99")
ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylim(min(0.97, min(coss) - 0.005), 1.0)
ax.set_ylabel("cos(monomer, dimer) per gold")
ax.set_title("basis collinearity (exact, gap 3.5 nm)", fontsize=8)
ax.legend(frameon=False, fontsize=7)
fig.suptitle("CTAC validation: exact T-matrix speciation fit, D/poly pinned "
             f"to TEM, A_surf={A_SURF}, gap={GAP} nm, jc", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig12_ctac_speciation.png"))
print(f"\nwrote outputs/fig12_ctac_speciation.png")
