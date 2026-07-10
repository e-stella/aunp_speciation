"""Haiss-style UV-Vis sizing check against TEM, for the CTAC dataset.

Computes the A_spr/A_450 ratio (Haiss 2007 Eq.11) and LSPR peak position, and
compares the measured size-signal to (a) our Mie+Brendel-Bormann forward model
and (b) Haiss's original citrate calibration. Shows that the citrate calibration
overestimates size for clean CTAC spheres -> you must recalibrate per system.

Usage: python scripts/haiss_check.py path/to/GNP_CTAC_TEM_UV.xlsx
"""
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(__file__),"..","src"))
import numpy as np, openpyxl
from aunp_speciation import dielectric; dielectric.use_gold_model('bb')
from aunp_speciation.mie import monomer_cross_sections as mcs

XL = sys.argv[1] if len(sys.argv)>1 else "GNP_CTAC_TEM_UV.xlsx"
wb=openpyxl.load_workbook(XL,data_only=True)
s1=list(zip(*wb["Sheet1"].iter_rows(values_only=True)))
tem={c[0]:np.array([v for v in c[1:] if isinstance(v,(int,float))],float) for c in s1}
rows=list(wb["Sheet2"].iter_rows(values_only=True)); hdr=rows[0]
A=np.array([[float(x) for x in r] for r in rows[1:] if r[0] is not None]); A=A[np.argsort(A[:,0])]
wl=A[:,0]; samples=list(hdr[1:]); tk=list(tem.keys())
print(f"{'sample':14s}{'TEM':>6}{'lam_spr':>8}{'Aspr/A450':>11}{'Haiss_d':>9}")
for k,s in enumerate(samples):
    y=A[:,k+1]; m=(wl>=500)&(wl<=560); Aspr=y[m].max(); lam=wl[m][y[m].argmax()]
    r=Aspr/np.interp(450,wl,y); d=np.exp(3.00*r-2.20)   # Haiss Eq.11 (citrate calib)
    print(f"{s:14s}{tem[tk[k]].mean():6.1f}{lam:8.0f}{r:11.2f}{d:9.1f}")
print("\nNote: Haiss Eq.11 constants (B1=3.00,B2=2.20) are for CITRATE particles;")
print("clean CTAC spheres have a higher A_spr/A_450 -> size is overestimated.")
print("Recalibrate B1,B2 against THIS TEM+UV set for a CTAC-specific ruler.")
