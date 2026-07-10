"""Validate the optics on real data: simulate each sample from its MEASURED
TEM size distribution (no free size parameter) and overlay on the experimental
UV-Vis. Uses the accurate Brendel-Bormann gold dielectric. Writes fig9.

Usage: python scripts/fit_ctac_realdata.py path/to/GNP_CTAC_TEM_UV.xlsx
"""
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(__file__),"..","src"))
import numpy as np, openpyxl, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from aunp_speciation import dielectric; dielectric.use_gold_model('bb')
from aunp_speciation.mie import monomer_cross_sections as mcs

XL=sys.argv[1] if len(sys.argv)>1 else "GNP_CTAC_TEM_UV.xlsx"  # xlsx: Sheet1=TEM diam cols, Sheet2=WV_nm+Abs cols
wb=openpyxl.load_workbook(XL,data_only=True)
s1=list(zip(*wb["Sheet1"].iter_rows(values_only=True)))
tem={c[0]:np.array([v for v in c[1:] if isinstance(v,(int,float))],float) for c in s1}
rows=list(wb["Sheet2"].iter_rows(values_only=True)); hdr=rows[0]
A=np.array([[float(x) for x in r] for r in rows[1:] if r[0] is not None]); A=A[np.argsort(A[:,0])]
wl=A[:,0]; samples=list(hdr[1:]); tkeys=list(tem.keys())

def spectrum_from_hist(diams,wlw):
    # coarse-grain the empirical diameters to speed Mie (bin to 1 nm)
    d=np.clip(diams,2,80); edges=np.arange(d.min()-0.5,d.max()+1.5,1.0)
    cnt,_=np.histogram(d,edges); cen=0.5*(edges[:-1]+edges[1:]); cen=cen[cnt>0]; cnt=cnt[cnt>0]
    w=cnt/cnt.sum()
    return sum(wi*mcs(dc,wlw,"water",size_correction=False)["ext"] for dc,wi in zip(cen,w))

WIN=(460,740); fig,axes=plt.subplots(2,3,figsize=(14,8)); peaks=[]
for k,(ax,scol,tk) in enumerate(zip(axes.ravel(),samples,tkeys)):
    y=A[:,k+1]; m=(wl>=WIN[0])&(wl<=WIN[1]); wlw=wl[m]; yw=y[m]
    Dt=tem[tk].mean(); pt=100*tem[tk].std(ddof=1)/Dt
    M=spectrum_from_hist(tem[tk],wlw); G=np.column_stack([M,np.ones_like(M)])
    c=np.linalg.lstsq(G,yw,rcond=None)[0]; model=G@c
    rms=np.sqrt(np.mean((model-yw)**2))/yw.max()
    ep=wlw[yw.argmax()]; mp=wlw[model.argmax()]; peaks.append((Dt,ep,mp))
    ax.plot(wlw,yw/yw.max(),'k.',ms=3,label='experiment')
    ax.plot(wlw,model/yw.max(),'-',color="#D55E00",lw=2,label='simulated (TEM dist.)')
    ax.set_title(f'{scol}: TEM {Dt:.0f} nm ({pt:.0f}%)\nexp peak {ep:.0f} / sim {mp:.0f} nm, rmsN {rms:.3f}',fontsize=9)
    ax.set_xlabel('wavelength (nm)'); ax.set_ylim(0,1.12)
    if k==0: ax.legend(fontsize=8,frameon=False)
axl=axes.ravel()[-1]
Dt=[p[0] for p in peaks]; ep=[p[1] for p in peaks]; mp=[p[2] for p in peaks]
axl.plot(Dt,ep,'o-',color='k',label='experiment')
axl.plot(Dt,mp,'s--',color="#D55E00",label='simulated (BB)')
axl.set_xlabel('TEM diameter (nm)'); axl.set_ylabel('LSPR peak (nm)')
axl.set_title('peak position vs size'); axl.legend(fontsize=8)
fig.suptitle('CTAC gold: experimental UV-Vis vs simulation from the measured TEM size distribution\n(Mie + Brendel-Bormann gold, no free size parameter)',fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(os.path.dirname(__file__),'..','outputs','fig9_ctac_realdata.png'),dpi=120)
print(f"{'sample':16s}{'TEM_nm':>8}{'exp_pk':>8}{'sim_pk':>8}{'rmsN':>7}")
for (Dt_,ep_,mp_),s in zip(peaks,samples):
    pass
for k,s in enumerate(samples):
    print(f"{s:16s}{peaks[k][0]:8.1f}{peaks[k][1]:8.0f}{peaks[k][2]:8.0f}")
print("wrote outputs/fig9_ctac_realdata.png")
