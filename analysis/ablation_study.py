# -*- coding: utf-8 -*-
"""ablation_study.py — remove each physics term, RE-CALIBRATE the rest fairly,
and measure what breaks: (i) Table-1 shape fit, (ii) zero-refit dyad
validation at the confirmed ~72 h test duration. Also builds the SI
response-surface / interaction figure from the exhaustive grid.
Ablations: no_closure (f0=0) · no_floor (f_res=0) · no_crack (f_c0=0) ·
no_tortuous (tau^2=1) · no_lag (steady-state only) · no_calibration
(generic 1 pinhole/mm^2, no fit at all)."""
import numpy as np, csv, os, sys
from scipy.optimize import least_squares
sys.path.insert(0, "/home/claude")
from tfe_physics_engine import ImplicitDiffusion1D

OUT = "/mnt/user-data/outputs/ablation"; os.makedirs(OUT, exist_ok=True)
RG, TREF, T = 8.314, 298.15, 311.15
DA = 1.0
G2 = lambda J: J*1e3*86400; K2 = lambda W: W/1e3/86400
arr = lambda x, Ea: x*np.exp(-Ea/RG*(1/T-1/TREF))
d_tab = np.array([15.,20.,30.,50.,60.]); W_tab = np.array([6.7e-3,7e-4,8e-4,1.3e-3,4.7e-3])
MEAS = {1:1.7e-4, 2:3.6e-5, 3:7.7e-6}
D_par,S_par,Ea_par = 5e-13,1.5,40e3; P_par = arr(D_par,Ea_par)*S_par
P_lat = arr(1e-21,60e3)*0.10; r_pin=50e-9
R_PET = DA/K2(3.0); L_PET,S_PET = 125e-6,5.0; P_PET = L_PET/R_PET
DIN, DORG = 50., 500.  # corrected anchor

SHAPES = {  # name -> (shape(d, p), n_params, p0, bounds)
 "baseline":     (lambda d,p: 10**p[0]*np.exp(-d/p[1]) + 1 + 10**p[2]*np.maximum(0,(d-p[3])/p[3])**2,
                  4, [np.log10(300),4,np.log10(3),40], ([0,1,-2,20],[7,15,4,60])),
 "no_closure":   (lambda d,p: 1 + 10**p[0]*np.maximum(0,(d-p[1])/p[1])**2,
                  2, [np.log10(3),40], ([-2,20],[4,60])),
 "no_floor":     (lambda d,p: 10**p[0]*np.exp(-d/p[1]) + 10**p[2]*np.maximum(0,(d-p[3])/p[3])**2,
                  4, [np.log10(300),4,np.log10(3),40], ([0,1,-2,20],[7,15,4,60])),
 "no_crack":     (lambda d,p: 10**p[0]*np.exp(-d/p[1]) + 1,
                  2, [np.log10(300),4], ([0,1],[7,15])),
 "no_tortuous":  None, "no_lag": None,   # share baseline shape
 "no_calibration": None,
}

def fit_shape(name):
    shape, npar, p0, bnd = SHAPES[name]
    def res(q):
        lK, p = q[0], q[1:]
        return np.log10(np.maximum(10**lK*shape(d_tab,p),1e-30)) - np.log10(W_tab)
    f = least_squares(res, [np.log10(7e-4)]+p0, bounds=([-6]+list(bnd[0]),[0]+list(bnd[1])))
    return (lambda d: shape(np.asarray(d,float), f.x[1:])), float(np.max(np.abs(res(f.x))))

def forward(ffun, tortuous=True):
    f30 = None
    # Stage B: anchor absolute scale on dyad-1 (same protocol as paper)
    R_top = DORG*1e-9/P_par
    R1 = 1.0/K2(MEAS[1])
    f30 = (DIN*1e-9/(R1-R_PET-R_top) - P_lat)/P_par
    scale = f30/ffun(DIN)
    f = lambda d: scale*ffun(d)
    def tau2(fv):
        if not tortuous: return 1.0
        s = r_pin*np.sqrt(np.pi/fv)
        return 1 + s**2*np.log(max(s/r_pin,np.e))/(2*np.pi*(DORG*1e-9)**2)
    def stack(n):
        fv=f(DIN); P_in=P_lat+fv*P_par; t2=tau2(fv)
        dx=[L_PET/24]*24; K=[P_PET]*24; C=[S_PET]*24
        for i in range(n):
            dx+=[DIN*1e-9/8]*8; K+=[P_in]*8; C+=[0.10]*8
            tt=1.0 if i==n-1 else t2
            dx+=[DORG*1e-9/8]*8; K+=[P_par/tt]*8; C+=[S_par]*8
        return map(np.array,(dx,K,C))
    def steady(n):
        fv=f(DIN); R = R_PET + n*DIN*1e-9/(P_lat+fv*P_par) + (n-1)*DORG*1e-9*tau2(fv)/P_par + R_top
        return G2(DA/R)
    def app72(n):
        dx,K,C = stack(n); sol=ImplicitDiffusion1D(dx)
        a=np.zeros(len(dx)); t,dt=0.,1.; ts,Js=[],[]
        while t < 10*86400:
            dtn=min(dt,10*86400-t)
            a,Jl,_=sol.step(a,dtn,K,C,("dirichlet",0.),("dirichlet",DA))
            t+=dtn; ts.append(t); Js.append(Jl); dt=min(dt*1.2,900.)
        return G2(np.interp(72*3600,ts,Js))
    return steady, app72, f30

rows=[]
def record(name, maxres, e2, e3, note):
    rows.append(dict(ablation=name, table1_maxres_dec=round(maxres,3),
                     dyad2_err_dec=round(e2,2), dyad3_err_dec=round(e3,2), note=note))
    print(f"{name:16s} fit={maxres:5.3f}  d2={e2:+5.2f}  d3={e3:+5.2f}  {note}")

for name in ("baseline","no_closure","no_floor","no_crack"):
    ffun, mres = fit_shape(name)
    st, ap, _ = forward(ffun, tortuous=True)
    e2, e3 = np.log10(ap(2)/MEAS[2]), np.log10(ap(3)/MEAS[3])
    record(name, mres, e2, e3, "refit + 72h transient")

fb,_ = fit_shape("baseline"); mres_b = fit_shape("baseline")[1]
st, ap, _ = forward(fb, tortuous=False)
record("no_tortuous", mres_b, np.log10(ap(2)/MEAS[2]), np.log10(ap(3)/MEAS[3]),
       "tau^2=1; shape refit unaffected")
st, ap, _ = forward(fb, tortuous=True)
record("no_lag", mres_b, np.log10(st(2)/MEAS[2]), np.log10(st(3)/MEAS[3]),
       "steady-state prediction only")
# uncalibrated: generic 1 pinhole/mm2, no fitting anywhere
fgen = lambda d: np.pi*r_pin**2*1e6*np.ones_like(np.asarray(d,float))
def fwd_nocal():
    fv = float(fgen(DIN)); P_in=P_lat+fv*P_par
    s=r_pin*np.sqrt(np.pi/fv); t2=1+s**2*np.log(s/r_pin)/(2*np.pi*(DORG*1e-9)**2)
    def steady(n):
        R=R_PET+n*DIN*1e-9/P_in+(n-1)*DORG*1e-9*t2/P_par+DORG*1e-9/P_par
        return G2(DA/R)
    return steady
stg = fwd_nocal()
e1 = np.log10(stg(1)/MEAS[1])
record("no_calibration", float("nan"), np.log10(stg(2)/MEAS[2]), np.log10(stg(3)/MEAS[3]),
       f"generic 1/mm2, no fit; dyad-1 err {e1:+.2f} dec")

with open(f"{OUT}/ablation_table.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

# ---------------- figure: ablation bars + response surface + interactions
import pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,
                     "legend.frameon":False,"savefig.bbox":"tight"})
fig,ax=plt.subplots(1,3,figsize=(11.5,3.2))
names=[r["ablation"] for r in rows]; x=np.arange(len(names))
a2=[abs(r["dyad2_err_dec"]) for r in rows]; a3=[abs(r["dyad3_err_dec"]) for r in rows]
ax[0].bar(x-0.18,a2,0.36,label="|err| 2-dyad",color="#0072B2")
ax[0].bar(x+0.18,a3,0.36,label="|err| 3-dyad",color="#D55E00")
ax[0].axhline(0.4,ls="--",c="k",lw=.8); ax[0].text(0.05,0.42,"0.4 dec",fontsize=7)
ax[0].set_xticks(x); ax[0].set_xticklabels([n.replace("_","\n") for n in names],fontsize=7)
ax[0].set(ylabel="zero-refit error [decades]",title="(a) Ablation: validation degradation")
ax[0].legend(fontsize=7)
df=pd.read_csv("/mnt/user-data/outputs/step6_optimization/grid_all.csv")
sl=df[df.d_org==100].pivot_table(index="n",columns="d_in",values="T80_h")
im=ax[1].imshow(np.log10(sl.values),aspect="auto",origin="lower",cmap="viridis",
   extent=[sl.columns.min(),sl.columns.max(),sl.index.min()-0.5,sl.index.max()+0.5])
ax[1].axvline(22.5,ls="--",c="w",lw=1); ax[1].axvline(44,ls="--",c="r",lw=1)
plt.colorbar(im,ax=ax[1],label="log10 T80 [h]")
ax[1].set(xlabel="d_inorg [nm]",ylabel="n_pairs",title="(b) Response surface (d_org=100 nm)")
for n in range(1,7):
    s2=df[(df.d_org==100)&(df.n==n)].sort_values("d_in")
    ax[2].semilogy(s2.d_in,s2.T80_h,label=f"n={n}",color=plt.cm.viridis((n-1)/5))
ax[2].axvspan(22.5,44,color="grey",alpha=.12)
ax[2].set(xlabel="d_inorg [nm]",ylabel="T80 [h]",title="(c) Interaction: d_in effect scales with n")
ax[2].legend(fontsize=6.5,ncol=2)
for ext,kw in (("png",dict(dpi=300)),("pdf",{})):
    fig.savefig(f"{OUT}/FigS4_ablation_sensitivity.{ext}",**kw)
print("saved:", OUT)
