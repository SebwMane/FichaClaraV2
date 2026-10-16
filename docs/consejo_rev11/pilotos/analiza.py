import json, numpy as np
R=json.load(open("piloto_curv.json")); ref=R["ref"]
geos=sorted({k.split("|")[0] for k in R["geo"]})
rows=[]
for k,v in R["geo"].items():
    g,N,m=k.split("|"); N=int(N); m=int(m)
    r=np.concatenate([x["r"] for x in v]); l=np.concatenate([x["l"] for x in v]); vv=np.concatenate([x["v"] for x in v])
    mr,ml=ref[str(m)][0]; sr,sl=ref[str(m)][1]
    dl=(l.mean()-ml)/ml; sdl=np.hypot(l.std()/np.sqrt(len(l)),sl)/ml
    dr=(np.nanmean(r)-mr)/mr; sdr=np.hypot(np.nanstd(r)/np.sqrt(len(r)),sr)/mr
    # dependencia con el boost: diferencia entre tercil rapido y lento
    q=np.quantile(vv,[1/3,2/3]); slow=l[vv<=q[0]].mean(); fast=l[vv>=q[1]].mean()
    rows.append((g,N,m,len(l),int(np.mean([x["n"] for x in v])),dl,sdl,dr,sdr,(fast-slow)/ml))
rows.sort()
print("geo N m nint nintervalos_m  dlinks_rel(se)  dr_rel(se)  boost(rapido-lento)/ref")
for x in rows: print(f"{x[0]:11s} {x[1]:5d} {x[2]:3d} {x[3]:4d} {x[4]:7d}  {x[5]:+.3f}({x[6]:.3f})  {x[7]:+.3f}({x[8]:.3f})  {x[9]:+.3f}")
print("\nAjuste por geometria: y = b + A*(m/N)^a  (a libre y a=1/2 fijo), y = dlinks_rel y dr_rel")
for g in geos:
    sel=[x for x in rows if x[0]==g]
    X=np.array([x[2]/x[1] for x in sel])
    for j,name in ((5,"dlinks"),(7,"dr")):
        y=np.array([x[j] for x in sel]); w=1/np.array([x[j+1] for x in sel])**2
        A=np.vstack([np.ones_like(X),np.sqrt(X)]).T
        W=np.sqrt(w)[:,None]; coef,*_=np.linalg.lstsq(A*W,y*np.sqrt(w),rcond=None)
        cov=np.linalg.inv((A*W).T@(A*W)); chi=np.sum(w*(y-A@coef)**2)
        # exponente libre (sin b) por rejilla
        best=min(((np.sum(w*(y-(np.sum(w*y*X**a)/np.sum(w*X**(2*a)))*X**a)**2),a) for a in np.linspace(0.1,1.5,141)))
        print(f"{g:11s} {name:6s} b={coef[0]:+.4f}({np.sqrt(cov[0,0]):.4f}) A={coef[1]:+.3f}({np.sqrt(cov[1,1]):.3f}) chi2/gl={chi/(len(y)-2):.2f}  a_libre(b=0)={best[1]:.2f}")
