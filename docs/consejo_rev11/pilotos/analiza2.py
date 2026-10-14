import json, numpy as np
R=json.load(open("piloto2.json"))
G=["plano","dS_fuerte","AdS_fuerte","FRW_rad","FRW_polvo"]
for j,name in ((1,"dlinks_rel"),(0,"dr_rel (frac. orden)")):
    print("\n==",name,"media(SE entre semillas); filas m, columnas N=1000,2000,4000")
    for g in G:
        print(g)
        for m in (10,20,40,80):
            cells=[]
            for N in (1000,2000,4000):
                v=np.array(R.get(f"{g}|{N}|{m}",[]))
                cells.append(f"{v[:,j].mean():+.4f}({v[:,j].std(ddof=1)/np.sqrt(len(v)):.4f})" if len(v)>2 else "   --   ")
            print(f"  m={m:3d} "+"  ".join(cells))
print("\nAjuste y = b + A x^a, x = m/N (SE entre semillas); a=1/2 fijo y a libre (rejilla, con b)")
for g in G:
  for j,name in ((1,"dl"),(0,"dr")):
    X=[];Y=[];W=[]
    for k,v in R.items():
        gg,N,m=k.split("|")
        if gg!=g: continue
        v=np.array(v)
        if len(v)<3: continue
        X.append(int(m)/int(N)); Y.append(v[:,j].mean()); W.append(1/max(v[:,j].std(ddof=1)**2/len(v),1e-8))
    X,Y,W=map(np.array,(X,Y,W))
    def fit(a):
        A=np.vstack([np.ones_like(X),X**a]).T; sw=np.sqrt(W)[:,None]
        c,*_=np.linalg.lstsq(A*sw,Y*sw[:,0],rcond=None); cov=np.linalg.inv((A*sw).T@(A*sw))
        return c,cov,np.sum(W*(Y-A@c)**2)
    c,cov,chi=fit(0.5)
    best=min((fit(a)[2],a) for a in np.linspace(0.2,1.2,51))
    print(f"{g:10s} {name}: b={c[0]:+.4f}({np.sqrt(cov[0,0]):.4f}) A={c[1]:+.3f}({np.sqrt(cov[1,1]):.3f}) chi2={chi:.1f}/{len(Y)-2}  a_mejor={best[1]:.2f} (chi2 {best[0]:.1f})")
print("\ncorrelacion media (por semilla) entre |dx|/dt del intervalo y links, m=40, N=2000:")
for g in G:
    v=np.array(R[f"{g}|2000|40"]); print(f"  {g:10s} {np.nanmean(v[:,2]):+.3f}({np.nanstd(v[:,2],ddof=1)/np.sqrt(len(v)):.3f})")
