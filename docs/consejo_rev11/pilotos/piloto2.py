"""Piloto 2 declarado: medias por semilla (8 semillas independientes), SE entre semillas.
Geometrias positivas: plano, dS_fuerte, AdS_fuerte, FRW_rad (R=0), FRW_polvo (a ~ eta^2, R>0). Semillas SeedSequence(900101, spawn_key=(2,...))."""
import numpy as np, json, sys, time
sys.path.insert(0,".")
from piloto_curv import rng, sub_obs, ref, GEOS
from f7 import frw_power_points
GEOS=dict(GEOS); GEOS["FRW_polvo"]=lambda n,r: frw_power_points(n,4,8.0,r)
G=["plano","dS_fuerte","AdS_fuerte","FRW_rad","FRW_polvo"]
MS=[10,20,40,80]; NS=[1000,2000,4000]; S=8
refs={m:ref(m,20000,rng(2,0,m)) for m in MS}
out={}; t0=time.time()
for gi,g in enumerate(G):
    for N in NS:
        for s in range(S):
            r=rng(2,1,gi,N,s); X=GEOS[g](N,r)
            from omega.causality.sprinkling import causal_matrix
            C=causal_matrix(X)
            for m in MS:
                a,_=sub_obs(C,X,m,r,cap=1500)
                if len(a)<5: continue
                mu=refs[m][0]
                out.setdefault(f"{g}|{N}|{m}",[]).append([float((np.nanmean(a[:,0])-mu[0])/mu[0]),float((a[:,1].mean()-mu[1])/mu[1]),
                    float(np.corrcoef(a[:,2],a[:,1])[0,1]) if len(a)>10 else float("nan")])
        print(g,N,round(time.time()-t0),flush=True)
json.dump(out,open("piloto2.json","w"))
