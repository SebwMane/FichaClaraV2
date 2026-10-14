# Control del sesgo en plano: medias por semilla (SE honesto entre semillas) y chequeo de uniformidad
import numpy as np, sys
sys.path.insert(0,".")
from piloto_curv import rng, sub_obs, ref
from omega.causality.sprinkling import sprinkle_diamond, causal_matrix
print("diag", causal_matrix(sprinkle_diamond(5,4,rng(9,9))).diagonal())
m=40; mu,se=ref(m,20000,rng(5,0,m)); print("ref r,ell",mu,se)
for N in (1000,2000):
    per=[]
    for s in range(12):
        r=rng(5,1,N,s); X=sprinkle_diamond(N,4,r); C=causal_matrix(X)
        a,_=sub_obs(C,X,m,r,cap=100000)
        per.append([np.nanmean(a[:,0]),a[:,1].mean()])
    per=np.array(per); d=(per-mu)/mu
    print(N, "dr,dl rel media", d.mean(0), "SE entre semillas", d.std(0,ddof=1)/np.sqrt(len(d)))
