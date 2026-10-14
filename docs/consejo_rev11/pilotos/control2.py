import numpy as np, sys
sys.path.insert(0,".")
from piloto_curv import rng, sub_obs, ref
from omega.causality.sprinkling import sprinkle_diamond, causal_matrix
m=10
for rs in range(2):
    mu,se=ref(m,40000,rng(6,rs,m)); print("ref",mu,se)
per=[]
for s in range(40):
    r=rng(6,2,s); X=sprinkle_diamond(1000,4,r); C=causal_matrix(X)
    a,_=sub_obs(C,X,m,r,cap=100000); per.append([a[:,0].mean(),a[:,1].mean()])
per=np.array(per); print("muestreo plano N=1000 m=10:",per.mean(0),per.std(0,ddof=1)/np.sqrt(40))
