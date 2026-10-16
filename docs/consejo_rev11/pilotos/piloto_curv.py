"""Piloto declarado (asesor curvatura). Solo geometrias positivas. SeedSequence(900101).
Observables por sub-intervalo [x,y] con exactamente m elementos interiores (m+2 con las puntas):
  r   = fraccion de pares relacionados entre los m interiores (Myrheim-Meyer)
  ell = links del sub-orden (m+2 elementos) por elemento
Desviacion = media(geom) - media(referencia exacta: m uniformes en diamante plano 4D + puntas).
"""
import sys, json, time
import numpy as np
sys.path.insert(0, "/home/user/FichaClaraV2"); sys.path.insert(0, "/home/user/FichaClaraV2/tools")
from omega.causality.sprinkling import sprinkle_diamond, causal_matrix
from f7 import de_sitter_eta_points, frw_power_points

SS = np.random.SeedSequence(900101)
def rng(*k): return np.random.Generator(np.random.PCG64(np.random.SeedSequence(900101, spawn_key=k)))

def ads_points(n, r, z0=0.7, zmin=0.2):
    # AdS Poincare conforme: ds^2=(1/z^2)(-dt^2+dz^2+dx^2), z = x1 + z0 in [0.2,1.2]; densidad prop. z^-4
    pts=[]; tot=0
    while tot<n:
        c=sprinkle_diamond(4*n,4,r); z=c[:,1]+z0
        k=c[r.random(len(c))<(zmin/z)**4]; pts.append(k); tot+=len(k)
    return np.concatenate(pts)[:n]

GEOS = {
 "plano": lambda n,r: sprinkle_diamond(n,4,r),
 "dS": lambda n,r: de_sitter_eta_points(n,4,r),
 "dS_fuerte": lambda n,r: de_sitter_eta_points(n,4,r,-1.2,-0.2),
 "AdS_fuerte": ads_points,
 "FRW_rad": lambda n,r: frw_power_points(n,4,4.0,r),
}

def obs(sub):  # sub: bool matrix (m+2)x(m+2) in causal order (transitively closed)
    m2=sub.shape[0]; m=m2-2
    inner=sub[1:-1,1:-1]
    r=inner.sum()/(m*(m-1)/2) if m>1 else np.nan
    s=sub.astype(np.float32); two=(s@s)>0
    links=(sub & ~two).sum()
    return r, links/m2

def sub_obs(C, coords, m, r, cap=400):
    Cf=C.astype(np.float32); B=np.rint(Cf@Cf).astype(np.int32); B[~C]=-1
    I,J=np.nonzero(B==m)
    if len(I)>cap:
        sel=r.choice(len(I),cap,replace=False); I,J=I[sel],J[sel]
    out=[]
    for i,j in zip(I,J):
        mid=np.nonzero(C[i]&C[:,j])[0]
        idx=np.concatenate([[i],mid,[j]])
        o=np.argsort(coords[idx,0]); idx=idx[o]
        rr,ll=obs(C[np.ix_(idx,idx)])
        dt=coords[j,0]-coords[i,0]; dx=np.linalg.norm(coords[j,1:]-coords[i,1:])
        out.append((rr,ll,dx/dt))
    return np.array(out), len(np.nonzero(B==m)[0])

def ref(m, reps, r):
    tips=np.array([[-0.5,0,0,0],[0.5,0,0,0]])
    v=[]
    for _ in range(reps):
        p=np.vstack([tips[:1],sprinkle_diamond(m,4,r),tips[1:]])
        p=p[np.argsort(p[:,0])]
        v.append(obs(causal_matrix(p)))
    v=np.array(v); return v.mean(0), v.std(0)/np.sqrt(reps)

if __name__=="__main__":
    MS=[5,10,20,40,80]; res={"ref":{}, "geo":{}}
    t0=time.time()
    for m in MS:
        mu,se=ref(m,4000,rng(0,m)); res["ref"][m]=(mu.tolist(),se.tolist())
    print("ref", round(time.time()-t0,1), res["ref"], flush=True)
    NS=[1000,2000,4000,8000]
    for gi,(g,f) in enumerate(GEOS.items()):
        for N in NS:
            if g in ("dS","FRW_rad","AdS_fuerte") and N==8000: continue
            for s in range(3 if N<8000 else 1):
                r=rng(1,gi,N,s); X=f(N,r); C=causal_matrix(X)
                for m in MS:
                    a,cnt=sub_obs(C,X,m,r)
                    if len(a)==0: continue
                    mu=np.array(res["ref"][m][0])
                    res["geo"].setdefault(f"{g}|{N}|{m}",[]).append(dict(
                        n=cnt, dr=float(np.nanmean(a[:,0])-mu[0]), dl=float(a[:,1].mean()-mu[1]),
                        sr=float(np.nanstd(a[:,0])/np.sqrt(len(a))), sl=float(a[:,1].std()/np.sqrt(len(a))),
                        v=a[:,2].tolist(), r=a[:,0].tolist(), l=a[:,1].tolist()))
                print(g,N,s,round(time.time()-t0,1),flush=True)
    json.dump(res,open("piloto_curv.json","w"))
