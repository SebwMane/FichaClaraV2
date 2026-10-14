import json, numpy as np
R=json.load(open("piloto2.json"))
def data(g,j,ms=(20,40,80)):
    out=[]
    for m in ms:
        for N in (1000,2000,4000):
            v=np.array(R[f"{g}|{N}|{m}"])[:,j]; out.append((m,N,v.mean(),v.std(ddof=1)/np.sqrt(len(v))))
    return np.array(out)
print("dl: amplitud propia por m, exponente comun en N: y = A_m N^-a (sin intercepto), m=20,40,80")
for g in ("dS_fuerte","AdS_fuerte"):
    D=data(g,1); best=None
    for a in np.linspace(0.0,1.5,151):
        chi=0
        for m in (20,40,80):
            s=D[D[:,0]==m]; x=s[:,1]**-a; w=1/s[:,3]**2; A=np.sum(w*x*s[:,2])/np.sum(w*x*x); chi+=np.sum(w*(s[:,2]-A*x)**2)
        if best is None or chi<best[0]: best=(chi,a)
    grid=np.linspace(0,1.5,151); 
    print(f"  {g}: a_mejor={best[1]:.2f} chi2={best[0]:.1f}/6")
print("dr: y = A (m/N)^a, sin intercepto, m=10..80; perfil chi2 para IC aproximado")
for g in ("dS_fuerte","AdS_fuerte"):
    D=data(g,0,(10,20,40,80)); x0=D[:,0]/D[:,1]; w=1/D[:,3]**2
    chis=[]
    for a in np.linspace(0.05,1.5,146):
        x=x0**a; A=np.sum(w*x*D[:,2])/np.sum(w*x*x); chis.append((np.sum(w*(D[:,2]-A*x)**2),a,A))
    c0=min(chis); ok=[c[1] for c in chis if c[0]<=c0[0]+1]
    print(f"  {g}: a={c0[1]:.2f} [{min(ok):.2f},{max(ok):.2f}] A={c0[2]:+.3f} chi2={c0[0]:.1f}/10")
