"""PILOTO DECLARADO (asesora de tamaño finito). Solo geometrías positivas. Semillas SeedSequence(900001, spawn_key=...).
No es evidencia; no usa entropías oficiales."""
import sys, json, time
sys.path.insert(0, "/home/user/FichaClaraV2"); sys.path.insert(0, "/home/user/FichaClaraV2/tools")
import numpy as np
from omega.causality.sprinkling import sprinkle_diamond, causal_matrix
from omega.causality.torus import sprinkle_slab, slab_causal_matrix
from omega.causality.estimators import between_counts
from f7 import de_sitter_eta_points

ENT = 900001
N = 4000
MS = (5, 10, 20, 40, 80, 160)
KMAX = 8
NSAMP = 300
T = 0.45
def rng(*k): return np.random.Generator(np.random.PCG64(np.random.SeedSequence(ENT, spawn_key=k)))

def frw_slab(n, r):  # a ∝ η, η = t + 0.5 ∈ [0.5, 0.95], densidad ∝ η^4, toro
    out = []; tot = 0
    while tot < n:
        c = sprinkle_slab(4 * n, r, T); eta = c[:, 0] + 0.5
        k = c[r.random(len(c)) < (eta / 0.95) ** 4]; out.append(k); tot += len(k)
    return np.concatenate(out)[:n]

def gen(name, s):
    r = rng(1, list(STR).index(name), s)
    if name == "diamante": x = sprinkle_diamond(N, 4, r); return x, causal_matrix(x)
    if name == "bloque": x = sprinkle_slab(N, r, T); return x, slab_causal_matrix(x)
    if name == "caja": x = sprinkle_slab(N, r, T); return x, causal_matrix(x)
    if name == "de_sitter": x = de_sitter_eta_points(N, 4, r); return x, causal_matrix(x)
    if name == "de_sitter_fuerte": x = de_sitter_eta_points(N, 4, r, -1.2, -0.2); return x, causal_matrix(x)
    if name == "frw_bloque": x = frw_slab(N, r); return x, slab_causal_matrix(x)
STR = {"diamante": 6, "bloque": 6, "caja": 4, "de_sitter": 4, "de_sitter_fuerte": 4, "frw_bloque": 4}

def tlim(name):
    return {"diamante": (-.5, .5), "de_sitter": (-.5, .5), "de_sitter_fuerte": (-.5, .5)}.get(name, (0, T))

def local_abund(B, idx):
    sub = B[np.ix_(idx, idx)]
    v = sub[sub >= 0]
    return np.bincount(v, minlength=KMAX + 1)[:KMAX + 1]

def truncated(x, i, j):  # ¿el diamante geométrico de (i,j) sale de la caja [0,1]^3?
    a, b = x[i, 1:], x[j, 1:]; A = (x[j, 0] - x[i, 0]) / 2; c = (a + b) / 2
    f = np.linalg.norm(b - a) / 2; u = (b - a) / (2 * f) if f > 0 else np.zeros(3)
    Bm = np.sqrt(max(A * A - f * f, 0)); ext = np.sqrt(A * A * u * u + Bm * Bm * (1 - u * u))
    return bool(np.any(c - ext < 0) or np.any(c + ext > 1))

def density_ratio(name, x, i, j):
    if name == "de_sitter": e0, e1 = x[i, 0] - 1, x[j, 0] - 1; return (e0 / e1) ** 4
    if name == "de_sitter_fuerte": e0, e1 = x[i, 0] - 0.7, x[j, 0] - 0.7; return (e0 / e1) ** 4
    if name == "frw_bloque": return ((x[j, 0] + .5) / (x[i, 0] + .5)) ** 4
    return 1.0

# referencia exacta: m uniformes en diamante 4D + dos puntas
REF = {}
def ref(m, reps=3000):
    if m in REF: return REF[m]
    r = rng(0, m); acc = np.zeros(KMAX + 1); fr = 0.0
    tips = np.array([[-.5, 0, 0, 0], [.5, 0, 0, 0]])
    for _ in range(reps):
        x = np.vstack([tips, sprinkle_diamond(m, 4, r)]); c = causal_matrix(x); B = between_counts(c)
        v = B[B >= 0]; acc += np.bincount(v, minlength=KMAX + 1)[:KMAX + 1]; fr += c.sum() / ((m + 2) * (m + 1) / 2)
    REF[m] = (acc / reps, fr / reps); return REF[m]

def dforma(nk, rf):
    v = np.log((nk + 1) / (rf + 1)); return float(np.max(np.abs(v - v.mean()))), float(np.max(np.abs(v)))

rows = []
t0 = time.time()
for name, ns in STR.items():
    for s in range(ns):
        x, C = gen(name, s); B = between_counts(C); lo, hi = tlim(name)
        for m in MS:
            ii, jj = np.nonzero(B == m); cnt = len(ii)
            row = {"estr": name, "seed": s, "m": m, "n_int": int(cnt)}
            if cnt == 0: rows.append(row); continue
            mid = ((x[ii, 0] + x[jj, 0]) / 2 - lo) / (hi - lo)
            row["tmid_cuartiles"] = np.percentile(mid, [10, 50, 90]).round(3).tolist()
            r = rng(2, list(STR).index(name), s, m); sel = r.choice(cnt, size=min(NSAMP, cnt), replace=False)
            profs = []; members = []; trunc = []; dens = []; frac = []
            for q in sel:
                i, j = ii[q], jj[q]; inner = np.flatnonzero(C[i] & C[:, j])
                idx = np.concatenate([[i], inner, [j]]); members.append(set(inner.tolist()))
                profs.append(local_abund(B, idx)); frac.append(C[np.ix_(idx, idx)].sum() / ((m + 2) * (m + 1) / 2))
                if name == "caja": trunc.append(truncated(x, i, j))
                dens.append(density_ratio(name, x, i, j))
            P = np.array(profs, float); rf, rfr = ref(m)
            mean = P.mean(0); row["n_sel"] = len(sel)
            row["dforma"], row["D"] = dforma(mean, rf)
            row["v"] = np.log((mean + 1) / (rf + 1)).round(4).tolist()
            row["links_mean"] = float(mean[0]); row["links_se_ingenuo"] = float(P[:, 0].std(ddof=1) / np.sqrt(len(sel))) if len(sel) > 1 else None
            row["frac_rel"] = float(np.mean(frac)); row["frac_rel_ref"] = float(rfr)
            # solapamiento: fracción de pares de intervalos muestreados que comparten algún elemento interior
            k = min(len(members), 100); sh = 0; tot = 0
            for a in range(k):
                for b in range(a + 1, k):
                    tot += 1; sh += bool(members[a] & members[b])
            row["frac_pares_solapados"] = sh / tot if tot else None
            allm = [e for st in members for e in st]; row["multiplicidad_media"] = len(allm) / max(1, len(set(allm)))
            if trunc:
                tr = np.array(trunc); row["frac_recortados"] = float(tr.mean())
                if tr.sum() > 5 and (~tr).sum() > 5:
                    row["dforma_recort"] = dforma(P[tr].mean(0), rf)[0]; row["dforma_no_recort"] = dforma(P[~tr].mean(0), rf)[0]
                    row["links_rel_recort"] = float(P[tr, 0].mean() / rf[0]); row["links_rel_no_recort"] = float(P[~tr, 0].mean() / rf[0])
            if name in ("de_sitter", "de_sitter_fuerte", "frw_bloque"):
                row["dens_ratio_mediana"] = float(np.median(dens)); row["dens_ratio_p90"] = float(np.percentile(dens, 90))
            rows.append(row)
        print(name, s, round(time.time() - t0, 1), flush=True)
json.dump(rows, open("/tmp/claude-0/-home-user-FichaClaraV2/90c7523c-b9d2-5eee-afba-260b3848fae0/scratchpad/asesores/piloto_tamano.json", "w"), indent=1)
