"""OMEGA E6-S v1.1 (docs/OMEGA_E6_1_PRERREGISTRO.md, congelado): panel cerrado por transporte + controles frescos.

Reutiliza tools/e6_static.py (v1, congelado) sin modificarlo. Salida: results/e6_1/{static_table.json, summary.json}.
Uso: python tools/e6_static_v11.py

Decisiones de implementacion:
 * Operadores de transporte sobre Z^d (torus_lattice(side, d)): truncacion T (nodo (v,j), j<d = +e_j, j>=d = -e_{j-d}; ciclo
   j~j+1 mod 2d; arista v->v+e_i une (v,i) con (v+e_i, d+i)), grafo de lineas L (arista (v,i) = v*d+i; adyacentes si comparten
   extremo), producto Z^d x K2 (nodo (v,s) = s*n+v). Geometrias deterministas (semilla 0).
 * Si N > 20000: reglas de vertice sobre 3000 vertices sin reposicion con rng_from_key((MASTER_E6, id, 0, 8)), ordenados;
   cantidades locales calculadas exactamente solo para esos vertices (sub-matriz de vecinos). Si N <= 20000: todos los vertices
   con e6_static.local_quantities (identico a v1).
 * C-5CIC: arista uv en un 5-ciclo u-v-x-y-z-u con x in N(v)\\{u}, z in N(u)\\{v}, x != z, y in N(x)&N(z) \\ {u,v}.
 * C-PAR: grado par. P0-S2(c): |S2(v)| = #vertices a distancia exacta 2 = c, c in {8,18,32,50}.
 * Las reglas de arista usan la muestra de aristas de v1 (misma clave (MASTER_E6, id, semilla, 7)).
"""

from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path
from typing import Any, Callable

import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import e6_static as e6  # noqa: E402
from e6_static import MASTER_E6, GEOM_ID as GEOM_ID_V1, rng_from_key  # noqa: E402
from omega.dynamics.sq import s_edge, t_edge  # noqa: E402
from p1d3_panel import _from_edges, torus_lattice  # noqa: E402

OUT = e6.ROOT / "results" / "e6_1"
N_VERT_SAMPLE = 3000
N_VERT_THRESHOLD = 20000

SIDES_T = {2: 50, 3: 12, 4: 6, 5: 5, 6: 5}
SIDES_L = {1: 10000, 2: 100, 3: 20, 4: 10, 5: 6, 6: 5}
SIDES_K2 = {1: 5000, 2: 70, 3: 16, 4: 8, 5: 6, 6: 5}

GEOM_ID: dict[str, int] = dict(GEOM_ID_V1)
for _d, _s in SIDES_T.items():
    GEOM_ID[f"T_Z{_d}_{_s}"] = 100 + _d
for _d, _s in SIDES_L.items():
    GEOM_ID[f"L_Z{_d}_{_s}"] = 110 + _d
for _d, _s in SIDES_K2.items():
    GEOM_ID[f"Z{_d}K2_{_s}"] = 120 + _d

GRID_S2 = (8, 18, 32, 50)

# ------------------------------------------------------------------ operadores de transporte


def _zd_edges(side: int, d: int) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    """(n, tails, heads, axis) de las aristas de Z^d torus; la arista (v,i) tiene indice v*d+i."""
    n = side**d
    idx = np.arange(n, dtype=np.int64)
    coords = np.stack(np.unravel_index(idx, (side,) * d), axis=1)
    heads = np.empty((n, d), dtype=np.int64)
    for ax in range(d):
        c = coords.copy()
        c[:, ax] = (c[:, ax] + 1) % side
        heads[:, ax] = np.ravel_multi_index(tuple(c.T), (side,) * d)
    return n, idx, heads, coords


def truncation(side: int, d: int) -> sparse.csr_array:
    assert d >= 2
    n, idx, heads, _ = _zd_edges(side, d)
    k = 2 * d
    us, vs = [], []
    for j in range(k):  # ciclo del gadget
        us.append(idx * k + j)
        vs.append(idx * k + (j + 1) % k)
    for i in range(d):  # aristas de Z^d
        us.append(idx * k + i)
        vs.append(heads[:, i] * k + d + i)
    return _from_edges(n * k, np.concatenate(us), np.concatenate(vs))


def line_graph(side: int, d: int) -> sparse.csr_array:
    n, idx, heads, _ = _zd_edges(side, d)
    inc = np.empty((n, 2 * d), dtype=np.int64)  # aristas incidentes a v
    for i in range(d):
        inc[:, i] = idx * d + i
        inc[heads[:, i], d + i] = idx * d + i  # v es la cabeza de (v,i): arista entrante de heads[v,i]
    us, vs = [], []
    for a, b in itertools.combinations(range(2 * d), 2):
        us.append(inc[:, a])
        vs.append(inc[:, b])
    return _from_edges(n * d, np.concatenate(us), np.concatenate(vs))


def prod_k2(side: int, d: int) -> sparse.csr_array:
    A = torus_lattice(side, d)
    n = A.shape[0]
    t = sparse.triu(A, k=1, format="coo")
    u, v = np.asarray(t.row, dtype=np.int64), np.asarray(t.col, dtype=np.int64)
    idx = np.arange(n, dtype=np.int64)
    return _from_edges(2 * n, np.concatenate([u, u + n, idx]), np.concatenate([v, v + n, idx + n]))


def panel_v11() -> dict[str, tuple[int, Callable[[Any], sparse.csr_array], bool]]:
    P = dict(e6._panel())
    for d, s in SIDES_T.items():
        P[f"T_Z{d}_{s}"] = (d, (lambda r, d=d, s=s: truncation(s, d)), False)
    for d, s in SIDES_L.items():
        P[f"L_Z{d}_{s}"] = (d, (lambda r, d=d, s=s: line_graph(s, d)), False)
    for d, s in SIDES_K2.items():
        P[f"Z{d}K2_{s}"] = (d, (lambda r, d=d, s=s: prod_k2(s, d)), False)
    assert set(P) == set(GEOM_ID)
    return P


# ------------------------------------------------------------------ cantidades locales en un subconjunto de vertices


def local_quantities_subset(A: sparse.csr_array, verts: np.ndarray) -> dict[str, Any]:
    """deg, tri, c4, f4 exactos para los vertices dados (misma semantica que e6.local_quantities)."""
    A = sparse.csr_array(A)
    n = verts.size
    deg = np.diff(A.indptr).astype(np.int64)[verts]
    tri = np.zeros(n, dtype=bool)
    c4 = np.zeros(n, dtype=np.int64)
    f4 = np.zeros(n, dtype=np.float64)
    for k, v in enumerate(verts):
        dv = int(deg[k])
        if dv < 2:
            continue
        nb = A.indices[A.indptr[v]:A.indptr[v + 1]]
        Anb = A[nb]
        tri[k] = Anb[:, nb].sum() > 0
        sub = (Anb @ Anb.T).toarray()
        np.fill_diagonal(sub, 0)
        npairs = dv * (dv - 1) // 2
        c4[k] = int(round(sub.sum() / 2.0)) - npairs
        f4[k] = float((sub >= 2).sum() / 2.0) / npairs
    return {"n": n, "deg": deg, "tri": tri, "c4": c4, "f4": f4}


def s2_sizes(A: sparse.csr_array, verts: np.ndarray) -> np.ndarray:
    out = np.zeros(verts.size, dtype=np.int64)
    ip, ix = A.indptr, A.indices
    for k, v in enumerate(verts):
        nb = ix[ip[v]:ip[v + 1]]
        if nb.size == 0:
            continue
        cand = np.unique(np.concatenate([ix[ip[a]:ip[a + 1]] for a in nb]))
        out[k] = np.setdiff1d(cand, np.concatenate([nb, [v]]), assume_unique=False).size
    return out


# ------------------------------------------------------------------ reglas de arista frescas


def in_5cycle(adj: list[set[int]], u: int, v: int) -> bool:
    """Existe ciclo simple u-v-x-y-z-u (5 vertices distintos)."""
    for x in adj[v]:
        if x == u:
            continue
        for z in adj[u]:
            if z == v or z == x:
                continue
            if (adj[x] & adj[z]) - {u, v}:
                return True
    return False


def edge_rest_v11(A: sparse.csr_array, u: np.ndarray, v: np.ndarray) -> dict[str, np.ndarray]:
    adj = [set(A.indices[A.indptr[i]:A.indptr[i + 1]].tolist()) for i in range(A.shape[0])]
    f = lambda fn: np.array([fn(adj, int(a), int(b)) for a, b in zip(u, v)], dtype=bool)  # noqa: E731
    return {"C-SQ": f(s_edge), "C-TRI": f(t_edge), "C-5CIC": f(in_5cycle)}


# ------------------------------------------------------------------ reglas / clases

NEW_RULES = ["C-5CIC", "C-PAR"] + [f"P0-S2({c})" for c in GRID_S2]


def rule_names() -> list[str]:
    return e6.rule_names() + NEW_RULES


def rule_activity(loc: dict[str, Any], edges: dict[str, np.ndarray], s2: np.ndarray) -> dict[str, float]:
    a = e6.rule_activity(loc, edges)
    a["C-5CIC"] = float(edges["C-5CIC"].mean())
    a["C-PAR"] = float((loc["deg"] % 2 == 0).mean())
    for c in GRID_S2:
        a[f"P0-S2({c})"] = float((s2 == c).mean())
    return a


def families() -> dict[str, list[str]]:
    f = e6.families()
    f["C-5CIC"] = ["C-5CIC"]
    f["C-PAR"] = ["C-PAR"]
    f["P0-S2"] = [f"P0-S2({c})" for c in GRID_S2]
    return f


def verdict(cls: dict[str, dict[str, Any]], dial: dict[str, bool]) -> dict[str, Any]:
    v1 = dial["P0-DEG"] and dial["P0-SQV"] and dial["P0-S2"]
    v2 = (cls["P0-SQV(12)"]["class"] == "SELECTOR" and 3 in cls["P0-SQV(12)"]["R"]
          and cls["P0-SQV(24)"]["class"] == "SELECTOR" and 4 in cls["P0-SQV(24)"]["R"]
          and cls["P0-S2(18)"]["class"] == "SELECTOR" and 3 in cls["P0-S2(18)"]["R"]
          and cls["P0-S2(32)"]["class"] == "SELECTOR" and 4 in cls["P0-S2(32)"]["R"])
    v3 = (all(cls[r]["class"] != "SELECTOR" for r in ("C-TRIV", "C-SQ", "C-TRI", "P0-LAT(2)", "C-5CIC", "C-PAR"))
          and not dial["P0-LAT"])
    return {"V1p": v1, "V2p": v2, "V3p": v3, "verdict": "VALIDO" if (v1 and v2 and v3) else "INVALIDO"}


# ------------------------------------------------------------------ ejecucion


def run_geometry(name: str, dim: int, ctor: Callable[[Any], sparse.csr_array], rand: bool) -> dict[str, Any]:
    gid = GEOM_ID[name]
    per_seed = []
    info: dict[str, Any] = {}
    for seed in (e6.SEEDS_RANDOM if rand else (0,)):
        rng = rng_from_key((MASTER_E6, gid, seed)) if rand else None
        adj = ctor(rng)
        A = sparse.csr_array(adj)
        A.setdiag(0)
        A.eliminate_zeros()
        A.data[:] = 1.0
        N = A.shape[0]
        deg_all = np.diff(A.indptr)
        sampled = N > N_VERT_THRESHOLD
        if sampled:
            verts = np.sort(rng_from_key((MASTER_E6, gid, seed, 8)).choice(N, size=N_VERT_SAMPLE, replace=False))
            loc = local_quantities_subset(A, verts)
        else:
            verts = np.arange(N)
            loc = e6.local_quantities(A)
        s2 = s2_sizes(A, verts)
        u, v = e6.sample_edges(A, (MASTER_E6, gid, seed, 7))
        act = rule_activity(loc, edge_rest_v11(A, u, v), s2)
        per_seed.append(act)
        info = {"N": int(N), "E": int(deg_all.sum() // 2), "mean_deg": float(deg_all.mean()),
                "max_deg": int(deg_all.max()), "vertex_sampled": bool(sampled), "n_vertices_eval": int(verts.size)}
    a = {r: float(np.mean([s[r] for s in per_seed])) for r in per_seed[0]}
    return {"name": name, "dim": dim, "seeds": len(per_seed), "a": a, "a_per_seed": per_seed, **info}


def main() -> int:
    t0 = time.time()
    geoms = []
    for name, (dim, ctor, rand) in panel_v11().items():
        t1 = time.time()
        g = run_geometry(name, dim, ctor, rand)
        geoms.append(g)
        print(f"  {name:16s} d={dim} N={g['N']:6d} <k>={g['mean_deg']:.2f} kmax={g['max_deg']:3d}"
              f"{' [muestreo]' if g['vertex_sampled'] else ''} ({time.time() - t1:.1f}s)", flush=True)
    dims = sorted({g["dim"] for g in geoms})
    table: dict[str, Any] = {}
    cls: dict[str, dict[str, Any]] = {}
    for r in rule_names():
        m = {d: max(g["a"][r] for g in geoms if g["dim"] == d) for d in dims}
        abar = {d: float(np.mean([g["a"][r] for g in geoms if g["dim"] == d])) for d in dims}
        c = e6.m_profile_classes(m)
        cls[r] = c
        table[r] = {"a_geom": {g["name"]: g["a"][r] for g in geoms}, "m": m, "abar": abar, **c}
    fam_out: dict[str, Any] = {}
    dial: dict[str, bool] = {}
    for fam, members in families().items():
        dial[fam] = e6.dial_flag({r: cls[r] for r in members}) if len(members) > 1 else False
        fam_out[fam] = {"members": {r: cls[r]["class"] for r in members}, "DIAL": dial[fam]}
    ver = verdict(cls, dial)
    runtime = time.time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    from omega.experiments.v11.gate import head_commit
    (OUT / "static_table.json").write_text(json.dumps(
        {"master": MASTER_E6, "geometries": [{k: v for k, v in g.items() if k != "a_per_seed"} for g in geoms],
         "rules": table}, indent=1))
    summary = {"master": MASTER_E6, "head_commit": head_commit(), "families": fam_out,
               "rule_classes": {r: {"class": c["class"], "R": c["R"], "Q": c["Q"]} for r, c in cls.items()},
               **ver, "runtime_s": runtime}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1))
    print("\nm(d) por regla (R = m>=0.9, Q = m<=0.5):")
    print(f"{'regla':14s} " + " ".join(f"d={d:<5d}" for d in dims) + "  clase")
    for r in rule_names():
        print(f"{r:14s} " + " ".join(f"{table[r]['m'][d]:7.3f}" for d in dims) + f"  {cls[r]['class']} R={cls[r]['R']}")
    print("\nFamilias:", {f: o["DIAL"] for f, o in fam_out.items()})
    print(f"V1'={ver['V1p']} V2'={ver['V2p']} V3'={ver['V3p']} -> E6-S v1.1 {ver['verdict']}  ({runtime:.0f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
