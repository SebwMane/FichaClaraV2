"""OMEGA E6-S: calibracion del brazo estatico (docs/OMEGA_E6_PRERREGISTRO.md §3-§4, congelado).

Salida: results/e6/{static_table.json, summary.json}. Uso: python tools/e6_static.py

Decisiones de implementacion:
 * Cantidades locales solo desde la estructura del grafo (matriz de adyacencia): deg(v); tri(v) = existe arista entre dos
   vecinos de v; c4(v) = sum_{a<b en N(v)} (|N(a) & N(b)| - 1); f4(v) = fraccion de pares a<b en N(v) con |N(a) & N(b)| >= 2.
   f4 := 0 si deg < 2 (no hay pares). Los pares con |N(a)&N(b)| = 1 solo comparten v.
 * Arista en 4-ciclo: omega.dynamics.sq.s_edge (misma semantica que el reposo de SQ, B-s). Arista en triangulo: t_edge.
 * P0-LAT(q): reposo sii no hay triangulo en v y 2*c4(v) == deg*(deg-q) (aritmetica entera exacta); si deg*(deg-q) es impar
   la igualdad es imposible, luego no hay reposo (deg(deg-q)/2 no entero).
 * P0-SQF: |f4 - f| <= 0.01 inclusivo, con holgura 1e-12 contra error de coma flotante. CIEGO: max m - min m <= 0.10 (+1e-12).
 * Aristas: se muestrean min(3000, |E|) aristas u<v (orden lexicografico por codigo) uniformes sin reposicion con
   rng_from_key((MASTER_E6, id, semilla, 7)); todas las reglas de arista usan la misma muestra de la instancia.
 * Geometrias deterministas: semilla 0 (una instancia). Aleatorias (RGG k=8,16, RGG1 anillo, toros+20 % diagonales): semillas 0 y 1,
   generador con rng_from_key((MASTER_E6, id, semilla)). Los ids son una tabla fija (GEOM_ID).
 * Toro con diagonales: rc3.torus_with_diagonals(side, dim, 0.2, rng). Panal: honeycomb_torus(100, 100). RGG periodico N=8000;
   RGG1 anillo k=10 N=10000 via rc3.random_ring. Circulante C_n(1,2,3) via _from_edges.
 * Diamante: caja periodica entera de lado 4L (unidades a/4); subred A = (x,y,z) pares con x+y+z = 0 mod 4; subred B = A+(1,1,1);
   aristas A -> A+(1,1,1), A+(1,-1,-1), A+(-1,1,-1), A+(-1,-1,1). N = 8 L^3 (L=10: 8000), grado 4, cintura 6.
 * DIAL: conjuntos R no vacios con >= 2 valores distintos (como conjuntos) y al menos un c con clase SELECTOR.
 * MONOTONO: R "intervalo" = dimensiones consecutivas; debe contener 1 o 6; toda d fuera de R pertenece a Q.
 * Para reglas sin rejilla (C-*), la "familia" es la propia regla (sin DIAL).
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import itertools  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any, Callable  # noqa: E402

import numpy as np  # noqa: E402
from scipy import sparse  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from a0_rc2 import cycle, honeycomb_torus, triangular_torus  # noqa: E402
from omega.c0.references import rng_from_key  # noqa: E402
from omega.dynamics.sq import s_edge, t_edge  # noqa: E402
from p1d3_panel import _from_edges, rgg, torus_lattice  # noqa: E402
from rc3 import bcc_lattice, fcc_lattice, random_ring, torus_with_diagonals  # noqa: E402

MASTER_E6 = 20261020
OUT = ROOT / "results" / "e6"
N_EDGES = 3000
EPS = 1e-12

# ------------------------------------------------------------------ tabla fija de ids de geometria

GEOM_ID: dict[str, int] = {
    "C10000": 1, "circ10000_123": 2, "RGG1_k10": 3,
    "Z2_100": 10, "tri_90": 11, "panal_100x100": 12, "RGG2_k8": 13, "RGG2_k16": 14, "Z2diag_90": 15,
    "Z3_20": 20, "FCC_L12": 21, "BCC_L10": 22, "diamante_L10": 23, "RGG3_k8": 24, "RGG3_k16": 25, "Z3diag_20": 26,
    "Z4_10": 30, "RGG4_k8": 31, "RGG4_k16": 32, "Z4diag_10": 33,
    "Z5_6": 40, "RGG5_k8": 41, "RGG5_k16": 42, "Z5diag_6": 43,
    "Z6_5": 50, "RGG6_k8": 51, "RGG6_k16": 52, "Z6diag_5": 53,
}
SEEDS_RANDOM = (0, 1)

# ------------------------------------------------------------------ generadores nuevos


def circulant(n: int, offsets: tuple[int, ...]) -> sparse.csr_array:
    i = np.arange(n, dtype=np.int64)
    return _from_edges(n, np.concatenate([i for _ in offsets]), np.concatenate([(i + o) % n for o in offsets]))


def diamond_lattice(L: int) -> sparse.csr_array:
    """Diamante periodico: 8 L^3 nodos, grado 4, cintura 6 (ver docstring del modulo)."""
    S = 4 * L
    g = np.stack(np.meshgrid(*(np.arange(S),) * 3, indexing="ij"), axis=-1).reshape(-1, 3).astype(np.int64)
    isA = ((g % 2 == 0).all(axis=1)) & ((g.sum(axis=1) % 4) == 0)
    cA = g[isA]
    # sitios B = A + (1,1,1)
    shift = np.array([1, 1, 1], dtype=np.int64)
    sites_B = (cA + shift) % S

    def key(c: np.ndarray) -> np.ndarray:
        return (c[:, 0] * S + c[:, 1]) * S + c[:, 2]

    nA = cA.shape[0]
    lab = -np.ones(S**3, dtype=np.int64)
    lab[key(cA)] = np.arange(nA)
    lab[key(sites_B)] = nA + np.arange(nA)
    us, vs = [], []
    for o in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)):
        t = (cA + np.array(o, dtype=np.int64)) % S
        dst = lab[key(t)]
        assert (dst >= nA).all()
        us.append(np.arange(nA, dtype=np.int64))
        vs.append(dst)
    return _from_edges(2 * nA, np.concatenate(us), np.concatenate(vs))


# ------------------------------------------------------------------ panel

# nombre -> (dim, constructor(rng|None), aleatoria)
def _panel() -> dict[str, tuple[int, Callable[[Any], sparse.csr_array], bool]]:
    P: dict[str, tuple[int, Callable[[Any], sparse.csr_array], bool]] = {
        "C10000": (1, lambda r: cycle(10000), False),
        "circ10000_123": (1, lambda r: circulant(10000, (1, 2, 3)), False),
        "RGG1_k10": (1, lambda r: random_ring(10000, 10, r), True),
        "Z2_100": (2, lambda r: torus_lattice(100, 2), False),
        "tri_90": (2, lambda r: triangular_torus(90), False),
        "panal_100x100": (2, lambda r: honeycomb_torus(100, 100), False),
        "Z3_20": (3, lambda r: torus_lattice(20, 3), False),
        "FCC_L12": (3, lambda r: fcc_lattice(12), False),
        "BCC_L10": (3, lambda r: bcc_lattice(10), False),
        "diamante_L10": (3, lambda r: diamond_lattice(10), False),
        "Z4_10": (4, lambda r: torus_lattice(10, 4), False),
        "Z5_6": (5, lambda r: torus_lattice(6, 5), False),
        "Z6_5": (6, lambda r: torus_lattice(5, 6), False),
    }
    for d in range(2, 7):
        for k in (8, 16):
            P[f"RGG{d}_k{k}"] = (d, (lambda r, d=d, k=k: rgg(8000, d, k, r, True)), True)
    for d, side in ((2, 90), (3, 20), (4, 10), (5, 6), (6, 5)):
        P[f"Z{d}diag_{side}"] = (d, (lambda r, d=d, side=side: torus_with_diagonals(side, d, 0.2, r)), True)
    assert set(P) == set(GEOM_ID)
    return P


# ------------------------------------------------------------------ cantidades locales


def local_quantities(adj: sparse.csr_array) -> dict[str, Any]:
    """deg, tri(v), c4(v), f4(v) por vertice (solo estructura del grafo)."""
    A = sparse.csr_array(adj)
    A.setdiag(0)
    A.eliminate_zeros()
    A.data[:] = 1.0
    n = A.shape[0]
    deg = np.diff(A.indptr).astype(np.int64)
    A2 = sparse.csr_array(A @ A)
    A2.setdiag(0)
    A2.eliminate_zeros()
    tri = np.asarray(A2.multiply(A).sum(axis=1)).ravel() > 0
    c4 = np.zeros(n, dtype=np.int64)
    f4 = np.zeros(n, dtype=np.float64)
    for v in range(n):
        dv = int(deg[v])
        if dv < 2:
            continue
        nb = A.indices[A.indptr[v]:A.indptr[v + 1]]
        sub = A2[nb][:, nb]
        npairs = dv * (dv - 1) // 2
        # sub simetrica, diagonal 0: suma/2 = sum_{a<b} |N(a)&N(b)| ; v siempre cuenta 1 en cada par
        c4[v] = int(round(sub.sum() / 2.0)) - npairs
        f4[v] = float((sub.data >= 2).sum() / 2.0) / npairs
    return {"n": n, "deg": deg, "tri": tri, "c4": c4, "f4": f4, "A": A}


def sample_edges(A: sparse.csr_array, key: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
    t = sparse.triu(A, k=1, format="coo")
    u, v = np.asarray(t.row, dtype=np.int64), np.asarray(t.col, dtype=np.int64)
    order = np.argsort(u * A.shape[0] + v)
    u, v = u[order], v[order]
    m = min(N_EDGES, u.size)
    pick = np.sort(rng_from_key(key).choice(u.size, size=m, replace=False))
    return u[pick], v[pick]


def edge_rest(A: sparse.csr_array, u: np.ndarray, v: np.ndarray) -> dict[str, np.ndarray]:
    adj = [set(A.indices[A.indptr[i]:A.indptr[i + 1]].tolist()) for i in range(A.shape[0])]
    sq = np.array([s_edge(adj, int(a), int(b)) for a, b in zip(u, v)], dtype=bool)
    tr = np.array([t_edge(adj, int(a), int(b)) for a, b in zip(u, v)], dtype=bool)
    return {"C-SQ": sq, "C-TRI": tr}


# ------------------------------------------------------------------ reglas

GRID_LAT = (1, 2, 3)
GRID_DEG = (2, 3, 4, 6, 8, 10, 12)
GRID_SQV = (4, 12, 24, 40)
GRID_SQF = ((2, 3), (4, 5), (6, 7), (8, 9))  # f = num/den


def rule_names() -> list[str]:
    return (["C-TRIV", "C-SQ", "C-TRI"] + [f"P0-LAT({q})" for q in GRID_LAT] + [f"P0-DEG({k})" for k in GRID_DEG]
            + [f"P0-SQV({c})" for c in GRID_SQV] + [f"P0-SQF({a}/{b})" for a, b in GRID_SQF])


def rule_activity(loc: dict[str, Any], edges: dict[str, np.ndarray]) -> dict[str, float]:
    deg, tri, c4, f4 = loc["deg"], loc["tri"], loc["c4"], loc["f4"]
    a: dict[str, float] = {"C-TRIV": 1.0, "C-SQ": float(edges["C-SQ"].mean()), "C-TRI": float(edges["C-TRI"].mean())}
    for q in GRID_LAT:
        a[f"P0-LAT({q})"] = float((~tri & (2 * c4 == deg * (deg - q))).mean())
    for k in GRID_DEG:
        a[f"P0-DEG({k})"] = float((deg == k).mean())
    for c in GRID_SQV:
        a[f"P0-SQV({c})"] = float((c4 == c).mean())
    for p, q in GRID_SQF:
        a[f"P0-SQF({p}/{q})"] = float((np.abs(f4 - p / q) <= 0.01 + EPS).mean())
    return a


# ------------------------------------------------------------------ clases (§3.3)


def m_profile_classes(m: dict[int, float]) -> dict[str, Any]:
    dims = sorted(m)
    vals = [m[d] for d in dims]
    R = [d for d in dims if m[d] >= 0.9]
    Q = [d for d in dims if m[d] <= 0.5]
    if max(vals) - min(vals) <= 0.10 + EPS:
        cls = "CIEGO"
    elif any(any(dp < d for dp in Q) and any(dpp > d for dpp in Q) for d in R):
        cls = "SELECTOR"
    elif (R and R == list(range(R[0], R[-1] + 1)) and (1 in R or 6 in R)
          and all(d in Q for d in dims if d not in R)):
        cls = "MONOTONO"
    else:
        cls = "INDETERMINADO"
    return {"class": cls, "R": R, "Q": Q}


def dial_flag(per_c: dict[str, dict[str, Any]]) -> bool:
    sets = {tuple(v["R"]) for v in per_c.values() if v["R"]}
    return len(sets) >= 2 and any(v["class"] == "SELECTOR" for v in per_c.values())


def families() -> dict[str, list[str]]:
    return {
        "C-TRIV": ["C-TRIV"], "C-SQ": ["C-SQ"], "C-TRI": ["C-TRI"],
        "P0-LAT": [f"P0-LAT({q})" for q in GRID_LAT], "P0-DEG": [f"P0-DEG({k})" for k in GRID_DEG],
        "P0-SQV": [f"P0-SQV({c})" for c in GRID_SQV], "P0-SQF": [f"P0-SQF({a}/{b})" for a, b in GRID_SQF],
    }


def verdict(cls: dict[str, dict[str, Any]], dial: dict[str, bool]) -> dict[str, Any]:
    v1 = dial["P0-DEG"] and dial["P0-SQV"]
    v2 = (cls["P0-SQV(12)"]["class"] == "SELECTOR" and 3 in cls["P0-SQV(12)"]["R"]
          and cls["P0-SQV(24)"]["class"] == "SELECTOR" and 4 in cls["P0-SQV(24)"]["R"])
    v3 = (all(cls[r]["class"] != "SELECTOR" for r in ("C-TRIV", "C-SQ", "C-TRI", "P0-LAT(2)")) and not dial["P0-LAT"])
    return {"V1": v1, "V2": v2, "V3": v3, "verdict": "VALIDO" if (v1 and v2 and v3) else "INVALIDO"}


# ------------------------------------------------------------------ ejecucion


def run_geometry(name: str, dim: int, ctor: Callable[[Any], sparse.csr_array], rand: bool) -> dict[str, Any]:
    gid = GEOM_ID[name]
    per_seed = []
    info: dict[str, Any] = {}
    for seed in (SEEDS_RANDOM if rand else (0,)):
        rng = rng_from_key((MASTER_E6, gid, seed)) if rand else None
        adj = ctor(rng)
        loc = local_quantities(adj)
        u, v = sample_edges(loc["A"], (MASTER_E6, gid, seed, 7))
        act = rule_activity(loc, edge_rest(loc["A"], u, v))
        per_seed.append(act)
        info = {"N": int(loc["n"]), "E": int(loc["deg"].sum() // 2), "mean_deg": float(loc["deg"].mean()),
                "max_deg": int(loc["deg"].max())}
    a = {r: float(np.mean([s[r] for s in per_seed])) for r in per_seed[0]}
    return {"name": name, "dim": dim, "seeds": len(per_seed), "a": a, "a_per_seed": per_seed, **info}


def main() -> int:
    t0 = time.time()
    geoms = []
    for name, (dim, ctor, rand) in _panel().items():
        t1 = time.time()
        g = run_geometry(name, dim, ctor, rand)
        geoms.append(g)
        print(f"  {name:16s} d={dim} N={g['N']:6d} <k>={g['mean_deg']:.2f} kmax={g['max_deg']:3d} ({time.time() - t1:.1f}s)", flush=True)
    dims = sorted({g["dim"] for g in geoms})
    table: dict[str, Any] = {}
    cls: dict[str, dict[str, Any]] = {}
    for r in rule_names():
        m = {d: max(g["a"][r] for g in geoms if g["dim"] == d) for d in dims}
        abar = {d: float(np.mean([g["a"][r] for g in geoms if g["dim"] == d])) for d in dims}
        c = m_profile_classes(m)
        cls[r] = c
        table[r] = {"a_geom": {g["name"]: g["a"][r] for g in geoms}, "m": m, "abar": abar, **c}
    fam_out: dict[str, Any] = {}
    dial: dict[str, bool] = {}
    for fam, members in families().items():
        dial[fam] = dial_flag({r: cls[r] for r in members}) if len(members) > 1 else False
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
    print(f"V1={ver['V1']} V2={ver['V2']} V3={ver['V3']} -> E6-S {ver['verdict']}  ({runtime:.0f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
