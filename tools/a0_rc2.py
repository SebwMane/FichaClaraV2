"""OMEGA-A0 Parte II: RC-2, juez de A, validacion fuera de muestra (docs/OMEGA_A0_PRERREGISTRO.md §II).

Salida: results/a0_rc2/{graphs.jsonl, summary.json}. Uso:
  python tools/a0_rc2.py [--procs 3] [--smoke]
Instrumentos y umbrales congelados de tools/cic_l0b.py (mismas subclaves de RNG: anillos key+(1,), localidad key+(2,),
coherencia key+(3,), curvatura key+(4,)); no hay Nivel II. Heisenberg (grupo E) se reporta aparte y no entra en veredictos.

Decisiones de lectura (ver informe): ver docstrings de cada generador.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from collections import defaultdict  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402
from scipy import sparse  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from c0_dynamics import _atomic_write  # noqa: E402
from cic_l0b import FNEG_MAX, KAPPA_HI, KAPPA_LO, LOCALITY_MIN, curvature_summary, locality_fraction  # noqa: E402
from coh_l0b import heisenberg  # noqa: E402
from p1d3_panel import _codes, _from_edges, _radius, patchwork, rgg, ring, torus_lattice  # noqa: E402

from omega.c0.references import rng_from_key  # noqa: E402
from omega.diagnostics.annulus import annulus_profile, annulus_status  # noqa: E402
from omega.diagnostics.coherence import coherence_status, edge_coherence_profile  # noqa: E402
from omega.experiments.v11.gate import head_commit  # noqa: E402

MASTER = 20261015
# Guarda de viabilidad (NO es un umbral de decision): el LP de W1 de ollivier_edge_sparse es denso en |sx|*|sy|;
# una arista hub-hub del 3-arbol s0 (1090 x 993 -> 1.0e6 variables) agota los 16 GB (OOM) y cuelga el Pool.
# Si el producto de soportes de alguna arista de la muestra supera MAX_LP_PRODUCT, la curvatura se declara
# no evaluable (mediana/f_neg = NaN => la condicion (iv) falla) en lugar de ejecutarse.
MAX_LP_PRODUCT = 400_000
NS = 20_000
SEEDS = (0, 1, 2)
OUT_FULL = ROOT / "results" / "a0_rc2"
OUT_SMOKE = ROOT / "results" / "a0_rc2_smoke"

FAMILY_ID: dict[str, int] = {
    "toro_triangular_141": 1, "panal_hexagonal": 2, "T4_12": 3, "RGG4_k16": 4, "RGG2_k20": 5,
    "RGG3_k12_del15": 6, "T3_27_diag20": 7, "cubo_abierto_27": 8, "RGG2_cilindro_k12": 9,
    "3arbol": 20, "apoloniana": 21, "arbol_de_rejillas": 22, "BA_m6": 23, "cuadrado_141_atajos1": 24,
    "T3_27_atajos0.1": 25, "retazos_b3": 26,
    "escalera_C5000xP4": 40, "cilindro_C2000xC10": 41, "anillo_k6": 42,
    "Heisenberg_Z27": 60,
}
SMOKE_FAMILIES = tuple(FAMILY_ID)  # humo: todas las familias, N~2000, semilla 0

# ------------------------------------------------------------------ generadores


def triangular_torus(side: int) -> sparse.csr_array:
    """Toro triangular: (i,j) ~ (i+1,j), (i,j+1), (i+1,j-1) mod side (y sus inversos): grado 6."""
    n = side * side
    i, j = np.meshgrid(np.arange(side), np.arange(side), indexing="ij")
    i, j = i.ravel().astype(np.int64), j.ravel().astype(np.int64)
    idx = i * side + j
    us, vs = [], []
    for di, dj in ((1, 0), (0, 1), (1, -1)):
        us.append(idx)
        vs.append(((i + di) % side) * side + (j + dj) % side)
    return _from_edges(n, np.concatenate(us), np.concatenate(vs))


def honeycomb_torus(rows: int, cols: int) -> sparse.csr_array:
    """Panal toroidal (ladrillo): (i,j) ~ (i+1,j) mod rows; vertical (i,j)~(i,j+1) si (i+j) par. `rows`,`cols` pares: grado 3, cintura 6."""
    if rows % 2 or cols % 2:
        raise ValueError("rows y cols deben ser pares")
    n = rows * cols
    i, j = np.meshgrid(np.arange(rows), np.arange(cols), indexing="ij")
    i, j = i.ravel().astype(np.int64), j.ravel().astype(np.int64)
    idx = i * cols + j
    even = (i + j) % 2 == 0
    u = np.concatenate([idx, idx[even]])
    v = np.concatenate([((i + 1) % rows) * cols + j, i[even] * cols + (j[even] + 1) % cols])
    return _from_edges(n, u, v)


def lattice_open(side: int, dim: int) -> sparse.csr_array:
    """Red cubica 27^3 sin periodicidad (vecinos a distancia 1)."""
    n = side**dim
    idx = np.arange(n, dtype=np.int64)
    coords = np.stack(np.unravel_index(idx, (side,) * dim), axis=1)
    us, vs = [], []
    for ax in range(dim):
        ok = coords[:, ax] < side - 1
        c = coords[ok].copy()
        c[:, ax] += 1
        us.append(idx[ok])
        vs.append(np.ravel_multi_index(tuple(c.T), (side,) * dim).astype(np.int64))
    return _from_edges(n, np.concatenate(us), np.concatenate(vs))


def giant(adj: sparse.csr_array) -> sparse.csr_array:
    _, lab = connected_components(adj, directed=False)
    big = np.argmax(np.bincount(lab))
    keep = np.flatnonzero(lab == big)
    return sparse.csr_array(adj[keep][:, keep])


def rgg_deleted(n: int, d: int, k: float, frac: float, rng: np.random.Generator) -> sparse.csr_array:
    """RGG toroidal; cada arista se elimina con prob. exacta-en-numero `frac` (muestra uniforme sin reemplazo); componente gigante."""
    a = rgg(n, d, k, rng)
    t = sparse.triu(a, k=1, format="coo")
    u, v = np.asarray(t.row, dtype=np.int64), np.asarray(t.col, dtype=np.int64)
    keep = np.sort(rng.permutation(u.size)[int(round(frac * u.size)):])
    return giant(_from_edges(n, u[keep], v[keep]))


def _edge_codes(adj: sparse.csr_array) -> np.ndarray:
    t = sparse.triu(adj, k=1, format="coo")
    return _codes(adj.shape[0], np.asarray(t.row, dtype=np.int64), np.asarray(t.col, dtype=np.int64))


def torus_with_diagonals(side: int, frac: float, rng: np.random.Generator) -> sparse.csr_array:
    """T3 periodico + round(frac*|E|) diagonales de cara (offsets con dos coordenadas +-1 de 3: 12 offsets) distintas y nuevas."""
    base = torus_lattice(side, 3)
    n = side**3
    have = _edge_codes(base)
    target = int(round(frac * have.size))
    offs = np.array([[a, b, 0] for a in (-1, 1) for b in (-1, 1)] + [[a, 0, b] for a in (-1, 1) for b in (-1, 1)]
                    + [[0, a, b] for a in (-1, 1) for b in (-1, 1)], dtype=np.int64)
    new = np.empty(0, dtype=np.int64)
    while new.size < target:
        m = 2 * (target - new.size) + 16
        x = rng.integers(n, size=m)
        o = offs[rng.integers(offs.shape[0], size=m)]
        c = np.stack(np.unravel_index(x, (side,) * 3), axis=1)
        y = np.ravel_multi_index(tuple(((c + o) % side).T), (side,) * 3).astype(np.int64)
        cc = np.unique(_codes(n, x, y))
        cc = cc[~np.isin(cc, have)]
        new = np.unique(np.concatenate([new, cc]))
    new = rng.permutation(new)[:target]
    allc = np.concatenate([have, new])
    return _from_edges(n, allc // n, allc % n)


def add_random_shortcuts(adj: sparse.csr_array, frac: float, rng: np.random.Generator) -> sparse.csr_array:
    """Anade round(frac*|E|) aristas entre pares de vertices uniformes (distintos; duplicados colapsan en la matriz binaria)."""
    n = adj.shape[0]
    t = sparse.triu(adj, k=1, format="coo")
    m = int(round(frac * t.nnz))
    a = rng.integers(n, size=m)
    b = (a + 1 + rng.integers(n - 1, size=m)) % n
    return _from_edges(n, np.concatenate([np.asarray(t.row, dtype=np.int64), a]), np.concatenate([np.asarray(t.col, dtype=np.int64), b]))


def rgg_cylinder(n: int, d_k: float, rng: np.random.Generator) -> sparse.csr_array:
    """RGG2 en cuadrado unidad periodico solo en x. cKDTree acepta boxsize por eje: boxsize=[1, 1000] hace que el
    envoltorio en y (distancia >= 999) no ocurra nunca, es decir y no periodico. Radio como rgg (sin correccion de borde)."""
    x = rng.random((n, 2))
    tree = cKDTree(x, boxsize=np.array([1.0, 1000.0]))
    p = tree.query_pairs(_radius(n, 2, d_k), output_type="ndarray")
    return _from_edges(n, p[:, 0].astype(np.int64), p[:, 1].astype(np.int64))


def random_three_tree(n: int, rng: np.random.Generator) -> sparse.csr_array:
    """3-arbol aleatorio: K4 inicial; cada nodo nuevo elige uniformemente un triangulo (3-clique) de la lista de
    triangulos creados (el original sigue disponible) y se une a sus 3 vertices; se anaden (v,a,b),(v,a,c),(v,b,c)."""
    tris = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    eu = [0, 0, 0, 1, 1, 2]
    ev = [1, 2, 3, 2, 3, 3]
    r = rng.random(n)
    for x in range(4, n):
        a, b, c = tris[int(r[x] * len(tris))]
        eu += [a, b, c]
        ev += [x, x, x]
        tris += [(x, a, b), (x, a, c), (x, b, c)]
    return _from_edges(n, np.asarray(eu, dtype=np.int64), np.asarray(ev, dtype=np.int64))


def apollonian(n: int, rng: np.random.Generator) -> sparse.csr_array:
    """Apoloniana aleatoria: K4 (4 caras); se elige una cara uniforme (a,b,c), se anade v unido a a,b,c y la cara se ELIMINA,
    sustituida por (v,a,b),(v,a,c),(v,b,c). Aristas = 3N-6."""
    faces = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    eu = [0, 0, 0, 1, 1, 2]
    ev = [1, 2, 3, 2, 3, 3]
    r = rng.random(n)
    for x in range(4, n):
        f = int(r[x] * len(faces))
        a, b, c = faces[f]
        eu += [a, b, c]
        ev += [x, x, x]
        faces[f] = (x, a, b)
        faces += [(x, a, c), (x, b, c)]
    return _from_edges(n, np.asarray(eu, dtype=np.int64), np.asarray(ev, dtype=np.int64))


def _grid_edges(p: int, side: int) -> tuple[list[int], list[int]]:
    u, v = [], []
    base = p * side * side
    for r in range(side):
        for c in range(side):
            i = base + r * side + c
            if c + 1 < side:
                u.append(i), v.append(i + 1)
            if r + 1 < side:
                u.append(i), v.append(i + side)
    return u, v


def tree_of_grids(patches: int, side: int = 10) -> sparse.csr_array:
    """Arbol binario completo (indexacion en monton: hijos 2p+1, 2p+2) de parches side x side abiertos (4 vecinos).
    Hijo izquierdo: su fila superior (fila 0) se une uno a uno (side aristas) con la fila inferior del padre (fila side-1);
    hijo derecho: su fila superior se une con la columna derecha del padre (columna side-1), de arriba abajo."""
    us, vs = [], []
    for p in range(patches):
        u, v = _grid_edges(p, side)
        us += u
        vs += v
    for ch in range(1, patches):
        par = (ch - 1) // 2
        cb, pb = ch * side * side, par * side * side
        for t in range(side):
            top = cb + t  # fila 0 del hijo
            if ch % 2 == 1:
                pv = pb + (side - 1) * side + t  # fila inferior del padre
            else:
                pv = pb + t * side + (side - 1)  # columna derecha del padre
            us.append(top), vs.append(pv)
    return _from_edges(patches * side * side, np.asarray(us, dtype=np.int64), np.asarray(vs, dtype=np.int64))


def barabasi_albert(n: int, m: int, rng: np.random.Generator) -> sparse.csr_array:
    """BA por lista de nodos repetidos: arranque K_{m+1}; cada nodo nuevo elige m destinos distintos de la lista
    (probabilidad proporcional al grado) y anade a la lista los destinos y el nodo m veces."""
    eu, ev, rep = [], [], []
    for a in range(m + 1):
        for b in range(a + 1, m + 1):
            eu.append(a), ev.append(b)
        rep += [a] * m
    for x in range(m + 1, n):
        tg: set[int] = set()
        while len(tg) < m:
            tg.add(rep[int(rng.random() * len(rep))])
        for t in tg:
            eu.append(x), ev.append(t)
            rep.append(t)
        rep += [x] * m
    return _from_edges(n, np.asarray(eu, dtype=np.int64), np.asarray(ev, dtype=np.int64))


def cartesian(a: sparse.csr_array, b: sparse.csr_array) -> sparse.csr_array:
    na, nb = a.shape[0], b.shape[0]
    m = sparse.kron(a, sparse.identity(nb), format="csr") + sparse.kron(sparse.identity(na), b, format="csr")
    m = sparse.csr_array(m)
    m.data[:] = 1.0
    return m


def cycle(n: int) -> sparse.csr_array:
    i = np.arange(n, dtype=np.int64)
    return _from_edges(n, i, (i + 1) % n)


def path(n: int) -> sparse.csr_array:
    i = np.arange(n - 1, dtype=np.int64)
    return _from_edges(n, i, i + 1)


# ------------------------------------------------------------------ panel


def build(spec: dict[str, Any], rng: np.random.Generator) -> sparse.csr_array:
    g, a, n = spec["gen"], spec.get("args", {}), int(spec["N"])
    if g == "tri_torus":
        return triangular_torus(a["side"])
    if g == "honey":
        return honeycomb_torus(a["rows"], a["cols"])
    if g == "lattice":
        return torus_lattice(a["side"], a["dim"])
    if g == "open":
        return lattice_open(a["side"], a["dim"])
    if g == "rgg":
        return rgg(n, a["d"], a["k"], rng)
    if g == "rgg_del":
        return rgg_deleted(n, 3, 12.0, 0.15, rng)
    if g == "diag":
        return torus_with_diagonals(a["side"], 0.2, rng)
    if g == "cyl":
        return rgg_cylinder(n, 12.0, rng)
    if g == "3tree":
        return random_three_tree(n, rng)
    if g == "apol":
        return apollonian(n, rng)
    if g == "grids":
        return tree_of_grids(a["patches"], 10)
    if g == "ba":
        return barabasi_albert(n, 6, rng)
    if g == "sq_short":
        return add_random_shortcuts(torus_lattice(a["side"], 2), 0.01, rng)
    if g == "t3_short":
        return add_random_shortcuts(torus_lattice(a["side"], 3), 0.001, rng)
    if g == "patch":
        return patchwork(n, 3, rng)
    if g == "ladder":
        return cartesian(cycle(a["nc"]), path(4))
    if g == "cylinder":
        return cartesian(cycle(a["nc"]), cycle(10))
    if g == "ring":
        return ring(n, 6)
    if g == "heis":
        return heisenberg(a["hn"])
    raise ValueError(g)


def specs(smoke: bool = False) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    sd = (0,) if smoke else SEEDS
    ns = 2000 if smoke else NS
    det = (0,)

    def add(fam: str, grp: str, gen: str, n: int, seeds: tuple[int, ...] = sd, **args: Any) -> None:
        if smoke and fam not in SMOKE_FAMILIES:
            return
        for s in seeds:
            out.append({"id": f"{fam}_s{s}", "family": fam, "group": grp, "gen": gen, "N": n, "args": args,
                        "key": (MASTER, FAMILY_ID[fam], s)})

    # V
    add("toro_triangular_141", "V", "tri_torus", 45**2 if smoke else 141**2, det, side=45 if smoke else 141)
    add("panal_hexagonal", "V", "honey", 44 * 46 if smoke else 140 * 142, det, rows=44 if smoke else 140, cols=46 if smoke else 142)
    add("T4_12", "V", "lattice", 6**4 if smoke else 12**4, det, side=6 if smoke else 12, dim=4)
    add("RGG4_k16", "V", "rgg", ns, d=4, k=16.0)
    add("RGG2_k20", "V", "rgg", ns, d=2, k=20.0)
    add("RGG3_k12_del15", "V", "rgg_del", ns)
    add("T3_27_diag20", "V", "diag", 13**3 if smoke else 27**3, side=13 if smoke else 27)
    add("cubo_abierto_27", "V", "open", 13**3 if smoke else 27**3, det, side=13 if smoke else 27, dim=3)
    add("RGG2_cilindro_k12", "V", "cyl", ns)
    # X
    add("3arbol", "X", "3tree", ns)
    add("apoloniana", "X", "apol", ns)
    add("arbol_de_rejillas", "X", "grids", 20 * 100 if smoke else 200 * 100, det, patches=20 if smoke else 200)
    add("BA_m6", "X", "ba", ns)
    add("cuadrado_141_atajos1", "X", "sq_short", 45**2 if smoke else 141**2, sd, side=45 if smoke else 141)
    add("T3_27_atajos0.1", "X", "t3_short", 13**3 if smoke else 27**3, sd, side=13 if smoke else 27)
    add("retazos_b3", "X", "patch", ns)
    # U
    add("escalera_C5000xP4", "U", "ladder", ns, det, nc=ns // 4)
    add("cilindro_C2000xC10", "U", "cylinder", ns, det, nc=ns // 10)
    add("anillo_k6", "U", "ring", ns, det)
    # E
    add("Heisenberg_Z27", "E", "heis", 12**3 if smoke else 27**3, det, hn=12 if smoke else 27)
    return out


# ------------------------------------------------------------------ medidas por grafo


def rc2_conditions(row: dict[str, Any]) -> dict[str, bool]:
    k = row["curvature"]
    c = {
        "i_locality": bool(row["locality"] >= LOCALITY_MIN),
        "ii_one_end": bool(row["annulus_status"] == "CONEXO" and len(row["annulus"]["scales"]) >= 2),
        "iii_coherence": bool(row["coherence_status"] in ("COHERENTE", "INTERMEDIO")),
        "iv_curvature": bool(KAPPA_LO <= k["median"] <= KAPPA_HI and k["f_neg"] <= FNEG_MAX),
    }
    c["nondegenerate"] = all(c.values())
    return c


def curvature_guarded(adj: sparse.csr_array, key: tuple[int, ...]) -> dict[str, Any]:
    """curvature_summary con la guarda de viabilidad. Replica el muestreo de aristas de cic_l0b.curvature_summary
    (misma rng, mismo rng.choice) solo para medir el coste; si es viable llama a curvature_summary tal cual."""
    tri = sparse.triu(adj, k=1, format="coo")
    xs, ys = np.asarray(tri.row, dtype=np.int64), np.asarray(tri.col, dtype=np.int64)
    n_edges = 300  # cic_l0b.N_KAPPA_EDGES
    if xs.size > n_edges:
        pick = np.sort(rng_from_key(key).choice(xs.size, size=n_edges, replace=False))
        xs, ys = xs[pick], ys[pick]
    deg = np.diff(adj.indptr)
    worst = int(((deg[xs] + 1) * (deg[ys] + 1)).max())
    if worst > MAX_LP_PRODUCT:
        nan = float("nan")
        return {"n_edges": int(xs.size), "median": nan, "mean": nan, "f_neg": nan, "infeasible": True, "max_lp_product": worst}
    out = curvature_summary(adj, rng_from_key(key))
    out["max_lp_product"] = worst
    return out


_SPECS: list[dict[str, Any]] = []


def run_spec(spec: dict[str, Any]) -> dict[str, Any]:
    t0 = time.perf_counter()
    key = tuple(int(k) for k in spec["key"])
    adj = build(spec, rng_from_key(key))
    t1 = time.perf_counter()
    ann = annulus_profile(adj, rng_from_key(key + (1,)))
    t2 = time.perf_counter()
    loc = locality_fraction(adj, rng_from_key(key + (2,)))
    coh = edge_coherence_profile(adj, rng_from_key(key + (3,)))
    t3 = time.perf_counter()
    kap = curvature_guarded(adj, key + (4,))
    t4 = time.perf_counter()
    row = {k: v for k, v in spec.items() if k not in ("gen", "args", "key")}
    row["key"] = list(key)
    row["n"] = int(adj.shape[0])
    row["mean_degree"] = float(adj.nnz / adj.shape[0])
    row["annulus"] = ann
    row["annulus_status"] = annulus_status(ann)
    row["locality"] = loc
    row["coherence"] = coh
    row["coherence_status"] = coherence_status(coh)
    row["curvature"] = kap
    row["rc2"] = rc2_conditions(row)
    row["rejected_by"] = [c for c in ("i_locality", "ii_one_end", "iii_coherence", "iv_curvature") if not row["rc2"][c]]
    row["seconds"] = {"gen": round(t1 - t0, 2), "annulus": round(t2 - t1, 2), "loc_coh": round(t3 - t2, 2),
                      "kappa": round(t4 - t3, 2)}
    return row


def _run_idx(i: int) -> dict[str, Any]:
    row = run_spec(_SPECS[i])
    print(f"  {row['id']:<30} n={row['n']:<6} deg={row['mean_degree']:.1f} ann={row['annulus_status']:<12} "
          f"sc={len(row['annulus']['scales'])} loc={row['locality']:.2f} coh={row['coherence_status']:<11} "
          f"kmed={row['curvature']['median']:.2f} fneg={row['curvature']['f_neg']:.2f} "
          f"RC2={int(row['rc2']['nondegenerate'])} rej={row['rejected_by']} {row['seconds']}", flush=True)
    return row


# ------------------------------------------------------------------ veredicto


def _frac(rows: list[dict[str, Any]], pred: Any) -> float | None:
    return float(sum(bool(pred(r)) for r in rows) / len(rows)) if rows else None


def verdict(sens: float | None, spec: float | None) -> str:
    if sens is None or spec is None:
        return "NO_EVALUABLE"
    if sens >= 0.8 and spec >= 0.95:
        return "RC2-VALIDA"
    if sens >= 0.7 and spec >= 0.9:
        return "RC2-PARCIAL"
    return "RC2-INVALIDA"


def evaluate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    v = [r for r in rows if r["group"] == "V"]
    neg = [r for r in rows if r["group"] in ("X", "U")]
    sens = _frac(v, lambda r: r["rc2"]["nondegenerate"])
    spec = _frac(neg, lambda r: not r["rc2"]["nondegenerate"])
    e = [r for r in rows if r["group"] == "E"]
    return {"verdict": verdict(sens, spec), "sensitivity": sens, "specificity": spec, "n_V": len(v), "n_X_U": len(neg),
            "false_negatives": [r["id"] for r in v if not r["rc2"]["nondegenerate"]],
            "false_positives": [r["id"] for r in neg if r["rc2"]["nondegenerate"]],
            "E_informative": [{"id": r["id"], "nondegenerate": r["rc2"]["nondegenerate"], "rejected_by": r["rejected_by"]} for r in e]}


def family_table(rows: list[dict[str, Any]]) -> dict[str, Any]:
    fam: dict[str, dict[str, Any]] = {}
    for r in rows:
        f = fam.setdefault(r["family"], {"group": r["group"], "n_instances": 0, "n_nondegenerate": 0, "annulus_status": [],
                                         "n_scales": [], "coherence_status": [], "kappa_median": [], "f_neg": [],
                                         "locality": [], "rejected_by": []})
        f["n_instances"] += 1
        f["n_nondegenerate"] += int(r["rc2"]["nondegenerate"])
        f["annulus_status"].append(r["annulus_status"])
        f["n_scales"].append(len(r["annulus"]["scales"]))
        f["coherence_status"].append(r["coherence_status"])
        f["kappa_median"].append(round(r["curvature"]["median"], 3))
        f["f_neg"].append(round(r["curvature"]["f_neg"], 3))
        f["locality"].append(round(r["locality"], 3))
        f["rejected_by"].append(r["rejected_by"])
    return fam


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--smoke", action="store_true", help="N~2000, semilla 0 -> results/a0_rc2_smoke/")
    args = ap.parse_args(argv)
    out = OUT_SMOKE if args.smoke else OUT_FULL
    _SPECS[:] = specs(args.smoke)
    order = sorted(range(len(_SPECS)), key=lambda i: -int(_SPECS[i]["N"]))
    t0 = time.perf_counter()
    print(f"{len(_SPECS)} grafos, procs={args.procs}", flush=True)
    done: dict[int, dict[str, Any]] = {}
    with mp.get_context("fork").Pool(args.procs) as pool:
        for i, row in zip(order, pool.imap(_run_idx, order, chunksize=1)):
            done[i] = row
    rows = [done[i] for i in range(len(_SPECS))]
    out.mkdir(parents=True, exist_ok=True)
    _atomic_write(out / "graphs.jsonl", "".join(json.dumps(r, default=str) + "\n" for r in rows))
    ev = evaluate(rows)
    summary: dict[str, Any] = {"step": "a0_rc2", "code_commit": head_commit(), "smoke": args.smoke, "n_graphs": len(rows),
                               "seconds": round(time.perf_counter() - t0, 1), "evaluation": ev,
                               "by_family": family_table(rows), "heisenberg_excluded_from_verdicts": True}
    _atomic_write(out / "summary.json", json.dumps(summary, indent=2, default=str))
    print(f"RC-2 {ev['verdict']} sens={ev['sensitivity']} spec={ev['specificity']} FN={ev['false_negatives']} FP={ev['false_positives']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
