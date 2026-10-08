"""Clasificacion del estado final de C0 (prerregistro §4, C0-L4) y validacion del criterio de localidad."""

from __future__ import annotations

from typing import Any

import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, shortest_path

from omega.c0 import references as R
from omega.landscape.structure import (
    CLASS_MULTI_CLIQUE,
    CLASS_OVERLAPPING,
    CLASS_SINGLE_CLIQUE,
    CLASS_UNIFORM,
    classify_structure,
)
from omega.types import BoolArray, FloatArray, IntArray

__all__ = [
    "strong_support",
    "giant_component",
    "degree_preserving_rewire",
    "mean_hop_connected",
    "h_null",
    "short_cycle_fraction",
    "classify_c0",
    "validate_locality",
]

EMPTY_MAX = 1e-6
DENSE_CLASSES = (CLASS_UNIFORM, CLASS_SINGLE_CLIQUE, CLASS_MULTI_CLIQUE, CLASS_OVERLAPPING)
DENSITY_TRIVIAL = 0.25
GIANT_MIN_FRAC = 0.5
H_NULL_MIN = 1.15
CYCLES_MIN = 0.5
N_NULL = 5


def _as_bool(adj: BoolArray | FloatArray) -> BoolArray:
    return np.asarray(np.asarray(adj, dtype=np.float64) > 0.0, dtype=np.bool_)


def strong_support(w: FloatArray, rel: float = 0.1) -> BoolArray:
    """A = {W > rel * max W} sin diagonal."""
    mx = float(np.max(w)) if w.size else 0.0
    a = w > rel * mx
    np.fill_diagonal(a, False)
    return np.asarray(a, dtype=np.bool_)


def giant_component(adj: BoolArray | FloatArray) -> IntArray:
    """Indices de la componente conexa mayor (la de menor indice minimo si hay empate)."""
    _, labels = connected_components(csr_matrix(_as_bool(adj)), directed=False)
    cnt = np.bincount(labels)
    big = int(np.argmax(cnt))
    return np.asarray(np.flatnonzero(labels == big), dtype=np.int64)


def degree_preserving_rewire(adj: BoolArray | FloatArray, n_swaps: int, rng: np.random.Generator) -> BoolArray:
    """Intercambios dobles de aristas (sin lazos ni multiaristas); determinista dado rng."""
    a = _as_bool(adj)
    n = a.shape[0]
    g = nx.from_numpy_array(a.astype(np.int64))
    if n_swaps > 0 and g.number_of_edges() >= 2 and n >= 4:
        try:
            nx.double_edge_swap(g, nswap=n_swaps, max_tries=100 * n_swaps, seed=int(rng.integers(2**31)))
        except nx.NetworkXError:
            pass  # grafo casi rigido: se queda con los intercambios logrados
    out = np.asarray(nx.to_numpy_array(g, nodelist=range(n), dtype=np.float64)) > 0
    return np.asarray(out, dtype=np.bool_)


def mean_hop_connected(adj: BoolArray) -> float:
    """Distancia media de saltos entre pares i<j conectados (nan si no hay)."""
    n = adj.shape[0]
    if n < 2:
        return float("nan")
    d = shortest_path(csr_matrix(adj), unweighted=True, directed=False)
    v = d[np.triu_indices(n, k=1)]
    v = v[np.isfinite(v)]
    return float(v.mean()) if v.size else float("nan")


def h_null(adj: BoolArray | FloatArray, rng: np.random.Generator, n_null: int = N_NULL) -> float:
    """Salto medio de la gigante / media sobre n_null recableados (10|E| swaps) de la gigante inducida."""
    a = _as_bool(adj)
    idx = giant_component(a)
    sub = a[np.ix_(idx, idx)]
    n_e = int(np.triu(sub, k=1).sum())
    if idx.size < 4 or n_e < 2:
        return float("nan")
    h0 = mean_hop_connected(sub)
    nulls = [mean_hop_connected(degree_preserving_rewire(sub, 10 * n_e, rng)) for _ in range(n_null)]
    hn = float(np.nanmean(nulls))
    return h0 / hn if hn > 0.0 and np.isfinite(hn) else float("nan")


def short_cycle_fraction(adj: BoolArray | FloatArray) -> float:
    """Fraccion de aristas de la gigante contenidas en un triangulo o un 4-ciclo."""
    a = _as_bool(adj)
    idx = giant_component(a)
    s = a[np.ix_(idx, idx)]
    if s.shape[0] < 3:
        return 0.0
    f = s.astype(np.float64)
    iu = np.triu_indices(s.shape[0], k=1)
    edge = s[iu]
    if not edge.any():
        return 0.0
    f2 = f @ f
    deg = f.sum(axis=1)
    # caminos i-u-v-j de longitud 3 menos los degenerados (u=j o v=i): deg_i + deg_j - 1 en una arista
    walks3 = (f2 @ f)[iu]
    tri = f2[iu] > 0.5
    c4 = (walks3 - (deg[iu[0]] + deg[iu[1]] - 1.0)) > 0.5
    return float(np.mean((tri | c4)[edge]))


def _mean_degree(a: BoolArray, nodes: IntArray) -> float:
    return float(a[np.ix_(nodes, nodes)].sum(axis=1).mean()) if nodes.size else 0.0


def classify_c0(w: FloatArray, rng: np.random.Generator) -> dict[str, Any]:
    """Clase del estado final segun el orden de decision de C0-L4 mas diagnosticos."""
    n = int(w.shape[0])
    out: dict[str, Any] = {
        "class": "VACIO", "structure": None, "density_A": None, "giant_frac": None, "mean_degree_A": None,
        "kmax_over_kmean": None, "h_null": None, "short_cycles": None, "n_edges_A": 0,
    }
    if float(np.max(w)) <= EMPTY_MAX:
        return out
    a = strong_support(w)
    deg = a.sum(axis=1)
    live = deg > 0
    n_live = int(live.sum())
    n_edges = int(a.sum() // 2)
    giant = giant_component(a)
    out["n_edges_A"] = n_edges
    out["giant_frac"] = float(giant.size) / n
    out["density_A"] = float(a.sum()) / (n_live * (n_live - 1)) if n_live > 1 else 0.0
    out["mean_degree_A"] = float(deg[live].mean()) if n_live else 0.0
    out["kmax_over_kmean"] = float(deg.max() / deg[live].mean()) if n_live else None
    st = classify_structure(w).state_class
    out["structure"] = st
    if st in DENSE_CLASSES or out["density_A"] >= DENSITY_TRIVIAL:
        out["class"] = "DENSO_TRIVIAL"
        return out
    if giant.size < GIANT_MIN_FRAC * n:
        out["class"] = "FRAGMENTADO"
        return out
    hn = h_null(a, rng)
    sc = short_cycle_fraction(a)
    out["h_null"] = hn
    out["short_cycles"] = sc
    out["class"] = "DISPERSO_LOCAL" if (np.isfinite(hn) and hn >= H_NULL_MIN and sc >= CYCLES_MIN) else "DISPERSO_NO_LOCAL"
    return out


def validate_locality(seed: int = 0) -> dict[str, Any]:
    """Validacion preregistrada: 4 positivos deben dar LOCAL, 5 negativos no."""
    rng = np.random.Generator(np.random.PCG64(seed))
    n = 216
    graphs: list[tuple[str, bool, FloatArray]] = [
        ("T3_6x6x6", True, R.torus_lattice_3d(6)),
        ("RGG3_k12", True, R.rgg3_torus_binary(n, 12.0, rng)),
        ("triangular_12x18", True, R.triangular_torus(12, 18)),
        ("ring_k6", True, R.ring_lattice(n, 6)),
        ("ER_k6", False, R.erdos_renyi_m(n, n * 6 // 2, rng)),
        ("RR_k6", False, R.random_regular(n, 6, rng)),
        ("K_108_108", False, R.complete_bipartite_half(n)),
        ("K8_union", False, R.clique_union(n, 8)),
        ("random_tree", False, R.random_tree(n, rng)),
    ]
    rows: dict[str, Any] = {}
    ok = True
    for name, positive, adj in graphs:
        r = classify_c0(adj, rng)
        good = (r["class"] == "DISPERSO_LOCAL") == positive
        ok = ok and good
        rows[name] = {**r, "expected_local": positive, "ok": good}
    return {"seed": seed, "ok": ok, "graphs": rows}
