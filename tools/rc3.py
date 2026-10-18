"""OMEGA RC-3: juez de exclusion A- (docs/OMEGA_RC3_PRERREGISTRO.md §1-§3, congelado).

Salida: results/rc3/{pairs.jsonl, summary.json} (smoke: results/rc3_smoke). Uso:
  python tools/rc3.py [--procs 3] [--smoke] [--fresh]
Cada instancia es un PAR (G_N, G_8N) del mismo proceso; clave (MASTER, FAMILY_ID[fam], semilla, indice_de_tamano).
RNG: generador = key; anillos = key+(1,); coherencia = key+(3,); perfil de bolas = key+(5,).
Es REANUDABLE por par: pairs.jsonl se anade por par terminado; al relanzar se saltan los pares ya presentes
(--fresh los borra). Guardas de viabilidad: un instrumento que lance MemoryError se declara NO_EVALUABLE
(abstencion: no activa ninguna exclusion) y se anota en `guards`.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import itertools  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402
from scipy import sparse  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from a0_rc2 import (  # noqa: E402
    add_random_shortcuts,
    apollonian,
    barabasi_albert,
    cartesian,
    cycle,
    giant,
    honeycomb_torus,
    lattice_open,
    path,
)
from c0_dynamics import _atomic_write  # noqa: E402
from cic_l0b import random_two_tree  # noqa: E402
from coh_l0b import cycle_times_rr, heisenberg, husimi_cactus, prufer_tree, subdivided_binary_tree  # noqa: E402
from p1d3_panel import (  # noqa: E402
    _codes,
    _from_edges,
    _rewire_one_end,
    erdos_renyi,
    patchwork,
    random_regular,
    rgg,
    torus_lattice,
)

from omega.c0.references import rng_from_key  # noqa: E402
from omega.diagnostics.annulus import annulus_profile, annulus_status  # noqa: E402
from omega.diagnostics.coherence import coherence_status, edge_coherence_profile  # noqa: E402
from omega.diagnostics.sampled_growth import sampled_ball_profile  # noqa: E402
from omega.experiments.v11.gate import head_commit  # noqa: E402

MASTER = 20261016
OUT_FULL = ROOT / "results" / "rc3"
OUT_SMOKE = ROOT / "results" / "rc3_smoke"
DELTA_LOW, DELTA_HIGH = 0.125, 0.75  # X1 (delta < 0.125), X3 (delta > 0.75): congelados

# (id, grupo, d, nueva, aleatoria)
FAMS: list[tuple[str, str, int | None, bool, bool]] = [
    ("toro_cuadrado", "V", 2, False, False), ("RGG2_k8", "V", 2, False, True),
    ("triangular_caja", "V", 2, True, False), ("panal_hexagonal", "V", 2, False, False),
    ("T3", "V", 3, False, False), ("RGG3_k8", "V", 3, False, True), ("FCC", "V", 3, True, False),
    ("BCC", "V", 3, True, False), ("RGG3_k12_caja", "V", 3, False, True), ("T3_diag20", "V", 3, False, True),
    ("T4", "V", 4, False, False), ("RGG4_k8", "V", 4, False, True), ("T4_diag20", "V", 4, True, True),
    ("RGG4_k16", "V", 4, False, True),
    ("T5", "V", 5, True, False), ("RGG5_k16", "V", 5, True, True),
    ("T6", "V", 6, True, False),
    ("Heisenberg", "E", None, False, False),
    ("clique", "X", None, True, False), ("4arbol", "X", None, True, True), ("apoloniana", "X", None, False, True),
    ("BA_m3", "X", None, True, True), ("RR_k3", "X", None, True, True), ("ER_k4", "X", None, True, True),
    ("arbol_prufer", "X", None, False, True), ("arbol_binario_sub10", "X", None, False, False),
    ("arbol_cubos_5", "X", None, True, False), ("arbol_binario_x_C10", "X", None, True, False),
    ("RGG3_k12_atajos0.1pc", "X", None, True, True), ("cuadrado_WS_b0.01", "X", None, True, True),
    ("retazos_b3", "X", None, False, True), ("cactus", "X", None, False, False), ("2arbol", "X", None, False, True),
    ("cilindro_CxC5", "U", None, True, False), ("RGG1_anillo_k10", "U", None, True, True),
    ("escalera_CxP4", "U", None, False, False), ("tubo_CxRR", "U", None, False, True),
]
FAMILY_ID = {f[0]: i + 1 for i, f in enumerate(FAMS)}
INFO = {f[0]: {"group": f[1], "d": f[2], "new": f[3], "random": f[4]} for f in FAMS}

# ------------------------------------------------------------------ generadores nuevos


def k_tree(n: int, k: int, rng: np.random.Generator) -> sparse.csr_array:
    """k-arbol aleatorio: K_{k+1} inicial; cada nodo nuevo elige uniformemente una k-clique de la lista (la original sigue
    disponible), se une a sus k vertices y se anaden las k cliques (c con un vertice sustituido por el nodo nuevo)."""
    cl = [tuple(c) for c in itertools.combinations(range(k + 1), k)]
    eu = [a for a, _ in itertools.combinations(range(k + 1), 2)]
    ev = [b for _, b in itertools.combinations(range(k + 1), 2)]
    r = rng.random(n)
    for x in range(k + 1, n):
        c = cl[int(r[x] * len(cl))]
        for a in c:
            eu.append(a), ev.append(x)
        for i in range(k):
            cl.append(c[:i] + (x,) + c[i + 1:])
    return _from_edges(n, np.asarray(eu, dtype=np.int64), np.asarray(ev, dtype=np.int64))


def triangular_open(side: int) -> sparse.csr_array:
    """Triangular en caja abierta: (i,j)~(i+1,j),(i,j+1),(i+1,j-1) SIN envoltura."""
    i, j = np.meshgrid(np.arange(side), np.arange(side), indexing="ij")
    i, j = i.ravel().astype(np.int64), j.ravel().astype(np.int64)
    idx = i * side + j
    us, vs = [], []
    for di, dj in ((1, 0), (0, 1), (1, -1)):
        ok = (i + di < side) & (j + dj >= 0) & (j + dj < side)
        us.append(idx[ok])
        vs.append((i[ok] + di) * side + j[ok] + dj)
    return _from_edges(side * side, np.concatenate(us), np.concatenate(vs))


def _cubic_offsets_lattice(S: int, keep: Any, offs: np.ndarray) -> sparse.csr_array:
    """Sitios de la caja periodica S^3 con keep(x,y,z); aristas por `offs` (mod S); reetiquetado 0..n-1."""
    g = np.stack(np.meshgrid(*(np.arange(S),) * 3, indexing="ij"), axis=-1).reshape(-1, 3).astype(np.int64)
    sel = keep(g[:, 0], g[:, 1], g[:, 2])
    lab = -np.ones(S**3, dtype=np.int64)
    lab[np.flatnonzero(sel)] = np.arange(int(sel.sum()))
    c = g[sel]
    src = lab[(c[:, 0] * S + c[:, 1]) * S + c[:, 2]]
    us, vs = [], []
    for o in offs:
        t = (c + o) % S
        dst = lab[(t[:, 0] * S + t[:, 1]) * S + t[:, 2]]
        assert (dst >= 0).all()
        us.append(src), vs.append(dst)
    return _from_edges(int(sel.sum()), np.concatenate(us), np.concatenate(vs))


def fcc_lattice(L: int) -> sparse.csr_array:
    """FCC periodica: puntos de la caja cubica de lado 2L con suma de coordenadas par (4 L^3 nodos); 12 vecinos (±1,±1,0) y perms."""
    offs = np.array([p for p in itertools.product((-1, 0, 1), repeat=3) if sum(abs(t) for t in p) == 2], dtype=np.int64)
    assert offs.shape[0] == 12
    return _cubic_offsets_lattice(2 * L, lambda x, y, z: (x + y + z) % 2 == 0, offs)


def bcc_lattice(L: int) -> sparse.csr_array:
    """BCC periodica: caja entera de lado 2L, sitios con las tres coordenadas de la misma paridad (2 L^3 nodos);
    cada sitio unido a los 8 vecinos (±1,±1,±1) (la subred opuesta; = (±1/2)^3 en unidades de la celda)."""
    offs = np.array(list(itertools.product((-1, 1), repeat=3)), dtype=np.int64)
    return _cubic_offsets_lattice(2 * L, lambda x, y, z: ((x % 2) == (y % 2)) & ((y % 2) == (z % 2)), offs)


def torus_with_diagonals(side: int, dim: int, frac: float, rng: np.random.Generator) -> sparse.csr_array:
    """T^dim periodico + round(frac*|E|) diagonales (offsets con exactamente dos coordenadas ±1: 4*C(dim,2)) distintas y nuevas."""
    base = torus_lattice(side, dim)
    n = side**dim
    t0 = sparse.triu(base, k=1, format="coo")
    have = _codes(n, np.asarray(t0.row, dtype=np.int64), np.asarray(t0.col, dtype=np.int64))
    target = int(round(frac * have.size))
    offl = []
    for i, j in itertools.combinations(range(dim), 2):
        for a in (-1, 1):
            for b in (-1, 1):
                o = [0] * dim
                o[i], o[j] = a, b
                offl.append(o)
    offs = np.array(offl, dtype=np.int64)
    new = np.empty(0, dtype=np.int64)
    while new.size < target:
        m = 2 * (target - new.size) + 16
        x = rng.integers(n, size=m)
        o = offs[rng.integers(offs.shape[0], size=m)]
        c = np.stack(np.unravel_index(x, (side,) * dim), axis=1)
        y = np.ravel_multi_index(tuple(((c + o) % side).T), (side,) * dim).astype(np.int64)
        cc = np.unique(_codes(n, x, y))
        cc = cc[~np.isin(cc, have)]
        new = np.unique(np.concatenate([new, cc]))
    new = rng.permutation(new)[:target]
    allc = np.concatenate([have, new])
    return _from_edges(n, allc // n, allc % n)


def clique_graph(n: int) -> sparse.csr_array:
    a = sparse.csr_array(np.ones((n, n)) - np.eye(n))
    return a


def tree_of_cubes(patches: int, s: int = 5) -> sparse.csr_array:
    """Arbol binario completo (monton: hijos 2p+1, 2p+2) de parches s^3 abiertos. Cara del hijo x=0 (s^2 vertices (a,b)) pegada
    vertice a vertice (s^2 aristas) al padre: hijo izquierdo -> cara x=s-1 del padre ((s-1,a,b)); derecho -> cara y=s-1 ((a,s-1,b))."""
    loc = sparse.triu(lattice_open(s, 3), k=1, format="coo")
    lu, lv = np.asarray(loc.row, dtype=np.int64), np.asarray(loc.col, dtype=np.int64)
    base = (np.arange(patches, dtype=np.int64) * s**3)[:, None]
    us = [(base + lu[None, :]).ravel()]
    vs = [(base + lv[None, :]).ravel()]
    a, b = (t.ravel().astype(np.int64) for t in np.meshgrid(np.arange(s), np.arange(s), indexing="ij"))
    ch = np.arange(1, patches, dtype=np.int64)
    par = (ch - 1) // 2
    cidx = ch[:, None] * s**3 + (0 * s + a[None, :]) * s + b[None, :]
    pl = par[:, None] * s**3 + ((s - 1) * s + a[None, :]) * s + b[None, :]
    pr = par[:, None] * s**3 + (a[None, :] * s + (s - 1)) * s + b[None, :]
    pidx = np.where((ch % 2 == 1)[:, None], pl, pr)
    us.append(cidx.ravel()), vs.append(pidx.ravel())
    return _from_edges(patches * s**3, np.concatenate(us), np.concatenate(vs))


def binary_tree(nt: int) -> sparse.csr_array:
    c = np.arange(1, nt, dtype=np.int64)
    return _from_edges(nt, (c - 1) // 2, c)


def watts_strogatz_square(side: int, beta: float, rng: np.random.Generator) -> sparse.csr_array:
    """Toro cuadrado; cada arista (u<v) reconecta su extremo v a w uniforme con prob beta (w != u y (u,w) inexistente)."""
    n = side * side
    t = sparse.triu(torus_lattice(side, 2), k=1, format="coo")
    u, v = np.asarray(t.row, dtype=np.int64), np.asarray(t.col, dtype=np.int64)
    return _rewire_one_end(n, u, v, rng.random(u.size) < beta, rng)


def random_ring(n: int, k: float, rng: np.random.Generator) -> sparse.csr_array:
    return rgg(n, 1, k, rng, True)


# ------------------------------------------------------------------ panel (tamanos documentados)
# fam -> (args_N, args_8N) completos; (args_N, args_8N) smoke. Cada args = dict con "gen" y parametros.


def _sizes(smoke: bool) -> dict[str, tuple[dict[str, Any], dict[str, Any]]]:
    def rg(d: int, k: float, n0: int, **kw: Any) -> tuple[dict[str, Any], dict[str, Any]]:
        return ({"gen": "rgg", "N": n0, "d": d, "k": k, **kw}, {"gen": "rgg", "N": 8 * n0, "d": d, "k": k, **kw})

    def lat(g: str, a: int, b: int, **kw: Any) -> tuple[dict[str, Any], dict[str, Any]]:
        return ({"gen": g, "side": a, **kw}, {"gen": g, "side": b, **kw})

    s = smoke
    n0 = 1000 if s else 10_000
    S: dict[str, tuple[dict[str, Any], dict[str, Any]]] = {
        "toro_cuadrado": lat("lattice", 32 if s else 100, 90 if s else 283, dim=2),
        "RGG2_k8": rg(2, 8.0, n0),
        "triangular_caja": lat("tri_open", 32 if s else 100, 90 if s else 283),
        "panal_hexagonal": ({"gen": "honey", "rows": 32 if s else 100, "cols": 32 if s else 100},
                            {"gen": "honey", "rows": 90 if s else 282, "cols": 90 if s else 284}),
        "T3": lat("lattice", 10 if s else 22, 20 if s else 44, dim=3),
        "RGG3_k8": rg(3, 8.0, n0),
        "FCC": ({"gen": "fcc", "L": 6 if s else 14}, {"gen": "fcc", "L": 12 if s else 27}),
        "BCC": ({"gen": "bcc", "L": 8 if s else 17}, {"gen": "bcc", "L": 16 if s else 34}),
        "RGG3_k12_caja": rg(3, 12.0, n0, periodic=False),
        "T3_diag20": lat("diag", 10 if s else 22, 20 if s else 44, dim=3),
        "T4": lat("lattice", 6 if s else 10, 10 if s else 17, dim=4),
        "RGG4_k8": rg(4, 8.0, n0),
        "T4_diag20": lat("diag", 6 if s else 10, 10 if s else 17, dim=4),
        "RGG4_k16": rg(4, 16.0, n0),
        "T5": lat("lattice", 4 if s else 6, 6 if s else 11, dim=5),
        "RGG5_k16": rg(5, 16.0, n0),
        "T6": lat("lattice", 3 if s else 4, 4 if s else 7, dim=6),
        "Heisenberg": ({"gen": "heis", "hn": 10 if s else 22}, {"gen": "heis", "hn": 20 if s else 44}),
        "clique": ({"gen": "clique", "N": 100 if s else 200}, {"gen": "clique", "N": 800 if s else 1600}),
        "4arbol": ({"gen": "ktree", "N": n0, "k": 4}, {"gen": "ktree", "N": 8 * n0, "k": 4}),
        "apoloniana": ({"gen": "apol", "N": n0}, {"gen": "apol", "N": 8 * n0}),
        "BA_m3": ({"gen": "ba", "N": n0}, {"gen": "ba", "N": 8 * n0}),
        "RR_k3": ({"gen": "rr", "N": n0}, {"gen": "rr", "N": 8 * n0}),
        "ER_k4": ({"gen": "er", "N": n0}, {"gen": "er", "N": 8 * n0}),
        "arbol_prufer": ({"gen": "prufer", "N": n0}, {"gen": "prufer", "N": 8 * n0}),
        "arbol_binario_sub10": ({"gen": "subtree", "levels": 7 if s else 10}, {"gen": "subtree", "levels": 10 if s else 13}),
        "arbol_cubos_5": ({"gen": "cubes", "patches": 8 if s else 80}, {"gen": "cubes", "patches": 64 if s else 640}),
        "arbol_binario_x_C10": ({"gen": "treecyc", "nt": n0 // 10}, {"gen": "treecyc", "nt": 8 * n0 // 10}),
        "RGG3_k12_atajos0.1pc": rg(3, 12.0, n0, shortcut=0.001),
        "cuadrado_WS_b0.01": lat("ws_sq", 32 if s else 100, 90 if s else 283),
        "retazos_b3": ({"gen": "patch", "N": n0}, {"gen": "patch", "N": 8 * n0}),
        "cactus": ({"gen": "cactus", "N": n0}, {"gen": "cactus", "N": 8 * n0}),
        "2arbol": ({"gen": "2tree", "N": n0}, {"gen": "2tree", "N": 8 * n0}),
        "cilindro_CxC5": ({"gen": "cyl", "nc": n0 // 5}, {"gen": "cyl", "nc": 8 * n0 // 5}),
        "RGG1_anillo_k10": ({"gen": "ring_rgg", "N": n0}, {"gen": "ring_rgg", "N": 8 * n0}),
        "escalera_CxP4": ({"gen": "ladder", "nc": n0 // 4}, {"gen": "ladder", "nc": 8 * n0 // 4}),
        "tubo_CxRR": ({"gen": "tube", "nc": 10 if s else 50, "nr": 100 if s else 200},
                      {"gen": "tube", "nc": 80 if s else 400, "nr": 100 if s else 200}),
    }
    assert set(S) == set(FAMILY_ID)
    return S


def build(a: dict[str, Any], rng: np.random.Generator) -> sparse.csr_array:
    g = a["gen"]
    n = int(a.get("N", 0))
    if g == "lattice":
        return torus_lattice(a["side"], a["dim"])
    if g == "rgg":
        adj = rgg(n, a["d"], a["k"], rng, a.get("periodic", True))
        return add_random_shortcuts(adj, a["shortcut"], rng) if "shortcut" in a else adj
    if g == "tri_open":
        return triangular_open(a["side"])
    if g == "honey":
        return honeycomb_torus(a["rows"], a["cols"])
    if g == "fcc":
        return fcc_lattice(a["L"])
    if g == "bcc":
        return bcc_lattice(a["L"])
    if g == "diag":
        return torus_with_diagonals(a["side"], a["dim"], 0.2, rng)
    if g == "heis":
        return heisenberg(a["hn"])
    if g == "clique":
        return clique_graph(n)
    if g == "ktree":
        return k_tree(n, a["k"], rng)
    if g == "apol":
        return apollonian(n, rng)
    if g == "ba":
        return barabasi_albert(n, 3, rng)
    if g == "rr":
        return random_regular(n, 3, rng)
    if g == "er":
        return erdos_renyi(n, 2 * n, rng)  # grado medio 4; se mide en la gigante
    if g == "prufer":
        return prufer_tree(n, rng)
    if g == "subtree":
        return subdivided_binary_tree(a["levels"], 10)
    if g == "cubes":
        return tree_of_cubes(a["patches"], 5)
    if g == "treecyc":
        return cartesian(binary_tree(a["nt"]), cycle(10))
    if g == "ws_sq":
        return watts_strogatz_square(a["side"], 0.01, rng)
    if g == "patch":
        return patchwork(n, 3, rng)
    if g == "cactus":
        return husimi_cactus(n)
    if g == "2tree":
        return random_two_tree(n, rng)
    if g == "cyl":
        return cartesian(cycle(a["nc"]), cycle(5))
    if g == "ring_rgg":
        return random_ring(n, 10.0, rng)
    if g == "ladder":
        return cartesian(cycle(a["nc"]), path(4))
    if g == "tube":
        return cycle_times_rr(a["nc"], a["nr"], 4, rng)
    raise ValueError(g)


def specs(smoke: bool = False) -> list[dict[str, Any]]:
    """Un spec por PAR (familia, semilla): aleatorias semillas 0,1; deterministas 0."""
    S = _sizes(smoke)
    out = []
    for fam, grp, d, new, rnd in FAMS:
        for seed in ((0, 1) if rnd else (0,)):
            out.append({"id": f"{fam}_s{seed}", "family": fam, "group": grp, "d": d, "new": new, "seed": seed,
                        "args": [S[fam][0], S[fam][1]]})
    return out


# ------------------------------------------------------------------ estadisticas


def half_mass_radius(profile: dict[str, Any]) -> float:
    """R (§1.1): m(0)=1, m(r)=profile['mean'][r-1]; primer r con m(r) >= n_giant/2; interpolacion lineal entre r-1 y r.
    Si no cruza dentro del perfil devuelto, R = r_cap."""
    n = float(profile["n_giant"])
    half = n / 2.0
    m = [1.0] + [float(x) for x in profile["mean"]]
    if m[0] >= half:
        return 0.0
    for r in range(1, len(m)):
        if m[r] >= half:
            return float((r - 1) + (half - m[r - 1]) / (m[r] - m[r - 1]))
    return float(profile["r_cap"])


def scaling_delta(r1: float | None, r8: float | None, n1: int, n8: int) -> float | None:
    if r1 is None or r8 is None or r1 <= 0 or r8 <= 0 or n1 <= 0 or n8 <= 1 or n8 == n1:
        return None
    return float(math.log(r8 / r1) / math.log(n8 / n1))


def pair_exclusions(delta: float | None, ann: tuple[str, str], coh: tuple[str, str]) -> list[str]:
    """Reglas X1-X4 (§1). delta None (no evaluable) o estados NO_EVALUABLE/SIN_VENTANA no activan nada (abstencion)."""
    ex = []
    ok = delta is not None and not math.isnan(delta)
    if ok and delta < DELTA_LOW:
        ex.append("X1")
    if ann[0] == "RAMIFICADO" and ann[1] == "RAMIFICADO":
        ex.append("X2")
    if (ann[0] == "DOS_EXTREMOS" and ann[1] == "DOS_EXTREMOS") or (ok and delta > DELTA_HIGH):
        ex.append("X3")
    if coh[0] == "INCOHERENTE" and coh[1] == "INCOHERENTE":
        ex.append("X4")
    return ex


def _guarded(name: str, guards: list[str], fn: Any, *a: Any) -> Any:
    try:
        return fn(*a)
    except MemoryError:
        guards.append(f"{name}:MemoryError")
        return None


def measure(adj: sparse.csr_array, key: tuple[int, ...]) -> dict[str, Any]:
    t0 = time.perf_counter()
    guards: list[str] = []
    prof = _guarded("ball", guards, sampled_ball_profile, adj, rng_from_key(key + (5,)))
    t1 = time.perf_counter()
    ann = _guarded("annulus", guards, annulus_profile, adj, rng_from_key(key + (1,)))
    t2 = time.perf_counter()
    coh = _guarded("coherence", guards, edge_coherence_profile, adj, rng_from_key(key + (3,)))
    t3 = time.perf_counter()
    return {
        "n": int(adj.shape[0]), "mean_degree": float(adj.nnz / adj.shape[0]),
        "n_giant": None if prof is None else int(prof["n_giant"]),
        "R": None if prof is None else half_mass_radius(prof), "r_cap": None if prof is None else int(prof["r_cap"]),
        "annulus": ann, "annulus_status": "NO_EVALUABLE" if ann is None else annulus_status(ann),
        "n_scales": None if ann is None else len(ann["scales"]), "r_w": None if ann is None else ann["r_w"],
        "coherence": coh, "coherence_status": "NO_EVALUABLE" if coh is None else coherence_status(coh),
        "guards": guards,
        "seconds": {"ball": round(t1 - t0, 2), "annulus": round(t2 - t1, 2), "coherence": round(t3 - t2, 2)},
    }


def run_pair(spec: dict[str, Any]) -> dict[str, Any]:
    t0 = time.perf_counter()
    fid = FAMILY_ID[spec["family"]]
    sizes = []
    for si in (0, 1):
        key = (MASTER, fid, int(spec["seed"]), si)
        tg = time.perf_counter()
        adj = build(spec["args"][si], rng_from_key(key))
        gen_s = round(time.perf_counter() - tg, 2)
        m = measure(adj, key)
        m["seconds"]["gen"] = gen_s
        m["key"] = list(key)
        sizes.append(m)
        del adj
    a, b = sizes
    delta = scaling_delta(a["R"], b["R"], a["n_giant"] or 0, b["n_giant"] or 0)
    ex = pair_exclusions(delta, (a["annulus_status"], b["annulus_status"]), (a["coherence_status"], b["coherence_status"]))
    row = {k: v for k, v in spec.items() if k != "args"}
    row.update({"sizes": sizes, "delta": delta, "exclusions": ex,
                "verdict": ("EXCLUIDO(" + ",".join(ex) + ")") if ex else "NO-EXCLUIDO",
                "excluded": bool(ex), "guards": a["guards"] + b["guards"], "seconds": round(time.perf_counter() - t0, 1)})
    return row


_SPECS: list[dict[str, Any]] = []


def _run_idx(i: int) -> dict[str, Any]:
    return run_pair(_SPECS[i])


def fmt(r: dict[str, Any]) -> str:
    a, b = r["sizes"]
    d = "None" if r["delta"] is None else f"{r['delta']:.3f}"
    return (f"  {r['id']:<30} n={a['n_giant']}/{b['n_giant']} R={a['R']:.2f}/{b['R']:.2f} delta={d} "
            f"ann={a['annulus_status']}/{b['annulus_status']} sc8={b['n_scales']} coh={a['coherence_status']}/{b['coherence_status']} "
            f"{r['verdict']} guards={r['guards']} {r['seconds']}s")


# ------------------------------------------------------------------ metricas y veredicto (§3)


def metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    v = [r for r in rows if r["group"] == "V"]  # Heisenberg (E) no entra en FE
    neg = [r for r in rows if r["group"] in ("X", "U")]
    by_d: dict[str, dict[str, Any]] = {}
    for r in v:
        e = by_d.setdefault(str(r["d"]), {"n": 0, "excluded": 0})
        e["n"] += 1
        e["excluded"] += int(r["excluded"])
    for e in by_d.values():
        e["FE"] = e["excluded"] / e["n"]
    by_fam: dict[str, dict[str, Any]] = {}
    for r in neg:
        e = by_fam.setdefault(r["family"], {"n": 0, "excluded": 0, "non_excluded": []})
        e["n"] += 1
        e["excluded"] += int(r["excluded"])
        if not r["excluded"]:
            e["non_excluded"].append(r["id"])
    for e in by_fam.values():
        e["EP"] = e["excluded"] / e["n"]
    fe = sum(r["excluded"] for r in v) / len(v) if v else None
    ep = sum(r["excluded"] for r in neg) / len(neg) if neg else None
    cond = None
    if fe is not None and ep is not None:
        cond = {"FE<=0.05": fe <= 0.05, "per_d_max1": all(e["excluded"] <= 1 for e in by_d.values()),
                "EP>=0.90": ep >= 0.90, "families_all_but_one": all(len(e["non_excluded"]) <= 1 for e in by_fam.values())}
    return {"n_V": len(v), "n_X_U": len(neg), "FE": fe, "FE_by_d": by_d, "EP": ep, "EP_by_family": by_fam,
            "valid_conditions": cond, "verdict": rc3_verdict(fe, ep, cond),
            "V_excluded": [r["id"] for r in v if r["excluded"]], "XU_not_excluded": [r["id"] for r in neg if not r["excluded"]]}


def rc3_verdict(fe: float | None, ep: float | None, cond: dict[str, bool] | None) -> str:
    if fe is None or ep is None or cond is None:
        return "NO_EVALUABLE"
    if all(cond.values()):
        return "RC3-VALIDA"
    if fe <= 0.10 and ep >= 0.75:
        return "RC3-PARCIAL"
    return "RC3-INVALIDA"


def evaluate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    heis = [r for r in rows if r["group"] == "E"]
    return {"new": metrics([r for r in rows if r["new"]]), "repeated": metrics([r for r in rows if not r["new"]]),
            "all": metrics(rows),
            "Heisenberg_V_ext": [{"id": r["id"], "verdict": r["verdict"], "delta": r["delta"],
                                  "coherence": [s["coherence_status"] for s in r["sizes"]]} for r in heis]}


def family_table(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for r in rows:
        a, b = r["sizes"]
        out.append({"id": r["id"], "family": r["family"], "group": r["group"], "d": r["d"], "new": r["new"],
                    "n_giant": [a["n_giant"], b["n_giant"]], "R": [a["R"], b["R"]], "delta": r["delta"],
                    "annulus": [a["annulus_status"], b["annulus_status"]], "n_scales_8N": b["n_scales"],
                    "r_w": [a["r_w"], b["r_w"]], "coherence": [a["coherence_status"], b["coherence_status"]],
                    "exclusions": r["exclusions"], "verdict": r["verdict"]})
    return out


def _load(path: Path) -> dict[str, dict[str, Any]]:
    done: dict[str, dict[str, Any]] = {}
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                try:
                    r = json.loads(line)
                    done[r["id"]] = r
                except json.JSONDecodeError:
                    pass  # linea truncada por una interrupcion: se recalcula
    return done


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--smoke", action="store_true", help="tamanos pequenos -> results/rc3_smoke/")
    ap.add_argument("--fresh", action="store_true", help="ignora/borra pairs.jsonl previo")
    args = ap.parse_args(argv)
    out = OUT_SMOKE if args.smoke else OUT_FULL
    out.mkdir(parents=True, exist_ok=True)
    pj = out / "pairs.jsonl"
    if args.fresh and pj.exists():
        pj.unlink()
    done = _load(pj)
    _SPECS[:] = specs(args.smoke)
    todo = [i for i, s in enumerate(_SPECS) if s["id"] not in done]
    # los pares mas grandes primero
    todo.sort(key=lambda i: -max(int(_SPECS[i]["args"][1].get("N", 0)), 1))
    t0 = time.perf_counter()
    print(f"{len(_SPECS)} pares ({len(done)} ya hechos), procs={args.procs}", flush=True)
    with mp.get_context("fork").Pool(args.procs) as pool, open(pj, "a") as fh:
        for row in pool.imap_unordered(_run_idx, todo, chunksize=1):
            fh.write(json.dumps(row, default=str) + "\n")
            fh.flush()
            done[row["id"]] = row
            print(fmt(row), flush=True)
    rows = [done[s["id"]] for s in _SPECS]
    _atomic_write(pj, "".join(json.dumps(r, default=str) + "\n" for r in rows))
    ev = evaluate(rows)
    guards = {r["id"]: r["guards"] for r in rows if r["guards"]}
    summary = {"step": "rc3", "code_commit": head_commit(), "smoke": args.smoke, "n_pairs": len(rows),
               "seconds_this_run": round(time.perf_counter() - t0, 1),
               "pair_seconds_total": round(sum(r["seconds"] for r in rows), 1), "evaluation": ev,
               "verdict_new": ev["new"]["verdict"], "guards": guards, "by_pair": family_table(rows),
               "heisenberg_excluded_from_FE": True}
    _atomic_write(out / "summary.json", json.dumps(summary, indent=2, default=str))
    for part in ("new", "repeated", "all"):
        m = ev[part]
        print(f"[{part}] {m['verdict']} FE={m['FE']} EP={m['EP']} FE_by_d={ {k: v['excluded'] for k, v in m['FE_by_d'].items()} }")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
