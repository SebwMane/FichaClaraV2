"""P1-D.3 (docs/OMEGA_P1D3_PRERREGISTRO.md): calibracion del instrumento geometrico a N grande (BFS muestreado).

Salida: results/p1d3/{graphs.jsonl, summary.json, level3.json}. Uso:
  python tools/p1d3_panel.py [--procs 3] [--skip-level3] [--smoke]
Todos los grafos son dispersos; solo los analogos de Nivel III (N ~ 729) se densifican.

Decisiones de lectura literal (ver informe): radio del gradiente analitico (grado medio = N V r^3 * 13/12, con 13/12 =
int f^2 para f = (1+2x)/2, ignorando la discontinuidad de densidad en la costura periodica); la mitad superior de t del
D_s es la mitad superior de los puntos de la rejilla log; la sensibilidad exige |D_conv - d| <= 0.25.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from collections import defaultdict  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402
from scipy import sparse  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from c0_dynamics import _atomic_write  # noqa: E402
from l3b_stability import certificate_report  # noqa: E402

from omega.c0.references import rng_from_key  # noqa: E402
from omega.config.seeds import SeedKey  # noqa: E402
from omega.diagnostics.sampled_growth import (  # noqa: E402
    LEVEL2_THRESHOLDS,
    report_category,
    sampled_ball_profile,
    spectral_return,
    variety_category,
)
from omega.experiments.v11.gate import head_commit  # noqa: E402

MASTER = 20261009
NS = 20_000
NBIG = 50_000
SEEDS = (0, 1, 2)
OUT_FULL = ROOT / "results" / "p1d3"
OUT_SMOKE = ROOT / "results" / "p1d3_smoke"

FAMILY_ID: dict[str, int] = {
    "RGG2_k12": 1, "RGG3_k12": 2, "RGG3_k8": 3, "RGG_S2_k12": 4, "anillo_k12": 5, "cuadrado_141": 6, "T3_27": 7,
    "RGG3_caja_k12": 10, "RGG2_caja_k12": 11, "RGG3_gradiente": 12,
    "caveman_K4": 24, "caveman_K5": 25, "caveman_K6": 26, "caveman_K8": 28, "caveman_K10": 30, "caveman_K12": 32,
    "caveman_K16": 36, "cliques3D_K6_L15": 40, "cliques3D_K8_L13": 41, "cliques2D_K6_L58": 42,
    "retazos_b2": 50, "retazos_b4": 51, "retazos_b8": 52, "RGG3_atajos_0.1": 53, "RGG3_atajos_1": 54,
    "WS_b0.001": 60, "WS_b0.002": 61, "WS_b0.005": 62, "WS_b0.01": 63, "ER_k12": 70, "RR_k12": 71,
    "T3_37": 80, "cliques3D_K6_L7": 90,
}
TRUE_DIM = {"RGG2_k12": 2, "RGG3_k12": 3, "RGG3_k8": 3, "RGG_S2_k12": 2, "anillo_k12": 1, "cuadrado_141": 2, "T3_27": 3}
TRUE_DIM_SCALE = {"RGG3_k12": 3, "RGG2_k12": 2}

# ------------------------------------------------------------------ generadores dispersos


def _from_edges(n: int, u: np.ndarray, v: np.ndarray) -> sparse.csr_array:
    keep = u != v
    u, v = u[keep], v[keep]
    m = sparse.coo_array(
        (np.ones(2 * u.size), (np.concatenate([u, v]), np.concatenate([v, u]))), shape=(n, n)
    )
    a = sparse.csr_array(m)
    a.data[:] = 1.0
    return a


def _codes(n: int, u: np.ndarray, v: np.ndarray) -> np.ndarray:
    lo, hi = np.minimum(u, v).astype(np.int64), np.maximum(u, v).astype(np.int64)
    return lo * n + hi


def _unit_ball(d: int) -> float:
    return math.pi ** (d / 2.0) / math.gamma(d / 2.0 + 1.0)


def _radius(n: int, d: int, k: float) -> float:
    return float((k / ((n - 1) * _unit_ball(d))) ** (1.0 / d))


def _pairs(x: np.ndarray, r: float, periodic: bool) -> tuple[np.ndarray, np.ndarray]:
    tree = cKDTree(x, boxsize=1.0) if periodic else cKDTree(x)
    p = tree.query_pairs(r, output_type="ndarray")
    return p[:, 0].astype(np.int64), p[:, 1].astype(np.int64)


def rgg_points(n: int, d: int, k: float, rng: np.random.Generator, periodic: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = rng.random((n, d))
    u, v = _pairs(x, _radius(n, d, k), periodic)
    return x, u, v


def rgg(n: int, d: int, k: float, rng: np.random.Generator, periodic: bool = True) -> sparse.csr_array:
    _, u, v = rgg_points(n, d, k, rng, periodic)
    return _from_edges(n, u, v)


def rgg_sphere(n: int, k: float, rng: np.random.Generator) -> sparse.csr_array:
    x = rng.normal(size=(n, 3))
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    cos_t = 1.0 - 2.0 * k / n  # 2*pi*(1-cos t)*N/(4 pi) = k
    r = math.sqrt(2.0 * (1.0 - cos_t))  # cuerda = 2 sin(t/2)
    u, v = _pairs(x, r, False)
    return _from_edges(n, u, v)


def rgg_gradient(n: int, k: float, rng: np.random.Generator) -> sparse.csr_array:
    """Toro 3D; x con densidad (1+2x)/2 (CDF inversa); radio analitico: grado medio = N V r^3 * (13/12)."""
    x = rng.random((n, 3))
    x[:, 0] = (-1.0 + np.sqrt(1.0 + 8.0 * x[:, 0])) / 2.0
    x[:, 0] = np.minimum(x[:, 0], np.nextafter(1.0, 0.0))
    r = (k / ((n - 1) * _unit_ball(3) * (13.0 / 12.0))) ** (1.0 / 3.0)
    u, v = _pairs(x, float(r), True)
    return _from_edges(n, u, v)


def ring(n: int, k: int) -> sparse.csr_array:
    i = np.arange(n, dtype=np.int64)
    u = np.concatenate([i for _ in range(k // 2)])
    v = np.concatenate([(i + d) % n for d in range(1, k // 2 + 1)])
    return _from_edges(n, u, v)


def torus_lattice(side: int, dim: int) -> sparse.csr_array:
    n = side**dim
    idx = np.arange(n, dtype=np.int64)
    coords = np.stack(np.unravel_index(idx, (side,) * dim), axis=1)
    us, vs = [], []
    for ax in range(dim):
        c = coords.copy()
        c[:, ax] = (c[:, ax] + 1) % side
        us.append(idx)
        vs.append(np.ravel_multi_index(tuple(c.T), (side,) * dim).astype(np.int64))
    return _from_edges(n, np.concatenate(us), np.concatenate(vs))


def caveman(n: int, size: int) -> sparse.csr_array:
    """Misma regla que connected_caveman: cliques consecutivas (sobrantes a la ultima); quitar (a,a+1), anadir (a,a-1 mod n)."""
    nb = n // size
    bounds = [(b * size, (b + 1) * size) for b in range(nb)]
    bounds[-1] = (bounds[-1][0], n)
    us, vs = [], []
    for a, b in bounds:
        iu, ju = np.triu_indices(b - a, 1)
        us.append(iu + a)
        vs.append(ju + a)
    u, v = np.concatenate(us).astype(np.int64), np.concatenate(vs).astype(np.int64)
    starts = np.array([a for a, _ in bounds], dtype=np.int64)
    drop = np.isin(_codes(n, u, v), _codes(n, starts, starts + 1))
    u, v = u[~drop], v[~drop]
    return _from_edges(n, np.concatenate([u, starts]), np.concatenate([v, (starts - 1) % n]))


def clique_lattice(side: int, m: int, dim: int) -> sparse.csr_array:
    """Toro side^dim de sitios, cada uno K_m; direccion j enlaza miembro j%m con miembro (j^1)%m del vecino en esa direccion."""
    ns = side**dim
    n = ns * m
    sites = np.arange(ns, dtype=np.int64)
    coords = np.stack(np.unravel_index(sites, (side,) * dim), axis=1)
    us, vs = [], []
    ia, ib = np.triu_indices(m, 1)
    us.append((sites[:, None] * m + ia[None, :]).ravel())
    vs.append((sites[:, None] * m + ib[None, :]).ravel())
    j = 0
    for ax in range(dim):
        for sgn in (1, -1):
            c = coords.copy()
            c[:, ax] = (c[:, ax] + sgn) % side
            t = np.ravel_multi_index(tuple(c.T), (side,) * dim).astype(np.int64)
            us.append(sites * m + j % m)
            vs.append(t * m + (j ^ 1) % m)
            j += 1
    return _from_edges(n, np.concatenate(us), np.concatenate(vs))


def _rewire_one_end(n: int, u: np.ndarray, v: np.ndarray, pick: np.ndarray, rng: np.random.Generator) -> sparse.csr_array:
    """Para cada arista elegida (u,v): w uniforme; si w != u y (u,w) no existe, cambia (u,v) por (u,w)."""
    w = rng.integers(n, size=u.size)
    exists = np.isin(_codes(n, u, w), _codes(n, u, v))
    ok = pick & (w != u) & ~exists
    return _from_edges(n, np.concatenate([u[~ok], u[ok]]), np.concatenate([v[~ok], w[ok]]))


def shortcuts(n: int, frac: float, rng: np.random.Generator) -> sparse.csr_array:
    _, u, v = rgg_points(n, 3, 12.0, rng)
    return _rewire_one_end(n, u, v, rng.random(u.size) < frac, rng)


def watts_strogatz(n: int, k: int, beta: float, rng: np.random.Generator) -> sparse.csr_array:
    i = np.arange(n, dtype=np.int64)
    u = np.concatenate([i for _ in range(k // 2)])
    v = np.concatenate([(i + d) % n for d in range(1, k // 2 + 1)])
    return _rewire_one_end(n, u, v, rng.random(u.size) < beta, rng)


def patchwork(n: int, blocks: int, rng: np.random.Generator) -> sparse.csr_array:
    """RGG3 k12 toroidal; extremos de aristas entre bloques distintos reemparejados al azar (grados preservados)."""
    x, u, v = rgg_points(n, 3, 12.0, rng)
    bid = np.floor(x * blocks).astype(np.int64) @ np.array([blocks * blocks, blocks, 1], dtype=np.int64)
    cross = bid[u] != bid[v]
    stubs = np.concatenate([u[cross], v[cross]])
    rng.shuffle(stubs)
    half = stubs.size // 2
    return _from_edges(n, np.concatenate([u[~cross], stubs[0 : 2 * half : 2]]), np.concatenate([v[~cross], stubs[1 : 2 * half : 2]]))


def erdos_renyi(n: int, m: int, rng: np.random.Generator) -> sparse.csr_array:
    have = np.empty(0, dtype=np.int64)
    while have.size < m:
        a = rng.integers(n, size=2 * (m - have.size) + 16)
        b = rng.integers(n, size=a.size)
        keep = a != b
        have = np.unique(np.concatenate([have, _codes(n, a[keep], b[keep])]))
    sel = rng.permutation(have)[:m] if have.size > m else have
    return _from_edges(n, sel // n, sel % n)


def random_regular(n: int, k: int, rng: np.random.Generator) -> sparse.csr_array:
    import networkx as nx

    g = nx.random_regular_graph(k, n, seed=int(rng.integers(2**31)))
    return sparse.csr_array(nx.to_scipy_sparse_array(g, format="csr", dtype=np.float64))


# ------------------------------------------------------------------ panel


def build(spec: dict[str, Any], rng: np.random.Generator) -> sparse.csr_array:
    g, a, n = spec["gen"], spec.get("args", {}), int(spec["N"])
    if g == "rgg":
        return rgg(n, a["d"], a["k"], rng, a.get("periodic", True))
    if g == "sphere":
        return rgg_sphere(n, 12.0, rng)
    if g == "gradient":
        return rgg_gradient(n, 12.0, rng)
    if g == "ring":
        return ring(n, 12)
    if g == "lattice":
        return torus_lattice(a["side"], a["dim"])
    if g == "caveman":
        return caveman(n, a["size"])
    if g == "cliques":
        return clique_lattice(a["side"], a["m"], a["dim"])
    if g == "patch":
        return patchwork(n, a["blocks"], rng)
    if g == "shortcuts":
        return shortcuts(n, a["frac"], rng)
    if g == "ws":
        return watts_strogatz(n, 12, a["beta"], rng)
    if g == "er":
        return erdos_renyi(n, 6 * n, rng)
    if g == "rr":
        return random_regular(n, 12, rng)
    raise ValueError(g)


def specs(smoke: bool = False) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []

    def add(fam: str, cls: str, gen: str, n: int, seeds: tuple[int, ...] = SEEDS, panel: str = "main", suffix: str = "",
            fid_off: int = 0, **args: Any) -> None:
        for s in seeds:
            out.append({"id": f"{fam}{suffix}_s{s}", "family": fam + suffix, "cls": cls, "panel": panel, "gen": gen, "N": n,
                        "args": args, "key": (MASTER, FAMILY_ID[fam] + fid_off, s),
                        "true_dim": (TRUE_DIM.get(fam) if panel == "main" else TRUE_DIM_SCALE.get(fam))})

    if smoke:
        sd = (0,)
        n = 2000
        add("RGG2_k12", "G-hom", "rgg", n, sd, d=2, k=12.0)
        add("RGG3_k12", "G-hom", "rgg", n, sd, d=3, k=12.0)
        add("RGG_S2_k12", "G-hom", "sphere", n, sd)
        add("anillo_k12", "G-hom", "ring", n, sd)
        add("RGG3_caja_k12", "G-inh", "rgg", n, sd, d=3, k=12.0, periodic=False)
        add("RGG3_gradiente", "G-inh", "gradient", n, sd)
        add("caveman_K8", "R", "caveman", n, sd, size=8)
        add("cliques3D_K6_L7", "R", "cliques", 7**3 * 6, sd, side=7, m=6, dim=3)
        add("retazos_b2", "L", "patch", n, sd, blocks=2)
        add("RGG3_atajos_1", "L", "shortcuts", n, sd, frac=0.01)
        add("WS_b0.005", "L", "ws", n, sd, beta=0.005)
        add("ER_k12", "NL", "er", n, sd)
        add("RR_k12", "NL", "rr", n, sd)
        return out
    add("RGG2_k12", "G-hom", "rgg", NS, d=2, k=12.0)
    add("RGG3_k12", "G-hom", "rgg", NS, d=3, k=12.0)
    add("RGG3_k8", "G-hom", "rgg", NS, d=3, k=8.0)
    add("RGG_S2_k12", "G-hom", "sphere", NS)
    add("anillo_k12", "G-hom", "ring", NS, (0,))
    add("cuadrado_141", "G-hom", "lattice", 141**2, (0,), side=141, dim=2)
    add("T3_27", "G-hom", "lattice", 27**3, (0,), side=27, dim=3)
    add("RGG3_caja_k12", "G-inh", "rgg", NS, d=3, k=12.0, periodic=False)
    add("RGG2_caja_k12", "G-inh", "rgg", NS, d=2, k=12.0, periodic=False)
    add("RGG3_gradiente", "G-inh", "gradient", NS)
    for size in (4, 5, 6, 8, 10, 12, 16):
        add(f"caveman_K{size}", "R", "caveman", NS, (0,), size=size)
    add("cliques3D_K6_L15", "R", "cliques", 15**3 * 6, (0,), side=15, m=6, dim=3)
    add("cliques3D_K8_L13", "R", "cliques", 13**3 * 8, (0,), side=13, m=8, dim=3)
    add("cliques2D_K6_L58", "R", "cliques", 58**2 * 6, (0,), side=58, m=6, dim=2)
    for b in (2, 4, 8):
        add(f"retazos_b{b}", "L", "patch", NS, blocks=b)
    add("RGG3_atajos_0.1", "L", "shortcuts", NS, frac=0.001)
    add("RGG3_atajos_1", "L", "shortcuts", NS, frac=0.01)
    for beta in (0.001, 0.002, 0.005, 0.01):
        add(f"WS_b{beta}", "L", "ws", NS, beta=beta)
    add("ER_k12", "NL", "er", NS)
    add("RR_k12", "NL", "rr", NS)
    # subpanel de escala N ~ 5e4 (semillas 0 y 1; el reticulo T3_37 es determinista: solo semilla 0)
    sc = dict(panel="escala", suffix="_N5e4", fid_off=100)
    s2 = (0, 1)
    add("RGG3_k12", "G-hom", "rgg", NBIG, s2, d=3, k=12.0, **sc)
    add("RGG2_k12", "G-hom", "rgg", NBIG, s2, d=2, k=12.0, **sc)
    add("retazos_b4", "L", "patch", NBIG, s2, blocks=4, **sc)
    add("WS_b0.002", "L", "ws", NBIG, s2, beta=0.002, **sc)
    add("caveman_K8", "R", "caveman", NBIG, (0,), size=8, **sc)
    add("T3_37", "G-hom", "lattice", 37**3, (0,), side=37, dim=3, **sc)
    return out


_SPECS: list[dict[str, Any]] = []


def run_spec(spec: dict[str, Any]) -> dict[str, Any]:
    t0 = time.perf_counter()
    key = tuple(spec["key"])
    adj = build(spec, rng_from_key(key))
    t1 = time.perf_counter()
    prof = sampled_ball_profile(adj, rng_from_key(key + (1,)))
    t2 = time.perf_counter()
    t_max = int(min(4000, 4 * max(prof["max_window"], 2) ** 2))
    sp = spectral_return(adj, rng_from_key(key + (2,)), 64, t_max)
    t3 = time.perf_counter()
    row = {k: v for k, v in spec.items() if k not in ("gen", "args", "key")}
    row["key"] = list(key)
    row["n"] = int(adj.shape[0])
    row["mean_degree"] = float(adj.nnz / adj.shape[0])
    row.update(prof)
    cat: dict[str, str] = {}
    for t in LEVEL2_THRESHOLDS:
        l2 = prof["level2"][str(t)]
        cat[str(t)] = report_category(bool(prof["level1"]["pass"]), l2["status"], l2["class"])
    row["category_by_threshold"] = cat
    row["category"] = cat["0.05"]
    row["variety_category"] = None
    row["spectral"] = sp
    row["seconds"] = {"gen": round(t1 - t0, 2), "profile": round(t2 - t1, 2), "spectral": round(t3 - t2, 2)}
    return row


def _run_idx(i: int) -> dict[str, Any]:
    row = run_spec(_SPECS[i])
    print(f"  {row['id']:<28} n={row['n']:<6} deg={row['mean_degree']:.1f} I={row['level1']['status']:<15} "
          f"II={row['level2']['0.05']['status']:<13} D={row['level2']['0.05']['D_conv']} {row['category']} {row['seconds']}", flush=True)
    return row


# ------------------------------------------------------------------ veredicto


def _joint(r: dict[str, Any], t: float) -> bool:
    return bool(r["level1"]["pass"]) and r["level2"][str(t)]["status"] == "CONVERGE"


def _correct(r: dict[str, Any], t: float) -> bool:
    d = r["level2"][str(t)]["D_conv"]
    return _joint(r, t) and d is not None and r["true_dim"] is not None and abs(d - r["true_dim"]) <= 0.25


def family_table(rows: list[dict[str, Any]], t: float) -> dict[str, Any]:
    fam: dict[str, dict[str, Any]] = defaultdict(lambda: {"cls": None, "panel": None, "n": 0, "I": 0, "II": 0, "joint": 0,
                                                          "I_status": [], "II_status": [], "D_conv": [], "rho": [], "category": []})
    for r in rows:
        f = fam[r["family"]]
        f["cls"], f["panel"] = r["cls"], r["panel"]
        f["n"] += 1
        f["I"] += int(r["level1"]["pass"])
        f["II"] += int(r["level2"][str(t)]["status"] == "CONVERGE")
        f["joint"] += int(_joint(r, t))
        f["I_status"].append(r["level1"]["status"])
        f["II_status"].append(r["level2"][str(t)]["status"])
        d = r["level2"][str(t)]["D_conv"]
        f["D_conv"].append(None if d is None else round(d, 3))
        f["rho"].append(None if r["level1"]["rho"] is None else round(r["level1"]["rho"], 3))
        f["category"].append(r["category_by_threshold"][str(t)])
        f.setdefault("D_s", []).append(None if r["spectral"]["D_s_median_upper"] is None else round(r["spectral"]["D_s_median_upper"], 3))
    return dict(fam)


def verdict(rows: list[dict[str, Any]], t: float, level3_cond: Any) -> dict[str, Any]:
    main = [r for r in rows if r["panel"] == "main"]
    ghom = [r for r in main if r["cls"] == "G-hom"]
    nongeo = [r for r in main if r["cls"] in ("L", "NL")]
    sens = sum(_correct(r, t) for r in ghom) / len(ghom)
    spec = sum(not _joint(r, t) for r in nongeo) / len(nongeo)
    if sens < 0.7 or spec < 0.7:
        v = "INSTRUMENTO-INVALIDO"
    elif sens >= 0.9 and spec >= 0.9 and level3_cond is not False and level3_cond != "NO_EVALUADO":
        v = "INSTRUMENTO-VALIDO"
    else:
        v = "INSTRUMENTO-PARCIAL"
    return {"threshold": t, "verdict": v, "sensitivity": sens, "specificity": spec, "n_ghom": len(ghom), "n_nongeo": len(nongeo),
            "level3_condition": level3_cond,
            "level3_note": "NO_EVALUADO: verdict provisional (sin Nivel III)" if level3_cond == "NO_EVALUADO" else None,
            "by_family": family_table(rows, t)}


# ------------------------------------------------------------------ Nivel III (§5)

ANALOGS: dict[str, tuple[str, dict[str, Any]]] = {
    "cliques3D_K6_L15": ("cliques3D_K6_L5", {"gen": "cliques", "side": 5, "m": 6, "dim": 3}),
    "cliques3D_K8_L13": ("cliques3D_K8_L4", {"gen": "cliques", "side": 4, "m": 8, "dim": 3}),
    "RGG3_caja_k12": ("RGG3_caja_n729", {"gen": "box"}),
    "RGG3_gradiente": ("RGG3_gradiente_n729", {"gen": "gradient"}),
}


def _analog(name: str, spec: dict[str, Any], rng: np.random.Generator) -> np.ndarray:
    g = spec["gen"]
    if g == "cliques":
        a = clique_lattice(spec["side"], spec["m"], spec["dim"])
    elif g == "box":
        a = rgg(729, 3, 12.0, rng, periodic=False)
    elif g == "gradient":
        a = rgg_gradient(729, 12.0, rng)
    elif g == "ref":
        a = rgg(729, 3, 12.0, rng, periodic=True)
    else:
        raise ValueError(name)
    return np.asarray(a.toarray(), dtype=np.float64)


def _cert_rejects(codes: Any) -> bool:
    return isinstance(codes, list) and any(str(c).endswith(("F9", "F4")) for c in codes)


def level3(rows: list[dict[str, Any]], smoke: bool = False) -> dict[str, Any]:
    main = [r for r in rows if r["panel"] == "main"]
    qualifying: list[str] = []
    by_fam: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for r in main:
        by_fam[r["family"]].append(r)
    for fam, rs in by_fam.items():
        if rs[0]["cls"] not in ("R", "G-inh"):
            continue
        ok = sum(1 for r in rs if r["category"].startswith("GEOMETRIA_GRUESA(") and r["category"][17:-1].isdigit()
                 and int(r["category"][17:-1]) >= 3)
        if ok >= 2.0 / 3.0 * len(rs):
            qualifying.append(fam)
    res: dict[str, Any] = {"qualifying_families": sorted(qualifying), "unmapped_families": [], "analogs": {}, "references": {}}
    jobs: list[tuple[str, str, str, dict[str, Any], int]] = []
    for i, fam in enumerate(sorted(qualifying)):
        if fam in ANALOGS:
            name, sp = ANALOGS[fam]
            jobs.append(("analogs", name, fam, sp, 200 + i))
        else:
            res["unmapped_families"].append(fam)
    jobs.append(("references", "RGG3_k12_n729", "RGG3_k12", {"gen": "ref"}, 299))
    for bucket, name, fam, sp, fid in jobs:
        adj = _analog(name, sp, rng_from_key((MASTER, fid, 0)))
        reps = [certificate_report(adj, "CONVERGED", SeedKey(MASTER, (fid, s))) for s in range(3)]
        res[bucket][name] = {"family": fam, "n": int(adj.shape[0]), "codes": [r.get("codes") for r in reps],
                             "errors": [r.get("error") for r in reps],
                             "passes_all": all(r.get("codes") == [] and "error" not in r for r in reps),
                             "rejects_F9_or_F4_all": all(_cert_rejects(r.get("codes")) for r in reps)}
        print(f"  nivel III {name}: {res[bucket][name]['codes']}", flush=True)
    r3 = [v for v in res["analogs"].values() if str(v["family"]).startswith("cliques3D")]
    res["R_class3_analogs"] = sorted(k for k, v in res["analogs"].items() if str(v["family"]).startswith("cliques3D"))
    res["R_class3_all_rejected"] = all(v["rejects_F9_or_F4_all"] for v in r3) if r3 else None
    return res


def apply_variety(rows: list[dict[str, Any]], l3: dict[str, Any]) -> None:
    cert: dict[str, bool] = {}
    for bucket in ("analogs", "references"):
        for v in l3[bucket].values():
            cert[v["family"]] = bool(v["passes_all"])
    for r in rows:
        cls = r["level2"]["0.05"]["class"]
        r["variety_category"] = variety_category(r["category"], cls, cert.get(r["family"]) if r["panel"] == "main" else None)


# ------------------------------------------------------------------ main


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--skip-level3", action="store_true")
    ap.add_argument("--smoke", action="store_true", help="N~2000, 1 semilla, pocas familias -> results/p1d3_smoke/")
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
    l3: dict[str, Any] | None = None
    cond: Any = "NO_EVALUADO"
    if not args.skip_level3:
        l3 = level3(rows, args.smoke)
        apply_variety(rows, l3)
        cond = True if l3["R_class3_all_rejected"] is None else bool(l3["R_class3_all_rejected"])
        _atomic_write(out / "level3.json", json.dumps(l3, indent=2, default=str))
    _atomic_write(out / "graphs.jsonl", "".join(json.dumps(r, default=str) + "\n" for r in rows))
    res = {str(t): verdict(rows, t, cond) for t in LEVEL2_THRESHOLDS}
    scope = {str(t): {"G-inh": {f: v for f, v in family_table(rows, t).items() if v["cls"] == "G-inh"},
                      "R": {f: v for f, v in family_table(rows, t).items() if v["cls"] == "R"}} for t in LEVEL2_THRESHOLDS}
    summary: dict[str, Any] = {"step": "p1d3_panel", "code_commit": head_commit(), "smoke": args.smoke, "n_graphs": len(rows),
                               "seconds": round(time.perf_counter() - t0, 1), "level3_condition": cond,
                               "final_threshold": 0.05, "result": res, "scope_tables": scope}
    _atomic_write(out / "summary.json", json.dumps(summary, indent=2, default=str))
    for t, v in res.items():
        print(f"umbral={t}: {v['verdict']} sens={v['sensitivity']:.2f} spec={v['specificity']:.2f} L3={v['level3_condition']}")
    for f, d in res["0.05"]["by_family"].items():
        print(f"  {f:<22} {d['cls']:<6} I={d['I']}/{d['n']} II={d['II']}/{d['n']} D={d['D_conv']} rho={d['rho']} Ds={d['D_s']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
