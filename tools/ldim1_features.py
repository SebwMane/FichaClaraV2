"""Omega L-DIM-1 (docs/OMEGA_LDIM1_PRERREGISTRO.md §3, §6, A1): generacion del panel + extraccion CIEGA de caracteristicas.

Salida: results/ldim1/{features.jsonl, truth.jsonl, meta.json}. Uso:
  python tools/ldim1_features.py [--procs 3] [--smoke]
`extract(adj, id, rng)` solo recibe la adyacencia, el ID opaco y un rng: no ve familia, d, grado nominal ni estrato (NC-3).
La verdad (familia, grupo, d, estrato, escala, N, semilla, grado medio calculado por el generador) va a truth.jsonl, que solo
lee tools/ldim1_eval.py. Smoke: tamanos reducidos, semilla 0, resultados en results/ldim1_smoke/.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
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

from c0_dynamics import _atomic_write  # noqa: E402
from cic_l0b import random_two_tree  # noqa: E402
from coh_l0b import husimi_cactus, prufer_tree, subdivided_binary_tree  # noqa: E402
from p1d3_panel import (  # noqa: E402
    _from_edges,
    caveman,
    erdos_renyi,
    patchwork,
    random_regular,
    rgg,
    ring,
    torus_lattice,
    watts_strogatz,
)

from omega.c0.references import rng_from_key  # noqa: E402
from omega.diagnostics.sampled_growth import sampled_ball_profile  # noqa: E402
from omega.experiments.v11.gate import head_commit  # noqa: E402

MASTER = 20261014
SEEDS = (0, 1, 2)
N_SOURCES = 400
OUT_FULL = ROOT / "results" / "ldim1"
OUT_SMOKE = ROOT / "results" / "ldim1_smoke"

FAMILY_ID: dict[str, int] = {
    "anillo_k2": 1, "anillo_k8": 2, "cuadrado": 3, "RGG2_k8": 4, "RGG2_k16": 5, "T3": 6, "RGG3_k8": 7, "RGG3_k16": 8,
    "T4": 9, "RGG4_k8": 10, "RGG4_k16": 11, "T5": 12,
    "anillo_k2_4N": 21, "cuadrado_4N": 22, "T3_4N": 23, "T4_4N": 24, "T5_4N": 25, "RGG3_k16_4N": 26,
    "clique_K200": 40, "arbol_prufer": 41, "arbol_binario_sub8": 42, "cactus_triangulos": 43, "2arbol": 44, "RR_k12": 45,
    "ER_k12": 46, "caveman_K8": 47, "WS_b0.01": 48, "retazos_b2": 49,
}

# fam, grupo, d, estrato, escala, gen, (N, args) completo, (N, args) smoke, deterministica
_G, _X = "G", "X"
_Row = tuple[str, str, "int | None", "str | None", str, str, tuple[int, dict[str, Any]], tuple[int, dict[str, Any]], bool]


def _lat(side: int, dim: int) -> tuple[int, dict[str, Any]]:
    return side**dim, {"side": side, "dim": dim}


def _rgg(n: int, d: int, k: float) -> tuple[int, dict[str, Any]]:
    return n, {"d": d, "k": k}


_PANEL: list[_Row] = [
    # geometrias, escala N
    ("anillo_k2", _G, 1, "reticulo", "N", "ring", (20000, {"k": 2}), (2000, {"k": 2}), True),
    ("anillo_k8", _G, 1, "k8", "N", "ring", (20000, {"k": 8}), (2000, {"k": 8}), True),
    ("cuadrado", _G, 2, "reticulo", "N", "lattice", _lat(141, 2), _lat(45, 2), True),
    ("RGG2_k8", _G, 2, "k8", "N", "rgg", _rgg(20000, 2, 8.0), _rgg(2000, 2, 8.0), False),
    ("RGG2_k16", _G, 2, "k16", "N", "rgg", _rgg(20000, 2, 16.0), _rgg(2000, 2, 16.0), False),
    ("T3", _G, 3, "reticulo", "N", "lattice", _lat(27, 3), _lat(13, 3), True),
    ("RGG3_k8", _G, 3, "k8", "N", "rgg", _rgg(20000, 3, 8.0), _rgg(2000, 3, 8.0), False),
    ("RGG3_k16", _G, 3, "k16", "N", "rgg", _rgg(20000, 3, 16.0), _rgg(2000, 3, 16.0), False),
    ("T4", _G, 4, "reticulo", "N", "lattice", _lat(12, 4), _lat(7, 4), True),
    ("RGG4_k8", _G, 4, "k8", "N", "rgg", _rgg(20000, 4, 8.0), _rgg(2000, 4, 8.0), False),
    ("RGG4_k16", _G, 4, "k16", "N", "rgg", _rgg(20000, 4, 16.0), _rgg(2000, 4, 16.0), False),
    ("T5", _G, 5, "reticulo", "N", "lattice", _lat(7, 5), _lat(5, 5), True),
    # geometrias, escala ~4N (L7)
    ("anillo_k2_4N", _G, 1, "reticulo", "4N", "ring", (80000, {"k": 2}), (8000, {"k": 2}), True),
    ("cuadrado_4N", _G, 2, "reticulo", "4N", "lattice", _lat(283, 2), _lat(90, 2), True),
    ("T3_4N", _G, 3, "reticulo", "4N", "lattice", _lat(43, 3), _lat(21, 3), True),
    ("T4_4N", _G, 4, "reticulo", "4N", "lattice", _lat(17, 4), _lat(10, 4), True),
    ("T5_4N", _G, 5, "reticulo", "4N", "lattice", _lat(9, 5), _lat(7, 5), True),
    ("RGG3_k16_4N", _G, 3, "k16", "4N", "rgg", _rgg(80000, 3, 16.0), _rgg(8000, 3, 16.0), False),
    # degenerados
    ("clique_K200", _X, None, None, "N", "clique", (200, {}), (200, {}), True),
    ("arbol_prufer", _X, None, None, "N", "prufer", (20000, {}), (2000, {}), False),
    ("arbol_binario_sub8", _X, None, None, "N", "subtree", (0, {"levels": 11, "seg": 8}), (0, {"levels": 7, "seg": 8}), True),
    ("cactus_triangulos", _X, None, None, "N", "cactus", (20000, {}), (2000, {}), True),
    ("2arbol", _X, None, None, "N", "2tree", (20000, {}), (2000, {}), False),
    ("RR_k12", _X, None, None, "N", "rr", (20000, {}), (2000, {}), False),
    ("ER_k12", _X, None, None, "N", "er", (20000, {}), (2000, {}), False),
    ("caveman_K8", _X, None, None, "N", "caveman", (20000, {"size": 8}), (2000, {"size": 8}), True),
    ("WS_b0.01", _X, None, None, "N", "ws", (20000, {"beta": 0.01}), (2000, {"beta": 0.01}), False),
    ("retazos_b2", _X, None, None, "N", "patch", (20000, {"blocks": 2}), (2000, {"blocks": 2}), False),
]


def clique(n: int) -> sparse.csr_array:
    iu, iv = np.triu_indices(n, 1)
    return _from_edges(n, iu.astype(np.int64), iv.astype(np.int64))


def build(gen: str, n: int, a: dict[str, Any], rng: np.random.Generator) -> sparse.csr_array:
    if gen == "ring":
        return ring(n, a["k"])
    if gen == "lattice":
        return torus_lattice(a["side"], a["dim"])
    if gen == "rgg":
        return rgg(n, a["d"], a["k"], rng, True)
    if gen == "clique":
        return clique(200)
    if gen == "prufer":
        return prufer_tree(n, rng)
    if gen == "subtree":
        return subdivided_binary_tree(a["levels"], a["seg"])
    if gen == "cactus":
        return husimi_cactus(n)
    if gen == "2tree":
        return random_two_tree(n, rng)
    if gen == "rr":
        return random_regular(n, 12, rng)
    if gen == "er":
        return erdos_renyi(n, 6 * n, rng)
    if gen == "caveman":
        return caveman(n, a["size"])
    if gen == "ws":
        return watts_strogatz(n, 12, a["beta"], rng)
    if gen == "patch":
        return patchwork(n, a["blocks"], rng)
    raise ValueError(gen)


def graph_id(key: tuple[int, ...]) -> str:
    return hashlib.sha256(repr(tuple(key)).encode()).hexdigest()[:12]


def specs(smoke: bool = False) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for fam, grp, d, stratum, scale, gen, full, small, det in _PANEL:
        n, args = small if smoke else full
        seeds = (0,) if (det or smoke) else SEEDS
        for s in seeds:
            key = (MASTER, FAMILY_ID[fam], s)
            out.append({"id": graph_id(key), "family": fam, "group": grp, "d": d, "stratum": stratum, "scale": scale,
                        "N_nominal": n, "seed": s, "gen": gen, "args": args, "key": key})
    return out


# ------------------------------------------------------------------ extraccion ciega


def extract(adj: sparse.csr_array, gid: str, rng: np.random.Generator) -> dict[str, Any]:
    """Caracteristicas ciegas (§1, A1.1): solo adyacencia + ID + rng. r* = max(1, max_window // 2)."""
    prof = sampled_ball_profile(adj, rng, n_sources=N_SOURCES)
    mean = [float(x) for x in prof["mean"]]  # mean[r-1] = m(r), r = 1..r_cap
    w = int(prof["max_window"])
    r_star = max(1, w // 2)
    m = [1.0] + mean  # m[r] = m(r), m(0) = 1
    s_star = m[r_star] - m[r_star - 1]
    phi = (m[r_star + 1] - m[r_star]) / m[r_star] if r_star + 1 <= len(mean) else 0.0
    psi = 1.0 / s_star
    return {"id": gid, "r_star": r_star, "max_window": w, "n_giant": int(prof["n_giant"]), "m": mean, "phi": float(phi),
            "psi": float(psi)}


_SPECS: list[dict[str, Any]] = []


def run_spec(spec: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], dict[str, float]]:
    t0 = time.perf_counter()
    key = tuple(int(k) for k in spec["key"])
    adj = build(spec["gen"], int(spec["N_nominal"]), spec["args"], rng_from_key(key))
    n = int(adj.shape[0])
    mean_degree = float(adj.nnz / n)  # = 2 * aristas / n (lado del generador)
    t1 = time.perf_counter()
    feat = extract(adj, spec["id"], rng_from_key(key + (1,)))  # el extractor no recibe spec
    t2 = time.perf_counter()
    truth = {"id": spec["id"], "family": spec["family"], "group": spec["group"], "d": spec["d"], "stratum": spec["stratum"],
             "scale": spec["scale"], "N": n, "seed": spec["seed"], "mean_degree": mean_degree}
    return feat, truth, {"gen": round(t1 - t0, 2), "extract": round(t2 - t1, 2)}


def _run_idx(i: int) -> tuple[dict[str, Any], dict[str, Any], dict[str, float]]:
    feat, truth, sec = run_spec(_SPECS[i])
    print(f"  {truth['family']:<20} s{truth['seed']} n={truth['N']:<6} deg={truth['mean_degree']:.1f} r*={feat['r_star']:<3} "
          f"phi={feat['phi']:.4f} psi={feat['psi']:.4f} {sec}", flush=True)
    return feat, truth, sec


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--smoke", action="store_true", help="tamanos reducidos, semilla 0 -> results/ldim1_smoke/")
    args = ap.parse_args(argv)
    out = OUT_SMOKE if args.smoke else OUT_FULL
    _SPECS[:] = specs(args.smoke)
    ids = [s["id"] for s in _SPECS]
    assert len(set(ids)) == len(ids), "colision de IDs"
    order = sorted(range(len(_SPECS)), key=lambda i: -int(_SPECS[i]["N_nominal"]))
    t0 = time.perf_counter()
    print(f"{len(_SPECS)} grafos, procs={args.procs}", flush=True)
    done: dict[int, tuple[dict[str, Any], dict[str, Any], dict[str, float]]] = {}
    with mp.get_context("fork").Pool(args.procs) as pool:
        for i, res in zip(order, pool.imap(_run_idx, order, chunksize=1)):
            done[i] = res
    rows = [done[i] for i in range(len(_SPECS))]
    out.mkdir(parents=True, exist_ok=True)
    _atomic_write(out / "features.jsonl", "".join(json.dumps(f) + "\n" for f, _, _ in rows))
    _atomic_write(out / "truth.jsonl", "".join(json.dumps(t) + "\n" for _, t, _ in rows))
    meta = {"step": "ldim1_features", "code_commit": head_commit(), "smoke": args.smoke, "master": MASTER, "n_graphs": len(rows),
            "n_sources": N_SOURCES, "seconds": round(time.perf_counter() - t0, 1),
            "seconds_by_graph": {t["id"]: s for _, t, s in rows}}
    _atomic_write(out / "meta.json", json.dumps(meta, indent=2))
    print(f"listo: {len(rows)} grafos en {meta['seconds']} s -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
