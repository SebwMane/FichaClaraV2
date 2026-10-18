"""Omega-C0: bateria geometrica DESCRIPTIVA sobre las celdas CANDIDATO-C0 (prerregistro §4, enmienda C0-A6).

Solo describe; no hay umbral nuevo ni veredicto de geometria (los umbrales del certificado no se tocan).
Por cada final DISPERSO-LOCAL de celdas candidatas (N = 216 de C0-L4; N = 125 y 343 de R6) sobre el soporte
fuerte A = {W > 0.1 max W}, componente gigante:
  D_eff (saltos, estimador del bloque 1.1), D_s (paseo lazy sobre A binario), homogeneidad, isotropia (k = round D_eff),
  salto medio, diametro, descriptores de cliques (tamano de la clique maxima, cliques maximales >= 4 por nodo).
Escalado con el tamano (C0-A6): L(N) ~ N^(1/D_L) con el salto medio en N = 125, 216, 343 por (celda, inicio).
Controles (R5): T3 6^3, RGG3 k12, toro triangular 12x18, anillo k6, caveman conectado K8, arbol aleatorio,
y sus versiones a N = 125 / 343 cuando existen (toros 5^3, 7^3) para calibrar D_L.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
from collections import Counter, defaultdict  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from c0_dynamics import _atomic_write, read_rows  # noqa: E402

from omega.c0 import references as R  # noqa: E402
from omega.c0.locality import giant_component, strong_support  # noqa: E402
from omega.config.settings import DimensionConfig, SpectralConfig  # noqa: E402
from omega.config.settings11 import LocalStructureConfig  # noqa: E402
from omega.experiments.v11.gate import head_commit  # noqa: E402
from omega.geometry.dimension import effective_dimension  # noqa: E402
from omega.geometry.distances import hop_distance_matrix  # noqa: E402
from omega.geometry.local_structure import homogeneity, isotropy  # noqa: E402
from omega.geometry.spectral import spectral_dimension  # noqa: E402
from omega.types import FloatArray  # noqa: E402

OUT = ROOT / "results" / "c0_battery"
MASTER = 20261005


def describe(w: FloatArray, seed: int) -> dict[str, Any]:
    a = strong_support(w)
    g = giant_component(a)
    ag = np.asarray(a[np.ix_(g, g)], dtype=np.bool_)
    d = hop_distance_matrix(ag)
    fin = d[np.isfinite(d) & (d > 0)]
    est = effective_dimension(d, DimensionConfig())
    ds = spectral_dimension(ag.astype(np.float64), SpectralConfig())
    lcfg = LocalStructureConfig()
    hom = homogeneity(d, est, lcfg)
    k = int(round(est.value)) if np.isfinite(est.value) and est.value >= 1 else None
    iso = isotropy(d, est, k, lcfg, np.random.Generator(np.random.PCG64(seed)))
    gx = nx.from_numpy_array(ag.astype(int))
    cl = [c for c in nx.find_cliques(gx) if len(c) >= 4]
    memb = Counter(v for c in cl for v in c)
    deg = ag.sum(axis=1)
    return {
        "n_giant": int(g.size), "mean_deg": float(deg.mean()), "mean_hop": float(fin.mean()), "diameter": float(fin.max()),
        "D_eff": est.value, "D_eff_status": est.status, "D_s": ds.value, "D_s_status": ds.status,
        "hom_cv": hom.cv, "iso_median": iso.median_ratio, "iso_k": k,
        "max_clique": max((len(c) for c in nx.find_cliques(gx)), default=0), "n_cliques_ge4": len(cl),
        "cliques_per_node": float(np.mean([memb.get(i, 0) for i in range(g.size)])),
        "frac_nodes_in_clique_ge4": float(np.mean([memb.get(i, 0) > 0 for i in range(g.size)])),
    }


def _task(args: tuple[str, str, int, dict[str, Any]]) -> dict[str, Any]:
    tag, path, seed, meta = args
    z = np.load(path)
    w = np.asarray(z[z.files[0]], dtype=np.float64)
    return {"tag": tag, **meta, **describe(w, seed)}


def controls() -> list[dict[str, Any]]:
    rng = np.random.Generator(np.random.PCG64(MASTER))
    out = []
    for n, side in ((125, 5), (216, 6), (343, 7)):
        cands = {
            "T3": R.torus_lattice_3d(side),
            "RGG3_k12": R.rgg3_torus_binary(n, 12.0, rng),
            "ring_k6": R.ring_lattice(n, 6),
            "caveman_K8": __import__("omega.landscape.references", fromlist=["x"]).connected_caveman(n - n % 8, 8),
            "random_tree": R.random_tree(n, rng),
        }
        if n == 216:
            cands["tri_12x18"] = R.triangular_torus(12, 18)
        for name, adj in cands.items():
            out.append({"tag": "control", "name": name, "n": int(adj.shape[0]), **describe(adj, 0)})
    return out


def size_exponent(points: dict[int, list[float]]) -> float | None:
    ns = sorted(n for n in points if points[n])
    if len(ns) < 3:
        return None
    x = np.log(np.array(ns, dtype=float))
    y = np.log(np.array([np.median(points[n]) for n in ns]))
    slope = float(np.polyfit(x, y, 1)[0])
    return 1.0 / slope if slope > 1e-9 else float("inf")


def main() -> int:
    rt = json.loads((ROOT / "results" / "c0_redteam" / "summary.json").read_text())
    cand = set(rt["result"]["candidato_c0_cells"])
    l4 = read_rows(ROOT / "results" / "c0_dynamics" / "runs.jsonl")
    r6 = read_rows(ROOT / "results" / "c0_redteam" / "runs.jsonl")
    tasks: list[tuple[str, str, int, dict[str, Any]]] = []
    for x in l4:
        if x["kind"] == "generic" and x["cell_idx"] in cand and x["cls"]["class"] == "DISPERSO_LOCAL":
            meta = {"id": x["id"], "cell_idx": x["cell_idx"], "c_star": x["c_star"], "k_star": x["k_star"], "a": x["a"],
                    "init": x["init"], "seed": x["seed"], "n": 216, "status": x["status"]}
            tasks.append(("final", str(ROOT / "runs" / "c0" / x["npz"]), x["seed"], meta))
    for x in r6:
        if x["kind"] == "R6" and x["cell_idx"] in cand and x["cls"]["class"] == "DISPERSO_LOCAL":
            meta = {"id": x["id"], "cell_idx": x["cell_idx"], "c_star": x["c_star"], "k_star": x["k_star"], "a": x["a"],
                    "init": x["init"], "seed": x["seed"], "n": x["n"], "status": x["status"]}
            tasks.append(("final", str(ROOT / "runs" / "c0_redteam" / "out" / f"{x['id']}.npz"), x["seed"], meta))
    with mp.get_context("fork").Pool(4) as pool:
        rows = pool.map(_task, tasks)
    ctrl = controls()
    by: dict[tuple[int, str], dict[int, list[float]]] = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by[(r["cell_idx"], r["init"])][r["n"]].append(r["mean_hop"])
    scaling = [{"cell_idx": c, "init": i, "D_L": size_exponent(p), "mean_hop_by_n": {n: float(np.median(v)) for n, v in sorted(p.items())}}
               for (c, i), p in sorted(by.items())]
    cby: dict[str, dict[int, list[float]]] = defaultdict(lambda: defaultdict(list))
    for r in ctrl:
        cby[r["name"]][r["n"]].append(r["mean_hop"])
    ctrl_scaling = {name: size_exponent(p) for name, p in cby.items()}

    def med(key: str, sel: list[dict[str, Any]]) -> float | None:
        v = [float(r[key]) for r in sel if r.get(key) is not None and np.isfinite(r[key])]
        return float(np.median(v)) if v else None

    finals216 = [r for r in rows if r["n"] == 216]
    summary = {
        "step": "c0_battery", "code_commit": head_commit(), "candidate_cells": sorted(cand), "n_finals": len(rows),
        "finals_216": {k: med(k, finals216) for k in ("D_eff", "D_s", "hom_cv", "iso_median", "mean_hop", "diameter",
                                                      "max_clique", "cliques_per_node", "frac_nodes_in_clique_ge4", "mean_deg")},
        "D_eff_status_216": dict(Counter(r["D_eff_status"] for r in finals216)),
        "D_s_status_216": dict(Counter(r["D_s_status"] for r in finals216)),
        "size_scaling": scaling,
        "D_L_distribution": [s["D_L"] for s in scaling if s["D_L"] is not None],
        "controls": ctrl, "controls_D_L": ctrl_scaling,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    _atomic_write(OUT / "summary.json", json.dumps(summary, indent=2, default=str))
    with (OUT / "finals.jsonl").open("w", encoding="utf-8") as fh:
        for r in rows:
            fh.write(json.dumps(r, default=str) + "\n")
    print(json.dumps({k: summary[k] for k in ("n_finals", "finals_216", "D_eff_status_216", "controls_D_L")}, indent=1, default=str))
    print("D_L:", sorted(round(x, 2) for x in summary["D_L_distribution"] if x is not None and np.isfinite(x)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
