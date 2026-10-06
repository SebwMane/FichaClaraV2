"""P1-D (docs/OMEGA_P1D_PRERREGISTRO.md): ¿la firma CV(|B_r|) decreciente distingue geometria de localidad?

Salida: results/p1_diagnostic/{summary.json, graphs.jsonl}. Uso: python tools/p1_diagnostic.py [--procs 4]
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
from collections import Counter, defaultdict  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any, Callable  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from c0_dynamics import _atomic_write, read_rows  # noqa: E402

from omega.c0 import references as R  # noqa: E402
from omega.c0.locality import strong_support  # noqa: E402
from omega.controls.random_geometric import rgg_torus  # noqa: E402
from omega.diagnostics.ball_growth import ball_profile  # noqa: E402
from omega.experiments.v11.gate import head_commit  # noqa: E402

MASTER = 20261007
SEEDS = (0, 1, 2)
OUT = ROOT / "results" / "p1_diagnostic"


def rgg_points(n: int, dim: int, k: float, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    vol = math.pi ** (dim / 2.0) / math.gamma(dim / 2.0 + 1.0)
    r = (k / ((n - 1) * vol)) ** (1.0 / dim)
    x = rng.random((n, dim))
    diff = np.abs(x[:, None, :] - x[None, :, :])
    diff = np.minimum(diff, 1.0 - diff)
    a = np.sqrt((diff**2).sum(axis=2)) < r
    np.fill_diagonal(a, False)
    return x, a.astype(np.float64)


def patchwork(n: int, blocks: int, rng: np.random.Generator) -> np.ndarray:
    """RGG3 k12; los extremos de aristas entre bloques distintos se reemparejan al azar (grados preservados)."""
    x, a = rgg_points(n, 3, 12.0, rng)
    bid = np.floor(x * blocks).astype(int) @ np.array([blocks * blocks, blocks, 1])
    iu, ju = np.nonzero(np.triu(a, 1))
    cross = bid[iu] != bid[ju]
    out = a.copy()
    out[iu[cross], ju[cross]] = 0.0
    out[ju[cross], iu[cross]] = 0.0
    stubs = np.concatenate([iu[cross], ju[cross]])
    rng.shuffle(stubs)
    for u, v in zip(stubs[0::2], stubs[1::2]):
        if u != v and out[u, v] == 0.0:
            out[u, v] = out[v, u] = 1.0
    return out


def shortcuts(n: int, frac: float, rng: np.random.Generator) -> np.ndarray:
    _, a = rgg_points(n, 3, 12.0, rng)
    iu, ju = np.nonzero(np.triu(a, 1))
    pick = rng.random(iu.size) < frac
    out = a.copy()
    for u, v in zip(iu[pick], ju[pick]):
        w = int(rng.integers(n))
        if w != u and out[u, w] == 0.0:
            out[u, v] = out[v, u] = 0.0
            out[u, w] = out[w, u] = 1.0
    return out


def diluted(n: int, rng: np.random.Generator) -> np.ndarray:
    a = rgg_torus(n, 3, 16.0, rng)
    iu, ju = np.nonzero(np.triu(a, 1))
    drop = rng.random(iu.size) < 0.25
    a[iu[drop], ju[drop]] = 0.0
    a[ju[drop], iu[drop]] = 0.0
    return a


def square_torus(side: int) -> np.ndarray:
    n = side * side
    a = np.zeros((n, n))
    for i in range(side):
        for j in range(side):
            v = i * side + j
            for di, dj in ((1, 0), (0, 1)):
                u = ((i + di) % side) * side + (j + dj) % side
                a[v, u] = a[u, v] = 1.0
    return a


def specs() -> list[dict[str, Any]]:
    """(id, etiqueta, panel, familia, generador sin argumentos -> adyacencia)."""
    out: list[dict[str, Any]] = []

    def add(fam: str, label: str, gen: Callable[[np.random.Generator], np.ndarray], fid: int, seeds: tuple[int, ...] = SEEDS,
            panel: str = "principal") -> None:
        for s in seeds:
            out.append({"id": f"{fam}_s{s}", "family": fam, "label": label, "panel": panel, "gen": gen, "key": (MASTER, fid, s)})

    add("RGG3_k8", "G", lambda g: rgg_torus(729, 3, 8.0, g), 1)
    add("RGG3_k12", "G", lambda g: rgg_torus(729, 3, 12.0, g), 2)
    add("RGG3_k16", "G", lambda g: rgg_torus(729, 3, 16.0, g), 3)
    add("RGG3_k12_n1728", "G", lambda g: rgg_torus(1728, 3, 12.0, g), 4)
    add("RGG2_k12", "G", lambda g: rgg_torus(729, 2, 12.0, g), 5)
    add("RGG4_k12", "G", lambda g: rgg_torus(729, 4, 12.0, g), 6)
    add("RGG3_diluido", "G", lambda g: diluted(729, g), 7)
    add("T3_9", "G", lambda g: R.torus_lattice_3d(9), 8, (0,))
    add("cuadrado_27", "G", lambda g: square_torus(27), 9, (0,))
    add("triangular_27", "G", lambda g: R.triangular_torus(27, 27), 10, (0,))
    add("anillo_k12", "G", lambda g: R.ring_lattice(729, 12), 11, (0,))
    add("retazos_27", "L", lambda g: patchwork(729, 3, g), 12)
    add("retazos_8", "L", lambda g: patchwork(729, 2, g), 13)
    add("retazos_64_n1728", "L", lambda g: patchwork(1728, 4, g), 14)
    add("RGG3_atajos10", "L", lambda g: shortcuts(729, 0.10, g), 15)
    add("WS_b005", "L", lambda g: R.watts_strogatz(729, 12, 0.05, g), 16)
    add("WS_b02", "L", lambda g: R.watts_strogatz(729, 12, 0.2, g), 17)
    add("caveman_K8", "L", lambda g: __import__("omega.landscape.references", fromlist=["x"]).connected_caveman(728, 8), 18, (0,))
    add("ER_k12", "NL", lambda g: R.erdos_renyi_m(729, 729 * 6, g), 19)
    add("RR_k12", "NL", lambda g: R.random_regular(729, 12, g), 20)
    add("RGG3_k12_n343", "G", lambda g: rgg_torus(343, 3, 12.0, g), 21, panel="secundario_C0")
    add("RGG2_k12_n343", "G", lambda g: rgg_torus(343, 2, 12.0, g), 22, panel="secundario_C0")
    for x in read_rows(ROOT / "results" / "c0_redteam" / "runs.jsonl"):
        if x["kind"] == "R6" and x["n"] == 343 and x["seed"] == 0 and x["cls"]["class"] == "DISPERSO_LOCAL":
            path = ROOT / "runs" / "c0_redteam" / "out" / f"{x['id']}.npz"
            out.append({"id": x["id"], "family": f"C0_{x['init']}_n343", "label": "L", "panel": "secundario_C0",
                        "npz": str(path), "key": None})
    return out


_SPECS: list[dict[str, Any]] = []


def run_idx(i: int) -> dict[str, Any]:
    return run(_SPECS[i])


def run(spec: dict[str, Any]) -> dict[str, Any]:
    if "npz" in spec:
        z = np.load(spec["npz"])
        adj = strong_support(np.asarray(z[z.files[0]], dtype=np.float64)).astype(np.float64)
    else:
        adj = spec["gen"](R.rng_from_key(spec["key"]))
    prof = ball_profile(adj)
    return {k: v for k, v in spec.items() if k not in ("gen", "npz")} | {"mean_degree": float(np.asarray(adj > 0).sum() / adj.shape[0]),
                                                                         **prof}


def verdict(rows: list[dict[str, Any]], tag: str) -> dict[str, Any]:
    g = [r for r in rows if r["label"] == "G" and r[tag]["status"] in ("HOMOGENEIZA", "NO_HOMOGENEIZA")]
    nl = [r for r in rows if r["label"] in ("L", "NL") and r[tag]["status"] in ("HOMOGENEIZA", "NO_HOMOGENEIZA")]
    sens = sum(r[tag]["status"] == "HOMOGENEIZA" for r in g) / len(g) if g else None
    spec = sum(r[tag]["status"] == "NO_HOMOGENEIZA" for r in nl) / len(nl) if nl else None
    if len(g) < 5 or len(nl) < 5 or sens is None or spec is None:
        v = "INDETERMINADO"
    elif sens >= 0.9 and spec >= 0.9:
        v = "SOBREVIVE"
    elif sens < 0.7 or spec < 0.7:
        v = "MUERE"
    else:
        v = "PARCIAL"
    fam: dict[str, Counter[str]] = defaultdict(Counter)
    rho: dict[str, list[float | None]] = defaultdict(list)
    for r in rows:
        fam[r["family"]][r[tag]["status"]] += 1
        rho[r["family"]].append(None if r[tag]["rho"] is None else round(r[tag]["rho"], 3))
    return {"verdict": v, "sensitivity": sens, "specificity": spec, "n_G_eval": len(g), "n_L_NL_eval": len(nl),
            "by_family": {f: {"label": next(r["label"] for r in rows if r["family"] == f), "status": dict(c), "rho": rho[f]}
                          for f, c in fam.items()}}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--procs", type=int, default=4)
    args = ap.parse_args(argv)
    _SPECS[:] = specs()
    with mp.get_context("fork").Pool(args.procs) as pool:
        rows = pool.map(run_idx, range(len(_SPECS)), chunksize=1)
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "graphs.jsonl").open("w", encoding="utf-8") as fh:
        for r in rows:
            fh.write(json.dumps(r, default=str) + "\n")
    res = {}
    for panel in ("principal", "secundario_C0"):
        sel = [r for r in rows if r["panel"] == panel]
        res[panel] = {"primario_q4": verdict(sel, "q4"), "secundario_q2": verdict(sel, "q2")}
    _atomic_write(OUT / "summary.json", json.dumps({"step": "p1_diagnostic", "code_commit": head_commit(),
                                                    "n_graphs": len(rows), "result": res}, indent=2, default=str))
    for panel, rr in res.items():
        for tag, v in rr.items():
            print(panel, tag, v["verdict"], "sens", v["sensitivity"], "spec", v["specificity"], v["n_G_eval"], v["n_L_NL_eval"])
        for f, d in rr["primario_q4"]["by_family"].items():
            print(f"   {f:<22} {d['label']:<3} {d['status']} rho={d['rho'][:6]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
