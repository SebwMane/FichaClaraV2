"""Fase 2, §3 (informacion de apoyo, DESCRIPTIVA): conteos locales por arista (k, c, q) en referencias y finales C0.

c_ij = vecinos comunes; q_ij = 4-ciclos que contienen la arista = (A^3)_ij - k_i - k_j + 1.
No es una prueba de C1: solo dice que conteos tendria que distinguir cualquier C1 local de segundo orden.
Salida: results/c1_local_stats/summary.json. Uso: python tools/c1_local_stats.py
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import itertools  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from c0_dynamics import _atomic_write  # noqa: E402

from omega.c0 import references as R  # noqa: E402
from omega.c0.locality import giant_component, strong_support  # noqa: E402
from omega.experiments.v11.gate import head_commit  # noqa: E402

MASTER = 20261006


def hypercube(d: int) -> np.ndarray:
    n = 2**d
    a = np.zeros((n, n))
    for v in range(n):
        for b in range(d):
            a[v, v ^ (1 << b)] = 1.0
    return a


def complete_bipartite(m: int) -> np.ndarray:
    a = np.zeros((2 * m, 2 * m))
    a[:m, m:] = 1.0
    a[m:, :m] = 1.0
    return a


def edge_stats(adj: np.ndarray) -> dict[str, Any]:
    a = (np.asarray(adj) > 0).astype(np.float64)
    np.fill_diagonal(a, 0.0)
    k = a.sum(axis=1)
    a2 = a @ a
    a3 = a2 @ a
    iu = np.triu_indices(a.shape[0], 1)
    e = a[iu] > 0
    i, j = iu[0][e], iu[1][e]
    c = a2[i, j]
    q = a3[i, j] - k[i] - k[j] + 1.0

    def d(x: np.ndarray) -> dict[str, float]:
        vals, cnt = np.unique(np.round(x, 6), return_counts=True)
        return {"mean": float(x.mean()), "cv": float(x.std() / x.mean()) if x.mean() > 0 else 0.0,
                "mode": float(vals[np.argmax(cnt)]), "frac_mode": float(cnt.max() / x.size)}

    pairs = np.stack([c, q], axis=1)
    _, pc = np.unique(pairs, axis=0, return_counts=True)
    return {"n": int(a.shape[0]), "m": int(e.sum()), "k": d(k), "c": d(c), "q": d(q),
            "frac_edges_modal_cq": float(pc.max() / e.sum())}


def main() -> int:
    rng = R.rng_from_key((MASTER, 90))
    refs: dict[str, np.ndarray] = {
        "T3_9": R.torus_lattice_3d(9),
        "T3_8": R.torus_lattice_3d(8),
        "Q6": hypercube(6),
        "Q9": hypercube(9),
        "K_6,6": complete_bipartite(6),
        "tri_12x18": R.triangular_torus(12, 18),
        "RGG3_k12_n729": R.rgg3_torus_binary(729, 12.0, rng),
        "ER_k12_n729": R.erdos_renyi_m(729, 729 * 6, rng),
        "RR_k6_n729": R.random_regular(729, 6, rng),
    }
    out: dict[str, Any] = {"references": {name: edge_stats(a) for name, a in refs.items()}, "c0_R_n729": {}}
    for c, s in itertools.product((18, 19, 20, 36, 37, 38, 55, 56, 73, 74), (0, 1, 2)):
        path = ROOT / "runs" / "c0_f1" / "out" / f"F1_c{c}_R_n729_s{s}.npz"
        if not path.exists():
            continue
        z = np.load(path)
        a = strong_support(np.asarray(z[z.files[0]], dtype=np.float64))
        g = giant_component(a)
        out["c0_R_n729"][f"c{c}_s{s}"] = edge_stats(np.asarray(a[np.ix_(g, g)], dtype=np.float64))
    summary = {"step": "c1_local_stats", "code_commit": head_commit(), "descriptive_only": True, **out}
    dest = ROOT / "results" / "c1_local_stats"
    dest.mkdir(parents=True, exist_ok=True)
    _atomic_write(dest / "summary.json", json.dumps(summary, indent=2))
    fmt = "{:<16} n={:<5} k={:5.2f} (cv {:.2f})  c={:5.2f} (cv {:.2f})  q={:6.2f} (cv {:.2f})  modal(c,q)={:.2f}"
    for name, st in list(out["references"].items()) + list(out["c0_R_n729"].items()):
        print(fmt.format(name, st["n"], st["k"]["mean"], st["k"]["cv"], st["c"]["mean"], st["c"]["cv"],
                         st["q"]["mean"], st["q"]["cv"], st["frac_edges_modal_cq"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
