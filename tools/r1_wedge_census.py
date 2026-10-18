"""Omega R1-0: censo exacto de vecinos comunes t(u,w) en aristas y cunas (pares a distancia 2) para el analisis de puntos
fijos de reglas de reescritura por umbral entero (docs/OMEGA_R1_0.md §3). No es dinamica: verifica afirmaciones estaticas.
Uso: python tools/r1_wedge_census.py -> results/r1_0/wedge_census.json
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from a0_rc2 import triangular_torus  # noqa: E402
from c0_dynamics import _atomic_write  # noqa: E402
from cic_l0b import random_two_tree  # noqa: E402
from coh_l0b import prufer_tree  # noqa: E402
from p1d3_panel import erdos_renyi, random_regular, rgg, torus_lattice  # noqa: E402

from omega.c0.references import rng_from_key  # noqa: E402


def census(adj: sparse.csr_array) -> dict[str, dict[str, float]]:
    a = sparse.csr_array(adj, dtype=np.int64)
    a.data[:] = 1
    common = (a @ a).tocsr()  # common[u, w] = t(u, w)
    common.setdiag(0)
    common.eliminate_zeros()
    edges = sparse.triu(a, k=1).tocoo()
    t_edge = np.asarray(common[edges.row, edges.col]).ravel()
    w = sparse.triu(common, k=1).tocoo()
    is_edge = np.asarray(a[w.row, w.col]).ravel() > 0
    t_wedge = w.data[~is_edge]  # pares no adyacentes con >= 1 vecino comun (distancia 2)

    def dist(x: np.ndarray) -> dict[str, float]:
        c = Counter(np.minimum(x, 9).tolist())
        n = max(1, x.size)
        return {("9+" if k == 9 else str(k)): round(v / n, 4) for k, v in sorted(c.items())}

    return {"edges": dist(t_edge), "wedges": dist(t_wedge)}


def main() -> int:
    key = (20261017, 0)
    g = {
        "Z2 (cuadrado 30^2)": torus_lattice(30, 2),
        "Z3 (T3 12^3)": torus_lattice(12, 3),
        "Z4 (T4 7^4)": torus_lattice(7, 4),
        "triangular 30^2": triangular_torus(30),
        "RGG3 k12": rgg(3000, 3, 12.0, rng_from_key(key + (1,))),
        "arbol Prufer": prufer_tree(3000, rng_from_key(key + (2,))),
        "ER k4": erdos_renyi(3000, 6000, rng_from_key(key + (3,))),
        "RR k6": random_regular(3000, 6, rng_from_key(key + (4,))),
        "2-arbol": random_two_tree(3000, rng_from_key(key + (5,))),
    }
    out = {name: census(a) for name, a in g.items()}
    for name, c in out.items():
        print(f"{name:<20} aristas t: {c['edges']}   cunas t: {c['wedges']}")
    _atomic_write(ROOT / "results" / "r1_0" / "wedge_census.json", json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
