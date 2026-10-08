"""Fase 2 (DESCRIPTIVO): crecimiento de bolas |B_r(i)| en escala intermedia.

Pregunta: si los conteos locales (k, c, q) no separan los finales C0 de N=729 de una RGG3 genuina
(results/c1_local_stats), ¿lo hace la regularidad del crecimiento de bolas?
Por nodo: |B_r| para r = 1..R_MAX y exponente local d_i = log(|B_4|/|B_2|)/log 2.
Por grafo: mediana y CV entre nodos de |B_r| y de d_i. No es una prueba de C1 ni un umbral nuevo.
Salida: results/c1_ball_growth/summary.json. Uso: python tools/c1_ball_growth.py
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
from omega.geometry.distances import hop_distance_matrix  # noqa: E402

MASTER = 20261006
R_MAX = 5


def growth(adj: np.ndarray) -> dict[str, Any]:
    a = np.asarray(adj) > 0
    g = giant_component(a)
    a = np.asarray(a[np.ix_(g, g)], dtype=np.bool_)
    d = hop_distance_matrix(a)
    balls = np.stack([(d <= r).sum(axis=1) for r in range(1, R_MAX + 1)], axis=1).astype(float)
    dloc = np.log(balls[:, 3] / balls[:, 1]) / np.log(2.0)

    def st(x: np.ndarray) -> dict[str, float]:
        return {"median": float(np.median(x)), "cv": float(x.std() / x.mean())}

    return {"n": int(a.shape[0]), "B_r": {str(r + 1): st(balls[:, r]) for r in range(R_MAX)}, "d_local_2_4": st(dloc)}


def main() -> int:
    rng = R.rng_from_key((MASTER, 91))
    graphs: dict[str, np.ndarray] = {"T3_9": R.torus_lattice_3d(9)}
    for s in range(3):
        graphs[f"RGG3_k12_n729_s{s}"] = R.rgg3_torus_binary(729, 12.0, rng)
    graphs["ER_k12_n729"] = R.erdos_renyi_m(729, 729 * 6, rng)
    for c, s in itertools.product((18, 19, 20, 36, 37, 38, 55, 56, 73, 74), (0, 1, 2)):
        path = ROOT / "runs" / "c0_f1" / "out" / f"F1_c{c}_R_n729_s{s}.npz"
        if path.exists():
            z = np.load(path)
            graphs[f"C0_c{c}_s{s}"] = strong_support(np.asarray(z[z.files[0]], dtype=np.float64))
    res = {name: growth(a) for name, a in graphs.items()}
    dest = ROOT / "results" / "c1_ball_growth"
    dest.mkdir(parents=True, exist_ok=True)
    _atomic_write(dest / "summary.json", json.dumps({"step": "c1_ball_growth", "code_commit": head_commit(),
                                                      "descriptive_only": True, "graphs": res}, indent=2))
    for name, st in res.items():
        b = st["B_r"]
        print(f"{name:<20} " + " ".join(f"B{r}={b[str(r)]['median']:6.0f}(cv {b[str(r)]['cv']:.2f})" for r in (2, 3, 4))
              + f"  d24={st['d_local_2_4']['median']:.2f}(cv {st['d_local_2_4']['cv']:.2f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
