"""Omega-D L-OD-0 (docs/OMEGA_D_L0_PRERREGISTRO.md, seccion 4): solapamiento de difusion sobre un panel de grafos fijos, sin dinamica.

Panel 4.1 (N ~ 1000, clave PCG64 (20261010, familia, semilla)), medidas 4.2 por grafo, criterios K1-K3 de 4.3 y regla de decision.
Los finales C0 (N = 729) son informativos y no entran en K1-K3.
Salida: results/omega_d_l0/{graphs.jsonl, summary.json}. `--smoke`: 3 grafos pequenos -> results/omega_d_l0_smoke/.
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
from typing import Any, Callable  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from c0_dynamics import _atomic_write  # noqa: E402
from p1_diagnostic import patchwork, shortcuts, square_torus  # noqa: E402
from p1d2_panel import clique_lattice_3d  # noqa: E402

from omega.c0 import references as R  # noqa: E402
from omega.c0.locality import strong_support  # noqa: E402
from omega.controls.random_geometric import rgg_torus  # noqa: E402
from omega.diffusion.overlap import (  # noqa: E402
    DEFAULT_TAUS,
    distance2_contrast,
    edge_scaling,
    giant_dense,
    lazy_sym_operator,
    overlap_powers,
)
from omega.experiments.v11.gate import head_commit  # noqa: E402
from omega.landscape.references import clique_union, connected_caveman  # noqa: E402

MASTER = 20261010
SEEDS = (0, 1, 2)
N = 1000
OUT = ROOT / "results" / "omega_d_l0"
OUT_SMOKE = ROOT / "results" / "omega_d_l0_smoke"
_SPECS: list[dict[str, Any]] = []


def specs() -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []

    def add(fam: str, cls: str, gen: Callable[[np.random.Generator], np.ndarray], fid: int, seeds: tuple[int, ...] = SEEDS) -> None:
        for s in seeds:
            out.append({"id": f"{fam}_s{s}", "family": fam, "cls": cls, "gen": gen, "key": (MASTER, fid, s)})

    add("RGG3_k12", "G", lambda g: rgg_torus(N, 3, 12.0, g), 1)
    add("RGG2_k12", "G", lambda g: rgg_torus(N, 2, 12.0, g), 2)
    add("T3_10", "G", lambda g: R.torus_lattice_3d(10), 3, (0,))
    add("cuadrado_32", "G", lambda g: square_torus(32), 4, (0,))
    add("anillo_k12", "G", lambda g: R.ring_lattice(N, 12), 5, (0,))
    add("retazos_2x3", "A", lambda g: patchwork(N, 2, g), 6)
    add("retazos_3x3", "A", lambda g: patchwork(N, 3, g), 7)
    add("RGG3_atajos1", "A", lambda g: shortcuts(N, 0.01, g), 8)
    add("WS_b001", "A", lambda g: R.watts_strogatz(N, 12, 0.01, g), 9)
    add("cliques_K10", "R", lambda g: clique_union(N, 10), 10, (0,))
    add("caveman_K8", "R", lambda g: connected_caveman(N, 8), 11, (0,))
    add("reticulo_cliques_K8_L5", "R", lambda g: clique_lattice_3d(5, 8), 12, (0,))
    add("ER_k12", "NL", lambda g: R.erdos_renyi_m(N, N * 6, g), 13)
    add("RR_k12", "NL", lambda g: R.random_regular(N, 12, g), 14)
    for p in sorted((ROOT / "runs" / "c0_f1" / "out").glob("F1_c*_R_n729_s*.npz")):
        out.append({"id": "C0_" + p.stem, "family": "C0_R_n729", "cls": "INFO", "npz": str(p), "key": None})
    return out


def specs_smoke() -> list[dict[str, Any]]:
    n = 216
    return [
        {"id": "RGG3_k12_s0", "family": "RGG3_k12", "cls": "G", "gen": lambda g: rgg_torus(n, 3, 12.0, g), "key": (MASTER, 1, 0)},
        {"id": "retazos_2x3_s0", "family": "retazos_2x3", "cls": "A", "gen": lambda g: patchwork(n, 2, g), "key": (MASTER, 6, 0)},
        {"id": "cliques_K10_s0", "family": "cliques_K10", "cls": "R", "gen": lambda g: clique_union(200, 10), "key": (MASTER, 10, 0)},
    ]


def measure(adj: np.ndarray, key: tuple[int, ...] | None) -> dict[str, Any]:
    a = giant_dense(adj)
    al = overlap_powers(lazy_sym_operator(a), DEFAULT_TAUS)
    iu, ju = np.nonzero(np.triu(a, 1))
    es = edge_scaling(al, iu, ju, (2, 4, 8))
    rng = R.rng_from_key((*key, 2) if key is not None else (MASTER, 99, 99))  # flujo propio para el muestreo a distancia 2
    d2 = distance2_contrast(a, al, rng, 5000)
    n = int(a.shape[0])
    return {"n_giant": n, "mean_degree": float(a.sum() / n), **es["summary"], "dist2": d2}


def run_idx(i: int) -> dict[str, Any]:
    spec = _SPECS[i]
    if "npz" in spec:
        z = np.load(spec["npz"])
        adj = strong_support(np.asarray(z[z.files[0]], dtype=np.float64)).astype(np.float64)
    else:
        adj = spec["gen"](R.rng_from_key(spec["key"]))
    t0 = time.time()
    row = {k: v for k, v in spec.items() if k not in ("gen", "npz")}
    row.update(measure(adj, spec["key"]))
    row["seconds"] = time.time() - t0
    return row


def _med(xs: list[float]) -> float:
    return float(np.median(xs)) if xs else float("nan")


def evaluate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    by: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        by[r["family"]].append(r)
    table = {
        f: {
            "cls": rs[0]["cls"], "n_graphs": len(rs), "n_giant": [r["n_giant"] for r in rs],
            "median_s": [r["median_s"] for r in rs], "f_nd": [r["f_nd"] for r in rs],
            "frac_degenerate": [r["frac_degenerate"] for r in rs],
            "median_Q8": [r["median_Q"]["8"] for r in rs],
            "median_alpha_dist2_t8": [r["dist2"]["median_alpha"].get("8") for r in rs],
        }
        for f, rs in by.items()
    }
    res: dict[str, Any] = {}
    # K1
    rgg = by.get("RGG3_k12", [])
    if rgg and "retazos_2x3" in by and "retazos_3x3" in by:
        base = _med([r["f_nd"] for r in rgg])
        k1fam = {f: [r["f_nd"] for r in by[f]] for f in ("retazos_2x3", "retazos_3x3")}
        ok = all(len(v) == 3 and all(x >= base + 0.05 for x in v) for v in k1fam.values())
        res["K1"] = {"pass": bool(ok), "f_nd_RGG3_median": base, "threshold": base + 0.05, "f_nd_retazos": k1fam,
                     "informative": {f: [r["f_nd"] for r in by.get(f, [])] for f in ("RGG3_atajos1", "WS_b001")}}
    else:
        res["K1"] = None
    # K2
    g_rows = [r for r in rows if r["cls"] == "G"]
    nl_rows = [r for r in rows if r["cls"] == "NL"]
    if g_rows and nl_rows:
        g_s = {r["id"]: r["median_s"] for r in g_rows}
        nl_s = {r["id"]: r["median_s"] for r in nl_rows}
        g_ok = all(abs(v) <= 0.15 for v in g_s.values())
        nl_ok = all(v <= -0.5 for v in nl_s.values())
        res["K2"] = {"pass": bool(g_ok and nl_ok), "G_pass": bool(g_ok), "NL_pass": bool(nl_ok), "G_median_s": g_s, "NL_median_s": nl_s,
                     "NL_family_median_s": {f: _med([r["median_s"] for r in by[f]]) for f in ("ER_k12", "RR_k12") if f in by}}
    else:
        res["K2"] = None
    # K3
    if rgg and "cliques_K10" in by and "caveman_K8" in by:
        q_ref = _med([r["median_Q"]["8"] for r in rgg])
        k3 = {f: by[f][0]["median_Q"]["8"] for f in ("cliques_K10", "caveman_K8")}
        res["K3"] = {"pass": bool(all(v <= 0.1 * q_ref for v in k3.values())), "Q8_RGG3_median": q_ref, "threshold": 0.1 * q_ref,
                     "Q8": k3}
    else:
        res["K3"] = None
    k1, k2 = res["K1"], res["K2"]
    if k1 is None or k2 is None:
        dec = "INCOMPLETO (panel parcial)"
    elif not k1["pass"]:
        dec = "K1 FALLA: Omega-D muere en su premisa; pasar a P3"
    elif not k2["pass"]:
        dec = "K1 PASA, K2 FALLA: Omega-D en suspenso; decide el Consejo"
    else:
        dec = "K1 y K2 PASAN: definir la accion candidata en L-OD-1 (otro prerregistro)"
    res["decision"] = dec
    res["family_table"] = table
    res["lattice_vs_half_d"] = {}
    for fam, d in (("T3_10", 3), ("cuadrado_32", 2), ("anillo_k12", 1)):
        if fam in by:
            r = by[fam][0]
            q = {t: v for t, v in r["median_Q"].items()}
            res["lattice_vs_half_d"][fam] = {
                "d": d, "d_over_2": d / 2.0, "median_Q": q,
                "rel_err_tau8": q["8"] / (d / 2.0) - 1.0,
                "note": "anillo k12 no es vecino-mas-proximo: d/2 no aplica, solo informativo" if fam == "anillo_k12" else "",
            }
    return res


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args(argv)
    out = OUT_SMOKE if args.smoke else OUT
    _SPECS[:] = specs_smoke() if args.smoke else specs()
    t0 = time.time()
    with mp.get_context("fork").Pool(args.procs) as pool:
        rows = pool.map(run_idx, range(len(_SPECS)), chunksize=1)
    out.mkdir(parents=True, exist_ok=True)
    _atomic_write(out / "graphs.jsonl", "".join(json.dumps(r, allow_nan=True) + "\n" for r in rows))
    res = evaluate([r for r in rows if r["cls"] != "INFO"])
    c0 = [r for r in rows if r["cls"] == "INFO"]
    res["C0_informative"] = {"n_graphs": len(c0), "f_nd": [r["f_nd"] for r in c0], "median_s": [r["median_s"] for r in c0],
                             "ids": [r["id"] for r in c0]}
    summary = {"step": "omega_d_l0", "smoke": args.smoke, "code_commit": head_commit(), "n_graphs": len(rows),
               "seconds": time.time() - t0, "result": res}
    _atomic_write(out / "summary.json", json.dumps(summary, indent=2, allow_nan=True) + "\n")
    for k in ("K1", "K2", "K3"):
        print(k, None if res[k] is None else res[k]["pass"])
    print("decision:", res["decision"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
