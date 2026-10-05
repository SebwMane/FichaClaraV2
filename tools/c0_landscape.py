"""C0-L1, C0-L2 y C0-L3 (prerregistro Omega-C0 §2-§4): paisaje extremal de S_C0 sobre referencias, N = 216.

Malla: c* {0,1,2,4,8} x k* {4,6,8,12,16} x a {0.25,1,4} (75 celdas). Cada referencia se evalua con amplitud optima
(t en (0,1], rejilla 2000 + refinamiento), salvo la union de cliques BINARIAS (amplitud 1 exacta).
Triviales: vacio, uniforme, clique_binaria, clique_difusa, bipartita. El resto es no trivial.

  L1: MUERTE-L1 si la mejor clique binaria es la mejor referencia en >= 50% de las celdas.
  L2: celda rota si S_nt < S_triv - 0.01|LB|; empate si |S_nt - S_triv| <= 0.01|LB| sin mejorar. Global: ROTA >= 10% rotas;
      EMPATE si no y >= 10% rotas o empatadas; NO ROTA en otro caso.
  L3: KKT en la caja de T3 y de RGG3 (k 6,8,12,16; semillas 0-4) con amplitud 1 y optima; INDICIO NEGATIVO si ninguna lo es.

Salida: <out>/summary.json, points.csv, log.txt. Modo full exige arbol git limpio (--allow-dirty solo con --smoke).
"""

from __future__ import annotations

import os

# Un hilo BLAS por proceso: el paralelismo es por corrida (4 procesos).
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import csv
import hashlib
import json
import multiprocessing as mp
import os
import sys
import time
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from omega.c0 import references as R  # noqa: E402
from omega.c0.functional import C0Params, action_c0, kkt_box, params_from_targets  # noqa: E402
from omega.experiments.v11.gate import head_commit, tree_dirty  # noqa: E402
from omega.io.provenance import dependency_lock  # noqa: E402
from omega.phases.scan import to_jsonable  # noqa: E402
from omega.types import FloatArray  # noqa: E402

N = 216
STEP = "c0_landscape"
TRIVIAL = ("vacio", "uniforme", "clique_binaria", "clique_difusa", "bipartita")
FAMILY_ID = {"rgg3": 1, "rr": 2, "er": 3, "ws": 4, "tree": 5}

FULL: dict[str, Any] = {
    "n": N, "c_stars": [0, 1, 2, 4, 8], "k_stars": [4, 6, 8, 12, 16], "a_values": [0.25, 1.0, 4.0],
    "seeds": [0, 1, 2, 3, 4], "rgg_k": [6, 8, 12, 16], "rr_k": [4, 6, 8, 12, 16], "er_k": [4, 6, 8, 12, 16],
    "ws": {"k": 6, "beta": 0.1}, "ring_k": 6, "tri_torus": [12, 18], "amp_grid": R.AMP_GRID,
    "l2_tol_rel": 0.01, "l2_global_frac": 0.10, "l1_global_frac": 0.50, "kkt_tol_rel": 1e-9, "smoke_cells": 3,
    "master_entropy": 20261005,
}


def canonical_hash(cfg: dict[str, Any]) -> str:
    text = json.dumps(to_jsonable(cfg), sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def divisors(n: int) -> list[int]:
    return [s for s in range(2, n + 1) if n % s == 0]


def build_config(smoke: bool) -> dict[str, Any]:
    cfg = dict(FULL)
    if smoke:
        cfg["seeds"] = [0]
    cfg["smoke"] = smoke
    return cfg


def cell_grid(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    cells = [
        {"cell_idx": i, "c_star": c, "k_star": k, "a": a}
        for i, (c, k, a) in enumerate(
            (c, k, a) for c in cfg["c_stars"] for k in cfg["k_stars"] for a in cfg["a_values"]
        )
    ]
    if cfg["smoke"]:
        # c*=0,k*=6,a=1 (favorable a T3), una celda con c*>0 y una con k* grande
        pick = [(0, 6, 1.0), (2, 8, 1.0), (8, 16, 4.0)]
        cells = [c for c in cells if (c["c_star"], c["k_star"], c["a"]) in pick][: cfg["smoke_cells"]]
    return cells


# ------------------------------------------------------------------ grafos (cache por proceso)


@lru_cache(maxsize=None)
def _adj(name: str, par: int, seed: int) -> FloatArray:
    """Adyacencia 0/1 de una referencia; (name, par, seed) la identifican."""
    if name == "t3":
        return R.torus_lattice_3d(6)
    if name == "ring":
        return R.ring_lattice(N, par)
    if name == "tri":
        return R.triangular_torus(12, 18)
    if name == "bip":
        return R.complete_bipartite_half(N)
    if name == "uni":
        return R.uniform_complete(N)
    if name == "clique":
        return R.clique_union(N, par)
    rng = R.rng_from_key((20261005, FAMILY_ID[name.split("_")[0]], par, seed))
    if name == "rgg3":
        return R.rgg3_torus_binary(N, float(par), rng)
    if name == "rr":
        return R.random_regular(N, par, rng)
    if name == "er":
        return R.erdos_renyi_m(N, int(round(N * par / 2)), rng)
    if name == "ws":
        return R.watts_strogatz(N, par, 0.1, rng)
    if name == "tree":
        return R.random_tree(N, rng)
    raise KeyError(name)


@lru_cache(maxsize=None)
def _prof(name: str, par: int, seed: int) -> R.AmpProfile:
    return R.amplitude_profile(_adj(name, par, seed))


def reference_list(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    """Todas las referencias: nombre, familia, grupo, (name, par, seed) y si usa amplitud optima."""
    refs: list[dict[str, Any]] = []

    def add(ref: str, fam: str, group: str, key: tuple[str, int, int], opt: bool = True) -> None:
        refs.append({"ref": ref, "family": fam, "group": group, "key": key, "opt": opt})

    for s in divisors(N):
        add(f"clique_binaria_s{s}", "clique_binaria", "trivial", ("clique", s, 0), opt=False)
        add(f"clique_difusa_s{s}", "clique_difusa", "trivial", ("clique", s, 0))
    add("bipartita", "bipartita", "trivial", ("bip", 0, 0))
    add("uniforme", "uniforme", "trivial", ("uni", 0, 0))
    add("t3", "t3", "no_trivial", ("t3", 0, 0))
    add("tree", "tree", "no_trivial", ("tree", 0, 0))
    add(f"ring_k{cfg['ring_k']}", "ring", "no_trivial", ("ring", cfg["ring_k"], 0))
    add(f"ws_k{cfg['ws']['k']}_b0.1", "ws", "no_trivial", ("ws", cfg["ws"]["k"], 0))
    add("tri_torus_12x18", "tri_torus", "no_trivial", ("tri", 0, 0))
    for fam, ks in (("rgg3", cfg["rgg_k"]), ("rr", cfg["rr_k"]), ("er", cfg["er_k"])):
        for k in ks:
            for s in cfg["seeds"]:
                add(f"{fam}_k{k}_s{s}", f"{fam}_k{k}", "no_trivial", (fam, k, s))
    return refs


# ------------------------------------------------------------------ una celda


def eval_cell(cell: dict[str, Any], cfg: dict[str, Any]) -> dict[str, Any]:
    t0 = time.perf_counter()
    p = params_from_targets(cell["c_star"], cell["k_star"], cell["a"])
    lb = p.lower_bound(N)
    rows: list[dict[str, Any]] = [
        {"ref": "vacio", "family": "vacio", "group": "trivial", "t": 0.0, "S": 0.0, "seed": None}
    ]
    for r in reference_list(cfg):
        name, par, seed = r["key"]
        if r["opt"]:
            t, s = R.optimal_amplitude_from_profile(_prof(name, par, seed), p)
        else:
            t, s = 1.0, action_c0(_adj(name, par, seed), p)
        rows.append({"ref": r["ref"], "family": r["family"], "group": r["group"], "t": t, "S": s, "seed": seed if name in FAMILY_ID_SEEDED else None})
    for row in rows:
        row["S_over_LB"] = row["S"] / lb
    best = min(rows, key=lambda r: r["S"])
    triv = min((r for r in rows if r["group"] == "trivial"), key=lambda r: r["S"])
    nt = min((r for r in rows if r["group"] == "no_trivial"), key=lambda r: r["S"])
    best_cb = min((r for r in rows if r["family"] == "clique_binaria"), key=lambda r: r["S"])
    others = min(r["S"] for r in rows if r["family"] != "clique_binaria")
    tol = cfg["l2_tol_rel"] * abs(lb)
    l1 = best_cb["S"] <= others + 1e-9 * abs(lb)
    l1_strict = best_cb["S"] < others - 1e-9 * abs(lb)
    if nt["S"] < triv["S"] - tol:
        l2 = "rota"
    elif abs(nt["S"] - triv["S"]) <= tol:
        l2 = "empate"
    else:
        l2 = "no"
    fam_min: dict[str, dict[str, Any]] = {}
    for r in rows:
        f = fam_min.get(r["family"])
        if f is None or r["S"] < f["S"]:
            fam_min[r["family"]] = {"S": r["S"], "ref": r["ref"], "t": r["t"]}
    l3 = kkt_cell(p, cfg)
    return {
        "cell": cell, "lambda": p.lam, "kappa": p.kappa, "psi_star": p.psi_star, "LB": lb,
        "best": {"ref": best["ref"], "S": best["S"], "S_over_LB": best["S"] / lb},
        "best_trivial": {"ref": triv["ref"], "S": triv["S"], "S_over_LB": triv["S"] / lb},
        "best_no_trivial": {"ref": nt["ref"], "S": nt["S"], "S_over_LB": nt["S"] / lb},
        "best_clique_binaria": {"ref": best_cb["ref"], "S": best_cb["S"], "S_over_LB": best_cb["S"] / lb},
        "L1": bool(l1), "L1_strict": bool(l1_strict), "L2": l2,
        "gap_nt_minus_triv_over_LB": (nt["S"] - triv["S"]) / abs(lb),
        "family_min": fam_min, "L3": l3, "points": rows, "elapsed_s": time.perf_counter() - t0,
    }


FAMILY_ID_SEEDED = ("rgg3", "rr", "er")


def _kkt_entry(w: FloatArray, p: C0Params, tol: float) -> dict[str, Any]:
    r = kkt_box(w, p, tol)
    return {"is_kkt": r["is_kkt"], "strict": r["strict"], "max_violation": r["max_violation"], "tol": r["tol"]}


def kkt_cell(p: C0Params, cfg: dict[str, Any]) -> dict[str, Any]:
    """KKT de T3 y de cada RGG3 (k, semilla) con amplitud 1 y con la optima."""
    tol = cfg["kkt_tol_rel"]
    out: dict[str, Any] = {}
    items: list[tuple[str, tuple[str, int, int]]] = [("t3", ("t3", 0, 0))]
    items += [(f"rgg3_k{k}_s{s}", ("rgg3", k, s)) for k in cfg["rgg_k"] for s in cfg["seeds"]]
    for label, key in items:
        a = _adj(*key)
        t_opt, _ = R.optimal_amplitude_from_profile(_prof(*key), p)
        out[label] = {"t_opt": t_opt, "amp1": _kkt_entry(a, p, tol), "opt": _kkt_entry(t_opt * a, p, tol)}
    return out


def run_cell(args: tuple[dict[str, Any], dict[str, Any]]) -> dict[str, Any]:
    cell, cfg = args
    return eval_cell(cell, cfg)


# ------------------------------------------------------------------ agregacion


def aggregate(results: list[dict[str, Any]], cfg: dict[str, Any]) -> dict[str, Any]:
    n = len(results)
    if n == 0:
        return {"decision": "SIN DATOS"}
    n_l1 = sum(1 for r in results if r["L1"])
    n_rota = sum(1 for r in results if r["L2"] == "rota")
    n_emp = sum(1 for r in results if r["L2"] == "empate")
    l2_global = (
        "ROTA" if n_rota / n >= cfg["l2_global_frac"]
        else "EMPATE" if (n_rota + n_emp) / n >= cfg["l2_global_frac"] else "NO ROTA"
    )
    t3_cells: list[dict[str, Any]] = []
    rgg_cells: list[dict[str, Any]] = []
    for r in results:
        c = r["cell"]
        tag = {"cell_idx": c["cell_idx"], "c_star": c["c_star"], "k_star": c["k_star"], "a": c["a"]}
        t3 = r["L3"]["t3"]
        if t3["amp1"]["is_kkt"] or t3["opt"]["is_kkt"]:
            t3_cells.append({**tag, "amp1": t3["amp1"]["is_kkt"], "amp1_strict": t3["amp1"]["strict"],
                             "opt": t3["opt"]["is_kkt"], "opt_strict": t3["opt"]["strict"]})
        hits = [
            {"label": k, "amp1": v["amp1"]["is_kkt"], "amp1_strict": v["amp1"]["strict"],
             "opt": v["opt"]["is_kkt"], "opt_strict": v["opt"]["strict"]}
            for k, v in r["L3"].items()
            if k.startswith("rgg3") and (v["amp1"]["is_kkt"] or v["opt"]["is_kkt"])
        ]
        if hits:
            rgg_cells.append({**tag, "hits": hits})
    best_counts: dict[str, int] = {}
    for r in results:
        fam = r["best"]["ref"].split("_s")[0] if r["best"]["ref"].startswith("clique") else r["best"]["ref"]
        best_counts[fam] = best_counts.get(fam, 0) + 1
    return {
        "n_cells": n,
        "L1": {"n_clique_binaria_best": n_l1, "frac": n_l1 / n, "n_strict": sum(1 for r in results if r["L1_strict"]),
               "decision": "MUERTE-L1" if n_l1 / n >= cfg["l1_global_frac"] else "SOBREVIVE"},
        "L2": {"n_rota": n_rota, "n_empate": n_emp, "n_no": n - n_rota - n_emp,
               "frac_rota": n_rota / n, "frac_rota_o_empate": (n_rota + n_emp) / n, "decision": l2_global},
        "L3": {"t3_kkt_cells": t3_cells, "rgg3_kkt_cells": rgg_cells,
               "t3_kkt_strict_cells": sum(1 for c in t3_cells if c["amp1_strict"] or c["opt_strict"]),
               "decision": "INDICIO NEGATIVO" if not t3_cells and not rgg_cells else "KKT PRESENTE"},
        "best_overall_counts": best_counts,
    }


def cell_table(results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "cell_idx": r["cell"]["cell_idx"], "c_star": r["cell"]["c_star"], "k_star": r["cell"]["k_star"], "a": r["cell"]["a"],
            "LB": r["LB"], "best": r["best"], "best_trivial": r["best_trivial"], "best_no_trivial": r["best_no_trivial"],
            "L1": r["L1"], "L2": r["L2"], "gap_nt_minus_triv_over_LB": r["gap_nt_minus_triv_over_LB"], "family_min": r["family_min"],
        }
        for r in results
    ]


def write_points(path: Path, results: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["cell_idx", "c_star", "k_star", "a", "ref", "family", "group", "seed", "t", "S", "S_over_LB", "LB"])
        for r in results:
            c = r["cell"]
            for pt in r["points"]:
                wr.writerow([c["cell_idx"], c["c_star"], c["k_star"], c["a"], pt["ref"], pt["family"], pt["group"],
                             "" if pt["seed"] is None else pt["seed"], f"{pt['t']:.12g}", f"{pt['S']:.12g}",
                             f"{pt['S_over_LB']:.12g}", f"{r['LB']:.12g}"])
    os.replace(tmp, path)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0] if __doc__ else "")
    ap.add_argument("--smoke", action="store_true", help="3 celdas, 1 semilla (no valida la decision)")
    ap.add_argument("--allow-dirty", action="store_true", help="solo smoke: permite arbol sucio (queda marcado)")
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)
    out = args.out or (ROOT / "results" / ("c0_landscape_smoke" if args.smoke else "c0_landscape"))
    dirty = tree_dirty()
    if dirty and not (args.smoke and args.allow_dirty):
        print("ERROR: arbol git sucio; el modo full exige commit limpio (--allow-dirty solo con --smoke).", file=sys.stderr)
        return 2
    cfg = build_config(args.smoke)
    chash = canonical_hash(cfg)
    meta: dict[str, Any] = {
        "step": STEP, "mode": "smoke" if args.smoke else "full", "code_commit": head_commit(), "git_dirty": dirty,
        "allow_dirty_used": bool(args.allow_dirty and dirty), "config_hash": chash, "config": cfg,
        "dependency_lock": dependency_lock(),
    }
    out.mkdir(parents=True, exist_ok=True)
    log_path = out / "log.txt"
    log_fh = log_path.open("a", encoding="utf-8")

    def log(s: str) -> None:
        print(s, flush=True)
        log_fh.write(s + "\n")
        log_fh.flush()

    cells = cell_grid(cfg)
    _atomic_write(out / "summary.json", json.dumps(to_jsonable({**meta, "phase": "preregistered"}), indent=2, allow_nan=False) + "\n")
    log(f"[C0-L1..3] {len(cells)} celdas, config_hash {chash[:12]}, commit {meta['code_commit'][:10]}, dirty={dirty}")
    t0 = time.perf_counter()
    payload = [(c, cfg) for c in cells]
    if args.procs <= 1:
        results = [run_cell(x) for x in payload]
    else:
        with mp.get_context("fork").Pool(args.procs) as pool:
            results = []
            for i, r in enumerate(pool.imap_unordered(run_cell, payload, chunksize=1), 1):
                results.append(r)
                log(f"  celda {i}/{len(cells)} (idx {r['cell']['cell_idx']}, {r['elapsed_s']:.1f}s)")
    results.sort(key=lambda r: int(r["cell"]["cell_idx"]))
    agg = aggregate(results, cfg)
    final = {**meta, "phase": "final", "total_seconds": time.perf_counter() - t0, "result": agg, "cells": cell_table(results)}
    _atomic_write(out / "summary.json", json.dumps(to_jsonable(final), indent=2, allow_nan=False) + "\n")
    write_points(out / "points.csv", results)
    log(f"[C0-L1] clique binaria mejor en {agg['L1']['n_clique_binaria_best']}/{agg['n_cells']} -> {agg['L1']['decision']}")
    log(f"[C0-L2] rotas {agg['L2']['n_rota']}, empates {agg['L2']['n_empate']} -> {agg['L2']['decision']}")
    log(f"[C0-L3] T3 KKT en {len(agg['L3']['t3_kkt_cells'])} celdas, RGG3 KKT en {len(agg['L3']['rgg3_kkt_cells'])} -> {agg['L3']['decision']}")
    log_fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
