"""Omega-C0: R6 (tamano N = 125 y 343) y diagnostico D-1 (continuacion a 100 000 pasos), enmienda C0-A5.

Lee results/c0_dynamics/runs.jsonl (resultado oficial de C0-L4). Salida: results/c0_redteam/{summary.json, runs.jsonl,
log.txt}; W finales en runs/c0_redteam/out/ (ignorado por git). Uso: python tools/c0_redteam.py [--resume] [--summarize-only]
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
from collections import Counter, defaultdict  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from c0_dynamics import _atomic_write, _params, _sym_from_upper, read_rows  # noqa: E402

from omega.c0 import references as R  # noqa: E402
from omega.c0.dynamics import evolve_c0  # noqa: E402
from omega.c0.locality import classify_c0  # noqa: E402
from omega.experiments.v11.gate import head_commit, tree_dirty  # noqa: E402
from omega.types import FloatArray  # noqa: E402

MASTER = 20261005
SIZES = (125, 343)
SEEDS = (0, 1, 2)
MAX_STEPS = 20000
D1_TOTAL = 100000
TOL, PATIENCE = 1e-10, 50
SRC = ROOT / "results" / "c0_dynamics"
OUT = ROOT / "results" / "c0_redteam"
RUNS = ROOT / "runs" / "c0_redteam"
SRC_RUNS = ROOT / "runs" / "c0"


def r6_plan(src: list[dict[str, Any]]) -> dict[int, list[str]]:
    """Celdas que superan R3 y R8 -> inicios que votaron LOCAL (>= 2/3 semillas)."""
    votes: dict[tuple[int, str], int] = Counter()
    kmax: dict[int, float] = defaultdict(float)
    for x in src:
        if x["kind"] == "generic" and x["cls"]["class"] == "DISPERSO_LOCAL":
            votes[(x["cell_idx"], x["init"])] += 1
            kmax[x["cell_idx"]] = max(kmax[x["cell_idx"]], float(x["cls"]["kmax_over_kmean"]))
    voted: dict[int, list[str]] = defaultdict(list)
    for (c, i), v in sorted(votes.items()):
        if v >= math.ceil(2 / 3 * len(SEEDS)):
            voted[c].append(i)
    return {c: v for c, v in sorted(voted.items()) if len(v) >= 2 and kmax[c] <= 3.0}


def cell_of(src: list[dict[str, Any]], idx: int) -> dict[str, Any]:
    x = next(r for r in src if r["cell_idx"] == idx)
    return {"cell_idx": idx, "c_star": x["c_star"], "k_star": x["k_star"], "a": x["a"]}


def build_input(init: str, n: int, k: float, rng: np.random.Generator) -> FloatArray:
    """Mismos inicios genericos que C0-L4, intensivos en N."""
    n_pairs = n * (n - 1) // 2
    if init == "U":
        return _sym_from_upper(np.clip(k / (n - 1) * (1.0 + 0.1 * rng.standard_normal(n_pairs)), 0.0, 1.0), n)
    if init == "E":
        return R.erdos_renyi_m(n, int(round(n * k / 2.0)), rng)
    return _sym_from_upper(rng.random(n_pairs) * 2.0 * k / (n - 1), n)


def build_tasks(src: list[dict[str, Any]]) -> list[dict[str, Any]]:
    tasks: list[dict[str, Any]] = []
    for x in src:  # D-1 primero
        if x["kind"] == "generic" and x["cls"]["class"] == "DISPERSO_LOCAL" and x["status"] != "converged":
            tasks.append({"id": f"D1_{x['id']}", "kind": "D1", "src_id": x["id"], "npz": x["npz"],
                          "cell": cell_of(src, x["cell_idx"]), "init": x["init"], "seed": x["seed"], "n": 216,
                          "steps_done": x["steps"]})
    for c, inits in r6_plan(src).items():
        for init in inits:
            for n in SIZES:
                for s in SEEDS:
                    tasks.append({"id": f"R6_c{c:02d}_{init}_n{n}_s{s}", "kind": "R6", "cell": cell_of(src, c),
                                  "init": init, "seed": s, "n": n})
    return tasks


def run_task(task: dict[str, Any]) -> dict[str, Any]:
    t0 = time.perf_counter()
    cell = task["cell"]
    p = _params(cell)
    n = int(task["n"])
    if task["kind"] == "D1":
        z = np.load(SRC_RUNS / task["npz"])
        w0 = np.asarray(z[z.files[0]], dtype=np.float64)
        steps = D1_TOTAL - int(task["steps_done"])
        rng_c = R.rng_from_key((MASTER, 61, cell["cell_idx"], "UER".index(task["init"]), task["seed"]))
    else:
        rng = R.rng_from_key((MASTER, 60, n, cell["cell_idx"], "UER".index(task["init"]), task["seed"]))
        w0 = build_input(task["init"], n, float(cell["k_star"]), rng)
        steps = MAX_STEPS
        rng_c = R.rng_from_key((MASTER, 62, n, cell["cell_idx"], "UER".index(task["init"]), task["seed"]))
    ev = evolve_c0(w0, p, max_steps=steps, tol=TOL, patience=PATIENCE)
    w = ev["w"]
    cls = classify_c0(w, rng_c)
    lb = p.lower_bound(n)
    out = RUNS / "out" / f"{task['id']}.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, w=w)
    return {"id": task["id"], "kind": task["kind"], "src_id": task.get("src_id"), "cell_idx": cell["cell_idx"],
            "c_star": cell["c_star"], "k_star": cell["k_star"], "a": cell["a"], "init": task["init"],
            "seed": task["seed"], "n": n, "status": ev["status"],
            "steps_total": ev["steps"] + int(task.get("steps_done", 0)), "S_over_LB": ev["s_final"] / lb,
            "kkt_residual": ev["kkt_residual"], "cls": cls, "elapsed_s": time.perf_counter() - t0}


def summarize(src: list[dict[str, Any]], rows: list[dict[str, Any]]) -> dict[str, Any]:
    srcmap = {x["id"]: x for x in src}
    d1 = [r for r in rows if r["kind"] == "D1"]
    d1_trans = Counter(f"{srcmap[r['src_id']]['cls']['class']}->{r['cls']['class']}" for r in d1)
    d1_status = Counter(r["status"] for r in d1)
    still = [r for r in d1 if r["cls"]["class"] == "DISPERSO_LOCAL"]
    r6 = [r for r in rows if r["kind"] == "R6"]
    plan = r6_plan(src)
    need = math.ceil(2 / 3 * len(SEEDS))
    cells_out: list[dict[str, Any]] = []
    passed: list[int] = []
    for c, inits in plan.items():
        per: dict[str, dict[int, int]] = {}
        ok_any = False
        for init in inits:
            per[init] = {n: sum(1 for r in r6 if r["cell_idx"] == c and r["init"] == init and r["n"] == n
                                and r["cls"]["class"] == "DISPERSO_LOCAL") for n in SIZES}
            ok_any = ok_any or all(per[init][n] >= need for n in SIZES)
        cells_out.append({"cell_idx": c, "inits": inits, "local_counts": per, "R6_pass": ok_any})
        if ok_any:
            passed.append(c)
    return {
        "D1": {"n": len(d1), "transitions": dict(d1_trans), "status": dict(d1_status),
               "still_local": len(still), "still_local_converged": sum(1 for r in still if r["status"] == "converged")},
        "R6": {"n_cells": len(plan), "n_runs": len(r6), "cells": cells_out, "cells_pass": passed,
               "class_counts_by_n": {str(n): dict(Counter(r["cls"]["class"] for r in r6 if r["n"] == n)) for n in SIZES}},
        "candidato_c0_cells": passed,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--summarize-only", action="store_true")
    ap.add_argument("--procs", type=int, default=4)
    args = ap.parse_args(argv)
    src = read_rows(SRC / "runs.jsonl")
    tasks = build_tasks(src)
    OUT.mkdir(parents=True, exist_ok=True)
    jsonl = OUT / "runs.jsonl"
    if not args.summarize_only:
        if tree_dirty() and not args.resume:
            print("arbol sucio: commitear antes de la corrida oficial", file=sys.stderr)
            return 2
        done = {r["id"] for r in read_rows(jsonl)} if args.resume else set()
        if jsonl.exists() and not args.resume:
            print("runs.jsonl existe: usar --resume", file=sys.stderr)
            return 2
        todo = [t for t in tasks if t["id"] not in done]
        t0 = time.perf_counter()
        with (OUT / "log.txt").open("a", encoding="utf-8") as log, jsonl.open("a", encoding="utf-8") as fh:
            log.write(f"commit {head_commit()} tareas {len(tasks)} pendientes {len(todo)}\n")
            with mp.get_context("fork").Pool(args.procs) as pool:
                for i, row in enumerate(pool.imap_unordered(run_task, todo), 1):
                    fh.write(json.dumps(row) + "\n")
                    fh.flush()
                    if i % 10 == 0 or i == len(todo):
                        msg = f"  {i}/{len(todo)} ({time.perf_counter() - t0:.0f}s)"
                        print(msg, flush=True)
                        log.write(msg + "\n")
                        log.flush()
    rows = read_rows(jsonl)
    summary = {"step": "c0_redteam", "code_commit": head_commit(), "n_tasks": len(tasks), "n_rows": len(rows),
               "complete": len({r["id"] for r in rows}) == len(tasks), "result": summarize(src, rows)}
    _atomic_write(OUT / "summary.json", json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary["result"]["D1"]), "\nR6 pasa:", summary["result"]["R6"]["cells_pass"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
