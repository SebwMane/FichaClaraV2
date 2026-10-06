"""Fase 2, protocolo C0-Theta (docs/OMEGA_FASE2_PRERREGISTRO.md §2): supervivencia (Theta-S) y formacion (Theta-F).

Salida: results/c0_theta/{runs.jsonl, summary.json, log.txt}; W finales en runs/c0_theta/out/.
Uso: python tools/c0_theta.py [--resume] [--summarize-only] [--procs 4] [--dt 0.02] [--smoke]
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

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from c0_dynamics import _atomic_write, _params, read_rows  # noqa: E402
from c0_redteam import build_input  # noqa: E402

from omega.c0 import references as R  # noqa: E402
from omega.c0.functional import action_c0  # noqa: E402
from omega.c0.langevin import background_fraction, langevin_c0  # noqa: E402
from omega.c0.locality import classify_c0, strong_support  # noqa: E402
from omega.experiments.v11.gate import head_commit, tree_dirty  # noqa: E402

MASTER = 20261006
CELLS = (19, 37, 55, 73)
GRID_S = (0.0, 0.1, 0.3, 1.0, 3.0)
GRID_F = (0.3, 1.0, 3.0)
SEEDS = (0, 1, 2)
N_STEPS = 40000
OUT = ROOT / "results" / "c0_theta"
RUNS = ROOT / "runs" / "c0_theta"


def cell_info(idx: int) -> dict[str, Any]:
    cs, ks, av = [0, 1, 2, 4, 8], [4, 6, 8, 12, 16], [0.25, 1.0, 4.0]
    c, rem = divmod(idx, 15)
    k, a = divmod(rem, 3)
    return {"cell_idx": idx, "c_star": cs[c], "k_star": ks[k], "a": av[a]}


def theta_n(cell: dict[str, Any], n: int) -> float:
    p = _params(cell)
    g_bg = p.psi_star - p.a
    return float(g_bg * cell["k_star"] / n)


def build_tasks() -> list[dict[str, Any]]:
    l4 = {(x["cell_idx"], x["init"], x["seed"]): x for x in read_rows(ROOT / "results" / "c0_dynamics" / "runs.jsonl")
          if x.get("kind") == "generic"}
    tasks: list[dict[str, Any]] = []
    for c in CELLS:
        for ti, mult in enumerate(GRID_S):
            for s in SEEDS:
                tasks.append({"id": f"S_c{c}_R_n343_t{ti}_s{s}", "proto": "S", "cell": c, "init": "R", "n": 343, "seed": s,
                              "ti": ti, "mult": mult, "src": f"runs/c0_redteam/out/R6_c{c:02d}_R_n343_s{s}.npz"})
                for init in ("U", "R"):
                    tasks.append({"id": f"S_c{c}_{init}_n216_t{ti}_s{s}", "proto": "S", "cell": c, "init": init, "n": 216,
                                  "seed": s, "ti": ti, "mult": mult, "src": "runs/c0/" + l4[(c, init, s)]["npz"]})
        for ti, mult in enumerate(GRID_F):
            for s in SEEDS:
                tasks.append({"id": f"F_c{c}_R_n216_t{ti}_s{s}", "proto": "F", "cell": c, "init": "R", "n": 216, "seed": s,
                              "ti": ti, "mult": mult, "src": None})
    tasks.sort(key=lambda t: -t["n"])  # las largas primero
    return tasks


def jaccard(a: np.ndarray, b: np.ndarray) -> float:
    iu = np.triu_indices(a.shape[0], 1)
    x, y = a[iu], b[iu]
    union = int(np.sum(x | y))
    return float(np.sum(x & y) / union) if union else 1.0


def run_task(task: dict[str, Any], dt: float, n_steps: int) -> dict[str, Any]:
    t0 = time.perf_counter()
    cell = cell_info(task["cell"])
    p = _params(cell)
    n = int(task["n"])
    ii = "UER".index(task["init"])
    blk = 70 if task["proto"] == "S" else 71
    if task["src"]:
        z = np.load(ROOT / task["src"])
        w0 = np.asarray(z[z.files[0]], dtype=np.float64)
    else:
        w0 = build_input("R", n, float(cell["k_star"]), R.rng_from_key((MASTER, 72, n, task["cell"], ii, task["seed"])))
    assert w0.shape == (n, n)
    th_n = theta_n(cell, n)
    theta = task["mult"] * th_n
    rng = R.rng_from_key((MASTER, blk, n, task["cell"], ii, task["seed"], task["ti"]))
    out = langevin_c0(w0, p, theta, dt=dt, n_steps=n_steps, rng=rng)
    w = out["w"]
    cls = classify_c0(w, R.rng_from_key((MASTER, 73, n, task["cell"], ii, task["seed"], task["ti"], blk)))
    a_fin = strong_support(w)
    row: dict[str, Any] = {k: v for k, v in task.items()}
    row.update(cell)
    lb = p.lower_bound(n)
    s_vals = [x["S"] for x in out["trace"]]
    q = max(1, len(s_vals) // 4)
    row.update({
        "dt": dt, "n_steps": n_steps, "theta_N": th_n, "theta": theta, "cls": cls,
        "f_bg": background_fraction(w), "m_rel": float(a_fin.sum() / 2) / (n * cell["k_star"] / 2.0),
        "S0_over_LB": action_c0(w0, p) / lb, "S_over_LB": out["s_final"] / lb,
        "stationarity": abs(np.mean(s_vals[-q:]) - np.mean(s_vals[-2 * q:-q])) / abs(np.mean(s_vals[-q:])) if len(s_vals) >= 2 * q else None,
        "trace": [{"step": x["step"], "S_over_LB": x["S"] / lb, "f_bg": x["f_bg"]} for x in out["trace"]],
    })
    row["J"] = jaccard(strong_support(w0), a_fin) if task["proto"] == "S" else None
    dest = RUNS / "out" / f"{task['id']}.npz"
    dest.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(dest, w=w)
    row["elapsed_s"] = time.perf_counter() - t0
    return row


def _verdict(sel: list[dict[str, Any]], proto: str) -> str:
    need = math.ceil(2 / 3 * len(SEEDS))
    ok = [r["cls"]["class"] == "DISPERSO_LOCAL" and r["f_bg"] <= 0.5 and (proto == "F" or r["J"] >= 0.5) for r in sel]
    bad = [r["cls"]["class"] != "DISPERSO_LOCAL" or r["f_bg"] > 0.5 for r in sel]
    if sum(ok) >= need:
        return "SOBREVIVE"
    if sum(bad) >= need:
        return "DESTRUIDA"
    return "MIXTA"


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    groups: dict[tuple[str, int, int, str, float], list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        groups[(r["proto"], r["cell"], r["n"], r["init"], r["mult"])].append(r)
    table = []
    for (proto, c, n, init, mult), sel in sorted(groups.items()):
        table.append({"proto": proto, "cell": c, "n": n, "init": init, "mult": mult, "verdict": _verdict(sel, proto),
                      "classes": [r["cls"]["class"] for r in sel], "f_bg": [round(r["f_bg"], 4) for r in sel],
                      "J": [None if r["J"] is None else round(r["J"], 4) for r in sel],
                      "m_rel": [round(r["m_rel"], 3) for r in sel], "S_over_LB": [round(r["S_over_LB"], 4) for r in sel]})
    vmap = {(t["proto"], t["cell"], t["n"], t["init"], t["mult"]): t["verdict"] for t in table}
    control_ok = all(vmap.get(("S", c, n, i, 0.0)) == "SOBREVIVE" for c in CELLS for n, i in ((216, "U"), (216, "R"), (343, "R")))
    s0 = [r for r in rows if r["proto"] == "S" and r["mult"] == 0.0]
    control_s_ok = all(r["S_over_LB"] >= r["S0_over_LB"] * (1 - 0.01) for r in s0)  # S/LB no empeora > 1 %
    bound: dict[str, float] = {}
    for c in CELLS:
        for n, i in ((216, "U"), (216, "R"), (343, "R")):
            tc = 0.0
            for m in GRID_S[1:]:
                if vmap.get(("S", c, n, i, m)) == "SOBREVIVE":
                    tc = m
                else:
                    break
            bound[f"c{c}_{i}_n{n}"] = tc
    idx = {m: k for k, m in enumerate(GRID_S)}
    r_cells = [(bound[f"c{c}_R_n216"], bound[f"c{c}_R_n343"]) for c in CELLS]
    if all(a in (0.1, 0.3, 1.0) and b in (0.1, 0.3, 1.0) and abs(idx[a] - idx[b]) <= 1 for a, b in r_cells):
        glob = "FRAGIL-1/N"
    elif sum(1 for c in CELLS if vmap.get(("S", c, 216, "R", 3.0)) == "SOBREVIVE" and vmap.get(("S", c, 343, "R", 3.0)) == "SOBREVIVE") >= 3:
        glob = "ROBUSTA"
    elif sum(1 for a, b in r_cells if a == 0.0 and b == 0.0) >= 3:
        glob = "INESTABLE"
    else:
        glob = "NO_CONCLUYENTE"
    form = {f"c{c}": ("FORMA" if vmap.get(("F", c, 216, "R", 0.3)) == "SOBREVIVE" else "NO_FORMA") for c in CELLS}
    return {"theta_S_global": glob, "theta_c_over_theta_N": bound, "theta_F": form,
            "control_theta0_class_ok": control_ok, "control_theta0_S_ok": control_s_ok, "table": table}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--summarize-only", action="store_true")
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--dt", type=float, default=0.0025)  # F2-A1
    ap.add_argument("--smoke", action="store_true", help="4 tareas, 400 pasos, salida en results/c0_theta_smoke")
    args = ap.parse_args(argv)
    out = ROOT / "results" / ("c0_theta_smoke" if args.smoke else "c0_theta")
    tasks = build_tasks()
    n_steps = N_STEPS
    if args.smoke:
        tasks = [t for t in tasks if t["seed"] == 0 and t["cell"] == 37 and t["mult"] in (0.0, 1.0) and t["n"] == 216][:4]
        n_steps = 400
    out.mkdir(parents=True, exist_ok=True)
    jsonl = out / "runs.jsonl"
    if not args.summarize_only:
        if tree_dirty() and not (args.resume or args.smoke):
            print("arbol sucio: commitear antes de la corrida oficial", file=sys.stderr)
            return 2
        if jsonl.exists() and not args.resume:
            print("runs.jsonl existe: usar --resume", file=sys.stderr)
            return 2
        done = {r["id"] for r in read_rows(jsonl)}
        todo = [t for t in tasks if t["id"] not in done]
        t0 = time.perf_counter()
        with (out / "log.txt").open("a", encoding="utf-8") as log, jsonl.open("a", encoding="utf-8") as fh:
            log.write(f"commit {head_commit()} dt {args.dt} tareas {len(tasks)} pendientes {len(todo)}\n")
            with mp.get_context("fork").Pool(args.procs) as pool:
                for i, row in enumerate(pool.imap_unordered(_run_star, [(t, args.dt, n_steps) for t in todo]), 1):
                    fh.write(json.dumps(row, default=str) + "\n")
                    fh.flush()
                    msg = f"  {i}/{len(todo)} {row['id']} {row['cls']['class']} f_bg={row['f_bg']:.3f} ({time.perf_counter() - t0:.0f}s)"
                    print(msg, flush=True)
                    log.write(msg + "\n")
                    log.flush()
    rows = read_rows(jsonl)
    res = summarize(rows) if not args.smoke else {"rows": len(rows)}
    summary = {"step": "c0_theta", "code_commit": head_commit(), "n_tasks": len(tasks), "n_rows": len(rows),
               "complete": len({r["id"] for r in rows}) == len(tasks), "result": res}
    _atomic_write(out / "summary.json", json.dumps(summary, indent=2, default=str))
    print(json.dumps({k: v for k, v in res.items() if k != "table"}, indent=1))
    return 0


def _run_star(args: tuple[dict[str, Any], float, int]) -> dict[str, Any]:
    return run_task(*args)


if __name__ == "__main__":
    raise SystemExit(main())
