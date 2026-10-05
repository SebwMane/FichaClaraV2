"""Omega-C0, prueba F1 (prerregistro §7): subfamilia R-3D con N = 512 y 729, certificado Omega-1.1 sin cambios.

Salida: results/c0_f1/{runs.jsonl, summary.json, log.txt}; W finales en runs/c0_f1/out/.
Uso: python tools/c0_f1.py [--resume] [--summarize-only] [--procs 4]
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

from c0_battery import describe, size_exponent  # noqa: E402
from c0_dynamics import _atomic_write, _params, read_rows  # noqa: E402
from c0_redteam import build_input  # noqa: E402
from l3b_stability import certificate_report  # noqa: E402

from omega.c0 import references as R  # noqa: E402
from omega.c0.dynamics import evolve_c0  # noqa: E402
from omega.c0.locality import classify_c0, degree_preserving_rewire, giant_component, mean_hop_connected, strong_support  # noqa: E402
from omega.config.seeds import SeedKey  # noqa: E402
from omega.experiments.v11.gate import head_commit, tree_dirty  # noqa: E402

MASTER = 20261005
CELLS_512 = (18, 19, 20, 36, 37, 38, 55, 56, 73, 74)
CELLS_729 = (19, 37, 55, 73)
SEEDS = (0, 1, 2)
MAX_STEPS = 40000
TOL, PATIENCE = 1e-10, 50
OUT = ROOT / "results" / "c0_f1"
RUNS = ROOT / "runs" / "c0_f1"


def cell_info(idx: int) -> dict[str, Any]:
    cs, ks, av = [0, 1, 2, 4, 8], [4, 6, 8, 12, 16], [0.25, 1.0, 4.0]
    c, rem = divmod(idx, len(ks) * len(av))
    k, a = divmod(rem, len(av))
    return {"cell_idx": idx, "c_star": cs[c], "k_star": ks[k], "a": av[a]}


def build_tasks() -> list[dict[str, Any]]:
    tasks: list[dict[str, Any]] = []
    for n, cells in ((729, CELLS_729), (512, CELLS_512)):  # las largas primero
        for c in cells:
            for s in SEEDS:
                tasks.append({"id": f"F1_c{c:02d}_R_n{n}_s{s}", "kind": "c0", "cell": cell_info(c), "n": n, "seed": s})
    for n, side in ((512, 8), (729, 9)):
        tasks.append({"id": f"REF_T3_n{n}", "kind": "ref", "ref": "T3", "side": side, "n": n, "seed": 0})
        for s in SEEDS:
            tasks.append({"id": f"REF_RGG3_n{n}_s{s}", "kind": "ref", "ref": "RGG3", "n": n, "seed": s})
    return tasks


def null_hop(w: np.ndarray, rng: np.random.Generator) -> float:
    a = strong_support(w)
    g = giant_component(a)
    ag = np.asarray(a[np.ix_(g, g)], dtype=np.bool_)
    an = degree_preserving_rewire(ag, 10 * int(ag.sum() // 2), rng)
    gn = giant_component(an)
    return mean_hop_connected(np.asarray(an[np.ix_(gn, gn)], dtype=np.bool_))


def run_task(task: dict[str, Any]) -> dict[str, Any]:
    t0 = time.perf_counter()
    n = int(task["n"])
    key = SeedKey(MASTER, (81, n, task.get("cell", {}).get("cell_idx", 99), task["seed"]))
    row: dict[str, Any] = {k: v for k, v in task.items() if k != "cell"}
    if task["kind"] == "ref":
        if task["ref"] == "T3":
            w = R.torus_lattice_3d(int(task["side"]))
        else:
            w = R.rgg3_torus_binary(n, 12.0, R.rng_from_key((MASTER, 82, n, task["seed"])))
        status = "CONVERGED"
    else:
        cell = task["cell"]
        row.update(cell)
        p = _params(cell)
        rng = R.rng_from_key((MASTER, 60, n, cell["cell_idx"], "UER".index("R"), task["seed"]))
        w0 = build_input("R", n, float(cell["k_star"]), rng)
        ev = evolve_c0(w0, p, max_steps=MAX_STEPS, tol=TOL, patience=PATIENCE)
        w = ev["w"]
        status = "CONVERGED" if ev["status"] == "converged" else "MAX_STEPS"
        row.update({"status": ev["status"], "steps": ev["steps"], "S_over_LB": ev["s_final"] / p.lower_bound(n),
                    "kkt_residual": ev["kkt_residual"]})
        row["cls"] = classify_c0(w, R.rng_from_key((MASTER, 62, n, cell["cell_idx"], 2, task["seed"])))
        out = RUNS / "out" / f"{task['id']}.npz"
        out.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(out, w=w)
    row["battery"] = describe(w, task["seed"])
    row["null_mean_hop"] = null_hop(w, R.rng_from_key((MASTER, 90, n, task["seed"], len(task["id"]))))
    row["certificate"] = certificate_report(w, status, key)
    row["elapsed_s"] = time.perf_counter() - t0
    return row


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    prev216 = [r for r in read_rows(ROOT / "results" / "c0_battery" / "finals.jsonl") if r["init"] == "R"]
    nullc = {(t["cell_idx"], t["init"]): t for t in json.loads((ROOT / "results" / "c0_battery" / "null_control.json").read_text())["table"]}
    refs = [r for r in rows if r["kind"] == "ref"]
    ref_codes: dict[int, set[str]] = defaultdict(set)
    ref_table = []
    for r in refs:
        codes = r["certificate"].get("codes") or []
        if r["ref"] == "RGG3":
            ref_codes[r["n"]].update(codes)
        ref_table.append({"id": r["id"], "codes": codes, "primary": r["certificate"].get("primary"),
                          "D_eff": r["battery"]["D_eff"], "D_s": r["battery"]["D_s"], "mean_hop": r["battery"]["mean_hop"]})
    c0 = [r for r in rows if r["kind"] == "c0"]
    cells = []
    need = math.ceil(2 / 3 * len(SEEDS))
    for c in CELLS_512:
        nmax = 729 if c in CELLS_729 else 512
        sel = [r for r in c0 if r["cell_idx"] == c]
        top = [r for r in sel if r["n"] == nmax]
        n_local = sum(1 for r in top if r["cls"]["class"] == "DISPERSO_LOCAL")
        n_cert = sum(1 for r in top if r["certificate"].get("codes") == [])
        hop: dict[int, list[float]] = defaultdict(list)
        nhop: dict[int, list[float]] = defaultdict(list)
        for r in prev216:
            if r["cell_idx"] == c and r["n"] in (216, 343):
                hop[r["n"]].append(r["mean_hop"])
        for r in sel:
            hop[r["n"]].append(r["battery"]["mean_hop"])
            nhop[r["n"]].append(r["null_mean_hop"])
        dl = size_exponent({n: v for n, v in hop.items() if n >= 216})
        dl_null_hi = size_exponent(nhop) if len(nhop) >= 3 else None
        dl_null_prev = nullc.get((c, "R"), {}).get("D_L_null")
        extra_codes = sorted({code for r in top for code in (r["certificate"].get("codes") or [])
                              if code != "Ω-F10" and code not in ref_codes[nmax]})
        sep = dl is not None and dl_null_prev is not None and dl < dl_null_prev - 0.5
        in_band = dl is not None and 2.5 <= dl <= 3.5
        if n_local >= need and n_cert >= need:
            verdict = "F1-POSITIVO"
        elif n_local >= need and in_band and sep and not extra_codes:
            verdict = "F1-INDETERMINADO"
        else:
            verdict = "F1-NEGATIVO"
        cells.append({"cell_idx": c, **cell_info(c), "n_max": nmax, "local_at_nmax": n_local, "cert_pass_at_nmax": n_cert,
                      "D_L": dl, "D_L_null_prev": dl_null_prev, "D_L_null_hi": dl_null_hi,
                      "mean_hop_by_n": {n: float(np.median(v)) for n, v in sorted(hop.items())},
                      "codes_at_nmax": [r["certificate"].get("codes") for r in top],
                      "extra_codes_vs_rgg3": extra_codes, "status_at_nmax": [r["status"] for r in top],
                      "D_s_at_nmax": [r["battery"]["D_s"] for r in top], "D_eff_at_nmax": [r["battery"]["D_eff"] for r in top],
                      "verdict": verdict})
    order = {"F1-POSITIVO": 2, "F1-INDETERMINADO": 1, "F1-NEGATIVO": 0}
    best = max((c["verdict"] for c in cells), key=lambda v: order[v]) if cells else "F1-NEGATIVO"
    return {"global": best, "verdicts": dict(Counter(c["verdict"] for c in cells)), "cells": cells,
            "references": ref_table, "rgg3_codes_by_n": {n: sorted(v) for n, v in ref_codes.items()}}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--summarize-only", action="store_true")
    ap.add_argument("--procs", type=int, default=4)
    args = ap.parse_args(argv)
    tasks = build_tasks()
    OUT.mkdir(parents=True, exist_ok=True)
    jsonl = OUT / "runs.jsonl"
    if not args.summarize_only:
        if tree_dirty() and not args.resume:
            print("arbol sucio: commitear antes de la corrida oficial", file=sys.stderr)
            return 2
        if jsonl.exists() and not args.resume:
            print("runs.jsonl existe: usar --resume", file=sys.stderr)
            return 2
        done = {r["id"] for r in read_rows(jsonl)}
        todo = [t for t in tasks if t["id"] not in done]
        t0 = time.perf_counter()
        with (OUT / "log.txt").open("a", encoding="utf-8") as log, jsonl.open("a", encoding="utf-8") as fh:
            log.write(f"commit {head_commit()} tareas {len(tasks)} pendientes {len(todo)}\n")
            with mp.get_context("fork").Pool(args.procs) as pool:
                for i, row in enumerate(pool.imap_unordered(run_task, todo), 1):
                    fh.write(json.dumps(row, default=str) + "\n")
                    fh.flush()
                    msg = f"  {i}/{len(todo)} {row['id']} ({time.perf_counter() - t0:.0f}s)"
                    print(msg, flush=True)
                    log.write(msg + "\n")
                    log.flush()
    rows = read_rows(jsonl)
    res = summarize(rows)
    summary = {"step": "c0_f1", "code_commit": head_commit(), "n_tasks": len(tasks), "n_rows": len(rows),
               "complete": len({r["id"] for r in rows}) == len(tasks), "result": res}
    _atomic_write(OUT / "summary.json", json.dumps(summary, indent=2, default=str))
    print("F1 global:", res["global"], res["verdicts"])
    for c in res["cells"]:
        print(c["cell_idx"], c["verdict"], "D_L", c["D_L"], "local", c["local_at_nmax"], "cert", c["cert_pass_at_nmax"], c["codes_at_nmax"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
