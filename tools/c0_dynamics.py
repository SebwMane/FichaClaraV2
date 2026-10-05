"""C0-L4 (prerregistro Omega-C0 §4): dinamica de S_C0 desde inicios genericos; clasificacion del estado final, N = 216.

Antes de cualquier dinamica se corre `validate_locality` (si falla: codigo de salida 3, hace falta enmienda).
Tareas: 75 celdas x inicios {U,E,R} x semillas {0,1,2} (genericos, votan) + 75 celdas x {T3, RGG3 k12} con ruido iid 1e-2
(semilla 0, solo descriptivos). Decision: MUERTE-A (todo VACIO/DENSO_TRIVIAL), MUERTE-B (hay dispersos/fragmentados pero ninguna
celda con DISPERSO_LOCAL en >= 2/3 semillas del mismo inicio), CONTINUA (alguna celda lo logra; solo se listan candidatas).
Red-team R1-R4, R7, R8 sobre toda la malla (R5 y R6 no se ejecutan aqui). Salida: <out>/summary.json, runs.jsonl, log.txt,
locality_validation.json, r7.json; W finales en <runs-dir>/out/<id>.npz (runs/ esta en .gitignore).
Modo full exige arbol git limpio (--allow-dirty solo con --smoke).
"""

from __future__ import annotations

import os

# Un hilo BLAS por proceso: el paralelismo es por corrida (4 procesos).
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import hashlib
import json
import math
import multiprocessing as mp
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from omega.c0 import references as R  # noqa: E402
from omega.c0.dynamics import evolve_c0  # noqa: E402
from omega.c0.functional import C0Params, params_from_targets  # noqa: E402
from omega.c0.locality import classify_c0, giant_component, mean_hop_connected, strong_support, validate_locality  # noqa: E402
from omega.experiments.v11.gate import head_commit, tree_dirty  # noqa: E402
from omega.io.provenance import dependency_lock  # noqa: E402
from omega.phases.scan import to_jsonable  # noqa: E402
from omega.types import FloatArray  # noqa: E402

N = 216
STEP = "c0_dynamics"
MASTER = 20261005
GENERIC = ("U", "E", "R")
GEO = ("T3", "RGG3k12")
CLASSES = ("VACIO", "DENSO_TRIVIAL", "FRAGMENTADO", "DISPERSO_LOCAL", "DISPERSO_NO_LOCAL")
CLIQUE_STRUCT = ("clique_única", "multi_clique", "cliques_solapadas")

FULL: dict[str, Any] = {
    "n": N, "c_stars": [0, 1, 2, 4, 8], "k_stars": [4, 6, 8, 12, 16], "a_values": [0.25, 1.0, 4.0],
    "seeds": [0, 1, 2], "max_steps": 20000, "tol": 1e-10, "patience": 50, "geo_noise": 1e-2,
    "vote_frac": 2 / 3, "r7_n": 10, "r7_seed": 7, "r3_max": 3.0, "master_entropy": MASTER,
}


def canonical_hash(cfg: dict[str, Any]) -> str:
    text = json.dumps(to_jsonable(cfg), sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_config(smoke: bool) -> dict[str, Any]:
    cfg = dict(FULL)
    cfg["smoke"] = smoke
    if smoke:
        cfg["seeds"] = [0]
        cfg["max_steps"] = 2000
        cfg["smoke_cells"] = [[0, 6, 1.0], [2, 8, 1.0]]
        cfg["r7_n"] = 2
    return cfg


def cell_grid(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    cells = [
        {"cell_idx": i, "c_star": c, "k_star": k, "a": a}
        for i, (c, k, a) in enumerate((c, k, a) for c in cfg["c_stars"] for k in cfg["k_stars"] for a in cfg["a_values"])
    ]
    if cfg["smoke"]:
        pick = {tuple(x) for x in cfg["smoke_cells"]}
        cells = [c for c in cells if (c["c_star"], c["k_star"], c["a"]) in pick]
    return cells


def build_tasks(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    tasks: list[dict[str, Any]] = []
    for cell in cell_grid(cfg):
        for ii, init in enumerate(GENERIC):
            for s in cfg["seeds"]:
                tasks.append({"id": f"c{cell['cell_idx']:02d}_{init}_s{s}", "cell": cell, "init": init,
                              "init_idx": ii, "seed": s, "kind": "generic"})
        for gi, init in enumerate(GEO):
            tasks.append({"id": f"c{cell['cell_idx']:02d}_{init}_noise", "cell": cell, "init": init,
                          "init_idx": 10 + gi, "seed": 0, "kind": "geo"})
    return tasks


# ------------------------------------------------------------------ inicios


def _sym_from_upper(v: FloatArray, n: int) -> FloatArray:
    w = np.zeros((n, n), dtype=np.float64)
    iu = np.triu_indices(n, k=1)
    w[iu] = v
    return np.asarray(w + w.T, dtype=np.float64)


def task_rng(task: dict[str, Any], stream: int = 0) -> np.random.Generator:
    """PCG64 determinista por (celda, inicio, semilla, flujo)."""
    return R.rng_from_key((MASTER, 10, task["cell"]["cell_idx"], task["init_idx"], task["seed"], stream))


def build_input(task: dict[str, Any], cfg: dict[str, Any]) -> tuple[FloatArray, FloatArray | None]:
    """(W0, adyacencia del input geometrico o None)."""
    k = float(task["cell"]["k_star"])
    rng = task_rng(task)
    n_pairs = N * (N - 1) // 2
    init = task["init"]
    if init == "U":
        v = np.clip(k / (N - 1) * (1.0 + 0.1 * rng.standard_normal(n_pairs)), 0.0, 1.0)
        return _sym_from_upper(v, N), None
    if init == "E":
        return R.erdos_renyi_m(N, int(round(N * k / 2.0)), rng), None
    if init == "R":
        return _sym_from_upper(rng.random(n_pairs) * 2.0 * k / (N - 1), N), None
    adj = R.torus_lattice_3d(6) if init == "T3" else R.rgg3_torus_binary(N, 12.0, rng)
    noise = _sym_from_upper(rng.standard_normal(n_pairs), N)
    return np.asarray(np.clip(adj + cfg["geo_noise"] * noise, 0.0, 1.0), dtype=np.float64), adj


def giant_hop(a: FloatArray) -> float:
    idx = giant_component(a)
    if idx.size < 2:
        return 0.0
    return mean_hop_connected(np.asarray(a[np.ix_(idx, idx)] > 0, dtype=np.bool_))


def geo_metrics(w: FloatArray, adj: FloatArray) -> dict[str, float]:
    """J (Jaccard aristas input vs soporte fuerte final) y H (salto medio final / input), como L-A5."""
    if float(np.max(w)) <= 1e-6:
        return {"J": 0.0, "H": 0.0}
    iu = np.triu_indices(N, k=1)
    s = strong_support(w)[iu]
    e = adj[iu] > 0
    union = float(np.sum(s | e))
    h0 = giant_hop(adj)
    return {"J": float(np.sum(s & e)) / union if union > 0 else 0.0,
            "H": giant_hop(strong_support(w).astype(np.float64)) / h0 if h0 > 0 else 0.0}


# ------------------------------------------------------------------ una corrida


def _params(cell: dict[str, Any]) -> C0Params:
    return params_from_targets(cell["c_star"], cell["k_star"], cell["a"])


def run_task(args: tuple[dict[str, Any], dict[str, Any]]) -> dict[str, Any]:
    task, run_cfg = args
    t0 = time.perf_counter()
    cfg = run_cfg["cfg"]
    cell = task["cell"]
    p = _params(cell)
    w0, adj = build_input(task, cfg)
    ev = evolve_c0(w0, p, max_steps=cfg["max_steps"], tol=cfg["tol"], patience=cfg["patience"])
    w = ev["w"]
    cls = classify_c0(w, task_rng(task, 1))
    lb = p.lower_bound(N)
    row: dict[str, Any] = {
        "id": task["id"], "kind": task["kind"], "cell_idx": cell["cell_idx"], "c_star": cell["c_star"],
        "k_star": cell["k_star"], "a": cell["a"], "init": task["init"], "seed": task["seed"],
        "lam": p.lam, "kappa": p.kappa, "psi_star": p.psi_star, "LB": lb,
        "status": ev["status"], "steps": ev["steps"], "rejections": ev["rejections"], "S": ev["s_final"],
        "S_over_LB": ev["s_final"] / lb, "kkt_residual": ev["kkt_residual"], "dt_final": ev["dt_final"],
        "w_max": float(w.max()), "w_mean_offdiag": float(w.sum() / (N * (N - 1))), "cls": cls,
        "config_hash": run_cfg["config_hash"],
    }
    if adj is not None:
        row.update(geo_metrics(w, adj))
    npz = Path(run_cfg["runs_dir"]) / "out" / f"{task['id']}.npz"
    npz.parent.mkdir(parents=True, exist_ok=True)
    tmp = npz.with_name(npz.name + ".tmp.npz")
    np.savez_compressed(tmp, w=w)
    os.replace(tmp, npz)
    row["npz"] = f"out/{npz.name}"
    row["npz_sha256"] = file_sha256(npz)
    row["elapsed_s"] = time.perf_counter() - t0
    return row


def r7_task(args: tuple[dict[str, Any], dict[str, Any]]) -> dict[str, Any]:
    """Reevolucion con el estado inicial reetiquetado: max|P W_final P^T - W_final_perm|."""
    task, run_cfg = args
    cfg = run_cfg["cfg"]
    p = _params(task["cell"])
    w0, _ = build_input(task, cfg)
    perm = R.rng_from_key((MASTER, 77, task["cell"]["cell_idx"], task["init_idx"], task["seed"])).permutation(N)
    kw = {"max_steps": cfg["max_steps"], "tol": cfg["tol"], "patience": cfg["patience"]}
    # C0-A1: (i) codigo, horizonte corto 60 pasos, max|dW| <= 1e-8
    s1 = evolve_c0(w0, p, max_steps=60, tol=0.0, patience=10**9)
    s2 = evolve_c0(w0[np.ix_(perm, perm)], p, max_steps=60, tol=0.0, patience=10**9)
    d_short = float(np.max(np.abs(s1["w"][np.ix_(perm, perm)] - s2["w"])))
    # (ii) resultado: misma clase y |dS|/|S| <= 1e-6 en finales convergidos; max|dW| solo se reporta
    r1 = evolve_c0(w0, p, **kw)
    r2 = evolve_c0(w0[np.ix_(perm, perm)], p, **kw)
    d = float(np.max(np.abs(r1["w"][np.ix_(perm, perm)] - r2["w"])))
    c1 = classify_c0(r1["w"], R.rng_from_key((MASTER, 78, task["cell"]["cell_idx"])))["class"]
    c2 = classify_c0(r2["w"], R.rng_from_key((MASTER, 78, task["cell"]["cell_idx"])))["class"]
    ds = abs(r1["s_final"] - r2["s_final"]) / max(abs(r1["s_final"]), 1e-300)
    both_conv = r1["status"] == "converged" and r2["status"] == "converged"
    ok_result = (c1 == c2 and ds <= 1e-6) if both_conv else None
    return {"id": task["id"], "max_abs_diff_short": d_short, "ok_code": bool(d_short <= 1e-8),
            "max_abs_diff": d, "rel_dS": ds, "class": c1, "class_perm": c2, "both_converged": both_conv,
            "ok_result": ok_result, "status": r1["status"], "status_perm": r2["status"],
            "steps": r1["steps"], "steps_perm": r2["steps"],
            "ok": bool(d_short <= 1e-8) and ok_result is not False}


# ------------------------------------------------------------------ agregacion


def read_rows(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return rows


def _dist(v: list[float]) -> dict[str, Any]:
    if not v:
        return {"n": 0}
    a = np.asarray(v, dtype=np.float64)
    q = np.quantile(a, [0.0, 0.1, 0.5, 0.9, 1.0])
    return {"n": int(a.size), "min": float(q[0]), "p10": float(q[1]), "median": float(q[2]), "p90": float(q[3]), "max": float(q[4])}


def summarize(rows: list[dict[str, Any]], cfg: dict[str, Any], n_expected: int) -> dict[str, Any]:
    gen = [r for r in rows if r["kind"] == "generic"]
    geo = [r for r in rows if r["kind"] == "geo"]
    out: dict[str, Any] = {"n_rows": len(rows), "n_expected": n_expected, "n_generic": len(gen), "n_geo": len(geo)}
    if not gen:
        out["decision"] = "SIN DATOS"
        return out
    out["status_counts"] = dict(Counter(r["status"] for r in rows))
    out["class_counts_generic"] = {c: sum(1 for r in gen if r["cls"]["class"] == c) for c in CLASSES}
    by_init = {i: {c: sum(1 for r in gen if r["init"] == i and r["cls"]["class"] == c) for c in CLASSES} for i in GENERIC}
    by_cell: dict[int, dict[str, Any]] = {}
    for r in gen:
        d = by_cell.setdefault(r["cell_idx"], {"c_star": r["c_star"], "k_star": r["k_star"], "a": r["a"], "classes": Counter(),
                                               "by_init": {i: [] for i in GENERIC}})
        d["classes"][r["cls"]["class"]] += 1
        d["by_init"][r["init"]].append(r["cls"]["class"])
    cell_rows = []
    candidates = []
    for ci in sorted(by_cell):
        d = by_cell[ci]
        votes = {}
        for i in GENERIC:
            lst = d["by_init"][i]
            need = math.ceil(cfg["vote_frac"] * len(lst) - 1e-9) if lst else 1
            votes[i] = bool(lst) and sum(1 for c in lst if c == "DISPERSO_LOCAL") >= need
        cell_rows.append({"cell_idx": ci, "c_star": d["c_star"], "k_star": d["k_star"], "a": d["a"],
                          "classes": dict(d["classes"]), "by_init": d["by_init"], "vote_local": votes})
        if any(votes.values()):
            kr = [r["cls"]["kmax_over_kmean"] for r in gen if r["cell_idx"] == ci and r["cls"]["class"] == "DISPERSO_LOCAL"
                  and r["cls"]["kmax_over_kmean"] is not None]
            candidates.append({"cell_idx": ci, "c_star": d["c_star"], "k_star": d["k_star"], "a": d["a"],
                               "inits_voting": [i for i, v in votes.items() if v], "R8_ge2_inits": sum(votes.values()) >= 2,
                               "R3_median_kmax_over_kmean": float(np.median(kr)) if kr else None,
                               "R3_ok": bool(kr) and float(np.max(kr)) <= cfg["r3_max"]})
    if all(r["cls"]["class"] in ("VACIO", "DENSO_TRIVIAL") for r in gen):
        decision = "MUERTE-A"
    elif not candidates:
        decision = "MUERTE-B"
    else:
        decision = "CONTINUA"
    nonempty = [r for r in gen if r["cls"]["class"] != "VACIO"]
    dispersed = [r for r in gen if r["cls"]["class"] in ("DISPERSO_LOCAL", "DISPERSO_NO_LOCAL")]
    out.update({
        "decision": decision, "by_init": by_init, "cells": cell_rows, "candidate_cells": candidates,
        "S_over_LB_generic": _dist([r["S_over_LB"] for r in gen]),
        "R1_cells_with_vacio": sorted({r["cell_idx"] for r in gen if r["cls"]["class"] == "VACIO"}),
        "R2_cells_with_clique_classes": sorted({r["cell_idx"] for r in gen if r["cls"]["structure"] in CLIQUE_STRUCT}),
        "R3_kmax_over_kmean_nonempty": _dist([r["cls"]["kmax_over_kmean"] for r in nonempty if r["cls"]["kmax_over_kmean"] is not None]),
        "R3_n_above_limit": sum(1 for r in nonempty if (r["cls"]["kmax_over_kmean"] or 0.0) > cfg["r3_max"]),
        "R4_h_null_dispersed": _dist([r["cls"]["h_null"] for r in dispersed if r["cls"]["h_null"] is not None]),
        "R4_short_cycles_dispersed": _dist([r["cls"]["short_cycles"] for r in dispersed if r["cls"]["short_cycles"] is not None]),
        "R8_class_by_init": by_init,
        "geo_descriptive": [
            {"id": r["id"], "cell_idx": r["cell_idx"], "init": r["init"], "class": r["cls"]["class"], "J": r.get("J"), "H": r.get("H"),
             "h_null": r["cls"]["h_null"], "short_cycles": r["cls"]["short_cycles"]} for r in geo
        ],
        "geo_class_counts": {i: dict(Counter(r["cls"]["class"] for r in geo if r["init"] == i)) for i in GEO},
        "total_run_seconds": sum(float(r.get("elapsed_s") or 0.0) for r in rows),
    })
    return out


# ------------------------------------------------------------------ ejecucion


def execute(tasks: list[dict[str, Any]], run_cfg: dict[str, Any], jsonl: Path, procs: int, log: Any) -> None:
    jsonl.parent.mkdir(parents=True, exist_ok=True)
    total = len(tasks)
    t0 = time.perf_counter()
    done = 0
    with jsonl.open("a", encoding="utf-8") as fh:
        payload = [(t, run_cfg) for t in tasks]
        pool = None
        if procs <= 1:
            it: Any = map(run_task, payload)
        else:
            pool = mp.get_context("fork").Pool(procs)
            it = pool.imap_unordered(run_task, payload, chunksize=1)
        try:
            for row in it:
                fh.write(json.dumps(to_jsonable(row), allow_nan=False, sort_keys=True) + "\n")
                fh.flush()
                os.fsync(fh.fileno())
                done += 1
                if done % 10 == 0 or done == total:
                    log(f"  {done}/{total} corridas ({time.perf_counter() - t0:.0f}s)")
        finally:
            if pool is not None:
                pool.close()
                pool.join()


def run_r7(tasks: list[dict[str, Any]], run_cfg: dict[str, Any], procs: int) -> dict[str, Any]:
    cfg = run_cfg["cfg"]
    generic = [t for t in tasks if t["kind"] == "generic"]
    rng = np.random.Generator(np.random.PCG64(cfg["r7_seed"]))
    pick = rng.choice(len(generic), size=min(cfg["r7_n"], len(generic)), replace=False)
    chosen = [generic[int(i)] for i in sorted(pick)]
    payload = [(t, run_cfg) for t in chosen]
    if procs <= 1:
        res = [r7_task(x) for x in payload]
    else:
        with mp.get_context("fork").Pool(procs) as pool:
            res = pool.map(r7_task, payload)
    return {"n": len(res), "n_ok": sum(1 for r in res if r["ok"]), "max_abs_diff": max(r["max_abs_diff"] for r in res),
            "threshold": 1e-8, "runs": res}


def time_probe(cfg: dict[str, Any], procs: int, log: Any) -> None:
    """Mide una corrida (inicio U, semilla 0) en celdas dificiles y estima el total (fuera del registro)."""
    cells = {(c["c_star"], c["k_star"], c["a"]): c for c in cell_grid(build_config(False))}
    n_tasks = len(build_tasks(build_config(False)))
    times = []
    tmp = ROOT / "runs" / "c0_probe"
    for key in ((8, 16, 4.0), (4, 12, 0.25), (0, 6, 1.0)):
        task = {"id": "probe", "cell": cells[key], "init": "U", "init_idx": 0, "seed": 0, "kind": "generic"}
        row = run_task((task, {"cfg": cfg, "config_hash": "probe", "runs_dir": str(tmp)}))
        times.append(row["elapsed_s"])
        log(f"[sonda] celda c*={key[0]} k*={key[1]} a={key[2]}: {row['elapsed_s']:.1f}s, {row['steps']} pasos, {row['status']}, "
            f"clase={row['cls']['class']}, S/LB={row['S_over_LB']:.4f}")
    mean, worst = float(np.mean(times)), float(np.max(times))
    log(f"[sonda] total estimado: {n_tasks} corridas x {mean:.0f}s / {procs} procs = {n_tasks * mean / procs / 3600:.2f} h "
        f"(cota alta con la peor: {n_tasks * worst / procs / 3600:.2f} h)")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0] if __doc__ else "")
    ap.add_argument("--smoke", action="store_true", help="2 celdas x 3 inicios x 1 semilla, max_steps 2000")
    ap.add_argument("--allow-dirty", action="store_true", help="solo smoke: permite arbol sucio (queda marcado)")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--summarize-only", action="store_true")
    ap.add_argument("--time-probe", action="store_true")
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--runs-dir", type=Path, default=None)
    args = ap.parse_args(argv)
    log0 = lambda s: print(s, flush=True)  # noqa: E731
    if args.time_probe:
        time_probe(build_config(False), args.procs, log0)
        return 0
    out = args.out or (ROOT / "results" / ("c0_dynamics_smoke" if args.smoke else "c0_dynamics"))
    runs_dir = args.runs_dir or (ROOT / "runs" / ("c0_smoke" if args.smoke else "c0"))
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
    log_fh = (out / "log.txt").open("a", encoding="utf-8")

    def log(s: str) -> None:
        print(s, flush=True)
        log_fh.write(s + "\n")
        log_fh.flush()

    tasks = build_tasks(cfg)
    jsonl = out / "runs.jsonl"
    summary_path = out / "summary.json"
    run_cfg = {"cfg": cfg, "config_hash": chash, "runs_dir": str(runs_dir)}

    if args.summarize_only:
        rows = [r for r in read_rows(jsonl) if r.get("config_hash") == chash]
        r7p = out / "r7.json"
        r7 = json.loads(r7p.read_text(encoding="utf-8")) if r7p.exists() else None
        final = {**meta, "phase": "final", "complete": len(rows) == len(tasks), "result": summarize(rows, cfg, len(tasks)), "R7": r7}
        _atomic_write(summary_path, json.dumps(to_jsonable(final), indent=2, allow_nan=False) + "\n")
        log(f"[C0-L4] re-agregado: {len(rows)}/{len(tasks)} filas; decision {final['result'].get('decision')}")
        return 0

    existing = read_rows(jsonl)
    if existing and not args.resume:
        print(f"ERROR: {jsonl} ya existe con {len(existing)} filas; usa --resume o borra el directorio.", file=sys.stderr)
        return 2
    if any(r.get("config_hash") != chash for r in existing):
        print("ERROR: runs.jsonl contiene filas con otro config_hash; no se reanuda.", file=sys.stderr)
        return 2

    val = validate_locality(0)
    _atomic_write(out / "locality_validation.json", json.dumps(to_jsonable(val), indent=2, allow_nan=False) + "\n")
    for name, r in val["graphs"].items():
        log(f"  validacion {name}: {r['class']} h_null={r['h_null']} ciclos={r['short_cycles']} ok={r['ok']}")
    if not val["ok"]:
        log("[C0-L4] VALIDACION DE LOCALIDAD FALLIDA: hace falta una enmienda antes de la dinamica; no se corre nada.")
        _atomic_write(summary_path, json.dumps(to_jsonable({**meta, "phase": "final", "aborted": True,
                      "result": {"decision": "LOCALITY_VALIDATION_FAILED", "validation": val}}), indent=2, allow_nan=False) + "\n")
        return 3
    _atomic_write(summary_path, json.dumps(to_jsonable({**meta, "phase": "preregistered", "complete": False}), indent=2, allow_nan=False) + "\n")
    log(f"[C0-L4] preregistrado (config_hash {chash[:12]}, commit {meta['code_commit'][:10]}, dirty={dirty}, {len(tasks)} corridas)")

    done_ids = {r["id"] for r in existing if (runs_dir / str(r.get("npz", "?"))).exists()}
    pending = [t for t in tasks if t["id"] not in done_ids]
    log(f"[C0-L4] {len(done_ids)} ya hechas, {len(pending)} pendientes")
    t0 = time.perf_counter()
    execute(pending, run_cfg, jsonl, args.procs, log)
    log(f"[C0-L4] corridas listas en {time.perf_counter() - t0:.0f}s; R7")
    r7 = run_r7(tasks, run_cfg, args.procs)
    _atomic_write(out / "r7.json", json.dumps(to_jsonable(r7), indent=2, allow_nan=False) + "\n")
    log(f"[R7] {r7['n_ok']}/{r7['n']} con max|diff| <= 1e-8 (maximo {r7['max_abs_diff']:.3g})")
    rows = [r for r in read_rows(jsonl) if r.get("config_hash") == chash]
    res = summarize(rows, cfg, len(tasks))
    final = {**meta, "phase": "final", "complete": len(rows) == len(tasks), "locality_validation": val, "result": res, "R7": r7}
    _atomic_write(summary_path, json.dumps(to_jsonable(final), indent=2, allow_nan=False) + "\n")
    log(f"[C0-L4] decision: {res['decision']}; clases genericos: {res['class_counts_generic']}")
    if res["decision"] == "CONTINUA":
        log(f"[C0-L4] celdas candidatas (NO se corre N=343/125 automaticamente): {[c['cell_idx'] for c in res['candidate_cells']]}")
    log_fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
