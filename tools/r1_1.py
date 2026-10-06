"""Omega R1-1, etapa A (y etapa B condicional): dinamica SQ (docs/OMEGA_R1_1_PRERREGISTRO.md, congelado).

Salida: results/r1_1/{runs.jsonl, pairs.jsonl, summary.json} (smoke: results/r1_1_smoke). Uso:
  python tools/r1_1.py [--procs 3] [--smoke] [--fresh] [--no-stage-b]
Unidad de trabajo = PAR (variante, inicio, semilla) con sus dos tamanos (N0, 8N0; etapa B: 2N0, 4N0). REANUDABLE: pairs.jsonl
se anade por par terminado (runs.jsonl se anade antes); al relanzar se saltan los pares presentes (--fresh borra).
Claves RNG: run = (20261018, id_variante, id_inicio, semilla, indice_de_tamano) (id_variante 1..18, id_inicio 1..5, tamanos
0..3 = N0, 8N0, 2N0, 4N0); grafo inicial = key; dinamica = key+(7,); RC-3 (tools/rc3.measure) = key + sus subclaves (5,1,3).
I1: 'completo con pesos' umbralizado W>1-p con p=20/(N-1) es EXACTAMENTE G(N,p) (cada W_ij>1-p independiente con prob. p), asi que se
genera como G(N,p): m ~ Binomial(N(N-1)/2, p) y un subconjunto uniforme de m aristas (erdos_renyi).
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
from collections import Counter  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402
from scipy import sparse  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from c0_dynamics import _atomic_write  # noqa: E402
from coh_l0b import prufer_tree  # noqa: E402
from p1d3_panel import erdos_renyi, ring, torus_lattice  # noqa: E402
from r1_wedge_census import census  # noqa: E402
from rc3 import measure, pair_exclusions, scaling_delta  # noqa: E402

from omega.c0.references import rng_from_key  # noqa: E402
from omega.dynamics.sq import (  # noqa: E402
    VARIANT_B,
    VARIANT_O,
    VARIANT_R,
    Graph,
    edge_flags,
    run_dynamics,
)
from omega.experiments.v11.gate import head_commit  # noqa: E402

MASTER = 20261018
OUT_FULL = ROOT / "results" / "r1_1"
OUT_SMOKE = ROOT / "results" / "r1_1_smoke"
VARIANTS = [(b, r, o) for b in VARIANT_B for r in VARIANT_R for o in VARIANT_O]  # 18; id = indice + 1
INITS = {1: "I1", 2: "I2", 3: "I3", 4: "I4", 5: "I5"}
AMORPHOUS = (1, 2, 3, 4)
T_MAX = 100
ATTRACTOR_FREQ = 0.8


def sizes(smoke: bool) -> dict[str, Any]:
    n0 = 300 if smoke else 5000
    return {"n": [n0, 8 * n0, 2 * n0, 4 * n0], "side": [7 if smoke else 17, 14 if smoke else 34, 8, 10]}


def build_init(init: int, si: int, sz: dict[str, Any], rng: np.random.Generator) -> sparse.csr_array:
    n = sz["n"][si]
    if init == 1:
        p = 20.0 / (n - 1)
        m = int(rng.binomial(n * (n - 1) // 2, p))
        return erdos_renyi(n, m, rng)
    if init == 2:
        return erdos_renyi(n, 2 * n, rng)
    if init == 3:
        return prufer_tree(n, rng)
    if init == 4:
        return ring(n, 2)
    if init == 5:
        return torus_lattice(sz["side"][si], 3)
    raise ValueError(init)


def specs(smoke: bool, stage: str = "A", triggered: list[tuple[str, str]] | None = None) -> list[dict[str, Any]]:
    """Etapa A: 18 variantes x 5 inicios x semillas 0-2 (par N0/8N0). Etapa B: para cada (B,R) disparado, 3 ordenes x I1-I4:
    semillas 3-4 en el par N0/8N0 y semillas 0-4 en el par 2N0/4N0 (I5 no se repite)."""
    out = []
    if stage == "A":
        for vi, (b, r, o) in enumerate(VARIANTS):
            for init in INITS:
                for seed in ((0,) if smoke else (0, 1, 2)):
                    out.append(_spec("A", vi, b, r, o, init, seed, "N0_8N0"))
    else:
        for vi, (b, r, o) in enumerate(VARIANTS):
            if (b, r) not in (triggered or []):
                continue
            for init in AMORPHOUS:
                for seed in ((3,) if smoke else (3, 4)):
                    out.append(_spec("B", vi, b, r, o, init, seed, "N0_8N0"))
                for seed in ((0,) if smoke else (0, 1, 2, 3, 4)):
                    out.append(_spec("B", vi, b, r, o, init, seed, "2N0_4N0"))
    return out


def _spec(stage: str, vi: int, b: str, r: str, o: str, init: int, seed: int, sp: str) -> dict[str, Any]:
    return {"id": f"{stage}|v{vi + 1}|I{init}|s{seed}|{sp}", "stage": stage, "variant_id": vi + 1, "B": b, "R": r, "O": o,
            "init": init, "seed": seed, "size_pair": sp, "sizes_idx": (0, 1) if sp == "N0_8N0" else (2, 3)}


# ------------------------------------------------------------------ una corrida


def _deg_summary(deg: np.ndarray) -> dict[str, Any]:
    c = Counter(np.minimum(deg, 12).tolist())
    return {"mean": float(deg.mean()), "max": int(deg.max()), "q": [float(x) for x in np.quantile(deg, [0.1, 0.5, 0.9])],
            "hist": {("12+" if k == 12 else str(k)): round(v / deg.size, 4) for k, v in sorted(c.items())}}


def run_one(spec: dict[str, Any], si: int, smoke: bool) -> tuple[dict[str, Any], sparse.csr_array]:
    t0 = time.perf_counter()
    sz = sizes(smoke)
    key = (MASTER, spec["variant_id"], spec["init"], spec["seed"], si)
    a0 = build_init(spec["init"], si, sz, rng_from_key(key))
    g = Graph.from_csr(a0)
    n = g.n
    init_edges = g.edge_set()
    res = run_dynamics(g, spec["B"], spec["R"], spec["O"], rng_from_key(key + (7,)), T=10 if smoke else T_MAX)
    fin = g.to_csr()
    nc, lab = connected_components(fin, directed=False)
    giant = int(np.bincount(lab).max())
    row = {"id": spec["id"] + f"|z{si}", "pair_id": spec["id"], "stage": spec["stage"], "variant_id": spec["variant_id"],
           "B": spec["B"], "R": spec["R"], "O": spec["O"], "init": spec["init"], "seed": spec["seed"], "size_idx": si,
           "N": n, "key": list(key), "E_initial": len(init_edges), "E_final": g.n_edges, "sweeps": res["sweeps"],
           "absorbed": res["absorbed"], "final_equals_initial": g.edge_set() == init_edges, "giant": giant,
           "giant_frac": giant / n, "n_components": int(nc), "final_deg": _deg_summary(g.degrees()),
           "trajectory": res["trajectory"], "seconds": round(time.perf_counter() - t0, 1)}
    return row, fin


def census_extra(fin: sparse.csr_array) -> dict[str, Any]:
    g = Graph.from_csr(fin)
    out: dict[str, Any] = {"t_census": census(fin), "degree": _deg_summary(g.degrees())}
    if g.n_edges:
        _, _, s, t = edge_flags(fin)
        out.update({"frac_edges_s": float(s.mean()), "frac_edges_t": float(t.mean())})
    return out


def run_pair(spec: dict[str, Any], smoke: bool) -> dict[str, Any]:
    t0 = time.perf_counter()
    runs, fins = [], []
    for si in spec["sizes_idx"]:
        r, f = run_one(spec, si, smoke)
        runs.append(r)
        fins.append(f)
    base = {k: spec[k] for k in ("id", "stage", "variant_id", "B", "R", "O", "init", "seed", "size_pair")}
    pair: dict[str, Any] = dict(base)
    pair.update({"N": [r["N"] for r in runs], "giant_frac": [r["giant_frac"] for r in runs],
                 "absorbed": [r["absorbed"] for r in runs], "sweeps": [r["sweeps"] for r in runs],
                 "E_over_N_final": [r["E_final"] / r["N"] for r in runs],
                 "final_equals_initial": [r["final_equals_initial"] for r in runs]})
    if min(pair["giant_frac"]) < 0.5:
        pair.update({"verdict": "FRAGMENTADO", "exclusions": [], "rc3": None})
    else:
        ms = []
        for r, f in zip(runs, fins):
            try:
                ms.append(measure(f, tuple(r["key"])))
            except Exception as e:  # abstencion declarada: el instrumento no es evaluable
                ms.append({"n": r["N"], "n_giant": None, "R": None, "annulus_status": "NO_EVALUABLE", "coherence_status": "NO_EVALUABLE",
                           "annulus": None, "coherence": None, "guards": [f"{type(e).__name__}:{e}"], "n_scales": None, "r_w": None})
        a, b = ms
        delta = scaling_delta(a["R"], b["R"], a["n_giant"] or 0, b["n_giant"] or 0)
        ex = pair_exclusions(delta, (a["annulus_status"], b["annulus_status"]), (a["coherence_status"], b["coherence_status"]))
        if not ex:
            verdict = "NO-EXCLUIDO"
        elif ex == ["X4"]:
            verdict = "EXCLUIDO-X4-marginal"
        else:
            verdict = "EXCLUIDO(" + ",".join(ex) + ")"

        def per_size(m: dict[str, Any]) -> dict[str, Any]:
            return {"n_giant": m["n_giant"], "R": m["R"], "annulus_status": m["annulus_status"],
                    "coherence_status": m["coherence_status"], "f_med": (m["annulus"] or {}).get("f_med"),
                    "c_med": (m["annulus"] or {}).get("c_med"), "n_scales": m["n_scales"], "r_w": m["r_w"],
                    "gamma": (m["coherence"] or {}).get("median_gamma"), "f_deg": (m["coherence"] or {}).get("f_deg"),
                    "f_low": (m["coherence"] or {}).get("f_low"), "guards": m["guards"]}

        pair.update({"verdict": verdict, "exclusions": ex, "delta": delta, "rc3": [per_size(a), per_size(b)]})
        if verdict in ("NO-EXCLUIDO", "EXCLUIDO-X4-marginal"):
            try:
                pair["inspection"] = [census_extra(f) for f in fins]
            except MemoryError:
                pair["inspection"] = "MemoryError"
    pair["seconds"] = round(time.perf_counter() - t0, 1)
    return {"runs": runs, "pair": pair}


_SPECS: list[dict[str, Any]] = []
_SMOKE = False


def _run_idx(i: int) -> dict[str, Any]:
    return run_pair(_SPECS[i], _SMOKE)


# ------------------------------------------------------------------ agregacion (§7)


def _dominant(classes: list[str]) -> dict[str, Any]:
    c = Counter(classes)
    top = sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))[0]
    freq = top[1] / len(classes) if classes else 0.0
    return {"n": len(classes), "dominant": top[0] if classes else None, "freq": freq,
            "attractor": bool(classes) and freq >= ATTRACTOR_FREQ,
            "result": top[0] if classes and freq >= ATTRACTOR_FREQ else "MULTIESTABLE", "distribution": dict(c)}


def aggregate(pairs: list[dict[str, Any]], stage: str) -> dict[str, Any]:
    sel = [p for p in pairs if p["stage"] == stage and p["init"] in AMORPHOUS]
    table = {}
    for b in VARIANT_B:
        for r in VARIANT_R:
            table[f"{b}/{r}"] = _dominant([p["verdict"] for p in sel if p["B"] == b and p["R"] == r])
    return table


def stage_a_report(pairs: list[dict[str, Any]], runs: list[dict[str, Any]], smoke: bool) -> dict[str, Any]:
    A = [p for p in pairs if p["stage"] == "A"]
    table = aggregate(pairs, "A")
    by_init: dict[str, Any] = {}
    for key in table:
        b, r = key.split("/")
        for init in INITS:
            by_init[f"{key}|{INITS[init]}"] = dict(Counter(p["verdict"] for p in A if p["B"] == b and p["R"] == r and p["init"] == init))
    by_order: dict[str, Any] = {}
    for key in table:
        b, r = key.split("/")
        for o in VARIANT_O:
            by_order[f"{key}|{o}"] = dict(Counter(p["verdict"] for p in A if p["B"] == b and p["R"] == r and p["O"] == o and p["init"] in AMORPHOUS))
    i5 = {}
    for p in A:
        if p["init"] == 5:
            e = i5.setdefault(f"{p['B']}/{p['R']}/{p['O']}", {"n": 0, "fixed_both_sizes": 0, "verdicts": Counter()})
            e["n"] += 1
            e["fixed_both_sizes"] += int(all(p["final_equals_initial"]))
            e["verdicts"][p["verdict"]] += 1
    for e in i5.values():
        e["verdicts"] = dict(e["verdicts"])
    A_runs = [r for r in runs if r["stage"] == "A"]
    absorb = {"runs": len(A_runs), "absorbed": sum(r["absorbed"] for r in A_runs),
              "absorbed_by_B_init": {f"{b}|{INITS[i]}": [sum(r["absorbed"] for r in A_runs if r["B"] == b and r["init"] == i),
                                                         sum(1 for r in A_runs if r["B"] == b and r["init"] == i)]
                                     for b in VARIANT_B for i in INITS}}
    non_excl = [p for p in A if p["verdict"] in ("NO-EXCLUIDO", "EXCLUIDO-X4-marginal")]
    trig = sorted({(p["B"], p["R"]) for p in A if p["verdict"] == "NO-EXCLUIDO"})
    # trayectorias tipicas: una corrida por (B, orden) (R2, I2, semilla 0, N0)
    traj = {}
    for b in VARIANT_B:
        for o in VARIANT_O:
            c = [r for r in A_runs if r["B"] == b and r["O"] == o and r["R"] == "R2" and r["init"] == 2 and r["seed"] == 0 and r["size_idx"] == 0]
            if c:
                traj[f"{b}|{o}"] = {"id": c[0]["id"], "sweeps": c[0]["sweeps"], "absorbed": c[0]["absorbed"], "trajectory": c[0]["trajectory"]}
    return {"n_pairs": len(A), "table_BR": table, "by_init": by_init, "by_order": by_order, "I5_control": i5, "absorption": absorb,
            "verdict_counts": dict(Counter(p["verdict"] for p in A)),
            "inspected": [{k: p[k] for k in ("id", "verdict", "B", "R", "O", "init", "seed", "delta", "rc3", "inspection", "giant_frac")} for p in non_excl],
            "triggered_BR": [list(t) for t in trig], "typical_trajectories": traj}


def final_summary(pairs: list[dict[str, Any]], runs: list[dict[str, Any]], smoke: bool, secs: float) -> dict[str, Any]:
    rep = stage_a_report(pairs, runs, smoke)
    trig = [tuple(t) for t in rep["triggered_BR"]]
    out: dict[str, Any] = {"step": "r1_1", "code_commit": head_commit(), "smoke": smoke, "seconds_this_run": round(secs, 1),
                           "stage_A": rep, "stage_B_triggered": bool(trig)}
    positives = []
    if trig:
        tb = aggregate(pairs, "B")
        out["stage_B"] = {"table_BR": {k: tb[k] for k in tb if tuple(k.split("/")) in trig},
                          "n_pairs": sum(1 for p in pairs if p["stage"] == "B")}
        for b, r in trig:
            ka, kb = rep["table_BR"][f"{b}/{r}"], tb.get(f"{b}/{r}")
            if ka["result"] == "NO-EXCLUIDO" and kb and kb["result"] == "NO-EXCLUIDO":
                positives.append(f"{b}/{r}")
    out["positive_BR"] = positives
    out["omega6_null_rate"] = len(positives) / 6
    out["omega6_note"] = "positivo = atractor NO-EXCLUIDO en A y mantenido en B (§7); sin disparo de B, 0/6"
    return out


# ------------------------------------------------------------------ ejecucion


def _load(path: Path) -> dict[str, dict[str, Any]]:
    done: dict[str, dict[str, Any]] = {}
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                try:
                    r = json.loads(line)
                    done[r["id"]] = r
                except json.JSONDecodeError:
                    pass
    return done


def fmt(p: dict[str, Any]) -> str:
    return (f"  {p['id']:<28} {p['B']}/{p['R']}/{p['O']:<5} N={p['N']} giant={[round(x, 2) for x in p['giant_frac']]} "
            f"E/N={[round(x, 2) for x in p['E_over_N_final']]} abs={p['absorbed']} sw={p['sweeps']} {p['verdict']} {p['seconds']}s")


def execute(todo_specs: list[dict[str, Any]], done: dict[str, dict[str, Any]], out: Path, procs: int) -> None:
    global _SPECS
    _SPECS = todo_specs
    todo = [i for i, s in enumerate(todo_specs) if s["id"] not in done]
    todo.sort(key=lambda i: -todo_specs[i]["sizes_idx"][1])
    print(f"{len(todo_specs)} pares ({len(todo_specs) - len(todo)} ya hechos), procs={procs}", flush=True)
    if not todo:
        return
    with mp.get_context("fork").Pool(procs) as pool, open(out / "pairs.jsonl", "a") as fp, open(out / "runs.jsonl", "a") as fr:
        for res in pool.imap_unordered(_run_idx, todo, chunksize=1):
            for r in res["runs"]:
                fr.write(json.dumps(r, default=str) + "\n")
            fr.flush()
            fp.write(json.dumps(res["pair"], default=str) + "\n")
            fp.flush()
            done[res["pair"]["id"]] = res["pair"]
            print(fmt(res["pair"]), flush=True)


def main(argv: list[str] | None = None) -> int:
    global _SMOKE
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--smoke", action="store_true", help="N0=300/2400, 1 semilla, 10 barridos -> results/r1_1_smoke/")
    ap.add_argument("--fresh", action="store_true")
    ap.add_argument("--no-stage-b", action="store_true", help="no ejecuta la etapa B aunque se dispare")
    args = ap.parse_args(argv)
    _SMOKE = args.smoke
    out = OUT_SMOKE if args.smoke else OUT_FULL
    out.mkdir(parents=True, exist_ok=True)
    if args.fresh:
        for f in ("pairs.jsonl", "runs.jsonl"):
            (out / f).unlink(missing_ok=True)
    done = _load(out / "pairs.jsonl")
    t0 = time.perf_counter()
    sa = specs(args.smoke, "A")
    execute(sa, done, out, args.procs)
    pairs_a = [done[s["id"]] for s in sa]
    trig = sorted({(p["B"], p["R"]) for p in pairs_a if p["verdict"] == "NO-EXCLUIDO"})
    print(f"Etapa B disparada por: {trig if trig else 'nadie (no se ejecuta)'}", flush=True)
    if trig and not args.no_stage_b:
        execute(specs(args.smoke, "B", trig), done, out, args.procs)
    # reescritura atomica ordenada y deduplicada
    all_specs = sa + (specs(args.smoke, "B", trig) if trig and not args.no_stage_b else [])
    pairs = [done[s["id"]] for s in all_specs if s["id"] in done]
    runs_d = _load(out / "runs.jsonl")
    runs = [runs_d[k] for k in sorted(runs_d)]
    _atomic_write(out / "pairs.jsonl", "".join(json.dumps(p, default=str) + "\n" for p in pairs))
    _atomic_write(out / "runs.jsonl", "".join(json.dumps(r, default=str) + "\n" for r in runs))
    summ = final_summary(pairs, runs, args.smoke, time.perf_counter() - t0)
    _atomic_write(out / "summary.json", json.dumps(summ, indent=2, default=str))
    for k, v in summ["stage_A"]["table_BR"].items():
        print(f"{k:<10} {v['result']:<24} freq={v['freq']:.2f} {v['distribution']}")
    print("Etapa B:", summ["stage_B_triggered"], "positivos:", summ["positive_BR"], "Omega6:", summ["omega6_null_rate"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
