"""ES-C0b: nucleo mixto (movimiento R de reinsercion + local) y rejilla recalibrada (docs/OMEGA_ESC0B_PREREGISTRO.md).

Diagnostico de mezcla; no es evidencia. Ejecutor reanudable por partes (una parte = una cadena):
  run [--stage 1|2] [--N 64] [--eps-idx 0|1] [--workers 4] [--sweeps-override k --out DIR --allow-dirty]
  summarize [--out DIR]   -> <out>/summary.json (+ tabla, con comparacion con ES-C0 b = 0)
  figure [--out DIR] [--fig PNG]

Semillas: SeedSequence(20261107, spawn_key=(330, etapa, eps_idx, N, kappa_idx, cadena)); en la etapa 1 kappa_idx = 0.
Etapa 1: beta = 0, eps in {1, 0.25}, N in {64, 128}, 4 inicios de ES-C0 (tools/esc0_piloto.start_order), 1000 / 500
barridos, p_r = 0.1. Etapa 2: solo celdas (eps, N) cuya etapa 1 mezclo (criterio de ES-C §4, mismo codigo de
diagnostico que ES-C0); sigma0 = desviacion (ddof = 1) de S_eps en las mitades finales de las 4 cadenas de la etapa
1; beta = kappa / sigma0, kappa in {-4, -1, 1, 4}. Clasificador: los archivos oficiales de results/omega_esc0/.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # matrices chicas; evita sobresuscripcion con --workers

import argparse
import json
import math
import multiprocessing as mp
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

_TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(_TOOLS.parent))
sys.path.insert(0, str(_TOOLS))

from es0_piloto import _clean, git_info
from esc0_piloto import (
    ESS_OK,
    OBS_MIX,
    RHAT_LOST,
    RHAT_OK,
    START_NAMES,
    STAY_UNFLAGGED,
    TRACE_KEYS,
    cell_summary,
    measure,
    start_order,
)
from esc0_piloto import OUT as ESC0_OUT
from esc0_piloto import clf_path as esc0_clf_path
from esc0_piloto import part_path as esc0_part_path

from omega.causality.mcmc_bd_eps import action_eps, init_state
from omega.causality.mcmc_reinsert import sweep_mixed

EXPERIMENT_ID = 330
ENT_CHAIN = 20261107
EPS = (1.0, 0.25)
NS = (64, 128)
SWEEPS = {64: 1000, 128: 500}
KAPPAS = (-4.0, -1.0, 1.0, 4.0)
P_R = 0.1
N_CHAINS = 4
ESC0_B0_IDX = 2  # b = 0 en la rejilla de ES-C0 (BS = (-12, -4, 0, 4, 12))
SPEC = Path("docs") / "OMEGA_ESC0B_PREREGISTRO.md"
OUT = Path("results") / "omega_esc0b"
FIG = Path("docs") / "figuras" / "esc0b_piloto.png"
COUNT_KEYS = ("local_proposed", "local_valid", "local_accepted", "r_proposed", "r_applicable", "r_accepted",
              "r_changed")


def rng_for(*key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(ENT_CHAIN, spawn_key=(EXPERIMENT_ID, *key))))


def part_path(out: Path, stage: int, ei: int, n: int, ki: int, chain: int) -> Path:
    tag = f"s1_e{ei}_N{n}" if stage == 1 else f"s2_e{ei}_N{n}_k{ki}"
    return out / "parts" / f"{tag}_c{chain}.json"


def load_json(p: Path) -> dict[str, Any] | None:
    if not p.exists():
        return None
    rec: dict[str, Any] = json.loads(p.read_text())
    return rec


def load_classifier(n: int) -> dict[str, Any]:
    """Clasificador oficial de ES-C0 (no se recalibra)."""
    p = esc0_clf_path(ESC0_OUT, n)
    rec = load_json(p)
    if rec is None:
        sys.exit(f"falta el clasificador oficial {p}")
    if not rec.get("official"):
        sys.exit(f"{p} no es oficial")
    return rec


# --- una cadena -------------------------------------------------------------------------------------


def run_part(task: dict[str, Any]) -> str:
    stage, ei, n, ki, chain = task["stage"], task["eps_idx"], task["N"], task["kappa_idx"], task["chain"]
    out = Path(task["out"])
    p = part_path(out, stage, ei, n, ki, chain)
    if p.exists():
        return f"{p.name} ya existe"
    clf = load_classifier(n)
    ref4, delta4 = np.asarray(clf["ref4"], dtype=np.float64), float(clf["delta4"])
    eps, beta = EPS[ei], float(task["beta"])
    rng = rng_for(stage, ei, n, ki, chain)
    t0 = time.perf_counter()
    st = init_state(start_order(n, chain, rng), eps)
    t_init = time.perf_counter() - t0
    traces: dict[str, list[Any]] = {k: [] for k in TRACE_KEYS}
    counts: dict[str, list[int]] = {k: [] for k in (*COUNT_KEYS, "valid", "accepted")}
    wall: list[float] = []
    t_meas = 0.0
    initial = measure(st, ref4, delta4, True)
    for _ in range(task["sweeps"]):
        ts = time.perf_counter()
        cnt = sweep_mixed(st, beta, rng, P_R)
        tm = time.perf_counter()
        wall.append(tm - ts)
        for k in COUNT_KEYS:
            counts[k].append(cnt[k])
        # para esc0_piloto.cell_summary: propuestas validas = locales validas + R aplicables
        counts["valid"].append(cnt["local_valid"] + cnt["r_applicable"])
        counts["accepted"].append(cnt["local_accepted"] + cnt["r_accepted"])
        row = measure(st, ref4, delta4, True)
        t_meas += time.perf_counter() - tm
        for k in TRACE_KEYS:
            traces[k].append(row[k])
    drift = float(st.s - action_eps(st.c, st.b, st.g, st.const))
    rec = {
        "spec": {"doc": str(SPEC), "experiment_id": EXPERIMENT_ID, "entropy": ENT_CHAIN, "stage": stage,
                 "eps_idx": ei, "eps": eps, "N": n, "kappa_idx": ki,
                 "kappa": KAPPAS[ki] if stage == 2 else None, "sigma0": task.get("sigma0"), "beta": beta,
                 "chain": chain, "start": START_NAMES[chain], "sweeps": task["sweeps"],
                 "proposals_per_sweep": n * (n - 1) // 2, "p_r": P_R, "flag_every": 1,
                 "spawn_key": [EXPERIMENT_ID, stage, ei, n, ki, chain], "official": task["official"],
                 "ignore_mixing": task.get("ignore_mixing", False)},
        "commit": task["commit"], "dirty": task["dirty"], "clasificador_estado": clf["estado"],
        "clasificador_commit": clf["commit"],
        "initial": initial, "traces": traces, "counts": counts, "wall_s": wall, "init_s": t_init,
        "measure_s": t_meas, "s_drift_final": drift,
    }
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(_clean(rec)))
    tmp.replace(p)
    return f"{p.name} ok {sum(wall):.1f} s barridos + {t_meas:.1f} s medicion"


def stage_parts(out: Path, stage: int, ei: int, n: int, ki: int = 0) -> list[dict[str, Any]]:
    return [p for ch in range(N_CHAINS) if (p := load_json(part_path(out, stage, ei, n, ki, ch))) is not None]


def stage1_cell(out: Path, ei: int, n: int) -> dict[str, Any] | None:
    """Mezcla y sigma0 de la celda de etapa 1 (None si faltan cadenas)."""
    parts = stage_parts(out, 1, ei, n)
    if len(parts) < N_CHAINS:
        return None
    cs = cell_summary(parts)
    fin = np.concatenate([np.asarray(p["traces"]["S"], dtype=np.float64)[len(p["traces"]["S"]) // 2:]
                          for p in parts])
    return {"mezcla": bool(cs.get("mezcla")), "sigma0": float(np.std(fin, ddof=1)),
            "official": all(p["spec"]["official"] for p in parts)}


def cmd_run(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    if args.sweeps_override is not None and not args.out:
        sys.exit("--sweeps-override exige --out (no se mezcla con la corrida oficial)")
    head, dirty = git_info()
    if dirty and not args.allow_dirty:
        sys.exit("exige arbol git limpio fuera de results/ (usar --allow-dirty solo para pruebas)")
    official = not dirty and args.sweeps_override is None and out == OUT and not args.ignore_mixing
    if args.ignore_mixing and out == OUT:
        sys.exit("--ignore-mixing es solo para pruebas (exige --out)")
    base = {"out": str(out), "commit": head, "dirty": dirty, "official": official}
    if official:  # una corrida oficial no reutiliza partes no oficiales
        stale = [p.name for p in out.glob("parts/*.json") if not json.loads(p.read_text())["spec"]["official"]]
        if stale:
            sys.exit(f"{out} contiene partes no oficiales ({len(stale)}, p. ej. {stale[0]}); moverlas antes")
    ns = [args.N] if args.N else list(NS)
    eis = [args.eps_idx] if args.eps_idx is not None else list(range(len(EPS)))
    for n in ns:
        clf = load_classifier(n)
        print(f"clasificador N={n} (ES-C0, commit {clf['commit'][:8]}): {clf['estado']}")
    tasks: list[dict[str, Any]] = []
    plan: dict[str, Any] = {}
    for n in ns:
        sw = args.sweeps_override if args.sweeps_override is not None else SWEEPS[n]
        for ei in eis:
            if args.stage == 1:
                for ch in range(N_CHAINS):
                    tasks.append({**base, "stage": 1, "eps_idx": ei, "N": n, "kappa_idx": 0, "chain": ch,
                                  "sweeps": sw, "beta": 0.0})
                continue
            cell = stage1_cell(out, ei, n)
            key = f"e{ei}_N{n}"
            if cell is None:
                sys.exit(f"etapa 1 incompleta en {key}: correr --stage 1 antes")
            if official and not cell["official"]:
                sys.exit(f"etapa 1 de {key} no es oficial")
            if not cell["mezcla"] and not args.ignore_mixing:
                plan[key] = {"eps": EPS[ei], "N": n, "etapa1_mezcla": False, "sigma0": cell["sigma0"],
                             "omitida": "la etapa 1 no mezcló (§3)"}
                print(f"{key}: etapa 1 no mezcló -> se omite la etapa 2")
                continue
            sigma0 = cell["sigma0"]
            if not (math.isfinite(sigma0) and sigma0 > 0):
                plan[key] = {"eps": EPS[ei], "N": n, "etapa1_mezcla": cell["mezcla"], "sigma0": sigma0,
                             "omitida": "sigma0 no positivo"}
                continue
            betas = [k / sigma0 for k in KAPPAS]
            plan[key] = {"eps": EPS[ei], "N": n, "etapa1_mezcla": cell["mezcla"], "sigma0": sigma0,
                         "kappas": list(KAPPAS), "betas": betas, "ignore_mixing": bool(args.ignore_mixing)}
            for ki, beta in enumerate(betas):
                for ch in range(N_CHAINS):
                    tasks.append({**base, "stage": 2, "eps_idx": ei, "N": n, "kappa_idx": ki, "chain": ch,
                                  "sweeps": sw, "beta": beta, "sigma0": sigma0,
                                  "ignore_mixing": bool(args.ignore_mixing)})
    if args.stage == 2:
        out.mkdir(parents=True, exist_ok=True)
        pp = out / "etapa2_plan.json"
        old = load_json(pp) or {}
        (out / "etapa2_plan.json").write_text(json.dumps(_clean({**old, **plan}), indent=2, allow_nan=False))
    tasks = [t for t in tasks
             if not part_path(out, t["stage"], t["eps_idx"], t["N"], t["kappa_idx"], t["chain"]).exists()]
    tasks.sort(key=lambda t: -t["N"] * t["N"] * t["sweeps"])  # las mas largas primero
    print(f"etapa {args.stage}: {len(tasks)} partes por correr en {out} "
          f"(commit {head[:8]}, dirty {dirty}, oficial {official})")
    t0 = time.time()
    if args.workers <= 1:
        for t in tasks:
            print(run_part(t), flush=True)
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            for msg in pool.imap_unordered(run_part, tasks):
                print(msg, flush=True)
    print(f"pared total {time.time() - t0:.1f} s")


# --- resumen ----------------------------------------------------------------------------------------


def r_stats(parts: list[dict[str, Any]]) -> dict[str, Any]:
    tot = {k: sum(sum(p["counts"][k]) for p in parts) for k in COUNT_KEYS}
    rp, lp = max(tot["r_proposed"], 1), max(tot["local_proposed"], 1)
    per_chain = {}
    for p in parts:
        c = {k: sum(p["counts"][k]) for k in COUNT_KEYS}
        per_chain[str(p["spec"]["chain"])] = {
            "r_aplicable_frac": c["r_applicable"] / max(c["r_proposed"], 1),
            "r_aceptada_por_aplicable": c["r_accepted"] / max(c["r_applicable"], 1),
            "r_cambia_frac": c["r_changed"] / max(c["r_proposed"], 1)}
    return {**tot, "r_aplicable_frac": tot["r_applicable"] / rp,
            "r_aceptada_por_propuesta": tot["r_accepted"] / rp,
            "r_aceptada_por_aplicable": tot["r_accepted"] / max(tot["r_applicable"], 1),
            "r_cambia_por_propuesta": tot["r_changed"] / rp,
            "local_valida_frac": tot["local_valid"] / lp, "local_aceptada_frac": tot["local_accepted"] / lp,
            "por_cadena": per_chain}


def esc0_b0(ei: int, n: int) -> dict[str, Any] | None:
    s = load_json(ESC0_OUT / "summary.json")
    if s is None:
        return None
    c = s["cells"].get(f"e{ei}_N{n}_b{ESC0_B0_IDX}")
    if c is None:
        return None
    mx = c.get("mezcla_detalle", {})
    return {"mezcla": c.get("mezcla"), "official": c.get("official"),
            "rhat": {o: mx[o]["rhat"] for o in OBS_MIX if o in mx},
            "ess": {o: mx[o]["ess"] for o in OBS_MIX if o in mx},
            "rel_frac_mitad_final": {k: v["rel_frac_mitad_final"] for k, v in c["por_cadena"].items()},
            "s_per_sweep": c.get("s_per_sweep")}


def full_cell(parts: list[dict[str, Any]]) -> dict[str, Any]:
    cs = cell_summary(parts)
    fin = [np.asarray(p["traces"]["S"], dtype=np.float64)[len(p["traces"]["S"]) // 2:] for p in parts]
    sd = float(np.std(np.concatenate(fin), ddof=1)) if sum(x.size for x in fin) > 1 else float("nan")
    return {**cs, "sigma_S_mitad_final": sd, "R": r_stats(parts)}


def cmd_summarize(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    cells: dict[str, Any] = {}
    for ei, eps in enumerate(EPS):
        for n in NS:
            parts = stage_parts(out, 1, ei, n)
            if parts:
                key = f"s1_e{ei}_N{n}"
                cells[key] = {"etapa": 1, "eps": eps, "N": n, "kappa": None, "beta": 0.0, **full_cell(parts),
                              "esc0_b0": esc0_b0(ei, n)}
                if len(parts) < N_CHAINS:
                    cells[key]["incompleta"] = True
            for ki, kap in enumerate(KAPPAS):
                parts = stage_parts(out, 2, ei, n, ki)
                if parts:
                    key = f"s2_e{ei}_N{n}_k{ki}"
                    cells[key] = {"etapa": 2, "eps": eps, "N": n, "kappa": kap, "beta": parts[0]["spec"]["beta"],
                                  "sigma0": parts[0]["spec"]["sigma0"], **full_cell(parts)}
                    if len(parts) < N_CHAINS:
                        cells[key]["incompleta"] = True
    res = {"spec": str(SPEC), "experiment_id": EXPERIMENT_ID, "entropy": ENT_CHAIN, "p_r": P_R,
           "criterio": {"rhat": RHAT_OK, "ess": ESS_OK, "rhat_no_estimable": RHAT_LOST, "observables": list(OBS_MIX),
                        "mitad": "final", "flag_every": 1, "salida_sprinkling_mediciones": STAY_UNFLAGGED},
           "clasificador": {str(n): {k: c[k] for k in ("estado", "validacion", "delta4", "commit", "official")}
                            for n in NS if (c := load_json(esc0_clf_path(ESC0_OUT, n))) is not None},
           "etapa2_plan": load_json(out / "etapa2_plan.json"), "cells": cells}
    (out / "summary.json").write_text(json.dumps(_clean(res), indent=2, allow_nan=False))
    print(f"{'celda':<16}{'eps':>5}{'N':>5}{'kappa':>6} {'s/barr':>7} {'R aplic':>8} {'R acep':>7} "
          f"{'mezcla':>7} {'Rhat S/rel/h':>20} {'ESS min':>8} {'salida':>7} | ES-C0 b=0: {'mezcla':>6} "
          f"{'Rhat S/rel/h':>20} {'ESS min':>8}")
    for k, c in cells.items():
        mx = c.get("mezcla_detalle", {})
        rh = "/".join(f"{mx[o]['rhat']:.2f}" for o in OBS_MIX) if mx else "-"
        ess = f"{min(mx[o]['ess'] for o in OBS_MIX):.1f}" if mx else "-"
        ref = c.get("esc0_b0")
        rr, re_, rm = "-", "-", "-"
        if ref and ref["rhat"]:
            rr = "/".join(f"{ref['rhat'][o]:.2f}" for o in OBS_MIX)
            re_ = f"{min(ref['ess'].values()):.1f}"
            rm = str(ref["mezcla"])
        kap_s = "-" if c["kappa"] is None else f"{c['kappa']:+g}"
        print(f"{k:<16}{c['eps']:>5}{c['N']:>5}{kap_s:>6} {c['s_per_sweep']:>7.3f} "
              f"{c['R']['r_aplicable_frac']:>8.3f} {c['R']['r_aceptada_por_aplicable']:>7.3f} "
              f"{c.get('mezcla')!s:>7} {rh:>20} {ess:>8} {c.get('sprinkling_4d_barrido_salida')!s:>7} | "
              f"{'':>11}{rm:>6} {rr:>20} {re_:>8}")
    plan = load_json(out / "etapa2_plan.json") or {}
    for k, v in plan.items():
        extra = (f", omitida ({v['omitida']})" if "omitida" in v
                 else ", betas " + ", ".join(f"{b:.4g}" for b in v["betas"]))
        print(f"etapa 2 {k}: sigma0 = {v['sigma0']:.4g}{extra}")


# --- figura -----------------------------------------------------------------------------------------


def cmd_figure(args: argparse.Namespace) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(args.out) if args.out else OUT
    fig_path = Path(args.fig) if args.fig else (out / "esc0b_piloto.png" if args.out else FIG)
    rows = [(ei, n) for ei in range(len(EPS)) for n in NS if stage_parts(out, 1, ei, n)]
    if not rows:
        sys.exit("no hay partes de la etapa 1")
    has2 = any(stage_parts(out, 2, ei, n, ki) for ei, n in rows for ki in range(len(KAPPAS)))
    ncol = 2 + (len(KAPPAS) if has2 else 0)
    colors = ("#1f77b4", "#d62728", "#2ca02c", "#9467bd")
    fig, axes = plt.subplots(len(rows), ncol, figsize=(3.3 * ncol + 0.6, 2.4 * len(rows) + 1.0), squeeze=False)

    def draw(ax: Any, parts: list[dict[str, Any]], title: str) -> None:
        for p in parts:
            ch = p["spec"]["chain"]
            y = [p["initial"]["rel_frac"], *p["traces"]["rel_frac"]]  # barrido 0 = inicio
            ax.plot(np.arange(len(y)), y, lw=0.7, color=colors[ch], label=START_NAMES[ch])
        if not parts:
            ax.text(0.5, 0.5, "sin datos", ha="center", va="center", fontsize=8, transform=ax.transAxes)
            ax.set_xticks([])
        ax.set_ylim(-0.03, 1.03)
        ax.tick_params(labelsize=7)
        ax.set_title(title, fontsize=8)

    for r, (ei, n) in enumerate(rows):
        e = f"ε = {EPS[ei]:g}, N = {n}"
        old = [p for ch in range(N_CHAINS)
               if (p := load_json(esc0_part_path(ESC0_OUT, ei, n, ESC0_B0_IDX, ch))) is not None]
        draw(axes[r, 0], old, f"ES-C0, solo local, b = 0\n{e}")
        draw(axes[r, 1], stage_parts(out, 1, ei, n), f"ES-C0b etapa 1, mixto, b = 0\n{e}")
        if has2:
            for ki, kap in enumerate(KAPPAS):
                parts = stage_parts(out, 2, ei, n, ki)
                ks = f"{kap:+g}".replace("-", "−")
                draw(axes[r, 2 + ki], parts, f"ES-C0b etapa 2, κ = {ks}\n{e}")
        axes[r, 0].set_ylabel("fracción de relaciones", fontsize=7)
    for c in range(ncol):
        axes[-1, c].set_xlabel("barrido", fontsize=7)
    h, lab = axes[0, 1].get_legend_handles_labels()
    leg = fig.legend(h, lab, loc="lower center", ncol=4, fontsize=8, frameon=False)
    for line in leg.get_lines():
        line.set_linewidth(2.0)
    fig.suptitle("ES-C0b: fracción de relaciones de las 4 cadenas, solo local (ES-C0) frente a núcleo mixto con R\n"
                 "(diagnóstico de mezcla, no es evidencia)", fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=110)
    print(fig_path)


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--stage", type=int, choices=(1, 2), default=1)
    r.add_argument("--N", type=int, choices=NS)
    r.add_argument("--eps-idx", type=int, choices=range(len(EPS)))
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--sweeps-override", type=int)
    r.add_argument("--out")
    r.add_argument("--allow-dirty", action="store_true")
    r.add_argument("--ignore-mixing", action="store_true",
                   help="solo pruebas: corre la etapa 2 aunque la etapa 1 no haya mezclado (exige --out)")
    s = sub.add_parser("summarize")
    s.add_argument("--out")
    f = sub.add_parser("figure")
    f.add_argument("--out")
    f.add_argument("--fig")
    args = ap.parse_args()
    {"run": cmd_run, "summarize": cmd_summarize, "figure": cmd_figure}[args.cmd](args)


if __name__ == "__main__":
    main()
