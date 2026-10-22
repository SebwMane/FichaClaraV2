"""ES-C0: diagnostico de mezcla de la accion BD 4D suavizada (docs/OMEGA_ESC_PREREGISTRO.md). No es evidencia.

Ejecutor reanudable por partes (una parte = una cadena):
  run [--N 32] [--eps-idx 0|1] [--workers 4] [--sweeps-override k --out DIR --allow-dirty]
  summarize [--out DIR]   -> <out>/summary.json
  figure [--out DIR] [--fig PNG]

Semillas: cadenas SeedSequence(20261104, spawn_key=(329, eps_idx, N, b_idx, cadena)); clasificador por N con
SeedSequence(20261105 | 20261106, spawn_key=(329, N, s)), misma estructura que tools/p4.py:classifier.
Rejilla: eps in {1, 0.25} (idx 0, 1); b in {-12, -4, 0, 4, 12} (idx 0..4); beta = b/N; N in {32, 64, 128}.
Cadenas: 0 orden vacio, 1 cadena total, 2 sprinkling 4D natural (p4.natural_sprinkle), 3 KR reetiquetado natural.
Se mide en cada barrido: S_eps, fraccion de relaciones, altura, links/N, d_MM y la bandera 4D (cada FLAG_EVERY = 1).
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
import numpy.typing as npt

_TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(_TOOLS.parent))
sys.path.insert(0, str(_TOOLS))

from es0_piloto import _clean, git_info, split_rhat, tau_sokal
from p4 import natural_sprinkle

from omega.causality.estimators import interval_abundances, mm_dimension
from omega.causality.mcmc_bd_eps import (
    EpsState,
    action_eps,
    init_state,
    sweep,
)
from omega.causality.mcmc_general import height
from omega.causality.orders import (
    is_partial_order,
    kr_order,
    random_bipartite,
)
from omega.causality.profile import (
    K_MAX,
    d_ref,
    profile_distance,
    reference_abundances,
)

EXPERIMENT_ID = 329
ENT_CHAIN, ENT_CAL, ENT_VAL = 20261104, 20261105, 20261106
EPS = (1.0, 0.25)
BS = (-12, -4, 0, 4, 12)
NS = (32, 64, 128)
SWEEPS = {32: 2000, 64: 1000, 128: 500}
N_CHAINS = 4
FLAG_EVERY = 1  # la bandera cuesta ~1% de un barrido: se mide en cada barrido
STAY_UNFLAGGED = 10  # mediciones seguidas sin bandera para fijar el barrido de salida de la cadena 2
RHAT_OK, ESS_OK, RHAT_LOST = 1.05, 200.0, 1.5
OBS_MIX = ("S", "rel_frac", "height")
TRACE_KEYS = ("S", "rel_frac", "height", "links_per_N", "d_mm", "flag")
SPEC = Path("docs") / "OMEGA_ESC_PREREGISTRO.md"
OUT = Path("results") / "omega_esc0"
FIG = Path("docs") / "figuras" / "esc0_piloto.png"
START_NAMES = ("orden vacío", "cadena total", "sprinkling 4D", "KR natural")
Bool = npt.NDArray[np.bool_]
Float = npt.NDArray[np.float64]


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


def part_path(out: Path, eps_idx: int, n: int, b_idx: int, chain: int) -> Path:
    return out / "parts" / f"e{eps_idx}_N{n}_b{b_idx}_c{chain}.json"


def clf_path(out: Path, n: int) -> Path:
    return out / f"clasificador_N{n}.json"


# --- clasificador (tools/p4.py:classifier, entropias propias) ---------------------------------------


def flag_from(nk: npt.NDArray[np.int64], dmm: float, ref4: Float, delta4: float) -> bool:
    return d_ref(dmm) == 4 and profile_distance(nk, ref4) <= delta4


def calibrate(n: int) -> dict[str, Any]:
    ref4 = reference_abundances(4, n, K_MAX, [rng_for(ENT_CAL, n, 1000 + s) for s in range(40)])
    dist = []
    for s in range(40):
        c = natural_sprinkle(n, rng_for(ENT_CAL, n, s))
        dist.append(profile_distance(interval_abundances(c, K_MAX), ref4))
    delta4 = 1.5 * float(np.percentile(dist, 95))

    def flag(c: Bool) -> bool:
        return flag_from(interval_abundances(c, K_MAX), mm_dimension(c), ref4, delta4)

    pos = float(np.mean([flag(natural_sprinkle(n, rng_for(ENT_VAL, n, s))) for s in range(100)]))
    kr = float(np.mean([flag(kr_order(n, rng_for(ENT_VAL, n, 500 + s))) for s in range(20)]))
    bi = float(np.mean([flag(random_bipartite(n, 0.5, rng_for(ENT_VAL, n, 700 + s))) for s in range(20)]))
    ok = pos >= 0.9 and kr <= 0.1 and bi <= 0.1
    return {"N": n, "ref4": ref4.tolist(), "delta4": delta4, "k_max": K_MAX,
            "validacion": {"pos": pos, "kr": kr, "bicapa": bi, "ok": ok},
            "estado": "VALIDADO" if ok else "INDETERMINADO"}


def ensure_classifier(out: Path, n: int, base: dict[str, Any]) -> dict[str, Any]:
    p = clf_path(out, n)
    if p.exists():
        rec: dict[str, Any] = json.loads(p.read_text())
        return rec
    t0 = time.perf_counter()
    rec = {**calibrate(n), "experiment_id": EXPERIMENT_ID, "entropias": [ENT_CAL, ENT_VAL],
           "commit": base["commit"], "dirty": base["dirty"], "official": base["official"]}
    rec["seconds"] = time.perf_counter() - t0
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(rec))
    tmp.replace(p)
    return rec


# --- inicios (§2) -----------------------------------------------------------------------------------


def natural_relabel(c: Bool) -> Bool:
    """Extension lineal por numero de ancestros (si a < b, a tiene estrictamente menos)."""
    perm = np.argsort(c.sum(axis=0), kind="stable")
    out: Bool = c[np.ix_(perm, perm)]
    return out


def start_order(n: int, chain: int, rng: np.random.Generator) -> Bool:
    if chain == 0:
        c = np.zeros((n, n), dtype=bool)
    elif chain == 1:
        c = np.triu(np.ones((n, n), dtype=bool), k=1)
    elif chain == 2:
        c = natural_sprinkle(n, rng)
    else:
        c = natural_relabel(kr_order(n, rng))
    if np.any(np.tril(c)) or not is_partial_order(c):
        raise RuntimeError(f"inicio {chain} sin etiquetado natural")
    return c


# --- una cadena -------------------------------------------------------------------------------------


def measure(st: EpsState, ref4: Float, delta4: float, with_flag: bool) -> dict[str, Any]:
    n = st.c.shape[0]
    rel = st.b[st.c]  # elementos entre cada par relacionado (= between_counts en los pares relacionados)
    nrel = int(rel.size)
    dmm = mm_dimension(st.c)
    row: dict[str, Any] = {"S": st.s, "rel_frac": nrel / (n * (n - 1) / 2), "height": height(st.c),
                           "links_per_N": int(np.count_nonzero(rel == 0)) / n,
                           "d_mm": dmm if math.isfinite(dmm) else None, "flag": None}
    if with_flag:
        nk = np.bincount(rel, minlength=K_MAX + 1)[: K_MAX + 1].astype(np.int64)
        row["flag"] = int(flag_from(nk, dmm, ref4, delta4))
    return row


def run_part(task: dict[str, Any]) -> str:
    ei, n, bi, chain = task["eps_idx"], task["N"], task["b_idx"], task["chain"]
    out = Path(task["out"])
    p = part_path(out, ei, n, bi, chain)
    if p.exists():
        return f"{p.name} ya existe"
    clf = json.loads(clf_path(out, n).read_text())
    ref4, delta4 = np.asarray(clf["ref4"], dtype=np.float64), float(clf["delta4"])
    eps, b = EPS[ei], BS[bi]
    beta = b / n
    rng = rng_for(ENT_CHAIN, ei, n, bi, chain)
    t0 = time.perf_counter()
    st = init_state(start_order(n, chain, rng), eps)
    t_init = time.perf_counter() - t0
    traces: dict[str, list[Any]] = {k: [] for k in TRACE_KEYS}
    counts: dict[str, list[int]] = {"valid": [], "accepted": []}
    wall: list[float] = []
    t_meas = 0.0
    initial = measure(st, ref4, delta4, True)
    for sw in range(task["sweeps"]):
        ts = time.perf_counter()
        a, v = sweep(st, beta, rng)
        tm = time.perf_counter()
        wall.append(tm - ts)
        counts["accepted"].append(a)
        counts["valid"].append(v)
        row = measure(st, ref4, delta4, (sw + 1) % FLAG_EVERY == 0)
        t_meas += time.perf_counter() - tm
        for k in TRACE_KEYS:
            traces[k].append(row[k])
    drift = float(st.s - action_eps(st.c, st.b, st.g, st.const))
    rec = {
        "spec": {"doc": str(SPEC), "experiment_id": EXPERIMENT_ID, "entropy": ENT_CHAIN, "eps_idx": ei, "eps": eps,
                 "N": n, "b_idx": bi, "b": b, "beta": beta, "chain": chain, "start": START_NAMES[chain],
                 "sweeps": task["sweeps"], "proposals_per_sweep": n * (n - 1) // 2, "flag_every": FLAG_EVERY,
                 "spawn_key": [EXPERIMENT_ID, ei, n, bi, chain], "official": task["official"]},
        "commit": task["commit"], "dirty": task["dirty"], "clasificador_estado": clf["estado"],
        "initial": initial, "traces": traces, "counts": counts, "wall_s": wall, "init_s": t_init,
        "measure_s": t_meas, "s_drift_final": drift,
    }
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(_clean(rec)))
    tmp.replace(p)
    return f"{p.name} ok {sum(wall):.1f} s barridos + {t_meas:.1f} s medicion"


def load_part(out: Path, ei: int, n: int, bi: int, chain: int) -> dict[str, Any] | None:
    p = part_path(out, ei, n, bi, chain)
    if not p.exists():
        return None
    rec: dict[str, Any] = json.loads(p.read_text())
    return rec


def cmd_run(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    if args.sweeps_override is not None and not args.out:
        sys.exit("--sweeps-override exige --out (no se mezcla con la corrida oficial)")
    head, dirty = git_info()
    if dirty and not args.allow_dirty:
        sys.exit("exige arbol git limpio fuera de results/ (usar --allow-dirty solo para pruebas)")
    official = not dirty and args.sweeps_override is None and out == OUT
    base = {"out": str(out), "commit": head, "dirty": dirty, "official": official}
    if official:  # una corrida oficial no reutiliza clasificadores ni partes no oficiales
        stale = [p.name for p in out.glob("clasificador_N*.json") if not json.loads(p.read_text())["official"]]
        stale += [p.name for p in out.glob("parts/*.json") if not json.loads(p.read_text())["spec"]["official"]]
        if stale:
            sys.exit(f"{out} contiene archivos no oficiales ({len(stale)}, p. ej. {stale[0]}); moverlos antes")
    ns = [args.N] if args.N else list(NS)
    eis = [args.eps_idx] if args.eps_idx is not None else list(range(len(EPS)))
    tasks: list[dict[str, Any]] = []
    for n in ns:
        clf = ensure_classifier(out, n, base)  # antes de mirar cualquier cadena
        v = clf["validacion"]
        print(f"clasificador N={n}: pos {v['pos']:.2f}, kr {v['kr']:.2f}, bicapa {v['bicapa']:.2f} -> {clf['estado']}")
        if clf["estado"] != "VALIDADO":  # aclaracion A1 del preregistro: la mezcla no depende de la bandera
            print(f"N={n}: bandera INDETERMINADA; las cadenas corren solo para el diagnostico de mezcla")
        sw = args.sweeps_override if args.sweeps_override is not None else SWEEPS[n]
        for ei in eis:
            for bi in range(len(BS)):
                for ch in range(N_CHAINS):
                    tasks.append({**base, "eps_idx": ei, "N": n, "b_idx": bi, "chain": ch, "sweeps": sw})
    tasks = [t for t in tasks if not part_path(out, t["eps_idx"], t["N"], t["b_idx"], t["chain"]).exists()]
    tasks.sort(key=lambda t: -t["N"] * t["N"] * t["sweeps"])  # las mas largas primero
    print(f"{len(tasks)} partes por correr en {out} (commit {head[:8]}, dirty {dirty}, oficial {official})")
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


def first_exit(flags: list[Any]) -> int | None:
    """Primer barrido (1..T) desde el que la bandera queda en 0 durante >= STAY_UNFLAGGED mediciones seguidas.

    Las trazas guardan None en los barridos sin medicion de bandera; el barrido devuelto es el de la primera."""
    meas = [(t + 1, f) for t, f in enumerate(flags) if f is not None]
    run = 0
    for k, (_, f) in enumerate(meas):
        if f == 0:
            run += 1
            if run >= STAY_UNFLAGGED:
                return meas[k - STAY_UNFLAGGED + 1][0]
        else:
            run = 0
    return None


def chain_summary(p: dict[str, Any]) -> dict[str, Any]:
    tr = p["traces"]
    t = len(tr["S"])
    h = t // 2
    fl = [f for f in tr["flag"][h:] if f is not None]
    dmm = [x for x in tr["d_mm"][h:] if x is not None]
    nprop = p["spec"]["proposals_per_sweep"] * t
    return {
        "start": p["spec"]["start"], "sweeps": t, "flag_inicial": p["initial"]["flag"],
        "flag_frac_mitad_final": float(np.mean(fl)) if fl else None,
        "d_mm_mitad_final": float(np.mean(dmm)) if dmm else None,
        "d_mm_nan_frac": 1.0 - len(dmm) / max(t - h, 1),
        **{f"{k}_mitad_final": float(np.mean(tr[k][h:])) for k in ("S", "rel_frac", "height", "links_per_N")},
        "accept_rate": sum(p["counts"]["accepted"]) / nprop, "valid_frac": sum(p["counts"]["valid"]) / nprop,
        "accept_per_valid": sum(p["counts"]["accepted"]) / max(sum(p["counts"]["valid"]), 1),
        "s_per_sweep": float(np.mean(p["wall_s"])), "measure_s_per_sweep": p["measure_s"] / t,
        "s_drift_final": p["s_drift_final"],
    }


def cell_summary(parts: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {
        "chains": len(parts), "commits": sorted({p["commit"] for p in parts}),
        "dirty": any(p["dirty"] for p in parts), "official": all(p["spec"]["official"] for p in parts),
        "por_cadena": {str(p["spec"]["chain"]): chain_summary(p) for p in parts},
    }
    nprop = sum(p["spec"]["proposals_per_sweep"] * len(p["wall_s"]) for p in parts)
    out["accept_rate"] = sum(sum(p["counts"]["accepted"]) for p in parts) / nprop
    out["valid_frac"] = sum(sum(p["counts"]["valid"]) for p in parts) / nprop
    out["s_per_sweep"] = float(np.mean([w for p in parts for w in p["wall_s"]]))
    c2 = [p for p in parts if p["spec"]["chain"] == 2]
    if c2:
        out["sprinkling_4d_barrido_salida"] = first_exit(c2[0]["traces"]["flag"])
        out["sprinkling_4d_flag_inicial"] = c2[0]["initial"]["flag"]
    if len(parts) < 2:
        return out
    half = [len(p["traces"]["S"]) // 2 for p in parts]
    fin = {k: [np.asarray(p["traces"][k], dtype=np.float64)[h:] for p, h in zip(parts, half, strict=True)]
           for k in OBS_MIX}
    mix: dict[str, Any] = {}
    ok = True
    lost = False
    tau_max = 0.0
    for k in OBS_MIX:
        rh = split_rhat(fin[k])
        tw = [tau_sokal(x) for x in fin[k]]
        taus = [tt for tt, _ in tw]
        ess = float(sum(x.size / tt for x, tt in zip(fin[k], taus, strict=True)))
        mix[k] = {"rhat": rh, "tau": taus, "ventana_sokal": [w for _, w in tw], "ess": ess}
        ok = ok and math.isfinite(rh) and rh < RHAT_OK and ess >= ESS_OK
        lost = lost or not math.isfinite(rh) or rh >= RHAT_LOST
        tau_max = max(tau_max, *taus)
    out["mezcla_detalle"] = mix
    out["mezcla"] = ok
    if not ok:
        nc = len(parts)
        need = 2 * math.ceil(ESS_OK / nc * tau_max)  # mitad final con ESS_OK/nc tau por cadena; largo total 2x
        out["sweeps_needed"] = "no estimable" if lost else need
        out["horas_cpu_estimadas"] = None if lost else need * out["s_per_sweep"] * nc / 3600
    return out


def cmd_summarize(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    clfs: dict[str, Any] = {}
    for n in NS:
        if clf_path(out, n).exists():
            c = json.loads(clf_path(out, n).read_text())
            clfs[str(n)] = {k: c[k] for k in ("estado", "validacion", "delta4", "commit", "dirty", "official")}
    cells: dict[str, Any] = {}
    for ei, eps in enumerate(EPS):
        for n in NS:
            for bi, b in enumerate(BS):
                parts = [p for ch in range(N_CHAINS) if (p := load_part(out, ei, n, bi, ch)) is not None]
                if not parts:
                    continue
                key = f"e{ei}_N{n}_b{bi}"
                cells[key] = {"eps": eps, "N": n, "b": b, "beta": b / n, **cell_summary(parts)}
                if len(parts) < N_CHAINS:
                    cells[key]["incompleta"] = True
    res = {"spec": str(SPEC), "experiment_id": EXPERIMENT_ID,
           "criterio": {"rhat": RHAT_OK, "ess": ESS_OK, "rhat_no_estimable": RHAT_LOST, "observables": list(OBS_MIX),
                        "mitad": "final", "flag_every": FLAG_EVERY, "salida_sprinkling_mediciones": STAY_UNFLAGGED},
           "clasificador": clfs, "cells": cells}
    (out / "summary.json").write_text(json.dumps(_clean(res), indent=2, allow_nan=False))
    print(f"{'celda':<14}{'eps':>5}{'N':>5}{'b':>5} {'s/barr':>7} {'valid':>6} {'acept':>6} {'mezcla':>7} "
          f"{'Rhat S/rel/h':>20} {'ESS min':>8} {'salida':>7} flag(c0..c3)")
    for k, c in cells.items():
        mx = c.get("mezcla_detalle", {})
        rh = "/".join(f"{mx[o]['rhat']:.2f}" for o in OBS_MIX) if mx else "-"
        ess = min(mx[o]["ess"] for o in OBS_MIX) if mx else float("nan")
        fl = " ".join(f"{v['flag_frac_mitad_final']:.2f}" if v["flag_frac_mitad_final"] is not None else "-"
                      for v in c["por_cadena"].values())
        print(f"{k:<14}{c['eps']:>5}{c['N']:>5}{c['b']:>5} {c['s_per_sweep']:>7.3f} {c['valid_frac']:>6.3f} "
              f"{c['accept_rate']:>6.3f} {c.get('mezcla')!s:>7} {rh:>20} {ess:>8.1f} "
              f"{c.get('sprinkling_4d_barrido_salida')!s:>7} {fl}")


# --- figura -----------------------------------------------------------------------------------------


def cmd_figure(args: argparse.Namespace) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(args.out) if args.out else OUT
    fig_path = Path(args.fig) if args.fig else (out / "esc0_piloto.png" if args.out else FIG)
    rows = [(ei, n) for ei in range(len(EPS)) for n in NS
            if any(part_path(out, ei, n, bi, ch).exists() for bi in range(len(BS)) for ch in range(N_CHAINS))]
    if not rows:
        sys.exit("no hay partes")
    colors = ("#1f77b4", "#d62728", "#2ca02c", "#9467bd")
    fig, axes = plt.subplots(len(rows), len(BS), figsize=(15, 2.25 * len(rows) + 0.9), squeeze=False,
                             sharey="row")
    for r, (ei, n) in enumerate(rows):
        for bi, b in enumerate(BS):
            ax = axes[r, bi]
            for ch in range(N_CHAINS):
                p = load_part(out, ei, n, bi, ch)
                if p is not None:
                    y = [p["initial"]["rel_frac"], *p["traces"]["rel_frac"]]  # barrido 0 = inicio
                    ax.plot(np.arange(len(y)), y, lw=0.7, color=colors[ch], label=START_NAMES[ch])
            ax.set_ylim(-0.03, 1.03)
            ax.tick_params(labelsize=7)
            bs = f"{b:+d}" if b else "0"
            ax.set_title(f"ε = {EPS[ei]:g}, N = {n}, b = {bs}".replace("-", "−"), fontsize=8)
            if bi == 0:
                ax.set_ylabel("fracción de relaciones", fontsize=7)
            if r == len(rows) - 1:
                ax.set_xlabel("barrido", fontsize=7)
    h, lab = axes[0, 0].get_legend_handles_labels()
    leg = fig.legend(h, lab, loc="lower center", ncol=4, fontsize=8, frameon=False)
    for line in leg.get_lines():
        line.set_linewidth(2.0)
    fig.suptitle("ES-C0: fracción de relaciones de las 4 cadenas por celda (diagnóstico de mezcla, no es evidencia)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0.035, 1, 0.975))
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=110)
    print(fig_path)


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--N", type=int, choices=NS)
    r.add_argument("--eps-idx", type=int, choices=range(len(EPS)))
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--sweeps-override", type=int)
    r.add_argument("--out")
    r.add_argument("--allow-dirty", action="store_true")
    s = sub.add_parser("summarize")
    s.add_argument("--out")
    f = sub.add_parser("figure")
    f.add_argument("--out")
    f.add_argument("--fig")
    args = ap.parse_args()
    {"run": cmd_run, "summarize": cmd_summarize, "figure": cmd_figure}[args.cmd](args)


if __name__ == "__main__":
    main()
