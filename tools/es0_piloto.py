"""ES-0: piloto de factibilidad del ensamble entropico de Omega-E (docs/OMEGA_ES0_PILOTO.md).

Ejecutor reanudable por partes (una parte = una cadena):
  run --stage 1|2|cost [--N 32] [--space Ea] [--workers 4] [--sweeps-override k --out DIR --allow-dirty]
  summarize [--out DIR]   -> <out>/summary.json
  figure [--out DIR] [--fig PNG]

Semillas: SeedSequence(20261100, spawn_key=(325, espacio, N, indice de lam, cadena)); espacio Ea = 0, Eb = 1;
indice de lam: 0 -> lam = 0, 1 -> -2/sigma0, 2 -> +2/sigma0. sigma0 = desviacion estandar (ddof = 1) de S_SJ en la
mitad final de las 4 cadenas de la etapa 1 del mismo espacio y N, juntas. Si sigma0 == 0 se registra y no hay etapa 2.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # matrices chicas; evita sobresuscripcion con --workers

import argparse
import json
import math
import multiprocessing as mp
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omega.causality.mcmc import order_from_v
from omega.causality.mcmc_sj import (
    GeneralOrderState,
    TwoOrderState,
    init_general,
    init_two_order,
    measure,
    sweep_general,
    sweep_two_order,
)

EXPERIMENT_ID = 325
ENTROPY = 20261100
SPACES = {"Ea": 0, "Eb": 1}
NS = (32, 64)
SWEEPS = {32: 400, 64: 200}
N_COST, SWEEPS_COST = 128, 5
N_CHAINS = 4
LAM_UNIT = 2.0
RHAT_OK, ESS_OK, RHAT_LOST, SOKAL_C = 1.05, 200.0, 1.5, 5.0
OBS_MIX = ("S", "rel_frac", "height")
TRACE_KEYS = ("S", "nR", "has_region", "rel_frac", "links", "height")
SPEC = Path("docs") / "OMEGA_ES0_PILOTO.md"
OUT = Path("results") / "omega_es0"
FIG = Path("docs") / "figuras" / "es0_piloto.png"
START_NAMES = {
    "Ea": ("cadena total", "anticadena", "2-orden al azar", "2-orden al azar"),
    "Eb": ("orden vacio", "cadena total", "sprinkling 2D", "KR 3 capas"),
}
Bool = npt.NDArray[np.bool_]


def rng_for(space: str, n: int, lam_idx: int, chain: int) -> np.random.Generator:
    key = (EXPERIMENT_ID, SPACES[space], n, lam_idx, chain)
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(ENTROPY, spawn_key=key)))


def git_info() -> tuple[str, bool]:
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=False).stdout.strip()
    st = subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True, check=False)
    return head, bool(st.stdout.strip())


def part_path(out: Path, space: str, n: int, lam_idx: int, chain: int) -> Path:
    return out / "parts" / f"{space}_N{n}_l{lam_idx}_c{chain}.json"


# --- inicios (§2) -----------------------------------------------------------------------------------


def kr3(n: int, rng: np.random.Generator) -> Bool:
    """KR de 3 capas en orden de etiquetas (N/4, N - 2(N/4), N/4): capas 1 y 3 todas; contiguas con p = 1/2."""
    a = n // 4
    l1, l2, l3 = np.arange(a), np.arange(a, n - a), np.arange(n - a, n)
    c = np.zeros((n, n), dtype=bool)
    c[np.ix_(l1, l3)] = True
    c[np.ix_(l1, l2)] = rng.random((l1.size, l2.size)) < 0.5
    c[np.ix_(l2, l3)] = rng.random((l2.size, l3.size)) < 0.5
    return c


def start_state(space: str, n: int, chain: int, cost: bool, rng: np.random.Generator) -> TwoOrderState | GeneralOrderState:
    if space == "Ea":
        if cost or chain >= 2:
            v = rng.permutation(n).astype(np.int64)
        elif chain == 0:
            v = np.arange(n, dtype=np.int64)
        else:
            v = np.arange(n - 1, -1, -1, dtype=np.int64)
        return init_two_order(v)
    if cost or chain == 2:
        c = order_from_v(rng.permutation(n).astype(np.int64))  # 2-orden con u = identidad: ya etiquetado natural
    elif chain == 0:
        c = np.zeros((n, n), dtype=bool)
    elif chain == 1:
        c = np.triu(np.ones((n, n), dtype=bool), k=1)
    else:
        c = kr3(n, rng)
    return init_general(c)


# --- una cadena -------------------------------------------------------------------------------------


def run_part(task: dict[str, Any]) -> str:
    space, n, lam_idx, chain = task["space"], task["N"], task["lam_idx"], task["chain"]
    out = Path(task["out"])
    p = part_path(out, space, n, lam_idx, chain)
    if p.exists():
        return f"{p.name} ya existe"
    rng = rng_for(space, n, lam_idx, chain)
    t0 = time.perf_counter()
    st = start_state(space, n, chain, task["cost"], rng)
    t_init = time.perf_counter() - t0
    lam = float(task["lam"])
    traces: dict[str, list[Any]] = {k: [] for k in TRACE_KEYS}
    counts: dict[str, list[int]] = {"proposals": [], "valid": [], "accepted": []}
    wall: list[float] = []
    initial = measure(st)
    for _ in range(task["sweeps"]):
        ts = time.perf_counter()
        if isinstance(st, TwoOrderState):
            cnt = sweep_two_order(st, lam, rng)
        else:
            cnt = sweep_general(st, lam, rng)
        wall.append(time.perf_counter() - ts)
        for k, x in cnt._asdict().items():
            counts[k].append(int(x))
        row = measure(st)
        for k in TRACE_KEYS:
            traces[k].append(row[k])
    rec = {
        "spec": {"doc": str(SPEC), "experiment_id": EXPERIMENT_ID, "entropy": ENTROPY, "space": space, "N": n,
                 "lam_idx": lam_idx, "chain": chain, "sweeps": task["sweeps"], "cost_cell": task["cost"],
                 "start": "2-orden al azar" if task["cost"] else START_NAMES[space][chain],
                 "spawn_key": [EXPERIMENT_ID, SPACES[space], n, lam_idx, chain], "official": task["official"]},
        "commit": task["commit"], "dirty": task["dirty"], "lam": lam, "sigma0": task.get("sigma0"),
        "initial": initial, "traces": traces, "counts": counts, "wall_s": wall, "init_s": t_init,
        "n_eval": st.n_eval, "t_eval_s": st.t_eval,
    }
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(rec))
    tmp.replace(p)
    return f"{p.name} ok {sum(wall):.1f} s ({1e3 * st.t_eval / max(st.n_eval, 1):.2f} ms/eval)"


def load_part(out: Path, space: str, n: int, lam_idx: int, chain: int) -> dict[str, Any] | None:
    p = part_path(out, space, n, lam_idx, chain)
    return json.loads(p.read_text()) if p.exists() else None


def sigma0_for(out: Path, space: str, n: int) -> float | None:
    """sigma0 de la etapa 1 (None si faltan partes); se registra en <out>/sigma0_<espacio>_N<N>.json."""
    parts = [load_part(out, space, n, 0, ch) for ch in range(N_CHAINS)]
    if any(p is None for p in parts):
        return None
    s = np.concatenate([np.asarray(p["traces"]["S"], dtype=np.float64)[len(p["traces"]["S"]) // 2:]
                        for p in parts if p is not None])
    sig = float(np.std(s, ddof=1)) if s.size > 1 else 0.0
    (out / f"sigma0_{space}_N{n}.json").write_text(json.dumps(
        {"space": space, "N": n, "sigma0": sig, "n_samples": int(s.size),
         "etapa2": "omitida (sigma0 = 0)" if sig == 0 else "lam = -/+ 2/sigma0"}))
    return sig


def cmd_run(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    if args.sweeps_override is not None and not args.out:
        sys.exit("--sweeps-override exige --out (no se mezcla con la corrida oficial)")
    head, dirty = git_info()
    if dirty and not args.allow_dirty:
        sys.exit("exige arbol git limpio (usar --allow-dirty solo para pruebas)")
    official = not dirty and args.sweeps_override is None and out == OUT
    spaces = [args.space] if args.space else list(SPACES)
    tasks: list[dict[str, Any]] = []
    base = {"out": str(out), "commit": head, "dirty": dirty, "official": official}
    if args.stage == "cost":
        sw = args.sweeps_override if args.sweeps_override is not None else SWEEPS_COST
        for sp in spaces:
            tasks.append({**base, "space": sp, "N": N_COST, "lam_idx": 0, "chain": 0, "lam": 0.0, "sweeps": sw,
                          "cost": True})
    else:
        ns = [args.N] if args.N else list(NS)
        for sp in spaces:
            for n in ns:
                sw = args.sweeps_override if args.sweeps_override is not None else SWEEPS[n]
                if args.stage == "1":
                    lams: list[tuple[int, float, float | None]] = [(0, 0.0, None)]
                else:
                    sig = sigma0_for(out, sp, n)
                    if sig is None:
                        print(f"{sp} N={n}: falta la etapa 1, se omite")
                        continue
                    if sig == 0:
                        print(f"{sp} N={n}: sigma0 = 0, sin etapa 2")
                        continue
                    lams = [(1, -LAM_UNIT / sig, sig), (2, LAM_UNIT / sig, sig)]
                for li, lam, sig0 in lams:
                    for ch in range(N_CHAINS):
                        tasks.append({**base, "space": sp, "N": n, "lam_idx": li, "chain": ch, "lam": lam,
                                      "sigma0": sig0, "sweeps": sw, "cost": False})
    tasks = [t for t in tasks if not part_path(out, t["space"], t["N"], t["lam_idx"], t["chain"]).exists()]
    print(f"{len(tasks)} partes por correr en {out}")
    if args.workers <= 1:
        for t in tasks:
            print(run_part(t), flush=True)
        return
    with mp.get_context("spawn").Pool(args.workers) as pool:
        for msg in pool.imap_unordered(run_part, tasks):
            print(msg, flush=True)


# --- diagnosticos de mezcla -------------------------------------------------------------------------


def split_rhat(chains: list[npt.NDArray[np.float64]]) -> float:
    """R-hat dividido clasico (Gelman-Rubin): cada cadena se parte en dos mitades."""
    seqs = []
    for x in chains:
        h = x.size // 2
        seqs += [x[:h], x[h:2 * h]]
    m, n = len(seqs), min(s.size for s in seqs)
    if n < 2:
        return float("nan")
    arr = np.array([s[:n] for s in seqs])
    w = float(np.mean(arr.var(axis=1, ddof=1)))
    b = n * float(arr.mean(axis=1).var(ddof=1)) if m > 1 else 0.0
    if w == 0:
        return 1.0 if b == 0 else float("inf")
    var = (n - 1) / n * w + b / n
    return math.sqrt(var / w)


def tau_sokal(x: npt.NDArray[np.float64], c: float = SOKAL_C) -> tuple[float, bool]:
    """Tiempo de autocorrelacion integrado con ventana automatica de Sokal (M >= c tau). (tau, ventana hallada).

    Serie constante: tau = n (conservador: aporta 1 al ESS)."""
    n = x.size
    y = x - x.mean()
    if n < 2 or not np.any(y):
        return float(n), False
    f = np.fft.rfft(y, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    rho = acf / acf[0]
    taus = 2.0 * np.cumsum(rho) - 1.0  # tau(M) = 1 + 2 sum_{t=1..M} rho(t)
    for m in range(1, n):
        if m >= c * taus[m]:
            return max(float(taus[m]), 1.0), True
    return max(float(taus[-1]), 1.0), False


def eta2(s: npt.NDArray[np.float64], g: npt.NDArray[np.int64]) -> float:
    """Fraccion de la varianza de S explicada por |R| (medias por grupo)."""
    tot = float(np.sum((s - s.mean()) ** 2))
    if tot == 0:
        return float("nan")
    within = sum(float(np.sum((s[g == k] - s[g == k].mean()) ** 2)) for k in np.unique(g))
    return 1.0 - within / tot


def cell_summary(parts: list[dict[str, Any]]) -> dict[str, Any]:
    t_eval = sum(p["t_eval_s"] for p in parts)
    n_eval = sum(p["n_eval"] for p in parts)
    tot = {k: sum(sum(p["counts"][k]) for p in parts) for k in ("proposals", "valid", "accepted")}
    s_sweep = float(np.mean([w for p in parts for w in p["wall_s"]]))
    half = [len(p["traces"]["S"]) // 2 for p in parts]
    fin = {k: [np.asarray(p["traces"][k], dtype=np.float64)[h:] for p, h in zip(parts, half)] for k in TRACE_KEYS}
    allv = {k: np.concatenate(fin[k]) for k in TRACE_KEYS}
    nr = allv["nR"].astype(np.int64)
    vals, cnts = np.unique(nr, return_counts=True)
    out: dict[str, Any] = {
        "chains": len(parts), "sweeps": [len(p["wall_s"]) for p in parts], "lam": parts[0]["lam"],
        "sigma0": parts[0].get("sigma0"), "commits": sorted({p["commit"] for p in parts}),
        "dirty": any(p["dirty"] for p in parts), "official": all(p["spec"]["official"] for p in parts),
        "ms_per_eval": 1e3 * t_eval / n_eval if n_eval else None, "s_per_sweep": s_sweep,
        "valid_frac": tot["valid"] / tot["proposals"] if tot["proposals"] else None,
        "accept_rate": tot["accepted"] / tot["proposals"] if tot["proposals"] else None,
        "accept_per_valid": tot["accepted"] / tot["valid"] if tot["valid"] else None,
        "frac_sin_region_mitad_final": float(1.0 - allv["has_region"].mean()) if nr.size else None,
        "frac_sin_region_todo": float(1.0 - np.mean([x for p in parts for x in p["traces"]["has_region"]])),
        "nR_mitad_final": {"mean": float(nr.mean()) if nr.size else None,
                           "hist": {str(int(v)): int(k) for v, k in zip(vals, cnts)}},
        "S_mitad_final": {"mean": float(allv["S"].mean()) if nr.size else None,
                          "std": float(allv["S"].std(ddof=1)) if nr.size > 1 else None},
        "rel_frac_mitad_final_mean": float(allv["rel_frac"].mean()) if nr.size else None,
    }
    reg = nr > 0
    s_reg = allv["S"][reg]
    out["S_vs_nR"] = {
        "pearson_con_region": float(np.corrcoef(s_reg, nr[reg])[0, 1])
        if reg.sum() > 2 and np.std(s_reg) > 0 and np.std(nr[reg]) > 0 else None,
        "eta2_todo": eta2(allv["S"], nr) if nr.size > 1 else None,
        "eta2_con_region": eta2(s_reg, nr[reg]) if reg.sum() > 1 else None,
    }
    if len(parts) < 2:
        return out
    mix: dict[str, Any] = {}
    ok = True
    lost = False
    tau_max = 0.0
    for k in OBS_MIX:
        rh = split_rhat(fin[k])
        tw = [tau_sokal(x) for x in fin[k]]
        taus = [t for t, _ in tw]
        ess = float(sum(x.size / t for x, t in zip(fin[k], taus)))
        mix[k] = {"rhat": rh, "tau": taus, "ventana_sokal": [w for _, w in tw], "ess": ess}
        ok = ok and math.isfinite(rh) and rh < RHAT_OK and ess >= ESS_OK
        lost = lost or not math.isfinite(rh) or rh >= RHAT_LOST
        tau_max = max(tau_max, *taus)
    out["mezcla_detalle"] = mix
    out["mezcla"] = ok
    if not ok:
        nc = len(parts)
        need = 2 * math.ceil(ESS_OK / nc * tau_max)  # mitad final de 50 tau por cadena (4 cadenas); largo total 2x
        out["sweeps_needed"] = "no estimable" if lost else need
        out["horas_cpu_estimadas"] = None if lost else need * s_sweep * nc / 3600
    return out


def _clean(x: Any) -> Any:
    """JSON estricto: NaN -> None, infinito -> "inf"."""
    if isinstance(x, float) and not math.isfinite(x):
        return None if math.isnan(x) else ("inf" if x > 0 else "-inf")
    if isinstance(x, dict):
        return {k: _clean(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_clean(v) for v in x]
    return x


def cmd_summarize(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    cells: dict[str, Any] = {}
    for sp in SPACES:
        for n in (*NS, N_COST):
            for li in (0, 1, 2):
                parts = [p for ch in range(N_CHAINS) if (p := load_part(out, sp, n, li, ch)) is not None]
                if not parts:
                    continue
                key = f"{sp}_N{n}_l{li}"
                cells[key] = {"space": sp, "N": n, "lam_idx": li, **cell_summary(parts)}
                if n != N_COST and len(parts) < N_CHAINS:
                    cells[key]["incompleta"] = True
    sig = {p.stem: json.loads(p.read_text()) for p in sorted(out.glob("sigma0_*.json"))}
    res = {"spec": str(SPEC), "criterio": {"rhat": RHAT_OK, "ess": ESS_OK, "rhat_no_estimable": RHAT_LOST,
                                           "sokal_c": SOKAL_C, "observables": list(OBS_MIX)},
           "sigma0": sig, "cells": cells}
    (out / "summary.json").write_text(json.dumps(_clean(res), indent=2, allow_nan=False))
    for k, c in cells.items():
        print(f"{k}: {c['ms_per_eval'] or 0:.2f} ms/eval, {c['s_per_sweep']:.2f} s/barrido, "
              f"validas {c['valid_frac']:.3f}, acept {c['accept_rate']:.3f}, mezcla {c.get('mezcla')}")


def cmd_figure(args: argparse.Namespace) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(args.out) if args.out else OUT
    fig_path = Path(args.fig) if args.fig else (out / "es0_piloto.png" if args.out else FIG)
    cells = []
    for li in (0, 1, 2):
        for sp in SPACES:
            for n in NS:
                parts = [load_part(out, sp, n, li, ch) for ch in range(N_CHAINS)]
                if any(p is not None for p in parts):
                    cells.append((sp, n, li, parts))
    if not cells:
        sys.exit("no hay partes")
    lam_name = {0: "λ = 0", 1: "λ = −2/σ₀", 2: "λ = +2/σ₀"}
    fig, axes = plt.subplots(len(cells), 2, figsize=(11, 2.3 * len(cells)), squeeze=False)
    for row, (sp, n, li, parts) in enumerate(cells):
        for col, (k, lab) in enumerate((("S", "S_SJ"), ("rel_frac", "fracción de relaciones"))):
            ax = axes[row, col]
            for ch, p in enumerate(parts):
                if p is not None:
                    ax.plot(np.arange(1, len(p["traces"][k]) + 1), p["traces"][k], lw=0.8,
                            label=START_NAMES[sp][ch])
            ax.set_title(f"{sp}, N = {n}, {lam_name[li]}", fontsize=9)
            ax.set_ylabel(lab, fontsize=8)
            ax.tick_params(labelsize=7)
            if row == len(cells) - 1:
                ax.set_xlabel("barrido", fontsize=8)
            if col == 1:  # los inicios difieren entre Ea y Eb
                ax.legend(fontsize=6, loc="best")
    fig.suptitle("ES-0: trazas de las 4 cadenas por celda (piloto, no es evidencia)", fontsize=10)
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=110)
    print(fig_path)


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--stage", choices=("1", "2", "cost"), required=True)
    r.add_argument("--N", type=int, choices=NS)
    r.add_argument("--space", choices=tuple(SPACES))
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
