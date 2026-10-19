"""R17-3 (docs/OMEGA_CONSEJO_REV17_ES0.md §5.2): perfil de s_sj, exactitud y tiempos de omega.causality.sj_fast.

  python tools/r17_3_bench.py profile [--ns 64,128,256,512]
  python tools/r17_3_bench.py exact --space Ea|Eb [--ns 32,64,128,256,512] [--states 2000]
  python tools/r17_3_bench.py summarize

Exactitud: cadenas de Metropolis en E-a y E-b (mismos movimientos que mcmc_sj) con el evaluador incremental; en
CADA propuesta valida (aceptada o rechazada) se compara con el calculo completo de s_sj (mismas operaciones que
mcmc_sj.s_sj, ademas de los margenes a los umbrales). Arranques como ES-0 (tools/es0_piloto.py); lam en
{0, -L, +L} con L = 2/sigma0 de ES-0 a N = 64 (solo valores de prueba). Semillas:
SeedSequence(20261102, spawn_key=(327, espacio, N, indice de lam, arranque)). Salida en --out (JSON).
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # como ES-0 (un hilo por cadena)

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omega.causality.entropy import interval_balanced, pauli_jordan, sj_from_eig, ssee
from omega.causality.mcmc import order_from_v
from omega.causality.mcmc_sj import s_sj
from omega.causality.sj_fast import SJEvaluator, s_sj_fast

ENTROPY = 20261102
EXPERIMENT_ID = 327
SPACES = {"Ea": 0, "Eb": 1}
LAM_TEST = {"Ea": 2 / 0.2876803676752637, "Eb": 2 / 0.11037570105977527}  # 2/sigma0 de ES-0, N = 64
OUT = Path("/tmp/claude-0/-home-user-FichaClaraV2/90c7523c-b9d2-5eee-afba-260b3848fae0/scratchpad/r17_3")
Bool = npt.NDArray[np.bool_]


def rng_for(space: int, n: int, lam_idx: int, start: int) -> np.random.Generator:
    ss = np.random.SeedSequence(ENTROPY, spawn_key=(EXPERIMENT_ID, space, n, lam_idx, start))
    return np.random.Generator(np.random.PCG64(ss))


def kr3(n: int, rng: np.random.Generator) -> Bool:
    a = n // 4
    l1, l2, l3 = np.arange(a), np.arange(a, n - a), np.arange(n - a, n)
    c = np.zeros((n, n), dtype=bool)
    c[np.ix_(l1, l3)] = True
    c[np.ix_(l1, l2)] = rng.random((l1.size, l2.size)) < 0.5
    c[np.ix_(l2, l3)] = rng.random((l2.size, l3.size)) < 0.5
    return c


def reference(c: Bool) -> tuple[float, int, float, float]:
    """Mismas operaciones que mcmc_sj.s_sj, mas los margenes a los umbrales global y de la region."""
    n = c.shape[0]
    r = interval_balanced(c, n // 8)
    if r.size == 0:
        return 0.0, 0, math.inf, math.inf
    cut = math.sqrt(n) / (4 * math.pi)
    mu, u = np.linalg.eigh(pauli_jordan(c))
    idt, w = sj_from_eig(mu, u, cut)
    cut_r = math.sqrt(r.size) / (4 * math.pi)
    s = ssee(idt, w, r, True, cut_r)
    nu = np.linalg.eigvalsh(idt[np.ix_(r, r)])
    return s, int(r.size), float(np.min(np.abs(np.abs(mu) - cut))), float(np.min(np.abs(np.abs(nu) - cut_r)))


def _tm(f: Any, reps: int) -> float:
    f()
    t0 = time.perf_counter()
    for _ in range(reps):
        f()
    return (time.perf_counter() - t0) / reps


# --- perfil -----------------------------------------------------------------------------------------------


def _profile_one(n: int) -> dict[str, Any]:
    rng = rng_for(9, n, 0, 0)
    c = order_from_v(rng.permutation(n).astype(np.int64))
    reps = max(3, 2000 // n)
    cut = math.sqrt(n) / (4 * math.pi)
    r = interval_balanced(c, n // 8)
    h = pauli_jordan(c)
    mu, u = np.linalg.eigh(h)
    idt, w = sj_from_eig(mu, u, cut)
    k = c.astype(np.float64) - c.T
    m = k.T @ k
    ur = np.random.default_rng(0).random((n, n))
    row: dict[str, Any] = {
        "N": n,
        "nR": int(r.size),
        "kept": int(np.sum(np.abs(mu) >= cut)),
        "interval_balanced_ms": 1e3 * _tm(lambda: interval_balanced(c, n // 8), reps),
        "pauli_jordan_ms": 1e3 * _tm(lambda: pauli_jordan(c), reps),
        "eigh_complex_ms": 1e3 * _tm(lambda: np.linalg.eigh(h), reps),
        "sj_from_eig_ms": 1e3 * _tm(lambda: sj_from_eig(mu, u, cut), reps),
        "ssee_ms": 1e3 * _tm(lambda: ssee(idt, w, r, True, math.sqrt(r.size) / (4 * math.pi)), reps),
        "s_sj_total_ms": 1e3 * _tm(lambda: s_sj(c), reps),
        "s_sj_fast_ms": 1e3 * _tm(lambda: s_sj_fast(c), reps),
        "eigh_real_M_ms": 1e3 * _tm(lambda: np.linalg.eigh(m), reps),
        "eigvalsh_real_M_ms": 1e3 * _tm(lambda: np.linalg.eigvalsh(m), reps),
        # cota inferior del costo de aceptar en una actualizacion de bajo rango: U (complejo) por U~ (real), rango 1
        "lowrank_commit_one_rank1_ms": 1e3 * _tm(lambda: u @ ur, reps),
    }
    ev = SJEvaluator.from_two_order(rng.permutation(n).astype(np.int64))
    pairs = rng.choice(n, size=(reps + 1, 2), replace=True)
    it = [0]

    def swap_once() -> None:
        a, b = (int(x) for x in pairs[it[0] % (reps + 1)])
        it[0] += 1
        ev.propose_swap(a, b if a != b else (a + 1) % n)
        ev.rollback()

    row["evaluator_swap_ms"] = 1e3 * _tm(swap_once, reps)
    return row


def profile(ns: list[int], out: Path) -> None:
    rows = []
    for n in ns:
        row = _profile_one(n)
        print(json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
        rows.append(row)
    out.mkdir(parents=True, exist_ok=True)
    (out / "profile.json").write_text(json.dumps(rows, indent=1))


# --- exactitud y tiempos a lo largo de cadenas ---------------------------------------------------------------


def _accept(lam: float, ds: float, u: float) -> bool:
    x = lam * ds
    return x >= 0 or u < math.exp(x)


def start(space: str, n: int, idx: int, rng: np.random.Generator) -> SJEvaluator:
    if space == "Ea":
        if idx == 0:
            v = np.arange(n, dtype=np.int64)
        elif idx == 1:
            v = np.arange(n - 1, -1, -1, dtype=np.int64)
        else:
            v = rng.permutation(n).astype(np.int64)
        return SJEvaluator.from_two_order(v)
    if idx == 0:
        c = np.zeros((n, n), dtype=bool)
    elif idx == 1:
        c = np.triu(np.ones((n, n), dtype=bool), k=1)
    elif idx == 2:
        c = order_from_v(rng.permutation(n).astype(np.int64))
    else:
        c = kr3(n, rng)
    return SJEvaluator(c)


def run_chain(space: str, n: int, lam_idx: int, idx: int, n_checked: int, burn: int) -> dict[str, Any]:
    rng = rng_for(SPACES[space], n, lam_idx, idx)
    lam = (0.0, -LAM_TEST[space], LAM_TEST[space])[lam_idx]
    ev = start(space, n, idx, rng)
    c = ev.c
    st: dict[str, Any] = {"space": space, "N": n, "lam": lam, "start": idx, "checked": 0, "proposals": 0,
                          "valid": 0, "accepted": 0, "max_dS": 0.0, "nR_mismatch": 0, "empty_region": 0,
                          "min_margin_ref_global": math.inf, "min_margin_ref_region": math.inf,
                          "t_fast": 0.0, "t_full": 0.0, "t_check": 0.0}
    t = 0
    while st["checked"] < n_checked:
        a = int(rng.integers(0, n))
        b = int(rng.integers(0, n - 1))
        b += b >= a
        if space == "Eb" and a > b:
            a, b = b, a
        u = float(rng.random())
        t += 1
        measured = t > burn
        t0 = time.perf_counter()
        if space == "Eb":
            if c[a, b]:
                ok = not np.any(c[a] & c[:, b])
            else:
                ok = not (np.any(c[b] & ~c[a]) or np.any(c[:, a] & ~c[:, b]))
            if measured:
                st["t_check"] += time.perf_counter() - t0
                st["proposals"] += 1
            if not ok:
                continue
        elif measured:
            st["proposals"] += 1
        s_old = ev.s
        t0 = time.perf_counter()
        r = ev.propose_toggle(a, b) if space == "Eb" else ev.propose_swap(a, b)
        t1 = time.perf_counter()
        if measured:
            st["valid"] += 1
            ref = reference(ev.c)
            t2 = time.perf_counter()
            st["t_full"] += t2 - t1
            st["checked"] += 1
            st["max_dS"] = max(st["max_dS"], abs(r.s - ref[0]))
            st["nR_mismatch"] += int(r.n_r != ref[1])
            st["empty_region"] += int(ref[1] == 0)
            st["min_margin_ref_global"] = min(st["min_margin_ref_global"], ref[2])
            st["min_margin_ref_region"] = min(st["min_margin_ref_region"], ref[3])
        t3 = time.perf_counter()
        if _accept(lam, r.s - s_old, u):
            ev.commit()
            st["accepted"] += int(measured)
        else:
            ev.rollback()
        if measured:
            st["t_fast"] += (t1 - t0) + (time.perf_counter() - t3)
    st["state_ok"] = ev.check_state()
    st["fallbacks"] = ev.n_fallback
    st["min_margin_fast_global"] = ev.min_margin_global
    st["min_margin_fast_region"] = ev.min_margin_region
    return st


def exact(space: str, ns: list[int], n_states: int, out: Path) -> None:
    n_starts = 4
    per = math.ceil(n_states / (n_starts * 3))
    res = []
    for n in ns:
        burn = {32: 1000, 64: 1000, 128: 600, 256: 300}.get(n, 100)
        for lam_idx in range(3):
            for idx in range(n_starts):
                st = run_chain(space, n, lam_idx, idx, per, burn)
                print(json.dumps({k: (f"{v:.3g}" if isinstance(v, float) else v) for k, v in st.items()}), flush=True)
                res.append(st)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"exact_{space}.json").write_text(json.dumps(res, indent=1))


def summarize(out: Path) -> None:
    summ: dict[str, Any] = {}
    for space in SPACES:
        p = out / f"exact_{space}.json"
        if not p.exists():
            continue
        rows = json.loads(p.read_text())
        for n in sorted({r["N"] for r in rows}):
            rs = [r for r in rows if r["N"] == n]
            chk = sum(r["checked"] for r in rs)
            prop = sum(r["proposals"] for r in rs)
            val = sum(r["valid"] for r in rs)
            t_fast_eval = sum(r["t_fast"] for r in rs) / val
            t_full_eval = sum(r["t_full"] for r in rs) / val
            t_check = sum(r["t_check"] for r in rs) / prop
            vf = val / prop
            m = n * (n - 1) // 2
            summ[f"{space}_N{n}"] = {
                "states_checked": chk,
                "max_dS": max(r["max_dS"] for r in rs),
                "nR_mismatch": sum(r["nR_mismatch"] for r in rs),
                "empty_region_states": sum(r["empty_region"] for r in rs),
                "fallbacks": sum(r["fallbacks"] for r in rs),
                "state_ok": all(r["state_ok"] for r in rs),
                "min_margin_ref_global": min(r["min_margin_ref_global"] for r in rs),
                "min_margin_ref_region": min(r["min_margin_ref_region"] for r in rs),
                "min_margin_fast_global": min(r["min_margin_fast_global"] for r in rs),
                "valid_frac": vf,
                "accept_per_valid": sum(r["accepted"] for r in rs) / val,
                "ms_per_valid_full": 1e3 * t_full_eval,
                "ms_per_valid_fast": 1e3 * t_fast_eval,
                "speedup": t_full_eval / t_fast_eval,
                "s_per_sweep_full": m * (t_check + vf * t_full_eval),
                "s_per_sweep_fast": m * (t_check + vf * t_fast_eval),
            }
    (out / "summary.json").write_text(json.dumps(summ, indent=1))
    for k, v in summ.items():
        print(k, json.dumps({a: (f"{b:.3g}" if isinstance(b, float) else b) for a, b in v.items()}))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("profile", "exact", "summarize"))
    ap.add_argument("--ns", default="")
    ap.add_argument("--space", choices=tuple(SPACES), default="Ea")
    ap.add_argument("--states", type=int, default=2000)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()
    if a.cmd == "profile":
        profile([int(x) for x in (a.ns or "64,128,256,512").split(",")], a.out)
    elif a.cmd == "exact":
        exact(a.space, [int(x) for x in (a.ns or "32,64,128,256,512").split(",")], a.states, a.out)
    else:
        summarize(a.out)


if __name__ == "__main__":
    main()
