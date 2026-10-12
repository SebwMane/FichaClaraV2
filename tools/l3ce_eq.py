"""L-3c-EQ (Rev. 3 §16): transicion de equilibrio por Wang-Landau + multicanonica. Constantes congeladas.

Tareas:
  --task validate                       (N = 8, enumeracion exacta)
  --task wl   --n {32,64} --chunk K     (continua desde el punto de control; K multiplo de BATCH)
  --task prod --n {32,64}               (un bloque de produccion por llamada; 4 bloques)
  --task final
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from l3ce_prime import STEP, classifier, flag, git_state  # noqa: E402

from omega.causality.wanglandau import BATCH, WLState, exact_dos, new_state, normalized_lng, prod_run, wl_run  # noqa: E402

EXPERIMENT_ID = 306
ENT_WL, ENT_PROD, ENT_VAL = 20261032, 20261033, 20261035  # EQ-A2: entropia nueva de validacion (20261034 fallo)
SIGMA0 = {32: 31, 64: 54}
WL_BUDGET = {32: 6_000_000, 64: 24_000_000}
PROD_BUDGET = {32: 2_000_000, 64: 6_000_000}
PROD_EVERY, BLOCKS = 20, 4
VAL_N, VAL_STEPS, VAL_SMAX_EXTRA, VAL_TOL, VAL_MIN_STATES = 8, 2_000_000, 40, 0.15, 10  # EQ-A1 smax 48; EQ-A2 2e6 pasos
B_SCAN = np.round(np.arange(0.5, 12.0001, 0.01), 2)
OUT = Path("results") / "l3ce_eq"


def seeded(entropy: int, n: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, n))))


def window(n: int) -> tuple[int, int]:
    # AUDITORIA A-20 / enmienda EQ-A1: la bicapa NO es el minimo; ventana inferior en -N^2/2 - 4N
    smin = -(n * n // 2) - 4 * n
    smax = n + 2 * ((4 * SIGMA0[n]) // 2)
    return smin, smax


def ckpt(n: int) -> Path:
    return OUT / f"state_{n}.json"


def log_chunk(name: str, head: str, info: dict[str, Any]) -> None:
    with (OUT / "chunks.jsonl").open("a") as fh:
        fh.write(json.dumps({"task": name, "commit": head, **info}) + "\n")


def task_validate(head: str) -> None:
    t0 = time.time()
    exact = exact_dos(VAL_N)
    st = new_state(VAL_N, -(VAL_N * VAL_N // 2) - 4 * VAL_N, VAL_N + VAL_SMAX_EXTRA, seeded(ENT_VAL, VAL_N))
    st = wl_run(st, VAL_STEPS)
    est = normalized_lng(st, float(math.factorial(VAL_N)))
    rows = []
    for s, cnt in sorted(exact.items()):
        if cnt >= VAL_MIN_STATES:
            d = abs(est.get(s, float("-inf")) - math.log(cnt))
            rows.append({"S": s, "exact": cnt, "abs_dlng": d})
    worst = max(r["abs_dlng"] for r in rows)
    res = {"rows": rows, "max_abs_dlng": worst, "pass": worst <= VAL_TOL, "missing_exact_bins": sorted(set(exact) - set(est))}
    prev = OUT / "validate.json"
    if prev.exists():  # EQ-A2: el intento fallido se conserva, no se sobrescribe
        prev.rename(OUT / "validate_intento1_fallido.json")
    (OUT / "validate.json").write_text(json.dumps({"commit": head, "data": res}, indent=2))
    log_chunk("validate", head, {"seconds": round(time.time() - t0, 1)})
    print("validate", res["pass"], round(worst, 4))


def task_wl(n: int, chunk: int, head: str) -> None:
    p = ckpt(n)
    if p.exists():
        st = WLState.from_json(json.loads(p.read_text()))
    else:
        smin, smax = window(n)
        st = new_state(n, smin, smax, seeded(ENT_WL, n))
    left = WL_BUDGET[n] - st.t
    if left <= 0:
        print("presupuesto WL agotado")
        return
    steps = min(chunk, left)
    t0 = time.time()
    st = wl_run(st, steps)
    p.write_text(json.dumps(st.to_json()))
    log_chunk(f"wl{n}", head, {"t": st.t, "lnf": st.lnf, "one_over_t": st.one_over_t, "seconds": round(time.time() - t0, 1)})
    print(f"wl{n}", st.t, st.lnf, round(time.time() - t0, 1), "s")


def task_prod(n: int, head: str) -> None:
    p = ckpt(n)
    st = WLState.from_json(json.loads(p.read_text()))
    if st.t < WL_BUDGET[n]:
        sys.exit("WL incompleto")
    per_block = PROD_BUDGET[n] // BLOCKS
    if st.prod_t >= PROD_BUDGET[n]:
        print("produccion completa")
        return
    if st.prod_t == 0:
        st.rng_state = dict(seeded(ENT_PROD, n).bit_generator.state)
    delta2, ref2 = classifier(n)
    t0 = time.time()
    st = prod_run(st, per_block, PROD_EVERY, lambda c: flag(c, delta2, ref2))
    block = st.prod_t // per_block - 1
    (OUT / f"block_{n}_{block}.json").write_text(json.dumps({"visits": st.prod_visits, "flags": st.prod_flags}))
    p.write_text(json.dumps(st.to_json()))
    log_chunk(f"prod{n}", head, {"prod_t": st.prod_t, "block": block, "seconds": round(time.time() - t0, 1)})
    print(f"prod{n} bloque {block}", round(time.time() - t0, 1), "s")


def flag_curve(lng: dict[int, float], visits: np.ndarray[Any, Any], flags: np.ndarray[Any, Any], smin: int, n: int) -> tuple[list[float], list[float]]:
    ss = np.array([s for s in sorted(lng) if visits[(s - smin) // 2] > 0])
    lg = np.array([lng[int(s)] for s in ss])
    phi = np.array([flags[(s - smin) // 2] / visits[(s - smin) // 2] for s in ss])
    all_s = np.array(sorted(lng))
    all_lg = np.array([lng[int(s)] for s in all_s])
    curve, missing = [], []
    for b in B_SCAN:
        w = lg - b / n * ss
        m = w.max()
        p = np.exp(w - m)
        curve.append(float(np.sum(p * phi) / np.sum(p)))
        wa = all_lg - b / n * all_s
        ma = wa.max()
        pa = np.exp(wa - ma)
        missing.append(float(1 - np.sum(p) * math.exp(m - ma) / np.sum(pa)))
    return curve, missing


def b_star(curve: list[float]) -> float:
    for b, f in zip(B_SCAN, curve):
        if f < 0.5:
            return float(b)
    return float("inf")


def task_final() -> None:
    up = {int(k): float(v) for k, v in json.loads((Path("results") / "l3ce_prime" / "summary.json").read_text())["b_c"].items()}
    down = {int(k): float(v) for k, v in json.loads((Path("results") / "l3ce_hyst" / "summary.json").read_text())["b_c_down"].items()}
    val = json.loads((OUT / "validate.json").read_text())["data"]
    commits = sorted({json.loads(line)["commit"] for line in (OUT / "chunks.jsonl").read_text().splitlines()})
    res: dict[str, Any] = {"phase": "final", "commits": commits, "V": val["pass"], "V_max_abs_dlng": val["max_abs_dlng"]}
    bstars: dict[int, float] = {}
    for n in (32, 64):
        st = WLState.from_json(json.loads(ckpt(n).read_text()))
        lng = normalized_lng(st, float(math.factorial(n)))
        visits, flags = np.array(st.prod_visits), np.array(st.prod_flags)
        curve, missing = flag_curve(lng, visits, flags, st.smin, n)
        bs = b_star(curve)
        blocks = []
        prev_v, prev_f = np.zeros_like(visits), np.zeros_like(flags)
        for k in range(BLOCKS):
            d = json.loads((OUT / f"block_{n}_{k}.json").read_text())
            cv, cf = np.array(d["visits"]), np.array(d["flags"])
            blocks.append(b_star(flag_curve(lng, cv - prev_v, cf - prev_f, st.smin, n)[0]))
            prev_v, prev_f = cv, cf
        sens = []  # EQ-A2: sensibilidad de b* a ruido 0.2 en ln g (sin voto)
        rs = np.random.Generator(np.random.PCG64(np.random.SeedSequence(20261036, spawn_key=(n,))))
        for _ in range(10):
            noisy = {k: val_ + float(rs.normal(0, 0.2)) for k, val_ in lng.items()}
            sens.append(b_star(flag_curve(noisy, visits, flags, st.smin, n)[0]))
        i = int(np.argmin(np.abs(B_SCAN - bs))) if math.isfinite(bs) else -1
        bstars[n] = bs
        res[f"N{n}"] = {
            "b_star": bs, "b_star_blocks": blocks, "b_star_block_sd": float(np.std(blocks, ddof=1)),
            # EQ-A1: cobertura = alcanzar al menos el mejor estado conocido (sandwich, S = -N^2/2 - N + 6)
            "ground_state_visited": min(lng) <= -(n * n // 2) - n + 6, "S_min_visitado": min(lng),
            "S_sandwich": -(n * n // 2) - n + 6, "toca_borde_inferior": min(lng) == st.smin, "bins_visitados": len(lng),
            "peso_sin_phi_en_b_star": missing[i] if i >= 0 else None,
            "consistente": bool(bs <= up[n] + 1e-9 and bs > down[n] / STEP),
            "b_c_up": up[n], "b_c_down": down[n], "b_star_sensibilidad_lng02": [min(sens), max(sens)],
            "curve_sample": {str(b): round(c, 3) for b, c in zip(B_SCAN[::50], curve[::50])},
        }
    ok_pre = val["pass"] and all(res[f"N{n}"]["ground_state_visited"] and res[f"N{n}"]["consistente"] for n in (32, 64))
    if not ok_pre:
        res["L-3c-EQ"] = "INDETERMINADO"
    elif bstars[64] >= bstars[32] / STEP - 1e-9:
        res["L-3c-EQ"] = "VENTANA DE EQUILIBRIO ESTABLE con beta = b/N"
    else:
        res["L-3c-EQ"] = "SE ANGOSTA EN EQUILIBRIO"
    (OUT / "summary.json").write_text(json.dumps(res, indent=2, default=str))
    print(json.dumps(res, indent=2, default=str))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=["validate", "wl", "prod", "final"], required=True)
    ap.add_argument("--n", type=int, default=0)
    ap.add_argument("--chunk", type=int, default=BATCH)
    args = ap.parse_args()
    head, dirty = git_state()
    if dirty:
        sys.exit("exige arbol git limpio (fuera de results/)")
    OUT.mkdir(parents=True, exist_ok=True)
    if args.task == "validate":
        task_validate(head)
    elif args.task == "wl":
        task_wl(args.n, args.chunk, head)
    elif args.task == "prod":
        task_prod(args.n, head)
    else:
        task_final()


if __name__ == "__main__":
    main()
