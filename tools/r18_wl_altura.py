"""R18-1: piloto tecnico multicanonico en la altura (docs/OMEGA_R18_WL_ALTURA_PREREGISTRO.md). No es evidencia.

Ejecutor reanudable:
  run [--phase A|B|AB] [--N 32] [--workers 4] [--sweeps-a k --sweeps-b k --out DIR --allow-dirty]
  summarize [--out DIR] [--esc0 DIR]   -> <out>/summary.json y tabla con los 5 criterios de §4
  figure [--out DIR] [--fig PNG]

Fase A (adaptacion): una cadena Wang-Landau por N desde el orden vacio (inicio 0 de ES-C0), generador
SeedSequence(20261108, spawn_key=(331, 0, N)); W congelado en <out>/W_N{N}.json con su registro.
Fase B (produccion): 4 cadenas por N desde los inicios de ES-C0 (tools/esc0_piloto.py:start_order), W congelado,
generador SeedSequence(20261109, spawn_key=(331, 1, N, cadena)); una parte por cadena en <out>/parts/N{N}_c{k}.json.
La fase B de un N exige su archivo W.

Trazas por barrido (al final del barrido): h, fraccion de relaciones, S4, links/N, sector (0 bajo, 1 alto,
-1 intermedio) y los contadores acumulados de cruces y viajes. Los cruces y viajes se siguen en cada movimiento
aceptado (resolucion de propuesta), no solo al final del barrido.

Errores (criterios 4 y 5): medias por lotes sobre la mitad final de cada cadena, NB = 20 lotes por cadena (se
descartan las primeras (largo mod NB) muestras de la mitad final). Para la media reponderada r = sum w x / sum w se
usa el estimador de razon: Z_b = sum_{t en b} w_t (x_t - r), SE = sqrt(sum_cadenas NB var(Z_b)) / sum w (cadenas
independientes). Con w = 1 es el error por lotes ordinario (asi se trata la referencia ES-C0).
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # matrices chicas; evita sobresuscripcion con --workers

import argparse
import hashlib
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
from esc0_piloto import START_NAMES, start_order

from omega.causality import mcmc_height as mh
from omega.causality.mcmc_general import action4, height

EXPERIMENT_ID = 331
ENT_A, ENT_B = 20261108, 20261109
NS = (32, 64, 128)
HMAX = {32: 6, 64: 7, 128: 8}
SWEEPS_A = {32: 3000, 64: 3000, 128: 1500}
SWEEPS_B = {32: 2000, 64: 2000, 128: 1000}
N_CHAINS = 4
LNF0, LNF_MIN = 1.0, 1e-4
WL_START = 0  # la cadena de adaptacion arranca del orden vacio
RHAT_OK, ESS_OK, TRIPS_OK, ESS_RW_OK, Z_OK = 1.05, 200.0, 10, 200.0, 2.0
NB = 20
OBS_MIX = ("rel_frac", "height", "S4")
TRACE_KEYS = ("height", "rel_frac", "S4", "links_per_N", "sector", "crossings", "trips")
ESC0_KEY = {"eps_idx": 0, "eps": 1.0, "N": 32, "b_idx": 2, "b": 0}  # ES-C0 N = 32, b = 0, eps = 1
SPEC = Path("docs") / "OMEGA_R18_WL_ALTURA_PREREGISTRO.md"
OUT = Path("results") / "r18_wl"
ESC0 = Path("results") / "omega_esc0"
FIG = Path("docs") / "figuras" / "r18_wl_altura.png"
Float = npt.NDArray[np.float64]


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


def w_path(out: Path, n: int) -> Path:
    return out / f"W_N{n}.json"


def part_path(out: Path, n: int, chain: int) -> Path:
    return out / "parts" / f"N{n}_c{chain}.json"


def _write(p: Path, rec: dict[str, Any]) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(_clean(rec)))
    tmp.replace(p)


def _load(p: Path) -> dict[str, Any] | None:
    if not p.exists():
        return None
    rec: dict[str, Any] = json.loads(p.read_text())
    return rec


# --- fase A -----------------------------------------------------------------------------------------


def run_phase_a(task: dict[str, Any]) -> str:
    n = task["N"]
    out = Path(task["out"])
    p = w_path(out, n)
    if p.exists():
        return f"{p.name} ya existe"
    hmax = HMAX[n]
    rng = rng_for(ENT_A, 0, n)
    st = mh.init_state(start_order(n, WL_START, rng))
    bias = mh.Bias(np.zeros(hmax - mh.HMIN + 1), hmax)
    lnf_trace: list[float] = []
    wall: list[float] = []
    t_last = [time.perf_counter()]

    def on_sweep(sw: int, s: mh.HeightState, lnf: float) -> None:
        now = time.perf_counter()
        wall.append(now - t_last[0])
        t_last[0] = now
        lnf_trace.append(lnf)

    t0 = time.perf_counter()
    log = mh.wang_landau(st, bias, rng, task["sweeps"], LNF0, LNF_MIN, mh.FLATNESS, on_sweep=on_sweep)
    w_norm = bias.w - bias.w[0]
    rec = {
        "spec": {"doc": str(SPEC), "experiment_id": EXPERIMENT_ID, "entropy": ENT_A, "spawn_key": [EXPERIMENT_ID, 0, n],
                 "N": n, "hmax": hmax, "hmin": mh.HMIN, "budget_sweeps": task["sweeps"], "lnf0": LNF0,
                 "lnf_min": LNF_MIN, "flatness": mh.FLATNESS, "flat_check_every_sweeps": 1,
                 "wl_update": "despues de cada propuesta, solo con h en la ventana", "start": START_NAMES[WL_START],
                 "proposals_per_sweep": n * (n - 1) // 2, "official": task["official"]},
        "commit": task["commit"], "dirty": task["dirty"],
        "W": w_norm.tolist(), "W_h": list(range(mh.HMIN, hmax + 1)),
        "frozen_without_reaching_lnf_min": not log["reached_lnf_min"],
        "log": {**{k: v for k, v in log.items() if k != "h_trace"}, "lnf_por_barrido": lnf_trace,
                "h_por_barrido": log["h_trace"]},
        "wall_s": sum(wall), "s_per_sweep": float(np.mean(wall)) if wall else None, "total_s": time.perf_counter() - t0,
    }
    _write(p, rec)
    return (f"{p.name} ok: {log['sweeps_used']} barridos, ln f final {log['lnf_final']:.2e}, "
            f"{len(log['halvings'])} mitades, {sum(wall):.1f} s")


# --- fase B -----------------------------------------------------------------------------------------


def run_phase_b(task: dict[str, Any]) -> str:
    n, chain = task["N"], task["chain"]
    out = Path(task["out"])
    p = part_path(out, n, chain)
    if p.exists():
        return f"{p.name} ya existe"
    wtxt = w_path(out, n).read_text()
    wrec = json.loads(wtxt)
    hmax = int(wrec["spec"]["hmax"])
    bias = mh.Bias(np.asarray(wrec["W"], dtype=np.float64), hmax)
    rng = rng_for(ENT_B, 1, n, chain)
    t0 = time.perf_counter()
    st = mh.init_state(start_order(n, chain, rng))
    t_init = time.perf_counter() - t0
    npairs = n * (n - 1) / 2
    initial = {"height": st.h, "rel_frac": st.nrel / npairs, "S4": st.s4}
    traces: dict[str, list[Any]] = {k: [] for k in TRACE_KEYS}
    counts: dict[str, list[int]] = {"valid": [], "accepted": []}
    wall: list[float] = []
    t_meas = 0.0
    for _ in range(task["sweeps"]):
        ts = time.perf_counter()
        a, v = mh.sweep(st, bias, rng)
        tm = time.perf_counter()
        wall.append(tm - ts)
        counts["accepted"].append(a)
        counts["valid"].append(v)
        row = {"height": st.h, "rel_frac": st.nrel / npairs, "S4": st.s4,
               "links_per_N": int(np.count_nonzero(st.b[st.c] == 0)) / n, "sector": mh.sector(st.h),
               "crossings": st.crossings, "trips": st.trips}
        for k in TRACE_KEYS:
            traces[k].append(row[k])
        t_meas += time.perf_counter() - tm
    rec = {
        "spec": {"doc": str(SPEC), "experiment_id": EXPERIMENT_ID, "entropy": ENT_B,
                 "spawn_key": [EXPERIMENT_ID, 1, n, chain], "N": n, "chain": chain, "start": START_NAMES[chain],
                 "sweeps": task["sweeps"], "proposals_per_sweep": n * (n - 1) // 2, "hmax": hmax, "hmin": mh.HMIN,
                 "sectores": {"bajo": f"h <= {mh.H_LOW}", "alto": f"h >= {mh.H_HIGH}"},
                 "official": task["official"]},
        "commit": task["commit"], "dirty": task["dirty"], "W": bias.w.tolist(),
        "W_file_sha256": hashlib.sha256(wtxt.encode()).hexdigest(), "W_file_official": wrec["spec"]["official"],
        "initial": initial, "traces": traces, "counts": counts, "wall_s": wall, "init_s": t_init, "measure_s": t_meas,
        "final_check": {"h_exact": st.h == height(st.c), "s4_drift": float(st.s4 - action4(st.c, st.b))},
    }
    _write(p, rec)
    return f"{p.name} ok {sum(wall):.1f} s barridos ({float(np.mean(wall)):.4f} s/barrido), viajes {st.trips}"


def _pool(fn: Any, tasks: list[dict[str, Any]], workers: int) -> None:
    if workers <= 1 or len(tasks) <= 1:
        for t in tasks:
            print(fn(t), flush=True)
    else:
        with mp.get_context("spawn").Pool(min(workers, len(tasks))) as pool:
            for msg in pool.imap_unordered(fn, tasks):
                print(msg, flush=True)


def cmd_run(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    override = args.sweeps_a is not None or args.sweeps_b is not None
    if override and not args.out:
        sys.exit("--sweeps-a/--sweeps-b exigen --out (no se mezcla con la corrida oficial)")
    head, dirty = git_info()
    if dirty and not args.allow_dirty:
        sys.exit("exige arbol git limpio fuera de results/ (usar --allow-dirty solo para pruebas)")
    official = not dirty and not override and out == OUT
    base: dict[str, Any] = {"out": str(out), "commit": head, "dirty": dirty, "official": official}
    if official:
        stale = [p.name for p in out.glob("W_N*.json") if not json.loads(p.read_text())["spec"]["official"]]
        stale += [p.name for p in out.glob("parts/*.json") if not json.loads(p.read_text())["spec"]["official"]]
        if stale:
            sys.exit(f"{out} contiene archivos no oficiales ({len(stale)}, p. ej. {stale[0]}); moverlos antes")
    ns = [args.N] if args.N else list(NS)
    print(f"{out} (commit {head[:8]}, dirty {dirty}, oficial {official})")
    t0 = time.time()
    if "A" in args.phase:
        ta: list[dict[str, Any]] = [{**base, "N": n, "sweeps": args.sweeps_a if args.sweeps_a is not None else SWEEPS_A[n]}
              for n in ns if not w_path(out, n).exists()]
        ta.sort(key=lambda t: -t["N"] * t["N"] * t["sweeps"])
        print(f"fase A: {len(ta)} cadenas Wang-Landau por correr")
        _pool(run_phase_a, ta, args.workers)
    if "B" in args.phase:
        tb: list[dict[str, Any]] = []
        for n in ns:
            if not w_path(out, n).exists():
                print(f"N={n}: falta {w_path(out, n).name}; la fase B de este N no corre")
                continue
            if official and not json.loads(w_path(out, n).read_text())["spec"]["official"]:
                sys.exit(f"{w_path(out, n)} no es oficial")
            sw = args.sweeps_b if args.sweeps_b is not None else SWEEPS_B[n]
            tb += [{**base, "N": n, "chain": ch, "sweeps": sw} for ch in range(N_CHAINS)
                   if not part_path(out, n, ch).exists()]
        tb.sort(key=lambda t: -t["N"] * t["N"] * t["sweeps"])
        print(f"fase B: {len(tb)} cadenas de produccion por correr")
        _pool(run_phase_b, tb, args.workers)
    print(f"pared total {time.time() - t0:.1f} s")


# --- resumen ----------------------------------------------------------------------------------------


def _batch_sums(x: Float, w: Float, r: float) -> Float:
    m = x.size - x.size % NB
    z = (w[x.size - m:] * (x[x.size - m:] - r)).reshape(NB, -1).sum(axis=1)
    out: Float = z
    return out


def ratio_mean_se(xs: list[Float], ws: list[Float]) -> tuple[float, float]:
    """Media de razon sum w x / sum w sobre las cadenas dadas y su error por lotes (cadenas independientes)."""
    num = sum(float(np.sum(w * x)) for x, w in zip(xs, ws, strict=True))
    den = sum(float(np.sum(w)) for w in ws)
    r = num / den
    var = sum(NB * float(_batch_sums(x, w, r).var(ddof=1)) for x, w in zip(xs, ws, strict=True))
    return r, math.sqrt(var) / den


def esc0_reference(esc0: Path) -> dict[str, Any]:
    k = ESC0_KEY
    xs = []
    for ch in range(N_CHAINS):
        p = esc0 / "parts" / f"e{k['eps_idx']}_N{k['N']}_b{k['b_idx']}_c{ch}.json"
        rec = json.loads(p.read_text())
        sp = rec["spec"]
        if (sp["eps"], sp["N"], sp["b"]) != (k["eps"], k["N"], k["b"]):
            raise RuntimeError(f"{p} no es la celda ES-C0 N=32, b=0, eps=1")
        x = np.asarray(rec["traces"]["rel_frac"], dtype=np.float64)
        xs.append(x[x.size // 2:])
    m, se = ratio_mean_se(xs, [np.ones_like(x) for x in xs])
    return {"celda": "ES-C0 e0_N32_b2 (eps = 1, b = 0)", "dir": str(esc0), "rel_frac_media": m, "se_lotes": se,
            "por_cadena": [float(x.mean()) for x in xs], "muestras_por_cadena": [int(x.size) for x in xs]}


def n_summary(out: Path, n: int, ref: dict[str, Any] | None) -> dict[str, Any] | None:
    wrec = _load(w_path(out, n))
    parts = [p for ch in range(N_CHAINS) if (p := _load(part_path(out, n, ch))) is not None]
    if wrec is None and not parts:
        return None
    res: dict[str, Any] = {"N": n, "hmax": HMAX[n]}
    if wrec is not None:
        lg = wrec["log"]
        res["fase_A"] = {"W": wrec["W"], "W_h": wrec["W_h"], "sweeps_used": lg["sweeps_used"],
                         "budget": wrec["spec"]["budget_sweeps"], "lnf_final": lg["lnf_final"],
                         "reached_lnf_min": lg["reached_lnf_min"], "halvings": len(lg["halvings"]),
                         "halving_sweeps": [e["sweep"] for e in lg["halvings"]],
                         "congelado_sin_llegar": wrec["frozen_without_reaching_lnf_min"],
                         "s_per_sweep": wrec["s_per_sweep"], "official": wrec["spec"]["official"],
                         "commit": wrec["commit"], "dirty": wrec["dirty"]}
    if not parts:
        res["veredicto"] = "incompleto"
        return res
    hmax = int(parts[0]["spec"]["hmax"])
    wfull = mh.Bias(np.asarray(parts[0]["W"], dtype=np.float64), hmax).full(n)
    fin = {k: [np.asarray(p["traces"][k], dtype=np.float64)[len(p["traces"][k]) // 2:] for p in parts]
           for k in OBS_MIX}
    hfin = [x.astype(np.int64) for x in fin["height"]]
    logw = [wfull[h] for h in hfin]
    lmax = max(float(np.max(lw)) for lw in logw)
    ws = [np.exp(lw - lmax) for lw in logw]
    per: dict[str, Any] = {}
    for p, w, x, h in zip(parts, ws, fin["rel_frac"], hfin, strict=True):
        ch = p["spec"]["chain"]
        nprop = p["spec"]["proposals_per_sweep"] * len(p["wall_s"])
        lw_c = wfull[h]
        w_c = np.exp(lw_c - np.max(lw_c))  # normalizacion propia: con W muy desparejo los pesos conjuntos se anulan
        r, se = ratio_mean_se([x], [w_c])
        per[str(ch)] = {"start": p["spec"]["start"], "sweeps": len(p["wall_s"]), "viajes": p["traces"]["trips"][-1],
                        "cruces": p["traces"]["crossings"][-1],
                        "accept_rate": sum(p["counts"]["accepted"]) / nprop,
                        "valid_frac": sum(p["counts"]["valid"]) / nprop,
                        "s_per_sweep": float(np.mean(p["wall_s"])), "measure_s_per_sweep": p["measure_s"] / len(p["wall_s"]),
                        "rel_frac_sesgada_mitad_final": float(x.mean()), "h_sesgada_mitad_final": float(h.mean()),
                        "rel_frac_reponderada": r, "se_lotes": se,
                        "ess_reponderacion": float(w_c.sum() ** 2 / np.sum(w_c * w_c)),
                        "peso_conjunto_frac": float(w.sum()) / sum(float(v.sum()) for v in ws),
                        "final_check": p["final_check"], "official": p["spec"]["official"],
                        "W_file_sha256": p["W_file_sha256"]}
    res["por_cadena"] = per
    res["s_per_sweep"] = float(np.mean([s for p in parts for s in p["wall_s"]]))
    res["chains"] = len(parts)
    res["official"] = all(p["spec"]["official"] for p in parts) and bool(wrec and wrec["spec"]["official"])
    res["commits"] = sorted({p["commit"] for p in parts} | ({wrec["commit"]} if wrec else set()))
    if len({p["W_file_sha256"] for p in parts}) != 1:
        res["aviso"] = "las cadenas usaron archivos W distintos"
    crit: dict[str, Any] = {}
    # 1. R-hat dividido y ESS de las cadenas sesgadas (mitad final)
    mix: dict[str, Any] = {}
    ok1 = len(parts) >= 2
    for k in OBS_MIX:
        rh = split_rhat(fin[k]) if len(parts) >= 2 else float("nan")
        taus = [tau_sokal(x)[0] for x in fin[k]]
        ess = float(sum(x.size / t for x, t in zip(fin[k], taus, strict=True)))
        mix[k] = {"rhat": rh, "tau": taus, "ess": ess}
        ok1 = ok1 and math.isfinite(rh) and rh < RHAT_OK and ess >= ESS_OK
    crit["1_rhat_ess"] = {"ok": ok1, "detalle": mix}
    # 2. viajes completos por cadena (toda la produccion)
    trips = [v["viajes"] for v in per.values()]
    crit["2_viajes"] = {"ok": min(trips) >= TRIPS_OK, "por_cadena": trips, "minimo": TRIPS_OK}
    # 3. ESS de reponderacion, 4 cadenas juntas (mitades finales)
    wall_ = np.concatenate(ws)
    ess_rw = float(wall_.sum() ** 2 / np.sum(wall_ * wall_))
    crit["3_ess_reponderacion"] = {"ok": ess_rw >= ESS_RW_OK, "ess": ess_rw}
    # 4. concordancia entre cadenas: |r_c - r| <= 2 SE_c
    r_pool, se_pool = ratio_mean_se(fin["rel_frac"], ws)
    devs = {c: (v["rel_frac_reponderada"] - r_pool) / v["se_lotes"] if v["se_lotes"] > 0 else
            (0.0 if v["rel_frac_reponderada"] == r_pool else float("inf")) for c, v in per.items()}
    crit["4_concordancia"] = {"ok": all(abs(z) <= Z_OK for z in devs.values()), "media_conjunta": r_pool,
                              "se_conjunta": se_pool, "z_por_cadena": devs}
    # 5. validacion externa (solo N = 32)
    if n == ESC0_KEY["N"] and ref is not None:
        z = (r_pool - ref["rel_frac_media"]) / math.hypot(se_pool, ref["se_lotes"])
        crit["5_validacion_externa"] = {"ok": abs(z) <= Z_OK, "z": z, "media_r18": r_pool, "se_r18": se_pool,
                                        "media_esc0": ref["rel_frac_media"], "se_esc0": ref["se_lotes"]}
    elif n == ESC0_KEY["N"]:
        crit["5_validacion_externa"] = {"ok": False, "nota": "falta la referencia ES-C0"}
    else:
        crit["5_validacion_externa"] = {"ok": None, "nota": "solo aplica con N = 32"}
    res["criterios"] = crit
    res["h_reponderada"] = float(sum(float(np.sum(w * h)) for w, h in zip(ws, hfin, strict=True)) / wall_.sum())
    complete = len(parts) == N_CHAINS
    allok = all(c["ok"] is not False for c in crit.values())
    res["veredicto"] = ("mezcla" if allok else "no convergido") if complete else "incompleto"
    return res


def cmd_summarize(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    esc0 = Path(args.esc0) if args.esc0 else ESC0
    try:
        ref: dict[str, Any] | None = esc0_reference(esc0)
    except (OSError, RuntimeError) as e:
        print(f"referencia ES-C0 no disponible: {e}")
        ref = None
    ns = {str(n): s for n in NS if (s := n_summary(out, n, ref)) is not None}
    res = {"spec": str(SPEC), "experiment_id": EXPERIMENT_ID, "no_es_evidencia": True,
           "criterio": {"rhat": RHAT_OK, "ess": ESS_OK, "viajes_por_cadena": TRIPS_OK, "ess_reponderacion": ESS_RW_OK,
                        "z_max": Z_OK, "observables_mezcla": list(OBS_MIX), "mitad": "final", "lotes_por_cadena": NB,
                        "sectores": {"bajo": mh.H_LOW, "alto": mh.H_HIGH},
                        "criterio_4": "|r_c - r_conjunta| <= 2 SE_c (SE por lotes de la cadena c)",
                        "criterio_5": "|r_conjunta - m_ESC0| <= 2 sqrt(SE^2 + SE_ESC0^2); solo N = 32",
                        "veredicto": "mezcla si se cumplen todos los criterios aplicables con 4 cadenas"},
           "referencia_esc0": ref, "N": ns}
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(_clean(res), indent=2, allow_nan=False))
    if ref is not None:
        print(f"ES-C0 N=32 b=0 eps=1: rel_frac {ref['rel_frac_media']:.4f} +- {ref['se_lotes']:.4f} (lotes)")
    print(f"{'N':>4} {'WL':>12} {'s/barr':>7} {'Rhat rel/h/S4':>16} {'ESS min':>8} {'viajes c0..c3':>14} "
          f"{'ESSrw':>7} {'r ± SE':>17} {'z max':>6} {'z ext':>6}  veredicto")
    for s in ns.values():
        fa = s.get("fase_A", {})
        wl = (f"{'ok' if fa.get('reached_lnf_min') else 'lnf'}{fa.get('lnf_final', float('nan')):.0e}/"
              f"{fa.get('sweeps_used', '-')}") if fa else "-"
        cr = s.get("criterios")
        if cr is None:
            print(f"{s['N']:>4} {wl:>12}  {s['veredicto']}")
            continue
        mx = cr["1_rhat_ess"]["detalle"]
        rh = "/".join(f"{mx[o]['rhat']:.2f}" for o in OBS_MIX)
        ess = min(mx[o]["ess"] for o in OBS_MIX)
        tr = " ".join(str(t) for t in cr["2_viajes"]["por_cadena"])
        c4 = cr["4_concordancia"]
        zmax = max(abs(z) for z in c4["z_por_cadena"].values())
        c5 = cr["5_validacion_externa"]
        zext = f"{c5['z']:.2f}" if "z" in c5 else "-"
        flags = "".join(("+" if c["ok"] else ("-" if c["ok"] is False else "·")) for c in cr.values())
        print(f"{s['N']:>4} {wl:>12} {s['s_per_sweep']:>7.4f} {rh:>16} {ess:>8.1f} {tr:>14} "
              f"{cr['3_ess_reponderacion']['ess']:>7.1f} {c4['media_conjunta']:>8.4f}±{c4['se_conjunta']:.4f} "
              f"{zmax:>6.2f} {zext:>6}  {s['veredicto']} [{flags}]")


# --- figura -----------------------------------------------------------------------------------------


def cmd_figure(args: argparse.Namespace) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(args.out) if args.out else OUT
    fig_path = Path(args.fig) if args.fig else (out / "r18_wl_altura.png" if args.out else FIG)
    rows = [n for n in NS if w_path(out, n).exists() or any(part_path(out, n, ch).exists() for ch in range(N_CHAINS))]
    if not rows:
        sys.exit("no hay resultados")
    colors = ("#1f77b4", "#d62728", "#2ca02c", "#9467bd")
    fig, axes = plt.subplots(len(rows), 3, figsize=(14, 3.1 * len(rows) + 0.8), squeeze=False,
                             gridspec_kw={"width_ratios": [1, 2.2, 1.3]})
    for r, n in enumerate(rows):
        wrec = _load(w_path(out, n))
        parts = [p for ch in range(N_CHAINS) if (p := _load(part_path(out, n, ch))) is not None]
        hmax = HMAX[n]
        ax = axes[r, 0]
        if wrec is not None:
            hw = np.asarray(wrec["W_h"])
            ax.plot(hw, wrec["W"], "o-", color="#333333")
            lg = wrec["log"]
            est = "ln f < 1e-4" if lg["reached_lnf_min"] else f"congelado con ln f = {lg['lnf_final']:.1e}"
            ax.set_title(f"N = {n}: W(h) aprendido ({lg['sweeps_used']} barridos, {est})", fontsize=8)
        ax.set_xlabel("altura h", fontsize=8)
        ax.set_ylabel("W(h)  (W(2) = 0)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax = axes[r, 1]
        for p in parts:
            ch = p["spec"]["chain"]
            y = [p["initial"]["height"], *p["traces"]["height"]]
            ax.plot(np.arange(len(y)), np.minimum(y, hmax + 2), lw=0.6, color=colors[ch], label=START_NAMES[ch],
                    drawstyle="steps-post")
        ax.axhspan(0.5, mh.H_LOW + 0.5, color="#cce5ff", alpha=0.5, lw=0)
        ax.axhspan(mh.H_HIGH - 0.5, hmax + 2.5, color="#ffd6cc", alpha=0.5, lw=0)
        ax.set_ylim(0.5, hmax + 2.5)
        ax.set_title(f"N = {n}: altura de las 4 cadenas de producción (azul: sector bajo, rosa: sector alto; "
                     f"recortada en {hmax + 2})", fontsize=8)
        ax.set_xlabel("barrido", fontsize=8)
        ax.set_ylabel("altura h", fontsize=8)
        ax.tick_params(labelsize=7)
        ax = axes[r, 2]
        if parts:
            wfull = mh.Bias(np.asarray(parts[0]["W"], dtype=np.float64), hmax).full(n)
            hs = np.concatenate([np.asarray(p["traces"]["height"], dtype=np.int64)[len(p["traces"]["height"]) // 2:]
                                 for p in parts])
            lw = wfull[hs]
            w = np.exp(lw - np.max(lw))
            top = int(min(np.max(hs), n))
            edges = np.arange(0.5, top + 1.5)
            pb, _ = np.histogram(hs, bins=edges)
            pr, _ = np.histogram(hs, bins=edges, weights=w)
            x = np.arange(1, top + 1)
            ax.bar(x - 0.2, pb / pb.sum(), width=0.4, color="#7f7f7f", label="sesgada (muestreada)")
            ax.bar(x + 0.2, pr / pr.sum(), width=0.4, color="#ff7f0e", label="reponderada por e^W(h)")
            ax.legend(fontsize=7, frameon=False)
        ax.set_title(f"N = {n}: distribución de h, mitades finales juntas", fontsize=8)
        ax.set_xlabel("altura h", fontsize=8)
        ax.set_ylabel("fracción", fontsize=8)
        ax.tick_params(labelsize=7)
    h, lab = axes[0, 1].get_legend_handles_labels()
    if h:
        leg = fig.legend(h, lab, loc="lower center", ncol=4, fontsize=8, frameon=False)
        for line in leg.get_lines():
            line.set_linewidth(2.0)
    fig.suptitle("R18-1: muestreo multicanónico en la altura, β = 0 (piloto técnico, no es evidencia)", fontsize=10)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=110)
    print(fig_path)


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--phase", choices=("A", "B", "AB"), default="AB")
    r.add_argument("--N", type=int, choices=NS)
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--sweeps-a", type=int)
    r.add_argument("--sweeps-b", type=int)
    r.add_argument("--out")
    r.add_argument("--allow-dirty", action="store_true")
    s = sub.add_parser("summarize")
    s.add_argument("--out")
    s.add_argument("--esc0")
    f = sub.add_parser("figure")
    f.add_argument("--out")
    f.add_argument("--fig")
    args = ap.parse_args()
    {"run": cmd_run, "summarize": cmd_summarize, "figure": cmd_figure}[args.cmd](args)


if __name__ == "__main__":
    main()
