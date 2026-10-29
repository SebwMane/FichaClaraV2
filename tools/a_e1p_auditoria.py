"""A-E1': auditoria descriptiva de E1' en su dominio original y paso de escala (docs/OMEGA_CONSEJO_REV17_ES0.md §7.2,
firmado en §7.3). Sin pesos ni veredicto.

Ejecutor reanudable por partes (una parte = una corrida = un orden):
  run [--fam geo2d] [--N 512 1024] [--reps k] [--workers 4] [--out DIR --allow-dirty] [--no-perp]
  summarize [--out DIR]  -> <out>/summary.json y tablas
  figure [--out DIR] [--fig PNG]  -> docs/figuras/a_e1p_auditoria.png

Semillas: SeedSequence(20261103, spawn_key=(328, indice de familia, N, repeticion)); indices en el orden de FAMS.
Instrumento (igual que E1'/E1''): region = interval_balanced(c, N // 8) (N/8 +- 5%), corte global
c sqrt(N)/(4 pi), corte de region c sqrt(|R|)/(4 pi), ssee truncado, c en {0.5, 1, 2}; una sola diagonalizacion
real de K^T K (omega.causality.sj_multicut, exacta; respaldo por la ruta de referencia cerca de un umbral).
KR: omega.causality.orders.kr_order (la familia "kr" de E1' y E1'', tools/e0_e1.py), capas N/4, N/2, N/4.
Tiempos (un hilo BLAS por proceso): t_fast = sj_fast.evaluate (evaluador exacto R17-3 sin estado, c = 1);
t_core = sj_fast._s_region con M ya calculada (nucleo por propuesta del evaluador incremental, sin la
actualizacion O(N^2) de M y B); t_multicut = los tres umbrales y diagnosticos.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # un hilo BLAS por proceso (--workers)

import argparse
import json
import math
import multiprocessing as mp
import subprocess
import sys
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omega.causality.orders import kr_order, random_bipartite
from omega.causality.sj_fast import GUARD, _s_region, evaluate
from omega.causality.sj_multicut import multicut
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond

EXPERIMENT_ID = 328
ENTROPY = 20261103
NS = (64, 128, 256, 512, 1024, 2048)
CS = (0.5, 1.0, 2.0)
WINDOWS = {"W1": (64, 128, 256), "W2": (256, 512, 1024), "W3": (512, 1024, 2048), "ALL": NS}
A_REF = 1 / 3
# criterio de N minima para el paso de escala (§7.2): geo2d, c = 1
CRIT_DA, CRIT_SE, CRIT_MODES = 0.1, 0.05, 10
VALID_FRAC_EB = 0.15  # SUPUESTO: fraccion de propuestas validas en E-b (E-a: todas validas)
CHAINS, SWEEPS = 4, 1000
SPEC = Path("docs") / "OMEGA_CONSEJO_REV17_ES0.md"
OUT = Path("results") / "omega_a_e1p"
FIG = Path("docs") / "figuras" / "a_e1p_auditoria.png"
Bool = npt.NDArray[np.bool_]


def _chain(n: int, rng: np.random.Generator) -> Bool:
    out: Bool = np.triu(np.ones((n, n), dtype=bool), k=1)
    return out


FAMS: dict[str, tuple[Callable[[int, np.random.Generator], Bool], int, str]] = {
    "geo2d": (lambda n, r: causal_matrix(sprinkle_diamond(n, 2, r)), 30, "sprinkling 2D"),
    "altura2": (lambda n, r: random_bipartite(n, 0.5, r), 10, "altura 2"),
    "kr3": (lambda n, r: kr_order(n, r), 10, "KR 3 capas"),
    "cadena": (_chain, 1, "cadena total"),
    "geo3d": (lambda n, r: causal_matrix(sprinkle_diamond(n, 3, r)), 10, "sprinkling 3D"),
    "geo4d": (lambda n, r: causal_matrix(sprinkle_diamond(n, 4, r)), 10, "sprinkling 4D"),
}


def rng_for(fam: str, n: int, rep: int) -> np.random.Generator:
    key = (EXPERIMENT_ID, list(FAMS).index(fam), n, rep)
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(ENTROPY, spawn_key=key)))


def git_info() -> tuple[str, bool]:
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=False).stdout.strip()
    st = subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True,
                        check=False)
    return head, bool(st.stdout.strip())


def part_path(out: Path, fam: str, n: int, rep: int) -> Path:
    return out / "parts" / f"{fam}_N{n}_r{rep:02d}.json"


def _num(x: float) -> float | str | None:
    if math.isnan(x):
        return None
    return x if math.isfinite(x) else ("inf" if x > 0 else "-inf")


def measure(c: Bool, perp: bool) -> dict[str, Any]:
    """Una corrida: instrumento con los tres umbrales, diagnosticos y tiempos."""
    n = c.shape[0]
    t0 = time.perf_counter()
    res = multicut(c, CS, GUARD, 1.0 if perp else None)
    t_multi = time.perf_counter() - t0
    t0 = time.perf_counter()
    fast = evaluate(c)
    t_fast = time.perf_counter() - t0
    r = res.region
    row: dict[str, Any] = {"n_R": int(r.size), "has_region": bool(r.size), "t_multicut_s": t_multi,
                           "t_fast_s": t_fast, "fast_fallback": bool(fast.fallback)}
    if r.size == 0:
        return row
    cf = c.astype(np.float64)
    k = cf - cf.T
    m = k.T @ k
    t0 = time.perf_counter()
    _s_region(c, k, m, r, GUARD)
    row["t_core_s"] = time.perf_counter() - t0
    past = c.sum(axis=0)
    fut = c.sum(axis=1)
    sub = c[np.ix_(r, r)]
    row["frac_R"] = float(r.size / n)
    row["frac_pasado"] = float(past[int(r[int(np.argmin(sub.sum(axis=0)))])] / n)  # como E1'
    row["frac_futuro"] = float(fut[int(r[int(np.argmin(sub.sum(axis=1)))])] / n)
    for cr in res.cuts:
        tag = f"c{cr.c_mult}"
        row[f"S_{tag}"] = cr.s
        row[f"kept_global_{tag}"] = cr.n_kept_global
        row[f"kept_region_{tag}"] = cr.n_kept_region
        row[f"kept_region_frac_{tag}"] = cr.n_kept_region / r.size
        row[f"margin_global_{tag}"] = _num(cr.margin_global)
        row[f"margin_region_{tag}"] = _num(cr.margin_region)
        row[f"fallback_{tag}"] = cr.fallback
    row["dS_fast_c1"] = abs(fast.s - row["S_c1.0"])
    if perp:
        row["n_Rperp"] = res.n_perp
        row["S_perp_c1.0"] = res.s_perp
        row["kept_region_perp_c1.0"] = res.n_kept_perp
        row["margin_region_perp_c1.0"] = _num(res.margin_perp)
    return row


def run_part(task: dict[str, Any]) -> str:
    fam, n, rep = task["fam"], task["N"], task["rep"]
    p = part_path(Path(task["out"]), fam, n, rep)
    if p.exists():
        return f"{p.name} ya existe"
    t0 = time.perf_counter()
    c = FAMS[fam][0](n, rng_for(fam, n, rep))
    t_gen = time.perf_counter() - t0
    row = measure(c, task["perp"])
    rec = {"spec": {"doc": str(SPEC), "experiment_id": EXPERIMENT_ID, "entropy": ENTROPY, "fam": fam, "N": n,
                    "rep": rep, "spawn_key": [EXPERIMENT_ID, list(FAMS).index(fam), n, rep], "cs": list(CS),
                    "guard": GUARD, "official": task["official"]},
           "commit": task["commit"], "dirty": task["dirty"], "fam": fam, "N": n, "rep": rep,
           "n_relations": int(c.sum()), "t_gen_s": t_gen, **row}
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(rec, allow_nan=False))
    tmp.replace(p)
    return f"{p.name} ok |R|={row['n_R']} t_fast={1e3 * row['t_fast_s']:.1f} ms total={time.perf_counter() - t0:.1f} s"


def cmd_run(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    head, dirty = git_info()
    if dirty and not args.allow_dirty:
        sys.exit("exige arbol git limpio (excluido results/); --allow-dirty solo con --out fuera de la ruta oficial")
    if args.allow_dirty and out.resolve() == OUT.resolve():
        sys.exit("--allow-dirty exige --out distinto de la ruta oficial")
    if args.reps is not None and out.resolve() == OUT.resolve():
        sys.exit("--reps exige --out (no se mezcla con la corrida oficial)")
    official = not dirty and out.resolve() == OUT.resolve()
    fams = [args.fam] if args.fam else list(FAMS)
    ns = list(args.N) if args.N else list(NS)
    base: dict[str, Any] = {"out": str(out), "commit": head, "dirty": dirty, "official": official, "perp": not args.no_perp}
    tasks = [{**base, "fam": f, "N": n, "rep": rep} for f in fams for n in ns
             for rep in range(min(FAMS[f][1], args.reps) if args.reps is not None else FAMS[f][1])]
    tasks = [t for t in tasks if not part_path(out, t["fam"], t["N"], t["rep"]).exists()]
    tasks.sort(key=lambda t: -t["N"])  # las grandes primero: mejor reparto entre procesos
    print(f"{len(tasks)} partes por correr en {out} (commit {head[:10]}, sucio={dirty}, oficial={official})")
    if args.workers <= 1:
        for t in tasks:
            print(run_part(t), flush=True)
        return
    with mp.get_context("spawn").Pool(args.workers) as pool:
        for msg in pool.imap_unordered(run_part, tasks):
            print(msg, flush=True)


# --- resumen ----------------------------------------------------------------------------------------


def load_rows(out: Path) -> list[dict[str, Any]]:
    return [json.loads(p.read_text()) for p in sorted((out / "parts").glob("*.json"))]


def ols(x: npt.NDArray[np.float64], y: npt.NDArray[np.float64]) -> dict[str, Any]:
    """Pendiente por MCO con su error estandar (como e1prima.criterion)."""
    if x.size < 3 or float(np.ptp(x)) == 0.0:
        return {"n": int(x.size), "a": None, "se": None, "r2": None}
    a = np.vstack([x, np.ones_like(x)]).T
    coef = np.linalg.lstsq(a, y, rcond=None)[0]
    res = y - a @ coef
    tot = float(np.sum((y - y.mean()) ** 2)) or 1e-300
    se = math.sqrt(float(res @ res) / (x.size - 2) * float(np.linalg.inv(a.T @ a)[0, 0]))
    return {"n": int(x.size), "a": float(coef[0]), "b": float(coef[1]), "se": se, "r2": 1 - float(res @ res) / tot}


def xlog(n_r: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    out: npt.NDArray[np.float64] = np.log(np.sqrt(n_r) / (4 * math.pi))
    return out


def area_fit(rows: list[dict[str, Any]], ns: tuple[int, ...], cc: float) -> dict[str, Any]:
    ok = [r for r in rows if r["N"] in ns and r["has_region"]]
    x = xlog(np.array([r["n_R"] for r in ok], dtype=float))
    y = np.array([r[f"S_c{cc}"] for r in ok], dtype=float)
    f = ols(x, y)
    f["da"] = abs(f["a"] - A_REF) if f["a"] is not None else None
    return f


def _f(v: Any) -> float:
    return math.inf if v == "inf" else float(v)


def _stats(vals: list[float]) -> dict[str, float | None]:
    if not vals:
        return {"media": None, "de": None}
    a = np.array(vals, dtype=float)
    return {"media": float(a.mean()), "de": float(a.std(ddof=1)) if a.size > 1 else 0.0}


def cell(rows: list[dict[str, Any]]) -> dict[str, Any]:
    ok = [r for r in rows if r["has_region"]]
    d: dict[str, Any] = {"corridas": len(rows), "frac_con_region": len(ok) / len(rows) if rows else None}
    if not ok:
        return d

    def med(k: str) -> float:
        return float(np.median([r[k] for r in ok]))

    d.update({"R_mediana": med("n_R"), "R_sobre_N_mediana": med("frac_R"), "frac_pasado_mediana": med("frac_pasado"),
              "frac_futuro_mediana": med("frac_futuro")})
    for cc in CS:
        t = f"c{cc}"
        d[t] = {"modos_globales_mediana": med(f"kept_global_{t}"), "modos_region_mediana": med(f"kept_region_{t}"),
                "modos_region_min": int(min(r[f"kept_region_{t}"] for r in ok)),
                "frac_modos_region_mediana": med(f"kept_region_frac_{t}"),
                "margen_global_min": min(_f(r[f"margin_global_{t}"]) for r in ok),
                "margen_region_min": min(_f(r[f"margin_region_{t}"]) for r in ok),
                "respaldos": int(sum(r[f"fallback_{t}"] for r in ok)),
                "S": _stats([r[f"S_{t}"] for r in ok])}
    perp = [r for r in ok if r.get("S_perp_c1.0") is not None]
    if perp:
        d["S_perp_c1"] = _stats([r["S_perp_c1.0"] for r in perp])
        d["Rperp_mediana"] = float(np.median([r["n_Rperp"] for r in perp]))
        d["dif_rel_S_Sperp_media"] = float(np.mean([abs(r["S_c1.0"] - r["S_perp_c1.0"]) / abs(r["S_c1.0"])
                                                    for r in perp if r["S_c1.0"] != 0])) if any(
            r["S_c1.0"] != 0 for r in perp) else None
    d["dS_fast_c1_max"] = max(r["dS_fast_c1"] for r in ok)
    return d


def power_fit(rows: list[dict[str, Any]], key: str) -> dict[str, Any]:
    pts = [(r["N"], r[key]) for r in rows if r.get(key)]
    if len({n for n, _ in pts}) < 2:
        return {"exponente": None}
    x = np.log(np.array([p[0] for p in pts], dtype=float))
    y = np.log(1e3 * np.array([p[1] for p in pts], dtype=float))
    f = ols(x, y) if x.size >= 3 else {"a": float(np.polyfit(x, y, 1)[0]), "b": float(np.polyfit(x, y, 1)[1]),
                                       "se": None}
    return {"exponente": f["a"], "se": f["se"], "log_ms_intercepto": f["b"], "puntos": len(pts)}


def ms_at(fit: dict[str, Any], n: float) -> float | None:
    if fit.get("exponente") is None:
        return None
    return float(math.exp(fit["log_ms_intercepto"] + fit["exponente"] * math.log(n)))


def scale_step(rows: list[dict[str, Any]], fits: dict[str, Any], cells: dict[str, Any]) -> dict[str, Any]:
    geo = [r for r in rows if r["fam"] == "geo2d"]
    pf = {k: power_fit(geo, k) for k in ("t_fast_s", "t_core_s", "t_multicut_s")}
    med_ms = {}
    for n in NS:
        g = [r for r in geo if r["N"] == n and r["has_region"]]
        if g:
            med_ms[str(n)] = {k: 1e3 * float(np.median([r[k] for r in g])) for k in ("t_fast_s", "t_core_s")}
    # N minima: menor N de una ventana cuyo ajuste de geo2d (c = 1) cumple |a - 1/3| <= 0.1 y SE < 0.05, tal que esa N
    # y todas las N mayores de la ventana tengan mediana de modos de region (c = 1) >= 10
    cand: list[tuple[int, str]] = []
    for w, ns in WINDOWS.items():
        if w == "ALL":
            continue
        f = fits.get("geo2d", {}).get(w, {}).get("1.0", {})
        if f.get("a") is None or not (f["da"] <= CRIT_DA and f["se"] < CRIT_SE):
            continue
        for i, n in enumerate(ns):
            if all((cells.get("geo2d", {}).get(str(m), {}).get("c1.0", {}).get("modos_region_mediana") or 0) >= CRIT_MODES
                   for m in ns[i:]):
                cand.append((n, w))
                break
    n_star = min(n for n, _ in cand) if cand else None

    def project(n: int) -> dict[str, Any]:
        ms_meas = med_ms.get(str(n), {})
        ms_fit = {k: ms_at(pf[k], n) for k in ("t_fast_s", "t_core_s")}
        props = n * (n - 1) / 2
        d: dict[str, Any] = {"N": n, "ms_medido": ms_meas, "ms_ajuste": ms_fit}
        for k in ("t_fast_s", "t_core_s"):
            ms = ms_meas.get(k) or ms_fit[k]
            if ms is None:
                continue
            src = "medido (mediana)" if ms_meas.get(k) else "ajuste de potencia"
            s_eb = VALID_FRAC_EB * props * ms / 1e3
            s_ea = props * ms / 1e3
            d[k] = {"ms_por_propuesta": ms, "fuente": src, "s_por_barrido_Eb": s_eb, "s_por_barrido_Ea": s_ea,
                    "horas_cpu_Eb": CHAINS * SWEEPS * s_eb / 3600, "horas_cpu_Ea": CHAINS * SWEEPS * s_ea / 3600}
        return d

    return {"ajuste_potencia_geo2d": pf, "ms_mediana_geo2d": med_ms,
            "criterio_N_minima": {"da_max": CRIT_DA, "se_max": CRIT_SE, "modos_region_min": CRIT_MODES, "c": 1.0},
            "candidatas": [{"N": n, "ventana": w} for n, w in cand], "N_minima": n_star,
            "supuestos": {"frac_validas_Eb": VALID_FRAC_EB, "frac_validas_Ea": 1.0, "cadenas": CHAINS,
                          "barridos": SWEEPS, "propuestas_por_barrido": "N(N-1)/2",
                          "nota": "frac_validas_Eb = 0.15 es un supuesto; un hilo BLAS por evaluacion"},
            "proyeccion": {str(n): project(n) for n in NS},
            "proyeccion_N_minima": project(n_star) if n_star else "ninguna N cumple el criterio"}


def _clean(x: Any) -> Any:
    if isinstance(x, float) and not math.isfinite(x):
        return None if math.isnan(x) else ("inf" if x > 0 else "-inf")
    if isinstance(x, dict):
        return {k: _clean(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_clean(v) for v in x]
    return x


def _fmt(v: Any, p: int = 3) -> str:
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:.{p}f}" if math.isfinite(v) else "inf"
    return str(v)


def cmd_summarize(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    rows = load_rows(out)
    if not rows:
        sys.exit("no hay partes")
    by: dict[str, list[dict[str, Any]]] = {f: [r for r in rows if r["fam"] == f] for f in FAMS}
    cells = {f: {str(n): cell([r for r in rs if r["N"] == n]) for n in NS if any(r["N"] == n for r in rs)}
             for f, rs in by.items() if rs}
    fits = {f: {w: {str(cc): area_fit(rs, ns, cc) for cc in CS} for w, ns in WINDOWS.items()} for f, rs in by.items() if rs}
    diffs: dict[str, Any] = {}
    for f in ("kr3", "geo3d", "geo4d", "altura2", "cadena"):
        if f not in fits or "geo2d" not in fits:
            continue
        diffs[f] = {}
        for w in WINDOWS:
            diffs[f][w] = {}
            for cc in CS:
                a, g = fits[f][w][str(cc)], fits["geo2d"][w][str(cc)]
                if a["a"] is None or g["a"] is None:
                    diffs[f][w][str(cc)] = None
                    continue
                se = math.hypot(a["se"], g["se"])
                diffs[f][w][str(cc)] = {"dif": a["a"] - g["a"], "se_comb": se,
                                        "z": (a["a"] - g["a"]) / se if se > 0 else None,
                                        "dentro_2se": bool(abs(a["a"] - g["a"]) <= 2 * se)}
    expected = {f: FAMS[f][1] * len(NS) for f in FAMS}
    complete = {f: len(by[f]) == expected[f] for f in FAMS}
    commits = sorted({r["commit"] for r in rows})
    res: dict[str, Any] = {"spec": str(SPEC), "experiment_id": EXPERIMENT_ID, "entropy": ENTROPY, "commits": commits,
           "alguna_parte_sucia": any(r["dirty"] for r in rows), "oficial": all(r["spec"]["official"] for r in rows),
           "corridas": {f: len(by[f]) for f in FAMS}, "completa": complete, "kr": "orders.kr_order (como E1'/E1'')",
           "celdas": cells, "ley_area": fits, "dif_respecto_geo2d": diffs,
           "escala": scale_step(rows, fits, cells)}
    (out / "summary.json").write_text(json.dumps(_clean(res), indent=2, allow_nan=False))
    print(f"{len(rows)} corridas; commits {commits}; completa {complete}")
    print("\n== Por familia y N (c = 0.5 / 1 / 2) ==")
    print(f"{'fam':8} {'N':>5} {'n':>3} {'reg':>5} {'|R|':>6} {'pas':>5} {'fut':>5} | modos region (frac)       "
          f"| margen min glob/reg (c=1)   | S media +- de")
    for f, cs in cells.items():
        for n, d in cs.items():
            if "c1.0" not in d:
                print(f"{f:8} {n:>5} {d['corridas']:>3} {_fmt(d['frac_con_region'], 2):>5}   (sin region)")
                continue
            modes = " ".join(f"{d[f'c{cc}']['modos_region_mediana']:.0f}({d[f'c{cc}']['frac_modos_region_mediana']:.2f})"
                             for cc in CS)
            mg = " ".join(f"{d[f'c{cc}']['margen_global_min']:.1e}/{d[f'c{cc}']['margen_region_min']:.1e}" for cc in (CS[1],))
            ss = " ".join(f"{d[f'c{cc}']['S']['media']:.3f}+-{d[f'c{cc}']['S']['de']:.3f}" for cc in CS)
            print(f"{f:8} {n:>5} {d['corridas']:>3} {d['frac_con_region']:>5.2f} {d['R_mediana']:>6.0f} "
                  f"{d['frac_pasado_mediana']:>5.2f} {d['frac_futuro_mediana']:>5.2f} | {modes:25} | {mg:27} | {ss}")
    print("\n== Coeficiente de ley de area a (S vs log(sqrt|R|/4pi)), por ventana; |a - 1/3| ==")
    for f, fw in fits.items():
        for w in WINDOWS:
            txt = []
            for cc in CS:
                q = fw[w][str(cc)]
                txt.append(f"c={cc}: " + ("-" if q["a"] is None else f"{q['a']:.3f}+-{q['se']:.3f} (|da|={q['da']:.3f})"))
            print(f"{f:8} {w:4} n={fw[w]['1.0']['n']:>3}  " + "   ".join(txt))
    print("\n== a_fam - a_geo2d (SE combinada), por ventana ==")
    for f, dw in diffs.items():
        for w in WINDOWS:
            txt = [f"c={cc}: " + ("-" if dw[w][str(cc)] is None else
                                  f"{dw[w][str(cc)]['dif']:+.3f}+-{dw[w][str(cc)]['se_comb']:.3f}") for cc in CS]
            print(f"{f:8} {w:4}  " + "   ".join(txt))
    sc: dict[str, Any] = res["escala"]
    print("\n== Escala (geo2d, un hilo BLAS) ==")
    for k, v in sc["ajuste_potencia_geo2d"].items():
        print(f"{k}: ms ~ N^{_fmt(v.get('exponente'))} (+-{_fmt(v.get('se'))})")
    for n, d in sc["proyeccion"].items():
        q = d.get("t_fast_s")
        if q:
            print(f"N={n:>5}: {q['ms_por_propuesta']:.2f} ms/eval ({q['fuente']}); s/barrido E-b {q['s_por_barrido_Eb']:.1f}, "
                  f"E-a {q['s_por_barrido_Ea']:.1f}; horas CPU 4x1000 E-b {q['horas_cpu_Eb']:.1f}, E-a {q['horas_cpu_Ea']:.1f}")
    print(f"N minima (criterio §7.2 interpretado): {sc['N_minima']} (candidatas {sc['candidatas']})")
    print(f"resumen: {out / 'summary.json'}")


# --- figura -----------------------------------------------------------------------------------------


def cmd_figure(args: argparse.Namespace) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(args.out) if args.out else OUT
    fig_path = Path(args.fig) if args.fig else (out / "a_e1p_auditoria.png" if args.out else FIG)
    rows = load_rows(out)
    if not rows:
        sys.exit("no hay partes")
    summ_p = out / "summary.json"
    summ: dict[str, Any] | None = json.loads(summ_p.read_text()) if summ_p.exists() else None
    colors = dict(zip(FAMS, ("#1f77b4", "#7f7f7f", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e")))
    fig, axes = plt.subplots(2, 2, figsize=(12, 9.5))
    ax = axes[0, 0]
    for f in FAMS:
        ok = [r for r in rows if r["fam"] == f and r["has_region"]]
        if ok:
            ax.scatter(xlog(np.array([r["n_R"] for r in ok], dtype=float)), [r["S_c1.0"] for r in ok], s=12,
                       alpha=0.7, color=colors[f], label=FAMS[f][2])
    geo = [r for r in rows if r["fam"] == "geo2d" and r["has_region"]]
    if geo:
        x = xlog(np.array([r["n_R"] for r in geo], dtype=float))
        y = np.array([r["S_c1.0"] for r in geo])
        xx = np.linspace(float(np.min(x)), float(np.max(x)), 2)
        ax.plot(xx, float(y.mean()) + A_REF * (xx - float(x.mean())), "k--", lw=1.2, label="pendiente 1/3 (referencia)")
    ax.set_xlabel("log(√|R| / 4π)")
    ax.set_ylabel("S (c = 1)")
    ax.set_title("(1) S frente a log(√|R| / 4π), c = 1", fontsize=10)
    ax.legend(fontsize=7)
    ax = axes[0, 1]
    for f in FAMS:
        pts = []
        for n in NS:
            kept = [r["kept_region_c1.0"] for r in rows if r["fam"] == f and r["N"] == n and r["has_region"]]
            if kept:
                pts.append((n, float(np.median(kept))))
        if pts:
            ax.plot(*zip(*pts), "o-", color=colors[f], label=FAMS[f][2])
    ax.axhline(CRIT_MODES, color="k", ls=":", lw=1, label=f"{CRIT_MODES} modos")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("N")
    ax.set_ylabel("modos de región conservados (mediana)")
    ax.set_title("(2) Modos de la región que sobreviven al corte, c = 1", fontsize=10)
    ax.legend(fontsize=7)
    ax = axes[1, 0]
    wins = list(WINDOWS)
    if summ:
        fams_fit = [f for f in FAMS if f in summ["ley_area"]]
        width = 0.8 / max(len(fams_fit), 1)
        for i, f in enumerate(fams_fit):
            xs, ys, es = [], [], []
            for j, w in enumerate(wins):
                q = summ["ley_area"][f][w]["1.0"]
                if q["a"] is not None:
                    xs.append(j + (i - (len(fams_fit) - 1) / 2) * width)
                    ys.append(q["a"])
                    es.append(q["se"])
            if xs:
                ax.errorbar(xs, ys, yerr=es, fmt="o", ms=4, capsize=3, color=colors[f], label=FAMS[f][2])
    ax.axhline(A_REF, color="k", ls="--", lw=1, label="1/3")
    ax.set_xticks(range(len(wins)))
    ax.set_xticklabels([f"{w}\n{'-'.join(map(str, WINDOWS[w]))}" if w != "ALL" else "todas" for w in wins], fontsize=8)
    ax.set_ylabel("coeficiente a (c = 1)")
    ax.set_title("(3) Coeficiente de ley de área por tramo de N (± EE)", fontsize=10)
    ax.legend(fontsize=7)
    ax = axes[1, 1]
    for key, lab, col in (("t_fast_s", "evaluate (R17-3)", "#1f77b4"), ("t_core_s", "núcleo por propuesta", "#d62728"),
                          ("t_multicut_s", "tres cortes + diagnóstico", "#7f7f7f")):
        pts = [(r["N"], 1e3 * r[key]) for r in geo if r.get(key)]
        if not pts:
            continue
        ax.scatter(*zip(*pts), s=10, alpha=0.5, color=col)
        pf = power_fit(geo, key)
        if pf["exponente"] is not None:
            nn = np.array([min(p[0] for p in pts), max(p[0] for p in pts)], dtype=float)
            ax.plot(nn, np.exp(pf["log_ms_intercepto"]) * nn ** pf["exponente"], color=col,
                    label=f"{lab}: ∝ N^{pf['exponente']:.2f}")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("N")
    ax.set_ylabel("ms por evaluación (1 hilo)")
    ax.set_title("(4) Costo por evaluación, sprinkling 2D", fontsize=10)
    ax.legend(fontsize=7)
    tag = "" if all(r["spec"]["official"] for r in rows) else " (prueba, no oficial)"
    fig.suptitle(f"A-E1′: auditoría de E1′ en su dominio y costo de escala{tag}", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=110)
    print(fig_path)


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--fam", choices=tuple(FAMS))
    r.add_argument("--N", type=int, nargs="+", choices=NS)
    r.add_argument("--reps", type=int, help="tope de repeticiones por (familia, N); solo con --out")
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--no-perp", action="store_true", help="omite S del complemento causal")
    r.add_argument("--out")
    r.add_argument("--allow-dirty", action="store_true")
    s = sub.add_parser("summarize")
    s.add_argument("--out")
    g = sub.add_parser("figure")
    g.add_argument("--out")
    g.add_argument("--fig")
    args = ap.parse_args()
    {"run": cmd_run, "summarize": cmd_summarize, "figure": cmd_figure}[args.cmd](args)


if __name__ == "__main__":
    main()
