"""R23: piloto (a), descendientes de órdenes no geométricos frente a un blanco 3+1. Ejecutor del preregistro
exploratorio docs/OMEGA_R23_PILOTO_A_PREREGISTRO.md (sin voto; no es evidencia confirmatoria de H_E).

Ejecutor reanudable:
  run [--part null|parents|all] [--family a,b] [--parents a-b] [--workers 4] [--budget-min 110]
      [--n0 N --nparents P --nnull K --out DIR --allow-dirty]
  summarize [--out DIR] [--refs DIR] [--nref K]  -> <out>/summary.json, tabla.csv
  figure [--out DIR] [--fig PNG]                 -> docs/figuras/r23.png
Cada `run` y `summarize` deja su log en <out>/logs/ y las versiones (python, numpy, scipy) en <out>/entorno.json.

Lo que se reutiliza (importado, sin copiar ni modificar): de `tools/r14b.py`, Φ_q (`dilute`), observables (`measure`),
el vector (`obs_vector`), la estandarización con Ledoit-Wolf (`Band`, solo media, escala, pisos y covarianza contraída),
`comp_floors`, `COMPONENTS`, `QS` y el diamante plano (`GEN["plana"]`); de `tools/r20.py`, el radio conforme
(`conformal_radius`, rango 96 de 100), `log_ratio`, `boot_indices`, `git_info`, `entorno` y el log (`_logged`); de
`tools/r19.py`, la escritura atómica (`_write`), `_check_matrix`, `_load`, `_write_csv` y `sha256`. Generadores:
`f6_intento2.d5_gemelos` y `d4_agregadas`, `f6_confirmacion.c5_links_borrados`, `omega.causality.orders.kr_order`,
`transitive_percolation`, `thick_layer_order` y `f1f2.random_dorder`.

Decisiones declaradas (todas siguen el preregistro; donde el texto no fija un detalle, se fija aquí):

Blanco (§2). Para cada q, banda de la plana con las referencias **guardadas** de R20 (`results/omega_r20/refs`, clave
  «row»): 100 de ajuste (media, escala con pisos y Ledoit-Wolf de `r14b.Band`) y 100 de calibración (radio conforme).
  Las 100 de validación de R20 no se usan. Vector común de 13 componentes en los 4 q (`COMUNES`: todas menos r_40).
  Base de cada q = las de `COMUNES` finitas en las 200 referencias (con las de R20, las 13 en los 4 q; se comprueba y
  se informa). Por par (Y, Z) en un q: S = componentes de la base finitas en los dos vectores; si |S| < 10, el par no
  es evaluable en ese q (D, u nan). Con S: `r14b.Band(ajuste, ajuste[:0], S, comp_floors())`, radio conforme con las
  distancias de calibración en S (caché por (q, S)) y D = distancia/radio de Y, de Z y del conjunto nulo.
  Un nulo al que le falte alguna componente de S queda fuera de la nula de ese par (se cuenta y se informa; con
  diamantes planos de N₀ = 8000 no se espera que ocurra).
Conjunto nulo (§2). 100 diamantes planos nuevos por q, como `r20.run_ref`: N = N₀ con q = 1 y Binomial(N₀, q) si no;
  el mismo generador da N y los puntos, otro aparte elige los intervalos de la medida. Nunca ajusta nada.
Familias (§3), con índice fijo f (posición en FAMILIAS) para el spawn_key; gen(n, rng_puntos, rng_aux) -> matriz
  booleana n × n (comprobada con `r19._check_matrix`). Las que no usan azar auxiliar ignoran rng_aux.
Brazos (§4). Padre X = gen(N₀). Para cada q: Y = Φ_q(X) (q = 1 sin azar: Y = X); Z = gen(|Y|) nuevo, con generadores
  propios (con q = 1, un segundo orden independiente de N₀). Se mide cada brazo con su propio generador de intervalos.
  Un archivo por (familia, padre) con una fila por q: {qi, q, N, obs (Y), ctrl (Z)}; se comprueba N(Z) = N(Y).
Estadísticos (§5). u_B = fracción de la nula (del par) con D < D_B, empates por la mitad. AUC_B = media de u_B sobre los
  padres evaluables en ese q. Δ = mediana sobre padres evaluables de ln(D_Y/D_Z) (`r20.log_ratio`). Bootstrap: 2000
  réplicas de padres completos (los mismos índices en los 4 q y los dos brazos); IC 95% percentil (cuantiles lineales,
  como numpy). Si en alguna réplica no queda ningún padre evaluable, el IC es nan (no clasifica).
Clasificación (§6), márgenes provisionales: visible al inicio (IC de AUC_Y con q = 1 entero ≥ 0.95); cerca (IC de AUC_Y
  con q = 1/8 entero ≤ 0.75); lejos (entero ≥ 0.95); mejora (IC de Δ con q = 1/8 entero < −ln 1.3); igual (IC de Δ con
  q = 1/8 dentro de [−ln 1.3, ln 1.3]). Lectura, en este orden de precedencia:
    1. falta de resolución: pares no evaluables en más del 20% de los padres en q = 1 o en q = 1/8 (los dos q que
       entran en la lectura) -> INDETERMINADO;
    2. no visible al inicio -> no aplicable;
    3. visible, cerca y mejora -> compatible con acercamiento;
    4. visible e igual -> efecto de tamaño o familia (se informa si queda cerca, lejos o entre ambos). Si además queda
       lejos, también cumple «sin acercamiento» y se informa en «tambien_cumple»;
    5. visible, lejos y sin mejora -> sin acercamiento;
    6. cualquier otra combinación -> INDETERMINADO.

Unidades (un archivo por unidad; las que existen se saltan, así que repetir el comando reanuda):
  <out>/nulos/q<qi>_s<sss>.json        un diamante plano del conjunto nulo
  <out>/padres/<familia>/p<pp>.json    un padre: descendiente y control nuevo en los cuatro q

Semillas: SeedSequence(entropía, spawn_key=(334, tipo, ...)).
  20261118 (conjunto nulo):
    0 = N (Binomial) y puntos              (334, 0, qi, s)
    1 = selección de intervalos            (334, 1, qi, s)
  20261119 (padres, diluciones, controles, intervalos y bootstrap), f = índice de la familia en FAMILIAS:
    0 = puntos del padre                   (334, 0, f, N0, p)
    1 = azar auxiliar del padre            (334, 1, f, N0, p)       gemelos de D5, links de C5, agregadas de D4
    2 = dilución                           (334, 2, f, N0, p, qi)   q = 1 no lo usa
    3 = intervalos del descendiente Y      (334, 3, f, N0, p, qi)
    4 = puntos del control Z               (334, 4, f, N0, p, qi)
    5 = azar auxiliar del control Z        (334, 5, f, N0, p, qi)
    6 = intervalos del control Z           (334, 6, f, N0, p, qi)
    7 = bootstrap del resumen              (334, 7, f)
  Próximos libres: EXPERIMENT_ID 335 y entropía 20261120.

Árbol limpio: `git status` limpio fuera de results/. Se registran commit, sucio, oficial (solo con árbol limpio, sin
parámetros alterados y con la salida por defecto) y sha256 del preregistro.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # un hilo por proceso; el paralelismo viene de --workers

import argparse
import hashlib
import json
import math
import multiprocessing as mp
import sys
import time
import warnings
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

_TOOLS = Path(__file__).resolve().parent
_ROOT = _TOOLS.parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_TOOLS))

import r14b  # noqa: E402
import r19  # noqa: E402
import r20  # noqa: E402
from es0_piloto import _clean  # noqa: E402
from f1f2 import random_dorder  # noqa: E402
from f6_confirmacion import c5_links_borrados  # noqa: E402
from f6_intento2 import d4_agregadas, d5_gemelos  # noqa: E402

from omega.causality.orders import kr_order, thick_layer_order, transitive_percolation  # noqa: E402

EXPERIMENT_ID = 334
ENT_NULO = 20261118
ENT_PAD = 20261119
(K_PTS, K_AUX, K_DIL, K_SEL_Y, K_CTRL_PTS, K_CTRL_AUX, K_SEL_Z, K_BOOT) = range(8)
K_NULO_PTS, K_NULO_SEL = 0, 1
N0 = 8000
QS = r14b.QS  # (1, 1/2, 1/4, 1/8)
QI_INI, QI_FIN = 0, 3  # q = 1 y q = 1/8
N_PADRES = 30
N_NULO = 100
N_REF = 100  # referencias de R20 por conjunto y por q
N_BOOT = 2000
MIN_COMPS = 10
COMUNES = tuple(i for i, c in enumerate(r14b.COMPONENTS) if c != "r_40")  # 13 componentes
VISIBLE, CERCA, LEJOS = 0.95, 0.75, 0.95
MARGEN = math.log(1.3)
MAX_NO_EVAL = 0.20
SPEC = Path("docs") / "OMEGA_R23_PILOTO_A_PREREGISTRO.md"
OUT = Path("results") / "omega_r23"
REFS = Path("results") / "omega_r20" / "refs"
FIG = Path("docs") / "figuras" / "r23.png"
Bool = npt.NDArray[np.bool_]
Float = npt.NDArray[np.float64]
FamGen = Callable[[int, np.random.Generator, np.random.Generator], Bool]

COMPATIBLE = "compatible con acercamiento"
TAMANO = "efecto de tamaño o familia"
SIN = "sin acercamiento"
NO_APLICA = "no aplicable"
INDET = "INDETERMINADO"


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


# --- Familias (§3, parámetros congelados) -------------------------------------------------------------


def _plana(n: int, r: np.random.Generator, aux: np.random.Generator) -> Bool:
    c: Bool = r14b.GEN["plana"](n, r, aux)[0]
    return c


def _sellada(g: Callable[..., Any]) -> FamGen:
    """Generador sellado de F6 (n, gen, aux) -> (C, x, info): se queda con C."""

    def gen(n: int, r: np.random.Generator, aux: np.random.Generator) -> Bool:
        c: Bool = g(n, r, aux)[0]
        return c

    return gen


_PROD4 = random_dorder(4)

# (nombre, papel, generador, origen). El orden es el índice f del spawn_key.
FAMILIAS: list[tuple[str, str, FamGen, str]] = [
    ("plana", "control nulo", _plana, 'r14b.GEN["plana"]'),
    ("D5_gemelos", "control positivo", _sellada(d5_gemelos), "f6_intento2.d5_gemelos (10%)"),
    ("C5_links_borrados", "pregunta principal", _sellada(c5_links_borrados), "f6_confirmacion.c5_links_borrados (13%)"),
    ("D4_relaciones_agregadas", "pregunta secundaria", _sellada(d4_agregadas), "f6_intento2.d4_agregadas (10%)"),
    ("KR", "no geométrico", lambda n, r, a: kr_order(n, r, 0.5), "orders.kr_order(n, rng, 0.5)"),
    ("percolacion_005", "no geométrico", lambda n, r, a: transitive_percolation(n, 0.05, r),
     "orders.transitive_percolation(n, 0.05, rng)"),
    ("producto_4", "no geométrico, autosemejante", lambda n, r, a: _PROD4(n, r), "f1f2.random_dorder(4)"),
    ("capas_6", "no geométrico", lambda n, r, a: thick_layer_order(n, 6, 0.5, r), "orders.thick_layer_order(n, 6, 0.5)"),
]
NOMBRES = [f[0] for f in FAMILIAS]
FIDX = {f[0]: i for i, f in enumerate(FAMILIAS)}
GEN: dict[str, FamGen] = {f[0]: f[2] for f in FAMILIAS}
PAPEL = {f[0]: f[1] for f in FAMILIAS}


def generar(fam: str, n: int, r: np.random.Generator, aux: np.random.Generator) -> Bool:
    """Orden de la familia con exactamente n elementos (matriz booleana n × n)."""
    return r19._check_matrix(GEN[fam](n, r, aux), n)


# --- Unidades de trabajo ------------------------------------------------------------------------------


def null_path(out: Path, qi: int, s: int) -> Path:
    return out / "nulos" / f"q{qi}_s{s:03d}.json"


def parent_path(out: Path, fam: str, p: int) -> Path:
    return out / "padres" / fam / f"p{p:02d}.json"


def run_null(task: dict[str, Any]) -> str:
    """Un diamante plano del conjunto nulo, como `r20.run_ref` (mismo generador para N y puntos)."""
    t0 = time.time()
    qi, s, n0 = task["qi"], task["s"], task["n0"]
    q = QS[qi]
    r = rng_for(ENT_NULO, K_NULO_PTS, qi, s)
    n = n0 if q >= 1.0 else int(r.binomial(n0, q))
    c, _ = r14b.GEN["plana"](n, r, r)
    row = {"qi": qi, "q": q, "s": s, **r14b.measure(c, rng_for(ENT_NULO, K_NULO_SEL, qi, s))}
    r19._write(null_path(Path(task["out"]), qi, s), {"spec": task["spec"], "row": row})
    return f"nulo q{qi} s{s} N={n} {time.time() - t0:.1f} s"


def parent_rows(fam: str, p: int, n0: int) -> tuple[list[dict[str, Any]], float]:
    """Mide el padre p de la familia: (filas por q, segundos de generación del padre)."""
    f = FIDX[fam]
    t0 = time.time()
    x = generar(fam, n0, rng_for(ENT_PAD, K_PTS, f, n0, p), rng_for(ENT_PAD, K_AUX, f, n0, p))
    t_x = time.time() - t0
    rows: list[dict[str, Any]] = []
    for qi, q in enumerate(QS):
        _, y = r14b.dilute(x, q, rng_for(ENT_PAD, K_DIL, f, n0, p, qi) if q < 1.0 else None)
        n = int(y.shape[0])
        obs = r14b.measure(y, rng_for(ENT_PAD, K_SEL_Y, f, n0, p, qi))
        del y
        z = generar(fam, n, rng_for(ENT_PAD, K_CTRL_PTS, f, n0, p, qi), rng_for(ENT_PAD, K_CTRL_AUX, f, n0, p, qi))
        if z.shape[0] != n:
            raise RuntimeError(f"control con N = {z.shape[0]} distinto del descendiente {n}")
        ctrl = r14b.measure(z, rng_for(ENT_PAD, K_SEL_Z, f, n0, p, qi))
        del z
        rows.append({"qi": qi, "q": q, "N": n, "obs": obs, "ctrl": ctrl})
    return rows, t_x


def run_parent(task: dict[str, Any]) -> str:
    t0 = time.time()
    fam, p, n0 = task["fam"], task["p"], task["n0"]
    rows, t_x = parent_rows(fam, p, n0)
    r19._write(parent_path(Path(task["out"]), fam, p),
               {"spec": task["spec"], "familia": fam, "padre": p, "t_X_s": round(t_x, 2), "rows": rows})
    return f"padre {fam} p{p} {time.time() - t0:.1f} s (X {t_x:.1f} s)"


def _work(task: dict[str, Any]) -> str:
    return run_null(task) if task["kind"] == "nulo" else run_parent(task)


def _cost(task: dict[str, Any]) -> float:
    return float(QS[task["qi"]] ** 2) if task["kind"] == "nulo" else 3.0


def _pool(tasks: list[dict[str, Any]], workers: int, budget_s: float | None) -> int:
    """Como `r20._pool`, con fork y `_work` de este módulo; no lanza unidades nuevas pasado el presupuesto."""
    t0 = time.time()
    done = 0

    def over() -> bool:
        return budget_s is not None and time.time() - t0 > budget_s

    if workers <= 1:
        for t in tasks:
            if over():
                break
            print(_work(t), flush=True)
            done += 1
        return done
    with mp.get_context("fork").Pool(min(workers, max(1, len(tasks))), maxtasksperchild=4) as pool:
        pending: list[Any] = []
        it = iter(tasks)
        while True:
            while len(pending) < workers and not over():
                nxt: dict[str, Any] | None = next(it, None)
                if nxt is None:
                    break
                pending.append(pool.apply_async(_work, (nxt,)))
            if not pending:
                break
            ready = [p for p in pending if p.ready()]
            for p in ready:
                print(p.get(), flush=True)
                pending.remove(p)
                done += 1
            if not ready:
                time.sleep(0.2)
    return done


def cmd_run(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    override = any(v is not None for v in (args.n0, args.nparents, args.nnull))
    if override and not args.out:
        sys.exit("--n0/--nparents/--nnull exigen --out (no se mezcla con la corrida oficial)")
    with r20._logged(out, "run"):
        _run(args, out, override)


def _run(args: argparse.Namespace, out: Path, override: bool) -> None:
    head, dirty = r20.git_info()
    if dirty and not args.allow_dirty:
        sys.exit("exige arbol git limpio fuera de results/ (usar --allow-dirty solo para pruebas)")
    if not (_ROOT / SPEC).exists():
        sys.exit(f"falta {SPEC}")
    n0 = args.n0 or N0
    npad = args.nparents or N_PADRES
    nnull = args.nnull or N_NULO
    sp = r19.sha256(_ROOT / SPEC)
    official = not dirty and not override and out == OUT
    spec = {"commit": head, "dirty": dirty, "official": official, "n0": n0, "nparents": npad, "nnull": nnull,
            "preregistro_sha256": sp, "experiment_id": EXPERIMENT_ID}
    keys = ("n0", "nparents", "nnull", "official", "preregistro_sha256", "experiment_id")
    manifest = out / "manifest.json"
    if manifest.exists():
        old = json.loads(manifest.read_text())
        if any(old.get(k) != spec[k] for k in keys):
            sys.exit(f"{manifest} tiene otra configuración: {old}")
    else:
        r19._write(manifest, {**{k: spec[k] for k in keys}, "commit_inicial": head,
                              "entropias": {"nulo": ENT_NULO, "padres": ENT_PAD},
                              "familias": {nm: {"indice": i, "papel": pa, "origen": o}
                                           for i, (nm, pa, _, o) in enumerate(FAMILIAS)}})
    tasks: list[dict[str, Any]] = []
    if args.part in ("null", "all"):
        for qi in range(len(QS)):
            for s in range(nnull):
                if not null_path(out, qi, s).exists():
                    tasks.append({"kind": "nulo", "qi": qi, "s": s, "n0": n0, "out": str(out), "spec": spec})
    if args.part in ("parents", "all"):
        for fam in r14b._names(args.family, NOMBRES):
            for p in r14b._range(args.parents, npad):
                if not parent_path(out, fam, p).exists():
                    tasks.append({"kind": "padre", "fam": fam, "p": p, "n0": n0, "out": str(out), "spec": spec})
    tasks.sort(key=lambda t: -_cost(t))
    print(f"{out} (commit {head[:8]}, dirty {dirty}, oficial {official}): {len(tasks)} unidades por correr")
    t0 = time.time()
    done = _pool(tasks, args.workers, args.budget_min * 60 if args.budget_min else None)
    print(f"{done} de {len(tasks)} unidades; pared {time.time() - t0:.1f} s"
          + ("" if done == len(tasks) else "; repetir el comando para seguir"))


# --- Instrumento: banda por par (§2) ----------------------------------------------------------------------


def u_rank(d: float, d_null: Float) -> float:
    """Fracción de la nula con D menor que d; empates por la mitad. nan si d es nan o la nula está vacía."""
    a = np.asarray(d_null, dtype=np.float64)
    a = a[np.isfinite(a)]
    if math.isnan(d) or a.size == 0:
        return math.nan
    return float((np.count_nonzero(a < d) + 0.5 * np.count_nonzero(a == d)) / a.size)


class PairBand:
    """Blanco plano de un q: ajuste (Band), calibración (radio conforme) y conjunto nulo. Banda por subconjunto de
    componentes S (en caché): (Band, radio, D de los nulos en S con nan para los nulos sin alguna componente de S)."""

    def __init__(self, fit: Float, cal: Float, null: Float, comunes: tuple[int, ...] = COMUNES) -> None:
        self.fit, self.cal, self.null = fit, cal, null
        refs = np.vstack([fit, cal])
        self.base = tuple(i for i in comunes if bool(np.all(np.isfinite(refs[:, i]))))
        self.floors = r14b.comp_floors()
        self._cache: dict[tuple[int, ...], tuple[r14b.Band, float, Float]] = {}

    def band(self, comps: tuple[int, ...]) -> tuple[r14b.Band, float, Float]:
        if comps not in self._cache:
            b = r14b.Band(self.fit, self.fit[:0], comps, self.floors)
            rad, _, _ = r20.conformal_radius(b.dist(self.cal[:, list(comps)], sub=True))
            nv = self.null[:, list(comps)]
            ok = np.all(np.isfinite(nv), axis=1)
            dn = np.full(nv.shape[0], np.nan)
            if np.any(ok):
                dn[ok] = b.dist(nv[ok], sub=True) / rad
            self._cache[comps] = (b, rad, dn)
        return self._cache[comps]

    def comps_par(self, vy: Float, vz: Float) -> tuple[int, ...]:
        return tuple(i for i in self.base if math.isfinite(vy[i]) and math.isfinite(vz[i]))

    def pair(self, vy: Float, vz: Float) -> dict[str, Any]:
        """D_Y, D_Z, u_Y, u_Z con las componentes de la base finitas en los dos vectores; nan si |S| < 10."""
        comps = self.comps_par(vy, vz)
        rec: dict[str, Any] = {"S": comps, "dim": len(comps), "evaluable": len(comps) >= MIN_COMPS}
        if not rec["evaluable"]:
            return {**rec, "D_Y": math.nan, "D_Z": math.nan, "u_Y": math.nan, "u_Z": math.nan, "n_nulo": 0, "radio": math.nan}
        b, rad, dn = self.band(comps)
        dy = float(b.dist(vy)[0]) / rad
        dz = float(b.dist(vz)[0]) / rad
        return {**rec, "D_Y": dy, "D_Z": dz, "u_Y": u_rank(dy, dn), "u_Z": u_rank(dz, dn),
                "n_nulo": int(np.count_nonzero(np.isfinite(dn))), "radio": rad}


# --- Estadísticos y clasificación (§5, §6) ----------------------------------------------------------------


def _nanmean_rows(a: Float) -> Float:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.asarray(np.nanmean(a, axis=-1), dtype=np.float64)


def _nanmedian_rows(a: Float) -> Float:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.asarray(np.nanmedian(a, axis=-1), dtype=np.float64)


def boot_ci(x: Float, idx: npt.NDArray[np.int64], stat: str) -> tuple[float, float, float]:
    """(estimación, lo, hi): media o mediana sobre los padres con valor (nan = no evaluable) e IC 95% percentil por
    bootstrap con los índices dados. Si alguna réplica no tiene valores, el IC es nan."""
    a = np.asarray(x, dtype=np.float64)
    f = _nanmean_rows if stat == "media" else _nanmedian_rows
    est = float(f(a[None, :])[0]) if np.any(np.isfinite(a)) else math.nan
    reps = f(a[idx])
    if np.any(np.isnan(reps)) or math.isnan(est):
        return est, math.nan, math.nan
    lo, hi = np.percentile(reps, [2.5, 97.5])
    return est, float(lo), float(hi)


def etiquetas(auc_q1: tuple[float, float], auc_q8: tuple[float, float], delta_q8: tuple[float, float],
              m: float = MARGEN) -> dict[str, bool]:
    """Etiquetas de §6 según los IC (lo, hi); un IC nan no cumple ninguna."""
    def ok(ci: tuple[float, float]) -> bool:
        return not (math.isnan(ci[0]) or math.isnan(ci[1]))

    return {"visible_al_inicio": ok(auc_q1) and auc_q1[0] >= VISIBLE,
            "cerca": ok(auc_q8) and auc_q8[1] <= CERCA,
            "lejos": ok(auc_q8) and auc_q8[0] >= LEJOS,
            "mejora": ok(delta_q8) and delta_q8[1] < -m,
            "igual": ok(delta_q8) and delta_q8[0] >= -m and delta_q8[1] <= m}


def lectura(et: dict[str, bool], frac_no_eval: dict[str, float]) -> dict[str, Any]:
    """Lectura por familia (§6) con la precedencia declarada en el docstring del módulo."""
    tambien: list[str] = []
    if any(v > MAX_NO_EVAL for v in frac_no_eval.values()):
        return {"lectura": INDET, "motivo": f"falta de resolución: no evaluables > {MAX_NO_EVAL:.0%} de los padres "
                f"({frac_no_eval})", "tambien_cumple": tambien}
    if not et["visible_al_inicio"]:
        return {"lectura": NO_APLICA, "motivo": "no visible al inicio", "tambien_cumple": tambien}
    if et["cerca"] and et["mejora"]:
        return {"lectura": COMPATIBLE, "motivo": "visible, cerca y mejora", "tambien_cumple": tambien}
    if et["igual"]:
        donde = "cerca" if et["cerca"] else ("lejos" if et["lejos"] else "ni cerca ni lejos")
        if et["lejos"]:
            tambien.append(SIN)
        return {"lectura": TAMANO, "motivo": f"visible e igual que el control; queda {donde}", "tambien_cumple": tambien}
    if et["lejos"] and not et["mejora"]:
        return {"lectura": SIN, "motivo": "visible, lejos y sin mejora", "tambien_cumple": tambien}
    return {"lectura": INDET, "motivo": "otra combinación de etiquetas", "tambien_cumple": tambien}


# --- Resumen ------------------------------------------------------------------------------------------


def _vecs(rows: list[dict[str, Any]]) -> Float:
    return np.array([r14b.obs_vector(r) for r in rows], dtype=np.float64)


def refs_sha(refs: Path, nref: int) -> str:
    """sha256 de los archivos de referencia usados (ajuste y calibración, en orden fijo)."""
    h = hashlib.sha256()
    for conj in ("ajuste", "calibracion"):
        for qi in range(len(QS)):
            for s in range(nref):
                h.update((refs / conj / f"q{qi}_s{s:03d}.json").read_bytes())
    return h.hexdigest()


def load_bands(out: Path, refs: Path, nref: int, nnull: int) -> tuple[dict[int, PairBand], dict[str, Any], set[tuple[Any, ...]]]:
    bands: dict[int, PairBand] = {}
    info: dict[str, Any] = {}
    prov: set[tuple[Any, ...]] = set()
    for qi in range(len(QS)):
        frows = [r19._load(refs / "ajuste" / f"q{qi}_s{s:03d}.json")["row"] for s in range(nref)]
        fit = _vecs(frows)
        cal = _vecs([r19._load(refs / "calibracion" / f"q{qi}_s{s:03d}.json")["row"] for s in range(nref)])
        nrows = []
        for s in range(nnull):
            d = r19._load(null_path(out, qi, s))
            prov.add(_prov(d))
            nrows.append(d["row"])
        pb = PairBand(fit, cal, _vecs(nrows))
        bands[qi] = pb
        b, rad, dn = pb.band(pb.base)
        info[str(qi)] = {"q": QS[qi], "base": [r14b.COMPONENTS[i] for i in pb.base], "dim_base": len(pb.base),
                         "comunes_no_finitas_en_refs": [r14b.COMPONENTS[i] for i in COMUNES if i not in pb.base],
                         "radio_base": rad, "lambda_lw": b.lam, "con_piso": b.floored,
                         "N_media_ajuste": float(np.mean([r["N"] for r in frows])),
                         "N_media_nulo": float(np.mean([r["N"] for r in nrows])),
                         "nulo_sin_base_completa": int(np.count_nonzero(~np.isfinite(dn))),
                         "nulo_D_mediana": float(np.median(dn[np.isfinite(dn)])) if np.any(np.isfinite(dn)) else math.nan,
                         "nulo_en_banda": float(np.mean(dn[np.isfinite(dn)] <= 1.0)) if np.any(np.isfinite(dn)) else math.nan}
    return bands, info, prov


def _prov(d: dict[str, Any]) -> tuple[str, bool, bool, str]:
    s = d["spec"]
    return s["commit"], s["dirty"], s["official"], s["preregistro_sha256"]


def summarize(out: Path, refs: Path, nref: int = N_REF) -> dict[str, Any]:
    man = r19._load(out / "manifest.json")
    npad, nnull = int(man["nparents"]), int(man["nnull"])
    bands, binfo, prov = load_bands(out, refs, nref, nnull)
    res: dict[str, Any] = {}
    tabla: list[list[Any]] = []
    for fam in NOMBRES:
        dy = np.full((npad, len(QS)), np.nan)
        dz, uy, uz, dims = (np.full_like(dy, np.nan) for _ in range(4))
        nq = np.full_like(dy, np.nan)
        faltan: dict[str, int] = {}
        nnul_excl: dict[str, int] = {}
        for p in range(npad):
            d = r19._load(parent_path(out, fam, p))
            prov.add(_prov(d))
            for row in d["rows"]:
                qi = int(row["qi"])
                pb = bands[qi]
                vy, vz = r14b.obs_vector(row["obs"]), r14b.obs_vector(row["ctrl"])
                for arm, v in (("Y", vy), ("Z", vz)):
                    for j in COMUNES:
                        if not math.isfinite(v[j]):
                            kk = f"q{qi}:{arm}:{r14b.COMPONENTS[j]}"
                            faltan[kk] = faltan.get(kk, 0) + 1
                pr = pb.pair(vy, vz)
                dy[p, qi], dz[p, qi], uy[p, qi], uz[p, qi] = pr["D_Y"], pr["D_Z"], pr["u_Y"], pr["u_Z"]
                dims[p, qi], nq[p, qi] = pr["dim"], float(row["N"])
                if pr["evaluable"] and pr["n_nulo"] < nnull:
                    nnul_excl[str(qi)] = nnul_excl.get(str(qi), 0) + (nnull - pr["n_nulo"])
        lr = r20.log_ratio(dy, dz)
        idx = r20.boot_indices(npad, rng_for(ENT_PAD, K_BOOT, FIDX[fam]))
        por_q: list[dict[str, Any]] = []
        for qi, q in enumerate(QS):
            ev = np.isfinite(uy[:, qi]) & np.isfinite(uz[:, qi]) & np.isfinite(lr[:, qi])
            ay = boot_ci(np.where(ev, uy[:, qi], np.nan), idx, "media")
            az = boot_ci(np.where(ev, uz[:, qi], np.nan), idx, "media")
            de = boot_ci(np.where(ev, lr[:, qi], np.nan), idx, "mediana")
            cell: dict[str, Any] = {"q": q, "N_media": float(np.mean(nq[:, qi])),
                    "AUC_Y": ay[0], "AUC_Y_ic95": [ay[1], ay[2]], "AUC_Z": az[0], "AUC_Z_ic95": [az[1], az[2]],
                    "Delta": de[0], "Delta_ic95": [de[1], de[2]],
                    "D_Y_mediana": float(np.nanmedian(dy[ev, qi])) if ev.any() else math.nan,
                    "D_Z_mediana": float(np.nanmedian(dz[ev, qi])) if ev.any() else math.nan,
                    "Y_en_banda": float(np.mean(dy[ev, qi] <= 1.0)) if ev.any() else math.nan,
                    "Z_en_banda": float(np.mean(dz[ev, qi] <= 1.0)) if ev.any() else math.nan,
                    "no_evaluables": int(np.count_nonzero(~ev)), "frac_no_evaluables": float(np.mean(~ev)),
                    "dim_S": sorted({int(v) for v in dims[:, qi] if math.isfinite(v)}),
                    "nulos_excluidos_total": nnul_excl.get(str(qi), 0)}
            por_q.append(cell)
            tabla.append([fam, PAPEL[fam], q, cell["N_media"], cell["AUC_Y"], *cell["AUC_Y_ic95"], cell["AUC_Z"],
                          *cell["AUC_Z_ic95"], cell["Delta"], *cell["Delta_ic95"], cell["D_Y_mediana"],
                          cell["D_Z_mediana"], cell["Y_en_banda"], cell["Z_en_banda"], cell["no_evaluables"],
                          ";".join(str(v) for v in cell["dim_S"])])
        a1, a8 = por_q[QI_INI], por_q[QI_FIN]
        et = etiquetas((a1["AUC_Y_ic95"][0], a1["AUC_Y_ic95"][1]), (a8["AUC_Y_ic95"][0], a8["AUC_Y_ic95"][1]),
                       (a8["Delta_ic95"][0], a8["Delta_ic95"][1]))
        fne = {"q=1": a1["frac_no_evaluables"], "q=1/8": a8["frac_no_evaluables"]}
        lec = lectura(et, fne)
        res[fam] = {"indice": FIDX[fam], "papel": PAPEL[fam], "etiquetas": et, **lec, "por_q": por_q,
                    "componentes_faltantes": faltan, "D_Y": dy.tolist(), "D_Z": dz.tolist(), "u_Y": uy.tolist(),
                    "u_Z": uz.tolist(), "ln_DY_DZ": lr.tolist()}
    r19._write_csv(out / "tabla.csv", ["familia", "papel", "q", "N_media", "AUC_Y", "AUC_Y_lo", "AUC_Y_hi", "AUC_Z",
                                       "AUC_Z_lo", "AUC_Z_hi", "Delta", "Delta_lo", "Delta_hi", "D_Y_mediana",
                                       "D_Z_mediana", "Y_en_banda", "Z_en_banda", "no_evaluables", "dim_S"], tabla)
    commits = sorted({p[0] for p in prov})
    return {"estado": "PILOTO EXPLORATORIO R23 (sin voto)", "preregistro": str(SPEC),
            "preregistro_sha256_actual": r19.sha256(_ROOT / SPEC), "experiment_id": EXPERIMENT_ID,
            "entropias": {"nulo": ENT_NULO, "padres": ENT_PAD}, "manifest": man, "refs": str(refs), "nref": nref,
            "refs_sha256": refs_sha(refs, nref), "n_boot": N_BOOT, "margen": MARGEN,
            "umbrales": {"visible": VISIBLE, "cerca": CERCA, "lejos": LEJOS, "max_no_evaluables": MAX_NO_EVAL,
                         "min_componentes": MIN_COMPS},
            "procedencia": {"commits": commits, "un_solo_commit": len(commits) == 1, "alguna_sucia": any(p[1] for p in prov),
                            "todas_oficiales": all(p[2] for p in prov), "sha_preregistro": sorted({p[3] for p in prov})},
            "componentes_comunes": [r14b.COMPONENTS[i] for i in COMUNES], "q": list(QS), "blanco": binfo,
            "familias": res, "alcance": "exploratorio; IC condicionales al instrumento y al conjunto nulo; "
                                        "no se infiere velocidad de decaimiento"}


def cmd_summarize(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    refs = Path(args.refs) if args.refs else REFS
    with r20._logged(out, "summarize"):
        s = summarize(out, refs, args.nref or N_REF)
        (out / "summary.json").write_text(json.dumps(_clean(s), indent=2, default=str, ensure_ascii=False))
        _print_summary(s)


def _f(x: float, w: int = 5) -> str:
    return f"{x:{w}.2f}" if math.isfinite(x) else f"{x!s:>{w}s}"


def _print_summary(s: dict[str, Any]) -> None:
    pv = s["procedencia"]
    print(f"commits {pv['commits']} oficial {pv['todas_oficiales']} sucia {pv['alguna_sucia']}")
    print("\nBlanco por q (base, radio, nulos en banda):")
    for qi, b in s["blanco"].items():
        print(f"  q={QS[int(qi)]:<6} d{b['dim_base']} R {b['radio_base']:.3f}  nulo en banda {_f(b['nulo_en_banda'])}  "
              f"faltan en refs {b['comunes_no_finitas_en_refs']}")
    print("\nPor familia y q: AUC_Y [IC]  AUC_Z [IC]  Δ [IC]  medD_Y medD_Z  no_eval")
    for fam, r in s["familias"].items():
        print(f"{fam} ({r['papel']}): {r['lectura']} ({r['motivo']})"
              + (f"; también cumple {r['tambien_cumple']}" if r["tambien_cumple"] else ""))
        for c in r["por_q"]:
            print(f"   q={c['q']:<6} N {c['N_media']:7.0f}  {_f(c['AUC_Y'])} [{_f(c['AUC_Y_ic95'][0])},{_f(c['AUC_Y_ic95'][1])}]"
                  f"  {_f(c['AUC_Z'])} [{_f(c['AUC_Z_ic95'][0])},{_f(c['AUC_Z_ic95'][1])}]"
                  f"  {_f(c['Delta'])} [{_f(c['Delta_ic95'][0])},{_f(c['Delta_ic95'][1])}]"
                  f"  {_f(c['D_Y_mediana'], 7)} {_f(c['D_Z_mediana'], 7)}  {c['no_evaluables']}")
        et = r["etiquetas"]
        print("   etiquetas: " + (", ".join(k for k, v in et.items() if v) or "ninguna"))


# --- Figura -------------------------------------------------------------------------------------------


def cmd_figure(args: argparse.Namespace) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(args.out) if args.out else OUT
    fig_path = Path(args.fig) if args.fig else FIG
    s = r19._load(out / "summary.json")
    fams = list(s["familias"])
    ncol = 4
    nrow = math.ceil(len(fams) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.8 * ncol, 3.2 * nrow + 0.6), squeeze=False)
    xs = np.log2(np.array(QS))
    ticks = [f"1/{round(1 / q)}" if q < 1 else "1" for q in QS]

    def val(x: Any) -> float:
        return math.nan if x is None else float(x)

    for ax, fam in zip(axes.ravel(), fams):
        r = s["familias"][fam]
        for key, lab, col, dx in (("AUC_Y", "descendiente Y", "C0", -0.06), ("AUC_Z", "control nuevo Z", "C1", 0.06)):
            est = np.array([val(c[key]) for c in r["por_q"]])
            lo = np.array([val(c[f"{key}_ic95"][0]) for c in r["por_q"]])
            hi = np.array([val(c[f"{key}_ic95"][1]) for c in r["por_q"]])
            ax.errorbar(xs + dx, est, yerr=[np.abs(est - lo), np.abs(hi - est)], marker="o", ms=3.5, color=col,
                        capsize=2.5, label=lab)
        ax.axhline(LEJOS, color="C3", ls="--", lw=0.8)
        ax.axhline(CERCA, color="C2", ls="--", lw=0.8)
        ax.axhline(0.5, color="0.5", ls=":", lw=0.8)
        ax.set_ylim(0.0, 1.03)
        ax.set_xticks(xs)
        ax.set_xticklabels(ticks)
        ax.set_xlim(xs[0] + 0.3, xs[-1] - 0.3)
        ax.set_title(f"{fam} ({r['papel']})\n{r['lectura']}", fontsize=7.5)
        ax.set_xlabel("q", fontsize=7)
        ax.set_ylabel("AUC frente al nulo plano", fontsize=7)
        ax.grid(alpha=0.3)
        ax.tick_params(labelsize=6.5)
    axes.ravel()[0].legend(fontsize=6.5, loc="lower left")
    for extra in axes.ravel()[len(fams):]:
        extra.axis("off")
    man = s["manifest"]
    fig.suptitle(f"R23 (piloto exploratorio, sin voto): AUC del descendiente y del control nuevo frente a {man['nnull']} "
                 f"nulos planos por q (N0 = {man['n0']}, {man['nparents']} padres; IC 95% bootstrap por padre; "
                 f"líneas 0.75 y 0.95)", fontsize=9)
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=110)
    plt.close(fig)
    print(fig_path)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--part", choices=("null", "parents", "all"), default="all")
    r.add_argument("--family", help="familias separadas por coma (por defecto, todas)")
    r.add_argument("--parents", help="padres a-b")
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--budget-min", type=float, default=110.0, help="no lanza unidades nuevas pasado este tiempo")
    r.add_argument("--n0", type=int)
    r.add_argument("--nparents", type=int)
    r.add_argument("--nnull", type=int)
    r.add_argument("--out")
    r.add_argument("--allow-dirty", action="store_true")
    s = sub.add_parser("summarize")
    s.add_argument("--out")
    s.add_argument("--refs", help=f"referencias de R20 (por defecto {REFS})")
    s.add_argument("--nref", type=int, help=f"referencias por conjunto y por q (por defecto {N_REF})")
    f = sub.add_parser("figure")
    f.add_argument("--out")
    f.add_argument("--fig")
    args = ap.parse_args()
    {"run": cmd_run, "summarize": cmd_summarize, "figure": cmd_figure}[args.cmd](args)


if __name__ == "__main__":
    main()
