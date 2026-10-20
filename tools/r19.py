"""R19-(b): convergencia de defectos bajo dilución con control de defecto reinyectado. Ejecutor del preregistro
confirmatorio docs/OMEGA_R19_GEMELOS_PREREGISTRO.md con la interfaz docs/R19_INTERFAZ_SELLADO.md.

Ejecutor reanudable:
  run [--part refs|structs|all] [--struct a,b] [--real a-b] [--s a-b] [--sellado PATH]
      [--workers 4] [--budget-min 110] [--n0 N --reals R --nref K --nnuevas M --out DIR --allow-dirty]
  summarize [--out DIR] [--sellado PATH]  -> <out>/summary.json, tabla_L.csv, tabla_bandas.csv y el veredicto
  figure [--out DIR] [--fig PNG]          -> docs/figuras/r19.png

Instrumento: el de R14-B v2, importado de `tools/r14b.py` sin copiar: Φ_q (`dilute`), observables (`measure`,
vector p = 14 con √gemelos, `obs_vector`), banda de Mahalanobis con Ledoit-Wolf y radio de la aclaración A1
(`Reference`/`Band`: máx(p95 de las de ajuste, 19.ª de 20 de validación)). D = distancia / radio frente a la
referencia **plana** del mismo q (§2), para toda estructura.

Unidades (un archivo por unidad; las que ya existen se saltan, así que repetir el comando reanuda):
- referencia: <out>/refs/q<qi>_s<ss>.json, diamante plano 4D con N ~ Binomial(N0, q) (N0 en q = 1), como R14-B.
  60 por q: s = 0..39 fijan la banda y 40..59 la validan (regla A1).
- validación nueva: <out>/refs_nuevas/q<qi>_s<ss>.json, 20 diamantes planos más por q, con otra clave de semilla,
  para la condición de calibración de §4 («cubre al menos 16 de 20 sprinklings nuevos»). Se declara así porque las 20
  de validación de la banda entran en el radio (lo cubren con 19/20 por construcción) y no sirven de prueba.
- estructura: <out>/structs/<estructura>/r<rr>.json, una realización: observado y control en los cuatro q.

Estructuras: las selladas (`ESTRUCTURAS` del módulo `--sellado`, votan) y cinco vistas que no votan:
  plana, de_sitter (calibración nula, «region»), D7_cono_huso y D8_cono_suma_normas (conos vistos, «region») y
  D5_gemelos («copias» con k = 1, frac = 0.10: 800 gemelos sobre 7200 de base con N0 = 8000, como
  `f6_intento2.d5_gemelos`, pero construido con el mismo generador genérico de copias que las selladas).

Generador de copias (interfaz, tipo 1, con la regla fijada por el orquestador): con un total N (tamaño del orden
que se construye), n_originales = round(frac·N/k) (round de Python), n_copias = k·n_originales y
n_base = N − n_copias (con k = 1, n_copias = round(frac·N)); base = `sprinkle_diamond(n_base, 4, ·)` con relación de
Minkowski (`causal_blocks`); se eligen n_originales originales distintos de la base (sin reemplazo, azar auxiliar) y
cada uno recibe exactamente k copias con su mismo pasado y futuro; luego las etiquetas se permutan con el azar
auxiliar, como `d5_build`.

Control (§2), declarado:
- region: observado = Φ_q(generar(N0)); control = generar(N_c) con N_c ~ Binomial(N0, q), puntos nuevos.
- copias: observado = Φ_q(X) con X = base de N0 − round(frac·N0) más copias. Control: primero
  N_c ~ Binomial(N0, q) (el total del control tiene la misma distribución que el total de Φ_q(X)); luego
  n_copias = k·round(frac·N_c/k), n_base = N_c − n_copias; la base del control es un subconjunto uniforme de tamaño
  n_base de la base de X (orden inducido), es decir, Φ_q de la base condicionada a su tamaño, que es otro sprinkling
  plano de n_base puntos; por último se reinyectan copias con la misma k y la misma frac (generador de arriba, azar
  auxiliar propio). Así la fracción de copias del control es k·round(frac·N_c/k)/N_c ≈ frac en cada q.
- q = 1: observado y control tienen la misma distribución; se calculan igual, con sorteos independientes (en copias,
  la base es la misma de X y las copias se sortean de nuevo; en region, un sprinkling nuevo).
- D_obs y D_ctrl de un par (realización, q) usan la misma banda: las componentes de la referencia que son finitas en
  los dos vectores (aclaración (c) de R14-B aplicada al par).
- L_q = ln(D_obs(q) / D_ctrl(q)); L = L_{1/8}.

Semillas: SeedSequence(entropía, spawn_key=(332, tipo, ...)).
  20261110 (estructuras y diluciones), con k = índice de la estructura (vistas 0..4; selladas 100 + posición en
  ESTRUCTURAS):
    0 = puntos de la estructura (vistas)          (332, 0, k, N0, real)
    1 = azar auxiliar del observado (copias)      (332, 1, k, N0, real)
    2 = dilución del observado                    (332, 2, k, N0, real, qi)  q = 1 no lo usa
    3 = selección de intervalos del observado     (332, 3, k, N0, real, qi)
    4 = control: Binomial y subconjunto de base   (332, 4, k, N0, real, qi)
    5 = control: puntos nuevos (region, vistas)   (332, 5, k, N0, real, qi)
    6 = control: azar auxiliar de las copias      (332, 6, k, N0, real, qi)
    7 = control: selección de intervalos          (332, 7, k, N0, real, qi)
    8 = bootstrap del resumen                     (332, 8, j)   j = k (estructuras), 9000 + familia
  20261111 (generador sellado): los **puntos** de las estructuras selladas, es decir, el generador que se pasa a
    `generar` (region) o a `sprinkle_diamond` de la base (copias): tipo 0 (observado) y tipo 5 (control region),
    con las mismas claves. Todo lo demás de las selladas (diluciones, copias, Binomial, intervalos) sale de 20261110.
  20261112 (referencias y validación):
    0 = puntos de la referencia (Binomial y puntos)  (332, 0, qi, s)
    1 = selección de intervalos de la referencia     (332, 1, qi, s)
    2 = puntos de la validación nueva                (332, 2, qi, s)
    3 = selección de intervalos de la validación     (332, 3, qi, s)

Árbol limpio: se exige `git status` limpio fuera de results/ y del propio archivo sellado (que el protocolo §8 copia
al repositorio después del commit del ejecutor); su sha256 se registra en cada unidad y se compara con
docs/SELLADO_R19.md si existe.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # un hilo por proceso; el paralelismo viene de --workers

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import multiprocessing as mp
import re
import subprocess
import sys
import time
from collections.abc import Callable
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import numpy.typing as npt

_TOOLS = Path(__file__).resolve().parent
_ROOT = _TOOLS.parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_TOOLS))

import r14b
from es0_piloto import _clean
from f6_exploracion import causal_blocks

from omega.causality.sprinkling import sprinkle_diamond

EXPERIMENT_ID = 332
ENT_EST = 20261110
ENT_SELLADO = 20261111
ENT_REF = 20261112
(K_PTS, K_AUX, K_DIL, K_SEL, K_CTRL, K_CTRL_PTS, K_CTRL_AUX, K_CTRL_SEL, K_BOOT) = range(9)
K_REF_PTS, K_REF_SEL, K_NEW_PTS, K_NEW_SEL = 0, 1, 2, 3
N0 = 8000
QS = r14b.QS  # (1, 1/2, 1/4, 1/8)
N_REAL = 20
N_REFS = 60
N_FIT = 40
N_NUEVAS = 20
MIN_COBERTURA = 16  # de 20 nuevas (§4); con otro número, la misma proporción redondeada hacia arriba
N_BOOT = 2000
MARGEN = math.log(1.3)
SPEC = Path("docs") / "OMEGA_R19_GEMELOS_PREREGISTRO.md"
INTERFAZ = Path("docs") / "R19_INTERFAZ_SELLADO.md"
SELLO_DOC = Path("docs") / "SELLADO_R19.md"
SELLADO = Path("tools") / "sellado_r19.py"
OUT = Path("results") / "omega_r19"
FIG = Path("docs") / "figuras" / "r19.png"
Bool = npt.NDArray[np.bool_]
Int = npt.NDArray[np.int64]
Float = npt.NDArray[np.float64]
Generar = Callable[[int, np.random.Generator], Bool]
TIPOS_FAMILIA = ("se_desvanece", "autosemejante")
FRAC_RANGO = (0.05, 0.25)


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


# --- Generador de copias y control reinyectado ---------------------------------------------------------


def n_copias_de(n_total: int, frac: float, k: int = 1) -> tuple[int, int]:
    """(n_copias, n_base) con un total N: n_originales = round(frac·N/k), n_copias = k·n_originales,
    n_base = N − n_copias (regla de la interfaz fijada por el orquestador; con k = 1, round(frac·N))."""
    nc = k * round(frac * n_total / k)
    return nc, n_total - nc


def add_copies(c: Bool, n_c: int, k: int, aux: np.random.Generator) -> tuple[Bool, dict[str, Any]]:
    """Agrega n_c copias (múltiplo de k) a la base C (n_b × n_b): n_orig = n_c/k originales distintos (sin
    reemplazo), k copias cada uno, con el pasado y el futuro del original; etiquetas permutadas con `aux` (como
    `d5_build`). Devuelve (C', info), con info["origen"][e] = índice en la base del elemento (o de su original) con
    etiqueta e en C', e info["es_copia"][e]. Verificación barata obligatoria: cada copia tiene exactamente las filas y
    columnas de su original y no está relacionada con él."""
    n_b = int(c.shape[0])
    if k < 1:
        raise ValueError("k debe ser >= 1")
    if n_c < 0 or n_c % k:
        raise ValueError("n_c debe ser >= 0 y múltiplo de k")
    n_orig = n_c // k
    if n_orig > n_b:
        raise ValueError(f"faltan originales: {n_orig} > {n_b}")
    orig = np.asarray(aux.choice(n_b, size=n_orig, replace=False), dtype=np.int64)
    src = np.repeat(orig, k)
    n = n_b + n_c
    phi = np.concatenate([np.arange(n_b, dtype=np.int64), src])
    perm = np.asarray(aux.permutation(n), dtype=np.int64)
    origen = phi[perm]
    cp = c[np.ix_(origen, origen)]
    is_copy = perm >= n_b
    pos = np.empty(n, dtype=np.int64)
    pos[perm] = np.arange(n)
    tw = pos[n_b + np.arange(n_c)]
    og = pos[src]
    ok = (np.all(cp[tw] == cp[og], axis=1) & np.all(cp[:, tw] == cp[:, og], axis=0) & ~cp[tw, og] & ~cp[og, tw]) if n_c else np.ones(0, bool)
    built = int(np.count_nonzero(ok))
    info: dict[str, Any] = {"n_base": n_b, "n_copias": n_c, "k": k, "n_originales": n_orig, "copias_verificadas": built,
                            "frac_copias": n_c / n if n else math.nan}
    if built != n_c:
        raise RuntimeError(f"copias: verificación fallida {info}")
    info["origen"] = origen
    info["es_copia"] = is_copy
    return cp, info


def copias_base(n_base: int, rng: np.random.Generator) -> Bool:
    """Base de sprinkling 4D del diamante plano con la relación de Minkowski."""
    return causal_blocks(sprinkle_diamond(n_base, 4, rng), "mink")


def copias_generar(n: int, k: int, frac: float, rng: np.random.Generator, aux: np.random.Generator) -> tuple[Bool, Bool, dict[str, Any]]:
    """X con N = n: (C', base, info)."""
    nc, nb = n_copias_de(n, frac, k)
    base = copias_base(nb, rng)
    cp, info = add_copies(base, nc, k, aux)
    return cp, base, info


def copias_control(base: Bool, n0: int, q: float, k: int, frac: float, rng: np.random.Generator,
                   aux: np.random.Generator) -> tuple[Bool, dict[str, Any]]:
    """Control reinyectado en q: N_c ~ Binomial(N0, q) (N0 en q = 1), n_copias = k·round(frac·N_c/k), base = subconjunto
    uniforme de n_base = N_c − n_copias elementos de la base de X (orden inducido) y copias nuevas con la misma k."""
    n_c = n0 if q >= 1.0 else int(rng.binomial(n0, q))
    nc, nb = n_copias_de(n_c, frac, k)
    if nb > base.shape[0]:
        raise RuntimeError(f"control: n_base {nb} > base de X {base.shape[0]}")
    sub = np.sort(rng.choice(base.shape[0], size=nb, replace=False)).astype(np.int64)
    cp, info = add_copies(base[np.ix_(sub, sub)], nc, k, aux)
    return cp, info


def region_control(generar: Generar, n0: int, q: float, rng_size: np.random.Generator, rng_pts: np.random.Generator) -> Bool:
    """Control de una región: generar(N_c) con N_c ~ Binomial(N0, q) (N0 en q = 1)."""
    n = n0 if q >= 1.0 else int(rng_size.binomial(n0, q))
    return _check_matrix(generar(n, rng_pts), n)


def _check_matrix(c: Any, n: int) -> Bool:
    a = np.asarray(c)
    if a.shape != (n, n) or a.dtype != np.bool_:
        raise TypeError(f"generar debe devolver bool {n}×{n}; devolvió {a.dtype} {a.shape}")
    return a


# --- Estructuras ----------------------------------------------------------------------------------------


def _r14b_region(name: str) -> Generar:
    g = r14b.GEN[name]

    def gen(n: int, rng: np.random.Generator) -> Bool:
        return g(n, rng, rng)[0]

    return gen


# (nombre, grupo, rol, estructura) de las vistas; el índice es k.
VISTAS: list[tuple[str, str, str, dict[str, Any]]] = [
    ("plana", "calibracion", "nula", {"tipo": "region", "generar": _r14b_region("plana"), "origen": "r14b.GEN['plana']"}),
    ("de_sitter", "calibracion", "nula", {"tipo": "region", "generar": _r14b_region("de_sitter"), "origen": "r14b.GEN['de_sitter']"}),
    ("D7_cono_huso", "vista", "autosemejante", {"tipo": "region", "generar": _r14b_region("D7_cono_huso"),
                                               "origen": "f6_intento2._gauge_gen(d7_gauge, BOX_D7)"}),
    ("D8_cono_suma_normas", "vista", "autosemejante", {"tipo": "region", "generar": _r14b_region("D8_cono_suma_normas"),
                                                      "origen": "f6_intento2._gauge_gen(d8_gauge, BOX_STD)"}),
    ("D5_gemelos", "vista", "se_desvanece", {"tipo": "copias", "k": 1, "frac": 0.10,
                                             "origen": "generador genérico de copias con k = 1, frac = 0.10 (≡ D5)"}),
]
CALIBRACION = ("plana", "de_sitter")


def load_sellado(path: Path) -> ModuleType:
    if not path.exists():
        sys.exit(f"falta el archivo sellado {path}")
    spec = importlib.util.spec_from_file_location("sellado_r19", path)
    if spec is None or spec.loader is None:
        sys.exit(f"no se puede cargar {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_SELLADO_CACHE: dict[str, ModuleType] = {}


def _sellado_cached(path: str) -> ModuleType:
    if path not in _SELLADO_CACHE:
        _SELLADO_CACHE[path] = load_sellado(Path(path))
    return _SELLADO_CACHE[path]


def validate_sellado(mod: ModuleType) -> list[str]:
    """Errores de forma -> SystemExit; reglas del guardián que no se cumplen -> advertencias (se registran)."""
    fams = getattr(mod, "FAMILIAS", None)
    ests = getattr(mod, "ESTRUCTURAS", None)
    if not isinstance(fams, list) or not isinstance(ests, dict):
        sys.exit("el sellado debe definir FAMILIAS (list) y ESTRUCTURAS (dict)")
    vistas = {v[0] for v in VISTAS}
    warn: list[str] = []
    for name, e in ests.items():
        if name in vistas:
            sys.exit(f"nombre sellado {name} choca con una vista")
        if not isinstance(e, dict) or e.get("tipo") not in ("copias", "region"):
            sys.exit(f"estructura {name}: tipo inválido")
        if e["tipo"] == "copias":
            if not isinstance(e.get("k"), int) or e["k"] < 1 or not isinstance(e.get("frac"), (int, float)):
                sys.exit(f"estructura {name}: k o frac inválidos")
            if not FRAC_RANGO[0] <= float(e["frac"]) <= FRAC_RANGO[1]:
                warn.append(f"{name}: frac {e['frac']} fuera de [0.05, 0.25]")
        elif not callable(e.get("generar")):
            sys.exit(f"estructura {name}: generar no es invocable")
    seen: dict[str, str] = {}
    for f in fams:
        if not isinstance(f, dict) or f.get("tipo") not in TIPOS_FAMILIA or not isinstance(f.get("estructuras"), list):
            sys.exit(f"familia inválida: {f}")
        for s in f["estructuras"]:
            if s not in ests:
                sys.exit(f"familia {f.get('nombre')}: estructura desconocida {s}")
            if s in seen:
                sys.exit(f"estructura {s} en dos familias")
            seen[s] = f["tipo"]
        if not 2 <= len(f["estructuras"]) <= 3:
            warn.append(f"familia {f.get('nombre')}: {len(f['estructuras'])} estructuras (la regla pide 2 o 3)")
    for t in TIPOS_FAMILIA:
        nt = sum(1 for f in fams if f["tipo"] == t)
        if nt < 2:
            warn.append(f"{nt} familias de tipo {t} (la regla pide al menos 2)")
    for s in ests:
        if s not in seen:
            warn.append(f"estructura {s} sin familia: no vota")
    return warn


def structures(mod: ModuleType) -> list[dict[str, Any]]:
    """Todas las estructuras con su clave k, grupo, rol, familia, tipo y la entropía de sus puntos."""
    out: list[dict[str, Any]] = []
    for k, (name, grupo, rol, e) in enumerate(VISTAS):
        out.append({"nombre": name, "k": k, "grupo": grupo, "rol": rol, "familia": None, "ent_puntos": ENT_EST, **e})
    fam_of = {s: f for f in mod.FAMILIAS for s in f["estructuras"]}
    for j, (name, e) in enumerate(mod.ESTRUCTURAS.items()):
        f = fam_of.get(name)
        out.append({"nombre": name, "k": 100 + j, "grupo": "sellada", "rol": f["tipo"] if f else None,
                    "familia": f["nombre"] if f else None, "ent_puntos": ENT_SELLADO, **e})
    return out


# --- Unidades de trabajo --------------------------------------------------------------------------------


def ref_path(out: Path, qi: int, s: int, nueva: bool = False) -> Path:
    return out / ("refs_nuevas" if nueva else "refs") / f"q{qi}_s{s:02d}.json"


def struct_path(out: Path, name: str, real: int) -> Path:
    return out / "structs" / name / f"r{real:02d}.json"


def _write(p: Path, rec: dict[str, Any]) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(_clean(rec), default=str))
    tmp.replace(p)


def run_ref(task: dict[str, Any]) -> str:
    t0 = time.time()
    qi, s, n0, nueva = task["qi"], task["s"], task["n0"], task["nueva"]
    q = QS[qi]
    kp, ks = (K_NEW_PTS, K_NEW_SEL) if nueva else (K_REF_PTS, K_REF_SEL)
    r = rng_for(ENT_REF, kp, qi, s)
    n = n0 if q >= 1.0 else int(r.binomial(n0, q))
    c, _ = r14b.GEN["plana"](n, r, r)
    row = {"qi": qi, "q": q, "s": s, "nueva": nueva, **r14b.measure(c, rng_for(ENT_REF, ks, qi, s))}
    _write(ref_path(Path(task["out"]), qi, s, nueva), {"spec": task["spec"], "row": row})
    return f"ref{'-nueva' if nueva else ''} q{qi} s{s} N={n} {time.time() - t0:.1f} s"


def _defect_info(info: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in info.items() if k not in ("origen", "es_copia")}


def run_struct(task: dict[str, Any]) -> str:
    t0 = time.time()
    name, real, n0 = task["name"], task["real"], task["n0"]
    st = {s["nombre"]: s for s in structures(_sellado_cached(task["sellado"]))}[name]
    k, ent = st["k"], st["ent_puntos"]
    pts = rng_for(ent, K_PTS, k, n0, real)
    base: Bool | None = None
    extra: dict[str, Any] = {}
    if st["tipo"] == "copias":
        c, base, info = copias_generar(n0, st["k"], float(st["frac"]), pts, rng_for(ENT_EST, K_AUX, k, n0, real))
        extra = _defect_info(info)
    else:
        c = _check_matrix(st["generar"](n0, pts), n0)
    t_c = time.time() - t0
    rows: list[dict[str, Any]] = []
    for qi, q in enumerate(QS):
        _, cd = r14b.dilute(c, q, rng_for(ENT_EST, K_DIL, k, n0, real, qi) if q < 1.0 else None)
        obs = r14b.measure(cd, rng_for(ENT_EST, K_SEL, k, n0, real, qi))
        del cd
        rc = rng_for(ENT_EST, K_CTRL, k, n0, real, qi)
        cinfo: dict[str, Any] = {}
        if st["tipo"] == "copias":
            assert base is not None
            cc, ci = copias_control(base, n0, q, st["k"], float(st["frac"]), rc, rng_for(ENT_EST, K_CTRL_AUX, k, n0, real, qi))
            cinfo = _defect_info(ci)
        else:
            cc = region_control(st["generar"], n0, q, rc, rng_for(ent, K_CTRL_PTS, k, n0, real, qi))
        ctrl = r14b.measure(cc, rng_for(ENT_EST, K_CTRL_SEL, k, n0, real, qi))
        del cc
        rows.append({"qi": qi, "q": q, "obs": obs, "ctrl": ctrl, "ctrl_defecto": cinfo})
    _write(struct_path(Path(task["out"]), name, real), {"spec": task["spec"], "t_C_s": round(t_c, 2), "extra": extra, "rows": rows})
    return f"est {name} r{real} {time.time() - t0:.1f} s (C {t_c:.1f} s)"


def _work(task: dict[str, Any]) -> str:
    return run_ref(task) if task["kind"] == "ref" else run_struct(task)


def _cost(task: dict[str, Any]) -> float:
    if task["kind"] == "ref":
        return float(QS[task["qi"]] ** 2)
    return 3.0


def _pool(tasks: list[dict[str, Any]], workers: int, budget_s: float | None) -> int:
    """Como `r14b._pool`, pero con `_work` de este módulo."""
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
    with mp.get_context("spawn").Pool(min(workers, max(1, len(tasks)))) as pool:
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


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else ""


def git_info(sellado: Path) -> tuple[str, bool]:
    """(HEAD, sucio) con `git status` fuera de results/ y del archivo sellado (si está dentro del repositorio)."""
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=False, cwd=_ROOT).stdout.strip()
    excl = [":!results"]
    try:
        rel = sellado.resolve().relative_to(_ROOT)
        excl.append(f":!{rel.as_posix()}")
    except ValueError:
        pass
    st = subprocess.run(["git", "status", "--porcelain", "--", ".", *excl], capture_output=True, text=True, check=False, cwd=_ROOT)
    return head, bool(st.stdout.strip())


def sello_registrado(sha: str) -> bool | None:
    """¿El sha256 del sellado figura en docs/SELLADO_R19.md? None si ese documento no existe."""
    if not SELLO_DOC.exists():
        return None
    return sha in set(re.findall(r"[0-9a-f]{64}", SELLO_DOC.read_text()))


def cmd_run(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    sellado = Path(args.sellado) if args.sellado else SELLADO
    override = any(v is not None for v in (args.n0, args.reals, args.nref, args.nnuevas))
    if override and not args.out:
        sys.exit("--n0/--reals/--nref/--nnuevas exigen --out (no se mezcla con la corrida oficial)")
    head, dirty = git_info(sellado)
    if dirty and not args.allow_dirty:
        sys.exit("exige arbol git limpio fuera de results/ y del sellado (usar --allow-dirty solo para pruebas)")
    if not SPEC.exists():
        sys.exit(f"falta {SPEC}")
    mod = load_sellado(sellado)
    warn = validate_sellado(mod)
    for w in warn:
        print(f"advertencia: {w}")
    sts = structures(mod)
    names = [s["nombre"] for s in sts]
    n0 = args.n0 or N0
    reals = args.reals or N_REAL
    nref = args.nref or N_REFS
    nnew = args.nnuevas or N_NUEVAS
    sh = sha256(sellado)
    official = not dirty and not override and out == OUT and sellado == SELLADO
    spec = {"commit": head, "dirty": dirty, "official": official, "n0": n0, "reals": reals, "nref": nref, "nnuevas": nnew,
            "preregistro_sha256": sha256(SPEC), "sellado_sha256": sh, "experiment_id": EXPERIMENT_ID}
    keys = ("n0", "reals", "nref", "nnuevas", "official", "sellado_sha256", "experiment_id")
    manifest = out / "manifest.json"
    if manifest.exists():
        old = json.loads(manifest.read_text())
        if any(old.get(k2) != spec[k2] for k2 in keys):
            sys.exit(f"{manifest} tiene otra configuración: {old}")
    else:
        _write(manifest, {**{k2: spec[k2] for k2 in keys}, "sellado": str(sellado),
                          "entropias": {"estructuras": ENT_EST, "sellado": ENT_SELLADO, "referencias": ENT_REF}})
    tasks: list[dict[str, Any]] = []
    if args.part in ("refs", "all"):
        for qi in range(len(QS)):
            for nueva, total in ((False, nref), (True, nnew)):
                for s in _clip(args.s, total):
                    if not ref_path(out, qi, s, nueva).exists():
                        tasks.append({"kind": "ref", "qi": qi, "s": s, "nueva": nueva, "n0": n0, "out": str(out), "spec": spec})
    if args.part in ("structs", "all"):
        for name in r14b._names(args.struct, names):
            for real in r14b._range(args.real, reals):
                if not struct_path(out, name, real).exists():
                    tasks.append({"kind": "est", "name": name, "real": real, "n0": n0, "out": str(out), "spec": spec,
                                  "sellado": str(sellado.resolve())})
    tasks.sort(key=lambda t: -_cost(t))
    print(f"{out} (commit {head[:8]}, dirty {dirty}, oficial {official}, sellado {sh[:12]}): {len(tasks)} unidades por correr")
    t0 = time.time()
    done = _pool(tasks, args.workers, args.budget_min * 60 if args.budget_min else None)
    print(f"{done} de {len(tasks)} unidades; pared {time.time() - t0:.1f} s"
          + ("" if done == len(tasks) else "; repetir el comando para seguir"))


def _clip(spec: str | None, total: int) -> list[int]:
    """--s: rango a-b recortado a 0..total−1 (se aplica a las referencias y a las nuevas)."""
    if spec is None:
        return list(range(total))
    a, _, b = spec.partition("-")
    lo, hi = int(a), int(b) if b else int(a)
    return list(range(lo, min(hi, total - 1) + 1))


# --- Estadístico y reglas de decisión (§2, §4) -------------------------------------------------------


def pair_D(ref: r14b.Reference, vo: Float, vc: Float) -> tuple[float, float, tuple[int, ...]]:
    """(D_obs/radio, D_ctrl/radio, componentes) con la banda de las componentes finitas en ambos vectores."""
    comps = tuple(i for i in ref.base if math.isfinite(vo[i]) and math.isfinite(vc[i]))
    b = ref.band(comps)
    return float(b.dist(vo)[0]) / b.radio, float(b.dist(vc)[0]) / b.radio, comps


def median_ci(x: Float, rng: np.random.Generator, n_boot: int = N_BOOT) -> tuple[float, float, float]:
    """Mediana y IC 95% percentil por bootstrap de realizaciones."""
    a = np.asarray(x, dtype=np.float64)
    a = a[np.isfinite(a)]
    n = a.size
    if n == 0:
        return math.nan, math.nan, math.nan
    boots = np.median(a[rng.integers(0, n, size=(n_boot, n))], axis=1)
    return float(np.median(a)), float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def clase(rol: str | None, lo: float, hi: float, m: float = MARGEN) -> str:
    """Clase de una estructura según el IC de L. se_desvanece: «separa» (IC < −m), «contrario» (IC entero por encima
    de −m) o «ambiguo». autosemejante y nula: «equivalente» (IC ⊂ [−m, m]), «fuera» (IC entero fuera) o «ambiguo»."""
    if not (math.isfinite(lo) and math.isfinite(hi)):
        return "sin datos"
    if rol == "se_desvanece":
        if hi < -m:
            return "separa"
        if lo > -m:
            return "contrario"
        return "ambiguo"
    if lo >= -m and hi <= m:
        return "equivalente"
    if hi < -m or lo > m:
        return "fuera"
    return "ambiguo"


def verdict(calibracion_ok: bool, selladas: dict[str, tuple[str, float, float]], m: float = MARGEN) -> str:
    """§4. selladas: nombre -> (tipo, lo, hi) del IC de L. INDETERMINADO si falla la calibración o no hay selladas;
    PASA si todo «se desvanece» separa y todo «autosemejante» es equivalente; FALLA si alguno queda del lado contrario
    (se desvanece con IC entero por encima de −m, o autosemejante con IC entero fuera); si no, INDETERMINADO."""
    if not calibracion_ok or not selladas:
        return "INDETERMINADO"
    cl = {n: (t, clase(t, lo, hi, m)) for n, (t, lo, hi) in selladas.items()}
    if all((t == "se_desvanece" and c == "separa") or (t == "autosemejante" and c == "equivalente") for t, c in cl.values()):
        return "PASA"
    if any((t == "se_desvanece" and c == "contrario") or (t == "autosemejante" and c == "fuera") for t, c in cl.values()):
        return "FALLA"
    return "INDETERMINADO"


def min_cobertura(nnew: int) -> int:
    return MIN_COBERTURA if nnew == N_NUEVAS else math.ceil(MIN_COBERTURA * nnew / N_NUEVAS)


# --- Resumen ----------------------------------------------------------------------------------------


def _load(p: Path) -> dict[str, Any]:
    if not p.exists():
        raise SystemExit(f"falta {p}")
    d: dict[str, Any] = json.loads(p.read_text())
    return d


def _prov(d: dict[str, Any]) -> tuple[str, bool, bool, str, str]:
    s = d["spec"]
    return s["commit"], s["dirty"], s["official"], s["preregistro_sha256"], s["sellado_sha256"]


def summarize(out: Path, sellado: Path) -> dict[str, Any]:
    man = _load(out / "manifest.json")
    mod = load_sellado(sellado)
    warn = validate_sellado(mod)
    sh = sha256(sellado)
    if sh != man["sellado_sha256"]:
        raise SystemExit(f"el sellado {sellado} ({sh[:12]}) no es el de la corrida ({man['sellado_sha256'][:12]})")
    sts = structures(mod)
    prov: set[tuple[str, bool, bool, str, str]] = set()
    nf = r14b.n_fit_for(man["nref"])  # 40 de 60 en la corrida oficial
    refs: dict[int, r14b.Reference] = {}
    bandas: dict[str, Any] = {}
    t_band: list[list[Any]] = []
    cob_ok = True
    mc = min_cobertura(man["nnuevas"])
    for qi in range(len(QS)):
        rows, news = [], []
        for s in range(man["nref"]):
            d = _load(ref_path(out, qi, s))
            prov.add(_prov(d))
            rows.append(d["row"])
        for s in range(man["nnuevas"]):
            d = _load(ref_path(out, qi, s, True))
            prov.add(_prov(d))
            news.append(d["row"])
        ref = r14b.Reference(np.array([r14b.obs_vector(r) for r in rows]), nf)
        refs[qi] = ref
        b = ref.band(ref.base)
        dn = [ref.normalized(r14b.obs_vector(r))[0] for r in news]
        cub = int(sum(v <= 1.0 for v in dn))
        cob_ok &= cub >= mc
        info = {**b.info(), "q": QS[qi], "N_media": float(np.mean([r["N"] for r in rows])),
                "omitidas": [r14b.COMPONENTS[i] for i in range(len(r14b.COMPONENTS)) if i not in ref.base],
                "nuevas": {"n": len(news), "cubiertas": cub, "minimo": mc, "ok": cub >= mc, "D_sobre_radio": dn}}
        bandas[str(qi)] = info
        t_band.append([QS[qi], info["N_media"], info["dim"], ";".join(info["omitidas"]), info["lambda_lw"], info["radio_ajuste_p95"],
                       info["cobertura_validacion"], info["radio_banda"], info["cobertura_final"], cub, len(news)])
    res: dict[str, Any] = {}
    t_L: list[list[Any]] = []
    for st in sts:
        name = st["nombre"]
        dobs = np.full((man["reals"], len(QS)), np.nan)
        dctl = np.full_like(dobs, np.nan)
        omit: dict[str, int] = {}
        nobs: list[list[int]] = [[] for _ in QS]
        nctl: list[list[int]] = [[] for _ in QS]
        fctl: list[list[float]] = [[] for _ in QS]
        tw_o: list[list[int]] = [[] for _ in QS]
        tw_c: list[list[int]] = [[] for _ in QS]
        for real in range(man["reals"]):
            d = _load(struct_path(out, name, real))
            prov.add(_prov(d))
            for row in sorted(d["rows"], key=lambda r: int(r["qi"])):
                qi = int(row["qi"])
                ref = refs[qi]
                a, c, comps = pair_D(ref, r14b.obs_vector(row["obs"]), r14b.obs_vector(row["ctrl"]))
                dobs[real, qi], dctl[real, qi] = a, c
                for j in ref.base:
                    if j not in comps:
                        key = f"{qi}:{r14b.COMPONENTS[j]}"
                        omit[key] = omit.get(key, 0) + 1
                nobs[qi].append(row["obs"]["N"])
                nctl[qi].append(row["ctrl"]["N"])
                tw_o[qi].append(row["obs"]["pares_gemelos"])
                tw_c[qi].append(row["ctrl"]["pares_gemelos"])
                if row["ctrl_defecto"]:
                    fctl[qi].append(row["ctrl_defecto"]["frac_copias"])
        lq = np.log(dobs / dctl)
        rng = rng_for(ENT_EST, K_BOOT, st["k"])
        por_q: list[dict[str, Any]] = []
        for qi, q in enumerate(QS):
            med, lo, hi = median_ci(lq[:, qi], rng)
            mo, mcx = float(np.nanmedian(dobs[:, qi])), float(np.nanmedian(dctl[:, qi]))
            cell = {"q": q, "L_mediana": med, "L_ic95": [lo, hi], "D_obs_mediana": mo, "D_ctrl_mediana": mcx,
                    "N_obs_media": float(np.mean(nobs[qi])), "N_ctrl_media": float(np.mean(nctl[qi])),
                    "pares_gemelos_obs_media": float(np.mean(tw_o[qi])), "pares_gemelos_ctrl_media": float(np.mean(tw_c[qi])),
                    "frac_copias_ctrl_media": float(np.mean(fctl[qi])) if fctl[qi] else None,
                    "obs_en_banda": float(np.mean(dobs[:, qi] <= 1.0)), "ctrl_en_banda": float(np.mean(dctl[:, qi] <= 1.0))}
            por_q.append(cell)
            t_L.append([name, st["grupo"], st["familia"] or "", st["rol"] or "", st["tipo"], q, cell["N_obs_media"], cell["N_ctrl_media"],
                        mo, mcx, med, lo, hi, cell["frac_copias_ctrl_media"]])
        lo8, hi8 = por_q[3]["L_ic95"]
        res[name] = {"grupo": st["grupo"], "rol": st["rol"], "familia": st["familia"], "tipo": st["tipo"],
                     "k": st.get("k") if st["tipo"] == "copias" else None, "frac": st.get("frac"),
                     "L": {"mediana": por_q[3]["L_mediana"], "ic95": [lo8, hi8]},
                     "clase": clase(st["rol"], lo8, hi8), "por_q": por_q,
                     "D_obs": dobs.tolist(), "D_ctrl": dctl.tolist(), "L_q": lq.tolist(), "componentes_omitidas": omit}
    cal = {n: res[n]["clase"] == "equivalente" for n in CALIBRACION}
    calibracion: dict[str, Any] = {"L_dentro_del_margen": cal, "cobertura_nuevas": {str(qi): bandas[str(qi)]["nuevas"] for qi in range(len(QS))},
                   "cobertura_ok": cob_ok, "ok": all(cal.values()) and cob_ok}
    selladas = {st["nombre"]: (st["rol"], *res[st["nombre"]]["L"]["ic95"]) for st in sts if st["grupo"] == "sellada" and st["rol"]}
    ver = verdict(calibracion["ok"], selladas)
    familias: dict[str, Any] = {}
    for fi, f in enumerate(mod.FAMILIAS):
        mats = [np.array(res[s]["L_q"])[:, 3] for s in f["estructuras"]]
        r = rng_for(ENT_EST, K_BOOT, 9000 + fi)
        point = float(np.median([np.nanmedian(m) for m in mats]))
        boots = np.empty(N_BOOT)
        for ib in range(N_BOOT):
            pick = r.integers(0, len(mats), size=len(mats))
            boots[ib] = np.median([np.nanmedian(mats[j][r.integers(0, mats[j].size, size=mats[j].size)]) for j in pick])
        familias[f["nombre"]] = {"tipo": f["tipo"], "estructuras": f["estructuras"],
                                 "clases": {s: res[s]["clase"] for s in f["estructuras"]},
                                 "L_mediana_de_medianas": point,
                                 "L_ic95_jerarquico": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
                                 "nota": "informativo; el veredicto es por estructura (§4)"}
    _write_csv(out / "tabla_bandas.csv", ["q", "N_media", "dim", "omitidas", "lambda_lw", "radio_ajuste_p95", "cobertura_validacion",
                                          "radio_banda", "cobertura_final", "nuevas_cubiertas", "nuevas"], t_band)
    _write_csv(out / "tabla_L.csv", ["estructura", "grupo", "familia", "rol", "tipo", "q", "N_obs_media", "N_ctrl_media",
                                     "D_obs_mediana", "D_ctrl_mediana", "L_mediana", "L_ic95_lo", "L_ic95_hi", "frac_copias_ctrl"], t_L)
    commits = sorted({p[0] for p in prov})
    registrado = sello_registrado(sh)
    return {"estado": "CONFIRMATORIA R19-(b)", "preregistro": str(SPEC), "preregistro_sha256_actual": sha256(SPEC),
            "interfaz": str(INTERFAZ), "sellado": str(sellado), "sellado_sha256": sh, "sello_registrado_en_doc": registrado,
            "experiment_id": EXPERIMENT_ID, "entropias": {"estructuras": ENT_EST, "sellado": ENT_SELLADO, "referencias": ENT_REF},
            "manifest": man, "n_ajuste": nf, "margen": MARGEN, "n_boot": N_BOOT, "advertencias_interfaz": warn,
            "procedencia": {"commits": commits, "un_solo_commit": len(commits) == 1, "alguna_sucia": any(p[1] for p in prov),
                            "todas_oficiales": all(p[2] for p in prov), "sha_preregistro": sorted({p[3] for p in prov}),
                            "sha_sellado": sorted({p[4] for p in prov})},
            "componentes": list(r14b.COMPONENTS), "q": list(QS), "bandas": bandas, "calibracion": calibracion,
            "estructuras": res, "familias": familias, "veredicto": ver,
            "votan": sorted(selladas), "no_votan": [s["nombre"] for s in sts if s["nombre"] not in selladas]}


def _write_csv(p: Path, header: list[str], rows: list[list[Any]]) -> None:
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def cmd_summarize(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    man = _load(out / "manifest.json")
    sellado = Path(args.sellado) if args.sellado else Path(man["sellado"])
    s = summarize(out, sellado)
    (out / "summary.json").write_text(json.dumps(_clean(s), indent=2, default=str, ensure_ascii=False))
    pv = s["procedencia"]
    print(f"commits {pv['commits']} oficial {pv['todas_oficiales']} sucia {pv['alguna_sucia']} "
          f"sellado {s['sellado_sha256'][:12]} registrado {s['sello_registrado_en_doc']}")
    for w in s["advertencias_interfaz"]:
        print(f"advertencia: {w}")
    print("\nBanda plana por q (dim, radio, cobertura de las nuevas):")
    for qi, b in s["bandas"].items():
        nv = b["nuevas"]
        print(f"  q={QS[int(qi)]:<6} N≈{b['N_media']:7.0f} d{b['dim']} R {b['radio_banda']:.2f}  nuevas {nv['cubiertas']}/{nv['n']}"
              f" (mín {nv['minimo']}) {'ok' if nv['ok'] else 'NO'}")
    print(f"\nL_q = ln(D_obs/D_ctrl), mediana [IC 95%]; margen ±{MARGEN:.3f}")
    print(f"  {'estructura':24s} {'grupo':12s} {'rol':14s} " + " ".join(f"{'q=' + str(q):>22s}" for q in QS) + "  clase(L_1/8)")
    for name, info in s["estructuras"].items():
        cells = " ".join(f"{c['L_mediana']:6.2f} [{c['L_ic95'][0]:5.2f},{c['L_ic95'][1]:5.2f}]" for c in info["por_q"])
        print(f"  {name:24s} {info['grupo']:12s} {info['rol']!s:14s} {cells}  {info['clase']}")
    print("\nD_obs / D_ctrl (medianas) por q:")
    for name, info in s["estructuras"].items():
        print(f"  {name:24s} " + "  ".join(f"{c['D_obs_mediana']:6.2f}/{c['D_ctrl_mediana']:6.2f}" for c in info["por_q"]))
    print("\nFamilias:")
    for fn, f in s["familias"].items():
        print(f"  {fn} ({f['tipo']}): {f['clases']}  L med {f['L_mediana_de_medianas']:.2f} IC {f['L_ic95_jerarquico']}")
    cal = s["calibracion"]
    print(f"\nCalibración: L plana/de Sitter dentro {cal['L_dentro_del_margen']}; cobertura nuevas ok {cal['cobertura_ok']}")
    print(f"VEREDICTO (solo selladas: {', '.join(s['votan']) or 'ninguna'}): {s['veredicto']}")


# --- Figura -------------------------------------------------------------------------------------------


def cmd_figure(args: argparse.Namespace) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(args.out) if args.out else OUT
    fig_path = Path(args.fig) if args.fig else FIG
    s = _load(out / "summary.json")
    names = list(s["estructuras"])
    npan = len(names) + 1
    ncol = 4
    nrow = math.ceil(npan / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 3.3 * nrow), squeeze=False)
    xs = np.log2(np.array(QS))
    n0 = s["manifest"]["n0"]
    ticks = [f"1/{round(1 / q)}" if q < 1 else "1" for q in QS]
    for ax, name in zip(axes.ravel(), names):
        info = s["estructuras"][name]
        do, dc = np.array(info["D_obs"], dtype=float), np.array(info["D_ctrl"], dtype=float)
        for arr, lab, col in ((do, "observado Φ_q(X)", "C0"), (dc, "control", "C1")):
            med = np.nanmedian(arr, axis=0)
            ax.plot(xs, med, marker="o", color=col, label=lab)
            ax.fill_between(xs, np.nanpercentile(arr, 25, axis=0), np.nanpercentile(arr, 75, axis=0), color=col, alpha=0.18)
        ax.axhline(1.0, color="k", ls="--", lw=0.8)
        ax.set_yscale("log")
        ax.set_xticks(xs)
        ax.set_xticklabels(ticks)
        ax.set_xlim(xs[0] + 0.3, xs[-1] - 0.3)
        tag = f"{info['grupo']}, {info['rol']}" + (f", {info['familia']}" if info["familia"] else "")
        ax.set_title(f"{name}\n({tag})", fontsize=8)
        ax.set_xlabel("q", fontsize=8)
        ax.set_ylabel("D / radio (plana)", fontsize=8)
        ax.grid(alpha=0.3, which="both")
        ax.tick_params(labelsize=7)
    axes.ravel()[0].legend(fontsize=7, loc="best")
    ax = axes.ravel()[len(names)]
    m = s["margen"]
    ys = np.arange(len(names))[::-1]
    ax.axvspan(-m, m, color="0.85", label="margen ±ln 1.3")
    for y, name in zip(ys, names):
        info = s["estructuras"][name]
        med, (lo, hi) = info["L"]["mediana"], info["L"]["ic95"]
        col = {"sellada": "C3", "calibracion": "C2", "vista": "C7"}[info["grupo"]]
        ax.errorbar([med], [y], xerr=[[med - lo], [hi - med]], fmt="o", color=col, capsize=3, ms=4)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks(ys)
    ax.set_yticklabels(names, fontsize=6)
    ax.set_xlabel("L = ln(D_obs/D_ctrl) en q = 1/8", fontsize=8)
    ax.set_title(f"L_1/8: mediana e IC 95% ({s['n_boot']} bootstrap)\nrojo sellada, verde calibración, gris vista", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=7, loc="best")
    for extra in axes.ravel()[npan:]:
        extra.axis("off")
    nr = s["manifest"]["reals"]
    fig.suptitle(f"R19-(b): dilución con control reinyectado (N0 = {n0}, {nr} realizaciones; mediana y rango intercuartil). "
                 f"Veredicto: {s['veredicto']}", fontsize=10)
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=120)
    print(fig_path)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--part", choices=("refs", "structs", "all"), default="all")
    r.add_argument("--struct", help="estructuras separadas por coma (por defecto, todas)")
    r.add_argument("--real", help="realizaciones a-b de las estructuras")
    r.add_argument("--s", help="sprinklings a-b de las referencias (y de las nuevas, recortado)")
    r.add_argument("--sellado", help=f"módulo sellado (por defecto {SELLADO})")
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--budget-min", type=float, help="no lanza unidades nuevas pasado este tiempo")
    r.add_argument("--n0", type=int)
    r.add_argument("--reals", type=int)
    r.add_argument("--nref", type=int)
    r.add_argument("--nnuevas", type=int)
    r.add_argument("--out")
    r.add_argument("--allow-dirty", action="store_true")
    s = sub.add_parser("summarize")
    s.add_argument("--out")
    s.add_argument("--sellado", help="por defecto, el registrado en el manifest")
    f = sub.add_parser("figure")
    f.add_argument("--out")
    f.add_argument("--fig")
    args = ap.parse_args()
    {"run": cmd_run, "summarize": cmd_summarize, "figure": cmd_figure}[args.cmd](args)


if __name__ == "__main__":
    main()
