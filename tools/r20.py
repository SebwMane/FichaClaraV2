"""R20: repetición corregida de R19-(b), alcance acotado (opción A). Ejecutor del preregistro confirmatorio
docs/OMEGA_R20_REPETICION_PREREGISTRO.md (versión 2, firmada), con la interfaz docs/R19_INTERFAZ_SELLADO.md más
`IDENTIFICABLES` (§8).

Ejecutor reanudable:
  run [--part refs|structs|all] [--struct a,b] [--real a-b] [--s a-b] [--sellado PATH]
      [--workers 4] [--budget-min 110] [--n0 N --reals R --nref K --out DIR --allow-dirty]
  summarize [--out DIR] [--sellado PATH]  -> <out>/summary.json, tabla_L.csv, tabla_bandas.csv, tabla_familias.csv
  figure [--out DIR] [--fig PNG]          -> docs/figuras/r20.png
Cada `run` y `summarize` deja su log en <out>/logs/ y las versiones (python, numpy, scipy) en <out>/entorno.json.

Lo que se reutiliza de `tools/r19.py` (importado, sin copiar ni modificar): el generador de copias (`n_copias_de`,
`add_copies`, `copias_generar`), la comprobación de matrices, las cinco estructuras vistas (`VISTAS`) y la escritura
atómica. Del instrumento de R14-B v2 (`tools/r14b.py`): Φ_q (`dilute`), observables (`measure`), el vector p = 14
(`obs_vector`) y la estandarización con Ledoit-Wolf (`Band`, solo media, escala, pisos y covarianza contraída).

Decisiones declaradas (todas siguen el preregistro; donde el texto no fija un detalle, se fija aquí):

Banda por q (§2 v2). Tres conjuntos disjuntos de diamantes planos 4D con N ~ Binomial(N0, q) (N0 con q = 1):
  ajuste (100): media, escala con pisos y covarianza de Ledoit-Wolf de `r14b.Band` (la regla A1 de R14-B **no** se usa);
  calibración (100): radio conforme = estadístico de orden k = ⌈0.95·(n+1)⌉ (96 con n = 100) de las distancias de
    calibración; los empates se ordenan por el índice de la referencia (el valor del radio no cambia);
  validación (100): solo la regla de cobertura: al menos 87 de 100 con D = distancia/radio ≤ 1, en cada uno de los 4 q
    (con otro n, ⌈0.87·n⌉). Nunca recalibra.
Componentes (como R14-B (c) y R19): en cada q, la banda usa las componentes finitas en **las 300** referencias de ese q
  (en R19, r_40 faltaba con q = 1/8). Las mismas componentes valen para ajuste, calibración, validación y estructuras.
  Si un brazo de una estructura no tiene medible alguna de esas componentes, su D no está definido (nan): §4 lo trata
  como «L no definido» y la estructura queda no evaluable en ese q. No se reduce la banda por par (a diferencia de R19).

Control emparejado por N (§1 v2). N = tamaño de Φ_q(X) ya diluido (N0 con q = 1).
  copias: n_copias = k·round(frac·N/k), base = subconjunto uniforme de tamaño N − n_copias de la base de la misma
    realización (orden inducido), copias nuevas con la misma k y frac y un generador auxiliar aparte. Con q = 1 la base
    es entera (mismo tamaño) y las copias se sortean de nuevo.
  region: generar(N, rng nuevo).

Estadístico (§4 v2). L_q = ln(D_obs/D_ctrl), en reales extendidos: D_obs = 0 < D_ctrl da −∞, D_ctrl = 0 < D_obs da +∞.
  D_obs = D_ctrl = 0 o algún D nan da L nan. Mediana e IC 95% percentil (2000 réplicas, remuestreo de realizaciones
  completas: los mismos índices para los 4 q) con cuantiles lineales (como numpy) extendidos: si el interpolante
  tiene peso positivo sobre un infinito, el cuantil es ese infinito; −∞ y +∞ con pesos positivos dan nan.
  Una estructura es **no evaluable** en q si algún L_q es nan o si su mediana o su IC quedan nan; la evaluabilidad
  que cuenta para el veredicto es la de q = 1/8 (el estadístico primario). Los IC son condicionales al instrumento.

Población (§3): votan solo las selladas con familia e `IDENTIFICABLES[nombre]["identificable"] is True`. Las demás
  selladas se miden igual y se informan como «no identificable por este instrumento». Mínimo: al menos 2 familias de
  cada tipo con al menos 2 identificables cada una; si no, INDETERMINADO (población insuficiente).
Veredicto (§4): calibración = regla de cobertura en los 4 q y plana y de Sitter evaluables con IC de L_1/8 dentro de
  ±ln 1.3; si falla, INDETERMINADO. Luego población. PASA si todas las votantes son evaluables y cumplen; FALLA
  (contradicción observada, sin control del error familiar) si alguna votante evaluable queda del lado contrario;
  si no, INDETERMINADO. Una votante no evaluable impide PASA.
Familias: informe por familia (identificables, excluidas, clases, mediana de las medianas de L_1/8 de las votantes
  evaluables); informativo, sin IC propio.

Unidades (un archivo por unidad; las que existen se saltan, así que repetir el comando reanuda):
  <out>/refs/<conjunto>/q<qi>_s<sss>.json  conjunto ∈ ajuste, calibracion, validacion; s = 0..99
  <out>/structs/<estructura>/r<rr>.json    una realización: observado y control en los cuatro q

Semillas: SeedSequence(entropía, spawn_key=(333, tipo, ...)).
  20261113 (vistas, diluciones, controles y bootstrap), k = índice de semilla de la estructura (vistas 0..4; selladas
  100 + posición en ESTRUCTURAS). Se guarda en la clave «idx», separada de la multiplicidad k de las copias: en
  `r19.structures` la k de una sellada de copias pisaba el índice (`{"k": 100 + j, **e}`); aquí no.
    0 = puntos de la estructura (vistas)          (333, 0, k, N0, real)
    1 = azar auxiliar del observado (copias)      (333, 1, k, N0, real)
    2 = dilución del observado                    (333, 2, k, N0, real, qi)  q = 1 no lo usa
    3 = selección de intervalos del observado     (333, 3, k, N0, real, qi)
    4 = control copias: subconjunto de la base    (333, 4, k, N0, real, qi)
    5 = control region: puntos nuevos (vistas)    (333, 5, k, N0, real, qi)
    6 = control copias: azar auxiliar             (333, 6, k, N0, real, qi)
    7 = control: selección de intervalos          (333, 7, k, N0, real, qi)
    8 = bootstrap                                 (333, 8, k)
  20261115 (selladas): los puntos de las selladas, tipos 0 (base de copias o generar(N0)) y 5 (generar(N) del control),
    con las mismas claves. Todo lo demás de las selladas sale de 20261113.
  20261116 (referencias planas): (333, t, c, qi, s) con t = 0 puntos (Binomial y puntos), t = 1 selección de
    intervalos, c = 0 ajuste, 1 calibración, 2 validación.
  20261114 es del guardián y aquí no se usa.

Árbol limpio: `git status` limpio fuera de results/ (el sellado ya está en un commit propio según §6). Se registran
commit, sucio, oficial, sha256 del preregistro y del sellado (y si figura en docs/SELLADO_R20.md).
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # un hilo por proceso; el paralelismo viene de --workers

import argparse
import contextlib
import importlib.util
import json
import math
import multiprocessing as mp
import platform
import re
import subprocess
import sys
import time
import warnings
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any, TextIO

import numpy as np
import numpy.typing as npt

_TOOLS = Path(__file__).resolve().parent
_ROOT = _TOOLS.parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_TOOLS))

import r14b
import r19
from es0_piloto import _clean

EXPERIMENT_ID = 333
ENT_EST = 20261113
ENT_SELLADO = 20261115
ENT_REF = 20261116
(K_PTS, K_AUX, K_DIL, K_SEL, K_CTRL, K_CTRL_PTS, K_CTRL_AUX, K_CTRL_SEL, K_BOOT) = range(9)
K_REF_PTS, K_REF_SEL = 0, 1
CONJUNTOS = ("ajuste", "calibracion", "validacion")
N0 = 8000
QS = r14b.QS  # (1, 1/2, 1/4, 1/8)
QI_PRIM = 3  # q = 1/8
N_REAL = 30
N_REF = 100  # por conjunto y por q
MIN_VALIDAS = 87  # de 100 (§2 v2)
NIVEL = 0.95
N_BOOT = 2000
MARGEN = math.log(1.3)
SPEC = Path("docs") / "OMEGA_R20_REPETICION_PREREGISTRO.md"
INTERFAZ = Path("docs") / "R19_INTERFAZ_SELLADO.md"
SELLO_DOC = Path("docs") / "SELLADO_R20.md"
SELLADO = Path("tools") / "sellado_r20.py"
OUT = Path("results") / "omega_r20"
FIG = Path("docs") / "figuras" / "r20.png"
Bool = npt.NDArray[np.bool_]
Float = npt.NDArray[np.float64]
TIPOS_FAMILIA = r19.TIPOS_FAMILIA  # ("se_desvanece", "autosemejante")
TIPO_ESPERADO = {"copias": "se_desvanece", "region": "autosemejante"}
MAX_UNIVERSO = 16
MIN_FAMILIAS_UNIVERSO = 4
MIN_CANDIDATAS = 3
MIN_FAM_VOTO = 2
MIN_IDENT_FAM = 2
CALIBRACION = r19.CALIBRACION  # ("plana", "de_sitter")
NO_IDENT = "no identificable por este instrumento"


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


# --- Control emparejado por N (§1 v2) -----------------------------------------------------------------


def copias_control_n(base: Bool, n: int, k: int, frac: float, rng: np.random.Generator,
                     aux: np.random.Generator) -> tuple[Bool, dict[str, Any]]:
    """Control de copias con total N: n_copias = k·round(frac·N/k); base = subconjunto uniforme de tamaño N − n_copias de
    la base de la realización (orden inducido); copias nuevas con la misma k (azar `aux`)."""
    nc, nb = r19.n_copias_de(n, frac, k)
    if nb > base.shape[0]:
        raise RuntimeError(f"control: n_base {nb} > base de X {base.shape[0]}")
    sub = np.sort(rng.choice(base.shape[0], size=nb, replace=False)).astype(np.int64)
    return r19.add_copies(base[np.ix_(sub, sub)], nc, k, aux)


def region_control_n(generar: r19.Generar, n: int, rng_pts: np.random.Generator) -> Bool:
    """Control de una región con total N: generar(N, rng nuevo)."""
    return r19._check_matrix(generar(n, rng_pts), n)


# --- Banda conforme (§2 v2) -------------------------------------------------------------------------


def conformal_k(n_cal: int, nivel: float = NIVEL) -> int:
    """Rango del radio: ⌈nivel·(n+1)⌉ (96 con n = 100). Se calcula en enteros para evitar errores de redondeo."""
    num = round(nivel * 1000)
    return -(-num * (n_cal + 1) // 1000)


def conformal_radius(d_cal: Float, nivel: float = NIVEL) -> tuple[float, int, int]:
    """(radio, k, índice de la referencia en el corte): el k-ésimo menor de las distancias de calibración, con los
    empates ordenados por índice de referencia (determinista)."""
    d = np.asarray(d_cal, dtype=np.float64)
    n = int(d.size)
    k = conformal_k(n, nivel)
    if k > n or not np.all(np.isfinite(d)):
        raise ValueError(f"radio conforme: k = {k} con n = {n} (o distancias no finitas)")
    order = np.lexsort((np.arange(n), d))
    j = int(order[k - 1])
    return float(d[j]), k, j


def min_validas(n_val: int) -> int:
    return MIN_VALIDAS if n_val == N_REF else math.ceil(MIN_VALIDAS * n_val / N_REF)


class ConformalBand:
    """Banda de la plana en un q: estandarización y Ledoit-Wolf con las de ajuste (`r14b.Band` sin validación), radio
    conforme con las de calibración y regla de cobertura con las de validación. Componentes = finitas en todas."""

    def __init__(self, fit: Float, cal: Float, val: Float) -> None:
        allv = np.vstack([fit, cal, val])
        self.comps = tuple(int(i) for i in np.flatnonzero(np.all(np.isfinite(allv), axis=0)))
        self.band = r14b.Band(fit, fit[:0], self.comps, r14b.comp_floors())
        self.d_cal = self.band.dist(cal[:, list(self.comps)], sub=True)
        self.radio, self.k, self.j_corte = conformal_radius(self.d_cal)
        self.D_val = self.band.dist(val[:, list(self.comps)], sub=True) / self.radio
        self.n_val = int(self.D_val.size)
        self.cubiertas = int(np.count_nonzero(self.D_val <= 1.0))
        self.minimo = min_validas(self.n_val)
        self.ok = self.cubiertas >= self.minimo

    def D(self, vec: Float) -> float:
        """distancia/radio; nan si el vector no tiene medible alguna componente de la banda."""
        v = np.asarray(vec, dtype=np.float64)
        if not np.all(np.isfinite(v[list(self.comps)])):
            return math.nan
        return float(self.band.dist(v)[0]) / self.radio

    def info(self) -> dict[str, Any]:
        return {"componentes": [r14b.COMPONENTS[i] for i in self.comps], "dim": len(self.comps),
                "omitidas": [c for i, c in enumerate(r14b.COMPONENTS) if i not in self.comps],
                "lambda_lw": self.band.lam, "con_piso": self.band.floored, "radio": self.radio, "k": self.k,
                "n_cal": int(self.d_cal.size), "referencia_en_el_corte": self.j_corte,
                "validacion": {"n": self.n_val, "cubiertas": self.cubiertas, "minimo": self.minimo, "ok": self.ok,
                               "D_sobre_radio": self.D_val.tolist()}}


# --- Estadístico en reales extendidos (§4 v2) -----------------------------------------------------------


def log_ratio(d_obs: Float, d_ctrl: Float) -> Float:
    """L = ln(D_obs/D_ctrl) con −∞ (D_obs = 0 < D_ctrl), +∞ (D_ctrl = 0 < D_obs) y nan (0/0 o algún nan)."""
    a = np.asarray(d_obs, dtype=np.float64)
    b = np.asarray(d_ctrl, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        out: Float = np.log(a) - np.log(b)
    out = np.where((a == 0) & (b == 0), np.nan, out)
    return out


def ext_quantile_sorted(s: Float, p: float) -> Float:
    """Cuantil lineal (el de numpy) de cada fila de `s` (ya ordenada en la última dimensión, sin nan) en reales
    extendidos: con peso positivo sobre un infinito el resultado es ese infinito; −∞ y +∞ con pesos positivos, nan."""
    a = np.atleast_2d(np.asarray(s, dtype=np.float64))
    n = a.shape[-1]
    h = (n - 1) * p
    lo = math.floor(h)
    hi = min(lo + 1, n - 1)
    t = h - lo
    x, y = a[:, lo], a[:, hi]
    if t == 0.0:
        return np.array(x, dtype=np.float64)
    fin = np.isfinite(x) & np.isfinite(y)
    with np.errstate(invalid="ignore"):
        lin = x + (y - x) * t
    out = np.where(fin, lin, np.where(~np.isfinite(x) & ~np.isfinite(y), np.where(x == y, x, np.nan),
                                      np.where(np.isfinite(x), y, x)))
    return np.asarray(out, dtype=np.float64)


def ext_median(x: Float) -> float:
    a = np.asarray(x, dtype=np.float64)
    if a.size == 0 or np.any(np.isnan(a)):
        return math.nan
    return float(ext_quantile_sorted(np.sort(a), 0.5)[0])


def boot_indices(n: int, rng: np.random.Generator, n_boot: int = N_BOOT) -> npt.NDArray[np.int64]:
    """Índices de realizaciones para el bootstrap (los mismos para los 4 q: realizaciones completas)."""
    return np.asarray(rng.integers(0, n, size=(n_boot, n)), dtype=np.int64)


def median_ci_ext(x: Float, idx: npt.NDArray[np.int64]) -> tuple[float, float, float]:
    """Mediana y IC 95% percentil por bootstrap con los índices dados; nan si algún valor es nan."""
    a = np.asarray(x, dtype=np.float64)
    if a.size == 0 or np.any(np.isnan(a)):
        return math.nan, math.nan, math.nan
    meds = ext_quantile_sorted(np.sort(a[idx], axis=1), 0.5)
    if np.any(np.isnan(meds)):
        return ext_median(a), math.nan, math.nan
    sm = np.sort(meds)[None, :]
    return ext_median(a), float(ext_quantile_sorted(sm, 0.025)[0]), float(ext_quantile_sorted(sm, 0.975)[0])


def evaluable(med: float, lo: float, hi: float) -> bool:
    return not (math.isnan(med) or math.isnan(lo) or math.isnan(hi))


def clase(rol: str | None, lo: float, hi: float, m: float = MARGEN) -> str:
    """Clase según el IC de L (admite ±∞). se_desvanece: «separa» (hi < −m), «contrario» (lo > −m) o «ambiguo».
    autosemejante y nula: «equivalente» (IC ⊂ [−m, m]), «fuera» (IC entero fuera) o «ambiguo». nan: «no evaluable»."""
    if math.isnan(lo) or math.isnan(hi):
        return "no evaluable"
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


def poblacion(familias: list[dict[str, Any]], identificables: dict[str, bool]) -> dict[str, Any]:
    """§3: por tipo, familias con al menos 2 identificables; ok si hay al menos 2 de cada tipo."""
    por_tipo: dict[str, list[str]] = {t: [] for t in TIPOS_FAMILIA}
    for f in familias:
        n_id = sum(1 for s in f["estructuras"] if identificables.get(s, False))
        if n_id >= MIN_IDENT_FAM:
            por_tipo[f["tipo"]].append(f["nombre"])
    ok = all(len(v) >= MIN_FAM_VOTO for v in por_tipo.values())
    return {"familias_suficientes": por_tipo, "ok": ok}


def verdict(calibracion_ok: bool, poblacion_ok: bool, votantes: dict[str, tuple[str, float, float, bool]],
            m: float = MARGEN) -> tuple[str, str]:
    """§4. votantes: nombre -> (tipo, lo, hi, evaluable). Devuelve (veredicto, motivo)."""
    if not calibracion_ok:
        return "INDETERMINADO", "calibración"
    if not poblacion_ok or not votantes:
        return "INDETERMINADO", "población insuficiente"
    cl = {n: (t, clase(t, lo, hi, m) if ev else "no evaluable") for n, (t, lo, hi, ev) in votantes.items()}
    if any((t == "se_desvanece" and c == "contrario") or (t == "autosemejante" and c == "fuera") for t, c in cl.values()):
        return "FALLA", "contradicción observada en alguna votante (sin control del error familiar)"
    if all((t == "se_desvanece" and c == "separa") or (t == "autosemejante" and c == "equivalente") for t, c in cl.values()):
        return "PASA", "todas las votantes cumplen (intersección-unión)"
    if any(c == "no evaluable" for _, c in cl.values()):
        return "INDETERMINADO", "alguna votante no evaluable"
    return "INDETERMINADO", "alguna votante ambigua"


# --- Sellado ------------------------------------------------------------------------------------------


def load_sellado(path: Path) -> ModuleType:
    if not path.exists():
        sys.exit(f"falta el archivo sellado {path}")
    spec = importlib.util.spec_from_file_location("sellado_r20", path)
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
    """Errores de forma -> SystemExit; reglas del guardián (§3) que no se cumplen -> advertencias registradas."""
    fams = getattr(mod, "FAMILIAS", None)
    ests = getattr(mod, "ESTRUCTURAS", None)
    ident = getattr(mod, "IDENTIFICABLES", None)
    if not isinstance(fams, list) or not isinstance(ests, dict) or not isinstance(ident, dict):
        sys.exit("el sellado debe definir FAMILIAS (list), ESTRUCTURAS (dict) e IDENTIFICABLES (dict)")
    vistas = {v[0] for v in r19.VISTAS}
    warn: list[str] = []
    for name, e in ests.items():
        if name in vistas:
            sys.exit(f"nombre sellado {name} choca con una vista")
        if not isinstance(e, dict) or e.get("tipo") not in ("copias", "region"):
            sys.exit(f"estructura {name}: tipo inválido")
        if e["tipo"] == "copias":
            if not isinstance(e.get("k"), int) or e["k"] < 1 or not isinstance(e.get("frac"), (int, float)):
                sys.exit(f"estructura {name}: k o frac inválidos")
            if not r19.FRAC_RANGO[0] <= float(e["frac"]) <= r19.FRAC_RANGO[1]:
                warn.append(f"{name}: frac {e['frac']} fuera de [0.05, 0.25]")
        elif not callable(e.get("generar")):
            sys.exit(f"estructura {name}: generar no es invocable")
        rec = ident.get(name)
        if not isinstance(rec, dict) or not isinstance(rec.get("identificable"), bool):
            sys.exit(f"IDENTIFICABLES[{name}] falta o no tiene 'identificable' bool")
    extra = sorted(set(ident) - set(ests))
    if extra:
        sys.exit(f"IDENTIFICABLES nombra estructuras que no están en ESTRUCTURAS: {extra}")
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
            if TIPO_ESPERADO[ests[s]["tipo"]] != f["tipo"]:
                warn.append(f"{s}: tipo {ests[s]['tipo']} en una familia {f['tipo']}")
        if len(f["estructuras"]) < MIN_CANDIDATAS:
            warn.append(f"familia {f.get('nombre')}: {len(f['estructuras'])} candidatas (§3 pide al menos {MIN_CANDIDATAS})")
    if len(ests) > MAX_UNIVERSO:
        warn.append(f"universo de {len(ests)} estructuras (§3: a lo sumo {MAX_UNIVERSO})")
    if len(fams) < MIN_FAMILIAS_UNIVERSO:
        warn.append(f"{len(fams)} familias (§3: al menos {MIN_FAMILIAS_UNIVERSO})")
    for s in ests:
        if s not in seen:
            warn.append(f"estructura {s} sin familia: no vota")
    return warn


def structures(mod: ModuleType) -> list[dict[str, Any]]:
    """Todas las estructuras con su clave k, grupo, rol, familia, tipo, entropía de sus puntos y estado de voto."""
    out: list[dict[str, Any]] = []
    for k, (name, grupo, rol, e) in enumerate(r19.VISTAS):
        out.append({"nombre": name, "idx": k, "grupo": grupo, "rol": rol, "familia": None, "ent_puntos": ENT_EST,
                    "identificable": None, "vota": False, "estado": grupo, **e})
    fam_of = {s: f for f in mod.FAMILIAS for s in f["estructuras"]}
    for j, (name, e) in enumerate(mod.ESTRUCTURAS.items()):
        f = fam_of.get(name)
        ident = bool(mod.IDENTIFICABLES[name]["identificable"])
        vota = ident and f is not None
        estado = "vota" if vota else (NO_IDENT if not ident else "sin familia")
        out.append({"nombre": name, "idx": 100 + j, "grupo": "sellada", "rol": f["tipo"] if f else None,
                    "familia": f["nombre"] if f else None, "ent_puntos": ENT_SELLADO, "identificable": ident,
                    "vota": vota, "estado": estado, **e})
    return out


# --- Unidades de trabajo --------------------------------------------------------------------------------


def ref_path(out: Path, conj: int, qi: int, s: int) -> Path:
    return out / "refs" / CONJUNTOS[conj] / f"q{qi}_s{s:03d}.json"


def struct_path(out: Path, name: str, real: int) -> Path:
    return out / "structs" / name / f"r{real:02d}.json"


def run_ref(task: dict[str, Any]) -> str:
    t0 = time.time()
    conj, qi, s, n0 = task["conj"], task["qi"], task["s"], task["n0"]
    q = QS[qi]
    r = rng_for(ENT_REF, K_REF_PTS, conj, qi, s)
    n = n0 if q >= 1.0 else int(r.binomial(n0, q))
    c, _ = r14b.GEN["plana"](n, r, r)
    row = {"conjunto": CONJUNTOS[conj], "qi": qi, "q": q, "s": s, **r14b.measure(c, rng_for(ENT_REF, K_REF_SEL, conj, qi, s))}
    r19._write(ref_path(Path(task["out"]), conj, qi, s), {"spec": task["spec"], "row": row})
    return f"ref {CONJUNTOS[conj]} q{qi} s{s} N={n} {time.time() - t0:.1f} s"


def run_struct(task: dict[str, Any]) -> str:
    t0 = time.time()
    name, real, n0 = task["name"], task["real"], task["n0"]
    st = {s["nombre"]: s for s in structures(_sellado_cached(task["sellado"]))}[name]
    k, ent = st["idx"], st["ent_puntos"]
    pts = rng_for(ent, K_PTS, k, n0, real)
    base: Bool | None = None
    extra: dict[str, Any] = {}
    if st["tipo"] == "copias":
        c, base, info = r19.copias_generar(n0, st["k"], float(st["frac"]), pts, rng_for(ENT_EST, K_AUX, k, n0, real))
        extra = r19._defect_info(info)
    else:
        c = r19._check_matrix(st["generar"](n0, pts), n0)
    t_c = time.time() - t0
    rows: list[dict[str, Any]] = []
    for qi, q in enumerate(QS):
        _, cd = r14b.dilute(c, q, rng_for(ENT_EST, K_DIL, k, n0, real, qi) if q < 1.0 else None)
        n = int(cd.shape[0])
        obs = r14b.measure(cd, rng_for(ENT_EST, K_SEL, k, n0, real, qi))
        del cd
        cinfo: dict[str, Any] = {}
        if st["tipo"] == "copias":
            assert base is not None
            cc, ci = copias_control_n(base, n, st["k"], float(st["frac"]), rng_for(ENT_EST, K_CTRL, k, n0, real, qi),
                                      rng_for(ENT_EST, K_CTRL_AUX, k, n0, real, qi))
            cinfo = r19._defect_info(ci)
        else:
            cc = region_control_n(st["generar"], n, rng_for(ent, K_CTRL_PTS, k, n0, real, qi))
        if cc.shape[0] != n:
            raise RuntimeError(f"control con N = {cc.shape[0]} distinto del observado {n}")
        ctrl = r14b.measure(cc, rng_for(ENT_EST, K_CTRL_SEL, k, n0, real, qi))
        del cc
        rows.append({"qi": qi, "q": q, "N": n, "obs": obs, "ctrl": ctrl, "ctrl_defecto": cinfo})
    r19._write(struct_path(Path(task["out"]), name, real), {"spec": task["spec"], "t_C_s": round(t_c, 2), "extra": extra, "rows": rows})
    return f"est {name} r{real} {time.time() - t0:.1f} s (C {t_c:.1f} s)"


def _work(task: dict[str, Any]) -> str:
    return run_ref(task) if task["kind"] == "ref" else run_struct(task)


def _cost(task: dict[str, Any]) -> float:
    return float(QS[task["qi"]] ** 2) if task["kind"] == "ref" else 3.0


def _pool(tasks: list[dict[str, Any]], workers: int, budget_s: float | None) -> int:
    """Como `r19._pool`, con `_work` de este módulo; no lanza unidades nuevas pasado el presupuesto."""
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


# --- Procedencia y logs -----------------------------------------------------------------------------------


def git_info() -> tuple[str, bool]:
    """(HEAD, sucio) con `git status` fuera de results/."""
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=False, cwd=_ROOT).stdout.strip()
    st = subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True,
                        check=False, cwd=_ROOT)
    return head, bool(st.stdout.strip())


def sello_registrado(sha: str) -> bool | None:
    """¿El sha256 del sellado figura en docs/SELLADO_R20.md? None si ese documento no existe."""
    p = _ROOT / SELLO_DOC
    if not p.exists():
        return None
    return sha in set(re.findall(r"[0-9a-f]{64}", p.read_text()))


def entorno() -> dict[str, Any]:
    try:
        import scipy

        sv: str | None = scipy.__version__
    except ImportError:  # pragma: no cover
        sv = None
    return {"python": sys.version, "numpy": np.__version__, "scipy": sv, "plataforma": platform.platform(),
            "maquina": platform.machine(), "cpus": os.cpu_count(), "argv": sys.argv}


class _Tee:
    def __init__(self, *streams: TextIO) -> None:
        self.streams = streams

    def write(self, s: str) -> int:
        for st in self.streams:
            st.write(s)
        return len(s)

    def flush(self) -> None:
        for st in self.streams:
            st.flush()


@contextlib.contextmanager
def _logged(out: Path, cmd: str) -> Iterator[None]:
    """Copia la salida a <out>/logs/<cmd>_<UTC>.log y registra el entorno (también en <out>/entorno.json)."""
    logs = out / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    head, dirty = git_info()
    env = {**entorno(), "cmd": cmd, "utc": ts, "commit": head, "dirty": dirty}
    (logs / f"entorno_{cmd}_{ts}.json").write_text(json.dumps(env, indent=2))
    (out / "entorno.json").write_text(json.dumps(env, indent=2))
    old = sys.stdout
    with (logs / f"{cmd}_{ts}.log").open("a") as fh:
        sys.stdout = _Tee(old, fh)
        try:
            yield
        finally:
            sys.stdout = old


def cmd_run(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    override = any(v is not None for v in (args.n0, args.reals, args.nref))
    if override and not args.out:
        sys.exit("--n0/--reals/--nref exigen --out (no se mezcla con la corrida oficial)")
    with _logged(out, "run"):
        _run(args, out, override)


def _run(args: argparse.Namespace, out: Path, override: bool) -> None:
    sellado = Path(args.sellado) if args.sellado else SELLADO
    head, dirty = git_info()
    if dirty and not args.allow_dirty:
        sys.exit("exige arbol git limpio fuera de results/ (usar --allow-dirty solo para pruebas)")
    if not (_ROOT / SPEC).exists():
        sys.exit(f"falta {SPEC}")
    mod = load_sellado(sellado)
    warn = validate_sellado(mod)
    for w in warn:
        print(f"advertencia: {w}")
    sts = structures(mod)
    names = [s["nombre"] for s in sts]
    n0 = args.n0 or N0
    reals = args.reals or N_REAL
    nref = args.nref or N_REF
    sh = r19.sha256(sellado)
    sp = r19.sha256(_ROOT / SPEC)
    official = not dirty and not override and out == OUT and sellado == SELLADO
    spec = {"commit": head, "dirty": dirty, "official": official, "n0": n0, "reals": reals, "nref": nref,
            "preregistro_sha256": sp, "sellado_sha256": sh, "experiment_id": EXPERIMENT_ID}
    keys = ("n0", "reals", "nref", "official", "preregistro_sha256", "sellado_sha256", "experiment_id")
    manifest = out / "manifest.json"
    if manifest.exists():
        old = json.loads(manifest.read_text())
        if any(old.get(k2) != spec[k2] for k2 in keys):
            sys.exit(f"{manifest} tiene otra configuración: {old}")
    else:
        r19._write(manifest, {**{k2: spec[k2] for k2 in keys}, "sellado": str(sellado), "commit_inicial": head,
                              "entropias": {"estructuras": ENT_EST, "sellado": ENT_SELLADO, "referencias": ENT_REF},
                              "sello_registrado_en_doc": sello_registrado(sh), "advertencias_interfaz": warn})
    tasks: list[dict[str, Any]] = []
    if args.part in ("refs", "all"):
        for conj in range(len(CONJUNTOS)):
            for qi in range(len(QS)):
                for s in r19._clip(args.s, nref):
                    if not ref_path(out, conj, qi, s).exists():
                        tasks.append({"kind": "ref", "conj": conj, "qi": qi, "s": s, "n0": n0, "out": str(out), "spec": spec})
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


# --- Resumen ----------------------------------------------------------------------------------------


def _prov(d: dict[str, Any]) -> tuple[str, bool, bool, str, str]:
    s = d["spec"]
    return s["commit"], s["dirty"], s["official"], s["preregistro_sha256"], s["sellado_sha256"]


def _fmean(v: list[float]) -> float | None:
    return float(np.mean(v)) if v else None


def summarize(out: Path, sellado: Path) -> dict[str, Any]:
    man = r19._load(out / "manifest.json")
    mod = load_sellado(sellado)
    warn = validate_sellado(mod)
    sh = r19.sha256(sellado)
    if sh != man["sellado_sha256"]:
        raise SystemExit(f"el sellado {sellado} ({sh[:12]}) no es el de la corrida ({man['sellado_sha256'][:12]})")
    sts = structures(mod)
    prov: set[tuple[str, bool, bool, str, str]] = set()
    nref, reals = int(man["nref"]), int(man["reals"])
    bands: dict[int, ConformalBand] = {}
    bandas: dict[str, Any] = {}
    t_band: list[list[Any]] = []
    for qi in range(len(QS)):
        sets: list[Float] = []
        nmed: list[float] = []
        for conj in range(len(CONJUNTOS)):
            rows = []
            for s in range(nref):
                d = r19._load(ref_path(out, conj, qi, s))
                prov.add(_prov(d))
                rows.append(d["row"])
            sets.append(np.array([r14b.obs_vector(r) for r in rows]))
            nmed.append(float(np.mean([r["N"] for r in rows])))
        cb = ConformalBand(*sets)
        bands[qi] = cb
        info = {**cb.info(), "q": QS[qi], "N_media": dict(zip(CONJUNTOS, nmed))}
        bandas[str(qi)] = info
        t_band.append([QS[qi], nmed[0], info["dim"], ";".join(info["omitidas"]), info["lambda_lw"], info["k"],
                       info["referencia_en_el_corte"], info["radio"], cb.cubiertas, cb.n_val, cb.minimo, cb.ok])
    cobertura_ok = all(b.ok for b in bands.values())
    res: dict[str, Any] = {}
    t_L: list[list[Any]] = []
    for st in sts:
        name = st["nombre"]
        dobs = np.full((reals, len(QS)), np.nan)
        dctl = np.full_like(dobs, np.nan)
        nobs: list[list[float]] = [[] for _ in QS]
        fctl: list[list[float]] = [[] for _ in QS]
        tw_o: list[list[float]] = [[] for _ in QS]
        tw_c: list[list[float]] = [[] for _ in QS]
        no_medible: dict[str, int] = {}
        for real in range(reals):
            d = r19._load(struct_path(out, name, real))
            prov.add(_prov(d))
            for row in sorted(d["rows"], key=lambda r: int(r["qi"])):
                qi = int(row["qi"])
                cb = bands[qi]
                for arm, key in (("obs", "observado"), ("ctrl", "control")):
                    v = r14b.obs_vector(row[arm])
                    for j in cb.comps:
                        if not math.isfinite(v[j]):
                            kk = f"q{qi}:{key}:{r14b.COMPONENTS[j]}"
                            no_medible[kk] = no_medible.get(kk, 0) + 1
                dobs[real, qi] = cb.D(r14b.obs_vector(row["obs"]))
                dctl[real, qi] = cb.D(r14b.obs_vector(row["ctrl"]))
                nobs[qi].append(float(row["N"]))
                tw_o[qi].append(float(row["obs"]["pares_gemelos"]))
                tw_c[qi].append(float(row["ctrl"]["pares_gemelos"]))
                if row["ctrl_defecto"]:
                    fctl[qi].append(float(row["ctrl_defecto"]["frac_copias"]))
        lq = log_ratio(dobs, dctl)
        idx = boot_indices(reals, rng_for(ENT_EST, K_BOOT, st["idx"]))
        por_q: list[dict[str, Any]] = []
        for qi, q in enumerate(QS):
            med, lo, hi = median_ci_ext(lq[:, qi], idx)
            ev = evaluable(med, lo, hi)
            cell = {"q": q, "L_mediana": med, "L_ic95": [lo, hi], "evaluable": ev,
                    "clase": clase(st["rol"], lo, hi) if ev else "no evaluable",
                    "D_obs_mediana": ext_median(dobs[:, qi]), "D_ctrl_mediana": ext_median(dctl[:, qi]),
                    "n_L_menos_inf": int(np.sum(lq[:, qi] == -np.inf)), "n_L_mas_inf": int(np.sum(lq[:, qi] == np.inf)),
                    "n_L_nan": int(np.sum(np.isnan(lq[:, qi]))), "N_media": _fmean(nobs[qi]),
                    "pares_gemelos_obs_media": _fmean(tw_o[qi]), "pares_gemelos_ctrl_media": _fmean(tw_c[qi]),
                    "frac_copias_ctrl_media": _fmean(fctl[qi]),
                    "obs_en_banda": float(np.mean(dobs[:, qi] <= 1.0)), "ctrl_en_banda": float(np.mean(dctl[:, qi] <= 1.0))}
            por_q.append(cell)
            t_L.append([name, st["grupo"], st["estado"], st["familia"] or "", st["rol"] or "", st["tipo"], q, cell["N_media"],
                        cell["D_obs_mediana"], cell["D_ctrl_mediana"], med, lo, hi, ev, cell["clase"],
                        cell["frac_copias_ctrl_media"]])
        prim = por_q[QI_PRIM]
        ident_rec = mod.IDENTIFICABLES.get(name) if st["grupo"] == "sellada" else None
        res[name] = {"grupo": st["grupo"], "estado": st["estado"], "vota": st["vota"], "rol": st["rol"],
                     "familia": st["familia"], "tipo": st["tipo"], "identificable": st["identificable"],
                     "requisitos_guardian": ident_rec,
                     "k": st.get("k") if st["tipo"] == "copias" else None, "frac": st.get("frac"),
                     "L": {"mediana": prim["L_mediana"], "ic95": prim["L_ic95"], "evaluable": prim["evaluable"]},
                     "clase": prim["clase"], "por_q": por_q, "componentes_no_medibles": no_medible,
                     "D_obs": dobs.tolist(), "D_ctrl": dctl.tolist(), "L_q": lq.tolist()}
    cal_L = {n: res[n]["L"]["evaluable"] and res[n]["clase"] == "equivalente" for n in CALIBRACION}
    calibracion: dict[str, Any] = {"cobertura": {str(qi): bandas[str(qi)]["validacion"] | {"D_sobre_radio": None} for qi in range(len(QS))},
                   "cobertura_ok": cobertura_ok, "L_dentro_del_margen": cal_L,
                   "ok": cobertura_ok and all(cal_L.values())}
    identificables = {s["nombre"]: bool(s["identificable"]) for s in sts if s["grupo"] == "sellada"}
    pob = poblacion(mod.FAMILIAS, identificables)
    votantes = {s["nombre"]: (s["rol"], *res[s["nombre"]]["L"]["ic95"], res[s["nombre"]]["L"]["evaluable"])
                for s in sts if s["vota"]}
    ver, motivo = verdict(calibracion["ok"], pob["ok"], votantes)
    familias: dict[str, Any] = {}
    t_f: list[list[Any]] = []
    for f in mod.FAMILIAS:
        ss = f["estructuras"]
        ids = [s for s in ss if identificables.get(s)]
        exc = [s for s in ss if not identificables.get(s)]
        meds = [res[s]["L"]["mediana"] for s in ids if res[s]["L"]["evaluable"]]
        mm = ext_median(np.array(meds)) if meds else math.nan
        fam = {"tipo": f["tipo"], "estructuras": ss, "identificables": ids, "n_excluidas": len(exc),
               "excluidas_no_identificables": exc, "clases": {s: res[s]["clase"] for s in ss},
               "votantes_no_evaluables": [s for s in ids if not res[s]["L"]["evaluable"]],
               "L_mediana_de_medianas_votantes": mm, "nota": "informativo; el veredicto es por estructura (§4)"}
        familias[f["nombre"]] = fam
        t_f.append([f["nombre"], f["tipo"], len(ss), len(ids), len(exc), ";".join(exc), mm,
                    ";".join(f"{s}:{res[s]['clase']}" for s in ss)])
    r19._write_csv(out / "tabla_bandas.csv", ["q", "N_media_ajuste", "dim", "omitidas", "lambda_lw", "k", "referencia_en_el_corte",
                                              "radio", "validas_cubiertas", "validas", "minimo", "ok"], t_band)
    r19._write_csv(out / "tabla_L.csv", ["estructura", "grupo", "estado", "familia", "rol", "tipo", "q", "N_media", "D_obs_mediana",
                                         "D_ctrl_mediana", "L_mediana", "L_ic95_lo", "L_ic95_hi", "evaluable", "clase",
                                         "frac_copias_ctrl"], t_L)
    r19._write_csv(out / "tabla_familias.csv", ["familia", "tipo", "candidatas", "identificables", "excluidas", "lista_excluidas",
                                                "L_mediana_de_medianas_votantes", "clases"], t_f)
    commits = sorted({p[0] for p in prov})
    return {"estado": "CONFIRMATORIA R20", "preregistro": str(SPEC), "preregistro_sha256_actual": r19.sha256(_ROOT / SPEC),
            "interfaz": str(INTERFAZ), "sellado": str(sellado), "sellado_sha256": sh, "sello_registrado_en_doc": sello_registrado(sh),
            "experiment_id": EXPERIMENT_ID, "entropias": {"estructuras": ENT_EST, "sellado": ENT_SELLADO, "referencias": ENT_REF},
            "manifest": man, "margen": MARGEN, "n_boot": N_BOOT, "nivel_conforme": NIVEL, "advertencias_interfaz": warn,
            "procedencia": {"commits": commits, "un_solo_commit": len(commits) == 1, "alguna_sucia": any(p[1] for p in prov),
                            "todas_oficiales": all(p[2] for p in prov), "sha_preregistro": sorted({p[3] for p in prov}),
                            "sha_sellado": sorted({p[4] for p in prov})},
            "componentes": list(r14b.COMPONENTS), "q": list(QS), "bandas": bandas, "calibracion": calibracion,
            "poblacion": pob, "estructuras": res, "familias": familias, "veredicto": ver, "motivo": motivo,
            "votan": sorted(votantes), "n_votan": len(votantes),
            "no_identificables": sorted(n for n, v in identificables.items() if not v),
            "no_votan": [s["nombre"] for s in sts if not s["vota"]],
            "alcance": "condicional al instrumento ajustado; solo estructuras identificables; sin control del error familiar en FALLA"}


def cmd_summarize(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    man = r19._load(out / "manifest.json")
    sellado = Path(args.sellado) if args.sellado else Path(man["sellado"])
    with _logged(out, "summarize"):
        s = summarize(out, sellado)
        (out / "summary.json").write_text(json.dumps(_clean(s), indent=2, default=str, ensure_ascii=False))
        _print_summary(s)


def _f(x: float) -> str:
    return f"{x:6.2f}" if math.isfinite(x) else f"{x!s:>6s}"


def _print_summary(s: dict[str, Any]) -> None:
    pv = s["procedencia"]
    print(f"commits {pv['commits']} oficial {pv['todas_oficiales']} sucia {pv['alguna_sucia']} "
          f"sellado {s['sellado_sha256'][:12]} registrado {s['sello_registrado_en_doc']}")
    for w in s["advertencias_interfaz"]:
        print(f"advertencia: {w}")
    print("\nBanda conforme por q (dim, radio, validación):")
    for qi, b in s["bandas"].items():
        v = b["validacion"]
        print(f"  q={QS[int(qi)]:<6} d{b['dim']} R {b['radio']:.3f} (k {b['k']}, ref {b['referencia_en_el_corte']})  "
              f"validación {v['cubiertas']}/{v['n']} (mín {v['minimo']}) {'ok' if v['ok'] else 'NO'}")
    print(f"\nL_q = ln(D_obs/D_ctrl), mediana [IC 95%]; margen ±{MARGEN:.3f}")
    for name, info in s["estructuras"].items():
        cells = " ".join(f"{_f(c['L_mediana'])} [{_f(c['L_ic95'][0])},{_f(c['L_ic95'][1])}]" for c in info["por_q"])
        print(f"  {name:26s} {info['estado'][:16]:16s} {info['rol']!s:14s} {cells}  {info['clase']}")
    print("\nFamilias:")
    for fn, f in s["familias"].items():
        print(f"  {fn} ({f['tipo']}): identificables {f['identificables']}; excluidas {f['n_excluidas']} "
              f"{f['excluidas_no_identificables']}; clases {f['clases']}")
    cal = s["calibracion"]
    print(f"\nCalibración: cobertura ok {cal['cobertura_ok']}; L plana/de Sitter dentro {cal['L_dentro_del_margen']}")
    print(f"Población: {s['poblacion']}")
    print(f"VEREDICTO ({s['n_votan']} votantes: {', '.join(s['votan']) or 'ninguna'}): {s['veredicto']} ({s['motivo']})")


# --- Figura -------------------------------------------------------------------------------------------


def _ext(x: Any) -> float:
    if x is None:
        return math.nan
    if isinstance(x, str):
        return math.inf if x == "inf" else -math.inf
    return float(x)


def _arr(a: Any) -> Float:
    return np.array([[_ext(v) for v in row] for row in a], dtype=np.float64)


ESTILO = {"vota": ("C3", "o", "vota"), NO_IDENT: ("C1", "s", "no identificable"), "sin familia": ("C1", "s", "sin familia"),
          "calibracion": ("C2", "D", "calibración"), "vista": ("C7", "^", "vista, no vota")}


def cmd_figure(args: argparse.Namespace) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(args.out) if args.out else OUT
    fig_path = Path(args.fig) if args.fig else FIG
    s = r19._load(out / "summary.json")
    names = list(s["estructuras"])
    npan = len(names) + 1
    ncol = 5
    nrow = math.ceil(npan / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.7 * ncol, 3.1 * nrow + 0.6), squeeze=False)
    xs = np.log2(np.array(QS))
    ticks = [f"1/{round(1 / q)}" if q < 1 else "1" for q in QS]
    for ax, name in zip(axes.ravel(), names):
        info = s["estructuras"][name]
        for arr, lab, col in ((_arr(info["D_obs"]), "observado", "C0"), (_arr(info["D_ctrl"]), "control (mismo N)", "C1")):
            a = np.where(arr > 0, arr, np.nan)  # D = 0 no cabe en escala log
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)  # columnas sin D definido (no evaluable en ese q)
                ax.plot(xs, np.nanmedian(a, axis=0), marker="o", ms=3, color=col, label=lab)
                ax.fill_between(xs, np.nanpercentile(a, 25, axis=0), np.nanpercentile(a, 75, axis=0), color=col, alpha=0.18)
        ax.axhline(1.0, color="k", ls="--", lw=0.8)
        ax.set_yscale("log")
        ax.set_xticks(xs)
        ax.set_xticklabels(ticks)
        ax.set_xlim(xs[0] + 0.3, xs[-1] - 0.3)
        col, _, etiqueta = ESTILO[info["estado"]]
        fam = f", {info['familia']}" if info["familia"] else ""
        ev = "" if info["L"]["evaluable"] else ", NO EVALUABLE"
        ax.set_title(f"{name}\n{etiqueta}{fam}{ev}", fontsize=7.5, color=col if info["estado"] != "vista" else "k")
        ax.set_xlabel("q", fontsize=7)
        ax.set_ylabel("D / radio", fontsize=7)
        ax.grid(alpha=0.3, which="both")
        ax.tick_params(labelsize=6.5)
    axes.ravel()[0].legend(fontsize=6.5, loc="best")
    ax = axes.ravel()[len(names)]
    m = s["margen"]
    ys = np.arange(len(names))[::-1]
    ax.axvspan(-m, m, color="0.88")
    ax.axvline(-m, color="0.5", lw=0.6, ls=":")
    fin = [v for n in names for v in (_ext(s["estructuras"][n]["L"]["ic95"][0]), _ext(s["estructuras"][n]["L"]["ic95"][1]))
           if math.isfinite(v)]
    lim = max([3 * m, *[abs(v) for v in fin]]) * 1.1
    for y, name in zip(ys, names):
        info = s["estructuras"][name]
        col, mk, _ = ESTILO[info["estado"]]
        med, lo, hi = _ext(info["L"]["mediana"]), _ext(info["L"]["ic95"][0]), _ext(info["L"]["ic95"][1])
        if not info["L"]["evaluable"]:
            ax.plot([0], [y], marker="x", color=col, ms=6)
            continue
        cm, cl, ch = (float(np.clip(v, -lim, lim)) for v in (med, lo, hi))
        ax.errorbar([cm], [y], xerr=[[cm - cl], [ch - cm]], fmt=mk, color=col, capsize=2.5, ms=4,
                    mfc=col if info["estado"] in ("vota", "calibracion") else "white")
        for v, cv in ((lo, cl), (hi, ch)):
            if not math.isfinite(v):
                ax.annotate("∞" if v > 0 else "−∞", (cv, y), fontsize=6, ha="center", va="bottom")
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlim(-lim, lim)
    ax.set_yticks(ys)
    ax.set_yticklabels(names, fontsize=5.5)
    ax.set_xlabel("L = ln(D_obs/D_ctrl) en q = 1/8", fontsize=7)
    ax.set_title(f"L_1/8: mediana e IC 95% ({s['n_boot']} bootstrap)\nrojo vota, naranja no identificable, verde calibración,\n"
                 "gris vista, x no evaluable; franja ±ln 1.3", fontsize=6.5)
    ax.tick_params(labelsize=6.5)
    for extra in axes.ravel()[npan:]:
        extra.axis("off")
    man = s["manifest"]
    fig.suptitle(f"R20: dilución con control emparejado por N (N0 = {man['n0']}, {man['reals']} realizaciones; mediana y "
                 f"rango intercuartil; banda conforme). Veredicto: {s['veredicto']} ({s['n_votan']} votantes; {s['motivo']})",
                 fontsize=9)
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=110)
    plt.close(fig)
    print(fig_path)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--part", choices=("refs", "structs", "all"), default="all")
    r.add_argument("--struct", help="estructuras separadas por coma (por defecto, todas)")
    r.add_argument("--real", help="realizaciones a-b de las estructuras")
    r.add_argument("--s", help="sprinklings a-b de las referencias (en los tres conjuntos)")
    r.add_argument("--sellado", help=f"módulo sellado (por defecto {SELLADO})")
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--budget-min", type=float, help="no lanza unidades nuevas pasado este tiempo")
    r.add_argument("--n0", type=int)
    r.add_argument("--reals", type=int)
    r.add_argument("--nref", type=int, help="referencias por conjunto y por q (oficial 100; mínimo 19)")
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
