"""R14-B v2: respuesta de órdenes dados a la dilución aleatoria (docs/OMEGA_R14B_V2_PREREGISTRO.md). Exploratoria,
sin voto.

Ejecutor reanudable:
  run [--part refs|structs|all] [--struct a,b] [--geom a,b] [--q 0,1,2,3] [--real a-b] [--s a-b]
      [--workers 4] [--budget-min 110] [--n0 N --reals R --nref K --out DIR --allow-dirty]
  summarize [--out DIR]          -> <out>/summary.json, tabla_trayectorias.csv, tabla_bandas.csv
  figure [--out DIR] [--fig PNG] -> docs/figuras/r14b_v2.png

Unidades (un archivo por unidad; las que ya existen se saltan, así que repetir el comando reanuda):
- referencia: <out>/refs/<geometría>/q<qi>_s<ss>.json, un sprinkling de la geometría con N ~ Binomial(N0, q)
  elementos (N0 en q = 1). Así la referencia de q tiene exactamente la distribución de Φ_q aplicado a un sprinkling
  de N0 elementos (aclaración (e) del handoff). 60 por (geometría, q): s = 0..39 fijan la banda, 40..59 la validan.
- estructura: <out>/structs/<estructura>/r<rr>.json, una realización con N0 elementos medida en los cuatro q.

Semillas: SeedSequence(20261088, spawn_key=(323, tipo, ...)), con
  0 = puntos de referencia      (323, 0, g, qi, s)       g = índice en GEOMS (Binomial y luego los puntos)
  1 = selección de intervalos   (323, 1, g, qi, s)       de la referencia
  2 = puntos de la estructura   (323, 2, k, N0, real)    k = índice en STRUCTS
  3 = azar auxiliar             (323, 3, k, N0, real)    links de C5, agregadas de D4, gemelos de D5
  4 = dilución                  (323, 4, k, N0, real, qi) un generador independiente por q (q = 1 no lo usa)
  5 = selección de intervalos   (323, 5, k, N0, real, qi)
  6 = bootstrap del resumen     (323, 6, j)               j = 100·k + g (estructuras), 9000 + familia
Las estructuras selladas en F6 (C0, C4, C5, D0, D2, D4 a D8) son aquí vistas: se regeneran con esta entropía.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # un hilo por proceso; el paralelismo viene de --workers

import argparse
import csv
import hashlib
import json
import math
import multiprocessing as mp
import sys
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

_TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(_TOOLS.parent))
sys.path.insert(0, str(_TOOLS))

from es0_piloto import _clean, git_info
from f6_confirmacion import SELLADAS as SELLADAS_I1
from f6_exploracion import STRUCT, between_blocks, interval_obs
from f6_intento2 import SELLADAS as SELLADAS_I2
from f6_intento2 import twin_pairs

from omega.causality.estimators import mm_dimension

EXPERIMENT_ID = 323
ENTROPY = 20261088
K_REF_GEN, K_REF_SEL, K_EST_GEN, K_EST_AUX, K_DIL, K_EST_SEL, K_BOOT = 0, 1, 2, 3, 4, 5, 6
N0 = 8000
QS = (1.0, 0.5, 0.25, 0.125)
N_REAL = 10
N_REFS = 60
N_FIT = 40  # s < N_FIT fijan la banda; el resto la valida
M_R = (10, 20, 40)
LOCAL_FRAC = 0.1
MIN_INTERVALS, MAX_INTERVALS = 30, 300
K_PROF = 8
N_BOOT = 1000
TWIN_SD_FLOOR = 0.5  # sd mínima de √(pares gemelos) (Poisson estabilizado)
COMPONENTS = ("frac_rel", *(f"r_{m}" for m in M_R), *(f"perfil_{k}" for k in range(K_PROF + 1)), "raiz_gemelos")
SPEC = Path("docs") / "OMEGA_R14B_V2_PREREGISTRO.md"
OUT = Path("results") / "omega_r14b"
FIG = Path("docs") / "figuras" / "r14b_v2.png"
Bool = npt.NDArray[np.bool_]
Float = npt.NDArray[np.float64]
Gen = Callable[[int, np.random.Generator, np.random.Generator], tuple[Bool, dict[str, Any]]]


def rng_for(*key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(ENTROPY, spawn_key=(EXPERIMENT_ID, *key))))


# --- Estructuras (§5) --------------------------------------------------------------------------------


def _from_struct(name: str) -> Gen:
    """Generador de la exploración (`f6_exploracion.STRUCT[name]`): (n, rng) -> (C, x); no usa aux."""
    g = STRUCT[name][2]

    def gen(n: int, r: np.random.Generator, aux: np.random.Generator) -> tuple[Bool, dict[str, Any]]:
        return g(n, r)[0], {}

    return gen


def _from_sealed(table: list[tuple[str, str, Any]], name: str) -> Gen:
    """Generador de una tabla sellada (C0 a C7 o D0 a D8): (n, gen, aux) -> (C, x, extra)."""
    g = {nm: f for nm, _, f in table}[name]

    def gen(n: int, r: np.random.Generator, aux: np.random.Generator) -> tuple[Bool, dict[str, Any]]:
        c, _, extra = g(n, r, aux)
        return c, extra

    return gen


# (nombre, familia, generador, origen). El orden es el índice k del spawn_key; las 8 primeras son GEOMS (índice g).
STRUCTS: list[tuple[str, str, Gen, str]] = [
    ("plana", "geometria", _from_struct("diamante_4d"), "f6_exploracion.STRUCT['diamante_4d']"),
    ("de_sitter", "geometria", _from_struct("de_sitter_4d"), "f6_exploracion.STRUCT['de_sitter_4d']"),
    ("de_sitter_fuerte", "geometria", _from_struct("de_sitter_4d_fuerte"), "f6_exploracion.STRUCT['de_sitter_4d_fuerte']"),
    ("anti_de_sitter", "geometria", _from_struct("anti_de_sitter_4d"), "f6_exploracion.STRUCT['anti_de_sitter_4d']"),
    ("frw_polvo", "geometria", _from_struct("frw_polvo_4d"), "f6_exploracion.STRUCT['frw_polvo_4d']"),
    ("D2_frw_exponencial", "geometria", _from_sealed(SELLADAS_I2, "D2_frw_exponencial"), "f6_intento2.d2_frw_exp"),
    ("C0_kasner", "geometria", _from_sealed(SELLADAS_I1, "C0_kasner_vacio"), "f6_confirmacion.c0_kasner"),
    ("D0_kasner", "geometria", _from_sealed(SELLADAS_I2, "D0_kasner_u2"), "f6_intento2.d0_kasner"),
    ("C5_links_borrados", "defecto", _from_sealed(SELLADAS_I1, "C5_links_borrados_13"), "f6_confirmacion.c5_links_borrados"),
    ("D4_relaciones_agregadas", "defecto", _from_sealed(SELLADAS_I2, "D4_relaciones_agregadas"), "f6_intento2.d4_agregadas"),
    ("D5_gemelos", "defecto", _from_sealed(SELLADAS_I2, "D5_gemelos_10"), "f6_intento2.d5_gemelos"),
    ("D6_cono_tapa_plana", "cono", _from_sealed(SELLADAS_I2, "D6_cono_tapa_plana"), "f6_intento2._gauge_gen(d6_gauge)"),
    ("D7_cono_huso", "cono", _from_sealed(SELLADAS_I2, "D7_cono_huso"), "f6_intento2._gauge_gen(d7_gauge, BOX_D7)"),
    ("D8_cono_suma_normas", "cono", _from_sealed(SELLADAS_I2, "D8_cono_suma_normas"), "f6_intento2._gauge_gen(d8_gauge)"),
    ("C4_cono_esfera_cubo", "cono", _from_sealed(SELLADAS_I1, "C4_cono_esfera_cubo"), "f6_confirmacion._cone_gen(_within_esfera_cubo)"),
    ("cono_octaedrico", "cono", _from_struct("cono_octaedrico_4d"), "f6_exploracion.STRUCT['cono_octaedrico_4d']"),
    ("producto_4", "producto", _from_struct("producto_4"), "f6_exploracion.STRUCT['producto_4']"),
]
NAMES = [s[0] for s in STRUCTS]
GEOMS = [s[0] for s in STRUCTS if s[1] == "geometria"]
FAMILY = {s[0]: s[1] for s in STRUCTS}
GEN = {s[0]: s[2] for s in STRUCTS}
FAMILIES = ("geometria", "defecto", "cono", "producto")
assert NAMES[: len(GEOMS)] == GEOMS


def primary_target(name: str) -> str:
    """Objetivo principal: la propia geometría; para defectos, conos y producto, la plana (aclaración (e))."""
    return name if FAMILY[name] == "geometria" else "plana"


# --- Mapa de grano grueso y observables (§2, §4) -----------------------------------------------------


def dilute(c: Bool, q: float, rng: np.random.Generator | None) -> tuple[npt.NDArray[np.int64], Bool]:
    """Φ_q: conserva cada elemento con probabilidad q (un uniforme por elemento, en orden de etiqueta) y devuelve
    (índices conservados, orden inducido). q = 1 no consume azar."""
    n = c.shape[0]
    if q >= 1.0:
        idx = np.arange(n, dtype=np.int64)
        return idx, c
    if rng is None:
        raise ValueError("q < 1 requiere generador")
    idx = np.flatnonzero(rng.random(n) < q).astype(np.int64)
    return idx, c[np.ix_(idx, idx)]


def measure(c: Bool, rng: np.random.Generator) -> dict[str, Any]:
    """Observables intensivos de §4 sobre un orden: fracción de relaciones y d_MM; r por sub-intervalo (m ∈ M_R,
    m ≤ 0.1·m_max, entre 30 y 300 intervalos con exactamente m interiores, como `f6_exploracion.measure_order`);
    perfil N_k/N_pares (k = 0..8; N_pares = N(N−1)/2); pares gemelos (igual pasado y futuro)."""
    t0 = time.time()
    n = int(c.shape[0])
    pairs = n * (n - 1) / 2
    n_rel = int(np.count_nonzero(c))
    b = between_blocks(c)
    m_max = int(np.max(b)) if n_rel else -1
    sel = b[(b >= 0) & (b <= K_PROF)]
    prof = np.bincount(sel.astype(np.int64), minlength=K_PROF + 1)[: K_PROF + 1]
    del sel
    por_m: dict[str, Any] = {}
    for m in M_R:
        local = m <= LOCAL_FRAC * m_max
        cell: dict[str, Any] = {"local": bool(local)}
        if local:
            ii, jj = np.nonzero(b == m)
            cnt = int(ii.size)
            cell.update({"disponibles": cnt, "medible": cnt >= MIN_INTERVALS})
            if cnt >= MIN_INTERVALS:
                pick = rng.choice(cnt, size=min(MAX_INTERVALS, cnt), replace=False)
                rs = np.array([interval_obs(c, b, int(ii[p]), int(jj[p]))[0] for p in pick])
                cell.update({"n_sel": int(pick.size), "r_media": float(rs.mean()), "r_sd": float(rs.std(ddof=1))})
        else:
            cell["medible"] = False
        por_m[str(m)] = cell
    del b
    tw = twin_pairs(c)
    return {"N": n, "relaciones": n_rel, "frac_rel": n_rel / pairs, "d_mm": mm_dimension(c), "m_max": m_max,
            "N_k": [int(v) for v in prof], "perfil": [float(v) / pairs for v in prof], "por_m": por_m,
            "pares_gemelos": tw, "frac_gemelos": tw / pairs, "t_medida_s": round(time.time() - t0, 2)}


def obs_vector(row: dict[str, Any]) -> Float:
    """Vector O de la distancia de Mahalanobis (aclaración (a)): frac_rel, r_m (nan si no medible), perfil N_k/N_pares
    y √(pares gemelos). d_MM queda fuera (función monótona de frac_rel) y se informa aparte."""
    r = [row["por_m"][str(m)]["r_media"] if row["por_m"][str(m)]["medible"] else math.nan for m in M_R]
    return np.array([row["frac_rel"], *r, *row["perfil"], math.sqrt(row["pares_gemelos"])], dtype=np.float64)


# --- Bandas de Mahalanobis (§3) -----------------------------------------------------------------------


def ledoit_wolf(z: Float) -> tuple[Float, float]:
    """Covarianza con contracción de Ledoit-Wolf (2004) hacia μ·I; z (n × p) ya estandarizada. (Σ, λ)."""
    n, p = z.shape
    x = z - z.mean(axis=0)
    s = x.T @ x / n
    mu = float(np.trace(s)) / p
    d2 = float(np.sum((s - mu * np.eye(p)) ** 2))
    b2 = float(sum(np.sum((np.outer(v, v) - s) ** 2) for v in x)) / (n * n)
    lam = 1.0 if d2 <= 0 else min(b2, d2) / d2
    out: Float = lam * mu * np.eye(p) + (1.0 - lam) * s
    return out, lam


class Band:
    """Banda de una referencia en un subconjunto de componentes: media, escala y covarianza contraída de las filas de
    ajuste; radio = max(p95 de las distancias de ajuste, menor radio que contiene el 95% de las de validación)."""

    def __init__(self, fit: Float, val: Float, comps: tuple[int, ...], floors: Float) -> None:
        self.comps = comps
        f, v = fit[:, comps], val[:, comps]
        self.mu = f.mean(axis=0)
        raw = f.std(axis=0, ddof=1)
        fl = floors[list(comps)]
        self.sd = np.maximum(raw, fl)
        cov, self.lam = ledoit_wolf((f - self.mu) / self.sd)
        # Componente con piso (en la práctica, √gemelos casi siempre 0 en la referencia): su varianza estandarizada
        # pasa a ser al menos 1, es decir, al menos piso² en unidades originales.
        floored = np.flatnonzero(raw < fl)
        cov[floored, floored] = np.maximum(cov[floored, floored], 1.0)
        self.floored = [COMPONENTS[comps[i]] for i in floored]
        self.prec = np.linalg.inv(cov)
        d_fit = self.dist(f, sub=True)
        d_val = self.dist(v, sub=True)
        self.r_fit = float(np.percentile(d_fit, 95))
        self.cobertura = float(np.mean(d_val <= self.r_fit)) if d_val.size else math.nan
        k = math.ceil(0.95 * d_val.size)
        self.r_val = float(np.sort(d_val)[k - 1]) if d_val.size else math.nan
        self.radio = max(self.r_fit, self.r_val) if d_val.size else self.r_fit
        self.cobertura_final = float(np.mean(d_val <= self.radio)) if d_val.size else math.nan

    def dist(self, x: Float, sub: bool = False) -> Float:
        y = np.atleast_2d(x if sub else x[..., list(self.comps)])
        z = (y - self.mu) / self.sd
        out: Float = np.sqrt(np.einsum("ij,jk,ik->i", z, self.prec, z))
        return out

    def info(self) -> dict[str, Any]:
        return {"componentes": [COMPONENTS[i] for i in self.comps], "dim": len(self.comps), "lambda_lw": self.lam,
                "con_piso": self.floored,
                "radio_ajuste_p95": self.r_fit, "cobertura_validacion": self.cobertura, "radio_validacion_95": self.r_val,
                "radio_banda": self.radio, "margen_delta": self.radio - self.r_fit, "cobertura_final": self.cobertura_final}


def comp_floors() -> Float:
    f = np.full(len(COMPONENTS), 1e-12)
    f[COMPONENTS.index("raiz_gemelos")] = TWIN_SD_FLOOR
    return f


class Reference:
    """Referencia de una geometría en un q: filas de ajuste y de validación, bandas por subconjunto (en caché)."""

    def __init__(self, vecs: Float, n_fit: int) -> None:
        self.fit, self.val = vecs[:n_fit], vecs[n_fit:]
        self.base = tuple(int(i) for i in np.flatnonzero(np.all(np.isfinite(vecs), axis=0)))
        self.floors = comp_floors()
        self._cache: dict[tuple[int, ...], Band] = {}

    def band(self, comps: tuple[int, ...]) -> Band:
        if comps not in self._cache:
            self._cache[comps] = Band(self.fit, self.val, comps, self.floors)
        return self._cache[comps]

    def normalized(self, vec: Float) -> tuple[float, tuple[int, ...]]:
        """D/radio de un vector, con las componentes de la referencia que el vector tiene medibles (aclaración (c))."""
        comps = tuple(i for i in self.base if math.isfinite(vec[i]))
        b = self.band(comps)
        return float(b.dist(vec)[0]) / b.radio, comps


# --- Etiquetas (§2, §6) -------------------------------------------------------------------------------


def state(lo: float, hi: float) -> str:
    """Estado de una escala según el IC del estadístico D/radio: dentro, fuera o ambiguo."""
    if hi <= 1.0:
        return "dentro"
    if lo > 1.0:
        return "fuera"
    return "ambiguo"


def label(states: list[str], grows: bool) -> str:
    """Etiqueta descriptiva (q en el orden de QS: 1, 1/2, 1/4, 1/8).

    - en banda microscópica: dentro en q = 1.
    - converge a banda: dentro en q = 1/4 y en q = 1/8, y la distancia no crece de q = 1/2 a q = 1/8.
    - no converge: fuera en q = 1/4 o en q = 1/8, o dentro en ambos pero la distancia crece de q = 1/2 a q = 1/8.
    - indeterminada: el resto (algún estado ambiguo en q = 1/4 o 1/8 sin ninguno fuera).

    «Crece» (aclaración (f)): el IC 95% emparejado por bootstrap de med(1/8) − med(1/2), con med la mediana de
    D/radio entre realizaciones, queda entero por encima de 0.
    """
    if states[0] == "dentro":
        return "en banda microscópica"
    s4, s8 = states[2], states[3]
    if s4 == "dentro" and s8 == "dentro":
        return "no converge" if grows else "converge a banda"
    if "fuera" in (s4, s8):
        return "no converge"
    return "indeterminada"


def trajectory(dn: Float, rng: np.random.Generator, n_boot: int = N_BOOT) -> dict[str, Any]:
    """dn: (realizaciones × q) de D/radio. Mediana por q e IC 95% por bootstrap de realizaciones completas (emparejado
    entre q). nan (componentes sin banda) se ignoran en la mediana."""
    n = dn.shape[0]
    med = np.nanmedian(dn, axis=0)
    boots = np.empty((n_boot, dn.shape[1]))
    for b in range(n_boot):
        boots[b] = np.nanmedian(dn[rng.integers(0, n, size=n)], axis=0)
    lo, hi = np.percentile(boots, 2.5, axis=0), np.percentile(boots, 97.5, axis=0)
    states = [state(float(a), float(c)) for a, c in zip(lo, hi)]
    diff = boots[:, 3] - boots[:, 1]
    d_ic = [float(np.percentile(diff, 2.5)), float(np.percentile(diff, 97.5))]
    grows = d_ic[0] > 0
    return {"mediana": med.tolist(), "ic95": [[float(a), float(c)] for a, c in zip(lo, hi)],
            "frac_en_banda": np.mean(dn <= 1.0, axis=0).tolist(), "estados": states,
            "cambio_1_8_menos_1_2": {"valor": float(med[3] - med[1]), "ic95": d_ic, "crece": grows},
            "etiqueta": label(states, grows), "en_banda_todo_q": all(s == "dentro" for s in states)}


# --- Unidades de trabajo ------------------------------------------------------------------------------


def ref_path(out: Path, geom: str, qi: int, s: int) -> Path:
    return out / "refs" / geom / f"q{qi}_s{s:02d}.json"


def struct_path(out: Path, name: str, real: int) -> Path:
    return out / "structs" / name / f"r{real:02d}.json"


def _write(p: Path, rec: dict[str, Any]) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(_clean(rec), default=str))
    tmp.replace(p)


def run_ref(task: dict[str, Any]) -> str:
    t0 = time.time()
    geom, qi, s, n0 = task["geom"], task["qi"], task["s"], task["n0"]
    g = GEOMS.index(geom)
    q = QS[qi]
    r = rng_for(K_REF_GEN, g, qi, s)
    n = n0 if q >= 1.0 else int(r.binomial(n0, q))
    c, extra = GEN[geom](n, r, r)  # las geometrías no usan aux
    t_c = time.time() - t0
    row = {"geometria": geom, "qi": qi, "q": q, "s": s, "t_C_s": round(t_c, 2), **extra,
           **measure(c, rng_for(K_REF_SEL, g, qi, s))}
    p = ref_path(Path(task["out"]), geom, qi, s)
    _write(p, {"spec": task["spec"], "row": row})
    return f"ref {geom} q{qi} s{s} N={n} {time.time() - t0:.1f} s"


def run_struct(task: dict[str, Any]) -> str:
    t0 = time.time()
    name, real, n0 = task["name"], task["real"], task["n0"]
    k = NAMES.index(name)
    c, extra = GEN[name](n0, rng_for(K_EST_GEN, k, n0, real), rng_for(K_EST_AUX, k, n0, real))
    t_c = time.time() - t0
    rows: list[dict[str, Any]] = []
    for qi, q in enumerate(QS):
        idx, cd = dilute(c, q, rng_for(K_DIL, k, n0, real, qi) if q < 1.0 else None)
        rows.append({"estructura": name, "realizacion": real, "qi": qi, "q": q, **measure(cd, rng_for(K_EST_SEL, k, n0, real, qi))})
        del cd, idx
    _write(struct_path(Path(task["out"]), name, real), {"spec": task["spec"], "t_C_s": round(t_c, 2), "extra": extra, "rows": rows})
    return f"est {name} r{real} {time.time() - t0:.1f} s (C {t_c:.1f} s)"


def _work(task: dict[str, Any]) -> str:
    return run_ref(task) if task["kind"] == "ref" else run_struct(task)


# Costo relativo para ordenar (los más largos primero): Kasner ≈ 4× el resto.
def _cost(task: dict[str, Any]) -> float:
    name = task.get("geom") or task.get("name") or ""
    w = 4.0 if "kasner" in name else (2.0 if name in ("C5_links_borrados", "D4_relaciones_agregadas") else 1.0)
    if task["kind"] == "ref":
        return w * QS[task["qi"]] ** 2
    return w * 1.5


def _pool(tasks: list[dict[str, Any]], workers: int, budget_s: float | None) -> int:
    """Corre las tareas con hasta `workers` procesos; con presupuesto deja de lanzar tareas nuevas al agotarlo."""
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


def _range(spec: str | None, total: int) -> list[int]:
    if spec is None:
        return list(range(total))
    a, _, b = spec.partition("-")
    lo, hi = int(a), int(b) if b else int(a)
    if not 0 <= lo <= hi < total:
        sys.exit(f"rango fuera de 0..{total - 1}: {spec}")
    return list(range(lo, hi + 1))


def _names(spec: str | None, valid: list[str]) -> list[str]:
    if spec is None:
        return list(valid)
    out = spec.split(",")
    bad = [x for x in out if x not in valid]
    if bad:
        sys.exit(f"desconocidos: {bad}; válidos: {valid}")
    return out


def spec_sha() -> str:
    return hashlib.sha256(SPEC.read_bytes()).hexdigest() if SPEC.exists() else ""


def cmd_run(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    override = any(v is not None for v in (args.n0, args.reals, args.nref))
    if override and not args.out:
        sys.exit("--n0/--reals/--nref exigen --out (no se mezcla con la corrida oficial)")
    head, dirty = git_info()
    if dirty and not args.allow_dirty:
        sys.exit("exige arbol git limpio fuera de results/ (usar --allow-dirty solo para pruebas)")
    if not SPEC.exists():
        sys.exit(f"falta {SPEC}")
    n0 = args.n0 or N0
    reals = args.reals or N_REAL
    nref = args.nref or N_REFS
    official = not dirty and not override and out == OUT
    spec = {"commit": head, "dirty": dirty, "official": official, "n0": n0, "reals": reals, "nref": nref,
            "preregistro_sha256": spec_sha(), "experiment_id": EXPERIMENT_ID, "entropia": ENTROPY}
    manifest = out / "manifest.json"
    if manifest.exists():
        old = json.loads(manifest.read_text())
        if (old["n0"], old["reals"], old["nref"], old["official"]) != (n0, reals, nref, official):
            sys.exit(f"{manifest} tiene otra configuración: {old}")
    else:
        _write(manifest, {k: spec[k] for k in ("n0", "reals", "nref", "official", "experiment_id", "entropia")})
    qis = [int(v) for v in args.q.split(",")] if args.q else list(range(len(QS)))
    tasks: list[dict[str, Any]] = []
    if args.part in ("refs", "all"):
        for geom in _names(args.geom, GEOMS):
            for qi in qis:
                for s in _range(args.s, nref):
                    if not ref_path(out, geom, qi, s).exists():
                        tasks.append({"kind": "ref", "geom": geom, "qi": qi, "s": s, "n0": n0, "out": str(out), "spec": spec})
    if args.part in ("structs", "all"):
        for name in _names(args.struct, NAMES):
            for real in _range(args.real, reals):
                if not struct_path(out, name, real).exists():
                    tasks.append({"kind": "est", "name": name, "real": real, "n0": n0, "out": str(out), "spec": spec})
    tasks.sort(key=lambda t: -_cost(t))
    print(f"{out} (commit {head[:8]}, dirty {dirty}, oficial {official}): {len(tasks)} unidades por correr")
    t0 = time.time()
    done = _pool(tasks, args.workers, args.budget_min * 60 if args.budget_min else None)
    print(f"{done} de {len(tasks)} unidades; pared {time.time() - t0:.1f} s"
          + ("" if done == len(tasks) else "; repetir el comando para seguir"))


# --- Resumen ------------------------------------------------------------------------------------------


def _load_json(p: Path) -> dict[str, Any]:
    data: dict[str, Any] = json.loads(p.read_text())
    return data


def load_all(out: Path) -> tuple[dict[str, Any], dict[tuple[str, int], list[dict[str, Any]]], dict[str, list[list[dict[str, Any]]]], set[tuple[str, bool, bool, str]]]:
    man = _load_json(out / "manifest.json")
    prov: set[tuple[str, bool, bool, str]] = set()
    refs: dict[tuple[str, int], list[dict[str, Any]]] = {}
    for geom in GEOMS:
        for qi in range(len(QS)):
            rows = []
            for s in range(man["nref"]):
                p = ref_path(out, geom, qi, s)
                if not p.exists():
                    raise SystemExit(f"falta {p}")
                d = _load_json(p)
                prov.add((d["spec"]["commit"], d["spec"]["dirty"], d["spec"]["official"], d["spec"]["preregistro_sha256"]))
                rows.append(d["row"])
            refs[(geom, qi)] = rows
    structs: dict[str, list[list[dict[str, Any]]]] = {}
    for name in NAMES:
        reals = []
        for real in range(man["reals"]):
            p = struct_path(out, name, real)
            if not p.exists():
                raise SystemExit(f"falta {p}")
            d = _load_json(p)
            prov.add((d["spec"]["commit"], d["spec"]["dirty"], d["spec"]["official"], d["spec"]["preregistro_sha256"]))
            reals.append(sorted(d["rows"], key=lambda r: int(r["qi"])))
        structs[name] = reals
    return man, refs, structs, prov


def n_fit_for(nref: int) -> int:
    """40 de 60 en la corrida oficial; en pruebas con otro nref, la misma proporción (2/3)."""
    return N_FIT if nref == N_REFS else max(2, math.ceil(2 * nref / 3))


def _stat(v: list[float]) -> dict[str, Any]:
    a = np.array([x for x in v if x is not None and math.isfinite(x)], dtype=np.float64)
    if a.size == 0:
        return {"n": 0, "media": None, "sd": None}
    return {"n": int(a.size), "media": float(a.mean()), "sd": float(a.std(ddof=1)) if a.size > 1 else None}


def summarize(out: Path) -> dict[str, Any]:
    man, refs, structs, prov = load_all(out)
    nf = n_fit_for(man["nref"])
    ref_obj = {key: Reference(np.array([obs_vector(r) for r in rows]), nf) for key, rows in refs.items()}
    bandas: dict[str, Any] = {}
    t_band: list[list[Any]] = []
    for (geom, qi), ref in ref_obj.items():
        b = ref.band(ref.base)
        rows = refs[(geom, qi)]
        info = {**b.info(), "q": QS[qi], "N_media": float(np.mean([r["N"] for r in rows])),
                "omitidas": [COMPONENTS[i] for i in range(len(COMPONENTS)) if i not in ref.base],
                "observables": {"frac_rel": _stat([r["frac_rel"] for r in rows]), "d_mm": _stat([r["d_mm"] for r in rows]),
                                **{f"r_{m}": _stat([r["por_m"][str(m)].get("r_media", math.nan) for r in rows]) for m in M_R},
                                "frac_gemelos": _stat([r["frac_gemelos"] for r in rows])}}
        bandas.setdefault(geom, {})[str(qi)] = info
        t_band.append([geom, QS[qi], info["N_media"], info["dim"], ";".join(info["omitidas"]), info["lambda_lw"],
                       info["radio_ajuste_p95"], info["cobertura_validacion"], info["radio_validacion_95"], info["radio_banda"],
                       info["cobertura_final"]])
    res: dict[str, Any] = {}
    t_traj: list[list[Any]] = []
    for k, name in enumerate(NAMES):
        reals = structs[name]
        est: dict[str, Any] = {"familia": FAMILY[name], "objetivo_principal": primary_target(name), "por_objetivo": {},
                                "observables": {}}
        for qi in range(len(QS)):
            rows = [rr[qi] for rr in reals]
            est["observables"][str(qi)] = {
                "q": QS[qi], "N": _stat([r["N"] for r in rows]), "frac_rel": _stat([r["frac_rel"] for r in rows]),
                "d_mm": _stat([r["d_mm"] for r in rows]), "frac_gemelos": _stat([r["frac_gemelos"] for r in rows]),
                **{f"r_{m}": _stat([r["por_m"][str(m)].get("r_media", math.nan) for r in rows]) for m in M_R},
                "r_medible_en": {str(m): sum(bool(r["por_m"][str(m)]["medible"]) for r in rows) for m in M_R}}
        for g, geom in enumerate(GEOMS):
            dn = np.full((len(reals), len(QS)), np.nan)
            omit: dict[str, int] = {}
            for i, rr in enumerate(reals):
                for qi in range(len(QS)):
                    ref = ref_obj[(geom, qi)]
                    val, comps = ref.normalized(obs_vector(rr[qi]))
                    dn[i, qi] = val
                    for j in ref.base:
                        if j not in comps:
                            key = f"{qi}:{COMPONENTS[j]}"
                            omit[key] = omit.get(key, 0) + 1
            tr = trajectory(dn, rng_for(K_BOOT, 100 * k + g))
            tr["D_sobre_radio"] = dn.tolist()
            tr["componentes_omitidas"] = omit
            est["por_objetivo"][geom] = tr
            if geom == primary_target(name) or FAMILY[name] != "geometria":
                for qi in range(len(QS)):
                    t_traj.append([name, FAMILY[name], geom, geom == primary_target(name), QS[qi],
                                   est["observables"][str(qi)]["N"]["media"], tr["mediana"][qi], tr["ic95"][qi][0],
                                   tr["ic95"][qi][1], tr["frac_en_banda"][qi], tr["estados"][qi], tr["etiqueta"]])
        prim = est["por_objetivo"][primary_target(name)]
        est["etiqueta"] = prim["etiqueta"]
        est["objetivos_en_banda_o_convergen"] = [g for g, t in est["por_objetivo"].items()
                                                  if t["etiqueta"] in ("en banda microscópica", "converge a banda")]
        res[name] = est
    fam = family_summary(res)
    pred = {
        "geometrias_en_banda_todo_q": {g: res[g]["por_objetivo"][g]["en_banda_todo_q"] for g in GEOMS},
        "D5_converge_a_banda": res["D5_gemelos"]["etiqueta"] == "converge a banda",
        "conos_no_convergen": {n: res[n]["etiqueta"] == "no converge" for n in NAMES if FAMILY[n] == "cono"},
    }
    _write_csv(out / "tabla_bandas.csv", ["geometria", "q", "N_media", "dim", "omitidas", "lambda_lw", "radio_ajuste_p95",
                                          "cobertura_validacion", "radio_validacion_95", "radio_banda", "cobertura_final"], t_band)
    _write_csv(out / "tabla_trayectorias.csv", ["estructura", "familia", "objetivo", "principal", "q", "N_media", "mediana_D_radio",
                                                "ic95_lo", "ic95_hi", "frac_en_banda", "estado", "etiqueta"], t_traj)
    commits = sorted({p[0] for p in prov})
    return {"estado": "EXPLORATORIA (R14-B v2), sin voto", "preregistro": str(SPEC), "preregistro_sha256_actual": spec_sha(),
            "experiment_id": EXPERIMENT_ID, "entropia": ENTROPY, "manifest": man, "n_ajuste": nf,
            "procedencia": {"commits": commits, "un_solo_commit": len(commits) == 1, "alguna_sucia": any(p[1] for p in prov),
                            "todas_oficiales": all(p[2] for p in prov), "sha_preregistro": sorted({p[3] for p in prov})},
            "componentes": list(COMPONENTS), "q": list(QS), "bandas": bandas, "estructuras": res, "familias": fam,
            "predicciones_congeladas": pred}


def family_summary(res: dict[str, Any]) -> dict[str, Any]:
    """Por familia y q: mediana entre estructuras de la mediana D/radio (objetivo principal), con IC por bootstrap
    jerárquico (estructuras y luego realizaciones), y conteo de etiquetas."""
    out: dict[str, Any] = {}
    for fi, fam in enumerate(FAMILIES):
        names = [n for n in NAMES if FAMILY[n] == fam]
        mats = [np.array(res[n]["por_objetivo"][res[n]["objetivo_principal"]]["D_sobre_radio"]) for n in names]
        r = rng_for(K_BOOT, 9000 + fi)
        point = np.median([np.nanmedian(m, axis=0) for m in mats], axis=0)
        boots = np.empty((N_BOOT, len(QS)))
        for b in range(N_BOOT):
            pick = r.integers(0, len(mats), size=len(mats))
            meds = [np.nanmedian(mats[j][r.integers(0, mats[j].shape[0], size=mats[j].shape[0])], axis=0) for j in pick]
            boots[b] = np.median(meds, axis=0)
        labels: dict[str, int] = {}
        for n in names:
            labels[res[n]["etiqueta"]] = labels.get(res[n]["etiqueta"], 0) + 1
        out[fam] = {"estructuras": names, "mediana_D_radio": point.tolist(),
                    "ic95": [[float(a), float(c)] for a, c in zip(np.percentile(boots, 2.5, axis=0), np.percentile(boots, 97.5, axis=0))],
                    "etiquetas": labels}
    return out


def _write_csv(p: Path, header: list[str], rows: list[list[Any]]) -> None:
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def cmd_summarize(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    s = summarize(out)
    (out / "summary.json").write_text(json.dumps(_clean(s), indent=2, default=str))
    pv = s["procedencia"]
    print(f"commits {pv['commits']} oficial {pv['todas_oficiales']} sucia {pv['alguna_sucia']}")
    print("\nBandas (dim, cobertura de validación con el p95 de ajuste, radio final):")
    for geom, by_q in s["bandas"].items():
        print(f"  {geom:22s} " + "  ".join(f"q={QS[int(qi)]:<5} d{b['dim']} cob {b['cobertura_validacion']:.2f} R {b['radio_banda']:.2f}"
                                          for qi, b in by_q.items()))
    print("\nTrayectorias (mediana D/radio frente al objetivo principal; q = 1, 1/2, 1/4, 1/8):")
    print(f"  {'estructura':24s} {'objetivo':18s} {'q=1':>7s} {'1/2':>7s} {'1/4':>7s} {'1/8':>7s}  etiqueta")
    for name, info in s["estructuras"].items():
        tr = info["por_objetivo"][info["objetivo_principal"]]
        print(f"  {name:24s} {info['objetivo_principal']:18s} " + " ".join(f"{v:7.2f}" for v in tr["mediana"]) + f"  {info['etiqueta']}")
    print("\nFamilias:", json.dumps({f: v["etiquetas"] for f, v in s["familias"].items()}, ensure_ascii=False))
    print("Predicciones congeladas:", json.dumps(s["predicciones_congeladas"], ensure_ascii=False))


# --- Figura -------------------------------------------------------------------------------------------

TITLES = {"geometria": "Geometrías (frente a su propia referencia)", "defecto": "Defectos (frente a la plana)",
          "cono": "Conos (frente a la plana)", "producto": "Producto 4 (frente a la plana)"}


def cmd_figure(args: argparse.Namespace) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(args.out) if args.out else OUT
    fig_path = Path(args.fig) if args.fig else FIG
    s = _load_json(out / "summary.json")
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    xs = np.log2(np.array(QS))
    n0 = s["manifest"]["n0"]
    for ax, fam in zip(axes.ravel(), FAMILIES):
        names = [n for n in NAMES if FAMILY[n] == fam]
        for j, name in enumerate(names):
            info = s["estructuras"][name]
            tr = info["por_objetivo"][info["objetivo_principal"]]
            med = np.array(tr["mediana"], dtype=float)
            ci = np.array(tr["ic95"], dtype=float)
            off = (j - (len(names) - 1) / 2) * 0.04
            ax.errorbar(xs + off, med, yerr=[med - ci[:, 0], ci[:, 1] - med], marker="o", ms=4, capsize=3,
                        label=f"{name}: {info['etiqueta']}")
        ax.axhline(1.0, color="k", ls="--", lw=1)
        ax.text(0.99, 1.02, "borde de la banda", ha="right", va="bottom", fontsize=8, transform=ax.get_yaxis_transform())
        ax.set_yscale("log")
        lo, hi = ax.get_ylim()
        ax.set_ylim(min(lo, 0.5), max(hi, 2.0))
        ax.set_title(TITLES[fam], fontsize=10)
        ax.set_ylabel("D de Mahalanobis / radio de la banda")
        ax.set_xticks(xs)
        ax.set_xticklabels([f"q = 1/{round(1 / q)}\nN ef. {int(n0 * q)}" if q < 1 else f"q = 1\nN ef. {n0}" for q in QS])
        ax.set_xlim(xs[0] + 0.3, xs[-1] - 0.3)  # q = 1 a la izquierda: el grano engrosa hacia la derecha
        ax.legend(fontsize=7, loc="best", framealpha=0.85)
        ax.grid(alpha=0.3, which="both")
    nr = s["manifest"]["reals"]
    fig.suptitle(f"R14-B v2: respuesta a la dilución aleatoria (mediana de {nr} realizaciones, IC 95% por bootstrap)", fontsize=11)
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=130)
    print(fig_path)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--part", choices=("refs", "structs", "all"), default="all")
    r.add_argument("--struct", help="estructuras separadas por coma (por defecto, las 17)")
    r.add_argument("--geom", help="geometrías de referencia separadas por coma (por defecto, las 8)")
    r.add_argument("--q", help="índices de q separados por coma (0 = 1, 1 = 1/2, 2 = 1/4, 3 = 1/8)")
    r.add_argument("--real", help="realizaciones a-b de las estructuras")
    r.add_argument("--s", help="sprinklings a-b de las referencias")
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--budget-min", type=float, help="no lanza unidades nuevas pasado este tiempo (tandas < 2 h)")
    r.add_argument("--n0", type=int)
    r.add_argument("--reals", type=int)
    r.add_argument("--nref", type=int)
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
