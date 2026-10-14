"""F6 exploratoria (docs/OMEGA_F6_EXPLORACION.md, R11-3). Partes: ref | nula_<N> | est_<nombre>_<N>[_r<s>] | final.

EXPLORATORIA, sin voto: nada de lo que produzca cuenta como evidencia.

- ``ref``: referencia exacta por m (4000 diamantes planos 4D de m puntos más las puntas): r_plano(m) y el
  perfil plano N_k(m), k = 0..⌊m/2⌋.
- ``nula_<N>`` (N ∈ {4000, 8000}): 30 diamantes planos 4D completos medidos con el instrumento (§1).
- ``est_<nombre>_<N>``: las realizaciones de una estructura en un N (10 en 4000 y 8000; 5 en 16000 para
  las marcadas ★). En N = 8000 cada realización se mide además diluida al azar a q = 0.5 (4000 elementos).
  ``est_<nombre>_<N>_r<s>`` es la misma parte restringida a la realización s (mismas semillas, mismas filas),
  por si la parte completa no cabe en un comando; ``final`` usa ``est_<nombre>_<N>.json`` si existe y, si no,
  une las partes por realización.
- ``final``: Δr, ajustes b + A·x (y b + A·x + B·x² aparte) con bootstrap de realizaciones, colapso de la
  dilución, D_loc(m) en sus dos versiones frente a la nula, tablas y ``summary.json``.

Semillas: SeedSequence(ENTROPIA, spawn_key=(320, familia, índice, N, realización)); familia 0 = referencias
(índice 0, N ← m, realización ← diamante de referencia), 1 = nula (índice 0), 2 = estructuras (índice = orden
en ``STRUCT``), 3 = dilución (índice = orden en ``STRUCT``, N = 8000). El bootstrap de ``final`` usa la
familia 4 (índice = orden en ``STRUCT``; no fijada por el protocolo).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Callable, Literal

import numpy as np
import numpy.typing as npt

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from f1f2 import random_dorder  # noqa: E402
from f5prima import sprinkle_polyhedral  # noqa: E402
from f7 import _check_p, de_sitter_eta_points, frw_power_points, sprinkle_lp  # noqa: E402

from omega.causality.orders import transitive_percolation  # noqa: E402
from omega.causality.sprinkling import DIAMOND_HALF_HEIGHT, causal_matrix, sprinkle_diamond  # noqa: E402
from omega.causality.torus import sprinkle_slab  # noqa: E402

EXPERIMENT_ID = 320
ENTROPY = 20261075
FAM_REF, FAM_NULA, FAM_EST, FAM_DIL, FAM_BOOT = 0, 1, 2, 3, 4
M_GRID = (10, 20, 40, 80, 160)
LOCAL_FRAC = 0.1
MIN_INTERVALS, MAX_INTERVALS = 30, 300
N_REF = 4000
NS = (4000, 8000)
N_BIG = 16000
N_REAL, N_REAL_BIG, N_NULA = 10, 5, 30
N_DIL, Q_DIL = 8000, 0.5
N_BOOT = 1000
SLAB_LEN = 0.45
ADS_Z0, ADS_ZMIN = 0.7, 0.2
BLOCK_ROWS = 1024
DOC = Path("docs") / "OMEGA_CONSEJO_REV11_F6.md"
SPEC = Path("docs") / "OMEGA_F6_EXPLORACION.md"
OUT = Path("results") / "omega_f6_exploracion"
Bool = npt.NDArray[np.bool_]
Float = npt.NDArray[np.float64]
Int = npt.NDArray[np.int32]
Kind = Literal["mink", "slab", "max", "l1", "lp"]
StructGen = Callable[[int, np.random.Generator], tuple[Bool, Float | None]]


def rng_for(*key: int) -> np.random.Generator:
    """Generador de la exploración: spawn_key = (320, familia, índice, N, realización)."""
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(ENTROPY, spawn_key=(EXPERIMENT_ID, *key))))


# --- Matriz causal por bloques y conteos «entre» ----------------------------------------------------


def _relation_rows(coords: Float, lo: int, hi: int, kind: Kind, p: float) -> Bool:
    """Filas lo..hi de la relación; mismas operaciones elemento a elemento que las funciones originales
    (`causal_matrix`, `slab_causal_matrix`, `polyhedral_matrix`, `lp_matrix`), así que el resultado es idéntico."""
    t = coords[:, 0]
    dt = t[None, :] - t[lo:hi, None]
    if kind in ("mink", "slab"):
        dx2 = np.zeros_like(dt)
        for k in range(1, coords.shape[1]):
            diff = coords[None, :, k] - coords[lo:hi, None, k]
            if kind == "slab":
                diff = (diff + 0.5) % 1.0 - 0.5
            dx2 += diff * diff
        rel: Bool = (dt > 0) & (dt * dt >= dx2)
        return rel
    acc = np.zeros_like(dt)
    for k in range(1, coords.shape[1]):
        diff = np.abs(coords[None, :, k] - coords[lo:hi, None, k])
        if kind == "max":
            np.maximum(acc, diff, out=acc)
        elif kind == "l1":
            acc += diff
        else:
            acc += diff**p
    if kind == "lp":
        acc = acc ** (1.0 / p)
    rel2: Bool = (dt > 0) & (dt >= acc)
    return rel2


def causal_blocks(coords: Float, kind: Kind = "mink", p: float = 2.0, block: int = BLOCK_ROWS) -> Bool:
    """Matriz causal construida por bloques de filas (sin matrices float64 n × n completas)."""
    if kind == "slab" and float(np.max(coords[:, 0]) - np.min(coords[:, 0])) >= 0.5:
        raise ValueError("la duracion debe ser < 1/2 para que la imagen minima sea la unica")
    if kind == "lp":
        _check_p(p)
    n = coords.shape[0]
    c = np.empty((n, n), dtype=np.bool_)
    for lo in range(0, n, block):
        hi = min(n, lo + block)
        c[lo:hi] = _relation_rows(coords, lo, hi, kind, p)
    return c


def between_blocks(c: Bool, block: int = 2048) -> Int:
    """B[i, j] = elementos estrictamente entre i y j si i precede a j; −1 si no (como `between_counts`).

    float32 es exacto hasta 2^24; el producto se hace por bloques de filas y se guarda en int32.
    """
    n = c.shape[0]
    if n >= 2**24:
        raise ValueError("n demasiado grande para float32 exacto")
    cf = c.astype(np.float32)
    b = np.empty((n, n), dtype=np.int32)
    for lo in range(0, n, block):
        hi = min(n, lo + block)
        blk = np.rint(cf[lo:hi] @ cf).astype(np.int32)
        blk[~c[lo:hi]] = -1
        b[lo:hi] = blk
    return b


# --- Generadores -----------------------------------------------------------------------------------


def ads_points(n: int, r: np.random.Generator, z0: float = ADS_Z0, zmin: float = ADS_ZMIN) -> Float:
    """Anti de Sitter de Poincaré: diamante de `sprinkle_diamond` con z = x₃ + z0 ∈ [0.2, 1.2] y aceptación
    (zmin/z)⁴ (densidad ∝ z⁻⁴, volumen propio); la relación es la de Minkowski (invariancia conforme)."""
    if not math.isclose(z0 - DIAMOND_HALF_HEIGHT, zmin) or zmin <= 0:
        raise ValueError("se requiere z0 − T = zmin > 0")
    pts: list[Float] = []
    total = 0
    while total < n:
        cand = sprinkle_diamond(4 * n, 4, r)
        z = cand[:, 3] + z0
        keep = cand[r.random(cand.shape[0]) < (zmin / z) ** 4]
        pts.append(keep)
        total += keep.shape[0]
    out: Float = np.concatenate(pts)[:n]
    return out


def _coord_gen(points: Callable[[int, np.random.Generator], Float], kind: Kind, p: float = 2.0) -> StructGen:
    def gen(n: int, r: np.random.Generator) -> tuple[Bool, Float | None]:
        x = points(n, r)
        return causal_blocks(x, kind, p), x

    return gen


def _order_gen(f: Callable[[int, np.random.Generator], Bool]) -> StructGen:
    def gen(n: int, r: np.random.Generator) -> tuple[Bool, Float | None]:
        return f(n, r), None

    return gen


# El orden del diccionario es el índice de estructura del spawn_key (familias 2, 3 y 4). Tabla de §2.
# (papel, ★, generador)
STRUCT: dict[str, tuple[str, bool, StructGen]] = {
    "diamante_4d": ("geometria", True, _coord_gen(lambda n, r: sprinkle_diamond(n, 4, r), "mink")),
    "de_sitter_4d": ("geometria", False, _coord_gen(lambda n, r: de_sitter_eta_points(n, 4, r), "mink")),
    "de_sitter_4d_fuerte": ("geometria", True, _coord_gen(lambda n, r: de_sitter_eta_points(n, 4, r, -1.2, -0.2), "mink")),
    "anti_de_sitter_4d": ("geometria", True, _coord_gen(ads_points, "mink")),
    "frw_radiacion_4d": ("geometria", False, _coord_gen(lambda n, r: frw_power_points(n, 4, 4.0, r), "mink")),
    "frw_polvo_4d": ("geometria", False, _coord_gen(lambda n, r: frw_power_points(n, 4, 8.0, r), "mink")),
    "bloque": ("geometria", True, _coord_gen(lambda n, r: sprinkle_slab(n, r, SLAB_LEN), "slab")),
    "caja": ("geometria", False, _coord_gen(lambda n, r: sprinkle_slab(n, r, SLAB_LEN), "mink")),
    "cono_cubico_4d": ("impostor", False, _coord_gen(lambda n, r: sprinkle_polyhedral(n, "max", r), "max")),
    "cono_octaedrico_4d": ("impostor", True, _coord_gen(lambda n, r: sprinkle_polyhedral(n, "l1", r), "l1")),
    "cono_p6_4d": ("impostor", False, _coord_gen(lambda n, r: sprinkle_lp(n, 6.0, r), "lp", 6.0)),
    "cono_p1_5_4d": ("zona_gris", False, _coord_gen(lambda n, r: sprinkle_lp(n, 1.5, r), "lp", 1.5)),
    "cono_p3_4d": ("zona_gris", False, _coord_gen(lambda n, r: sprinkle_lp(n, 3.0, r), "lp", 3.0)),
    "producto_4": ("impostor", False, _order_gen(random_dorder(4))),
    "percolacion_005": ("impostor", False, _order_gen(lambda n, r: transitive_percolation(n, 0.05, r))),
}


def struct_ns(name: str) -> tuple[int, ...]:
    return (*NS, N_BIG) if STRUCT[name][1] else NS


def n_real(n: int) -> int:
    return N_REAL_BIG if n == N_BIG else N_REAL


# --- Instrumento -----------------------------------------------------------------------------------


def local_ms(m_max: int) -> list[int]:
    """Regla de localidad S1: m de la rejilla con m ≤ 0.1·m_max."""
    return [m for m in M_GRID if m <= LOCAL_FRAC * m_max]


def interval_suborder(c: Bool, i: int, j: int) -> npt.NDArray[np.int64]:
    """Índices del sub-orden del intervalo [i, j]: la punta i, los interiores y la punta j."""
    inner = np.flatnonzero(c[i] & c[:, j])
    return np.concatenate([[i], inner, [j]]).astype(np.int64)


def interval_obs(c: Bool, b: npt.NDArray[np.integer[Any]], i: int, j: int) -> tuple[float, npt.NDArray[np.int64]]:
    """r sobre los m interiores (sin puntas) y perfil N_k, k = 0..⌊m/2⌋, del sub-orden con puntas.

    Un intervalo es convexo: los elementos entre dos de sus elementos están en él, así que los conteos del
    sub-orden coinciden con los de B global.
    """
    idx = interval_suborder(c, i, j)
    m = idx.size - 2
    inner = idx[1:-1]
    r = float(np.count_nonzero(c[np.ix_(inner, inner)])) / (m * (m - 1) / 2) if m > 1 else float("nan")
    sub = b[np.ix_(idx, idx)]
    kmax = m // 2
    prof = np.bincount(sub[sub >= 0].astype(np.int64), minlength=kmax + 1)[: kmax + 1].astype(np.int64)
    return r, prof


def box_truncated(x: Float, i: int, j: int) -> bool:
    """¿El diamante geométrico de (i, j) se sale de la caja espacial [0, 1]³? (solo para el informe)."""
    a, bb = x[i, 1:], x[j, 1:]
    half_t = (x[j, 0] - x[i, 0]) / 2
    centre = (a + bb) / 2
    f = float(np.linalg.norm(bb - a)) / 2
    u = (bb - a) / (2 * f) if f > 0 else np.zeros(3)
    semi = math.sqrt(max(half_t * half_t - f * f, 0.0))
    ext = np.sqrt(half_t * half_t * u * u + semi * semi * (1 - u * u))
    return bool(np.any(centre - ext < 0) or np.any(centre + ext > 1))


def measure_order(c: Bool, r: np.random.Generator, coords: Float | None = None, box: bool = False) -> dict[str, Any]:
    """Instrumento común (§1) sobre una realización: B una sola vez, m_max, y por m local hasta 300 intervalos."""
    t0 = time.time()
    b = between_blocks(c)
    m_max = int(np.max(b))
    t_b = time.time() - t0
    loc = set(local_ms(m_max))
    por_m: dict[str, Any] = {}
    for m in M_GRID:
        if m not in loc:
            por_m[str(m)] = {"local": False, "disponibles": int(np.count_nonzero(b == m)), "medible": False}
            continue
        ii, jj = np.nonzero(b == m)
        cnt = int(ii.size)
        row: dict[str, Any] = {"local": True, "disponibles": cnt, "medible": cnt >= MIN_INTERVALS}
        if cnt >= MIN_INTERVALS:
            sel = r.choice(cnt, size=min(MAX_INTERVALS, cnt), replace=False)
            rs: list[float] = []
            profs: list[npt.NDArray[np.int64]] = []
            trunc: list[bool] = []
            for q in sel:
                rr, pr = interval_obs(c, b, int(ii[q]), int(jj[q]))
                rs.append(rr)
                profs.append(pr)
                if box and coords is not None:
                    trunc.append(box_truncated(coords, int(ii[q]), int(jj[q])))
            ra = np.array(rs)
            row.update({"n_sel": int(sel.size), "r_media": float(ra.mean()), "r_sd": float(ra.std(ddof=1)),
                        "perfil_medio": np.mean(np.array(profs, dtype=np.float64), axis=0).tolist()})
            if trunc:
                row["frac_recortados"] = float(np.mean(trunc))
        por_m[str(m)] = row
    return {"N": int(c.shape[0]), "m_max": m_max, "por_m": por_m, "t_B_s": round(t_b, 2), "t_medida_s": round(time.time() - t0, 2)}


# --- Partes ----------------------------------------------------------------------------------------


def ref_rows(reps: int = N_REF, ms: tuple[int, ...] = M_GRID) -> list[dict[str, Any]]:
    """Referencia exacta por m: m puntos uniformes en el diamante plano 4D más las puntas (±T, 0, 0, 0)."""
    tips = np.array([[-DIAMOND_HALF_HEIGHT, 0, 0, 0], [DIAMOND_HALF_HEIGHT, 0, 0, 0]], dtype=np.float64)
    rows: list[dict[str, Any]] = []
    for m in ms:
        rs, profs = [], []
        for s in range(reps):
            x = np.vstack([tips[:1], sprinkle_diamond(m, 4, rng_for(FAM_REF, 0, m, s)), tips[1:]])
            c = causal_matrix(x)
            b = between_blocks(c)
            if np.count_nonzero(c[0] & c[:, m + 1]) != m:
                raise RuntimeError("las puntas no contienen los m puntos")
            rr, pr = interval_obs(c, b, 0, m + 1)
            rs.append(rr)
            profs.append(pr)
        ra, pa = np.array(rs), np.array(profs, dtype=np.float64)
        rows.append({"m": m, "reps": reps, "r_plano": float(ra.mean()), "r_sd": float(ra.std(ddof=1)),
                     "perfil_plano": pa.mean(axis=0).tolist(), "perfil_sd": pa.std(axis=0, ddof=1).tolist()})
    return rows


def nula_rows(n: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for s in range(N_NULA):
        r = rng_for(FAM_NULA, 0, n, s)
        c = causal_blocks(sprinkle_diamond(n, 4, r), "mink")
        rows.append({"realizacion": s, "diluida": False, **measure_order(c, r)})
        print(f"nula {n} {s}", rows[-1]["t_medida_s"], "s", file=sys.stderr, flush=True)
    return rows


def est_rows(name: str, n: int, reals: list[int] | None = None) -> list[dict[str, Any]]:
    role, _, gen = STRUCT[name]
    k = list(STRUCT).index(name)
    box = name == "caja"
    rows: list[dict[str, Any]] = []
    for s in reals if reals is not None else range(n_real(n)):
        t0 = time.time()
        r = rng_for(FAM_EST, k, n, s)
        c, x = gen(n, r)
        t_c = time.time() - t0
        rows.append({"estructura": name, "papel": role, "realizacion": s, "diluida": False, "q": 1.0,
                     "t_C_s": round(t_c, 2), **measure_order(c, r, x, box)})
        print(f"{name} {n} {s}", round(time.time() - t0, 1), "s", file=sys.stderr, flush=True)
        if n == N_DIL:
            t0 = time.time()
            rd = rng_for(FAM_DIL, k, n, s)
            keep = np.sort(rd.choice(n, size=int(round(Q_DIL * n)), replace=False))
            cd = c[np.ix_(keep, keep)]
            xd = x[keep] if x is not None else None
            del c
            rows.append({"estructura": name, "papel": role, "realizacion": s, "diluida": True, "q": Q_DIL,
                         **measure_order(cd, rd, xd, box)})
            print(f"{name} {n} {s} diluida", round(time.time() - t0, 1), "s", file=sys.stderr, flush=True)
    return rows


# --- Final -----------------------------------------------------------------------------------------


def log_ratio(nk: Float, ref: Float) -> Float:
    out: Float = np.log((nk + 1) / (ref + 1))
    return out


def d_loc(nk: Float, ref: Float) -> tuple[float, float]:
    """(forma, máximo): max_k |v_k − media_k(v)| (F7) y max_k |v_k| (F5′)."""
    v = log_ratio(nk, ref)
    return float(np.max(np.abs(v - v.mean()))), float(np.max(np.abs(v)))


def fit_poly(x: Float, y: Float, deg: int) -> list[float] | None:
    """Mínimos cuadrados y = Σ c_i x^i (i = 0..deg); None si no hay puntos suficientes."""
    if np.unique(x).size < deg + 1:
        return None
    a = np.vander(x, deg + 1, increasing=True)
    coef = np.linalg.lstsq(a, y, rcond=None)[0]
    return [float(v) for v in coef]


def _cells(dr: dict[int, Float], keep: dict[int, npt.NDArray[np.bool_]], idx: dict[int, npt.NDArray[np.int64]]) -> tuple[Float, Float]:
    """Puntos (x, Δr medio) de las celdas incluidas, con las realizaciones `idx` por N."""
    xs, ys = [], []
    for n, arr in dr.items():
        sub = arr[idx[n]]
        for jm, m in enumerate(M_GRID):
            if not keep[n][jm]:
                continue
            col = sub[:, jm]
            col = col[np.isfinite(col)]
            if col.size:
                xs.append(math.sqrt(m / n))
                ys.append(float(col.mean()))
    return np.array(xs), np.array(ys)


def fit_structure(dr: dict[int, Float], k: int) -> dict[str, Any]:
    """Δr = b + A·x (y aparte b + A·x + B·x²), con bootstrap de realizaciones completas dentro de cada N."""
    keep = {n: np.sum(np.isfinite(a), axis=0) >= a.shape[0] / 2 for n, a in dr.items()}
    full = {n: np.arange(a.shape[0]) for n, a in dr.items()}
    x, y = _cells(dr, keep, full)
    out: dict[str, Any] = {"puntos": [[float(a), float(b)] for a, b in zip(x, y)]}
    r = rng_for(FAM_BOOT, k)
    boots: dict[int, list[list[float]]] = {1: [], 2: []}
    for _ in range(N_BOOT):
        idx = {n: np.asarray(r.integers(0, a.shape[0], size=a.shape[0]), dtype=np.int64) for n, a in dr.items()}
        xb, yb = _cells(dr, keep, idx)
        for deg in (1, 2):
            f = fit_poly(xb, yb, deg)
            if f is not None:
                boots[deg].append(f)
    for deg, names, label in ((1, ("b", "A"), "lineal"), (2, ("b", "A", "B"), "cuadratico")):
        est = fit_poly(x, y, deg)
        if est is None:
            out[label] = None
            continue
        bt = np.array(boots[deg]) if boots[deg] else np.full((1, deg + 1), np.nan)
        out[label] = {nm: {"valor": est[i], "se_boot": float(np.std(bt[:, i], ddof=1)) if bt.shape[0] > 1 else None,
                           "ic95": [float(np.percentile(bt[:, i], 2.5)), float(np.percentile(bt[:, i], 97.5))]}
                      for i, nm in enumerate(names)}
        out[label]["replicas_validas"] = len(boots[deg])
    return out


def load(part: str) -> dict[str, Any]:
    data: dict[str, Any] = json.loads((OUT / f"{part}.json").read_text())
    return data


def load_est(name: str, n: int) -> list[dict[str, Any]]:
    if (OUT / f"est_{name}_{n}.json").exists():
        rows: list[dict[str, Any]] = load(f"est_{name}_{n}")["rows"]
        return rows
    rows = [row for s in range(n_real(n)) for row in load(f"est_{name}_{n}_r{s}")["rows"]]
    return rows


def _stat(v: list[float]) -> dict[str, Any]:
    a = np.array([u for u in v if math.isfinite(u)])
    if a.size == 0:
        return {"n": 0, "media": None, "se": None}
    return {"n": int(a.size), "media": float(a.mean()), "se": float(a.std(ddof=1) / math.sqrt(a.size)) if a.size > 1 else None}


def _dr_matrix(rows: list[dict[str, Any]], rplano: dict[int, float]) -> Float:
    """Δr por realización (filas) y m (columnas); nan si la (realización, m) no es medible."""
    out = np.full((len(rows), len(M_GRID)), np.nan)
    for i, row in enumerate(rows):
        for jm, m in enumerate(M_GRID):
            cell = row["por_m"][str(m)]
            if cell["medible"]:
                out[i, jm] = cell["r_media"] - rplano[m]
    return out


def _dloc_rows(rows: list[dict[str, Any]], pref: dict[int, Float]) -> dict[int, list[tuple[float, float]]]:
    out: dict[int, list[tuple[float, float]]] = {m: [] for m in M_GRID}
    for row in rows:
        for m in M_GRID:
            cell = row["por_m"][str(m)]
            if cell["medible"]:
                out[m].append(d_loc(np.array(cell["perfil_medio"]), pref[m]))
    return out


def final(head: str) -> dict[str, Any]:
    ref = {row["m"]: row for row in load("ref")["rows"]}
    rplano = {m: float(ref[m]["r_plano"]) for m in M_GRID}
    pref = {m: np.array(ref[m]["perfil_plano"], dtype=np.float64) for m in M_GRID}
    # Nula por (N, m): distribución de D_loc (forma y máximo) y Δr.
    nula: dict[str, Any] = {}
    p95: dict[tuple[int, int], tuple[float, float]] = {}
    nula_maxt: dict[int, list[float]] = {}
    for n in NS:
        rows = load(f"nula_{n}")["rows"]
        dl = _dloc_rows(rows, pref)
        drm = _dr_matrix(rows, rplano)
        nula[str(n)] = {}
        for jm, m in enumerate(M_GRID):
            if len(dl[m]) < 2:
                nula[str(n)][str(m)] = {"n": len(dl[m])}
                continue
            f, mx = np.array([v[0] for v in dl[m]]), np.array([v[1] for v in dl[m]])
            p95[(n, m)] = (float(np.percentile(f, 95)), float(np.percentile(mx, 95)))
            nula[str(n)][str(m)] = {"n": len(dl[m]), "dloc_forma_media": float(f.mean()), "dloc_forma_p95": p95[(n, m)][0],
                                    "dloc_max_media": float(mx.mean()), "dloc_max_p95": p95[(n, m)][1],
                                    "delta_r": _stat(drm[:, jm].tolist())}
        nula_maxt[n] = [_max_t(row, pref, p95, n) for row in rows]
    n_half = int(round(Q_DIL * N_DIL))
    res: dict[str, Any] = {}
    t_dr: list[list[Any]] = []
    t_dl: list[list[Any]] = []
    t_fit: list[list[Any]] = []
    t_col: list[list[Any]] = []
    for k, (name, (role, star, _)) in enumerate(STRUCT.items()):
        by_n = {n: load_est(name, n) for n in struct_ns(name)}
        plain = {n: [r for r in rows if not r["diluida"]] for n, rows in by_n.items()}
        dil = [r for r in by_n[N_DIL] if r["diluida"]]
        dr = {n: _dr_matrix(rows, rplano) for n, rows in plain.items()}
        info: dict[str, Any] = {"papel": role, "estrella": star, "por_N": {}}
        groups: list[tuple[str, int, list[dict[str, Any]], Float]] = [(str(n), n, plain[n], dr[n]) for n in plain]
        groups.append((f"{N_DIL}_diluida_q{Q_DIL}", n_half, dil, _dr_matrix(dil, rplano)))
        for label, n_eff, rows, drm in groups:
            mm = [r["m_max"] for r in rows]
            dl = _dloc_rows(rows, pref)
            per_m: dict[str, Any] = {}
            for jm, m in enumerate(M_GRID):
                cells = [r["por_m"][str(m)] for r in rows]
                st = _stat(drm[:, jm].tolist())
                e: dict[str, Any] = {"x": math.sqrt(m / n_eff), "local_en": sum(c["local"] for c in cells),
                                     "medible_en": sum(c["medible"] for c in cells), "realizaciones": len(rows),
                                     "disponibles_mediana": float(np.median([c["disponibles"] for c in cells])),
                                     "delta_r": st}
                tr = [c["frac_recortados"] for c in cells if "frac_recortados" in c]
                if tr:
                    e["frac_recortados_media"] = float(np.mean(tr))
                if dl[m]:
                    f, mx = np.array([v[0] for v in dl[m]]), np.array([v[1] for v in dl[m]])
                    nu = p95.get((n_eff, m)) if label.isdigit() else None
                    e["dloc_forma_media"], e["dloc_max_media"] = float(f.mean()), float(mx.mean())
                    e["nula_p95"] = {"forma": nu[0], "maximo": nu[1]} if nu else None
                    if nu:
                        e["frac_sobre_p95"] = {"forma": float(np.mean(f > nu[0])), "maximo": float(np.mean(mx > nu[1]))}
                    t_dl.append([name, label, m, len(dl[m]), e["dloc_forma_media"], nu[0] if nu else "", e["dloc_max_media"], nu[1] if nu else ""])
                per_m[str(m)] = e
                t_dr.append([name, label, m, e["x"], e["medible_en"], len(rows), st["media"], st["se"], e.get("frac_recortados_media", "")])
            maxt = [_max_t(r, pref, p95, n_eff) for r in rows] if label.isdigit() and n_eff in NS else []
            info["por_N"][label] = {"m_max": {"media": float(np.mean(mm)), "min": int(min(mm)), "max": int(max(mm))},
                                    "por_m": per_m,
                                    "max_T": ({"mediana": float(np.nanmedian(maxt)) if np.any(np.isfinite(maxt)) else None,
                                               "nula_p95": float(np.nanpercentile(nula_maxt[n_eff], 95))
                                               if np.any(np.isfinite(nula_maxt[n_eff])) else None} if maxt else None)}
        fit = fit_structure(dr, k)
        info["ajuste"] = fit
        for label in ("lineal", "cuadratico"):
            fv = fit[label]
            if fv is None:
                t_fit.append([name, label, len(fit["puntos"]), "", "", "", "", "", ""])
                continue
            vals = [fv[nm] for nm in ("b", "A", "B") if nm in fv]
            row_f: list[Any] = [name, label, len(fit["puntos"])]
            for v in vals:
                row_f += [v["valor"], f"[{v['ic95'][0]:.5f}, {v['ic95'][1]:.5f}]"]
            t_fit.append(row_f + [""] * (9 - len(row_f)))
        # Colapso: Δr diluida (q = 0.5, N = 8000) frente a Δr en N = 4000 y el mismo m.
        drd = _dr_matrix(dil, rplano)
        col: dict[str, Any] = {}
        for jm, m in enumerate(M_GRID):
            a, b4, b8 = _stat(drd[:, jm].tolist()), _stat(dr[n_half][:, jm].tolist()), _stat(dr[N_DIL][:, jm].tolist())
            d = None
            if a["media"] is not None and b4["media"] is not None:
                se = math.sqrt((a["se"] or 0.0) ** 2 + (b4["se"] or 0.0) ** 2)
                d = {"diluida_menos_4000": a["media"] - b4["media"], "se": se}
            col[str(m)] = {"diluida": a, "N4000": b4, "N8000": b8, "diferencia": d}
            t_col.append([name, m, a["media"], b4["media"], b8["media"], d["diluida_menos_4000"] if d else "", d["se"] if d else ""])
        info["colapso"] = col
        res[name] = info
    _write_csv("tabla_delta_r.csv", ["estructura", "N", "m", "x", "medible_en", "realizaciones", "delta_r", "se", "frac_recortados"], t_dr)
    _write_csv("tabla_dloc.csv", ["estructura", "N", "m", "n", "dloc_forma", "nula_p95_forma", "dloc_max", "nula_p95_max"], t_dl)
    _write_csv("tabla_ajustes.csv", ["estructura", "modelo", "puntos", "b", "b_ic95", "A", "A_ic95", "B", "B_ic95"], t_fit)
    _write_csv("tabla_colapso.csv", ["estructura", "m", "delta_r_diluida", "delta_r_4000", "delta_r_8000", "diluida_menos_4000", "se"], t_col)
    return {"commit": head, "estado": "EXPLORATORIA (R11-3), sin voto", "experiment_id": EXPERIMENT_ID,
            "referencia": {str(m): {"r_plano": rplano[m], "reps": ref[m]["reps"]} for m in M_GRID},
            "nula": nula, "nula_max_T_p95": {str(n): float(np.nanpercentile(v, 95)) if np.any(np.isfinite(v)) else None
                                             for n, v in nula_maxt.items()},
            "estructuras": res}


def _max_t(row: dict[str, Any], pref: dict[int, Float], p95: dict[tuple[int, int], tuple[float, float]], n: int) -> float:
    """Descriptor estilo max-T: max_m D_loc,forma(m)/p95_nula(N, m) sobre los m medibles de la realización."""
    vals = []
    for m in M_GRID:
        cell = row["por_m"][str(m)]
        if cell["medible"] and (n, m) in p95 and p95[(n, m)][0] > 0:
            vals.append(d_loc(np.array(cell["perfil_medio"]), pref[m])[0] / p95[(n, m)][0])
    return max(vals) if vals else float("nan")


def _write_csv(name: str, header: list[str], rows: list[list[Any]]) -> None:
    with (OUT / name).open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


# --- CLI -------------------------------------------------------------------------------------------


def parse_est(part: str) -> tuple[str, int, list[int] | None] | None:
    mt = re.fullmatch(r"est_(?P<name>.+)_(?P<n>\d+)(?:_r(?P<s>\d+))?", part)
    if mt is None or mt["name"] not in STRUCT:
        return None
    n = int(mt["n"])
    if n not in struct_ns(mt["name"]):
        return None
    if mt["s"] is None:
        return mt["name"], n, None
    s = int(mt["s"])
    return (mt["name"], n, [s]) if s < n_real(n) else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", required=True)
    args = ap.parse_args()
    part: str = args.part
    if "R11-3 | **Firmada" not in DOC.read_text() or not SPEC.exists():
        sys.exit("F6 exploratoria requiere la firma R11-3 y el protocolo congelado")
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip():
        sys.exit("exige arbol git limpio")
    t0 = time.time()
    rows: list[dict[str, Any]]
    est = parse_est(part) if part.startswith("est_") else None
    if part == "ref":
        rows = ref_rows()
    elif part.startswith("nula_") and part[5:].isdigit() and int(part[5:]) in NS:
        rows = nula_rows(int(part[5:]))
    elif est is not None:
        rows = est_rows(*est)
    elif part == "final":
        out = final(head)
        (OUT / "summary.json").write_text(json.dumps(out, indent=2, default=str))
        print("final ok", round(time.time() - t0, 1), "s")
        return
    else:
        sys.exit("parte desconocida")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{part}.json").write_text(json.dumps({"commit": head, "tiempo_s": round(time.time() - t0, 1), "rows": rows}, default=str))
    print(part, "ok", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
