"""R14-C v2: observables causales multivariables (docs/OMEGA_R14C_V2_PREREGISTRO.md). Exploratoria, sin voto.

Ejecutor reanudable:
  run [--part mc|refs|structs|all] [--struct a,b] [--geom a,b] [--n 4000,8000] [--real a-b] [--s a-b]
      [--workers 4] [--blas-threads 1] [--budget-min 110] [--ns N1,N2 --reals R --nref K --nmc M --out DIR --allow-dirty]
  summarize [--out DIR]          -> <out>/summary.json, <out>/tabla_r14c.csv (la tabla única de §3)
  figure [--out DIR] [--fig PNG] -> docs/figuras/r14c_v2.png
  cost [--n 8000] [--struct plana]  -> mide una realización completa (tiempos de C, B, sub-intervalos, Laplaciano y
        eigh); hilos BLAS con OPENBLAS_NUM_THREADS=k delante del comando

Unidades (un archivo por unidad; las que ya existen se saltan, así que repetir el comando reanuda):
- mc:   <out>/mc/m<m>.json, 200 diamantes planos 4D de m puntos (más las puntas para el punto medio): referencia
        de los observables por sub-intervalo (§1, «calibración de referencias»).
- refs: <out>/refs/<geometría>/N<N>_s<ss>.json, un sprinkling de la geometría (60 por N y por geometría). Cada
        geometría es su propia referencia (R14-B v2 §3); los impostores usan la plana.
- structs: <out>/structs/<estructura>/N<N>_r<rr>.json, una realización (10 por estructura y N).

Semillas: SeedSequence(20261089, spawn_key=(324, tipo, ...)), con
  0 = puntos de la referencia global      (324, 0, g, N, s)      g = índice en GEOMS
  1 = azar auxiliar de la referencia      (324, 1, g, N, s)      solo C2 (modos)
  2 = selección de intervalos (referencia) (324, 2, g, N, s)
  3 = puntos de la estructura             (324, 3, k, N, real)   k = índice en STRUCTS
  4 = azar auxiliar de la estructura      (324, 4, k, N, real)   modos de C2, links de C5, agregadas de D4, gemelos de D5
  5 = selección de intervalos (estructura) (324, 5, k, N, real)
  6 = diamantes de Monte Carlo            (324, 6, m, s)
  7 = bootstrap del resumen               (324, 7, j)            j = 10000·iN + 100·familia + observable
Todas las estructuras son vistas (v1 §3): se regeneran con esta entropía; sus semillas no coinciden con F6 ni R14-B.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # un hilo por proceso salvo --blas-threads (se fija antes de lanzar los procesos)

import argparse
import csv
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
import scipy.linalg as sla
from scipy.stats import energy_distance, norm

_TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(_TOOLS.parent))
sys.path.insert(0, str(_TOOLS))

from es0_piloto import _clean, git_info
from f6_confirmacion import SELLADAS as SELLADAS_I1
from f6_exploracion import between_blocks, interval_suborder
from f6_intento2 import SELLADAS as SELLADAS_I2
from r14b import Gen, _from_sealed, _from_struct

from omega.causality.estimators import midpoint_dimension, mm_dimension
from omega.causality.sprinkling import (
    DIAMOND_HALF_HEIGHT,
    causal_matrix,
    sprinkle_diamond,
)

EXPERIMENT_ID = 324
ENTROPY = 20261089
K_REF_GEN, K_REF_AUX, K_REF_SEL, K_EST_GEN, K_EST_AUX, K_EST_SEL, K_MC, K_BOOT = 0, 1, 2, 3, 4, 5, 6, 7
NS = (4000, 8000)
N_REAL = 10
N_REFS = 60
N_MC = 200
MS = (20, 40, 80, 160)
LOCAL_FRAC = 0.1
MIN_INTERVALS, MAX_INTERVALS = 30, 300
KS = (2, 3, 4, 5)
# χ_k en 4D (VERIFICACION_R15_2 §B1): χ_k = 12^(k−1)·6 / [k (2k−1)! (2k+1)!]
CHI = {2: 1 / 20, 3: 1 / 2100, 4: 1 / 705600, 5: 1 / 582120000}
R_PLANO = 0.1
LAP_BINS = 1000  # histograma del espectro en [0, 2]
N_BOOT = 2000
MIN_REF_VALID = 20  # celdas de referencia con menos de 20 de 60 valores: sin comparación
MIN_REAL_VALID = 5  # estructura con menos de 5 de 10 realizaciones válidas en una celda: fuera de esa celda
ALPHA = 0.05
SPEC = Path("docs") / "OMEGA_R14C_V2_PREREGISTRO.md"
OUT = Path("results") / "omega_r14c"
FIG = Path("docs") / "figuras" / "r14c_v2.png"
Bool = npt.NDArray[np.bool_]
Float = npt.NDArray[np.float64]

# --- Estructuras (v1 §3) -------------------------------------------------------------------------------
# (nombre, familia, generador, origen). Orden = índice k del spawn_key. Las de geometrías van primero (índice g).
STRUCTS: list[tuple[str, str, Gen, str]] = [
    ("plana", "conforme", _from_struct("diamante_4d"), "f6_exploracion.STRUCT['diamante_4d']"),
    ("de_sitter", "conforme", _from_struct("de_sitter_4d"), "f6_exploracion.STRUCT['de_sitter_4d']"),
    ("de_sitter_fuerte", "conforme", _from_struct("de_sitter_4d_fuerte"), "f6_exploracion.STRUCT['de_sitter_4d_fuerte']"),
    ("anti_de_sitter", "conforme", _from_struct("anti_de_sitter_4d"), "f6_exploracion.STRUCT['anti_de_sitter_4d']"),
    ("frw_polvo", "conforme", _from_struct("frw_polvo_4d"), "f6_exploracion.STRUCT['frw_polvo_4d']"),
    ("D2_frw_exponencial", "conforme", _from_sealed(SELLADAS_I2, "D2_frw_exponencial"), "f6_intento2.d2_frw_exp"),
    ("D3_frw_rigida", "conforme", _from_sealed(SELLADAS_I2, "D3_frw_rigida"), "f6_intento2.d3_frw_rigida"),
    ("C2_densidad_suave", "conforme", _from_sealed(SELLADAS_I1, "C2_minkowski_densidad_suave"), "f6_confirmacion.c2_conforme (modos con aux)"),
    ("C0_kasner", "kasner", _from_sealed(SELLADAS_I1, "C0_kasner_vacio"), "f6_confirmacion.c0_kasner"),
    ("D0_kasner", "kasner", _from_sealed(SELLADAS_I2, "D0_kasner_u2"), "f6_intento2.d0_kasner"),
    ("C1_elipsoidal_rotado", "lineal", _from_sealed(SELLADAS_I1, "C1_cono_elipsoidal_rotado"), "f6_confirmacion.c1_elipsoidal"),
    ("D1_cizalla", "lineal", _from_sealed(SELLADAS_I2, "D1_minkowski_cizalla"), "f6_intento2.d1_cizalla"),
    ("D6_cono_tapa_plana", "cono_suave", _from_sealed(SELLADAS_I2, "D6_cono_tapa_plana"), "f6_intento2._gauge_gen(d6_gauge, BOX_STD)"),
    ("D7_cono_huso", "cono_suave", _from_sealed(SELLADAS_I2, "D7_cono_huso"), "f6_intento2._gauge_gen(d7_gauge, BOX_D7)"),
    ("D8_cono_suma_normas", "cono_suave", _from_sealed(SELLADAS_I2, "D8_cono_suma_normas"), "f6_intento2._gauge_gen(d8_gauge, BOX_STD)"),
    ("C3_cono_cilindrico", "cono_suave", _from_sealed(SELLADAS_I1, "C3_cono_cilindrico"), "f6_confirmacion._cone_gen(_within_cilindro)"),
    ("C4_cono_esfera_cubo", "cono_suave", _from_sealed(SELLADAS_I1, "C4_cono_esfera_cubo"), "f6_confirmacion._cone_gen(_within_esfera_cubo)"),
    ("C6_cono_p_por_eje", "cono_suave", _from_sealed(SELLADAS_I1, "C6_cono_p_por_eje"), "f6_confirmacion._cone_gen(_within_p_por_eje)"),
    ("cono_p1_5", "cono_suave", _from_struct("cono_p1_5_4d"), "f6_exploracion.STRUCT['cono_p1_5_4d'] (f7.sprinkle_lp, p = 1.5)"),
    ("cono_p3", "cono_suave", _from_struct("cono_p3_4d"), "f6_exploracion.STRUCT['cono_p3_4d'] (f7.sprinkle_lp, p = 3)"),
    ("cono_p6", "cono_suave", _from_struct("cono_p6_4d"), "f6_exploracion.STRUCT['cono_p6_4d'] (f7.sprinkle_lp, p = 6)"),
    ("cono_cubico", "poliedrico", _from_struct("cono_cubico_4d"), "f6_exploracion.STRUCT['cono_cubico_4d'] (f5prima.sprinkle_polyhedral, max)"),
    ("cono_octaedrico", "poliedrico", _from_struct("cono_octaedrico_4d"), "f6_exploracion.STRUCT['cono_octaedrico_4d'] (f5prima.sprinkle_polyhedral, l1)"),
    ("C5_links_borrados", "defecto", _from_sealed(SELLADAS_I1, "C5_links_borrados_13"), "f6_confirmacion.c5_links_borrados (aux)"),
    ("D4_relaciones_agregadas", "defecto", _from_sealed(SELLADAS_I2, "D4_relaciones_agregadas"), "f6_intento2.d4_agregadas (aux)"),
    ("D5_gemelos", "defecto", _from_sealed(SELLADAS_I2, "D5_gemelos_10"), "f6_intento2.d5_gemelos (aux)"),
    ("producto_4", "producto", _from_struct("producto_4"), "f6_exploracion.STRUCT['producto_4'] (f1f2.random_dorder(4))"),
    ("percolacion_005", "percolacion", _from_struct("percolacion_005"), "f6_exploracion.STRUCT['percolacion_005'] (transitive_percolation, p = 0.05)"),
]
NAMES = [s[0] for s in STRUCTS]
FAMILY = {s[0]: s[1] for s in STRUCTS}
GEN = {s[0]: s[2] for s in STRUCTS}
GEOM_FAMILIES = ("conforme", "kasner", "lineal")
FAMILIES = (*GEOM_FAMILIES, "cono_suave", "poliedrico", "defecto", "producto", "percolacion")
GEOMS = [s[0] for s in STRUCTS if s[1] in GEOM_FAMILIES]
assert NAMES[: len(GEOMS)] == GEOMS


def rng_for(*key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(ENTROPY, spawn_key=(EXPERIMENT_ID, *key))))


def target(name: str) -> str:
    """Referencia global: la propia geometría; impostores, la plana (R14-B v2 §3 y aclaración A1)."""
    return name if FAMILY[name] in GEOM_FAMILIES else "plana"


# --- Observables por sub-intervalo (§1, 1 a 5 y 8) ------------------------------------------------------


def chain_counts(c: Bool, kmax: int | None = None) -> tuple[list[int], int]:
    """(C_1..C_K, L): C_k = número de k-cadenas (v_1 = 1, v_k = Cᵀ v_{k−1}, C_k = Σ v_k) y L = cadena más larga
    (el mayor k con C_k > 0). Vale con cualquier etiquetado (no exige orden natural). Con kmax se corta en C_kmax
    (lista de largo kmax, con ceros) y L vale −1."""
    n = c.shape[0]
    if n == 0:
        return ([0] * kmax if kmax else []), 0
    a = c.astype(np.float64)
    v = np.ones(n)
    out = [n]
    while kmax is None or len(out) < kmax:
        v = a.T @ v
        tot = float(v.sum())
        if tot <= 0:
            break
        out.append(round(tot))
    if kmax is not None:
        return out + [0] * (kmax - len(out)), -1
    return out, len(out)


def longest_chain(c: Bool) -> int:
    return chain_counts(c)[1]


def expected_chains(n: int, k: int) -> float:
    """⟨C_k⟩ = n(n−1)···(n−k+1)·χ_k (n puntos uniformes en el diamante 4D; factorial descendente, VERIFICACION B1)."""
    ff = 1.0
    for i in range(k):
        ff *= n - i
    return ff * CHI[k]


def interval_measures(c_in: Bool, b_in: npt.NDArray[np.integer[Any]], c_tips: Bool) -> dict[str, Any]:
    """Observables de un sub-intervalo. `c_in` y `b_in`: los m elementos interiores (sin puntas) y sus conteos «entre»
    (B global restringido: el intervalo es convexo). `c_tips`: el sub-orden con las dos puntas (para el punto medio).

    - r = relaciones interiores / C(m, 2) (como F6).
    - abundancias N_k (k = 0..m−2) de los pares relacionados del interior; S_{k+2}/n = N_k / R (Eichhorn ec. 16).
    - C_k (k = 2..5) y la cadena más larga L del interior (n = m; ⟨C_k⟩ con n = m).
    - d_mid (C0-4: `midpoint_dimension` sobre el sub-orden con puntas) − d_r (`mm_dimension` del interior).
    - s̄ = (m − N_0 + 9N_1 − 16N_2 + 8N_3)/m: media de s(x) = 1 − N_0(x) + 9N_1(x) − 16N_2(x) + 8N_3(x) sobre el
      interior (sin suavizar; igual a `bd_action_4d(interior)/m`).
    """
    m = c_in.shape[0]
    rel = int(np.count_nonzero(c_in))
    r = rel / (m * (m - 1) / 2)
    vals = b_in[b_in >= 0].astype(np.int64)
    nk = np.bincount(vals, minlength=max(m - 1, 1))[: max(m - 1, 1)]
    counts, longest = chain_counts(c_in)
    ck = [counts[k - 1] if k <= len(counts) else 0 for k in KS]
    n4 = [int(nk[i]) if i < nk.size else 0 for i in range(4)]
    s_mean = (m - n4[0] + 9 * n4[1] - 16 * n4[2] + 8 * n4[3]) / m
    dmid = midpoint_dimension(c_tips)
    dr = mm_dimension(c_in)
    return {"r": r, "R": rel, "N_k": nk, "C_k": ck, "L": longest, "dmid_menos_dr": dmid - dr, "s_media": s_mean}


def _cell_summary(ms: list[dict[str, Any]], m: int) -> dict[str, Any]:
    rs = np.array([x["r"] for x in ms])
    prof = np.array([x["N_k"] / x["R"] for x in ms if x["R"] > 0], dtype=np.float64).reshape(-1, m - 1)
    ratios = np.array([[x["C_k"][i] / expected_chains(m, k) for i, k in enumerate(KS)] for x in ms])
    dd = np.array([x["dmid_menos_dr"] for x in ms], dtype=np.float64)
    return {"r": [round(float(v), 6) for v in rs], "r_media": float(rs.mean()), "r_sd": float(rs.std(ddof=1)),
            "abund_media": [float(v) for v in prof.mean(axis=0)] if prof.shape[0] else None,
            "abund_sin_relaciones": int(len(ms) - prof.shape[0]),
            "cad_media": [float(v) for v in ratios.mean(axis=0)], "cad_sd": [float(v) for v in ratios.std(axis=0, ddof=1)],
            "L_media": float(np.mean([x["L"] for x in ms])),
            "dmid_menos_dr": float(np.nanmean(dd)) if np.any(np.isfinite(dd)) else math.nan,
            "dmid_nan": int(np.count_nonzero(~np.isfinite(dd))),
            "s_media": float(np.mean([x["s_media"] for x in ms]))}


# --- Observables globales (§1, 6 y 7) ----------------------------------------------------------------------


def link_degrees(links: Bool) -> npt.NDArray[np.int64]:
    """d_in + d_out por elemento (= 2·d_j de Eichhorn ec. 14)."""
    out: npt.NDArray[np.int64] = links.sum(axis=0, dtype=np.int64) + links.sum(axis=1, dtype=np.int64)
    return out


def normalized_laplacian(links: Bool) -> Float:
    """Δ = 𝟙 − 𝒟^(−1/2) L_S 𝒟^(−1/2), L_S = L + Lᵀ, 𝒟 = diag(sumas de filas de L_S); 𝒟^(−1/2) := 0 en vértices
    aislados (pseudoinversa de Moore-Penrose), que quedan con autovalor 1."""
    deg = link_degrees(links).astype(np.float64)
    inv = np.zeros_like(deg)
    nz = deg > 0
    inv[nz] = 1.0 / np.sqrt(deg[nz])
    a = links.astype(np.float64)
    a += a.T
    a *= inv[:, None]
    a *= inv[None, :]
    np.negative(a, out=a)
    a[np.diag_indices_from(a)] += 1.0
    return a


def laplacian_spectrum(links: Bool) -> Float:
    """Espectro completo (ascendente) con LAPACK denso (syevd, solo autovalores; misma llamada que `measure`)."""
    w: Float = sla.eigh(normalized_laplacian(links), eigvals_only=True, overwrite_a=True, check_finite=False, driver="evd")
    return w


def spectrum_summary(w: Float) -> dict[str, Any]:
    hist = np.histogram(np.clip(w, 0.0, 2.0), bins=LAP_BINS, range=(0.0, 2.0))[0]
    return {"lambda2": float(w[1]) if w.size > 1 else math.nan, "lambda_max": float(w[-1]), "lambda_min": float(w[0]),
            "fuera_0_2": int(np.count_nonzero((w < -1e-9) | (w > 2 + 1e-9))), "hist": [int(v) for v in hist]}


# --- Medida de una realización ----------------------------------------------------------------------------


def measure(c: Bool, rng: np.random.Generator, ms: tuple[int, ...] = MS, spectrum: bool = True) -> dict[str, Any]:
    """Todos los observables de §1 sobre un orden: B una sola vez; por m local (m ≤ 0.1·m_max, al menos 30
    intervalos con exactamente m interiores) hasta 300 intervalos elegidos con `rng.choice` (como F6 y R14-B); luego
    los links (B == 0), los grados y el espectro del Laplaciano."""
    t0 = time.time()
    n = int(c.shape[0])
    b = between_blocks(c)
    m_max = int(np.max(b)) if np.any(c) else -1
    t_b = time.time() - t0
    por_m: dict[str, Any] = {}
    for m in ms:
        cell: dict[str, Any] = {"local": bool(m <= LOCAL_FRAC * m_max)}
        if cell["local"]:
            ii, jj = np.nonzero(b == m)
            cnt = int(ii.size)
            cell.update({"disponibles": cnt, "medible": cnt >= MIN_INTERVALS})
            if cnt >= MIN_INTERVALS:
                pick = rng.choice(cnt, size=min(MAX_INTERVALS, cnt), replace=False)
                res = []
                for p in pick:
                    idx = interval_suborder(c, int(ii[p]), int(jj[p]))
                    inner = idx[1:-1]
                    res.append(interval_measures(c[np.ix_(inner, inner)], b[np.ix_(inner, inner)], c[np.ix_(idx, idx)]))
                cell.update({"n_sel": int(pick.size), **_cell_summary(res, m)})
        else:
            cell["medible"] = False
        por_m[str(m)] = cell
    t_sub = time.time() - t0 - t_b
    links = b == 0
    del b
    deg = link_degrees(links)
    out: dict[str, Any] = {"N": n, "m_max": m_max, "relaciones": int(np.count_nonzero(c)), "links": int(np.count_nonzero(links)),
                           "por_m": por_m, "grado_bincount": [int(v) for v in np.bincount(deg)],
                           "aislados": int(np.count_nonzero(deg == 0)), "t_B_s": round(t_b, 2), "t_subintervalos_s": round(t_sub, 2)}
    if spectrum:
        t1 = time.time()
        a = normalized_laplacian(links)
        del links
        t2 = time.time()
        w = sla.eigh(a, eigvals_only=True, overwrite_a=True, check_finite=False, driver="evd")
        del a
        out["laplaciano"] = spectrum_summary(w)
        out["t_laplaciano_s"] = round(time.time() - t1, 2)
        out["t_eigh_s"] = round(time.time() - t2, 2)
    out["t_medida_s"] = round(time.time() - t0, 2)
    return out


def mc_rows(m: int, nmc: int) -> dict[str, Any]:
    """Referencia de Monte Carlo por m: nmc diamantes planos 4D de m puntos (más las puntas (±T, 0, 0, 0) para el
    punto medio); mismos observables que un sub-intervalo."""
    tips = np.array([[-DIAMOND_HALF_HEIGHT, 0, 0, 0], [DIAMOND_HALF_HEIGHT, 0, 0, 0]], dtype=np.float64)
    res = []
    for s in range(nmc):
        x = np.vstack([tips[:1], sprinkle_diamond(m, 4, rng_for(K_MC, m, s)), tips[1:]])
        ct = causal_matrix(x)
        if np.count_nonzero(ct[0] & ct[:, m + 1]) != m:
            raise RuntimeError("las puntas no contienen los m puntos")
        b = between_blocks(ct)
        inner = np.arange(1, m + 1)
        res.append(interval_measures(ct[np.ix_(inner, inner)], b[np.ix_(inner, inner)], ct))
    summ = _cell_summary(res, m)
    return {"m": m, "reps": nmc, **summ, "L": [int(x["L"]) for x in res], "C_k": [x["C_k"] for x in res],
            "s": [float(x["s_media"]) for x in res], "dmid_menos_dr_todos": [float(x["dmid_menos_dr"]) for x in res]}


# --- Unidades de trabajo ------------------------------------------------------------------------------------


def mc_path(out: Path, m: int) -> Path:
    return out / "mc" / f"m{m}.json"


def ref_path(out: Path, geom: str, n: int, s: int) -> Path:
    return out / "refs" / geom / f"N{n}_s{s:02d}.json"


def struct_path(out: Path, name: str, n: int, real: int) -> Path:
    return out / "structs" / name / f"N{n}_r{real:02d}.json"


def _write(p: Path, rec: dict[str, Any]) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(_clean(rec), default=str))
    tmp.replace(p)


def run_mc(task: dict[str, Any]) -> str:
    t0 = time.time()
    row = mc_rows(task["m"], task["nmc"])
    _write(mc_path(Path(task["out"]), task["m"]), {"spec": task["spec"], "row": row})
    return f"mc m={task['m']} {time.time() - t0:.1f} s"


def run_ref(task: dict[str, Any]) -> str:
    t0 = time.time()
    geom, n, s = task["geom"], task["n"], task["s"]
    g = GEOMS.index(geom)
    c, extra = GEN[geom](n, rng_for(K_REF_GEN, g, n, s), rng_for(K_REF_AUX, g, n, s))
    t_c = time.time() - t0
    row = {"geometria": geom, "s": s, "t_C_s": round(t_c, 2), **extra, **measure(c, rng_for(K_REF_SEL, g, n, s))}
    _write(ref_path(Path(task["out"]), geom, n, s), {"spec": task["spec"], "row": row})
    return f"ref {geom} N={n} s{s} {time.time() - t0:.1f} s (C {t_c:.1f} s, lap {row.get('t_laplaciano_s')} s)"


def run_struct(task: dict[str, Any]) -> str:
    t0 = time.time()
    name, n, real = task["name"], task["n"], task["real"]
    k = NAMES.index(name)
    c, extra = GEN[name](n, rng_for(K_EST_GEN, k, n, real), rng_for(K_EST_AUX, k, n, real))
    t_c = time.time() - t0
    row = {"estructura": name, "realizacion": real, "t_C_s": round(t_c, 2), "extra": extra,
           **measure(c, rng_for(K_EST_SEL, k, n, real))}
    _write(struct_path(Path(task["out"]), name, n, real), {"spec": task["spec"], "row": row})
    return f"est {name} N={n} r{real} {time.time() - t0:.1f} s (C {t_c:.1f} s, lap {row.get('t_laplaciano_s')} s)"


def _work(task: dict[str, Any]) -> str:
    kind = task["kind"]
    if kind == "mc":
        return run_mc(task)
    return run_ref(task) if kind == "ref" else run_struct(task)


def _cost(task: dict[str, Any]) -> float:
    """Costo relativo para lanzar primero lo más largo: ∝ N³ (espectro), Kasner y defectos más caros."""
    if task["kind"] == "mc":
        return 1e9  # antes que nada: el resumen las necesita y son baratas
    name = task.get("geom") or task.get("name") or ""
    w = 1.5 if "kasner" in name else (1.2 if FAMILY.get(name) == "defecto" else 1.0)
    return w * (float(task["n"]) / 1000.0) ** 3


def _pool(tasks: list[dict[str, Any]], workers: int, budget_s: float | None) -> int:
    """Corre las tareas con `workers` procesos (spawn: heredan OPENBLAS_NUM_THREADS); con presupuesto deja de lanzar
    tareas nuevas al agotarlo."""
    t0 = time.time()
    done = 0

    def over() -> bool:
        return budget_s is not None and time.time() - t0 > budget_s

    with mp.get_context("spawn").Pool(max(1, min(workers, len(tasks)))) as pool:
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
    override = any(v is not None for v in (args.ns, args.reals, args.nref, args.nmc))
    if override and not args.out:
        sys.exit("--ns/--reals/--nref/--nmc exigen --out (no se mezcla con la corrida oficial)")
    head, dirty = git_info()
    if dirty and not args.allow_dirty:
        sys.exit("exige arbol git limpio fuera de results/ (usar --allow-dirty solo para pruebas)")
    if not SPEC.exists():
        sys.exit(f"falta {SPEC}")
    ns = tuple(int(v) for v in args.ns.split(",")) if args.ns else NS
    reals = args.reals or N_REAL
    nref = args.nref or N_REFS
    nmc = args.nmc or N_MC
    official = not dirty and not override and out == OUT
    spec = {"commit": head, "dirty": dirty, "official": official, "ns": list(ns), "reals": reals, "nref": nref, "nmc": nmc,
            "preregistro_sha256": spec_sha(), "experiment_id": EXPERIMENT_ID, "entropia": ENTROPY,
            "blas_threads": args.blas_threads}
    manifest = out / "manifest.json"
    keys = ("ns", "reals", "nref", "nmc", "official", "experiment_id", "entropia")
    if manifest.exists():
        old = json.loads(manifest.read_text())
        if any(old[k] != spec[k] for k in keys):
            sys.exit(f"{manifest} tiene otra configuración: {old}")
    else:
        _write(manifest, {k: spec[k] for k in keys})
    sel_n = [int(v) for v in args.n.split(",")] if args.n else list(ns)
    if any(v not in ns for v in sel_n):
        sys.exit(f"--n fuera de {ns}")
    tasks: list[dict[str, Any]] = []
    base = {"out": str(out), "spec": spec}
    if args.part in ("mc", "all"):
        tasks += [{"kind": "mc", "m": m, "nmc": nmc, **base} for m in MS if not mc_path(out, m).exists()]
    if args.part in ("refs", "all"):
        for geom in _names(args.geom, GEOMS):
            for n in sel_n:
                for s in _range(args.s, nref):
                    if not ref_path(out, geom, n, s).exists():
                        tasks.append({"kind": "ref", "geom": geom, "n": n, "s": s, **base})
    if args.part in ("structs", "all"):
        for name in _names(args.struct, NAMES):
            for n in sel_n:
                for real in _range(args.real, reals):
                    if not struct_path(out, name, n, real).exists():
                        tasks.append({"kind": "est", "name": name, "n": n, "real": real, **base})
    tasks.sort(key=lambda t: -_cost(t))
    os.environ["OPENBLAS_NUM_THREADS"] = os.environ["OMP_NUM_THREADS"] = os.environ["MKL_NUM_THREADS"] = str(args.blas_threads)
    print(f"{out} (commit {head[:8]}, dirty {dirty}, oficial {official}, hilos BLAS {args.blas_threads}): "
          f"{len(tasks)} unidades por correr")
    t0 = time.time()
    done = _pool(tasks, args.workers, args.budget_min * 60 if args.budget_min else None) if tasks else 0
    print(f"{done} de {len(tasks)} unidades; pared {time.time() - t0:.1f} s"
          + ("" if done == len(tasks) else "; repetir el comando para seguir"))


def cmd_cost(args: argparse.Namespace) -> None:
    """Mide (sin escribir resultados) una realización completa y el Laplaciano aparte, con la entropía de una semilla
    de prueba (no la del experimento)."""
    n = int(args.n or 8000)
    name = args.struct or "plana"
    r = np.random.default_rng(np.random.SeedSequence(990324, spawn_key=(n,)))
    t0 = time.time()
    c, _ = GEN[name](n, r, r)
    t_c = time.time() - t0
    row = measure(c, r)
    t_all = time.time() - t0
    print(json.dumps({"estructura": name, "N": n, "hilos_blas": os.environ.get("OPENBLAS_NUM_THREADS"), "t_generar_s": round(t_c, 1),
                      "t_B_s": row["t_B_s"], "t_subintervalos_s": row["t_subintervalos_s"],
                      "t_laplaciano_s": row.get("t_laplaciano_s"), "t_eigh_s": row.get("t_eigh_s"),
                      "t_total_realizacion_s": round(t_all, 1), "lambda2": row["laplaciano"]["lambda2"],
                      "lambda_max": row["laplaciano"]["lambda_max"], "links": row["links"], "aislados": row["aislados"],
                      "medibles": {m: v["medible"] for m, v in row["por_m"].items()}}, indent=1))


# --- Resumen ---------------------------------------------------------------------------------------------------

# Observables escalares por realización: (nombre, m o None). Orden = índice de observable del bootstrap.
SUB_OBS = ("r_media", "r_sd", "r_ED", "abund_ED", "cad2", "cad3", "cad4", "cad5", "L_ratio", "dmid_menos_dr", "bd_s_menos_ref")
GLOB_OBS = ("pendiente_A", "grado_ED", "lap_lambda2", "lap_lambda_max", "lap_ED")
OBS: list[tuple[str, int | None]] = [(o, m) for m in MS for o in SUB_OBS] + [(o, None) for o in GLOB_OBS]


def obs_label(o: str, m: int | None) -> str:
    return o if m is None else f"{o}@m{m}"


def _f(v: Any) -> float:
    return math.nan if v is None or (isinstance(v, str)) else float(v)


def _load_json(p: Path) -> dict[str, Any]:
    data: dict[str, Any] = json.loads(p.read_text())
    return data


def ed_weighted(u: Float, uw: Float, v: Float, vw: Float) -> float:
    """Distancia de energía (Székely) entre dos distribuciones 1D con pesos: √(2∫(F − G)² dx) (`scipy.stats`)."""
    if uw.sum() <= 0 or vw.sum() <= 0:
        return math.nan
    return float(energy_distance(u, v, uw, vw))


def slope_a(row: dict[str, Any], n: int) -> float:
    """Pendiente A de F6 por realización: MCO de Δr(m) = r_media(m) − 1/10 frente a x = √(m/N), con los m medibles
    (al menos 2)."""
    xs, ys = [], []
    for m in MS:
        cell = row["por_m"][str(m)]
        if cell.get("medible"):
            xs.append(math.sqrt(m / n))
            ys.append(float(cell["r_media"]) - R_PLANO)
    if len(xs) < 2:
        return math.nan
    return float(np.polyfit(np.array(xs), np.array(ys), 1)[0])


class MCRef:
    def __init__(self, rows: dict[int, dict[str, Any]]) -> None:
        self.rows = rows
        self.r = {m: np.array(rows[m]["r"]) for m in rows}
        self.abund = {m: np.array(rows[m]["abund_media"]) for m in rows}
        self.L = {m: float(np.mean(rows[m]["L"])) for m in rows}
        self.s = {m: float(np.mean(rows[m]["s"])) for m in rows}


def scalar_obs(row: dict[str, Any], n: int, mc: MCRef) -> dict[str, float]:
    """Observables escalares de una realización (salvo los globales de distancia de energía, que dependen del
    conjunto de referencia y se calculan aparte)."""
    out: dict[str, float] = {}
    for m in MS:
        cell = row["por_m"][str(m)]
        keys = [obs_label(o, m) for o in SUB_OBS]
        if not cell.get("medible"):
            out.update({k: math.nan for k in keys})
            continue
        r = np.array([_f(v) for v in cell["r"]])
        prof = np.array([_f(v) for v in cell["abund_media"] or []])
        ks: Float = np.arange(prof.size, dtype=np.float64)
        out[obs_label("r_media", m)] = float(cell["r_media"])
        out[obs_label("r_sd", m)] = float(cell["r_sd"])
        out[obs_label("r_ED", m)] = float(energy_distance(r, mc.r[m]))
        out[obs_label("abund_ED", m)] = ed_weighted(ks, prof, ks, mc.abund[m]) if prof.size else math.nan
        for i, k in enumerate(KS):
            out[obs_label(f"cad{k}", m)] = float(cell["cad_media"][i])
        out[obs_label("L_ratio", m)] = float(cell["L_media"]) / mc.L[m]
        out[obs_label("dmid_menos_dr", m)] = _f(cell["dmid_menos_dr"])
        out[obs_label("bd_s_menos_ref", m)] = float(cell["s_media"]) - mc.s[m]
    out["pendiente_A"] = slope_a(row, n)
    lap = row.get("laplaciano")
    out["lap_lambda2"] = _f(lap["lambda2"]) if lap else math.nan
    out["lap_lambda_max"] = _f(lap["lambda_max"]) if lap else math.nan
    return out


def _vec(counts: list[int], size: int) -> Float:
    v = np.zeros(size)
    v[: len(counts)] = counts
    return v


class GlobalRef:
    """Distribuciones agrupadas (grado de links y espectro) de las referencias de una geometría en un N; para una
    referencia se excluye ella misma (dejar uno fuera)."""

    def __init__(self, rows: list[dict[str, Any]]) -> None:
        self.dmax = max(len(r["grado_bincount"]) for r in rows) + 1
        self.deg = np.array([_vec(r["grado_bincount"], self.dmax) for r in rows])
        self.has_lap = all("laplaciano" in r for r in rows)
        self.hist = np.array([r["laplaciano"]["hist"] for r in rows], dtype=np.float64) if self.has_lap else None
        self.deg_tot = self.deg.sum(axis=0)
        self.hist_tot = self.hist.sum(axis=0) if self.hist is not None else None

    def ed(self, row: dict[str, Any], leave_out: int | None = None) -> tuple[float, float]:
        size = max(self.dmax, len(row["grado_bincount"]))
        u = _vec(row["grado_bincount"], size)
        pool = _vec(list(self.deg_tot), size)
        if leave_out is not None:
            pool -= _vec(list(self.deg[leave_out]), size)
        x: Float = np.arange(size, dtype=np.float64) / 2.0  # d_j = ½(d_in + d_out) (Eichhorn ec. 14)
        g = ed_weighted(x, u, x, pool)
        lap = row.get("laplaciano")
        if lap is None or self.hist_tot is None or self.hist is None:
            return g, math.nan
        centers: Float = (np.arange(LAP_BINS, dtype=np.float64) + 0.5) * (2.0 / LAP_BINS)
        hp = self.hist_tot - (self.hist[leave_out] if leave_out is not None else 0.0)
        return g, ed_weighted(centers, np.array(lap["hist"], dtype=np.float64), centers, hp)


def holm(p: list[float]) -> list[float]:
    """p ajustados de Holm (paso descendente, monótonos, ≤ 1); nan queda nan y no cuenta en el número de pruebas."""
    idx = [i for i, v in enumerate(p) if math.isfinite(v)]
    order = sorted(idx, key=lambda i: p[i])
    k = len(order)
    out = [math.nan] * len(p)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (k - rank) * p[i]))
        out[i] = run
    return out


def family_effect(zs: list[Float], rng: np.random.Generator, n_boot: int = N_BOOT) -> dict[str, Any]:
    """Δ de la familia = media (entre estructuras) de la media de z por realización. Bootstrap jerárquico:
    estructuras con reposición y, dentro de cada una, realizaciones con reposición. p bilateral de Wald con el error
    estándar max(se del bootstrap, se bajo H0) (permite p pequeños para Holm sin fiarse de un bootstrap con pocas
    unidades); IC 95% por percentiles del bootstrap."""
    zs = [z for z in zs if z.size]
    if not zs:
        return {"delta": math.nan, "ic95": [math.nan, math.nan], "se": math.nan, "se_boot": math.nan, "se_h0": math.nan,
                "p": math.nan, "n_estructuras": 0}
    s = len(zs)
    point = float(np.mean([z.mean() for z in zs]))
    per = np.empty((s, n_boot))
    for j, z in enumerate(zs):
        per[j] = z[rng.integers(0, z.size, size=(n_boot, z.size))].mean(axis=1)
    pick = rng.integers(0, s, size=(n_boot, s))
    boots = per[pick, np.arange(n_boot)[:, None]].mean(axis=1)
    se_boot = float(boots.std(ddof=1))
    # piso: desviación de Δ bajo H0 (cada z ~ N(0, 1) independiente, en unidades de la referencia): evita p diminutos
    # cuando el bootstrap subestima la dispersión (pocas realizaciones o pocas estructuras)
    se_h0 = math.sqrt(sum(1.0 / z.size for z in zs)) / s
    se = max(se_boot, se_h0)
    p = float(2 * norm.sf(abs(point) / se))
    return {"delta": point, "ic95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))], "se": se,
            "se_boot": se_boot, "se_h0": se_h0, "p": p, "n_estructuras": s}


def load_all(out: Path) -> tuple[dict[str, Any], dict[int, dict[str, Any]], dict[tuple[str, int], list[dict[str, Any]]],
                                  dict[tuple[str, int], list[dict[str, Any]]], set[tuple[str, bool, bool, str]]]:
    man = _load_json(out / "manifest.json")
    prov: set[tuple[str, bool, bool, str]] = set()

    def get(p: Path) -> dict[str, Any]:
        if not p.exists():
            raise SystemExit(f"falta {p}")
        d = _load_json(p)
        prov.add((d["spec"]["commit"], d["spec"]["dirty"], d["spec"]["official"], d["spec"]["preregistro_sha256"]))
        row: dict[str, Any] = d["row"]
        return row

    mc = {m: get(mc_path(out, m)) for m in MS}
    refs = {(g, n): [get(ref_path(out, g, n, s)) for s in range(man["nref"])] for g in GEOMS for n in man["ns"]}
    structs = {(k, n): [get(struct_path(out, k, n, r)) for r in range(man["reals"])] for k in NAMES for n in man["ns"]}
    return man, mc, refs, structs, prov


def chain_check(mc: dict[int, dict[str, Any]]) -> dict[str, Any]:
    """⟨C_k⟩ de Monte Carlo frente a n(n−1)···(n−k+1)·χ_k (z con el error estándar de la media)."""
    out: dict[str, Any] = {}
    for m, row in mc.items():
        ck = np.array(row["C_k"], dtype=np.float64)
        for i, k in enumerate(KS):
            mean, se = float(ck[:, i].mean()), float(ck[:, i].std(ddof=1) / math.sqrt(ck.shape[0]))
            th = expected_chains(m, k)
            out[f"m{m}_k{k}"] = {"teoria": th, "mc": mean, "se": se, "z": (mean - th) / se if se > 0 else math.nan}
    return out


def summarize(out: Path) -> dict[str, Any]:
    man, mc_rows_, refs, structs, prov = load_all(out)
    mc = MCRef(mc_rows_)
    # mínimos de validez: 20 de 60 y 5 de 10 en la corrida oficial; en pruebas, la misma proporción
    min_ref = MIN_REF_VALID if man["nref"] == N_REFS else max(2, math.ceil(man["nref"] / 3))
    min_real = MIN_REAL_VALID if man["reals"] == N_REAL else max(1, math.ceil(man["reals"] / 2))
    labels = [obs_label(o, m) for o, m in OBS]
    # valores escalares por realización: vals[(estructura o geometría-ref, N)] -> (realizaciones × observables)
    ref_vals: dict[tuple[str, int], Float] = {}
    gref: dict[tuple[str, int], GlobalRef] = {}
    for (g, n), rows in refs.items():
        gr = GlobalRef(rows)
        gref[(g, n)] = gr
        mat = np.full((len(rows), len(labels)), np.nan)
        for i, row in enumerate(rows):
            sc = scalar_obs(row, n, mc)
            sc["grado_ED"], sc["lap_ED"] = gr.ed(row, leave_out=i)
            mat[i] = [sc[lb] for lb in labels]
        ref_vals[(g, n)] = mat
    # valores de las estructuras frente a cada referencia de geometría (los ED globales dependen de la referencia)
    st_vals: dict[tuple[str, int, str], Float] = {}
    for (name, n), rows in structs.items():
        base = [scalar_obs(row, n, mc) for row in rows]
        for g in GEOMS:
            mat = np.full((len(rows), len(labels)), np.nan)
            for i, row in enumerate(rows):
                sc = dict(base[i])
                sc["grado_ED"], sc["lap_ED"] = gref[(g, n)].ed(row)
                mat[i] = [sc[lb] for lb in labels]
            st_vals[(name, n, g)] = mat

    zcache: dict[tuple[str, int, str], Float] = {}

    def zscores(name: str, n: int, g: str) -> Float:
        if (name, n, g) in zcache:
            return zcache[(name, n, g)]
        ref = ref_vals[(g, n)]
        x = st_vals[(name, n, g)]
        z = np.full(x.shape, np.nan)
        for j in range(len(labels)):
            col = ref[:, j][np.isfinite(ref[:, j])]
            if col.size >= min_ref and float(col.std(ddof=1)) > 0:
                z[:, j] = (x[:, j] - col.mean()) / col.std(ddof=1)
        zcache[(name, n, g)] = z
        return z

    z_target: dict[tuple[str, int], Float] = {}
    per_struct: dict[str, Any] = {}
    for name in NAMES:
        per_struct[name] = {"familia": FAMILY[name], "referencia": target(name), "por_N": {}}
        for n in man["ns"]:
            z = zscores(name, n, target(name))
            z_target[(name, n)] = z
            raw = st_vals[(name, n, target(name))]
            per_struct[name]["por_N"][str(n)] = {
                lb: {"media_z": float(np.nanmean(z[:, j])) if np.any(np.isfinite(z[:, j])) else None,
                     "validas": int(np.count_nonzero(np.isfinite(z[:, j]))),
                     "media_valor": float(np.nanmean(raw[:, j])) if np.any(np.isfinite(raw[:, j])) else None}
                for j, lb in enumerate(labels)}
            rows = structs[(name, n)]
            lap = [r["laplaciano"]["hist"] for r in rows if "laplaciano" in r]
            if lap:
                h = np.array(lap, dtype=np.float64).mean(axis=0).reshape(50, -1).sum(axis=1)
                per_struct[name]["por_N"][str(n)]["histograma_laplaciano_50"] = (h / h.sum()).tolist()
    table: list[dict[str, Any]] = []
    for ni, n in enumerate(man["ns"]):
        for fi, fam in enumerate(FAMILIES):
            names = [k for k in NAMES if FAMILY[k] == fam]
            for oi, lb in enumerate(labels):
                zs: list[Float] = []
                for k in names:
                    col = z_target[(k, n)][:, oi]
                    col = col[np.isfinite(col)]
                    if col.size >= min_real:
                        zs.append(col)
                eff = family_effect(zs, rng_for(K_BOOT, 10000 * ni + 100 * fi + oi))
                rec: dict[str, Any] = {"N": n, "familia": fam,
                                       "comparacion": "propia_geometria" if fam in GEOM_FAMILIES else "plana",
                                       "observable": lb, **eff, "estructuras": len(names)}
                if fam not in GEOM_FAMILIES:
                    # informativo: Δ frente a cada geometría (punto, sin bootstrap); la más cercana en |Δ|
                    best = (math.inf, "")
                    for g in GEOMS:
                        means = []
                        for k in names:
                            zg = zscores(k, n, g)
                            col = zg[:, oi][np.isfinite(zg[:, oi])]
                            if col.size >= min_real:
                                means.append(float(col.mean()))
                        if means and abs(float(np.mean(means))) < best[0]:
                            best = (abs(float(np.mean(means))), g)
                    rec["min_abs_delta_geometrias"] = best[0] if best[1] else math.nan
                    rec["geometria_mas_cercana"] = best[1]
                table.append(rec)
    padj = holm([r["p"] for r in table])
    for r, pa in zip(table, padj):
        r["p_holm"] = pa
        sig = math.isfinite(pa) and pa < ALPHA
        r["holm_significativo"] = sig
        r["separa"] = bool(sig and r["familia"] not in GEOM_FAMILIES)
    header = ["N", "familia", "comparacion", "observable", "estructuras", "n_estructuras", "delta", "ic95_lo", "ic95_hi", "se",
              "se_boot", "se_h0", "p", "p_holm", "holm_significativo", "separa", "min_abs_delta_geometrias", "geometria_mas_cercana"]
    with (out / "tabla_r14c.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in table:
            w.writerow([r["N"], r["familia"], r["comparacion"], r["observable"], r["estructuras"], r["n_estructuras"], r["delta"],
                        r["ic95"][0], r["ic95"][1], r["se"], r["se_boot"], r["se_h0"], r["p"], r["p_holm"], r["holm_significativo"], r["separa"],
                        r.get("min_abs_delta_geometrias", ""), r.get("geometria_mas_cercana", "")])
    commits = sorted({p[0] for p in prov})
    refsum = {f"{g}:{n}": {lb: {"media": float(np.nanmean(v[:, j])) if np.any(np.isfinite(v[:, j])) else None,
                                "sd": float(np.nanstd(v[:, j], ddof=1)) if np.count_nonzero(np.isfinite(v[:, j])) > 1 else None,
                                "validas": int(np.count_nonzero(np.isfinite(v[:, j])))} for j, lb in enumerate(labels)}
              for (g, n), v in ref_vals.items()}
    geo_sig = [r for r in table if r["familia"] in GEOM_FAMILIES and r["holm_significativo"]]
    return {"estado": "EXPLORATORIA (R14-C v2), sin voto", "preregistro": str(SPEC), "preregistro_sha256_actual": spec_sha(),
            "experiment_id": EXPERIMENT_ID, "entropia": ENTROPY, "manifest": man,
            "minimos_validez": {"referencias": min_ref, "realizaciones": min_real}, "n_boot": N_BOOT, "alfa_holm": ALPHA,
            "procedencia": {"commits": commits, "un_solo_commit": len(commits) == 1, "alguna_sucia": any(p[1] for p in prov),
                            "todas_oficiales": all(p[2] for p in prov), "sha_preregistro": sorted({p[3] for p in prov})},
            "observables": labels, "familias": {f: [k for k in NAMES if FAMILY[k] == f] for f in FAMILIES},
            "generadores": {s[0]: s[3] for s in STRUCTS},
            "referencia_mc": {"L_ref": mc.L, "s_ref": mc.s, "r_media": {m: float(mc.r[m].mean()) for m in MS},
                              "chequeo_cadenas": chain_check(mc_rows_)},
            "referencias_globales": refsum, "estructuras": per_struct,
            "tabla": table, "pruebas_holm": int(sum(math.isfinite(r["p"]) for r in table)),
            "separa": [f"{r['familia']} N={r['N']} {r['observable']}" for r in table if r["separa"]],
            "geometrias_significativas_contra_su_referencia": [f"{r['familia']} N={r['N']} {r['observable']}" for r in geo_sig]}


def cmd_summarize(args: argparse.Namespace) -> None:
    out = Path(args.out) if args.out else OUT
    s = summarize(out)
    (out / "summary.json").write_text(json.dumps(_clean(s), indent=1, default=str))
    pv = s["procedencia"]
    print(f"commits {pv['commits']} oficial {pv['todas_oficiales']} sucia {pv['alguna_sucia']}")
    print(f"filas de la tabla: {len(s['tabla'])}; pruebas en Holm: {s['pruebas_holm']}")
    print("chequeo de cadenas (z):", {k: round(v["z"], 2) for k, v in s["referencia_mc"]["chequeo_cadenas"].items()})
    print(f"separa (Holm < {ALPHA}, familias de impostores): {len(s['separa'])}")
    for x in s["separa"][:40]:
        print("  ", x)
    print("geometrías significativas contra su propia referencia (deberían ser ninguna):", s["geometrias_significativas_contra_su_referencia"])


# --- Figura ----------------------------------------------------------------------------------------------------

FAM_LABEL = {"conforme": "conformes", "kasner": "Kasner", "lineal": "lineales", "cono_suave": "conos suaves",
             "poliedrico": "poliédricos", "defecto": "defectos", "producto": "producto 4", "percolacion": "percolación"}


def cmd_figure(args: argparse.Namespace) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm

    out = Path(args.out) if args.out else OUT
    fig_path = Path(args.fig) if args.fig else FIG
    s = _load_json(out / "summary.json")
    ns = s["manifest"]["ns"]
    labels = s["observables"]
    fig, axes = plt.subplots(1, len(ns), figsize=(4.2 * len(ns) + 1.5, 0.19 * len(labels) + 2.2), sharey=True, squeeze=False)
    norm_ = SymLogNorm(linthresh=1.0, vmin=-30, vmax=30, base=10)
    im = None
    for ax, n in zip(axes[0], ns):
        mat = np.full((len(labels), len(FAMILIES)), np.nan)
        sig = np.zeros_like(mat, dtype=bool)
        for r in s["tabla"]:
            if r["N"] != n:
                continue
            i, j = labels.index(r["observable"]), FAMILIES.index(r["familia"])
            mat[i, j] = np.nan if r["delta"] is None else r["delta"]
            sig[i, j] = bool(r["holm_significativo"])
        cmap = matplotlib.colormaps["RdBu_r"].copy()
        cmap.set_bad("0.82")  # sin dato (m no local o celdas sin validez)
        im = ax.imshow(np.ma.masked_invalid(np.clip(mat, -30, 30)), aspect="auto", cmap=cmap, norm=norm_, interpolation="nearest")
        yy, xx = np.nonzero(sig)
        ax.scatter(xx, yy, s=7, c="k", marker="o", linewidths=0)
        ax.set_xticks(range(len(FAMILIES)))
        ax.set_xticklabels([FAM_LABEL[f] for f in FAMILIES], rotation=60, ha="right", fontsize=7)
        ax.axvline(len(GEOM_FAMILIES) - 0.5, color="k", lw=0.8)
        for k in range(1, len(MS) + 1):
            ax.axhline(k * len(SUB_OBS) - 0.5, color="0.5", lw=0.4)
        ax.set_title(f"N = {n}", fontsize=9)
    axes[0][0].set_yticks(range(len(labels)))
    axes[0][0].set_yticklabels(labels, fontsize=5.5)
    if im is not None:
        cb = fig.colorbar(im, ax=axes[0].tolist(), shrink=0.6, pad=0.02)
        cb.set_label("Δ de la familia (desv. estándar de la referencia; recortado a ±30)", fontsize=7)
        cb.ax.tick_params(labelsize=6)
    fig.suptitle("R14-C v2 (exploratoria, sin voto): tamaño de efecto por familia y observable.\n"
                 "Geometrías frente a su propia referencia; impostores frente a la plana. Punto: Holm < 0.05 en toda la tabla. Gris: sin dato.",
                 fontsize=8)
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=140, bbox_inches="tight")
    print(fig_path)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--part", choices=("mc", "refs", "structs", "all"), default="all")
    r.add_argument("--struct", help="estructuras separadas por coma (por defecto, las 28)")
    r.add_argument("--geom", help="geometrías de referencia separadas por coma (por defecto, las 12)")
    r.add_argument("--n", help="N separados por coma (subconjunto de los del manifiesto)")
    r.add_argument("--real", help="realizaciones a-b de las estructuras")
    r.add_argument("--s", help="sprinklings a-b de las referencias")
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--blas-threads", type=int, default=1, help="hilos BLAS/LAPACK por proceso (workers × hilos ≤ núcleos)")
    r.add_argument("--budget-min", type=float, help="no lanza unidades nuevas pasado este tiempo (tandas < 2 h)")
    r.add_argument("--ns", help="N de prueba separados por coma (exige --out)")
    r.add_argument("--reals", type=int)
    r.add_argument("--nref", type=int)
    r.add_argument("--nmc", type=int)
    r.add_argument("--out")
    r.add_argument("--allow-dirty", action="store_true")
    s = sub.add_parser("summarize")
    s.add_argument("--out")
    f = sub.add_parser("figure")
    f.add_argument("--out")
    f.add_argument("--fig")
    c = sub.add_parser("cost")
    c.add_argument("--n", type=int)
    c.add_argument("--struct")
    args = ap.parse_args()
    {"run": cmd_run, "summarize": cmd_summarize, "figure": cmd_figure, "cost": cmd_cost}[args.cmd](args)


if __name__ == "__main__":
    main()
