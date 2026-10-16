"""F6 confirmatoria (docs/OMEGA_F6_CONFIRMACION_PREREGISTRO.md, R12-1). Partes: ref | nula_<N> | cal_<g> |
vista_<nombre> | sellada_<k> | final.

Instrumento idéntico al de la exploración (``tools/f6_exploracion.py``, que se importa): matriz causal por
bloques, B, m_max, m ∈ {10, 20, 40, 80, 160} con m ≤ 0.1·m_max, hasta 300 intervalos con m interiores.
Sin dilución ni N = 16000.

- ``ref``: referencia exacta por m (4000 diamantes planos 4D de m puntos más las puntas), entropía 20261078.
- ``nula_<N>``: 30 diamantes planos 4D por N (entropía 20261080, familia 1).
- ``cal_<g>``, g = 0..14: un grupo de calibración (10 realizaciones en N = 4000 y 10 en N = 8000; entropía
  20261080, familia 2). Grupos 0..9 diamantes planos, 10..12 de Sitter (η ∈ [−1.5, −0.5]), 13..14 FRW polvo.
- ``vista_<nombre>``: estructura ya vista con semillas nuevas (entropía 20261081, familia 3).
- ``sellada_<k>``, k = 0..7: estructura del archivo sellado (``docs/SELLADO_F6_CONFIRMACION_CIEGA.md``, hash
  verificado). Puntos y relación con la entropía 20261076, azar auxiliar (modos de C2, links borrados de C5)
  con la 20261077, ambas con spawn_key (321, k, N, s).
- ``final``: b y A con IC por bootstrap, b_plano, τ_b, clasificación, veredicto y lo informativo (§5).

Partición: cualquier parte de estructura (``cal_<g>``, ``vista_<nombre>``, ``sellada_<k>``) admite los sufijos
``_<N>`` (un solo N) y ``_<N>_r<s>`` (una sola realización): mismas semillas, mismas filas. ``final`` usa,
por estructura y N, la parte completa si existe; si no, la de ese N; si no, las de cada realización.
``--listar`` imprime la partición recomendada (cada parte por debajo de unos 9 minutos).

Semillas: SeedSequence(entropía, spawn_key=(321, familia, índice, N, realización)); familia 0 = referencia
(índice 0, N ← m), 1 = nula (índice 0), 2 = calibración (índice = grupo), 3 = vistas (índice = orden en
``VISTAS``), 4 = bootstrap (spawn_key (321, 4, j) con j = g para la calibración, 100 + índice para las vistas y
200 + k para las selladas; entropía 20261080). Las selladas usan (321, k, N, s) sin familia, como fija su archivo.
La selección de intervalos de cada realización sigue con el generador de puntos (como en la exploración); en las
selladas, el de la entropía 20261076. La 20261079 queda reservada y no se usa (§7, A1).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable

import numpy as np
import numpy.typing as npt

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from f6_exploracion import (  # noqa: E402
    STRUCT,
    _cells,
    _dloc_rows,
    _dr_matrix,
    _max_t,
    _stat,
    between_blocks,
    causal_blocks,
    fit_poly,
    interval_obs,
    measure_order,
)
from f7 import de_sitter_eta_points, frw_power_points  # noqa: E402

from omega.causality.kasner import kasner_causal_matrix, sprinkle_kasner_interval, validate_kasner_exponents  # noqa: E402
from omega.causality.orders import is_partial_order  # noqa: E402
from omega.causality.sprinkling import DIAMOND_HALF_HEIGHT, causal_matrix, sprinkle_diamond  # noqa: E402

EXPERIMENT_ID = 321
ENT_SELLADA_GEN = 20261076
ENT_SELLADA_AUX = 20261077
ENT_REF = 20261078
ENT_ANALISIS = 20261080  # calibración, nula y bootstrap
ENT_VISTAS = 20261081
FAM_REF, FAM_NULA, FAM_CAL, FAM_VISTA, FAM_BOOT = 0, 1, 2, 3, 4
BOOT_VISTA, BOOT_SELLADA = 100, 200
M_GRID = (10, 20, 40, 80, 160)
NS = (4000, 8000)
N_REAL = 10
N_REF = 4000
N_NULA = 30
N_BOOT = 1000
N_CAL = 15
TAU_FACTOR = 1.5
BLOCK_ROWS = 1024
T_HALF = DIAMOND_HALF_HEIGHT
PREREG = Path("docs") / "OMEGA_F6_CONFIRMACION_PREREGISTRO.md"
SELLADO = Path("docs") / "SELLADO_F6_CONFIRMACION_CIEGA.md"
SELLADO_SHA256 = "b8439d85fcf2ea3c6ea76b92cce8caced98483d85dac64bcb7ca43f3fbc91c79"
FIRMA = "R12-1 | **Firmada"
OUT = Path("results") / "omega_f6_confirmacion"
Bool = npt.NDArray[np.bool_]
Float = npt.NDArray[np.float64]
Within = Callable[[Float, Float, Float, Float], Bool]
SealedGen = Callable[[int, np.random.Generator, np.random.Generator], tuple[Bool, Float | None, dict[str, Any]]]


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    """SeedSequence(entropía, spawn_key=(321, *key))."""
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


# --- Generadores sellados ---------------------------------------------------------------------------

KASNER_P = (-6.0 / 19.0, 10.0 / 19.0, 15.0 / 19.0)
KASNER_T_LO, KASNER_T_HI, KASNER_K = 1.0, 1.9, 800


def c0_kasner(n: int, r: np.random.Generator, aux: np.random.Generator, k: int = KASNER_K) -> tuple[Bool, Float | None, dict[str, Any]]:
    """C0: Kasner de vacío p = (−6/19, 10/19, 15/19), intervalo t ∈ [1.0, 1.9], K = 800 direcciones."""
    p = validate_kasner_exponents(KASNER_P)
    x = sprinkle_kasner_interval(n, p, KASNER_T_LO, KASNER_T_HI, r, k=k)
    return kasner_causal_matrix(x, p, k=k), x, {}


def _rot_x(a: float) -> Float:
    c, s = math.cos(a), math.sin(a)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def _rot_y(a: float) -> Float:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def _rot_z(a: float) -> Float:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


C1_C = (0.70, 1.00, 1.45)
C1_R: Float = _rot_z(0.6) @ _rot_y(0.9) @ _rot_x(0.3)
C1_A: Float = np.diag([1.0 / c for c in C1_C]) @ C1_R.T


def c1_transform(x: Float) -> Float:
    """(t, x) → (t, A x): en estas coordenadas la relación de C1 es la de Minkowski."""
    out = x.copy()
    out[:, 1:] = x[:, 1:] @ C1_A.T
    return out


def c1_elipsoidal(n: int, r: np.random.Generator, aux: np.random.Generator) -> tuple[Bool, Float | None, dict[str, Any]]:
    """C1: diamante estándar de coordenadas; Δt > 0 y Δt ≥ ‖AΔx‖₂, construida como Minkowski tras x' = Ax
    (A lineal: AΔx = Ax_j − Ax_i salvo redondeo)."""
    x = sprinkle_diamond(n, 4, r)
    return causal_blocks(c1_transform(x), "mink"), x, {}


C2_MODES, C2_AMP, C2_KMAX, C2_SHIFT = 5, 0.07, 1.5, 0.35


def c2_modes(aux: np.random.Generator) -> tuple[Float, Float]:
    """kⱼ (5 × 4) y θⱼ (5,), en el orden j = 1..5: primero las 4 componentes de kⱼ y luego θⱼ."""
    ks = np.empty((C2_MODES, 4))
    th = np.empty(C2_MODES)
    for j in range(C2_MODES):
        ks[j] = aux.uniform(-C2_KMAX, C2_KMAX, size=4)
        th[j] = aux.uniform(0.0, 2.0 * math.pi)
    return ks, th


def c2_phi(x: Float, ks: Float, th: Float) -> Float:
    """φ(x) = Σⱼ 0.07·cos(2π kⱼ·x + θⱼ), x = (t, x, y, z)."""
    out: Float = np.sum(C2_AMP * np.cos(2.0 * math.pi * (x @ ks.T) + th[None, :]), axis=1)
    return out


def c2_points(n: int, r: np.random.Generator, ks: Float, th: Float) -> Float:
    """Candidatos uniformes del diamante estándar aceptados con probabilidad exp(φ − 0.35); los n primeros."""
    pts: list[Float] = []
    total = 0
    while total < n:
        cand = sprinkle_diamond(4 * n, 4, r)
        keep = cand[r.random(cand.shape[0]) < np.exp(c2_phi(cand, ks, th) - C2_SHIFT)]
        pts.append(keep)
        total += keep.shape[0]
    out: Float = np.concatenate(pts)[:n]
    return out


def c2_conforme(n: int, r: np.random.Generator, aux: np.random.Generator) -> tuple[Bool, Float | None, dict[str, Any]]:
    """C2: Minkowski con densidad exp(φ); los modos salen de `aux` antes de generar los puntos."""
    ks, th = c2_modes(aux)
    x = c2_points(n, r, ks, th)
    return causal_blocks(x, "mink"), x, {"c2_k": ks.tolist(), "c2_theta": th.tolist()}


# Conos: within(|x|, |y|, |z|, s) ⇔ ‖(x, y, z)‖_cono ≤ s (s ≥ 0). Región {|t| ≤ T − ‖x‖} y relación
# Δt > 0 y Δt ≥ ‖Δx‖_cono (convenciones del archivo sellado).


def _within_cilindro(ax: Float, ay: Float, az: Float, s: Float) -> Bool:
    out: Bool = np.maximum(np.sqrt(ax * ax + ay * ay), az) <= s
    return out


C4_S = 0.79


def _within_esfera_cubo(ax: Float, ay: Float, az: Float, s: Float) -> Bool:
    out: Bool = np.maximum(np.sqrt(ax * ax + ay * ay + az * az), np.maximum(np.maximum(ax, ay), az) / C4_S) <= s
    return out


C6_P = (1.0, 8.0, 32.0)


def _within_p_por_eje(ax: Float, ay: Float, az: Float, s: Float) -> Bool:
    """‖x‖ ≤ s ⇔ Σᵢ |xᵢ/s|^{pᵢ} ≤ 1 (gauge de Minkowski de la bola de C6), sin bisección.

    Cada término es ≤ 1 dentro de la bola, así que primero se exige |xᵢ| ≤ s (evita desbordes con p = 32) y
    solo ahí se divide; s = 0 (frontera de medida nula) se resuelve como x = 0, sin 0/0.
    """
    ax, ay, az, s = np.broadcast_arrays(ax, ay, az, s)
    out = np.zeros(ax.shape, dtype=np.bool_)
    pos = s > 0
    cand = pos & (ax <= s) & (ay <= s) & (az <= s)
    sc = s[cand]
    acc = (ax[cand] / sc) ** C6_P[0] + (ay[cand] / sc) ** C6_P[1] + (az[cand] / sc) ** C6_P[2]
    out[cand] = acc <= 1.0
    zero = ~pos
    out[zero] = (ax[zero] == 0) & (ay[zero] == 0) & (az[zero] == 0)
    return out


def c6_norm_bisect(x: Float, rtol: float = 1e-12) -> Float:
    """Norma de C6 por bisección (definición literal): λ > 0 con Σᵢ |xᵢ/λ|^{pᵢ} = 1; 0 en x = 0."""
    a = np.abs(np.atleast_2d(np.asarray(x, dtype=np.float64)))
    out = np.zeros(a.shape[0])
    for i, row in enumerate(a):
        inf = float(row.max())
        if inf == 0.0:
            continue
        lo, hi = inf, 3.0 * inf  # g(‖x‖∞) ≥ 1 y g(3‖x‖∞) ≤ 1/3 + 3⁻⁸ + 3⁻³² < 1

        def g(lam: float) -> float:
            return sum(float((row[k] / lam) ** C6_P[k]) for k in range(3))

        while hi - lo > rtol * hi:
            mid = 0.5 * (lo + hi)
            if g(mid) > 1.0:
                lo = mid
            else:
                hi = mid
        out[i] = 0.5 * (lo + hi)
    return out


def _within_anidado(ax: Float, ay: Float, az: Float, s: Float) -> Bool:
    q = ax * ax + ay * ay
    out: Bool = (q * q + az**4) ** 0.25 <= s
    return out


def cone_blocks(coords: Float, within: Within, block: int = BLOCK_ROWS) -> Bool:
    """Relación Δt > 0 y Δt ≥ ‖Δx‖_cono por bloques de filas."""
    n = coords.shape[0]
    t = coords[:, 0]
    c = np.empty((n, n), dtype=np.bool_)
    for lo in range(0, n, block):
        hi = min(n, lo + block)
        dt = t[None, :] - t[lo:hi, None]
        d = [np.abs(coords[None, :, k] - coords[lo:hi, None, k]) for k in (1, 2, 3)]
        c[lo:hi] = (dt > 0) & within(d[0], d[1], d[2], dt)
    return c


def in_cone_region(x: Float, within: Within) -> Bool:
    """{|t| ≤ T − ‖x‖_cono}."""
    a = np.abs(x)
    out: Bool = (a[:, 0] <= T_HALF) & within(a[:, 1], a[:, 2], a[:, 3], T_HALF - a[:, 0])
    return out


def sprinkle_cone(n: int, within: Within, rng: np.random.Generator) -> Float:
    """n puntos uniformes en {|t| ≤ T − ‖x‖_cono} por rechazo en [−T, T]⁴ (como `f7.sprinkle_lp`)."""
    if n < 0:
        raise ValueError("n debe ser >= 0")
    out = np.empty((n, 4), dtype=np.float64)
    filled = 0
    while filled < n:
        batch = max(64, 32 * (n - filled))
        cand = rng.uniform(-T_HALF, T_HALF, size=(batch, 4))
        ok = cand[in_cone_region(cand, within)]
        take = min(len(ok), n - filled)
        out[filled : filled + take] = ok[:take]
        filled += take
    return out


def _cone_gen(within: Within) -> SealedGen:
    def gen(n: int, r: np.random.Generator, aux: np.random.Generator) -> tuple[Bool, Float | None, dict[str, Any]]:
        x = sprinkle_cone(n, within, r)
        return cone_blocks(x, within), x, {}

    return gen


C5_FRAC = Fraction(13, 100)


def link_mask(c: Bool) -> Bool:
    """Links: C[i, j] y ningún k con C[i, k] y C[k, j] (B == 0; B = −1 fuera de C)."""
    out: Bool = between_blocks(c) == 0
    return out


def transitive_step(c: Bool, block: int = 2048) -> Bool:
    """C ∪ (C ∘ C) por bloques (float32 exacto); igual a C si y solo si el cierre transitivo de C es C."""
    n = c.shape[0]
    cf = c.astype(np.float32)
    out = c.copy()
    for lo in range(0, n, block):
        hi = min(n, lo + block)
        out[lo:hi] |= (cf[lo:hi] @ cf) > 0
    return out


def c5_delete_links(c: Bool, aux: np.random.Generator) -> tuple[Bool, dict[str, Any]]:
    """C' = C sin exactamente round(0.13·|L|) links elegidos con `aux.choice` sobre L en orden lexicográfico."""
    lmask = link_mask(c)
    ii, jj = np.nonzero(lmask)  # np.nonzero recorre por filas: orden lexicográfico de (i, j)
    n_links = int(ii.size)
    n_del = round(C5_FRAC * n_links)  # redondeo exacto (racional) al entero más cercano, empates a par
    sel = aux.choice(n_links, size=n_del, replace=False)
    cp = c.copy()
    cp[ii[sel], jj[sel]] = False
    lmask2 = link_mask(cp)
    n_links2 = int(np.count_nonzero(lmask2))
    new_links = int(np.count_nonzero(lmask2 & ~lmask))
    n_rel, n_rel2 = int(np.count_nonzero(c)), int(np.count_nonzero(cp))
    po = bool(is_partial_order(cp))
    closed = bool(np.array_equal(transitive_step(cp), cp))
    info = {"links_C": n_links, "borrados": int(n_del), "relaciones_C": n_rel, "relaciones_Cp": n_rel2,
            "links_Cp": n_links2, "links_nuevos": new_links,
            "frac_links_C": n_links / n_rel if n_rel else None, "frac_links_Cp": n_links2 / n_rel2 if n_rel2 else None,
            "orden_parcial": po, "cierre_igual": closed}
    if not (po and closed):
        raise RuntimeError(f"C5: verificación obligatoria fallida {info}")
    return cp, info


def c5_links_borrados(n: int, r: np.random.Generator, aux: np.random.Generator) -> tuple[Bool, Float | None, dict[str, Any]]:
    """C5: diamante estándar, C de Minkowski y 13% de sus links borrados; la prueba usa C'."""
    x = sprinkle_diamond(n, 4, r)
    cp, info = c5_delete_links(causal_blocks(x, "mink"), aux)
    return cp, x, {"c5": info}


# Orden = índice k del archivo sellado. (nombre, papel, generador)
SELLADAS: list[tuple[str, str, SealedGen]] = [
    ("C0_kasner_vacio", "positiva", c0_kasner),
    ("C1_cono_elipsoidal_rotado", "positiva", c1_elipsoidal),
    ("C2_minkowski_densidad_suave", "positiva", c2_conforme),
    ("C3_cono_cilindrico", "negativa", _cone_gen(_within_cilindro)),
    ("C4_cono_esfera_cubo", "negativa", _cone_gen(_within_esfera_cubo)),
    ("C5_links_borrados_13", "negativa", c5_links_borrados),
    ("C6_cono_p_por_eje", "negativa", _cone_gen(_within_p_por_eje)),
    ("C7_cono_anidado", "zona_gris", _cone_gen(_within_anidado)),
]
CONE_WITHIN: dict[int, Within] = {3: _within_cilindro, 4: _within_esfera_cubo, 6: _within_p_por_eje, 7: _within_anidado}

# Ya vistas (§4), con el generador de la exploración. El orden es el índice de la familia 3.
VISTAS: dict[str, str] = {
    "diamante_4d": "positiva",
    "de_sitter_4d": "positiva",
    "de_sitter_4d_fuerte": "positiva",
    "anti_de_sitter_4d": "positiva",
    "frw_radiacion_4d": "positiva",
    "frw_polvo_4d": "positiva",
    "cono_cubico_4d": "negativa",
    "cono_octaedrico_4d": "negativa",
    "cono_p6_4d": "negativa",
    "producto_4": "negativa",
    "percolacion_005": "negativa",
    "cono_p1_5_4d": "zona_gris",
    "cono_p3_4d": "zona_gris",
}


def cal_kind(g: int) -> str:
    if not 0 <= g < N_CAL:
        raise ValueError("grupo de calibración fuera de 0..14")
    return "diamante_plano" if g < 10 else ("de_sitter" if g < 13 else "frw_polvo")


def cal_points(g: int, n: int, r: np.random.Generator) -> Float:
    kind = cal_kind(g)
    if kind == "diamante_plano":
        return sprinkle_diamond(n, 4, r)
    if kind == "de_sitter":
        return de_sitter_eta_points(n, 4, r, -1.5, -0.5)
    return frw_power_points(n, 4, 8.0, r)  # a ∝ η²


# --- Partes ----------------------------------------------------------------------------------------


def ref_rows(reps: int = N_REF, ms: tuple[int, ...] = M_GRID, entropy: int = ENT_REF) -> list[dict[str, Any]]:
    """Como `f6_exploracion.ref_rows` con la entropía de la confirmación."""
    tips = np.array([[-T_HALF, 0, 0, 0], [T_HALF, 0, 0, 0]], dtype=np.float64)
    rows: list[dict[str, Any]] = []
    for m in ms:
        rs, profs = [], []
        for s in range(reps):
            x = np.vstack([tips[:1], sprinkle_diamond(m, 4, rng_for(entropy, FAM_REF, 0, m, s)), tips[1:]])
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


def nula_rows(n: int, reals: list[int] | None = None, entropy: int = ENT_ANALISIS) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for s in reals if reals is not None else range(N_NULA):
        r = rng_for(entropy, FAM_NULA, 0, n, s)
        c = causal_blocks(sprinkle_diamond(n, 4, r), "mink")
        rows.append({"realizacion": s, **measure_order(c, r)})
        print(f"nula {n} {s}", rows[-1]["t_medida_s"], "s", file=sys.stderr, flush=True)
    return rows


def _measure(label: str, s: int, gen: Callable[[], tuple[Bool, dict[str, Any]]], r: np.random.Generator) -> dict[str, Any]:
    t0 = time.time()
    c, extra = gen()
    t_c = time.time() - t0
    row = {"realizacion": s, "t_C_s": round(t_c, 2), **extra, **measure_order(c, r)}
    print(f"{label} {c.shape[0]} {s}", round(time.time() - t0, 1), "s", file=sys.stderr, flush=True)
    return row


def cal_rows(g: int, n: int, reals: list[int] | None = None, entropy: int = ENT_ANALISIS) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for s in reals if reals is not None else range(N_REAL):
        r = rng_for(entropy, FAM_CAL, g, n, s)
        rows.append(_measure(f"cal_{g}", s, lambda: (causal_blocks(cal_points(g, n, r), "mink"), {}), r))
    return rows


def vista_rows(name: str, n: int, reals: list[int] | None = None, entropy: int = ENT_VISTAS) -> list[dict[str, Any]]:
    gen = STRUCT[name][2]
    k = list(VISTAS).index(name)
    rows: list[dict[str, Any]] = []
    for s in reals if reals is not None else range(N_REAL):
        r = rng_for(entropy, FAM_VISTA, k, n, s)
        rows.append(_measure(name, s, lambda: (gen(n, r)[0], {}), r))
    return rows


def sellada_rows(k: int, n: int, reals: list[int] | None = None,
                 entropies: tuple[int, int] = (ENT_SELLADA_GEN, ENT_SELLADA_AUX)) -> list[dict[str, Any]]:
    name, _, gen = SELLADAS[k]
    rows: list[dict[str, Any]] = []
    for s in reals if reals is not None else range(N_REAL):
        r = rng_for(entropies[0], k, n, s)
        aux = rng_for(entropies[1], k, n, s)

        def build() -> tuple[Bool, dict[str, Any]]:
            c, _, extra = gen(n, r, aux)
            return c, extra

        rows.append(_measure(name, s, build, r))
    return rows


# --- Nombres de partes -------------------------------------------------------------------------------


def struct_bases() -> list[str]:
    return [f"cal_{g}" for g in range(N_CAL)] + [f"vista_{v}" for v in VISTAS] + [f"sellada_{k}" for k in range(len(SELLADAS))]


def parse_struct(part: str) -> tuple[str, list[int], list[int] | None] | None:
    """`<base>[_<N>[_r<s>]]` → (base, Ns, realizaciones o None)."""
    mt = re.fullmatch(r"(?P<base>(?:cal|vista|sellada)_.+?)(?:_(?P<n>4000|8000)(?:_r(?P<s>\d+))?)?", part)
    if mt is None or mt["base"] not in struct_bases():
        return None
    ns = [int(mt["n"])] if mt["n"] else list(NS)
    if mt["s"] is None:
        return mt["base"], ns, None
    s = int(mt["s"])
    return (mt["base"], ns, [s]) if s < N_REAL else None


def parse_nula(part: str) -> tuple[int, list[int] | None] | None:
    mt = re.fullmatch(r"nula_(?P<n>4000|8000)(?:_r(?P<s>\d+))?", part)
    if mt is None:
        return None
    if mt["s"] is None:
        return int(mt["n"]), None
    s = int(mt["s"])
    return (int(mt["n"]), [s]) if s < N_NULA else None


def struct_rows(base: str, ns: list[int], reals: list[int] | None) -> list[dict[str, Any]]:
    kind, _, rest = base.partition("_")
    rows: list[dict[str, Any]] = []
    for n in ns:
        if kind == "cal":
            rows += cal_rows(int(rest), n, reals)
        elif kind == "vista":
            rows += vista_rows(rest, n, reals)
        else:
            rows += sellada_rows(int(rest), n, reals)
    return rows


# Partición recomendada (tiempos medidos por realización con semillas de prueba, 4 núcleos: N = 8000 ≈ 10 a 18 s
# en casi todas; percolación ≈ 30 s, C5 ≈ 32 s y Kasner (K = 800) ≈ 78 s; N = 4000 ≈ 3 a 5 veces menos).
# Kasner va por N y, en N = 8000, por realización; C5 y percolación, por N; el resto, una parte por estructura.
SPLIT_POR_N = {"sellada_0", "sellada_5", "vista_percolacion_005"}
SPLIT_POR_REAL_8000 = {"sellada_0"}


def recommended_parts() -> list[str]:
    parts = ["ref", "nula_4000", "nula_8000"]
    for base in struct_bases():
        if base in SPLIT_POR_N:
            parts.append(f"{base}_4000")
            if base in SPLIT_POR_REAL_8000:
                parts += [f"{base}_8000_r{s}" for s in range(N_REAL)]
            else:
                parts.append(f"{base}_8000")
        else:
            parts.append(base)
    return parts + ["final"]


# --- Final -----------------------------------------------------------------------------------------


def load(part: str) -> dict[str, Any]:
    data: dict[str, Any] = json.loads((OUT / f"{part}.json").read_text())
    return data


def load_struct(base: str, n: int) -> list[dict[str, Any]]:
    """Filas de una estructura en un N: parte completa, parte de ese N o partes por realización."""
    if (OUT / f"{base}.json").exists():
        rows = [row for row in load(base)["rows"] if row["N"] == n]
    elif (OUT / f"{base}_{n}.json").exists():
        rows = load(f"{base}_{n}")["rows"]
    else:
        rows = [row for s in range(N_REAL) for row in load(f"{base}_{n}_r{s}")["rows"]]
    rows = sorted(rows, key=lambda row: int(row["realizacion"]))
    if [int(row["realizacion"]) for row in rows] != list(range(N_REAL)):
        raise RuntimeError(f"{base} N = {n}: realizaciones incompletas")
    return rows


def load_nula(n: int) -> list[dict[str, Any]]:
    if (OUT / f"nula_{n}.json").exists():
        rows: list[dict[str, Any]] = load(f"nula_{n}")["rows"]
    else:
        rows = [row for s in range(N_NULA) for row in load(f"nula_{n}_r{s}")["rows"]]
    if sorted(int(row["realizacion"]) for row in rows) != list(range(N_NULA)):
        raise RuntimeError(f"nula N = {n}: realizaciones incompletas")
    return rows


def fit_linear_boot(dr: dict[int, Float], rng: np.random.Generator, n_boot: int = N_BOOT) -> dict[str, Any]:
    """Δr = b + A·x sobre las celdas (m, N) con al menos la mitad de realizaciones medibles (como
    `f6_exploracion.fit_structure`, solo el modelo lineal); bootstrap de realizaciones completas dentro de cada N."""
    keep = {n: np.sum(np.isfinite(a), axis=0) >= a.shape[0] / 2 for n, a in dr.items()}
    full = {n: np.arange(a.shape[0]) for n, a in dr.items()}
    x, y = _cells(dr, keep, full)
    out: dict[str, Any] = {"puntos": [[float(a), float(b)] for a, b in zip(x, y)],
                           "x_distintos": int(np.unique(x).size)}
    est = fit_poly(x, y, 1)
    if est is None:
        out.update({"medible": False, "b": None, "A": None, "replicas_validas": 0})
        return out
    boots: list[list[float]] = []
    for _ in range(n_boot):
        idx = {n: np.asarray(rng.integers(0, a.shape[0], size=a.shape[0]), dtype=np.int64) for n, a in dr.items()}
        xb, yb = _cells(dr, keep, idx)
        f = fit_poly(xb, yb, 1)
        if f is not None:
            boots.append(f)
    bt = np.array(boots) if boots else np.full((1, 2), np.nan)
    out["medible"] = True
    for i, nm in enumerate(("b", "A")):
        out[nm] = {"valor": est[i], "se_boot": float(np.std(bt[:, i], ddof=1)) if bt.shape[0] > 1 else None,
                   "ic95": [float(np.percentile(bt[:, i], 2.5)), float(np.percentile(bt[:, i], 97.5))]}
    out["replicas_validas"] = len(boots)
    return out


PRED_A = {"de_sitter_4d": "+", "de_sitter_4d_fuerte": "+", "anti_de_sitter_4d": "-", "frw_radiacion_4d": "0",
          "C0_kasner_vacio": "0"}


def _pred_a(name: str, fit: dict[str, Any]) -> dict[str, Any] | None:
    pred = PRED_A.get(name)
    if pred is None or fit["A"] is None:
        return None if pred is None else {"prediccion": pred, "A": None, "concuerda": None}
    a, (lo, hi) = fit["A"]["valor"], fit["A"]["ic95"]
    excl = lo > 0 or hi < 0
    ok = (not excl) if pred == "0" else (excl and (a > 0) == (pred == "+"))
    return {"prediccion": pred, "A": a, "ic95": [lo, hi], "signo": "+" if a > 0 else "-", "ic_excluye_0": excl, "concuerda": ok}


def final(head: str) -> dict[str, Any]:
    ref = {row["m"]: row for row in load("ref")["rows"]}
    rplano = {m: float(ref[m]["r_plano"]) for m in M_GRID}
    pref = {m: np.array(ref[m]["perfil_plano"], dtype=np.float64) for m in M_GRID}
    # Rama 2: nula por (N, m) y max-T (informativo).
    nula: dict[str, Any] = {}
    p95: dict[tuple[int, int], tuple[float, float]] = {}
    nula_maxt: dict[int, list[float]] = {}
    for n in NS:
        rows = load_nula(n)
        dl = _dloc_rows(rows, pref)
        drm = _dr_matrix(rows, rplano)
        nula[str(n)] = {}
        for jm, m in enumerate(M_GRID):
            if len(dl[m]) < 2:
                nula[str(n)][str(m)] = {"n": len(dl[m])}
                continue
            f, mx = np.array([v[0] for v in dl[m]]), np.array([v[1] for v in dl[m]])
            p95[(n, m)] = (float(np.percentile(f, 95)), float(np.percentile(mx, 95)))
            nula[str(n)][str(m)] = {"n": len(dl[m]), "dloc_forma_p95": p95[(n, m)][0], "dloc_max_p95": p95[(n, m)][1],
                                    "delta_r": _stat(drm[:, jm].tolist())}
        nula_maxt[n] = [_max_t(row, pref, p95, n) for row in rows]
    nula_maxt_p95 = {n: float(np.nanpercentile(v, 95)) if np.any(np.isfinite(v)) else None for n, v in nula_maxt.items()}

    # (nombre, conjunto, papel, base de la parte, índice de bootstrap)
    items: list[tuple[str, str, str, str, int]] = [(f"cal_{g}_{cal_kind(g)}", "calibracion", cal_kind(g), f"cal_{g}", g)
                                                   for g in range(N_CAL)]
    items += [(v, "vista", role, f"vista_{v}", BOOT_VISTA + i) for i, (v, role) in enumerate(VISTAS.items())]
    items += [(nm, "sellada", role, f"sellada_{k}", BOOT_SELLADA + k) for k, (nm, role, _) in enumerate(SELLADAS)]
    res: dict[str, Any] = {}
    t_dr: list[list[Any]] = []
    t_dl: list[list[Any]] = []
    c5_rows: list[list[Any]] = []
    for name, conj, role, base, jb in items:
        by_n = {n: load_struct(base, n) for n in NS}
        dr = {n: _dr_matrix(rows, rplano) for n, rows in by_n.items()}
        info: dict[str, Any] = {"conjunto": conj, "papel": role, "por_N": {}}
        for n, rows in by_n.items():
            dl = _dloc_rows(rows, pref)
            per_m: dict[str, Any] = {}
            for jm, m in enumerate(M_GRID):
                cells = [r["por_m"][str(m)] for r in rows]
                st = _stat(dr[n][:, jm].tolist())
                e: dict[str, Any] = {"x": math.sqrt(m / n), "local_en": sum(c["local"] for c in cells),
                                     "medible_en": sum(c["medible"] for c in cells), "delta_r": st}
                if dl[m]:
                    f, mx = np.array([v[0] for v in dl[m]]), np.array([v[1] for v in dl[m]])
                    nu = p95.get((n, m))
                    e["dloc_forma_media"], e["dloc_max_media"] = float(f.mean()), float(mx.mean())
                    e["nula_p95"] = {"forma": nu[0], "maximo": nu[1]} if nu else None
                    t_dl.append([name, n, m, len(dl[m]), e["dloc_forma_media"], nu[0] if nu else "", e["dloc_max_media"], nu[1] if nu else ""])
                per_m[str(m)] = e
                t_dr.append([name, conj, n, m, e["x"], e["medible_en"], len(rows), st["media"], st["se"]])
            maxt = [_max_t(r, pref, p95, n) for r in rows]
            mm = [r["m_max"] for r in rows]
            info["por_N"][str(n)] = {"m_max": {"media": float(np.mean(mm)), "min": int(min(mm)), "max": int(max(mm))},
                                     "por_m": per_m,
                                     "max_T": {"mediana": float(np.nanmedian(maxt)) if np.any(np.isfinite(maxt)) else None,
                                               "frac_sobre_nula_p95": (float(np.mean(np.array(maxt) > nula_maxt_p95[n]))
                                                                       if nula_maxt_p95[n] is not None else None),
                                               "nula_p95": nula_maxt_p95[n]}}
            if base == "sellada_5":
                for r in rows:
                    c5 = r["c5"]
                    c5_rows.append([n, r["realizacion"], c5["links_C"], c5["borrados"], c5["links_Cp"], c5["links_nuevos"],
                                    c5["frac_links_C"], c5["frac_links_Cp"], c5["orden_parcial"], c5["cierre_igual"]])
        info["ajuste"] = fit_linear_boot(dr, rng_for(ENT_ANALISIS, FAM_BOOT, jb))
        res[name] = info

    cal = [res[f"cal_{g}_{cal_kind(g)}"]["ajuste"] for g in range(N_CAL)]
    cal_unmeasurable = [g for g in range(N_CAL) if cal[g]["b"] is None]
    flat_b = [cal[g]["b"]["valor"] for g in range(10) if cal[g]["b"] is not None]
    b_plano = float(np.median(flat_b)) if flat_b else float("nan")
    devs = [abs(c["b"]["valor"] - b_plano) for c in cal if c["b"] is not None]
    tau_b = TAU_FACTOR * max(devs) if devs else float("nan")
    calibracion_valida = not cal_unmeasurable and math.isfinite(tau_b)

    t_fit: list[list[Any]] = []
    votos: dict[str, Any] = {}
    for name, info in res.items():
        fit = info["ajuste"]
        if fit["b"] is None:
            dev, margin, et = None, None, False  # «no medible»: no se clasifica como espacio-tiempo (§1, A2)
        else:
            dev = abs(fit["b"]["valor"] - b_plano)
            margin = dev / tau_b if tau_b > 0 else None
            et = bool(dev <= tau_b)
        info["clasificacion"] = {"medible": fit["medible"], "dif_b": dev, "margen": margin, "espacio_tiempo": et}
        pa = _pred_a(name, fit)
        if pa is not None:
            info["prediccion_A"] = pa
        role = info["papel"]
        if info["conjunto"] != "calibracion" and role in ("positiva", "negativa"):
            ok = et if role == "positiva" else not et
            votos[name] = {"papel": role, "conjunto": info["conjunto"], "medible": fit["medible"], "espacio_tiempo": et, "acierto": ok}
        b, a = fit["b"], fit["A"]
        t_fit.append([name, info["conjunto"], role, len(fit["puntos"]), fit["x_distintos"],
                      b["valor"] if b else "", f"[{b['ic95'][0]:.6f}, {b['ic95'][1]:.6f}]" if b else "",
                      a["valor"] if a else "", f"[{a['ic95'][0]:.6f}, {a['ic95'][1]:.6f}]" if a else "",
                      dev if dev is not None else "", margin if margin is not None else "",
                      "no medible" if not fit["medible"] else ("espacio-tiempo" if et else "rechazada"),
                      ("acierto" if votos[name]["acierto"] else "fallo") if name in votos else "sin voto"])
    pasa = calibracion_valida and all(v["acierto"] for v in votos.values())
    # Informativo §5: de Sitter fuerte con A mayor que de Sitter.
    a_ds, a_dsf = res["de_sitter_4d"]["ajuste"]["A"], res["de_sitter_4d_fuerte"]["ajuste"]["A"]
    c5_all = [bool(r[8]) and bool(r[9]) for r in c5_rows]

    OUT.mkdir(parents=True, exist_ok=True)
    _write_csv("tabla_delta_r.csv", ["estructura", "conjunto", "N", "m", "x", "medible_en", "realizaciones", "delta_r", "se"], t_dr)
    _write_csv("tabla_dloc.csv", ["estructura", "N", "m", "n", "dloc_forma", "nula_p95_forma", "dloc_max", "nula_p95_max"], t_dl)
    _write_csv("tabla_ajustes.csv", ["estructura", "conjunto", "papel", "celdas", "x_distintos", "b", "b_ic95", "A", "A_ic95",
                                     "abs_b_menos_b_plano", "margen_sobre_tau_b", "clasificacion", "voto"], t_fit)
    _write_csv("tabla_c5.csv", ["N", "realizacion", "links_C", "borrados", "links_Cp", "links_nuevos", "frac_links_C",
                                "frac_links_Cp", "orden_parcial", "cierre_igual"], c5_rows)
    return {"commit": head, "estado": "CONFIRMATORIA (R12-1), intento 1 de 2", "experiment_id": EXPERIMENT_ID,
            "sellado_sha256": SELLADO_SHA256,
            "referencia": {str(m): {"r_plano": rplano[m], "reps": ref[m]["reps"]} for m in M_GRID},
            "calibracion": {"b_plano": b_plano, "tau_b": tau_b, "regla": "b_plano = mediana de b en los 10 grupos planos; "
                            "tau_b = 1.5 × max_g |b_g − b_plano| sobre los 15 grupos", "grupos_no_medibles": cal_unmeasurable,
                            "valida": calibracion_valida},
            "veredicto": "PASA" if pasa else "FALLA",
            "votos": votos,
            "no_medibles": [nm for nm, info in res.items() if not info["ajuste"]["medible"]],
            "informativo": {"de_sitter_fuerte_A_mayor": (a_dsf["valor"] > a_ds["valor"]) if a_ds and a_dsf else None,
                            "c5_verificaciones_ok": all(c5_all) if c5_all else None,
                            "nula": nula, "nula_max_T_p95": {str(n): v for n, v in nula_maxt_p95.items()}},
            "estructuras": res}


def _write_csv(name: str, header: list[str], rows: list[list[Any]]) -> None:
    with (OUT / name).open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


# --- CLI -------------------------------------------------------------------------------------------


def sealed_hash_ok() -> bool:
    return hashlib.sha256(SELLADO.read_bytes()).hexdigest() == SELLADO_SHA256


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--part")
    ap.add_argument("--listar", action="store_true", help="imprime la partición recomendada y sale")
    args = ap.parse_args()
    if args.listar:
        print("\n".join(recommended_parts()))
        return
    part: str | None = args.part
    if part is None:
        sys.exit("falta --part")
    if FIRMA not in PREREG.read_text():
        sys.exit("F6 confirmatoria requiere la firma R12-1 en el preregistro")
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip():
        sys.exit("exige arbol git limpio")
    if (part.startswith("sellada_") or part == "final") and not sealed_hash_ok():
        sys.exit("el SHA-256 del archivo sellado no coincide: la confirmación queda anulada")
    t0 = time.time()
    rows: list[dict[str, Any]]
    st = parse_struct(part)
    nu = parse_nula(part)
    if part == "ref":
        rows = ref_rows()
    elif nu is not None:
        rows = nula_rows(*nu)
    elif st is not None:
        rows = struct_rows(*st)
    elif part == "final":
        out = final(head)
        (OUT / "summary.json").write_text(json.dumps(out, indent=2, default=str))
        print("final ok", out["veredicto"], round(time.time() - t0, 1), "s")
        return
    else:
        sys.exit("parte desconocida")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{part}.json").write_text(json.dumps({"commit": head, "tiempo_s": round(time.time() - t0, 1), "rows": rows}, default=str))
    print(part, "ok", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
