"""F6 intento 2 (docs/OMEGA_F6_INTENTO2_PREREGISTRO.md, R13-1; sello docs/SELLADO_F6_INTENTO2_CIEGA.md).

Partes: refperfil | nula_<N> | cal_<g> | vista_<nombre> | sellada_<k> | final.

Instrumento idéntico al del intento 1 (``tools/f6_confirmacion.py``, que importa de ``tools/f6_exploracion.py``;
ambos se importan, nada se copia), salvo la referencia: **Δr = media de r − 1/10** (exacta, sin parte ``ref``).

- ``refperfil``: perfil plano N_k(m) de referencia, solo para la rama 2 (informativa): 4000 diamantes planos 4D de
  m puntos más las puntas por m; entropía 20261086, familia 5 (spawn_key (322, 5, 0, m, s)).
- ``nula_<N>``: 30 diamantes planos 4D por N (entropía 20261084, familia 1), para la rama 2.
- ``cal_<g>``, g = 0..14: grupo de calibración, 30 realizaciones en N = 4000 y 30 en N = 8000 (entropía 20261084,
  familia 2). Grupos 0..9 diamantes planos, 10..12 de Sitter (η ∈ [−1.5, −0.5]), 13..14 FRW polvo.
- ``vista_<nombre>``: las 21 estructuras vistas (13 geometrías y conos de la exploración y C0 a C7 del intento 1),
  entropía 20261085, familia 3, índice = orden en ``VISTAS``.
- ``sellada_<k>``, k = 0..8: D0 a D8 del sello; puntos con la entropía 20261082 y azar auxiliar (D4, D5) con la
  20261083, ambas con spawn_key (322, k, N, s).
- ``final``: b y A con bootstrap (familia 4), K, b_plano, τ_b, τ_K, clasificación, veredicto y lo informativo.

Partición: toda parte de estructura y ``nula_<N>`` admiten ``_<N>`` y ``_<N>_r<a>-<b>`` (realizaciones a..b,
inclusive; ``_r<s>`` = una sola): mismas semillas, mismas filas. ``final`` une, por estructura y N, todas las partes
presentes y exige exactamente las realizaciones 0..29 sin repetir, y que todas las partes usadas tengan el mismo
commit. ``--listar`` imprime la partición recomendada (cada parte por debajo de 8 minutos con los tiempos medidos).

Semillas: SeedSequence(entropía, spawn_key=(322, familia, índice, N, realización)). Bootstrap: spawn_key
(322, 4, j), j = g en la calibración, 100 + índice en las vistas y 200 + k en las selladas (entropía 20261084).
La selección de intervalos de cada realización sigue con el generador de puntos (como en los intentos previos);
en las vistas C0 a C7, el azar auxiliar de C2 y C5 sale del mismo generador de la realización (un solo spawn_key).
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

from f6_confirmacion import (  # noqa: E402
    SELLADAS as SELLADAS_I1,
    VISTAS as VISTAS_I1,
    _measure,
    cal_kind,
    cal_points,
    fit_linear_boot,
    link_mask,
    transitive_step,
)
from f6_exploracion import (  # noqa: E402
    STRUCT,
    _dloc_rows,
    _dr_matrix,
    _max_t,
    _stat,
    between_blocks,
    causal_blocks,
    interval_obs,
    measure_order,
)
from f7 import frw_power_points  # noqa: E402

from omega.causality.kasner import kasner_causal_matrix, sprinkle_kasner_interval, validate_kasner_exponents  # noqa: E402
from omega.causality.orders import is_partial_order  # noqa: E402
from omega.causality.sprinkling import DIAMOND_HALF_HEIGHT, causal_matrix, sprinkle_diamond  # noqa: E402

EXPERIMENT_ID = 322
ENT_SELLADA_GEN = 20261082
ENT_SELLADA_AUX = 20261083
ENT_ANALISIS = 20261084  # calibración, nula y bootstrap
ENT_VISTAS = 20261085
ENT_REFPERFIL = 20261086
FAM_NULA, FAM_CAL, FAM_VISTA, FAM_BOOT, FAM_REFPERFIL = 1, 2, 3, 4, 5
BOOT_VISTA, BOOT_SELLADA = 100, 200
M_GRID = (10, 20, 40, 80, 160)
K_MS = (10, 20, 40, 80)  # pares (m, 4000) ↔ (2m, 8000)
NS = (4000, 8000)
N_REAL = 30
N_NULA = 30
N_REFPERFIL = 4000
N_BOOT = 1000
N_CAL = 15
TAU_FACTOR = 1.5
R_PLANO = Fraction(1, 10)
BLOCK_ROWS = 1024
T_HALF = DIAMOND_HALF_HEIGHT
PREREG = Path("docs") / "OMEGA_F6_INTENTO2_PREREGISTRO.md"
SELLADO = Path("docs") / "SELLADO_F6_INTENTO2_CIEGA.md"
SELLADO_SHA256 = "b64503ffb0d5f5772a0ab8c8dcf80e689371d75ec6f93a6d93f9912c930866cf"
FIRMA = "R13-1 | **Firmada"
OUT = Path("results") / "omega_f6_intento2"
Bool = npt.NDArray[np.bool_]
Float = npt.NDArray[np.float64]
Gauge = Callable[[Float, Float, Float], Float]
SealedGen = Callable[[int, np.random.Generator, np.random.Generator], tuple[Bool, Float | None, dict[str, Any]]]


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    """SeedSequence(entropía, spawn_key=(322, *key))."""
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


# --- D0 a D3 (positivas) ----------------------------------------------------------------------------

D0_P = (6.0 / 7.0, -2.0 / 7.0, 3.0 / 7.0)
D0_T_LO, D0_T_HI, D0_K = 0.70, 1.35, 800


def d0_kasner(n: int, gen: np.random.Generator, aux: np.random.Generator, k: int = D0_K) -> tuple[Bool, Float | None, dict[str, Any]]:
    """D0: Kasner de vacío p = (6/7, −2/7, 3/7), t ∈ [0.70, 1.35], K = 800 direcciones (solo `gen`)."""
    p = validate_kasner_exponents(D0_P)
    x = sprinkle_kasner_interval(n, p, D0_T_LO, D0_T_HI, gen, k=k)
    return kasner_causal_matrix(x, p, k=k), x, {}


D1_B: Float = np.array([0.25, -0.15, 0.10])


def d1_shear(x: Float) -> Float:
    """x' = x − b·t por fila: en (t, x') la relación de D1 es la de Minkowski."""
    out = x.copy()
    out[:, 1:] = x[:, 1:] - x[:, :1] * D1_B[None, :]
    return out


def d1_cizalla(n: int, gen: np.random.Generator, aux: np.random.Generator) -> tuple[Bool, Float | None, dict[str, Any]]:
    """D1: diamante estándar de coordenadas; Δt > 0 y ‖Δx − bΔt‖₂ ≤ Δt, como Minkowski tras la cizalla."""
    x = sprinkle_diamond(n, 4, gen)
    return causal_blocks(d1_shear(x), "mink"), x, {}


def d2_points(n: int, gen: np.random.Generator) -> Float:
    """Candidatos `sprinkle_diamond(4n)` aceptados con probabilidad exp(4·(t − 0.5)) (densidad ∝ e^{4η}); los n primeros."""
    pts: list[Float] = []
    total = 0
    while total < n:
        cand = sprinkle_diamond(4 * n, 4, gen)
        keep = cand[gen.random(cand.shape[0]) < np.exp(4.0 * (cand[:, 0] - 0.5))]
        pts.append(keep)
        total += keep.shape[0]
    out: Float = np.concatenate(pts)[:n]
    return out


def d2_frw_exp(n: int, gen: np.random.Generator, aux: np.random.Generator) -> tuple[Bool, Float | None, dict[str, Any]]:
    """D2: FRW con a(η) = e^η, η = t ∈ [−0.5, 0.5]; relación de Minkowski."""
    x = d2_points(n, gen)
    return causal_blocks(x, "mink"), x, {}


def d3_frw_rigida(n: int, gen: np.random.Generator, aux: np.random.Generator) -> tuple[Bool, Float | None, dict[str, Any]]:
    """D3: FRW de materia rígida a ∝ η^{1/2}, η = t + 1: `frw_power_points(n, 4, 2.0)`; relación de Minkowski."""
    x = frw_power_points(n, 4, 2.0, gen)
    return causal_blocks(x, "mink"), x, {}


# --- D4: relaciones agregadas y cierre transitivo -----------------------------------------------------

D4_Q = 0.10
D4_V = math.pi / 24.0


def d4_delta2(n: int) -> float:
    """δ² = ℓ_N², ℓ_N = (V/n)^{1/4}, V = π/24."""
    return float((D4_V / n) ** 0.5)


def d4_candidates(x: Float, c: Bool, delta2: float, block: int = BLOCK_ROWS) -> Bool:
    """Pares con t_j − t_i > 0, C[i, j] = False y ‖x_j − x_i‖₂² − (t_j − t_i)² ≤ δ²."""
    n = x.shape[0]
    t = x[:, 0]
    out = np.empty((n, n), dtype=np.bool_)
    for lo in range(0, n, block):
        hi = min(n, lo + block)
        dt = t[None, :] - t[lo:hi, None]
        dx2 = np.zeros_like(dt)
        for k in (1, 2, 3):
            diff = x[None, :, k] - x[lo:hi, None, k]
            dx2 += diff * diff
        out[lo:hi] = (dt > 0) & ~c[lo:hi] & (dx2 - dt * dt <= delta2)
    return out


def transitive_closure(s: Bool) -> tuple[Bool, int, bool]:
    """Itera S ← S ∪ (S∘S) (`transitive_step`) hasta que no cambie. Devuelve (C'', pasos, cierre_igual), donde
    cierre_igual es la comparación `transitive_step(C'') == C''` del último paso (siempre True al salir)."""
    steps = 0
    while True:
        nxt = transitive_step(s)
        steps += 1
        if np.array_equal(nxt, s):
            return s, steps, True
        s = nxt


def d4_build(x: Float, c: Bool, aux: np.random.Generator) -> tuple[Bool, dict[str, Any]]:
    """Pasos 2 a 5 de D4 y la verificación obligatoria (B2: si falla, error)."""
    n = x.shape[0]
    cand = d4_candidates(x, c, d4_delta2(n))
    ii, jj = np.nonzero(cand)  # por filas: orden lexicográfico de (i, j)
    del cand
    u = aux.random(ii.size)
    add = u < D4_Q
    s = c.copy()
    s[ii[add], jj[add]] = True
    n_c, n_s = int(np.count_nonzero(c)), int(np.count_nonzero(s))
    cpp, steps, closed = transitive_closure(s)
    del s
    n_cpp = int(np.count_nonzero(cpp))
    po = bool(is_partial_order(cpp))
    info = {"candidatos": int(ii.size), "agregados": int(np.count_nonzero(add)), "relaciones_C": n_c, "relaciones_S": n_s,
            "relaciones_Cpp": n_cpp, "frac_nueva": (n_cpp - n_c) / n_c if n_c else None,
            "links_C": int(np.count_nonzero(link_mask(c))), "links_Cpp": int(np.count_nonzero(link_mask(cpp))),
            "pasos_cierre": steps, "orden_parcial": po, "cierre_igual": closed}
    if not (po and closed):
        raise RuntimeError(f"D4: verificación obligatoria fallida (B2) {info}")
    return cpp, info


def d4_agregadas(n: int, gen: np.random.Generator, aux: np.random.Generator) -> tuple[Bool, Float | None, dict[str, Any]]:
    """D4: diamante estándar, C de Minkowski, candidatos a menos de ℓ_N del cono, 10% agregados con `aux` y cierre."""
    x = sprinkle_diamond(n, 4, gen)
    cpp, info = d4_build(x, causal_blocks(x, "mink"), aux)
    return cpp, x, {"d4": info}


# --- D5: gemelos -------------------------------------------------------------------------------------

D5_FRAC = Fraction(1, 10)


def _row_keys(c: Bool) -> npt.NDArray[np.uint8]:
    """Clave por elemento: (pasado, futuro) empaquetados."""
    out: npt.NDArray[np.uint8] = np.concatenate([np.packbits(c, axis=0).T, np.packbits(c, axis=1)], axis=1)
    return out


def twin_pairs(c: Bool) -> int:
    """Número de pares de elementos con el mismo pasado y el mismo futuro."""
    _, counts = np.unique(_row_keys(c), axis=0, return_counts=True)
    return int(sum(int(k) * (int(k) - 1) // 2 for k in counts))


def d5_expected_pairs(c: Bool, orig: npt.NDArray[np.int64]) -> tuple[int, int]:
    """(pares ya presentes en C, pares esperados tras duplicar `orig`). Si C tiene grupos de elementos con el mismo
    pasado y futuro (en la práctica, elementos aislados de un sprinkling finito) con tamaños g y d duplicados cada
    uno, el total es Σ C(g + d, 2); sin tales grupos vale n_dup."""
    _, inv, counts = np.unique(_row_keys(c), axis=0, return_inverse=True, return_counts=True)
    inv = np.asarray(inv).reshape(-1)
    dups = np.bincount(inv[orig], minlength=counts.size)
    tot = counts + dups
    return int(np.sum(counts * (counts - 1) // 2)), int(np.sum(tot * (tot - 1) // 2))


def d5_build(c: Bool, n: int, aux: np.random.Generator) -> tuple[Bool, npt.NDArray[np.int64], npt.NDArray[np.int64], dict[str, Any]]:
    """Pasos 3 a 5 de D5 sobre C (n_base × n_base) y la verificación obligatoria (B2: si falla, error).
    Devuelve (C', φ, perm, info).

    Verificación: orden parcial y cierre (sello) y, como coherencia de la implementación, que cada uno de los n_dup
    gemelos tenga exactamente el pasado y el futuro de su original en C' y que el número de pares con el mismo
    pasado y futuro en C' coincida con el esperado a partir de C (n_dup más los pares que ya traen los elementos
    aislados del sprinkling; ver `d5_expected_pairs`)."""
    n_base = c.shape[0]
    n_dup = n - n_base
    orig = np.asarray(aux.choice(n_base, size=n_dup, replace=False), dtype=np.int64)
    perm = np.asarray(aux.permutation(n), dtype=np.int64)
    phi = np.concatenate([np.arange(n_base), orig]).astype(np.int64)
    cp = c[np.ix_(phi, phi)][np.ix_(perm, perm)]
    po = bool(is_partial_order(cp))
    closed = bool(np.array_equal(transitive_step(cp), cp))
    pos = np.empty(n, dtype=np.int64)
    pos[perm] = np.arange(n)  # etiqueta en C' del elemento u de C₀
    tw, og = pos[n_base + np.arange(n_dup)], pos[orig]
    built = int(np.count_nonzero(np.all(cp[tw] == cp[og], axis=1) & np.all(cp[:, tw] == cp[:, og], axis=0) & ~cp[tw, og] & ~cp[og, tw]))
    pairs = twin_pairs(cp)
    prev, expected = d5_expected_pairs(c, orig)
    info = {"n_dup": n_dup, "n_base": n_base, "gemelos_construidos": built, "pares_gemelos": pairs,
            "pares_previos_en_C": prev, "pares_esperados": expected,
            "aislados_C": int(np.count_nonzero(~c.any(axis=0) & ~c.any(axis=1))),
            "relaciones_C": int(np.count_nonzero(c)), "relaciones_Cp": int(np.count_nonzero(cp)),
            "orden_parcial": po, "cierre_igual": closed}
    if not (po and closed and built == n_dup and pairs == expected):
        raise RuntimeError(f"D5: verificación obligatoria fallida (B2) {info}")
    return cp, phi, perm, info


def d5_gemelos(n: int, gen: np.random.Generator, aux: np.random.Generator) -> tuple[Bool, Float | None, dict[str, Any]]:
    """D5: n_dup = round(0.10·n) gemelos de originales distintos, etiquetas permutadas con `aux`."""
    n_dup = round(D5_FRAC * n)
    x = sprinkle_diamond(n - n_dup, 4, gen)
    cp, phi, perm, info = d5_build(causal_blocks(x, "mink"), n, aux)
    return cp, x[phi][perm], {"d5": info}


# --- D6 a D8: conos con gauge F -------------------------------------------------------------------------
# Relación: Δt > 0 y Δt ≥ F(Δx) (Δx con signo). Región: {t + T ≥ F(x) y T − t ≥ F(−x)}.

D6_N: Float = np.array([2.0, -1.0, 2.0]) / 3.0
D6_S = 0.80


def d6_gauge(dx: Float, dy: Float, dz: Float) -> Float:
    """F(x) = max(‖x‖₂, (x·n)/s), n = (2, −1, 2)/3, s = 0.80."""
    out: Float = np.maximum(np.sqrt(dx * dx + dy * dy + dz * dz), (dx * D6_N[0] + dy * D6_N[1] + dz * D6_N[2]) / D6_S)
    return out


D7_L = 1.45
D7_SQ = math.sqrt(D7_L * D7_L - 1.0)


def d7_gauge(dx: Float, dy: Float, dz: Float) -> Float:
    """Huso: r si θ ≥ φ*; r(cos θ + sin θ·√(L² − 1))/L si θ < φ*, con r cos θ = |z|, r sin θ = √(x² + y²) y
    θ < φ* = arccos(1/L) ⇔ |z|·L > r (forma algebraica de la fórmula cerrada del sello)."""
    rho = np.sqrt(dx * dx + dy * dy)
    az = np.abs(dz)
    r = np.sqrt(rho * rho + az * az)
    out: Float = np.where(az * D7_L > r, (az + rho * D7_SQ) / D7_L, r)
    return out


def d7_gauge_trig(x: Float) -> Float:
    """Fórmula cerrada del sello, literal (con θ y φ*), para filas (n, 3); solo para pruebas."""
    r = np.linalg.norm(x, axis=1)
    out = np.zeros(x.shape[0])
    nz = r > 0
    theta = np.arccos(np.clip(np.abs(x[nz, 2]) / r[nz], 0.0, 1.0))
    phis = math.acos(1.0 / D7_L)
    val = np.where(theta >= phis, r[nz], r[nz] * (np.cos(theta) + np.sin(theta) * D7_SQ) / D7_L)
    out[nz] = val
    return out


def d8_gauge(dx: Float, dy: Float, dz: Float) -> Float:
    """F(x) = 0.68‖x‖₂ + 0.32‖x‖₁."""
    out: Float = 0.68 * np.sqrt(dx * dx + dy * dy + dz * dz) + 0.32 * (np.abs(dx) + np.abs(dy) + np.abs(dz))
    return out


def gauge_blocks(coords: Float, gauge: Gauge, block: int = BLOCK_ROWS) -> Bool:
    """C[i, j] ⇔ Δt > 0 y Δt ≥ F(x_j − x_i), por bloques de filas."""
    n = coords.shape[0]
    t = coords[:, 0]
    c = np.empty((n, n), dtype=np.bool_)
    for lo in range(0, n, block):
        hi = min(n, lo + block)
        dt = t[None, :] - t[lo:hi, None]
        d = [coords[None, :, k] - coords[lo:hi, None, k] for k in (1, 2, 3)]
        c[lo:hi] = (dt > 0) & (dt >= gauge(d[0], d[1], d[2]))
    return c


def in_gauge_region(x: Float, gauge: Gauge) -> Bool:
    """{t + T ≥ F(x) y T − t ≥ F(−x)}."""
    out: Bool = (x[:, 0] + T_HALF >= gauge(x[:, 1], x[:, 2], x[:, 3])) & (T_HALF - x[:, 0] >= gauge(-x[:, 1], -x[:, 2], -x[:, 3]))
    return out


BOX_STD: tuple[Float, Float] = (np.full(4, -T_HALF), np.full(4, T_HALF))
BOX_D7: tuple[Float, Float] = (np.array([-T_HALF, -T_HALF, -T_HALF, -D7_L * T_HALF]),
                               np.array([T_HALF, T_HALF, T_HALF, D7_L * T_HALF]))


def sprinkle_gauge(n: int, gauge: Gauge, box: tuple[Float, Float], rng: np.random.Generator) -> Float:
    """n puntos uniformes en la región por rechazo en la caja (como `sprinkle_cone`: lotes de max(64, 32·faltantes),
    se aceptan en orden y se toman los n primeros)."""
    if n < 0:
        raise ValueError("n debe ser >= 0")
    out = np.empty((n, 4), dtype=np.float64)
    filled = 0
    while filled < n:
        batch = max(64, 32 * (n - filled))
        cand = rng.uniform(box[0], box[1], size=(batch, 4))
        ok = cand[in_gauge_region(cand, gauge)]
        take = min(len(ok), n - filled)
        out[filled : filled + take] = ok[:take]
        filled += take
    return out


def _gauge_gen(gauge: Gauge, box: tuple[Float, Float]) -> SealedGen:
    def gen(n: int, r: np.random.Generator, aux: np.random.Generator) -> tuple[Bool, Float | None, dict[str, Any]]:
        x = sprinkle_gauge(n, gauge, box, r)
        return gauge_blocks(x, gauge), x, {}

    return gen


# Orden = índice k del sello. (nombre, papel, generador)
SELLADAS: list[tuple[str, str, SealedGen]] = [
    ("D0_kasner_u2", "positiva", d0_kasner),
    ("D1_minkowski_cizalla", "positiva", d1_cizalla),
    ("D2_frw_exponencial", "positiva", d2_frw_exp),
    ("D3_frw_rigida", "positiva", d3_frw_rigida),
    ("D4_relaciones_agregadas", "negativa", d4_agregadas),
    ("D5_gemelos_10", "negativa", d5_gemelos),
    ("D6_cono_tapa_plana", "negativa", _gauge_gen(d6_gauge, BOX_STD)),
    ("D7_cono_huso", "negativa", _gauge_gen(d7_gauge, BOX_D7)),
    ("D8_cono_suma_normas", "negativa", _gauge_gen(d8_gauge, BOX_STD)),
]
GAUGES: dict[int, Gauge] = {6: d6_gauge, 7: d7_gauge, 8: d8_gauge}

# Vistas (§5): 13 de la exploración y C0 a C7 del intento 1. El orden es el índice de la familia 3.
VISTAS: dict[str, str] = {**VISTAS_I1, **{nm: role for nm, role, _ in SELLADAS_I1}}
_C_GEN: dict[str, SealedGen] = {nm: g for nm, _, g in SELLADAS_I1}


def vista_build(name: str, n: int, r: np.random.Generator) -> tuple[Bool, dict[str, Any]]:
    """Matriz de una vista. C0 a C7 reciben el mismo generador como `r` y como `aux` (consumo secuencial)."""
    if name in _C_GEN:
        c, _, extra = _C_GEN[name](n, r, r)
        return c, extra
    return STRUCT[name][2](n, r)[0], {}


# --- Partes ----------------------------------------------------------------------------------------


def refperfil_rows(reps: int = N_REFPERFIL, ms: tuple[int, ...] = M_GRID, entropy: int = ENT_REFPERFIL) -> list[dict[str, Any]]:
    """Perfil plano por m (m puntos del diamante plano 4D más las puntas), spawn_key (322, 5, 0, m, s)."""
    tips = np.array([[-T_HALF, 0, 0, 0], [T_HALF, 0, 0, 0]], dtype=np.float64)
    rows: list[dict[str, Any]] = []
    for m in ms:
        rs, profs = [], []
        for s in range(reps):
            x = np.vstack([tips[:1], sprinkle_diamond(m, 4, rng_for(entropy, FAM_REFPERFIL, 0, m, s)), tips[1:]])
            c = causal_matrix(x)
            b = between_blocks(c)
            if np.count_nonzero(c[0] & c[:, m + 1]) != m:
                raise RuntimeError("las puntas no contienen los m puntos")
            rr, pr = interval_obs(c, b, 0, m + 1)
            rs.append(rr)
            profs.append(pr)
        ra, pa = np.array(rs), np.array(profs, dtype=np.float64)
        rows.append({"m": m, "reps": reps, "r_media": float(ra.mean()), "r_sd": float(ra.std(ddof=1)),
                     "perfil_plano": pa.mean(axis=0).tolist(), "perfil_sd": pa.std(axis=0, ddof=1).tolist()})
    return rows


def nula_rows(n: int, reals: list[int] | None = None, entropy: int = ENT_ANALISIS) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for s in reals if reals is not None else range(N_NULA):
        r = rng_for(entropy, FAM_NULA, 0, n, s)
        rows.append(_measure("nula", s, lambda: (causal_blocks(sprinkle_diamond(n, 4, r), "mink"), {}), r))
    return rows


def cal_rows(g: int, n: int, reals: list[int] | None = None, entropy: int = ENT_ANALISIS) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for s in reals if reals is not None else range(N_REAL):
        r = rng_for(entropy, FAM_CAL, g, n, s)
        rows.append(_measure(f"cal_{g}", s, lambda: (causal_blocks(cal_points(g, n, r), "mink"), {}), r))
    return rows


def vista_rows(name: str, n: int, reals: list[int] | None = None, entropy: int = ENT_VISTAS) -> list[dict[str, Any]]:
    k = list(VISTAS).index(name)
    rows: list[dict[str, Any]] = []
    for s in reals if reals is not None else range(N_REAL):
        r = rng_for(entropy, FAM_VISTA, k, n, s)
        rows.append(_measure(name, s, lambda: vista_build(name, n, r), r))
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


def _reals(a: str | None, b: str | None, total: int) -> list[int] | None | bool:
    """None = todas; lista a..b inclusive; False = rango inválido."""
    if a is None:
        return None
    lo, hi = int(a), int(b) if b is not None else int(a)
    return list(range(lo, hi + 1)) if lo <= hi < total else False


def parse_struct(part: str) -> tuple[str, list[int], list[int] | None] | None:
    """`<base>[_<N>[_r<a>[-<b>]]]` → (base, Ns, realizaciones o None)."""
    mt = re.fullmatch(r"(?P<base>(?:cal|vista|sellada)_.+?)(?:_(?P<n>4000|8000)(?:_r(?P<a>\d+)(?:-(?P<b>\d+))?)?)?", part)
    if mt is None or mt["base"] not in struct_bases():
        return None
    reals = _reals(mt["a"], mt["b"], N_REAL)
    if reals is False:
        return None
    ns = [int(mt["n"])] if mt["n"] else list(NS)
    return mt["base"], ns, reals if isinstance(reals, list) else None


def parse_nula(part: str) -> tuple[int, list[int] | None] | None:
    mt = re.fullmatch(r"nula_(?P<n>4000|8000)(?:_r(?P<a>\d+)(?:-(?P<b>\d+))?)?", part)
    if mt is None:
        return None
    reals = _reals(mt["a"], mt["b"], N_NULA)
    if reals is False:
        return None
    return int(mt["n"]), reals if isinstance(reals, list) else None


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


# Partición recomendada. Segundos por realización (generación + medida), N = 4000 / 8000. Medidos con
# SeedSequence(840001) en una máquina de 4 núcleos (D0 10.9/47.2, C0 11.4/48.1, D4 5.8/34.2, D5 2.7/17.1,
# C5 4.8/21.0, percolación 3.3/21.5, conos D6-D8 ≈ 2.1/11-13, Minkowski y FRW ≈ 1.7/8.3, nula 2.3/8.2,
# refperfil 10 s en total) y escalados ×1.6 para cubrir la máquina de la corrida del intento 1, que fue más lenta
# (C0 78 s, Minkowski 12.5 s en N = 8000). Sobre eso, `SAFETY` y el límite de 8 minutos por parte.
_STD, _CONE, _KASNER = (3.0, 14.0), (4.0, 21.0), (18.5, 78.0)
SECS: dict[str, tuple[float, float]] = {
    "nula": (3.7, 13.5),
    "sellada_0": _KASNER, "vista_C0_kasner_vacio": _KASNER,
    "sellada_4": (9.5, 55.0), "sellada_5": (4.5, 27.5),
    "vista_C5_links_borrados_13": (8.0, 34.0), "vista_percolacion_005": (5.5, 34.0),
    **{b: _CONE for b in ("sellada_6", "sellada_7", "sellada_8", "vista_cono_p6_4d", "vista_cono_p1_5_4d", "vista_cono_p3_4d",
                          "vista_C3_cono_cilindrico", "vista_C4_cono_esfera_cubo", "vista_C6_cono_p_por_eje", "vista_C7_cono_anidado")},
}
DEFAULT_SECS = _STD  # Minkowski, de Sitter, AdS, FRW, conos poliédricos, producto, C1, C2, D1 a D3, calibración
LIMIT_S = 8 * 60
SAFETY = 1.35  # margen sobre lo medido


def _chunks(total: int, per_real: float) -> list[tuple[int, int]]:
    size = max(1, min(total, int(LIMIT_S / (SAFETY * per_real))))
    k = math.ceil(total / size)
    size = math.ceil(total / k)
    return [(a, min(total, a + size) - 1) for a in range(0, total, size)]


def _split(base: str, per: tuple[float, float], total: int) -> list[str]:
    if SAFETY * total * sum(per) <= LIMIT_S:
        return [base]
    parts: list[str] = []
    for n, sec in zip(NS, per):
        ch = _chunks(total, sec)
        parts += [f"{base}_{n}"] if len(ch) == 1 else [f"{base}_{n}_r{a}-{b}" for a, b in ch]
    return parts


def recommended_parts() -> list[str]:
    parts = ["refperfil"]
    for n, sec in zip(NS, SECS.get("nula", DEFAULT_SECS)):
        ch = _chunks(N_NULA, sec)
        parts += [f"nula_{n}"] if len(ch) == 1 else [f"nula_{n}_r{a}-{b}" for a, b in ch]
    for base in struct_bases():
        parts += _split(base, SECS.get(base, DEFAULT_SECS), N_REAL)
    return parts + ["final"]


def estimated_seconds() -> dict[str, float]:
    """Estimación (sin margen) por parte recomendada, salvo refperfil y final."""
    out: dict[str, float] = {}
    for p in recommended_parts():
        st, nu = parse_struct(p), parse_nula(p)
        if st is not None:
            per = SECS.get(st[0], DEFAULT_SECS)
            nreal = len(st[2]) if st[2] is not None else N_REAL
            out[p] = sum(per[NS.index(n)] * nreal for n in st[1])
        elif nu is not None:
            per = SECS.get("nula", DEFAULT_SECS)
            out[p] = per[NS.index(nu[0])] * (len(nu[1]) if nu[1] is not None else N_NULA)
    return out


# --- Final -----------------------------------------------------------------------------------------


def _part_files() -> list[tuple[str, dict[str, Any]]]:
    return [(f.stem, json.loads(f.read_text())) for f in sorted(OUT.glob("*.json")) if f.stem != "summary"]


def collect_rows(files: list[tuple[str, dict[str, Any]]], base: str, n: int, total: int) -> tuple[list[dict[str, Any]], list[str]]:
    """Filas de `base` (estructura o «nula») en N, de todas las partes presentes; exige 0..total−1 sin repetir."""
    rows: list[dict[str, Any]] = []
    used: list[str] = []
    for part, data in files:
        if base == "nula":
            nu = parse_nula(part)
            match = nu is not None and nu[0] == n
        else:
            st = parse_struct(part)
            match = st is not None and st[0] == base and n in st[1]
        if not match:
            continue
        sel = [row for row in data["rows"] if int(row["N"]) == n]
        if sel:
            rows += sel
            used.append(part)
    rows.sort(key=lambda row: int(row["realizacion"]))
    if [int(row["realizacion"]) for row in rows] != list(range(total)):
        raise RuntimeError(f"{base} N = {n}: realizaciones incompletas o repetidas")
    return rows, used


def cell_stats(dr: dict[int, Float]) -> dict[int, dict[int, tuple[float, float, int]]]:
    """Celdas (m, N) presentes (al menos la mitad de realizaciones medibles): (media de Δr, se, medibles),
    se = desviación estándar entre realizaciones (ddof = 1) / √medibles."""
    out: dict[int, dict[int, tuple[float, float, int]]] = {}
    for n, arr in dr.items():
        out[n] = {}
        for jm, m in enumerate(M_GRID):
            col = arr[:, jm]
            col = col[np.isfinite(col)]
            if col.size >= arr.shape[0] / 2 and col.size > 0:
                se = float(col.std(ddof=1) / math.sqrt(col.size)) if col.size > 1 else float("nan")
                out[n][m] = (float(col.mean()), se, int(col.size))
    return out


def collapse_k(dr: dict[int, Float]) -> dict[str, Any]:
    """K = Σⱼ (dⱼ/seⱼ)², dⱼ = Δr(2m, 8000) − Δr(m, 4000), m ∈ {10, 20, 40, 80}, seⱼ = √(se₁² + se₂²); solo
    pares con las dos celdas presentes. No medible con menos de dos pares. seⱼ = 0: (d/se)² = 0 si dⱼ = 0 e ∞ si no."""
    cs = cell_stats(dr)
    pares: list[dict[str, Any]] = []
    for m in K_MS:
        a, b = cs.get(4000, {}).get(m), cs.get(8000, {}).get(2 * m)
        if a is None or b is None:
            continue
        d = b[0] - a[0]
        se = math.sqrt(a[1] ** 2 + b[1] ** 2)
        if se > 0:
            z = d / se
        else:
            z = 0.0 if d == 0 else math.inf
        pares.append({"m": m, "d": d, "se": se, "z": z, "dr_4000": a[0], "se_4000": a[1], "dr_8000": b[0], "se_8000": b[1]})
    k = float(sum(p["z"] ** 2 for p in pares))
    return {"pares": pares, "n_pares": len(pares), "medible": len(pares) >= 2, "K": k if len(pares) >= 2 else None}


def final(head: str) -> dict[str, Any]:
    files = _part_files()
    commits: dict[str, str] = {}
    rp = {m: float(R_PLANO) for m in M_GRID}
    ref = {row["m"]: row for row in dict(files)["refperfil"]["rows"]}
    commits["refperfil"] = dict(files)["refperfil"]["commit"]
    pref = {m: np.array(ref[m]["perfil_plano"], dtype=np.float64) for m in M_GRID}
    # Rama 2 (informativa): nula por (N, m) y max-T.
    nula: dict[str, Any] = {}
    p95: dict[tuple[int, int], tuple[float, float]] = {}
    nula_maxt: dict[int, list[float]] = {}
    for n in NS:
        rows, used = collect_rows(files, "nula", n, N_NULA)
        commits.update({p: dict(files)[p]["commit"] for p in used})
        dl = _dloc_rows(rows, pref)
        drm = _dr_matrix(rows, rp)
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
    t_col: list[list[Any]] = []
    t_d4: list[list[Any]] = []
    t_d5: list[list[Any]] = []
    t_c5: list[list[Any]] = []
    for name, conj, role, base, jb in items:
        by_n: dict[int, list[dict[str, Any]]] = {}
        for n in NS:
            by_n[n], used = collect_rows(files, base, n, N_REAL)
            commits.update({p: dict(files)[p]["commit"] for p in used})
        dr = {n: _dr_matrix(rows, rp) for n, rows in by_n.items()}
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
            for r in rows:
                if "d4" in r:
                    d = r["d4"]
                    t_d4.append([n, r["realizacion"], d["candidatos"], d["agregados"], d["relaciones_C"], d["relaciones_S"],
                                 d["relaciones_Cpp"], d["frac_nueva"], d["links_C"], d["links_Cpp"], d["pasos_cierre"],
                                 d["orden_parcial"], d["cierre_igual"]])
                if "d5" in r:
                    d = r["d5"]
                    t_d5.append([n, r["realizacion"], d["n_dup"], d["gemelos_construidos"], d["pares_gemelos"],
                                 d["pares_previos_en_C"], d["pares_esperados"], d["aislados_C"], d["relaciones_C"],
                                 d["relaciones_Cp"], d["orden_parcial"], d["cierre_igual"]])
                if "c5" in r:
                    c5 = r["c5"]
                    t_c5.append([n, r["realizacion"], c5["links_C"], c5["borrados"], c5["links_Cp"], c5["links_nuevos"],
                                 c5["orden_parcial"], c5["cierre_igual"]])
        info["ajuste"] = fit_linear_boot(dr, rng_for(ENT_ANALISIS, FAM_BOOT, jb))
        info["colapso"] = collapse_k(dr)
        for p in info["colapso"]["pares"]:
            t_col.append([name, p["m"], p["dr_4000"], p["se_4000"], p["dr_8000"], p["se_8000"], p["d"], p["se"], p["z"]])
        res[name] = info

    if len(set(commits.values())) != 1:
        raise RuntimeError(f"las partes usadas no tienen el mismo commit: {sorted(set(commits.values()))}")

    cal = [res[f"cal_{g}_{cal_kind(g)}"] for g in range(N_CAL)]
    cal_unmeasurable = [g for g in range(N_CAL) if cal[g]["ajuste"]["b"] is None or not cal[g]["colapso"]["medible"]]
    flat_b = [cal[g]["ajuste"]["b"]["valor"] for g in range(10) if cal[g]["ajuste"]["b"] is not None]
    b_plano = float(np.median(flat_b)) if flat_b else float("nan")
    devs = [abs(c["ajuste"]["b"]["valor"] - b_plano) for c in cal if c["ajuste"]["b"] is not None]
    tau_b = TAU_FACTOR * max(devs) if devs else float("nan")
    ks = [c["colapso"]["K"] for c in cal if c["colapso"]["K"] is not None]
    tau_k = TAU_FACTOR * max(ks) if ks else float("nan")
    calibracion_valida = not cal_unmeasurable and math.isfinite(tau_b) and math.isfinite(tau_k)

    t_fit: list[list[Any]] = []
    votos: dict[str, Any] = {}
    for name, info in res.items():
        fit, col = info["ajuste"], info["colapso"]
        medible = bool(fit["medible"]) and bool(col["medible"])
        dev = abs(fit["b"]["valor"] - b_plano) if fit["b"] is not None else None
        m_b = dev / tau_b if dev is not None and tau_b > 0 else None
        m_k = col["K"] / tau_k if col["K"] is not None and tau_k > 0 else None
        c1 = dev is not None and dev <= tau_b
        c2 = col["K"] is not None and col["K"] <= tau_k
        et = bool(medible and c1 and c2)  # no medible = no espacio-tiempo (§3; B1)
        info["clasificacion"] = {"medible": medible, "intercepto_medible": bool(fit["medible"]), "colapso_medible": bool(col["medible"]),
                                 "dif_b": dev, "margen_b": m_b, "K": col["K"], "margen_K": m_k,
                                 "cond_intercepto": bool(c1), "cond_colapso": bool(c2), "espacio_tiempo": et}
        a = fit["A"]
        if a is not None:
            lo, hi = a["ic95"]
            info["pendiente_A"] = {"A": a["valor"], "ic95": [lo, hi], "signo": "+" if a["valor"] > 0 else "-",
                                   "ic_excluye_0": bool(lo > 0 or hi < 0)}
        role = info["papel"]
        if info["conjunto"] != "calibracion" and role in ("positiva", "negativa"):
            ok = et if role == "positiva" else not et
            votos[name] = {"papel": role, "conjunto": info["conjunto"], "medible": medible, "espacio_tiempo": et, "acierto": ok}
        b = fit["b"]
        t_fit.append([name, info["conjunto"], role, len(fit["puntos"]), fit["x_distintos"],
                      b["valor"] if b else "", f"[{b['ic95'][0]:.6f}, {b['ic95'][1]:.6f}]" if b else "",
                      a["valor"] if a else "", f"[{a['ic95'][0]:.6f}, {a['ic95'][1]:.6f}]" if a else "",
                      dev if dev is not None else "", m_b if m_b is not None else "",
                      col["n_pares"], col["K"] if col["K"] is not None else "", m_k if m_k is not None else "",
                      "no medible" if not medible else ("espacio-tiempo" if et else "rechazada"),
                      ("acierto" if votos[name]["acierto"] else "fallo") if name in votos else "sin voto"])
    pasa = calibracion_valida and all(v["acierto"] for v in votos.values())
    d4_ok = [bool(r[11]) and bool(r[12]) for r in t_d4]
    d5_ok = [bool(r[10]) and bool(r[11]) and int(r[3]) == int(r[2]) and int(r[4]) == int(r[6]) for r in t_d5]
    d5_pairs_ok = [int(r[4]) == int(r[2]) for r in t_d5]

    OUT.mkdir(parents=True, exist_ok=True)
    _write_csv("tabla_delta_r.csv", ["estructura", "conjunto", "N", "m", "x", "medible_en", "realizaciones", "delta_r", "se"], t_dr)
    _write_csv("tabla_dloc.csv", ["estructura", "N", "m", "n", "dloc_forma", "nula_p95_forma", "dloc_max", "nula_p95_max"], t_dl)
    _write_csv("tabla_colapso.csv", ["estructura", "m", "delta_r_m_4000", "se_4000", "delta_r_2m_8000", "se_8000", "d", "se", "d_sobre_se"], t_col)
    _write_csv("tabla_ajustes.csv", ["estructura", "conjunto", "papel", "celdas", "x_distintos", "b", "b_ic95", "A", "A_ic95",
                                     "abs_b_menos_b_plano", "margen_sobre_tau_b", "pares_K", "K", "margen_sobre_tau_K",
                                     "clasificacion", "voto"], t_fit)
    _write_csv("tabla_d4.csv", ["N", "realizacion", "candidatos", "agregados", "relaciones_C", "relaciones_S", "relaciones_Cpp",
                                "frac_nueva", "links_C", "links_Cpp", "pasos_cierre", "orden_parcial", "cierre_igual"], t_d4)
    _write_csv("tabla_d5.csv", ["N", "realizacion", "n_dup", "gemelos_construidos", "pares_gemelos", "pares_previos_en_C",
                                "pares_esperados", "aislados_C", "relaciones_C", "relaciones_Cp", "orden_parcial", "cierre_igual"], t_d5)
    _write_csv("tabla_c5.csv", ["N", "realizacion", "links_C", "borrados", "links_Cp", "links_nuevos", "orden_parcial", "cierre_igual"], t_c5)
    _write_csv("tabla_partes.csv", ["parte", "commit"], [[p, c] for p, c in sorted(commits.items())])
    return {"commit_final": head, "commit_partes": next(iter(commits.values())), "partes": commits,
            "estado": "F6 intento 2 (R13-1), último intento", "experiment_id": EXPERIMENT_ID, "sellado_sha256": SELLADO_SHA256,
            "referencia_r": {"r_plano": "1/10", "valor": float(R_PLANO)},
            "calibracion": {"b_plano": b_plano, "tau_b": tau_b, "tau_K": tau_k,
                            "regla": "b_plano = mediana de b en los 10 grupos planos; tau_b = 1.5 × max_g |b_g − b_plano|; "
                                     "tau_K = 1.5 × max_g K_g (15 grupos)",
                            "grupos_no_medibles": cal_unmeasurable, "valida": calibracion_valida},
            "veredicto": "PASA" if pasa else "FALLA",
            "votos": votos,
            "no_medibles": [nm for nm, info in res.items() if not info["clasificacion"]["medible"]],
            "informativo": {"d4_verificaciones_ok": all(d4_ok) if d4_ok else None,
                            "d5_verificaciones_ok": all(d5_ok) if d5_ok else None,
                            "d5_pares_gemelos_igual_n_dup": all(d5_pairs_ok) if d5_pairs_ok else None,
                            "d5_nota": "pares_gemelos cuenta también los pares de elementos aislados del sprinkling (mismo "
                                       "pasado y futuro vacíos); pares_esperados los descuenta desde C",
                            "nula": nula, "nula_max_T_p95": {str(n): v for n, v in nula_maxt_p95.items()},
                            "refperfil": {str(m): {"reps": ref[m]["reps"], "r_media": ref[m]["r_media"]} for m in M_GRID}},
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
    st = parse_struct(part)
    nu = parse_nula(part)
    if part not in ("refperfil", "final") and st is None and nu is None:
        sys.exit("parte desconocida")
    if FIRMA not in PREREG.read_text():
        sys.exit("F6 intento 2 requiere la firma R13-1 en el preregistro")
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip():
        sys.exit("exige arbol git limpio")
    if (part.startswith("sellada_") or part == "final") and not sealed_hash_ok():
        sys.exit("el SHA-256 del archivo sellado no coincide: el intento queda anulado")
    t0 = time.time()
    rows: list[dict[str, Any]]
    if part == "refperfil":
        rows = refperfil_rows()
    elif nu is not None:
        rows = nula_rows(*nu)
    elif st is not None:
        rows = struct_rows(*st)
    else:
        out = final(head)
        (OUT / "summary.json").write_text(json.dumps(out, indent=2, default=str))
        print("final ok", out["veredicto"], round(time.time() - t0, 1), "s")
        return
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{part}.json").write_text(json.dumps({"commit": head, "tiempo_s": round(time.time() - t0, 1), "rows": rows}, default=str))
    print(part, "ok", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
