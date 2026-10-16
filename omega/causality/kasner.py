"""Relacion causal y sprinkling en Kasner de vacio (R11-4): curvatura no conformemente plana.

Metrica: ds^2 = -dt^2 + sum_i t^(2 p_i) dx_i^2, con sum p_i = 1 y sum p_i^2 = 1; sqrt(-g) = t.

Relacion causal. El espacio-tiempo es homogeneo en x, asi que el conjunto alcanzable desde
(t1, x1) hasta t2 es x1 + R(t1, t2), con R la integral de Aumann de los elipsoides de
velocidades E(t) = {v : sum_i a_i(t)^2 v_i^2 <= 1}, a_i = t^p_i. R es convexo y compacto y su
funcion soporte es h(u; t1, t2) = int_{t1}^{t2} sqrt(sum_i u_i^2 / a_i(t)^2) dt. Luego
y esta en J^+(x) si y solo si t2 >= t1 y u.dx <= h(u; t1, t2) para todo u.

Discretizacion. Se usa un conjunto fijo de K direcciones simetrico por +-u (K/2 pares, con
los tres ejes incluidos) y H_u(t) acumulada sobre los tiempos de los puntos, h = H_u(t2) - H_u(t1).
La condicion i < j se escribe como dominancia componente a componente
    A[j] <= A[i] y B[j] <= B[i],  A = P - H,  B = -P - H,  P = u.x,
mas t_j > t_i. Comparar los mismos numeros almacenados hace la relacion exactamente
transitiva (tambien en coma flotante) para cualquier K. El poliedro resultante contiene al
conjunto alcanzable verdadero: la relacion discreta solo puede anadir pares de borde.

Las coordenadas solo existen para construir la relacion; ningun estimador las recibe.
"""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt

__all__ = [
    "validate_kasner_exponents",
    "kasner_exponents",
    "kasner_directions",
    "support_table",
    "kasner_causal_matrix",
    "in_kasner_interval",
    "sprinkle_kasner_interval",
]

FloatArray = npt.NDArray[np.float64]
BoolMatrix = npt.NDArray[np.bool_]

_GL_NODES = 8  # nodos de Gauss-Legendre por subintervalo en s = ln t
_MAX_DS = 0.02  # paso maximo en ln t de la integracion de H_u


def validate_kasner_exponents(p: tuple[float, float, float] | FloatArray, tol: float = 1e-9) -> FloatArray:
    """Devuelve p como arreglo (3,) si cumple sum p = 1 y sum p^2 = 1 con tolerancia tol."""
    arr = np.asarray(p, dtype=np.float64)
    if arr.shape != (3,) or not np.all(np.isfinite(arr)):
        raise ValueError("p debe tener tres exponentes finitos")
    if abs(float(arr.sum()) - 1.0) > tol:
        raise ValueError(f"sum p = {arr.sum():.12g} != 1: no es Kasner de vacio")
    if abs(float(np.sum(arr * arr)) - 1.0) > tol:
        raise ValueError(f"sum p^2 = {np.sum(arr * arr):.12g} != 1: no es Kasner de vacio")
    return arr


def kasner_exponents(u: float) -> tuple[float, float, float]:
    """Parametrizacion estandar de Lifshitz-Khalatnikov: (-u, 1+u, u(1+u)) / (1+u+u^2).

    Cubre todas las soluciones de vacio (hasta permutacion) con u >= 1. No fija ningun valor
    para una prueba: el parametro lo elige quien llama.
    """
    if not math.isfinite(u):
        raise ValueError("u debe ser finito")
    den = 1.0 + u + u * u
    out = (-u / den, (1.0 + u) / den, u * (1.0 + u) / den)
    validate_kasner_exponents(out)
    return out


def kasner_directions(k: int) -> FloatArray:
    """K/2 direcciones unitarias de un hemisferio (Fibonacci) con los tres ejes al frente.

    El conjunto efectivo es {+-u}, de tamano K. K debe ser par y >= 6.
    """
    if k < 6 or k % 2:
        raise ValueError("K debe ser par y >= 6")
    m = k // 2
    axes = np.eye(3)
    rest = m - 3
    if rest == 0:
        return axes.copy()
    idx = np.arange(rest, dtype=np.float64)
    z = (idx + 0.5) / rest
    r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    phi = idx * math.pi * (3.0 - math.sqrt(5.0))
    fib = np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)
    return np.concatenate([axes, fib], axis=0)


def support_table(t: FloatArray, p: FloatArray, dirs: FloatArray) -> FloatArray:
    """H[i, k] = int_{t_ref}^{t_i} sqrt(sum_j u_kj^2 t'^(-2 p_j)) dt', con t_ref = min(t).

    Integra en s = ln t con Gauss-Legendre por subintervalos de ancho <= _MAX_DS, cuyos
    extremos incluyen todos los t_i; la suma acumulada da h = H(t2) - H(t1) exactamente aditiva.
    """
    t = np.asarray(t, dtype=np.float64)
    if t.ndim != 1:
        raise ValueError("t debe ser 1D")
    if t.size == 0:
        return np.zeros((0, dirs.shape[0]))
    if np.any(t <= 0) or not np.all(np.isfinite(t)):
        raise ValueError("Kasner requiere t > 0 finito")
    s_pts = np.log(t)
    s_lo, s_hi = float(s_pts.min()), float(s_pts.max())
    n_grid = max(2, int(math.ceil((s_hi - s_lo) / _MAX_DS)) + 1)
    knots = np.unique(np.concatenate([s_pts, np.linspace(s_lo, s_hi, n_grid)]))
    u2 = dirs * dirs  # (m, 3)
    if knots.size == 1:
        return np.zeros((t.size, dirs.shape[0]))
    x_gl, w_gl = np.polynomial.legendre.leggauss(_GL_NODES)
    left, right = knots[:-1], knots[1:]
    half = 0.5 * (right - left)
    mid = 0.5 * (right + left)
    m = dirs.shape[0]
    seg = np.empty((left.size, m))
    for c0 in range(0, left.size, 1024):  # por trozos: memoria O(1024 * g * m)
        c1 = min(left.size, c0 + 1024)
        s_nodes = mid[c0:c1, None] + half[c0:c1, None] * x_gl[None, :]  # (c, g)
        # integrando en s: e^s sqrt(sum_j u_j^2 e^{-2 p_j s})
        scale = np.exp(-2.0 * s_nodes[:, :, None] * p[None, None, :])  # (c, g, 3)
        integrand = np.sqrt(np.einsum("igj,kj->igk", scale, u2)) * np.exp(s_nodes)[:, :, None]
        seg[c0:c1] = np.einsum("igk,g->ik", integrand, w_gl) * half[c0:c1, None]
    cum = np.concatenate([np.zeros((1, dirs.shape[0])), np.cumsum(seg, axis=0)], axis=0)
    pos = np.searchsorted(knots, s_pts)
    out: FloatArray = cum[pos]
    return out


def _dominance_keys(coords: FloatArray, p: FloatArray, dirs: FloatArray) -> tuple[FloatArray, FloatArray]:
    proj = coords[:, 1:] @ dirs.T
    h = support_table(coords[:, 0], p, dirs)
    return proj - h, -proj - h


def _check_coords(coords: FloatArray) -> None:
    if coords.ndim != 2 or coords.shape[1] != 4:
        raise ValueError("coords debe ser (n, 4) con columna 0 = t")


def kasner_causal_matrix(
    coords: FloatArray,
    p: tuple[float, float, float] | FloatArray,
    k: int = 400,
    block: int = 256,
    chunk: int = 8,
) -> BoolMatrix:
    """C[i, j] = True si i precede estrictamente a j en Kasner(p), con K direcciones.

    Por bloques de `block` filas: se parte de los pares con t_j > t_i y se filtra direccion a
    direccion (de `chunk` en `chunk`, ejes primero), compactando los pares candidatos.
    """
    _check_coords(coords)
    pv = validate_kasner_exponents(p)
    dirs = kasner_directions(k)
    a_key, b_key = _dominance_keys(np.asarray(coords, dtype=np.float64), pv, dirs)
    t = coords[:, 0]
    n = coords.shape[0]
    m = dirs.shape[0]
    rel = np.zeros((n, n), dtype=bool)
    for r0 in range(0, n, block):
        r1 = min(n, r0 + block)
        later = t[None, :] > t[r0:r1, None]
        ii, jj = np.nonzero(later)
        ii = ii + r0
        for c0 in range(0, m, chunk):
            if ii.size == 0:
                break
            c1 = min(m, c0 + chunk)
            ok = np.all(a_key[jj, c0:c1] <= a_key[ii, c0:c1], axis=1)
            ok &= np.all(b_key[jj, c0:c1] <= b_key[ii, c0:c1], axis=1)
            ii, jj = ii[ok], jj[ok]
        rel[ii, jj] = True
    return rel


def in_kasner_interval(
    coords: FloatArray,
    p: tuple[float, float, float] | FloatArray,
    t_lo: float,
    t_hi: float,
    k: int = 400,
) -> npt.NDArray[np.bool_]:
    """Pertenencia a J^+(a) cap J^-(b), a = (t_lo, 0), b = (t_hi, 0), con la relacion discreta."""
    _check_coords(coords)
    pv = validate_kasner_exponents(p)
    if not 0.0 < t_lo < t_hi:
        raise ValueError("se requiere 0 < t_lo < t_hi")
    dirs = kasner_directions(k)
    ends = np.zeros((2, 4))
    ends[0, 0], ends[1, 0] = t_lo, t_hi
    allc = np.concatenate([ends, np.asarray(coords, dtype=np.float64)], axis=0)
    a_key, b_key = _dominance_keys(allc, pv, dirs)
    a0, b0 = a_key[0], b_key[0]
    a1, b1 = a_key[1], b_key[1]
    ak, bk = a_key[2:], b_key[2:]
    after_a = np.all(ak <= a0[None, :], axis=1) & np.all(bk <= b0[None, :], axis=1)
    before_b = np.all(a1[None, :] <= ak, axis=1) & np.all(b1[None, :] <= bk, axis=1)
    tt = coords[:, 0]
    inside: npt.NDArray[np.bool_] = after_a & before_b & (tt >= t_lo) & (tt <= t_hi)
    return inside


def _box_half_widths(pv: FloatArray, t_lo: float, t_hi: float) -> FloatArray:
    """Semianchos |x_i| del intervalo: max_t min(h(e_i; t_lo, t), h(e_i; t, t_hi)) = h(e_i; t_lo, t_hi)/2."""
    h = support_table(np.array([t_lo, t_hi]), pv, np.eye(3))
    total = h[1] - h[0]
    out: FloatArray = 0.5 * total * (1.0 + 1e-9) + 1e-12
    return out


def sprinkle_kasner_interval(
    n: int,
    p: tuple[float, float, float] | FloatArray,
    t_lo: float,
    t_hi: float,
    rng: np.random.Generator,
    k: int = 400,
) -> FloatArray:
    """n puntos uniformes en volumen propio (densidad coordenada prop. a t) en J^+(a) cap J^-(b).

    Rechazo: caja [t_lo, t_hi] x prod_i [-w_i, w_i] que contiene al poliedro discreto (los ejes
    estan entre las direcciones), aceptacion con probabilidad t / t_hi y pertenencia con la
    misma relacion de K direcciones. Columna 0 = t.
    """
    if n < 0:
        raise ValueError("n debe ser >= 0")
    pv = validate_kasner_exponents(p)
    if not 0.0 < t_lo < t_hi:
        raise ValueError("se requiere 0 < t_lo < t_hi")
    w = _box_half_widths(pv, t_lo, t_hi)
    out = np.empty((n, 4), dtype=np.float64)
    filled = 0
    while filled < n:
        batch = min(8192, max(256, 4 * (n - filled)))
        cand = np.empty((batch, 4))
        cand[:, 0] = rng.uniform(t_lo, t_hi, size=batch)
        cand[:, 1:] = rng.uniform(-1.0, 1.0, size=(batch, 3)) * w[None, :]
        keep = rng.random(batch) * t_hi < cand[:, 0]
        cand = cand[keep]
        ok = cand[in_kasner_interval(cand, pv, t_lo, t_hi, k)]
        take = min(len(ok), n - filled)
        out[filled : filled + take] = ok[:take]
        filled += take
    return out
