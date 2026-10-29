"""Bloque T^3 x [0, T] con espacio periodico de lado 1 (L4b, Rev. 3 §20)."""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt

from omega.causality.estimators import between_counts, smearing_weight_4d
from omega.causality.sprinkling import interval_volume

__all__ = ["SLAB_T", "sprinkle_slab", "slab_causal_matrix", "retarded_action_4d", "slab_pair_fractions"]

SLAB_T = 0.4
FloatArr = npt.NDArray[np.float64]
BoolArr = npt.NDArray[np.bool_]


def sprinkle_slab(n: int, rng: np.random.Generator, t_len: float = SLAB_T) -> FloatArr:
    out = np.empty((n, 4))
    out[:, 0] = rng.uniform(0, t_len, n)
    out[:, 1:] = rng.uniform(0, 1, (n, 3))
    return out


def _min_image(d: FloatArr) -> FloatArr:
    out: FloatArr = (d + 0.5) % 1.0 - 0.5
    return out


def slab_causal_matrix(x: FloatArr) -> BoolArr:
    if x[:, 0].max() - x[:, 0].min() >= 0.5:
        raise ValueError("la duracion debe ser < 1/2 para que la imagen minima sea la unica")
    dt = x[None, :, 0] - x[:, None, 0]
    dx2 = np.zeros_like(dt)
    for k in range(1, 4):
        d = _min_image(x[None, :, k] - x[:, None, k])
        dx2 += d * d
    rel: BoolArr = (dt > 0) & (dt * dt >= dx2)
    return rel


def retarded_action_4d(c: BoolArr, eps: float) -> FloatArr:
    """s(x) = 1 - eps sum_{y<x} f4(n(y,x), eps); con eps = 1 usa los coeficientes (1, -9, 16, -8).

    Densidad local: la accion es sqrt(eps) * sum_x s(x) (bd_action_4d_smeared).
    """
    b = between_counts(c)
    if eps == 1.0:
        coef = np.zeros(int(max(b.max(), 0)) + 1)
        coef[: min(4, coef.size)] = np.array([1.0, -9.0, 16.0, -8.0])[: min(4, coef.size)]
        w = np.where(c, coef[np.clip(b, 0, None)], 0.0)
        out1: FloatArr = 1.0 - w.sum(axis=0)
        return out1
    w = np.where(c, smearing_weight_4d(np.clip(b, 0, None).astype(np.float64), eps), 0.0)
    out: FloatArr = 1.0 - eps * w.sum(axis=0)
    return out


def slab_pair_fractions(pairs: int, rng: np.random.Generator, t_len: float = SLAB_T, chunk: int = 1_000_000) -> FloatArr:
    """Referencia por pares en el bloque: v = c4 tau^4 / V_bloque si relacionados, 0 si no."""
    parts = []
    for k in range(0, pairs, chunk):
        m = min(chunk, pairs - k)
        p, q = sprinkle_slab(m, rng, t_len), sprinkle_slab(m, rng, t_len)
        dt = np.abs(q[:, 0] - p[:, 0])
        dx2 = np.sum(_min_image(q[:, 1:] - p[:, 1:]) ** 2, axis=1)
        rel = dt * dt >= dx2
        tau = np.sqrt(np.maximum(dt * dt - dx2, 0.0))
        parts.append(np.where(rel, interval_volume(1.0, 4) * tau**4 / t_len, 0.0))
    out: FloatArr = np.concatenate(parts)
    return out


_ = math  # reservado
