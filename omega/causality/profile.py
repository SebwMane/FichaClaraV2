"""Perfil de abundancias de intervalos (C0-4b, Rev. 3 §9)."""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
import numpy.typing as npt

from omega.causality.estimators import interval_abundances
from omega.causality.exact import expected_abundance_2order
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond

__all__ = ["K_MAX", "chain_abundances", "reference_abundances", "profile_distance", "d_ref"]

K_MAX = 8
FloatArray = npt.NDArray[np.float64]


def chain_abundances(n: int, k_max: int = K_MAX) -> FloatArray:
    """Cadena de n elementos: N_k = n - 1 - k (0 si es negativo)."""
    return np.array([max(n - 1 - k, 0) for k in range(k_max + 1)], dtype=np.float64)


def reference_abundances(
    d: int, n: int, k_max: int = K_MAX, rngs: Sequence[np.random.Generator] = ()
) -> FloatArray:
    """d = 1 cadena exacta; d = 2 exacta (2-ordenes); d = 3, 4 media sobre sprinklings con los rngs dados."""
    if d == 1:
        return chain_abundances(n, k_max)
    if d == 2:
        return np.array([expected_abundance_2order(k, n) for k in range(k_max + 1)], dtype=np.float64)
    if d in (3, 4):
        if not rngs:
            raise ValueError("d = 3, 4 requiere generadores para la referencia")
        rows = [interval_abundances(causal_matrix(sprinkle_diamond(n, d, r)), k_max) for r in rngs]
        out: FloatArray = np.mean(np.array(rows, dtype=np.float64), axis=0)
        return out
    raise ValueError("d debe estar en {1, 2, 3, 4}")


def profile_distance(nk: npt.NDArray[np.int64] | FloatArray, ref: FloatArray) -> float:
    """max_k |ln((N_k + 1) / (ref_k + 1))|."""
    a = np.asarray(nk, dtype=np.float64)
    if a.shape != ref.shape:
        raise ValueError("perfiles de distinta longitud")
    return float(np.max(np.abs(np.log((a + 1) / (ref + 1)))))


def d_ref(d_hat: float) -> int | None:
    """Entero mas cercano recortado a {1..4}; None si d_hat no es finito."""
    if not math.isfinite(d_hat):
        return None
    return int(min(4, max(1, round(d_hat))))
