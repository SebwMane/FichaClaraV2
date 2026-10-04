"""Utilidades puras sobre matrices de pesos W (ANALYSIS §5, M§4-M§8)."""

from __future__ import annotations

import numpy as np

from omega.types import FloatArray, IntArray

__all__ = [
    "validate_weight_matrix",
    "upper_triangle",
    "from_upper_triangle",
    "clip_unit",
    "permute",
]


def validate_weight_matrix(w: FloatArray, *, atol: float = 1e-12) -> None:
    """Exige (N,N) float64, N>=2, finita, simetrica (atol), diagonal 0 y valores en [0,1]."""
    if not isinstance(w, np.ndarray):
        raise TypeError("w debe ser numpy.ndarray")
    if w.dtype != np.float64:
        raise TypeError(f"w debe ser float64, recibido {w.dtype}")
    if w.ndim != 2 or w.shape[0] != w.shape[1]:
        raise ValueError(f"w debe ser cuadrada (N,N), forma {w.shape}")
    if w.shape[0] < 2:
        raise ValueError("w requiere N >= 2")
    if not np.all(np.isfinite(w)):
        raise ValueError("w contiene valores no finitos")
    if not np.allclose(w, w.T, rtol=0.0, atol=atol):
        raise ValueError("w no es simetrica")
    if np.any(np.diag(w) != 0.0):
        raise ValueError("la diagonal de w debe ser 0")
    if w.min() < 0.0 or w.max() > 1.0:
        raise ValueError("los valores de w deben estar en [0,1]")


def upper_triangle(w: FloatArray) -> FloatArray:
    """Vector (M,) de aristas w_ab con a<b (orden por filas), M=N(N-1)/2."""
    validate_weight_matrix(w)
    iu = np.triu_indices(w.shape[0], k=1)
    return np.asarray(w[iu], dtype=np.float64).copy()


def from_upper_triangle(v: FloatArray, n: int) -> FloatArray:
    """Matriz simetrica (n,n), diag 0, desde el vector de aristas a<b; len(v)=n(n-1)/2."""
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise ValueError(f"n debe ser entero >= 2, recibido {n!r}")
    v = np.asarray(v, dtype=np.float64)
    if v.ndim != 1 or v.shape[0] != n * (n - 1) // 2:
        raise ValueError(f"v debe tener longitud {n * (n - 1) // 2}, forma {v.shape}")
    w = np.zeros((n, n), dtype=np.float64)
    iu = np.triu_indices(n, k=1)
    w[iu] = v
    w[(iu[1], iu[0])] = v
    return w


def clip_unit(w: FloatArray) -> FloatArray:
    """Copia de w recortada a [0,1]; no muta la entrada."""
    return np.clip(np.asarray(w, dtype=np.float64), 0.0, 1.0)


def permute(w: FloatArray, perm: IntArray) -> FloatArray:
    """P W P^T con (PWP^T)_ij = w[perm[i], perm[j]]; perm debe ser permutacion de 0..N-1."""
    validate_weight_matrix(w)
    p = np.asarray(perm)
    n = w.shape[0]
    if p.ndim != 1 or p.shape[0] != n or not np.issubdtype(p.dtype, np.integer):
        raise ValueError(f"perm debe ser vector entero de longitud {n}")
    if not np.array_equal(np.sort(p), np.arange(n)):
        raise ValueError("perm no es una permutacion de 0..N-1")
    return np.asarray(w[np.ix_(p, p)], dtype=np.float64).copy()
