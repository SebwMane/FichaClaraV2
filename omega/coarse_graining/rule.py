"""Regla de coarse-graining preregistrada `heavy_edge_matching/max/v1` (OMEGA_1_1_DESIGN §1.7).

1. Aristas con W > w_min ordenadas por (-W, u), u~U(0,1) del rng (desempate aleatorio).
2. Emparejamiento greedy maximal; los no emparejados quedan como bloques unitarios.
3. W'_AB = max_{i in A, j in B} W_ij ("mean" es la alternativa informada).
4. Bloques canónicos: tuplas ordenadas, ordenadas por su mínimo.
5. Niveles hasta n' < n_min o max_levels.

Con pesos sin empates el resultado es equivariante exacto bajo permutación de etiquetas;
con empates lo es en distribución. Referencias: Karypis & Kumar (1998), heavy-edge matching.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Literal

import numpy as np

from omega.config.settings11 import CoarseGrainConfig
from omega.contracts import CoarseLevel
from omega.network.weights import validate_weight_matrix
from omega.types import FloatArray

__all__ = ["heavy_edge_matching", "aggregate", "coarse_grain", "coarse_grain_hierarchy"]


def _check_rng(rng: np.random.Generator) -> None:
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng debe ser numpy.random.Generator")


def heavy_edge_matching(
    w: FloatArray, w_min: float, rng: np.random.Generator
) -> tuple[tuple[int, ...], ...]:
    """Partición en bloques de tamaño <= 2 por emparejamiento greedy de aristas pesadas."""
    validate_weight_matrix(w)
    if not 0.0 <= w_min <= 1.0:
        raise ValueError("w_min debe estar en [0,1]")
    _check_rng(rng)
    n = w.shape[0]
    iu, ju = np.triu_indices(n, k=1)
    vals = w[iu, ju]
    sel = np.flatnonzero(vals > w_min)
    u = rng.random(sel.size)
    order = np.lexsort((u, -vals[sel]))  # clave primaria -W, secundaria u
    matched = np.zeros(n, dtype=bool)
    pairs: list[tuple[int, int]] = []
    for e in sel[order]:
        a, b = int(iu[e]), int(ju[e])
        if not matched[a] and not matched[b]:
            matched[a] = matched[b] = True
            pairs.append((a, b))
    blocks: list[tuple[int, ...]] = [tuple(sorted(p)) for p in pairs]
    blocks.extend((int(i),) for i in np.flatnonzero(~matched))
    blocks.sort(key=lambda b: b[0])
    return tuple(blocks)


def aggregate(
    w: FloatArray, blocks: Sequence[Sequence[int]], how: Literal["max", "mean"]
) -> FloatArray:
    """W' entre bloques: máximo (o media) de W_ij con i en A, j en B; diagonal 0."""
    validate_weight_matrix(w)
    if how not in ("max", "mean"):
        raise ValueError(f"how debe ser 'max' o 'mean', recibido {how!r}")
    n = w.shape[0]
    flat = [int(i) for b in blocks for i in b]
    if any(len(b) == 0 for b in blocks) or sorted(flat) != list(range(n)):
        raise ValueError("blocks debe ser una partición de range(n) sin bloques vacíos")
    sizes = np.array([len(b) for b in blocks], dtype=np.int64)
    starts = np.concatenate(([0], np.cumsum(sizes)[:-1]))
    perm = np.array(flat, dtype=np.int64)
    ws = w[np.ix_(perm, perm)]
    if how == "max":
        out = np.maximum.reduceat(np.maximum.reduceat(ws, starts, axis=0), starts, axis=1)
    else:
        s = np.add.reduceat(np.add.reduceat(ws, starts, axis=0), starts, axis=1)
        out = s / np.outer(sizes, sizes)
        out = 0.5 * (out + out.T)  # elimina asimetría por orden de suma en coma flotante
    out = np.array(out, dtype=np.float64)
    np.fill_diagonal(out, 0.0)
    return out


def coarse_grain(
    w: FloatArray, w_min: float, cfg: CoarseGrainConfig, rng: np.random.Generator, level: int = 1
) -> CoarseLevel:
    """Un paso de coarse-graining; `blocks` indexa los nodos de `w` (nivel anterior)."""
    if not isinstance(cfg, CoarseGrainConfig):
        raise TypeError("cfg debe ser CoarseGrainConfig")
    blocks = heavy_edge_matching(w, w_min, rng)
    w2 = aggregate(w, blocks, cfg.aggregation)
    return CoarseLevel(level=level, w=w2, blocks=blocks, n=len(blocks))


def coarse_grain_hierarchy(
    w: FloatArray, w_min: float, cfg: CoarseGrainConfig, rng: np.random.Generator
) -> tuple[CoarseLevel, ...]:
    """Niveles 1..max_levels; se detiene al producir n' < n_min (ese nivel se descarta por no
    medible) o si el emparejamiento no reduce n (sin aristas > w_min)."""
    if not isinstance(cfg, CoarseGrainConfig):
        raise TypeError("cfg debe ser CoarseGrainConfig")
    levels: list[CoarseLevel] = []
    cur = w
    for lv in range(1, cfg.max_levels + 1):
        nxt = coarse_grain(cur, w_min, cfg, rng, level=lv)
        if nxt.n >= cur.shape[0] or nxt.n < cfg.n_min:
            break
        levels.append(nxt)
        cur = nxt.w
    return tuple(levels)
