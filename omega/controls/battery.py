"""Batería de nulos contra un candidato (OMEGA_1_1_DESIGN §1.10, §3.5)."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from omega.config.seeds import SeedKey, make_rng
from omega.config.settings11 import NullModel
from omega.controls.configuration_model import matched_configuration_model
from omega.controls.degree_preserving_rewire import degree_preserving_rewire
from omega.controls.erdos_renyi import matched_erdos_renyi
from omega.controls.random_geometric import rgg_torus
from omega.controls.small_world import matched_small_world
from omega.network.weights import validate_weight_matrix
from omega.phases.stability import shuffled_weight_null
from omega.types import FloatArray

__all__ = ["null_seed_key", "make_null", "null_battery"]

_WS_P_DEFAULT = 0.1
_RGG_DIM_DEFAULT = 2


def null_seed_key(base: SeedKey, kind: NullModel) -> SeedKey:
    """SeedKey independiente por tipo: spawn_key + (1000 + índice del enum,)."""
    if not isinstance(base, SeedKey):
        raise TypeError("base debe ser SeedKey")
    if not isinstance(kind, NullModel):
        raise TypeError("kind debe ser NullModel")
    idx = list(NullModel).index(kind)
    return SeedKey(base.entropy, base.spawn_key + (1000 + idx,))


def make_null(w: FloatArray, w_min: float, kind: NullModel, rng: np.random.Generator) -> FloatArray:
    """Genera un nulo del tipo `kind` emparejado con `w`.

    SMALL_WORLD usa p=0.1; RANDOM_GEOMETRIC usa el toro T^2 binario con el k̄ de A.
    """
    validate_weight_matrix(w)
    if not isinstance(kind, NullModel):
        raise TypeError("kind debe ser NullModel")
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng debe ser numpy.random.Generator")
    if kind is NullModel.ERDOS_RENYI:
        return matched_erdos_renyi(w, w_min, rng, "binary")
    if kind is NullModel.CONFIGURATION_MODEL:
        return matched_configuration_model(w, w_min, rng)
    if kind is NullModel.DEGREE_PRESERVING_REWIRE:
        return degree_preserving_rewire(w, w_min, rng)
    if kind is NullModel.SHUFFLED_WEIGHTS:
        return shuffled_weight_null(w, w_min, rng)
    if kind is NullModel.SMALL_WORLD:
        return matched_small_world(w, w_min, _WS_P_DEFAULT, rng)
    n = w.shape[0]
    kbar = float((w > w_min).sum()) / n
    return rgg_torus(n, _RGG_DIM_DEFAULT, max(kbar, 1e-9), rng, "binary")


def null_battery(
    w: FloatArray, w_min: float, kinds: Sequence[NullModel], base: SeedKey
) -> dict[NullModel, FloatArray]:
    """Un nulo por tipo, cada uno con su propio flujo `null_seed_key(base, kind)`."""
    out: dict[NullModel, FloatArray] = {}
    for kind in kinds:
        out[kind] = make_null(w, w_min, kind, make_rng(null_seed_key(base, kind)))
    return out
