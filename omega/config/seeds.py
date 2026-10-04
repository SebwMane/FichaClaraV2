"""Semillas reproducibles e independientes del orden (D-28)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from omega.config.settings import SeedConfig

__all__ = ["SeedKey", "seed_key", "make_rng"]


@dataclass(frozen=True, slots=True)
class SeedKey:
    """Identifica un flujo aleatorio: entropia y spawn_key de SeedSequence."""

    entropy: int
    spawn_key: tuple[int, ...]

    def __post_init__(self) -> None:
        if isinstance(self.entropy, bool) or not isinstance(self.entropy, int) or self.entropy < 0:
            raise ValueError(f"entropy debe ser entero >= 0, recibido {self.entropy!r}")
        if not isinstance(self.spawn_key, tuple):
            raise TypeError("spawn_key debe ser tuple")
        for k in self.spawn_key:
            if isinstance(k, bool) or not isinstance(k, int) or k < 0:
                raise ValueError(f"spawn_key debe contener enteros >= 0, recibido {k!r}")


def seed_key(cfg: SeedConfig, point_index: int, replicate: int) -> SeedKey:
    """SeedKey(entropy, (experiment_id, point_index, replicate)); indices >= 0 y replicate < replicates."""
    for name, v in (("point_index", point_index), ("replicate", replicate)):
        if isinstance(v, bool) or not isinstance(v, int):
            raise TypeError(f"{name} debe ser entero")
        if v < 0:
            raise ValueError(f"{name} debe ser >= 0, recibido {v}")
    if replicate >= cfg.replicates:
        raise ValueError(f"replicate={replicate} fuera de rango (replicates={cfg.replicates})")
    return SeedKey(cfg.master_entropy, (cfg.experiment_id, point_index, replicate))


def make_rng(key: SeedKey) -> np.random.Generator:
    """Generator(PCG64(SeedSequence(entropy, spawn_key))): determinista por clave."""
    if not isinstance(key, SeedKey):
        raise TypeError("key debe ser SeedKey")
    ss = np.random.SeedSequence(key.entropy, spawn_key=key.spawn_key)
    return np.random.Generator(np.random.PCG64(ss))
