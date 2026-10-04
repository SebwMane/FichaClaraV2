"""Fixtures compartidas (solo dependen de WP1)."""

from __future__ import annotations

import numpy as np
import pytest

from omega.config.settings import (
    DimensionConfig,
    DynamicsConfig,
    FunctionalParams,
    GraphConfig,
    InitConfig,
    OmegaConfig,
    PhaseThresholds,
    ScanConfig,
    SeedConfig,
    SpectralConfig,
)


def make_config(n: int = 12, scan: bool = True) -> OmegaConfig:
    """Configuracion por defecto de prueba."""
    return OmegaConfig(
        init=InitConfig(n=n),
        functional=FunctionalParams(alpha=0.3, beta=1.0, gamma=0.05),
        dynamics=DynamicsConfig(),
        graph=GraphConfig(),
        dimension=DimensionConfig(),
        spectral=SpectralConfig(),
        phases=PhaseThresholds(),
        seeds=SeedConfig(master_entropy=12345, experiment_id=7, replicates=3),
        scan=ScanConfig(alpha_hat=(0.0, 1.0, 2.0), gamma_hat=(0.0, 1.0)) if scan else None,
    )


@pytest.fixture
def cfg() -> OmegaConfig:
    return make_config()


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(2024)))
