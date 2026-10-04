"""Fabricas de RunEvidence sinteticos y validos para tests de Omega-1.1 (sin E/S, sin aleatoriedad global)."""

from __future__ import annotations

import dataclasses
from typing import Any

import numpy as np

from omega.config.settings11 import Omega11Config
from omega.contracts import (
    AnnulusReport,
    BettiCurves,
    BettiResult,
    CurvatureSummary,
    DistanceSuiteResult,
    HomogeneityReport,
    IsotropyReport,
    LocalityReport,
    PersistenceH0,
    RunEvidence,
    ShortCycleBetti,
    TopologySummary,
)
from omega.config.settings11 import DistanceMode
from omega.types import DimensionEstimate, GeometryObservables, RunStatus, TopologyObservables
from tests.conftest import make_config

__all__ = [
    "make_cfg11",
    "make_dim",
    "make_topology",
    "make_evidence",
    "with_dimensions",
    "make_series",
]


def make_cfg11() -> Omega11Config:
    return Omega11Config(base=make_config())


def make_dim(value: float = 3.0, *, status: str = "ok", plateau: bool = True, method: str = "synthetic") -> DimensionEstimate:
    return DimensionEstimate(
        value=value,
        stderr=0.01,
        status=status,  # type: ignore[arg-type]
        window=(1, 5),
        scales=np.arange(1.0, 6.0),
        profile=np.arange(1.0, 6.0),
        local_slopes=np.full(4, value),
        plateau=plateau,
        method=method,
    )


def make_topology(**kw: Any) -> TopologyObservables:
    base: dict[str, Any] = dict(
        n=800, mean_strength=8.0, std_strength=1.0, cv_strength=0.1, mean_binary_degree=12.0, binary_density=0.015,
        mean_weight=0.015, giant_size=800, giant_fraction=1.0, n_components=1, n_large_components=1,
        clustering_weighted=0.4, clustering_binary=0.4, frac_at_zero=0.9, frac_at_one=0.0,
    )
    base.update(kw)
    return TopologyObservables(**base)


def _topo_summary(stable: bool = True) -> TopologySummary:
    thetas = np.array([0.5, 0.1], dtype=np.float64)
    curves = BettiCurves(
        thetas=thetas,
        beta0=np.array([1, 1], dtype=np.int64),
        beta1_short=np.array([0, 0], dtype=np.int64),
        b1_density=np.array([0.0, 0.0], dtype=np.float64),
        giant_fraction=np.array([1.0, 1.0], dtype=np.float64),
        status=("ok", "ok"),
    )
    h0 = PersistenceH0(births=np.ones(3), deaths=np.array([0.5, 0.4, 0.3]), n_essential=1)
    return TopologySummary(
        curves=curves,
        h0=h0,
        at_w_min=ShortCycleBetti(b0=1, b1=0, n_edges=100, n_faces=50, b1_density=0.0, max_length=4, status="ok"),
        clique_at_w_min=BettiResult(betti=(1, 0, 0), counts=(10, 20, 5), euler_betti=1, euler_counts=-5, status="ok"),
        stable_flags=(stable,) * 5,
        stable=stable,
    )


def make_evidence(**overrides: Any) -> RunEvidence:
    """Estado "pasa todo": D_vol = D_s = D_W = 3.0, clase 3, conexo, local, estable, no trivial."""
    ok_dim = make_dim(3.0)
    base: dict[str, Any] = dict(
        status=RunStatus.CONVERGED,
        topology=make_topology(),
        geometry=GeometryObservables(
            path_length=5.0, path_length_hops=5.0, diameter=12.0,
            d_eff=ok_dim, d_eff_ball=ok_dim, d_eff_hops=ok_dim, d_s=ok_dim,
        ),
        suite=DistanceSuiteResult(
            estimates={DistanceMode.HOP: ok_dim},
            fallback_used={DistanceMode.HOP: False},
            resistance_exponent=0.2,
            metric_spread=0.05,
            distance_sensitive=False,
            reason="synthetic",
        ),
        d_weyl=make_dim(3.0, method="weyl"),
        homogeneity=HomogeneityReport(mu=2.9, sigma=0.1, cv=0.05, n_valid=700, ok=True),
        isotropy=IsotropyReport(median_ratio=0.7, p10_ratio=0.5, radius=4, k=3, n_sources=64, ok=True),
        locality=LocalityReport(detour_fraction=0.99, ok=True),
        annulus=AnnulusReport(fractions={2: 0.97, 3: 0.98}, ok=True),
        topo=_topo_summary(True),
        curvature=CurvatureSummary(
            mean=0.0, std=0.05, se=0.002, tail_fraction=0.0, n_edges=1000, sampled=True,
            edge_values=np.zeros(4), ok=True,
        ),
        qrc=None,
        weight_mean=0.015,
        weight_cv=0.3,
        consensus_dimension=3.0,
        dimension_class=3,
        clustering_ratio=1.5,
        fiedler_length=4.0,
        n=800,
    )
    base.update(overrides)
    return RunEvidence(**base)


def with_dimensions(ev: RunEvidence, d_star: float, *, cls: int | None = None) -> RunEvidence:
    """Copia con los tres estimadores y D* iguales a `d_star` (clase = round si no se da y esta a <= 0.35)."""
    dim = make_dim(d_star)
    geo = dataclasses.replace(ev.geometry, d_eff=dim, d_s=dim)
    if cls is None:
        r = round(d_star)
        cls = r if abs(d_star - r) <= 0.35 and r >= 1 else None
    return dataclasses.replace(
        ev, geometry=geo, d_weyl=make_dim(d_star, method="weyl"), consensus_dimension=d_star, dimension_class=cls
    )


def make_series(d_star: float, count: int, *, cv: float = 0.05, **overrides: Any) -> list[RunEvidence]:
    """`count` corridas "pasa todo" con D* = d_star y homogeneidad cv dada."""
    ev = with_dimensions(make_evidence(**overrides), d_star)
    hom = HomogeneityReport(mu=2.9, sigma=0.1 * cv, cv=cv, n_valid=700, ok=True)
    ev = dataclasses.replace(ev, homogeneity=hom)
    return [ev for _ in range(count)]
