"""Taxonomia de fallos Omega-F0..Omega-F10 y banderas por corrida (OMEGA_1_1_DESIGN §1.11, §1.12, §3.6).

Logica pura: recibe `RunEvidence` ya construida; no calcula geometria. Evidencia ausente o no finita cuenta como
False en toda bandera que deba cumplirse.
"""

from __future__ import annotations

import math
import statistics
from collections.abc import Iterable, Mapping, Sequence

from omega.config.settings11 import DistanceSuiteConfig, Omega11Config
from omega.contracts import (
    FAILURE_PRECEDENCE,
    FailureCode,
    RunAssessment11,
    RunEvidence,
)
from omega.types import DimensionEstimate, RunStatus

__all__ = [
    "RUN_FLAG_KEYS",
    "consensus_dimension",
    "resistance_consistent",
    "run_flags",
    "run_failure_codes",
    "primary_code",
    "assess_run",
]

RUN_FLAG_KEYS: tuple[str, ...] = (
    "converged",
    "empty",
    "dense_or_uniform",
    "nontrivial",
    "connected",
    "small_world",
    "locality",
    "d_volume",
    "d_spectral",
    "d_weyl",
    "metric_robust",
    "isotropy_ok",
    "homogeneity_ok",
    "topology_stable",
    "manifold_proxy_ok",
)


def _fin(*xs: float) -> bool:
    return all(math.isfinite(x) for x in xs)


def consensus_dimension(estimates: Sequence[DimensionEstimate], class_tol: float) -> tuple[float, int | None]:
    """D* = mediana de las estimaciones, solo si TODAS son "ok" y finitas; si no, (NaN, None).

    La clase es round(D*) si |D* - round(D*)| <= class_tol; si no (o si round < 1), None.
    """
    if len(estimates) == 0 or any(e.status != "ok" or not math.isfinite(e.value) for e in estimates):
        return math.nan, None
    d_star = float(statistics.median(e.value for e in estimates))
    cls = round(d_star)
    if cls < 1 or abs(d_star - cls) > class_tol:
        return d_star, None
    return d_star, int(cls)


def resistance_consistent(zeta: float, dimension_class: int | None, suite: DistanceSuiteConfig) -> bool:
    """Exponente de resistencia coherente con la clase (§1.1). Clase None o zeta no finito: False."""
    if dimension_class is None or not math.isfinite(zeta):
        return False
    if dimension_class == 1:
        return zeta >= suite.zeta_1d_min
    if dimension_class == 2:
        return suite.zeta_2d_min <= zeta <= suite.zeta_2d_max
    if dimension_class >= 3:
        return zeta < suite.zeta_3d_max
    return False


def _dim_ok(est: DimensionEstimate, d_star: float, tol: float) -> bool:
    return (
        est.status == "ok"
        and est.plateau is True
        and _fin(est.value, d_star)
        and abs(est.value - d_star) <= tol
    )


def run_flags(ev: RunEvidence, cfg: Omega11Config) -> dict[str, bool]:
    """Las 15 banderas de una corrida (claves en `RUN_FLAG_KEYS`). Todas son `bool` estrictos."""
    c = cfg.certificate
    top = ev.topology
    d_star = ev.consensus_dimension

    converged = ev.status is RunStatus.CONVERGED
    empty = bool(
        _fin(top.giant_fraction, top.mean_binary_degree)
        and top.giant_fraction < c.empty_g
        and top.mean_binary_degree < c.empty_kbin
    )
    dense = bool(
        (_fin(top.binary_density) and top.binary_density >= c.dense_rho)
        or (_fin(ev.weight_mean) and ev.weight_mean >= c.dense_meanw)
        or (_fin(top.mean_weight) and top.mean_weight >= c.dense_meanw)
        or (_fin(ev.weight_cv, ev.weight_mean) and ev.weight_cv <= c.uniform_cv_max and ev.weight_mean > 0.0)
    )
    # Conservador: evidencia no finita en lo que decide trivialidad impide declarar "nontrivial".
    trivial_evidence_ok = _fin(
        top.giant_fraction, top.mean_binary_degree, top.binary_density, top.mean_weight, ev.weight_mean, ev.weight_cv
    )
    nontrivial = (not empty) and (not dense) and trivial_evidence_ok
    connected = bool(_fin(top.giant_fraction) and top.giant_fraction >= c.g_connected)

    cv_ok = _fin(ev.clustering_ratio, ev.curvature.tail_fraction)
    small_world = bool(
        cv_ok
        and ev.clustering_ratio >= cfg.local.small_world_clustering_ratio
        and ev.curvature.tail_fraction > cfg.curvature.tail_frac_max
    )
    locality = bool(
        cv_ok
        and _fin(ev.locality.detour_fraction)
        and ev.locality.detour_fraction >= cfg.local.locality_min
        and not small_world
    )

    d_volume = _dim_ok(ev.geometry.d_eff, d_star, c.dim_tol)
    d_spectral = _dim_ok(ev.geometry.d_s, d_star, c.dim_tol)
    # Para D_Weyl, `plateau` es R^2 >= r2_min (§1.3).
    d_weyl = _dim_ok(ev.d_weyl, d_star, c.dim_tol)

    metric_robust = bool(
        (not ev.suite.distance_sensitive)
        and resistance_consistent(ev.suite.resistance_exponent, ev.dimension_class, cfg.distance)
    )

    isotropy_ok = bool(ev.isotropy.ok)
    h = ev.homogeneity
    homogeneity_ok = bool(
        h.ok and _fin(h.cv) and h.cv <= cfg.local.homogeneity_cv_max and h.n_valid >= cfg.local.homogeneity_min_nodes
    )

    flags_theta = ev.topo.stable_flags
    topology_stable = bool(
        ev.topo.stable
        and len(flags_theta) > 0
        and sum(flags_theta) / len(flags_theta) >= cfg.topology.stable_fraction
    )

    cur = ev.curvature
    w = ev.topo.at_w_min
    manifold_proxy_ok = bool(
        ev.annulus.ok
        and cur.ok
        and _fin(cur.mean, cur.tail_fraction, w.b1_density)
        and cur.mean >= cfg.curvature.mean_min
        and cur.tail_fraction <= cfg.curvature.tail_frac_max
        and w.status == "ok"
        and w.b1_density <= cfg.topology.b1_density_max
    )

    return {
        "converged": converged,
        "empty": empty,
        "dense_or_uniform": dense,
        "nontrivial": nontrivial,
        "connected": connected,
        "small_world": small_world,
        "locality": locality,
        "d_volume": d_volume,
        "d_spectral": d_spectral,
        "d_weyl": d_weyl,
        "metric_robust": metric_robust,
        "isotropy_ok": isotropy_ok,
        "homogeneity_ok": homogeneity_ok,
        "topology_stable": topology_stable,
        "manifold_proxy_ok": manifold_proxy_ok,
    }


def _flag(flags: Mapping[str, bool], key: str) -> bool | None:
    """Valor estricto de la bandera; None si falta. np.bool_ u otros tipos: TypeError."""
    if key not in flags:
        return None
    v = flags[key]
    if type(v) is not bool:
        raise TypeError(f"la bandera {key!r} debe ser bool estricto, recibido {type(v).__name__}")
    return v


def run_failure_codes(flags: Mapping[str, bool]) -> tuple[FailureCode, ...]:
    """Conjunto completo de codigos de una corrida, en orden de precedencia. Bandera ausente = no cumplida.

    F0/F1 se activan cuando `empty`/`dense_or_uniform` son True. Si `nontrivial` no se cumple sin que haya
    causa F0/F1 explicita (evidencia no finita o bandera ausente), se informa F1 (opcion conservadora).
    """
    empty = _flag(flags, "empty") is True
    dense = _flag(flags, "dense_or_uniform") is True
    nontrivial = _flag(flags, "nontrivial") is True

    def fails(*keys: str) -> bool:
        return any(_flag(flags, k) is not True for k in keys)

    codes: set[FailureCode] = set()
    if fails("converged"):
        codes.add(FailureCode.F10)
    if empty:
        codes.add(FailureCode.F0)
    if dense or (not nontrivial and not empty):
        codes.add(FailureCode.F1)
    if fails("connected"):
        codes.add(FailureCode.F2)
    if fails("locality"):
        codes.add(FailureCode.F3)
    if fails("d_volume", "d_spectral", "d_weyl"):
        codes.add(FailureCode.F4)
    if fails("metric_robust"):
        codes.add(FailureCode.F5)
    if fails("isotropy_ok", "homogeneity_ok", "topology_stable", "manifold_proxy_ok"):
        codes.add(FailureCode.F9)
    return tuple(c for c in FAILURE_PRECEDENCE if c in codes)


def primary_code(codes: Iterable[FailureCode]) -> FailureCode | None:
    """Codigo primario segun `FAILURE_PRECEDENCE` (F10 > F0 > F1 > F2 > F3 > F4 > F5 > F7 > F6 > F8 > F9)."""
    present = set(codes)
    for c in FAILURE_PRECEDENCE:
        if c in present:
            return c
    return None


def assess_run(ev: RunEvidence, cfg: Omega11Config) -> RunAssessment11:
    flags = run_flags(ev, cfg)
    codes = run_failure_codes(flags)
    return RunAssessment11(codes=codes, primary=primary_code(codes), flags=flags, passes=len(codes) == 0)
