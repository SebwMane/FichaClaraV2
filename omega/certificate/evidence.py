"""Recoleccion de evidencia por corrida (OMEGA_1_1_DESIGN §3.7, §1.12).

Une los modulos de geometria, topologia y curvatura de la oleada 1 en un `RunEvidence`.
`evidence.geometry` es exactamente `geometry_observables` de Omega-1.0 (igual que `scan.observe`).
Solo calcula; no decide (la decision es de `omega.certificate.taxonomy`).
"""

from __future__ import annotations

import hashlib
import math
from typing import Any

import numpy as np

from omega.certificate.taxonomy import consensus_dimension
from omega.config.seeds import SeedKey, make_rng
from omega.config.settings11 import DistanceMode, Omega11Config, TopologyConfig
from omega.contracts import BettiCurves, RunEvidence, ShortCycleBetti, TopologySummary
from omega.curvature.ollivier import ollivier_curvature
from omega.curvature.qrc import quantum_ricci_profile
from omega.geometry.distance_suite import distance_suite
from omega.geometry.distances import hop_distance_matrix, threshold_adjacency
from omega.geometry.local_structure import annulus_connectivity, homogeneity, isotropy, locality
from omega.geometry.observables import geometry_observables
from omega.geometry.weyl import fiedler_length, weyl_dimension
from omega.network.topology import component_labels, giant_component_nodes, topology_observables
from omega.network.weights import upper_triangle, validate_weight_matrix
from omega.topology.betti import clique_complex_betti, short_cycle_betti1
from omega.topology.filtration import threshold_graph
from omega.topology.persistence import h0_persistence
from omega.types import BoolArray, FloatArray, RunStatus

__all__ = ["evidence_rng", "weight_stats", "collect_run_evidence", "evidence_summary"]

_EVIDENCE_STREAM = 2000


def evidence_rng(key: SeedKey) -> np.random.Generator:
    """Flujo propio de la evidencia: SeedKey(entropy, spawn_key + (2000,)); independiente de la dinamica y los nulos."""
    if not isinstance(key, SeedKey):
        raise TypeError("key debe ser SeedKey")
    return make_rng(SeedKey(key.entropy, key.spawn_key + (_EVIDENCE_STREAM,)))


def _cached_topology_summary(w: FloatArray, w_min: float, sensitivity: tuple[float, ...], g_connected: float, cfg: TopologyConfig) -> TopologySummary:
    """Mismo resultado que `topology_summary` de WP-C, pero memoizando beta1^(4) por grafo G_theta.

    En estados binarios todos los niveles theta<1 son el mismo grafo; `short_cycle_betti1` (el coste dominante,
    ~1 s con N=800) se evalua una vez por grafo distinto. La malla de sensibilidad reutiliza los niveles ya
    calculados cuando theta pertenece a la malla de Betti.
    """
    n = w.shape[0]
    cache: dict[bytes, ShortCycleBetti] = {}
    levels: dict[float, tuple[int, int, float, float, str]] = {}

    def short(a: BoolArray) -> ShortCycleBetti:
        key = hashlib.sha1(a.tobytes()).digest()
        hit = cache.get(key)
        if hit is None:
            hit = short_cycle_betti1(a, cfg.short_cycle_length, cfg.max_faces)
            cache[key] = hit
        return hit

    def level(theta: float) -> tuple[int, int, float, float, str]:
        got = levels.get(theta)
        if got is None:
            a = threshold_graph(w, theta)
            n_comp, labels = component_labels(a)
            giant = float(np.bincount(labels).max()) / n
            sc = short(a)
            got = (n_comp, sc.b1, sc.b1_density, giant, sc.status)
            levels[theta] = got
        return got

    rows = [level(float(t)) for t in cfg.thetas]
    curves = BettiCurves(
        thetas=np.asarray(cfg.thetas, dtype=np.float64),
        beta0=np.asarray([r[0] for r in rows], dtype=np.int64),
        beta1_short=np.asarray([r[1] for r in rows], dtype=np.int64),
        b1_density=np.asarray([r[2] for r in rows], dtype=np.float64),
        giant_fraction=np.asarray([r[3] for r in rows], dtype=np.float64),
        status=tuple(r[4] for r in rows),
    )
    a_w = threshold_graph(w, w_min)
    flags = []
    for t in sensitivity:
        _, _, dens, giant, status = level(float(t))
        flags.append(bool(giant >= g_connected and dens <= cfg.b1_density_max and status == "ok"))
    stable = bool(flags) and sum(flags) / len(flags) >= cfg.stable_fraction
    return TopologySummary(
        curves=curves,
        h0=h0_persistence(w),
        at_w_min=short(a_w),
        clique_at_w_min=clique_complex_betti(a_w, cfg.clique_max_dim, cfg.max_simplices),
        stable_flags=tuple(flags),
        stable=stable,
    )


def weight_stats(w: FloatArray) -> tuple[float, float]:
    """(media, cv) del triangulo superior de W. Con media 0 el cv es 0.0 (sin dispersion)."""
    validate_weight_matrix(w)
    v = upper_triangle(w)
    mean = float(v.mean())
    if mean <= 0.0:
        return mean, 0.0
    return mean, float(v.std() / mean)


def collect_run_evidence(
    w: FloatArray,
    status: RunStatus,
    cfg: Omega11Config,
    rng: np.random.Generator,
    *,
    with_qrc: bool = False,
) -> RunEvidence:
    """Toda la evidencia de una corrida sobre el estado final `w`.

    Orden fijo de consumo del `rng` (determinismo): isotropia, anillos, curvatura, QRC. Las medidas locales
    (homogeneidad, isotropia, localidad, anillos) se toman sobre la componente gigante de A=(W>w_min) y la
    ventana de saltos es la de `geometry.d_eff_hops` (identica a la de la suite).
    """
    validate_weight_matrix(w)
    if not isinstance(status, RunStatus):
        raise TypeError("status debe ser RunStatus")
    if not isinstance(cfg, Omega11Config):
        raise TypeError("cfg debe ser Omega11Config")
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng debe ser numpy.random.Generator")
    base = cfg.base
    w_min = base.graph.w_min

    top = topology_observables(w, w_min, base.phases.large_component_frac)
    geo = geometry_observables(w, base.graph, base.dimension, base.spectral)
    suite = distance_suite(w, base.graph, base.dimension, cfg.distance)
    d_weyl = weyl_dimension(w, w_min, cfg.weyl)
    d_star, cls = consensus_dimension([geo.d_eff, geo.d_s, d_weyl], cfg.certificate.class_tol)

    a = threshold_adjacency(w, w_min)
    _, labels = component_labels(a)
    nodes = giant_component_nodes(labels)
    if nodes.size >= 2:
        asub = np.ascontiguousarray(a[np.ix_(nodes, nodes)])
        d_hop = hop_distance_matrix(asub)
    else:  # gigante trivial: mismas convenciones que geometry_observables
        asub = np.zeros((1, 1), dtype=np.bool_)
        d_hop = np.zeros((1, 1), dtype=np.float64)

    hop_est = geo.d_eff_hops
    homo = homogeneity(d_hop, hop_est, cfg.local)
    iso = isotropy(d_hop, hop_est, cls, cfg.local, rng)
    loc = locality(asub, cfg.local)
    ann = annulus_connectivity(asub, d_hop, hop_est, cfg.local, rng)
    topo = _cached_topology_summary(
        w, w_min, tuple(base.graph.w_min_sensitivity), cfg.certificate.g_connected, cfg.topology
    )
    cur = ollivier_curvature(w, w_min, base.graph, cfg.distance, cfg.curvature, rng)
    qrc = (
        quantum_ricci_profile(d_hop, cfg.curvature.qrc_deltas, cfg.curvature.qrc_max_pairs, rng)
        if with_qrc and nodes.size >= 2
        else None
    )

    wmean, wcv = weight_stats(w)
    ratio = (
        top.clustering_binary / top.binary_density
        if math.isfinite(top.clustering_binary) and top.binary_density > 0.0
        else math.nan
    )
    return RunEvidence(
        status=status,
        topology=top,
        geometry=geo,
        suite=suite,
        d_weyl=d_weyl,
        homogeneity=homo,
        isotropy=iso,
        locality=loc,
        annulus=ann,
        topo=topo,
        curvature=cur,
        qrc=qrc,
        weight_mean=wmean,
        weight_cv=wcv,
        consensus_dimension=d_star,
        dimension_class=cls,
        clustering_ratio=float(ratio),
        fiedler_length=fiedler_length(w, w_min),
        n=int(w.shape[0]),
    )


def _fin(x: float) -> float | None:
    return float(x) if math.isfinite(x) else None


def evidence_summary(ev: RunEvidence) -> dict[str, Any]:
    """Resumen JSON-compatible (sin matrices; no finitos -> None)."""
    top, geo = ev.topology, ev.geometry
    hop = ev.suite.estimates.get(DistanceMode.HOP)
    return {
        "n": ev.n,
        "status": ev.status.value,
        "giant_fraction": _fin(top.giant_fraction),
        "mean_binary_degree": _fin(top.mean_binary_degree),
        "binary_density": _fin(top.binary_density),
        "mean_weight": _fin(ev.weight_mean),
        "weight_cv": _fin(ev.weight_cv),
        "clustering_ratio": _fin(ev.clustering_ratio),
        "fiedler_length": _fin(ev.fiedler_length),
        "d_vol": _fin(geo.d_eff.value),
        "d_vol_status": geo.d_eff.status,
        "d_s": _fin(geo.d_s.value),
        "d_s_status": geo.d_s.status,
        "d_weyl": _fin(ev.d_weyl.value),
        "d_weyl_status": ev.d_weyl.status,
        "d_weyl_r2_plateau": bool(ev.d_weyl.plateau),
        "consensus_dimension": _fin(ev.consensus_dimension),
        "dimension_class": ev.dimension_class,
        "d_hop": _fin(hop.value) if hop is not None else None,
        "distance_sensitive": bool(ev.suite.distance_sensitive),
        "metric_spread": _fin(ev.suite.metric_spread),
        "resistance_exponent": _fin(ev.suite.resistance_exponent),
        "homogeneity_cv": _fin(ev.homogeneity.cv),
        "isotropy_median": _fin(ev.isotropy.median_ratio),
        "isotropy_p10": _fin(ev.isotropy.p10_ratio),
        "detour_fraction": _fin(ev.locality.detour_fraction),
        "annulus_ok": bool(ev.annulus.ok),
        "topology_stable": bool(ev.topo.stable),
        "b1_density_w_min": _fin(ev.topo.at_w_min.b1_density),
        "kappa_mean": _fin(ev.curvature.mean),
        "kappa_tail_fraction": _fin(ev.curvature.tail_fraction),
        "kappa_n_edges": ev.curvature.n_edges,
    }
