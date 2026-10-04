"""Persistencia H0 exacta, curvas de Betti y estabilidad topologica (OMEGA_1_1_DESIGN §1.2, §3.3).

H0 por Kruskal sobre aristas W>0 ordenadas por (-W, i, j): barras (nacimiento 1.0, muerte W_e) y
una barra esencial por componente de {W>0}. No hay codigo de barras H1 (coste no justificado a
este N). Referencia: Edelsbrunner-Letscher-Zomorodian (2002), persistencia topologica.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from omega.config.settings11 import TopologyConfig
from omega.contracts import BettiCurves, PersistenceH0, TopologySummary
from omega.network.weights import validate_weight_matrix
from omega.topology.betti import clique_complex_betti, short_cycle_betti1
from omega.topology.filtration import component_labels, descending_edges, threshold_graph
from omega.types import FloatArray

__all__ = ["h0_persistence", "betti_curves", "topology_stable", "topology_summary"]


def h0_persistence(w: FloatArray) -> PersistenceH0:
    """Barras H0 exactas (arbol generador maximo). #barras finitas = N - #componentes."""
    validate_weight_matrix(w)
    n = w.shape[0]
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    ii, jj, vv = descending_edges(w)
    deaths: list[float] = []
    for i, j, v in zip(ii.tolist(), jj.tolist(), vv.tolist(), strict=True):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
            deaths.append(v)
    d = np.asarray(deaths, dtype=np.float64)
    return PersistenceH0(births=np.ones_like(d), deaths=d, n_essential=n - len(deaths))


def _level_stats(w: FloatArray, theta: float, cfg: TopologyConfig) -> tuple[int, int, float, float, str]:
    a = threshold_graph(w, theta)
    n_comp, labels = component_labels(a)
    giant = float(np.bincount(labels).max()) / w.shape[0]
    sc = short_cycle_betti1(a, cfg.short_cycle_length, cfg.max_faces)
    return n_comp, sc.b1, sc.b1_density, giant, sc.status


def betti_curves(w: FloatArray, cfg: TopologyConfig) -> BettiCurves:
    """beta0(theta), beta1^(L)(theta), b1/E y fraccion gigante sobre la malla cfg.thetas."""
    validate_weight_matrix(w)
    rows = [_level_stats(w, t, cfg) for t in cfg.thetas]
    return BettiCurves(
        thetas=np.asarray(cfg.thetas, dtype=np.float64),
        beta0=np.asarray([r[0] for r in rows], dtype=np.int64),
        beta1_short=np.asarray([r[1] for r in rows], dtype=np.int64),
        b1_density=np.asarray([r[2] for r in rows], dtype=np.float64),
        giant_fraction=np.asarray([r[3] for r in rows], dtype=np.float64),
        status=tuple(r[4] for r in rows),
    )


def topology_stable(
    w: FloatArray, thetas: Sequence[float], g_connected: float, cfg: TopologyConfig
) -> tuple[bool, tuple[bool, ...]]:
    """Estabilidad topologica sobre una malla de umbrales (§1.12).

    s(theta) = (giant_fraction >= g_connected) y (b1_density <= cfg.b1_density_max) y status "ok".
    Estable si s(theta) se cumple en >= cfg.stable_fraction de los theta. Malla vacia: False.
    """
    validate_weight_matrix(w)
    flags: list[bool] = []
    for t in thetas:
        _, _, dens, giant, status = _level_stats(w, float(t), cfg)
        flags.append(bool(giant >= g_connected and dens <= cfg.b1_density_max and status == "ok"))
    stable = bool(flags) and sum(flags) / len(flags) >= cfg.stable_fraction
    return stable, tuple(flags)


def topology_summary(
    w: FloatArray, w_min: float, sensitivity: Sequence[float], g_connected: float, cfg: TopologyConfig
) -> TopologySummary:
    """Resumen topologico: curvas, H0, beta1^(L) y complejo de cliques en w_min, y estabilidad."""
    curves = betti_curves(w, cfg)
    h0 = h0_persistence(w)
    a = threshold_graph(w, w_min)
    at_w_min = short_cycle_betti1(a, cfg.short_cycle_length, cfg.max_faces)
    clique = clique_complex_betti(a, cfg.clique_max_dim, cfg.max_simplices)
    stable, flags = topology_stable(w, sensitivity, g_connected, cfg)
    return TopologySummary(curves=curves, h0=h0, at_w_min=at_w_min, clique_at_w_min=clique,
                           stable_flags=flags, stable=stable)
