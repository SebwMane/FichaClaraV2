"""Tests de omega.topology.persistence (WP-C)."""

from __future__ import annotations

import networkx as nx
import numpy as np

from omega.config.settings11 import TopologyConfig
from omega.experiments.reference_graphs import complete_graph, periodic_lattice, random_geometric_torus
from omega.topology.persistence import betti_curves, h0_persistence, topology_stable, topology_summary


def _random_w(n: int, p: float, seed: int) -> np.ndarray:
    rng = np.random.Generator(np.random.PCG64(seed))
    m = np.triu(rng.random((n, n)) * (rng.random((n, n)) < p), 1)
    return np.asarray(m + m.T, dtype=np.float64)


def test_h0_bars_match_maximum_spanning_tree() -> None:
    w = _random_w(30, 0.4, 1)
    g = nx.from_numpy_array(w)
    assert nx.is_connected(g)
    h = h0_persistence(w)
    assert h.deaths.size == 29 and h.n_essential == 1
    assert np.all(h.births == 1.0)
    mst = nx.maximum_spanning_tree(g)
    expected = sorted(d["weight"] for _, _, d in mst.edges(data=True))
    assert np.allclose(sorted(h.deaths.tolist()), expected)


def test_h0_disconnected_essential() -> None:
    w = np.zeros((6, 6))
    w[0, 1] = w[1, 0] = 0.5
    w[2, 3] = w[3, 2] = 0.7
    h = h0_persistence(w)
    assert h.deaths.size == 2 and h.n_essential == 4


def test_betti_curves_monotone() -> None:
    w = _random_w(40, 0.15, 2)
    cfg = TopologyConfig()
    c = betti_curves(w, cfg)
    assert c.thetas.shape == (len(cfg.thetas),)
    assert np.all(np.diff(c.beta0) <= 0)  # theta desciende: beta0 no crece
    assert np.all(np.diff(c.giant_fraction) >= -1e-15)
    assert all(s == "ok" for s in c.status)


def test_topology_stable_lattice_and_dense_random() -> None:
    cfg = TopologyConfig()
    ok, flags = topology_stable(periodic_lattice((6, 6)), (0.01, 0.05, 0.1, 0.2, 0.5), 0.95, cfg)
    assert ok and all(flags)
    ok, flags = topology_stable(_random_w(60, 0.1, 5), (0.01, 0.05, 0.1), 0.95, cfg)
    assert not ok and not any(flags)  # b1/E por encima de 0.03
    assert topology_stable(complete_graph(5), (), 0.95, cfg) == (False, ())


def test_topology_summary_torus() -> None:
    cfg = TopologyConfig()
    s = topology_summary(periodic_lattice((6, 6)), 0.1, (0.1, 0.5), 0.95, cfg)
    assert s.at_w_min.b1 == 2 and s.stable
    assert s.clique_at_w_min.betti[0] == 1 and s.h0.n_essential == 1


def test_curves_over_budget_flagged() -> None:
    cfg = TopologyConfig(thetas=(0.5,), max_faces=10)
    c = betti_curves(complete_graph(12), cfg)
    assert c.status == ("over_budget",)
    assert topology_stable(complete_graph(12), (0.5,), 0.95, cfg)[0] is False


def test_rgg_curves_run() -> None:
    w = random_geometric_torus(80, 2, 8, np.random.Generator(np.random.PCG64(0)))
    c = betti_curves(w, TopologyConfig(thetas=(0.5, 0.1)))
    assert c.beta0.shape == (2,)
