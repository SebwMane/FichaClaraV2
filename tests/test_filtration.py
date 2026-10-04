"""Tests de omega.topology.filtration (WP-C)."""

from __future__ import annotations

import numpy as np
import pytest

from omega.experiments.reference_graphs import complete_graph, periodic_lattice
from omega.topology.filtration import descending_edges, filtration_levels, threshold_graph


def test_threshold_graph_strict_and_symmetric() -> None:
    w = np.array([[0.0, 0.5, 0.2], [0.5, 0.0, 0.5], [0.2, 0.5, 0.0]])
    a = threshold_graph(w, 0.5)
    assert not a.any()  # umbral estricto
    a = threshold_graph(w, 0.2)
    assert a[0, 1] and a[1, 2] and not a[0, 2]
    assert np.array_equal(a, a.T) and not a.diagonal().any()


def test_threshold_graph_validates() -> None:
    with pytest.raises(ValueError):
        threshold_graph(complete_graph(3), 1.5)
    with pytest.raises(ValueError):
        threshold_graph(np.array([[0.0, 1.0], [0.5, 0.0]]), 0.1)


def test_descending_edges_order_and_ties() -> None:
    w = np.zeros((4, 4))
    for i, j, v in [(0, 1, 0.5), (2, 3, 0.9), (0, 2, 0.5), (1, 3, 0.1)]:
        w[i, j] = w[j, i] = v
    ei, ej, ev = descending_edges(w)
    assert list(zip(ei.tolist(), ej.tolist(), ev.tolist())) == [(2, 3, 0.9), (0, 1, 0.5), (0, 2, 0.5), (1, 3, 0.1)]


def test_filtration_levels_ring() -> None:
    w = periodic_lattice((10,))
    lv = filtration_levels(w, (0.5, 0.0))
    assert lv[0].n_edges == 10 and lv[0].n_components == 1 and lv[0].giant_fraction == 1.0
    assert lv[0].mean_degree == 2.0
    w2 = w.copy()
    w2[0, 1] = w2[1, 0] = 0.3
    lv2 = filtration_levels(w2, (0.5,))
    assert lv2[0].n_edges == 9 and lv2[0].n_components == 1


def test_filtration_levels_isolated_nodes() -> None:
    w = np.zeros((5, 5))
    w[0, 1] = w[1, 0] = 0.8
    lv = filtration_levels(w, (0.5,))[0]
    assert lv.n_components == 4 and lv.giant_fraction == pytest.approx(0.4)
