"""Tests de omega.curvature.qrc (WP-C)."""

from __future__ import annotations

import numpy as np

from omega.experiments.reference_graphs import periodic_lattice
from omega.curvature.qrc import average_sphere_distance, quantum_ricci_profile
from omega.geometry.distances import hop_distance_matrix


def test_ring_by_hand() -> None:
    # C_12, x=0, y=2, delta=2: S(0)={2,10}, S(2)={0,4}; d: (2,0)=2,(2,4)=2,(10,0)=2,(10,4)=6 -> media 3
    d = hop_distance_matrix(periodic_lattice((12,)) > 0)
    assert average_sphere_distance(d, 0, 2, 2) == 3.0
    # delta=1, y=1: S(0)={1,11}, S(1)={0,2}: d=(1,0)=1,(1,2)=1,(11,0)=1,(11,2)=3 -> 1.5
    assert average_sphere_distance(d, 0, 1, 1) == 1.5


def test_empty_sphere_nan() -> None:
    d = hop_distance_matrix(periodic_lattice((6,)) > 0)
    assert np.isnan(average_sphere_distance(d, 0, 1, 5))


def test_profile_ring_and_torus_artifacts() -> None:
    rng = np.random.Generator(np.random.PCG64(0))
    d = hop_distance_matrix(periodic_lattice((12,)) > 0)
    p = quantum_ricci_profile(d, (1, 2), 200, rng)
    assert p.ratio_mean[0] == 1.5 and p.ratio_mean[1] == 1.5  # d_bar(delta=2)/2 = 3/2
    assert p.n_pairs[0] == 12 and p.ratio_se[0] == 0.0
    d2 = hop_distance_matrix(periodic_lattice((12, 12)) > 0)
    p2 = quantum_ricci_profile(d2, (1, 2, 3), 200, np.random.Generator(np.random.PCG64(1)))
    assert p2.ratio_mean[0] == 1.875  # artefacto de reticula (hallazgo 6: ~1.88)
    assert np.all(p2.n_pairs == np.array([200, 200, 200]) ) or p2.n_pairs[0] == 200


def test_profile_deterministic_and_sampled() -> None:
    d = hop_distance_matrix(periodic_lattice((10, 10)) > 0)
    a = quantum_ricci_profile(d, (2,), 20, np.random.Generator(np.random.PCG64(7)))
    b = quantum_ricci_profile(d, (2,), 20, np.random.Generator(np.random.PCG64(7)))
    assert a.n_pairs[0] == 20 and np.array_equal(a.ratio_mean, b.ratio_mean)
