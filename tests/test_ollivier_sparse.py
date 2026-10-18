"""Equivalencia dispersa/densa de la curvatura de Ollivier (auditoria K)."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse
from scipy.spatial import cKDTree

from omega.c0.references import rng_from_key
from omega.curvature.ollivier import ollivier_edge
from omega.curvature.ollivier_sparse import ollivier_edge_sparse
from omega.geometry.distances import hop_distance_matrix


def _csr(n: int, u: np.ndarray, v: np.ndarray) -> sparse.csr_array:
    keep = u != v
    u, v = u[keep], v[keep]
    m = sparse.coo_array((np.ones(2 * u.size), (np.concatenate([u, v]), np.concatenate([v, u]))), shape=(n, n))
    a = sparse.csr_array(m)
    a.data[:] = 1.0
    return a


def _torus() -> sparse.csr_array:
    s = 5
    idx = np.arange(s**3)
    c = np.stack(np.unravel_index(idx, (s,) * 3), axis=1)
    us, vs = [], []
    for ax in range(3):
        cc = c.copy()
        cc[:, ax] = (cc[:, ax] + 1) % s
        us.append(idx)
        vs.append(np.ravel_multi_index(tuple(cc.T), (s,) * 3))
    return _csr(s**3, np.concatenate(us), np.concatenate(vs))


def _rgg() -> sparse.csr_array:
    rng = rng_from_key((7, 1))
    x = rng.random((300, 2))
    p = cKDTree(x).query_pairs(0.09, output_type="ndarray")
    return _csr(300, p[:, 0], p[:, 1])


def _tree() -> sparse.csr_array:
    rng = rng_from_key((7, 2))
    v = np.arange(1, 200)
    u = np.asarray([rng.integers(0, i) for i in v])
    return _csr(200, u, v)


def _caveman() -> sparse.csr_array:
    size, nb = 6, 10
    us, vs = [], []
    for b in range(nb):
        iu, ju = np.triu_indices(size, 1)
        us.append(iu + b * size)
        vs.append(ju + b * size)
    u, v = np.concatenate(us), np.concatenate(vs)
    starts = np.arange(nb) * size
    drop = np.isin(u * 1000 + v, starts * 1000 + starts + 1)
    u, v = u[~drop], v[~drop]
    return _csr(size * nb, np.concatenate([u, starts]), np.concatenate([v, (starts - 1) % (size * nb)]))


@pytest.mark.parametrize("make", [_torus, _rgg, _tree, _caveman])
def test_sparse_equals_dense(make: object) -> None:
    adj = make()  # type: ignore[operator]
    a = adj.toarray() != 0
    d = hop_distance_matrix(a)
    w = a.astype(np.float64)
    xs, ys = np.nonzero(np.triu(a, 1))
    rng = rng_from_key((7, 3))
    pick = rng.choice(xs.size, size=min(25, xs.size), replace=False)
    for i in pick:
        x, y = int(xs[i]), int(ys[i])
        kd = ollivier_edge(w, a, d, x, y, 0.5)
        ks = ollivier_edge_sparse(adj, x, y, 0.5)
        assert abs(kd - ks) < 1e-9


def test_rejects_non_edge() -> None:
    adj = _tree()
    with pytest.raises(ValueError):
        ollivier_edge_sparse(adj, 0, 0, 0.5)
