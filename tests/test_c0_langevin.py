import numpy as np

from omega.c0.functional import action_c0, params_from_targets
from omega.c0.langevin import background_fraction, langevin_c0, reflect_unit


def _rng(s: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(s))


def _rand_w(n: int, scale: float, s: int) -> np.ndarray:
    w = np.triu(_rng(s).random((n, n)) * scale, 1)
    return w + w.T


def test_reflect_unit() -> None:
    x = np.array([-0.2, 0.3, 1.25, 2.1, -1.5])
    np.testing.assert_allclose(reflect_unit(x), [0.2, 0.3, 0.75, 0.1, 0.5])


def test_background_fraction() -> None:
    w = np.zeros((3, 3))
    w[0, 1] = w[1, 0] = 1.0
    w[0, 2] = w[2, 0] = 0.05
    assert abs(background_fraction(w) - 0.1 / 2.1) < 1e-12


def test_theta_zero_is_deterministic_descent() -> None:
    p = params_from_targets(2, 8, 1.0)
    w0 = _rand_w(40, 0.4, 0)
    out = langevin_c0(w0, p, 0.0, dt=0.01, n_steps=200, rng=_rng(1), record_every=50)
    assert out["s_final"] < action_c0(w0, p)
    w = out["w"]
    assert np.allclose(w, w.T) and np.all(np.diag(w) == 0) and w.min() >= 0 and w.max() <= 1


def test_noise_reproducible_and_raises_background() -> None:
    p = params_from_targets(2, 8, 1.0)
    w0 = _rand_w(40, 0.4, 0)
    a = langevin_c0(w0, p, 0.05, dt=0.01, n_steps=100, rng=_rng(2))["w"]
    b = langevin_c0(w0, p, 0.05, dt=0.01, n_steps=100, rng=_rng(2))["w"]
    c = langevin_c0(w0, p, 0.0, dt=0.01, n_steps=100, rng=_rng(2))["w"]
    assert np.array_equal(a, b)
    assert background_fraction(a) > background_fraction(c)
