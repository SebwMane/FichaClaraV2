import numpy as np

from omega.c0 import references as R
from omega.diagnostics.ball_growth import ball_profile, giant


def test_lattice_is_cv_zero() -> None:
    p = ball_profile(R.torus_lattice_3d(7))
    assert p["q4"]["status"] == "CV_CERO"
    assert p["mean"][0] == 7.0 and p["mean"][1] == 25.0


def test_giant_component_only() -> None:
    a = np.zeros((6, 6))
    for i, j in ((0, 1), (1, 2), (2, 3), (4, 5)):
        a[i, j] = a[j, i] = 1.0
    assert giant(a).shape == (4, 4)


def test_path_has_no_window_or_profile() -> None:
    n = 200
    a = np.zeros((n, n))
    for i in range(n - 1):
        a[i, i + 1] = a[i + 1, i] = 1.0
    p = ball_profile(a)
    assert p["q4"]["status"] in ("HOMOGENEIZA", "NO_HOMOGENEIZA")
    assert abs(p["D_B"][3] - 1.0) < 0.2


def test_expander_has_no_window() -> None:
    rng = np.random.Generator(np.random.PCG64(0))
    p = ball_profile(R.random_regular(300, 12, rng))
    assert p["q4"]["status"] == "SIN_VENTANA"


def test_dimension_plateau_synthetic() -> None:
    from omega.diagnostics.ball_growth import dimension_plateau

    prof = {"q4": {"window": [2, 3, 4, 5, 6]}, "D_B": [1.0, 2.9, 3.0, 3.05, 3.0, 2.9]}
    out = dimension_plateau(prof, 0.15)
    assert out["status"] == "PLATEAU" and abs(out["D_plat"] - 3.0) < 0.1 and out["run"] == [2, 3, 4, 5]
    osc = {"q4": {"window": [2, 3, 4, 5, 6, 7]}, "D_B": [0.5, 1.66, 0.7, 1.51, 0.67, 1.46]}
    assert dimension_plateau(osc, 0.15)["status"] == "SIN_PLATEAU"
    short = {"q4": {"window": [2, 3, 4]}, "D_B": [1.0, 2.0]}
    assert dimension_plateau(short)["status"] == "SIN_VENTANA_D"


def test_dimension_plateau_lattice_and_rgg() -> None:
    from omega.controls.random_geometric import rgg_torus
    from omega.diagnostics.ball_growth import dimension_plateau

    ring = ball_profile(R.ring_lattice(600, 12))
    p = dimension_plateau(ring, 0.15)
    assert p["status"] == "PLATEAU" and abs(p["D_plat"] - 1.0) < 0.1
    g = ball_profile(rgg_torus(1500, 2, 12.0, np.random.Generator(np.random.PCG64(0))))
    assert dimension_plateau(g, 0.25)["status"] in ("PLATEAU", "SIN_VENTANA_D")
