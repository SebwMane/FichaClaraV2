"""omega.phases.finite_size: size_config, fss_observables y fss_fit con controles de pendiente conocida."""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest

from omega.certificate.evidence import collect_run_evidence
from omega.config.seeds import SeedKey, make_rng
from omega.config.settings11 import Omega11Config
from omega.experiments.reference_graphs import periodic_lattice, random_regular
from omega.geometry.weyl import fiedler_length
from omega.phases.finite_size import DENSE_CURVATURE_EDGES, FSS_KEYS, evidence_cfg, fss_fit, fss_observables, size_config
from omega.phases.scan import default_config
from omega.types import RunStatus
from tests.factories_v11 import make_evidence

KEY = SeedKey(20240901, (11, 0, 0))


def _cfg(n: int = 100) -> Omega11Config:
    return Omega11Config(base=default_config(n, master_entropy=20240901, replicates=3))


def _xi(shape: tuple[int, ...]) -> tuple[int, float]:
    w = periodic_lattice(shape)
    ev = collect_run_evidence(w, RunStatus.CONVERGED, _cfg(w.shape[0]), make_rng(KEY))
    return int(w.shape[0]), fss_observables(ev)["xi_f"]


def test_size_config_only_changes_n() -> None:
    cfg = _cfg(100)
    c2 = size_config(cfg, 300)
    assert c2.base.init.n == 300 and cfg.base.init.n == 100
    assert dataclasses.replace(c2, base=cfg.base) == cfg  # nada mas cambia (umbrales intactos)
    with pytest.raises(ValueError):
        size_config(cfg, 2)
    with pytest.raises(TypeError):
        size_config(cfg, 10.0)  # type: ignore[arg-type]


def test_fss_observables_keys_and_values() -> None:
    ev = make_evidence()
    o = fss_observables(ev)
    assert tuple(o) == FSS_KEYS
    assert o["d_star"] == 3.0 and o["d_vol"] == 3.0 and o["d_s"] == 3.0 and o["d_weyl"] == 3.0
    assert o["c_bin"] == 0.4 and o["giant_fraction"] == 1.0 and o["rho"] == 0.015
    assert o["xi_f"] == 4.0 and o["l_hop"] == 5.0
    ev_nan = dataclasses.replace(ev, consensus_dimension=math.nan)
    assert math.isnan(fss_observables(ev_nan)["d_star"])


def test_fss_fit_recovers_known_power_law() -> None:
    n = [64, 100, 200, 400, 800]
    v = [2.5 * x**0.5 for x in n]
    f = fss_fit(n, v)
    assert f["slope_loglog"] == pytest.approx(0.5, abs=1e-12)
    assert f["intercept_loglog"] == pytest.approx(math.log(2.5), abs=1e-12)
    assert f["r2_loglog"] == pytest.approx(1.0, abs=1e-12)
    assert f["se_loglog"] == pytest.approx(0.0, abs=1e-9)
    assert f["n_points"] == 5.0


def test_fss_fit_recovers_log_law_and_noise() -> None:
    n = [64, 100, 150, 200, 300, 500, 800]
    f = fss_fit(n, [1.0 + 0.7 * math.log(x) for x in n])
    assert f["slope_ln"] == pytest.approx(0.7, abs=1e-12)
    assert f["intercept_ln"] == pytest.approx(1.0, abs=1e-12)
    rng = np.random.Generator(np.random.PCG64(7))
    noisy = [x**-0.3 * math.exp(0.01 * float(rng.standard_normal())) for x in n]
    g = fss_fit(n, noisy)
    assert abs(g["slope_loglog"] + 0.3) < 3 * g["se_loglog"] + 0.01
    assert g["se_loglog"] > 0.0


def test_fss_fit_handles_nan_nonpositive_and_errors() -> None:
    f = fss_fit([10, 20, 40, 80], [1.0, math.nan, 4.0, 8.0])
    assert f["n_points"] == 3.0 and f["slope_loglog"] == pytest.approx(1.0, abs=1e-9)
    z = fss_fit([10, 20, 40], [0.0, 0.0, 0.0])  # v <= 0: sin ajuste log-log, constante en ln N
    assert math.isnan(z["slope_loglog"]) and z["slope_ln"] == pytest.approx(0.0, abs=1e-12)
    two = fss_fit([10, 20], [1.0, 2.0])
    assert two["slope_loglog"] == pytest.approx(1.0) and math.isnan(two["se_loglog"])
    with pytest.raises(ValueError):
        fss_fit([10, 20, 30], [1.0, 2.0])
    with pytest.raises(ValueError):
        fss_fit([10], [1.0])
    with pytest.raises(ValueError):
        fss_fit([10, 0], [1.0, 2.0])


def test_xi_f_scales_as_n_to_one_over_d_on_tori() -> None:
    # anillo (D=1): xi_F = N / (2 pi) exacto en pendiente
    ns, xs = zip(*[_xi((n,)) for n in (50, 100, 200)], strict=True)
    assert fss_fit(ns, xs)["slope_loglog"] == pytest.approx(1.0, abs=0.02)
    # toro 2D (D=2)
    ns, xs = zip(*[_xi((el, el)) for el in (8, 12, 16)], strict=True)
    assert fss_fit(ns, xs)["slope_loglog"] == pytest.approx(0.5, abs=0.03)
    # toro 3D (D=3): mismo observable (fiedler_length) sin la evidencia completa, por coste
    pts = [(el**3, fiedler_length(periodic_lattice((el, el, el)), 0.5)) for el in (6, 8, 10, 12)]
    assert fss_fit([p[0] for p in pts], [p[1] for p in pts])["slope_loglog"] == pytest.approx(1.0 / 3.0, abs=0.04)


def test_xi_f_roughly_constant_on_expanders() -> None:
    ns = [100, 200, 400, 800]
    xs = []
    for n in ns:
        w = random_regular(n, 6, make_rng(SeedKey(5, (n,))))
        xs.append(fiedler_length(w, 0.5))
    assert abs(fss_fit(ns, xs)["slope_loglog"]) < 0.05


def test_evidence_cfg_budgets_only_dense_states() -> None:
    cfg = _cfg(60)
    sparse = periodic_lattice((60,))
    dense = np.ones((60, 60)) - np.eye(60)
    assert evidence_cfg(cfg, sparse) is cfg
    c2 = evidence_cfg(cfg, dense)
    assert c2.curvature.max_edges == DENSE_CURVATURE_EDGES < cfg.curvature.max_edges
    assert dataclasses.replace(c2, curvature=cfg.curvature) == cfg  # ningun umbral cambia
    with pytest.raises(TypeError):
        evidence_cfg("cfg", dense)  # type: ignore[arg-type]


def test_evidence_cfg_keeps_full_budget_in_intermediate_density() -> None:
    """Enmienda A-6 (auditoria B2): 0.1 < rho < 0.5 no es F1 garantizado, no se recorta."""
    cfg = _cfg(60)
    ii, jj = np.indices((60, 60))
    a = np.triu(((ii * 7 + jj * 3) % 10 < 3).astype(float), 1)  # densidad ~0.3 determinista (sin RNG)
    mid = a + a.T
    assert 0.1 < float(np.mean(mid[np.triu_indices(60, 1)] > 0)) < cfg.certificate.dense_rho
    assert evidence_cfg(cfg, mid) is cfg
    # umbral inclusivo (>=): rho exactamente dense_rho si se recorta
    half = np.zeros((60, 60))
    iu = np.triu_indices(60, 1)
    idx = np.arange(iu[0].shape[0])
    sel = idx[: idx.shape[0] // 2]
    half[iu[0][sel], iu[1][sel]] = 1.0
    half = half + half.T
    assert evidence_cfg(cfg, half).curvature.max_edges == DENSE_CURVATURE_EDGES
