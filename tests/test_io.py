"""Pruebas de io/passport (M§49): ida y vuelta, digests, campos obligatorios."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from omega.config.convert import config_from_dict
from omega.config.seeds import SeedKey, make_rng
from omega.config.settings import OmegaConfig
from omega.io.passport import (
    PassportIntegrityError,
    array_digest,
    build_passport,
    load_run,
    save_run,
    software_versions,
)
from omega.network.initialization import random_uniform_weights
from omega.network.topology import topology_observables
from omega.network.weights import upper_triangle
from omega.types import (
    DimensionEstimate,
    GeometryObservables,
    Observables,
    PhaseAssessment,
    PhaseLabel,
    RunResult,
    RunStatus,
    Trajectory,
)
from tests.conftest import make_config


def _dim(value: float, status: str, method: str) -> DimensionEstimate:
    x = np.linspace(1.0, 5.0, 7)
    return DimensionEstimate(
        value=value, stderr=0.01, status=status,  # type: ignore[arg-type]
        window=(1, 5) if status == "ok" else None,
        scales=x, profile=x**2, local_slopes=np.full(6, 2.0), plateau=status == "ok", method=method,
    )


def make_result(cfg: OmegaConfig) -> RunResult:
    n = cfg.init.n
    key = SeedKey(cfg.seeds.master_entropy, (cfg.seeds.experiment_id, 0, 0))
    w0 = random_uniform_weights(n, make_rng(key))
    wf = np.clip(w0 * 1.1, 0.0, 1.0)
    np.fill_diagonal(wf, 0.0)
    m = n * (n - 1) // 2
    traj = Trajectory(
        w_final=wf, status=RunStatus.CONVERGED, steps=3, dt=0.01, tau=0.03,
        scalars={"action": np.array([3.0, 2.0, 1.0, 0.5]), "max_dw": np.array([0.1, 0.01, 0.0])},
        snapshot_steps=np.array([0, 1, 2], dtype=np.int64),
        snapshots=np.vstack([upper_triangle(w0)] * 3).reshape(3, m),
    )
    geo = GeometryObservables(
        path_length=1.5, path_length_hops=1.2, diameter=3.0,
        d_eff=_dim(2.9, "ok", "shell"), d_eff_ball=_dim(float("nan"), "no_window", "ball"),
        d_eff_hops=_dim(2.8, "ok", "shell"), d_s=_dim(3.0, "ok", "lazy_walk"),
    )
    obs = Observables(topology_observables(wf, cfg.graph.w_min, cfg.phases.large_component_frac), geo)
    return RunResult(
        seed=key, params=cfg.functional, alpha_hat=1.0, gamma_hat=0.5, w0=w0, trajectory=traj,
        observables=obs, assessment=PhaseAssessment(PhaseLabel.F, {"geometric": True, "dense": False}),
    )


@pytest.fixture
def saved(tmp_path: Path, cfg: OmegaConfig) -> tuple[RunResult, Path, Path]:
    res = make_result(cfg)
    j, z = save_run(res, cfg, tmp_path, "run_001", "abc123")
    return res, j, z


def test_array_digest_properties() -> None:
    a = np.arange(6, dtype=np.float64)
    assert array_digest(a) == array_digest(a.copy())
    assert array_digest(a) != array_digest(a.reshape(2, 3))
    b = a.copy()
    b[0] += 1e-15
    assert array_digest(a) != array_digest(b)
    assert len(array_digest(a)) == 64


def test_roundtrip_identical(saved: tuple[RunResult, Path, Path], cfg: OmegaConfig) -> None:
    res, j, z = saved
    assert j.exists() and z.exists()
    passport, arrays = load_run(j)
    assert config_from_dict(passport["config"]) == cfg
    assert np.array_equal(arrays["w0_upper"], upper_triangle(res.w0))
    assert np.array_equal(arrays["w_final_upper"], upper_triangle(res.trajectory.w_final))
    assert np.array_equal(arrays["snapshots_upper"], res.trajectory.snapshots)
    assert np.array_equal(arrays["snapshot_steps"], res.trajectory.snapshot_steps)
    assert np.array_equal(arrays["scalars_action"], res.trajectory.scalars["action"])
    assert np.array_equal(arrays["deff_scales"], res.observables.geometry.d_eff.scales)
    assert np.array_equal(arrays["ds_profile"], res.observables.geometry.d_s.profile)
    assert passport["omega0"]["sha256"] == array_digest(upper_triangle(res.w0))
    # el JSON en disco es JSON estricto y coincide con lo cargado
    assert json.loads(j.read_text()) == passport


def test_required_fields_m49(saved: tuple[RunResult, Path, Path]) -> None:
    passport, _ = load_run(saved[1])
    for key in (
        "schema_version", "run_id", "created_utc", "git_commit", "versions", "config", "seed",
        "params_raw", "params_reduced", "N", "dt", "tau_max", "steps", "status", "omega0",
        "observables", "phase", "curvature", "npz_file", "npz_sha256",
    ):
        assert key in passport, key
    assert passport["seed"]["bit_generator"] == "PCG64"
    assert passport["seed"]["spawn_key"] == [7, 0, 0]
    assert passport["curvature"] is None and passport["git_commit"] == "abc123"
    assert passport["run_id"] == "run_001" and passport["status"] == "CONVERGED"
    for key in ("alpha", "beta", "gamma", "eta", "mu"):
        assert key in passport["params_raw"]
    assert set(passport["params_reduced"]) == {"alpha_hat", "gamma_hat"}
    obs = passport["observables"]
    for key in ("mean_strength", "std_strength", "clustering_weighted", "clustering_binary",
                "path_length", "giant_size", "giant_fraction", "d_eff", "d_s"):
        assert key in obs, key
    assert obs["d_eff"]["window"] == [1, 5] and obs["d_eff"]["status"] == "ok"
    assert obs["d_eff_ball"]["value"] is None  # NaN -> null (JSON estricto)
    assert passport["phase"] == {"label": "F", "flags": {"geometric": True, "dense": False}}


def test_altered_npz_raises(saved: tuple[RunResult, Path, Path]) -> None:
    _, j, z = saved
    data = bytearray(z.read_bytes())
    data[len(data) // 2] ^= 0xFF
    z.write_bytes(bytes(data))
    with pytest.raises(PassportIntegrityError):
        load_run(j)


def test_altered_array_with_consistent_file_hash_raises(
    saved: tuple[RunResult, Path, Path]
) -> None:
    _, j, z = saved
    passport = json.loads(j.read_text())
    with np.load(z) as d:
        arrays = {k: d[k] for k in d.files}
    arrays["w_final_upper"] = arrays["w_final_upper"] + 1e-9
    np.savez_compressed(z, **arrays)
    import hashlib

    passport["npz_sha256"] = hashlib.sha256(z.read_bytes()).hexdigest()  # solo el hash del archivo
    j.write_text(json.dumps(passport))
    with pytest.raises(ValueError, match="w_final_upper"):
        load_run(j)


def test_altered_omega0_digest_raises(saved: tuple[RunResult, Path, Path]) -> None:
    _, j, _ = saved
    passport = json.loads(j.read_text())
    passport["omega0"]["sha256"] = "0" * 64
    j.write_text(json.dumps(passport))
    with pytest.raises(PassportIntegrityError):
        load_run(j)


def test_missing_npz_digest_field_raises(saved: tuple[RunResult, Path, Path]) -> None:
    _, j, _ = saved
    passport = json.loads(j.read_text())
    passport["npz_sha256"] = None
    j.write_text(json.dumps(passport))
    with pytest.raises(ValueError):
        load_run(j)


def test_invalid_run_id_and_no_mutation(tmp_path: Path, cfg: OmegaConfig) -> None:
    res = make_result(cfg)
    for bad in ("../x", "a/b", "", ".hidden"):
        with pytest.raises(ValueError):
            save_run(res, cfg, tmp_path, bad, "x")
    w0 = res.w0.copy()
    save_run(res, cfg, tmp_path / "sub", "ok", "x")
    assert np.array_equal(res.w0, w0)


def test_build_passport_is_json_serializable(cfg: OmegaConfig) -> None:
    p = build_passport(cfg, make_result(cfg), software_versions(), "deadbeef")
    json.dumps(p, allow_nan=False)
    assert {"python", "numpy", "scipy", "omega"} <= set(p["versions"])
