"""Pruebas de io/provenance: pasaporte Ω-1.1, IDs atómicos y reconstrucción (WP-E)."""

from __future__ import annotations

import json
import threading
from pathlib import Path

import numpy as np
import pytest

from omega.baseline import BASELINE_COMMIT
from omega.config.convert11 import config_hash
from omega.config.seeds import SeedKey, make_rng
from omega.config.settings11 import NullModel, Omega11Config
from omega.dynamics.evolution import evolve
from omega.io.passport import PassportIntegrityError, array_digest
from omega.io.provenance import (
    allocate_experiment_id,
    build_passport_v11,
    dependency_lock,
    dependency_lock_hash,
    experiment_file_stem,
    format_experiment_id,
    git_info,
    load_passport_v11,
    parse_experiment_id,
    reconstruct,
    save_passport_v11,
)
from omega.network.initialization import random_uniform_weights
from omega.network.weights import upper_triangle
from tests.conftest import make_config

FIELDS = (
    "model_version", "baseline", "code_commit", "git_branch", "python_version", "dependency_lock_hash",
    "config_hash", "initialization_distribution", "distance_definition", "threshold_rule",
    "dimension_algorithm", "coarse_graining_algorithm", "null_model", "random_seed",
    "termination_reason", "experiment_id", "reconstruct", "omega10",
)


def _simulate(cfg: Omega11Config, seed: SeedKey) -> tuple[np.ndarray, np.ndarray]:
    w0 = random_uniform_weights(cfg.base.init.n, make_rng(seed))
    traj = evolve(w0, cfg.base.functional, cfg.base.dynamics)
    return w0, traj.w_final


def _passport(cfg: Omega11Config, seed: SeedKey, w0: np.ndarray, number: int = 4821) -> dict:  # type: ignore[type-arg]
    return build_passport_v11(
        cfg, experiment_number=number, seed=seed, termination_reason="CONVERGED",
        initialization_distribution="uniform/upper_mirror", null_model=NullModel.ERDOS_RENYI,
        coarse_graining=True, results={"x": float("nan"), "y": np.float64(1.5)},
        code_commit="abc", git_branch="main", entrypoint="omega.dynamics.evolution:evolve",
        w0_sha256=array_digest(upper_triangle(w0)), omega10_passport=None,
    )


def test_fields_present() -> None:
    cfg = Omega11Config(base=make_config(12))
    seed = SeedKey(1, (2, 0, 0))
    w0, _ = _simulate(cfg, seed)
    p = _passport(cfg, seed, w0)
    for f in FIELDS:
        assert f in p, f
    assert p["model_version"] == "Ω-1.1"
    assert p["baseline"]["commit"] == BASELINE_COMMIT
    assert p["experiment_id"] == "Ω-EXP-004821"
    assert p["coarse_graining_algorithm"] == "heavy_edge_matching/max/v1"
    assert p["null_model"] == "erdos_renyi"
    assert p["config_hash"] == config_hash(cfg)
    assert p["results"]["x"] is None
    assert set(p["random_seed"]) == {"entropy", "spawn_key", "bit_generator"}
    assert set(p["reconstruct"]) == {"entrypoint", "config", "seed", "w0_sha256"}
    with pytest.raises(ValueError):
        build_passport_v11(
            cfg, experiment_number=1, seed=seed, termination_reason="BAD",
            initialization_distribution="u", null_model=None, coarse_graining=False, results={},
            code_commit="a", git_branch="b", entrypoint="m:f", w0_sha256=None, omega10_passport=None,
        )


def test_dependency_lock_hash() -> None:
    lock = dependency_lock()
    assert {"python", "numpy", "scipy", "networkx"} <= set(lock)
    assert dependency_lock_hash(lock) == dependency_lock_hash(dict(reversed(list(lock.items()))))
    assert dependency_lock_hash({"a": "1"}) != dependency_lock_hash({"a": "2"})


def test_git_info_unknown_on_failure(tmp_path: Path) -> None:
    assert git_info(tmp_path / "no_existe") == ("unknown", "unknown")
    c, b = git_info(Path(__file__).resolve().parent)
    assert isinstance(c, str) and isinstance(b, str)


def test_id_format_parse() -> None:
    assert format_experiment_id(4821) == "Ω-EXP-004821"
    assert experiment_file_stem(4821) == "OMEGA-EXP-004821"
    assert parse_experiment_id("Ω-EXP-004821") == 4821
    assert experiment_file_stem(4821).isascii()
    with pytest.raises(ValueError):
        parse_experiment_id("EXP-1")
    with pytest.raises(ValueError):
        format_experiment_id(-1)


def test_allocate_sequential_and_concurrent(tmp_path: Path) -> None:
    assert allocate_experiment_id(tmp_path) == 1
    assert allocate_experiment_id(tmp_path) == 2
    got: list[int] = []
    lock = threading.Lock()

    def work() -> None:
        for _ in range(5):
            n = allocate_experiment_id(tmp_path)
            with lock:
                got.append(n)

    threads = [threading.Thread(target=work) for _ in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert len(got) == 40 and len(set(got)) == 40
    assert len(list((tmp_path / "ids").glob("*.claim"))) == 42


def test_save_load_tamper_and_reconstruct(tmp_path: Path) -> None:
    cfg = Omega11Config(base=make_config(12))
    seed = SeedKey(cfg.base.seeds.master_entropy, (cfg.base.seeds.experiment_id, 0, 0))
    w0, wf = _simulate(cfg, seed)
    arrays = {"w0_upper": upper_triangle(w0), "w_final_upper": upper_triangle(wf)}
    p = _passport(cfg, seed, w0, number=7)
    jp, npz = save_passport_v11(p, arrays, tmp_path)
    assert jp.name == "OMEGA-EXP-000007.json" and npz.name == "OMEGA-EXP-000007.npz"
    p2, arrs = load_passport_v11(jp)
    assert p2["experiment_id"] == "Ω-EXP-000007"
    cfg2, seed2 = reconstruct(p2)
    assert cfg2 == cfg and seed2 == seed
    w0b, wfb = _simulate(cfg2, seed2)
    assert array_digest(upper_triangle(w0b)) == p2["array_digests"]["w0_upper"] == p2["reconstruct"]["w0_sha256"]
    assert array_digest(upper_triangle(wfb)) == p2["array_digests"]["w_final_upper"]
    assert set(arrs) == set(arrays)
    # digest alterado en el JSON
    bad = json.loads(jp.read_text(encoding="utf-8"))
    bad["array_digests"]["w0_upper"] = "0" * 64
    jp_bad = tmp_path / "bad.json"
    jp_bad.write_text(json.dumps(bad, ensure_ascii=False), encoding="utf-8")
    with pytest.raises(PassportIntegrityError):
        load_passport_v11(jp_bad)
    # NPZ alterado
    npz.write_bytes(npz.read_bytes() + b"x")
    with pytest.raises(PassportIntegrityError):
        load_passport_v11(jp)
    # config alterada
    bad2 = dict(p2)
    rec = dict(bad2["reconstruct"])
    cfgd = json.loads(json.dumps(rec["config"]))
    cfgd["topology"]["b1_density_max"] = 0.04
    rec["config"] = cfgd
    bad2["reconstruct"] = rec
    with pytest.raises(PassportIntegrityError):
        reconstruct(bad2)
