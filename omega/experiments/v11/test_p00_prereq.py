"""Pasos 1-2 de P§22 (Ω-1.1): suite analitica y golden de Ω-1.0 como prerrequisitos de la compuerta (auditoria B7).

Enmienda A-13 (auditoria B7): `STEP_ORDER` empieza con `p01_analytical` y `p02_golden`. `run()` ejecuta UNA vez
`pytest.main(["-q", "-p", "no:cacheprovider", "-m", "slow or not slow", tests/analytical, tests/test_baseline_golden.py])`
(incluidos los tests lentos) y escribe ambos `summary.json` con `complete = (rc == 0)`. En modo full exige un arbol de git
limpio (el summary registra `code_commit`, `git_dirty` y `config_hash`, que `require_prerequisites` verifica).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal

import pytest

from omega.config.settings11 import Omega11Config
from omega.experiments.v11.gate import (
    PrerequisiteError,
    require_clean_tree,
    require_prerequisites,
    write_summary,
)
from omega.phases.scan import default_config

STEP_ANALYTICAL = "p01_analytical"
STEP_GOLDEN = "p02_golden"
TARGETS = {STEP_ANALYTICAL: "tests/analytical", STEP_GOLDEN: "tests/test_baseline_golden.py"}
REPO_ROOT = Path(__file__).resolve().parents[3]


def pytest_args() -> list[str]:
    """Argumentos exactos de `pytest.main` (incluye los tests `slow`)."""
    return ["-q", "-p", "no:cacheprovider", "-m", "slow or not slow", *(str(REPO_ROOT / t) for t in TARGETS.values())]


def run(cfg: Omega11Config, out_root: Path, *, mode: Literal["smoke", "full"]) -> Path:
    """Ejecuta pasos 1-2 y devuelve el `summary.json` de `p02_golden`. En full rechaza un arbol sucio."""
    if mode == "full":
        require_clean_tree(STEP_ANALYTICAL)
    args = pytest_args()
    rc = int(pytest.main(args))
    last = out_root
    for step, target in TARGETS.items():
        data: dict[str, Any] = {
            "description": f"prerrequisito P§22: {target}", "target": target, "pytest_args": args,
            "returncode": rc, "complete": rc == 0, "expectations_met": rc == 0,
        }
        last = write_summary(out_root, step, data, mode=mode, cfg=cfg)
    return last


def test_pytest_args_cover_both_targets_with_slow() -> None:
    args = pytest_args()
    assert args[:5] == ["-q", "-p", "no:cacheprovider", "-m", "slow or not slow"]
    assert args[5].endswith("tests/analytical") and args[6].endswith("tests/test_baseline_golden.py")
    assert Path(args[5]).is_dir() and Path(args[6]).is_file()


def test_smoke_writes_both_summaries(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[list[str]] = []

    def fake_main(args: list[str]) -> int:
        calls.append(args)
        return 0

    monkeypatch.setattr(pytest, "main", fake_main)
    monkeypatch.setattr("omega.experiments.v11.gate.head_commit", lambda: "abc")
    monkeypatch.setattr("omega.experiments.v11.gate.tree_dirty", lambda: False)
    cfg = Omega11Config(base=default_config(24, master_entropy=20240901, replicates=2))
    run(cfg, tmp_path, mode="smoke")
    assert len(calls) == 1 and calls[0] == pytest_args()
    for step in TARGETS:
        d = json.loads((tmp_path / step / "summary.json").read_text(encoding="utf-8"))
        assert d["complete"] is True and d["mode"] == "smoke" and d["step"] == step
        assert d["code_commit"] == "abc" and d["git_dirty"] is False and len(d["config_hash"]) == 64


def test_failure_marks_incomplete(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(pytest, "main", lambda args: 1)
    monkeypatch.setattr("omega.experiments.v11.gate.tree_dirty", lambda: False)
    cfg = Omega11Config(base=default_config(24, master_entropy=20240901, replicates=2))
    run(cfg, tmp_path, mode="full")
    for step in TARGETS:
        d = json.loads((tmp_path / step / "summary.json").read_text(encoding="utf-8"))
        assert d["complete"] is False and d["returncode"] == 1


def test_full_rejects_dirty_tree(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    called: list[int] = []
    def fake_main(args: list[str]) -> int:
        called.append(1)
        return 0

    monkeypatch.setattr(pytest, "main", fake_main)
    monkeypatch.setattr("omega.experiments.v11.gate.tree_dirty", lambda: True)
    cfg = Omega11Config(base=default_config(24, master_entropy=20240901, replicates=2))
    with pytest.raises(PrerequisiteError):
        run(cfg, tmp_path, mode="full")
    assert not called and not (tmp_path / STEP_ANALYTICAL).exists()
    with pytest.raises(PrerequisiteError):
        require_prerequisites(tmp_path, "o00_validation")


def test_gate_requires_commit_and_clean_tree(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("omega.experiments.v11.gate.head_commit", lambda: "c1")
    monkeypatch.setattr("omega.experiments.v11.gate.tree_dirty", lambda: False)
    for step in TARGETS:
        write_summary(tmp_path, step, {"complete": True}, mode="full")
    require_prerequisites(tmp_path, "o00_validation")
    monkeypatch.setattr("omega.experiments.v11.gate.head_commit", lambda: "c2")  # otro commit
    with pytest.raises(PrerequisiteError, match="otro commit"):
        require_prerequisites(tmp_path, "o00_validation")
    monkeypatch.setattr("omega.experiments.v11.gate.head_commit", lambda: "c1")
    monkeypatch.setattr("omega.experiments.v11.gate.tree_dirty", lambda: True)  # summary escrito con arbol sucio
    write_summary(tmp_path, STEP_GOLDEN, {"complete": True}, mode="full")
    monkeypatch.setattr("omega.experiments.v11.gate.tree_dirty", lambda: False)
    with pytest.raises(PrerequisiteError, match="sucio"):
        require_prerequisites(tmp_path, "o00_validation")
