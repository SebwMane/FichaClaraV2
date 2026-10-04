"""Golden baseline de Omega-1.0 (OMEGA_1_1_DESIGN §2): recalcula cada fixture congelado."""

from __future__ import annotations

import ast
import hashlib
import json
import re
import time
from pathlib import Path
from typing import Any

import pytest

from tools.freeze_baseline import (
    BASELINE_COMMIT,
    FIXTURE_NAMES,
    MANIFEST_NAME,
    ROOT,
    compute_fixtures,
    manifest_paths,
    sha256_file,
)

EXPECTED = ROOT / "tests" / "expected_baseline"
TOL = 1e-12
DOCSTRING_ONLY = ("omega/curvature/__init__.py", "omega/coarse_graining/__init__.py")


def _load(name: str) -> dict[str, Any]:
    doc = json.loads((EXPECTED / f"{name}.json").read_text(encoding="utf-8"))
    assert doc.pop("generated_from_commit") == BASELINE_COMMIT
    assert isinstance(doc, dict)
    return doc


def _compare(exp: Any, got: Any, where: str) -> None:
    """Exactos: int/str/bool/None/claves; floats <= 1e-12 absoluto."""
    if isinstance(exp, dict):
        assert isinstance(got, dict) and set(exp) == set(got), f"{where}: claves distintas"
        for k in exp:
            _compare(exp[k], got[k], f"{where}.{k}")
    elif isinstance(exp, list):
        assert isinstance(got, list) and len(exp) == len(got), f"{where}: longitud distinta"
        for i, (a, b) in enumerate(zip(exp, got)):
            _compare(a, b, f"{where}[{i}]")
    elif isinstance(exp, bool) or exp is None or isinstance(exp, str):
        assert type(exp) is type(got) and exp == got, f"{where}: {exp!r} != {got!r}"
    elif isinstance(exp, int):
        assert type(got) is int and exp == got, f"{where}: {exp!r} != {got!r}"
    else:
        assert isinstance(exp, float) and isinstance(got, float), f"{where}: tipos {exp!r} vs {got!r}"
        assert abs(exp - got) <= TOL, f"{where}: {exp!r} vs {got!r}"


@pytest.fixture(scope="module")
def recomputed() -> dict[str, dict[str, Any]]:
    t0 = time.perf_counter()
    out = compute_fixtures()
    assert time.perf_counter() - t0 < 20.0, "recalculo del golden > 20 s"
    return out


@pytest.mark.parametrize("name", FIXTURE_NAMES)
def test_baseline_golden(name: str, recomputed: dict[str, dict[str, Any]]) -> None:
    _compare(_load(name), recomputed[name], name)


def test_w0_digests_exact(recomputed: dict[str, dict[str, Any]]) -> None:
    for name in ("empty_phase", "complete_phase"):
        exp = [r["w0_digest"] for r in _load(name)["runs"]]
        got = [r["w0_digest"] for r in recomputed[name]["runs"]]
        assert exp == got
        assert all(re.fullmatch(r"[0-9a-f]{64}", d) for d in exp)


def _docstring_only(path: Path) -> bool:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return all(isinstance(s, ast.Expr) and isinstance(s.value, ast.Constant) for s in tree.body)


def _sha_modulo_version(rel: str) -> str:
    """sha256 del archivo con la version normalizada a 1.0.0 (solo metadato)."""
    path = ROOT / rel
    text = path.read_text(encoding="utf-8")
    if rel == "pyproject.toml":
        text = re.sub(r'(?m)^version = "[^"]*"$', 'version = "1.0.0"', text, count=1)
    else:
        text = re.sub(r'(?m)^__version__ = "[^"]*"$', '__version__ = "1.0.0"', text, count=1)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def test_omega10_manifest_unchanged() -> None:
    manifest: dict[str, str] = _load(MANIFEST_NAME)["sha256"]
    assert set(manifest) <= set(manifest_paths()), "archivos de Omega-1.0 eliminados"
    changed: list[str] = []
    for rel, digest in manifest.items():
        if sha256_file(ROOT / rel) == digest:
            continue
        if rel in ("omega/__init__.py", "pyproject.toml") and _sha_modulo_version(rel) == digest:
            continue
        if rel in DOCSTRING_ONLY and _docstring_only(ROOT / rel):
            continue
        changed.append(rel)
    assert not changed, f"archivos de Omega-1.0 modificados: {changed}"
