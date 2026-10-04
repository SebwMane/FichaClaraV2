"""Prueba estatica de arquitectura (ANALYSIS §6.2): RNG explicito y aislamiento de reference_graphs."""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "omega"
ALLOWED_NP_RANDOM = {"Generator", "PCG64", "SeedSequence", "BitGenerator"}
STDLIB_RNG = {"random", "secrets"}


def _py_files(base: Path) -> list[Path]:
    return sorted(p for p in base.rglob("*.py") if "__pycache__" not in p.parts)


def _parse(path: Path) -> ast.Module:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _is_numpy_random(node: ast.AST) -> bool:
    """True si node es la expresion `np.random` o `numpy.random`."""
    return (
        isinstance(node, ast.Attribute)
        and node.attr == "random"
        and isinstance(node.value, ast.Name)
        and node.value.id in {"np", "numpy"}
    )


def global_rng_violations(tree: ast.Module) -> list[str]:
    """Usos de np.random.<x> distintos de Generator/PCG64/SeedSequence y de modulos RNG de stdlib."""
    found: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute) and _is_numpy_random(node.value):
            if node.attr not in ALLOWED_NP_RANDOM:
                found.append(f"np.random.{node.attr} (linea {node.lineno})")
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            if node.module == "numpy.random":
                bad = [a.name for a in node.names if a.name not in ALLOWED_NP_RANDOM]
                if bad:
                    found.append(f"from numpy.random import {bad} (linea {node.lineno})")
            elif node.module.split(".")[0] in STDLIB_RNG:
                found.append(f"from {node.module} import ... (linea {node.lineno})")
        elif isinstance(node, ast.Import):
            for a in node.names:
                if a.name.split(".")[0] in STDLIB_RNG:
                    found.append(f"import {a.name} (linea {node.lineno})")
                elif a.name == "numpy.random" and a.asname is None:
                    found.append(f"import numpy.random (linea {node.lineno})")
    return found


def imported_modules(tree: ast.Module) -> list[str]:
    mods: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods += [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            base = ("." * node.level) + (node.module or "")
            mods.append(base)
            mods += [f"{base}.{a.name}" for a in node.names]
    return mods


def test_no_global_numpy_random_anywhere() -> None:
    offenders = {}
    for path in _py_files(PACKAGE) + _py_files(ROOT / "tests"):
        v = global_rng_violations(_parse(path))
        if v:
            offenders[str(path.relative_to(ROOT))] = v
    assert not offenders, f"RNG global prohibido: {offenders}"


def test_detector_catches_violations() -> None:
    bad = ast.parse("import numpy as np\nnp.random.seed(1)\nx = np.random.rand(3)\n")
    assert len(global_rng_violations(bad)) == 2
    assert global_rng_violations(ast.parse("import random\n"))
    assert global_rng_violations(ast.parse("from numpy.random import default_rng, PCG64\n"))
    ok = ast.parse("import numpy as np\ng = np.random.Generator(np.random.PCG64(1))\n")
    assert global_rng_violations(ok) == []


@pytest.mark.parametrize("subpackage", ["network", "dynamics"])
def test_core_does_not_import_reference_graphs(subpackage: str) -> None:
    base = PACKAGE / subpackage
    if not base.exists():
        pytest.skip(f"{subpackage}/ aun no existe")
    for path in _py_files(base):
        mods = imported_modules(_parse(path))
        bad = [m for m in mods if "reference_graphs" in m or "experiments" in m]
        assert not bad, f"{path.relative_to(ROOT)} importa {bad}"


def test_network_package_exists_and_checked() -> None:
    assert _py_files(PACKAGE / "network")


def test_io_uses_absolute_imports_only() -> None:
    for path in _py_files(PACKAGE / "io"):
        for node in ast.walk(_parse(path)):
            if isinstance(node, ast.ImportFrom):
                assert node.level == 0, f"{path.name}: import relativo prohibido en omega.io"


def test_no_mutable_module_state_in_wp1_modules() -> None:
    """Sin estado global mutable: ningun `global` ni asignacion de list/dict/set a nivel de modulo."""
    wp1 = [PACKAGE / "network", PACKAGE / "config", PACKAGE / "statistics", PACKAGE / "io"]
    for base in wp1:
        for path in _py_files(base):
            tree = _parse(path)
            assert not any(isinstance(n, ast.Global) for n in ast.walk(tree)), path
            for stmt in tree.body:
                value = getattr(stmt, "value", None)
                if isinstance(stmt, (ast.Assign, ast.AnnAssign)) and isinstance(
                    value, (ast.List, ast.Dict, ast.Set)
                ):
                    # constantes privadas inmutables por convencion (p. ej. _SECTIONS) se permiten
                    names = [getattr(t, "id", "") for t in getattr(stmt, "targets", [])]
                    names.append(getattr(getattr(stmt, "target", None), "id", ""))
                    assert all(n == "" or n.startswith("_") or n.isupper() for n in names), (path, names)
