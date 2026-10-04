"""Pruebas de arquitectura de Omega-1.1 (WP-F, OMEGA_1_1_DESIGN §1.12 y §3.6). Solo lectura de fuentes (AST)."""

from __future__ import annotations

import ast
import re
from collections.abc import Iterator
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
OMEGA = ROOT / "omega"
NEW_TESTS = ("factories_v11.py", "test_taxonomy.py", "test_certificate.py", "test_architecture_v11.py")

D3_RE = re.compile(r"(?i)\b(d3|is_?3d|three_?d)\b")
D3_ALLOWED = {"RULE_D3_NEVER_SUFFICIENT"}
NP_RANDOM_OK = {"Generator", "default_rng", "SeedSequence", "PCG64", "BitGenerator"}


def _py_files(base: Path) -> list[Path]:
    return sorted(p for p in base.rglob("*.py") if "__pycache__" not in p.parts)


def _parse(path: Path) -> ast.Module:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _identifiers(tree: ast.AST) -> Iterator[tuple[str, int]]:
    for node in ast.walk(tree):
        line = getattr(node, "lineno", 0)
        if isinstance(node, ast.Name):
            yield node.id, line
        elif isinstance(node, ast.Attribute):
            yield node.attr, line
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            yield node.name, line
        elif isinstance(node, ast.arg):
            yield node.arg, line
        elif isinstance(node, ast.keyword) and node.arg is not None:
            yield node.arg, line
        elif isinstance(node, ast.alias):
            yield (node.asname or node.name).split(".")[-1], line
        elif isinstance(node, (ast.Global, ast.Nonlocal)):
            for n in node.names:
                yield n, line
        elif isinstance(node, ast.ExceptHandler) and node.name:
            yield node.name, line


def _d3_hits(tree: ast.AST) -> list[tuple[str, int]]:
    return [(i, ln) for i, ln in _identifiers(tree) if i not in D3_ALLOWED and D3_RE.search(i)]


# ------------------------------------------------------------------ (a) sin booleano "3D"


def test_d3_regex_self_check() -> None:
    bad = ast.parse("d3 = 1\nD3 = 2\nis_3d = 3\nis3d = 4\nthree_d = 5\nthreeD = 6\ndef f(is_3D): ...\nx.d3\n")
    assert len(_d3_hits(bad)) == 8
    ok = ast.parse("RULE_D3_NEVER_SUFFICIENT = 'x'\nd_volume = 1\ndimension_class = 3\nd30 = 2\n")
    assert _d3_hits(ok) == []


def test_no_d3_identifiers_in_omega() -> None:
    offenders: list[str] = []
    for path in _py_files(OMEGA):
        for ident, line in _d3_hits(_parse(path)):
            offenders.append(f"{path.relative_to(ROOT)}:{line}: {ident}")
    assert not offenders, "\n".join(offenders)


def test_only_allowed_d3_identifier_is_the_rule_constant() -> None:
    seen = {i for p in _py_files(OMEGA) for i, _ in _identifiers(_parse(p)) if "D3" in i or "d3" in i}
    assert seen <= D3_ALLOWED | {i for i in seen if not D3_RE.search(i)}
    assert "RULE_D3_NEVER_SUFFICIENT" in seen


# ------------------------------------------------------------------ (b) certificado solo en omega/certificate


def _called_names(tree: ast.AST) -> set[str]:
    out: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            f = node.func
            if isinstance(f, ast.Name):
                out.add(f.id)
            elif isinstance(f, ast.Attribute):
                out.add(f.attr)
    return out


def test_geometry_certificate_instantiated_only_in_certificate_package() -> None:
    offenders = [
        str(p.relative_to(ROOT))
        for p in _py_files(OMEGA)
        if "GeometryCertificate" in _called_names(_parse(p)) and p.parent != OMEGA / "certificate"
    ]
    assert not offenders, offenders
    assert "GeometryCertificate" in _called_names(_parse(OMEGA / "certificate" / "certificate.py"))


def test_candidate_verdict_only_produced_by_contracts_and_certificate() -> None:
    """GEOMETRIC_CANDIDATE solo se construye en contracts.py y omega/certificate/; en el resto, solo se compara."""
    offenders: list[str] = []
    for path in _py_files(OMEGA):
        if path == OMEGA / "contracts.py" or path.parent == OMEGA / "certificate":
            continue
        tree = _parse(path)
        parents = {child: parent for parent in ast.walk(tree) for child in ast.iter_child_nodes(parent)}
        for node in ast.walk(tree):
            if isinstance(node, ast.Attribute) and node.attr == "GEOMETRIC_CANDIDATE":
                if not isinstance(parents.get(node), ast.Compare):
                    offenders.append(f"{path.relative_to(ROOT)}:{node.lineno}")
    assert not offenders, offenders


def test_certificate_package_gates_candidate_on_all_fields() -> None:
    """El unico return/uso de GEOMETRIC_CANDIDATE en certificate.py esta dentro de finalize_verdict, tras satisfied_all."""
    src = (OMEGA / "certificate" / "certificate.py").read_text(encoding="utf-8")
    tree = ast.parse(src)
    uses = [n for n in ast.walk(tree) if isinstance(n, ast.Attribute) and n.attr == "GEOMETRIC_CANDIDATE"]
    assert len(uses) == 1
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "finalize_verdict")
    assert uses[0] in list(ast.walk(fn))
    fn_src = ast.get_source_segment(src, fn) or ""
    assert "satisfied_all" in fn_src and "dimension_class is not None" in fn_src and "len(codes) == 0" in fn_src
    assert not (OMEGA / "certificate" / "taxonomy.py").read_text(encoding="utf-8").count("GEOMETRIC_CANDIDATE")


# ------------------------------------------------------------------ (c) sin np.random global


def _global_random_hits(tree: ast.AST) -> list[int]:
    hits: list[int] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute) and node.attr not in NP_RANDOM_OK:
            v = node.value
            if (
                isinstance(v, ast.Attribute)
                and v.attr == "random"
                and isinstance(v.value, ast.Name)
                and v.value.id in ("np", "numpy")
            ):
                hits.append(node.lineno)
        elif isinstance(node, ast.Import):
            hits += [node.lineno for a in node.names if a.name == "random"]
        elif isinstance(node, ast.ImportFrom):
            if node.module == "random":
                hits.append(node.lineno)
            elif node.module == "numpy.random":
                hits += [node.lineno for a in node.names if a.name not in NP_RANDOM_OK]
            elif node.module == "numpy":
                hits += [node.lineno for a in node.names if a.name == "random"]
    return hits


def test_global_random_detector_self_check() -> None:
    assert _global_random_hits(ast.parse("import numpy as np\nnp.random.seed(1)\n")) == [2]
    assert _global_random_hits(ast.parse("import random\n")) == [1]
    assert _global_random_hits(ast.parse("np.random.default_rng(1)\nx: np.random.Generator\n")) == []


def test_no_global_np_random_in_omega_and_new_tests() -> None:
    offenders: list[str] = []
    files = _py_files(OMEGA) + [ROOT / "tests" / n for n in NEW_TESTS]
    for path in files:
        offenders += [f"{path.relative_to(ROOT)}:{ln}" for ln in _global_random_hits(_parse(path))]
    assert not offenders, offenders


# ------------------------------------------------------------------ (d) imports de omega.experiments


def _module_of(path: Path) -> str:
    rel = path.relative_to(ROOT).with_suffix("")
    parts = list(rel.parts)
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def _imported_modules(path: Path, tree: ast.AST) -> list[tuple[str, int]]:
    out: list[tuple[str, int]] = []
    pkg = _module_of(path).split(".")
    if path.name != "__init__.py":
        pkg = pkg[:-1]
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            out += [(a.name, node.lineno) for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = pkg[: len(pkg) - (node.level - 1)]
                mod = ".".join(base + ([node.module] if node.module else []))
            else:
                mod = node.module or ""
            out.append((mod, node.lineno))
            out += [(f"{mod}.{a.name}", node.lineno) for a in node.names]
    return out


def test_modules_do_not_import_experiments_except_v11() -> None:
    """Ningun modulo de omega fuera de omega/experiments (incluido experiments/v11) importa omega.experiments."""
    offenders: list[str] = []
    exp_dir = OMEGA / "experiments"
    for path in _py_files(OMEGA):
        inside = exp_dir in path.parents
        for mod, line in _imported_modules(path, _parse(path)):
            if not inside and (mod == "omega.experiments" or mod.startswith("omega.experiments.")):
                offenders.append(f"{path.relative_to(ROOT)}:{line} importa {mod}")
    assert not offenders, "\n".join(offenders)


def test_wp_f_modules_are_pure() -> None:
    """taxonomy/certificate solo dependen de contracts, settings11, types y entre si (sin E/S ni experimentos)."""
    allowed_prefix = ("omega.contracts", "omega.config.settings11", "omega.types", "omega.certificate")
    for name in ("taxonomy.py", "certificate.py"):
        path = OMEGA / "certificate" / name
        for mod, _ in _imported_modules(path, _parse(path)):
            if mod.startswith("omega."):
                assert mod.startswith(allowed_prefix), f"{name} importa {mod}"


@pytest.mark.parametrize("name", NEW_TESTS)
def test_new_test_files_exist(name: str) -> None:
    assert (ROOT / "tests" / name).is_file()
