"""Orden de P§22 codificado como compuerta, y utilidades comunes de los experimentos Omega-1.1 (DESIGN §3.7, §4).

`require_prerequisites(out_root, step)` exige que TODOS los pasos anteriores de `STEP_ORDER` tengan
`out_root/<paso>/summary.json` con `mode == "full"` y `complete == True`. `complete` significa "toda la malla
preregistrada se ejecuto", no "las expectativas se cumplieron" (eso va en `expectations_met`).

Enmienda A-13 (auditoria B7): `STEP_ORDER` empieza con `p01_analytical` y `p02_golden` (pasos 1-2 de P§22, ejecutados por
`test_p00_prereq`). Cada summary registra `code_commit`, `git_dirty` y `config_hash`; `require_prerequisites` exige que
cada summary previo sea del commit HEAD actual, con arbol limpio, y rechaza ejecutar si el arbol actual esta sucio.

E/S permitida: este modulo vive en `omega/experiments/`.
"""

from __future__ import annotations

import dataclasses
import json
import os
import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Final, Literal

from omega.certificate.evidence import evidence_summary
from omega.certificate.taxonomy import assess_run
from omega.config.convert11 import config_hash
from omega.config.seeds import SeedKey
from omega.config.settings import FunctionalParams
from omega.config.settings11 import Engine, FixedDensityConfig, NullModel, Omega11Config
from omega.contracts import RunEvidence
from omega.io.passport import array_digest
from omega.io.provenance import allocate_experiment_id, build_passport_v11, git_info, save_passport_v11
from omega.network.weights import upper_triangle
from omega.phases.scan import to_jsonable
from omega.types import FloatArray

__all__ = [
    "STEP_ORDER",
    "PrerequisiteError",
    "runs_root",
    "write_summary",
    "head_commit",
    "tree_dirty",
    "require_clean_tree",
    "require_prerequisites",
    "derive_cfg",
    "PassportWriter",
    "assessment_row",
    "check_expectation",
]

STEP_ORDER: Final = (
    "p01_analytical",
    "p02_golden",
    "o00_validation",
    "o01_baseline",
    "o02_distance",
    "o03_nulls",
    "o04_ablation",
    "s06_phase_diagram",
    "o06_finite_size",
    "o05_ensemble",
    "o08_topology",
    "o10_curvature",
    "o11_coarse",
)

Mode = Literal["smoke", "full"]


class PrerequisiteError(RuntimeError):
    """Falta un paso anterior en modo full completo."""


def runs_root() -> Path:
    """Raiz de salida: `OMEGA_RUNS_DIR` o `runs/omega11` (ignorada por git)."""
    return Path(os.environ.get("OMEGA_RUNS_DIR", "runs/omega11"))


_REPO_ROOT: Final = Path(__file__).resolve().parents[3]


def head_commit() -> str:
    """Commit HEAD del repositorio ('unknown' si no se puede leer)."""
    return git_info(_REPO_ROOT)[0]


def tree_dirty() -> bool:
    """True si `git status --porcelain` no es vacio (o si git falla: no se puede garantizar un arbol limpio)."""
    try:
        res = subprocess.run(
            ["git", "-C", str(_REPO_ROOT), "status", "--porcelain"],
            capture_output=True, text=True, timeout=30, check=False,
        )
    except (OSError, subprocess.SubprocessError, ValueError):
        return True
    return res.returncode != 0 or bool(res.stdout.strip())


def require_clean_tree(step: str) -> None:
    """Levanta PrerequisiteError si el arbol de trabajo actual esta sucio (o git no responde)."""
    if tree_dirty():
        raise PrerequisiteError(f"{step}: el arbol de git esta sucio; el modo full exige un commit limpio")


def write_summary(
    out_root: Path, name: str, data: Mapping[str, Any], *, mode: Mode, cfg: Omega11Config | None = None
) -> Path:
    """Escribe `out_root/<name>/summary.json` (JSON estricto, atomico) y devuelve su ruta.

    Anade `step`, `mode`, `code_commit`, `git_dirty`, `config_hash` (None si no se pasa `cfg`) y, si falta,
    `complete=False`. `complete` debe ser bool.
    """
    if mode not in ("smoke", "full"):
        raise ValueError(f"mode invalido: {mode!r}")
    if not name or Path(name).name != name:
        raise ValueError(f"nombre de paso invalido: {name!r}")
    out = dict(data)
    out["step"] = name
    out["mode"] = mode
    out["code_commit"] = head_commit()
    out["git_dirty"] = tree_dirty()
    out["config_hash"] = None if cfg is None else config_hash(cfg)
    out.setdefault("complete", False)
    if type(out["complete"]) is not bool:
        raise TypeError("complete debe ser bool")
    directory = Path(out_root) / name
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "summary.json"
    tmp = directory / "summary.json.tmp"
    tmp.write_text(json.dumps(to_jsonable(out), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(tmp, path)
    return path


def require_prerequisites(out_root: Path, step: str) -> None:
    """Levanta PrerequisiteError si el arbol esta sucio o algun paso anterior a `step` no tiene summary full completo,
    del commit HEAD actual y de un arbol limpio (Enmienda A-13)."""
    if step not in STEP_ORDER:
        raise ValueError(f"paso desconocido: {step!r}")
    require_clean_tree(step)
    head = head_commit()
    if head == "unknown":
        raise PrerequisiteError(f"{step}: no se pudo leer el commit HEAD")
    for prev in STEP_ORDER[: STEP_ORDER.index(step)]:
        path = Path(out_root) / prev / "summary.json"
        if not path.is_file():
            raise PrerequisiteError(f"{step}: falta {prev} (no existe {path})")
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise PrerequisiteError(f"{step}: {prev} tiene un summary ilegible: {exc}") from exc
        if not isinstance(data, dict) or data.get("mode") != "full" or data.get("complete") is not True:
            raise PrerequisiteError(f"{step}: {prev} no esta en modo full completo")
        if data.get("code_commit") != head:
            raise PrerequisiteError(f"{step}: {prev} es de otro commit ({data.get('code_commit')!r} != HEAD {head})")
        if data.get("git_dirty") is not False:
            raise PrerequisiteError(f"{step}: {prev} se ejecuto con el arbol sucio")


def derive_cfg(
    cfg: Omega11Config,
    *,
    n: int | None = None,
    functional: FunctionalParams | None = None,
    experiment_id: int | None = None,
    replicates: int | None = None,
    w_min: float | None = None,
    engine: Engine | None = None,
    fixed_density: FixedDensityConfig | None = None,
) -> Omega11Config:
    """Copia de `cfg` con los campos pedidos cambiados (los umbrales nunca se tocan)."""
    base = cfg.base
    if n is not None:
        base = dataclasses.replace(base, init=dataclasses.replace(base.init, n=n))
    if functional is not None:
        base = dataclasses.replace(base, functional=functional)
    if experiment_id is not None or replicates is not None:
        seeds = base.seeds
        if experiment_id is not None:
            seeds = dataclasses.replace(seeds, experiment_id=experiment_id)
        if replicates is not None:
            seeds = dataclasses.replace(seeds, replicates=replicates)
        base = dataclasses.replace(base, seeds=seeds)
    if w_min is not None:
        base = dataclasses.replace(base, graph=dataclasses.replace(base.graph, w_min=w_min))
    extra: dict[str, Any] = {}
    if engine is not None:
        extra["engine"] = engine
    if fixed_density is not None:
        extra["fixed_density"] = fixed_density
    return dataclasses.replace(cfg, base=base, **extra)


class PassportWriter:
    """Guarda pasaportes v1.1 (`OMEGA-EXP-NNNNNN.json/.npz`) con ids Omega-EXP asignados atomicamente."""

    def __init__(self, out_root: Path, step: str, entrypoint: str) -> None:
        self.directory = Path(out_root) / step / "passports"
        self.registry = Path(out_root) / "registry"
        self.entrypoint = entrypoint
        self.commit, self.branch = git_info(Path(__file__).resolve().parents[3])
        self.labels: list[str] = []

    def save(
        self,
        cfg: Omega11Config,
        *,
        seed: SeedKey,
        w: FloatArray,
        termination: str,
        init_distribution: str,
        results: Mapping[str, Any],
        null_model: NullModel | None = None,
        w0: FloatArray | None = None,
    ) -> str:
        """Pasaporte de un estado; devuelve su etiqueta `Ω-EXP-%06d`."""
        number = allocate_experiment_id(self.registry)
        arrays: dict[str, FloatArray] = {"w_final": upper_triangle(w)}
        if w0 is not None:
            arrays["w0"] = upper_triangle(w0)
        passport = build_passport_v11(
            cfg,
            experiment_number=number,
            seed=seed,
            termination_reason=termination,
            initialization_distribution=init_distribution,
            null_model=null_model,
            coarse_graining=False,
            results=to_jsonable(dict(results)),
            code_commit=self.commit,
            git_branch=self.branch,
            entrypoint=self.entrypoint,
            w0_sha256=None if w0 is None else array_digest(arrays["w0"]),
            omega10_passport=None,
        )
        save_passport_v11(passport, arrays, self.directory)
        label = str(passport["experiment_id"])
        self.labels.append(label)
        return label


def assessment_row(ev: RunEvidence, cfg: Omega11Config, *, full: bool = False) -> dict[str, Any]:
    """Fila JSON de una corrida: codigos, primario, banderas fallidas y resumen de evidencia."""
    a = assess_run(ev, cfg)
    row: dict[str, Any] = {
        "codes": [c.value for c in a.codes],
        "primary": None if a.primary is None else a.primary.value,
        "passes": a.passes,
        "flags_false": sorted(k for k, v in a.flags.items() if not v and k not in ("empty", "dense_or_uniform", "small_world")),
        "flags_true_trivial": sorted(k for k in ("empty", "dense_or_uniform", "small_world") if a.flags[k]),
    }
    summ = evidence_summary(ev)
    if full:
        row["evidence"] = summ
    else:
        row["evidence"] = {
            k: summ[k]
            for k in (
                "status", "giant_fraction", "mean_weight", "binary_density", "d_vol", "d_s", "d_weyl",
                "consensus_dimension", "dimension_class",
            )
        }
    return row


def check_expectation(op: str, measured: Any, target: Any, tol: float | None = None) -> bool:
    """Evalua una expectativa preregistrada. `None` (medida no finita/ausente) nunca la cumple.

    ops: "abs_within" (|m-target|<=tol), "lt", "gt", "eq", "contains" (target in m), "not_contains".
    """
    if measured is None:
        return False
    if op == "abs_within":
        assert tol is not None
        return bool(abs(float(measured) - float(target)) <= tol)
    if op == "lt":
        return bool(float(measured) < float(target))
    if op == "gt":
        return bool(float(measured) > float(target))
    if op == "eq":
        return bool(measured == target)
    if op == "contains":
        return bool(target in measured)
    if op == "not_contains":
        return bool(target not in measured)
    raise ValueError(f"op desconocido: {op!r}")
