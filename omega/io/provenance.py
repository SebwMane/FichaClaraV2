"""Pasaporte ampliado Ω-1.1 (OMEGA_1_1_DESIGN §1.13, P§17). Único módulo de E/S de WP-E.

El pasaporte v1.0 queda intacto y se embebe en `omega10`. Imports absolutos siempre.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import re
import subprocess
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Any, Final

import numpy as np

from omega.baseline import BASELINE_COMMIT, BASELINE_NAME, MODEL_VERSION
from omega.config.convert11 import canonical_json, config11_from_dict, config11_to_dict, config_hash
from omega.config.seeds import SeedKey
from omega.config.settings11 import NullModel, Omega11Config
from omega.io.passport import PassportIntegrityError, array_digest
from omega.types import FloatArray

__all__ = [
    "TERMINATION_REASONS",
    "dependency_lock",
    "dependency_lock_hash",
    "git_info",
    "format_experiment_id",
    "experiment_file_stem",
    "parse_experiment_id",
    "allocate_experiment_id",
    "build_passport_v11",
    "save_passport_v11",
    "load_passport_v11",
    "reconstruct",
]

TERMINATION_REASONS: Final[frozenset[str]] = frozenset(
    {"CONVERGED", "MAX_STEPS", "NONFINITE", "STATIC", "CHAIN_COMPLETE", "CHAIN_NOT_EQUILIBRATED"}
)
_LABEL = re.compile(r"^Ω-EXP-(\d{6,})$")
_THRESHOLD_RULE: Final = (
    "distances/curvature/filtration: W>w_min (estricto); topology_observables Ω-1.0: W>=w_min"
)
_DIMENSION_ALGORITHM: Final = (
    "d_vol: shell log-grid Ω-1.0; d_s: lazy q=0.5; "
    "d_weyl: combinatorial binary staircase [10,0.2N] R2>=0.95"
)
_MAX_ID_RETRIES: Final = 10_000


def dependency_lock(packages: Sequence[str] = ("numpy", "scipy", "networkx")) -> dict[str, str]:
    """Versiones instaladas (incluye 'python'); 'unknown' si el paquete no está."""
    out = {"python": platform.python_version()}
    for pkg in packages:
        try:
            out[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            out[pkg] = "unknown"
    return out


def dependency_lock_hash(lock: Mapping[str, str]) -> str:
    """sha256 de '\\n'.join(sorted(f'{n}=={v}'))."""
    text = "\n".join(sorted(f"{n}=={v}" for n, v in lock.items()))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _git(repo: Path, *args: str) -> str:
    try:
        res = subprocess.run(
            ["git", "-C", str(repo), *args],
            capture_output=True, text=True, timeout=10, check=False,
        )
    except (OSError, subprocess.SubprocessError, ValueError):
        return "unknown"
    out = res.stdout.strip()
    return out if res.returncode == 0 and out else "unknown"


def git_info(repo: Path) -> tuple[str, str]:
    """(commit, rama) vía `git rev-parse HEAD` y `--abbrev-ref HEAD`; 'unknown' si falla."""
    return _git(Path(repo), "rev-parse", "HEAD"), _git(Path(repo), "rev-parse", "--abbrev-ref", "HEAD")


def format_experiment_id(number: int) -> str:
    """'Ω-EXP-%06d' (etiqueta JSON)."""
    if isinstance(number, bool) or not isinstance(number, int) or number < 0:
        raise ValueError(f"número de experimento inválido: {number!r}")
    return f"Ω-EXP-{number:06d}"


def experiment_file_stem(number: int) -> str:
    """'OMEGA-EXP-%06d' (nombre de archivo, solo ASCII)."""
    return format_experiment_id(number).replace("Ω", "OMEGA", 1)


def parse_experiment_id(label: str) -> int:
    """Inversa de format_experiment_id."""
    m = _LABEL.match(label) if isinstance(label, str) else None
    if m is None:
        raise ValueError(f"etiqueta de experimento inválida: {label!r}")
    return int(m.group(1))


def allocate_experiment_id(registry_dir: Path) -> int:
    """Reserva atómica: crea `registry/ids/NNNNNN.claim` con O_CREAT|O_EXCL (max+1, reintenta)."""
    ids = Path(registry_dir) / "ids"
    ids.mkdir(parents=True, exist_ok=True)
    for _ in range(_MAX_ID_RETRIES):
        existing = [int(p.stem) for p in ids.glob("*.claim") if p.stem.isdigit()]
        nxt = max(existing, default=0) + 1
        try:
            fd = os.open(ids / f"{nxt:06d}.claim", os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            continue
        os.close(fd)
        return nxt
    raise RuntimeError("no se pudo asignar un id de experimento")


def _distance_definition(cfg: Omega11Config) -> str:
    parts = {
        "hop": "hop: conteo de saltos BFS sobre A=(W>w_min)",
        "weighted_inverse": "weighted_inverse: d=1/W sobre aristas W>w_min",
        "weighted_log": "weighted_log: d=-log(max(W,log_floor)) sobre aristas W>w_min",
        "resistance": "resistance: resistencia efectiva sobre el grafo ponderado",
    }
    return "; ".join(parts[m.value] for m in cfg.distance.modes)


def _plain(obj: Any) -> Any:
    if isinstance(obj, Mapping):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _plain(obj.tolist())
    if isinstance(obj, np.generic):
        return _plain(obj.item())
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    return obj


def build_passport_v11(
    cfg: Omega11Config,
    *,
    experiment_number: int,
    seed: SeedKey,
    termination_reason: str,
    initialization_distribution: str,
    null_model: NullModel | None,
    coarse_graining: bool,
    results: Mapping[str, Any],
    code_commit: str,
    git_branch: str,
    entrypoint: str,
    w0_sha256: str | None,
    omega10_passport: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Diccionario JSON-compatible del pasaporte Ω-1.1 (todos los campos de §1.13)."""
    if not isinstance(cfg, Omega11Config):
        raise TypeError("cfg debe ser Omega11Config")
    if not isinstance(seed, SeedKey):
        raise TypeError("seed debe ser SeedKey")
    if termination_reason not in TERMINATION_REASONS:
        raise ValueError(f"termination_reason inválido: {termination_reason!r}")
    if null_model is not None and not isinstance(null_model, NullModel):
        raise TypeError("null_model debe ser NullModel o None")
    if ":" not in entrypoint:
        raise ValueError("entrypoint debe tener forma 'módulo:función'")
    lock = dependency_lock()
    cfg_dict = config11_to_dict(cfg)
    seed_d = {"entropy": seed.entropy, "spawn_key": list(seed.spawn_key), "bit_generator": "PCG64"}
    passport: dict[str, Any] = {
        "schema_version": cfg.schema_version,
        "model_version": MODEL_VERSION,
        "baseline": {"name": BASELINE_NAME, "commit": BASELINE_COMMIT},
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "code_commit": code_commit,
        "git_branch": git_branch,
        "python_version": lock["python"],
        "dependency_lock": lock,
        "dependency_lock_hash": dependency_lock_hash(lock),
        "config": cfg_dict,
        "config_hash": config_hash(cfg),
        "initialization_distribution": initialization_distribution,
        "distance_definition": _distance_definition(cfg),
        "threshold_rule": _THRESHOLD_RULE,
        "dimension_algorithm": _DIMENSION_ALGORITHM,
        "coarse_graining_algorithm": (
            f"{cfg.coarse.rule}/{cfg.coarse.aggregation}/v1" if coarse_graining else None
        ),
        "null_model": None if null_model is None else null_model.value,
        "random_seed": seed_d,
        "termination_reason": termination_reason,
        "experiment_id": format_experiment_id(experiment_number),
        "reconstruct": {
            "entrypoint": entrypoint,
            "config": cfg_dict,
            "seed": seed_d,
            "w0_sha256": w0_sha256,
        },
        "omega10": None if omega10_passport is None else dict(omega10_passport),
        "results": dict(results),
        "npz_file": None,
        "npz_sha256": None,
        "array_digests": {},
    }
    plain = _plain(passport)
    assert isinstance(plain, dict)
    canonical_json(plain)  # valida serializabilidad estricta
    return plain


def _file_digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def save_passport_v11(
    passport: Mapping[str, Any], arrays: Mapping[str, FloatArray], directory: Path
) -> tuple[Path, Path]:
    """Escribe `OMEGA-EXP-NNNNNN.json/.npz` en `directory`; devuelve (json, npz)."""
    number = parse_experiment_id(str(passport.get("experiment_id")))
    stem = experiment_file_stem(number)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    arrs = {k: np.asarray(v) for k, v in arrays.items()}
    npz_path = directory / f"{stem}.npz"
    json_path = directory / f"{stem}.json"
    np.savez_compressed(npz_path, **arrs)  # type: ignore[arg-type]
    out = dict(passport)
    out["npz_file"] = npz_path.name
    out["npz_sha256"] = _file_digest(npz_path)
    out["array_digests"] = {k: array_digest(v) for k, v in arrs.items()}
    json_path.write_text(
        json.dumps(_plain(out), indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return json_path, npz_path


def load_passport_v11(json_path: Path) -> tuple[dict[str, Any], dict[str, FloatArray]]:
    """Lee pasaporte y arreglos verificando sha256 del NPZ y de cada arreglo."""
    json_path = Path(json_path)
    passport = json.loads(json_path.read_text(encoding="utf-8"))
    for key in ("npz_file", "npz_sha256", "array_digests", "config", "reconstruct"):
        if passport.get(key) is None:
            raise ValueError(f"pasaporte incompleto: falta {key!r}")
    npz_name = str(passport["npz_file"])
    if Path(npz_name).name != npz_name:
        raise ValueError(f"npz_file no puede contener rutas: {npz_name!r}")
    npz_path = json_path.parent / npz_name
    if _file_digest(npz_path) != passport["npz_sha256"]:
        raise PassportIntegrityError(f"sha256 del NPZ no coincide: {npz_path}")
    with np.load(npz_path, allow_pickle=False) as data:
        arrays = {name: np.asarray(data[name]) for name in data.files}
    expected: dict[str, str] = passport["array_digests"]
    if set(expected) != set(arrays):
        raise PassportIntegrityError("los arreglos del NPZ no coinciden con los registrados")
    for name, arr in arrays.items():
        if array_digest(arr) != expected[name]:
            raise PassportIntegrityError(f"digest alterado en el arreglo {name!r}")
    return passport, arrays


def reconstruct(passport: Mapping[str, Any]) -> tuple[Omega11Config, SeedKey]:
    """(cfg, SeedKey) desde `passport['reconstruct']`; verifica config_hash si está presente."""
    rec = passport.get("reconstruct")
    if not isinstance(rec, Mapping):
        raise ValueError("pasaporte sin sección 'reconstruct'")
    cfg = config11_from_dict(rec["config"])
    seed_d = rec["seed"]
    seed = SeedKey(int(seed_d["entropy"]), tuple(int(x) for x in seed_d["spawn_key"]))
    expected = passport.get("config_hash")
    if expected is not None and config_hash(cfg) != expected:
        raise PassportIntegrityError("config_hash no coincide con la configuración reconstruida")
    return cfg, seed
