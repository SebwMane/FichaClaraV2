"""Pasaporte de reproducibilidad (M§49): JSON + NPZ comprimido con digests sha256.

Unico modulo del nucleo con E/S. Imports absolutos siempre.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import re
from dataclasses import asdict
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Any

import numpy as np

import omega
from omega.config.convert import config_from_dict, config_to_dict, raw_to_reduced
from omega.config.settings import OmegaConfig
from omega.network.weights import upper_triangle
from omega.types import DimensionEstimate, FloatArray, RunResult

__all__ = [
    "PassportIntegrityError",
    "array_digest",
    "software_versions",
    "build_passport",
    "save_run",
    "load_run",
]

_RUN_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


class PassportIntegrityError(ValueError):
    """El digest de un archivo o arreglo no coincide con el registrado en el pasaporte."""


def array_digest(a: FloatArray) -> str:
    """sha256 hexadecimal de dtype, forma y bytes C-contiguos del arreglo."""
    arr = np.ascontiguousarray(a)
    h = hashlib.sha256()
    h.update(f"{arr.dtype.str}|{arr.shape}|".encode())
    h.update(arr.tobytes())
    return h.hexdigest()


def _file_digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def software_versions() -> dict[str, str]:
    """Versiones de Python y bibliotecas relevantes."""
    out = {"python": platform.python_version(), "omega": omega.__version__}
    for pkg in ("numpy", "scipy", "networkx"):
        try:
            out[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            out[pkg] = "unknown"
    return out


def _plain(obj: Any) -> Any:
    """Convierte a tipos JSON puros; no finitos -> None."""
    if isinstance(obj, dict):
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


def _dim_summary(e: DimensionEstimate) -> dict[str, Any]:
    return {
        "value": e.value,
        "stderr": e.stderr,
        "status": e.status,
        "window": None if e.window is None else list(e.window),
        "plateau": e.plateau,
        "method": e.method,
    }


def _as_upper(w: FloatArray) -> FloatArray:
    """Acepta matriz (N,N) o vector triangular superior; devuelve vector float64."""
    a = np.asarray(w, dtype=np.float64)
    return upper_triangle(a) if a.ndim == 2 else a.copy()


def _arrays(result: RunResult) -> dict[str, FloatArray]:
    traj, geo = result.trajectory, result.observables.geometry
    arrays: dict[str, FloatArray] = {
        "w0_upper": _as_upper(result.w0),
        "w_final_upper": _as_upper(traj.w_final),
        "snapshot_steps": np.asarray(traj.snapshot_steps),
        "snapshots_upper": np.asarray(traj.snapshots, dtype=np.float64),
        "deff_scales": np.asarray(geo.d_eff.scales, dtype=np.float64),
        "deff_profile": np.asarray(geo.d_eff.profile, dtype=np.float64),
        "ds_times": np.asarray(geo.d_s.scales, dtype=np.float64),
        "ds_profile": np.asarray(geo.d_s.profile, dtype=np.float64),
    }
    for name, series in traj.scalars.items():
        arrays[f"scalars_{name}"] = np.asarray(series, dtype=np.float64)
    return arrays


def build_passport(
    cfg: OmegaConfig, result: RunResult, versions: dict[str, str], git_commit: str
) -> dict[str, Any]:
    """Diccionario JSON-compatible del pasaporte (M§49) sin campos de archivo.

    `run_id`, `npz_file`, `npz_sha256` y `array_digests` los completa save_run.
    """
    topo, geo = result.observables.topology, result.observables.geometry
    traj = result.trajectory
    n = cfg.init.n
    alpha_hat, gamma_hat = raw_to_reduced(result.params, n)
    passport: dict[str, Any] = {
        "schema_version": cfg.schema_version,
        "run_id": None,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit,
        "versions": dict(versions),
        "config": config_to_dict(cfg),
        "seed": {
            "entropy": result.seed.entropy,
            "spawn_key": list(result.seed.spawn_key),
            "bit_generator": "PCG64",
        },
        "params_raw": asdict(result.params),
        "params_reduced": {"alpha_hat": alpha_hat, "gamma_hat": gamma_hat},
        "N": n,
        "dt": traj.dt,
        "tau_max": traj.tau,
        "steps": traj.steps,
        "status": traj.status.value,
        "omega0": {
            "distribution": cfg.init.distribution,
            "symmetrization": cfg.init.symmetrization,
            "sha256": array_digest(_as_upper(result.w0)),
        },
        "observables": {
            "mean_strength": topo.mean_strength,
            "std_strength": topo.std_strength,
            "clustering_weighted": topo.clustering_weighted,
            "clustering_binary": topo.clustering_binary,
            "path_length": geo.path_length,
            "path_length_hops": geo.path_length_hops,
            "diameter": geo.diameter,
            "giant_size": topo.giant_size,
            "giant_fraction": topo.giant_fraction,
            "topology": asdict(topo),
            "d_eff": _dim_summary(geo.d_eff),
            "d_eff_ball": _dim_summary(geo.d_eff_ball),
            "d_eff_hops": _dim_summary(geo.d_eff_hops),
            "d_s": _dim_summary(geo.d_s),
        },
        "phase": {
            "label": result.assessment.label.value,
            "flags": {k: bool(v) for k, v in result.assessment.flags.items()},
        },
        "curvature": None,
        "npz_file": None,
        "npz_sha256": None,
        "array_digests": {},
    }
    plain = _plain(passport)
    assert isinstance(plain, dict)
    return plain


def save_run(
    result: RunResult, cfg: OmegaConfig, directory: Path, run_id: str, git_commit: str
) -> tuple[Path, Path]:
    """Escribe `<run_id>.json` y `<run_id>.npz` en `directory`; devuelve (json, npz)."""
    if not _RUN_ID.match(run_id):
        raise ValueError(f"run_id invalido: {run_id!r}")
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    arrays = _arrays(result)
    npz_path = directory / f"{run_id}.npz"
    json_path = directory / f"{run_id}.json"
    np.savez_compressed(npz_path, **arrays)  # type: ignore[arg-type]
    passport = build_passport(cfg, result, software_versions(), git_commit)
    passport["run_id"] = run_id
    passport["npz_file"] = npz_path.name
    passport["npz_sha256"] = _file_digest(npz_path)
    passport["array_digests"] = {k: array_digest(v) for k, v in arrays.items()}
    json_path.write_text(
        json.dumps(passport, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )
    return json_path, npz_path


def load_run(json_path: Path) -> tuple[dict[str, Any], dict[str, FloatArray]]:
    """Lee pasaporte y arreglos verificando todos los digests.

    Lanza PassportIntegrityError (ValueError) si algun digest no coincide.
    """
    json_path = Path(json_path)
    passport = json.loads(json_path.read_text(encoding="utf-8"))
    for key in ("npz_file", "npz_sha256", "array_digests", "omega0", "config"):
        if passport.get(key) is None:
            raise ValueError(f"pasaporte incompleto: falta {key!r}")
    config_from_dict(passport["config"])  # valida la configuracion registrada
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
    if array_digest(arrays["w0_upper"]) != passport["omega0"]["sha256"]:
        raise PassportIntegrityError("digest de omega0 no coincide con w0_upper")
    return passport, arrays
