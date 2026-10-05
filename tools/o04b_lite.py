"""O4B-1 enmendada (Rev. 2, R2.4/R2.5): clasificacion estructural y residuos KKT de los finales de O-04.

Entrada: `runs/omega11/o04_ablation/passports/OMEGA-EXP-*.json` + `.npz` (arrays `w0`, `w_final`: triangular superior).
Cada array se verifica con `omega.io.passport.array_digest` frente a `array_digests` del JSON (ademas, sha256 del
archivo .npz frente a `npz_sha256`). Si algo no coincide, el pasaporte se EXCLUYE y se informa.

Identificacion (inspeccion de un pasaporte real, schema 1.1):
  * Omega-B  : `config.fixed_density != null` (engine fixed_density). rho = config.fixed_density.rho.
  * S0       : `config.fixed_density == null` (engine gradient); la ablacion esta en `results.ablation`. Solo se informa.
  * N        : config.base.init.n;  parametros: config.base.functional (FunctionalParams crudos).
  * gamma_hat/alpha_hat se recuperan con `raw_to_reduced`; factor = alpha_hat / uniform_state_threshold(gamma_hat, rho, N).
  * Estado de la corrida: `termination_reason` ("CONVERGED" / ...).

Por final: `classify_structure` (R2.4), kappa_inj por componente (informe), P3 abiertos y, solo si
termination_reason == CONVERGED: residuos KKT (Omega-B: `kkt_residuals_fixed_density`; S0:
`kkt_residuals_s0`, SOLO para la ablacion "full": en las demas ablaciones el funcional dinamico no es el de
`config.base.functional`, asi que sus residuos no tendrian sentido; se informa "n/a").
Tolerancia de cero/uno de los residuos: 1e-6 (valor por defecto del modulo kkt).

Decision congelada (R2.5):
  EXITO   : 100% de los finales Omega-B verificados en {vacio, uniforme, clique_unica, multi_clique, cliques_solapadas}
            (con o sin halo) Y residuo KKT relativo <= 1e-7 en todos los convergidos.
  FRACASO : algun final Omega-B "otro" o algun convergido Omega-B con residuo relativo > 1e-7.
Se informa tambien por celda (rho, gamma_hat, factor) y el desglose de clases. Smoke: ~20 pasaportes (10 Omega-B y 10 S0
equiespaciados) y decision "N/A (smoke)". El modo full exige arbol git limpio (--allow-dirty solo con --smoke).
Salida: <out>/summary.json (fase preregistered -> final) y <out>/states.csv.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from omega.config.convert import raw_to_reduced  # noqa: E402
from omega.config.settings import FunctionalParams  # noqa: E402
from omega.dynamics.fixed_density import uniform_state_threshold  # noqa: E402
from omega.experiments.v11.gate import head_commit, tree_dirty  # noqa: E402
from omega.io.passport import array_digest  # noqa: E402
from omega.io.provenance import dependency_lock  # noqa: E402
from omega.landscape.kkt import kkt_residuals_fixed_density, kkt_residuals_s0  # noqa: E402
from omega.landscape.structure import classify_structure  # noqa: E402
from omega.network.weights import from_upper_triangle  # noqa: E402
from omega.phases.scan import to_jsonable  # noqa: E402
from omega.types import FloatArray  # noqa: E402

STEP = "o04b_lite"
KKT_REL_TOL = 1e-7
ALLOWED_CLASSES = ("vacío", "uniforme", "clique_única", "multi_clique", "cliques_solapadas")
DEFAULT_PASSPORTS = ROOT / "runs" / "omega11" / "o04_ablation" / "passports"
SMOKE_PER_BLOCK = 10

CSV_FIELDS = [
    "run_id", "block", "ablation", "n", "rho", "gamma_hat", "alpha_hat", "factor", "status", "steps",
    "state_class", "has_halo", "halo_is_clique_union", "n_components", "component_sizes", "kappa_inj", "open_p3",
    "kkt_rel", "kkt_max_violation", "kkt_multiplier", "kkt_ok",
]


def canonical_hash(cfg: dict[str, Any]) -> str:
    text = json.dumps(to_jsonable(cfg), sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_npz(npz_path: Path, expected: dict[str, str], npz_sha256: str | None = None) -> tuple[bool, str, dict[str, FloatArray]]:
    """Verifica cada array del .npz con `array_digest` (la misma funcion que usa el pasaporte) contra `expected`.

    Devuelve (ok, motivo, arrays). Si `npz_sha256` no es None tambien se comprueba el sha256 del archivo.
    """
    if not npz_path.is_file():
        return False, "npz ausente", {}
    if npz_sha256 is not None and file_sha256(npz_path) != npz_sha256:
        return False, "sha256 del archivo npz no coincide", {}
    with np.load(npz_path, allow_pickle=False) as data:
        arrays = {name: np.asarray(data[name]) for name in data.files}
    if set(arrays) != set(expected):
        return False, "los arreglos del npz no coinciden con array_digests", {}
    for name, arr in arrays.items():
        if array_digest(arr) != expected[name]:
            return False, f"digest alterado en {name!r}", {}
    return True, "", arrays


def describe_passport(p: dict[str, Any]) -> dict[str, Any]:
    """Metadatos de un pasaporte v1.1: bloque, ablacion, N, rho, parametros crudos y reducidos, factor, estado."""
    cfg = p["config"]
    base = cfg["base"]
    n = int(base["init"]["n"])
    fp = FunctionalParams(**{k: float(v) for k, v in base["functional"].items()})
    a_hat, g_hat = raw_to_reduced(fp, n)
    fd = cfg.get("fixed_density")
    out: dict[str, Any] = {
        "run_id": str(p["experiment_id"]), "n": n, "params": fp, "alpha_hat": round(float(a_hat), 9), "gamma_hat": round(float(g_hat), 9),
        "status": str(p["termination_reason"]), "steps": p.get("results", {}).get("steps"),
    }
    if fd is not None:
        rho = float(fd["rho"])
        out.update(block="omega_b", ablation="", rho=rho, factor=round(float(a_hat) / uniform_state_threshold(float(g_hat), rho, n), 6))
    else:
        out.update(block="s0", ablation=str(p.get("results", {}).get("ablation", "full")), rho=float("nan"), factor=float("nan"))
    return out


def analyse_state(w: FloatArray, meta: dict[str, Any]) -> dict[str, Any]:
    """Clase R2.4, diagnosticos y residuos KKT (solo convergidos; S0 solo ablacion full)."""
    cl = classify_structure(w)
    row: dict[str, Any] = {
        "state_class": cl.state_class, "has_halo": cl.has_halo, "halo_is_clique_union": cl.halo_is_clique_union,
        "n_components": len(cl.components), "component_sizes": list(cl.component_sizes),
        "kappa_inj": [None if not np.isfinite(k) else float(k) for k in cl.kappa_inj], "open_p3": cl.open_p3,
        "kkt_rel": None, "kkt_max_violation": None, "kkt_multiplier": None, "kkt_ok": None,
    }
    if meta["status"] != "CONVERGED":
        return row
    if meta["block"] == "omega_b":
        res = kkt_residuals_fixed_density(w, meta["params"])
    elif meta["ablation"] == "full":
        res = kkt_residuals_s0(w, meta["params"])
    else:
        return row
    row.update(kkt_rel=res.relative_violation, kkt_max_violation=res.max_violation, kkt_multiplier=res.multiplier,
               kkt_ok=bool(res.relative_violation <= KKT_REL_TOL))
    return row


def decide(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Decision congelada sobre las filas Omega-B (cada fila: state_class, status, kkt_ok)."""
    b = [r for r in rows if r["block"] == "omega_b"]
    bad_class = [r["run_id"] for r in b if r["state_class"] not in ALLOWED_CLASSES]
    bad_kkt = [r["run_id"] for r in b if r["status"] == "CONVERGED" and r["kkt_ok"] is not True]
    if not b:
        decision = "SIN DATOS Omega-B"
    elif bad_class or bad_kkt:
        decision = "FRACASO"
    else:
        decision = "EXITO"
    return {"decision": decision, "n_omega_b": len(b), "n_otro_or_disallowed": len(bad_class), "runs_disallowed_class": bad_class,
            "n_converged_not_kkt": len(bad_kkt), "runs_converged_not_kkt": bad_kkt}


def cell_report(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    cells: dict[tuple[float, float, float, int], list[dict[str, Any]]] = {}
    for r in rows:
        if r["block"] == "omega_b":
            cells.setdefault((r["rho"], r["gamma_hat"], r["factor"], r["n"]), []).append(r)
    out = []
    for (rho, g, f, n), rs in sorted(cells.items()):
        kk = [r["kkt_rel"] for r in rs if r["kkt_rel"] is not None]
        out.append({
            "rho": rho, "gamma_hat": g, "factor": f, "n": n, "n_runs": len(rs),
            "classes": dict(Counter(r["state_class"] for r in rs)), "n_halo": sum(1 for r in rs if r["has_halo"]),
            "n_converged": sum(1 for r in rs if r["status"] == "CONVERGED"),
            "n_converged_not_kkt": sum(1 for r in rs if r["status"] == "CONVERGED" and r["kkt_ok"] is not True),
            "max_kkt_rel": max(kk) if kk else None,
            "max_open_p3": max(r["open_p3"] for r in rs),
        })
    return out


def select(entries: list[tuple[Path, dict[str, Any]]], smoke: bool) -> list[tuple[Path, dict[str, Any]]]:
    if not smoke:
        return entries
    chosen: list[tuple[Path, dict[str, Any]]] = []
    for blk in ("omega_b", "s0"):
        sub = [e for e in entries if describe_passport(e[1])["block"] == blk]
        if sub:
            idx = np.unique(np.linspace(0, len(sub) - 1, min(SMOKE_PER_BLOCK, len(sub))).round().astype(int))
            chosen += [sub[int(i)] for i in idx]
    return chosen


def write_states(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        wr.writeheader()
        for r in rows:
            wr.writerow({k: (repr(v) if isinstance(v, float) else ("" if v is None else v)) for k, v in ((k, r.get(k)) for k in CSV_FIELDS)})
    os.replace(tmp, path)


def run_computation(passport_dir: Path, smoke: bool, log: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]], int]:
    jsons = sorted(passport_dir.glob("OMEGA-EXP-*.json"))
    entries = [(j, json.loads(j.read_text(encoding="utf-8"))) for j in jsons]
    entries = select(entries, smoke)
    rows: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    for k, (jpath, p) in enumerate(entries):
        ok, reason, arrays = verify_npz(jpath.parent / str(p["npz_file"]), dict(p["array_digests"]), p.get("npz_sha256"))
        if not ok:
            excluded.append({"passport": jpath.name, "reason": reason})
            continue
        meta = describe_passport(p)
        n = int(meta["n"])
        w = from_upper_triangle(arrays["w_final"], n)
        rows.append({k2: v for k2, v in meta.items() if k2 != "params"} | analyse_state(w, meta))
        if (k + 1) % 100 == 0:
            log(f"  {k + 1}/{len(entries)}")
    return rows, excluded, len(entries)


def summarize(rows: list[dict[str, Any]], excluded: list[dict[str, Any]], n_listed: int, smoke: bool) -> dict[str, Any]:
    res = decide(rows)
    if smoke:
        res["decision"] = f"N/A (smoke; valor provisional {res['decision']})"
    s0 = [r for r in rows if r["block"] == "s0"]
    return {
        "decision_result": res,
        "n_passports_considered": n_listed,
        "n_verified": len(rows),
        "n_excluded": len(excluded),
        "excluded": excluded,
        "class_breakdown_omega_b": dict(Counter(r["state_class"] for r in rows if r["block"] == "omega_b")),
        "halo_omega_b": dict(Counter(r["state_class"] + ("+halo" if r["has_halo"] else "") for r in rows if r["block"] == "omega_b")),
        "cells_omega_b": cell_report(rows),
        "s0_report_no_vote": {
            "n": len(s0),
            "class_breakdown": dict(Counter(r["state_class"] for r in s0)),
            "by_ablation": {a: dict(Counter(r["state_class"] for r in s0 if r["ablation"] == a)) for a in sorted({r["ablation"] for r in s0})},
            "n_converged_with_kkt": sum(1 for r in s0 if r["kkt_rel"] is not None),
            "n_converged_not_kkt": sum(1 for r in s0 if r["kkt_ok"] is False),
        },
        "kappa_inj_note": "kappa_inj por componente en states.csv (solo informe)",
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0] if __doc__ else "")
    ap.add_argument("--smoke", action="store_true", help="~20 pasaportes (no valida la decision)")
    ap.add_argument("--allow-dirty", action="store_true", help="solo smoke: permite arbol sucio (queda marcado)")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "o04b_lite")
    ap.add_argument("--passports", type=Path, default=DEFAULT_PASSPORTS)
    args = ap.parse_args(argv)

    dirty = tree_dirty()
    if dirty and not (args.smoke and args.allow_dirty):
        print("ERROR: arbol git sucio; el modo full exige commit limpio (--allow-dirty solo con --smoke).", file=sys.stderr)
        return 2
    cfg: dict[str, Any] = {
        "passport_dir": str(args.passports.relative_to(ROOT)) if args.passports.is_relative_to(ROOT) else str(args.passports),
        "smoke": bool(args.smoke), "kkt_rel_tol": KKT_REL_TOL, "allowed_classes": list(ALLOWED_CLASSES),
        "kkt_tol_zero_one": 1e-6, "smoke_per_block": SMOKE_PER_BLOCK,
    }
    chash = canonical_hash(cfg)
    commit = head_commit()
    meta: dict[str, Any] = {
        "step": STEP, "mode": "smoke" if args.smoke else "full", "code_commit": commit, "git_dirty": dirty,
        "allow_dirty_used": bool(args.allow_dirty and dirty), "config_hash": chash, "dependency_lock": dependency_lock(), "config": cfg,
    }
    summary_path = args.out / "summary.json"
    _atomic_write(summary_path, json.dumps(to_jsonable({**meta, "phase": "preregistered", "complete": False}), indent=2, allow_nan=False) + "\n")
    print(f"[O4B-1] preregistrado en {summary_path} (config_hash {chash[:12]}, commit {commit[:10]}, dirty={dirty})")

    t0 = time.perf_counter()
    rows, excluded, n_listed = run_computation(args.passports, args.smoke, lambda s: print(s, flush=True))
    elapsed = time.perf_counter() - t0
    write_states(args.out / "states.csv", rows)
    result = summarize(rows, excluded, n_listed, args.smoke)
    final = {**meta, "phase": "final", "complete": True, "elapsed_s": elapsed, "result": result}
    _atomic_write(summary_path, json.dumps(to_jsonable(final), indent=2, allow_nan=False) + "\n")
    d = result["decision_result"]
    print(f"[O4B-1] {result['n_verified']}/{n_listed} verificados, {result['n_excluded']} excluidos, {elapsed:.1f}s")
    print(f"[O4B-1] decision: {d['decision']}; Omega-B={d['n_omega_b']}; clase no permitida={d['n_otro_or_disallowed']}; convergidos no KKT={d['n_converged_not_kkt']}")
    print(f"[O4B-1] clases Omega-B: {result['class_breakdown_omega_b']}")
    print(f"[O4B-1] S0 (informe): {result['s0_report_no_vote']['class_breakdown']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
