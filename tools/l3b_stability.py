"""L-3b (Rev. 2, R2.5): estabilidad de arranques geometricos bajo S0 y Omega-B (N = 216).

Pregunta: los estados que S0 y Omega-B seleccionan cuando arrancan CERCA de una geometria, siguen siendo geometricos?

Todo lo de abajo es la lectura mas literal de R2.5 / R2.4; las decisiones de interpretacion estan marcadas con [I].

Inputs (N = 216), indice de input `ii`: 0 T3 = torus_lattice_3d(6); 1 RGG3 = rgg3_torus_binary(216, 12, rng del grafo s) con
  rng del grafo = SeedKey(20261005, (3, 0, s)); 2 T3_decorated = decorated_lattice_3d(3, 8); 3 ER = erdos_renyi_m con la misma m
  que el RGG3 de la semilla (rng [I]: SeedKey(20261005, (3, 0, 100 + s))); 4 K8_union = clique_union(216, 8).
  La semilla s usa el grafo s (T3, decorada y K8 son deterministas). La corrida eps = 0 (diagnostico, sin voto) usa el grafo 0.
Celdas: S0 = alpha_hat in {0.5, 1, 1.5, 1.9, 2.1, 2.5, 3} x gamma_hat in {0, 1, 10, 100} (28; cell_idx 0..27, alpha mayor-orden);
  Omega-B = f in {0.25, 0.5, 1, 1.5, 3} x gamma_hat in {0, 1, 10, 100} (20; cell_idx 28..47), rho = densidad del input
  = sum_{i<j} A / C(N,2), alpha_hat = f * uniform_state_threshold(gamma_hat, rho, N). [I] En RGG3 y ER, rho (y por tanto alpha_hat)
  se calcula con el grafo de CADA semilla (el input de la corrida es ese grafo); la celda agrupa las semillas.
Ruido: W0 = clip(A + eps * xi, 0, 1), xi = triangular superior de N(0,1) espejada con diagonal 0, rng PCG64 =
  make_rng(SeedKey(20261005, (3, 1 + input_idx, cell_idx, eps_idx, s))); eps_idx 0 -> 1e-2, 1 -> 1e-3, 2 -> 0 (s = 0 unica).
  En Omega-B ademas `project_fixed_density` a la masa del input (density_to_total(rho, N)).
Dinamica: `evolve` (S0) / `evolve_fixed_density` con FixedDensityConfig(rho = rho_input) (Omega-B); DynamicsConfig() por defecto
  salvo max_steps = 40000. [I] Para obtener el estado al 90% de los pasos se evoluciona en dos tramos: 36000 pasos y, SOLO si no
  convergio, 4000 mas desde ese estado (`evolve_two_stage`; con smoke max_steps = 2000 el tramo 1 es 1800). Estado y pasos
  finales son los del tramo 2; la historia de max_dw para el criterio de parada se reinicia (efecto despreciable).
  Deriva = R_orig(final) - R_orig(estado al 90%), solo si el tramo 1 no convergio.
Observables por final: clase R2.4; R_orig = sum_{(i,j) arista del input} W_ij / sum_{i<j} W_ij (0 si suma 0);
  rho_S = Spearman(W_ij, -d0_ij) sobre los pares i<j con d0 <= 4 (d0 = distancia de saltos del input; pares inalcanzables
  quedan fuera; scipy.stats.spearmanr; NaN si W o d0 son constantes -> cuenta como fallo); fraccion de componente gigante de
  {W > 0.1 max W} (0 si max = 0); P3 abiertos en {W > 1e-6}; residuos KKT relativos y lambda' (solo CONVERGED; tol cero/uno 1e-6);
  status y pasos.
Geometrico-persistente: R_orig >= 0.5 y rho_S >= 0.5 y gigante >= 0.5 (N/2 nodos) y clase fuera de {clique_unica, multi_clique,
  cliques_solapadas, vacio, uniforme} (es decir, clase "otro").
Estabilidad [I]: TODO final geometrico-persistente (cualquier eps) se re-perturba con eps = 1e-2 (xi como arriba con
  rng = make_rng(SeedKey(20261005, (3, 99, input_idx, cell_idx, eps_idx, s))), W0' = clip(W_final + 1e-2 xi, 0, 1), proyectado en
  Omega-B a la masa del input) y se re-evoluciona igual; "estable" si el nuevo final sigue siendo geometrico-persistente.
Validacion del clasificador (ANTES de las corridas, sobre los inputs sin evolucionar, grafo 0, R_orig = 1, rho_S sobre W = A):
  T3, RGG3 y decorada deben ser geometrico-persistentes; ER y K8_union no. Si falla: se escribe el summary con decision
  "CLASSIFIER_VALIDATION_FAILED" (con los valores medidos) y se sale con codigo 3, sin ajustar nada.
Voto por (input, celda): persistente si >= 2/3 de las semillas con eps = 1e-2 son geometrico-persistentes Y estables.
Decision global (en este orden): "NO" si alguna celda persistente con input T3 o RGG3; si no "NO-decorado" si solo hay celdas
  persistentes con T3_decorated; si no "METAESTABLE" si alguna celda de T3/RGG3/decorada [I: los controles no votan] tiene >= 2/3
  de semillas persistentes (mismo criterio geo+estable) con eps = 1e-3 pero no con 1e-2, o algun final MAX_STEPS con
  |deriva R_orig| > 0.01 en T3/RGG3/decorada; si no "SI (incompatible)". Informe sin voto: convergidos con residuo KKT relativo
  > 1e-7 ("no KKT"), controles ER/K8_union geometrico-persistentes (alarma del clasificador), eps = 0.
Certificado completo (collect_run_evidence + assess_run, `certificate_report`): solo sobre finales geometrico-persistentes
  (nivel 3, sin voto); desactivable con --no-certificate.
Paralelismo: Pool de 4 procesos, una corrida por tarea, cada tarea autocontenida y determinista; resultados anadidos a
  runs.jsonl (una fila por corrida, sin matrices) en cuanto terminan; --resume salta corridas ya hechas (config_hash igual).
  Matrices finales en <runs-dir>/<id>.npz con sha256 del archivo en la fila.
Modo full exige arbol git limpio (--allow-dirty solo con --smoke). Smoke: T3 y ER, 2 celdas S0 y 2 Omega-B, eps = 1e-2 una
  semilla, max_steps = 2000 (--skip-validation permitido solo en smoke, queda marcado).
"""

from __future__ import annotations

import os

# Un hilo BLAS por proceso: el paralelismo es por corrida (4 procesos).
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402
from collections import Counter  # noqa: E402
from dataclasses import dataclass  # noqa: E402
from functools import lru_cache  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402
from scipy.sparse import csr_matrix  # noqa: E402
from scipy.sparse.csgraph import connected_components, shortest_path  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from omega.certificate.evidence import collect_run_evidence, evidence_rng  # noqa: E402
from omega.config.convert import reduced_to_raw  # noqa: E402
from omega.config.seeds import SeedKey, make_rng  # noqa: E402
from omega.config.settings import DynamicsConfig, FunctionalParams  # noqa: E402
from omega.config.settings11 import FixedDensityConfig, Omega11Config  # noqa: E402
from omega.dynamics.evolution import evolve  # noqa: E402
from omega.dynamics.fixed_density import (  # noqa: E402
    density_to_total,
    evolve_fixed_density,
    project_fixed_density,
    uniform_state_threshold,
)
from omega.experiments.v11.gate import assessment_row, head_commit, tree_dirty  # noqa: E402
from omega.io.passport import array_digest  # noqa: E402
from omega.io.provenance import dependency_lock  # noqa: E402
from omega.landscape import references as R  # noqa: E402
from omega.landscape.kkt import kkt_residuals_fixed_density, kkt_residuals_s0  # noqa: E402
from omega.landscape.structure import classify_structure  # noqa: E402
from omega.network.weights import from_upper_triangle, upper_triangle  # noqa: E402
from omega.phases.finite_size import evidence_cfg  # noqa: E402
from omega.phases.scan import default_config, to_jsonable  # noqa: E402
from omega.types import FloatArray, RunStatus  # noqa: E402

STEP = "l3b_stability"
N = 216
MASTER_ENTROPY = 20261005
EXPERIMENT_ID = 3
INPUTS = ("T3", "RGG3", "T3_decorated", "ER", "K8_union")
GEOMETRIC_INPUTS = ("T3", "RGG3", "T3_decorated")
CONTROL_INPUTS = ("ER", "K8_union")
EPS_LIST = (1e-2, 1e-3, 0.0)  # eps_idx 0, 1, 2
STAB_STREAM = 99
ER_GRAPH_OFFSET = 100
MAX_STEPS_FULL = 40000
STAGE1_FRACTION = 0.9
HOP_MAX = 4
R_ORIG_MIN = 0.5
RHO_S_MIN = 0.5
GIANT_MIN = 0.5
GIANT_REL_THR = 0.1
NON_GEOMETRIC_CLASSES = ("clique_única", "multi_clique", "cliques_solapadas", "vacío", "uniforme")
KKT_REL_TOL = 1e-7
DRIFT_TOL = 0.01
VOTE_FRACTION = (2, 3)  # >= 2/3

FULL: dict[str, Any] = {
    "n": N,
    "inputs": list(INPUTS),
    "s0_alpha_hats": [0.5, 1.0, 1.5, 1.9, 2.1, 2.5, 3.0],
    "gamma_hats": [0.0, 1.0, 10.0, 100.0],
    "omega_b_factors": [0.25, 0.5, 1.0, 1.5, 3.0],
    "eps": [{"eps_idx": 0, "eps": 1e-2, "seeds": [0, 1, 2]}, {"eps_idx": 1, "eps": 1e-3, "seeds": [0, 1, 2]}, {"eps_idx": 2, "eps": 0.0, "seeds": [0]}],
    "max_steps": MAX_STEPS_FULL,
    "certificate": True,
}
SMOKE_OVERRIDE: dict[str, Any] = {
    "inputs": ["T3", "ER"],
    "smoke_s0_cells": [[1.9, 0.0], [3.0, 10.0]],
    "smoke_omega_b_cells": [[0.5, 10.0], [3.0, 0.0]],
    "eps": [{"eps_idx": 0, "eps": 1e-2, "seeds": [0]}],
    "max_steps": 2000,
    "certificate": False,
}


# ------------------------------------------------------------------ utilidades


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


def build_config(smoke: bool) -> dict[str, Any]:
    cfg = dict(FULL)
    if smoke:
        cfg.update(SMOKE_OVERRIDE)
    cfg.update(
        {
            "master_entropy": MASTER_ENTROPY,
            "experiment_id": EXPERIMENT_ID,
            "smoke": bool(smoke),
            "hop_max": HOP_MAX,
            "geometric_persistent": {"r_orig_min": R_ORIG_MIN, "rho_s_min": RHO_S_MIN, "giant_min_frac": GIANT_MIN,
                                      "giant_rel_threshold": GIANT_REL_THR, "excluded_classes": list(NON_GEOMETRIC_CLASSES)},
            "stage1_fraction": STAGE1_FRACTION,
            "kkt_rel_tol": KKT_REL_TOL,
            "drift_tol": DRIFT_TOL,
            "vote_fraction": list(VOTE_FRACTION),
            "stability_eps": EPS_LIST[0],
            "stability_stream": STAB_STREAM,
            "er_graph_offset": ER_GRAPH_OFFSET,
        }
    )
    return cfg


def cell_grid(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    """Celdas con cell_idx global: S0 0..27 (alpha mayor-orden), Omega-B 28..47 (f mayor-orden); el smoke filtra sin reindexar."""
    cells: list[dict[str, Any]] = []
    idx = 0
    for a in cfg["s0_alpha_hats"]:
        for g in cfg["gamma_hats"]:
            cells.append({"cell_idx": idx, "block": "s0", "alpha_hat": float(a), "gamma_hat": float(g), "factor": None})
            idx += 1
    for f in cfg["omega_b_factors"]:
        for g in cfg["gamma_hats"]:
            cells.append({"cell_idx": idx, "block": "omega_b", "alpha_hat": None, "gamma_hat": float(g), "factor": float(f)})
            idx += 1
    if cfg.get("smoke"):
        s0 = {tuple(c) for c in cfg["smoke_s0_cells"]}
        ob = {tuple(c) for c in cfg["smoke_omega_b_cells"]}
        cells = [c for c in cells if (c["block"] == "s0" and (c["alpha_hat"], c["gamma_hat"]) in s0)
                 or (c["block"] == "omega_b" and (c["factor"], c["gamma_hat"]) in ob)]
    return cells


def row_id(input_name: str, block: str, cell_idx: int, eps_idx: int, s: int) -> str:
    return f"{input_name}|{block}|c{cell_idx}|e{eps_idx}|s{s}"


def build_tasks(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    tasks: list[dict[str, Any]] = []
    cells = cell_grid(cfg)
    for name in cfg["inputs"]:
        ii = INPUTS.index(name)
        for c in cells:
            for e in cfg["eps"]:
                for s in e["seeds"]:
                    tasks.append({"id": row_id(name, c["block"], c["cell_idx"], e["eps_idx"], s), "input": name, "input_idx": ii,
                                  "cell": c, "eps_idx": int(e["eps_idx"]), "eps": float(e["eps"]), "s": int(s)})
    return tasks


# ------------------------------------------------------------------ inputs y geometria


def build_input(name: str, graph: int) -> FloatArray:
    """Matriz de adyacencia binaria del input `name` (grafo `graph` para RGG3/ER)."""
    if name == "T3":
        return R.torus_lattice_3d(6)
    if name == "RGG3":
        return R.rgg3_torus_binary(N, 12, make_rng(SeedKey(MASTER_ENTROPY, (EXPERIMENT_ID, 0, graph))))
    if name == "T3_decorated":
        return R.decorated_lattice_3d(3, 8)
    if name == "ER":
        m = int(round(float(build_input("RGG3", graph).sum()) / 2.0))
        return R.erdos_renyi_m(N, m, make_rng(SeedKey(MASTER_ENTROPY, (EXPERIMENT_ID, 0, ER_GRAPH_OFFSET + graph))))
    if name == "K8_union":
        return R.clique_union(N, 8)
    raise ValueError(f"input desconocido: {name}")


def graph_index(name: str, s: int, eps: float) -> int:
    """Indice de grafo: s para RGG3/ER con eps > 0; 0 si eps = 0 o el input es determinista."""
    return int(s) if (name in ("RGG3", "ER") and eps > 0.0) else 0


@dataclass(frozen=True, slots=True)
class InputGeometry:
    """Adyacencia del input y pares i<j con d0 <= 4 (distancia de saltos) para rho_S."""

    n: int
    adj: FloatArray  # binaria 0/1
    edge_upper: Any  # mascara booleana (M,) de aristas del input sobre triu_indices
    pair_mask: Any  # mascara booleana (M,) de pares con d0 <= HOP_MAX (finito)
    neg_d0: FloatArray  # -d0 sobre los pares de pair_mask


def make_geometry(adj: FloatArray, hop_max: int = HOP_MAX) -> InputGeometry:
    n = int(adj.shape[0])
    d = shortest_path(csr_matrix(adj), unweighted=True, directed=False)
    iu = np.triu_indices(n, k=1)
    d0 = np.asarray(d[iu], dtype=np.float64)
    pair_mask = np.isfinite(d0) & (d0 <= hop_max)
    return InputGeometry(n=n, adj=np.asarray(adj, dtype=np.float64), edge_upper=np.asarray(adj[iu] > 0.0),
                         pair_mask=pair_mask, neg_d0=-d0[pair_mask])


@lru_cache(maxsize=24)
def cached_input(name: str, graph: int) -> tuple[FloatArray, InputGeometry]:
    a = build_input(name, graph)
    return a, make_geometry(a)


# ------------------------------------------------------------------ observables puras


def r_orig(w: FloatArray, edge_upper: Any) -> float:
    """sum_{aristas del input} W_ij / sum_{i<j} W_ij (0 si la suma es 0)."""
    v = upper_triangle(w)
    tot = float(v.sum())
    return float(v[edge_upper].sum() / tot) if tot > 0.0 else 0.0


def rho_spearman(w: FloatArray, geom: InputGeometry) -> float:
    """Spearman(W_ij, -d0_ij) sobre los pares con d0 <= 4. NaN si algun vector es constante o no hay pares."""
    x = upper_triangle(w)[geom.pair_mask]
    y = geom.neg_d0
    if x.size < 3 or float(np.ptp(x)) == 0.0 or float(np.ptp(y)) == 0.0:
        return float("nan")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = float(spearmanr(x, y).statistic)
    return r


def giant_fraction(w: FloatArray, rel_thr: float = GIANT_REL_THR) -> float:
    """Tamano / N de la componente mas grande del grafo {W > rel_thr * max W} (0 si max W = 0)."""
    mx = float(np.max(w))
    if mx <= 0.0:
        return 0.0
    ncomp, labels = connected_components(csr_matrix(w > rel_thr * mx), directed=False)
    return float(np.max(np.bincount(labels))) / float(w.shape[0]) if ncomp > 0 else 0.0


def is_geometric_persistent(r: float, rho_s: float, giant: float, state_class: str) -> bool:
    """Regla congelada: R_orig >= 0.5, rho_S >= 0.5 (NaN falla), gigante >= 0.5 N y clase fuera de las 5 no geometricas."""
    if not math.isfinite(rho_s):
        return False
    return bool(r >= R_ORIG_MIN and rho_s >= RHO_S_MIN and giant >= GIANT_MIN and state_class not in NON_GEOMETRIC_CLASSES)


def observe(w: FloatArray, geom: InputGeometry, block: str, params: FunctionalParams | None, status: str) -> dict[str, Any]:
    """Observables de un final (todo salvo la dinamica). KKT solo si status == CONVERGED y hay `params`."""
    cl = classify_structure(w)
    ro = r_orig(w, geom.edge_upper)
    rs = rho_spearman(w, geom)
    gf = giant_fraction(w)
    out: dict[str, Any] = {
        "state_class": cl.state_class, "has_halo": cl.has_halo, "halo_is_clique_union": cl.halo_is_clique_union,
        "n_components": len(cl.components), "open_p3": cl.open_p3,
        "R_orig": ro, "rho_S": None if not math.isfinite(rs) else rs, "rho_S_nan": not math.isfinite(rs), "giant_frac": gf,
        "geo": is_geometric_persistent(ro, rs, gf, cl.state_class),
        "kkt_rel": None, "kkt_ok": None, "kkt_multiplier": None,
    }
    if status == RunStatus.CONVERGED.value and params is not None:
        res = kkt_residuals_fixed_density(w, params) if block == "omega_b" else kkt_residuals_s0(w, params)
        out.update(kkt_rel=res.relative_violation, kkt_ok=bool(res.relative_violation <= KKT_REL_TOL), kkt_multiplier=res.multiplier)
    return out


# ------------------------------------------------------------------ dinamica


def noisy_start(base: FloatArray, eps: float, rng: np.random.Generator) -> FloatArray:
    """clip(base + eps * xi, 0, 1), xi gaussiano simetrico (triangular superior espejada, diagonal 0)."""
    n = base.shape[0]
    xi = np.triu(rng.standard_normal((n, n)), k=1)
    xi = xi + xi.T
    w = np.clip(base + eps * xi, 0.0, 1.0)
    np.fill_diagonal(w, 0.0)
    return np.asarray(w, dtype=np.float64)


def evolve_two_stage(
    w0: FloatArray, block: str, params: FunctionalParams, rho: float, max_steps: int
) -> tuple[FloatArray, str, int, FloatArray | None]:
    """(w_final, status, pasos, w al 90% o None). Tramo 1 = int(0.9 max_steps); tramo 2 (max_steps - tramo 1) solo si el 1 no convergio."""
    n1 = int(STAGE1_FRACTION * max_steps)

    def run(w: FloatArray, steps: int) -> Any:
        cfg = DynamicsConfig(max_steps=steps)
        if block == "omega_b":
            return evolve_fixed_density(w, params, cfg, FixedDensityConfig(rho=rho))
        return evolve(w, params, cfg)

    t1 = run(w0, n1)
    if t1.status is not RunStatus.MAX_STEPS or max_steps - n1 < 1:
        return t1.w_final, t1.status.value, int(t1.steps), None
    t2 = run(t1.w_final, max_steps - n1)
    return t2.w_final, t2.status.value, n1 + int(t2.steps), t1.w_final


def cell_params(cell: dict[str, Any], rho: float) -> tuple[FunctionalParams, float]:
    """(FunctionalParams, alpha_hat) de la celda; en Omega-B alpha_hat = f * uniform_state_threshold(gamma_hat, rho, N)."""
    if cell["block"] == "omega_b":
        a_hat = float(cell["factor"]) * uniform_state_threshold(float(cell["gamma_hat"]), rho, N)
    else:
        a_hat = float(cell["alpha_hat"])
    return reduced_to_raw(a_hat, float(cell["gamma_hat"]), N), a_hat


def certificate_report(w: FloatArray, status: str, key: SeedKey) -> dict[str, Any]:
    """Certificado completo (collect_run_evidence + assess_run) sobre `w`; informe de nivel 3, sin voto."""
    try:
        base = Omega11Config(base=default_config(w.shape[0], master_entropy=MASTER_ENTROPY, experiment_id=EXPERIMENT_ID, replicates=1))
        c_ev = evidence_cfg(base, w)
        ev = collect_run_evidence(w, RunStatus(status), c_ev, evidence_rng(key))
        row = assessment_row(ev, c_ev)
        return {"codes": row["codes"], "primary": row["primary"], "passes": row["passes"], "flags_false": row["flags_false"],
                "evidence": row["evidence"]}
    except Exception as exc:  # informe: nunca tumba la corrida
        return {"error": f"{type(exc).__name__}: {exc}"}


# ------------------------------------------------------------------ una corrida


def _save_npz(path: Path, arrays: dict[str, FloatArray]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp.npz")
    np.savez_compressed(tmp, **arrays)  # type: ignore[arg-type]
    os.replace(tmp, path)
    return file_sha256(path)


def run_task(args: tuple[dict[str, Any], dict[str, Any]]) -> dict[str, Any]:
    """Corrida autocontenida y determinista: input -> ruido -> evolucion -> observables -> estabilidad -> (certificado)."""
    task, run_cfg = args
    t0 = time.perf_counter()
    name, ii, cell, eps_idx, eps, s = task["input"], int(task["input_idx"]), task["cell"], int(task["eps_idx"]), float(task["eps"]), int(task["s"])
    block, ci = str(cell["block"]), int(cell["cell_idx"])
    max_steps = int(run_cfg["max_steps"])
    a, geom = cached_input(name, graph_index(name, s, eps))
    total = float(a.sum()) / 2.0
    rho = total / (N * (N - 1) / 2.0)
    params, a_hat = cell_params(cell, rho)

    def start(base: FloatArray, key: SeedKey, e: float) -> FloatArray:
        w0 = noisy_start(base, e, make_rng(key))
        return project_fixed_density(w0, density_to_total(rho, N)) if block == "omega_b" else w0

    w0 = start(a, SeedKey(MASTER_ENTROPY, (EXPERIMENT_ID, 1 + ii, ci, eps_idx, s)), eps)
    wf, status, steps, w90 = evolve_two_stage(w0, block, params, rho, max_steps)
    obs = observe(wf, geom, block, params, status)
    drift = None if w90 is None else obs["R_orig"] - r_orig(w90, geom.edge_upper)
    row: dict[str, Any] = {
        "id": task["id"], "config_hash": run_cfg["config_hash"], "input": name, "input_idx": ii, "graph": graph_index(name, s, eps),
        "block": block, "cell_idx": ci, "alpha_hat": a_hat, "gamma_hat": float(cell["gamma_hat"]), "factor": cell["factor"],
        "rho_input": rho, "n_edges_input": total, "eps": eps, "eps_idx": eps_idx, "s": s, "status": status, "steps": steps,
        "max_steps_final": status == RunStatus.MAX_STEPS.value, "drift_R_orig": drift, **obs,
        "stable": None, "stab": None, "certificate": None,
    }
    arrays: dict[str, FloatArray] = {"w_final_upper": upper_triangle(wf)}
    if obs["geo"]:
        key = SeedKey(MASTER_ENTROPY, (EXPERIMENT_ID, STAB_STREAM, ii, ci, eps_idx, s))
        w0b = start(wf, key, EPS_LIST[0])
        wf2, st2, steps2, _ = evolve_two_stage(w0b, block, params, rho, max_steps)
        obs2 = observe(wf2, geom, block, params, st2)
        row["stable"] = bool(obs2["geo"])
        row["stab"] = {"status": st2, "steps": steps2, **obs2}
        arrays["w_stab_final_upper"] = upper_triangle(wf2)
        if run_cfg.get("certificate"):
            row["certificate"] = certificate_report(wf, status, SeedKey(MASTER_ENTROPY, (EXPERIMENT_ID, 98, ii, ci, eps_idx, s)))
    path = Path(run_cfg["runs_dir"]) / (str(task["id"]).replace("|", "_") + ".npz")
    row["npz"] = path.name
    row["npz_sha256"] = _save_npz(path, arrays)
    row["array_digests"] = {k: array_digest(v) for k, v in arrays.items()}
    row["elapsed_s"] = time.perf_counter() - t0
    return dict(to_jsonable(row))


# ------------------------------------------------------------------ validacion del clasificador


def classifier_validation() -> dict[str, Any]:
    """Inputs sin evolucionar (grafo 0): R_orig = 1, rho_S sobre W = A. Esperado: T3, RGG3, decorada geo; ER, K8 no."""
    out: dict[str, Any] = {}
    ok = True
    for name in INPUTS:
        a, geom = cached_input(name, 0)
        o = observe(a, geom, "s0", None, RunStatus.MAX_STEPS.value)
        want = name in GEOMETRIC_INPUTS
        out[name] = {"R_orig": o["R_orig"], "rho_S": o["rho_S"], "giant_frac": o["giant_frac"], "state_class": o["state_class"],
                     "geo": o["geo"], "expected_geo": want, "ok": o["geo"] == want}
        ok = ok and o["geo"] == want
    return {"passed": ok, "inputs": out}


# ------------------------------------------------------------------ decision


def _persistent_count(rows: list[dict[str, Any]]) -> tuple[int, int]:
    n = len(rows)
    return sum(1 for r in rows if r.get("geo") and r.get("stable")), n


def _vote(rows: list[dict[str, Any]]) -> bool:
    k, n = _persistent_count(rows)
    return n > 0 and k * VOTE_FRACTION[1] >= VOTE_FRACTION[0] * n


def group_cells(rows: list[dict[str, Any]]) -> dict[tuple[str, str, int], dict[int, list[dict[str, Any]]]]:
    out: dict[tuple[str, str, int], dict[int, list[dict[str, Any]]]] = {}
    for r in rows:
        out.setdefault((r["input"], r["block"], r["cell_idx"]), {}).setdefault(int(r["eps_idx"]), []).append(r)
    return out


def global_decision(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Decision global congelada sobre filas de corridas (ver docstring del modulo)."""
    cells = group_cells(rows)
    persistent_1e2: list[tuple[str, str, int]] = []
    persistent_1e3_only: list[tuple[str, str, int]] = []
    for key, by_eps in cells.items():
        p2 = _vote(by_eps.get(0, []))
        p3 = _vote(by_eps.get(1, []))
        if p2:
            persistent_1e2.append(key)
        if p3 and not p2 and key[0] in GEOMETRIC_INPUTS:
            persistent_1e3_only.append(key)
    drift_rows = [r["id"] for r in rows if r["input"] in GEOMETRIC_INPUTS and r.get("max_steps_final")
                  and r.get("drift_R_orig") is not None and abs(r["drift_R_orig"]) > DRIFT_TOL and int(r["eps_idx"]) in (0, 1)]
    inputs_persistent = {k[0] for k in persistent_1e2}
    if inputs_persistent & {"T3", "RGG3"}:
        decision = "NO"
    elif "T3_decorated" in inputs_persistent:
        decision = "NO-decorado"
    elif persistent_1e3_only or drift_rows:
        decision = "METAESTABLE"
    else:
        decision = "SÍ (incompatible)"
    conv = [r for r in rows if r["status"] == RunStatus.CONVERGED.value and r.get("kkt_ok") is not None]
    alarm = sorted({r["id"] for r in rows if r["input"] in CONTROL_INPUTS and r.get("geo")})
    return {
        "decision": decision,
        "persistent_cells_eps_1e-2": [{"input": k[0], "block": k[1], "cell_idx": k[2]} for k in sorted(persistent_1e2)],
        "metastable_cells_eps_1e-3_only": [{"input": k[0], "block": k[1], "cell_idx": k[2]} for k in sorted(persistent_1e3_only)],
        "max_steps_drift_rows": drift_rows,
        "n_converged_with_kkt": len(conv),
        "n_converged_not_kkt_no_vote": sum(1 for r in conv if r["kkt_ok"] is False),
        "control_geometric_persistent_alarm": {"alarm": bool(alarm), "ids": alarm},
    }


def per_input_report(rows: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for name in INPUTS:
        rs = [r for r in rows if r["input"] == name]
        if not rs:
            continue
        vote = [r for r in rs if int(r["eps_idx"]) in (0, 1)]
        out[name] = {
            "n_runs": len(rs),
            "classes": dict(Counter(r["state_class"] for r in rs)),
            "status": dict(Counter(r["status"] for r in rs)),
            "n_geo": sum(1 for r in rs if r["geo"]), "n_geo_stable": sum(1 for r in rs if r["geo"] and r["stable"]),
            "n_geo_votes_eps": sum(1 for r in vote if r["geo"]),
            "eps0_n_geo": sum(1 for r in rs if int(r["eps_idx"]) == 2 and r["geo"]),
            "n_max_steps": sum(1 for r in rs if r["max_steps_final"]),
            "n_rho_S_nan": sum(1 for r in rs if r["rho_S_nan"]),
        }
    return out


def cell_table(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for (name, block, ci), by_eps in sorted(group_cells(rows).items()):
        r0 = (by_eps.get(0) or by_eps.get(1) or by_eps.get(2))[0]  # type: ignore[index]
        out.append({
            "input": name, "block": block, "cell_idx": ci, "gamma_hat": r0["gamma_hat"], "alpha_hat": r0["alpha_hat"], "factor": r0["factor"],
            "classes_1e-2": dict(Counter(r["state_class"] for r in by_eps.get(0, []))),
            "geo_1e-2": sum(1 for r in by_eps.get(0, []) if r["geo"]), "persistent_1e-2": sum(1 for r in by_eps.get(0, []) if r["geo"] and r["stable"]),
            "persistent_1e-3": sum(1 for r in by_eps.get(1, []) if r["geo"] and r["stable"]),
            "n_1e-2": len(by_eps.get(0, [])), "n_1e-3": len(by_eps.get(1, [])),
            "cell_persistent_1e-2": _vote(by_eps.get(0, [])), "cell_persistent_1e-3": _vote(by_eps.get(1, [])),
        })
    return out


# ------------------------------------------------------------------ E/S


def read_rows(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                continue  # linea truncada por un corte: se recalcula
    return rows


def execute(tasks: list[dict[str, Any]], run_cfg: dict[str, Any], jsonl: Path, procs: int, log: Any) -> None:
    jsonl.parent.mkdir(parents=True, exist_ok=True)
    total = len(tasks)
    t0 = time.perf_counter()
    done = 0
    with jsonl.open("a", encoding="utf-8") as fh:
        payload = [(t, run_cfg) for t in tasks]
        if procs <= 1:
            it: Any = map(run_task, payload)
            pool = None
        else:
            pool = mp.get_context("fork").Pool(procs)
            it = pool.imap_unordered(run_task, payload, chunksize=1)
        try:
            for row in it:
                fh.write(json.dumps(row, allow_nan=False, sort_keys=True) + "\n")
                fh.flush()
                os.fsync(fh.fileno())
                done += 1
                if done % 10 == 0 or done == total:
                    log(f"  {done}/{total} corridas ({time.perf_counter() - t0:.0f}s)")
        finally:
            if pool is not None:
                pool.close()
                pool.join()


def time_probe(max_steps: int, log: Any) -> None:
    """Mide unas pocas corridas Omega-B reales con N = 216 (fuera del registro; no escribe nada)."""
    cfg = build_config(False)
    cells = {(c["factor"], c["gamma_hat"]): c for c in cell_grid(cfg) if c["block"] == "omega_b"}
    for name, key in (("T3", (0.5, 10.0)), ("RGG3", (1.5, 1.0)), ("T3_decorated", (3.0, 100.0)), ("RGG3", (0.25, 0.0))):
        c = cells[key]
        a, geom = cached_input(name, 0)
        rho = float(a.sum()) / 2.0 / (N * (N - 1) / 2.0)
        params, _ = cell_params(c, rho)
        w0 = project_fixed_density(noisy_start(a, 1e-2, make_rng(SeedKey(MASTER_ENTROPY, (EXPERIMENT_ID, 1000, 0, 0, 0)))), density_to_total(rho, N))
        t = time.perf_counter()
        wf, st, steps, _ = evolve_two_stage(w0, "omega_b", params, rho, max_steps)
        dt = time.perf_counter() - t
        o = observe(wf, geom, "omega_b", params, st)
        log(f"[sonda] {name} f={key[0]} g={key[1]}: {dt:.1f}s, {steps} pasos ({1e3 * dt / max(steps, 1):.2f} ms/paso), {st}, clase={o['state_class']}, geo={o['geo']}")


def summarize(rows: list[dict[str, Any]], tasks_expected: int) -> dict[str, Any]:
    return {
        "n_rows": len(rows), "n_expected": tasks_expected,
        "global": global_decision(rows) if rows else {"decision": "SIN DATOS"},
        "per_input": per_input_report(rows),
        "cells": cell_table(rows),
        "n_certificates": sum(1 for r in rows if r.get("certificate") is not None),
        "certificate_passes": sum(1 for r in rows if (r.get("certificate") or {}).get("passes") is True),
        "total_run_seconds": sum(float(r.get("elapsed_s") or 0.0) for r in rows),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0] if __doc__ else "")
    ap.add_argument("--smoke", action="store_true", help="malla minima (no valida la decision)")
    ap.add_argument("--allow-dirty", action="store_true", help="solo smoke: permite arbol sucio (queda marcado)")
    ap.add_argument("--skip-validation", action="store_true", help="solo smoke: omite la validacion del clasificador (queda marcado)")
    ap.add_argument("--resume", action="store_true", help="salta corridas ya presentes en runs.jsonl (mismo config_hash)")
    ap.add_argument("--summarize-only", action="store_true", help="solo re-agrega runs.jsonl en summary.json")
    ap.add_argument("--no-certificate", action="store_true", help="no calcula el certificado completo de los geometrico-persistentes")
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "l3b_stability")
    ap.add_argument("--runs-dir", type=Path, default=ROOT / "runs" / "l3b_stability")
    ap.add_argument("--time-probe", action="store_true", help="mide unas corridas Omega-B reales (N=216) y sale")
    args = ap.parse_args(argv)
    log = lambda s: print(s, flush=True)  # noqa: E731

    if args.time_probe:
        time_probe(MAX_STEPS_FULL, log)
        return 0
    if (args.skip_validation) and not args.smoke:
        print("ERROR: --skip-validation solo con --smoke.", file=sys.stderr)
        return 2
    dirty = tree_dirty()
    if dirty and not (args.smoke and args.allow_dirty):
        print("ERROR: arbol git sucio; el modo full exige commit limpio (--allow-dirty solo con --smoke).", file=sys.stderr)
        return 2
    cfg = build_config(args.smoke)
    if args.no_certificate:
        cfg["certificate"] = False
    chash = canonical_hash(cfg)
    commit = head_commit()
    meta: dict[str, Any] = {
        "step": STEP, "mode": "smoke" if args.smoke else "full", "code_commit": commit, "git_dirty": dirty,
        "allow_dirty_used": bool(args.allow_dirty and dirty), "validation_skipped": bool(args.skip_validation),
        "config_hash": chash, "dependency_lock": dependency_lock(), "config": cfg,
    }
    summary_path = args.out / "summary.json"
    jsonl = args.out / "runs.jsonl"
    tasks = build_tasks(cfg)

    if args.summarize_only:
        rows = [r for r in read_rows(jsonl) if r.get("config_hash") == chash]
        final = {**meta, "phase": "final", "complete": len(rows) == len(tasks), "result": summarize(rows, len(tasks))}
        _atomic_write(summary_path, json.dumps(to_jsonable(final), indent=2, allow_nan=False) + "\n")
        print(f"[L-3b] re-agregado: {len(rows)}/{len(tasks)} filas; decision {final['result']['global']['decision']}")
        return 0

    existing = read_rows(jsonl)
    if existing and not args.resume:
        print(f"ERROR: {jsonl} ya existe con {len(existing)} filas; usa --resume o borra el directorio.", file=sys.stderr)
        return 2
    if any(r.get("config_hash") != chash for r in existing):
        print("ERROR: runs.jsonl contiene filas con otro config_hash; no se reanuda.", file=sys.stderr)
        return 2
    _atomic_write(summary_path, json.dumps(to_jsonable({**meta, "phase": "preregistered", "complete": False}), indent=2, allow_nan=False) + "\n")
    print(f"[L-3b] preregistrado en {summary_path} (config_hash {chash[:12]}, commit {commit[:10]}, dirty={dirty}, {len(tasks)} corridas)")

    if not args.skip_validation:
        val = classifier_validation()
        if not val["passed"]:
            final = {**meta, "phase": "final", "complete": False, "aborted": True,
                     "result": {"decision": "CLASSIFIER_VALIDATION_FAILED", "classifier_validation": val}}
            _atomic_write(summary_path, json.dumps(to_jsonable(final), indent=2, allow_nan=False) + "\n")
            print("[L-3b] CLASSIFIER_VALIDATION_FAILED (sin ajustar):")
            for k, v in val["inputs"].items():
                print(f"    {k}: R_orig={v['R_orig']:.3f} rho_S={v['rho_S']} gigante={v['giant_frac']:.2f} clase={v['state_class']} geo={v['geo']} esperado={v['expected_geo']}")
            return 3
        meta["classifier_validation"] = val
        print("[L-3b] validacion del clasificador: OK")

    done_ids = {r["id"] for r in existing if (args.runs_dir / str(r.get("npz", "?"))).exists()}
    pending = [t for t in tasks if t["id"] not in done_ids]
    print(f"[L-3b] {len(done_ids)} ya hechas, {len(pending)} pendientes")
    run_cfg = {"max_steps": cfg["max_steps"], "config_hash": chash, "runs_dir": str(args.runs_dir), "certificate": bool(cfg["certificate"])}
    t0 = time.perf_counter()
    execute(pending, run_cfg, jsonl, args.procs, log)
    elapsed = time.perf_counter() - t0

    wanted = {t["id"] for t in tasks}
    rows = [r for r in read_rows(jsonl) if r.get("config_hash") == chash and r["id"] in wanted]
    uniq = {r["id"]: r for r in rows}
    rows = list(uniq.values())
    result = summarize(rows, len(tasks))
    final = {**meta, "phase": "final", "complete": len(rows) == len(tasks), "elapsed_s_this_invocation": elapsed, "result": result}
    if args.smoke:
        final["result"]["global"]["decision"] = f"N/A (smoke; valor provisional {result['global']['decision']})"
    _atomic_write(summary_path, json.dumps(to_jsonable(final), indent=2, allow_nan=False) + "\n")
    g = result["global"]
    print(f"[L-3b] {len(rows)}/{len(tasks)} corridas en {elapsed:.1f}s; decision: {g['decision']}")
    print(f"[L-3b] celdas persistentes (1e-2): {len(g['persistent_cells_eps_1e-2'])}; metaestables 1e-3: {len(g['metastable_cells_eps_1e-3_only'])}; "
          f"deriva MAX_STEPS: {len(g['max_steps_drift_rows'])}; convergidos no KKT: {g['n_converged_not_kkt_no_vote']}/{g['n_converged_with_kkt']}")
    print(f"[L-3b] alarma controles: {g['control_geometric_persistent_alarm']['alarm']}")
    for k, v in result["per_input"].items():
        print(f"[L-3b] {k}: clases={v['classes']} status={v['status']} geo={v['n_geo']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
