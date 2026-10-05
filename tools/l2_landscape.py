"""L-2 (Rev. 2, R2.5): paisaje extremal sobre referencias concretas, a masa fija (Omega-B), mu = 0.

Compara, por punto (rho, alpha_hat, gamma_hat, eta_hat) y semilla, la accion S0 minima de las referencias
GEOMETRICAS evaluadas con la de las NO geometricas evaluadas, ambas con sum_{i<j} W = W0 = rho C(N,2):

    Delta S = [min S(geometricas) - min S(no geometricas)] / (beta C(N,2)).

Decision congelada (R2.5):
  EXITO   : Delta S > 0 en el 100% de los puntos, en >= 3/5 semillas (semilla "buena" = todos sus puntos > 0).
  FRACASO : algun punto con Delta S < 0 en >= 3/5 semillas (esos puntos quedan marcados como prioritarios para L-3b).
  En otro caso: INDETERMINADO (se informa, no hay criterio congelado para ello).

Convenciones (documentadas, no ajustables a posteriori):
  * Parametros reducidos: alpha = 2 beta alpha_hat/(N-2), gamma = beta gamma_hat/(N-2) (omega.config.convert);
    eta = eta_hat * beta/(N-2)  (eta_hat = eta (N-2)/beta, ficha L-2). beta = 1, mu = 0.
  * Semillas: PCG64(SeedSequence(entropy=20261005, spawn_key=(EXPERIMENT_ID, stream, seed))) via omega.config.seeds.
    stream 0 = puntos/RGG3 de la semilla (compartidos entre rho); stream 1+i = G(n,m) del i-esimo rho.
  * Igualacion de masa: `match_mass` (arista parcial colex). Una referencia geometrica BINARIA cuya masa propia
    difiere de W0 en mas de MATCH_REL_TOL (relativo) se DESCARTA en ese rho (se informa): igualarla anadiria o
    quitaria una fraccion grande de aristas colex y dejaria de ser una referencia geometrica. Las no geometricas
    siempre se igualan. La RGG3 ponderada se ajusta por biseccion en `a` y se omite el par (perfil, r) si la masa
    no es alcanzable.
  * Union de cliques y caveman usan K_{k+1} para cada k de las referencias geometricas activas en ese rho
    (6 para T3 6^3, 12 para RGG3 binaria/ponderada, 13 para la decorada).
Salida: <out>/summary.json (fase "preregistered" antes de calcular, luego "final") y <out>/points.csv.
Modo full: exige arbol git limpio. --allow-dirty solo para smoke (queda marcado en el summary).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
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

from omega.config.convert import reduced_to_raw  # noqa: E402
from omega.config.seeds import SeedKey, make_rng  # noqa: E402
from omega.config.settings import FunctionalParams  # noqa: E402
from omega.dynamics.fixed_density import density_to_total, uniform_state_threshold  # noqa: E402
from omega.dynamics.functional import action, degree_irregularity, density, smoothness, triangles  # noqa: E402
from omega.experiments.v11.gate import head_commit, tree_dirty  # noqa: E402
from omega.io.provenance import dependency_lock  # noqa: E402
from omega.landscape import references as R  # noqa: E402
from omega.phases.scan import to_jsonable  # noqa: E402
from omega.types import FloatArray  # noqa: E402

N = 216
BETA = 1.0
MASTER_ENTROPY = 20261005
EXPERIMENT_ID = 2  # L-2
MATCH_REL_TOL = 0.10
# Diagnostico SIN voto (declarado antes del full, tras el smoke): con r grande la RGG ponderada degenera en casi
# uniforme; delta_S_local usa solo geometricas binarias y ponderadas con r <= LOCAL_R_MULT_MAX * r12.
LOCAL_R_MULT_MAX = 1.5
MIN_SEEDS = 3  # "3/5"
K_T3, K_RGG, K_DEC = 6, 12, 13
STEP = "l2_landscape"

FULL: dict[str, Any] = {
    "n": N,
    "rhos": [6 / 215, 12 / 215, 0.1],
    "alpha_hats": [0.5, 1.0, 1.5, 2.0, 3.0],
    "gamma_hats": [0.0, 1.0, 10.0, 100.0],
    "eta_hats": [0.0, 0.1, 1.0],
    "mu": 0.0,
    "beta": BETA,
    "n_seeds": 5,
    "r_multipliers": [float(x) for x in np.linspace(0.5, 3.0, 12)],
    "profiles": list(R.PROFILES),
}
# L-2b (Consejo Rev. 2, R2.8 L-A4): suplemento declarado tras L-2. Mismos criterios, referencias y semillas; la malla
# de alpha se expresa en la escala de Omega-B, alpha_hat = f * alpha_c(gamma_hat, rho, N) (la de L-3b), porque la malla
# absoluta de L-2 (alpha_hat <= 3) queda muy por debajo del umbral alpha_c ~ (1+gamma_hat)/rho de estas densidades.
L2B_OVERRIDE: dict[str, Any] = {
    "alpha_mode": "omega_b_factor",
    "factors": [0.25, 0.5, 1.0, 1.5, 3.0],
}
SMOKE_OVERRIDE: dict[str, Any] = {
    "rhos": [6 / 215, 0.1],
    "alpha_hats": [1.0, 3.0],
    "gamma_hats": [0.0, 10.0],
    "eta_hats": [0.0, 1.0],
    "n_seeds": 2,
}


def eta_from_reduced(eta_hat: float, n: int, beta: float = 1.0) -> float:
    """eta = eta_hat * beta/(n-2), es decir eta_hat = eta (n-2)/beta (ficha L-2)."""
    return eta_hat * beta / (n - 2)


def params_for(alpha_hat: float, gamma_hat: float, eta_hat: float, n: int) -> FunctionalParams:
    p = reduced_to_raw(alpha_hat, gamma_hat, n, BETA)
    return FunctionalParams(alpha=p.alpha, beta=p.beta, gamma=p.gamma, eta=eta_from_reduced(eta_hat, n, BETA), mu=0.0)


def canonical_hash(cfg: dict[str, Any]) -> str:
    text = json.dumps(to_jsonable(cfg), sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def build_config(smoke: bool, l2b: bool = False) -> dict[str, Any]:
    cfg = dict(FULL)
    cfg["alpha_mode"] = "absolute"
    if l2b:
        cfg.update(L2B_OVERRIDE)
        cfg.pop("alpha_hats")
    if smoke:
        cfg.update(SMOKE_OVERRIDE)
    cfg.update(
        {
            "master_entropy": MASTER_ENTROPY,
            "experiment_id": EXPERIMENT_ID,
            "match_rel_tol": MATCH_REL_TOL,
            "local_r_mult_max": LOCAL_R_MULT_MAX,
            "k_by_geometric": {"T3_6^3": K_T3, "RGG3_k12": K_RGG, "RGG3_weighted": K_RGG, "T3_decorated": K_DEC},
            "criteria": criteria(),
        }
    )
    return cfg


def criteria() -> dict[str, Any]:
    return {
        "observable": "DeltaS = [min S(geometricas) - min S(no geometricas)] / (beta*C(N,2))",
        "exito": "DeltaS>0 en 100% de los puntos, en >=3/5 semillas",
        "fracaso": "DeltaS<0 en algun punto, en >=3/5 semillas (puntos marcados prioritarios para L-3b)",
        "min_seeds": MIN_SEEDS,
        "eta_convention": "eta = eta_hat*beta/(N-2)",
        "geometricas": ["T3 6^3 binaria", "RGG3 k12 binaria", "RGG3 ponderada optimizada (perfil x r, a por biseccion)", "T3 decorada 3^3 x K8"],
        "no_geometricas": ["clique colex", "union K_{k+1}", "caveman conectado", "ER G(n,m)", "uniforme"],
    }


# ------------------------------------------------------------------ referencias


class Ref:
    """Referencia evaluada: nombre, detalle, matriz y componentes (T, S_dens, S_deg, S_smooth)."""

    __slots__ = ("name", "detail", "w", "comps", "local")

    def __init__(self, name: str, detail: str, w: FloatArray, local: bool = True) -> None:
        self.name = name
        self.local = local
        self.detail = detail
        self.w = w
        self.comps = np.array(
            [triangles(w), density(w), degree_irregularity(w), smoothness(w)], dtype=np.float64
        )


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def build_references(
    cfg: dict[str, Any], rho: float, rho_idx: int, seed: int
) -> tuple[list[Ref], list[Ref], list[dict[str, Any]]]:
    """(geometricas, no geometricas, informe de descartes/omisiones) para (rho, semilla)."""
    n = int(cfg["n"])
    total = density_to_total(rho, n)
    report: list[dict[str, Any]] = []
    geo: list[Ref] = []
    non: list[Ref] = []
    key0 = SeedKey(MASTER_ENTROPY, (EXPERIMENT_ID, 0, seed))

    def add_binary_geo(name: str, w: FloatArray) -> None:
        own = float(w.sum() / 2.0)
        rel = abs(own - total) / total
        if rel > MATCH_REL_TOL:
            report.append({"ref": name, "status": "descartada", "reason": f"masa propia {own:g} vs W0 {total:g} (rel {rel:.3f} > {MATCH_REL_TOL})"})
            return
        m = R.match_mass(w, total)
        if m is None:
            report.append({"ref": name, "status": "descartada", "reason": "no igualable"})
            return
        geo.append(Ref(name, f"own_mass={own:g}", m))

    add_binary_geo("T3_6^3", R.torus_lattice_3d(6))
    add_binary_geo("RGG3_k12", R.rgg3_torus_binary(n, K_RGG, make_rng(key0)))
    add_binary_geo("T3_decorated", R.decorated_lattice_3d(3, 8))

    pts = R.rgg3_torus_points(n, make_rng(key0))
    r12 = R.rgg3_radius(n, K_RGG)
    omitted = 0
    for profile in cfg["profiles"]:
        for mult in cfg["r_multipliers"]:
            r = float(mult) * r12
            try:
                a = R.fit_amplitude_for_mass(pts, r, profile, total)
            except ValueError:
                omitted += 1
                report.append({"ref": "RGG3_weighted", "status": "omitida", "reason": f"masa inalcanzable perfil={profile} r_mult={mult:.4f}"})
                continue
            w = R.weighted_rgg_from_points(pts, r, profile, a)
            geo.append(Ref("RGG3_weighted", f"{profile},r_mult={mult:.4f},a={a:.6g}", w, local=float(mult) <= LOCAL_R_MULT_MAX))

    active_k: set[int] = set()
    for g in geo:
        active_k.add({"T3_6^3": K_T3, "RGG3_k12": K_RGG, "RGG3_weighted": K_RGG, "T3_decorated": K_DEC}[g.name])

    def add_non(name: str, detail: str, w: FloatArray) -> None:
        m = R.match_mass(w, total)
        if m is None:
            report.append({"ref": name, "status": "descartada", "reason": "no igualable"})
        else:
            non.append(Ref(name, detail, m))

    add_non("clique_colex", "", R.colex_clique_with_mass(n, total))
    for k in sorted(active_k):
        add_non(f"clique_union_K{k + 1}", "", R.clique_union(n, k + 1))
        add_non(f"caveman_K{k + 1}", "", R.connected_caveman(n, k + 1))
    rng_er = make_rng(SeedKey(MASTER_ENTROPY, (EXPERIMENT_ID, 1 + rho_idx, seed)))
    add_non("ER_Gnm", f"m={int(math.floor(total + 1e-9))}", R.erdos_renyi_m(n, int(math.floor(total + 1e-9)), rng_er))
    add_non("uniform", "", R.uniform_with_mass(n, total))
    return geo, non, report


# ------------------------------------------------------------------ evaluacion


def evaluate(refs: list[Ref], coef: FloatArray) -> FloatArray:
    """S de cada referencia en cada punto de parametros: (n_ref, n_pts). coef columnas (-alpha, beta, gamma, eta)."""
    comps = np.stack([r.comps for r in refs])
    return np.asarray(comps @ coef.T, dtype=np.float64)


def self_check(refs: list[Ref], p: FunctionalParams) -> None:
    """La composicion por componentes debe coincidir con omega.dynamics.functional.action."""
    coef = np.array([[-p.alpha, p.beta, p.gamma, p.eta]])
    for r in refs:
        s = float(evaluate([r], coef)[0, 0])
        ref = action(r.w, p)
        if abs(s - ref) > 1e-8 * max(1.0, abs(ref)):
            raise AssertionError(f"S por componentes {s} != action {ref} en {r.name}")


def run_computation(cfg: dict[str, Any], log: Any) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    n = int(cfg["n"])
    by_factor = cfg.get("alpha_mode") == "omega_b_factor"
    first = cfg["factors"] if by_factor else cfg["alpha_hats"]
    base_grid = list(itertools.product(first, cfg["gamma_hats"], cfg["eta_hats"]))
    norm = BETA * n * (n - 1) / 2.0
    rows: list[dict[str, Any]] = []
    reports: dict[str, Any] = {}
    for rho_idx, rho in enumerate(cfg["rhos"]):
        if by_factor:
            grid = [(f * uniform_state_threshold(g, rho, n), g, e) for f, g, e in base_grid]
        else:
            grid = list(base_grid)
        params = [params_for(a, g, e, n) for a, g, e in grid]
        coef = np.array([[-p.alpha, p.beta, p.gamma, p.eta] for p in params], dtype=np.float64)
        for seed in range(int(cfg["n_seeds"])):
            t0 = time.perf_counter()
            geo, non, report = build_references(cfg, rho, rho_idx, seed)
            reports[f"rho={rho:.6f}|seed={seed}"] = {
                "geo_active": dict(Counter(g.name for g in geo)),
                "non_active": [r.name for r in non],
                "omissions": report,
            }
            if not geo:
                raise RuntimeError(f"sin referencias geometricas activas en rho={rho}, semilla={seed}")
            self_check(geo[:3] + non[:3], params[-1])
            s_geo = evaluate(geo, coef)
            s_non = evaluate(non, coef)
            ig = np.argmin(s_geo, axis=0)
            loc_idx = [i for i, g in enumerate(geo) if g.local]
            s_loc = s_geo[loc_idx] if loc_idx else None
            inn = np.argmin(s_non, axis=0)
            for p_idx, (a, g, e) in enumerate(grid):
                sg = float(s_geo[ig[p_idx], p_idx])
                sn = float(s_non[inn[p_idx], p_idx])
                gb, nb = geo[ig[p_idx]], non[inn[p_idx]]
                rows.append(
                    {
                        "rho": rho, "alpha_hat": a, "gamma_hat": g, "eta_hat": e, "seed": seed,
                        "factor": base_grid[p_idx][0] if by_factor else float("nan"),
                        "S_geo_min": sg, "geo_best": gb.name, "geo_best_detail": gb.detail,
                        "S_non_min": sn, "non_best": nb.name, "non_best_detail": nb.detail,
                        "delta_S": (sg - sn) / norm, "n_geo": len(geo), "n_non": len(non),
                        "delta_S_local": (float(s_loc[:, p_idx].min()) - sn) / norm if s_loc is not None else float("nan"),
                    }
                )
            log(f"  rho={rho:.5f} semilla={seed}: {len(geo)} geo, {len(non)} no-geo, {time.perf_counter() - t0:.1f}s")
    return rows, reports


def decide(rows: list[dict[str, Any]], n_seeds: int) -> dict[str, Any]:
    """Aplica la decision congelada sobre las filas de puntos."""
    by_point: dict[tuple[float, float, float, float], dict[int, float]] = {}
    by_seed: dict[int, list[float]] = {}
    for r in rows:
        key = (r["rho"], r["alpha_hat"], r["gamma_hat"], r["eta_hat"])
        by_point.setdefault(key, {})[r["seed"]] = r["delta_S"]
        by_seed.setdefault(r["seed"], []).append(r["delta_S"])
    good_seeds = [s for s, v in sorted(by_seed.items()) if all(x > 0.0 for x in v)]
    neg_points = [k for k, d in by_point.items() if sum(1 for x in d.values() if x < 0.0) >= MIN_SEEDS]
    pos_points_3 = [k for k, d in by_point.items() if sum(1 for x in d.values() if x > 0.0) >= MIN_SEEDS]
    full = n_seeds >= 5
    if not full:
        decision = "N/A (smoke: menos de 5 semillas)"
    elif len(good_seeds) >= MIN_SEEDS:
        decision = "EXITO"
    elif neg_points:
        decision = "FRACASO"
    else:
        decision = "INDETERMINADO"
    allv = [r["delta_S"] for r in rows]
    return {
        "decision": decision,
        "seeds_with_all_points_positive": good_seeds,
        "n_points": len(by_point),
        "n_points_negative_in_ge3_seeds": len(neg_points),
        "priority_points_for_L3b": [dict(zip(("rho", "alpha_hat", "gamma_hat", "eta_hat"), k)) for k in sorted(neg_points)],
        "n_points_positive_in_ge3_seeds": len(pos_points_3),
        "delta_S_min": min(allv),
        "delta_S_max": max(allv),
        "frac_rows_positive": sum(1 for x in allv if x > 0.0) / len(allv),
        "geo_best_counts": dict(Counter(r["geo_best"] for r in rows)),
        "diagnostic_local_no_vote": {
            "local_r_mult_max": LOCAL_R_MULT_MAX,
            "delta_S_local_min": min((r["delta_S_local"] for r in rows if math.isfinite(r["delta_S_local"])), default=float("nan")),
            "frac_rows_local_positive": sum(1 for r in rows if math.isfinite(r["delta_S_local"]) and r["delta_S_local"] > 0.0) / len(rows),
            "n_rows_local_negative": sum(1 for r in rows if math.isfinite(r["delta_S_local"]) and r["delta_S_local"] < 0.0),
        },
        "non_best_counts": dict(Counter(r["non_best"] for r in rows)),
    }


CSV_FIELDS = [
    "rho", "alpha_hat", "gamma_hat", "eta_hat", "seed", "S_geo_min", "geo_best", "geo_best_detail",
    "S_non_min", "non_best", "non_best_detail", "delta_S", "n_geo", "n_non", "delta_S_local", "factor",
]


def write_points(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        wr.writeheader()
        for r in rows:
            wr.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in r.items()})
    os.replace(tmp, path)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0] if __doc__ else "")
    ap.add_argument("--smoke", action="store_true", help="malla minima (no valida la decision)")
    ap.add_argument("--allow-dirty", action="store_true", help="solo smoke: permite arbol sucio (queda marcado)")
    ap.add_argument("--out", type=Path, default=None, help="directorio de salida")
    ap.add_argument("--l2b", action="store_true", help="suplemento L-2b: alpha_hat = f * alpha_c(gamma_hat, rho, N)")
    args = ap.parse_args(argv)

    dirty = tree_dirty()
    if dirty and not (args.smoke and args.allow_dirty):
        print("ERROR: arbol git sucio; el modo full exige commit limpio (--allow-dirty solo con --smoke).", file=sys.stderr)
        return 2
    cfg = build_config(args.smoke, args.l2b)
    if args.out is None:
        args.out = ROOT / "results" / ("landscape_l2b" if args.l2b else "landscape_l2")
    chash = canonical_hash(cfg)
    commit = head_commit()
    meta: dict[str, Any] = {
        "step": "l2b_landscape" if args.l2b else STEP,
        "mode": "smoke" if args.smoke else "full",
        "code_commit": commit,
        "git_dirty": dirty,
        "allow_dirty_used": bool(args.allow_dirty and dirty),
        "config_hash": chash,
        "dependency_lock": dependency_lock(),
        "config": cfg,
    }
    summary_path = args.out / "summary.json"
    _atomic_write(summary_path, json.dumps(to_jsonable({**meta, "phase": "preregistered", "complete": False}), indent=2, allow_nan=False) + "\n")
    print(f"[L-2] preregistrado en {summary_path} (config_hash {chash[:12]}, commit {commit[:10]}, dirty={dirty})")

    t0 = time.perf_counter()
    rows, reports = run_computation(cfg, lambda s: print(s, flush=True))
    elapsed = time.perf_counter() - t0
    write_points(args.out / "points.csv", rows)
    result = decide(rows, int(cfg["n_seeds"]))
    final = {
        **meta,
        "phase": "final",
        "complete": True,
        "elapsed_s": elapsed,
        "result": result,
        "reference_reports": reports,
        "scope_note": "Compara representantes evaluados; no minimiza globalmente. Redactar como 'ninguna de las referencias geometricas evaluadas...'.",
    }
    _atomic_write(summary_path, json.dumps(to_jsonable(final), indent=2, allow_nan=False) + "\n")
    print(f"[L-2] {result['n_points']} puntos x {cfg['n_seeds']} semillas en {elapsed:.1f}s")
    print(f"[L-2] decision: {result['decision']}; semillas con todos los puntos > 0: {result['seeds_with_all_points_positive']}")
    print(f"[L-2] puntos con DeltaS<0 en >=3 semillas: {result['n_points_negative_in_ge3_seeds']}; DeltaS in [{result['delta_S_min']:.4g}, {result['delta_S_max']:.4g}]")
    print(f"[L-2] mejor geometrica: {result['geo_best_counts']}")
    print(f"[L-2] mejor no geometrica: {result['non_best_counts']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
