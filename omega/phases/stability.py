"""Estabilidad y robustez de fases (ANALYSIS §5.8, §4.1; M§26, M§29, M§31, M§32).

- seed_summary, permutation_invariance (exacta <= 1e-9, D-21), nulos (R6), histeresis por
  continuacion y confirmacion de F por punto.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Callable, Mapping, Sequence

import numpy as np

from omega.config.seeds import seed_key
from omega.config.settings import OmegaConfig, PhaseThresholds
from omega.network.weights import from_upper_triangle, permute, upper_triangle, validate_weight_matrix
from omega.phases.classification import classify_run
from omega.phases.scan import observe, point_params, simulate_from
from omega.types import FloatArray, IntArray, Observables, PhaseAssessment, PhaseLabel, RunResult, RunStatus

__all__ = [
    "seed_summary",
    "permutation_invariance",
    "shuffled_weight_null",
    "degree_preserving_null",
    "hysteresis_sweep",
    "hysteresis_gap",
    "wmin_candidate_flags",
    "confirm_geometric",
    "random_permutation",
]

PERMUTATION_TOL = 1e-9  # D-21


def seed_summary(
    results: Sequence[RunResult], key: Callable[[RunResult], float]
) -> tuple[float, float, int]:
    """(media, desviacion estandar muestral ddof=1 [0 si n=1], n) de key(r) sobre las replicas
    con valor finito; (NaN, NaN, 0) si no hay ninguno."""
    vals = np.asarray([key(r) for r in results], dtype=np.float64)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return (math.nan, math.nan, 0)
    std = float(vals.std(ddof=1)) if vals.size > 1 else 0.0
    return (float(vals.mean()), std, int(vals.size))


def _scalar_map(o: Observables) -> dict[str, float]:
    out: dict[str, float] = {}
    for f in dataclasses.fields(o.topology):
        out[f"topology.{f.name}"] = float(getattr(o.topology, f.name))
    geo = o.geometry
    for name in ("path_length", "path_length_hops", "diameter"):
        out[f"geometry.{name}"] = float(getattr(geo, name))
    for name in ("d_eff", "d_eff_ball", "d_eff_hops", "d_s"):
        est = getattr(geo, name)
        out[f"geometry.{name}"] = float(est.value)
        out[f"geometry.{name}.stderr"] = float(est.stderr)
    return out


def _diff(a: float, b: float) -> float:
    if math.isnan(a) and math.isnan(b):
        return 0.0
    if math.isnan(a) or math.isnan(b):
        return math.inf
    if math.isinf(a) and math.isinf(b) and a == b:
        return 0.0
    return abs(a - b)


def permutation_invariance(w: FloatArray, perm: IntArray, cfg: OmegaConfig) -> dict[str, float]:
    """Diferencia absoluta de cada observable escalar entre W y P W P^T (M§32, D-21).

    Claves: cada escalar, 'status_mismatches' (n.o de estimaciones con estado distinto) y
    'max_abs_diff' (maximo; inf si hay estados distintos). La invariancia exige <= 1e-9.
    """
    validate_weight_matrix(w)
    wp = permute(w, perm)
    o1, o2 = observe(w, cfg), observe(wp, cfg)
    s1, s2 = _scalar_map(o1), _scalar_map(o2)
    out = {k: _diff(s1[k], s2[k]) for k in s1}
    statuses = [(o1.geometry.d_eff.status, o2.geometry.d_eff.status), (o1.geometry.d_s.status, o2.geometry.d_s.status),
                (o1.geometry.d_eff_ball.status, o2.geometry.d_eff_ball.status),
                (o1.geometry.d_eff_hops.status, o2.geometry.d_eff_hops.status)]
    mism = float(sum(1 for a, b in statuses if a != b))
    out["status_mismatches"] = mism
    out["max_abs_diff"] = math.inf if mism > 0 else float(max(out.values()))
    return out


def shuffled_weight_null(w: FloatArray, w_min: float, rng: np.random.Generator) -> FloatArray:
    """Nulo de pesos barajados (R6): permuta al azar los M pesos entre TODAS las aristas
    (conserva exactamente el multiconjunto de pesos y la densidad). `w_min` solo se valida
    (contrato fijo de §5.8): el barajado es global, no depende del umbral."""
    validate_weight_matrix(w)
    if not 0.0 <= w_min <= 1.0:
        raise ValueError("w_min debe estar en [0,1]")
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng debe ser numpy.random.Generator")
    v = upper_triangle(w)
    return from_upper_triangle(rng.permutation(v), w.shape[0])


def degree_preserving_null(
    w: FloatArray, w_min: float, n_swaps: int, rng: np.random.Generator
) -> FloatArray:
    """Recableado con grados binarios preservados (R6): `n_swaps` intentos de intercambio doble
    de aristas A=(W>w_min): (a,b),(c,d) -> (a,d),(c,b) si ambas posiciones destino no son aristas
    ni autolazos. Se intercambian los VALORES de las entradas, de modo que el multiconjunto de
    pesos se conserva exactamente y A conserva su secuencia de grados."""
    validate_weight_matrix(w)
    if not 0.0 <= w_min <= 1.0:
        raise ValueError("w_min debe estar en [0,1]")
    if isinstance(n_swaps, bool) or not isinstance(n_swaps, int) or n_swaps < 0:
        raise ValueError("n_swaps debe ser entero >= 0")
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng debe ser numpy.random.Generator")
    m = np.array(w, dtype=np.float64, copy=True)
    iu, ju = np.triu_indices(m.shape[0], k=1)
    sel = m[iu, ju] > w_min
    ei, ej = iu[sel].copy(), ju[sel].copy()
    if ei.size < 2:
        return m
    for _ in range(n_swaps):
        e1, e2 = rng.integers(0, ei.size, size=2)
        if e1 == e2:
            continue
        a, b = int(ei[e1]), int(ej[e1])
        c, d = int(ei[e2]), int(ej[e2])
        if rng.random() < 0.5:
            c, d = d, c
        if len({a, b, c, d}) < 4:
            continue
        if m[a, d] > w_min or m[c, b] > w_min:
            continue
        wab, wcd, wad, wcb = m[a, b], m[c, d], m[a, d], m[c, b]
        m[a, b] = m[b, a] = wad
        m[c, d] = m[d, c] = wcb
        m[a, d] = m[d, a] = wab
        m[c, b] = m[b, c] = wcd
        ei[e1], ej[e1] = min(a, d), max(a, d)
        ei[e2], ej[e2] = min(c, b), max(c, b)
    return m


def hysteresis_sweep(
    w_start: FloatArray, alpha_hat_path: Sequence[float], gamma_hat: float, cfg: OmegaConfig
) -> tuple[RunResult, ...]:
    """Continuacion (M§31): en cada alpha_hat de la ruta se evoluciona desde el estado final del
    anterior (el primero, desde w_start). Claves: seed_key(cfg.seeds, i, 0); no usa RNG."""
    validate_weight_matrix(w_start)
    if len(alpha_hat_path) == 0:
        raise ValueError("alpha_hat_path no puede estar vacia")
    out: list[RunResult] = []
    w = w_start
    for i, a_hat in enumerate(alpha_hat_path):
        p = point_params(float(a_hat), gamma_hat, cfg.init.n, cfg.scan)
        r = simulate_from(cfg, p, seed_key(cfg.seeds, i, 0), w)
        out.append(r)
        w = r.trajectory.w_final
    return tuple(out)


def hysteresis_gap(
    fwd: Sequence[RunResult], bwd: Sequence[RunResult], key: Callable[[RunResult], float]
) -> FloatArray:
    """Brecha key(ida) - key(vuelta) por valor del parametro. `bwd` esta en orden de recorrido
    (inverso): bwd[len-1-i] corresponde a fwd[i] y deben coincidir en alpha_hat (ValueError si no)."""
    if len(fwd) != len(bwd) or len(fwd) == 0:
        raise ValueError("fwd y bwd deben tener la misma longitud (>0)")
    n = len(fwd)
    gap = np.empty(n, dtype=np.float64)
    for i in range(n):
        b = bwd[n - 1 - i]
        if not math.isclose(fwd[i].alpha_hat, b.alpha_hat, rel_tol=0.0, abs_tol=1e-9):
            raise ValueError("las rutas ida/vuelta no coinciden en alpha_hat")
        gap[i] = key(fwd[i]) - key(b)
    return gap


def wmin_candidate_flags(w: FloatArray, cfg: OmegaConfig) -> tuple[bool, ...]:
    """Bandera F-candidata del estado `w` para cada w_min de cfg.graph.w_min_sensitivity (R2)."""
    flags: list[bool] = []
    for wm in cfg.graph.w_min_sensitivity:
        c = dataclasses.replace(cfg, graph=dataclasses.replace(cfg.graph, w_min=wm))
        flags.append(bool(classify_run(observe(w, c), RunStatus.CONVERGED, cfg.phases).flags["F_candidate"]))
    return tuple(flags)


def _f_fraction(rs: Sequence[RunResult]) -> float:
    if len(rs) == 0:
        return 0.0
    return sum(1 for r in rs if r.assessment.flags.get("F_candidate", False)) / len(rs)


def _deff(r: RunResult) -> float:
    return float(r.observables.geometry.d_eff.value)


def confirm_geometric(
    point: Sequence[RunResult],
    by_size: Mapping[int, Sequence[RunResult]],
    neighbors: Sequence[Sequence[RunResult]],
    t: PhaseThresholds,
    *,
    permutation_max_diff: float | None = None,
    wmin_flags: Sequence[bool] | None = None,
    null_candidates: Sequence[bool] | None = None,
) -> PhaseAssessment:
    """F-confirmada por punto (§4.1). Criterios (todos requeridos; evidencia ausente = no cumple):
    seed_fraction (>= f_seed_fraction de replicas F-cand), seed_std (std D_eff <= f_seed_std),
    size_stability (|dD| <= f_size_tol entre tamanos, al menos 2, sin deriva monotona estricta con
    >=3), param_stability (cada vecino +-param_perturbation con fraccion F-cand >= f_seed_fraction),
    permutation (<= 1e-9), wmin_stability (>= 4/5 de w_min F-cand), null_distinct (fraccion de F-cand
    del nulo de pesos barajados <= 1 - f_seed_fraction). Etiqueta F si se confirma; si no, la moda
    de las replicas, con F-candidata no confirmada rebajada a U."""
    if len(point) == 0:
        raise ValueError("point no puede estar vacio")
    flags: dict[str, bool] = {}
    flags["seed_fraction"] = _f_fraction(point) >= t.f_seed_fraction
    _, std, n_ok = seed_summary(point, _deff)
    flags["seed_std"] = n_ok >= 2 and std <= t.f_seed_std

    means: list[float] = []
    for n_size in sorted(by_size):
        m, _, k = seed_summary(by_size[n_size], _deff)
        if k > 0:
            means.append(m)
    size_ok = len(means) >= 2 and len(means) == len(by_size)
    if size_ok:
        diffs = np.diff(np.asarray(means))
        size_ok = bool(np.max(np.abs(np.subtract.outer(means, means))) <= t.f_size_tol)
        if len(means) >= 3 and (np.all(diffs > 0) or np.all(diffs < 0)):
            size_ok = False
    flags["size_stability"] = size_ok
    flags["param_stability"] = len(neighbors) > 0 and all(_f_fraction(nb) >= t.f_seed_fraction for nb in neighbors)
    flags["permutation"] = permutation_max_diff is not None and permutation_max_diff <= PERMUTATION_TOL
    flags["wmin_stability"] = wmin_flags is not None and len(wmin_flags) > 0 and (
        sum(bool(f) for f in wmin_flags) / len(wmin_flags) >= 0.8
    )
    flags["null_distinct"] = null_candidates is not None and len(null_candidates) > 0 and (
        sum(bool(f) for f in null_candidates) / len(null_candidates) <= 1.0 - t.f_seed_fraction
    )
    confirmed = all(flags.values())
    flags["F_confirmed"] = confirmed
    if confirmed:
        label = PhaseLabel.F
    else:
        counts: dict[PhaseLabel, int] = {}
        for r in point:
            counts[r.assessment.label] = counts.get(r.assessment.label, 0) + 1
        label = max(sorted(counts, key=lambda lab: lab.value), key=lambda lab: counts[lab])
        if label is PhaseLabel.F:
            label = PhaseLabel.U
    return PhaseAssessment(label=label, flags=flags)


def random_permutation(n: int, rng: np.random.Generator) -> IntArray:
    """Permutacion aleatoria de 0..n-1 con rng explicito (auxiliar para M§32)."""
    return np.asarray(rng.permutation(n), dtype=np.int64)
