"""GeometryCertificate y veredicto por punto (OMEGA_1_1_DESIGN §1.11, §1.12, §3.6). Logica pura.

Regla de oro: GEOMETRIC_CANDIDATE exige los 14 campos True, `dimension_class` no None y conjunto de codigos vacio.
D~3 jamas basta por si sola (`RULE_D3_NEVER_SUFFICIENT`).
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Mapping, Sequence
from typing import Final

import numpy as np

from omega.certificate.taxonomy import assess_run, primary_code, run_flags
from omega.config.settings11 import NullModel, Omega11Config
from omega.contracts import (
    CERTIFICATE_FIELDS,
    FAILURE_PRECEDENCE,
    FIELD_CODE,
    RULE_D3_NEVER_SUFFICIENT as _RULE,
    EnsembleSummary,
    FailureCode,
    GeometryCertificate,
    PointVerdict,
    RunEvidence,
    Verdict,
)

__all__ = [
    "RULE_D3_NEVER_SUFFICIENT",
    "RUN_LEVEL_FIELDS",
    "finalize_verdict",
    "seed_robust",
    "size_robust",
    "null_separated",
    "build_certificate",
    "point_verdict",
]

RULE_D3_NEVER_SUFFICIENT: Final = _RULE

RUN_LEVEL_FIELDS: Final[tuple[str, ...]] = (
    "connected",
    "nontrivial",
    "locality",
    "d_volume",
    "d_spectral",
    "d_weyl",
    "metric_robust",
    "isotropy_ok",
    "homogeneity_ok",
    "topology_stable",
    "manifold_proxy_ok",
)

_Evidence = dict[str, float | int | str | bool | None]


def _finite_dstar(runs: Sequence[RunEvidence]) -> list[float]:
    return [r.consensus_dimension for r in runs if math.isfinite(r.consensus_dimension)]


def _clean(runs: Sequence[RunEvidence], cfg: Omega11Config) -> list[bool]:
    return [assess_run(r, cfg).passes for r in runs]


def finalize_verdict(cert: GeometryCertificate, codes: Sequence[FailureCode]) -> Verdict:
    """Candidato solo con certificado completo (14 campos True y clase no None) y SIN ningun codigo."""
    if cert.satisfied_all and cert.dimension_class is not None and len(codes) == 0:
        return Verdict.GEOMETRIC_CANDIDATE
    return Verdict.NOT_CANDIDATE


def seed_robust(runs: Sequence[RunEvidence], cfg: Omega11Config) -> tuple[bool, dict[str, float]]:
    """>= 80% sin codigo, std(D*) <= 0.2 y >= 80% con la misma clase. Exige >= 2 corridas (conservador)."""
    c = cfg.certificate
    n = len(runs)
    ev: dict[str, float] = {"seed_n": float(n)}
    if n < 2:
        ev["seed_missing"] = 1.0
        return False, ev
    frac_clean = sum(_clean(runs, cfg)) / n
    ds = _finite_dstar(runs)
    std = float(np.std(ds)) if len(ds) >= 1 else math.nan
    classes = Counter(r.dimension_class for r in runs if r.dimension_class is not None)
    agree = (max(classes.values()) / n) if classes else 0.0
    ev.update(seed_frac_clean=frac_clean, seed_std_dstar=std, seed_class_agreement=agree)
    ok = (
        frac_clean >= c.seed_fraction
        and len(ds) / n >= c.seed_fraction
        and math.isfinite(std)
        and std <= c.seed_std
        and agree >= c.seed_fraction
    )
    return bool(ok), ev


def size_robust(by_size: Mapping[int, Sequence[RunEvidence]], cfg: Omega11Config) -> tuple[bool, dict[str, float]]:
    """Robustez en tamano (§1.12, B9). Se evalua el SUFIJO CONTIGUO de tamanos limpios (fraccion de corridas sin
    codigo >= seed_fraction) que incluye al mayor; debe tener >= `min_sizes` tamanos. Un tamano pequeno malo
    fuera del sufijo no invalida; un hueco malo dentro lo corta. Se registra `size_suffix_min` (menor N)."""
    c = cfg.certificate
    sizes_all = sorted(n for n, rs in by_size.items() if len(rs) > 0)
    ev: dict[str, float] = {"size_n_sizes": float(len(sizes_all))}
    if len(sizes_all) < c.min_sizes:
        ev["size_missing"] = 1.0
        return False, ev
    fracs = {n: sum(_clean(by_size[n], cfg)) / len(by_size[n]) for n in sizes_all}
    for n in sizes_all:
        ev[f"size_{n}_frac_clean"] = fracs[n]
    sizes: list[int] = []
    for n in reversed(sizes_all):
        if fracs[n] < c.seed_fraction:
            break
        sizes.append(n)
    sizes.reverse()
    ev["size_suffix_n"] = float(len(sizes))
    ev["size_suffix_min"] = float(sizes[0]) if sizes else 0.0
    if len(sizes) < c.min_sizes:
        return False, ev
    clean_ok = True
    means: list[float] = []
    for n in sizes:
        ds = _finite_dstar(by_size[n])
        means.append(float(np.mean(ds)) if ds else math.nan)
    if not all(math.isfinite(m) for m in means):
        ev["size_missing"] = 1.0
        return False, ev
    rng_ = max(means) - min(means)
    ln_n = np.log(np.asarray(sizes, dtype=np.float64))
    slope = float(np.polyfit(ln_n, np.asarray(means, dtype=np.float64), 1)[0]) if len(set(sizes)) > 1 else math.nan

    def cv_at(n: int) -> float:
        cvs = [r.homogeneity.cv for r in by_size[n] if math.isfinite(r.homogeneity.cv)]
        return float(np.mean(cvs)) if cvs else math.nan

    cv_small, cv_large = cv_at(sizes[0]), cv_at(sizes[-1])
    ev.update(size_range=rng_, size_slope=slope, size_cv_small=cv_small, size_cv_large=cv_large)
    ok = (
        clean_ok
        and rng_ <= c.size_range_tol
        and math.isfinite(slope)
        and abs(slope) <= c.size_slope_tol
        and math.isfinite(cv_small)
        and math.isfinite(cv_large)
        and cv_large <= cv_small + c.homogeneity_size_slack
    )
    return bool(ok), ev


def null_separated(
    cand: Sequence[RunEvidence],
    nulls: Mapping[NullModel, Sequence[RunEvidence]],
    cfg: Omega11Config,
) -> tuple[bool, dict[str, float]]:
    """Separacion de cada nulo de la bateria (`cfg.null_models`, mas los extra presentes). Nulo ausente o vacio: False."""
    c = cfg.certificate
    ev: dict[str, float] = {}
    cand_d = _finite_dstar(cand)
    if not cand_d:
        ev["null_missing"] = 1.0
        return False, ev
    mu_c = float(np.mean(cand_d))
    required = list(cfg.null_models) + [m for m in nulls if m not in cfg.null_models]
    ok = True
    if NullModel.RANDOM_GEOMETRIC in required:
        # RANDOM_GEOMETRIC es un control positivo (geometrico por construccion), nunca un nulo (desviacion 3).
        ev["null_random_geometric_invalid"] = 1.0
        ok = False
    for model in required:
        reps = nulls.get(model, ())
        key = f"null_{model.value}"
        if len(reps) == 0:
            ev[f"{key}_missing"] = 1.0
            ok = False
            continue
        frac_pass = sum(_clean(reps, cfg)) / len(reps)
        ev[f"{key}_frac_pass"] = frac_pass
        if frac_pass > c.null_pass_max:
            ok = False
        nd = _finite_dstar(reps)
        if nd:
            mu_n = float(np.mean(nd))
            sigma = float(np.std(nd, ddof=1)) if len(nd) >= 2 else 0.0
            gap = abs(mu_c - mu_n)
            need = max(c.null_sigma * sigma, c.dim_tol)
            ev[f"{key}_gap"] = gap
            ev[f"{key}_sigma"] = sigma
            if gap < need:
                ok = False
    return ok, ev


def _point_dimension(runs: Sequence[RunEvidence], cfg: Omega11Config) -> tuple[float, int | None]:
    ds = _finite_dstar(runs)
    if not ds:
        return math.nan, None
    d_star = float(np.mean(ds))
    cls = round(d_star)
    if cls < 1 or abs(d_star - cls) > cfg.certificate.class_tol:
        return d_star, None
    return d_star, int(cls)


def build_certificate(
    runs: Sequence[RunEvidence],
    by_size: Mapping[int, Sequence[RunEvidence]],
    nulls: Mapping[NullModel, Sequence[RunEvidence]],
    cfg: Omega11Config,
) -> GeometryCertificate:
    n = len(runs)
    frac = cfg.certificate.seed_fraction
    fields: dict[str, bool] = {}
    evidence: _Evidence = {"n_runs": n}
    if n > 0:
        per_run = [run_flags(r, cfg) for r in runs]
        for name in RUN_LEVEL_FIELDS:
            f = sum(1 for fl in per_run if fl[name]) / n
            evidence[f"frac_{name}"] = f
            fields[name] = bool(f >= frac)
        evidence["converged_fraction"] = sum(1 for fl in per_run if fl["converged"]) / n
    else:
        for name in RUN_LEVEL_FIELDS:
            fields[name] = False
        evidence["converged_fraction"] = 0.0
    seed_ok, seed_ev = seed_robust(runs, cfg)
    size_ok, size_ev = size_robust(by_size, cfg)
    null_ok, null_ev = null_separated(runs, nulls, cfg)
    fields["seed_robust"], fields["size_robust"], fields["null_separated"] = seed_ok, size_ok, null_ok
    evidence.update(seed_ev)
    evidence.update(size_ev)
    evidence.update(null_ev)

    missing: list[str] = []
    if n == 0:
        missing.append("runs")
    if "size_missing" in size_ev:
        missing.append("sizes")
    if "null_missing" in null_ev or any(k.endswith("_missing") and k.startswith("null_") for k in null_ev):
        missing.append("nulls")
    if "seed_missing" in seed_ev:
        missing.append("seeds")
    if missing:
        evidence["missing"] = ",".join(missing)

    d_star, cls = _point_dimension(runs, cfg)
    return GeometryCertificate(
        consensus_dimension=d_star,
        dimension_class=cls,
        n_runs=n,
        evidence=evidence,
        **{name: fields[name] for name in CERTIFICATE_FIELDS},
    )


def point_verdict(
    runs: Sequence[RunEvidence],
    by_size: Mapping[int, Sequence[RunEvidence]],
    nulls: Mapping[NullModel, Sequence[RunEvidence]],
    cfg: Omega11Config,
    ensemble: EnsembleSummary | None = None,
) -> PointVerdict:
    cert = build_certificate(runs, by_size, nulls, cfg)
    n = len(runs)
    frac = cfg.certificate.seed_fraction
    codes: set[FailureCode] = set()

    assessments = [assess_run(r, cfg) for r in runs]
    n_conv = sum(1 for a in assessments if a.flags["converged"])
    if n == 0 or n_conv / n < frac or (ensemble is not None and not ensemble.equilibrated):
        codes.add(FailureCode.F10)

    outcomes: Counter[str] = Counter()
    for a in assessments:
        outcomes["PASS" if a.primary is None else a.primary.value] += 1
    fractions = {k: v / n for k, v in outcomes.items()} if n else {}
    modal_fail = False
    if n:
        modal, count = max(outcomes.items(), key=lambda kv: (kv[1], kv[0]))
        if count / n < frac:
            codes.add(FailureCode.F7)
        elif modal != "PASS":
            codes.add(FailureCode(modal))
            modal_fail = True
    if not modal_fail:  # B5: con un fallo modal claro, F6/F7/F8 no desplazan al codigo primario
        if not cert.size_robust:
            codes.add(FailureCode.F6)
        if not cert.seed_robust:
            codes.add(FailureCode.F7)
        if not cert.null_separated:
            codes.add(FailureCode.F8)
    if not codes:
        # Salvaguarda: un certificado incompleto nunca queda sin codigo.
        if cert.dimension_class is None:
            codes.add(FailureCode.F4)
        for name in cert.failed_fields:
            codes.add(FIELD_CODE[name])

    ordered = tuple(c for c in FAILURE_PRECEDENCE if c in codes)
    return PointVerdict(
        certificate=cert,
        codes=ordered,
        primary=primary_code(ordered),
        verdict=finalize_verdict(cert, ordered),
        run_outcome_fractions=fractions,
    )
