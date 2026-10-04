"""Analisis de coarse-graining (OMEGA_1_1_DESIGN §1.7, §3.7): evidencia por nivel y consistencia de clase.

`hierarchy_evidence` aplica la regla preregistrada `heavy_edge_matching/max/v1` (`omega.coarse_graining.rule`) con
`cfg.coarse.replicates` replicas del rng y recolecta la evidencia de cada nivel medible. `class_consistency` aplica la
metrica de P§14: todos los niveles medibles con la misma `dimension_class` y sin nuevo codigo de fallo.

Solo calcula; no decide fases. Los niveles sin clase (p. ej. `no_window`) cuentan como indeterminados, no como fallo.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from omega.certificate.evidence import collect_run_evidence
from omega.certificate.taxonomy import assess_run, primary_code
from omega.coarse_graining.rule import coarse_grain_hierarchy
from omega.config.settings11 import Omega11Config
from omega.contracts import CoarseConsistency, FailureCode, RunEvidence
from omega.network.weights import validate_weight_matrix
from omega.types import FloatArray, RunStatus

__all__ = ["hierarchy_evidence", "class_consistency"]


def hierarchy_evidence(w: FloatArray, cfg: Omega11Config, rng: np.random.Generator) -> tuple[RunEvidence, ...]:
    """Evidencia del nivel 0 (el estado `w`) seguida, por cada una de las `cfg.coarse.replicates` replicas, de la
    de sus niveles 1..L (L <= max_levels, n' >= n_min). Orden: nivel 0, replica 0 (niveles 1..), replica 1, ...

    El consumo del rng es secuencial (emparejamiento y evidencia de cada nivel), asi que es determinista por rng.
    """
    validate_weight_matrix(w)
    if not isinstance(cfg, Omega11Config):
        raise TypeError("cfg debe ser Omega11Config")
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng debe ser numpy.random.Generator")
    w_min = cfg.base.graph.w_min
    out = [collect_run_evidence(w, RunStatus.CONVERGED, cfg, rng)]
    for _ in range(cfg.coarse.replicates):
        for lv in coarse_grain_hierarchy(w, w_min, cfg.coarse, rng):
            out.append(collect_run_evidence(lv.w, RunStatus.CONVERGED, cfg, rng))
    return tuple(out)


def class_consistency(levels: Sequence[RunEvidence], cfg: Omega11Config) -> CoarseConsistency:
    """Consistencia de clase entre niveles (`levels[0]` es el estado original).

    * `classes[i]` = `dimension_class` del nivel i (None = indeterminado: no cuenta como fallo).
    * `codes[i]` = codigo primario de los codigos del nivel i que NO estan en el nivel 0 (nuevo codigo de fallo); None si
      no hay ninguno o si el nivel es indeterminado. Para i = 0 siempre None.
    * `consistent` := hay >= 2 niveles medibles, todos con la misma clase y ninguno con codigo nuevo.
    """
    if len(levels) == 0:
        raise ValueError("levels vacio")
    base_codes = set(assess_run(levels[0], cfg).codes)
    classes: list[int | None] = []
    codes: list[FailureCode | None] = []
    for i, ev in enumerate(levels):
        classes.append(ev.dimension_class)
        if i == 0 or ev.dimension_class is None:
            codes.append(None)
            continue
        new = [c for c in assess_run(ev, cfg).codes if c not in base_codes]
        codes.append(primary_code(new))
    measured = [c for c in classes if c is not None]
    undetermined = len(classes) - len(measured)
    new_code = any(c is not None for c in codes)
    consistent = len(measured) >= 2 and len(set(measured)) == 1 and not new_code
    return CoarseConsistency(
        classes=tuple(classes),
        codes=tuple(codes),
        measured_levels=len(measured),
        undetermined_levels=undetermined,
        consistent=bool(consistent),
    )
