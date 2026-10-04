"""Clasificacion de fases por corrida (ANALYSIS §4.1, preregistrada; M§19, M§21).

Predicados puros sobre `Observables` y `PhaseThresholds`. Los umbrales no se ajustan tras ver
resultados (M§44).
"""

from __future__ import annotations

import math

from omega.config.settings import PhaseThresholds
from omega.types import Observables, PhaseAssessment, PhaseLabel, RunStatus

__all__ = [
    "is_hyperdense",
    "is_dispersed",
    "is_fragmented",
    "is_geometric_candidate",
    "is_clustered",
    "is_connected_homogeneous",
    "phase_flags",
    "classify_run",
]

# Orden de precedencia E > A > B > F > D > C (§4.1): (etiqueta, clave de bandera).
_PRECEDENCE: tuple[tuple[PhaseLabel, str], ...] = (
    (PhaseLabel.E, "E"),
    (PhaseLabel.A, "A"),
    (PhaseLabel.B, "B"),
    (PhaseLabel.F, "F_candidate"),
    (PhaseLabel.D, "D"),
    (PhaseLabel.C, "C"),
)


def is_hyperdense(o: Observables, t: PhaseThresholds) -> bool:
    """E: densidad binaria >= rho_hyperdense o peso medio >= meanw_hyperdense."""
    top = o.topology
    return bool(top.binary_density >= t.rho_hyperdense or top.mean_weight >= t.meanw_hyperdense)


def is_dispersed(o: Observables, t: PhaseThresholds) -> bool:
    """A: fraccion gigante < g_dispersed y grado binario medio < kbin_dispersed."""
    top = o.topology
    return bool(top.giant_fraction < t.g_dispersed and top.mean_binary_degree < t.kbin_dispersed)


def is_fragmented(o: Observables, t: PhaseThresholds) -> bool:
    """B: >= 2 componentes grandes (>= large_component_frac*N) y fraccion gigante < g_connected."""
    top = o.topology
    return bool(top.n_large_components >= 2 and top.giant_fraction < t.g_connected)


def _ok(status: str) -> bool:
    return status == "ok"


def is_geometric_candidate(o: Observables, t: PhaseThresholds) -> bool:
    """F-candidata (por corrida): G>=g_geometric, cv<=cv_homogeneous, D_eff ok con plateau,
    D_s ok, |D_s-D_eff|<=ds_deff_tol y no hiperdensa."""
    top, geo = o.topology, o.geometry
    if top.giant_fraction < t.g_geometric or top.cv_strength > t.cv_homogeneous:
        return False
    if is_hyperdense(o, t):
        return False
    if not (_ok(geo.d_eff.status) and geo.d_eff.plateau and _ok(geo.d_s.status)):
        return False
    if not (math.isfinite(geo.d_eff.value) and math.isfinite(geo.d_s.value)):
        return False
    return bool(abs(geo.d_s.value - geo.d_eff.value) <= t.ds_deff_tol)


def is_clustered(o: Observables, t: PhaseThresholds) -> bool:
    """D: G>=g_connected, C_bin>=c_clustered y C_bin/rho>=c_ratio_clustered."""
    top = o.topology
    if top.giant_fraction < t.g_connected or top.clustering_binary < t.c_clustered:
        return False
    if top.binary_density <= 0.0:
        return False
    return bool(top.clustering_binary / top.binary_density >= t.c_ratio_clustered)


def is_connected_homogeneous(o: Observables, t: PhaseThresholds) -> bool:
    """C: G>=g_connected y cv<=cv_homogeneous."""
    top = o.topology
    return bool(top.giant_fraction >= t.g_connected and top.cv_strength <= t.cv_homogeneous)


def phase_flags(o: Observables, t: PhaseThresholds) -> dict[str, bool]:
    """Banderas independientes de cada predicado (sin precedencia)."""
    return {
        "E": is_hyperdense(o, t),
        "A": is_dispersed(o, t),
        "B": is_fragmented(o, t),
        "F_candidate": is_geometric_candidate(o, t),
        "D": is_clustered(o, t),
        "C": is_connected_homogeneous(o, t),
    }


def classify_run(o: Observables, status: RunStatus, t: PhaseThresholds) -> PhaseAssessment:
    """Etiqueta por precedencia E>A>B>F>D>C>U; una corrida no convergida recibe U.

    Se devuelven siempre todas las banderas, mas `converged`.
    """
    if not isinstance(status, RunStatus):
        raise TypeError("status debe ser RunStatus")
    flags = phase_flags(o, t)
    converged = status is RunStatus.CONVERGED
    flags["converged"] = converged
    label = PhaseLabel.U
    if converged:
        for lab, key in _PRECEDENCE:
            if flags[key]:
                label = lab
                break
    return PhaseAssessment(label=label, flags=flags)
