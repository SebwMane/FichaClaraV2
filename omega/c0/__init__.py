"""Omega-C0: hipotesis de competencia relacional (prerregistro docs/OMEGA_C0_PRERREGISTRO.md)."""

from __future__ import annotations

from omega.c0.dynamics import evolve_c0
from omega.c0.functional import (
    C0Params,
    action_and_grad,
    action_c0,
    dpsi,
    grad_c0,
    kkt_box,
    params_from_targets,
    psi,
)
from omega.c0.locality import (
    classify_c0,
    degree_preserving_rewire,
    giant_component,
    h_null,
    short_cycle_fraction,
    strong_support,
    validate_locality,
)
from omega.c0.references import optimal_amplitude

__all__ = [
    "C0Params",
    "action_and_grad",
    "action_c0",
    "classify_c0",
    "degree_preserving_rewire",
    "dpsi",
    "evolve_c0",
    "giant_component",
    "grad_c0",
    "h_null",
    "kkt_box",
    "optimal_amplitude",
    "params_from_targets",
    "psi",
    "short_cycle_fraction",
    "strong_support",
    "validate_locality",
]
