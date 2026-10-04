"""Observables geometricos agregados sobre la componente gigante (D-11, D-22; M§13)."""

from __future__ import annotations

import dataclasses

import numpy as np

from omega.config.settings import DimensionConfig, GraphConfig, SpectralConfig
from omega.geometry.dimension import effective_dimension
from omega.geometry.distances import (
    diameter,
    distance_matrix,
    hop_distance_matrix,
    path_length,
    threshold_adjacency,
)
from omega.geometry.spectral import spectral_dimension
from omega.network.topology import component_labels, giant_component_nodes, submatrix
from omega.network.weights import validate_weight_matrix
from omega.types import FloatArray, GeometryObservables

__all__ = ["geometry_observables"]


def geometry_observables(
    w: FloatArray, graph: GraphConfig, dim: DimensionConfig, spec: SpectralConfig
) -> GeometryObservables:
    """L, L en saltos, diametro, D_eff (shell/ball segun dim, ball, saltos) y D_s sobre la
    componente gigante del grafo A=(W>w_min) (D-2, D-11, D-22). Con <2 nodos, L/diametro=NaN."""
    validate_weight_matrix(w)
    a = threshold_adjacency(w, graph.w_min)
    _, labels = component_labels(a)
    nodes = giant_component_nodes(labels)
    if nodes.size < 2:  # gigante trivial: sin pares; las estimaciones resultan insufficient_component
        d = np.zeros((1, 1), dtype=np.float64)
        return GeometryObservables(
            path_length=path_length(d),
            path_length_hops=path_length(d),
            diameter=diameter(d),
            d_eff=effective_dimension(d, dim),
            d_eff_ball=effective_dimension(d, dataclasses.replace(dim, estimator="ball")),
            d_eff_hops=effective_dimension(d, dim),
            d_s=spectral_dimension(np.zeros((2, 2), dtype=np.float64), spec, plateau_tol=dim.plateau_tol),
        )
    ws = submatrix(w, nodes)
    asub = a[np.ix_(nodes, nodes)]
    d = distance_matrix(ws, graph.w_min, graph.epsilon)
    dh = hop_distance_matrix(asub)
    if spec.graph == "thresholded_weighted":
        m = np.where(asub, ws, 0.0)
    elif spec.graph == "weighted":
        m = ws
    else:
        m = asub.astype(np.float64)
    return GeometryObservables(
        path_length=path_length(d),
        path_length_hops=path_length(dh),
        diameter=diameter(d),
        d_eff=effective_dimension(d, dim),
        d_eff_ball=effective_dimension(d, dataclasses.replace(dim, estimator="ball")),
        d_eff_hops=effective_dimension(dh, dim),
        d_s=spectral_dimension(m, spec, plateau_tol=dim.plateau_tol),
    )
