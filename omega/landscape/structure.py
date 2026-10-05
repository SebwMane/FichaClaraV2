"""Clasificacion estructural compartida de un estado W (Consejo Omega Rev. 2, R2.4; congelada).

Pasos: (1) vacio si max W <= 1e-6; (2) uniforme si cv de las entradas fuera de la diagonal < uniform_cv_max
(`CertificateThresholds`, 1e-6); (3) grafo fuerte A_s = {W > 1/2} y estrato debil {1e-6 < W <= 1/2}
(bandera halo); (4) componentes de A_s de tamano >= 2: `clique` (densidad interna 1), `solapada`
(no clique, todas sus aristas en cliques maximales de tamano >= 4 y a lo sumo 10 de esos cliques
maximales) u `otro`; (5) clase del estado. Diagnosticos sin voto: kappa_inj por componente,
tamanos, P3 abiertos del soporte {W > 1e-6}.

Interpretacion de "union de cliques" del estrato debil: union disjunta de cliques, es decir, cada
componente conexa del soporte debil es completa.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final, Literal, TypeAlias

import networkx as nx
import numpy as np

from omega.config.settings11 import CertificateThresholds
from omega.landscape.kkt import open_p3_count
from omega.network.weights import validate_weight_matrix
from omega.types import BoolArray, FloatArray

__all__ = [
    "StateClass",
    "ComponentKind",
    "CLASS_EMPTY",
    "CLASS_UNIFORM",
    "CLASS_SINGLE_CLIQUE",
    "CLASS_MULTI_CLIQUE",
    "CLASS_OVERLAPPING",
    "CLASS_OTHER",
    "ComponentInfo",
    "StructureClassification",
    "kappa_inj",
    "classify_structure",
]

StateClass: TypeAlias = Literal["vacío", "uniforme", "clique_única", "multi_clique", "cliques_solapadas", "otro"]
ComponentKind: TypeAlias = Literal["clique", "solapada", "otro"]

CLASS_EMPTY: Final = "vacío"
CLASS_UNIFORM: Final = "uniforme"
CLASS_SINGLE_CLIQUE: Final = "clique_única"
CLASS_MULTI_CLIQUE: Final = "multi_clique"
CLASS_OVERLAPPING: Final = "cliques_solapadas"
CLASS_OTHER: Final = "otro"

EMPTY_THRESHOLD: Final = 1e-6
STRONG_THRESHOLD: Final = 0.5
OVERLAP_MIN_CLIQUE: Final = 4
OVERLAP_MAX_CLIQUES: Final = 10


@dataclass(frozen=True, slots=True)
class ComponentInfo:
    """Componente de A_s con tamano >= 2: nodos, tipo, densidad interna, kappa_inj y n.o de cliques >= 4."""

    nodes: tuple[int, ...]
    size: int
    kind: ComponentKind
    density: float
    kappa_inj: float
    n_big_cliques: int


@dataclass(frozen=True, slots=True)
class StructureClassification:
    """Resultado de R2.4: clase, banderas halo y diagnosticos sin voto."""

    state_class: StateClass
    has_halo: bool
    halo_is_clique_union: bool | None
    components: tuple[ComponentInfo, ...]
    component_sizes: tuple[int, ...]
    kappa_inj: tuple[float, ...]
    open_p3: int


def kappa_inj(adj: BoolArray) -> float:
    """kappa = t_inj(K3)/t_inj(K2)^{3/2} del grafo binario `adj` (nan si n < 3 o sin aristas).

    t_inj(K2) = sum_{i!=j} A_ij / (n)_2 y t_inj(K3) = tr(A^3) / (n)_3 (A sin lazos: todo camino cerrado
    de longitud 3 es inyectivo). Vale 1 en una clique.
    """
    n = int(adj.shape[0])
    if n < 3:
        return float("nan")
    a = adj.astype(np.float64)
    t2 = float(a.sum()) / (n * (n - 1))
    if t2 <= 0.0:
        return float("nan")
    t3 = float(np.trace(a @ a @ a)) / (n * (n - 1) * (n - 2))
    return float(t3 / t2**1.5)


def _component_kind(sub: BoolArray) -> tuple[ComponentKind, int]:
    """(tipo, n.o de cliques maximales de tamano >= 4 contados; se corta al pasar de 10)."""
    n = int(sub.shape[0])
    n_edges = int(np.sum(np.triu(sub, k=1)))
    if n_edges == n * (n - 1) // 2:
        return "clique", 1
    g = nx.from_numpy_array(sub.astype(np.int64))
    covered = np.zeros_like(sub, dtype=bool)
    big = 0
    for clq in nx.find_cliques(g):
        if len(clq) >= OVERLAP_MIN_CLIQUE:
            big += 1
            if big > OVERLAP_MAX_CLIQUES:
                return "otro", big
            idx = np.array(clq)
            covered[np.ix_(idx, idx)] = True
    covered &= sub
    if big > 0 and int(np.sum(np.triu(covered, k=1))) == n_edges:
        return "solapada", big
    return "otro", big


def _is_clique_union(adj: BoolArray) -> bool:
    """True si cada componente conexa de `adj` es completa (union disjunta de cliques)."""
    s = adj.astype(np.float64)
    a2 = s @ s
    # i~k, k~j => i~j (transitividad) <=> union disjunta de cliques
    off = ~np.eye(adj.shape[0], dtype=bool)
    return bool(np.all(((a2 > 0.0) & off) <= adj))


def classify_structure(w: FloatArray) -> StructureClassification:
    """Clasifica W segun R2.4 (vacio, uniforme, clique_unica, multi_clique, cliques_solapadas, otro)."""
    validate_weight_matrix(w)
    n = int(w.shape[0])
    iu = np.triu_indices(n, k=1)
    off = w[iu]
    p3 = open_p3_count(w, EMPTY_THRESHOLD)
    weak = (w > EMPTY_THRESHOLD) & (w <= STRONG_THRESHOLD)
    has_halo = bool(np.any(weak))
    halo_union = _is_clique_union(weak) if has_halo else None
    strong = w > STRONG_THRESHOLD

    def result(cls: StateClass, comps: tuple[ComponentInfo, ...] = ()) -> StructureClassification:
        return StructureClassification(
            state_class=cls,
            has_halo=has_halo,
            halo_is_clique_union=halo_union,
            components=comps,
            component_sizes=tuple(c.size for c in comps),
            kappa_inj=tuple(c.kappa_inj for c in comps),
            open_p3=p3,
        )

    if float(np.max(off)) <= EMPTY_THRESHOLD:
        return result(CLASS_EMPTY)
    mean = float(np.mean(off))
    if mean > 0.0 and float(np.std(off)) / mean < CertificateThresholds().uniform_cv_max:
        return result(CLASS_UNIFORM)

    g = nx.from_numpy_array(strong.astype(np.int64))
    comps: list[ComponentInfo] = []
    for nodes in sorted((sorted(c) for c in nx.connected_components(g) if len(c) >= 2), key=lambda c: c[0]):
        idx = np.array(nodes)
        sub = strong[np.ix_(idx, idx)]
        kind, nbig = _component_kind(sub)
        m = len(nodes)
        dens = float(np.sum(np.triu(sub, k=1))) / (m * (m - 1) / 2.0)
        comps.append(ComponentInfo(tuple(int(i) for i in nodes), m, kind, dens, kappa_inj(sub), nbig))
    ctup = tuple(comps)
    if not ctup or any(c.kind == "otro" for c in ctup):
        return result(CLASS_OTHER, ctup)
    if all(c.kind == "clique" for c in ctup):
        return result(CLASS_SINGLE_CLIQUE if len(ctup) == 1 else CLASS_MULTI_CLIQUE, ctup)
    return result(CLASS_OVERLAPPING, ctup)
