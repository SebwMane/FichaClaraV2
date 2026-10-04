"""Estructura local de la red (DESIGN §1.4-§1.5): homogeneidad, isotropía, localidad, anillos.

Todas las entradas (`a`, `d_hop`) YA están restringidas a la componente gigante. Son diagnósticos
sobre una métrica ya medida: no alimentan la dinámica. Referencias: escalado de cáscaras (D_i
por nodo), MDS clásico de Torgerson (1952) para el proxy de isotropía local, y el criterio de
desvío corto (≤3 saltos en G−e) como proxy de localidad (Watts & Strogatz, 1998).
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components

from omega.config.settings11 import LocalStructureConfig
from omega.contracts import AnnulusReport, HomogeneityReport, IsotropyReport, LocalityReport
from omega.geometry.distances import validate_distance_matrix
from omega.types import BoolArray, DimensionEstimate, FloatArray

__all__ = [
    "node_dimensions",
    "homogeneity",
    "ball_mds_ratio",
    "isotropy",
    "edge_detour_fraction",
    "locality",
    "annulus_connectivity",
]

_NAN = float("nan")


def _window_radii(hop_est: DimensionEstimate) -> FloatArray | None:
    """Radios enteros de la ventana ajustada de HOP, o None si la estimación no es "ok"."""
    if hop_est.status != "ok" or hop_est.window is None:
        return None
    lo, hi = hop_est.window
    return np.asarray(hop_est.scales[lo : hi + 1], dtype=np.float64)


def node_dimensions(d_hop: FloatArray, radii: FloatArray) -> FloatArray:
    """D_i = 1 + pendiente de ln S_i(r) frente a ln r sobre los radios enteros `radii`.

    S_i(r) es el número de nodos a exactamente r saltos de i. Si algún S_i(r)=0, D_i es NaN.
    """
    validate_distance_matrix(d_hop)
    r = np.asarray(radii, dtype=np.float64)
    if r.ndim != 1 or r.size < 2 or np.any(r <= 0.0) or np.any(np.diff(r) <= 0.0):
        raise ValueError("radii debe ser vector 1D creciente, positivo y con >=2 radios")
    n = d_hop.shape[0]
    shells = np.stack([np.sum(d_hop == rk, axis=1) for rk in r], axis=1).astype(np.float64)
    out = np.full(n, _NAN, dtype=np.float64)
    lr = np.log(r)
    lrc = lr - lr.mean()
    sxx = float(np.sum(lrc**2))
    valid = np.all(shells > 0.0, axis=1)
    if valid.any():
        ly = np.log(shells[valid])
        slope = (ly - ly.mean(axis=1, keepdims=True)) @ lrc / sxx
        out[valid] = 1.0 + slope
    return out


def homogeneity(d_hop: FloatArray, hop_est: DimensionEstimate, cfg: LocalStructureConfig) -> HomogeneityReport:
    """cv = σ_D/|μ_D| de las dimensiones por nodo (σ poblacional, nodos válidos).

    `ok := cv ≤ homogeneity_cv_max ∧ n_valid ≥ homogeneity_min_nodes`. Medido: toros 0, RGG3 0.14.
    """
    radii = _window_radii(hop_est)
    if radii is None or radii.size < 2:
        return HomogeneityReport(_NAN, _NAN, _NAN, 0, False)
    di = node_dimensions(d_hop, radii)
    v = di[np.isfinite(di)]
    if v.size == 0:
        return HomogeneityReport(_NAN, _NAN, _NAN, 0, False)
    mu = float(v.mean())
    sigma = float(v.std())
    cv = sigma / abs(mu) if mu != 0.0 else _NAN
    ok = bool(v.size >= cfg.homogeneity_min_nodes and cv <= cfg.homogeneity_cv_max)
    return HomogeneityReport(mu, sigma, cv, int(v.size), ok)


def ball_mds_ratio(d_hop: FloatArray, center: int, radius: int, k: int) -> float:
    """Cociente λ_k/λ_1 del MDS clásico (G = −½ J D² J) de la bola B_center(radius) en saltos.

    k es 1-indexado (k=1 da 1). Devuelve 0.0 si la bola tiene ≤ k nodos o λ_1 ≤ 0 (sin estructura
    k-dimensional); el cociente se acota a [0, 1].
    """
    validate_distance_matrix(d_hop)
    n = d_hop.shape[0]
    if isinstance(center, bool) or not 0 <= int(center) < n:
        raise ValueError(f"center fuera de rango: {center!r}")
    if radius < 0 or k < 1:
        raise ValueError("radius >= 0 y k >= 1")
    idx = np.flatnonzero(d_hop[center] <= radius)
    if idx.size <= k:
        return 0.0
    d2 = d_hop[np.ix_(idx, idx)] ** 2
    m = idx.size
    j = np.eye(m) - np.full((m, m), 1.0 / m)
    g = -0.5 * j @ d2 @ j
    ev = np.sort(np.linalg.eigvalsh(0.5 * (g + g.T)))[::-1]
    if not ev[0] > 1e-12:
        return 0.0
    return float(min(1.0, max(0.0, ev[k - 1] / ev[0])))


def isotropy(
    d_hop: FloatArray,
    hop_est: DimensionEstimate,
    k: int | None,
    cfg: LocalStructureConfig,
    rng: np.random.Generator,
) -> IsotropyReport:
    """Proxy MDS local: mediana y percentil 10 de λ_k/λ_1 en `isotropy_sources` bolas.

    Radio R = r_hi(HOP) + isotropy_radius_offset; fuentes muestreadas con `rng` sin reemplazo.
    `ok := mediana ≥ isotropy_median_min ∧ p10 ≥ isotropy_p10_min` (clase 1: trivialmente True;
    k None o HOP sin ventana: False). Medido (k=3): toro 1.00, RGG3 k12 0.70/0.58, slab 32×5×5 0.27.
    """
    validate_distance_matrix(d_hop)
    radii = _window_radii(hop_est)
    if k is None or radii is None or radii.size == 0:
        return IsotropyReport(_NAN, _NAN, 0, 0 if k is None else int(k), 0, False)
    if k < 1:
        raise ValueError("k debe ser >= 1")
    radius = int(round(float(radii[-1]))) + cfg.isotropy_radius_offset
    n = d_hop.shape[0]
    n_src = min(cfg.isotropy_sources, n)
    src = np.sort(rng.choice(n, size=n_src, replace=False))
    if k == 1:
        return IsotropyReport(1.0, 1.0, radius, 1, int(n_src), True)
    ratios = np.array([ball_mds_ratio(d_hop, int(i), radius, k) for i in src], dtype=np.float64)
    med = float(np.median(ratios))
    p10 = float(np.percentile(ratios, 10))
    ok = bool(med >= cfg.isotropy_median_min and p10 >= cfg.isotropy_p10_min)
    return IsotropyReport(med, p10, radius, int(k), int(n_src), ok)


def _check_a(a: BoolArray) -> None:
    if not isinstance(a, np.ndarray) or a.dtype != np.bool_ or a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("a debe ser ndarray bool cuadrada")
    if not np.array_equal(a, a.T) or np.any(np.diag(a)):
        raise ValueError("a debe ser simetrica con diagonal False")


def edge_detour_fraction(a: BoolArray) -> float:
    """Fracción de aristas con desvío ≤ 3 saltos en G−e: (A²)_uv > 0 ∨ (A³)_uv − k_u − k_v + 1 > 0.

    NaN si no hay aristas. Medido: toros 1.0, RGG3 k12 0.998, WS p=0.05 0.949, ER k12 0.86.
    """
    _check_a(a)
    af = a.astype(np.float64)
    iu, iv = np.nonzero(np.triu(a, 1))
    if iu.size == 0:
        return _NAN
    a2 = af @ af
    a3 = a2 @ af
    deg = af.sum(axis=1)
    has2 = a2[iu, iv] > 0.5
    has3 = (a3[iu, iv] - deg[iu] - deg[iv] + 1.0) > 0.5
    return float(np.mean(has2 | has3))


def locality(a: BoolArray, cfg: LocalStructureConfig) -> LocalityReport:
    """`ok := fracción de desvío ≥ locality_min` (0.97)."""
    f = edge_detour_fraction(a)
    return LocalityReport(f, bool(f >= cfg.locality_min))


def annulus_connectivity(
    a: BoolArray,
    d_hop: FloatArray,
    hop_est: DimensionEstimate,
    cfg: LocalStructureConfig,
    rng: np.random.Generator,
) -> AnnulusReport:
    """Conectividad de anillos A_i(r)={j: r ≤ d_ij ≤ r+1}, r ∈ [annulus_r_min, r_hi(HOP)].

    fractions[r] = fracción de fuentes (annulus_sources, con rng) cuyo anillo induce un subgrafo
    conexo (anillo vacío: no conexo). `ok := min_r fracción ≥ annulus_min`; si r_hi < r_min, False
    ("escala insuficiente", fractions vacío). Medido: 9³ 1.0; RGG3 k12 0.97; árbol y anillo 0.
    """
    _check_a(a)
    validate_distance_matrix(d_hop)
    if a.shape != d_hop.shape:
        raise ValueError("a y d_hop deben tener la misma forma")
    radii = _window_radii(hop_est)
    if radii is None or radii.size == 0:
        return AnnulusReport({}, False)
    r_hi = int(round(float(radii[-1])))
    if r_hi < cfg.annulus_r_min:
        return AnnulusReport({}, False)
    n = a.shape[0]
    n_src = min(cfg.annulus_sources, n)
    src = np.sort(rng.choice(n, size=n_src, replace=False))
    fractions: dict[int, float] = {}
    for r in range(cfg.annulus_r_min, r_hi + 1):
        good = 0
        for i in src:
            idx = np.flatnonzero((d_hop[i] >= r) & (d_hop[i] <= r + 1))
            if idx.size == 0:
                continue
            nc, _ = connected_components(csr_array(a[np.ix_(idx, idx)].astype(np.float64)), directed=False)
            good += int(nc == 1)
        fractions[r] = good / float(n_src)
    return AnnulusReport(fractions, bool(min(fractions.values()) >= cfg.annulus_min))
