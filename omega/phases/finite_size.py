"""Tamano finito (OMEGA_1_1_DESIGN §1.14, §3.7): configuracion por tamano, observables FSS y ajustes.

Observables: D*, cada D, C_bin, G/N, rho binaria, xi_F = lambda2^{-1/2} (longitud de Fiedler) y L_hop.
`xi_F ~ N^{1/D}` en estructuras geometricas y ~ constante en expansores (discriminante F3/F6).
Solo calcula; la regla de decision (`size_robust`) vive en `omega.certificate.certificate`.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Sequence

import numpy as np

from omega.config.settings11 import Omega11Config
from omega.contracts import RunEvidence
from omega.types import FloatArray

__all__ = ["size_config", "fss_observables", "fss_fit", "evidence_cfg", "FSS_KEYS", "DENSE_CURVATURE_EDGES"]

DENSE_CURVATURE_EDGES = 25  # presupuesto de aristas de Ollivier en estados densos (solo coste, ver `evidence_cfg`)

FSS_KEYS = ("d_star", "d_vol", "d_s", "d_weyl", "c_bin", "giant_fraction", "rho", "xi_f", "l_hop")


def size_config(cfg: Omega11Config, n: int) -> Omega11Config:
    """Copia de `cfg` con `base.init.n = n` (los umbrales no se tocan)."""
    if isinstance(n, bool) or not isinstance(n, int):
        raise TypeError("n debe ser entero")
    if n < 3:
        raise ValueError("n debe ser >= 3")
    base = dataclasses.replace(cfg.base, init=dataclasses.replace(cfg.base.init, n=n))
    return dataclasses.replace(cfg, base=base)


def evidence_cfg(cfg: Omega11Config, w: FloatArray) -> Omega11Config:
    """Configuracion de evidencia con presupuesto de coste para estados F1 garantizados.

    Ollivier resuelve un LP de transporte por arista con soporte ~ grado: con grafos casi completos cuesta ~90 s con
    N=120 y es inviable con N=800. El recorte se aplica SOLO cuando F1 (denso o uniforme) esta garantizado por los mismos
    umbrales del certificado (`certificate.dense_rho` / `dense_meanw`, tests de `assess_run`): en esos estados la
    curvatura es (casi) constante y se muestrean `DENSE_CURVATURE_EDGES` aristas. Estados no triviales con
    0.1 < rho < 0.5 conservan el presupuesto completo (Enmienda A-6 (auditoria B1, B2)). Solo presupuesto de muestreo
    (`curvature.max_edges`); NINGUN umbral cambia.
    """
    if not isinstance(cfg, Omega11Config):
        raise TypeError("cfg debe ser Omega11Config")
    n = w.shape[0]
    iu = np.triu_indices(n, 1)
    v = w[iu]
    dense = bool(
        np.mean(v > cfg.base.graph.w_min) >= cfg.certificate.dense_rho or float(v.mean()) >= cfg.certificate.dense_meanw
    )
    if not dense or cfg.curvature.max_edges <= DENSE_CURVATURE_EDGES:
        return cfg
    return dataclasses.replace(cfg, curvature=dataclasses.replace(cfg.curvature, max_edges=DENSE_CURVATURE_EDGES))


def fss_observables(ev: RunEvidence) -> dict[str, float]:
    """Observables de tamano finito de una corrida; los no medibles son NaN.

    Claves: d_star (consenso), d_vol, d_s, d_weyl, c_bin (clustering binario), giant_fraction (G/N), rho (densidad
    binaria), xi_f (longitud de Fiedler) y l_hop (camino medio en saltos sobre la gigante).
    """
    top, geo = ev.topology, ev.geometry
    vals = {
        "d_star": ev.consensus_dimension,
        "d_vol": geo.d_eff.value,
        "d_s": geo.d_s.value,
        "d_weyl": ev.d_weyl.value,
        "c_bin": top.clustering_binary,
        "giant_fraction": top.giant_fraction,
        "rho": top.binary_density,
        "xi_f": ev.fiedler_length,
        "l_hop": geo.path_length_hops,
    }
    return {k: float(v) if math.isfinite(v) else math.nan for k, v in vals.items()}


def _linfit(x: FloatArray, y: FloatArray) -> tuple[float, float, float, float]:
    """(pendiente, ordenada, error estandar de la pendiente, R^2) por minimos cuadrados; se exigen >= 2 puntos."""
    n = x.shape[0]
    xm, ym = float(x.mean()), float(y.mean())
    sxx = float(np.sum((x - xm) ** 2))
    if sxx <= 0.0:
        return math.nan, math.nan, math.nan, math.nan
    slope = float(np.sum((x - xm) * (y - ym)) / sxx)
    icpt = ym - slope * xm
    res = y - (icpt + slope * x)
    sst = float(np.sum((y - ym) ** 2))
    ss_res = float(np.sum(res**2))
    r2 = 1.0 - ss_res / sst if sst > 0.0 else 1.0
    se = math.sqrt(ss_res / (n - 2) / sxx) if n > 2 else math.nan
    return slope, icpt, se, r2


def fss_fit(sizes: Sequence[int], values: Sequence[float]) -> dict[str, float]:
    """Ajustes de escalado: pendiente log-log (ln v vs ln N, solo v > 0) y ajuste en ln N (v vs ln N).

    Los pares con v no finito se descartan. Devuelve `slope_loglog`, `intercept_loglog`, `se_loglog`, `r2_loglog`,
    `slope_ln`, `intercept_ln`, `se_ln`, `r2_ln` y `n_points` (los campos no estimables son NaN).
    Se exigen >= 2 pares finitos; `se_*` requiere >= 3.
    """
    if len(sizes) != len(values):
        raise ValueError("sizes y values deben tener la misma longitud")
    n_arr = np.asarray(sizes, dtype=np.float64)
    v_arr = np.asarray(values, dtype=np.float64)
    if np.any(n_arr <= 0.0) or not np.all(np.isfinite(n_arr)):
        raise ValueError("sizes debe ser finito y > 0")
    ok = np.isfinite(v_arr)
    if int(ok.sum()) < 2:
        raise ValueError("se requieren >= 2 pares finitos")
    ln_n = np.log(n_arr[ok])
    v = v_arr[ok]
    out: dict[str, float] = {"n_points": float(v.shape[0])}
    pos = v > 0.0
    if int(pos.sum()) >= 2:
        s, i, se, r2 = _linfit(ln_n[pos], np.log(v[pos]))
    else:
        s = i = se = r2 = math.nan
    out.update(slope_loglog=s, intercept_loglog=i, se_loglog=se, r2_loglog=r2)
    s, i, se, r2 = _linfit(ln_n, v)
    out.update(slope_ln=s, intercept_ln=i, se_ln=se, r2_ln=r2)
    return out
