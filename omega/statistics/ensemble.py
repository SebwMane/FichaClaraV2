"""Diagnosticos de ensemble MCMC (DESIGN §1.8): tau_int de Sokal, ESS, split-R̂, Geweke.

Theta es temperatura estadistica; aqui solo se analizan series, sin tiempo fisico.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np

from omega.config.settings11 import EnsembleConfig
from omega.contracts import ChainResult, EnsembleSummary
from omega.types import FloatArray

__all__ = [
    "autocorrelation",
    "integrated_autocorr_time",
    "effective_sample_size",
    "split_rhat",
    "geweke_z",
    "post_burn_in",
    "ensemble_summary",
]

PRIMARY_KEYS = ("action", "mean_weight")


def _as1d(x: FloatArray) -> FloatArray:
    a = np.asarray(x, dtype=np.float64)
    if a.ndim != 1 or a.shape[0] < 1:
        raise ValueError("la serie debe ser 1-D no vacia")
    return a


def autocorrelation(x: FloatArray) -> FloatArray:
    """rho_t (t=0..n-1) por FFT, estimador sesgado (divide por n); rho_0 = 1.

    Serie constante: rho = (1, 0, ..., 0).
    """
    a = _as1d(x)
    n = a.shape[0]
    d = a - a.mean()
    nfft = 1 << (2 * n - 1).bit_length()
    f = np.fft.rfft(d, nfft)
    acov = np.fft.irfft(f * np.conj(f), nfft)[:n] / n
    if acov[0] <= 0.0 or not np.isfinite(acov[0]):
        out = np.zeros(n, dtype=np.float64)
        out[0] = 1.0
        return out
    return np.asarray(acov / acov[0], dtype=np.float64)


def integrated_autocorr_time(x: FloatArray, c: float = 5.0) -> float:
    """tau = 1/2 + sum_{t=1}^M rho_t con la ventana de Sokal: menor M >= c*tau(M)."""
    rho = autocorrelation(x)
    tau_m = 0.5 + np.cumsum(rho[1:]) if rho.shape[0] > 1 else np.zeros(0)
    if tau_m.shape[0] == 0:
        return 0.5
    m = np.arange(1, tau_m.shape[0] + 1)
    ok = np.nonzero(m >= c * tau_m)[0]
    j = int(ok[0]) if ok.shape[0] > 0 else tau_m.shape[0] - 1
    return float(max(tau_m[j], 0.5))


def effective_sample_size(x: FloatArray, c: float = 5.0) -> float:
    """ESS = n / (2 tau_int)."""
    a = _as1d(x)
    return float(a.shape[0] / (2.0 * integrated_autocorr_time(a, c)))


def split_rhat(chains: Sequence[FloatArray]) -> float:
    """R̂ dividido (Gelman-Rubin con cada cadena partida en dos mitades; trunca a la longitud minima)."""
    arrs = [_as1d(c) for c in chains]
    if len(arrs) < 1:
        raise ValueError("se requiere al menos una cadena")
    n = min(a.shape[0] for a in arrs) // 2
    if n < 2:
        raise ValueError("cadenas demasiado cortas para split-R̂")
    parts = np.array([h for a in arrs for h in (a[:n], a[n : 2 * n])])
    means = parts.mean(axis=1)
    w = float(parts.var(axis=1, ddof=1).mean())
    b = float(n * means.var(ddof=1))
    if w <= 0.0:
        return 1.0 if b <= 0.0 else math.inf
    var_plus = (n - 1) / n * w + b / n
    return float(math.sqrt(var_plus / w))


def _mean_var(seg: FloatArray, c: float) -> tuple[float, float]:
    tau = integrated_autocorr_time(seg, c)
    return float(seg.mean()), float(seg.var() * 2.0 * tau / seg.shape[0])


def geweke_z(x: FloatArray, first: float = 0.1, last: float = 0.5, c: float = 5.0) -> float:
    """z de Geweke: primer `first` frente al ultimo `last` (varianza de la media con tau_int)."""
    a = _as1d(x)
    n = a.shape[0]
    na, nb = int(first * n), int(last * n)
    if na < 2 or nb < 2 or na + nb > n:
        raise ValueError("serie demasiado corta o fracciones solapadas")
    ma, va = _mean_var(a[:na], c)
    mb, vb = _mean_var(a[n - nb :], c)
    den = math.sqrt(va + vb)
    if den == 0.0:
        return 0.0 if ma == mb else math.inf
    return float((ma - mb) / den)


def post_burn_in(x: FloatArray, fraction: float) -> FloatArray:
    """Copia de la serie sin la fraccion inicial de burn-in."""
    a = _as1d(x)
    if not (0.0 <= fraction < 1.0):
        raise ValueError("fraction debe estar en [0,1)")
    return np.array(a[int(fraction * a.shape[0]) :], dtype=np.float64, copy=True)


def ensemble_summary(
    chains: Sequence[ChainResult],
    keys: Sequence[str],
    cfg: EnsembleConfig,
    burn_in_fraction: float,
) -> EnsembleSummary:
    """Resume el ensemble por clave: media, SE (via ESS), R̂, tau_int medio y ESS total.

    equilibrated := R̂ <= rhat_max y ESS_total >= ess_min en las claves primarias
    (action y mean_weight si estan en `keys`; si no, todas) y |z_Geweke| <= geweke_z_max
    en cada cadena para esas claves. geweke_ok recoge solo el criterio de Geweke.
    """
    if len(chains) < 2:
        raise ValueError("se requieren al menos 2 cadenas")
    if len(keys) == 0:
        raise ValueError("keys vacio")
    means: dict[str, float] = {}
    ses: dict[str, float] = {}
    rhat: dict[str, float] = {}
    tau: dict[str, float] = {}
    ess: dict[str, float] = {}
    primary = [k for k in keys if k in PRIMARY_KEYS] or list(keys)
    geweke_ok = True
    rhat_ok = True
    ess_ok = True
    for key in keys:
        series = [post_burn_in(c.samples[key], burn_in_fraction) for c in chains]
        taus = [integrated_autocorr_time(s, cfg.autocorr_c) for s in series]
        ess_c = [s.shape[0] / (2.0 * t) for s, t in zip(series, taus, strict=True)]
        pooled = np.concatenate(series)
        ess_tot = float(sum(ess_c))
        means[key] = float(pooled.mean())
        ses[key] = float(pooled.std() / math.sqrt(ess_tot)) if ess_tot > 0 else math.inf
        r = split_rhat(series)
        rhat[key] = r
        tau[key] = float(np.mean(taus))
        ess[key] = ess_tot
        if key in primary:
            rhat_ok = rhat_ok and r <= cfg.rhat_max
            ess_ok = ess_ok and ess_tot >= cfg.ess_min
            for s in series:
                try:
                    z = geweke_z(s, c=cfg.autocorr_c)
                except ValueError:
                    z = math.inf
                geweke_ok = geweke_ok and abs(z) <= cfg.geweke_z_max
    return EnsembleSummary(
        means=means,
        ses=ses,
        rhat=rhat,
        tau_int=tau,
        ess=ess,
        geweke_ok=bool(geweke_ok),
        equilibrated=bool(rhat_ok and ess_ok and geweke_ok),
    )
