"""Metropolis por aristas con paso adaptativo solo en burn-in (DESIGN §1.8, §3.4).

Theta = Theta_hat*beta es temperatura estadistica, no tiempo fisico. Propuesta simetrica
w' = B(w + s*U(-1,1)) con reflexion B; aceptacion min(1, exp(-dH/Theta)).
"""

from __future__ import annotations

import numpy as np

from omega.config.settings import FunctionalParams
from omega.config.settings11 import Engine, MetropolisConfig
from omega.contracts import ChainResult
from omega.network.weights import validate_weight_matrix
from omega.statistics.langevin import (
    OBSERVABLE_KEYS,
    chain_observables,
    pick_state_indices,
    reflect_unit,
    temperature,
)
from omega.types import FloatArray

__all__ = ["delta_action_edge", "metropolis_sweep", "run_metropolis"]

STEP_MIN = 1e-4
STEP_MAX = 1.0
ADAPT_FACTOR = 1.1


def _require_eta0(p: FunctionalParams) -> None:
    if p.eta != 0.0:
        raise ValueError("Metropolis exige eta == 0")


def delta_action_edge(
    w: FloatArray, a: int, b: int, new_value: float, p: FunctionalParams, k: FloatArray
) -> float:
    """H(w') - H(w) al cambiar la arista (a,b) a new_value, en O(N). k = fortalezas de w."""
    _require_eta0(p)
    n = w.shape[0]
    old = float(w[a, b])
    d = float(new_value) - old
    w2_ab = float(w[a] @ w[b])
    big_k = float(k.sum())
    deg = (2.0 * d * float(k[a]) + d * d) + (2.0 * d * float(k[b]) + d * d)
    deg -= ((big_k + 2.0 * d) ** 2 - big_k**2) / n
    return float(
        -p.alpha * d * w2_ab + p.beta * (new_value * new_value - old * old) + p.gamma * deg - p.mu * d
    )


def metropolis_sweep(
    w: FloatArray, p: FunctionalParams, theta: float, step: float, rng: np.random.Generator
) -> tuple[FloatArray, int]:
    """Un barrido (una propuesta por arista, orden rng.permutation(M)). Devuelve (w', aceptadas)."""
    validate_weight_matrix(w)
    _require_eta0(p)
    if not (theta > 0.0) or not (step > 0.0):
        raise ValueError("theta y step deben ser > 0")
    n = w.shape[0]
    out = np.array(w, dtype=np.float64, copy=True)
    ia, ib = np.triu_indices(n, 1)
    m = ia.shape[0]
    order = rng.permutation(m)
    props = reflect_unit(out[ia, ib] + step * rng.uniform(-1.0, 1.0, size=m))
    logu = np.log(rng.uniform(size=m))
    k = out.sum(axis=1)
    acc = 0
    for j in order:
        a, b = int(ia[j]), int(ib[j])
        new = float(props[j])
        dh = delta_action_edge(out, a, b, new, p, k)
        if dh <= 0.0 or logu[j] < -dh / theta:
            d = new - out[a, b]
            out[a, b] = new
            out[b, a] = new
            k[a] += d
            k[b] += d
            acc += 1
    return out, acc


def run_metropolis(
    w0: FloatArray,
    p: FunctionalParams,
    cfg: MetropolisConfig,
    w_min: float,
    rng: np.random.Generator,
    n_states: int = 0,
) -> ChainResult:
    """Cadena Metropolis; el paso se adapta SOLO durante el burn-in y luego queda fijo."""
    validate_weight_matrix(w0)
    _require_eta0(p)
    n = w0.shape[0]
    m = n * (n - 1) // 2
    theta = temperature(cfg.theta_hat, p.beta)
    n_burn = int(cfg.burn_in_fraction * cfg.n_sweeps)
    step = cfg.step_init
    w = np.array(w0, dtype=np.float64, copy=True)
    sweeps = list(range(cfg.thin, cfg.n_sweeps + 1, cfg.thin))
    post = [i for i, s in enumerate(sweeps) if s > n_burn]
    keep = {post[i] for i in pick_state_indices(len(post), n_states)}
    cols: dict[str, list[float]] = {k: [] for k in OBSERVABLE_KEYS}
    states: list[FloatArray] = []
    win_acc = 0
    win_n = 0
    post_acc = 0
    post_n = 0
    idx = 0
    for s in range(1, cfg.n_sweeps + 1):
        w, acc = metropolis_sweep(w, p, theta, step, rng)
        if s <= n_burn:
            win_acc += acc
            win_n += m
            if s % cfg.adapt_every == 0:
                rate = win_acc / win_n
                if rate > cfg.acceptance_high:
                    step = min(STEP_MAX, step * ADAPT_FACTOR)
                elif rate < cfg.acceptance_low:
                    step = max(STEP_MIN, step / ADAPT_FACTOR)
                win_acc = 0
                win_n = 0
        else:
            post_acc += acc
            post_n += m
        if s % cfg.thin == 0:
            for key, val in chain_observables(w, p, w_min).items():
                cols[key].append(val)
            if idx in keep:
                states.append(w.copy())
            idx += 1
    if post_n > 0:
        rate_final = post_acc / post_n
    else:
        rate_final = win_acc / win_n if win_n > 0 else 0.0
    return ChainResult(
        engine=Engine.METROPOLIS,
        theta=theta,
        samples={k: np.asarray(v, dtype=np.float64) for k, v in cols.items()},
        sample_steps=np.asarray(sweeps, dtype=np.int64),
        w_final=w,
        states=tuple(states),
        acceptance=float(rate_final),
        step_size=float(step),
        n_steps=cfg.n_sweeps,
    )
