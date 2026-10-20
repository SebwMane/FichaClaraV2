"""R19-(b): archivo sellado del guardián (interfaz: docs/R19_INTERFAZ_SELLADO.md).

Define FAMILIAS y ESTRUCTURAS. Autocontenido: solo depende de numpy.

- «se_desvanece»: entradas {"tipo": "copias", "k", "frac"}; el ejecutor implementa generador y control reinyectado.
- «autosemejante»: entradas {"tipo": "region", "generar"}; generar(n, rng) -> C (bool n×n), un sprinkling de Poisson
  de n puntos en el diamante de gauge {t + T ≥ F(x) y T − t ≥ F(−x)} (T = 0.5) con la relación
  C[i, j] ⇔ Δt > 0 y Δt ≥ F(x_j − x_i). F es sublineal (positivamente homogénea y subaditiva) y positiva fuera
  del origen, así que la relación es un orden parcial estricto. La región y la relación no dependen de n.
  Muestreo por rechazo en una caja que contiene la región; la relación se construye por bloques de filas.
"""

from __future__ import annotations

import math
from typing import Any, Callable

import numpy as np
import numpy.typing as npt

Bool = npt.NDArray[np.bool_]
Float = npt.NDArray[np.float64]
Gauge = Callable[[Float, Float, Float], Float]

T_HALF = 0.5
BLOCK_ROWS = 1024


# --- Gauges -------------------------------------------------------------------------------------------

def _tetra(dx: Float, dy: Float, dz: Float) -> Float:
    """Unidad = tetraedro: F(x) = max(x+y+z, x−y−z, −x+y−z, −x−y+z)."""
    out: Float = np.maximum(np.maximum(dx + dy + dz, dx - dy - dz), np.maximum(-dx + dy - dz, -dx - dy + dz))
    return out


_PHI = (1.0 + math.sqrt(5.0)) / 2.0
_ICO = np.array([[0.0, 1.0, _PHI], [0.0, 1.0, -_PHI], [1.0, _PHI, 0.0],
                 [1.0, -_PHI, 0.0], [_PHI, 0.0, 1.0], [-_PHI, 0.0, 1.0]]) / math.sqrt(1.0 + _PHI * _PHI)


def _dodeca(dx: Float, dy: Float, dz: Float) -> Float:
    """Unidad = dodecaedro: F(x) = max sobre los 6 ejes del icosaedro de |a·x|."""
    out: Float = np.abs(_ICO[0, 0] * dx + _ICO[0, 1] * dy + _ICO[0, 2] * dz)
    for a in _ICO[1:]:
        out = np.maximum(out, np.abs(a[0] * dx + a[1] * dy + a[2] * dz))
    return out


def _rombo(dx: Float, dy: Float, dz: Float) -> Float:
    """F(x) = max(|x|+|y|, |y|+|z|, |z|+|x|)."""
    ax, ay, az = np.abs(dx), np.abs(dy), np.abs(dz)
    out: Float = np.maximum(np.maximum(ax + ay, ay + az), az + ax)
    return out


_Q1 = np.array([1.8, 1.0, 0.55])
_Q2 = np.array([0.55, 1.8, 1.0])


def _suma_elip(dx: Float, dy: Float, dz: Float) -> Float:
    """F(x) = ½(√(xᵀQ₁x) + √(xᵀQ₂x)), Q₁, Q₂ diagonales distintas."""
    a = np.sqrt(_Q1[0] * dx * dx + _Q1[1] * dy * dy + _Q1[2] * dz * dz)
    b = np.sqrt(_Q2[0] * dx * dx + _Q2[1] * dy * dy + _Q2[2] * dz * dz)
    out: Float = 0.5 * (a + b)
    return out


_V = np.array([1.0, 2.0, -2.0]) / 3.0
_W = 0.45


def _sesgo(dx: Float, dy: Float, dz: Float) -> Float:
    """F(x) = ‖x‖₂ + w·max(0, v·x), v unitario, w = 0.45."""
    r = np.sqrt(dx * dx + dy * dy + dz * dz)
    out: Float = r + _W * np.maximum(0.0, _V[0] * dx + _V[1] * dy + _V[2] * dz)
    return out


# --- Sprinkling y relación ------------------------------------------------------------------------------

def _in_region(x: Float, g: Gauge) -> Bool:
    out: Bool = (x[:, 0] + T_HALF >= g(x[:, 1], x[:, 2], x[:, 3])) & (T_HALF - x[:, 0] >= g(-x[:, 1], -x[:, 2], -x[:, 3]))
    return out


def _box_for(g: Gauge, n_dir: int = 400_000) -> tuple[Float, Float]:
    """Caja que contiene la región: para cada eje, |x_i| ≤ 2T·max_u |u_i|/(F(u)+F(−u)) (cota exacta del máximo
    sobre direcciones), estimada en n_dir direcciones fijas y ampliada un 15%. Determinista (semilla fija interna,
    no relacionada con el análisis)."""
    u = np.random.default_rng(0).standard_normal((n_dir, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    s = g(u[:, 0], u[:, 1], u[:, 2]) + g(-u[:, 0], -u[:, 1], -u[:, 2])
    h = 1.15 * 2.0 * T_HALF * np.max(np.abs(u) / s[:, None], axis=0)
    lo = np.array([-T_HALF, -h[0], -h[1], -h[2]])
    return lo, -lo


def _sprinkle(n: int, g: Gauge, box: tuple[Float, Float], rng: np.random.Generator) -> Float:
    if n < 0:
        raise ValueError("n debe ser >= 0")
    out = np.empty((n, 4), dtype=np.float64)
    filled = 0
    while filled < n:
        batch = max(64, 32 * (n - filled))
        cand = rng.uniform(box[0], box[1], size=(batch, 4))
        ok = cand[_in_region(cand, g)]
        take = min(len(ok), n - filled)
        out[filled: filled + take] = ok[:take]
        filled += take
    return out


def _relation(x: Float, g: Gauge, block: int = BLOCK_ROWS) -> Bool:
    n = x.shape[0]
    t = x[:, 0]
    c = np.empty((n, n), dtype=np.bool_)
    for lo in range(0, n, block):
        hi = min(n, lo + block)
        dt = t[None, :] - t[lo:hi, None]
        d = [x[None, :, k] - x[lo:hi, None, k] for k in (1, 2, 3)]
        c[lo:hi] = (dt > 0) & (dt >= g(d[0], d[1], d[2]))
    return c


def _region_gen(g: Gauge) -> Callable[[int, np.random.Generator], Bool]:
    box = _box_for(g)

    def generar(n: int, rng: np.random.Generator) -> Bool:
        return _relation(_sprinkle(int(n), g, box, rng), g)

    generar.gauge = g  # type: ignore[attr-defined]
    generar.box = box  # type: ignore[attr-defined]
    return generar


# --- Tablas selladas --------------------------------------------------------------------------------------

ESTRUCTURAS: dict[str, dict[str, Any]] = {
    "S1a_gemelos_14": {"tipo": "copias", "k": 1, "frac": 0.14},
    "S1b_gemelos_22": {"tipo": "copias", "k": 1, "frac": 0.22},
    "S2a_tripletes_12": {"tipo": "copias", "k": 2, "frac": 0.12},
    "S2b_tripletes_18": {"tipo": "copias", "k": 2, "frac": 0.18},
    "A1a_cono_tetraedrico": {"tipo": "region", "generar": _region_gen(_tetra)},
    "A1b_cono_dodecaedrico": {"tipo": "region", "generar": _region_gen(_dodeca)},
    "A1c_cono_rombododecaedrico": {"tipo": "region", "generar": _region_gen(_rombo)},
    "A2a_cono_suma_elipsoides": {"tipo": "region", "generar": _region_gen(_suma_elip)},
    "A2b_cono_esfera_sesgada": {"tipo": "region", "generar": _region_gen(_sesgo)},
}

FAMILIAS: list[dict[str, Any]] = [
    {"nombre": "S1_gemelos", "tipo": "se_desvanece", "estructuras": ["S1a_gemelos_14", "S1b_gemelos_22"]},
    {"nombre": "S2_tripletes", "tipo": "se_desvanece", "estructuras": ["S2a_tripletes_12", "S2b_tripletes_18"]},
    {"nombre": "A1_conos_poliedricos", "tipo": "autosemejante",
     "estructuras": ["A1a_cono_tetraedrico", "A1b_cono_dodecaedrico", "A1c_cono_rombododecaedrico"]},
    {"nombre": "A2_conos_suaves", "tipo": "autosemejante",
     "estructuras": ["A2a_cono_suma_elipsoides", "A2b_cono_esfera_sesgada"]},
]
