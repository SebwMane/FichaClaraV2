"""R20: archivo sellado del guardián (interfaz: docs/R19_INTERFAZ_SELLADO.md más el agregado IDENTIFICABLES de
docs/OMEGA_R20_REPETICION_PREREGISTRO.md §8).

Autocontenido: solo depende de numpy.

- «se_desvanece»: entradas {"tipo": "copias", "k", "frac"}; el ejecutor implementa el generador y el control.
- «autosemejante»: entradas {"tipo": "region", "generar"}; generar(n, rng) -> C (bool n×n), un sprinkling de Poisson
  de n puntos en el diamante de gauge {t + T ≥ F(x) y T − t ≥ F(−x)} (T = 0.5) con la relación
  C[i, j] ⇔ Δt > 0 y Δt ≥ F(x_j − x_i). F es el gauge (funcional de Minkowski) de un cuerpo convexo con el origen
  en su interior: positivamente homogéneo, subaditivo y positivo fuera del origen, así que la relación es un orden
  parcial estricto. La región y la relación no dependen de n. Muestreo por rechazo en una caja fija que contiene la
  región; la relación se construye por bloques de filas.
- IDENTIFICABLES: requisitos de visibilidad de §3, calculados por el guardián con la entropía 20261114 (ver el
  recibo). Solo votan las estructuras con "identificable": True.
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


# --- Gauges poliédricos: F(x) = max_i a_i·x (caras {a_i·x ≤ 1}) ------------------------------------------

def _poly(normals: Float) -> Gauge:
    a = np.asarray(normals, dtype=np.float64)

    def g(dx: Float, dy: Float, dz: Float) -> Float:
        out: Float = a[0, 0] * dx + a[0, 1] * dy + a[0, 2] * dz
        for v in a[1:]:
            out = np.maximum(out, v[0] * dx + v[1] * dy + v[2] * dz)
        return out

    return g


def _ngon_normals(m: int) -> Float:
    """Normales unitarias exteriores de las aristas del m-gono regular en el plano xy (circunradio 1)."""
    ang = 2.0 * math.pi * (np.arange(m) + 0.5) / m
    return np.stack([np.cos(ang), np.sin(ang)], axis=1)


_TRI = _ngon_normals(3)          # apotema 1/2
_HEX = _ngon_normals(6)          # apotema √3/2

# Prisma triangular: triángulo equilátero (circunradio 1) × [−1, 1].
_N_PRISMA3 = np.vstack([np.column_stack([_TRI / 0.5, np.zeros(3)]), [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]])

# Prisma hexagonal: hexágono regular (circunradio 1) × [−1, 1].
_N_PRISMA6 = np.vstack([np.column_stack([_HEX / (math.sqrt(3.0) / 2.0), np.zeros(6)]),
                        [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]])

# Pirámide de base cuadrada: base [−1, 1]² en z = −h0, ápice (0, 0, h1); h0 = 0.5, h1 = 1.5.
_H0, _H1 = 0.5, 1.5
_B = (_H1 + _H0) / _H1
_N_PIRAMIDE4 = np.array([[0.0, 0.0, -1.0 / _H0],
                         [_B, 0.0, 1.0 / _H1], [-_B, 0.0, 1.0 / _H1],
                         [0.0, _B, 1.0 / _H1], [0.0, -_B, 1.0 / _H1]])

# Bipirámide triangular: triángulo equilátero (circunradio 1) en z = 0, ápices (0, 0, ±1.2).
_HB = 1.2
_N_BIPIR3 = np.vstack([np.column_stack([_TRI / 0.5, np.full(3, 1.0 / _HB)]),
                       np.column_stack([_TRI / 0.5, np.full(3, -1.0 / _HB)])])

_prisma3 = _poly(_N_PRISMA3)
_prisma6 = _poly(_N_PRISMA6)
_piramide4 = _poly(_N_PIRAMIDE4)
_bipir3 = _poly(_N_BIPIR3)


# --- Gauges curvos (cuerpos con frontera curva, no elipsoidales) ---------------------------------------------

def _disco_desplazado(px: Float, py: Float, cx: float, cy: float, rad: float) -> Float:
    """Gauge 2D del disco de radio rad con centro (cx, cy), |c| < rad."""
    s = rad * rad - (cx * cx + cy * cy)
    pc = px * cx + py * cy
    out: Float = (pc + np.sqrt(pc * pc + (px * px + py * py) * s)) / s
    return out


# Cono circular: base disco de radio 1 en z = −h0 y ápice (0, 0, h1); h0 = 0.5, h1 = 1.5.
def _cono_circular(dx: Float, dy: Float, dz: Float) -> Float:
    rho = np.sqrt(dx * dx + dy * dy)
    out: Float = np.maximum(-dz / _H0, (rho * (_H0 + _H1) + dz) / _H1)
    return out


# Prisma de Reuleaux: triángulo de Reuleaux (discos de radio √3 centrados en los vértices del triángulo equilátero
# de circunradio 1) × [−1, 1].
_RV = np.stack([np.cos(2.0 * math.pi * np.arange(3) / 3 + math.pi / 2), np.sin(2.0 * math.pi * np.arange(3) / 3 + math.pi / 2)], axis=1)
_RR = math.sqrt(3.0)


def _prisma_reuleaux(dx: Float, dy: Float, dz: Float) -> Float:
    out: Float = np.abs(dz)
    for v in _RV:
        out = np.maximum(out, _disco_desplazado(dx, dy, float(v[0]), float(v[1]), _RR))
    return out


# Prisma semicircular: semidisco {x² + y² ≤ 1, y ≥ 0} trasladado por (0, −0.42) × [−1, 1].
_YS = 0.42


def _prisma_semicircular(dx: Float, dy: Float, dz: Float) -> Float:
    out: Float = np.maximum(np.maximum(_disco_desplazado(dx, dy, 0.0, -_YS, 1.0), -dy / _YS), np.abs(dz))
    return out


# Prisma lenticular: intersección de los discos de radio 1 con centros (±0.6, 0) × [−1, 1].
_CL = 0.6


def _prisma_lenticular(dx: Float, dy: Float, dz: Float) -> Float:
    out: Float = np.maximum(np.maximum(_disco_desplazado(dx, dy, _CL, 0.0, 1.0), _disco_desplazado(dx, dy, -_CL, 0.0, 1.0)),
                            np.abs(dz))
    return out


# --- Sprinkling y relación ------------------------------------------------------------------------------------

def _in_region(x: Float, g: Gauge) -> Bool:
    out: Bool = (x[:, 0] + T_HALF >= g(x[:, 1], x[:, 2], x[:, 3])) & (T_HALF - x[:, 0] >= g(-x[:, 1], -x[:, 2], -x[:, 3]))
    return out


def _box_for(g: Gauge, n_dir: int = 400_000) -> tuple[Float, Float]:
    """Caja que contiene la región. Un punto (t, x) de la región cumple F(x) ≤ t + T y F(−x) ≤ T − t, luego
    F(x) + F(−x) ≤ 2T y |x_i| ≤ 2T·max_u |u_i|/(F(u) + F(−u)). El máximo se estima en n_dir direcciones fijas
    (semilla interna fija, ajena al análisis) y se amplía un 15%."""
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


# --- Tablas selladas ---------------------------------------------------------------------------------------------

ESTRUCTURAS: dict[str, dict[str, Any]] = {
    "G1a_gemelos_07": {"tipo": "copias", "k": 1, "frac": 0.07},
    "G1b_gemelos_12": {"tipo": "copias", "k": 1, "frac": 0.12},
    "G1c_gemelos_17": {"tipo": "copias", "k": 1, "frac": 0.17},
    "G1d_gemelos_23": {"tipo": "copias", "k": 1, "frac": 0.23},
    "G2a_tripletes_07": {"tipo": "copias", "k": 2, "frac": 0.07},
    "G2b_tripletes_12": {"tipo": "copias", "k": 2, "frac": 0.12},
    "G2c_tripletes_17": {"tipo": "copias", "k": 2, "frac": 0.17},
    "G2d_tripletes_23": {"tipo": "copias", "k": 2, "frac": 0.23},
    "P1a_cono_prisma_triangular": {"tipo": "region", "generar": _region_gen(_prisma3)},
    "P1b_cono_prisma_hexagonal": {"tipo": "region", "generar": _region_gen(_prisma6)},
    "P1c_cono_piramide_cuadrada": {"tipo": "region", "generar": _region_gen(_piramide4)},
    "P1d_cono_bipiramide_triangular": {"tipo": "region", "generar": _region_gen(_bipir3)},
    "S2a_cono_circular": {"tipo": "region", "generar": _region_gen(_cono_circular)},
    "S2b_cono_prisma_reuleaux": {"tipo": "region", "generar": _region_gen(_prisma_reuleaux)},
    "S2c_cono_prisma_semicircular": {"tipo": "region", "generar": _region_gen(_prisma_semicircular)},
    "S2d_cono_prisma_lenticular": {"tipo": "region", "generar": _region_gen(_prisma_lenticular)},
}

FAMILIAS: list[dict[str, Any]] = [
    {"nombre": "G1_gemelos", "tipo": "se_desvanece",
     "estructuras": ["G1a_gemelos_07", "G1b_gemelos_12", "G1c_gemelos_17", "G1d_gemelos_23"]},
    {"nombre": "G2_tripletes", "tipo": "se_desvanece",
     "estructuras": ["G2a_tripletes_07", "G2b_tripletes_12", "G2c_tripletes_17", "G2d_tripletes_23"]},
    {"nombre": "P1_conos_poliedricos", "tipo": "autosemejante",
     "estructuras": ["P1a_cono_prisma_triangular", "P1b_cono_prisma_hexagonal", "P1c_cono_piramide_cuadrada",
                     "P1d_cono_bipiramide_triangular"]},
    {"nombre": "S2_conos_suaves", "tipo": "autosemejante",
     "estructuras": ["S2a_cono_circular", "S2b_cono_prisma_reuleaux", "S2c_cono_prisma_semicircular",
                     "S2d_cono_prisma_lenticular"]},
]

# --- IDENTIFICABLES (BEGIN) ---
# Bandas de selección (§2, entropía 20261114): radio conforme (96.º de 100 de calibración) q=1: 4.8232; q=1/8: 4.4397.
IDENTIFICABLES: dict[str, dict[str, Any]] = {'G1a_gemelos_07': {'identificable': False,
                    'familia': 'G1_gemelos',
                    'regla': 'copias: med D_obs/radio(q=1) >= 2.0 y med D_ctrl/radio(q=1/8) >= 1.5',
                    'medidas': {'med_Dobs_radio_q1': 3.5655,
                                'med_Dctrl_radio_q1': 3.5471,
                                'med_Dobs_radio_q1_8': 0.7239,
                                'med_Dctrl_radio_q1_8': 1.0553},
                    'radio': {'q1': 4.8232, 'q1_8': 4.4397},
                    'semillas': 'SeedSequence(20261114, spawn_key=(333, 9, 0, real, sub[, qi])), real = 0..9; sub 0 '
                                'puntos, 1 copias del observado, 2 dilución, 3 intervalos del observado, 4 '
                                'subconjunto de base del control, 5 puntos del control (región), 6 copias del '
                                'control, 7 intervalos del control; qi = 0 (q=1), 3 (q=1/8). Banda: '
                                'SeedSequence(20261114, spawn_key=(333, 8, qi, s, 0|1)), s = 0..99 ajuste, 100..199 '
                                'calibración',
                    'N': {'N0': 8000, 'q1': 8000, 'q1_8_mediana': 1011.5, 'q1_8_rango': [974, 1048]}},
 'G1b_gemelos_12': {'identificable': True,
                    'familia': 'G1_gemelos',
                    'regla': 'copias: med D_obs/radio(q=1) >= 2.0 y med D_ctrl/radio(q=1/8) >= 1.5',
                    'medidas': {'med_Dobs_radio_q1': 6.2243,
                                'med_Dctrl_radio_q1': 6.0896,
                                'med_Dobs_radio_q1_8': 0.7217,
                                'med_Dctrl_radio_q1_8': 1.9714},
                    'radio': {'q1': 4.8232, 'q1_8': 4.4397},
                    'semillas': 'SeedSequence(20261114, spawn_key=(333, 9, 1, real, sub[, qi])), real = 0..9; sub 0 '
                                'puntos, 1 copias del observado, 2 dilución, 3 intervalos del observado, 4 '
                                'subconjunto de base del control, 5 puntos del control (región), 6 copias del '
                                'control, 7 intervalos del control; qi = 0 (q=1), 3 (q=1/8). Banda: '
                                'SeedSequence(20261114, spawn_key=(333, 8, qi, s, 0|1)), s = 0..99 ajuste, 100..199 '
                                'calibración',
                    'N': {'N0': 8000, 'q1': 8000, 'q1_8_mediana': 1005.0, 'q1_8_rango': [975, 1034]}},
 'G1c_gemelos_17': {'identificable': True,
                    'familia': 'G1_gemelos',
                    'regla': 'copias: med D_obs/radio(q=1) >= 2.0 y med D_ctrl/radio(q=1/8) >= 1.5',
                    'medidas': {'med_Dobs_radio_q1': 9.1244,
                                'med_Dctrl_radio_q1': 9.3471,
                                'med_Dobs_radio_q1_8': 0.8658,
                                'med_Dctrl_radio_q1_8': 2.7817},
                    'radio': {'q1': 4.8232, 'q1_8': 4.4397},
                    'semillas': 'SeedSequence(20261114, spawn_key=(333, 9, 2, real, sub[, qi])), real = 0..9; sub 0 '
                                'puntos, 1 copias del observado, 2 dilución, 3 intervalos del observado, 4 '
                                'subconjunto de base del control, 5 puntos del control (región), 6 copias del '
                                'control, 7 intervalos del control; qi = 0 (q=1), 3 (q=1/8). Banda: '
                                'SeedSequence(20261114, spawn_key=(333, 8, qi, s, 0|1)), s = 0..99 ajuste, 100..199 '
                                'calibración',
                    'N': {'N0': 8000, 'q1': 8000, 'q1_8_mediana': 1002.0, 'q1_8_rango': [960, 1051]}},
 'G1d_gemelos_23': {'identificable': True,
                    'familia': 'G1_gemelos',
                    'regla': 'copias: med D_obs/radio(q=1) >= 2.0 y med D_ctrl/radio(q=1/8) >= 1.5',
                    'medidas': {'med_Dobs_radio_q1': 14.064,
                                'med_Dctrl_radio_q1': 13.8746,
                                'med_Dobs_radio_q1_8': 0.9012,
                                'med_Dctrl_radio_q1_8': 3.8533},
                    'radio': {'q1': 4.8232, 'q1_8': 4.4397},
                    'semillas': 'SeedSequence(20261114, spawn_key=(333, 9, 3, real, sub[, qi])), real = 0..9; sub 0 '
                                'puntos, 1 copias del observado, 2 dilución, 3 intervalos del observado, 4 '
                                'subconjunto de base del control, 5 puntos del control (región), 6 copias del '
                                'control, 7 intervalos del control; qi = 0 (q=1), 3 (q=1/8). Banda: '
                                'SeedSequence(20261114, spawn_key=(333, 8, qi, s, 0|1)), s = 0..99 ajuste, 100..199 '
                                'calibración',
                    'N': {'N0': 8000, 'q1': 8000, 'q1_8_mediana': 996.5, 'q1_8_rango': [956, 1052]}},
 'G2a_tripletes_07': {'identificable': False,
                      'familia': 'G2_tripletes',
                      'regla': 'copias: med D_obs/radio(q=1) >= 2.0 y med D_ctrl/radio(q=1/8) >= 1.5',
                      'medidas': {'med_Dobs_radio_q1': 3.5376,
                                  'med_Dctrl_radio_q1': 3.4392,
                                  'med_Dobs_radio_q1_8': 0.7546,
                                  'med_Dctrl_radio_q1_8': 1.4403},
                      'radio': {'q1': 4.8232, 'q1_8': 4.4397},
                      'semillas': 'SeedSequence(20261114, spawn_key=(333, 9, 4, real, sub[, qi])), real = 0..9; sub '
                                  '0 puntos, 1 copias del observado, 2 dilución, 3 intervalos del observado, 4 '
                                  'subconjunto de base del control, 5 puntos del control (región), 6 copias del '
                                  'control, 7 intervalos del control; qi = 0 (q=1), 3 (q=1/8). Banda: '
                                  'SeedSequence(20261114, spawn_key=(333, 8, qi, s, 0|1)), s = 0..99 ajuste, '
                                  '100..199 calibración',
                      'N': {'N0': 8000, 'q1': 8000, 'q1_8_mediana': 993.0, 'q1_8_rango': [947, 1038]}},
 'G2b_tripletes_12': {'identificable': True,
                      'familia': 'G2_tripletes',
                      'regla': 'copias: med D_obs/radio(q=1) >= 2.0 y med D_ctrl/radio(q=1/8) >= 1.5',
                      'medidas': {'med_Dobs_radio_q1': 5.8477,
                                  'med_Dctrl_radio_q1': 5.7785,
                                  'med_Dobs_radio_q1_8': 0.8046,
                                  'med_Dctrl_radio_q1_8': 1.7049},
                      'radio': {'q1': 4.8232, 'q1_8': 4.4397},
                      'semillas': 'SeedSequence(20261114, spawn_key=(333, 9, 5, real, sub[, qi])), real = 0..9; sub '
                                  '0 puntos, 1 copias del observado, 2 dilución, 3 intervalos del observado, 4 '
                                  'subconjunto de base del control, 5 puntos del control (región), 6 copias del '
                                  'control, 7 intervalos del control; qi = 0 (q=1), 3 (q=1/8). Banda: '
                                  'SeedSequence(20261114, spawn_key=(333, 8, qi, s, 0|1)), s = 0..99 ajuste, '
                                  '100..199 calibración',
                      'N': {'N0': 8000, 'q1': 8000, 'q1_8_mediana': 993.5, 'q1_8_rango': [971, 1030]}},
 'G2c_tripletes_17': {'identificable': True,
                      'familia': 'G2_tripletes',
                      'regla': 'copias: med D_obs/radio(q=1) >= 2.0 y med D_ctrl/radio(q=1/8) >= 1.5',
                      'medidas': {'med_Dobs_radio_q1': 8.4554,
                                  'med_Dctrl_radio_q1': 8.135,
                                  'med_Dobs_radio_q1_8': 0.8829,
                                  'med_Dctrl_radio_q1_8': 2.7013},
                      'radio': {'q1': 4.8232, 'q1_8': 4.4397},
                      'semillas': 'SeedSequence(20261114, spawn_key=(333, 9, 6, real, sub[, qi])), real = 0..9; sub '
                                  '0 puntos, 1 copias del observado, 2 dilución, 3 intervalos del observado, 4 '
                                  'subconjunto de base del control, 5 puntos del control (región), 6 copias del '
                                  'control, 7 intervalos del control; qi = 0 (q=1), 3 (q=1/8). Banda: '
                                  'SeedSequence(20261114, spawn_key=(333, 8, qi, s, 0|1)), s = 0..99 ajuste, '
                                  '100..199 calibración',
                      'N': {'N0': 8000, 'q1': 8000, 'q1_8_mediana': 993.0, 'q1_8_rango': [947, 1018]}},
 'G2d_tripletes_23': {'identificable': True,
                      'familia': 'G2_tripletes',
                      'regla': 'copias: med D_obs/radio(q=1) >= 2.0 y med D_ctrl/radio(q=1/8) >= 1.5',
                      'medidas': {'med_Dobs_radio_q1': 11.6121,
                                  'med_Dctrl_radio_q1': 11.7423,
                                  'med_Dobs_radio_q1_8': 0.9708,
                                  'med_Dctrl_radio_q1_8': 3.7573},
                      'radio': {'q1': 4.8232, 'q1_8': 4.4397},
                      'semillas': 'SeedSequence(20261114, spawn_key=(333, 9, 7, real, sub[, qi])), real = 0..9; sub '
                                  '0 puntos, 1 copias del observado, 2 dilución, 3 intervalos del observado, 4 '
                                  'subconjunto de base del control, 5 puntos del control (región), 6 copias del '
                                  'control, 7 intervalos del control; qi = 0 (q=1), 3 (q=1/8). Banda: '
                                  'SeedSequence(20261114, spawn_key=(333, 8, qi, s, 0|1)), s = 0..99 ajuste, '
                                  '100..199 calibración',
                      'N': {'N0': 8000, 'q1': 8000, 'q1_8_mediana': 1001.0, 'q1_8_rango': [948, 1065]}},
 'P1a_cono_prisma_triangular': {'identificable': True,
                                'familia': 'P1_conos_poliedricos',
                                'regla': 'region: med D_obs/radio(q=1) >= 2.0 y, con q=1/8, med D_obs/radio >= 1.25 '
                                         'y med D_ctrl/radio >= 1.25',
                                'medidas': {'med_Dobs_radio_q1': 17.3606,
                                            'med_Dctrl_radio_q1': 17.39,
                                            'med_Dobs_radio_q1_8': 2.9665,
                                            'med_Dctrl_radio_q1_8': 2.9974},
                                'radio': {'q1': 4.8232, 'q1_8': 4.4397},
                                'semillas': 'SeedSequence(20261114, spawn_key=(333, 9, 8, real, sub[, qi])), real = '
                                            '0..9; sub 0 puntos, 1 copias del observado, 2 dilución, 3 intervalos '
                                            'del observado, 4 subconjunto de base del control, 5 puntos del control '
                                            '(región), 6 copias del control, 7 intervalos del control; qi = 0 (q=1), '
                                            '3 (q=1/8). Banda: SeedSequence(20261114, spawn_key=(333, 8, qi, s, '
                                            '0|1)), s = 0..99 ajuste, 100..199 calibración',
                                'N': {'N0': 8000, 'q1': 8000, 'q1_8_mediana': 1010.5, 'q1_8_rango': [984, 1057]}},
 'P1b_cono_prisma_hexagonal': {'identificable': True,
                               'familia': 'P1_conos_poliedricos',
                               'regla': 'region: med D_obs/radio(q=1) >= 2.0 y, con q=1/8, med D_obs/radio >= 1.25 y '
                                        'med D_ctrl/radio >= 1.25',
                               'medidas': {'med_Dobs_radio_q1': 10.8972,
                                           'med_Dctrl_radio_q1': 10.7839,
                                           'med_Dobs_radio_q1_8': 1.5957,
                                           'med_Dctrl_radio_q1_8': 1.6548},
                               'radio': {'q1': 4.8232, 'q1_8': 4.4397},
                               'semillas': 'SeedSequence(20261114, spawn_key=(333, 9, 9, real, sub[, qi])), real = '
                                           '0..9; sub 0 puntos, 1 copias del observado, 2 dilución, 3 intervalos del '
                                           'observado, 4 subconjunto de base del control, 5 puntos del control '
                                           '(región), 6 copias del control, 7 intervalos del control; qi = 0 (q=1), '
                                           '3 (q=1/8). Banda: SeedSequence(20261114, spawn_key=(333, 8, qi, s, '
                                           '0|1)), s = 0..99 ajuste, 100..199 calibración',
                               'N': {'N0': 8000, 'q1': 8000, 'q1_8_mediana': 1002.5, 'q1_8_rango': [970, 1041]}},
 'P1c_cono_piramide_cuadrada': {'identificable': True,
                                'familia': 'P1_conos_poliedricos',
                                'regla': 'region: med D_obs/radio(q=1) >= 2.0 y, con q=1/8, med D_obs/radio >= 1.25 '
                                         'y med D_ctrl/radio >= 1.25',
                                'medidas': {'med_Dobs_radio_q1': 18.3047,
                                            'med_Dctrl_radio_q1': 18.0707,
                                            'med_Dobs_radio_q1_8': 2.9735,
                                            'med_Dctrl_radio_q1_8': 2.9791},
                                'radio': {'q1': 4.8232, 'q1_8': 4.4397},
                                'semillas': 'SeedSequence(20261114, spawn_key=(333, 9, 10, real, sub[, qi])), real = '
                                            '0..9; sub 0 puntos, 1 copias del observado, 2 dilución, 3 intervalos '
                                            'del observado, 4 subconjunto de base del control, 5 puntos del control '
                                            '(región), 6 copias del control, 7 intervalos del control; qi = 0 (q=1), '
                                            '3 (q=1/8). Banda: SeedSequence(20261114, spawn_key=(333, 8, qi, s, '
                                            '0|1)), s = 0..99 ajuste, 100..199 calibración',
                                'N': {'N0': 8000, 'q1': 8000, 'q1_8_mediana': 997.0, 'q1_8_rango': [976, 1051]}},
 'P1d_cono_bipiramide_triangular': {'identificable': True,
                                    'familia': 'P1_conos_poliedricos',
                                    'regla': 'region: med D_obs/radio(q=1) >= 2.0 y, con q=1/8, med D_obs/radio >= '
                                             '1.25 y med D_ctrl/radio >= 1.25',
                                    'medidas': {'med_Dobs_radio_q1': 16.3005,
                                                'med_Dctrl_radio_q1': 16.1111,
                                                'med_Dobs_radio_q1_8': 2.6783,
                                                'med_Dctrl_radio_q1_8': 2.5813},
                                    'radio': {'q1': 4.8232, 'q1_8': 4.4397},
                                    'semillas': 'SeedSequence(20261114, spawn_key=(333, 9, 11, real, sub[, qi])), '
                                                'real = 0..9; sub 0 puntos, 1 copias del observado, 2 dilución, 3 '
                                                'intervalos del observado, 4 subconjunto de base del control, 5 '
                                                'puntos del control (región), 6 copias del control, 7 intervalos del '
                                                'control; qi = 0 (q=1), 3 (q=1/8). Banda: SeedSequence(20261114, '
                                                'spawn_key=(333, 8, qi, s, 0|1)), s = 0..99 ajuste, 100..199 '
                                                'calibración',
                                    'N': {'N0': 8000, 'q1': 8000, 'q1_8_mediana': 989.0, 'q1_8_rango': [949, 1025]}},
 'S2a_cono_circular': {'identificable': True,
                       'familia': 'S2_conos_suaves',
                       'regla': 'region: med D_obs/radio(q=1) >= 2.0 y, con q=1/8, med D_obs/radio >= 1.25 y med '
                                'D_ctrl/radio >= 1.25',
                       'medidas': {'med_Dobs_radio_q1': 13.0784,
                                   'med_Dctrl_radio_q1': 13.1095,
                                   'med_Dobs_radio_q1_8': 2.373,
                                   'med_Dctrl_radio_q1_8': 2.259},
                       'radio': {'q1': 4.8232, 'q1_8': 4.4397},
                       'semillas': 'SeedSequence(20261114, spawn_key=(333, 9, 12, real, sub[, qi])), real = 0..9; '
                                   'sub 0 puntos, 1 copias del observado, 2 dilución, 3 intervalos del observado, 4 '
                                   'subconjunto de base del control, 5 puntos del control (región), 6 copias del '
                                   'control, 7 intervalos del control; qi = 0 (q=1), 3 (q=1/8). Banda: '
                                   'SeedSequence(20261114, spawn_key=(333, 8, qi, s, 0|1)), s = 0..99 ajuste, '
                                   '100..199 calibración',
                       'N': {'N0': 8000, 'q1': 8000, 'q1_8_mediana': 1006.5, 'q1_8_rango': [968, 1036]}},
 'S2b_cono_prisma_reuleaux': {'identificable': True,
                              'familia': 'S2_conos_suaves',
                              'regla': 'region: med D_obs/radio(q=1) >= 2.0 y, con q=1/8, med D_obs/radio >= 1.25 y '
                                       'med D_ctrl/radio >= 1.25',
                              'medidas': {'med_Dobs_radio_q1': 9.734,
                                          'med_Dctrl_radio_q1': 9.8693,
                                          'med_Dobs_radio_q1_8': 1.3793,
                                          'med_Dctrl_radio_q1_8': 1.7277},
                              'radio': {'q1': 4.8232, 'q1_8': 4.4397},
                              'semillas': 'SeedSequence(20261114, spawn_key=(333, 9, 13, real, sub[, qi])), real = '
                                          '0..9; sub 0 puntos, 1 copias del observado, 2 dilución, 3 intervalos del '
                                          'observado, 4 subconjunto de base del control, 5 puntos del control '
                                          '(región), 6 copias del control, 7 intervalos del control; qi = 0 (q=1), 3 '
                                          '(q=1/8). Banda: SeedSequence(20261114, spawn_key=(333, 8, qi, s, 0|1)), s '
                                          '= 0..99 ajuste, 100..199 calibración',
                              'N': {'N0': 8000, 'q1': 8000, 'q1_8_mediana': 1022.0, 'q1_8_rango': [936, 1041]}},
 'S2c_cono_prisma_semicircular': {'identificable': True,
                                  'familia': 'S2_conos_suaves',
                                  'regla': 'region: med D_obs/radio(q=1) >= 2.0 y, con q=1/8, med D_obs/radio >= '
                                           '1.25 y med D_ctrl/radio >= 1.25',
                                  'medidas': {'med_Dobs_radio_q1': 10.7931,
                                              'med_Dctrl_radio_q1': 10.7395,
                                              'med_Dobs_radio_q1_8': 1.8821,
                                              'med_Dctrl_radio_q1_8': 1.7662},
                                  'radio': {'q1': 4.8232, 'q1_8': 4.4397},
                                  'semillas': 'SeedSequence(20261114, spawn_key=(333, 9, 14, real, sub[, qi])), real '
                                              '= 0..9; sub 0 puntos, 1 copias del observado, 2 dilución, 3 '
                                              'intervalos del observado, 4 subconjunto de base del control, 5 puntos '
                                              'del control (región), 6 copias del control, 7 intervalos del control; '
                                              'qi = 0 (q=1), 3 (q=1/8). Banda: SeedSequence(20261114, '
                                              'spawn_key=(333, 8, qi, s, 0|1)), s = 0..99 ajuste, 100..199 '
                                              'calibración',
                                  'N': {'N0': 8000, 'q1': 8000, 'q1_8_mediana': 1005.5, 'q1_8_rango': [971, 1043]}},
 'S2d_cono_prisma_lenticular': {'identificable': True,
                                'familia': 'S2_conos_suaves',
                                'regla': 'region: med D_obs/radio(q=1) >= 2.0 y, con q=1/8, med D_obs/radio >= 1.25 '
                                         'y med D_ctrl/radio >= 1.25',
                                'medidas': {'med_Dobs_radio_q1': 8.9065,
                                            'med_Dctrl_radio_q1': 9.0431,
                                            'med_Dobs_radio_q1_8': 1.3093,
                                            'med_Dctrl_radio_q1_8': 1.383},
                                'radio': {'q1': 4.8232, 'q1_8': 4.4397},
                                'semillas': 'SeedSequence(20261114, spawn_key=(333, 9, 15, real, sub[, qi])), real = '
                                            '0..9; sub 0 puntos, 1 copias del observado, 2 dilución, 3 intervalos '
                                            'del observado, 4 subconjunto de base del control, 5 puntos del control '
                                            '(región), 6 copias del control, 7 intervalos del control; qi = 0 (q=1), '
                                            '3 (q=1/8). Banda: SeedSequence(20261114, spawn_key=(333, 8, qi, s, '
                                            '0|1)), s = 0..99 ajuste, 100..199 calibración',
                                'N': {'N0': 8000, 'q1': 8000, 'q1_8_mediana': 976.5, 'q1_8_rango': [947, 1020]}}}
# --- IDENTIFICABLES (END) ---
