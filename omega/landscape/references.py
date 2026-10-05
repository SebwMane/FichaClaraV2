"""Referencias concretas para el paisaje extremal L-2/L-3b (Rev. 2, R2.5); puro, sin E/S.

Todas las funciones devuelven matrices W simetricas (n,n), float64, diagonal 0 y valores en [0,1],
validas para `validate_weight_matrix`. No se importa `omega.experiments` (prohibido fuera de
experiments/): la reticula periodica se reimplementa aqui y los tests la comparan con
`reference_graphs.periodic_lattice`.

Orden colex de aristas (i<j): ordenadas por j y luego por i: (0,1),(0,2),(1,2),(0,3),(1,3),(2,3),...
Llenar en orden colex produce el clique sobre los primeros nodos (la configuracion con mas
triangulos para una masa dada, cota de Kruskal-Katona/Lovasz).
"""

from __future__ import annotations

import math
from typing import Final, Literal

import numpy as np

from omega.controls.erdos_renyi import erdos_renyi_gnm
from omega.controls.random_geometric import rgg_torus
from omega.network.weights import from_upper_triangle, validate_weight_matrix
from omega.types import FloatArray, IntArray

__all__ = [
    "Profile",
    "PROFILES",
    "torus_lattice_3d",
    "rgg3_radius",
    "rgg3_torus_binary",
    "rgg3_torus_points",
    "torus_distances",
    "weighted_rgg_from_points",
    "fit_amplitude_for_mass",
    "decorated_lattice_3d",
    "colex_pairs",
    "colex_clique_with_mass",
    "clique_union",
    "connected_caveman",
    "erdos_renyi_m",
    "uniform_with_mass",
    "match_mass",
]

Profile = Literal["step", "linear", "smooth"]
PROFILES: Final[tuple[Profile, ...]] = ("step", "linear", "smooth")
_SNAP: Final = 1e-9


def _check_int(name: str, v: int, minimum: int) -> None:
    if isinstance(v, bool) or not isinstance(v, int) or v < minimum:
        raise ValueError(f"{name} debe ser entero >= {minimum}, recibido {v!r}")


def _check_rng(rng: np.random.Generator) -> None:
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng debe ser numpy.random.Generator")


def _n_pairs(n: int) -> int:
    return n * (n - 1) // 2


def colex_pairs(n: int) -> tuple[IntArray, IntArray]:
    """Pares (i<j) en orden colex (por j, luego i): devuelve (I, J) de longitud C(n,2)."""
    _check_int("n", n, 2)
    i, j = np.triu_indices(n, k=1)
    order = np.lexsort((i, j))
    return np.asarray(i[order], dtype=np.int64), np.asarray(j[order], dtype=np.int64)


def _from_pairs(n: int, i: IntArray, j: IntArray, v: FloatArray) -> FloatArray:
    w = np.zeros((n, n), dtype=np.float64)
    w[i, j] = v
    w[j, i] = v
    return w


def _snap_total(total: float, n: int) -> float:
    t = float(total)
    if not math.isfinite(t):
        raise ValueError("total debe ser finito")
    if t < -_SNAP or t > _n_pairs(n) + _SNAP:
        raise ValueError(f"total debe estar en [0, {_n_pairs(n)}], recibido {t}")
    if abs(t - round(t)) < _SNAP:
        t = float(round(t))
    return min(max(t, 0.0), float(_n_pairs(n)))


# ---------------------------------------------------------------- geometricas


def torus_lattice_3d(side: int) -> FloatArray:
    """Toro cubico side^3 con W=1 entre vecinos mas proximos (grado 6); side >= 3. Indice = x*s^2+y*s+z."""
    _check_int("side", side, 3)
    n = side**3
    idx = np.arange(n)
    coords = np.stack(np.unravel_index(idx, (side, side, side)), axis=1)
    w = np.zeros((n, n), dtype=np.float64)
    for axis in range(3):
        nxt = coords.copy()
        nxt[:, axis] = (nxt[:, axis] + 1) % side
        dst = np.ravel_multi_index(tuple(nxt.T), (side, side, side))
        w[idx, dst] = 1.0
        w[dst, idx] = 1.0
    return w


def rgg3_radius(n: int, mean_degree: float) -> float:
    """Radio r del RGG3 en el toro unidad con (n-1)(4 pi/3) r^3 = mean_degree (igual que rgg_torus)."""
    _check_int("n", n, 3)
    if not mean_degree > 0.0:
        raise ValueError("mean_degree debe ser > 0")
    return float((mean_degree / ((n - 1) * (4.0 * math.pi / 3.0))) ** (1.0 / 3.0))


def rgg3_torus_binary(n: int, k: float, rng: np.random.Generator) -> FloatArray:
    """RGG3 binario en el toro con grado medio esperado k (reutiliza `rgg_torus`)."""
    _check_rng(rng)
    return rgg_torus(n, 3, float(k), rng, weights="binary")


def rgg3_torus_points(n: int, rng: np.random.Generator) -> FloatArray:
    """n puntos U[0,1)^3. Mismo consumo de `rng` que rgg_torus: mismo rng/semilla => mismos puntos."""
    _check_int("n", n, 3)
    _check_rng(rng)
    return np.asarray(rng.random((n, 3)), dtype=np.float64)


def torus_distances(points: FloatArray) -> FloatArray:
    """Matriz (n,n) de distancias euclidianas periodicas en el toro unidad."""
    p = np.asarray(points, dtype=np.float64)
    if p.ndim != 2 or p.shape[0] < 2:
        raise ValueError("points debe tener forma (n,dim), n >= 2")
    diff = np.abs(p[:, None, :] - p[None, :, :])
    diff = np.minimum(diff, 1.0 - diff)
    return np.asarray(np.sqrt(np.sum(diff * diff, axis=2)), dtype=np.float64)


def _profile_values(x: FloatArray, profile: Profile) -> FloatArray:
    """phi(x) sobre x = d/r >= 0: escalon 1[x<1], lineal (1-x)+, suave (1-x^2)+."""
    if profile == "step":
        out = (x < 1.0).astype(np.float64)
    elif profile == "linear":
        out = np.where(x < 1.0, 1.0 - x, 0.0)
    elif profile == "smooth":
        out = np.where(x < 1.0, 1.0 - x * x, 0.0)
    else:
        raise ValueError(f"profile invalido: {profile!r}")
    return np.asarray(out, dtype=np.float64)


def _upper_phi(points: FloatArray, r: float, profile: Profile) -> FloatArray:
    if not r > 0.0:
        raise ValueError("r debe ser > 0")
    dist = torus_distances(points)
    iu = np.triu_indices(dist.shape[0], k=1)
    return _profile_values(dist[iu] / r, profile)


def weighted_rgg_from_points(points: FloatArray, r: float, profile: Profile, a: float) -> FloatArray:
    """w_ij = min(1, a phi(d_ij/r)) con distancias periodicas; a >= 0. Con step y a=1 es el RGG binario."""
    if not a >= 0.0 or not math.isfinite(a):
        raise ValueError("a debe ser finito y >= 0")
    phi = _upper_phi(points, r, profile)
    n = int(np.asarray(points).shape[0])
    return from_upper_triangle(np.minimum(1.0, a * phi), n)


def fit_amplitude_for_mass(points: FloatArray, r: float, profile: Profile, total: float) -> float:
    """Amplitud a tal que sum_{i<j} min(1, a phi) = total, por biseccion; ValueError si inalcanzable.

    La masa es continua y creciente en a; su supremo (a -> inf) es el numero de pares con d < r.
    """
    t = float(total)
    if not math.isfinite(t) or t < 0.0:
        raise ValueError("total debe ser finito y >= 0")
    phi = _upper_phi(points, r, profile)
    tol = 1e-9 * max(1.0, t)

    def mass(a: float) -> float:
        return float(np.minimum(1.0, a * phi).sum())

    if t <= tol:
        return 0.0
    hi = 1.0
    while mass(hi) < t and hi < 1e12:
        hi *= 2.0
    if mass(hi) < t - tol:
        raise ValueError(f"masa {t} inalcanzable con r={r}, perfil={profile} (maximo ~ {mass(hi):.6f})")
    lo = 0.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mass(mid) < t:
            lo = mid
        else:
            hi = mid
    return hi


def decorated_lattice_3d(side: int = 3, block: int = 8) -> FloatArray:
    """Reticula decorada side^3 (x) K_block: un K_block por sitio del toro; sitios vecinos se unen vertice a vertice
    (vertice t del bloque a con vertice t del bloque b). Nodo = sitio*block + t. Grado = (block-1) + 6 (=13 para 3,8).

    side >= 3 garantiza que los vecinos +-1 de cada sitio son 6 sitios distintos (en lado 3, +1 y -1 difieren).
    """
    _check_int("side", side, 3)
    _check_int("block", block, 2)
    sites = torus_lattice_3d(side)
    inner = np.kron(np.eye(side**3), np.ones((block, block)) - np.eye(block))
    links = np.kron(sites, np.eye(block))
    w = np.asarray(inner + links, dtype=np.float64)
    validate_weight_matrix(w)
    return w


# ---------------------------------------------------------------- no geometricas


def colex_clique_with_mass(n: int, total: float) -> FloatArray:
    """Llena aristas en orden colex con peso 1 hasta masa `total`; la ultima arista es parcial."""
    _check_int("n", n, 2)
    t = _snap_total(total, n)
    i, j = colex_pairs(n)
    v = np.zeros(i.shape[0], dtype=np.float64)
    m = int(math.floor(t))
    v[:m] = 1.0
    rest = t - m
    if rest > 0.0 and m < v.shape[0]:
        v[m] = rest
    return _from_pairs(n, i, j, v)


def _block_bounds(n: int, size: int) -> list[tuple[int, int]]:
    _check_int("n", n, 2)
    _check_int("size", size, 2)
    if size > n:
        raise ValueError("size no puede superar n")
    nb = n // size
    bounds = [(b * size, (b + 1) * size) for b in range(nb)]
    bounds[-1] = (bounds[-1][0], n)
    return bounds


def clique_union(n: int, size: int) -> FloatArray:
    """floor(n/size) cliques K_size disjuntas sobre nodos consecutivos; los nodos sobrantes se incorporan a la ultima."""
    w = np.zeros((n, n), dtype=np.float64)
    for a, b in _block_bounds(n, size):
        w[a:b, a:b] = 1.0
    np.fill_diagonal(w, 0.0)
    return w


def connected_caveman(n: int, size: int) -> FloatArray:
    """Caveman conectado (estandar de networkx `connected_caveman_graph`): cliques en anillo; en cada clique se quita
    la arista (a, a+1) y se anade (a, a-1 mod n) hacia la clique anterior. Sobrantes en la ultima clique. size >= 3."""
    _check_int("size", size, 3)
    w = clique_union(n, size)
    bounds = _block_bounds(n, size)
    if len(bounds) < 2:
        raise ValueError("se requieren >= 2 cliques")
    for a, _ in bounds:
        w[a, a + 1] = w[a + 1, a] = 0.0
        prev = (a - 1) % n
        w[a, prev] = w[prev, a] = 1.0
    return w


def erdos_renyi_m(n: int, m: int, rng: np.random.Generator) -> FloatArray:
    """G(n,m) binario (reutiliza `erdos_renyi_gnm`)."""
    return erdos_renyi_gnm(n, m, rng)


def uniform_with_mass(n: int, total: float) -> FloatArray:
    """W_ij = total / C(n,2) para todo i != j."""
    _check_int("n", n, 2)
    t = _snap_total(total, n)
    w = np.full((n, n), t / _n_pairs(n), dtype=np.float64)
    np.fill_diagonal(w, 0.0)
    return w


def match_mass(w: FloatArray, total: float, *, tol: float = 1e-9) -> FloatArray | None:
    """Iguala sum_{i<j} W a `total` modificando aristas en orden colex con a lo sumo una arista parcial.

    Si falta masa: sube a 1 las aristas en orden colex (las que ya valen 1 no cuentan) y la ultima, parcialmente.
    Si sobra: baja a 0 las aristas existentes en orden colex INVERSO (de la ultima a la primera), la ultima parcialmente.
    Devuelve None si no es igualable (total fuera de [0, C(n,2)]). No muta `w`.
    """
    validate_weight_matrix(w)
    n = w.shape[0]
    t = float(total)
    if not math.isfinite(t) or t < -tol or t > _n_pairs(n) + tol:
        return None
    t = min(max(t, 0.0), float(_n_pairs(n)))
    i, j = colex_pairs(n)
    v = np.array(w[i, j], dtype=np.float64)
    d = t - float(v.sum())
    if abs(d) <= tol:
        return np.array(w, dtype=np.float64)
    if d > 0.0:
        cum = np.cumsum(1.0 - v)
        if cum[-1] < d - tol:
            return None
        k = int(np.searchsorted(cum, d, side="left"))
        k = min(k, v.shape[0] - 1)
        before = float(cum[k - 1]) if k > 0 else 0.0
        v[:k] = 1.0
        v[k] = min(1.0, v[k] + (d - before))
    else:
        e = -d
        rv = v[::-1]
        cum = np.cumsum(rv)
        if cum[-1] < e - tol:
            return None
        k = int(np.searchsorted(cum, e, side="left"))
        k = min(k, rv.shape[0] - 1)
        before = float(cum[k - 1]) if k > 0 else 0.0
        rv = rv.copy()
        rv[:k] = 0.0
        rv[k] = max(0.0, rv[k] - (e - before))
        v = rv[::-1].copy()
    out = _from_pairs(n, i, j, v)
    validate_weight_matrix(out)
    return out
