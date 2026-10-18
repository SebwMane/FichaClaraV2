"""Dinamica SQ (docs/OMEGA_R1_1_PRERREGISTRO.md §1, congelado). Solo N es parametro; E, grado, ciclos y componentes son salidas.

Reglas (3 B x 2 R x 3 O): borrado por B-s / B-t / B-st, reenganche R2 / R3, ordenes SYNC / ASYNC / RAND.
La dinamica NO lee ninguna metrica de `record` (E/N, fracciones s,t, componentes, gigante): son solo salida (G1).
Unico RNG: el `numpy.random.Generator` recibido (de `rng_from_key`), consumido en bloques por `Rnd`.

Decisiones de implementacion documentadas:
 * SYNC: disparadores evaluados con matrices dispersas sobre el grafo actual; reenganche con candidatos calculados sobre la
   instantanea posterior al borrado; si dos reenganches coinciden en la misma arista, el segundo se omite (no hay multiaristas).
 * RAND: E y N se fijan al inicio del barrido; si en un evento de arista el grafo ya no tiene aristas, el evento no hace nada.
 * Un evento de vertice (RAND) o la fase de vertices (ASYNC) solo actua si grado <= 1 en el grafo vigente.
 * Vecino no existente: si el reenganche no encuentra candidato a distancia r, cae al vertice uniforme no vecino.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

VARIANT_B = ("B-s", "B-t", "B-st")
VARIANT_R = ("R2", "R3")
VARIANT_O = ("SYNC", "ASYNC", "RAND")


# ------------------------------------------------------------------ disparadores locales (conjuntos)


def s_edge(adj: list[set[int]], u: int, v: int) -> bool:
    """True sii la arista uv esta en algun 4-ciclo: existe a in N(u)\\{v}, b in N(v)\\{u}, a != b, a ~ b.

    Para a in N(u), adj[a] & N(v) contiene siempre a u; basta otro elemento b != u (b != a pues a no es su propio vecino)."""
    nu, nv = adj[u], adj[v]
    if len(nu) > len(nv):
        nu, nv = nv, nu
    for a in nu:
        if a != u and a != v and len(adj[a] & nv) > 1:
            return True
    return False


def t_edge(adj: list[set[int]], u: int, v: int) -> bool:
    return bool(adj[u] & adj[v])


def should_delete(adj: list[set[int]], u: int, v: int, b: str) -> bool:
    if b == "B-s":
        return not s_edge(adj, u, v)
    if b == "B-t":
        return not t_edge(adj, u, v)
    if b == "B-st":
        return not t_edge(adj, u, v) and not s_edge(adj, u, v)
    raise ValueError(b)


# ------------------------------------------------------------------ grafo mutable


class Graph:
    """Adyacencia como lista de conjuntos + lista de aristas (codigos u*N+v, u<v) con borrado O(1) (swap-pop)."""

    def __init__(self, n: int, adj: list[set[int]] | None = None):
        self.n = n
        self.adj: list[set[int]] = adj if adj is not None else [set() for _ in range(n)]
        self.el: list[int] = []
        self.pos: dict[int, int] = {}
        for u in range(n):
            for v in self.adj[u]:
                if u < v:
                    self.pos[u * n + v] = len(self.el)
                    self.el.append(u * n + v)

    @classmethod
    def from_csr(cls, a: sparse.spmatrix | sparse.sparray) -> Graph:
        m = sparse.csr_array(a)
        n = m.shape[0]
        adj = [set(m.indices[m.indptr[i]:m.indptr[i + 1]].tolist()) - {i} for i in range(n)]
        return cls(n, adj)

    def code(self, u: int, v: int) -> int:
        return u * self.n + v if u < v else v * self.n + u

    def has(self, u: int, v: int) -> bool:
        return v in self.adj[u]

    def add(self, u: int, v: int) -> bool:
        if u == v or v in self.adj[u]:
            return False
        self.adj[u].add(v)
        self.adj[v].add(u)
        c = self.code(u, v)
        self.pos[c] = len(self.el)
        self.el.append(c)
        return True

    def remove(self, u: int, v: int) -> None:
        self.adj[u].discard(v)
        self.adj[v].discard(u)
        c = self.code(u, v)
        i = self.pos.pop(c)
        last = self.el.pop()
        if i < len(self.el):
            self.el[i] = last
            self.pos[last] = i

    @property
    def n_edges(self) -> int:
        return len(self.el)

    def edge_array(self) -> np.ndarray:
        return np.asarray(self.el, dtype=np.int64)

    def edge_set(self) -> set[int]:
        return set(self.el)

    def to_csr(self) -> sparse.csr_array:
        e = self.edge_array()
        u, v = e // self.n, e % self.n
        m = sparse.coo_array((np.ones(2 * e.size), (np.concatenate([u, v]), np.concatenate([v, u]))), shape=(self.n, self.n))
        a = sparse.csr_array(m)
        a.data[:] = 1.0
        return a

    def degrees(self) -> np.ndarray:
        return np.fromiter((len(s) for s in self.adj), dtype=np.int64, count=self.n)


class Rnd:
    """Aleatoriedad por bloques sobre un unico Generator (el de rng_from_key)."""

    def __init__(self, rng: np.random.Generator, block: int = 1 << 16):
        self.rng, self.block = rng, block
        self.buf: list[float] = []
        self.i = 0

    def u(self) -> float:
        if self.i >= len(self.buf):
            self.buf = self.rng.random(self.block).tolist()
            self.i = 0
        x = self.buf[self.i]
        self.i += 1
        return x

    def below(self, k: int) -> int:
        return min(int(self.u() * k), k - 1)

    def permutation(self, k: int) -> np.ndarray:
        return self.rng.permutation(k)


# ------------------------------------------------------------------ reenganche


def exact_distance_set(adj: list[set[int]], v: int, r: int) -> set[int]:
    """Vertices a distancia exactamente r (r = 2 o 3) de v (BFS de radio r)."""
    seen = {v} | adj[v]
    layer = set(adj[v])
    for _ in range(r - 1):
        nxt: set[int] = set()
        for x in layer:
            nxt |= adj[x]
        nxt -= seen
        seen |= nxt
        layer = nxt
        if not layer:
            break
    return layer


def choose_target(adj: list[set[int]], v: int, rmode: str, rnd: Rnd) -> int:
    """Destino del reenganche de v (grado 0 o 1) sobre `adj`; -1 si no hay (solo en grafos minusculos)."""
    n = len(adj)
    if len(adj[v]) == 1:
        cand = exact_distance_set(adj, v, 2 if rmode == "R2" else 3)
        if cand:
            lst = sorted(cand)
            return lst[rnd.below(len(lst))]
    for _ in range(10_000):  # grado 0, o sin candidato: vertice uniforme distinto de v y no vecino (forzamiento no local)
        t = rnd.below(n)
        if t != v and t not in adj[v]:
            return t
    return -1


# ------------------------------------------------------------------ SYNC: disparadores con matrices dispersas


def edge_flags(a: sparse.csr_array) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """(u, v, s, t) por arista u<v. s = (A^3)_uv - deg u - deg v + 1 > 0 ; t = (A^2)_uv > 0. A binaria simetrica."""
    n = a.shape[0]
    m = sparse.csr_array(a, dtype=np.int32)
    m.data[:] = 1
    m.setdiag(0)
    m.eliminate_zeros()
    m.sort_indices()
    deg = np.diff(m.indptr).astype(np.int64)
    rows = np.repeat(np.arange(n, dtype=np.int64), deg)
    cols = m.indices.astype(np.int64)
    if cols.size == 0:
        z = np.empty(0, dtype=np.int64)
        return z, z, np.empty(0, dtype=bool), np.empty(0, dtype=bool)
    a2 = sparse.csr_array(m @ m)
    # t: (A^2)_uv > 0 en las aristas
    ta = sparse.csr_array(m.multiply(a2))
    ta.eliminate_zeros()
    ta.sort_indices()
    trow = np.repeat(np.arange(n, dtype=np.int64), np.diff(ta.indptr))
    tcode = trow * n + ta.indices.astype(np.int64)
    code = rows * n + cols
    t = np.zeros(code.size, dtype=bool)
    if tcode.size:
        idx = np.searchsorted(code, tcode)
        t[idx] = True
    # s: A3_uv por bloques de filas
    s = np.zeros(code.size, dtype=bool)
    avg = max(1.0, a2.nnz / n) * max(1.0, m.nnz / n)
    chunk = int(max(64, min(n, 4e6 / avg)))
    for lo in range(0, n, chunk):
        hi = min(n, lo + chunk)
        sub = sparse.csr_array(a2[lo:hi] @ m)
        sub = sparse.csr_array(sub.multiply(m[lo:hi]))
        sub.sort_indices()
        srow = np.repeat(np.arange(lo, hi, dtype=np.int64), np.diff(sub.indptr))
        scode = srow * n + sub.indices.astype(np.int64)
        a3 = sub.data.astype(np.int64)
        val = a3 - deg[srow] - deg[sub.indices] + 1 > 0
        idx = np.searchsorted(code, scode)
        s[idx] = val
    keep = rows < cols
    return rows[keep], cols[keep], s[keep], t[keep]


def delete_mask(s: np.ndarray, t: np.ndarray, b: str) -> np.ndarray:
    if b == "B-s":
        return ~s
    if b == "B-t":
        return ~t
    if b == "B-st":
        return ~s & ~t
    raise ValueError(b)


# ------------------------------------------------------------------ barridos


def sweep_sync(g: Graph, b: str, r: str, rnd: Rnd) -> tuple[int, int]:
    n = g.n
    u, v, s, t = edge_flags(g.to_csr()) if g.n_edges else (np.empty(0, int),) * 2 + (np.empty(0, bool),) * 2
    dm = delete_mask(s, t, b)
    for x, y in zip(u[dm].tolist(), v[dm].tolist()):
        g.remove(x, y)
    dels = int(dm.sum())
    low = [x for x in range(n) if len(g.adj[x]) <= 1]  # grados sobre el grafo posterior al borrado
    order = rnd.permutation(len(low))
    targets = [(low[i], choose_target(g.adj, low[i], r, rnd)) for i in order.tolist()]  # candidatos sobre la instantanea
    adds = 0
    for x, y in targets:
        if y >= 0 and g.add(x, y):  # duplicado (coincidencia de dos reenganches) -> se omite
            adds += 1
    return dels, adds


def sweep_async(g: Graph, b: str, r: str, rnd: Rnd) -> tuple[int, int]:
    dels = adds = 0
    el = g.edge_array()
    n = g.n
    for i in rnd.permutation(el.size).tolist():
        c = int(el[i])
        x, y = divmod(c, n)
        if should_delete(g.adj, x, y, b):
            g.remove(x, y)
            dels += 1
    for x in rnd.permutation(n).tolist():
        if len(g.adj[x]) <= 1:
            y = choose_target(g.adj, x, r, rnd)
            if y >= 0 and g.add(x, y):
                adds += 1
    return dels, adds


def sweep_rand(g: Graph, b: str, r: str, rnd: Rnd) -> tuple[int, int]:
    dels = adds = 0
    n = g.n
    e0 = g.n_edges
    p_edge = e0 / (e0 + n)
    for _ in range(e0 + n):
        if rnd.u() < p_edge:
            if g.n_edges:
                c = g.el[rnd.below(g.n_edges)]
                x, y = divmod(c, n)
                if should_delete(g.adj, x, y, b):
                    g.remove(x, y)
                    dels += 1
        else:
            x = rnd.below(n)
            if len(g.adj[x]) <= 1:
                y = choose_target(g.adj, x, r, rnd)
                if y >= 0 and g.add(x, y):
                    adds += 1
    return dels, adds


SWEEPS = {"SYNC": sweep_sync, "ASYNC": sweep_async, "RAND": sweep_rand}


def is_absorbed(g: Graph, b: str) -> bool:
    """Ningun vertice de grado <= 1 y ningun disparador de borrado activo (con salida temprana)."""
    if any(len(s) <= 1 for s in g.adj):
        return False
    n = g.n
    for c in g.el:
        x, y = divmod(c, n)
        if should_delete(g.adj, x, y, b):
            return False
    return True


# ------------------------------------------------------------------ registro (solo salida)


def record(g: Graph, sweep: int, dels: int, adds: int) -> dict[str, Any]:
    a = g.to_csr()
    n = g.n
    if g.n_edges:
        _, _, s, t = edge_flags(a)
        fs, ft = float(s.mean()), float(t.mean())
    else:
        fs = ft = 0.0
    nc, lab = connected_components(a, directed=False)
    giant = int(np.bincount(lab).max())
    deg = g.degrees()
    return {"sweep": sweep, "E_over_N": g.n_edges / n, "frac_s": fs, "frac_t": ft, "n_components": int(nc),
            "giant_frac": giant / n, "frac_deg_le1": float((deg <= 1).mean()), "deletions": dels, "reattachments": adds}


def run_dynamics(g: Graph, b: str, r: str, o: str, rng: np.random.Generator, T: int = 100,
                 record_every: int = 5) -> dict[str, Any]:
    """Hasta T barridos; absorcion comprobada tras cada barrido. Devuelve trayectoria (cada 5 barridos y final)."""
    rnd = Rnd(rng)
    step = SWEEPS[o]
    traj = [record(g, 0, 0, 0)]
    absorbed = is_absorbed(g, b)
    sweeps = 0
    last_recorded = 0
    while sweeps < T and not absorbed:
        dels, adds = step(g, b, r, rnd)
        sweeps += 1
        absorbed = is_absorbed(g, b)
        if sweeps % record_every == 0 or sweeps == T or absorbed:
            traj.append(record(g, sweeps, dels, adds))
            last_recorded = sweeps
    if last_recorded != sweeps:
        traj.append(record(g, sweeps, 0, 0))
    return {"sweeps": sweeps, "absorbed": bool(absorbed), "trajectory": traj}
