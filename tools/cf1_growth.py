"""OMEGA CF-1: crecimiento y coalescencia por completion flag causal (docs/OMEGA_CF1_PRERREGISTRO.md, congelado).

Salida: results/cf1/{runs.jsonl, summary.json} (humo: results/cf1_smoke). Uso:
  python tools/cf1_growth.py --smoke                 # N <= 2000, 1 semilla, todas las variantes
  python tools/cf1_growth.py [--procs 3] [--fresh]   # corrida completa (N1=1e4, N2=8e4, N3=6.4e5)
runs.jsonl: un registro por (variante, w, semilla), determinista (los tiempos solo van a summary.json). REANUDABLE por id
(--fresh lo borra). Todo el azar via omega.c0.references.rng_from_key((20261025, id_variante, w, semilla, k)).

Decisiones de implementacion (lectura mas literal de las ambiguedades)
----------------------------------------------------------------------
1. Estructura. Orden parcial con listas up/down de cubiertas (ids crecientes; los elementos nuevos solo se anaden ENCIMA, asi
   que las relaciones entre elementos existentes nunca cambian). Regla C_flag = la de tools/cf01_completion.py con k = flag:
   instancia (z, S), S subconjunto de up(z), |S| >= 2, con J(S\\{s}) existente para todo s y J(S) inexistente; J(T) (|T|=1 el
   propio elemento; si no, el elemento de MENOR ID en la interseccion de up(J(T\\{t})), t en T). El codigo de deteccion
   (nivel a nivel sobre A) se COPIA de cf01 (no se importa: cf01 mantiene cache global y recalcula `below`, que no escala).
   El elemento nuevo cubre {J(S\\{s})}. Instancias degeneradas (los sub-supremos no son todos distintos) se descartan y se
   cuentan (`degenerate`); son permanentes (los J positivos son estables), asi que no se mantienen como pendientes.
2. Mantenimiento incremental. Al crear u con cubiertas L: se refrescan (se recalculan sus instancias) los elementos de L y todo z
   por debajo de u a distancia de cobertura d <= min(|up(z)|, maxm) con |up(z)| >= 2 y (sin CAP) con >= 2 cubiertas superiores
   estrictamente por debajo de u (criterio de cf01), donde maxm = mayor |S| visto en una instancia o conjunto A. Es suficiente:
   una instancia en z solo cambia si algun J(T), T subconjunto de up(z) con |T| = m, pasa a existir (entonces u = J(T), todo t
   de T esta bajo u y existe un camino de cobertura de longitud m de z a u) o, con CAP, si cambia |up| de algun l de L. Los J
   positivos se cachean globalmente (jpos; son estables). El test `test_incremental_matches_full_recompute` audita el
   mantenimiento (instancias, bloqueos por ancestros y candidatos ER) contra el recalculo completo.
3. CAUSAL-LOCAL. Z = conjunto de elementos con alguna instancia pendiente. z es elegible sii ningun z' de Z es antecesor
   ESTRICTO de z. Se mantiene blockers[z] = {z' en Z : z' < z}. Cuando z entra en Z se prueba el orden con cada w de Z
   (is_anc(y, w): DFS por `down` desde w podando rangos < rango(y); un antecesor estricto tiene rango estrictamente menor, asi que
   la poda es exacta). Como las relaciones entre elementos existentes no cambian, blockers solo se actualiza al entrar/salir de
   Z. NO se usa ningun atajo "por rango": la elegibilidad es exactamente la definicion de cono pasado. Se elige uniformemente entre
   las INSTANCIAS elegibles (z, S) con el RNG k = 0 (lista ordenada por (z, S)).
4. CAUSAL-RANK. rango(u) = 1 + max rango de sus cubiertas inferiores, rango 0 en los elementos minimales (para semillas con varios
   minimales: R20 y D2+2 antes de la coalescencia = cadena mas larga desde el minimal que corresponda, es decir 0 en
   TODOS los minimales). Clave de seleccion: (rango z, -aridad |S|, z, S) minimo (aridad grande primero, luego id).
5. ER. Si no hay instancias pendientes: candidatos = {p : |up(p)| < |down(p)|} (la raiz y todo minimal tienen down = 0, nunca
   son candidatos); se elige uniforme (RNG k = 1) y se anade una cubierta superior nueva (down = {p}). Si no hay candidatos se
   registra status = ER_SIN_CANDIDATOS y la corrida termina en ese tamano. EL: p uniforme sobre TODOS los elementos (k = 1).
   ASYNC: en cada paso se saca un numero del RNG k = 3 (siempre, para que el flujo sea determinista); con prob. 1/2 se
   extiende (ER) si hay candidatos aunque haya compleciones pendientes; si no, se completa; si no hay pendientes se extiende.
6. CAP(c). Una operacion se prohibe si dejaria |up(l)| > c en algun l que reciba una cubierta superior nueva (en
   compleciones: los l de L; en extensiones: p). Las instancias que violarian la cota se retiran de las pendientes
   (`cap_blocked`; el bloqueo es permanente porque |up| solo crece) y los candidatos ER exigen |up| < c. El maximo de la semilla
   (raiz con up = 6 > c) no es una violacion: ninguna operacion le anade cubiertas. PAIR: solo m = 2 (C2).
7. Semillas. S_w: raiz (id 0) + w atomos. R20: 20 puntos uniformes en [0,1]^2 (RNG k = 2), orden producto estricto
   (p < q sii x_p<x_q e y_p<y_q), ids por x creciente, cubiertas = diagrama de Hasse. D2+2 / S2uS3 (V2, w = 22 / 23 en
   la clave): cada componente se crece INDEPENDIENTEMENTE con las reglas de V1 hasta N1/2 elementos CADA UNA (lectura literal del
   prerregistro "cada componente hasta N1/2"; el total es N1), componente A y luego B con los mismos flujos RNG; despues UN
   elemento nuevo cubre a un maximal elegido uniformemente (RNG k = 1) de cada componente, y se sigue con V1 (estructura
   unica). La instantanea N1 se toma justo tras la coalescencia (n = N1 + 1).
8. Claves RNG: id_variante = V1:1, V1-R:2, V1-A:3, V2:4, V3:5, V4:60+c (c = 2,3,4), V5:7, V7:8, V0:0. w = ancho de semilla
   (V2: 22 y 23; V7: 20). semilla = 0,1,2. Medida (rc3.measure): clave (MASTER, id_variante, w, semilla, indice_de_tamano 0/1/2);
   rc3.measure le anade sufijos (1,3,5) internamente.
9. Instantaneas: primer n >= N_i (n crece de 1 en 1; en V2 N1 se toma en n = N1 + 1). Grafo medido = Hasse no dirigido del propio
   J. delta12 / delta23 = rc3.scaling_delta(R, R', n_giant, n_giant'); exclusiones X1-X4 por par = rc3.pair_exclusions;
   W5 = |E/N' - E/N| / (E/N) < 0.25 con E/N = nnz/(2n) del grafo completo; W6 = w6_calib.family_metrics + w6_rule con
   (tau_d, tau_k) congelados leidos de results/w6/summary.json (calibration.thresholds): W6.2 |d23 - d12| <= tau_d, W6.3
   |k3 - k2|/k2 <= tau_k (k = grado medio de la componente gigante = w6_calib.giant_degree), W6.4 = W5 en ambos pares. Una corrida
   que no llega a N3 (ER sin candidatos / coste) es 'NO_EVALUABLE a N3': w6_valid = None y cuenta como NO W6-valida en los
   criterios (W6 exige los tres tamanos).
10. Observables solo reportados. r_loc = mediana de |up| sobre elementos con |up| = |down|; r_loc*delta (r_loc(N2)*d12 y
    r_loc(N3)*d23); elementos por nivel de rango (lista completa si hay <= 3000 niveles; si no solo resumen); sigma_k =
    desviacion tipica del grado de Hasse (|up|+|down|); fraccion de esquinas ambiguas = #{p: |down|-|up| >= 2} / #{p: |up| <
    |down|}; supremos ambiguos = # de T distintos (|T| >= 2, subconjunto de up(z) para algun z, con J(T) existente segun el
    algoritmo nivel a nivel) con mas de un candidato en la interseccion; r_loc del futuro conjunto: V2 = mediana de |up| entre los
    elementos con |up|=|down| estrictamente por encima del elemento de coalescencia; V7 = idem para el futuro comun (elementos
    por encima de TODOS los minimales de R20); V1 = cono futuro de la raiz (todo). Dimension espectral: NO se calcula (no
    calibrada; se omitio por coste), campo null.
11. Criterios §5 (resumen.json): una semilla 'califica' si es W6-valida y no esta excluida por RC-3 en ninguno de los dos pares;
    A-positivo = al menos 2 semillas califican (se informa tambien la lectura estricta: >= 2 W6-validas y ninguna semilla
    excluida). Clase W6 = (n_w6_validas >= 2); clase RC-3 = (n_semillas_no_excluidas >= 2); 'robusto' = misma clase W6 y RC-3
    en V1, V1-R, V1-A y |d12(V1) - d12(V1-x)| <= 0.04 (medias sobre semillas con delta definido). Trazabilidad 'inicio': |d12(V1,w)
    - 1/w| <= 0.06 para w = 2,3,4. K1/K3: >= 2 de 3 semillas NO W6-validas (no evaluables incluidas); K2: |d12(V4,c) - 1/c| <= 0.08
    para c = 2,3,4. D1: d_w = 1/d12(V1,w); 'difiere de w0' = no 'inicio'; 'converge a un valor comun' = max d_w - min d_w <= 0.5 (umbral de
    esta implementacion, el prerregistro no da numero; se informa d_w para que el cerebro lo juzgue).
12. Guardas de coste: --max-seconds por corrida (defecto sin limite en humo, 6 h completo) y tope de pendientes; si se
    exceden, status = ABORTADO_COSTE y se registra el tamano alcanzado ('no evaluable a N3'); el diseno no cambia.
13. V0: C_flag causal (CAUSAL de cf01 reutilizando cf01.Poset + boolean_check, y los planificadores LOCAL y RANK de este
    modulo) desde S7 y S8; relabelado: se permutan los ids de S_w (w = 3..6) y se comprueba isomorfia (networkx, grafo
    dirigido de cubiertas) con el resultado sin permutar y con B_w.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import itertools  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import multiprocessing as mp  # noqa: E402
import statistics  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any, Callable  # noqa: E402

import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from p1d3_panel import _from_edges  # noqa: E402
from rc3 import measure  # noqa: E402
from w6_calib import family_metrics, giant_degree, w6_rule  # noqa: E402

import cf01_completion as cf01  # noqa: E402
from omega.c0.references import rng_from_key  # noqa: E402

MASTER_CF1 = 20261025
SIZES = (10_000, 80_000, 640_000)
SMOKE_SIZES = (200, 800, 2_000)
OUT_FULL = ROOT / "results" / "cf1"
OUT_SMOKE = ROOT / "results" / "cf1_smoke"
W6_SUMMARY = ROOT / "results" / "w6" / "summary.json"
VID = {"V0": 0, "V1": 1, "V1R": 2, "V1A": 3, "V2": 4, "V3": 5, "V4": 60, "V5": 7, "V7": 8}
MAX_PENDING = 50_000
LEVELS_FULL_MAX = 3000


# ====================================================================== estructura


class Grower:
    """Orden parcial creciente con C_flag (o PAIR), planificador, ER/EL, ASYNC y CAP."""

    def __init__(self, up: list[list[int]], down: list[list[int]], rngs: dict[int, np.random.Generator], *,
                 rule: str = "flag", sched: str = "LOCAL", ext: str = "ER", asynchronous: bool = False,
                 cap: int | None = None):
        assert rule in ("flag", "pair") and sched in ("LOCAL", "RANK") and ext in ("ER", "EL")
        self.rule, self.sched, self.ext, self.asynchronous, self.cap = rule, sched, ext, asynchronous, cap
        self.rngs = rngs
        self.up = [list(x) for x in up]
        self.down = [list(x) for x in down]
        n = len(self.up)
        self.rank = self._ranks()
        self.maxup = max((len(u) for u in self.up), default=0)
        self.n_edges = sum(len(u) for u in self.up)
        self.traj: list[dict] = []
        self.traj_next = 500
        self.maxm = 2  # mayor |S| visto en una instancia o en un conjunto A (cota de la profundidad de refresco)
        self.jpos: dict[tuple[int, ...], int] = {}
        self.pend: dict[int, list[tuple[int, ...]]] = {}
        self.blockers: dict[int, set[int]] = {}
        self.cand: list[int] = []
        self.cpos: dict[int, int] = {}
        self.degenerate: set[tuple] = set()
        self.cap_blocked: set[tuple] = set()
        self.n_complete = 0
        self.n_extend = 0
        for p in range(n):
            self._cand_update(p)
        for z in range(n):
            self._refresh(z)

    # ---- basicos
    @property
    def n(self) -> int:
        return len(self.up)

    def _ranks(self) -> list[int]:
        n = len(self.up)
        rank = [-1] * n
        indeg = [len(d) for d in self.down]
        stack = [i for i in range(n) if indeg[i] == 0]
        for i in stack:
            rank[i] = 0
        order = []
        while stack:
            v = stack.pop()
            order.append(v)
            for u in self.up[v]:
                indeg[u] -= 1
                if indeg[u] == 0:
                    stack.append(u)
        assert len(order) == n, "ciclo en la semilla"
        for v in order:
            if self.down[v]:
                rank[v] = 1 + max(rank[d] for d in self.down[v])
        return rank

    def _cand_update(self, p: int) -> None:
        if self.ext == "EL":
            return
        ok = len(self.up[p]) < len(self.down[p]) and (self.cap is None or len(self.up[p]) < self.cap)
        pos = self.cpos.get(p)
        if ok and pos is None:
            self.cpos[p] = len(self.cand)
            self.cand.append(p)
        elif not ok and pos is not None:
            last = self.cand.pop()
            if last != p:
                self.cand[pos] = last
                self.cpos[last] = pos
            del self.cpos[p]

    # ---- sub-supremos e instancias (copia de la logica de cf01_completion, sin caches globales)
    def _J(self, T: tuple[int, ...], memo: dict) -> int | None:
        if len(T) == 1:
            return T[0]
        r = self.jpos.get(T)  # los J positivos son estables (los ids nuevos son mayores): cache global
        if r is not None:
            return r
        if T in memo:
            return memo[T]
        if len(T) == 2:
            uy = self.up[T[1]]
            r = min((a for a in self.up[T[0]] if a in uy), default=None)
        else:
            subs = []
            for i in range(len(T)):
                j = self._J(T[:i] + T[i + 1:], memo)
                if j is None:
                    memo[T] = None
                    return None
                subs.append(j)
            c = set(self.up[subs[0]])
            for j in subs[1:]:
                c &= set(self.up[j])
                if not c:
                    break
            r = min(c) if c else None
        if r is None:
            memo[T] = None
        else:
            self.jpos[T] = r
        return r

    def _lower(self, S: tuple[int, ...], memo: dict) -> list[int] | None:
        L = {self._J(S[:i] + S[i + 1:], memo) for i in range(len(S))}
        if len(L) < len(S):
            return None
        return sorted(L)  # type: ignore[arg-type]

    def compute_z(self, z: int) -> list[tuple[int, ...]]:
        ups = sorted(self.up[z])
        if len(ups) < 2:
            return []
        memo: dict = {}
        kmax = len(ups) if self.rule == "flag" else 2
        res: list[tuple[int, ...]] = []
        A: set[tuple[int, ...]] = set()
        for x, y in itertools.combinations(ups, 2):
            if self._J((x, y), memo) is not None:
                A.add((x, y))
            else:
                res.append((x, y))
        m = 3
        last_a = 2 if A else 1
        while m <= kmax and A:
            A2: set[tuple[int, ...]] = set()
            for S in sorted(A):
                for e in ups:
                    if e <= S[-1]:
                        continue
                    S2 = S + (e,)
                    if all(S2[:i] + S2[i + 1:] in A for i in range(m - 1)):
                        if self._J(S2, memo) is not None:
                            A2.add(S2)
                        else:
                            res.append(S2)
            A = A2
            if A:
                last_a = m
            m += 1
        # cota de profundidad de refresco: un T de tamano t puede pasar a tener J solo si T era una instancia (res) o,
        # trivialmente, si sus (t-1)-subconjuntos tienen J (A nivel t-1)
        self.maxm = max(self.maxm, last_a + 1, max((len(S) for S in res), default=0))
        out = []
        for S in sorted(res):
            L = self._lower(S, memo)
            if L is None:
                self.degenerate.add((z, S))
                continue
            if self.cap is not None and any(len(self.up[l]) + 1 > self.cap for l in L):
                self.cap_blocked.add((z, S))
                continue
            out.append(S)
        return out

    # ---- conjunto de pendientes y bloqueos (cono pasado)
    def is_anc(self, y: int, w: int) -> bool:
        """y es antecesor estricto de w (poda exacta por rango)."""
        ry = self.rank[y]
        if self.rank[w] <= ry:
            return False
        seen = set()
        stack = [w]
        while stack:
            v = stack.pop()
            for d in self.down[v]:
                if d == y:
                    return True
                if d not in seen and self.rank[d] > ry:
                    seen.add(d)
                    stack.append(d)
        return False

    def _refresh(self, z: int) -> None:
        new = self.compute_z(z)
        had = z in self.pend
        if new:
            self.pend[z] = new
            if not had:
                bl: set[int] = set()
                for w in self.pend:
                    if w == z:
                        continue
                    if self.rank[w] < self.rank[z]:
                        if self.is_anc(w, z):
                            bl.add(w)
                    elif self.rank[w] > self.rank[z] and self.is_anc(z, w):
                        self.blockers[w].add(z)
                self.blockers[z] = bl
        elif had:
            del self.pend[z]
            del self.blockers[z]
            for w in self.blockers.values():
                w.discard(z)

    # ---- creacion de elementos
    def _create(self, L: list[int]) -> int:
        u = len(self.up)
        self.up.append([])
        self.down.append(list(L))
        self.rank.append(1 + max(self.rank[l] for l in L))
        self.n_edges += len(L)
        for l in L:
            self.up[l].append(u)
            if len(self.up[l]) > self.maxup:
                self.maxup = len(self.up[l])
        for l in L:
            self._cand_update(l)
        self._cand_update(u)
        # elementos a refrescar
        dirty = set(L)
        cand_dirty: list[tuple[int, int]] = []
        dist = {u: 0}
        frontier = [u]
        for d in range(1, min(self.maxup, self.maxm) + 1):
            nxt = []
            for v in frontier:
                for x in self.down[v]:
                    if x not in dist:
                        dist[x] = d
                        nxt.append(x)
                        if len(self.up[x]) >= 2 and d <= len(self.up[x]):
                            cand_dirty.append((x, d))
            frontier = nxt
            if not frontier:
                break
        for x, d in cand_dirty:
            # sin CAP basta que >= 2 cubiertas superiores de x esten estrictamente por debajo de u (criterio de cf01: un J(T)
            # nuevo = u exige que todos los t de T esten bajo u); con CAP un cambio de |up| de un l de L tambien cuenta
            if self.cap is not None or sum(1 for y in self.up[x] if y in dist) >= 2:
                dirty.add(x)
        for z in sorted(dirty):
            self._refresh(z)
        return u

    # ---- operaciones
    def eligible_instances(self) -> list[tuple[int, tuple[int, ...]]]:
        if self.sched == "LOCAL":
            return [(z, S) for z in sorted(self.pend) if not self.blockers[z] for S in self.pend[z]]
        return [(z, S) for z in sorted(self.pend) for S in self.pend[z]]

    def complete(self) -> None:
        if self.sched == "LOCAL":
            el = self.eligible_instances()
            z, S = el[int(self.rngs[0].integers(len(el)))]
        else:
            z, S = min(((self.rank[z], -len(S), z, S) for z in self.pend for S in self.pend[z]))[2:]
        L = self._lower(S, {})
        assert L is not None
        self._create(L)
        self.n_complete += 1

    def extend(self) -> None:
        if self.ext == "EL":
            p = int(self.rngs[1].integers(self.n))
        else:
            p = self.cand[int(self.rngs[1].integers(len(self.cand)))]
        self._create([p])
        self.n_extend += 1

    def step(self) -> str | None:
        draw_ext = self.asynchronous and self.rngs[3].random() < 0.5
        can_ext = self.ext == "EL" or bool(self.cand)
        if self.pend:
            if draw_ext and can_ext:
                self.extend()
            else:
                self.complete()
        else:
            if not can_ext:
                return "ER_SIN_CANDIDATOS"
            self.extend()
        return None

    def _traj_point(self) -> None:
        n = self.n
        ups = np.fromiter((len(u) for u in self.up), dtype=np.int64, count=n)
        dns = np.fromiter((len(d) for d in self.down), dtype=np.int64, count=n)
        inter = ups == dns
        self.traj.append({"n": n, "mean_hasse_degree": 2.0 * self.n_edges / n, "max_up": int(self.maxup),
                          "r_loc": float(np.median(ups[inter])) if inter.any() else None})

    def run(self, sizes: tuple[int, ...] | list[int], on_snap: Callable[[int], None] | None, max_n: int,
            deadline: float | None = None, si: int = 0, track: bool = False) -> tuple[str, int]:
        """Crece hasta max_n elementos. Devuelve (status, indice_de_la_siguiente_instantanea). El tiempo de las
        instantaneas (on_snap) NO cuenta para el limite de crecimiento. track: trayectoria (solo reporte, §8.2)."""
        ops = 0
        while True:
            while si < len(sizes) and self.n >= sizes[si]:
                if on_snap is not None:
                    t0 = time.perf_counter()
                    on_snap(si)
                    if deadline is not None:
                        deadline += time.perf_counter() - t0
                si += 1
            if self.n >= max_n:
                return "COMPLETO", si
            st = self.step()
            if st:
                return st, si
            ops += 1
            if track and self.n >= self.traj_next:
                self._traj_point()
                while self.traj_next <= self.n:
                    self.traj_next *= 2
            if len(self.pend) > MAX_PENDING:
                return "ABORTADO_COSTE", si
            if deadline is not None and ops % 64 == 0 and time.perf_counter() > deadline:
                return "ABORTADO_COSTE", si

    def closure(self, max_n: int = 100_000) -> bool:
        """Solo compleciones hasta que no queden instancias. True si termina."""
        while self.pend:
            if self.n >= max_n:
                return False
            self.complete()
        return True

    # ---- auditoria
    def full_pending(self) -> dict[int, list[tuple[int, ...]]]:
        out = {}
        for z in range(self.n):
            r = self.compute_z(z)
            if r:
                out[z] = r
        return out

    def below(self, u: int) -> set[int]:
        seen, st = set(), [u]
        while st:
            v = st.pop()
            for d in self.down[v]:
                if d not in seen:
                    seen.add(d)
                    st.append(d)
        return seen

    def audit(self) -> None:
        full = self.full_pending()
        assert full == self.pend, "pendientes incrementales != recalculo completo"
        for z in self.pend:
            true_bl = {y for y in self.pend if y != z and y in self.below(z)}
            assert true_bl == self.blockers[z], f"bloqueos de {z}: {true_bl} != {self.blockers[z]}"
        assert set(self.blockers) == set(self.pend)
        cand = {p for p in range(self.n) if self.ext == "ER" and len(self.up[p]) < len(self.down[p])
                and (self.cap is None or len(self.up[p]) < self.cap)}
        assert cand == set(self.cand) and len(self.cand) == len(self.cpos)

    def ambiguous_subjoins(self) -> int:
        """# de T (|T|>=2) con J(T) existente y mas de un candidato (nivel a nivel, como compute_z)."""
        counted: set[tuple[int, ...]] = set()
        amb = 0
        memo: dict = {}
        for z in range(self.n):
            ups = sorted(self.up[z])
            if len(ups) < 2 or len(ups) > 10:
                continue
            A: set[tuple[int, ...]] = set()
            for T in itertools.combinations(ups, 2):
                c = set(self.up[T[0]]) & set(self.up[T[1]])
                if c:
                    A.add(T)
                    if T not in counted:
                        counted.add(T)
                        amb += len(c) > 1
            m = 3
            while A and m <= len(ups) and self.rule == "flag":
                A2: set[tuple[int, ...]] = set()
                for S in sorted(A):
                    for e in ups:
                        if e <= S[-1]:
                            continue
                        S2 = S + (e,)
                        if all(S2[:i] + S2[i + 1:] in A for i in range(m - 1)):
                            subs = [self._J(S2[:i] + S2[i + 1:], memo) for i in range(m)]
                            c = set(self.up[subs[0]])
                            for j in subs[1:]:
                                c &= set(self.up[j])
                            if c:
                                A2.add(S2)
                                if S2 not in counted:
                                    counted.add(S2)
                                    amb += len(c) > 1
                A = A2
                m += 1
        return amb


# ====================================================================== semillas


def seed_S(w: int) -> tuple[list[list[int]], list[list[int]]]:
    up = [list(range(1, w + 1))] + [[] for _ in range(w)]
    down = [[]] + [[0] for _ in range(w)]
    return up, down


def seed_R20(rng: np.random.Generator, n: int = 20) -> tuple[list[list[int]], list[list[int]]]:
    pts = rng.random((n, 2))
    pts = pts[np.argsort(pts[:, 0], kind="stable")]
    less = [[bool(pts[i, 0] < pts[j, 0] and pts[i, 1] < pts[j, 1]) for j in range(n)] for i in range(n)]
    up: list[list[int]] = [[] for _ in range(n)]
    down: list[list[int]] = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if less[i][j] and not any(less[i][k] and less[k][j] for k in range(n)):
                up[i].append(j)
                down[j].append(i)
    return up, down


def relabel(up: list[list[int]], down: list[list[int]], perm: list[int]) -> tuple[list[list[int]], list[list[int]]]:
    n = len(up)
    nu: list[list[int]] = [[] for _ in range(n)]
    nd: list[list[int]] = [[] for _ in range(n)]
    for i in range(n):
        nu[perm[i]] = [perm[x] for x in up[i]]
        nd[perm[i]] = [perm[x] for x in down[i]]
    return nu, nd


def hasse_digraph(up: list[list[int]]) -> nx.DiGraph:
    g = nx.DiGraph()
    g.add_nodes_from(range(len(up)))
    g.add_edges_from((a, b) for a, bs in enumerate(up) for b in bs)
    return g


def boolean_digraph(w: int) -> nx.DiGraph:
    g = nx.DiGraph()
    g.add_nodes_from(range(2 ** w))
    g.add_edges_from((m, m | (1 << i)) for m in range(2 ** w) for i in range(w) if not m >> i & 1)
    return g


# ====================================================================== configuracion de variantes


def variant_config(variant: str, c: int | None = None) -> dict[str, Any]:
    base = {"rule": "flag", "sched": "LOCAL", "ext": "ER", "asynchronous": False, "cap": None}
    if variant == "V1R":
        base["sched"] = "RANK"
    elif variant == "V1A":
        base["asynchronous"] = True
    elif variant == "V3":
        base["ext"] = "EL"
    elif variant == "V4":
        base["cap"] = c
    elif variant == "V5":
        base["rule"] = "pair"
    return base


def variant_id(variant: str, c: int | None = None) -> int:
    return VID["V4"] + c if variant == "V4" else VID[variant]


def make_specs(seeds: tuple[int, ...] = (0, 1, 2)) -> list[dict[str, Any]]:
    out = []
    for variant in ("V1", "V1R", "V1A"):
        for w in (2, 3, 4):
            for s in seeds:
                out.append({"variant": variant, "w": w, "seed": s, "c": None, "seed_kind": f"S{w}"})
    for w, kind in ((22, "D2+2"), (23, "S2uS3")):
        for s in seeds:
            out.append({"variant": "V2", "w": w, "seed": s, "c": None, "seed_kind": kind})
    for s in seeds:
        out.append({"variant": "V3", "w": 2, "seed": s, "c": None, "seed_kind": "S2"})
    for c in (2, 3, 4):
        for s in seeds:
            out.append({"variant": "V4", "w": 6, "seed": s, "c": c, "seed_kind": "S6"})
    for s in seeds:
        out.append({"variant": "V5", "w": 3, "seed": s, "c": None, "seed_kind": "S3"})
    for s in seeds:
        out.append({"variant": "V7", "w": 20, "seed": s, "c": None, "seed_kind": "R20"})
    for sp in out:
        sp["id"] = f"{sp['variant']}" + (f"_c{sp['c']}" if sp["c"] else "") + f"_w{sp['w']}_s{sp['seed']}"
    return out


# ====================================================================== observables


def hasse_edges(g: Grower) -> tuple[np.ndarray, np.ndarray]:
    n = g.n
    cnt = np.fromiter((len(d) for d in g.down), dtype=np.int64, count=n)
    v = np.repeat(np.arange(n, dtype=np.int64), cnt)
    u = np.fromiter(itertools.chain.from_iterable(g.down), dtype=np.int64, count=int(cnt.sum()))
    return u, v


def future_set(g: Grower, sources: list[int]) -> set[int]:
    """Elementos estrictamente por encima de TODOS los sources (interseccion de los conos futuros)."""
    res: set[int] | None = None
    for s in sources:
        seen = set()
        st = [s]
        while st:
            v = st.pop()
            for u in g.up[v]:
                if u not in seen:
                    seen.add(u)
                    st.append(u)
        res = seen if res is None else (res & seen)
    return res or set()


def observables(g: Grower, future_sources: list[int] | None) -> dict[str, Any]:
    n = g.n
    ups = np.fromiter((len(u) for u in g.up), dtype=np.int64, count=n)
    dns = np.fromiter((len(d) for d in g.down), dtype=np.int64, count=n)
    interior = ups == dns
    r_loc = float(np.median(ups[interior])) if interior.any() else None
    front = ups < dns
    nfront = int(front.sum())
    namb = int(((dns - ups) >= 2).sum())
    deg = ups + dns
    ranks = np.asarray(g.rank, dtype=np.int64)
    lv = np.bincount(ranks)
    obs: dict[str, Any] = {
        "r_loc": r_loc, "n_interior": int(interior.sum()), "sigma_k": float(deg.std()), "mean_hasse_degree": float(deg.mean()),
        "max_up": int(ups.max()), "n_levels": int(lv.size), "levels_median": float(np.median(lv)), "levels_max": int(lv.max()),
        "levels": lv.tolist() if lv.size <= LEVELS_FULL_MAX else None,
        "n_frontier": nfront, "n_ambiguous_corners": namb,
        "ambiguous_corner_fraction": (namb / nfront) if nfront else None,
        "ambiguous_subjoins": g.ambiguous_subjoins(), "spectral_dimension": None,
    }
    if future_sources is not None:
        fut = future_set(g, future_sources)
        idx = np.fromiter(fut, dtype=np.int64, count=len(fut)) if fut else np.empty(0, dtype=np.int64)
        sel = idx[interior[idx]] if idx.size else idx
        obs["future_size"] = int(idx.size)
        obs["future_r_loc"] = float(np.median(ups[sel])) if sel.size else None
    return obs


# ====================================================================== corrida


def load_tau() -> tuple[float, float]:
    th = json.loads(W6_SUMMARY.read_text())["calibration"]["thresholds"]
    return float(th["tau_d"]), float(th["tau_k"])


def snapshot_record(g: Grower, key: tuple[int, ...], future_sources: list[int] | None, do_measure: bool) -> tuple[dict, dict]:
    t0 = time.perf_counter()
    obs = observables(g, future_sources)
    t1 = time.perf_counter()
    rec: dict[str, Any] = {"n": g.n, "obs": obs, "key": list(key)}
    secs: dict[str, Any] = {"obs": round(t1 - t0, 2)}
    if do_measure:
        u, v = hasse_edges(g)
        adj = _from_edges(g.n, u, v)
        kg, ng = giant_degree(adj)
        m = measure(adj, key)
        m.pop("annulus", None)
        m.pop("coherence", None)
        secs.update(m.pop("seconds"))
        m.update(k_giant=kg, n_giant_check=ng, mean_degree_graph=m["mean_degree"], E_over_N=m["mean_degree"] / 2.0)
        rec["meas"] = m
        secs["total_measure"] = round(time.perf_counter() - t1, 2)
    return rec, secs


def evaluate_run(snaps: list[dict], tau: tuple[float, float]) -> dict[str, Any]:
    """delta, W5, exclusiones por par y W6 con los umbrales congelados."""
    ev: dict[str, Any] = {"delta12": None, "delta23": None, "w6_valid": None, "rc3_excluded": None,
                          "exclusions_pair1": None, "exclusions_pair2": None, "w5_pair1": None, "w5_pair2": None}
    ms = [s["meas"] for s in snaps if "meas" in s]
    if len(ms) >= 2:
        from rc3 import pair_exclusions, scaling_delta
        a, b = ms[0], ms[1]
        ev["delta12"] = scaling_delta(a["R"], b["R"], a["n_giant"] or 0, b["n_giant"] or 0)
        ev["exclusions_pair1"] = pair_exclusions(ev["delta12"], (a["annulus_status"], b["annulus_status"]),
                                                 (a["coherence_status"], b["coherence_status"]))
        ev["w5_pair1"] = bool(abs(b["E_over_N"] - a["E_over_N"]) / a["E_over_N"] < 0.25)
    if len(ms) >= 3:
        from rc3 import pair_exclusions, scaling_delta
        b, c = ms[1], ms[2]
        ev["delta23"] = scaling_delta(b["R"], c["R"], b["n_giant"] or 0, c["n_giant"] or 0)
        ev["exclusions_pair2"] = pair_exclusions(ev["delta23"], (b["annulus_status"], c["annulus_status"]),
                                                 (b["coherence_status"], c["coherence_status"]))
        ev["w5_pair2"] = bool(abs(c["E_over_N"] - b["E_over_N"]) / b["E_over_N"] < 0.25)
        fm = family_metrics([{**m, "R": m["R"], "n_giant": m["n_giant"]} for m in ms[:3]])
        rule = w6_rule(fm, tau[0], tau[1])
        ev["w6"] = rule
        ev["w6_valid"] = bool(rule["valid"])
        ev["k"] = [fm["k1"], fm["k2"], fm["k3"]]
        ev["abs_ddelta"] = fm["abs_ddelta"]
        ev["rel_dk"] = fm["rel_dk"]
        ev["rc3_excluded"] = bool(ev["exclusions_pair1"] or ev["exclusions_pair2"])
        r2, r3 = snaps[1]["obs"]["r_loc"], snaps[2]["obs"]["r_loc"]
        ev["rloc_x_delta12"] = None if None in (r2, ev["delta12"]) else r2 * ev["delta12"]
        ev["rloc_x_delta23"] = None if None in (r3, ev["delta23"]) else r3 * ev["delta23"]
    elif len(ms) == 2:
        ev["rc3_excluded"] = bool(ev["exclusions_pair1"])
    return ev


def run_spec(spec: dict[str, Any], sizes: tuple[int, ...], tau: tuple[float, float], do_measure: bool = True,
             max_seconds: float | None = None) -> tuple[dict, dict]:
    """Una corrida (variante, w, semilla). Devuelve (registro determinista, tiempos)."""
    t_start = time.perf_counter()
    deadline = None if max_seconds is None else t_start + max_seconds
    variant, w, seed, c = spec["variant"], spec["w"], spec["seed"], spec["c"]
    vid = variant_id(variant, c)
    mk = lambda k: rng_from_key((MASTER_CF1, vid, w, seed, k))  # noqa: E731
    rngs = {0: mk(0), 1: mk(1), 3: mk(3)}
    cfg = variant_config(variant, c)
    snaps: list[dict] = []
    secs: dict[str, Any] = {"snapshots": []}
    fut_sources: list[int] | None = None
    coalescence: int | None = None
    coal_down: list[int] | None = None
    pre_status = "COMPLETO"

    def make_growers() -> Grower | None:
        nonlocal fut_sources, coalescence, coal_down, pre_status
        if variant == "V7":
            up, down = seed_R20(mk(2))
            g = Grower(up, down, rngs, **cfg)
            fut_sources = [i for i in range(g.n) if not g.down[i]]
            return g
        if variant == "V2":
            wa, wb = (2, 2) if w == 22 else (2, 3)
            half = sizes[0] // 2
            comps = []
            for ww in (wa, wb):
                up, down = seed_S(ww)
                gg = Grower(up, down, rngs, **cfg)
                st, _ = gg.run([], None, half, deadline)
                if st != "COMPLETO":
                    pre_status = st
                    return None
                comps.append(gg)
            A, B = comps
            off = A.n
            up = [list(x) for x in A.up] + [[y + off for y in x] for x in B.up]
            down = [list(x) for x in A.down] + [[y + off for y in x] for x in B.down]
            maxA = [i for i in range(A.n) if not A.up[i]]
            maxB = [i + off for i in range(B.n) if not B.up[i]]
            ia = maxA[int(rngs[1].integers(len(maxA)))]
            ib = maxB[int(rngs[1].integers(len(maxB)))]
            g = Grower(up, down, rngs, **cfg)
            coal_down = sorted([ia, ib])
            coalescence = g._create(coal_down)
            fut_sources = [coalescence]
            return g
        up, down = seed_S(w)
        g = Grower(up, down, rngs, **cfg)
        fut_sources = [0]
        return g

    g = make_growers()
    if g is None:
        status, n_final = pre_status, 0
    else:
        def on_snap(si: int) -> None:
            key = (MASTER_CF1, vid, w, seed, si)
            rec, sc = snapshot_record(g, key, fut_sources, do_measure)
            rec["size_index"] = si
            rec["target"] = sizes[si]
            snaps.append(rec)
            secs["snapshots"].append(sc)

        tg = time.perf_counter()
        # los snapshots se miden dentro del bucle: el reloj de crecimiento los descuenta al final
        status, _ = g.run(sizes, on_snap, sizes[-1], deadline, track=True)
        if variant == "V2" and not snaps and status == "COMPLETO":
            pass
        if not g.traj or g.traj[-1]["n"] != g.n:
            g._traj_point()
        n_final = g.n
        secs["growth_plus_measure_s"] = round(time.perf_counter() - tg, 1)
    # V2: la instantanea N1 cae justo tras la coalescencia (n = N1 + 1)
    ev = evaluate_run(snaps, tau)
    reached = [s["size_index"] for s in snaps]
    rec = {
        "id": spec["id"], "variant": variant, "w": w, "seed": seed, "c": c, "seed_kind": spec["seed_kind"],
        "key": [MASTER_CF1, vid, w, seed], "status": status, "n_final": n_final, "sizes_target": list(sizes),
        "sizes_reached": reached, "evaluable_N3": 2 in reached,
        "n_complete": None if g is None else g.n_complete, "n_extend": None if g is None else g.n_extend,
        "degenerate": None if g is None else len(g.degenerate), "cap_blocked": None if g is None else len(g.cap_blocked),
        "pending_final": None if g is None else len(g.pend),
        "trajectory": [] if g is None else g.traj,
        "final_state": None if g is None else {"mean_hasse_degree": 2.0 * g.n_edges / g.n, "max_up": g.maxup}, "coalescence_element": coalescence, "coalescence_down": coal_down,
        "snapshots": snaps, **ev,
    }
    secs["total_s"] = round(time.perf_counter() - t_start, 1)
    return rec, secs


# ====================================================================== V0


def v0_block() -> dict[str, Any]:
    out: dict[str, Any] = {"closure": {}, "relabel": {}}
    for w in (7, 8):
        P = cf01.Poset(w, cf01.FLAG, "CAUSAL", 0, cf01.CAP)
        term = P.run()
        b = cf01.boolean_check(P.up, P.down, w) if term else {"isomorphic": False}
        out["closure"][f"S{w}/cf01-CAUSAL"] = {"terminated": term, "n": P.n, "boolean": bool(b["isomorphic"])}
        for sched in ("LOCAL", "RANK"):
            rngs = {0: rng_from_key((MASTER_CF1, 0, w, 0, 0)), 1: rng_from_key((MASTER_CF1, 0, w, 0, 1)),
                    3: rng_from_key((MASTER_CF1, 0, w, 0, 3))}
            up, down = seed_S(w)
            g = Grower(up, down, rngs, sched=sched)
            term = g.closure(5000)
            iso = term and nx.is_isomorphic(hasse_digraph(g.up), boolean_digraph(w))
            out["closure"][f"S{w}/{sched}"] = {"terminated": term, "n": g.n, "boolean": bool(iso)}
    for w in (3, 4, 5, 6):
        rng = rng_from_key((MASTER_CF1, 0, w, 0, 2))
        perm = [int(x) for x in rng.permutation(w + 1)]
        res = {}
        for sched in ("LOCAL", "RANK"):
            rngs = {0: rng_from_key((MASTER_CF1, 0, w, 1, 0)), 1: rng_from_key((MASTER_CF1, 0, w, 1, 1)),
                    3: rng_from_key((MASTER_CF1, 0, w, 1, 3))}
            up, down = seed_S(w)
            g0 = Grower(up, down, rngs, sched=sched)
            g0.closure(5000)
            pu, pd = relabel(*seed_S(w), perm)
            g1 = Grower(pu, pd, rngs, sched=sched)
            g1.closure(5000)
            res[sched] = {"n": [g0.n, g1.n],
                          "iso_permuted_vs_original": bool(nx.is_isomorphic(hasse_digraph(g0.up), hasse_digraph(g1.up))),
                          "iso_permuted_vs_Bw": bool(nx.is_isomorphic(hasse_digraph(g1.up), boolean_digraph(w)))}
        out["relabel"][f"S{w}"] = {"perm": perm, **res}
    out["ok"] = (all(v["boolean"] for v in out["closure"].values())
                 and all(r["iso_permuted_vs_original"] and r["iso_permuted_vs_Bw"]
                         for s in out["relabel"].values() for r in (s["LOCAL"], s["RANK"])))
    return out


# ====================================================================== criterios §5


def _mean(xs: list[float]) -> float | None:
    xs = [x for x in xs if x is not None]
    return float(statistics.fmean(xs)) if xs else None


def group_stats(runs: list[dict]) -> dict[str, Any]:
    n = len(runs)
    valid = [r["w6_valid"] is True for r in runs]
    excl = [bool(r["rc3_excluded"]) if r["rc3_excluded"] is not None else None for r in runs]
    qual = [v and e is False for v, e in zip(valid, excl)]
    return {
        "n_seeds": n, "n_evaluable_N3": sum(bool(r["evaluable_N3"]) for r in runs), "n_w6_valid": sum(valid),
        "n_excluded": sum(e is True for e in excl), "n_not_excluded": sum(e is False for e in excl),
        "n_qualify": sum(qual), "delta12_mean": _mean([r["delta12"] for r in runs]),
        "delta23_mean": _mean([r["delta23"] for r in runs]), "status": [r["status"] for r in runs],
        "a_positive": sum(qual) >= 2, "a_positive_strict": sum(valid) >= 2 and not any(e is True for e in excl),
        "w6_class": sum(valid) >= 2, "rc3_class": sum(e is False for e in excl) >= 2,
        "seeds": [{"seed": r["seed"], "status": r["status"], "n_final": r["n_final"], "w6_valid": r["w6_valid"],
                   "delta12": r["delta12"], "delta23": r["delta23"], "excl1": r["exclusions_pair1"],
                   "excl2": r["exclusions_pair2"], "w5": [r["w5_pair1"], r["w5_pair2"]],
                   "k": r.get("k"), "rloc_N": [s["obs"]["r_loc"] for s in r["snapshots"]],
                   "rloc_x_delta": [r.get("rloc_x_delta12"), r.get("rloc_x_delta23")],
                   "corner_frac": [s["obs"]["ambiguous_corner_fraction"] for s in r["snapshots"]],
                   "future_rloc": [s["obs"].get("future_r_loc") for s in r["snapshots"]]} for r in runs],
    }


def evaluate_all(runs: list[dict]) -> dict[str, Any]:
    groups: dict[tuple, list[dict]] = {}
    for r in runs:
        groups.setdefault((r["variant"], r["w"], r["c"]), []).append(r)
    gs = {k: group_stats(v) for k, v in groups.items()}
    G = lambda v, w, c=None: gs.get((v, w, c))  # noqa: E731
    out: dict[str, Any] = {"groups": {f"{k[0]}|w={k[1]}|c={k[2]}": v for k, v in gs.items()}}
    # robustez / trazabilidad
    rob: dict[int, Any] = {}
    for w in (2, 3, 4):
        base = G("V1", w)
        if base is None:
            continue
        entry: dict[str, Any] = {}
        ok = True
        for x in ("V1R", "V1A"):
            o = G(x, w)
            if o is None:
                entry[x] = None
                ok = False
                continue
            dd = None if None in (base["delta12_mean"], o["delta12_mean"]) else abs(base["delta12_mean"] - o["delta12_mean"])
            same = base["w6_class"] == o["w6_class"] and base["rc3_class"] == o["rc3_class"]
            entry[x] = {"same_classes": same, "abs_diff_delta12": dd, "ok": bool(same and dd is not None and dd <= 0.04)}
            ok = ok and entry[x]["ok"]
        entry["robust"] = ok
        rob[w] = entry
    out["order_robustness"] = rob
    tr = {w: (None if (G("V1", w) is None or G("V1", w)["delta12_mean"] is None)
              else abs(G("V1", w)["delta12_mean"] - 1.0 / w) <= 0.06) for w in (2, 3, 4)}
    out["traceability_inicio"] = {"per_w": tr, "inicio": all(v is True for v in tr.values())}
    # calibraciones
    k1 = G("V3", 2)
    k3 = G("V5", 3)
    k2c = {c: G("V4", 6, c) for c in (2, 3, 4)}
    cal = {
        "K1": None if k1 is None else bool(k1["n_seeds"] - k1["n_w6_valid"] >= 2),
        "K2": None if any(v is None for v in k2c.values()) else bool(all(
            v["delta12_mean"] is not None and abs(v["delta12_mean"] - 1.0 / c) <= 0.08 for c, v in k2c.items())),
        "K3": None if k3 is None else bool(k3["n_seeds"] - k3["n_w6_valid"] >= 2),
    }
    cal["K2_detail"] = {str(c): None if v is None else v["delta12_mean"] for c, v in k2c.items()}
    out["calibration"] = cal
    cal_ok = all(cal[k] is True for k in ("K1", "K2", "K3"))
    # desenlace
    v1_pos = [w for w in (2, 3, 4) if G("V1", w) and G("V1", w)["a_positive"]]
    v1_pos_rob = [w for w in v1_pos if rob.get(w, {}).get("robust")]
    any_planner = [(v, w) for v in ("V1", "V1R", "V1A") for w in (2, 3, 4) if G(v, w) and G(v, w)["a_positive"]]
    dw = {w: (None if (G("V1", w) is None or not G("V1", w)["delta12_mean"]) else 1.0 / G("V1", w)["delta12_mean"]) for w in (2, 3, 4)}
    dvals = [v for v in dw.values() if v is not None]
    d1 = (not out["traceability_inicio"]["inicio"]) and len(dvals) == 3 and (max(dvals) - min(dvals) <= 0.5)
    if not cal_ok:
        outcome = "NO_INTERPRETABLE (K1-K3 no cumplidas o incompletas)"
    elif v1_pos_rob:
        outcome = "CF1-A+"
    elif any_planner:
        outcome = "CF1-A±"
    else:
        outcome = "CF1-neg"
    out["outcome"] = {"outcome": outcome, "V1_A_positive_w": v1_pos, "V1_A_positive_robust_w": v1_pos_rob,
                      "A_positive_any_planner": any_planner, "d_traced_by_w": dw, "D1": bool(d1),
                      "calibrations_ok": cal_ok,
                      "note": "con menos de 3 semillas/bloques incompletos el desenlace no es leible (humo)"}
    out["V2_V7_descriptive"] = {k: v for k, v in out["groups"].items() if k.startswith(("V2|", "V7|"))}
    return out


# ====================================================================== E/S


def _load(path: Path) -> dict[str, dict]:
    done: dict[str, dict] = {}
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                try:
                    r = json.loads(line)
                    done[r["id"]] = r
                except json.JSONDecodeError:
                    pass
    return done


_CTX: dict[str, Any] = {}


def _job(spec: dict) -> tuple[dict, dict]:
    return run_spec(spec, _CTX["sizes"], _CTX["tau"], True, _CTX["max_seconds"])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smoke", action="store_true", help="N<=2000, 1 semilla -> results/cf1_smoke")
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--fresh", action="store_true")
    ap.add_argument("--sizes", type=int, nargs=3, default=None)
    ap.add_argument("--seeds", type=int, default=None, help="n.o de semillas por bloque (defecto 3; humo 1)")
    ap.add_argument("--only", default=None, help="lista de variantes separadas por coma (p.ej. V1,V3)")
    ap.add_argument("--max-seconds", type=float, default=None)
    ap.add_argument("--no-v0", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    sizes = tuple(a.sizes) if a.sizes else (SMOKE_SIZES if a.smoke else SIZES)
    nseeds = a.seeds or (1 if a.smoke else 3)
    out = Path(a.out) if a.out else (OUT_SMOKE if a.smoke else OUT_FULL)
    out.mkdir(parents=True, exist_ok=True)
    tau = load_tau()
    rj = out / "runs.jsonl"
    if a.fresh and rj.exists():
        rj.unlink()
    done = _load(rj)
    specs = make_specs(tuple(range(nseeds)))
    if a.only:
        keep = set(a.only.split(","))
        specs = [s for s in specs if s["variant"] in keep]
    todo = [s for s in specs if s["id"] not in done]
    _CTX.update(sizes=sizes, tau=tau, max_seconds=a.max_seconds if a.max_seconds else (None if a.smoke else 6 * 3600))
    t0 = time.perf_counter()
    times: dict[str, Any] = {}
    print(f"{len(specs)} corridas ({len(done)} hechas), sizes={sizes}, tau={tau}", flush=True)
    sfile = out / "summary.json"
    prev_times = json.loads(sfile.read_text()).get("timing", {}) if sfile.exists() and not a.fresh else {}
    times.update(prev_times.get("runs", {}))
    with mp.get_context("fork").Pool(a.procs) as pool, open(rj, "a") as fh:
        for rec, sc in pool.imap_unordered(_job, todo, chunksize=1):
            fh.write(json.dumps(rec, default=str) + "\n")
            fh.flush()
            done[rec["id"]] = rec
            times[rec["id"]] = sc
            print(f"  {rec['id']:<18} {rec['status']:<18} n={rec['n_final']} d12={rec['delta12']} d23={rec['delta23']} "
                  f"w6={rec['w6_valid']} excl={rec['rc3_excluded']} {sc['total_s']}s", flush=True)
    runs = [done[s["id"]] for s in specs]
    rj.write_text("".join(json.dumps(r, default=str) + "\n" for r in runs))
    v0 = None if a.no_v0 else v0_block()
    summary = {"step": "cf1", "master": MASTER_CF1, "sizes": list(sizes), "tau_d_tau_k": list(tau), "n_runs": len(runs),
               "V0": v0, "evaluation": evaluate_all(runs),
               "timing": {"seconds_this_invocation": round(time.perf_counter() - t0, 1), "runs": times}}
    sfile.write_text(json.dumps(summary, indent=1, default=str))
    ev = summary["evaluation"]
    print("V0 ok:", None if v0 is None else v0["ok"])
    print("calibration:", ev["calibration"], "| outcome:", ev["outcome"]["outcome"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
