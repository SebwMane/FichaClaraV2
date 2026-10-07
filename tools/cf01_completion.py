"""CF-0.1: formalizacion de la completitud, calculo exacto sin dinamica de crecimiento (docs/OMEGA_CF0.md §8, prerregistro congelado).

Semilla: raiz (id 0) con w cubiertas superiores (ids 1..w, los atomos). Cada elemento nuevo se define por su conjunto de cubiertas
inferiores. Reglas: C2, C_k (aridades 3..k), C_flag (k = |up(z)|, sin cota; k_code = 99). Ordenes ALTA / BAJA / RONDA.
Barrido w=2..6, k in {2,3,4,5,flag}, 3 ordenes, perm 0 (desempate por id) + perm 1..3 (aleatorias), limite CAP elementos.
Salida: results/cf01/runs.jsonl, results/cf01/summary.json. Los tiempos solo van a summary.json (runs.jsonl es determinista).

Decisiones de implementacion (lectura mas literal de las ambiguedades)
1. Sub-supremo J(T), T subconjunto de up(z): |T|=1 -> el propio elemento; |T|>=2 -> el elemento u de menor id que cubre J(T\\{t}) para
   TODO t en T (u en interseccion de up(J(T\\{t}))). Si no hay ninguno, J(T) no existe. J es global (no depende de z). Los J positivos
   son estables (los ids nuevos son mayores), por eso se cachean; los negativos se invalidan en cada creacion.
2. Instancia de aridad m (2<=m<=k, o <=|up(z)| en flag): (z, S), S subconjunto de m elementos de up(z), con J(T) existente para todo T
   subconjunto de S de tamano m-1 (m=2: trivial) y J(S) inexistente. Esto equivale a "los m sub-supremos no tienen cubierta superior
   comun" porque J(S) es justamente una cubierta comun de los J(S\\{s}). Para m=2 es C2 (up(x) y up(y) disjuntos).
3. El elemento nuevo cubre exactamente {J(S\\{s}) : s en S}. Si esos m sub-supremos no son distintos entre si (coinciden), la
   instancia se descarta (no es un m-cubo genuino); se cuenta en 'degenerate_skipped'.
4. Ambiguedad: se informa 'ambiguous_joins' = numero de J(T) cacheados (|T|>=2) que en el estado FINAL tienen mas de un candidato.
5. Ordenes. ALTA: en cada paso una sola instancia, la de mayor aridad disponible. BAJA: la de menor aridad. Empates por prioridad:
   perm 0 -> clave (id z, ids de S ordenados); perm p>=1 -> cada elemento recibe, al crearse, una prioridad aleatoria de
   rng_from_key((20261024, w, k_code, p)) y la clave usa esas prioridades (permutacion de identificadores). La eleccion de J (menor
   id numerico) no se permuta.
   RONDA: se calculan todas las instancias del estado actual, se deduplican por conjunto de cubiertas inferiores (dos instancias que
   piden el mismo conjunto crean UN elemento) y se aplican todas sin revalidar (por construccion no interfieren: cada una
   crea un elemento propio). Los ids nuevos se asignan en orden de clave.
6. Limite: CAP=5000 elementos. En RONDA la ronda se trunca exactamente al llegar a CAP. Si al llegar a CAP aun hay instancias, no
   termina. (Si no queda ninguna con n==CAP, termina.)
7. Booleano: se verifica un isomorfismo real: mascara de atomos bajo cada elemento; 2^w elementos, mascaras biyectivas sobre todos los
   subconjuntos, relacion <= (cierre transitivo) coincide con inclusion de mascaras en ambos sentidos, y las aristas de cobertura
   guardadas son exactamente las de B_w (mascara v = mascara u + un atomo). Se informan ademas techo unico y rangos binomiales.
8. Rango = cadena mas larga desde la raiz. Certificado de periodicidad (solo no terminantes, perm 0): config(r) = subgrafo dirigido
   inducido (relaciones de cobertura) sobre los niveles r y r+1, con atributos de nodo (nivel relativo 0/1, n.o de cubiertas inferiores
   fuera de la ventana, n.o de cubiertas superiores fuera de la ventana); se compara config(r) con config(r+1) por isomorfismo (hash WL
   + networkx). Solo se prueban r <= rango_max - 4 para no usar niveles truncados por el limite. Es un certificado empirico (el
   limite trunca el proceso); se informa el primer r con repeticion, si lo hay.
9. Evaluacion de predicciones/lectura: 'cierra' = termina y es isomorfo a B_w. La lectura CF01-a/b/c se evalua con los ordenes
   ALTA y RONDA (los de la prediccion de cierre) y se informa la dependencia del orden por separado.
"""

from __future__ import annotations

import argparse
import itertools
import json
import multiprocessing as mp
import time
import warnings
from math import comb
from pathlib import Path

import networkx as nx

from omega.c0.references import rng_from_key

warnings.filterwarnings("ignore", message="The hashes produced for directed graphs")
MASTER_CF = 20261024
CAP = 5000
FLAG = 99
KS = (2, 3, 4, 5, FLAG)
ORDERS = ("ALTA", "BAJA", "RONDA")
WS = (2, 3, 4, 5, 6)
PERMS = (0, 1, 2, 3)


class Poset:
    def __init__(self, w: int, k: int, order: str, perm: int, cap: int = CAP):
        self.w, self.k, self.order, self.perm, self.cap = w, k, order, perm, cap
        self.up: list[set[int]] = [set()]
        self.down: list[set[int]] = [set()]
        self.rank: list[int] = [0]
        self.rng = rng_from_key((MASTER_CF, w, k, perm)) if perm > 0 else None
        self.prio: list[float] = [0.0]
        self.jc: dict[frozenset, int] = {}
        self.jneg: set[frozenset] = set()
        self.inst: dict[int, dict[tuple, tuple]] = {}  # z -> {S: key}
        self.by_arity: dict[int, dict[tuple, tuple]] = {}  # m -> {(z,S): key}
        self.degenerate = 0
        self.steps = 0
        for _ in range(w):
            self.new(frozenset({0}))

    @property
    def n(self) -> int:
        return len(self.up)

    def new(self, lower: frozenset) -> int:
        u = len(self.up)
        self.up.append(set())
        self.down.append(set(lower))
        self.rank.append(1 + max(self.rank[d] for d in lower))
        self.prio.append(float(u) if self.rng is None else float(self.rng.random()))
        for d in lower:
            self.up[d].add(u)
        self.jneg.clear()
        return u

    # ---- sub-supremos ----
    def cands(self, T: frozenset) -> set[int] | None:
        sub = []
        for t in T:
            j = self.J(T - {t})
            if j is None:
                return None
            sub.append(j)
        c = set(self.up[sub[0]])
        for j in sub[1:]:
            c &= self.up[j]
        return c

    def J(self, T: frozenset) -> int | None:
        if len(T) == 1:
            return next(iter(T))
        if T in self.jc:
            return self.jc[T]
        if T in self.jneg:
            return None
        c = self.cands(T)
        if c:
            r = min(c)
            self.jc[T] = r
            return r
        self.jneg.add(T)
        return None

    # ---- instancias ----
    def key(self, z: int, S: tuple) -> tuple:
        return (self.prio[z], tuple(sorted(self.prio[s] for s in S)))

    def compute_z(self, z: int) -> dict[tuple, tuple]:
        ups = sorted(self.up[z])
        res: dict[tuple, tuple] = {}
        if len(ups) < 2:
            return res
        kmax = len(ups) if self.k == FLAG else self.k
        A: set[tuple] = set()
        for x, y in itertools.combinations(ups, 2):
            if self.J(frozenset((x, y))) is not None:
                A.add((x, y))
            else:
                res[(x, y)] = self.key(z, (x, y))
        m = 3
        while m <= kmax and A:
            A2: set[tuple] = set()
            for S in A:
                for e in ups:
                    if e <= S[-1]:
                        continue
                    S2 = S + (e,)
                    if all(S2[:i] + S2[i + 1:] in A for i in range(m - 1)):
                        if self.J(frozenset(S2)) is not None:
                            A2.add(S2)
                        else:
                            res[S2] = self.key(z, S2)
            A = A2
            m += 1
        return res

    def set_z(self, z: int) -> None:
        for S in self.inst.get(z, {}):
            self.by_arity[len(S)].pop((z, S), None)
        new = self.compute_z(z)
        self.inst[z] = new
        for S, kk in new.items():
            self.by_arity.setdefault(len(S), {})[(z, S)] = kk

    def below(self, u: int) -> set[int]:
        seen, st = set(), [u]
        while st:
            v = st.pop()
            for d in self.down[v]:
                if d not in seen:
                    seen.add(d)
                    st.append(d)
        return seen

    def full_instances(self) -> dict[tuple, tuple]:
        """Recalculo completo (referencia para auditar el mantenimiento incremental)."""
        out = {}
        for z in range(self.n):
            for S, kk in self.compute_z(z).items():
                out[(z, S)] = kk
        return out

    def lower_of(self, S: tuple) -> frozenset | None:
        L = set()
        for s in S:
            j = self.J(frozenset(S) - {s})
            L.add(j)
        if len(L) < len(S):
            return None
        return frozenset(L)

    def apply(self, lowers: list[frozenset]) -> set[int]:
        dirty: set[int] = set()
        for L in lowers:
            if self.n >= self.cap:
                break
            u = self.new(L)
            bu = self.below(u)
            # instancias en z solo cambian si up(z) cambia (z en L) o si algun J(T), T subconjunto de up(z), pasa a existir;
            # entonces |T|>=2 y todo t en T esta estrictamente bajo u: basta |up(z) ∩ bajo(u)| >= 2.
            dirty |= set(L)
            dirty |= {z for z in bu if len(self.up[z] & bu) >= 2}
        return dirty

    def run(self) -> bool:
        for z in range(self.n):
            self.set_z(z)
        while True:
            arities = sorted(m for m, d in self.by_arity.items() if d)
            if not arities:
                return True
            if self.n >= self.cap:
                return False
            self.steps += 1
            if self.order == "RONDA":
                items = sorted((kk, zs) for m in arities for zs, kk in self.by_arity[m].items())
                seen, lowers = set(), []
                for _, (z, S) in items:
                    L = self.lower_of(S)
                    if L is None:
                        self.degenerate += 1
                    elif L not in seen:
                        seen.add(L)
                        lowers.append(L)
                n0 = self.n
                dirty = self.apply(lowers)
                if self.n == n0:  # solo instancias degeneradas: sin progreso posible
                    return True
            else:
                chosen = None
                for m in (reversed(arities) if self.order == "ALTA" else arities):
                    for zs, kk in sorted(self.by_arity[m].items(), key=lambda it: it[1]):
                        L = self.lower_of(zs[1])
                        if L is None:
                            self.degenerate += 1
                            continue
                        chosen = L
                        break
                    if chosen is not None:
                        break
                if chosen is None:
                    return True
                dirty = self.apply([chosen])
            dirty |= set(range(len(self.inst), self.n))
            for z in sorted(dirty):
                self.set_z(z)

    # ---- medidas ----
    def maximal(self) -> int:
        return sum(1 for u in self.up if not u)

    def levels(self) -> list[int]:
        c = [0] * (max(self.rank) + 1)
        for r in self.rank:
            c[r] += 1
        return c

    def ambiguous(self) -> int:
        return sum(1 for T in self.jc if len(T) >= 2 and len(self.cands(T) or ()) > 1)


def boolean_check(up: list[set[int]], down: list[set[int]], w: int) -> dict:
    n = len(up)
    out = {"n_ok": n == 2**w, "one_top": sum(1 for u in up if not u) == 1}
    atoms = sorted(up[0])
    out["atoms_ok"] = len(atoms) == w
    if n != 2**w or len(atoms) != w:
        out["isomorphic"] = False
        return out
    bit = {a: 1 << i for i, a in enumerate(atoms)}
    mask = [0] * n
    for v in range(1, n):  # ids crecientes: las cubiertas inferiores tienen id menor
        mask[v] = bit[v] if v in bit else 0
        if v not in bit:
            for d in down[v]:
                mask[v] |= mask[d]
    out["bijective"] = len(set(mask)) == n
    anc = [0] * n  # ancestros reflexivos como bitset de elementos
    for v in range(n):
        anc[v] = 1 << v
        for d in down[v]:
            anc[v] |= anc[d]
    order_ok = all(((anc[v] >> u) & 1) == (mask[u] & ~mask[v] == 0) for u in range(n) for v in range(n))
    cov = {(d, v) for v in range(n) for d in down[v]}
    cov_b = {(u, v) for u in range(n) for v in range(n) if (mask[v] & ~mask[u]) and bin(mask[v] ^ mask[u]).count("1") == 1 and mask[u] & ~mask[v] == 0}
    out["order_iso"] = bool(order_ok)
    out["covers_equal_Bw"] = cov == cov_b
    ranks = [0] * (w + 1)
    for m in mask:
        ranks[bin(m).count("1")] += 1
    out["binomial_ranks"] = ranks == [comb(w, i) for i in range(w + 1)]
    out["isomorphic"] = bool(out["bijective"] and order_ok and out["covers_equal_Bw"])
    return out


def periodicity(P: Poset) -> dict:
    R = max(P.rank)
    lv: dict[int, list[int]] = {}
    for v, r in enumerate(P.rank):
        lv.setdefault(r, []).append(v)

    def cfg(r: int) -> nx.DiGraph:
        nodes = lv.get(r, []) + lv.get(r + 1, [])
        S = set(nodes)
        g = nx.DiGraph()
        for v in nodes:
            g.add_node(v, a=f"{P.rank[v] - r}|{sum(1 for d in P.down[v] if d not in S)}|{sum(1 for u in P.up[v] if u not in S)}")
        for v in nodes:
            for d in P.down[v]:
                if d in S:
                    g.add_edge(d, v)
        return g

    first = None
    tested = 0
    for r in range(1, R - 3):
        g1, g2 = cfg(r), cfg(r + 1)
        tested = r
        if g1.number_of_nodes() != g2.number_of_nodes() or g1.number_of_edges() != g2.number_of_edges():
            continue
        h1 = nx.weisfeiler_lehman_graph_hash(g1, node_attr="a")
        h2 = nx.weisfeiler_lehman_graph_hash(g2, node_attr="a")
        if h1 == h2 and nx.is_isomorphic(g1, g2, node_match=lambda a, b: a["a"] == b["a"]):
            first = r
            break
    return {"max_rank": R, "last_r_tested": tested, "first_repeat_r": first, "periodic": first is not None}


def run_one(w: int, k: int, order: str, perm: int, cap: int = CAP) -> tuple[dict, float]:
    t0 = time.perf_counter()
    P = Poset(w, k, order, perm, cap)
    term = P.run()
    rec = {
        "w": w, "k_code": k, "order": order, "perm": perm, "terminated": term, "n": P.n, "maximal": P.maximal(),
        "steps": P.steps, "degenerate_skipped": P.degenerate, "ambiguous_joins": P.ambiguous(),
        "levels": P.levels(),
    }
    b = boolean_check(P.up, P.down, w) if term else {"isomorphic": False}
    rec["boolean"] = bool(b["isomorphic"])
    rec["boolean_detail"] = b
    if not term and perm == 0:
        rec["periodicity"] = periodicity(P)
    return rec, time.perf_counter() - t0


def evaluate(recs: list[dict]) -> dict:
    d = {(r["w"], r["k_code"], r["order"], r["perm"]): r for r in recs}
    ws = sorted({r["w"] for r in recs})
    ks = sorted({r["k_code"] for r in recs})
    closes = lambda w, k, o, p=0: d[(w, k, o, p)]["terminated"] and d[(w, k, o, p)]["boolean"]  # noqa: E731
    pred: dict = {}
    pred["C2_w2_cierra_B2"] = all(closes(2, 2, o) for o in ORDERS) if 2 in ws else None
    pred["C2_w>=3_no_termina"] = {o: [w for w in ws if w >= 3 and d[(w, 2, o, 0)]["terminated"]] for o in ORDERS}
    pred["C2_w3_certificado_periodicidad"] = {
        o: d[(3, 2, o, 0)].get("periodicity") for o in ORDERS if (3, 2, o, 0) in d}
    exp_k = {f"{k}/{o}": [(w, closes(w, k, o)) for w in ws if closes(w, k, o) != (w <= k)]
             for k in ks if k != FLAG for o in ("ALTA", "RONDA")}
    pred["Ck_ALTA_RONDA_cierra_sii_w<=k"] = {"violaciones_(w,cierra)": exp_k, "ok": all(not v for v in exp_k.values())}
    fl = {o: [w for w in ws if not closes(w, FLAG, o)] for o in ("ALTA", "RONDA")}
    pred["Cflag_ALTA_RONDA_cierra_todo_w"] = {"w_que_fallan": fl, "ok": all(not v for v in fl.values())}
    baja = {k: closes(3, k, "BAJA") for k in ks}
    pred["BAJA_falla_w3_todo_k"] = {"cierra_en_w3_por_k": baja, "ok": not any(baja.values()),
                                   "BAJA_cierra_otros": {f"w{w}/k{k}": closes(w, k, "BAJA") for w in ws for k in ks if w != 3 and closes(w, k, "BAJA")}}
    varies = []
    for w in ws:
        for k in ks:
            for o in ORDERS:
                vs = {(d[(w, k, o, p)]["terminated"], d[(w, k, o, p)]["boolean"]) for p in PERMS if (w, k, o, p) in d}
                if len(vs) > 1:
                    varies.append((w, k, o))
    detail_var = []
    for w in ws:
        for k in ks:
            for o in ORDERS:
                ns = {(d[(w, k, o, p)]["n"], d[(w, k, o, p)]["maximal"]) for p in PERMS if (w, k, o, p) in d}
                if len(ns) > 1:
                    detail_var.append((w, k, o, sorted(ns)))
    pred["Permutaciones_no_cambian_veredicto"] = {"ok": not varies, "cambian_veredicto": varies,
                                                 "cambian_n_o_maximales_(solo_informativo)": detail_var}
    # lecturas
    read: dict = {}
    for o in ("ALTA", "RONDA"):
        bounded_all = [k for k in ks if k != FLAG and all(closes(w, k, o) for w in ws)]
        read[o] = {
            "C_k_acotadas_que_cierran_todo_w": bounded_all,
            "C_flag_cierra_todo_w": all(closes(w, FLAG, o) for w in ws),
            "CF01-a": (not bounded_all) and all(closes(w, FLAG, o) for w in ws),
            "CF01-b": all(closes(w, 3, o) for w in ws),
            "CF01-c": any(not closes(w, FLAG, o) for w in ws),
        }
    dep = [(w, k) for w in ws for k in ks if len({closes(w, k, o) for o in ORDERS}) > 1]
    read["dependencia_de_orden"] = {"hay": bool(dep), "casos_(w,k)": dep}
    return {"predicciones": pred, "lectura": read}


def _job(args: tuple) -> tuple[dict, float]:
    return run_one(*args)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--wmax", type=int, default=6)
    ap.add_argument("--cap", type=int, default=CAP)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent.parent / "results" / "cf01"))
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    recs, times = [], {}
    t0 = time.perf_counter()
    jobs = [(w, k, o, p, a.cap) for w in range(2, a.wmax + 1) for k in KS for o in ORDERS for p in PERMS]
    with mp.Pool(a.jobs) as pool:  # cada corrida es determinista e independiente; imap conserva el orden
        for (w, k, o, p, _), (r, dt) in zip(jobs, pool.imap(_job, jobs)):
            recs.append(r)
            times[f"w{w}/k{k}/{o}/p{p}"] = round(dt, 3)
            print(w, k, o, p, r["terminated"], r["n"], r["maximal"], r["boolean"], f"{dt:.1f}s", flush=True)
    (out / "runs.jsonl").write_text("".join(json.dumps(r, sort_keys=True) + "\n" for r in recs))
    table = [{k: r[k] for k in ("w", "k_code", "order", "terminated", "n", "maximal", "boolean", "ambiguous_joins")}
             for r in recs if r["perm"] == 0]
    summ = {"cap": a.cap, "master": MASTER_CF, "table_perm0": table, **evaluate(recs),
            "ambiguous_joins_total": sum(r["ambiguous_joins"] for r in recs),
            "degenerate_skipped_total": sum(r["degenerate_skipped"] for r in recs),
            "runtime_total_s": round(time.perf_counter() - t0, 1), "runtime_per_run_s": times}
    (out / "summary.json").write_text(json.dumps(summ, indent=1, sort_keys=True, default=str))
    print(json.dumps({"predicciones": summ["predicciones"], "lectura": summ["lectura"]}, indent=1, default=str)[:6000])


if __name__ == "__main__":
    main()
