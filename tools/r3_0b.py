"""OMEGA R3-0b: verificacion estatica con el juez RC-3 sobre espacios de cortes L(J) (docs/OMEGA_R3_0.md §1 y §5, congelado).

Salida: results/r3_0b/{pairs.jsonl, summary.json}. Uso:
  python tools/r3_0b.py [--procs 3] [--fresh] [--cap 2000000]
Cada instancia es un PAR (L(J_N), L(J_8N)); clave RNG = (MASTER_R3, id_familia, semilla, indice_de_tamano), la MISMA clave
se pasa a rc3.measure (que deriva key+(1,), key+(3,), key+(5,) internamente). Es REANUDABLE por par (pairs.jsonl).

Decisiones de implementacion (lectura mas literal de las ambiguedades)
----------------------------------------------------------------------
* Convencion de n para delta: igual que rc3.run_pair, n = n_giant devuelto por measure (nodos de la componente gigante;
  L(J) es conexo, luego n_giant = |L|). R = rc3.half_mass_radius sobre el perfil de bolas. delta = rc3.scaling_delta.
* "n" de la tabla es el numero de elementos de J. Para J1/J2 el parametro buscado es la longitud de cadena n_c (n = w*n_c);
  para J3/J4 se busca n. Se elige el MENOR n con |L| >= objetivo (1e4 indice 0; 8e4 indice 1) por busqueda exponencial
  (n = 1, 2, 4, ...) y biseccion entre el ultimo n fallido y el primero logrado. Los dos tamanos (indices 0 y 1) se
  buscan de forma independiente, cada uno con su clave (por tanto J3/J4 usan posets independientes en los dos tamanos).
* Determinismo del poset para un n dado: el poset de J3/J4 se regenera por prueba con la MISMA clave; las variables
  aleatorias se consumen en orden (J3: rng.random((n,2)); J4: para j=0..n-1 una fila rng.random(j)), de modo que el
  poset de n elementos es el subposet inducido por los n primeros del de n+1. Consecuencia: |L| es monotono en n y la
  biseccion es exacta (|L(P')| <= |L(P)| para P' inducido: D -> down-cierre es inyectiva).
* Representacion para enumerar: predecesores como TUPLAS de indices (no bitmasks) con la relacion de COBERTURA
  (reduccion transitiva) para J1/J2/J4 y todos los predecesores para J3 (n pequeno). Enumerar down-sets con la cobertura es
  equivalente a hacerlo con la clausura (un down-set que contiene los predecesores inmediatos contiene toda la clausura).
  Los posets J2/J4 se construyen con su clausura transitiva en bitmasks (poset_j2/poset_j4, probados), y `covers`
  extrae la reduccion (requiere etiquetas topologicas: J2 id=k*w+i, J4 labels 0..n-1). Motivo: J4(0.5) necesita n ~ 2e4
  elementos y J1(1) 8e4; la clausura densa no escala en la enumeracion.
* J3: orden producto estricto (x<y sii ambas coordenadas estrictamente menores). J4: para i<j, i<j con prob. p y cierre
  transitivo. J2: cadenas x_{i,k} (id = k*w+i), relacion x_{i,k} < x_{j,k+2} para todo i!=j mas el orden de cada cadena,
  con cierre transitivo. J1: las w cadenas con SOLO la relacion de cobertura (predecesor inmediato): para enumerar
  down-sets es equivalente a la clausura (D es down-set => contiene la clausura) y evita n^2/2 bits (J1(1) llega a 8e4
  elementos).
* Enumeracion BFS desde el vacio, anadiendo elementos minimales del complemento (todos sus predecesores en D). Los
  conjuntos de addables se actualizan de forma incremental (sucesores de x ya sin predecesores pendientes); es el mismo
  conjunto de minimales. Down-sets = ints de Python, dict -> indice. Arista D ~ D U {x} una vez por (D,x).
* Tope duro CAP (defecto 2_000_000 down-sets): al superarlo se aborta limpiamente y se registra status CAP_EXCEEDED
  (explosion de tamano; dato, no error). Si la busqueda no alcanza el objetivo con n <= N_MAX = 200000: SEARCH_FAIL.
  En las pruebas de la busqueda se aborta en cuanto |L| >= objetivo (no hace falta el recuento completo).
* "Fuera de [objetivo, 1.5*objetivo]" se registra por tamano (out_of_band) y se mide igual.
* W5: validez si |E/N(8N) - E/N(N)| / (E/N(N)) < 0.25 con E/N = mean_degree/2 de measure (aristas/nodos del grafo L(J)).
* Prediccion J1(4) 'X4-marginal': match si X4 esta en exclusiones o exactamente un tamano es INCOHERENTE.
* Prediccion J3 'X1 o abstencion con senales de explosion': match si X1 en exclusiones, o (no excluido y abstencion
  [delta None o algun estado NO_EVALUABLE/SIN_VENTANA] y explosion). Explosion (|L| superpolinomico en n): el exponente
  local ln|L|/ln n crece mas de un 25 % entre el tamano N y el 8N (para |L| = n^k seria constante; para e^{cn} crece).
  Se registra k_local1, k_local8. J4(0.02): INCIERTO (match = None); si no excluido se marca residuo a inspeccionar.
* Condicion (i): |delta - 1/w| <= 0.08 con semilla 0. (iii): par con W5 valido, no excluido y delta en [0.2,0.6] (J3, J4).
  Pares sin delta (fallos de tamano / NO_EVALUABLE) no cuentan como no excluidos con delta en el rango.
* Familia 'seed': J1/J2 semilla 0 (deterministas); J3 y J4 semillas 0,1,2.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from collections import deque  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402
from scipy import sparse  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from c0_dynamics import _atomic_write  # noqa: E402
from rc3 import half_mass_radius, measure, pair_exclusions, scaling_delta  # noqa: E402

from omega.c0.references import rng_from_key  # noqa: E402
from omega.experiments.v11.gate import head_commit  # noqa: E402

MASTER_R3 = 20261021
OUT_FULL = ROOT / "results" / "r3_0b"
CAP = 2_000_000
N_MAX = 200_000
TARGETS = (10_000, 80_000)

# (id_str, id_num, generador, parametro, semillas)
FAMS: list[tuple[str, int, str, float | int, tuple[int, ...]]] = [
    ("J1_w1", 1, "J1", 1, (0,)), ("J1_w2", 2, "J1", 2, (0,)), ("J1_w3", 3, "J1", 3, (0,)), ("J1_w4", 4, "J1", 4, (0,)),
    ("J2_w2", 5, "J2", 2, (0,)), ("J2_w3", 6, "J2", 3, (0,)), ("J2_w4", 7, "J2", 4, (0,)),
    ("J3", 8, "J3", 0, (0, 1, 2)),
    ("J4_p0.02", 9, "J4", 0.02, (0, 1, 2)), ("J4_p0.1", 10, "J4", 0.1, (0, 1, 2)), ("J4_p0.5", 11, "J4", 0.5, (0, 1, 2)),
]
FAMILY_ID = {f[0]: f[1] for f in FAMS}

# ------------------------------------------------------------------ posets (predecesores como bitmasks de python)


def _row_to_int(row: np.ndarray) -> int:
    return int.from_bytes(np.packbits(row.astype(np.uint8), bitorder="little").tobytes(), "little")


def transitive_closure(direct: list[int]) -> list[int]:
    """Clausura transitiva de predecesores directos (bitmasks); orden arbitrario (DAG): topologico por Kahn."""
    n = len(direct)
    indeg = [0] * n
    succ: list[list[int]] = [[] for _ in range(n)]
    for j, m in enumerate(direct):
        i = 0
        while m:
            low = m & -m
            i = low.bit_length() - 1
            succ[i].append(j)
            indeg[j] += 1
            m ^= low
    out = list(direct)
    dq = deque(j for j in range(n) if indeg[j] == 0)
    seen = 0
    while dq:
        i = dq.popleft()
        seen += 1
        full_i = out[i] | (1 << i)
        for j in succ[i]:
            out[j] |= full_i
            indeg[j] -= 1
            if indeg[j] == 0:
                dq.append(j)
    if seen != n:
        raise ValueError("la relacion tiene un ciclo")
    return out


def covers(closure: list[int]) -> list[tuple[int, ...]]:
    """Reduccion transitiva (predecesores inmediatos) desde la clausura; etiquetas topologicas (pred < elemento)."""
    out = []
    for m in closure:
        cov = []
        while m:
            i = m.bit_length() - 1
            cov.append(i)
            m &= ~(closure[i] | (1 << i))
        out.append(tuple(cov))
    return out


def masks_to_tuples(masks: list[int]) -> list[tuple[int, ...]]:
    return [tuple(i for i in range(m.bit_length()) if m >> i & 1) for m in masks]


def poset_j1(w: int, nc: int) -> list[tuple[int, ...]]:
    """w cadenas disjuntas de longitud nc; id = k*w+i (cadena i, nivel k); SOLO predecesor inmediato (cobertura)."""
    return [((k - 1) * w + i,) if k > 0 else () for k in range(nc) for i in range(w)]


def poset_j2(w: int, nc: int) -> list[int]:
    """w cadenas, x_{i,k} < x_{j,k+2} (i!=j) mas orden de cadena; id = k*w+i; CIERRE transitivo (bitmasks)."""
    direct = []
    for k in range(nc):
        for i in range(w):
            m = 0
            if k >= 1:
                m |= 1 << ((k - 1) * w + i)
            if k >= 2:
                for j in range(w):
                    if j != i:
                        m |= 1 << ((k - 2) * w + j)
            direct.append(m)
    return transitive_closure(direct)


def poset_j3(n: int, rng: np.random.Generator) -> list[int]:
    """n puntos uniformes en [0,1]^2; x<y sii ambas coordenadas menores (ya transitivo; bitmasks)."""
    pts = rng.random((n, 2))
    xs, ys = pts[:, 0], pts[:, 1]
    less = (xs[None, :] < xs[:, None]) & (ys[None, :] < ys[:, None])  # less[j,i]: i<j
    return [_row_to_int(less[j]) for j in range(n)]


def poset_j4(n: int, p: float, rng: np.random.Generator) -> list[int]:
    """Percolacion transitiva: para i<j, i<j con prob p; cierre transitivo (bitmasks). Fila j = rng.random(j)
    (prefijo-consistente). Cierre incremental: se recorre el directo de mayor a menor indice, saltando lo ya alcanzado."""
    closure = [0]
    for j in range(1, n):
        rem = _row_to_int(rng.random(j) < p)
        acc = 0
        while rem:
            i = rem.bit_length() - 1
            acc |= closure[i] | (1 << i)
            rem &= ~acc
        closure.append(acc)
    return closure


def make_poset(gen: str, par: float | int, n_param: int, key: tuple[int, ...]) -> list[tuple[int, ...]]:
    """Predecesores (cobertura, o todos para J3) como tuplas, listos para enumerate_downsets."""
    if gen == "J1":
        return poset_j1(int(par), n_param)
    if gen == "J2":
        return covers(poset_j2(int(par), n_param))
    if gen == "J3":
        return masks_to_tuples(poset_j3(n_param, rng_from_key(key)))
    if gen == "J4":
        return covers(poset_j4(n_param, float(par), rng_from_key(key)))
    raise ValueError(gen)


# ------------------------------------------------------------------ enumeracion de down-sets


def enumerate_downsets(preds: list[tuple[int, ...]], limit: int, edges: bool = True) -> dict[str, Any]:
    """BFS de todos los down-sets desde el vacio (preds[x] = tupla de predecesores de x; cobertura o clausura).
    Aborta (complete=False) en cuanto |L| > limit. Devuelve count (= limit+1 si aborta), complete y, si edges,
    eu, ev: aristas D -> D U {x}, una por (D,x)."""
    n = len(preds)
    succ: list[list[int]] = [[] for _ in range(n)]
    for j, ps in enumerate(preds):
        for i in ps:
            succ[i].append(j)
    idx: dict[int, int] = {0: 0}
    queue: deque[tuple[int, int, list[int]]] = deque([(0, 0, [x for x in range(n) if not preds[x]])])
    eu: list[int] = []
    ev: list[int] = []
    while queue:
        d, i, add = queue.popleft()
        for x in add:
            d2 = d | (1 << x)
            j = idx.get(d2)
            if j is None:
                j = len(idx)
                idx[d2] = j
                if j + 1 > limit:
                    return {"count": j + 1, "complete": False, "eu": None, "ev": None}
                add2 = [y for y in add if y != x]
                for y in succ[x]:
                    if all(d2 >> q & 1 for q in preds[y]):
                        add2.append(y)
                queue.append((d2, j, add2))
            if edges:
                eu.append(i)
                ev.append(j)
    return {"count": len(idx), "complete": True,
            "eu": np.asarray(eu, dtype=np.int64) if edges else None, "ev": np.asarray(ev, dtype=np.int64) if edges else None}


def cover_graph(count: int, eu: np.ndarray, ev: np.ndarray) -> sparse.csr_array:
    """Adyacencia simetrica csr, data 1.0 (mismo formato que rc3/p1d3_panel._from_edges)."""
    m = sparse.coo_array((np.ones(2 * eu.size), (np.concatenate([eu, ev]), np.concatenate([ev, eu]))), shape=(count, count))
    a = sparse.csr_array(m)
    a.data[:] = 1.0
    return a


# ------------------------------------------------------------------ seleccion de tamano


def _elements(gen: str, par: float | int, n_param: int) -> int:
    return n_param * int(par) if gen in ("J1", "J2") else n_param


def choose_size(gen: str, par: float | int, key: tuple[int, ...], target: int) -> dict[str, Any]:
    """Menor n_param con |L| >= target: n=1,2,4,... hasta lograrlo (o N_MAX), luego biseccion. Poset regenerado por prueba."""
    def ge(n: int) -> bool:
        return not enumerate_downsets(make_poset(gen, par, n, key), target - 1, edges=False)["complete"]

    n, prev, trials = 1, 0, 0
    while True:
        trials += 1
        if ge(n):
            break
        prev = n
        if n >= N_MAX:
            return {"found": False, "trials": trials}
        n = min(2 * n, N_MAX)
    lo, hi = prev, n  # lo no cumple (0 = vacio), hi cumple
    while hi - lo > 1:
        mid = (lo + hi) // 2
        trials += 1
        if ge(mid):
            hi = mid
        else:
            lo = mid
    return {"found": True, "n_param": hi, "trials": trials}


def build_size(gen: str, par: float | int, seed: int, fid: int, si: int, cap: int) -> dict[str, Any]:
    """Elige tamano, enumera L(J) completo (con tope) y construye el grafo. Devuelve info (+ 'adj' si OK)."""
    key = (MASTER_R3, fid, seed, si)
    target = TARGETS[si]
    t0 = time.perf_counter()
    info: dict[str, Any] = {"key": list(key), "target": target}
    ch = choose_size(gen, par, key, target)
    info["search_trials"] = ch["trials"]
    if not ch["found"]:
        info.update(status="SEARCH_FAIL", n_param=None, n_elements=None, L=None)
        info["seconds_build"] = round(time.perf_counter() - t0, 2)
        return info
    n_param = ch["n_param"]
    info.update(n_param=n_param, n_elements=_elements(gen, par, n_param))
    en = enumerate_downsets(make_poset(gen, par, n_param, key), cap)
    if not en["complete"]:
        info.update(status="CAP_EXCEEDED", L=None, L_lower_bound=en["count"], cap=cap)
        info["seconds_build"] = round(time.perf_counter() - t0, 2)
        return info
    L = en["count"]
    info.update(status="OK", L=L, n_edges=int(en["eu"].size), out_of_band=bool(not (target <= L <= 1.5 * target)),
                L_over_target=L / target)
    info["adj"] = cover_graph(L, en["eu"], en["ev"])
    info["seconds_build"] = round(time.perf_counter() - t0, 2)
    return info


# ------------------------------------------------------------------ par


def run_pair(spec: dict[str, Any], cap: int = CAP) -> dict[str, Any]:
    t0 = time.perf_counter()
    fid = FAMILY_ID[spec["family"]]
    sizes = []
    for si in (0, 1):
        info = build_size(spec["gen"], spec["param"], spec["seed"], fid, si, cap)
        adj = info.pop("adj", None)
        if adj is not None:
            m = measure(adj, tuple(info["key"]))
            info.update(m)
            info["E_over_N"] = m["mean_degree"] / 2.0
            del adj
        sizes.append(info)
    a, b = sizes
    ok = all(s["status"] == "OK" for s in sizes)
    delta = None
    ex: list[str] = []
    w5 = None
    if ok:
        delta = scaling_delta(a["R"], b["R"], a["n_giant"] or 0, b["n_giant"] or 0)
        ex = pair_exclusions(delta, (a["annulus_status"], b["annulus_status"]), (a["coherence_status"], b["coherence_status"]))
        w5 = bool(abs(b["E_over_N"] - a["E_over_N"]) / a["E_over_N"] < 0.25)
    row = {k: v for k, v in spec.items()}
    kl = None
    if ok and a["n_elements"] > 1 and b["n_elements"] > 1:
        kl = [math.log(a["L"]) / math.log(a["n_elements"]), math.log(b["L"]) / math.log(b["n_elements"])]
    status = "OK" if ok else "+".join(sorted({s["status"] for s in sizes if s["status"] != "OK"}))
    row.update({"sizes": sizes, "status": status, "delta": delta, "exclusions": ex, "w5_valid": w5, "k_local": kl,
                "verdict": "SIN_MEDIDA" if not ok else (("EXCLUIDO(" + ",".join(ex) + ")") if ex else "NO-EXCLUIDO"),
                "excluded": bool(ex), "seconds": round(time.perf_counter() - t0, 1)})
    return row


def specs() -> list[dict[str, Any]]:
    return [{"id": f"{name}_s{seed}", "family": name, "gen": gen, "param": par, "seed": seed}
            for name, _fid, gen, par, seeds in FAMS for seed in seeds]


_SPECS: list[dict[str, Any]] = []
_CAP = [CAP]


def _run_idx(i: int) -> dict[str, Any]:
    return run_pair(_SPECS[i], _CAP[0])


def fmt(r: dict[str, Any]) -> str:
    a, b = r["sizes"]
    d = "None" if r["delta"] is None else f"{r['delta']:.3f}"
    return (f"  {r['id']:<14} n={a['n_elements']}/{b['n_elements']} L={a['L']}/{b['L']} delta={d} "
            f"ann={a.get('annulus_status')}/{b.get('annulus_status')} coh={a.get('coherence_status')}/{b.get('coherence_status')} "
            f"{r['verdict']} W5={r['w5_valid']} status={r['status']} {r['seconds']}s")


# ------------------------------------------------------------------ prediccion y corroboracion (§5)


def _abstains(r: dict[str, Any]) -> bool:
    bad = {"NO_EVALUABLE", "SIN_VENTANA"}
    return r["delta"] is None or any(s.get("annulus_status") in bad or s.get("coherence_status") in bad for s in r["sizes"])


def predict(r: dict[str, Any]) -> dict[str, Any]:
    """Fila de prediccion congelada (§5) para un par: match True/False/None (None = incierto/no aplicable)."""
    fam, ex = r["family"], r["exclusions"]
    out: dict[str, Any] = {"id": r["id"], "family": fam, "prediction": None, "match": None, "note": ""}
    if r["status"] != "OK":
        out.update(prediction="(ver familia)", match=False, note=f"sin medida: {r['status']}")
        return out
    d = r["delta"]
    if fam == "J1_w1":
        out.update(prediction="X3", match="X3" in ex)
    elif fam in ("J1_w2", "J1_w3"):
        w = int(fam[-1])
        out.update(prediction=f"no excluido, delta en 1/{w} +- 0.08",
                   match=bool((not r["excluded"]) and d is not None and abs(d - 1.0 / w) <= 0.08))
    elif fam == "J1_w4":
        n_inc = sum(s.get("coherence_status") == "INCOHERENTE" for s in r["sizes"])
        out.update(prediction="X4-marginal", match=bool("X4" in ex or n_inc == 1))
    elif fam.startswith("J2_"):
        out.update(prediction="X3", match="X3" in ex)
    elif fam == "J3":
        k = r["k_local"]
        explosion = bool(k is not None and k[1] > 1.25 * k[0])
        out.update(prediction="X1 o abstencion con explosion",
                   match=bool("X1" in ex or ((not r["excluded"]) and _abstains(r) and explosion)),
                   note=f"k_local={k} explosion={explosion} abstains={_abstains(r)}")
    elif fam in ("J4_p0.5", "J4_p0.1"):
        out.update(prediction="X3", match="X3" in ex)
    elif fam == "J4_p0.02":
        out.update(prediction="incierto", match=None,
                   note="RESIDUO a inspeccionar (W5=%s)" % r["w5_valid"] if not r["excluded"] else "excluido")
    return out


def corroboration(rows: list[dict[str, Any]]) -> dict[str, Any]:
    by = {r["id"]: r for r in rows}

    def one(i: str) -> dict[str, Any] | None:
        return by.get(i)

    ci = []
    for w in (2, 3):
        r = one(f"J1_w{w}_s0")
        ok = bool(r and r["status"] == "OK" and (not r["excluded"]) and r["delta"] is not None and abs(r["delta"] - 1.0 / w) <= 0.08)
        ci.append({"id": f"J1_w{w}_s0", "ok": ok, "delta": None if r is None else r["delta"],
                   "exclusions": None if r is None else r["exclusions"]})
    cii = []
    for w in (2, 3, 4):
        r = one(f"J2_w{w}_s0")
        ok = bool(r and r["status"] == "OK" and "X3" in r["exclusions"])
        cii.append({"id": f"J2_w{w}_s0", "ok": ok, "delta": None if r is None else r["delta"],
                    "exclusions": None if r is None else r["exclusions"]})
    viol = []
    for r in rows:
        if r["gen"] in ("J3", "J4") and r["status"] == "OK" and r["w5_valid"] and (not r["excluded"]) \
                and r["delta"] is not None and 0.2 <= r["delta"] <= 0.6:
            viol.append({"id": r["id"], "delta": r["delta"]})
    c = {"i": all(x["ok"] for x in ci), "ii": all(x["ok"] for x in cii), "iii": not viol}
    return {"conditions": {"i_J1(2),J1(3)": ci, "ii_J2_X3": cii, "iii_residuos": viol}, "i": c["i"], "ii": c["ii"], "iii": c["iii"],
            "R3_T1_T2_corroborados": all(c.values()), "residuo": bool(viol)}


def table(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for r in rows:
        a, b = r["sizes"]
        out.append({"id": r["id"], "status": r["status"], "n_elements": [a["n_elements"], b["n_elements"]],
                    "L": [a["L"], b["L"]], "out_of_band": [a.get("out_of_band"), b.get("out_of_band")],
                    "mean_degree": [a.get("mean_degree"), b.get("mean_degree")], "R": [a.get("R"), b.get("R")],
                    "delta": r["delta"], "annulus": [a.get("annulus_status"), b.get("annulus_status")],
                    "coherence": [a.get("coherence_status"), b.get("coherence_status")], "exclusions": r["exclusions"],
                    "w5_valid": r["w5_valid"], "k_local": r["k_local"], "size_status": [a["status"], b["status"]],
                    "seconds": r["seconds"]})
    return out


def _load(path: Path) -> dict[str, dict[str, Any]]:
    done: dict[str, dict[str, Any]] = {}
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                try:
                    r = json.loads(line)
                    done[r["id"]] = r
                except json.JSONDecodeError:
                    pass
    return done


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--fresh", action="store_true", help="ignora/borra pairs.jsonl previo")
    ap.add_argument("--cap", type=int, default=CAP, help="tope duro de down-sets")
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args(argv)
    out = Path(args.out) if args.out else OUT_FULL
    out.mkdir(parents=True, exist_ok=True)
    pj = out / "pairs.jsonl"
    if args.fresh and pj.exists():
        pj.unlink()
    done = _load(pj)
    _CAP[0] = args.cap
    _SPECS[:] = specs()
    todo = [i for i, s in enumerate(_SPECS) if s["id"] not in done]
    order = {"J3": 0, "J4": 1, "J2": 2, "J1": 3}
    todo.sort(key=lambda i: (order[_SPECS[i]["gen"]], -float(_SPECS[i]["param"])))
    t0 = time.perf_counter()
    print(f"{len(_SPECS)} pares ({len(done)} ya hechos), procs={args.procs}, cap={args.cap}", flush=True)
    with mp.get_context("fork").Pool(args.procs) as pool, open(pj, "a") as fh:
        for row in pool.imap_unordered(_run_idx, todo, chunksize=1):
            fh.write(json.dumps(row, default=str) + "\n")
            fh.flush()
            done[row["id"]] = row
            print(fmt(row), flush=True)
    rows = [done[s["id"]] for s in _SPECS]
    _atomic_write(pj, "".join(json.dumps(r, default=str) + "\n" for r in rows))
    preds = [predict(r) for r in rows]
    corr = corroboration(rows)
    summary = {"step": "r3_0b", "code_commit": head_commit(), "master": MASTER_R3, "cap": args.cap, "n_pairs": len(rows),
               "seconds_this_run": round(time.perf_counter() - t0, 1),
               "pair_seconds_total": round(sum(r["seconds"] for r in rows), 1),
               "corroboration": corr, "predictions": preds, "by_pair": table(rows),
               "failures": {r["id"]: r["status"] for r in rows if r["status"] != "OK"}}
    _atomic_write(out / "summary.json", json.dumps(summary, indent=2, default=str))
    for p in preds:
        print(f"  pred {p['id']:<14} {p['prediction']!s:<40} match={p['match']} {p['note']}")
    print(f"corroboracion: i={corr['i']} ii={corr['ii']} iii={corr['iii']} -> {corr['R3_T1_T2_corroborados']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
