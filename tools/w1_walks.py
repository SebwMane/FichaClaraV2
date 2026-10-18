"""OMEGA W-1: leyes de paseo (L1-L4), resolubilidad, estabilidad y universalidad (docs/OMEGA_W1_PRERREGISTRO.md, congelado).

Salida: results/w1/{cells.jsonl, summary.json} (smoke: results/w1_smoke). Uso:
  python tools/w1_walks.py [--procs 3] [--smoke] [--fresh]
Una celda = (geometria, tamano, semilla). Reanudable por celda (cells.jsonl se anade por celda terminada).

Decisiones de implementacion (ambiguedades del prerregistro -> lectura mas literal)
-----------------------------------------------------------------------------------
RNG: solo omega.c0.references.rng_from_key. Clave (MASTER_W1=20261023, geom_id, semilla, idx_tamano, k) con
  k=0 geometria, 1 origenes, 2 paseos (spec); k=3 bootstrap, k=4 perfil de bolas para R (anadidos; no previstos en el spec).
  idx_tamano 0 = N1, 1 = N2. Las geometrias deterministas usan semilla 0.
R: radio de media masa de RC-3 (rc3.half_mass_radius) sobre omega.diagnostics.sampled_growth.sampled_ball_profile
  (n_sources=400 por defecto, sobre la componente gigante) en vez de rc3.measure (mas pesado). OJO: ese perfil esta truncado
  en R_MAX=200 (r_cap); si la mediana de la bola no alcanza n/2 antes de r=200, R=200 y T_max=40000. Es la lectura literal
  de "R de RC-3" y se declara como anomalia en el informe (afecta a Z1, Z2, triangular, peines, cilindro...).
T: T_max=floor(R^2); Ts=(T_max//16, T_max//4, T_max). Si T_max//16<10 la celda es INSUFICIENTE (no se simulan paseos; a, SE = None).
Paseos: perezosos (quedarse con prob. 1/2; si no, vecino uniforme) sobre la componente gigante, vectorizados sobre CSR.
  Un unico uniforme u por paseante y paso: u<0.5 -> se queda; si no, indice de vecino = floor((u-0.5)*2*deg).
Origenes: M=256 uniformes sin reemplazo de la gigante. Por origen se simulan 3 paseos independientes; se reutilizan entre leyes:
  W1 para L1; (W1,W2) para L2 y L3; (W1,W2,W3) para L4. Las 3 T de una celda salen de los mismos prefijos de trayectoria.
L1: visitas de W1 al origen en pasos 1..t. L2: #t in 1..t con W1[t]==W2[t]. L3: |rango(W1) ∩ rango(W2)| hasta t (rangos incluyen
  el instante 0, luego el origen siempre esta en la interseccion). L4: |R1 ∩ R2 ∩ R3|.
a: pendiente OLS de ln(1+mean C(T)) frente a ln T en las 3 T. SE: bootstrap sobre origenes (500 remuestras, k=3); los mismos
  indices de remuestreo se usan para las 4 leyes de la celda. SE = desviacion tipica (ddof=1) de las 500 pendientes.
Clase: INSUFICIENTE si SE>0.03 o T_max//16<10 (primero); si no CRECE (a>=0.10), ACOTADA (a<=0.05), AMBIGUA.
Semillas aleatorias (RGG, Z3 diluido, retazos): 2 semillas (0,1). Los retazos tambien llevan 2 semillas (el spec no lo dice).
  Clase de geometria = clase comun si las semillas coinciden; si discrepan: INSUFICIENTE si alguna lo es (conservador, evita
  declarar determinada una geometria dudosa), si no AMBIGUA (incluye CRECE vs ACOTADA, que se anota como discordancia de semillas).
  No se agrupan origenes de las dos semillas porque T_max (y por tanto las T) difiere entre semillas. Se reportan ambas clases por semilla.
Resolubilidad: d* se explora en la rejilla {1, 1.5, ..., 6} (cortes dentro del rango del panel); d̂_c = [min, max] de los d*
  validos (intervalo cerrado). Condicion: d<d* -> CRECE, d>d* -> ACOTADA, d=d* cualquiera, y ninguna INSUFICIENTE con d!=d*.
  Sin d* valido: NO_RESUELTA si hay alguna INSUFICIENTE, si no NO_RESOLUBLE. Se usa la clase de geometria (arriba) de cada
  geometria homogenea; varias geometrias con el mismo d deben cumplir todas.
Estabilidad: RESOLUBLE tambien con las clases de N1 y d̂_c(N1) solapa (cerrado) con d̂_c(N2); si no, INESTABLE-N.
Universalidad (solo leyes RESOLUBLES y estables, capa 2 a N2): clase esperada de d_g = CRECE si d_g<min d̂_c, ACOTADA si
  d_g>max d̂_c, cualquiera si d_g in d̂_c. Discordancia determinada = observada CRECE vs esperada ACOTADA o viceversa -> NO UNIVERSAL.
  Observada AMBIGUA/INSUFICIENTE frente a esperada determinada no es discordancia determinada pero tampoco coincide: la ley
  queda INDETERMINADA (no cuenta como UNIVERSAL para W1-B). Retazos: sin d_g, solo se reportan.
Desenlace: leyes "validas" = RESOLUBLES y estables. W1-D si ninguna; W1-C si algun par con d̂_c disjuntos; W1-B si alguna es
  UNIVERSAL (los d̂_c son entonces mutuamente solapados = "comparten d̂_c"); W1-B- si ninguna es UNIVERSAL.
Tamanos (N real de la componente gigante en cells.jsonl). Se elige el parametro cuyo N nominal esta mas cerca de N1=1.25e5 y
  N2=1e6: Z1 N; Z2/triang lado 354/1000; Z3 50/100; Z4 19/32; Z5 10/16; Z6 7/10; FCC (4L^3) L=31/63; BCC (2L^3) L=40/79;
  RGG n=N; peine-1 (L vertices de ciclo, dientes con L vertices nuevos cada uno, N=L(L+1)) L=353/999; peine-2 (N=s^2(s+1)) s=100;
  cilindro C_n x C_10 con n=1e5; Z3 diluido lado 100 (N2); retazos n=1e6. Capa 2 solo N2.
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
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402
from scipy import sparse  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from a0_rc2 import cartesian, cycle, triangular_torus  # noqa: E402
from c0_dynamics import _atomic_write  # noqa: E402
from p1d3_panel import _from_edges, patchwork, rgg, torus_lattice  # noqa: E402
from rc3 import bcc_lattice, fcc_lattice, half_mass_radius  # noqa: E402

from omega.c0.references import rng_from_key  # noqa: E402
from omega.diagnostics.sampled_growth import giant_component, sampled_ball_profile  # noqa: E402
from omega.experiments.v11.gate import head_commit  # noqa: E402

MASTER_W1 = 20261023
OUT_FULL = ROOT / "results" / "w1"
OUT_SMOKE = ROOT / "results" / "w1_smoke"
M_ORIGINS = 256
N_BOOT = 500
LAWS = ("L1", "L2", "L3", "L4")
A_GROW, A_BOUND, SE_MAX, T_MIN = 0.10, 0.05, 0.03, 10
CRIT = {"L1": 2, "L2": 2, "L3": 4, "L4": 3}  # umbrales teoricos en Z^d (§1)
N_WALKERS = {"L1": 1, "L2": 2, "L3": 2, "L4": 3}

# ------------------------------------------------------------------ generadores nuevos


def comb(backbone: sparse.csr_array, tooth: int) -> sparse.csr_array:
    """Peine: a cada vertice v del grafo base se une un camino de `tooth` vertices nuevos (longitud tooth en aristas).
    N = n_b * (1 + tooth). Dientes: vertice (n_b + v*tooth + j), j=0..tooth-1; v - d_0 - d_1 - ..."""
    nb = backbone.shape[0]
    t = sparse.triu(backbone, k=1, format="coo")
    v = np.arange(nb, dtype=np.int64)[:, None]
    j = np.arange(tooth, dtype=np.int64)[None, :]
    node = nb + v * tooth + j
    prev = np.where(j == 0, v, node - 1)
    us = np.concatenate([np.asarray(t.row, dtype=np.int64), prev.ravel()])
    vs = np.concatenate([np.asarray(t.col, dtype=np.int64), node.ravel()])
    return _from_edges(nb * (1 + tooth), us, vs)


def comb1(L: int) -> sparse.csr_array:
    return comb(cycle(L), L)


def comb2(s: int) -> sparse.csr_array:
    return comb(torus_lattice(s, 2), s)


def diluted_z3(side: int, p: float, rng: np.random.Generator) -> sparse.csr_array:
    t = sparse.triu(torus_lattice(side, 3), k=1, format="coo")
    u, v = np.asarray(t.row, dtype=np.int64), np.asarray(t.col, dtype=np.int64)
    keep = rng.random(u.size) < p
    return _from_edges(side**3, u[keep], v[keep])


# (nombre, capa, d, d_g, aleatoria, gen, params N1, params N2, params N1 smoke, params N2 smoke)
def _g(name: str, layer: int, d: int | None, dg: int | None, rnd: bool, gen: str, p1: Any, p2: Any, s1: Any, s2: Any) -> dict[str, Any]:
    return {"name": name, "layer": layer, "d": d, "d_g": dg, "random": rnd, "gen": gen, "p": [p1, p2], "ps": [s1, s2]}


GEOMS: list[dict[str, Any]] = [
    _g("Z1", 1, 1, None, False, "torus1", 125_000, 1_000_000, 2000, 16_000),
    _g("Z2", 1, 2, None, False, "torus2", 354, 1000, 45, 126),
    _g("triangular", 1, 2, None, False, "tri", 354, 1000, 45, 126),
    _g("RGG2_k8", 1, 2, None, True, "rgg2", 125_000, 1_000_000, 2000, 16_000),
    _g("Z3", 1, 3, None, False, "torus3", 50, 100, 13, 25),
    _g("FCC", 1, 3, None, False, "fcc", 31, 63, 8, 16),
    _g("BCC", 1, 3, None, False, "bcc", 40, 79, 10, 20),
    _g("RGG3_k8", 1, 3, None, True, "rgg3", 125_000, 1_000_000, 2000, 16_000),
    _g("Z4", 1, 4, None, False, "torus4", 19, 32, 7, 11),
    _g("Z5", 1, 5, None, False, "torus5", 10, 16, 4, 6),
    _g("Z6", 1, 6, None, False, "torus6", 7, 10, 3, 5),
    _g("peine1", 2, None, 2, False, "comb1", None, 999, None, 126),
    _g("peine2", 2, None, 3, False, "comb2", None, 100, None, 25),
    _g("cilindro_CxC10", 2, None, 1, False, "cyl", None, 100_000, None, 1600),
    _g("Z3_diluido_p0.6", 2, None, 3, True, "dil", None, 100, None, 25),
    _g("retazos_b3", 2, None, None, True, "patch", None, 1_000_000, None, 16_000),
]
GID = {g["name"]: i + 1 for i, g in enumerate(GEOMS)}
BY_NAME = {g["name"]: g for g in GEOMS}


def build(gen: str, p: Any, rng: np.random.Generator) -> sparse.csr_array:
    if gen == "torus1":
        return torus_lattice(int(p), 1)
    if gen.startswith("torus"):
        return torus_lattice(int(p), int(gen[5:]))
    if gen == "tri":
        return triangular_torus(int(p))
    if gen == "rgg2":
        return rgg(int(p), 2, 8.0, rng, True)
    if gen == "rgg3":
        return rgg(int(p), 3, 8.0, rng, True)
    if gen == "fcc":
        return fcc_lattice(int(p))
    if gen == "bcc":
        return bcc_lattice(int(p))
    if gen == "comb1":
        return comb1(int(p))
    if gen == "comb2":
        return comb2(int(p))
    if gen == "cyl":
        return cartesian(cycle(int(p)), cycle(10))
    if gen == "dil":
        return diluted_z3(int(p), 0.6, rng)
    if gen == "patch":
        return patchwork(int(p), 3, rng)
    raise ValueError(gen)


def specs(smoke: bool) -> list[dict[str, Any]]:
    out = []
    for g in GEOMS:
        for si in (0, 1):
            par = g["ps" if smoke else "p"][si]
            if par is None:
                continue
            for seed in ((0, 1) if g["random"] else (0,)):
                out.append({"id": f"{g['name']}_n{si}_s{seed}", "geom": g["name"], "layer": g["layer"], "d": g["d"],
                            "d_g": g["d_g"], "size_idx": si, "seed": seed, "gen": g["gen"], "param": par})
    return out


# ------------------------------------------------------------------ paseos y leyes


def lazy_walks(indptr: np.ndarray, indices: np.ndarray, start: np.ndarray, T: int, rng: np.random.Generator,
               chunk: int = 512) -> np.ndarray:
    """Paseos perezosos vectorizados sobre CSR. Devuelve traj (T+1, K) int32; traj[0]=start."""
    K = start.size
    deg = np.diff(indptr).astype(np.int64)
    traj = np.empty((T + 1, K), dtype=np.int32)
    pos = start.astype(np.int32)
    traj[0] = pos
    t = 1
    while t <= T:
        c = min(chunk, T + 1 - t)
        u = rng.random((c, K))
        for i in range(c):
            ui = u[i]
            d = deg[pos]
            k = np.minimum(((ui - 0.5) * 2.0 * d).astype(np.int64), d - 1)
            nxt = indices[indptr[pos] + np.maximum(k, 0)]
            pos = np.where(ui < 0.5, pos, nxt).astype(np.int32)
            traj[t + i] = pos
        t += c
    return traj


def law_counts(traj: np.ndarray, Ts: tuple[int, ...]) -> dict[str, np.ndarray]:
    """traj (Tmax+1, 3, M): paseos 0,1,2 de cada origen. Devuelve C por ley (M, len(Ts))."""
    M = traj.shape[2]
    ti = np.asarray(Ts, dtype=np.int64)
    out: dict[str, np.ndarray] = {}
    c1 = np.cumsum(traj[1:, 0, :] == traj[0, 0, :][None, :], axis=0, dtype=np.int64)
    out["L1"] = c1[ti - 1].T.astype(np.float64)
    c2 = np.cumsum(traj[1:, 0, :] == traj[1:, 1, :], axis=0, dtype=np.int64)
    out["L2"] = c2[ti - 1].T.astype(np.float64)
    l3 = np.zeros((M, len(Ts)))
    l4 = np.zeros((M, len(Ts)))
    for m in range(M):
        cols = [np.ascontiguousarray(traj[:, w, m]) for w in range(3)]
        for j, t in enumerate(Ts):
            r1, r2, r3 = (np.unique(c[: t + 1]) for c in cols)
            i12 = np.intersect1d(r1, r2, assume_unique=True)
            l3[m, j] = i12.size
            l4[m, j] = np.intersect1d(i12, r3, assume_unique=True).size
    out["L3"], out["L4"] = l3, l4
    return out


def ols_slope(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Pendiente OLS por filas de y (..., n) frente a x (n,)."""
    xc = x - x.mean()
    return (y - y.mean(axis=-1, keepdims=True)) @ xc / float(xc @ xc)


def exponent(C: np.ndarray, Ts: tuple[int, ...], boot_idx: np.ndarray) -> tuple[float, float, list[float]]:
    x = np.log(np.asarray(Ts, dtype=np.float64))
    mean = C.mean(axis=0)
    a = float(ols_slope(x, np.log1p(mean)))
    bm = C[boot_idx].mean(axis=1)  # (B, nT)
    se = float(np.std(ols_slope(x, np.log1p(bm)), ddof=1))
    return a, se, [float(m) for m in mean]


def classify(a: float | None, se: float | None, t_max: int) -> str:
    if a is None or se is None or t_max // 16 < T_MIN or se > SE_MAX:
        return "INSUFICIENTE"
    if a >= A_GROW:
        return "CRECE"
    if a <= A_BOUND:
        return "ACOTADA"
    return "AMBIGUA"


def run_cell(spec: dict[str, Any], smoke: bool = False) -> dict[str, Any]:
    t0 = time.perf_counter()
    key = (MASTER_W1, GID[spec["geom"]], int(spec["seed"]), int(spec["size_idx"]))
    adj = build(spec["gen"], spec["param"], rng_from_key(key + (0,)))
    g = giant_component(adj)
    n_total, ng = int(adj.shape[0]), int(g.shape[0])
    del adj
    t1 = time.perf_counter()
    prof = sampled_ball_profile(g, rng_from_key(key + (4,)))
    R = half_mass_radius(prof)
    tb = time.perf_counter()
    t_max = int(math.floor(R * R))
    Ts = (t_max // 16, t_max // 4, t_max)
    row: dict[str, Any] = {k: spec[k] for k in ("id", "geom", "layer", "d", "d_g", "size_idx", "seed")}
    row.update({"key": list(key), "N_total": n_total, "N_giant": ng, "R": R, "r_cap": int(prof["r_cap"]), "T_max": t_max,
                "Ts": list(Ts), "mean_degree": float(g.nnz / ng)})
    laws: dict[str, Any] = {}
    if t_max // 16 < T_MIN:
        for L in LAWS:
            laws[L] = {"a": None, "se": None, "cls": "INSUFICIENTE", "means": None, "reason": "T_max//16<10"}
        t2 = t3 = time.perf_counter()
    else:
        orig = np.sort(rng_from_key(key + (1,)).choice(ng, size=min(M_ORIGINS, ng), replace=False)).astype(np.int64)
        M = orig.size
        start = np.tile(orig, 3)
        traj = lazy_walks(g.indptr, g.indices.astype(np.int32), start, t_max, rng_from_key(key + (2,)))
        t2 = time.perf_counter()
        C = law_counts(traj.reshape(t_max + 1, 3, M), Ts)
        del traj
        brng = rng_from_key(key + (3,))
        bidx = brng.integers(0, M, size=(N_BOOT, M))
        for L in LAWS:
            a, se, means = exponent(C[L], Ts, bidx)
            laws[L] = {"a": a, "se": se, "cls": classify(a, se, t_max), "means": means}
        t3 = time.perf_counter()
    row["laws"] = laws
    row["seconds"] = {"gen": round(t1 - t0, 1), "ball": round(tb - t1, 1), "walk": round(t2 - tb, 1),
                      "laws": round(t3 - t2, 1), "total": round(time.perf_counter() - t0, 1)}
    return row


# ------------------------------------------------------------------ analisis (§3)


def geometry_class(classes: list[str]) -> str:
    if len(set(classes)) == 1:
        return classes[0]
    if "INSUFICIENTE" in classes:
        return "INSUFICIENTE"
    return "AMBIGUA"


def dstar_grid(ds: list[float]) -> list[float]:
    lo, hi = int(min(ds)), int(max(ds))
    return [x / 2 for x in range(2 * lo, 2 * hi + 1)]


def resolve(items: list[tuple[float, str]]) -> dict[str, Any]:
    """items: (d, clase) de las geometrias homogeneas. Devuelve estado y d̂_c."""
    valid = []
    for ds in dstar_grid([d for d, _ in items]):
        ok = True
        for d, c in items:
            if d < ds and c != "CRECE":
                ok = False
            elif d > ds and c != "ACOTADA":
                ok = False
        if ok:
            valid.append(ds)
    if valid:
        return {"status": "RESOLUBLE", "dc": [min(valid), max(valid)]}
    return {"status": "NO_RESUELTA" if any(c == "INSUFICIENTE" for _, c in items) else "NO_RESOLUBLE", "dc": None}


def overlap(a: list[float], b: list[float]) -> bool:
    return max(a[0], b[0]) <= min(a[1], b[1])


def stability(res_n2: dict[str, Any], res_n1: dict[str, Any]) -> str:
    if res_n2["status"] != "RESOLUBLE":
        return "NA"
    if res_n1["status"] == "RESOLUBLE" and overlap(res_n1["dc"], res_n2["dc"]):
        return "ESTABLE"
    return "INESTABLE-N"


def expected_class(dg: float, dc: list[float]) -> str:
    if dc[0] <= dg <= dc[1]:
        return "CRITICO"
    return "CRECE" if dg < dc[0] else "ACOTADA"


def universality(dc: list[float], layer2: list[tuple[str, float | None, str]]) -> dict[str, Any]:
    """layer2: (nombre, d_g, clase observada). d_g None -> solo se reporta."""
    det: list[dict[str, Any]] = []
    status = "UNIVERSAL"
    for name, dg, obs in layer2:
        if dg is None:
            det.append({"geom": name, "d_g": None, "obs": obs, "expected": None, "verdict": "SOLO_REPORTE"})
            continue
        exp = expected_class(dg, dc)
        if exp == "CRITICO" or obs == exp:
            v = "OK"
        elif obs in ("CRECE", "ACOTADA"):
            v = "DISCORDANCIA"
            status = "NO_UNIVERSAL"
        else:
            v = "INDETERMINADA"
            if status == "UNIVERSAL":
                status = "INDETERMINADA"
        det.append({"geom": name, "d_g": dg, "obs": obs, "expected": exp, "verdict": v})
    return {"status": status, "detail": det}


def outcome(laws: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """laws[L] = {"res": .., "stab": .., "univ": ..}; devuelve desenlace y subclases."""
    valid = {L: v for L, v in laws.items() if v["res"]["status"] == "RESOLUBLE" and v["stab"] == "ESTABLE"}
    sub = [f"{L}: {v['res']['status']}/{v['stab']}/{v['univ']}" for L, v in laws.items()]
    if not valid:
        return {"outcome": "W1-D", "valid_laws": [], "subclasses": sub}
    names = sorted(valid)
    disjoint = [(a, b) for i, a in enumerate(names) for b in names[i + 1:]
                if not overlap(valid[a]["res"]["dc"], valid[b]["res"]["dc"])]
    if disjoint:
        o = "W1-C"
    elif any(v["univ"] == "UNIVERSAL" for v in valid.values()):
        o = "W1-B"
    else:
        o = "W1-B-"
    return {"outcome": o, "valid_laws": names, "disjoint_pairs": disjoint, "subclasses": sub}


def analyse(rows: list[dict[str, Any]]) -> dict[str, Any]:
    # clase por (geometria, tamano, ley) a partir de las semillas
    gc: dict[tuple[str, int], dict[str, Any]] = {}
    for r in rows:
        e = gc.setdefault((r["geom"], r["size_idx"]), {"seeds": {}})
        e["seeds"][r["seed"]] = {L: r["laws"][L]["cls"] for L in LAWS}
    anomalies: list[str] = []
    for (name, si), e in gc.items():
        e["cls"] = {L: geometry_class([s[L] for s in e["seeds"].values()]) for L in LAWS}
        for L in LAWS:
            cs = {s[L] for s in e["seeds"].values()}
            if len(cs) > 1:
                anomalies.append(f"semillas discrepan: {name} N{si + 1} {L}: {sorted(cs)}")
    laws: dict[str, Any] = {}
    for L in LAWS:
        l1 = [(g["d"], gc[(g["name"], si)]["cls"][L]) for g in GEOMS if g["layer"] == 1 for si in (1,) if (g["name"], si) in gc]
        l1b = [(g["d"], gc[(g["name"], 0)]["cls"][L]) for g in GEOMS if g["layer"] == 1 and (g["name"], 0) in gc]
        r2, r1 = (resolve(x) if x else {"status": "NA", "dc": None} for x in (l1, l1b))
        st = stability(r2, r1)
        lay2 = [(g["name"], g["d_g"], gc[(g["name"], 1)]["cls"][L]) for g in GEOMS if g["layer"] == 2 and (g["name"], 1) in gc]
        un = universality(r2["dc"], lay2) if st == "ESTABLE" else {"status": "NA", "detail": []}
        laws[L] = {"res": r2, "res_N1": r1, "stab": st, "univ": un["status"], "univ_detail": un["detail"]}
    out = outcome({L: {"res": v["res"], "stab": v["stab"], "univ": v["univ"]} for L, v in laws.items()})
    table5 = {}
    for L, v in laws.items():
        table5[L] = {
            "detecta_dimension": f"{v['res']['status']}, d_c={v['res']['dc']}",
            "robusta_heterogeneidad": v["univ"], "robusta_N": v["stab"],
            "caminantes": f"{N_WALKERS[L]} paseo(s); d_c={v['res']['dc']}",
            "red_vs_continuo": "No se mide (W-1 solo redes); nivel 3: R4 da 3 en continuo y 4 en red para k=2",
            "depende_de_la_ley": "comparar d_c entre leyes (ver 'dc_by_law')",
            "genera_dimension": "Solo medirla (Omega-1)"}
    out.update({"laws": laws, "table5": table5, "dc_by_law": {L: laws[L]["res"]["dc"] for L in LAWS}, "anomalies": anomalies})
    out["geometry_classes"] = {f"{n}_N{si + 1}": {"cls": e["cls"], "by_seed": e["seeds"]} for (n, si), e in gc.items()}
    out["predictions"] = predictions(out)
    return out


def predictions(s: dict[str, Any]) -> dict[str, Any]:
    gcl = s["geometry_classes"]
    laws = s["laws"]
    p: dict[str, Any] = {}
    for L in LAWS:
        crit = [gcl[f"{g['name']}_N2"]["cls"][L] for g in GEOMS if g["layer"] == 1 and g["d"] == CRIT[L] and f"{g['name']}_N2" in gcl]
        p[f"critico_AMBIGUO_{L}(d={CRIT[L]})"] = {"classes": crit, "match": bool(crit) and all(c == "AMBIGUA" for c in crit)}
    hi = {n: gcl[f"{n}_N2"]["cls"]["L1"] for n in ("Z5", "Z6") if f"{n}_N2" in gcl}
    p["Z5/Z6_INSUFICIENTE_N2_y_L3_NO_RESUELTA"] = {"Z5_Z6": hi, "L3": laws["L3"]["res"]["status"],
                                                   "match": any(c == "INSUFICIENTE" for c in hi.values()) and laws["L3"]["res"]["status"] == "NO_RESUELTA"}
    for L, d in (("L1", 2), ("L2", 2), ("L4", 3)):
        dc = laws[L]["res"]["dc"]
        p[f"{L}_RESOLUBLE_con_{d}_en_dc"] = {"res": laws[L]["res"]["status"], "dc": dc,
                                            "match": laws[L]["res"]["status"] == "RESOLUBLE" and dc is not None and dc[0] <= d <= dc[1]}
    c1 = {L: gcl.get("peine1_N2", {"cls": {}})["cls"].get(L) for L in LAWS}
    p["L2_NO_UNIVERSAL_y_peine1_ACOTADA"] = {"L2_univ": laws["L2"]["univ"], "peine1_L2": c1["L2"],
                                             "match": laws["L2"]["univ"] == "NO_UNIVERSAL" and c1["L2"] == "ACOTADA"}
    p["L1_UNIVERSAL_y_peine1_CRECE"] = {"L1_univ": laws["L1"]["univ"], "peine1_L1": c1["L1"],
                                        "match": laws["L1"]["univ"] == "UNIVERSAL" and c1["L1"] == "CRECE"}
    p["desenlace_W1-C(alt_W1-D)"] = {"outcome": s["outcome"], "match": s["outcome"] == "W1-C", "alt_match": s["outcome"] == "W1-D"}
    return p


# ------------------------------------------------------------------ ejecucion

_SPECS: list[dict[str, Any]] = []
_SMOKE = False


def _run_idx(i: int) -> dict[str, Any]:
    return run_cell(_SPECS[i], _SMOKE)


def fmt(r: dict[str, Any]) -> str:
    s = " ".join(f"{L}:{'-' if v['a'] is None else format(v['a'], '.3f')}±{'-' if v['se'] is None else format(v['se'], '.3f')}{v['cls'][:3]}"
                 for L, v in r["laws"].items())
    return f"  {r['id']:<28} N={r['N_giant']} R={r['R']:.1f} Tmax={r['T_max']} {s} {r['seconds']['total']}s"


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
    global _SMOKE
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--fresh", action="store_true")
    args = ap.parse_args(argv)
    _SMOKE = args.smoke
    out = OUT_SMOKE if args.smoke else OUT_FULL
    out.mkdir(parents=True, exist_ok=True)
    cj = out / "cells.jsonl"
    if args.fresh and cj.exists():
        cj.unlink()
    done = _load(cj)
    _SPECS[:] = specs(args.smoke)
    todo = [i for i, s in enumerate(_SPECS) if s["id"] not in done]
    todo.sort(key=lambda i: -(_SPECS[i]["size_idx"]))
    t0 = time.perf_counter()
    print(f"{len(_SPECS)} celdas ({len(done)} hechas), procs={args.procs}", flush=True)
    with mp.get_context("fork").Pool(args.procs) as pool, open(cj, "a") as fh:
        for row in pool.imap_unordered(_run_idx, todo, chunksize=1):
            fh.write(json.dumps(row, default=str) + "\n")
            fh.flush()
            done[row["id"]] = row
            print(fmt(row), flush=True)
    rows = [done[s["id"]] for s in _SPECS]
    _atomic_write(cj, "".join(json.dumps(r, default=str) + "\n" for r in rows))
    an = analyse(rows)
    summary = {"step": "w1", "code_commit": head_commit(), "smoke": args.smoke, "n_cells": len(rows),
               "seconds_this_run": round(time.perf_counter() - t0, 1), "cell_seconds_total": round(sum(r["seconds"]["total"] for r in rows), 1),
               "analysis": an,
               "cells": [{"id": r["id"], "geom": r["geom"], "size_idx": r["size_idx"], "seed": r["seed"], "N_giant": r["N_giant"],
                          "T_max": r["T_max"], "laws": {L: {k: r["laws"][L][k] for k in ("a", "se", "cls")} for L in LAWS}} for r in rows]}
    _atomic_write(out / "summary.json", json.dumps(summary, indent=2, default=str))
    print("OUTCOME", an["outcome"], {L: (v["res"]["status"], v["res"]["dc"], v["stab"], v["univ"]) for L, v in an["laws"].items()})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
