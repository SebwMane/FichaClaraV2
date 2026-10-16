"""L-CIC-0b + RC-1 (docs/OMEGA_CIC_L0.md §3-§6): observable de anillos y validacion fuera de muestra de la regla combinada.

Salida: results/cic_l0b/{graphs.jsonl, summary.json}. Uso:
  python tools/cic_l0b.py [--procs 3] [--smoke]
Grupos: V (geometrias d >= 2), X (adversarios), U (1D), E (informativos). Familias marcadas `new` = las estrella de §5.

Decisiones de lectura (ver informe): instancias deterministas con una sola semilla (0), aleatorias con semillas 0-2;
localidad (i) sobre 2000 aristas muestreadas (sin reemplazo, key + (2,)) del grafo completo; curvatura sobre 300 aristas
(key + (4,), igual que tools/k_audit.py); Nivel II = estado de D_B2 con umbral 0.05; "no geometricos" = X + U + retazos;
Heisenberg no entra en ningun veredicto; las familias nuevas son solo las marcadas con estrella en §5.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from collections import defaultdict  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402
from scipy import sparse  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from c0_dynamics import _atomic_write  # noqa: E402
from coh_l0b import cycle_times_rr, heisenberg, husimi_cactus, prufer_tree, subdivided_binary_tree  # noqa: E402
from p1d3_panel import (  # noqa: E402
    _codes,
    _from_edges,
    caveman,
    erdos_renyi,
    patchwork,
    random_regular,
    rgg,
    rgg_gradient,
    rgg_sphere,
    ring,
    torus_lattice,
    watts_strogatz,
)

from omega.c0.references import rng_from_key  # noqa: E402
from omega.curvature.ollivier_sparse import ollivier_edge_sparse  # noqa: E402
from omega.diagnostics.annulus import annulus_profile, annulus_status  # noqa: E402
from omega.diagnostics.coherence import coherence_status, edge_coherence_profile  # noqa: E402
from omega.diagnostics.sampled_growth import sampled_ball_profile  # noqa: E402
from omega.experiments.v11.gate import head_commit  # noqa: E402

MASTER = 20261013
NS = 20_000
SEEDS = (0, 1, 2)
OUT_FULL = ROOT / "results" / "cic_l0b"
OUT_SMOKE = ROOT / "results" / "cic_l0b_smoke"

# RC-1 (§4), umbrales congelados
LOCALITY_MIN = 0.5
KAPPA_LO, KAPPA_HI, FNEG_MAX, FNEG_CUT = -0.0566, 0.2700, 0.4283, -0.1
N_LOC_EDGES = 2000
N_KAPPA_EDGES = 300
IDLENESS = 0.5

FAMILY_ID: dict[str, int] = {
    "RGG2_k12": 1, "RGG3_k12": 2, "RGG3_k8": 3, "RGG3_k20": 4, "RGG_S2_k12": 5, "RGG3_caja_k12": 6, "RGG2_caja_k12": 7,
    "cuadrado_141": 8, "T3_27": 9,
    "arbol_prufer": 20, "arbol_binario_sub8": 21, "cactus_triangulos": 22, "arbol_cliques_K6": 23, "2arbol": 24,
    "grafo_triangulos_k12": 25, "ER_k12": 26, "RR_k12": 27, "cliques_K12_enlazadas": 28,
    "anillo_k12": 40, "caveman_K8": 41, "C100xRR200_k4": 42, "WS_b0.001": 43,
    "Heisenberg_Z27": 60, "retazos_b2": 61,
}
NEW_FAMILIES = frozenset({"RGG3_k8", "RGG3_k20", "arbol_cliques_K6", "2arbol", "grafo_triangulos_k12", "cliques_K12_enlazadas"})
DIM = {"RGG2_k12": 2, "RGG2_caja_k12": 2, "RGG_S2_k12": 2, "cuadrado_141": 2,
       "RGG3_k12": 3, "RGG3_k8": 3, "RGG3_k20": 3, "RGG3_caja_k12": 3, "T3_27": 3}
DENSITY = ("RGG3_k8", "RGG3_k12", "RGG3_k20")
SMOKE_FAMILIES = ("RGG3_k12", "cuadrado_141", "arbol_binario_sub8", "arbol_cliques_K6", "2arbol", "cliques_K12_enlazadas",
                  "anillo_k12", "grafo_triangulos_k12")

# ------------------------------------------------------------------ generadores dispersos nuevos


def clique_tree(n_target: int, m: int, rng: np.random.Generator) -> sparse.csr_array:
    """Arbol de cliques K_m: cada clique nueva comparte exactamente un nodo (miembro uniforme) con una clique uniforme."""
    nc = 1 + max(n_target - m, 0) // (m - 1)
    members = np.empty((nc, m), dtype=np.int64)
    members[0] = np.arange(m)
    pick = (rng.random(nc) * np.arange(nc)).astype(np.int64)  # pick[j] uniforme en 0..j-1 (j >= 1)
    col = rng.integers(m, size=nc)
    nxt = m
    for j in range(1, nc):
        members[j, 0] = members[pick[j], col[j]]
        members[j, 1:] = np.arange(nxt, nxt + m - 1)
        nxt += m - 1
    ia, ib = np.triu_indices(m, 1)
    return _from_edges(nxt, members[:, ia].ravel(), members[:, ib].ravel())


def random_two_tree(n: int, rng: np.random.Generator) -> sparse.csr_array:
    """2-arbol aleatorio: triangulo inicial; cada nodo nuevo se une a los dos extremos de una arista uniforme."""
    eu = np.empty(3 + 2 * (n - 3), dtype=np.int64)
    ev = np.empty_like(eu)
    eu[:3], ev[:3] = (0, 0, 1), (1, 2, 2)
    m = 3
    t = rng.random(n)
    for x in range(3, n):
        e = int(t[x] * m)
        a, b = int(eu[e]), int(ev[e])
        eu[m], ev[m] = a, x
        eu[m + 1], ev[m + 1] = b, x
        m += 2
    return _from_edges(n, eu[:m], ev[:m])


def triangles_plus_random(n: int, mean_deg: float, rng: np.random.Generator) -> sparse.csr_array:
    """N/3 triangulos sobre ternas uniformes (distintas) mas aristas aleatorias nuevas hasta grado medio `mean_deg`."""
    t = n // 3
    tri = rng.integers(n, size=(t, 3))
    bad = (tri[:, 0] == tri[:, 1]) | (tri[:, 0] == tri[:, 2]) | (tri[:, 1] == tri[:, 2])
    while bad.any():
        tri[bad] = rng.integers(n, size=(int(bad.sum()), 3))
        bad = (tri[:, 0] == tri[:, 1]) | (tri[:, 0] == tri[:, 2]) | (tri[:, 1] == tri[:, 2])
    tu = np.concatenate([tri[:, 0], tri[:, 0], tri[:, 1]])
    tv = np.concatenate([tri[:, 1], tri[:, 2], tri[:, 2]])
    tcodes = np.unique(_codes(n, tu, tv))
    m_target = int(round(mean_deg * n / 2.0))
    extra = np.empty(0, dtype=np.int64)
    need = m_target - tcodes.size
    while extra.size < need:
        a = rng.integers(n, size=2 * (need - extra.size) + 16)
        b = rng.integers(n, size=a.size)
        keep = a != b
        c = np.unique(_codes(n, a[keep], b[keep]))
        c = c[~np.isin(c, tcodes)]
        extra = np.unique(np.concatenate([extra, c]))
    extra = rng.permutation(extra)[: max(need, 0)]
    allc = np.concatenate([tcodes, extra])
    return _from_edges(n, allc // n, allc % n)


def linked_cliques(n: int, size: int, links: int, rng: np.random.Generator) -> sparse.csr_array:
    """Bloques K_size disjuntos; cada bloque emite `links` aristas a bloques uniformes distintos (miembros uniformes)."""
    nb = n // size
    ia, ib = np.triu_indices(size, 1)
    base = (np.arange(nb, dtype=np.int64) * size)[:, None]
    us = [(base + ia[None, :]).ravel()]
    vs = [(base + ib[None, :]).ravel()]
    src = np.repeat(np.arange(nb, dtype=np.int64), links)
    dst = (src + 1 + rng.integers(nb - 1, size=src.size)) % nb  # uniforme sobre los demas bloques
    us.append(src * size + rng.integers(size, size=src.size))
    vs.append(dst * size + rng.integers(size, size=src.size))
    return _from_edges(nb * size, np.concatenate(us), np.concatenate(vs))


# ------------------------------------------------------------------ panel


def build(spec: dict[str, Any], rng: np.random.Generator) -> sparse.csr_array:
    g, a, n = spec["gen"], spec.get("args", {}), int(spec["N"])
    if g == "rgg":
        return rgg(n, a["d"], a["k"], rng, a.get("periodic", True))
    if g == "sphere":
        return rgg_sphere(n, 12.0, rng)
    if g == "gradient":
        return rgg_gradient(n, 12.0, rng)
    if g == "ring":
        return ring(n, 12)
    if g == "lattice":
        return torus_lattice(a["side"], a["dim"])
    if g == "caveman":
        return caveman(n, a["size"])
    if g == "patch":
        return patchwork(n, a["blocks"], rng)
    if g == "ws":
        return watts_strogatz(n, 12, a["beta"], rng)
    if g == "er":
        return erdos_renyi(n, 6 * n, rng)
    if g == "rr":
        return random_regular(n, 12, rng)
    if g == "prufer":
        return prufer_tree(n, rng)
    if g == "subtree":
        return subdivided_binary_tree(a["levels"], a["seg"])
    if g == "cactus":
        return husimi_cactus(n)
    if g == "cyc_rr":
        return cycle_times_rr(a["nc"], a["nr"], 4, rng)
    if g == "heis":
        return heisenberg(a["hn"])
    if g == "ctree":
        return clique_tree(n, 6, rng)
    if g == "2tree":
        return random_two_tree(n, rng)
    if g == "tri_rand":
        return triangles_plus_random(n, 12.0, rng)
    if g == "linked":
        return linked_cliques(n, 12, 2, rng)
    raise ValueError(g)


def specs(smoke: bool = False) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    sd = (0,) if smoke else SEEDS
    ns = 2000 if smoke else NS
    det = (0,)

    def add(fam: str, grp: str, gen: str, n: int, seeds: tuple[int, ...] = sd, **args: Any) -> None:
        if smoke and fam not in SMOKE_FAMILIES:
            return
        for s in seeds:
            out.append({"id": f"{fam}_s{s}", "family": fam, "group": grp, "new": fam in NEW_FAMILIES, "gen": gen, "N": n,
                        "args": args, "key": (MASTER, FAMILY_ID[fam], s)})

    # V
    add("RGG2_k12", "V", "rgg", ns, d=2, k=12.0)
    add("RGG3_k12", "V", "rgg", ns, d=3, k=12.0)
    add("RGG3_k8", "V", "rgg", ns, d=3, k=8.0)
    add("RGG3_k20", "V", "rgg", ns, d=3, k=20.0)
    add("RGG_S2_k12", "V", "sphere", ns)
    add("RGG3_caja_k12", "V", "rgg", ns, d=3, k=12.0, periodic=False)
    add("RGG2_caja_k12", "V", "rgg", ns, d=2, k=12.0, periodic=False)
    add("cuadrado_141", "V", "lattice", 45**2 if smoke else 141**2, det, side=45 if smoke else 141, dim=2)
    add("T3_27", "V", "lattice", 13**3 if smoke else 27**3, det, side=13 if smoke else 27, dim=3)
    # X
    add("arbol_prufer", "X", "prufer", ns)
    add("arbol_binario_sub8", "X", "subtree", 0, det, levels=7 if smoke else 11, seg=8)
    add("cactus_triangulos", "X", "cactus", ns, det)
    add("arbol_cliques_K6", "X", "ctree", ns)
    add("2arbol", "X", "2tree", ns)
    add("grafo_triangulos_k12", "X", "tri_rand", ns)
    add("ER_k12", "X", "er", ns)
    add("RR_k12", "X", "rr", ns)
    add("cliques_K12_enlazadas", "X", "linked", ns)
    # U
    add("anillo_k12", "U", "ring", ns, det)
    add("caveman_K8", "U", "caveman", ns, det, size=8)
    add("C100xRR200_k4", "U", "cyc_rr", 0, nc=20 if smoke else 100, nr=100 if smoke else 200)
    add("WS_b0.001", "U", "ws", ns, beta=0.001)
    # E
    add("Heisenberg_Z27", "E", "heis", 12**3 if smoke else 27**3, det, hn=12 if smoke else 27)
    add("retazos_b2", "E", "patch", ns, blocks=2)
    return out


# ------------------------------------------------------------------ medidas por grafo


def locality_fraction(adj: sparse.csr_array, rng: np.random.Generator, n_edges: int = N_LOC_EDGES) -> float:
    """Fraccion de aristas (muestra de `n_edges` sin reemplazo, todas si hay menos) en algun triangulo o 4-ciclo."""
    a = sparse.csr_array(adj)
    tri = sparse.triu(a, k=1, format="coo")
    eu, ev = np.asarray(tri.row, dtype=np.int64), np.asarray(tri.col, dtype=np.int64)
    if eu.size == 0:
        return 0.0
    if eu.size > n_edges:
        sel = np.sort(rng.choice(eu.size, size=n_edges, replace=False))
        eu, ev = eu[sel], ev[sel]
    hit = 0
    ind, ptr = a.indices, a.indptr
    for u, v in zip(eu.tolist(), ev.tolist()):
        nu = ind[ptr[u] : ptr[u + 1]]
        nv = ind[ptr[v] : ptr[v + 1]]
        if np.intersect1d(nu, nv, assume_unique=True).size:
            hit += 1
            continue
        au, bv = nu[nu != v], nv[nv != u]
        if au.size and bv.size and a[au][:, bv].nnz:  # a ~ b con a != b (sin diagonal)
            hit += 1
    return hit / float(eu.size)


def curvature_summary(adj: sparse.csr_array, rng: np.random.Generator, n_edges: int = N_KAPPA_EDGES) -> dict[str, Any]:
    tri = sparse.triu(adj, k=1, format="coo")
    xs, ys = np.asarray(tri.row, dtype=np.int64), np.asarray(tri.col, dtype=np.int64)
    if xs.size > n_edges:
        pick = np.sort(rng.choice(xs.size, size=n_edges, replace=False))
        xs, ys = xs[pick], ys[pick]
    k = np.asarray([ollivier_edge_sparse(adj, int(x), int(y), IDLENESS) for x, y in zip(xs, ys, strict=True)], dtype=np.float64)
    return {"n_edges": int(k.size), "median": float(np.median(k)), "mean": float(k.mean()), "f_neg": float(np.mean(k < FNEG_CUT))}


_SPECS: list[dict[str, Any]] = []


def run_spec(spec: dict[str, Any]) -> dict[str, Any]:
    t0 = time.perf_counter()
    key = tuple(int(k) for k in spec["key"])
    adj = build(spec, rng_from_key(key))
    t1 = time.perf_counter()
    ann = annulus_profile(adj, rng_from_key(key + (1,)))
    t2 = time.perf_counter()
    loc = locality_fraction(adj, rng_from_key(key + (2,)))
    coh = edge_coherence_profile(adj, rng_from_key(key + (3,)))
    t3 = time.perf_counter()
    kap = curvature_summary(adj, rng_from_key(key + (4,)))
    t4 = time.perf_counter()
    prof = sampled_ball_profile(adj, rng_from_key(key + (5,)))
    t5 = time.perf_counter()
    row = {k: v for k, v in spec.items() if k not in ("gen", "args", "key")}
    row["key"] = list(key)
    row["n"] = int(adj.shape[0])
    row["mean_degree"] = float(adj.nnz / adj.shape[0])
    row["annulus"] = ann
    row["annulus_status"] = annulus_status(ann)
    row["locality"] = loc
    row["coherence"] = coh
    row["coherence_status"] = coherence_status(coh)
    row["curvature"] = kap
    l2 = prof["level2"]["0.05"]
    row["level2"] = {"status": l2["status"], "D_conv": l2["D_conv"], "max_window": prof["max_window"]}
    row["rc1"] = rc1_conditions(row)
    row["seconds"] = {"gen": round(t1 - t0, 2), "annulus": round(t2 - t1, 2), "loc_coh": round(t3 - t2, 2),
                      "kappa": round(t4 - t3, 2), "level2": round(t5 - t4, 2)}
    return row


def rc1_conditions(row: dict[str, Any]) -> dict[str, bool]:
    k = row["curvature"]
    c = {
        "i_locality": bool(row["locality"] >= LOCALITY_MIN),
        "ii_coherence": bool(row["coherence_status"] == "COHERENTE"),
        "iii_curvature": bool(KAPPA_LO <= k["median"] <= KAPPA_HI and k["f_neg"] <= FNEG_MAX),
        "iv_level2": bool(row["level2"]["status"] == "CONVERGE"),
    }
    c["geometric"] = all(c.values())
    return c


def _run_idx(i: int) -> dict[str, Any]:
    row = run_spec(_SPECS[i])
    c = row["rc1"]
    print(f"  {row['id']:<28} n={row['n']:<6} deg={row['mean_degree']:.1f} ann={row['annulus_status']:<12} "
          f"loc={row['locality']:.2f} coh={row['coherence_status']:<11} kmed={row['curvature']['median']:.2f} "
          f"fneg={row['curvature']['f_neg']:.2f} II={row['level2']['status']:<11} RC1={int(c['geometric'])} {row['seconds']}",
          flush=True)
    return row


# ------------------------------------------------------------------ veredictos


def _verdict(sens: float | None, spec: float | None, valid: str, partial: str, invalid: str) -> str:
    if sens is None or spec is None:
        return "NO_EVALUABLE"
    if sens < 0.7 or spec < 0.7:
        return invalid
    if sens >= 0.9 and spec >= 0.9:
        return valid
    return partial


def _frac(rows: list[dict[str, Any]], pred: Any) -> float | None:
    return float(sum(bool(pred(r)) for r in rows) / len(rows)) if rows else None


def cic_verdict(rows: list[dict[str, Any]]) -> dict[str, Any]:
    v = [r for r in rows if r["group"] == "V"]
    x = [r for r in rows if r["group"] == "X"]
    sens = _frac(v, lambda r: r["annulus_status"] == "CONEXO")
    excl = _frac(x, lambda r: r["annulus_status"] != "CONEXO")
    return {"verdict": _verdict(sens, excl, "REDUNDANCIA-VALIDA", "REDUNDANCIA-PARCIAL", "REDUNDANCIA-INVALIDA"),
            "sensitivity": sens, "exclusion": excl, "n_V": len(v), "n_X": len(x)}


def rc1_verdict(rows: list[dict[str, Any]]) -> dict[str, Any]:
    v = [r for r in rows if r["group"] == "V"]
    ng = [r for r in rows if r["group"] in ("X", "U") or r["family"].startswith("retazos")]
    sens = _frac(v, lambda r: r["rc1"]["geometric"])
    spec = _frac(ng, lambda r: not r["rc1"]["geometric"])
    return {"verdict": _verdict(sens, spec, "RC1-VALIDA", "RC1-PARCIAL", "RC1-INVALIDA"),
            "sensitivity": sens, "specificity": spec, "n_V": len(v), "n_nongeometric": len(ng)}


def family_table(rows: list[dict[str, Any]]) -> dict[str, Any]:
    fam: dict[str, dict[str, Any]] = defaultdict(lambda: {
        "group": None, "new": None, "n": 0, "annulus_status": [], "f_med": [], "c_med": [], "f2_med": [], "r_w": [],
        "locality": [], "coherence_status": [], "kappa_median": [], "f_neg": [], "level2": [], "rc1_geometric": 0,
        "rc1_fail": defaultdict(int)})
    for r in rows:
        f = fam[r["family"]]
        a = r["annulus"]
        f["group"], f["new"] = r["group"], r["new"]
        f["n"] += 1
        f["annulus_status"].append(r["annulus_status"])
        f["f_med"].append(None if a["f_med"] is None else round(a["f_med"], 3))
        f["c_med"].append(a["c_med"])
        f["f2_med"].append(None if a["f2_med"] is None else round(a["f2_med"], 3))
        f["r_w"].append(a["r_w"])
        f["locality"].append(round(r["locality"], 3))
        f["coherence_status"].append(r["coherence_status"])
        f["kappa_median"].append(round(r["curvature"]["median"], 3))
        f["f_neg"].append(round(r["curvature"]["f_neg"], 3))
        f["level2"].append(r["level2"]["status"])
        f["rc1_geometric"] += int(r["rc1"]["geometric"])
        for c in ("i_locality", "ii_coherence", "iii_curvature", "iv_level2"):
            if not r["rc1"][c]:
                f["rc1_fail"][c] += 1
    out = {k: dict(v) for k, v in fam.items()}
    for v in out.values():
        v["rc1_fail"] = dict(v["rc1_fail"])
    return out


def same_state(rows: list[dict[str, Any]], fams: tuple[str, ...]) -> dict[str, Any]:
    sel = [r for r in rows if r["family"] in fams]
    states = sorted({r["annulus_status"] for r in sel})
    return {"families": [f for f in fams if any(r["family"] == f for r in sel)], "states": states,
            "same_state": bool(len(states) == 1) if sel else None}


def evaluate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    evaluated = [r for r in rows if r["family"] != "Heisenberg_Z27"]  # Heisenberg no entra en ningun veredicto
    d2 = tuple(f for f, d in DIM.items() if d == 2)
    d3 = tuple(f for f, d in DIM.items() if d == 3)
    out: dict[str, Any] = {"L_CIC_0b": {}, "RC1": {}}
    for name, sel in (("all", evaluated), ("new", [r for r in evaluated if r["new"]]),
                      ("repeated", [r for r in evaluated if not r["new"]])):
        out["L_CIC_0b"][name] = cic_verdict(sel)
        out["RC1"][name] = rc1_verdict(sel)
    out["L_CIC_0b"]["density_independence"] = same_state(rows, DENSITY)
    d2s, d3s = same_state(rows, d2), same_state(rows, d3)
    both = sorted(set(d2s["states"]) | set(d3s["states"]))
    out["L_CIC_0b"]["dimension_independence"] = {"d2": d2s, "d3": d3s, "same_state": bool(len(both) == 1) if both else None}
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--smoke", action="store_true", help="N~2000, 1 semilla, pocas familias -> results/cic_l0b_smoke/")
    args = ap.parse_args(argv)
    out = OUT_SMOKE if args.smoke else OUT_FULL
    _SPECS[:] = specs(args.smoke)
    order = sorted(range(len(_SPECS)), key=lambda i: -int(_SPECS[i]["N"]))
    t0 = time.perf_counter()
    print(f"{len(_SPECS)} grafos, procs={args.procs}", flush=True)
    done: dict[int, dict[str, Any]] = {}
    with mp.get_context("fork").Pool(args.procs) as pool:
        for i, row in zip(order, pool.imap(_run_idx, order, chunksize=1)):
            done[i] = row
    rows = [done[i] for i in range(len(_SPECS))]
    out.mkdir(parents=True, exist_ok=True)
    _atomic_write(out / "graphs.jsonl", "".join(json.dumps(r, default=str) + "\n" for r in rows))
    ev = evaluate(rows)
    summary: dict[str, Any] = {"step": "cic_l0b", "code_commit": head_commit(), "smoke": args.smoke, "n_graphs": len(rows),
                               "seconds": round(time.perf_counter() - t0, 1), "evaluation": ev,
                               "by_family": family_table(rows),
                               "heisenberg_excluded_from_verdicts": True}
    _atomic_write(out / "summary.json", json.dumps(summary, indent=2, default=str))
    for part in ("all", "new", "repeated"):
        c, r = ev["L_CIC_0b"][part], ev["RC1"][part]
        print(f"[{part}] L-CIC-0b {c['verdict']} sens={c['sensitivity']} excl={c['exclusion']} | "
              f"RC-1 {r['verdict']} sens={r['sensitivity']} spec={r['specificity']}")
    print(f"densidad={ev['L_CIC_0b']['density_independence']} dimension={ev['L_CIC_0b']['dimension_independence']['same_state']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
