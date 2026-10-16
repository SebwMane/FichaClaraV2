"""L-COH-0b (docs/OMEGA_COH_L0.md §3): prueba de coherencia multiescala contra falsos positivos, sin dinamica.

Salida: results/coh_l0b/{graphs.jsonl, summary.json}. Uso:
  python tools/coh_l0b.py [--procs 3] [--smoke]
Grupos: V (geometrias), F (falsos positivos disenados), E (estatus especial; se informan, no deciden).
Decisiones de lectura (ver informe): arbol binario de 11 niveles (2^11-1 nodos), cada arista -> camino de 8 aristas;
cactus de Husimi de 19999 nodos (raiz + 9999 triangulos); las familias deterministas usan una sola instancia (semilla 0).
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
from p1d3_panel import (  # noqa: E402
    _from_edges,
    caveman,
    clique_lattice,
    erdos_renyi,
    patchwork,
    random_regular,
    rgg,
    rgg_gradient,
    rgg_sphere,
    ring,
    shortcuts,
    torus_lattice,
    watts_strogatz,
)

from omega.c0.references import rng_from_key  # noqa: E402
from omega.diagnostics.coherence import coherence_status, edge_coherence_profile  # noqa: E402
from omega.experiments.v11.gate import head_commit  # noqa: E402

MASTER = 20261012
NS = 20_000
SEEDS = (0, 1, 2)
OUT_FULL = ROOT / "results" / "coh_l0b"
OUT_SMOKE = ROOT / "results" / "coh_l0b_smoke"

FAMILY_ID: dict[str, int] = {
    "RGG2_k12": 1, "RGG3_k12": 2, "RGG_S2_k12": 3, "T3_27": 4, "cuadrado_141": 5, "RGG3_caja_k12": 6, "RGG2_caja_k12": 7,
    "RGG3_gradiente": 8,
    "arbol_prufer": 20, "arbol_binario_sub8": 21, "cactus_triangulos": 22, "ER_k12": 23, "RR_k12": 24, "retazos_b2": 25,
    "retazos_b4": 26, "retazos_b8": 27, "WS_b0.01": 28, "WS_b0.005": 29, "RGG3_atajos_1": 30, "C100xRR200_k4": 31,
    "anillo_k12": 40, "caveman_K8": 41, "cliques3D_K6_L15": 42, "Heisenberg_Z27": 43, "WS_b0.001": 44, "RGG3_atajos_0.1": 45,
}

# ------------------------------------------------------------------ generadores dispersos nuevos


def _edges_of(a: sparse.csr_array) -> tuple[np.ndarray, np.ndarray]:
    t = sparse.coo_array(sparse.triu(a, k=1))
    return np.asarray(t.row, dtype=np.int64), np.asarray(t.col, dtype=np.int64)


def prufer_tree(n: int, rng: np.random.Generator) -> sparse.csr_array:
    """Arbol uniforme sobre n nodos etiquetados, decodificando una secuencia de Pruefer (lineal)."""
    seq = rng.integers(n, size=n - 2)
    deg = np.ones(n, dtype=np.int64)
    np.add.at(deg, seq, 1)
    us: list[int] = []
    vs: list[int] = []
    ptr = int(np.flatnonzero(deg == 1)[0])
    leaf = ptr
    for x in seq.tolist():
        us.append(leaf)
        vs.append(x)
        deg[x] -= 1
        if deg[x] == 1 and x < ptr:
            leaf = x
        else:
            ptr += 1
            while deg[ptr] != 1:
                ptr += 1
            leaf = ptr
    us.append(leaf)
    vs.append(n - 1)
    return _from_edges(n, np.asarray(us, dtype=np.int64), np.asarray(vs, dtype=np.int64))


def subdivided_binary_tree(levels: int, seg: int) -> sparse.csr_array:
    """Arbol binario completo de `levels` niveles; cada arista sustituida por un camino de `seg` aristas."""
    nt = 2**levels - 1
    child = np.arange(1, nt, dtype=np.int64)
    parent = (child - 1) // 2
    ne = child.size
    n = nt + ne * (seg - 1)
    base = nt + np.arange(ne, dtype=np.int64) * (seg - 1)
    chain = [parent] + [base + j for j in range(seg - 1)] + [child]
    u = np.concatenate(chain[:-1])
    v = np.concatenate(chain[1:])
    return _from_edges(n, u, v)


def husimi_cactus(n_target: int) -> sparse.csr_array:
    """Arbol de triangulos: todo nodo esta en 2 triangulos salvo la frontera; crecimiento BFS hasta <= n_target nodos."""
    us: list[int] = []
    vs: list[int] = []
    queue: list[tuple[int, int]] = [(0, 2)]
    n = 1
    head = 0
    while head < len(queue):
        node, need = queue[head]
        head += 1
        for _ in range(need):
            if n + 2 > n_target:
                head = len(queue)
                break
            a, b = n, n + 1
            n += 2
            us += [node, node, a]
            vs += [a, b, b]
            queue.append((a, 1))
            queue.append((b, 1))
    return _from_edges(n, np.asarray(us, dtype=np.int64), np.asarray(vs, dtype=np.int64))


def cycle_times_rr(n_cycle: int, n_rr: int, k: int, rng: np.random.Generator) -> sparse.csr_array:
    """Producto cartesiano C_{n_cycle} x RR(n_rr, k); nodo (i, j) = i * n_rr + j."""
    ru, rv = _edges_of(random_regular(n_rr, k, rng))
    i = np.arange(n_cycle, dtype=np.int64)
    off = (i * n_rr)[:, None]
    u_rr = (off + ru[None, :]).ravel()
    v_rr = (off + rv[None, :]).ravel()
    j = np.arange(n_rr, dtype=np.int64)
    u_c = (i[:, None] * n_rr + j[None, :]).ravel()
    v_c = (((i + 1) % n_cycle)[:, None] * n_rr + j[None, :]).ravel()
    return _from_edges(n_cycle * n_rr, np.concatenate([u_rr, u_c]), np.concatenate([v_rr, v_c]))


def heisenberg(n: int) -> sparse.csr_array:
    """Cayley de H3(Z_n), generadores x=(1,0,0), y=(0,1,0) e inversos; (a,b,c)(a',b',c') = (a+a', b+b', c+c'+a b')."""
    a, b, c = (t.ravel().astype(np.int64) for t in np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing="ij"))
    idx = (a * n + b) * n + c
    # g*x = (a+1, b, c); g*y = (a, b+1, c+a)
    gx = (((a + 1) % n) * n + b) * n + c
    gy = (a * n + (b + 1) % n) * n + (c + a) % n
    return _from_edges(n**3, np.concatenate([idx, idx]), np.concatenate([gx, gy]))


# ------------------------------------------------------------------ panel


def build(spec: dict[str, Any], rng: np.random.Generator) -> sparse.csr_array:
    g, a, n = spec["gen"], spec.get("args", {}), int(spec["N"])
    if g == "rgg":
        return rgg(n, a["d"], 12.0, rng, a.get("periodic", True))
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
    if g == "cliques":
        return clique_lattice(a["side"], a["m"], a["dim"])
    if g == "patch":
        return patchwork(n, a["blocks"], rng)
    if g == "shortcuts":
        return shortcuts(n, a["frac"], rng)
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
    raise ValueError(g)


def specs(smoke: bool = False) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    sd = (0,) if smoke else SEEDS
    ns = 2000 if smoke else NS

    def add(fam: str, grp: str, gen: str, n: int, seeds: tuple[int, ...] = sd, **args: Any) -> None:
        for s in seeds:
            out.append({"id": f"{fam}_s{s}", "family": fam, "group": grp, "gen": gen, "N": n, "args": args,
                        "key": (MASTER, FAMILY_ID[fam], s)})

    det = (0,)
    # V
    add("RGG2_k12", "V", "rgg", ns, d=2)
    add("RGG3_k12", "V", "rgg", ns, d=3)
    add("RGG_S2_k12", "V", "sphere", ns)
    add("RGG3_caja_k12", "V", "rgg", ns, d=3, periodic=False)
    add("RGG2_caja_k12", "V", "rgg", ns, d=2, periodic=False)
    add("RGG3_gradiente", "V", "gradient", ns)
    if smoke:
        add("T3_27", "V", "lattice", 13**3, det, side=13, dim=3)
        add("cuadrado_141", "V", "lattice", 45**2, det, side=45, dim=2)
    else:
        add("T3_27", "V", "lattice", 27**3, det, side=27, dim=3)
        add("cuadrado_141", "V", "lattice", 141**2, det, side=141, dim=2)
    # F
    add("arbol_prufer", "F", "prufer", ns)
    add("arbol_binario_sub8", "F", "subtree", 0, det, levels=7 if smoke else 11, seg=8)
    add("cactus_triangulos", "F", "cactus", ns, det)
    add("ER_k12", "F", "er", ns)
    add("RR_k12", "F", "rr", ns)
    for b in (2, 4, 8):
        add(f"retazos_b{b}", "F", "patch", ns, blocks=b)
    add("WS_b0.01", "F", "ws", ns, beta=0.01)
    add("WS_b0.005", "F", "ws", ns, beta=0.005)
    add("RGG3_atajos_1", "F", "shortcuts", ns, frac=0.01)
    add("C100xRR200_k4", "F", "cyc_rr", 0, nc=20 if smoke else 100, nr=100 if smoke else 200)
    # E
    add("anillo_k12", "E", "ring", ns, det)
    add("caveman_K8", "E", "caveman", ns, det, size=8)
    if smoke:
        add("cliques3D_K6_L15", "E", "cliques", 7**3 * 6, det, side=7, m=6, dim=3)
        add("Heisenberg_Z27", "E", "heis", 12**3, det, hn=12)
    else:
        add("cliques3D_K6_L15", "E", "cliques", 15**3 * 6, det, side=15, m=6, dim=3)
        add("Heisenberg_Z27", "E", "heis", 27**3, det, hn=27)
    add("WS_b0.001", "E", "ws", ns, beta=0.001)
    add("RGG3_atajos_0.1", "E", "shortcuts", ns, frac=0.001)
    return out


_SPECS: list[dict[str, Any]] = []


def run_spec(spec: dict[str, Any]) -> dict[str, Any]:
    t0 = time.perf_counter()
    key = tuple(spec["key"])
    adj = build(spec, rng_from_key(key))
    t1 = time.perf_counter()
    summ = edge_coherence_profile(adj, rng_from_key(key + (1,)))
    t2 = time.perf_counter()
    row = {k: v for k, v in spec.items() if k not in ("gen", "args", "key")}
    row["key"] = list(key)
    row["n"] = int(adj.shape[0])
    row["mean_degree"] = float(adj.nnz / adj.shape[0])
    row["summary"] = summ
    row["status"] = coherence_status(summ)
    row["seconds"] = {"gen": round(t1 - t0, 2), "profile": round(t2 - t1, 2)}
    return row


def _run_idx(i: int) -> dict[str, Any]:
    row = run_spec(_SPECS[i])
    s = row["summary"]
    mg = s["median_gamma"]
    print(f"  {row['id']:<28} n={row['n']:<6} deg={row['mean_degree']:.1f} r_w={s['r_w']:<4} "
          f"gamma={'NA' if mg is None else f'{mg:.2f}'} f_deg={s['f_deg']:.2f} f_low={s['f_low']:.2f} "
          f"{row['status']} {row['seconds']}", flush=True)
    return row


# ------------------------------------------------------------------ veredicto


def family_table(rows: list[dict[str, Any]]) -> dict[str, Any]:
    fam: dict[str, dict[str, Any]] = defaultdict(lambda: {"group": None, "n": 0, "n_coherente": 0, "status": [],
                                                          "median_gamma": [], "f_deg": [], "f_low": [], "median_I1": [],
                                                          "r_w": []})
    for r in rows:
        f = fam[r["family"]]
        s = r["summary"]
        f["group"] = r["group"]
        f["n"] += 1
        f["n_coherente"] += int(r["status"] == "COHERENTE")
        f["status"].append(r["status"])
        f["median_gamma"].append(None if s["median_gamma"] is None else round(s["median_gamma"], 3))
        f["f_deg"].append(round(s["f_deg"], 3))
        f["f_low"].append(round(s["f_low"], 3))
        f["median_I1"].append(None if s["median_I1"] is None else round(s["median_I1"], 3))
        f["r_w"].append(s["r_w"])
    return dict(fam)


def verdict(rows: list[dict[str, Any]]) -> dict[str, Any]:
    v = [r for r in rows if r["group"] == "V"]
    f = [r for r in rows if r["group"] == "F"]
    sens = sum(r["status"] == "COHERENTE" for r in v) / len(v)
    spec = sum(r["status"] != "COHERENTE" for r in f) / len(f)
    if sens < 0.7 or spec < 0.7:
        out = "COHERENCIA-INVALIDA"
    elif sens >= 0.9 and spec >= 0.9:
        out = "COHERENCIA-VALIDA"
    else:
        out = "COHERENCIA-PARCIAL"
    return {"verdict": out, "sensitivity": sens, "specificity": spec, "n_V": len(v), "n_F": len(f)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--smoke", action="store_true", help="N~2000, 1 semilla -> results/coh_l0b_smoke/")
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
    _atomic_write(out / "graphs.jsonl", "".join(json.dumps(r, default=str) + "\n" for r in rows))
    table = family_table(rows)
    ver = verdict(rows)
    box = {f: table[f]["status"] for f in ("RGG3_caja_k12", "RGG2_caja_k12")}
    summary: dict[str, Any] = {
        "step": "coh_l0b", "code_commit": head_commit(), "smoke": args.smoke, "n_graphs": len(rows),
        "seconds": round(time.perf_counter() - t0, 1), "sensitivity": ver["sensitivity"], "specificity": ver["specificity"],
        "verdict": ver["verdict"], "n_V": ver["n_V"], "n_F": ver["n_F"], "box_families": box,
        "box_both_coherente": all(s == "COHERENTE" for st in box.values() for s in st),
        "E_group": {f: d["status"] for f, d in table.items() if d["group"] == "E"},
        "by_family": table,
        "graphs": [{"id": r["id"], "group": r["group"], "status": r["status"], "summary": r["summary"]} for r in rows],
    }
    _atomic_write(out / "summary.json", json.dumps(summary, indent=2, default=str))
    print(f"veredicto={ver['verdict']} sens={ver['sensitivity']:.2f} spec={ver['specificity']:.2f} cajas={box}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
