"""P1-D.2 (docs/OMEGA_P1D2_PRERREGISTRO.md): coherencia multiescala (Nivel I) + estabilidad dimensional (Nivel II).

Salida: results/p1d2/{graphs.jsonl, summary.json, level3.json}. Uso: python tools/p1d2_panel.py [--procs 3] [--skip-level3]
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
from collections import defaultdict  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any, Callable  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from c0_dynamics import _atomic_write  # noqa: E402
from l3b_stability import certificate_report  # noqa: E402
from p1_diagnostic import patchwork, rgg_points, shortcuts, square_torus  # noqa: E402

from omega.c0 import references as R  # noqa: E402
from omega.c0.locality import strong_support  # noqa: E402
from omega.config.seeds import SeedKey  # noqa: E402
from omega.controls.random_geometric import rgg_torus  # noqa: E402
from omega.diagnostics.ball_growth import ball_profile, dimension_plateau  # noqa: E402
from omega.experiments.v11.gate import head_commit  # noqa: E402
from omega.landscape.references import connected_caveman  # noqa: E402

MASTER = 20261008
N = 4096
DELTAS = (0.15, 0.10, 0.25)
OUT = ROOT / "results" / "p1d2"


def clique_lattice_3d(side: int, m: int) -> np.ndarray:
    """Toro side^3 de sitios; cada sitio es una K_m. Direccion j (0..5): miembro j%m <-> miembro (j^1)%m del vecino."""
    n = side**3 * m
    a = np.zeros((n, n))
    idx = lambda x, y, z: ((x % side) * side + (y % side)) * side + (z % side)  # noqa: E731
    dirs = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    for x in range(side):
        for y in range(side):
            for z in range(side):
                s = idx(x, y, z)
                blk = np.arange(s * m, (s + 1) * m)
                a[np.ix_(blk, blk)] = 1.0
                for j, (dx, dy, dz) in enumerate(dirs):
                    t = idx(x + dx, y + dy, z + dz)
                    u, v = s * m + j % m, t * m + (j ^ 1) % m
                    a[u, v] = a[v, u] = 1.0
    np.fill_diagonal(a, 0.0)
    return a


def watts_strogatz_fast(n: int, k: int, beta: float, rng: np.random.Generator) -> np.ndarray:
    a = R.ring_lattice(n, k)
    iu, ju = np.nonzero(np.triu(a, 1))
    pick = rng.random(iu.size) < beta
    for u, v in zip(iu[pick], ju[pick]):
        w = int(rng.integers(n))
        if w != u and a[u, w] == 0.0:
            a[u, v] = a[v, u] = 0.0
            a[u, w] = a[w, u] = 1.0
    return a


def triangular_big(side: int) -> np.ndarray:
    return R.triangular_torus(side, side)


SEEDS = (0, 1, 2)
_SPECS: list[dict[str, Any]] = []


def specs() -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []

    def add(fam: str, cls: str, gen: Callable[[np.random.Generator], np.ndarray], fid: int, seeds: tuple[int, ...] = SEEDS) -> None:
        for s in seeds:
            out.append({"id": f"{fam}_s{s}", "family": fam, "cls": cls, "gen": gen, "key": (MASTER, fid, s)})

    add("RGG2_k12", "G", lambda g: rgg_torus(N, 2, 12.0, g), 1)
    add("RGG3_k12", "G", lambda g: rgg_torus(N, 3, 12.0, g), 2)
    add("RGG3_k8", "G", lambda g: rgg_torus(N, 3, 8.0, g), 3)
    add("anillo_k12", "G", lambda g: R.ring_lattice(N, 12), 4, (0,))
    add("cuadrado_64", "G", lambda g: square_torus(64), 5, (0,))
    add("triangular_64", "G", lambda g: triangular_big(64), 6, (0,))
    add("T3_16", "G", lambda g: R.torus_lattice_3d(16), 7, (0,))
    for size in (4, 8, 16):
        add(f"caveman_K{size}", "R", lambda g, size=size: connected_caveman(N, size), 10 + size, (0,))
    add("reticulo_cliques_3D_K8", "R", lambda g: clique_lattice_3d(8, 8), 30, (0,))
    add("retazos_8", "L", lambda g: patchwork(N, 2, g), 40)
    add("retazos_27", "L", lambda g: patchwork(N, 3, g), 41)
    add("retazos_64", "L", lambda g: patchwork(N, 4, g), 42)
    add("RGG3_atajos1", "L", lambda g: shortcuts(N, 0.01, g), 43)
    for beta, fid in ((0.003, 44), (0.01, 45), (0.03, 46)):
        add(f"WS_b{beta}", "L", lambda g, beta=beta: watts_strogatz_fast(N, 12, beta, g), fid)
    add("ER_k12", "NL", lambda g: R.erdos_renyi_m(N, N * 6, g), 50, (0,))
    add("RR_k12", "NL", lambda g: R.random_regular(N, 12, g), 51, (0,))
    for c in (18, 19, 20, 36, 37, 38, 55, 56, 73, 74):
        for s in SEEDS:
            path = ROOT / "runs" / "c0_f1" / "out" / f"F1_c{c:02d}_R_n729_s{s}.npz"
            if path.exists():
                out.append({"id": f"C0_c{c}_s{s}", "family": "C0_R_n729", "cls": "INFO", "npz": str(path), "key": None})
    return out


def run_idx(i: int) -> dict[str, Any]:
    spec = _SPECS[i]
    if "npz" in spec:
        z = np.load(spec["npz"])
        adj = strong_support(np.asarray(z[z.files[0]], dtype=np.float64)).astype(np.float64)
    else:
        adj = spec["gen"](R.rng_from_key(spec["key"]))
    prof = ball_profile(adj)
    row = {k: v for k, v in spec.items() if k not in ("gen", "npz")}
    row.update({"n": int(adj.shape[0]), "mean_degree": float(np.asarray(adj > 0).sum() / adj.shape[0]), **prof})
    row["level2"] = {str(d): dimension_plateau(prof, d) for d in DELTAS}
    return row


def passes(row: dict[str, Any], delta: float) -> dict[str, Any]:
    st1 = row["q4"]["status"]
    lvl1 = st1 in ("HOMOGENEIZA", "CV_CERO")
    st2 = row["level2"][str(delta)]["status"]
    return {"I": lvl1, "II": st2 == "PLATEAU", "joint": lvl1 and st2 == "PLATEAU", "I_status": st1, "II_status": st2,
            "evaluable": st1 != "SIN_VENTANA" and st2 != "SIN_VENTANA_D"}


def verdict(rows: list[dict[str, Any]], delta: float) -> dict[str, Any]:
    dec = [r for r in rows if r["cls"] in ("G", "R", "L", "NL")]
    p = {r["id"]: passes(r, delta) for r in dec}
    g = [r for r in dec if r["cls"] == "G"]
    nongeo = [r for r in dec if r["cls"] in ("L", "NL")]
    sens = sum(p[r["id"]]["joint"] for r in g) / len(g)
    spec = sum(not p[r["id"]]["joint"] for r in nongeo) / len(nongeo)
    n_eval = sum(p[r["id"]]["evaluable"] for r in nongeo)
    fam: dict[str, dict[str, Any]] = defaultdict(lambda: {"cls": None, "n": 0, "I": 0, "II": 0, "joint": 0, "evaluable": 0, "D_plat": [], "rho": []})
    for r in dec:
        f = fam[r["family"]]
        f["cls"] = r["cls"]
        f["n"] += 1
        f["I"] += p[r["id"]]["I"]
        f["II"] += p[r["id"]]["II"]
        f["joint"] += p[r["id"]]["joint"]
        f["evaluable"] += p[r["id"]]["evaluable"]
        pl = r["level2"][str(delta)]["D_plat"]
        f["D_plat"].append(None if pl is None else round(pl, 2))
        f["rho"].append(None if r["q4"]["rho"] is None else round(r["q4"]["rho"], 3))
    r_pass = sorted(k for k, v in fam.items() if v["cls"] == "R" and v["joint"] >= 1)
    if n_eval < 5:
        v = "INDETERMINADO"
    elif sens < 0.7 or spec < 0.7:
        v = "NO_SEPARA"
    elif sens >= 0.9 and spec >= 0.9:
        v = "SEPARA-PARCIAL" if r_pass else "SEPARA"
    else:
        v = "INDETERMINADO"
    return {"delta": delta, "verdict": v, "sensitivity": sens, "specificity": spec, "nongeo_evaluable": n_eval,
            "nongeo_total": len(nongeo), "R_families_passing_joint": r_pass, "by_family": dict(fam)}


def level3(rows: list[dict[str, Any]]) -> dict[str, Any]:
    r_pass = verdict(rows, DELTAS[0])["R_families_passing_joint"]
    analogs: dict[str, np.ndarray] = {}
    for f in r_pass:
        if f.startswith("caveman_K"):
            size = int(f.split("K")[1])
            analogs[f"{f}_n728"] = connected_caveman(728, size)
    if "reticulo_cliques_3D_K8" in r_pass:
        analogs["reticulo_cliques_3D_K6_n750"] = clique_lattice_3d(5, 6)
    refs = {"RGG3_k12_n729": rgg_torus(729, 3, 12.0, R.rng_from_key((MASTER, 60, 0))),
            "RGG2_k12_n729": rgg_torus(729, 2, 12.0, R.rng_from_key((MASTER, 61, 0))),
            "anillo_k12_n729": R.ring_lattice(729, 12)}
    res: dict[str, Any] = {"R_pass_family_list": r_pass, "analogs": {}, "references": {}}
    for tag, bucket in (("analogs", analogs), ("references", refs)):
        for name, adj in bucket.items():
            res[tag][name] = [certificate_report(adj.astype(np.float64), "CONVERGED", SeedKey(MASTER, (70, s))).get("codes") for s in range(3)]
    return res


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--skip-level3", action="store_true")
    ap.add_argument("--only-level3", action="store_true", help="P1D2-A1: solo Nivel III sobre results/p1d2/graphs.jsonl")
    args = ap.parse_args(argv)
    if args.only_level3:
        rows0 = [json.loads(line) for line in (OUT / "graphs.jsonl").read_text().splitlines() if line.strip()]
        l3 = level3(rows0)
        _atomic_write(OUT / "level3.json", json.dumps(l3, indent=2, default=str))
        print(json.dumps(l3, default=str))
        return 0
    _SPECS[:] = specs()
    with mp.get_context("fork").Pool(args.procs) as pool:
        rows = pool.map(run_idx, range(len(_SPECS)), chunksize=1)
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "graphs.jsonl").open("w", encoding="utf-8") as fh:
        for r in rows:
            fh.write(json.dumps(r, default=str) + "\n")
    res = {str(d): verdict(rows, d) for d in DELTAS}
    summary: dict[str, Any] = {"step": "p1d2_panel", "code_commit": head_commit(), "n_graphs": len(rows), "result": res}
    _atomic_write(OUT / "summary.json", json.dumps(summary, indent=2, default=str))
    if not args.skip_level3:
        l3 = level3(rows)
        _atomic_write(OUT / "level3.json", json.dumps(l3, indent=2, default=str))
        print("Nivel III:", json.dumps(l3, default=str))
    for d, v in res.items():
        print(f"delta={d}: {v['verdict']} sens={v['sensitivity']:.2f} spec={v['specificity']:.2f} R_pasan={v['R_families_passing_joint']}")
    for f, d in res[str(DELTAS[0])]["by_family"].items():
        print(f"  {f:<24} {d['cls']:<3} I={d['I']}/{d['n']} II={d['II']}/{d['n']} joint={d['joint']}/{d['n']} Dplat={d['D_plat']} rho={d['rho']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
