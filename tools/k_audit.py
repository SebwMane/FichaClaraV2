"""Auditoria (b) de curvatura de Ollivier sobre el panel P1-D.3 (docs/OMEGA_K_AUDIT_PRERREGISTRO.md).

Salida: results/k_audit/{graphs.jsonl, summary.json}. Uso:
  python tools/k_audit.py [--procs 3] [--smoke] [--only-index I]
Regenera los grafos con la misma clave que tools/p1d3_panel.py::run_spec (rng_from_key(spec["key"])), muestrea 600
aristas con rng_from_key(key + (4,)) y calcula kappa perezosa (idleness 0.5, coste en saltos) con ollivier_edge_sparse.
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
from scipy.stats import spearmanr  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import p1d3_panel  # noqa: E402
from c0_dynamics import _atomic_write  # noqa: E402

from omega.c0.references import rng_from_key  # noqa: E402
from omega.curvature.ollivier_sparse import ollivier_edge_sparse  # noqa: E402
from omega.experiments.v11.gate import head_commit  # noqa: E402

IDLENESS = 0.5
N_EDGES = 600
NEG_CUT = -0.1
POS_CUT = 0.1
BAND_PAD = 0.02
FNEG_PAD = 0.05
KA1_THRESHOLD = 0.7
OUTSIDE_FAMILIES = ("RGG3_atajos_0.1", "ER_k12", "WS_b0.001", "WS_b0.002")
INSIDE_FAMILIES = ("RGG3_caja_k12", "RGG2_caja_k12")
HARD = OUTSIDE_FAMILIES + INSIDE_FAMILIES
OUT_FULL = ROOT / "results" / "k_audit"
OUT_SMOKE = ROOT / "results" / "k_audit_smoke"
SMOKE_FAMILIES = ("RGG3_k12", "caveman_K8", "ER_k12")

_SPECS: list[dict[str, Any]] = []


def kappa_sample(adj: sparse.csr_array, key: tuple[int, ...]) -> np.ndarray:
    """kappa en N_EDGES aristas (sin reemplazo, PCG64 de key + (4,)) del triangulo superior; todas si hay menos."""
    tri = sparse.triu(adj, k=1, format="coo")
    xs, ys = np.asarray(tri.row, dtype=np.int64), np.asarray(tri.col, dtype=np.int64)
    if xs.size > N_EDGES:
        pick = np.sort(rng_from_key(key + (4,)).choice(xs.size, size=N_EDGES, replace=False))
        xs, ys = xs[pick], ys[pick]
    return np.asarray([ollivier_edge_sparse(adj, int(x), int(y), IDLENESS) for x, y in zip(xs, ys, strict=True)],
                      dtype=np.float64)


def run_spec(spec: dict[str, Any]) -> dict[str, Any]:
    t0 = time.perf_counter()
    key = tuple(int(k) for k in spec["key"])
    adj = p1d3_panel.build(spec, rng_from_key(key))
    t1 = time.perf_counter()
    vals = kappa_sample(adj, key)
    t2 = time.perf_counter()
    return {
        "id": spec["id"], "family": spec["family"], "cls": spec["cls"], "panel": spec["panel"], "key": list(key),
        "n": int(adj.shape[0]), "mean_degree": float(adj.nnz / adj.shape[0]), "n_edges_sampled": int(vals.size),
        "kappa_median": float(np.median(vals)), "kappa_mean": float(vals.mean()),
        "f_neg": float(np.mean(vals < NEG_CUT)), "f_pos": float(np.mean(vals > POS_CUT)),
        "seconds": {"gen": round(t1 - t0, 2), "kappa": round(t2 - t1, 2)},
    }


def _run_idx(i: int) -> dict[str, Any]:
    row = run_spec(_SPECS[i])
    print(f"  {row['id']:<28} n={row['n']:<6} deg={row['mean_degree']:.2f} med={row['kappa_median']:+.3f} "
          f"fneg={row['f_neg']:.3f} fpos={row['f_pos']:.3f} {row['seconds']}", flush=True)
    return row


# ------------------------------------------------------------------ analisis


def _finite(x: Any) -> bool:
    return x is not None and math.isfinite(float(x))


def ka1(rows: list[dict[str, Any]]) -> dict[str, Any]:
    sel = [r for r in rows if _finite(r.get("rho")) and float(r["rho"]) > 0.0]
    if len(sel) < 3:
        return {"n": len(sel), "spearman": None, "p": None, "non_redundant": None}
    res = spearmanr([r["kappa_median"] for r in sel], [math.log(float(r["rho"])) for r in sel])
    c = float(res.statistic)
    return {"n": len(sel), "spearman": c, "p": float(res.pvalue), "non_redundant": bool(abs(c) < KA1_THRESHOLD)}


def ka2(rows: list[dict[str, Any]]) -> dict[str, Any]:
    hom = [r for r in rows if r["cls"] == "G-hom"]
    meds = [r["kappa_median"] for r in hom]
    lo, hi = min(meds) - BAND_PAD, max(meds) + BAND_PAD
    fneg_max = max(r["f_neg"] for r in hom) + FNEG_PAD

    def inside(r: dict[str, Any]) -> bool:
        return bool(lo <= r["kappa_median"] <= hi and r["f_neg"] <= fneg_max)

    fam: dict[str, Any] = {}
    for f in HARD:
        rs = [r for r in rows if r["family"] == f]
        want_inside = f in INSIDE_FAMILIES
        good = sum(int(inside(r) == want_inside) for r in rs)
        fam[f] = {"expected": "DENTRO" if want_inside else "FUERA", "n": len(rs), "n_correct_side": good,
                  "inside": [inside(r) for r in rs], "kappa_median": [round(r["kappa_median"], 4) for r in rs],
                  "f_neg": [round(r["f_neg"], 4) for r in rs],
                  "corrects": bool(rs and good * 3 >= 2 * len(rs))}
    hard_ids = set(OUTSIDE_FAMILIES)
    rest = [r for r in rows if r["cls"] in ("L", "NL") and r["family"] not in hard_ids]
    cr = [r for r in rows if r["cls"] == "R"]

    def frac_out(rs: list[dict[str, Any]]) -> float | None:
        return float(np.mean([not inside(r) for r in rs])) if rs else None

    return {
        "band": {"kappa_median": [lo, hi], "f_neg_max": fneg_max, "G_hom_n": len(hom)},
        "families": fam,
        "n_corrected": sum(int(v["corrects"]) for v in fam.values()),
        "control": {
            "G_hom_fraction_inside": float(np.mean([inside(r) for r in hom])),
            "R_fraction_outside": frac_out(cr), "R_n": len(cr),
            "LNL_rest_fraction_outside": frac_out(rest), "LNL_rest_n": len(rest),
            "LNL_rest_ids": [r["id"] for r in rest],
            "LNL_rest_outside_by_id": {r["id"]: (not inside(r)) for r in rest},
            "G_inh_other_inside_by_id": {r["id"]: inside(r) for r in rows
                                         if r["cls"] == "G-inh" and r["family"] not in INSIDE_FAMILIES},
        },
    }


def verdict(k1: dict[str, Any], k2: dict[str, Any]) -> str:
    nonred = k1["non_redundant"]
    corrected = int(k2["n_corrected"])
    if nonred and corrected >= 1:
        return "K APORTA"
    if nonred is False and corrected == 0:
        return "K REDUNDANTE"
    return "K MARGINAL"


def join_battery(rows: list[dict[str, Any]], path: Path) -> None:
    bat: dict[str, dict[str, Any]] = {}
    if path.exists():
        for line in path.read_text(encoding="utf-8").splitlines():
            if line.strip():
                b = json.loads(line)
                bat[b["id"]] = b
    for r in rows:
        b = bat.get(r["id"])
        if b is None:
            r["regen_ok"] = None
            continue
        r["regen_ok"] = bool(r["n"] == b["n"] and abs(r["mean_degree"] - b["mean_degree"]) < 1e-3)
        r["rho"] = b["level1"]["rho"]
        r["level1_status"] = b["level1"]["status"]
        r["level2_status"] = b["level2"]["0.05"]["status"]
        r["D_conv"] = b["level2"]["0.05"]["D_conv"]
        r["D_s"] = b["spectral"]["D_s_median_upper"]
        r["category"] = b["category"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--only-index", type=int, default=None, help="solo ese indice de specs() (cronometraje)")
    a = ap.parse_args()
    global _SPECS
    if a.smoke:
        _SPECS = [s for s in p1d3_panel.specs(smoke=True) if s["family"] in SMOKE_FAMILIES]
        out, bat = OUT_SMOKE, ROOT / "results" / "p1d3_smoke" / "graphs.jsonl"
    else:
        _SPECS = p1d3_panel.specs()
        out, bat = OUT_FULL, ROOT / "results" / "p1d3" / "graphs.jsonl"
    idx = list(range(len(_SPECS))) if a.only_index is None else [a.only_index]
    t0 = time.perf_counter()
    order = sorted(idx, key=lambda i: -int(_SPECS[i]["N"]))
    if a.procs > 1 and len(order) > 1:
        with mp.get_context("fork").Pool(a.procs) as pool:
            res = list(pool.imap_unordered(_run_idx, order))
    else:
        res = [_run_idx(i) for i in order]
    pos = {s["id"]: i for i, s in enumerate(_SPECS)}
    rows = sorted(res, key=lambda r: pos[r["id"]])
    if a.only_index is not None:
        print(f"tiempo total {time.perf_counter() - t0:.1f}s")
        return
    join_battery(rows, bat)
    _atomic_write(out / "graphs.jsonl", "".join(json.dumps(r) + "\n" for r in rows))
    summary: dict[str, Any] = {
        "commit": head_commit(), "smoke": a.smoke, "n_graphs": len(rows), "idleness": IDLENESS, "n_edges": N_EDGES,
        "regen_mismatches": [r["id"] for r in rows if r.get("regen_ok") is False],
        "seconds": round(time.perf_counter() - t0, 1),
    }
    if not a.smoke:
        k1, k2 = ka1(rows), ka2(rows)
        summary.update({"KA1": k1, "KA2": k2, "verdict": verdict(k1, k2)})
        print(f"KA1 {k1}\nKA2 corrected={k2['n_corrected']} control={k2['control']['G_hom_fraction_inside']}"
              f" -> {summary['verdict']}")
    _atomic_write(out / "summary.json", json.dumps(summary, indent=2, default=str) + "\n")


if __name__ == "__main__":
    main()
