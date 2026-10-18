"""OMEGA W6: calibracion de la condicion de uso de RC-3 a tres escalas (docs/OMEGA_W6_PRERREGISTRO.md, congelado).

Salida: results/w6/{runs.jsonl, summary.json}. Uso: python tools/w6_calib.py [--procs 3] [--fresh]
Unidad = (familia, indice de tamano); runs.jsonl se anade por unidad (REANUDABLE; --fresh lo borra).
RNG solo via omega.c0.references.rng_from_key; clave (MASTER_W6, id_familia, semilla, indice_de_tamano); la misma clave se
pasa a rc3.measure. Medida: rc3.measure sin cambios; R = rc3.half_mass_radius (dentro de measure) y delta = rc3.scaling_delta
sobre (R_i, R_j, n_giant_i, n_giant_j), como rc3.run_pair. Controles J3 (semillas 0-2): se REUTILIZAN los registros I2 de
results/r3_0b_inspect/runs.jsonl (no se vuelve a medir).

Decisiones de implementacion (lectura mas literal de las ambiguedades)
----------------------------------------------------------------------
* Tamanos objetivo 1e4, 8e4, 6.4e5 nodos. Se elige el parametro de tamano cuyo N real este mas cerca del objetivo en
  escala logaritmica (|ln(N/objetivo)|); se registra el N real.
* J1(w): L(J1(w)) con cadenas de n elementos es la caja de rejilla [0,n]^w, cuyo grafo de cobertura es la rejilla abierta de
  lado n+1 (identica a la construccion de tools/r3_0b.py; se construye directamente con a0_rc2.lattice_open). Las cajas
  quedan anidadas por longitud de cadena. (n+1) = lado mas cercano a objetivo^(1/w).
* Toros Z^d: p1d3_panel.torus_lattice(lado, d), lado mas cercano a objetivo^(1/d). Triangular: rc3.triangular_open(lado),
  N = lado^2. FCC: rc3.fcc_lattice(L), N = 4 L^3 (periodica, como en rc3).
* RGG_d (d=2,3): caja abierta [0,1]^d, grado medio objetivo k=8 en el bulk. Se generan UNA vez N3 = 640000 puntos uniformes
  con rng_from_key((MASTER_W6, id, 0, 2)) (semilla 0). El radio sale de la densidad de N3: r = (k / (N3 * V_d))^(1/d), V_d
  volumen de la bola unidad (se usa N3, no N3-1, por "radio desde la densidad"). Para el tamano i se toman los puntos en la
  subcaja [0, (N_i/N3)^(1/d)]^d (misma esquina; N_i nominal 1e4, 8e4, 6.4e5) y se construye el RGG con el mismo radio
  (cKDTree.query_pairs, sin periodicidad). El N real es el numero de puntos contenidos. La clave de MEDIDA es
  (MASTER_W6, id, 0, i) para cada tamano.
* k_i = grado medio (nnz/n) de la componente GIGANTE, calculado aqui (rc3.measure devuelve el grado medio del grafo
  completo, que se guarda aparte como mean_degree_graph). J3: L(J) es conexo (n_giant == n en los tres tamanos y semillas),
  luego k = mean_degree del registro I2.
* E/N para W5 = mean_degree_graph/2 del grafo completo, como r3_0b_inspect (para J3, el campo E_over_N del registro I2).
  W5 valido en un par si |E/N_b - E/N_a| / E/N_a < 0.25 (estricto, como r3_0b_inspect). W6.4 = ambos pares.
* delta_12 = scaling_delta(R1,R2,ng1,ng2), delta_23 = scaling_delta(R2,R3,ng2,ng3). Si algun delta es None (no evaluable),
  W6.2 NO se cumple (no se puede afirmar |d23-d12|<=tau) y la familia no aporta a max|Dd| de los controles.
* W6.1 se cumple por construccion para todas las familias (anidadas o deterministas); se registra True.
* W6.3: |k3-k2|/k2 <= tau_k (<=, como en el texto). W6.2: |d23-d12| <= tau_d.
* "Cuando haya mas de tres tamanos..." (pares de pares consecutivos y prohibicion de delta descendente en todos los pasos):
  con exactamente tres tamanos hay un solo par de pares y un solo paso de delta, por lo que no hay pares adicionales; la
  prohibicion de tendencia no se aplica (se registra delta_descends solo como dato).
* tau_d = max(0.02, 1.5*max|Dd|) y tau_k = max(0.05, 1.5*max|Dk|/k2) SOLO sobre los controles finitos (10 familias:
  J1(2), J1(3), J1(4), Z2, Z3, Z4, triangular, FCC, RGG2, RGG3) con valor definido.
* Validez W6 de cada familia = W6.1 y W6.2 y W6.3 y W6.4. Los controles finitos que no sean W6-validos se registran como
  anomalia (W5 puede fallarles; el umbral no se ajusta).
* P1: las 3 semillas J3 W6-invalidas. P2: el control "que fija" cada tau es el argmax (solo si 1.5*max supera el suelo;
  si manda el suelo no hay control que fije ese tau); se recalculan AMBOS umbrales sin ese control y si el veredicto de algun
  J3 pasa de invalido a valido -> 'potencia fragil' (se hace una prueba por cada control fijador distinto). P3: tau_d <= 0.10.
  Veredicto: VALIDA si P1; si no INVALIDA. P2/P3 solo condicionan la lectura.
* Exclusiones RC-3 por par con rc3.pair_exclusions(delta, (ann_a, ann_b), (coh_a, coh_b)) usando annulus_status y
  coherence_status de rc3.measure (informativas; no entran en la validez W6).
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
from scipy.sparse.csgraph import connected_components  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from a0_rc2 import lattice_open  # noqa: E402
from c0_dynamics import _atomic_write  # noqa: E402
from p1d3_panel import _from_edges, _unit_ball, torus_lattice  # noqa: E402
from rc3 import fcc_lattice, measure, pair_exclusions, scaling_delta, triangular_open  # noqa: E402

from omega.c0.references import rng_from_key  # noqa: E402
from omega.experiments.v11.gate import head_commit  # noqa: E402

MASTER_W6 = 20261022
OUT = ROOT / "results" / "w6"
INSPECT = ROOT / "results" / "r3_0b_inspect" / "runs.jsonl"
TARGETS = (10_000, 80_000, 640_000)
RGG_K = 8.0
N3_RGG = TARGETS[2]
TAU_D_FLOOR, TAU_K_FLOOR, MULT = 0.02, 0.05, 1.5
W5_TOL = 0.25
P3_MAX = 0.10

# (id_str, id_num, tipo, parametro)
FAMS: list[tuple[str, int, str, int]] = [
    ("J1_w2", 1, "box", 2), ("J1_w3", 2, "box", 3), ("J1_w4", 3, "box", 4),
    ("Z2", 4, "torus", 2), ("Z3", 5, "torus", 3), ("Z4", 6, "torus", 4),
    ("triangular_caja", 7, "tri", 2), ("FCC", 8, "fcc", 3),
    ("RGG2_k8", 9, "rgg", 2), ("RGG3_k8", 10, "rgg", 3),
]
FAMILY_ID = {f[0]: f[1] for f in FAMS}
FAMILY = {f[0]: f for f in FAMS}
CONTROLS = [f[0] for f in FAMS]
J3_IDS = ["J3_s0", "J3_s1", "J3_s2"]


# ------------------------------------------------------------------ construcciones


def closest_param(count_fn: Any, guess: float, target: int) -> int:
    """Parametro entero p (>=2) cuyo N = count_fn(p) esta mas cerca de target en |ln(N/target)|."""
    base = max(2, int(math.floor(guess)))
    cands = [p for p in (base - 1, base, base + 1, base + 2) if p >= 2]
    return min(cands, key=lambda p: abs(math.log(count_fn(p) / target)))


def size_param(kind: str, d: int, target: int) -> int:
    if kind in ("box", "torus"):
        return closest_param(lambda p: p**d, target ** (1.0 / d), target)
    if kind == "tri":
        return closest_param(lambda p: p * p, math.sqrt(target), target)
    if kind == "fcc":
        return closest_param(lambda p: 4 * p**3, (target / 4.0) ** (1.0 / 3.0), target)
    raise ValueError(kind)


def rgg_radius(d: int, k: float = RGG_K, n3: int = N3_RGG) -> float:
    return float((k / (n3 * _unit_ball(d))) ** (1.0 / d))


def rgg_points_all(d: int, fid: int, n3: int = N3_RGG) -> np.ndarray:
    """N3 puntos uniformes en [0,1]^d (semilla 0, tamano 2 = N3)."""
    return rng_from_key((MASTER_W6, fid, 0, 2)).random((n3, d))


def rgg_subbox(x: np.ndarray, frac: float, d: int) -> np.ndarray:
    """Puntos en la subcaja [0, frac^(1/d)]^d (misma esquina), en el orden original."""
    s = frac ** (1.0 / d)
    return x[np.all(x < s, axis=1)]


def rgg_open(pts: np.ndarray, r: float) -> sparse.csr_array:
    p = cKDTree(pts).query_pairs(r, output_type="ndarray")
    return _from_edges(pts.shape[0], p[:, 0].astype(np.int64), p[:, 1].astype(np.int64))


def build(name: str, si: int) -> tuple[sparse.csr_array, dict[str, Any]]:
    _, fid, kind, d = FAMILY[name]
    target = TARGETS[si]
    info: dict[str, Any] = {"target": target}
    if kind in ("box", "torus", "tri", "fcc"):
        p = size_param(kind, d, target)
        info["param"] = p
        if kind == "box":
            adj = lattice_open(p, d)
            info["n_chain"] = p - 1
        elif kind == "torus":
            adj = torus_lattice(p, d)
        elif kind == "tri":
            adj = triangular_open(p)
        else:
            adj = fcc_lattice(p)
        return adj, info
    x = rgg_points_all(d, fid)
    pts = rgg_subbox(x, target / N3_RGG, d)
    r = rgg_radius(d)
    info.update(radius=r, n_points=int(pts.shape[0]))
    return rgg_open(pts, r), info


def giant_degree(adj: sparse.csr_array) -> tuple[float, int]:
    """(grado medio nnz/n de la componente gigante, n_giant)."""
    _, lab = connected_components(adj, directed=False)
    big = int(np.argmax(np.bincount(lab)))
    keep = np.flatnonzero(lab == big)
    sub = sparse.csr_array(adj[keep][:, keep])
    return float(sub.nnz / keep.size), int(keep.size)


# ------------------------------------------------------------------ unidades


def run_unit(u: tuple[str, int]) -> dict[str, Any]:
    name, si = u
    fid = FAMILY_ID[name]
    key = (MASTER_W6, fid, 0, si)
    t0 = time.perf_counter()
    adj, info = build(name, si)
    gen_s = round(time.perf_counter() - t0, 2)
    kg, ng = giant_degree(adj)
    m = measure(adj, key)
    del adj
    m.pop("annulus", None)
    m.pop("coherence", None)
    m.update(info)
    m.update(id=f"{name}_i{si}", family=name, size_index=si, key=list(key), k_giant=kg, n_giant_check=ng,
             mean_degree_graph=m["mean_degree"], E_over_N=m["mean_degree"] / 2.0)
    m["seconds"]["gen"] = gen_s
    m["seconds_total"] = round(time.perf_counter() - t0, 1)
    return m


def load_j3() -> dict[str, list[dict[str, Any]]]:
    out: dict[str, list[dict[str, Any]]] = {}
    for line in INSPECT.read_text().splitlines():
        r = json.loads(line)
        if r.get("unit") != "I2" or r.get("status") != "OK":
            continue
        sz = []
        for i, s in enumerate(r["sizes"]):
            assert s["n_giant"] == s["n"], "J3 no conexo: k_giant != mean_degree"
            sz.append({"id": f"J3_s{r['seed']}_i{i}", "family": f"J3_s{r['seed']}", "size_index": i, "key": s["key"],
                       "n": s["n"], "n_giant": s["n_giant"], "R": s["R"], "k_giant": s["mean_degree"],
                       "mean_degree_graph": s["mean_degree"], "E_over_N": s["E_over_N"],
                       "annulus_status": s["annulus_status"], "coherence_status": s["coherence_status"],
                       "source": "results/r3_0b_inspect/runs.jsonl (I2)"})
        out[f"J3_s{r['seed']}"] = sz
    return out


# ------------------------------------------------------------------ analisis


def w5_valid(ea: float, eb: float) -> bool:
    return bool(abs(eb - ea) / ea < W5_TOL)


def family_metrics(sizes: list[dict[str, Any]]) -> dict[str, Any]:
    a, b, c = sizes
    d12 = scaling_delta(a["R"], b["R"], a["n_giant"] or 0, b["n_giant"] or 0)
    d23 = scaling_delta(b["R"], c["R"], b["n_giant"] or 0, c["n_giant"] or 0)
    k1, k2, k3 = a["k_giant"], b["k_giant"], c["k_giant"]
    pairs = []
    for x, y, dl in ((a, b, d12), (b, c, d23)):
        pairs.append({"w5": w5_valid(x["E_over_N"], y["E_over_N"]),
                      "exclusions": pair_exclusions(dl, (x["annulus_status"], y["annulus_status"]),
                                                    (x["coherence_status"], y["coherence_status"]))})
    dd = None if d12 is None or d23 is None else abs(d23 - d12)
    return {"N": [s["n"] for s in sizes], "n_giant": [s["n_giant"] for s in sizes], "R": [s["R"] for s in sizes],
            "delta12": d12, "delta23": d23, "abs_ddelta": dd,
            "delta_descends": None if dd is None else bool(d23 < d12),
            "k1": k1, "k2": k2, "k3": k3, "rel_dk": abs(k3 - k2) / k2,
            "rel_dk12": abs(k2 - k1) / k1, "E_over_N": [s["E_over_N"] for s in sizes],
            "w5_pair1": pairs[0]["w5"], "w5_pair2": pairs[1]["w5"],
            "exclusions_pair1": pairs[0]["exclusions"], "exclusions_pair2": pairs[1]["exclusions"],
            "status": [[s["annulus_status"], s["coherence_status"]] for s in sizes]}


def thresholds(ctrl: list[dict[str, Any]]) -> dict[str, Any]:
    """tau_d, tau_k por la regla congelada, sobre las metricas de los controles finitos dados."""
    dds = [(m["abs_ddelta"], m["family"]) for m in ctrl if m["abs_ddelta"] is not None]
    dks = [(m["rel_dk"], m["family"]) for m in ctrl]
    max_d, arg_d = max(dds) if dds else (0.0, None)
    max_k, arg_k = max(dks) if dks else (0.0, None)
    td, tk = max(TAU_D_FLOOR, MULT * max_d), max(TAU_K_FLOOR, MULT * max_k)
    return {"tau_d": td, "tau_k": tk, "max_ddelta": max_d, "max_rel_dk": max_k,
            "setter_d": arg_d if MULT * max_d > TAU_D_FLOOR else None,
            "setter_k": arg_k if MULT * max_k > TAU_K_FLOOR else None,
            "arg_max_d": arg_d, "arg_max_k": arg_k}


def w6_rule(m: dict[str, Any], tau_d: float, tau_k: float) -> dict[str, Any]:
    w61 = True
    w62 = m["abs_ddelta"] is not None and m["abs_ddelta"] <= tau_d
    w63 = m["rel_dk"] <= tau_k
    w64 = bool(m["w5_pair1"] and m["w5_pair2"])
    return {"W6.1": w61, "W6.2": bool(w62), "W6.3": bool(w63), "W6.4": w64, "valid": bool(w61 and w62 and w63 and w64)}


def calibrate(ctrl: list[dict[str, Any]], j3: list[dict[str, Any]]) -> dict[str, Any]:
    th = thresholds(ctrl)
    ctrl_rules = {m["family"]: w6_rule(m, th["tau_d"], th["tau_k"]) for m in ctrl}
    j3_rules = {m["family"]: w6_rule(m, th["tau_d"], th["tau_k"]) for m in j3}
    p1 = all(not r["valid"] for r in j3_rules.values()) and len(j3_rules) == 3
    # P2: leave-one-out sobre el/los control(es) que fijan tau
    setters = sorted({s for s in (th["setter_d"], th["setter_k"]) if s is not None})
    loo, fragile = [], False
    for s in setters:
        rest = [m for m in ctrl if m["family"] != s]
        t2 = thresholds(rest)
        rules = {m["family"]: w6_rule(m, t2["tau_d"], t2["tau_k"])["valid"] for m in j3}
        changed = [f for f, v in rules.items() if v != j3_rules[f]["valid"]]
        fragile = fragile or bool(changed)
        loo.append({"removed": s, "tau_d": t2["tau_d"], "tau_k": t2["tau_k"], "j3_valid": rules, "changed": changed})
    p2 = not fragile
    p3 = th["tau_d"] <= P3_MAX
    return {"thresholds": th, "controls_w6": ctrl_rules, "j3_w6": j3_rules,
            "P1": p1, "P2": p2, "P2_detail": loo, "P2_label": "OK" if p2 else "potencia fragil",
            "P3": p3, "P3_label": "OK" if p3 else "W6.2 no informativa",
            "controls_not_valid": [f for f, r in ctrl_rules.items() if not r["valid"]],
            "verdict": "W6 VALIDA" if p1 else "W6 INVALIDA"}


# ------------------------------------------------------------------ E/S


def _load(path: Path) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                out[r["id"]] = r
    return out


def fmt_table(fams: dict[str, dict[str, Any]], rules: dict[str, dict[str, Any]]) -> str:
    def f(x: Any, p: int = 4) -> str:
        return "None" if x is None else f"{x:.{p}f}"
    lines = ["family            N1/N2/N3             d12     d23     |Dd|    k1      k2      k3      |Dk|/k2 W5a W5b  2   3   valid"]
    for n, m in fams.items():
        r = rules[n]
        lines.append(f"{n:17s} {'/'.join(str(x) for x in m['N']):21s} {f(m['delta12'])}  {f(m['delta23'])}  {f(m['abs_ddelta'])}  "
                     f"{f(m['k1'], 3)}  {f(m['k2'], 3)}  {f(m['k3'], 3)}  {f(m['rel_dk'])}  {int(m['w5_pair1'])}   {int(m['w5_pair2'])}   "
                     f"{int(r['W6.2'])}   {int(r['W6.3'])}   {r['valid']}  X:{m['exclusions_pair1']}{m['exclusions_pair2']}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--fresh", action="store_true")
    args = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    rj = OUT / "runs.jsonl"
    if args.fresh and rj.exists():
        rj.unlink()
    done = _load(rj)
    units = [(n, si) for si in (2, 1, 0) for n in CONTROLS]  # grandes primero
    todo = [u for u in units if f"{u[0]}_i{u[1]}" not in done]
    t0 = time.perf_counter()
    print(f"{len(units)} unidades ({len(done)} hechas), procs={args.procs}", flush=True)
    with mp.get_context("fork").Pool(args.procs) as pool, open(rj, "a") as fh:
        for r in pool.imap_unordered(run_unit, todo, chunksize=1):
            fh.write(json.dumps(r, default=str) + "\n")
            fh.flush()
            done[r["id"]] = r
            print(f"  {r['id']} n={r['n']} ng={r['n_giant']} k={r['k_giant']:.3f} R={r['R']} {r['seconds_total']}s", flush=True)
    j3 = load_j3()
    # los registros J3 reutilizados se anaden a runs.jsonl una sola vez (marcados source)
    with open(rj, "a") as fh:
        for sz in j3.values():
            for s in sz:
                if s["id"] not in done:
                    fh.write(json.dumps(s, default=str) + "\n")
                    done[s["id"]] = s
    fams: dict[str, dict[str, Any]] = {}
    for n in CONTROLS:
        fams[n] = {"family": n, **family_metrics([done[f"{n}_i{i}"] for i in range(3)])}
    for n, sz in j3.items():
        fams[n] = {"family": n, **family_metrics(sz)}
    cal = calibrate([fams[n] for n in CONTROLS], [fams[n] for n in J3_IDS])
    rules = {**cal["controls_w6"], **cal["j3_w6"]}
    table = fmt_table(fams, rules)
    summary = {"step": "w6_calib", "code_commit": head_commit(), "master": MASTER_W6, "seconds_this_run": round(time.perf_counter() - t0, 1),
               "families": fams, "calibration": cal, "table": table.splitlines()}
    _atomic_write(OUT / "summary.json", json.dumps(summary, indent=2, default=str))
    print(table)
    print(json.dumps({k: cal[k] for k in ("thresholds", "P1", "P2", "P2_label", "P2_detail", "P3", "P3_label", "controls_not_valid", "verdict")},
                     indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
