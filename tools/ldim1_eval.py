"""Omega L-DIM-1 (docs/OMEGA_LDIM1_PRERREGISTRO.md §7 + A1): evaluador. Unico script que lee truth.jsonl.

Uso: python tools/ldim1_eval.py [--smoke]   (lee results/ldim1[_smoke]/{features,truth}.jsonl, escribe summary.json)

Decisiones de lectura (A1 las deja abiertas):
 - empates en argmin_d: gana la d menor; empates en la moda de NC-5: gana la d menor.
 - la meseta (L3) es la racha mas larga de indices consecutivos de la rejilla con L1&L2 y el mismo d*; empate -> la primera.
   Pasa si (longitud - 1) * 0.25 >= 1 decada. L4, L7 y NC-2 se evaluan en los theta de esa racha; sin meseta fallan (NC-2 sin bandera).
 - L4 usa F(d*) = media de F de las geometrias de dimension d* a tamano N (todos los estratos).
 - Nulo A: z = (log x - media)/desv sobre todas las instancias de tamano N (geometrias + degenerados); phi se usa como log(phi + 1e-6)
   porque phi = 0 en la clique (psi: log psi). Tasa = fraccion que cumple L1 + L2 + L4.
 - Nulos B y C: se cuenta aprobado si algun theta cumple L1 + L2 + L3 (el nulo B solo con p fijo, tasa por p).
 - L9 (ceguera) = True por construccion: extract() solo recibe (adj, id, rng) y truth.jsonl solo se lee aqui.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from c0_dynamics import _atomic_write  # noqa: E402

from omega.c0.references import rng_from_key  # noqa: E402
from omega.experiments.v11.gate import head_commit  # noqa: E402

MASTER = 20261014
THETAS = np.logspace(-3, 3, 25)
P_VALUES = (0.5, 1.0, 2.0)
DECADE_STEP = 0.25  # decadas por paso de la rejilla
STRATA = ("reticulo", "k8", "k16")
INTERIOR = (2, 3, 4)
EPS_PHI = 1e-6
M_A, N_B, N_C = 500, 200, 200
NULL_MAX = 0.05
RHO_FLAG = 0.9
NC4_TOL = 0.10
LAW = {0.5: 4, 1.0: 3, 2.0: 2}  # L-DIM-T1
OUT_FULL = ROOT / "results" / "ldim1"
OUT_SMOKE = ROOT / "results" / "ldim1_smoke"


# ------------------------------------------------------------------ nucleo vectorizado (filas = candidatos)


def dstar_rows(F: np.ndarray, d: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """argmin_d de la media de F (filas) sobre las instancias de `mask` con etiqueta d; 0 si `mask` es vacia."""
    labels = sorted({int(x) for x in d[mask]})
    if not labels:
        return np.zeros(F.shape[0], dtype=int)
    means = np.stack([F[:, mask & (d == lab)].mean(axis=1) for lab in labels], axis=1)
    return np.asarray(labels)[np.argmin(means, axis=1)]


def l1_l2(F: np.ndarray, d: np.ndarray, stratum: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(d* combinado, L1, L2) por fila. `d`, `stratum` etiquetan las columnas (geometrias de tamano N)."""
    everyone = np.ones(d.size, dtype=bool)
    pooled = dstar_rows(F, d, everyone)
    l1 = np.isin(pooled, INTERIOR)
    l2 = np.ones(F.shape[0], dtype=bool)
    for s in STRATA:
        m = stratum == s
        if not m.any():
            l2[:] = False
            continue
        l2 &= dstar_rows(F, d, m) == pooled
    return pooled, l1, l2


def longest_run(ok: np.ndarray, dst: np.ndarray) -> tuple[int, int, int]:
    """Racha mas larga de indices consecutivos con ok y el mismo d*. Devuelve (inicio, longitud, d*); (0, 0, 0) si no hay."""
    best = (0, 0, 0)
    i = 0
    while i < ok.size:
        if not ok[i]:
            i += 1
            continue
        j = i
        while j + 1 < ok.size and ok[j + 1] and dst[j + 1] == dst[i]:
            j += 1
        if j - i + 1 > best[1]:
            best = (i, j - i + 1, int(dst[i]))
        i = j + 1
    return best


def spans_decade(length: int) -> bool:
    return length >= 1 and (length - 1) * DECADE_STEP >= 1.0 - 1e-9


def l3_pass(ok: np.ndarray, dst: np.ndarray) -> bool:
    return spans_decade(longest_run(ok, dst)[1])


def l4_rows(F: np.ndarray, d: np.ndarray, degen: np.ndarray, pooled: np.ndarray) -> np.ndarray:
    """Cada instancia degenerada tiene F > media de F de las geometrias con d = d*(fila). `d` etiqueta las columnas geometricas."""
    out = np.zeros(F.shape[0], dtype=bool)
    for i in range(F.shape[0]):
        m = d == pooled[i]
        out[i] = bool(m.any()) and float(degen[i].min()) > float(F[i, m].mean())
    return out


# ------------------------------------------------------------------ carga


def load(out: Path) -> list[dict[str, Any]]:
    feats = {r["id"]: r for r in map(json.loads, (out / "features.jsonl").read_text().splitlines())}
    rows = []
    for t in map(json.loads, (out / "truth.jsonl").read_text().splitlines()):
        rows.append({**t, **{k: feats[t["id"]][k] for k in ("r_star", "max_window", "n_giant", "phi", "psi")}})
    return rows


def rel_close(a: float, b: float, tol: float = NC4_TOL) -> bool:
    mx = max(abs(a), abs(b))
    return True if mx == 0.0 else abs(a - b) / mx <= tol


# ------------------------------------------------------------------ evaluacion


def evaluate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    geo = [r for r in rows if r["group"] == "G" and r["scale"] == "N"]
    sc = [r for r in rows if r["group"] == "G" and r["scale"] == "4N"]
    deg = [r for r in rows if r["group"] == "X"]
    dg = np.array([r["d"] for r in geo], dtype=int)
    sg = np.array([r["stratum"] for r in geo])
    phi_g, psi_g = (np.array([r[k] for r in geo]) for k in ("phi", "psi"))
    phi_x, psi_x = (np.array([r[k] for r in deg]) for k in ("phi", "psi"))
    phi_s, psi_s = (np.array([r[k] for r in sc]) for k in ("phi", "psi"))
    d_s = np.array([r["d"] for r in sc], dtype=int)
    ret_s = np.array([r["stratum"] == "reticulo" for r in sc])
    deg_mean = np.array([r["mean_degree"] for r in geo])

    def F_of(p: float, phi: np.ndarray, psi: np.ndarray) -> np.ndarray:  # (T, n)
        return phi[None, :] + THETAS[:, None] * (psi[None, :] ** p)

    # NC-4 (independiente de p, theta)
    pairs = [(a, b) for i, a in enumerate(geo) for b in geo[i + 1 :]
             if a["d"] != b["d"] and rel_close(a["phi"], b["phi"]) and rel_close(a["psi"], b["psi"])]
    nc4 = {"n_pairs": len(pairs), "pairs": [[a["family"], a["seed"], a["d"], b["family"], b["seed"], b["d"]] for a, b in pairs[:50]]}

    per_p: dict[str, Any] = {}
    grids: dict[str, Any] = {}
    fam_f: dict[str, Any] = {}
    for p in P_VALUES:
        Fg, Fx, Fs = F_of(p, phi_g, psi_g), F_of(p, phi_x, psi_x), F_of(p, phi_s, psi_s)
        pooled, l1, l2 = l1_l2(Fg, dg, sg)
        ok12 = l1 & l2
        start, length, dval = longest_run(ok12, pooled)
        l3 = spans_decade(length)
        run = np.arange(start, start + length)
        by_stratum = {s: dstar_rows(Fg, dg, sg == s).tolist() for s in STRATA}
        l4_all = l4_rows(Fg, dg, Fx, pooled)  # por theta
        # L7: d* del estrato reticulo a 4N frente a N
        ds_scale = dstar_rows(Fs, d_s, ret_s) if ret_s.any() else np.zeros(len(THETAS), dtype=int)
        ds_n = dstar_rows(Fg, dg, sg == "reticulo")
        l7_all = ds_scale == ds_n
        # NC-2
        rho = np.array([_spearman(Fg[t], deg_mean) for t in range(len(THETAS))])
        max_rho = float(np.abs(rho).max())
        plateau_rho = np.abs(rho[run]) if length else np.array([])
        nc2_flag = bool((plateau_rho >= RHO_FLAG).any()) if length else False
        # NC-5
        inter = pooled[np.isin(pooled, INTERIOR)]
        mode = None
        if inter.size:
            cnt = Counter(inter.tolist())
            mode = min(cnt, key=lambda k: (-cnt[k], k))
        per_p[str(p)] = {"pooled_dstar": pooled.tolist(), "L1_thetas": int(l1.sum()), "L2_thetas": int((l1 & l2).sum()),
                         "plateau": {"start_index": int(start), "n_points": int(length), "decades": (length - 1) * DECADE_STEP if length else 0.0,
                                     "dstar": dval, "theta_lo": float(THETAS[start]) if length else None,
                                     "theta_hi": float(THETAS[start + length - 1]) if length else None},
                         "L1": bool(l1.any()), "L2": bool(ok12.any()), "L3": l3,
                         "L4": bool(l3 and l4_all[run].all()), "L4_thetas_ok": int(l4_all.sum()),
                         "L7": bool(l3 and l7_all[run].all()), "L7_thetas_ok": int(l7_all.sum()),
                         "NC2_flag": nc2_flag, "NC2_plateau_rho": [float(x) for x in rho[run]] if length else [],
                         "NC2_max_abs_rho_grid": max_rho,
                         "NC5_mode_dstar": mode}
        grids[str(p)] = {"theta": THETAS.tolist(), "pooled": pooled.tolist(), "stratum": by_stratum,
                         "L4_ok": l4_all.tolist(), "L7_ok": l7_all.tolist(), "reticulo_4N": ds_scale.tolist(), "rho": rho.tolist()}
        t1 = int(np.argmin(np.abs(np.log10(THETAS))))  # theta = 1
        fam_f[str(p)] = {}
        for fam in sorted({r["family"] for r in rows if r["scale"] == "N"}):
            vals = [float(F_of(p, np.array([r["phi"]]), np.array([r["psi"]]))[t1, 0]) for r in rows if r["family"] == fam]
            fam_f[str(p)][fam] = float(np.mean(vals))
    # NC-5 conjunto
    modes = [per_p[str(p)]["NC5_mode_dstar"] for p in P_VALUES]
    nc5_conf = bool(all(m is not None for m in modes) and modes[0] > modes[1] > modes[2])
    nc5 = {"mode_by_p": dict(zip(map(str, P_VALUES), modes)), "confirmed": nc5_conf,
           "exact_match_4_3_2": bool(modes == [LAW[p] for p in P_VALUES])}

    nulls = null_tests(geo, deg, dg, sg)

    verdict_p: dict[str, Any] = {}
    for p in P_VALUES:
        r = per_p[str(p)]
        crit = {"L1": r["L1"], "L2": r["L2"], "L3": r["L3"], "L4": r["L4"], "L7": r["L7"],
                "L8": bool((not r["NC2_flag"]) and (not nc5_conf)), "L9": True,
                "null_A": nulls["A"]["rate"] <= NULL_MAX, "null_B": nulls["B"][str(p)]["rate"] <= NULL_MAX,
                "null_C": nulls["C"][str(p)]["rate"] <= NULL_MAX}
        verdict_p[str(p)] = {"criteria": crit, "pass": all(crit.values()), "failing": [k for k, v in crit.items() if not v]}
    viable = any(v["pass"] for v in verdict_p.values())
    return {"verdict": "SELECCION-VIABLE" if viable else "SELECCION-NO-VIABLE", "by_p": verdict_p, "details": per_p, "NC4": nc4,
            "NC5": nc5, "nulls": nulls, "grid": grids, "F_theta1_by_family": fam_f,
            "L9_blind_by_construction": True, "n_instances": {"geo_N": len(geo), "geo_4N": len(sc), "degenerate": len(deg)},
            "instances": [{k: r[k] for k in ("id", "family", "group", "d", "stratum", "scale", "N", "seed", "mean_degree", "r_star",
                                              "max_window", "n_giant", "phi", "psi")} for r in rows]}


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return 0.0
    v = float(spearmanr(x, y).statistic)
    return 0.0 if np.isnan(v) else v


# ------------------------------------------------------------------ nulos


def _z(x: np.ndarray) -> np.ndarray:
    s = x.std()
    return (x - x.mean()) / (s if s > 0 else 1.0)


def null_tests(geo: list[dict[str, Any]], deg: list[dict[str, Any]], dg: np.ndarray, sg: np.ndarray) -> dict[str, Any]:
    phi_g, psi_g = (np.array([r[k] for r in geo]) for k in ("phi", "psi"))
    phi_x, psi_x = (np.array([r[k] for r in deg]) for k in ("phi", "psi"))
    # A
    lphi = np.log(np.concatenate([phi_g, phi_x]) + EPS_PHI)
    lpsi = np.log(np.concatenate([psi_g, psi_x]))
    zphi, zpsi = _z(lphi), _z(lpsi)
    ng = phi_g.size
    t = np.array([rng_from_key((MASTER, 99, i)).uniform(0.0, 2.0 * np.pi) for i in range(M_A)])
    a, b = np.cos(t), np.sin(t)
    FA = a[:, None] * zphi[None, :] + b[:, None] * zpsi[None, :]
    pooled, l1, l2 = l1_l2(FA[:, :ng], dg, sg)
    l4 = l4_rows(FA[:, :ng], dg, FA[:, ng:], pooled)
    passed = l1 & l2 & l4
    out: dict[str, Any] = {"A": {"M": M_A, "n_pass": int(passed.sum()), "rate": float(passed.mean()),
                                 "rate_L1": float(l1.mean()), "rate_L1_L2": float((l1 & l2).mean())}, "B": {}, "C": {}}
    # B
    perms = []
    for i in range(N_B):
        rng = rng_from_key((MASTER, 99, 1000 + i))
        d2 = dg.copy()
        for s in STRATA:
            idx = np.flatnonzero(sg == s)
            d2[idx] = dg[rng.permutation(idx)]
        perms.append(d2)
    # C
    reps = []
    for i in range(N_C):
        rng = rng_from_key((MASTER, 99, 2000 + i))
        reps.append((rng.integers(1, 6, size=len(deg)), np.array(STRATA)[rng.integers(0, 3, size=len(deg))]))
    for p in P_VALUES:
        Fg = phi_g[None, :] + THETAS[:, None] * psi_g[None, :] ** p
        Fx = phi_x[None, :] + THETAS[:, None] * psi_x[None, :] ** p
        nb = 0
        for d2 in perms:
            pl, a1, a2 = l1_l2(Fg, d2, sg)
            nb += int(l3_pass(a1 & a2, pl))
        nc = 0
        for d2, s2 in reps:
            pl, a1, a2 = l1_l2(Fx, d2, s2)
            nc += int(l3_pass(a1 & a2, pl))
        out["B"][str(p)] = {"n": N_B, "n_pass": nb, "rate": nb / N_B}
        out["C"][str(p)] = {"n": N_C, "n_pass": nc, "rate": nc / N_C}
    return out


# ------------------------------------------------------------------ main


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--smoke", action="store_true", help="lee results/ldim1_smoke/")
    args = ap.parse_args(argv)
    out = OUT_SMOKE if args.smoke else OUT_FULL
    rows = load(out)
    ev = evaluate(rows)
    summary = {"step": "ldim1_eval", "code_commit": head_commit(), "smoke": args.smoke, **ev}
    _atomic_write(out / "summary.json", json.dumps(summary, indent=2, default=str))
    print(f"VEREDICTO {ev['verdict']}  (instancias {ev['n_instances']})")
    for p in P_VALUES:
        v, d = ev["by_p"][str(p)], ev["details"][str(p)]
        pl = d["plateau"]
        print(f"p={p}: pass={v['pass']} failing={v['failing']} meseta={pl['n_points']} pts d*={pl['dstar']} "
              f"max|rho|={d['NC2_max_abs_rho_grid']:.2f} moda_d*={d['NC5_mode_dstar']}")
    print(f"NC-5 {ev['NC5']}  NC-4 pares={ev['NC4']['n_pairs']}")
    n = ev["nulls"]
    print(f"nulo A {n['A']['rate']:.3f} | B {[n['B'][str(p)]['rate'] for p in P_VALUES]} | C {[n['C'][str(p)]['rate'] for p in P_VALUES]}")
    print("d* combinado por theta:")
    for p in P_VALUES:
        print(f"  p={p}: {''.join(str(x) for x in ev['grid'][str(p)]['pooled'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
