"""Bloque A-S (Rev. 3 §17): familia sandwich (capas a, 2, b) en la rama lorentziana. Constantes congeladas.

Tipos de un elemento externo: indice 0 = (0,0), 1 = (0,1), 2 = (1,0), 3 = (1,1) segun sus relaciones con m1, m2.
"""

from __future__ import annotations

import argparse
import cmath
import itertools
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import mpmath  # type: ignore[import-untyped]
import numpy as np
from scipy.special import gammaln

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omega.causality.estimators import bd_action_2d_smeared  # noqa: E402
from omega.causality.orders import transitive_closure  # noqa: E402

EPS_AS2 = (0.5, 0.25, 0.1)
X_AS2 = (0.3, 1.2)
N_AS3 = (64, 256, 1024)
EPS_AS3 = (0.25, 0.1)
X_AS3 = (0.3, 1.2, math.pi / 2, 1.9)
X_WINDOW = (1.2, math.pi / 2, 1.9)
N_EXACT = (8, 12, 16, 20, 24)
N_TREND = (64, 128, 256, 512, 1024)
OUT = Path("results") / "bloque_as"
TYPES = ((0, 0), (0, 1), (1, 0), (1, 1))


def git_state() -> tuple[str, bool]:
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip())
    return head, dirty


def f12(eps: float) -> tuple[float, float]:
    return 1 - 3 * eps, 1 - 6 * eps + 6 * eps * eps


def action_from_types(n: int, eps: float, ta: list[int], tb: list[int]) -> float:
    f1, f2 = f12(eps)
    links = sum(sum(TYPES[t]) for t in ta) + sum(sum(TYPES[u]) for u in tb)
    n1 = n2 = 0
    for t in ta:
        for u in tb:
            k = TYPES[t][0] * TYPES[u][0] + TYPES[t][1] * TYPES[u][1]
            n1 += k == 1
            n2 += k == 2
    return eps * (n - 2 * eps * (links + f1 * n1 + f2 * n2))


def sandwich_order(a: int, b: int, ta: list[int], tb: list[int]) -> np.ndarray[Any, Any]:
    n = a + b + 2
    m1, m2 = a, a + 1
    adj = np.zeros((n, n), dtype=bool)
    for i, t in enumerate(ta):
        adj[i, m1], adj[i, m2] = bool(TYPES[t][0]), bool(TYPES[t][1])
    for j, u in enumerate(tb):
        adj[m1, a + 2 + j], adj[m2, a + 2 + j] = bool(TYPES[u][0]), bool(TYPES[u][1])
    return transitive_closure(adj)


def g_terms(theta: float, f1: float, f2: float, al: np.ndarray[Any, Any]) -> np.ndarray[Any, Any]:
    """G(alpha) para una matriz de conteos al[:, t]; theta = -2 lam eps^2."""
    a01, a10, a11 = al[:, 1], al[:, 2], al[:, 3]
    out: np.ndarray[Any, Any] = (1 + np.exp(1j * theta * (1 + f1 * (a01 + a11))) + np.exp(1j * theta * (1 + f1 * (a10 + a11)))
            + np.exp(1j * theta * (2 + f1 * (a01 + a10) + f2 * a11)))
    return out


def z_types_small(n: int, eps: float, x: float) -> complex:
    """Z por la formula de tipos (precision doble, solo N chico), sin factor de etiquetas de capa."""
    a = b = (n - 2) // 2
    lam = x / eps**2
    theta = -2 * x
    f1, f2 = f12(eps)
    total = 0j
    for al in itertools.product(range(a + 1), repeat=3):
        if sum(al) > a:
            continue
        alpha = np.array([[a - sum(al), *al]], dtype=float)
        mult = math.factorial(a) / math.prod(math.factorial(int(v)) for v in alpha[0])
        g = complex(g_terms(theta, f1, f2, alpha)[0])
        total += mult * cmath.exp(1j * theta * (alpha[0, 1] + alpha[0, 2] + 2 * alpha[0, 3])) * g**b
    return complex(cmath.exp(1j * lam * eps * n) * total)


def as2() -> dict[str, Any]:
    rows = []
    for n in (6, 8):
        a = b = (n - 2) // 2
        configs = list(itertools.product(range(4), repeat=a + b))
        for eps in EPS_AS2:
            maxdiff = 0.0
            acts = []
            for cfg in configs:
                ta, tb = list(cfg[:a]), list(cfg[a:])
                s_gen = bd_action_2d_smeared(sandwich_order(a, b, ta, tb), eps)
                s_typ = action_from_types(n, eps, ta, tb)
                maxdiff = max(maxdiff, abs(s_gen - s_typ))
                acts.append(s_gen)
            for x in X_AS2:
                lam = x / eps**2
                z_enum = complex(np.sum(np.exp(1j * lam * np.array(acts))))
                z_typ = z_types_small(n, eps, x)
                rows.append({"n": n, "eps": eps, "x": x, "max_action_diff": maxdiff, "z_diff": abs(z_enum - z_typ),
                             "ok": maxdiff < 1e-9 and abs(z_enum - z_typ) < 1e-9 * len(configs)})
    return {"rows": rows, "pass": all(r["ok"] for r in rows)}


def layer_factor(n: int) -> float:
    a = (n - 2) // 2
    return math.lgamma(n + 1) - 2 * math.lgamma(a + 1) - math.log(2)


def upper_bound(n: int, eps: float, x: float) -> float:
    a = b = (n - 2) // 2
    theta = -2 * x
    f1, f2 = f12(eps)
    lga = math.lgamma(a + 1)
    best = -math.inf
    acc = 0.0
    for a11 in range(a + 1):
        for a10 in range(a + 1 - a11):
            a01 = np.arange(a + 1 - a11 - a10)
            a00 = a - a11 - a10 - a01
            al = np.stack([a00, a01, np.full_like(a01, a10), np.full_like(a01, a11)], axis=1).astype(float)
            logm = lga - np.sum(gammaln(al + 1), axis=1)
            val = logm + b * np.log(np.abs(g_terms(theta, f1, f2, al)) + 1e-300)
            m = float(val.max())
            if m > best:
                acc = acc * math.exp(best - m) if best > -math.inf else 0.0
                best = m
            acc += float(np.sum(np.exp(val - best)))
    return layer_factor(n) + best + math.log(acc)


def exact_logabs(n: int, eps: float, x: float) -> float:
    mpmath.mp.dps = 80
    a = b = (n - 2) // 2
    theta = mpmath.mpf(-2 * x)
    f1, f2 = (mpmath.mpf(v) for v in f12(eps))
    total = mpmath.mpc(0)
    for al in itertools.product(range(a + 1), repeat=3):
        if sum(al) > a:
            continue
        a01, a10, a11 = al
        mult = mpmath.factorial(a) / (mpmath.factorial(a - sum(al)) * mpmath.factorial(a01) * mpmath.factorial(a10) * mpmath.factorial(a11))
        g = (1 + mpmath.expj(theta * (1 + f1 * (a01 + a11))) + mpmath.expj(theta * (1 + f1 * (a10 + a11)))
             + mpmath.expj(theta * (2 + f1 * (a01 + a10) + f2 * a11)))
        total += mult * mpmath.expj(theta * (a01 + a10 + 2 * a11)) * g**b
    return float(mpmath.log(abs(total))) + layer_factor(n)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["as2", "as3", "final"], required=True)
    args = ap.parse_args()
    head, dirty = git_state()
    if dirty:
        sys.exit("exige arbol git limpio (fuera de results/)")
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if args.part == "as2":
        data: dict[str, Any] = as2()
    elif args.part == "as3":
        data = {
            "U": {f"N{n}_eps{e}_x{x:.4f}": upper_bound(n, e, x) for n in N_TREND for e in EPS_AS3 for x in X_AS3},
            "exact_small": {f"N{n}_eps{e}_x{x:.4f}": {"exact": exact_logabs(n, e, x), "U": upper_bound(n, e, x)}
                            for n in N_EXACT for e in EPS_AS3 for x in X_AS3},
        }
    else:
        fits = json.loads((Path("results") / "bloque_a" / "a3.json").read_text())["data"]["exponent"]
        u = json.loads((OUT / "as3.json").read_text())["data"]["U"]
        as2_pass = json.loads((OUT / "as2.json").read_text())["data"]["pass"]
        geo = {}
        comp = {}
        for e in EPS_AS3:
            f = fits[str(e)]
            for x in X_WINDOW:
                for n in N_TREND:
                    g = math.lgamma(n + 1) - x * x * math.exp(f["c"]) * n ** f["gamma"] / 2
                    key = f"N{n}_eps{e}_x{x:.4f}"
                    geo[key] = g
                    comp[key] = {"U": u[key], "geo": g, "sandwich_suprimido": u[key] < g}
        win = [x for x in X_WINDOW if all(comp[f"N1024_eps{e}_x{x:.4f}"]["sandwich_suprimido"] for e in EPS_AS3)]
        data = {"AS-2": as2_pass, "comparacion": comp,
                "AS-4": ("SOBREVIVE (condicional)" if win and as2_pass else ("NO CONCLUYENTE" if as2_pass else "INDETERMINADO (AS-2 falla)")),
                "x_en_ventana_que_suprimen": win}
        (OUT / "summary.json").write_text(json.dumps({"commit": head, **data}, indent=2))
        print(json.dumps({k: v for k, v in data.items() if k != "comparacion"}, indent=2))
        return
    (OUT / f"{args.part}.json").write_text(json.dumps({"commit": head, "seconds": round(time.time() - t0, 1), "data": data}))
    print(args.part, "ok", round(time.time() - t0, 1), "s", "" if args.part != "as2" else data["pass"])


if __name__ == "__main__":
    main()
