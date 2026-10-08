"""Omega L-ARQ-0: comprobacion del lema del dial (L-ARQ-T1) sobre datos cerrados de L-DIM-1 (results/ldim1/summary.json).

No es un experimento nuevo: verifica una consecuencia de un teorema sobre datos ya publicados. Para cada p, toma los puntos
(B_d, A_d) = (media de psi^p, media de phi) de las geometrias de tamano N por dimension, calcula la envolvente convexa inferior
y las ventanas de theta en que cada d es argmin de A + theta*B, y las compara con la escalera observada.
Uso: python tools/arq_dial_check.py   -> results/arq/dial_check.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))

from c0_dynamics import _atomic_write  # noqa: E402


def lower_hull_windows(pts: dict[int, tuple[float, float]]) -> dict[int, tuple[float, float]]:
    """Ventanas exactas [theta_lo, theta_hi] en que d minimiza A + theta*B (theta > 0). Solo vertices de la envolvente inferior."""
    order = sorted(pts, key=lambda d: (-pts[d][0], pts[d][1]))  # B decreciente: theta pequeno favorece A minimo
    out: dict[int, tuple[float, float]] = {}
    for d in order:
        b, a = pts[d]
        lo, hi = 0.0, np.inf
        ok = True
        for e in pts:
            if e == d:
                continue
            be, ae = pts[e]
            # A_d + t B_d <= A_e + t B_e  <=>  t (B_d - B_e) <= A_e - A_d
            db, da = b - be, ae - a
            if db > 0:
                hi = min(hi, da / db)
            elif db < 0:
                lo = max(lo, da / db)
            elif da < 0:
                ok = False
        if ok and lo < hi:
            out[d] = (float(lo), float(hi))
    return out


def main() -> int:
    s = json.loads((ROOT / "results" / "ldim1" / "summary.json").read_text())
    geo = [r for r in s["instances"] if r["group"] == "G" and r["scale"] == "N"]
    res = {}
    for p in (0.5, 1.0, 2.0):
        pts = {d: (float(np.mean([r["psi"] ** p for r in geo if r["d"] == d])), float(np.mean([r["phi"] for r in geo if r["d"] == d])))
               for d in range(1, 6)}
        win = lower_hull_windows(pts)
        observed = s["grid"][str(p)]["pooled"]
        thetas = s["grid"][str(p)]["theta"]
        predicted = [min(pts, key=lambda d: pts[d][1] + t * pts[d][0]) for t in thetas]
        res[str(p)] = {"points_B_A": pts, "windows": win, "selectable": sorted(win),
                       "decades": {d: (float(np.log10(h / lo)) if lo > 0 and np.isfinite(h) else None) for d, (lo, h) in win.items()},
                       "staircase_matches_observed": predicted == observed}
        print(p, {d: f"[{lo:.3g}, {h:.3g}]" for d, (lo, h) in win.items()}, "coincide:", predicted == observed)
    _atomic_write(ROOT / "results" / "arq" / "dial_check.json", json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
