"""Piloto 2 DECLARADO (Consejo Rev. 4 §8): valida R-rodilla-2 solo en sprinklings 2D. Entropia 20261051. No es evidencia."""

from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omega.causality.entropy import interval_balanced, knee_powerlaw, pauli_jordan, sj_both, ssee  # noqa: E402
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond  # noqa: E402

ENT = 20261051
NS = (800, 1200, 1600)
SEEDS = 3
BAND = (0.5, 2.0)
OUT = Path("results") / "omega_e" / "piloto2_e1p.json"


def main() -> None:
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip():
        sys.exit("exige arbol git limpio")
    t0 = time.time()
    rows: list[dict[str, Any]] = []
    for i, n in enumerate(NS):
        for s in range(SEEDS):
            rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ENT, spawn_key=(i, s))))
            c = causal_matrix(sprinkle_diamond(n, 2, rng))
            r = interval_balanced(c, n // 4)
            idelta = pauli_jordan(c)
            fixed_big = math.sqrt(n) / (4 * math.pi)
            fixed_r = math.sqrt(r.size) / (4 * math.pi)
            kb, beta_b = knee_powerlaw(idelta, n)
            row: dict[str, Any] = {"N": n, "n_R": int(r.size), "rodilla_grande": kb, "beta_grande": beta_b,
                                   "fijo_grande": fixed_big, "fijo_R": fixed_r}
            _, idt_f, w_f = sj_both(idelta, fixed_big)
            row["S_fijo"] = ssee(idt_f, w_f, r, True)
            if kb is not None:
                _, idt_k, w_k = sj_both(idelta, kb)
                kr, beta_r = knee_powerlaw(idt_k[np.ix_(r, r)], int(r.size))
                row.update({"rodilla_R": kr, "beta_R": beta_r})
                row["S_rodilla"] = ssee(idt_k, w_k, r, True, kr if kr is not None else 0.0)
            else:
                row.update({"rodilla_R": None, "beta_R": float("nan"), "S_rodilla": None})
            ok_b = kb is not None and BAND[0] * fixed_big <= kb <= BAND[1] * fixed_big
            ok_r = row["rodilla_R"] is not None and BAND[0] * fixed_r <= row["rodilla_R"] <= BAND[1] * fixed_r
            row["ok"] = bool(ok_b and ok_r)
            rows.append(row)
    n_ok = sum(r["ok"] for r in rows)
    verdict = "ACEPTA" if n_ok >= 8 else "RECHAZA"
    OUT.write_text(json.dumps({"commit": head, "declarado": "piloto 2, no evidencia", "seconds": round(time.time() - t0, 1),
                               "n_ok": n_ok, "veredicto": verdict, "rows": rows}, indent=1, default=str))
    print(verdict, n_ok, "de", len(rows), round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
