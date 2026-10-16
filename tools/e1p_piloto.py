"""Paso 3 (Consejo Rev. 4): piloto DECLARADO de E1' solo con sprinklings 2D. No es evidencia. Entropia 20261050."""

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

from omega.causality.entropy import (  # noqa: E402
    causal_complement,
    interval_balanced,
    knee_cutoff,
    pauli_jordan,
    sj_both,
    ssee,
)
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond  # noqa: E402

ENT = 20261050
NS = (800, 1200, 1600)
SEEDS = 3
OUT = Path("results") / "omega_e" / "piloto_e1p.json"


def main() -> None:
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip():
        sys.exit("exige arbol git limpio")
    t0 = time.time()
    rows: list[dict[str, Any]] = []
    for i, n in enumerate(NS):
        for s in range(SEEDS):
            rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ENT, spawn_key=(i, s))))
            x = sprinkle_diamond(n, 2, rng)
            c = causal_matrix(x)
            r = interval_balanced(c, n // 4)
            rp = causal_complement(c, r)
            a, b = int(r[0]), int(r[-1])
            past_a = int(c[:, r].sum(axis=0).min())
            fut_b = int(c[r, :].sum(axis=1).min())
            idelta = pauli_jordan(c)
            fixed = math.sqrt(n) / (4 * math.pi)
            w_full, idt, w_t = sj_both(idelta, fixed)
            row: dict[str, Any] = {
                "N": n, "n_R": int(r.size), "n_Rperp": int(rp.size),
                "frac_pasado_min": past_a / n, "frac_futuro_min": fut_b / n,
                "S_R_tr": ssee(idt, w_t, r, True), "S_Rperp_tr": ssee(idt, w_t, rp, True),
                "S_R_un": ssee(idelta, w_full, r, False), "S_Rperp_un": ssee(idelta, w_full, rp, False),
                "umbral_fijo": fixed,
            }
            k = knee_cutoff(idelta)
            row["rodilla"] = k
            if k is not None:
                _, idt_k, w_k = sj_both(idelta, k)
                kr = knee_cutoff(idt_k[np.ix_(r, r)])
                row["rodilla_R"] = kr
                row["S_R_rodilla"] = ssee(idt_k, w_k, r, True, kr if kr is not None else 0.0)
            rows.append(row)
            _ = (a, b)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"commit": head, "declarado": "piloto, no evidencia", "seconds": round(time.time() - t0, 1), "rows": rows}, indent=1))
    print("ok", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
