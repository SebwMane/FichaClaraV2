"""Ejecuta los pasos 1–5 del panel (P§22) en modo full, en orden y sobre un commit limpio.

Pasos: p00 (analítica + golden) → O-00 → O-01 → O-02 → O-03 → O-04. No ejecuta s06 ni los
ensembles (decisión del usuario: congelar y analizar 1–5 primero). Usa exactamente las mismas
configuraciones que los `test_full_run` de cada experimento. La compuerta (`gate.py`) exige
árbol limpio y code_commit == HEAD; este script no modifica ningún parámetro.

Uso: python tools/run_steps_1_5.py [--from STEP]
Salida: $OMEGA_RUNS_DIR (por defecto runs/omega11) y un log JSON de tiempos en <out_root>/steps_1_5_log.json.
"""

from __future__ import annotations

import argparse
import importlib
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omega.config.settings11 import Omega11Config  # noqa: E402
from omega.experiments.v11.gate import runs_root  # noqa: E402
from omega.phases.scan import default_config  # noqa: E402

ENTROPY = 20240901

# (paso, módulo, N, réplicas) — idénticos a los test_full_run de cada módulo.
STEPS: tuple[tuple[str, str, int, int], ...] = (
    ("p00", "omega.experiments.v11.test_p00_prereq", 200, 3),
    ("o00", "omega.experiments.v11.test_o00_validation_v11", 200, 3),
    ("o01", "omega.experiments.v11.test_o01_baseline", 200, 10),
    ("o02", "omega.experiments.v11.test_o02_distance", 200, 5),
    ("o03", "omega.experiments.v11.test_o03_nulls", 800, 5),
    ("o04", "omega.experiments.v11.test_o04_ablation", 100, 10),
)


def main() -> int:
    """Corre los pasos en orden; se detiene en el primer fallo o resumen incompleto."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--from", dest="start", default="p00", choices=[s[0] for s in STEPS])
    args = parser.parse_args()
    out_root = runs_root()
    out_root.mkdir(parents=True, exist_ok=True)
    log_path = out_root / "steps_1_5_log.json"
    log: list[dict[str, object]] = json.loads(log_path.read_text()) if log_path.exists() else []
    names = [s[0] for s in STEPS]
    for step, module, n, reps in STEPS[names.index(args.start):]:
        cfg = Omega11Config(base=default_config(n, master_entropy=ENTROPY, replicates=reps))
        mod = importlib.import_module(module)
        t0 = time.time()
        print(f"[{time.strftime('%H:%M:%S')}] inicio {step} ({module}, N={n}, reps={reps})", flush=True)
        path = mod.run(cfg, out_root, mode="full")
        dt = time.time() - t0
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        entry = {"step": step, "seconds": round(dt, 1), "summary": str(path),
                 "complete": data.get("complete"), "expectations_met": data.get("expectations_met")}
        log.append(entry)
        log_path.write_text(json.dumps(log, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"[{time.strftime('%H:%M:%S')}] fin {step}: {entry}", flush=True)
        if data.get("complete") is not True:
            print(f"ABORTADO: {step} no está completo", flush=True)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
