"""R14-A: auditoría de inferencia y trazabilidad de F6 (docs/OMEGA_CONSEJO_REV14_EMERGENCIA.md, R14-2).

Solo lee resultados ya publicados (results/omega_f6_intento2, results/omega_f6_confirmacion); no simula nada y
no cambia ningún veredicto. Salida: results/omega_r14a/summary.json.

Análisis:
  A1 unidad independiente: varianza entre realizaciones de r frente a la varianza «ingenua» dentro de una realización
     (r_sd²/n_sel), por celda (m, N), en los 10 grupos planos de calibración.
  A2 covarianza entre escalas: correlación entre realizaciones de Δr en distintos m del mismo N (planos).
  A3 nula de K: los 15 valores de K de calibración frente a χ²₄ (lo que darían pares independientes y normales).
  A4 estabilidad de los umbrales: bootstrap de los 15 grupos (regla 1.5 × max) y frecuencia con que cambia la
     clasificación de cada estructura de prueba.
  A5 sensibilidad a la regla (descriptivo, NO rescata): errores con 1.0 × max, 1.5 × max y 2.0 × max.
  A6 forma del cono frente a b: ε y b de los conos con voto.
"""

from __future__ import annotations

import glob
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
from scipy import stats

I2 = Path("results") / "omega_f6_intento2"
OUT = Path("results") / "omega_r14a"
MS = (10, 20, 40, 80, 160)
NS = (4000, 8000)
SEED = 20261087  # entropía nueva, solo para el bootstrap de grupos (no genera estructuras)
EPS = {"cono_cubico_4d": 1 - 1 / math.sqrt(3), "cono_octaedrico_4d": math.sqrt(3) - 1, "cono_p6_4d": 1 - 3 ** (1 / 6 - 0.5),
       "C3_cono_cilindrico": 1 - 1 / math.sqrt(2), "C4_cono_esfera_cubo": 1 / 0.79 - 1, "C6_cono_p_por_eje": 0.3144,
       "D6_cono_tapa_plana": 0.25, "D7_cono_huso": 1 - 1 / 1.45, "D8_cono_suma_normas": 0.32 * (math.sqrt(3) - 1),
       "cono_p1_5_4d": 3 ** (1 / 6) - 1, "cono_p3_4d": 1 - 3 ** (-1 / 6), "C7_cono_anidado": 1 - 2 ** (-0.25)}


def realizations(base: str, n: int) -> list[dict[str, Any]]:
    rx = re.compile(rf"^{re.escape(base)}_{n}(_r\d+-\d+)?\.json$")
    rows: list[dict[str, Any]] = []
    for f in sorted(glob.glob(str(I2 / f"{base}_{n}*.json"))):
        if rx.match(Path(f).name):
            rows += json.loads(Path(f).read_text())["rows"]
    return sorted(rows, key=lambda r: r["realizacion"])


def a1_a2() -> dict[str, Any]:
    ratios: dict[str, list[float]] = {}
    corr_rows: dict[int, list[list[float]]] = {n: [] for n in NS}
    for g in range(10):
        for n in NS:
            rows = realizations(f"cal_{g}", n)
            mat = np.full((len(rows), len(MS)), np.nan)
            for i, r in enumerate(rows):
                for j, m in enumerate(MS):
                    c = r["por_m"][str(m)]
                    if c.get("medible"):
                        mat[i, j] = c["r_media"]
            for j, m in enumerate(MS):
                col = mat[:, j]
                ok = np.isfinite(col)
                if ok.sum() < 5:
                    continue
                between = float(np.var(col[ok], ddof=1))
                naive = float(np.mean([r["por_m"][str(m)]["r_sd"] ** 2 / r["por_m"][str(m)]["n_sel"]
                                       for r in rows if r["por_m"][str(m)].get("medible")]))
                ratios.setdefault(f"{n},{m}", []).append(between / naive)
            corr_rows[n].extend((mat - np.nanmean(mat, axis=0)).tolist())  # centrado por grupo
    a1 = {k: {"mediana": float(np.median(v)), "min": float(min(v)), "max": float(max(v))} for k, v in ratios.items()}
    a2: dict[str, Any] = {}
    for n in NS:
        x = np.array(corr_rows[n])
        cols = [j for j in range(len(MS)) if np.isfinite(x[:, j]).mean() > 0.9]
        sub = x[:, cols]
        sub = sub[np.all(np.isfinite(sub), axis=1)]
        a2[str(n)] = {"m": [MS[j] for j in cols], "corr": np.round(np.corrcoef(sub.T), 3).tolist(), "realizaciones": int(sub.shape[0])}
    return {"A1_varianza_entre_sobre_ingenua": a1, "A2_correlacion_entre_escalas": a2}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    s = json.loads((I2 / "summary.json").read_text())
    E = s["estructuras"]
    cal = {k: v for k, v in E.items() if v["conjunto"] == "calibracion"}
    test = {k: v for k, v in E.items() if v["conjunto"] != "calibracion"}

    def bk(v: dict[str, Any]) -> tuple[float, float]:
        return float(v["ajuste"]["b"]["valor"]), float(v["colapso"]["K"])

    flat = [bk(v)[0] for k, v in cal.items() if "diamante_plano" in k]
    cb = {k: bk(v) for k, v in cal.items()}
    b_plano = float(np.median(flat))
    out: dict[str, Any] = {"fuente": str(I2), "commit_partes": s["commit_partes"], "b_plano": b_plano,
                           "tau_b_oficial": s["calibracion"]["tau_b"], "tau_K_oficial": s["calibracion"]["tau_K"]}
    out.update(a1_a2())

    kc = np.array([v[1] for v in cb.values()])
    out["A3_nula_K"] = {"K_calibracion": np.round(np.sort(kc), 2).tolist(), "media": float(kc.mean()),
                        "media_chi2_4": 4.0, "p95_chi2_4": float(stats.chi2.ppf(0.95, 4)),
                        "max_esperado_15_chi2_4_mediana": float(stats.chi2.ppf(0.5 ** (1 / 15), 4)),
                        "ks_vs_chi2_4_p": float(stats.kstest(kc, "chi2", args=(4,)).pvalue)}

    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(SEED)))
    names = list(cb)
    flips = {k: 0 for k in test}
    taus = []
    nb = 4000
    for _ in range(nb):
        idx = rng.integers(0, len(names), len(names))
        sel = [cb[names[i]] for i in idx]
        fl = [b for (b, _), i in zip(sel, idx) if "diamante_plano" in names[i]]
        bp = float(np.median(fl)) if fl else b_plano
        tb = 1.5 * max(abs(b - bp) for b, _ in sel)
        tk = 1.5 * max(K for _, K in sel)
        taus.append((tb, tk))
        for k, v in test.items():
            b, K = bk(v)
            et = abs(b - bp) <= tb and K <= tk
            if et != v["clasificacion"]["espacio_tiempo"]:
                flips[k] += 1
    t = np.array(taus)
    out["A4_bootstrap_umbrales"] = {"replicas": nb, "tau_b_p5_p50_p95": np.round(np.percentile(t[:, 0], [5, 50, 95]), 5).tolist(),
                                    "tau_K_p5_p50_p95": np.round(np.percentile(t[:, 1], [5, 50, 95]), 2).tolist(),
                                    "frac_cambia_clasificacion": {k: round(c / nb, 3) for k, c in flips.items() if c}}

    sens: dict[str, Any] = {}
    devs = [abs(b - b_plano) for b, _ in cb.values()]
    for f in (1.0, 1.5, 2.0):
        tb, tk = f * max(devs), f * float(kc.max())
        err = []
        for k, v in test.items():
            if v["papel"] not in ("positiva", "negativa"):
                continue
            b, K = bk(v)
            et = abs(b - b_plano) <= tb and K <= tk
            if et != (v["papel"] == "positiva"):
                err.append(k)
        sens[f"{f}x_max"] = {"tau_b": tb, "tau_K": tk, "errores": err}
    out["A5_sensibilidad_regla_descriptiva"] = sens

    out["A6_cono_eps_frente_a_b"] = sorted(([k, round(EPS[k], 3), round(bk(test[k])[0] - b_plano, 5), test[k]["papel"]]
                                            for k in EPS if k in test), key=lambda r: r[1])
    (OUT / "summary.json").write_text(json.dumps(out, indent=2))
    print(json.dumps({k: out[k] for k in ("A3_nula_K", "A4_bootstrap_umbrales", "A5_sensibilidad_regla_descriptiva")}, indent=1)[:3000])


if __name__ == "__main__":
    main()
