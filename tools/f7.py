"""F7 (docs/OMEGA_CONSEJO_REV10_F6F7.md §8). Partes: ref_<d> | cal_<d> | cal_<d>_<N> | test_<estructura> | final.

Igual que F5′ (``tools/f5prima.py``) salvo lo que fija §8: EXPERIMENT_ID 319, entropías nuevas, 17
estructuras de prueba (§8.1) y la distancia de forma D_forma = max_k |v_k − media_k(v)| como criterio
que vota (D de F5′ queda informativa).

- ``ref_<d>``: 20 diamantes planos por N; referencia = media de N_k, k = 0..K_MAX.
- ``cal_<d>``: 40 diamantes planos por N (semillas 0..9 con espectro, d_MM y d_mid) y, para d ∈ {3, 4},
  20 FRW (a ∝ η²) por N, solo abundancias. ``cal_<d>_<N>`` es la misma parte restringida a un N.
- ``test_<estructura>``: 5 semillas por N de cada una de las 17 estructuras de §8.1.
- ``final``: umbrales δ_forma y δ, clasificación y veredicto único PASA/FALLA (§8).
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from f1f2 import random_dorder  # noqa: E402
from f3f4 import windows  # noqa: E402
from f5 import d_spec, measure  # noqa: E402
from f5prima import F3_WIN, Gen, polyhedral_cone  # noqa: E402

from omega.causality.estimators import interval_abundances  # noqa: E402
from omega.causality.orders import kr_order, thick_layer_order, transitive_percolation  # noqa: E402
from omega.causality.profile import K_MAX, profile_distance  # noqa: E402
from omega.causality.sprinkling import DIAMOND_HALF_HEIGHT, causal_matrix, sprinkle_diamond  # noqa: E402

EXPERIMENT_ID = 319
ENT_REF, ENT_CAL, ENT_TEST = 20261072, 20261073, 20261074
NS = (1600, 2400)
DIMS = (2, 3, 4, 5)
FRW_DIMS = (3, 4)
N_REF, N_CAL, N_CAL_SPEC, N_FRW, N_TEST = 20, 40, 10, 20, 5
FAM_PLANA, FAM_FRW = 0, 1
DOC = Path("docs") / "OMEGA_CONSEJO_REV10_F6F7.md"
OUT = Path("results") / "omega_f7"
Bool = npt.NDArray[np.bool_]
Float = npt.NDArray[np.float64]


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


# --- Generadores -----------------------------------------------------------------------------------


def de_sitter_eta_points(n: int, d: int, r: np.random.Generator, eta_lo: float = -1.5, eta_hi: float = -0.5) -> Float:
    """Diamante de `sprinkle_diamond` con η = t + (η_lo + η_hi)/2 ∈ [η_lo, η_hi], aceptación (|η_hi|/|η|)^d.

    El intervalo debe medir la altura del diamante (1) y quedar en η < 0. Con los valores por defecto
    (η = t − 1, aceptación (0.5/|η|)^d) consume el generador igual que `f5prima.de_sitter_points`.
    """
    if not (eta_lo < eta_hi < 0) or not math.isclose(eta_hi - eta_lo, 2 * DIAMOND_HALF_HEIGHT):
        raise ValueError("se requiere η_lo < η_hi < 0 y η_hi − η_lo = 1")
    shift = (eta_lo + eta_hi) / 2
    a = abs(eta_hi)
    pts: list[Float] = []
    total = 0
    while total < n:
        cand = sprinkle_diamond(4 * n, d, r)
        eta = cand[:, 0] + shift
        keep = cand[r.random(cand.shape[0]) < (a / np.abs(eta)) ** d]
        pts.append(keep)
        total += keep.shape[0]
    out: Float = np.concatenate(pts)[:n]
    return out


def de_sitter_eta(d: int, eta_lo: float = -1.5, eta_hi: float = -0.5) -> Gen:
    """Generador (n, rng) -> matriz causal de de Sitter con η ∈ [η_lo, η_hi] (relación de Minkowski)."""

    def gen(n: int, r: np.random.Generator) -> Bool:
        return causal_matrix(de_sitter_eta_points(n, d, r, eta_lo, eta_hi))

    return gen


def frw_power_points(n: int, d: int, power: float, r: np.random.Generator) -> Float:
    """FRW conforme con η = t + 1 ∈ [0.5, 1.5], aceptación (η/1.5)^power.

    power = 2d es a ∝ η² (`f5prima.frw_points`, misma secuencia aleatoria); power = d es radiación (a ∝ η).
    """
    pts: list[Float] = []
    total = 0
    while total < n:
        cand = sprinkle_diamond(4 * n, d, r)
        eta = cand[:, 0] + 1.0
        keep = cand[r.random(cand.shape[0]) < (eta / 1.5) ** power]
        pts.append(keep)
        total += keep.shape[0]
    out: Float = np.concatenate(pts)[:n]
    return out


def frw_radiation(d: int) -> Gen:
    """FRW de radiación (a ∝ η): densidad ∝ η^d."""

    def gen(n: int, r: np.random.Generator) -> Bool:
        return causal_matrix(frw_power_points(n, d, float(d), r))

    return gen


def _check_p(p: float) -> None:
    if not (math.isfinite(p) and p >= 1):
        raise ValueError("p debe ser finito y >= 1 (p = ∞ y p = 1 usan polyhedral_cone)")


def lp_norm(x: Float, p: float) -> Float:
    """Norma L^p por filas."""
    _check_p(p)
    out: Float = np.sum(np.abs(x) ** p, axis=1) ** (1.0 / p)
    return out


def sprinkle_lp(n: int, p: float, rng: np.random.Generator) -> Float:
    """n puntos uniformes en {|t| <= T, ‖x‖_p <= T − |t|} (1 + 3 dimensiones); rechazo en [−T, T]^4."""
    _check_p(p)
    if n < 0:
        raise ValueError("n debe ser >= 0")
    t_half = DIAMOND_HALF_HEIGHT
    out = np.empty((n, 4), dtype=np.float64)
    filled = 0
    while filled < n:
        batch = max(64, 32 * (n - filled))
        cand = rng.uniform(-t_half, t_half, size=(batch, 4))
        ok = cand[lp_norm(cand[:, 1:], p) <= t_half - np.abs(cand[:, 0])]
        take = min(len(ok), n - filled)
        out[filled : filled + take] = ok[:take]
        filled += take
    return out


def lp_matrix(coords: Float, p: float) -> Bool:
    """C[i, j] = True si Δt > 0 y Δt >= ‖Δx‖_p."""
    _check_p(p)
    t = coords[:, 0]
    dt = t[None, :] - t[:, None]
    acc = np.zeros_like(dt)
    for k in range(1, coords.shape[1]):
        acc += np.abs(coords[None, :, k] - coords[:, None, k]) ** p
    rel: Bool = (dt > 0) & (dt >= acc ** (1.0 / p))
    return rel


def lp_cone(p: float) -> Gen:
    _check_p(p)

    def gen(n: int, r: np.random.Generator) -> Bool:
        return lp_matrix(sprinkle_lp(n, p, r), p)

    return gen


# --- Distancias ------------------------------------------------------------------------------------


def log_ratio(nk: npt.NDArray[np.int64] | Float, ref: Float) -> Float:
    """v_k = ln((N_k + 1)/(ref_k + 1))."""
    a = np.asarray(nk, dtype=np.float64)
    if a.shape != ref.shape:
        raise ValueError("perfiles de distinta longitud")
    out: Float = np.log((a + 1) / (ref + 1))
    return out


def shape_distance(nk: npt.NDArray[np.int64] | Float, ref: Float) -> float:
    """D_forma = max_k |v_k − media_k(v)|."""
    v = log_ratio(nk, ref)
    return float(np.max(np.abs(v - np.mean(v))))


# El orden del diccionario es el índice de estructura de §8.1 (spawn_key).
TEST: dict[str, tuple[str, int, Gen]] = {
    "diamante_4d": ("positiva", 4, lambda n, r: causal_matrix(sprinkle_diamond(n, 4, r))),
    "diamante_3d": ("positiva", 3, lambda n, r: causal_matrix(sprinkle_diamond(n, 3, r))),
    "de_sitter_4d": ("positiva", 4, de_sitter_eta(4)),
    "de_sitter_3d": ("positiva", 3, de_sitter_eta(3)),
    "de_sitter_4d_fuerte": ("positiva", 4, de_sitter_eta(4, -1.2, -0.2)),
    "frw_radiacion_4d": ("positiva", 4, frw_radiation(4)),
    "cono_cubico_4d": ("negativa", 0, polyhedral_cone("max")),
    "cono_octaedrico_4d": ("negativa", 0, polyhedral_cone("l1")),
    "cono_p1_5_4d": ("negativa", 0, lp_cone(1.5)),
    "cono_p3_4d": ("negativa", 0, lp_cone(3.0)),
    "cono_p6_4d": ("informativa", 0, lp_cone(6.0)),
    "producto_3": ("negativa", 0, random_dorder(3)),
    "producto_4": ("negativa", 0, random_dorder(4)),
    "producto_5": ("negativa", 0, random_dorder(5)),
    "kr": ("negativa", 0, lambda n, r: kr_order(n, r)),
    "percolacion_005": ("negativa", 0, lambda n, r: transitive_percolation(n, 0.05, r)),
    "capas_gruesas": ("negativa", 0, lambda n, r: thick_layer_order(n, 8, 0.5, r)),
}


# --- Partes ----------------------------------------------------------------------------------------


def abund(c: Bool) -> list[int]:
    return [int(v) for v in interval_abundances(c, K_MAX)]


def ref_rows(d: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for n in NS:
        ab = [interval_abundances(causal_matrix(sprinkle_diamond(n, d, rng_for(ENT_REF, d, n, s))), K_MAX) for s in range(N_REF)]
        rows.append({"d": d, "N": n, "ref": np.mean(np.array(ab, dtype=float), axis=0).tolist()})
    return rows


def cal_rows(d: int, n: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for s in range(N_CAL):
        c = causal_matrix(sprinkle_diamond(n, d, rng_for(ENT_CAL, FAM_PLANA, d, n, s)))
        m = measure(c) if s < N_CAL_SPEC else {"N": n, "abund": abund(c)}
        rows.append({"familia": "plana", "d": d, "seed": s, **m})
    if d in FRW_DIMS:
        for s in range(N_FRW):
            c = causal_matrix(frw_power_points(n, d, 2.0 * d, rng_for(ENT_CAL, FAM_FRW, d, n, s)))
            rows.append({"familia": "frw", "d": d, "seed": s, "N": n, "abund": abund(c)})
    return rows


def test_rows(name: str) -> list[dict[str, Any]]:
    role, dtrue, gen = TEST[name]
    k = list(TEST).index(name)
    rows: list[dict[str, Any]] = []
    for n in NS:
        for s in range(N_TEST):
            rows.append({"estructura": name, "papel": role, "d_verdadera": dtrue, "seed": s, **measure(gen(n, rng_for(ENT_TEST, k, n, s)))})
    return rows


def load(part: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = json.loads((OUT / f"{part}.json").read_text())["rows"]
    return rows


def load_cal(d: int) -> list[dict[str, Any]]:
    if (OUT / f"cal_{d}.json").exists():
        return load(f"cal_{d}")
    return [r for n in NS for r in load(f"cal_{d}_{n}")]


def final(head: str) -> dict[str, Any]:
    refs = {(r["d"], r["N"]): np.array(r["ref"], dtype=float) for d in DIMS for r in load(f"ref_{d}")}
    cal = {d: load_cal(d) for d in DIMS}
    flat = {d: [r for r in cal[d] if r["familia"] == "plana"] for d in DIMS}
    spec = {d: [r for r in flat[d] if "beta_L" in r] for d in DIMS}
    curve = {n: {d: float(np.nanmedian([r["beta_L"] for r in spec[d] if r["N"] == n])) for d in DIMS} for n in NS}

    def dist(r: dict[str, Any], d: int) -> float:
        return profile_distance(np.array(r["abund"], dtype=float), refs[(d, r["N"])])

    def dform(r: dict[str, Any], d: int) -> float:
        return shape_distance(np.array(r["abund"], dtype=float), refs[(d, r["N"])])

    def ds(r: dict[str, Any]) -> float:
        return d_spec(r["beta_L"], curve[r["N"]])

    delta = {(d, n): 1.5 * float(np.percentile([dist(r, d) for r in cal[d] if r["N"] == n], 95)) for d in DIMS for n in NS}
    delta_f = {(d, n): 1.5 * float(np.percentile([dform(r, d) for r in cal[d] if r["N"] == n], 95)) for d in DIMS for n in NS}
    win = {d: (F3_WIN if d == 4 else {k: tuple(v) for k, v in windows(spec[d]).items()}) for d in DIMS}
    res: dict[str, Any] = {}
    for name, (role, dtrue, _) in TEST.items():
        rr = load(f"test_{name}")
        dsm = float(np.nanmedian([ds(r) for r in rr]))
        dstar = int(min(5, max(2, round(dsm)))) if math.isfinite(dsm) else None
        b = float(np.nanmedian([r["beta_L"] for r in rr]))
        ks = [r["rodilla_L"] for r in rr if r["rodilla_L"] is not None]
        km = float(np.median(ks)) if ks else None
        if dstar is None:
            firma, ratio_f, ratio, w = False, float("nan"), float("nan"), None
        else:
            w = win[dstar]
            firma = bool(w["beta"][0] <= b <= w["beta"][1] and km is not None and w["rodilla"][0] <= km <= w["rodilla"][1])
            ratio_f = float(np.median([dform(r, dstar) / delta_f[(dstar, r["N"])] for r in rr]))
            ratio = float(np.median([dist(r, dstar) / delta[(dstar, r["N"])] for r in rr]))
        cls = bool(firma and ratio_f <= 1.0)
        cls_d = bool(firma and ratio <= 1.0)
        ok: bool | None
        if role == "positiva":
            ok = cls and dstar == dtrue
        elif role == "negativa":
            ok = not cls
        else:
            ok = None
        res[name] = {"papel": role, "d_verdadera": dtrue, "d_spec_mediana": dsm, "d_estrella": dstar,
                     "beta_L": b, "rodilla_L": km, "ventana_d_estrella": w, "firma": firma,
                     "razon_Dforma_deltaforma": ratio_f, "razon_D_delta": ratio,
                     "clasificado": cls, "clasif_D_info": cls_d, "correcto": ok}
    verd = "PASA" if all(v["correcto"] for v in res.values() if v["correcto"] is not None) else "FALLA"
    return {"commit": head, "curva_beta_L": curve, "ventanas": win,
            "delta_forma": {f"{d},{n}": v for (d, n), v in delta_f.items()},
            "delta": {f"{d},{n}": v for (d, n), v in delta.items()},
            "veredicto": verd, "detalle": res}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", required=True)
    args = ap.parse_args()
    part: str = args.part
    if "R10-1 | **Firmada" not in DOC.read_text():
        sys.exit("F7 requiere la firma R10-1")
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip():
        sys.exit("exige arbol git limpio")
    t0 = time.time()
    rows: list[dict[str, Any]]
    tail = part.split("_")
    if part.startswith("ref_") and len(tail) == 2 and tail[1].isdigit() and int(tail[1]) in DIMS:
        OUT.mkdir(parents=True, exist_ok=True)
        rows = ref_rows(int(tail[1]))
    elif part.startswith("cal_") and len(tail) == 2 and tail[1].isdigit() and int(tail[1]) in DIMS:
        OUT.mkdir(parents=True, exist_ok=True)
        rows = [r for n in NS for r in cal_rows(int(tail[1]), n)]
    elif (part.startswith("cal_") and len(tail) == 3 and tail[1].isdigit() and tail[2].isdigit()
          and int(tail[1]) in DIMS and int(tail[2]) in NS):
        OUT.mkdir(parents=True, exist_ok=True)
        rows = cal_rows(int(tail[1]), int(tail[2]))
    elif part.startswith("test_") and part[5:] in TEST:
        OUT.mkdir(parents=True, exist_ok=True)
        rows = test_rows(part[5:])
    elif part == "final":
        out = final(head)
        (OUT / "summary.json").write_text(json.dumps(out, indent=2, default=str))
        print(json.dumps({"veredicto": out["veredicto"]}))
        return
    else:
        sys.exit("parte desconocida")
    (OUT / f"{part}.json").write_text(json.dumps({"commit": head, "rows": rows}, default=str))
    print(part, "ok", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
