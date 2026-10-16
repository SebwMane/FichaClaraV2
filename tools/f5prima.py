"""F5′ (docs/OMEGA_CONSEJO_REV9_F5PRIMA.md §6). Partes: ref_<d> | cal_<d> | cal_<d>_<N> | test_<estructura> | final.

Orden firmado (R9-1): referencias -> calibración -> prueba -> final.

- ``ref_<d>``: 20 diamantes planos por N; referencia = media de N_k, k = 0..K_MAX.
- ``cal_<d>``: 40 diamantes planos por N (``"familia": "plana"``); semillas 0..9 con espectro de links,
  d_MM y d_mid (``measure``), semillas 10..39 solo abundancias. Para d ∈ {3, 4}, además 20 diamantes FRW
  por N (``"familia": "frw"``), solo abundancias.
- ``cal_<d>_<N>``: la misma parte restringida a un solo N, por si ``cal_<d>`` no cabe en un comando
  (≈ 9 min). Las filas y las semillas son idénticas a las de ``cal_<d>``; ``final`` lee ``cal_<d>.json``
  si existe y, si no, une ``cal_<d>_1600.json`` y ``cal_<d>_2400.json``.
- ``test_<estructura>``: 5 semillas por N de cada una de las 12 estructuras de §6.3.
- ``final``: umbrales, clasificación y veredicto único PASA/FALLA (§6.2).
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Callable, Literal

import numpy as np
import numpy.typing as npt

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from f1f2 import random_dorder  # noqa: E402
from f3f4 import windows  # noqa: E402
from f5 import d_spec, mad, measure  # noqa: E402

from omega.causality.estimators import interval_abundances  # noqa: E402
from omega.causality.orders import kr_order, thick_layer_order, transitive_percolation  # noqa: E402
from omega.causality.profile import K_MAX, profile_distance  # noqa: E402
from omega.causality.sprinkling import DIAMOND_HALF_HEIGHT, causal_matrix, sprinkle_diamond  # noqa: E402

EXPERIMENT_ID = 318
ENT_REF, ENT_CAL, ENT_TEST = 20261069, 20261070, 20261071
NS = (1600, 2400)
DIMS = (2, 3, 4, 5)
FRW_DIMS = (3, 4)
N_REF, N_CAL, N_CAL_SPEC, N_FRW, N_TEST = 20, 40, 10, 20, 5
FAM_PLANA, FAM_FRW = 0, 1
F3_WIN = {"beta": (-0.363, -0.263), "rodilla": (1.511, 6.046)}
DOC = Path("docs") / "OMEGA_CONSEJO_REV9_F5PRIMA.md"
OUT = Path("results") / "omega_f5prima"
Bool = npt.NDArray[np.bool_]
Float = npt.NDArray[np.float64]
Gen = Callable[[int, np.random.Generator], Bool]
Norm = Literal["max", "l1"]


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


# --- Generadores -----------------------------------------------------------------------------------


def de_sitter_points(n: int, d: int, r: np.random.Generator) -> Float:
    """Diamante de `sprinkle_diamond` con η = t − 1 ∈ [−1.5, −0.5], aceptación (0.5/|η|)^d.

    Para d = 4 consume el generador igual que `f3f4.de_sitter_4d` (misma distribución y mismos puntos).
    """
    pts: list[Float] = []
    total = 0
    while total < n:
        cand = sprinkle_diamond(4 * n, d, r)
        eta = cand[:, 0] - 1.0
        keep = cand[r.random(cand.shape[0]) < (0.5 / np.abs(eta)) ** d]
        pts.append(keep)
        total += keep.shape[0]
    out: Float = np.concatenate(pts)[:n]
    return out


def frw_points(n: int, d: int, r: np.random.Generator) -> Float:
    """FRW con a(η) ∝ η², coordenadas conformes: η = t + 1 ∈ [0.5, 1.5], aceptación (η/1.5)^{2d}."""
    pts: list[Float] = []
    total = 0
    while total < n:
        cand = sprinkle_diamond(4 * n, d, r)
        eta = cand[:, 0] + 1.0
        keep = cand[r.random(cand.shape[0]) < (eta / 1.5) ** (2 * d)]
        pts.append(keep)
        total += keep.shape[0]
    out: Float = np.concatenate(pts)[:n]
    return out


def de_sitter(d: int) -> Gen:
    """Generador (n, rng) -> matriz causal; la relación es la de Minkowski (invariancia conforme)."""

    def gen(n: int, r: np.random.Generator) -> Bool:
        return causal_matrix(de_sitter_points(n, d, r))

    return gen


def frw(d: int) -> Gen:
    def gen(n: int, r: np.random.Generator) -> Bool:
        return causal_matrix(frw_points(n, d, r))

    return gen


def spatial_norm(x: Float, norm: Norm) -> Float:
    """Norma por filas de las coordenadas espaciales: "max" (cubo) o "l1" (octaedro)."""
    a = np.abs(x)
    out: Float = np.max(a, axis=1) if norm == "max" else np.sum(a, axis=1)
    return out


def sprinkle_polyhedral(n: int, norm: Norm, rng: np.random.Generator) -> Float:
    """n puntos uniformes en {|t| <= T, norm(x) <= T − |t|} (1 + 3 dimensiones); rechazo en [−T, T]^4."""
    if n < 0:
        raise ValueError("n debe ser >= 0")
    t_half = DIAMOND_HALF_HEIGHT
    out = np.empty((n, 4), dtype=np.float64)
    filled = 0
    while filled < n:
        batch = max(64, 32 * (n - filled))
        cand = rng.uniform(-t_half, t_half, size=(batch, 4))
        ok = cand[spatial_norm(cand[:, 1:], norm) <= t_half - np.abs(cand[:, 0])]
        take = min(len(ok), n - filled)
        out[filled : filled + take] = ok[:take]
        filled += take
    return out


def polyhedral_matrix(coords: Float, norm: Norm) -> Bool:
    """C[i, j] = True si Δt > 0 y Δt >= norm(Δx)."""
    t = coords[:, 0]
    dt = t[None, :] - t[:, None]
    dn = np.zeros_like(dt)
    for k in range(1, coords.shape[1]):
        diff = np.abs(coords[None, :, k] - coords[:, None, k])
        if norm == "max":
            np.maximum(dn, diff, out=dn)
        else:
            dn += diff
    rel: Bool = (dt > 0) & (dt >= dn)
    return rel


def polyhedral_cone(norm: Norm) -> Gen:
    if norm not in ("max", "l1"):
        raise ValueError("norm debe ser 'max' o 'l1'")

    def gen(n: int, r: np.random.Generator) -> Bool:
        return polyhedral_matrix(sprinkle_polyhedral(n, norm, r), norm)

    return gen


# El orden del diccionario es el índice de estructura de §6.3 (spawn_key).
TEST: dict[str, tuple[str, int, Gen]] = {
    "diamante_4d": ("positiva", 4, lambda n, r: causal_matrix(sprinkle_diamond(n, 4, r))),
    "diamante_3d": ("positiva", 3, lambda n, r: causal_matrix(sprinkle_diamond(n, 3, r))),
    "de_sitter_4d": ("positiva", 4, de_sitter(4)),
    "de_sitter_3d": ("positiva", 3, de_sitter(3)),
    "cono_cubico_4d": ("negativa", 0, polyhedral_cone("max")),
    "cono_octaedrico_4d": ("negativa", 0, polyhedral_cone("l1")),
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
            c = frw(d)(n, rng_for(ENT_CAL, FAM_FRW, d, n, s))
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

    def ds(r: dict[str, Any]) -> float:
        return d_spec(r["beta_L"], curve[r["N"]])

    delta = {(d, n): 1.5 * float(np.percentile([dist(r, d) for r in cal[d] if r["N"] == n], 95)) for d in DIMS for n in NS}
    win = {d: (F3_WIN if d == 4 else {k: tuple(v) for k, v in windows(spec[d]).items()}) for d in DIMS}
    c4 = spec[4]
    tau_a = max(0.15, 3 * 1.4826 * mad([r["d_mm"] - ds(r) for r in c4]))
    tau_b = max(0.15, 3 * 1.4826 * mad([r["d_mm"] - r["d_mid"] for r in c4 if math.isfinite(r["d_mid"])]))
    res: dict[str, Any] = {}
    for name, (role, dtrue, _) in TEST.items():
        rr = load(f"test_{name}")
        dsm = float(np.nanmedian([ds(r) for r in rr]))
        dstar = int(min(5, max(2, round(dsm)))) if math.isfinite(dsm) else None
        b = float(np.nanmedian([r["beta_L"] for r in rr]))
        ks = [r["rodilla_L"] for r in rr if r["rodilla_L"] is not None]
        km = float(np.median(ks)) if ks else None
        gap_a = float(np.nanmedian([abs(r["d_mm"] - ds(r)) for r in rr]))
        gap_b = float(np.nanmedian([abs(r["d_mm"] - r["d_mid"]) for r in rr]))
        if dstar is None:
            firma, ratio, w = False, float("nan"), None
        else:
            w = win[dstar]
            firma = bool(w["beta"][0] <= b <= w["beta"][1] and km is not None and w["rodilla"][0] <= km <= w["rodilla"][1])
            ratio = float(np.median([dist(r, dstar) / delta[(dstar, r["N"])] for r in rr]))
        cls = bool(firma and ratio <= 1.0)
        ok = (cls and dstar == dtrue) if role == "positiva" else not cls
        res[name] = {"papel": role, "d_verdadera": dtrue, "d_spec_mediana": dsm, "d_estrella": dstar,
                     "beta_L": b, "rodilla_L": km, "ventana_d_estrella": w, "firma": firma, "razon_D_delta": ratio,
                     "clasificado": cls, "correcto": ok,
                     "brecha_A": gap_a, "clasif_A_info": bool(firma and gap_a <= tau_a),
                     "brecha_B": gap_b, "clasif_B_info": bool(firma and gap_b <= tau_b)}
    verd = "PASA" if all(v["correcto"] for v in res.values()) else "FALLA"
    return {"commit": head, "curva_beta_L": curve, "ventanas": win,
            "delta": {f"{d},{n}": v for (d, n), v in delta.items()}, "tau_A": tau_a, "tau_B": tau_b,
            "veredicto": verd, "detalle": res}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", required=True)
    args = ap.parse_args()
    part: str = args.part
    if "R9-1 | **Firmada" not in DOC.read_text():
        sys.exit("F5′ requiere la firma R9-1")
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip():
        sys.exit("exige arbol git limpio")
    t0 = time.time()
    rows: list[dict[str, Any]]
    tail = part.split("_")
    if part.startswith("ref_") and len(tail) == 2 and int(tail[1]) in DIMS:
        OUT.mkdir(parents=True, exist_ok=True)
        rows = ref_rows(int(tail[1]))
    elif part.startswith("cal_") and len(tail) == 2 and int(tail[1]) in DIMS:
        OUT.mkdir(parents=True, exist_ok=True)
        rows = [r for n in NS for r in cal_rows(int(tail[1]), n)]
    elif part.startswith("cal_") and len(tail) == 3 and int(tail[1]) in DIMS and int(tail[2]) in NS:
        OUT.mkdir(parents=True, exist_ok=True)
        rows = cal_rows(int(tail[1]), int(tail[2]))
    elif part.startswith("test_") and part[5:] in TEST:
        OUT.mkdir(parents=True, exist_ok=True)
        rows = test_rows(part[5:])
    elif part == "final":
        out = final(head)
        (OUT / "summary.json").write_text(json.dumps(out, indent=2, default=str))
        print(json.dumps({"veredicto": out["veredicto"], "tau_A": round(out["tau_A"], 3), "tau_B": round(out["tau_B"], 3)}))
        return
    else:
        sys.exit("parte desconocida")
    (OUT / f"{part}.json").write_text(json.dumps({"commit": head, "rows": rows}, default=str))
    print(part, "ok", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
