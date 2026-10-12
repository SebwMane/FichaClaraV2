"""Piloto P4 (Rev. 3 §21). Partes: p40 | n64 | n128a | n128b | final. Constantes congeladas."""

from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import chisquare

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omega.causality.estimators import between_counts, interval_abundances, mm_dimension, ordering_fraction  # noqa: E402
from omega.causality.mcmc_general import init_state, sweep, height  # noqa: E402
from omega.causality.orders import is_partial_order, kr_order, random_bipartite  # noqa: E402
from omega.causality.profile import K_MAX, d_ref, profile_distance, reference_abundances  # noqa: E402
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond  # noqa: E402

EXPERIMENT_ID = 310
ENT_CHAIN, ENT_CAL, ENT_VAL = 20261044, 20261045, 20261046
BETAS = (-0.1, -0.03, -0.01, -0.003, 0.0, 0.003, 0.01, 0.03, 0.1)
THERM, MEAS, EVERY = 100, 100, 10
OUT = Path("results") / "p4"


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


def natural_sprinkle(n: int, rng: np.random.Generator) -> np.ndarray[Any, Any]:
    x = sprinkle_diamond(n, 4, rng)
    return causal_matrix(x[np.argsort(x[:, 0])])


def p40() -> dict[str, Any]:
    st = init_state(natural_sprinkle(32, rng_for(ENT_VAL, 1)))
    r = rng_for(ENT_VAL, 2)
    for _ in range(400):  # 400 barridos de 496 propuestas ~ 2e5 movimientos
        sweep(st, 0.01, r)
    full = between_counts(st.c)
    a_ok = is_partial_order(st.c) and bool(np.array_equal(st.b[st.c], full[st.c]))
    # (b) N = 4, beta = 0: uniforme sobre los ordenes etiquetados naturalmente
    pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    valid = []
    for bits in itertools.product((0, 1), repeat=6):
        c = np.zeros((4, 4), dtype=bool)
        for (i, j), on in zip(pairs, bits):
            c[i, j] = bool(on)
        if is_partial_order(c):
            valid.append(bits)
    index = {v: k for k, v in enumerate(valid)}
    st4 = init_state(np.zeros((4, 4), dtype=bool))
    r4 = rng_for(ENT_VAL, 3)
    counts = np.zeros(len(valid))
    for _ in range(200_000):
        sweep(st4, 0.0, r4)
        sweep(st4, 0.0, r4)
        counts[index[tuple(int(st4.c[i, j]) for i, j in pairs)]] += 1
    p = float(chisquare(counts).pvalue)
    # ENMIENDA P4-A1: con beta = -0.5 la cadena debe reproducir e^{-beta S4} exacto sobre los 40 ordenes
    from omega.causality.estimators import bd_action_4d
    beta_neg = -0.5
    w = np.array([np.exp(-beta_neg * bd_action_4d(_mat(v, pairs))) for v in valid])
    stn = init_state(np.zeros((4, 4), dtype=bool))
    rn = rng_for(ENT_VAL, 4)
    cn = np.zeros(len(valid))
    for _ in range(200_000):
        sweep(stn, beta_neg, rn)
        sweep(stn, beta_neg, rn)
        cn[index[tuple(int(stn.c[i, j]) for i, j in pairs)]] += 1
    pn = float(chisquare(cn, w / w.sum() * cn.sum()).pvalue)
    ok = a_ok and len(valid) == 40 and p > 0.001 and pn > 0.001
    return {"a_ok": a_ok, "n_posets_4": len(valid), "chi2_p": p, "chi2_p_beta_neg": pn, "pass": ok}


def _mat(bits: tuple[int, ...], pairs: list[tuple[int, int]]) -> np.ndarray[Any, Any]:
    c = np.zeros((4, 4), dtype=bool)
    for (i, j), on in zip(pairs, bits):
        c[i, j] = bool(on)
    return c


def classifier(n: int) -> tuple[float, np.ndarray[Any, Any], dict[str, Any]]:
    ref4 = reference_abundances(4, n, K_MAX, [rng_for(ENT_CAL, n, 1000 + s) for s in range(40)])
    dist = []
    for s in range(40):
        c = natural_sprinkle(n, rng_for(ENT_CAL, n, s))
        dist.append(profile_distance(interval_abundances(c, K_MAX), ref4))
    delta4 = 1.5 * float(np.percentile(dist, 95))

    def flag(c: np.ndarray[Any, Any]) -> bool:
        return d_ref(mm_dimension(c)) == 4 and profile_distance(interval_abundances(c, K_MAX), ref4) <= delta4

    pos = float(np.mean([flag(natural_sprinkle(n, rng_for(ENT_VAL, n, s))) for s in range(100)]))
    kr = float(np.mean([flag(kr_order(n, rng_for(ENT_VAL, n, 500 + s))) for s in range(20)]))
    bi = float(np.mean([flag(random_bipartite(n, 0.5, rng_for(ENT_VAL, n, 700 + s))) for s in range(20)]))
    val = {"pos": pos, "kr": kr, "bicapa": bi, "ok": pos >= 0.9 and kr <= 0.1 and bi <= 0.1, "delta4": delta4}
    return delta4, ref4, val


def run(n: int, betas: tuple[float, ...]) -> dict[str, Any]:
    delta4, ref4, val = classifier(n)

    def flag(c: np.ndarray[Any, Any]) -> bool:
        return d_ref(mm_dimension(c)) == 4 and profile_distance(interval_abundances(c, K_MAX), ref4) <= delta4

    rows = []
    for beta in betas:
        k = BETAS.index(beta)
        rng = rng_for(ENT_CHAIN, n, k)
        st = init_state(natural_sprinkle(n, rng))
        acc = val_n = 0
        samples = []
        for sw in range(THERM + MEAS):
            a, v = sweep(st, beta, rng)
            acc += a
            val_n += v
            if sw >= THERM and (sw - THERM + 1) % EVERY == 0:
                samples.append({"flag": flag(st.c), "d_mm": mm_dimension(st.c), "frac_rel": ordering_fraction(st.c),
                                "altura": height(st.c), "S4": st.s})
        rows.append({"beta": beta, "flag_fraction": float(np.mean([s["flag"] for s in samples])),
                     "d_mm": float(np.nanmean([s["d_mm"] for s in samples])),
                     "frac_rel": float(np.mean([s["frac_rel"] for s in samples])),
                     "altura": float(np.mean([s["altura"] for s in samples])),
                     "S4": float(np.mean([s["S4"] for s in samples])),
                     "aceptacion_sobre_validas": acc / max(val_n, 1)})
    return {"N": n, "validacion": val, "rows": rows}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["p40", "n64", "n128a", "n128b", "final"], required=True)
    args = ap.parse_args()
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip():
        sys.exit("exige arbol git limpio")
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if args.part == "final":
        p40d = json.loads((OUT / "p40.json").read_text())["data"]
        r64 = json.loads((OUT / "n64.json").read_text())["data"]
        ra = json.loads((OUT / "n128a.json").read_text())["data"]
        rb = json.loads((OUT / "n128b.json").read_text())["data"]
        r128 = {"N": 128, "validacion": ra["validacion"], "rows": ra["rows"] + rb["rows"]}
        out: dict[str, Any] = {"commit": head, "P4-0": p40d["pass"]}
        win = {}
        for r in (r64, r128):
            w = sorted(x["beta"] for x in r["rows"] if x["beta"] != 0 and x["flag_fraction"] >= 0.5)
            ok = r["validacion"]["ok"] and p40d["pass"]
            out[f"N{r['N']}"] = "INDETERMINADO" if not ok else ("VENTANA" if w else "SIN VENTANA")
            win[r["N"]] = w
        if out["N64"] == "VENTANA" and out["N128"] == "VENTANA":
            out["angosta"] = set(win[128]) < set(win[64])
        out["ventanas"] = {str(k): v for k, v in win.items()}
        (OUT / "summary.json").write_text(json.dumps(out, indent=2))
        print(json.dumps(out, indent=2))
        return
    data = (p40() if args.part == "p40" else run(64, BETAS) if args.part == "n64"
            else run(128, BETAS[:5]) if args.part == "n128a" else run(128, BETAS[5:]))
    (OUT / f"{args.part}.json").write_text(json.dumps({"commit": head, "seconds": round(time.time() - t0, 1), "data": data}))
    print(args.part, "ok", round(time.time() - t0, 1), "s", data.get("pass", ""))


if __name__ == "__main__":
    main()
