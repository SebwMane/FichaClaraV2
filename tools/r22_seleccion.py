"""R22, fase 3, punto 4: reproducción de la selección del guardián de R20 (auditoría, no confirmatoria).

Pedido: dictamen del Consejo, sesión 11 (docs/OMEGA_CONSEJO_REV22_R20.md, §3.4). Se regeneran las realizaciones de
selección del guardián con su esquema de semillas publicado (docs/sellado_r20/RECIBO.txt §2, entropía 20261114) y se
evalúan los requisitos de visibilidad de §3 del preregistro de R20 con dos métricas:
  guardian: banda propia del guardián (100 de ajuste + 100 de calibración, radio conforme 96.º) y D por par con las
            componentes de la banda finitas en los dos vectores (como r19.pair_D); debe reproducir el recibo;
  ejecutor: la banda conforme de la corrida confirmatoria (results/omega_r20/refs, componentes finitas en las 300
            referencias de cada q) y D nan si a un brazo le falta alguna (como tools/r20.py).
Universo: las 16 finales y las 4 S1 descartadas (docs/sellado_r20/universo_v1_descartado.py.txt; índices si = 12..15,
como en el recibo). No se usan las entropías 20261113, 20261115 ni 20261116 para generar nada.

  run  [--workers 4] [--out results/omega_r22/seleccion]
  summarize [--out ...]  -> <out>/summary.json, tabla.csv
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import multiprocessing as mp
import sys
import time
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
_TOOLS = _ROOT / "tools"
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_TOOLS))

import r14b  # noqa: E402
import r19  # noqa: E402
import r20  # noqa: E402

ENT_SEL = 20261114
EXP = 333
N0 = 8000
REALS = 10
QSEL = {0: 1.0, 3: 0.125}
FINAL = _ROOT / "tools" / "sellado_r20.py"
V1 = _ROOT / "docs" / "sellado_r20" / "universo_v1_descartado.py.txt"
OUT = _ROOT / "results" / "omega_r22" / "seleccion"
UMBRAL_Q1, UMBRAL_COPIAS, UMBRAL_REGION = 2.0, 1.5, 1.25


def rng(*key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(ENT_SEL, spawn_key=(EXP, *key))))


def _load(path: Path, name: str) -> ModuleType:
    spec = importlib.util.spec_from_loader(name, loader=None)
    assert spec is not None
    mod = importlib.util.module_from_spec(spec)
    exec(compile(path.read_text(), str(path), "exec"), mod.__dict__)
    return mod


_MODS: dict[str, ModuleType] = {}


def mods() -> dict[str, ModuleType]:
    if not _MODS:
        _MODS["final"] = _load(FINAL, "sellado_r20_final")
        _MODS["v1"] = _load(V1, "sellado_r20_v1")
    return _MODS


def universo() -> list[dict[str, Any]]:
    """(nombre, si, módulo, familia): índice si = posición en ESTRUCTURAS del universo de origen (recibo §2)."""
    out = []
    for key in ("final", "v1"):
        m = mods()[key]
        fam = {e: f["nombre"] for f in m.FAMILIAS for e in f["estructuras"]}
        tipo = {f["nombre"]: f["tipo"] for f in m.FAMILIAS}
        for si, (name, e) in enumerate(m.ESTRUCTURAS.items()):
            if key == "v1" and not name.startswith("S1"):
                continue
            out.append({"nombre": name, "si": si, "mod": key, "familia": fam[name], "rol": tipo[fam[name]],
                        "tipo": e["tipo"]})
    return out


def ref_task(t: dict[str, Any]) -> str:
    qi, s = t["qi"], t["s"]
    p = Path(t["out"]) / "refs_guardian" / f"q{qi}_s{s:03d}.json"
    if p.exists():
        return f"skip {p.name}"
    t0 = time.time()
    q = QSEL[qi]
    r = rng(8, qi, s, 0)
    n = N0 if q >= 1.0 else int(r.binomial(N0, q))
    c, _ = r14b.GEN["plana"](n, r, r)
    row = r14b.measure(c, rng(8, qi, s, 1))
    r19._write(p, {"qi": qi, "s": s, "row": row})
    return f"ref q{qi} s{s} N={n} {time.time() - t0:.1f} s"


def struct_task(t: dict[str, Any]) -> str:
    name, real = t["nombre"], t["real"]
    p = Path(t["out"]) / "estructuras" / name / f"r{real:02d}.json"
    if p.exists():
        return f"skip {name} r{real}"
    t0 = time.time()
    e = mods()[t["mod"]].ESTRUCTURAS[name]
    si = t["si"]
    base = None
    if e["tipo"] == "copias":
        c, base, _ = r19.copias_generar(N0, int(e["k"]), float(e["frac"]), rng(9, si, real, 0), rng(9, si, real, 1))
    else:
        c = r19._check_matrix(e["generar"](N0, rng(9, si, real, 0)), N0)
    rows = []
    for qi, q in QSEL.items():
        _, cd = r14b.dilute(c, q, rng(9, si, real, 2, qi) if q < 1.0 else None)
        n = int(cd.shape[0])
        obs = r14b.measure(cd, rng(9, si, real, 3, qi))
        del cd
        if e["tipo"] == "copias":
            assert base is not None
            cc, _ = r20.copias_control_n(base, n, int(e["k"]), float(e["frac"]), rng(9, si, real, 4, qi),
                                         rng(9, si, real, 6, qi))
        else:
            cc = r20.region_control_n(e["generar"], n, rng(9, si, real, 5, qi))
        ctrl = r14b.measure(cc, rng(9, si, real, 7, qi))
        del cc
        rows.append({"qi": qi, "q": q, "N": n, "obs": obs, "ctrl": ctrl})
    r19._write(p, {"nombre": name, "real": real, "rows": rows})
    return f"est {name} r{real} {time.time() - t0:.1f} s"


def _work(t: dict[str, Any]) -> str:
    return ref_task(t) if t["kind"] == "ref" else struct_task(t)


def cmd_run(a: argparse.Namespace) -> None:
    out = Path(a.out)
    tasks: list[dict[str, Any]] = [{"kind": "ref", "qi": qi, "s": s, "out": str(out)} for qi in QSEL for s in range(200)]
    tasks += [{"kind": "est", **{k: v for k, v in u.items() if k in ("nombre", "si", "mod")}, "real": r, "out": str(out)}
              for u in universo() for r in range(REALS)]
    if a.only:
        tasks = [t for t in tasks if t["kind"] == a.only]
    t0 = time.time()
    with mp.get_context("fork").Pool(a.workers, maxtasksperchild=4) as pool:
        for i, msg in enumerate(pool.imap_unordered(_work, tasks), 1):
            print(f"[{i}/{len(tasks)} {time.time() - t0:.0f} s] {msg}", flush=True)


class GuardianBand:
    """Banda del guardián: r14b.Band con las 100 de ajuste, radio conforme con las 100 de calibración; componentes
    finitas en las 200. D por par con las componentes finitas en los dos vectores (banda y radio recalculados)."""

    def __init__(self, fit: np.ndarray, cal: np.ndarray) -> None:
        self.fit, self.cal = fit, cal
        self.base = tuple(int(i) for i in np.flatnonzero(np.all(np.isfinite(np.vstack([fit, cal])), axis=0)))
        self._cache: dict[tuple[int, ...], tuple[r14b.Band, float]] = {}
        self.radio = self.band(self.base)[1]

    def band(self, comps: tuple[int, ...]) -> tuple[r14b.Band, float]:
        if comps not in self._cache:
            b = r14b.Band(self.fit, self.fit[:0], comps, r14b.comp_floors())
            rad, _, _ = r20.conformal_radius(b.dist(self.cal[:, list(comps)], sub=True))
            self._cache[comps] = (b, rad)
        return self._cache[comps]

    def pair(self, vo: np.ndarray, vc: np.ndarray) -> tuple[float, float, int]:
        comps = tuple(i for i in self.base if math.isfinite(vo[i]) and math.isfinite(vc[i]))
        b, rad = self.band(comps)
        return float(b.dist(vo)[0]) / rad, float(b.dist(vc)[0]) / rad, len(comps)


def _vecs(rows: list[dict[str, Any]]) -> np.ndarray:
    return np.array([r14b.obs_vector(r) for r in rows])


def executor_bands() -> dict[int, r20.ConformalBand]:
    src = _ROOT / "results" / "omega_r20" / "refs"
    out = {}
    for qi in QSEL:
        sets = [_vecs([json.loads((src / c / f"q{qi}_s{s:03d}.json").read_text())["row"] for s in range(100)])
                for c in r20.CONJUNTOS]
        out[qi] = r20.ConformalBand(*sets)
    return out


def _med(x: list[float]) -> float:
    a = np.asarray(x, dtype=np.float64)
    return float(np.median(a)) if a.size and not np.any(np.isnan(a)) else math.nan


def requisitos(tipo: str, m: dict[str, float]) -> bool:
    if any(math.isnan(v) for v in m.values()):
        return False
    if m["obs_q1"] < UMBRAL_Q1:
        return False
    if tipo == "copias":
        return m["ctrl_q8"] >= UMBRAL_COPIAS
    return m["obs_q8"] >= UMBRAL_REGION and m["ctrl_q8"] >= UMBRAL_REGION


def summarize(out: Path) -> dict[str, Any]:
    gb: dict[int, GuardianBand] = {}
    for qi in QSEL:
        rows = [json.loads((out / "refs_guardian" / f"q{qi}_s{s:03d}.json").read_text())["row"] for s in range(200)]
        v = _vecs(rows)
        gb[qi] = GuardianBand(v[:100], v[100:])
    eb = executor_bands()
    recibo = {**mods()["final"].IDENTIFICABLES}
    v1 = json.loads((_ROOT / "docs" / "sellado_r20" / "identificables_v1.json").read_text())["estructuras"]
    res: dict[str, Any] = {"bandas": {}, "estructuras": {}}
    for qi in QSEL:
        res["bandas"][str(qi)] = {
            "guardian": {"dim": len(gb[qi].base), "radio": gb[qi].radio,
                         "omitidas": [c for i, c in enumerate(r14b.COMPONENTS) if i not in gb[qi].base]},
            "ejecutor": {"dim": len(eb[qi].comps), "radio": eb[qi].radio,
                         "omitidas": [c for i, c in enumerate(r14b.COMPONENTS) if i not in eb[qi].comps]}}
    for u in universo():
        name = u["nombre"]
        d: dict[str, dict[str, list[float]]] = {m: {k: [] for k in ("obs_q1", "ctrl_q1", "obs_q8", "ctrl_q8")}
                                                for m in ("guardian", "ejecutor")}
        dims: list[int] = []
        nq8: list[int] = []
        for real in range(REALS):
            rec = json.loads((out / "estructuras" / name / f"r{real:02d}.json").read_text())
            for row in rec["rows"]:
                qi = int(row["qi"])
                tag = "q1" if qi == 0 else "q8"
                vo, vc = r14b.obs_vector(row["obs"]), r14b.obs_vector(row["ctrl"])
                go, gc, nd = gb[qi].pair(vo, vc)
                d["guardian"][f"obs_{tag}"].append(go)
                d["guardian"][f"ctrl_{tag}"].append(gc)
                d["ejecutor"][f"obs_{tag}"].append(eb[qi].D(vo))
                d["ejecutor"][f"ctrl_{tag}"].append(eb[qi].D(vc))
                if qi == 3:
                    dims.append(nd)
                    nq8.append(int(row["N"]))
        meds = {m: {k: _med(v) for k, v in dd.items()} for m, dd in d.items()}
        rec_g = recibo.get(name) or v1.get(name) or {}
        res["estructuras"][name] = {
            "familia": u["familia"], "universo": u["mod"], "tipo": u["tipo"],
            "medianas": meds,
            "identificable": {m: requisitos(u["tipo"], meds[m]) for m in meds},
            "recibo": {"identificable": rec_g.get("identificable"), "medidas": rec_g.get("medidas")},
            "dim_guardian_q8": sorted(set(dims)), "N_q8_mediana": float(np.median(nq8)),
            "n_nan_ejecutor": {k: int(np.sum(np.isnan(v))) for k, v in d["ejecutor"].items()}}
    (out / "summary.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float))
    with (out / "tabla.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["estructura", "familia", "metrica", "obs_q1", "ctrl_q1", "obs_q8", "ctrl_q8", "identificable",
                    "identificable_recibo"])
        for name, e in res["estructuras"].items():
            for m in ("guardian", "ejecutor"):
                mm = e["medianas"][m]
                w.writerow([name, e["familia"], m, *(f"{mm[k]:.3f}" for k in ("obs_q1", "ctrl_q1", "obs_q8", "ctrl_q8")),
                            e["identificable"][m], e["recibo"]["identificable"]])
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("run")
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--out", default=str(OUT))
    r.add_argument("--only", choices=["ref", "est"])
    s = sp.add_parser("summarize")
    s.add_argument("--out", default=str(OUT))
    a = ap.parse_args()
    if a.cmd == "run":
        cmd_run(a)
    else:
        res = summarize(Path(a.out))
        print(json.dumps(res["bandas"], indent=1))
        for name, e in res["estructuras"].items():
            g, x = e["medianas"]["guardian"], e["medianas"]["ejecutor"]
            print(f"{name:32s} G {g['obs_q1']:6.2f} {g['ctrl_q1']:6.2f} {g['obs_q8']:5.2f} {g['ctrl_q8']:5.2f} "
                  f"{e['identificable']['guardian']!s:5s} | E {x['obs_q1']:6.2f} {x['obs_q8']:5.2f} {x['ctrl_q8']:5.2f} "
                  f"{e['identificable']['ejecutor']!s:5s} | recibo {e['recibo']['identificable']}")


if __name__ == "__main__":
    main()
