"""R22, fase 3, puntos 1, 2, 3 y 5: auditoría de R20 sobre los datos ya guardados (no genera datos nuevos).

Pedido: dictamen del Consejo, sesión 11 (docs/OMEGA_CONSEJO_REV22_R20.md, §3). Lee results/omega_r20 (referencias y
estructuras de la corrida confirmatoria, commit c2e227d) y escribe results/omega_r22/auditoria.json:
  1. máscaras: componentes finitas por q y por conjunto de referencias, y componentes no medibles por brazo;
  2. sensibilidad a la máscara: (i) máscara fijada solo con las 100 de ajuste; (ii) vector común a los 4 q;
     cobertura de validación y L_q de todas las estructuras con cada máscara;
  5. borde cercano a 0.8: distribución de D/radio de la plana (validación) contra D_obs de las copias en cada q,
     con la referencia teórica de una distancia de Mahalanobis nula (raíz de chi² con p grados de libertad).
Ningún número de esta auditoría reemplaza el veredicto de R20; es descriptiva (nivel 2).
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy import stats

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "tools"))

import r14b  # noqa: E402
import r20  # noqa: E402

SRC = _ROOT / "results" / "omega_r20"
OUT = _ROOT / "results" / "omega_r22" / "auditoria.json"
COMP = r14b.COMPONENTS
QS = r20.QS


def load_refs() -> dict[int, list[np.ndarray]]:
    out = {}
    for qi in range(len(QS)):
        out[qi] = [np.array([r14b.obs_vector(json.loads((SRC / "refs" / c / f"q{qi}_s{s:03d}.json").read_text())["row"])
                             for s in range(100)]) for c in r20.CONJUNTOS]
    return out


def load_structs() -> dict[str, dict[int, dict[str, np.ndarray]]]:
    summ = json.loads((SRC / "summary.json").read_text())
    out: dict[str, Any] = {}
    for name in summ["estructuras"]:
        reals = sorted((SRC / "structs" / name).glob("r*.json"))
        per = {qi: {"obs": [], "ctrl": []} for qi in range(len(QS))}
        for p in reals:
            for row in json.loads(p.read_text())["rows"]:
                per[int(row["qi"])]["obs"].append(r14b.obs_vector(row["obs"]))
                per[int(row["qi"])]["ctrl"].append(r14b.obs_vector(row["ctrl"]))
        out[name] = {qi: {k: np.array(v) for k, v in d.items()} for qi, d in per.items()}
    return out


class MaskedBand(r20.ConformalBand):
    """ConformalBand con una máscara dada (en vez de «finitas en las 300»)."""

    def __init__(self, fit: np.ndarray, cal: np.ndarray, val: np.ndarray, comps: tuple[int, ...]) -> None:
        self.comps = comps
        self.band = r14b.Band(fit, fit[:0], comps, r14b.comp_floors())
        self.d_cal = self.band.dist(cal[:, list(comps)], sub=True)
        self.radio, self.k, self.j_corte = r20.conformal_radius(self.d_cal)
        self.D_val = self.band.dist(val[:, list(comps)], sub=True) / self.radio
        self.n_val = int(self.D_val.size)
        self.cubiertas = int(np.count_nonzero(self.D_val <= 1.0))
        self.minimo = r20.min_validas(self.n_val)
        self.ok = self.cubiertas >= self.minimo


def Dvec(cb: r20.ConformalBand, x: np.ndarray) -> np.ndarray:
    return np.array([cb.D(v) for v in x])


def finitos(x: np.ndarray) -> dict[str, int]:
    return {COMP[j]: int(np.sum(~np.isfinite(x[:, j]))) for j in range(x.shape[1]) if np.any(~np.isfinite(x[:, j]))}


def main() -> None:
    refs = load_refs()
    sts = load_structs()
    summ = json.loads((SRC / "summary.json").read_text())
    res: dict[str, Any] = {"fuente": "results/omega_r20 (commit c2e227d)", "q": list(QS)}

    # 1. Máscaras y faltantes
    m1: dict[str, Any] = {}
    for qi in range(len(QS)):
        f, c, v = refs[qi]
        m1[str(qi)] = {"q": QS[qi],
                       "no_finitas_por_conjunto": {n: finitos(x) for n, x in zip(r20.CONJUNTOS, (f, c, v))},
                       "mascara_300": [COMP[i] for i in np.flatnonzero(np.all(np.isfinite(np.vstack([f, c, v])), axis=0))],
                       "mascara_ajuste": [COMP[i] for i in np.flatnonzero(np.all(np.isfinite(f), axis=0))]}
        brazos: dict[str, Any] = {}
        for name, d in sts.items():
            for arm in ("obs", "ctrl"):
                nf = finitos(d[qi][arm])
                if nf:
                    brazos[f"{name}:{arm}"] = nf
        m1[str(qi)]["no_finitas_en_estructuras"] = brazos
    res["mascaras"] = m1

    # 2. Sensibilidad a la máscara
    masks: dict[str, dict[int, tuple[int, ...]]] = {"oficial": {}, "solo_ajuste": {}, "comun_4q": {}}
    for qi in range(len(QS)):
        f, c, v = refs[qi]
        masks["oficial"][qi] = tuple(int(i) for i in np.flatnonzero(np.all(np.isfinite(np.vstack([f, c, v])), axis=0)))
        masks["solo_ajuste"][qi] = tuple(int(i) for i in np.flatnonzero(np.all(np.isfinite(f), axis=0)))
    comun = tuple(sorted(set.intersection(*(set(m) for m in masks["oficial"].values()))))
    masks["comun_4q"] = {qi: comun for qi in range(len(QS))}
    sens: dict[str, Any] = {}
    for mname, mm in masks.items():
        bands = {}
        for qi in range(len(QS)):
            f, c, v = refs[qi]
            ok_all = all(np.all(np.isfinite(x[:, list(mm[qi])])) for x in (f, c, v))
            bands[qi] = MaskedBand(f, c, v, mm[qi]) if ok_all else None
        cell: dict[str, Any] = {"mascara": {str(qi): [COMP[i] for i in mm[qi]] for qi in mm},
                                "cobertura": {str(qi): (bands[qi].cubiertas if bands[qi] else None) for qi in bands},
                                "radio": {str(qi): (bands[qi].radio if bands[qi] else None) for qi in bands},
                                "L_mediana": {}, "D_obs_mediana": {}, "D_ctrl_mediana": {}}
        for name, d in sts.items():
            Lm, Do, Dc = [], [], []
            for qi in range(len(QS)):
                if bands[qi] is None:
                    Lm.append(None); Do.append(None); Dc.append(None)
                    continue
                do, dc = Dvec(bands[qi], d[qi]["obs"]), Dvec(bands[qi], d[qi]["ctrl"])
                Lm.append(r20.ext_median(r20.log_ratio(do, dc)))
                Do.append(r20.ext_median(do)); Dc.append(r20.ext_median(dc))
            cell["L_mediana"][name] = Lm
            cell["D_obs_mediana"][name] = Do
            cell["D_ctrl_mediana"][name] = Dc
        sens[mname] = cell
    # control: la máscara oficial debe reproducir summary.json
    rep = max(abs(sens["oficial"]["L_mediana"][n][qi] - summ["estructuras"][n]["por_q"][qi]["L_mediana"])
              for n in sts for qi in range(len(QS))
              if math.isfinite(summ["estructuras"][n]["por_q"][qi]["L_mediana"] or math.nan)
              and sens["oficial"]["L_mediana"][n][qi] is not None and math.isfinite(sens["oficial"]["L_mediana"][n][qi]))
    sens["reproduce_summary_max_dif_L"] = rep
    res["sensibilidad_mascara"] = sens

    # 5. Borde cercano a 0.8
    borde: dict[str, Any] = {}
    for qi in range(len(QS)):
        b = summ["bandas"][str(qi)]
        dval = np.array(b["validacion"]["D_sobre_radio"])
        p = int(b["dim"])
        radio = float(b["radio"])
        teor = math.sqrt(stats.chi2.ppf(0.5, p)) / math.sqrt(stats.chi2.ppf(96 / 101, p))
        cell = {"q": QS[qi], "dim": p,
                "plana_validacion": {"mediana": float(np.median(dval)), "p10": float(np.percentile(dval, 10)),
                                     "p90": float(np.percentile(dval, 90)), "min": float(dval.min()), "max": float(dval.max())},
                "referencia_chi_mediana_sobre_cuantil_96_101": teor,
                "estructuras": {}}
        for name, e in summ["estructuras"].items():
            if e["tipo"] != "copias" and name not in ("plana", "de_sitter"):
                continue
            do = np.array([row[qi] for row in e["D_obs"]], dtype=float)
            dc = np.array([row[qi] for row in e["D_ctrl"]], dtype=float)
            do_f = do[np.isfinite(do)]
            mw = stats.mannwhitneyu(do_f, dval, alternative="two-sided") if do_f.size else None
            ks = stats.ks_2samp(do_f, dval) if do_f.size else None
            Dc_med = float(np.nanmedian(dc))
            cell["estructuras"][name] = {
                "vota": e["vota"], "D_obs_mediana": float(np.median(do_f)) if do_f.size else None,
                "D_obs_p10_p90": [float(np.percentile(do_f, 10)), float(np.percentile(do_f, 90))] if do_f.size else None,
                "D_ctrl_mediana": Dc_med, "obs_en_banda": float(np.mean(do_f <= 1.0)) if do_f.size else None,
                "MW_p_vs_plana": float(mw.pvalue) if mw else None, "KS_p_vs_plana": float(ks.pvalue) if ks else None,
                "L_mediana": e["por_q"][qi]["L_mediana"],
                "L_si_obs_fuera_plano_tipico": math.log(float(np.median(dval)) / Dc_med) if Dc_med > 0 else None}
        borde[str(qi)] = cell
    res["borde"] = borde
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(f"escrito {OUT}")


if __name__ == "__main__":
    main()
