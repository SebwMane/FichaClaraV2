"""Pruebas rápidas de tools/f6_intento2.py (F6 intento 2). Semillas propias (860001), nunca las entropías del
intento."""

from __future__ import annotations

import importlib
import math
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import numpy.typing as npt

from omega.causality.orders import is_partial_order
from omega.causality.sprinkling import DIAMOND_HALF_HEIGHT, causal_matrix, sprinkle_diamond

Float = npt.NDArray[np.float64]


def _f6() -> ModuleType:
    tools = str(Path(__file__).resolve().parents[1] / "tools")
    if tools not in sys.path:
        sys.path.insert(0, tools)
    return importlib.import_module("f6_intento2")


def _rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence(860001, spawn_key=(seed,)))


def test_d4_partial_order_and_counts() -> None:
    f = _f6()
    n = 500
    x = sprinkle_diamond(n, 4, _rng(1))
    c = causal_matrix(x)
    cpp, info = f.d4_build(x, c, _rng(2))
    assert is_partial_order(cpp) and info["orden_parcial"] and info["cierre_igual"]
    assert np.array_equal(f.transitive_step(cpp, block=97), cpp)
    assert not np.any(c & ~cpp)  # C ⊆ C''
    # candidatos: Δt > 0, no relacionados en C y a menos de δ² del cono
    cand = f.d4_candidates(x, c, f.d4_delta2(n), block=64)
    dt = x[None, :, 0] - x[:, None, 0]
    dx2 = np.sum((x[None, :, 1:] - x[:, None, 1:]) ** 2, axis=2)
    assert np.array_equal(cand, (dt > 0) & ~c & (dx2 - dt * dt <= (math.pi / 24 / n) ** 0.5))
    ii, jj = np.nonzero(cand)
    add = _rng(2).random(ii.size) < 0.10
    s = c.copy()
    s[ii[add], jj[add]] = True
    assert info["candidatos"] == ii.size and info["agregados"] == int(add.sum()) > 0
    assert info["relaciones_C"] == np.count_nonzero(c) and info["relaciones_S"] == np.count_nonzero(c) + int(add.sum())
    assert info["relaciones_Cpp"] == np.count_nonzero(cpp) >= info["relaciones_S"]
    assert not np.any(s & ~cpp)  # S ⊆ C''
    assert math.isclose(info["frac_nueva"], (info["relaciones_Cpp"] - info["relaciones_C"]) / info["relaciones_C"])
    assert np.all(dt[cpp] > 0)
    assert info["links_C"] == np.count_nonzero(f.link_mask(c)) and info["pasos_cierre"] >= 1


def test_d5_partial_order_and_twins() -> None:
    f = _f6()
    n = 500
    c_full, x, extra = f.d5_gemelos(n, _rng(10), _rng(11))
    info = extra["d5"]
    assert c_full.shape == (n, n) and is_partial_order(c_full)
    assert info["n_dup"] == 50 and info["n_base"] == 450 and info["gemelos_construidos"] == 50
    assert info["pares_gemelos"] == info["pares_esperados"] == 50 + info["pares_previos_en_C"] + (
        info["pares_esperados"] - 50 - info["pares_previos_en_C"])
    assert info["pares_esperados"] >= 50 + info["pares_previos_en_C"]
    # reconstrucción: C' = C[φ][:, φ] permutada
    xb = sprinkle_diamond(450, 4, _rng(10))
    c = causal_matrix(xb)
    aux = _rng(11)
    orig = aux.choice(450, size=50, replace=False)
    perm = aux.permutation(n)
    phi = np.concatenate([np.arange(450), orig])
    assert np.array_equal(c_full, c[np.ix_(phi, phi)][np.ix_(perm, perm)])
    assert len(set(orig.tolist())) == 50
    # sin pares previos (se quitan los elementos que comparten pasado y futuro, p. ej. los aislados):
    # exactamente n_dup pares con el mismo pasado y futuro
    cs = c
    while f.twin_pairs(cs):
        _, inv, cnt = np.unique(f._row_keys(cs), axis=0, return_inverse=True, return_counts=True)
        keep = np.flatnonzero(cnt[np.asarray(inv).reshape(-1)] == 1)
        cs = cs[np.ix_(keep, keep)]
    assert cs.shape[0] > 300
    cp, _, _, info2 = f.d5_build(cs, cs.shape[0] + 40, _rng(12))
    assert info2["pares_gemelos"] == 40 == info2["pares_esperados"] and info2["gemelos_construidos"] == 40
    assert is_partial_order(cp)


def test_twin_pairs_counts_isolated() -> None:
    f = _f6()
    c = np.zeros((5, 5), dtype=bool)
    c[0, 1] = True
    # 2, 3, 4 aislados: 3 pares con mismo pasado y futuro
    assert f.twin_pairs(c) == 3
    prev, exp = f.d5_expected_pairs(c, np.array([0, 2]))
    assert prev == 3 and exp == 1 + 6


def _region_ok(f: ModuleType, x: Float, gauge: object) -> bool:
    t = DIAMOND_HALF_HEIGHT
    fx = gauge(x[:, 1], x[:, 2], x[:, 3])  # type: ignore[operator]
    fm = gauge(-x[:, 1], -x[:, 2], -x[:, 3])  # type: ignore[operator]
    return bool(np.all(x[:, 0] + t >= fx) and np.all(t - x[:, 0] >= fm))


def test_d6_d7_d8_partial_order_and_region() -> None:
    f = _f6()
    n = 300
    for k in (6, 7, 8):
        c, x, _ = f.SELLADAS[k][2](n, _rng(20 + k), _rng(30 + k))
        assert x is not None and x.shape == (n, 4)
        assert is_partial_order(c) and np.count_nonzero(c) > 0
        g = f.GAUGES[k]
        assert _region_ok(f, x, g) and np.all(f.in_gauge_region(x, g))
        dt = x[None, :, 0] - x[:, None, 0]
        d = (x[None, :, 1:] - x[:, None, 1:]).reshape(-1, 3)
        direct = (dt > 0) & (dt >= g(d[:, 0], d[:, 1], d[:, 2]).reshape(n, n))
        assert np.array_equal(direct, c)
        assert np.array_equal(f.gauge_blocks(x, g, block=37), c)
    # D7: puntos fuera de [−T, T] en z aparecen (la caja z llega a ±L·T)
    _, x7, _ = f.SELLADAS[7][2](3000, _rng(40), _rng(41))
    assert x7 is not None and np.abs(x7[:, 3]).max() > DIAMOND_HALF_HEIGHT


def test_gauges_values() -> None:
    f = _f6()
    r = _rng(50)
    u = r.normal(size=(200000, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    # D6: ε = 0.25, no par; D8: ε = 0.32(√3 − 1); D7: ε = 1 − 1/L
    f6 = f.d6_gauge(u[:, 0], u[:, 1], u[:, 2])
    assert 0.9999999 < f6.min() and 1.2499 < f6.max() <= 1.25 + 1e-12
    assert math.isclose(float(f.d6_gauge(*(f.D6_N[:, None]))[0]), 1.25)
    assert math.isclose(float(f.d6_gauge(*(-f.D6_N[:, None]))[0]), 1.0)
    f8 = f.d8_gauge(u[:, 0], u[:, 1], u[:, 2])
    assert 0.9999 < f8.min() and f8.max() <= 0.68 + 0.32 * math.sqrt(3) + 1e-12
    f7 = f.d7_gauge(u[:, 0], u[:, 1], u[:, 2])
    assert 1 / 1.45 <= float(f7.min()) < 1 / 1.45 * 1.01 and float(f7.max()) <= 1 + 1e-12
    assert math.isclose(float(f.d7_gauge(*np.array([[0.0], [0.0], [-1.0]]))[0]), 1 / 1.45)
    # D7: forma algebraica = fórmula cerrada literal del sello
    assert np.allclose(f7, f.d7_gauge_trig(u), rtol=1e-12, atol=1e-14)
    # D7: gauge = sup_ν x·ν / h(ν), h(ν) = max(1, L|ν₃|) (definición por función soporte)
    nu = r.normal(size=(400000, 3))
    nu /= np.linalg.norm(nu, axis=1)[:, None]
    h = np.maximum(1.0, 1.45 * np.abs(nu[:, 2]))
    dirs = np.array([[0, 0, 1.0], [1.0, 0, 0], [0.3, -0.2, 0.9], [0.5, 0.1, 0.4], [0.6, 0.0, 0.8], [-0.2, 0.7, -0.6]])
    for v in dirs:
        sup = float(np.max(nu @ v / h))
        val = float(f.d7_gauge(*(v[:, None]))[0])
        assert val >= sup - 1e-12 and val - sup < 2e-3 * np.linalg.norm(v)


def test_d1_equals_minkowski_after_shear() -> None:
    f = _f6()
    c, x, _ = f.d1_cizalla(300, _rng(60), _rng(61))
    assert x is not None
    xp = f.d1_shear(x)
    assert np.allclose(xp[:, 1:], x[:, 1:] - np.outer(x[:, 0], [0.25, -0.15, 0.10]))
    assert np.array_equal(c, causal_matrix(xp))
    dt = x[None, :, 0] - x[:, None, 0]
    dx = x[None, :, 1:] - x[:, None, 1:] - dt[:, :, None] * f.D1_B[None, None, :]
    direct = (dt > 0) & (np.linalg.norm(dx, axis=2) <= dt)
    assert np.count_nonzero(direct != c) <= 2  # solo pares en el borde por redondeo
    assert is_partial_order(c)


def test_d0_d2_d3_small() -> None:
    f = _f6()
    c0, x0, _ = f.d0_kasner(120, _rng(70), _rng(71), k=40)
    assert x0 is not None and is_partial_order(c0) and np.all((x0[:, 0] >= 0.70) & (x0[:, 0] <= 1.35))
    for gen, seed in ((f.d2_frw_exp, 72), (f.d3_frw_rigida, 73)):
        c, x, _ = gen(800, _rng(seed), _rng(seed + 10))
        assert x is not None and is_partial_order(c)
        assert np.all(np.linalg.norm(x[:, 1:], axis=1) <= 0.5 - np.abs(x[:, 0]) + 1e-12)
        assert x[:, 0].mean() > 0.02  # densidad creciente en t
    # D2 sigue el patrón de frw_power_points con aceptación exp(4(t − 0.5))
    r = _rng(80)
    cand = sprinkle_diamond(4 * 50, 4, r)
    keep = cand[r.random(cand.shape[0]) < np.exp(4.0 * (cand[:, 0] - 0.5))]
    assert np.array_equal(f.d2_points(50, _rng(80))[: min(50, len(keep))], keep[:50])


def test_collapse_k_fabricated() -> None:
    f = _f6()
    # 4 realizaciones; columnas m = 10..160. N = 4000: Δr = m-dependiente; N = 8000 lo mismo en 2m salvo m = 20.
    base4 = np.array([[0.01, 0.02, 0.03, 0.04, np.nan]] * 4)
    base4 += np.array([[0.0], [0.002], [-0.002], [0.0]])
    base8 = np.full((4, 5), np.nan)
    base8[:, 1:] = base4[:, :4]
    base8[:, 2] += 0.01  # (2m = 40, 8000) frente a (20, 4000): d = 0.01
    res = f.collapse_k({4000: base4, 8000: base8})
    se = 0.002 * math.sqrt(2 / 3) / 2  # sd (ddof 1) de (0, .002, −.002, 0) / √4
    sej = math.sqrt(2) * se
    assert res["n_pares"] == 4 and res["medible"]
    assert math.isclose(res["K"], (0.01 / sej) ** 2, rel_tol=1e-9)
    assert [p["m"] for p in res["pares"]] == [10, 20, 40, 80]
    # celdas ausentes (menos de la mitad medibles) quitan pares; con uno solo, no medible
    b8 = base8.copy()
    b8[:3, 1:4] = np.nan
    r2 = f.collapse_k({4000: base4, 8000: b8})
    assert r2["n_pares"] == 1 and not r2["medible"] and r2["K"] is None
    # se = 0: d = 0 aporta 0
    flat = np.zeros((3, 5))
    assert f.collapse_k({4000: flat, 8000: flat})["K"] == 0.0


def test_parts_seeds_and_listing() -> None:
    f = _f6()
    assert f.parse_struct("cal_14") == ("cal_14", [4000, 8000], None)
    assert f.parse_struct("cal_15") is None
    assert f.parse_struct("sellada_8_8000_r3-7") == ("sellada_8", [8000], [3, 4, 5, 6, 7])
    assert f.parse_struct("sellada_0_8000_r29") == ("sellada_0", [8000], [29])
    assert f.parse_struct("sellada_0_8000_r25-30") is None and f.parse_struct("sellada_0_8000_r7-3") is None
    assert f.parse_struct("sellada_9") is None
    assert f.parse_struct("vista_C5_links_borrados_13_4000") == ("vista_C5_links_borrados_13", [4000], None)
    assert f.parse_struct("vista_de_sitter_4d_8000_r0-14") == ("vista_de_sitter_4d", [8000], list(range(15)))
    assert f.parse_nula("nula_8000_r0-9") == (8000, list(range(10))) and f.parse_nula("nula_16000") is None
    assert len(f.VISTAS) == 21 and list(f.VISTAS)[13:] == [nm for nm, _, _ in f.SELLADAS_I1]
    roles = [role for _, role, _ in f.SELLADAS]
    assert roles.count("positiva") == 4 and roles.count("negativa") == 5 and len(roles) == 9
    vr = list(f.VISTAS.values())
    assert vr.count("positiva") == 9 and vr.count("negativa") == 9 and vr.count("zona_gris") == 3
    assert f.sealed_hash_ok() and f.FIRMA in f.PREREG.read_text()
    # semillas: spawn_key (322, ...)
    a = f.rng_for(1, 2, 3).random()
    b = np.random.Generator(np.random.PCG64(np.random.SeedSequence(1, spawn_key=(322, 2, 3)))).random()
    assert a == b
    parts = f.recommended_parts()
    assert parts[0] == "refperfil" and parts[-1] == "final" and len(set(parts)) == len(parts)
    est = f.estimated_seconds()
    assert all(v * f.SAFETY <= f.LIMIT_S for v in est.values())
    # la partición cubre cada (estructura, N, realización) exactamente una vez
    cover: dict[tuple[str, int], list[int]] = {}
    for p in parts[1:-1]:
        st, nu = f.parse_struct(p), f.parse_nula(p)
        assert (st is None) != (nu is None)
        base, ns, reals = st if st is not None else ("nula", [nu[0]], nu[1])
        for n in ns:
            cover.setdefault((base, n), []).extend(reals if reals is not None else range(30))
    assert len(cover) == 2 * (1 + len(f.struct_bases()))
    assert all(sorted(v) == list(range(30)) for v in cover.values())


def test_small_pipeline_and_collect(tmp_path: Path) -> None:
    f = _f6()
    rows = f.cal_rows(0, 1200, reals=[0, 1], entropy=860001)
    assert all(r["por_m"]["10"]["medible"] for r in rows)
    files = [("cal_0_4000_r0-0", {"commit": "x", "rows": [dict(rows[0], N=4000)]}),
             ("cal_0_4000_r1", {"commit": "x", "rows": [dict(rows[1], N=4000)]})]
    got, used = f.collect_rows(files, "cal_0", 4000, 2)
    assert [r["realizacion"] for r in got] == [0, 1] and used == ["cal_0_4000_r0-0", "cal_0_4000_r1"]
    try:
        f.collect_rows(files + [("cal_0_4000", {"commit": "x", "rows": [dict(rows[1], N=4000)]})], "cal_0", 4000, 2)
        raise AssertionError("debía rechazar realizaciones repetidas")
    except RuntimeError:
        pass
