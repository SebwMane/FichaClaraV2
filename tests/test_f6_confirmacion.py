"""Pruebas rápidas de tools/f6_confirmacion.py (F6 confirmatoria). Semillas propias (850001), nunca las
entropías de la confirmación."""

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


def _f6c() -> ModuleType:
    tools = str(Path(__file__).resolve().parents[1] / "tools")
    if tools not in sys.path:
        sys.path.insert(0, tools)
    return importlib.import_module("f6_confirmacion")


def _rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence(850001, spawn_key=(seed,)))


Float = npt.NDArray[np.float64]


def _norm_c3(a: Float) -> Float:
    out: Float = np.maximum(np.hypot(a[:, 0], a[:, 1]), a[:, 2])
    return out


def _norm_c4(a: Float) -> Float:
    out: Float = np.maximum(np.linalg.norm(a, axis=1), np.max(a, axis=1) / 0.79)
    return out


def _norm_c7(a: Float) -> Float:
    q = a[:, 0] ** 2 + a[:, 1] ** 2
    out: Float = (q * q + a[:, 2] ** 4) ** 0.25
    return out


def test_sealed_cones_partial_order_and_region() -> None:
    f = _f6c()
    t = DIAMOND_HALF_HEIGHT
    norms = {3: _norm_c3, 4: _norm_c4, 6: f.c6_norm_bisect, 7: _norm_c7}
    for k in (3, 4, 6, 7):
        c, x, _ = f.SELLADAS[k][2](300, _rng(10 + k), _rng(20 + k))
        assert x is not None and x.shape == (300, 4)
        assert is_partial_order(c) and np.count_nonzero(c) > 0
        nx = norms[k](np.abs(x[:, 1:]))
        assert np.all(np.abs(x[:, 0]) <= t - nx + 1e-9)
        assert np.all(f.in_cone_region(x, f.CONE_WITHIN[k]))
        # la relación por bloques coincide con la norma directa
        dt = x[None, :, 0] - x[:, None, 0]
        d = np.abs(x[None, :, 1:] - x[:, None, 1:]).reshape(-1, 3)
        direct = (dt > 0) & (dt >= norms[k](d).reshape(300, 300))
        assert np.count_nonzero(direct != c) == 0
    for k in (0, 1, 2, 5):
        c, _, _ = f.SELLADAS[k][2](250, _rng(30 + k), _rng(40 + k))
        assert is_partial_order(c) and np.count_nonzero(c) > 0


def test_c6_power_sum_matches_bisection() -> None:
    f = _f6c()
    r = _rng(50)
    a = np.abs(r.uniform(-0.6, 0.6, size=(400, 3)))
    s = r.uniform(0.05, 0.8, size=400)
    lam = f.c6_norm_bisect(a)
    # se descartan los casos a menos de 1e-9 relativo de la frontera (ahí manda el redondeo)
    clear = np.abs(lam - s) > 1e-9 * s
    got = f._within_p_por_eje(a[:, 0], a[:, 1], a[:, 2], s)
    assert clear.sum() > 350
    assert np.array_equal(got[clear], (lam <= s)[clear])
    assert 0.2 < got.mean() < 0.8
    # potencias altas y 0/0: s = 0 solo contiene el origen; |x| ≫ s no desborda
    z = np.zeros(3)
    assert f._within_p_por_eje(z, z, z, z).all()
    assert not f._within_p_por_eje(np.array([1e-300]), z[:1], z[:1], z[:1]).any()
    assert not f._within_p_por_eje(np.array([0.0]), np.array([0.0]), np.array([1e10]), np.array([1e-10])).any()
    # norma: eje z con p = 32 da λ = |z|; valores de la tabla del archivo sellado (mín 0.6856, máx 1.0088)
    assert math.isclose(f.c6_norm_bisect(np.array([0.0, 0.0, 0.3]))[0], 0.3, rel_tol=1e-11)
    u = r.normal(size=(20000, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    nu = f.c6_norm_bisect(u)
    assert 0.68 < nu.min() < 0.70 and 1.0 < nu.max() < 1.0089


def test_c1_equals_minkowski_after_transform() -> None:
    f = _f6c()
    x = sprinkle_diamond(300, 4, _rng(60))
    xp = f.c1_transform(x)
    assert np.allclose(xp[:, 1:], (f.C1_A @ x[:, 1:].T).T)
    c = f.causal_blocks(xp, "mink", block=64)
    assert np.array_equal(c, causal_matrix(xp))
    dt = x[None, :, 0] - x[:, None, 0]
    dx = (x[None, :, 1:] - x[:, None, 1:]) @ f.C1_A.T
    direct = (dt > 0) & (dt >= np.linalg.norm(dx, axis=2))
    assert np.array_equal(direct, c)
    # A = diag(1/c) Rᵀ con R ortogonal de determinante 1
    r = f.C1_R
    assert np.allclose(r @ r.T, np.eye(3)) and math.isclose(float(np.linalg.det(r)), 1.0)
    assert np.allclose(np.linalg.svd(f.C1_A, compute_uv=False), sorted(1 / np.array(f.C1_C), reverse=True))
    gen_c, gen_x, _ = f.SELLADAS[1][2](200, _rng(61), _rng(62))
    assert gen_x is not None and np.all(np.linalg.norm(gen_x[:, 1:], axis=1) <= 0.5 - np.abs(gen_x[:, 0]) + 1e-12)
    assert np.array_equal(gen_c, causal_matrix(f.c1_transform(gen_x)))


def test_c5_partial_order_closure_and_count() -> None:
    f = _f6c()
    c = causal_matrix(sprinkle_diamond(400, 4, _rng(70)))
    cp, info = f.c5_delete_links(c, _rng(71))
    assert is_partial_order(cp) and info["orden_parcial"] and info["cierre_igual"]
    assert np.array_equal(f.transitive_step(cp, block=97), cp)
    assert not np.any(cp & ~c)
    assert info["borrados"] == round(0.13 * info["links_C"]) == np.count_nonzero(c) - np.count_nonzero(cp)
    # los borrados son links de C y siguen sin intermedios en C'
    lm = f.link_mask(c)
    removed = c & ~cp
    assert np.all(lm[removed])
    # la elección es rng.choice sobre L en orden lexicográfico
    ii, jj = np.nonzero(lm)
    sel = _rng(71).choice(ii.size, size=info["borrados"], replace=False)
    exp = c.copy()
    exp[ii[sel], jj[sel]] = False
    assert np.array_equal(exp, cp)
    assert info["links_Cp"] >= info["links_C"] - info["borrados"]


def test_c2_density_range_and_modes() -> None:
    f = _f6c()
    ks, th = f.c2_modes(_rng(80))
    # orden de consumo: kⱼ (4 componentes) y luego θⱼ, j = 1..5
    r = _rng(80)
    for j in range(5):
        assert np.array_equal(ks[j], r.uniform(-1.5, 1.5, size=4))
        assert th[j] == r.uniform(0, 2 * math.pi)
    assert np.all(np.abs(ks) <= 1.5) and np.all((th >= 0) & (th < 2 * math.pi))
    x = sprinkle_diamond(20000, 4, _rng(81))
    phi = f.c2_phi(x, ks, th)
    dens = np.exp(phi)
    assert np.all(dens >= math.exp(-0.35)) and np.all(dens <= math.exp(0.35))
    assert dens.max() / dens.min() <= math.exp(0.7)
    pts = f.c2_points(3000, _rng(82), ks, th)
    assert pts.shape == (3000, 4) and np.all(np.linalg.norm(pts[:, 1:], axis=1) <= 0.5 - np.abs(pts[:, 0]) + 1e-12)
    # la densidad aceptada sigue a exp(φ): la media de φ en los puntos supera a la de puntos uniformes
    assert f.c2_phi(pts, ks, th).mean() > phi.mean()


def test_parts_and_seeds() -> None:
    f = _f6c()
    assert f.parse_struct("cal_14") == ("cal_14", [4000, 8000], None)
    assert f.parse_struct("cal_15") is None
    assert f.parse_struct("vista_cono_p1_5_4d_8000") == ("vista_cono_p1_5_4d", [8000], None)
    assert f.parse_struct("sellada_0_8000_r9") == ("sellada_0", [8000], [9])
    assert f.parse_struct("sellada_0_8000_r10") is None
    assert f.parse_struct("sellada_8") is None
    assert f.parse_struct("vista_bloque") is None
    assert f.parse_nula("nula_8000") == (8000, None) and f.parse_nula("nula_16000") is None
    assert [f.cal_kind(g) for g in (0, 9, 10, 12, 13, 14)] == ["diamante_plano"] * 2 + ["de_sitter"] * 2 + ["frw_polvo"] * 2
    parts = f.recommended_parts()
    assert parts[0] == "ref" and parts[-1] == "final" and len(set(parts)) == len(parts)
    for p in parts[3:-1]:
        assert f.parse_struct(p) is not None
    assert [nm[:2] for nm, _, _ in f.SELLADAS] == [f"C{k}" for k in range(8)]
    roles = [role for _, role, _ in f.SELLADAS]
    assert roles.count("positiva") == 3 and roles.count("negativa") == 4 and roles.count("zona_gris") == 1
    assert f.sealed_hash_ok()


def test_small_pipeline_fit() -> None:
    """Instrumento + ajuste lineal en miniatura con semillas propias."""
    f = _f6c()
    ref = f.ref_rows(reps=200, ms=(10,), entropy=850001)[0]
    assert abs(ref["r_plano"] - 0.1) < 0.012
    rows = f.cal_rows(0, 1200, reals=[0, 1], entropy=850001)
    assert all(r["por_m"]["10"]["medible"] for r in rows)
    rp = {m: 0.1 for m in f.M_GRID}
    dr = {1200: f._dr_matrix(rows, rp), 2400: f._dr_matrix(rows, rp)}
    fit = f.fit_linear_boot(dr, _rng(90), n_boot=20)
    assert fit["medible"] and fit["x_distintos"] >= 2 and fit["replicas_validas"] == 20
    lo, hi = fit["b"]["ic95"]
    assert lo <= hi and abs(fit["b"]["valor"]) < 0.05
    # «no medible»: un solo valor distinto de x (solo m = 10 en un solo N)
    one = np.full((4, len(f.M_GRID)), np.nan)
    one[:, 0] = 0.01
    nm = f.fit_linear_boot({4000: one}, _rng(91), n_boot=5)
    assert nm["medible"] is False and nm["b"] is None and nm["x_distintos"] == 1
