"""Pruebas rápidas de tools/f6_exploracion.py (F6 exploratoria). Semillas propias, no la entropía de F6."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import ModuleType

import numpy as np

from omega.causality.estimators import between_counts
from omega.causality.orders import is_partial_order
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond
from omega.causality.torus import slab_causal_matrix, sprinkle_slab


def _f6() -> ModuleType:
    tools = str(Path(__file__).resolve().parents[1] / "tools")
    if tools not in sys.path:
        sys.path.insert(0, tools)
    return importlib.import_module("f6_exploracion")


def _rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence(880001, spawn_key=(seed,)))


def test_blocks_equal_original_relations() -> None:
    f = _f6()
    f5p = importlib.import_module("f5prima")
    f7 = importlib.import_module("f7")
    x = sprinkle_diamond(300, 4, _rng(1))
    assert np.array_equal(f.causal_blocks(x, "mink", block=64), causal_matrix(x))
    assert np.array_equal(f.causal_blocks(x, "max", block=64), f5p.polyhedral_matrix(x, "max"))
    assert np.array_equal(f.causal_blocks(x, "l1", block=64), f5p.polyhedral_matrix(x, "l1"))
    assert np.array_equal(f.causal_blocks(x, "lp", 3.0, block=64), f7.lp_matrix(x, 3.0))
    s = sprinkle_slab(300, _rng(2), 0.45)
    assert np.array_equal(f.causal_blocks(s, "slab", block=64), slab_causal_matrix(s))
    c = causal_matrix(x)
    assert np.array_equal(f.between_blocks(c, block=70), between_counts(c))


def test_flat_diamond_mean_r_near_reference() -> None:
    f = _f6()
    ref = f.ref_rows(reps=300, ms=(10,))[0]
    assert abs(ref["r_plano"] - 0.1) < 0.01  # r(4) = 1/10 exacto para m puntos uniformes
    c = f.causal_blocks(sprinkle_diamond(1500, 4, _rng(3)), "mink")
    out = f.measure_order(c, _rng(4))
    cell = out["por_m"]["10"]
    assert cell["local"] and cell["medible"] and cell["n_sel"] == 300
    assert abs(cell["r_media"] - ref["r_plano"]) < 0.015
    assert len(cell["perfil_medio"]) == 10 // 2 + 1


def test_suborder_has_tips_as_min_and_max() -> None:
    f = _f6()
    c = causal_matrix(sprinkle_diamond(400, 4, _rng(5)))
    b = between_counts(c)
    ii, jj = np.nonzero(b == 12)
    assert ii.size > 0
    for i, j in zip(ii[:20], jj[:20]):
        idx = f.interval_suborder(c, int(i), int(j))
        assert idx.size == 12 + 2 and idx[0] == i and idx[-1] == j
        sub = c[np.ix_(idx, idx)]
        assert np.all(sub[0, 1:]) and np.all(sub[:-1, -1])  # mínimo y máximo
        assert not np.any(sub[:, 0]) and not np.any(sub[-1, :])
        r, prof = f.interval_obs(c, b, int(i), int(j))
        assert 0 <= r <= 1 and prof.size == 7  # k = 0..6; el par de puntas (k = 12) queda fuera


def test_ads_generator_partial_order_in_diamond() -> None:
    f = _f6()
    x = f.ads_points(400, _rng(6))
    assert x.shape == (400, 4)
    assert np.all(np.linalg.norm(x[:, 1:], axis=1) <= 0.5 - np.abs(x[:, 0]) + 1e-12)
    z = x[:, 3] + 0.7
    assert np.all((z >= 0.2) & (z <= 1.2))
    assert np.median(z) < 0.5  # densidad ∝ z⁻⁴: concentrada cerca de z = 0.2
    c, coords = f.STRUCT["anti_de_sitter_4d"][2](400, _rng(7))
    assert coords is not None and is_partial_order(c) and np.count_nonzero(c) > 0


def test_locality_rule_excludes_large_m() -> None:
    f = _f6()
    assert f.local_ms(100) == [10]
    assert f.local_ms(99) == []
    assert f.local_ms(800) == [10, 20, 40, 80]
    assert f.local_ms(1600) == [10, 20, 40, 80, 160]
    # en una realización chica de bloque solo se miden los m locales
    c = f.causal_blocks(sprinkle_slab(1200, _rng(8), 0.45), "slab")
    out = f.measure_order(c, _rng(9))
    for m, cell in out["por_m"].items():
        assert cell["local"] == (int(m) <= 0.1 * out["m_max"])
        if not cell["local"]:
            assert not cell["medible"] and "r_media" not in cell


def test_parse_parts() -> None:
    f = _f6()
    assert f.parse_est("est_cono_p1_5_4d_8000") == ("cono_p1_5_4d", 8000, None)
    assert f.parse_est("est_bloque_16000_r3") == ("bloque", 16000, [3])
    assert f.parse_est("est_bloque_16000_r5") is None
    assert f.parse_est("est_caja_16000") is None
    assert f.parse_est("est_kasner_4000") is None
