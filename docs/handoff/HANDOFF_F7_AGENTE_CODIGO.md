# Handoff F7 para el agente de código (2026-10-14)

**Especificación única:** `docs/OMEGA_CONSEJO_REV10_F6F7.md` §8 (y §2 para el contexto). Implementa solo eso; lo no especificado se anota en el informe, no se inventa.

## Qué entregar

1. **`tools/f7.py`**, construido sobre `tools/f5prima.py` (impórtalo y reutiliza lo que sirva: `measure`, `d_spec`, `mad`, `windows`, `frw`, `de_sitter_points`, `polyhedral_cone`, `random_dorder`, generadores de negativas). Partes `ref_<d>`, `cal_<d>`, `test_<nombre>` (17 estructuras de §8.1, índices exactos), `final`.
   - Compuerta: salir si el documento no contiene `R10-1 | **Firmada` o si el árbol git no está limpio fuera de `results/`.
   - Constantes: `EXPERIMENT_ID = 319`, `ENT_REF, ENT_CAL, ENT_TEST = 20261072, 20261073, 20261074`. Salida en `results/omega_f7/`.
   - Generadores nuevos tipados: de Sitter con intervalo de η parametrizado (la versión estándar debe reproducir exactamente a `de_sitter(d)` de f5prima con la misma semilla), FRW de radiación (a ∝ η, exponente d), y cono L^p con p finito (`lp_cone(p)`); p = ∞ y p = 1 usan `polyhedral_cone("max")` y `polyhedral_cone("l1")`.
   - `final`: `summary.json` con commit, curva β_L, ventanas, δ_forma y δ por (d, N), y por estructura: papel, d verdadera, d_spec mediana, d*, β_L y rodilla medianas, firma, razón D_forma/δ_forma mediana, razón D/δ mediana, clasificación que vota, clasificación informativa con D, correcto (None para la informativa). Veredicto único según §8.
2. **Tests** en `tests/test_causality.py` (rápidos): `lp_cone(p)` da órdenes parciales con puntos en la región; para p = 1.5, 3, 6 la relación L^p queda entre la octaédrica y la cúbica sobre los mismos puntos y es monótona en p (mayor p, más relaciones); el de Sitter fuerte y el FRW de radiación dan órdenes parciales dentro del diamante; `D_forma` es invariante si se suma una constante a todos los v_k.
3. **Verificación:** `python -m pytest tests/test_causality.py` sin fallos; `python -m mypy tools/f7.py tests/test_causality.py` sin errores nuevos (hay 6 previos conocidos en torus.py, mcmc_general.py y una línea antigua del test; ver auditoría §21).

## Qué NO hacer

No ejecutar partes oficiales, no escribir en `results/omega_f7/`, no usar las entropías oficiales fuera de las constantes, no cambiar criterios, no tocar otros archivos, no hacer commit.

## Informe

Archivos tocados, decisiones no explícitas, salida resumida de pytest y mypy, tiempos medidos con una semilla propia (en particular el rechazo del de Sitter fuerte, cuya aceptación mínima es (0.2/1.2)^4 ≈ 1/1300).
