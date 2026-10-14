# Handoff F5′ para el agente de código (2026-10-14)

**Especificación única:** `docs/OMEGA_CONSEJO_REV9_F5PRIMA.md` §6 (congelada en el commit 1d8c770). Implementa solo eso. Si algo no está especificado, no lo inventes: déjalo anotado en tu informe.

## Qué entregar

1. **`tools/f5prima.py`**, con la misma estructura que `tools/f5.py` (léelo primero y reutiliza sus funciones: `measure`, `d_spec`, `mad`, `windows`, `links`, `spectrum`, `de_sitter_4d` como modelo, `random_dorder`, etc.).
   - Partes (`--part`): `ref_<d>` (d ∈ 2..5), `cal_<d>` (40 planas por N; espectro y d_MM/d_mid solo en semillas 0..9; para d ∈ {3, 4} también las 20 FRW, marcadas con `"familia": "frw"`), `test_<nombre>` (las 12 estructuras de §6.3), `final`.
   - Si una parte `cal_<d>` resulta demasiado lenta para un comando (límite práctico ≈ 9 minutos), puede partirse por N (`cal_<d>_<N>`); documéntalo.
   - Compuerta: salir si el documento no contiene `R9-1 | **Firmada` o si el árbol git no está limpio (excluyendo `results/`), igual que F5.
   - Constantes: `EXPERIMENT_ID = 318`; `ENT_REF, ENT_CAL, ENT_TEST = 20261069, 20261070, 20261071`; `NS = (1600, 2400)`; `DIMS = (2, 3, 4, 5)`; ventanas de F3 para 4D copiadas de F5.
   - RNG: `Generator(PCG64(SeedSequence(entropía, spawn_key=(318, familia, d_o_índice, N, semilla))))`, con familia 0 = plana, 1 = FRW (en calibración); en referencias y prueba basta `(318, d_o_índice, N, semilla)` porque cada entropía es de una sola familia.
   - Salida en `results/omega_f5prima/`. `final` escribe `summary.json` con: commit, curva β_L, ventanas por d, δ(d, N), τ_A, τ_B y, por estructura: papel, d verdadera, d_spec mediana, d*, β_L mediana, rodilla mediana, firma, razón D/δ mediana, clasificación, correcto, brechas A y B y sus clasificaciones informativas. Y el veredicto único PASA/FALLA de §6.2.
2. **Generadores nuevos** (en `tools/f5prima.py`, tipados):
   - `de_sitter(d)` → generador `(n, rng) -> matriz causal`, generalizando `de_sitter_4d` de `tools/f3f4.py` al exponente d. Para d = 4 debe dar exactamente la misma distribución que `de_sitter_4d`.
   - `frw(d)` → η = t + 1, aceptación (η/1.5)^{2d}, misma estructura de rechazo.
   - `polyhedral_cone(norm)` con norm ∈ {"max", "l1"}: puntos uniformes en {|t| ≤ 0.5, norm(x) ≤ 0.5 − |t|} (4D: 1 tiempo + 3 espacio) por rechazo en la caja; relación x ≺ y si Δt > 0 y Δt ≥ norm(Δx).
3. **Tests** en `tests/test_causality.py` (rápidos, N pequeño): los conos poliédricos son órdenes parciales (`is_partial_order`) y sus puntos están dentro de la región; `de_sitter(4)` reproduce a `de_sitter_4d` con la misma semilla; `frw(d)` y `de_sitter(d)` devuelven órdenes parciales con todos los puntos en el diamante; la relación cúbica contiene a la de Minkowski y la octaédrica está contenida en ella (verificación de la definición).
4. **Verificación antes de entregar:** `mypy --strict tools/f5prima.py` (o el `mypy` del proyecto sobre los archivos tocados) limpio; `pytest tests/test_causality.py` sin fallos. Usa `bash -c 'set -euo pipefail; ...'` y no ocultes códigos de salida con `| tail`.

## Qué NO hacer

- No ejecutes ninguna parte oficial (`ref_`, `cal_`, `test_`, `final`) ni escribas en `results/omega_f5prima/`. Las corridas oficiales las hace el orquestador tras auditar.
- No uses las entropías oficiales en pruebas; en los tests usa semillas pequeñas propias.
- No cambies ningún criterio, umbral ni estructura. No toques archivos de Ω-1.1 ni otros tools.
- No hagas commit: deja los cambios en el árbol de trabajo. El orquestador audita y hace el commit.

## Informe

Al terminar, informa: archivos tocados, decisiones que tomaste donde el texto no era explícito, salida resumida de mypy y pytest, y una estimación de tiempo por parte (medida con una sola corrida pequeña fuera de las entropías oficiales, si la hiciste).
