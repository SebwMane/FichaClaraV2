# Ω-1.1 — Registro de enmiendas al preregistro (antes del paso 6)

Este registro existe para cumplir M§44: no ajustar el modelo tras ver resultados. Cada enmienda indica el texto original, el nuevo, qué la motivó y su justificación.

**Estado:** a la fecha de este registro **no se ha ejecutado ninguna corrida full**. Todas las enmiendas se originan en ejecuciones smoke o en la auditoría `docs/OMEGA_1_1_AUDIT.md`.

**Invariante:** ningún umbral de decisión de `omega/config/settings11.py` (de `CertificateThresholds`, `LocalStructureConfig`, `TopologyConfig`, `CurvatureConfig`, `WeylConfig` ni `DistanceSuiteConfig`) cambia por estas enmiendas. Única excepción: `resistance_max_nodes` (A-7), que es un tope de coste, no un umbral.

**Congelación:** el `config_hash` se recalcula y se congela con el commit que cierra estas enmiendas, antes de ejecutar O-00 en modo full.

- Fecha: 2026-10-04.
- Base: rama `claude/hopeful-galileo-88u1l1`, auditoría sobre `4b036d7`.

| ID | Experimento / módulo | Original | Nuevo | Motivación | Justificación |
|---|---|---|---|---|---|
| A-1 | O-08 | Identidad de Euler χ(Betti)=χ(conteos) siempre | Solo si `counts[-1]==0` | Smoke G2 | Matemática: con Betti truncado a dim 2 y tetraedros presentes, la identidad no aplica. |
| A-2 | O-05 | Validación Gibbs para Langevin y Metropolis | Afirmada solo para Metropolis; Langevin informa `bias` | Smoke G2 y medición de la auditoría | Langevin sin corrección tiene sesgo O(dt): +21.8 SE, y la mitad con dt/2. Metropolis es exacto (≤0.51 SE). MALA se pospone a Ω-1.2. |
| A-3 | O-06 smoke | N=(24,32,48) | N=(16,20,24); la expectativa de S0 se preregistra para N≥64 | Coste smoke | N≤24 está fuera de la malla preregistrada; el modo full no cambia. |
| A-4 | O-02 | Integrador sigmoide hasta convergencia | `max_steps` acotado (1500 smoke / 20000 full); se informa MAX_STEPS | Coste | Es un informe, no una decisión. |
| A-5 | s06 (B6) | χ = 0 exacto para S0 determinista | `CHI_ZERO_TOL = 1e-12` | Smoke: peso F0 ~1e-10 ⇒ χ~4e-21 y λ* espurio | Precisión numérica, no física. |
| A-6 | `evidence_cfg` (B1, B2) | Muestreo de Ollivier reducido a 25 aristas si ρ_bin>0.1 (aplicado en O-05..O-11) | Reducción solo si F1 está garantizado (ρ_bin ≥ dense_rho o ⟨W⟩ ≥ dense_meanw), aplicada uniformemente en todos los experimentos | Auditoría | En 0.1<ρ<0.5 la SE de la cola (~0.034) superaba el umbral 0.03. |
| A-7 | settings11 (B3) | `resistance_max_nodes = 2000` | 4000 | Auditoría: el control RGG3 con N=3000 daba F5 por omisión | Tope de coste, no umbral. |
| A-8 | certificate `size_robust` (B9) | Todos los tamaños limpios | Sufijo contiguo de tamaños limpios que incluye el mayor, con longitud ≥ min_sizes | Auditoría | Coherencia con R1: en 3D no hay ventana a N pequeño. No es cherry-picking: el sufijo debe ser contiguo e incluir el N mayor. |
| A-9 | certificate `null_separated` (B10) | Separación ≥ 3σ_null | ≥ max(3σ_null, dim_tol) | Auditoría | Evita "separar" con σ diminuta. Es más estricto. |
| A-10 | certificate códigos por punto (B5) | F6, F7 y F8 siempre que falte el campo | F6, F7 y F8 solo si el resultado modal es PASS | Auditoría | Diagnóstico correcto. El veredicto no cambia. |
| A-11 | O-05 (B17) | `GIBBS_SE_FLOOR = 0.005` | Eliminado | Auditoría | Era un ensanchamiento no autorizado de 3·SE. |
| A-12 | O-05 (B16) | Región dominante con cualquier motor | Exige confirmación con Metropolis; si no, `langevin_only_unconfirmed` | Auditoría | Langevin está sesgado (A-2). |
| A-13 | gate (B7) | Prerrequisitos desde O-00 | Añade `p01_analytical` y `p02_golden`, más commit == HEAD y árbol limpio | Auditoría y P§22 | Integridad del orden del panel. |
| A-14 | O-04 (desv. 7) | Full solo con N=100 | Añade N=200 para FULL y NO_TRIANGLES | Auditoría | La afirmación T=0 (P§3) debe valer al N de s06. |
| A-15 | Diseño §0/§1.2/§1.6 | β1^(4)(T³ 4³)=3 y κ=0 | Requiere lado ≥5 (se usan 6³ y 9³) | WP-C | Error del diseño: con lado 4, los ciclos de envoltura miden 4 y quedan rellenados. |
| A-16 | Diseño O-00 | ER k12 da F4 | Primario F3, con F4 ∈ codes | G1 | Precedencia F3>F4 del diseño. |

| A-17 | `annulus_connectivity` (D1/B4) | r ∈ [2, r_hi] | r ∈ [2, r_hi−1] (si r_hi > r_min); se registran `excluded_radius = r_hi`, `evaluated_radii` y el diagnóstico NO decisorio de sensibilidad [2,r_hi−1] vs [2,r_hi−2] | Auditoría: RGG3 N=3000 daba 0.80–0.875 en r=r_hi (F9) mientras T³ 15³ da 1.0 | **Aprobada por el usuario (D1).** r_hi = radio contaminado por saturación finita/topológica. Ninguna tolerancia ni criterio de certificación cambia; la regla se aplica por igual a positivos y nulos (tests). Medido: RGG3 N=3000 min 0.906/0.938/0.969 (3/3 ok); WS 0.25, árbol y anillo 0.0 (siguen fallando). |

**Nota de alcance de A-17:** se aplica a la prueba de conectividad de anillos (lo que proponía D1). El ajuste volumétrico de D_eff es código congelado de Ω-1.0 y ya excluye la saturación con su propio corte (N̄(r) ≤ 0.2·N_gc); no se modificó.

## Decisiones del usuario (2026-10-04)
- D1: **aprobada** (A-17).
- Ejecución: **aprobada** para los pasos 1–5 completos sobre un commit limpio, con manifiesto y hash actualizados, sin cambios de parámetros después de iniciar las corridas. **No** se ejecuta s06 ni los ensembles antes de congelar y analizar los pasos 1–5.
- **Criterio preregistrado:** un resultado F0/F1 de S0 es un resultado científico válido del baseline, no un fallo que justifique modificar retrospectivamente la funcional.

## Congelación para los pasos 1–5
Script: `tools/run_steps_1_5.py` (configuraciones idénticas a los `test_full_run`, `master_entropy=20240901`).

| Paso | N | Réplicas | config_hash (sha256 canónico) |
|---|---|---|---|
| p00 | 200 | 3 | `6a178e19b1b838936a18464876ecdcf62d8598dea4c38b9a70390f1032e67e11` |
| O-00 | 200 | 3 | `6a178e19b1b838936a18464876ecdcf62d8598dea4c38b9a70390f1032e67e11` |
| O-01 | 200 | 10 | `a72304f8ab551f38dfa9996ec6f8d254f2ced4240f31027add30d216fee21726` |
| O-02 | 200 | 5 | `a60f804833d5fb0e6e87937f4fd1e4f1b40fc920a3923e4a8a7b46ca5cf7d8ac` |
| O-03 | 800 | 5 | `474194734022c90796bb04e812d48a755e60217cd1fd19d797f3bee79bd3ac6b` |
| O-04 | 100 | 10 | `4090c3ab83638cf4029d7b29328c8f42d52d9b728bb131f9150631bb3487b88c` |

Cada experimento deriva su configuración (experiment_id, N y réplicas) y registra su propio `config_hash` en `summary.json`. El commit de ejecución queda en `code_commit` de cada resumen.

## Pendientes (no aplicadas)
- **B14 (homogeneidad, n_valid ≥ 0.9·n_gigante):** requiere un campo nuevo y una nueva versión del diseño.
- **D2:** confirmar R1, R3, R4 y R8 del §6 del diseño antes de s06.
