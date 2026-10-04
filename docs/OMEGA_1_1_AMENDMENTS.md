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

## Pendientes de decisión del usuario (no aplicadas)

- **D1 / B4: regla de anillos.**
  - Propuesta: evaluar r ∈ [2, r_hi−1] en lugar de [2, r_hi].
  - Motivo: con RGG3 N=3000, r = r_hi da 0.80–0.875 y produce F9 por saturación del toro, mientras que T³ 15³ da 1.0.
  - Sin el cambio, el certificado 3D es inalcanzable con tamaños XL.
- **B14: homogeneidad.**
  - Propuesta: exigir n_valid ≥ 0.9·n_gigante.
  - Requiere un campo nuevo de configuración y una nueva versión del diseño.
- **D2:** confirmar R1 (tamaños XL y coste), R3 (rigor de umbrales), R4 (isotropía por MDS local) y R8 (expectativa honesta) del §6 de `docs/OMEGA_1_1_DESIGN.md`.
