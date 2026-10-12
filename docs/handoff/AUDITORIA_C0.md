# Auditoría C0: orquestador sobre la entrega del agente de código

- **Auditado:** commit `e31b970` (entrega del agente de código).
- **Especificación:** `docs/OMEGA_2_PROPUESTA_CAUSAL.md` §4 y `docs/handoff/HANDOFF_C0_AGENTE_CODIGO.md`.
- **Método:** lectura del código, verificación independiente en el scratchpad (entropía de auditoría 20261014, no versionada) y ejecución de smoke y piloto declarados.

## 1. Veredicto

**Apto para la corrida oficial tras correcciones pequeñas.** Los estimadores son correctos. Hubo un test vacuo con una afirmación falsa, la herramienta no evaluaba los criterios congelados y la especificación tenía dos huecos (p de la percolación y tratamiento de nan) que el agente resolvió sin consultar.

## 2. Verificaciones que pasaron

| # | Verificación | Resultado |
|---|---|---|
| V1 | E[S_BD] exacto para 2-órdenes aleatorios por integral doble, independiente del código | 1.0000 en N = 64, 256 y 1024. Confirma la fórmula S = N − 2N₀ + 4N₁ − 2N₂ y el valor de literatura |
| V2 | Medias de N₀, N₁, N₂ del código frente a las exactas (N = 256, 300 réplicas) | 1063.6 / 811.9 / 691.2 frente a 1062.0 / 812.1 / 689.7, todas dentro de 1 SE |
| V3 | C0-2 en d = 4 con muestra mayor (1000 pares) | razón 0.997 |
| V4 | Volúmenes de intervalo c_d | d = 2: τ²/2; d = 3: πτ³/12; d = 4: πτ⁴/24 |
| V5 | RNG y arquitectura | `tests/test_architecture*.py` pasan; solo Generator(PCG64) |
| V6 | mypy --strict | limpio en `omega/causality`, la herramienta y los tests |

## 3. Hallazgos y correcciones

| ID | Gravedad | Hallazgo | Corrección |
|---|---|---|---|
| A-1 | Alta | La herramienta no calculaba ningún veredicto: el `summary.json` final no decía si C0 pasaba. | Se añadió `evaluate()` con los criterios congelados, sin cambiarlos. |
| A-2 | Alta (potencia) | sd(S_BD) ≈ 322 en N = 1024. Con 200 semillas, 3·SE ≈ 68: C0-3 casi no puede fallar. | Enmienda C-A2: criterio adicional, más estricto, sobre N₀, N₁, N₂ frente a sus esperanzas exactas (`omega/causality/exact.py`). |
| A-3 | Media | `test_controls_are_partial_orders_and_layered` terminaba en `or True` y su afirmación era falsa (un par capa 1 a capa 3 de KR tiene unos N/8 elementos entre ellos). | Se reemplazó por la propiedad correcta: altura ≤ 3, es decir, C³ = 0. |
| A-4 | Media | El handoff no fijaba p de la percolación; el agente eligió 0.05. | Enmienda C-A1: se congela y se declara. |
| A-5 | Media | Sin regla para nan en C0-4; en KR d_mid es siempre nan. | Enmienda C-A3, declarada como capaz de crear un positivo. |
| A-6 | Baja | En C0-1, N = 256 y N = 1024 usan el mismo flujo aleatorio por celda. Correlaciona celdas, no sesga. | Se deja; se anota. |
| A-7 | Hallazgo científico | El piloto muestra que la percolación con p = 0.05 tiene d_MM ≈ d_mid ≈ 1.2: C0-4 tal como está diseñado no la distingue. | No se toca el criterio. Se preregistra el fracaso esperado (nota de §8) y se planifica C0-4b. |

## 4. Lo que el orquestador NO hizo

- No cambió ninguna tolerancia de C0-1, C0-2 ni C0-3.
- No eligió τ para salvar la percolación (ningún τ razonable la salva).
- No ejecutó ninguna dinámica sobre órdenes.

## 5. Auditoría de C0-4b (entrega `e348899`)

**Veredicto: apto para la corrida oficial sin correcciones.**

| # | Verificación | Resultado |
|---|---|---|
| V7 | La calibración de δ no usa ningún control (inspección del código de `calibrate`) | cumple |
| V8 | Constantes de la herramienta frente a §9 (entropías, 40 réplicas, 20 semillas, 1.5 × p95, 90% y 10%, K = 8) | todas coinciden |
| V9 | Referencia d = 2 exacta coincide con el sprinkling (ya verificado en C0-3 exacto) | cumple |
| V10 | Flujos aleatorios de referencia, calibración y corrida oficial en entropías distintas | cumple |
| V11 | mypy --strict y tests | limpio; 15 tests pasan |

**Observación sin corrección.** En el smoke (N = 128, solo 4 réplicas de referencia) el sprinkling 4D quedó marcado en 2 de 3 semillas. Con N y réplicas tan pequeñas no es evidencia de nada y no motiva ningún cambio; se anota para comparar con la corrida oficial.

## 6. Auditoría de L-2c (entrega `4497916`)

**Veredicto: apto tras tres correcciones.**

| ID | Gravedad | Hallazgo | Corrección |
|---|---|---|---|
| A-9 | Baja | δ₂ de C0-4b copiado a mano en la herramienta (coincidía). | Se lee del `summary.json` oficial de C0-4b. |
| A-10 | Baja | Si φ̂ = 1, ln(−ln φ̂) lanza un error de dominio y aborta L-iii. | Enmienda L2c-A1. |
| A-11 | Media | Un ΔS_E = nan se leía como «geometría NO preferida». | Enmienda L2c-A3: INDETERMINADO. |
| A-12 | Sin corrección | Entropía de L-i elegida por el agente. | Enmienda L2c-A2, declarada. |
| A-13 | Observación | En el smoke (N = 128) ningún 2-orden aleatorio quedó marcado por C0-4b: δ₂ se calibró en N = 1024 y no vale en N = 128. La corrida oficial usa N = 1024, donde sí vale. Confirma que C0-4b depende del tamaño (ya anotado en sus límites). | — |

## 7. Auditoría de C1 y L-3c-E (entrega del agente tras `d3570d1`)

**Veredicto: apto tras una enmienda.**

| ID | Gravedad | Hallazgo | Corrección |
|---|---|---|---|
| A-14 | Media | §11.6 no pedía validar el clasificador recalibrado a N = 32 y 64. | Enmienda L3cE-A1 (solo puede bloquear). Scratchpad: N = 32, δ₂ = 1.59, 98% positivos y 0% capas; N = 64, δ₂ = 0.59, 99% y 0%. |
| A-15 | Baja | La herramienta tolera cambios en `results/` al exigir árbol limpio, porque escribe por partes. | Se acepta; `final` verifica que todas las partes provienen del mismo commit. |
| A-16 | Observación | Arrancar siempre desde geometría mide estabilidad y puede sobrestimar β_c si hay histéresis (transición de primer orden). Como ambos N se miden igual, la comparación sigue siendo justa; el valor absoluto de β_c no se interpreta. | — |
| — | Verificado | El agente detectó y corrigió por su cuenta una colisión de semillas entre las dos mitades de N = 64. | — |

## 8. Auditoría de L-3c-E′ (entrega `6abac2f`)

**Veredicto: apto tras una corrección.**

| ID | Gravedad | Hallazgo | Corrección |
|---|---|---|---|
| A-17 | Media | Efecto piso: si la geometría ya se pierde en el b más bajo, b_c queda fijado en 2.0 y dos tamaños con piso producirían un PASA falso. | Enmienda L3cE′-A1. |
| — | Verificado | Calibración y validación del clasificador idénticas a §11.6 con las entropías de §12; δ determinista entre tramos; comparación «un paso» con tolerancia numérica correcta. | — |

## 9. Auditoría de L-3c-E″ (entrega tras `122bf99`)

**Veredicto: apto sin correcciones.** δ del clasificador idéntico al de L-3c-E′ en N = 32 y 64 (mismas entropías); cálculo de brechas verificado con casos sintéticos (0/0/1 y 1/2/6 pasos); bicapa inicial correcta (N = 8: 16 relaciones, S = −24 = N − N²/2).
