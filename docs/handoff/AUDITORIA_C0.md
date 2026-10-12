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
