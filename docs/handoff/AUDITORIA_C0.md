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

## 10. Auditoría del Bloque A (entrega tras `4cce482`)

**Veredicto: apto tras una corrección menor.**

| ID | Gravedad | Hallazgo | Corrección |
|---|---|---|---|
| A-18 | Baja | `phase_average` declaraba acciones enteras y el agente silenció mypy con `type: ignore` para pasar acciones suavizadas (reales). | La firma acepta `Sequence[float]`; se quitan los `ignore`. |
| A-19 | Observación de alcance | A-3 usa ε fijo al crecer N (escala de no localidad fija en unidades de la discreción). La alternativa física, escala fija en unidades de volumen (ε ∝ 1/N), no se prueba aquí. | Se declara como límite en los resultados. |
| — | Verificado | φ_lead de A-2b coincide con la derivación de §15.1 (L ~ Bin(16, 1/2)); Δ se calcula como S_ε − (εN − 2ε²L) sobre el orden con cierre transitivo completo; la tolerancia del supuesto gaussiano sigue la fórmula de §15.3. | — |

## 11. Auditoría de L-3c-EQ (entrega tras `6894e3f`)

| ID | Gravedad | Hallazgo | Corrección |
|---|---|---|---|
| A-20 | **Alta** | La ventana de energía y el criterio de cobertura de §16 suponían que la bicapa es el mínimo de S. La enumeración exacta de N = 8 lo desmiente (S_min = −34, sándwich 3-2-3). Con la ventana original el WL habría excluido la verdadera fase de baja energía. | Enmienda EQ-A1. |
| — | Verificado | Determinismo por tramos (test), normalización Σg = N!, uso del clasificador de §12, criterio de consistencia con b_c↑ y b_c↓ oficiales. | — |

Smoke de velocidad (no evidencia): 0.04 ms por paso en N = 32.

## 12. Auditoría del Bloque A-S (entrega `1af6e5e`)

**Veredicto: apto sin correcciones.**

- f(1, ε) = 1 − 3ε y f(2, ε) = 1 − 6ε + 6ε² coinciden con `smearing_weight` del instrumento validado (C1) en ε = 0.5, 0.25 y 0.1.
- El factor de etiquetas N!/(a! b! 2!) **sobrecuenta** órdenes (una misma configuración puede surgir de asignaciones de capa distintas, por ejemplo con elementos aislados), así que U sobreestima la contribución del sándwich: la comparación queda del lado conservador, como pide §17.
- La cota aplica la desigualdad triangular solo sobre α, después de sumar exactamente sobre B; la fase del lado A tiene módulo 1. Correcto.
- Smoke (no evidencia): U ≥ exacto en N = 8 y 12; 11 s por cota en N = 1024.

## 13. Auditoría de L-ε (`435bbbe`)

Apto sin correcciones: el cálculo de S_ε de la herramienta coincide con `bd_action_2d_smeared` (validado en C1) a 10⁻⁹; ε = K/N < 1 en toda la malla; constantes iguales a §18.

## 14. Auditoría de C4 y L4 (entrega `832e841`)

**Veredicto: apto tras una corrección menor.**

| ID | Gravedad | Hallazgo | Corrección |
|---|---|---|---|
| A-21 | Baja | 4·10⁶ pares generados de una vez en un contenedor de 3 GB. | Generación en tramos de 10⁶ (mismo número de pares y misma secuencia determinista). |
| — | Verificado | El método de referencia por pares, aplicado a 2D (scratchpad, entropía 9090), reproduce E[S] = 1: 0.66 ± 0.63 (N = 64) y −4.0 ± 5.8 (N = 256). | — |
| — | Verificado | Cadena de 5 da S₄ = 4; ε → 1 recupera S₄; coeficientes de f₄ en ε = 1 dan (1, −9, 16, −8). | — |
| — | Observación | En ε = 1 y N = 512 el error de la referencia puede ser grande y restarle potencia a C4; se informa. | — |

## 15. Auditoría de L4b (entrega `3886c70`)

**Veredicto: apto sin correcciones.** Imagen mínima correcta (test con separación 0.9 → 0.1); con T = 0.4 < 1/2 la relación es transitiva (test de orden parcial); la suma de s(x) reproduce exactamente S₄ y S₄,ε del instrumento validado (tests); la referencia por pares usa v relativo al volumen del bloque, que es lo que requiere `expected_action_4d`; los pares no relacionados quedan enmascarados antes de aplicar coeficientes.

## 16. Auditoría de P4 (entrega `f683d1c`)

**Veredicto: apto sin correcciones.** Actualización incremental verificada por test (tras 20 barridos con β = 0.05, B y S₄ coinciden con el recálculo completo); validez de las propuestas (solo se quita un link; solo se agrega si todo a ≺ i ya cumple a ≺ j y todo b ≻ j ya cumple i ≺ b) implica transitividad; propuesta simétrica; presupuestos y entropías iguales a §21. Riesgo anotado: con N = 64 el redondeo de d_MM a 4 puede fallar en algunos sprinklings, y la validación del clasificador lo detectaría.

### 16.1 Corrección tras la corrida (A-22)

| ID | Gravedad | Hallazgo | Corrección |
|---|---|---|---|
| A-22 | **Alta** | Aceptación con atajo «ΔS ≤ 0 acepta», incorrecto para β < 0. La auditoría de §16 no lo vio: la validación P4-0 solo probaba β = 0 y β = 0.01. | Enmienda P4-A1: criterio corregido y validación exacta con β negativo. Lección: validar el muestreador en todo el rango de parámetros que se va a usar, no solo en un punto. |

## 17. Auditoría de Ω-E, E0/E1 (entrega tras `32aae6a`)

**Veredicto: apto sin correcciones.**
- Convenciones iguales al artículo: C_xy = 1 si x ≺ y, iΔ = (i/2)(C − Cᵀ), W = parte positiva; S = Σ λ ln|λ| con λ de (iΔ_R)⁻¹ W_R en la imagen de iΔ_R.
- Test de cordura: en la región completa el vacío es puro y S = 0 (< 10⁻⁶).
- Doble truncamiento: √N_L/(4π) en la región grande (en iΔ y W) y √n_R/(4π) en la restringida, como en el artículo.
- Regla de intervalo con extremos incluidos (+2) y desempate por índices; test con una cadena.
- Smoke (no evidencia, semilla 1): N = 1600, n ≈ 400: S sin truncar 123.6 (publicado ≈ 0.32·n − 6.6 ≈ 121); truncado 2.03 (publicado ≈ 2.04). 7 s por corrida.

## 18. Auditoría de la herramienta E1′ (paso 4)

**Veredicto: apto sin correcciones.** Comprobaciones en scratchpad (semilla 4242, no evidencia, nada escrito en `results/`):
- La decisión de congelamiento N/8 cumple su propósito: con N = 800 tienen región de 100 ± 5% la geometría (95), **KR (99)**, la percolación (95) y las capas gruesas (105); bicapa y sándwich no (0), como se esperaba.
- Corrida completa de geometría: fracciones de forma 0.114 y 0.116 (ventana de V′ [0.085, 0.125]); entropías 4.83 / 1.95 / 1.03 para c = 0.5 / 1 / 2; complemento 2.50; rodilla 0.48 y β = −0.999, coherentes con el piloto 2.
- `sj_from_eig` y `knee_from_spectrum` equivalen a las funciones originales (test).
- La guarda funciona: la etapa e1 se niega a correr sin «R4-1 FIRMADA» en el preregistro.
- Identificación de extremos: el elemento sin predecesores dentro de R es a, y el sin sucesores es b.

## 19. Auditoría de la herramienta E1″

**Veredicto: apto sin correcciones.** El control geométrico usa la condición más estricta de Rev. 5 (|a_g − 1/3| ≤ 2·SE_g, sin el +0.05 de E1′); la comparación con cada familia es simétrica en los errores; las ventanas espectrales son las de C5-A; la etiqueta sigue C5-B; la guarda exige «C5 DECISIÓN: EJECUTAR». Prueba con una semilla fuera del rango oficial (99): la corrida produce región, entropía, rodilla y β.
