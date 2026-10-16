# Ω — W-1: resultados (umbrales de paseo, estático)

- Prerregistro: `docs/OMEGA_W1_PRERREGISTRO.md` (commit a41e1ed), congelado antes del código.
- Código: `tools/w1_walks.py`, `tests/test_w1_walks.py` (13/13 pasan). Agente Sonnet; revisión del cerebro.
- Datos: `results/w1/{cells.jsonl, summary.json}`. 33 celdas; 15 min.

## 1. Veredicto literal

**W1-D.** «Las leyes de paseo examinadas no proporcionan un selector dimensional resoluble dentro del dominio computacional prerregistrado» (formulación del Consejo).

**Causa exacta (verificada en `cells.jsonl`):**
- En N₁ ≈ 1.25·10⁵, Z⁵ tiene R = 12 y T_max = 144. Como T_max/16 = 9 < 10, la celda es INSUFICIENTE para las 4 leyes, y lo mismo pasa con Z⁶.
- Por la regla §3, ninguna ley queda resuelta en N₁. La estabilidad en N no se puede establecer, y ninguna ley es «RESOLUBLE y estable».
- El agente etiquetó esto como INESTABLE-N. La etiqueta correcta según el §3 es **NO RESUELTA en N₁ por insuficiencia**. El desenlace es el mismo.

**Es un W1-D por insuficiencia, no por ausencia de corte.** Es justo el caso que el Consejo pidió no leer como imposibilidad.

## 2. Datos en N₂ ≈ 10⁶ (descriptivos; no deciden)

| Geometría | d | L1 retorno | L2 colisión | L3 cruce de 2 | L4 cruce de 3 |
|---|---|---|---|---|---|
| Z¹ | 1 | CRECE .49 | CRECE .47 | CRECE .46 | CRECE .46 |
| Z², triangular | 2 | AMBIGUA .09 | AMBIGUA .07 | CRECE .77–.79 | INS (Z²) / CRECE (triangular) |
| RGG₂ | 2 | CRECE .12 | AMBIGUA (las semillas discrepan) | CRECE | CRECE |
| Z³, FCC, BCC, RGG₃ | 3 | ACOTADA | ACOTADA | CRECE .41–.44 | AMBIGUA / ACOTADA / AMBIGUA / CRECE |
| Z⁴ | 4 | ACOTADA | ACOTADA | AMBIGUA .085 | ACOTADA |
| Z⁵, Z⁶ | 5, 6 | ACOTADA | ACOTADA | ACOTADA | ACOTADA |
| peine-1 (d_g = 2) | | CRECE | CRECE | CRECE | CRECE |
| **peine-2 (d_g = 3)** | | **CRECE .49** | **CRECE .43** | CRECE | CRECE |
| cilindro (d_g = 1) | | CRECE | CRECE | CRECE | CRECE |
| Z³ diluido (d_g = 3) | | ACOTADA | ACOTADA | CRECE | AMBIGUA |
| retazos | | INS | INS | INS | INS |

Lectura descriptiva (nivel 1, solo para estas instancias a N₂; **no** cambia el veredicto):

1. **El corte depende de la ley.** L1 y L2 cortan en 2 y L3 en 4, coincidiendo con la teoría (nivel 3). Es el patrón de W1-C: la «dimensión seleccionada» depende de qué ley se mire. Corrobora W-T4 descriptivamente.
2. **El umbral no es universal.**
   - El peine-2 tiene dimensión de crecimiento 3, y L1 y L2 lo clasifican CRECE, como un objeto de d ≤ 2. Los dientes atrapan al paseo.
   - **Recurrencia no es dimensión:** mide la conductancia o la isoperimetría, no el crecimiento de bolas.
   - Si la estabilidad en N se hubiera establecido, esto habría dado W1-B⁻ o W1-C con L1 y L2 no universales.
3. **L4 es sensible a la estructura local.** Las cuatro geometrías 3D se reparten entre ACOTADA, AMBIGUA y CRECE. El cruce de tres trayectorias no ve solo d.
4. **El tamaño necesario es grande.** Para resolver d = 5 hace falta T_max/16 ≥ 10, es decir R ≥ 13. En Z⁵ eso exige N ≳ 1.6·10⁵. Los umbrales de dimensión alta solo son accesibles cerca de 10⁶ nodos, lo que confirma cuantitativamente W-T5 (Ω-2).

## 3. Tabla del Consejo (§5)

| Pregunta | L1 | L2 | L3 | L4 |
|---|---|---|---|---|
| ¿Detecta dimensión? | En N₂ corta en 2; no resuelta en N₁ | En N₂ corta en 2; no resuelta en N₁ | En N₂ corta en 4; no resuelta en N₁ | No resuelta |
| ¿Robusta a la heterogeneidad? | No (peine-2), descriptivo | No (peine-2), descriptivo | Sin discordancias en la capa 2, descriptivo | No evaluable |
| ¿Robusta a N? | No establecida (insuficiencia en N₁) | ídem | ídem | — |
| ¿Depende del número de caminantes? | Sí: 1–2 caminantes con colisión → 2; 2 trayectorias → 4 | | | 3 trayectorias: sin corte limpio |
| ¿Depende de red o continuo? | Nivel 3: para k = 2, 4 en red y 3 en el continuo | | | |
| ¿Depende de la ley? | **Sí** (2 frente a 4) | | | |
| ¿Genera o solo mide? | **Solo mide** | Solo mide | Solo mide | Solo mide |

## 4. Contraste de predicciones del cerebro

| Predicción | Resultado |
|---|---|
| Desenlace más probable W1-C | **Fallida**: W1-D (la alternativa) |
| W1-D, si ocurre, por insuficiencia en d alto en N₂ | **Fallida en el dónde**: la insuficiencia fue en N₁ (Z⁵, Z⁶). En N₂, Z⁵ y Z⁶ sí se resuelven |
| L3 NO RESUELTA | Fallida: corta en 4 en N₂ |
| Valor crítico AMBIGUO en L2 y L3 | Acertada (Z², Z⁴) |
| Valor crítico AMBIGUO en L1 | Parcial (RGG₂ CRECE) |
| L4 con 3 ∈ d̂_c | Fallida: no resuelta, y sensible a la estructura local |
| L2 no universal por el peine-1 ACOTADA | **Fallida**: el peine-1 crece a T ≤ 4·10⁴. La finitud de Krishnapur y Peres es asintótica y no se ve aquí. La no universalidad apareció en el peine-2, que no predije |

Error de diseño reconocido:
- El umbral T_max/16 ≥ 10 junto con N₁ = 1.25·10⁵ hizo que **una sola celda** (Z⁵ en N₁, a 16 pasos del mínimo) decidiera el desenlace.
- No lo anticipé, aunque Ω-2 avisaba del problema de tamaño en d alto.
- No se corrige retroactivamente.

## 5. Decisión del cerebro y propuesta al Consejo

**¿Repetir W-1 con tamaños mayores (N₁ ≈ 10⁶, N₂ ≈ 8·10⁶)? No lo recomiendo, por un argumento de valor de la información.**
- Por diseño, W-1 no puede dar W1-A (Ω-1).
- Los mejores desenlaces alcanzables son W1-C y W1-B⁻, y los datos descriptivos de N₂ ya muestran ambos rasgos: corte dependiente de la ley (2 frente a 4) y no universalidad (peine-2).
- Ninguno de los dos abre una vía de selección: C traslada la selección a la elección de ley, y B⁻ es un detector dependiente de la homogeneidad, no un selector.
- Repetir con más tamaño solo podría convertir el W1-D formal en C o B⁻. **La decisión sobre el programa de selección sería la misma.**

**Propuesta: cerrar el programa de selección dimensional de Ω (problema B) con el resultado negativo estructurado.** Texto propuesto, que reúne las formulaciones ya aprobadas:

> Entre las familias examinadas, Ω no posee un selector de dimensión. Cuando d queda determinado de forma reproducible, aparece como parámetro, recuento de direcciones independientes (anchura, rango), condición inicial o dial (L-DIM-1, R1, R3, W-0). Las rutas sin parámetro continuo, los umbrales de paseo sobre el propio grafo, (i) no son resolubles de forma estable a tamaños ≤ 10⁶ en el dominio prerregistrado (W-1: W1-D); (ii) cuando se ven, dan un corte que depende de la ley elegida; y (iii) no son universales: la recurrencia no distingue un peine de dimensión 3 de una red de dimensión ≤ 2. Medir dimensión no es generarla. Que exista un principio pregeométrico que fije la ley de interacción sigue abierto y no es demostrable dentro de Ω. **Esto no es un teorema de imposibilidad.**

**Lo que queda en pie del programa:**
- El problema A (exclusión de degenerados) tiene instrumentos validados en sus dominios: RC-3 con W5 y W6, E6-S v1.1, y el mapa de exclusiones.
- Hay 27 lecciones registradas.
- La pregunta abierta queda formulada con precisión: ¿qué principio fijaría la ley?
