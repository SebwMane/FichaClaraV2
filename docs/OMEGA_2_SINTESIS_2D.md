# Ω-C en dos dimensiones: síntesis del ciclo (2026-10-11 a 2026-10-12)

Documento de cierre para quien retome el trabajo. Detalle en `OMEGA_2_PROPUESTA_CAUSAL.md` (§1 a §18) y en los `OMEGA_2_*_RESULTADOS.md`.

## Hipótesis

Ω-C: el estado fundamental es un orden parcial finito (sin coordenadas ni métrica) con la acción discreta de Benincasa-Dowker. Dos ramas con criterios separados (R3-2): euclídea (e^{−βS}, 2-órdenes, β = b/N por R3-5) y lorentziana (e^{iλS}, todos los órdenes, acción suavizada por R3-4 con ε = K/N por R3-6). Dimensión 2 por defecto (R3-3, sin firma).

## Instrumentos validados (bloque C0, C1)

| Instrumento | Estado |
|---|---|
| Dimensión de Myrheim-Meyer (d = 2, 3, 4) | Validado (error ≤ 0.013) |
| Número = volumen | Validado (razones 0.997 a 1.014) |
| Acción BD 2D y su media exacta 1 | Validado, también contra valores exactos de N_k |
| Discriminador de geometría por perfil de intervalos (C0-4b) | Validado frente a 4 familias no geométricas en N = 1024; válido en N = 32 a 128 tras recalibrar |
| Acción suavizada y su media exacta | Validada |
| Clasificador de dos estimadores (C0-4) | **Retirado**: la percolación lo engaña |

## Resultados de la compuerta

| Prueba | Resultado |
|---|---|
| L-2c-E | Geometría no es el mínimo de acción (esperado) |
| L-3c-E | Con β fijo la ventana geométrica se angosta (β_c·N ≈ 4.5) |
| L-3c-E′ | Con β = b/N la geometría es localmente estable: PASA |
| L-3c-E″ | Histéresis fuerte y creciente: INDETERMINADO |
| L-3c-EQ | Wang-Landau no alcanza la fase de baja energía: INDETERMINADO |
| L-2c-L | Acción sin suavizar: la geometría se cancela ya en N = 16 |
| Bloque A | La geometría supera a las familias por capas grandes: PASA (nivel 2) |
| Bloque A-S | El sándwich casi no se cancela: NO CONCLUYENTE, lectura adversa |
| L-ε | Ningún escalamiento sirve en 2D: INCOMPATIBLE |

## Lo aprendido

1. **Contar relaciones de orden basta para medir dimensión, volumen y curvatura media** en sprinklings, y un perfil completo de intervalos separa geometría de las familias no geométricas probadas.
2. **El mínimo de la acción en 2-órdenes es un sándwich** (capas a, 2, a; S = −N²/2 − N + 6), no la bicapa, y esa fase es vidriosa: la dinámica local no la alcanza.
3. **Rama euclídea:** la geometría es estable con β = b/N, pero no se pudo decidir si es la fase de equilibrio. Requiere movimientos no locales o métodos de población.
4. **Rama lorentziana en 2D: cerrada.** Cancelar a las familias por capas exige una fase por link de orden 1, y en esas condiciones la acción de la propia geometría se dispersa cada vez más con N.

## Estado de R3-1

Ω-C no se adopta ni se descarta. En d = 2: rama L cerrada, rama E abierta con un problema de método. La pregunta física natural siguiente es d = 4 (R3-3), donde la acción no es topológica y donde se formulan los resultados de supresión de la literatura, con un coste computacional mucho mayor.

---

# Anexo: fase 4D (2026-10-12)

| Prueba | Resultado |
|---|---|
| C4 acción 4D y suavizada (referencia independiente por pares) | ÉXITO |
| L4 compatibilidad lorentziana en diamantes | INCOMPATIBLE (7 esquemas) |
| C4b instrumento en bloque T³ × [0, T] | ÉXITO |
| L4b compatibilidad lorentziana sin borde espacial | INCOMPATIBLE; η_int de 0.84 a 1.48, no √N |

**Lección 5.** La dispersión de la acción discreta en geometrías planas 4D crece más rápido que √N incluso lejos de los bordes: las fluctuaciones están correlacionadas a larga distancia por la no localidad. En los esquemas probados, eso impide que la geometría sea coherente en la suma lorentziana donde las capas se cancelan, en 2D y en 4D.

**Estado global de Ω-C.** Rama lorentziana: cerrada en 2D e incompatible en 4D en todo lo probado. Rama euclídea: estable localmente en 2D, equilibrio sin resolver; no iniciada en 4D (exige Metropolis sobre órdenes generales).
