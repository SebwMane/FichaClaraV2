# Consejo Ω, Revisión 21: revisión del preregistro de R20 (sesión 10, 2026-10-21)

- **Fuente:** síntesis colectiva y propuesta de enmienda pegadas por el usuario. Las 8 notas individuales están en rutas locales del usuario, fuera de este repositorio.
- **Estado:** R20 sigue **PREPARADO** (commit 5068497). **No se sella ni se ejecuta** hasta que el usuario firme la enmienda.

## 1. Dictamen

- **R20 no se sella ni se ejecuta con su redacción actual.** Se conserva el alcance A.
- **R19-(b) sigue INDETERMINADO.**
- **Falta antes de sellar:**
  1. corregir el riesgo de cobertura;
  2. hacer operativa la regla de visibilidad, indicando para cada familia el brazo, q, la banda y la unidad de la mediana;
  3. fijar el emparejamiento del N entre observado y control;
  4. precisar las reglas de inferencia: qué significa FALLA, cómo se hace el bootstrap, qué pasa con los L no finitos, y que la inferencia es condicional al instrumento ajustado.

## 2. Corrección aceptada: mi cifra del 4.5% describía otro modelo

- El rango 96 de 100 da una cobertura **marginal** de 96/101. La cobertura de una banda ya calibrada es aleatoria y sigue Beta(96, 5). Las 100 validaciones comparten ese radio, así que su conteo es **Beta-Binomial(100, 96, 5)**.
- **Verificación del orquestador**, recalculada de forma independiente (coincide con el Consejo):

  | Modelo | Fallo por q (< 90/100) | Fallo en alguno de los 4 q |
  |---|---|---|
  | Conforme por lote, el que describe R20 | 0.0523 | 0.1933 con q independientes; 0.2092 con la cota de unión |
  | Cobertura fija p = 0.95, binomial | 0.0115 | 0.0451 |

- **Con la regla de 90/100:**
  - para bajar el riesgo global a cerca del 5% hacen falta n_cal = 1938 por q con independencia, o 2398 con la cota de unión;
  - verificado: 0.0500 y 0.0500.
- **Mi error:** escribí «con cobertura real de 0.95, la probabilidad de que falle algún q es 0.045». Eso supone una cobertura condicional fija, que el procedimiento no garantiza. El gate es **diagnóstico** y no certifica una cobertura mínima.

## 3. Una tercera opción que el Consejo no listó (cálculo del orquestador, para decidir)

Bajar el umbral del gate **antes de los datos**, sin cambiar las muestras. Riesgo con n_cal = 100 (Beta-Binomial) y sensibilidad (probabilidad de que el gate detecte una cobertura real de 0.85 o de 0.88 en un q):

| Gate | Fallo por q | Fallo en algún q | Detecta 0.85 | Detecta 0.88 |
|---|---|---|---|---|
| 90/100 | 0.052 | 0.193 | 0.90 | 0.67 |
| 88/100 | 0.020 | 0.077 | 0.75 | 0.42 |
| **87/100** | **0.012** | **0.047** | **0.65** | **0.31** |
| 86/100 | 0.007 | 0.028 | 0.54 | 0.22 |

No es un retoque retrospectivo, porque no hay datos de R20. El costo es una sensibilidad menor del diagnóstico.

## 4. Decisiones pendientes de firma

| # | Decisión | Opciones |
|---|---|---|
| R21-1 | Riesgo de cobertura | **(a)** mantener 100 de calibración y el gate 90/100, aceptando el 19.3% (cota 20.9%); **(b)** subir la calibración a unas 1938 a 2398 por q (las referencias cuestan de 7 a 9 veces más); **(c)** mantener 100 y bajar el gate a 87/100 (4.7%, sin costo extra y con menos sensibilidad) |
| R21-2 | Visibilidad (propuesta mayoritaria del Consejo) | **Copias:** D_obs/radio(q = 1) ≥ 2.0 y D_ctrl/radio(q = 1/8) ≥ 1.5. **Regiones:** D_obs/radio(q = 1) ≥ 2.0, y con q = 1/8, tanto D_obs/radio como D_ctrl/radio con mediana ≥ 1.25. Medianas por brazo sobre 10 realizaciones de selección, con la banda plana del mismo q. Se informan cada cantidad, la semilla, la banda y N. |
| R21-3 | N del control | **(a)** emparejado: el control usa el mismo N realizado que el observado diluido, en copias y en regiones; **(b)** sorteos Binomial independientes, como en R19 |
| R21-4 | Inferencia | FALLA = una contradicción observada en el conjunto, sin afirmar control del error familiar. Bootstrap percentil de la mediana, remuestreando realizaciones completas (observado y control emparejados, los 4 q juntos), 2000 réplicas. Un L no finito excluye esa realización y se cuenta; si más del 10% de una estructura son no finitos, la estructura no es evaluable. Los IC son condicionales al instrumento ajustado. |

## 5. Recomendación del usuario (2026-10-21) y propagación al preregistro

- **El usuario recomienda aprobar** R21-1(c), R21-2 y R21-3(a), y R21-4 con una enmienda para los L no finitos:
  - tratarlos como límites extendidos (−∞ / +∞);
  - si D_obs = D_ctrl = 0 o L no está definido, la estructura queda **no evaluable**;
  - no hay exclusiones selectivas.

  También pidió precisar R21-3(a): N es el tamaño total del orden ya diluido, y el control se genera condicionado a ese N, conservando la regla de copias.
- **Verificación:** con el gate de 87/100, el fallo es de 0.0119 por q; 0.0468 con independencia y 0.0476 por la cota de unión (4 × 0.0119). Coincide.
- **Propagado** al preregistro (versión 2): §1 control emparejado, §2 cobertura corregida, §3 visibilidad por brazo y §4 inferencia y no finitos. La cifra errónea de la versión 1 (0.045) queda corregida y explicada en el texto.
- **Estado:** el texto de la recomendación se lee como recomendación, no como firma. **Antes de lanzar al guardián se pide la firma explícita** de la versión 2.
