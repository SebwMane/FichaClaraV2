# Ω — E6: resultados de la calibración estática (E6-S v1 y v1.1)

Prerregistro: `docs/OMEGA_E6_PRERREGISTRO.md` (commit a1f741f), congelado antes del código.
Código: `tools/e6_static.py`, `tests/test_e6_static.py` (11/11 pasan). Agente Sonnet; revisión del cerebro.
Datos: `results/e6/{static_table.json, summary.json}`. Tiempo de ejecución: 50 s.

## 1. Veredicto

**E6-S v1: INVÁLIDO.** V1 se cumple, V2 se cumple y V3 falla.

Por la regla congelada (§4), no se ajusta ningún umbral y E6-S v1 no puede usarse como puerta. Su sustituto es una fase nueva, E6-S v1.1 (`docs/OMEGA_E6_1_PRERREGISTRO.md`).

## 2. Tabla m(d)

| Regla | d = 1 | 2 | 3 | 4 | 5 | 6 | Clase |
|---|---|---|---|---|---|---|---|
| C-TRIV | 1 | 1 | 1 | 1 | 1 | 1 | CIEGO |
| C-SQ | 1 | 1 | 1 | 1 | 1 | 1 | CIEGO |
| C-TRI | 1 | 1 | 1 | .996 | .987 | .974 | CIEGO |
| P0-LAT(1) | 0 | .004 | .004 | .004 | .005 | .004 | CIEGO (R vacío) |
| P0-LAT(2) | 1 | 1 | 1 | 1 | 1 | 1 | CIEGO |
| **P0-LAT(3)** | 0 | **1** | .003 | .004 | .006 | .010 | **SELECTOR {2}** |
| P0-DEG(2) | 1 | .010 | .011 | .012 | .010 | .010 | MONÓTONO {1} |
| P0-DEG(3) | .007 | 1 | .028 | .028 | .025 | .029 | SELECTOR {2} |
| P0-DEG(4) | .017 | 1 | 1 | .058 | .059 | .060 | SELECTOR {2, 3} |
| P0-DEG(6) | 1 | 1 | 1 | .118 | .124 | .125 | MONÓTONO {1, 2, 3} |
| P0-DEG(8) | .118 | .140 | 1 | 1 | .143 | .135 | SELECTOR {3, 4} |
| P0-DEG(10) | .127 | .099 | .106 | .266 | 1 | .100 | SELECTOR {5} |
| P0-DEG(12) | .095 | .065 | 1 | .061 | .278 | 1 | SELECTOR {3, 6} |
| P0-SQV(4) | .001 | 1 | .015 | .019 | .023 | .033 | SELECTOR {2} |
| P0-SQV(12) | .002 | 1 | 1 | .016 | .018 | .022 | SELECTOR {2, 3} |
| P0-SQV(24) | .003 | .007 | .015 | 1 | .013 | .013 | SELECTOR {4} |
| P0-SQV(40) | .003 | .006 | .007 | .015 | 1 | .008 | SELECTOR {5} |
| P0-SQF(2/3) | .007 | 1 | .112 | .078 | .054 | .086 | SELECTOR {2} |
| P0-SQF(4/5) | .004 | 1 | 1 | .062 | .066 | .106 | SELECTOR {2, 3} |
| P0-SQF(6/7) | .016 | .051 | 1 | 1 | .055 | .021 | SELECTOR {3, 4} |
| P0-SQF(8/9) | .013 | .034 | .070 | .073 | 1 | .012 | SELECTOR {5} |

Familias DIAL: P0-DEG, P0-SQV, P0-SQF y **P0-LAT**.

## 3. Por qué falla V3 (verificado por el cerebro)

**Hecho.** P0-LAT(3) descansa en el panal (grado 3, sin triángulos ni 4-ciclos: c₄ = 0 = 3·(3 − 3)/2), y solo en él. Lo comprobé de forma independiente: panal 100 × 100 con grado 3 en todos los vértices y traza(A³) = 0. No es un error del código.

**Diagnóstico (nivel 2).** El instrumento hizo exactamente lo que se le pidió.
- Dentro del panel, el único objeto de grado 3 sin triángulos ni cuadrados es el panal, que es 2D.
- El SELECTOR es **relativo al panel**. Hay grafos de grado 3 sin ciclos cortos en cualquier dimensión, por ejemplo la truncación de Z^d (cada vértice sustituido por un ciclo de longitud 2d) para d ≥ 3.
- Lo que E6-S v1 detectó es una **ausencia de cobertura del panel**, no una codificación de d.
- La misma causa afecta a P0-DEG(3) (SELECTOR {2}), y yo lo había predicho como SELECTOR genuino: mi predicción compartía el punto ciego del panel.

**Lección 23 (nivel 1 para el hecho, nivel 2 para la generalización).** «Codificar d» solo tiene sentido relativo a una clase de geometrías de referencia.
- Si un motivo local aparece en el panel en una sola dimensión, cualquier regla que descanse en ese motivo parece un selector.
- Es la lección 17 aplicada al juez: la ausencia de geometrías en otras d no demuestra que el motivo sea dimensional.
- Un panel de E6 debe estar **cerrado bajo operaciones genéricas que conservan d** (truncación, grafo de líneas, producto con K₂). Si no lo está, sus SELECTOR son provisionales.

**Sentido del error.** Es un falso positivo: rechaza una regla que no codifica d.
- En una puerta, ese error es el barato: un falso negativo dejaría pasar una codificación.
- Aun así, V3 se escribió para exigir especificidad y no se cumplió. El veredicto es INVÁLIDO sin atenuantes.

## 4. Lo que sí queda establecido (nivel 1)

- **Sensibilidad.** Las tres familias codificadoras conocidas (grado objetivo, cuadrados por vértice, fracción de cuadrados) se detectan como DIAL. Los 13 valores individuales coinciden con lo predicho.
- **C-SQ CIEGO.** El criterio de reposo de SQ no distingue d en el panel. El fracaso de R1-1 no se debió a codificación de d: queda confirmado R1-T3.
- **P0-LAT(2) CIEGO.** La regla que fija «red regular hipercúbica» sin fijar d descansa en Z^d para todo d. La reticularidad no es dimensión (A8).
- **C-TRI CIEGO**, siempre que el panel tenga geometrías ricas en triángulos en cada d (A9).

## 5. Contraste de predicciones del cerebro

| Predicción | Resultado |
|---|---|
| E6-S será VÁLIDO | **Falló** (V3) |
| P0-LAT con q = 1, 3: reposo casi nulo; familia no DIAL | **Falló** para q = 3 (el panal) |
| P0-DEG por valor de k (7 valores), P0-SQV y P0-SQF | Acertadas, pero P0-DEG(3) acertó por el mismo artefacto de panel |
| Riesgo de coincidencias en geometrías ajenas en P0-SQV | Ocurrió (c₄ = 12 en la triangular) sin efecto en V2 |
| C-SQ CIEGO | Acertada |

## 6. Decisión del cerebro

E6-S v1 queda inválido y congelado tal cual. El sustituto es una fase nueva, **E6-S v1.1**, con estas condiciones:
- Mismos umbrales y mismas reglas.
- Panel base intacto, más un **cierre por transporte**: truncación, grafo de líneas y producto con K₂ de cada Z^d.
- Validez con los controles anteriores **más controles nuevos nunca evaluados**, para que el rediseño no se ajuste al fallo observado.

Justificación: la causa se identificó con un mecanismo concreto, comprobable y general, el cierre del panel. No es una retirada de umbrales.

Riesgo de sesgo que se reconoce: v1.1 se diseña después de ver qué control falló. Por eso:
- (i) se añaden controles frescos, ciegos y codificadores;
- (ii) los codificadores antiguos deben seguir detectándose sobre el panel ampliado, aunque el cierre podría hacerlos fallar;
- (iii) si v1.1 pasa, su validez se registra como **más débil que una validación ciega** y se somete al Consejo.

---

## 7. E6-S v1.1: resultado

Prerregistro: `docs/OMEGA_E6_1_PRERREGISTRO.md` (commit be471f5).
Código: `tools/e6_static_v11.py`, que reutiliza v1 sin modificarlo, y `tests/test_e6_static_v11.py` (30/30 tests pasan entre v1 y v1.1).
Datos: `results/e6_1/`. Panel: 28 + 17 geometrías. Tiempo de ejecución: 97 s.

Verificaciones del cerebro:
- La re-ejecución de las reglas v1 sobre las geometrías v1 reproduce `results/e6/` con 0 discrepancias (comprobado por el agente con test).
- Comprobé de forma independiente que T(Z³) de lado 6 tiene 1296 nodos, grado 3, 0 triángulos y 0 cuadrados.

**Veredicto: E6-S v1.1 VÁLIDO** (V1′, V2′ y V3′ se cumplen). Por la forma de su diseño, la validación es **más débil que una ciega** (§6).

| Regla | d = 1 | 2 | 3 | 4 | 5 | 6 | Clase v1.1 | Clase v1 |
|---|---|---|---|---|---|---|---|---|
| C-TRIV / C-SQ / C-TRI | 1 | 1 | 1 | 1 | 1 | 1 | CIEGO | CIEGO |
| P0-LAT(2) | 1 | 1 | 1 | 1 | 1 | 1 | CIEGO | CIEGO |
| **P0-LAT(3)** | 0 | 1 | 1 | 1 | 1 | 1 | **MONÓTONO {2–6}** | SELECTOR {2} |
| **P0-DEG(3)** | 1 | 1 | 1 | 1 | 1 | 1 | **CIEGO** | SELECTOR {2} |
| P0-DEG(4) | .017 | 1 | 1 | .058 | .059 | .060 | SELECTOR {2, 3} | igual |
| P0-DEG(8) | .118 | .140 | 1 | 1 | .143 | .135 | SELECTOR {3, 4} | igual |
| P0-DEG(10) | .127 | .099 | **1** | .266 | 1 | .100 | SELECTOR {3, 5} | SELECTOR {5} |
| P0-DEG(12) | .095 | .065 | 1 | .061 | .278 | 1 | SELECTOR {3, 6} | igual |
| P0-SQV(4 / 12 / 24 / 40) | — | — | — | — | — | — | SELECTOR {2} / {2, 3} / {4} / {5} | igual |
| P0-SQF(2/3) | **1** | 1 | .112 | .078 | .054 | .086 | MONÓTONO {1, 2} | SELECTOR {2} |
| P0-SQF(8/9) | .013 | .034 | .070 | **1** | 1 | .012 | SELECTOR {4, 5} | SELECTOR {5} |
| C-5CIC (fresco) | 1 | 1 | 1 | 1 | 1 | 1 | CIEGO | — |
| C-PAR (fresco) | 1 | 1 | 1 | 1 | 1 | 1 | CIEGO | — |
| P0-S2(8 / 18 / 32 / 50) (fresco) | — | — | — | — | — | — | SELECTOR {2} / {3} / {4} / {5} | — |

Familias DIAL: P0-DEG, P0-SQV, P0-SQF y P0-S2. **P0-LAT ya no es DIAL.**

### 7.1 Lectura

1. **El cierre funcionó en la dirección de la especificidad.**
   - Los dos falsos SELECTOR de v1 (P0-LAT(3) y P0-DEG(3)) desaparecen.
   - La sensibilidad se conserva: las cuatro familias codificadoras siguen siendo DIAL.
   - El codificador fresco P0-S2, que nunca se había evaluado, se detecta con selección limpia, un d por valor.
   - Las dos reglas ciegas frescas salen CIEGO.
2. **El panel cerrado también corrige mi comprensión de los codificadores, a nivel 1:**
   - «Grado = 3» **no** codifica d: existen geometrías cúbicas en todo d, incluida la escalera para d = 1.
   - El grado codifica d solo dentro de la clase de redes rígidas, lo que coincide con el enunciado de Ω3 en R1-0.
   - «Grado = 10» ahora selecciona {3, 5}: L(Z³) tiene grado 10. Sigue siendo SELECTOR, pero de un conjunto no contiguo.
3. **Los SELECTOR siguen siendo relativos al panel.**
   - El cierre usa solo tres operaciones.
   - Es esperable que un panel más rico convierta otros SELECTOR en no selectores, como pasó con P0-DEG(10), que ganó d = 3.
   - Consecuencia de uso: **un FAIL de E6-S es un rechazo conservador.** Un candidato rechazado puede pedir un *transporte* explícito: construir, con operaciones genéricas que conservan d, su motivo de reposo en las dimensiones donde falta. Si lo consigue, el FAIL se revisa ante el Consejo.
   - **El PASS no tiene esa debilidad.** Si la regla no distingue d en un panel, tampoco lo distingue en ese subconjunto de un panel mayor. CIEGO es monótono bajo ampliación del panel (nivel 1, inmediato por definición de m(d) como máximo).

### 7.2 Contraste de predicciones (v1.1)

| Predicción | Resultado |
|---|---|
| E6-S v1.1 VÁLIDO | Acertada |
| P0-LAT(3) MONÓTONO {2–6} | Acertada |
| P0-DEG(3) MONÓTONO {2–6} | **Fallida** (sale CIEGO: olvidé que Z¹ □ K₂, la escalera, es cúbica) |
| Resto de reglas v1, misma clase | **Fallida en 3 reglas** (P0-DEG(10), P0-SQF(2/3), P0-SQF(8/9)); ninguna entra en V1′–V3′ |
| Sin colisiones en Z^d □ K₂ para c₄ y \|S₂\| | Acertada |

### 7.3 Estado de E6 tras esta fase

| Pieza | Estado |
|---|---|
| E6-S (brazo estático) | **VÁLIDO (v1.1)**, validación débil. Puerta utilizable: SELECTOR o DIAL ⇒ FAIL conservador, revisable por transporte |
| S0–S2 (auditoría escrita) | Definida (prerregistro E6 §3.4) |
| E6-D (brazo dinámico) | Diseño congelado (§5 del prerregistro E6). Su calibración dinámica se ejecuta con el primer candidato |
| Definición antigua de E6 en ARQ-0 («d no es función del tamaño del alfabeto») | **Sustituida** por E6 = S0–S2 + E6-S v1.1 + E6-D |

### 7.4 Decisión del cerebro y siguiente paso

- E6 queda definido y su brazo estático calibrado. Se congela como `claude/omega-e6-congelado`.
- **Siguiente paso: R3-0**, análisis sin dinámica, como R1-0. Se estudian las estructuras de orden capaces de crecimiento y coalescencia sin especificar d, y cada candidato pasa S0–S2.
- Si el candidato es una regla local sobre grafos sin etiquetas, pasa también E6-S v1.1.
- Para órdenes (posets), el panel debe ampliarse a geometrías de orden de dimensión conocida antes de usar E6-S: productos de cadenas, d-órdenes aleatorios (intersección de d órdenes lineales) y redes de Minkowski discretas como generador del panel. Esa ampliación es parte del prerregistro de R3-0, no de esta fase.
