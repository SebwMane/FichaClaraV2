# Ω — E6: resultados de la calibración estática (E6-S v1)

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
