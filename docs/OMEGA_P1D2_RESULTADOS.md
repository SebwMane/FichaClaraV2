# Ω — P1-D.2: resultados (coherencia multiescala y estabilidad dimensional)

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-p1d2`.
- **Prerregistro:** `docs/OMEGA_P1D2_PRERREGISTRO.md` (§0–§6 y enmienda P1D2-A1).
- **Datos:** `results/p1d2/` (`graphs.jsonl`, `summary.json`, `level3.json`). 64 grafos con N = 4096 (40 de decisión, 24 informativos de C0).

## 1. Veredicto oficial: SEPARA-PARCIAL

Con los tres valores de δ (0.15, 0.10 y 0.25), el veredicto es el mismo:

| Medida | Resultado |
|---|---|
| Sensibilidad (G: pasan I y II) | 13/13 = 1.00 |
| Especificidad (L ∪ NL: no pasan el conjunto) | 23/23 = 1.00 |
| Familias R que pasan el conjunto | **caveman K4** |

El criterio prerregistrado exige que ninguna familia R pase para dar SEPARA. Pasa una, de modo que el veredicto es **SEPARA-PARCIAL**.

## 2. Resultado por familia (δ = 0.15)

| Clase | Familia | Nivel I | Nivel II | D_plat | ρ |
|---|---|---|---|---|---|
| G | RGG2 k12 | 3/3 | 3/3 | 1.92–1.96 | 0.20–0.23 |
| G | RGG3 k12 / k8 | 6/6 | 6/6 | 2.84–2.87 | 0.36–0.48 |
| G | Anillo k12 / toro cuadrado 64² / triangular 64² | 3/3 (CV≡0) | 3/3 | 0.99 / 1.87 / 1.87 | — |
| G | T³ 16³ | 1/1 (CV≡0) | 1/1 | **2.66** | — |
| R | **Caveman K4** | 1/1 | **1/1** | **1.04** | 0.14 |
| R | Caveman K8 / K16 | 2/2 | 0/2 (SIN_PLATEAU) | — | 0.085 / 0.064 |
| R | Retículo de cliques 3D (K8, 8³) | 1/1 | 0/1 (SIN_PLATEAU) | — | 0.50 |
| L | Retazos-8 / 27 | 0/6 | 0/6 | — | 1.13–1.35 |
| L | Retazos-64 | 1/3 | sin ventana de D | — | 0.99–1.03 |
| L | RGG3 + 1 % de atajos | 1/3 | sin ventana de D | — | 1.00–1.08 |
| L | **WS β = 0.003** | **0/3** | **3/3 («meseta» 2.6–2.8)** | 2.62–2.77 | 5.0–6.1 |
| L | WS β = 0.01 / 0.03 | 0/6 | 0/6 | — | 1.4–2.8 |
| NL | ER k12, RR k12 | sin ventana | sin ventana | — | — |
| Informativo | C0-R, N = 729 (24 grafos) | NO_HOMOGENEIZA 18, sin ventana 6 | sin ventana de D 24/24 | — | 1.64–2.53 |

### Nivel III (certificado Ω-1.1, N ≈ 729; enmienda P1D2-A1)

| Grafo | Códigos |
|---|---|
| Caveman K4 (N = 728) | **F3, F9** (3/3 claves) |
| RGG3 k12 (referencia) | ninguno: pasa |
| RGG2 k12 (referencia) | F9 |
| Anillo k12 (referencia) | F5, F9 |

El certificado, tal como está calibrado, solo valida la clase 3: rechaza también el RGG2 y el anillo, que son geometrías auténticas. Su fallo en dimensiones 1 y 2 **no** es evidencia contra la geometría. Para K4, el resultado es coherente con «cliques locales, no vecindarios de variedad», pero tampoco aporta información que RGG2 y el anillo no compartan.

## 3. Contraste con las predicciones congeladas

| # | Predicción | Resultado |
|---|---|---|
| 1 | G pasa I y II; retículos con meseta en 1, 2, 3 (tolerancia 0.25) | ✔ salvo **T³ 16³ (2.66, 0.34 por debajo de 3)** y RGG3 (2.85, error 0.15). Hay un sesgo de tamaño finito de −5 % a −11 % |
| 2 | Caveman K4, K8 y K16 pasan I y fallan II por oscilación persistente | **Parcial.** K8 y K16 ✔; **K4 pasa I y II** (su oscilación es despreciable) |
| 3 | El retículo de cliques 3D pasa I y II (meseta ≈ 3) | ✘ **Pasa I pero no II.** No se corrió su Nivel III |
| 4 | Retazos, atajos y WS no pasan I | ✔ salvo 2 grafos en el límite (ρ ≈ 0.99–1.0) |
| 5 | ER y RR sin ventana | ✔ |
| 6 | El Nivel II aislado añade poco sobre el I; solo separa caveman | ✘ **Falso.** El Nivel II por sí solo es engañado por WS β = 0.003 (sensibilidad 1.0, especificidad 20/23 = 0.87), y la separación del caveman es un artefacto (§4) |
| Global | SEPARA-PARCIAL | ✔ la etiqueta, pero con el mecanismo equivocado: yo esperaba que pasara el retículo de cliques |

El Nivel I aislado da especificidad 21/23 = 0.91, con las dos excepciones en el límite (ρ ≈ 1).

## 4. Lo que revelan las trazas (descriptivo; análisis *post hoc* anunciado en P1D2-A1)

1. **El fallo del Nivel II para K8, K16 y el retículo de cliques es un artefacto de periodicidad.** D_B oscila con el período de la decoración (caveman K8: 0.52 ↔ 1.66 ↔ 0.70 ↔ 1.51…; retículo de cliques: 2.1 ↔ 3.1). Suavizado a dos pasos, D_B2 converge a **1.02–1.05** en K4, K8 y K16 (los tres son una línea 1D a escala gruesa), y a 2.5–2.84 en el retículo de cliques, cuya ventana es demasiado corta (el toro tiene lado 8).
2. **La definición de meseta es demasiado permisiva con rampas.** WS β = 0.003 recorre D_B de 0.99 a 2.81 sin saturarse; es un cruce de 1D a expansor, no una dimensión. Una racha de 3 valores con rango relativo ≤ 0.15 no distingue una rampa de una convergencia.
3. **Las geometrías convergen monótonamente.** En RGG2, anillo, cuadrado y triangular, las diferencias de D_B2 entre escalas consecutivas se reducen (0.10, 0.06, 0.03, 0.03, 0.01, 0.00). En WS β = 0.003 se mantienen en ≈ 0.2 durante cinco escalas.
4. **El RGG3 con N = 4096 solo tiene tres valores suavizados** (2.56, 2.78, 2.90). Todavía no demuestra convergencia: con este tamaño, el RGG3 y un cruce lento no se distinguen por D_B. Hace falta más N.
5. **Los grafos de C0 con N = 729 no son evaluables en el Nivel II.** El Nivel I da NO_HOMOGENEIZA en 18/18 evaluables (ρ = 1.64–2.53), como en Fase 2.

## 5. Dictamen del cerebro

**1. Respaldo a la formulación del Consejo, con una cota que hay que declarar.**
- «Geometría ⇒ homogeneización intermedia» está respaldado: todas las geometrías evaluadas (13 grafos) la cumplen.
- Pero todas son **toros planos**, es decir, espacios homogéneos con simetría de traslación. El Nivel I mide en realidad **invariancia estadística de traslación**.
- Una geometría real inhomogénea, como una esfera con curvatura, una caja con borde o un espacio hiperbólico, tendría nodos de borde o curvatura variable, y con ellas el CV no bajaría. Por tanto: **necesaria para geometrías homogéneas; no demostrada en general.**

**2. Homogeneidad ≠ dimensión ≠ variedad queda confirmada y se concreta:**

| Estructura | Nivel I | Nivel II | Nivel III |
|---|---|---|---|
| Geometría homogénea (RGG, retículos) | pasa | pasa (D_B converge) | pasa en clase 3 |
| **Geometría gruesa con decoración periódica** (caveman, retículo de cliques) | pasa (ρ ≪ 1) | falla por oscilación, o pasa si es pequeña (K4) | falla (F3, F9): vecindarios no variedad |
| Cruce o transitorio (WS β = 0.003) | falla (ρ ≈ 5–6) | **pasa de forma espuria** | — |
| Pegado incoherente (retazos, atajos) | falla | sin meseta | — |
| Localidad C0 | falla | sin ventana | falla |

**3. La decisión que no puedo tomar yo, y que el Consejo debe tomar:** qué significa «geometría» para Ω.
- **Geometría gruesa:** un grafo cuasi-isométrico a un espacio homogéneo. Caveman, el retículo de cliques y K4 son geometrías gruesas auténticas (1D, ≈3D, 1D). Los niveles I y II las detectan correctamente, y ahí no hay falso positivo.
- **Geometría tipo variedad:** vecindarios que se parecen a los de una variedad. Es lo que exige el certificado (F9), y por eso rechaza esas estructuras.
- Si Ω aspira a lo primero, los niveles I y II son una buena prueba. Si aspira a lo segundo, el Nivel III es imprescindible. Recomiendo **mantener ambas categorías** («geometría gruesa» y «geometría de variedad») en todos los informes.
- La frase «un falso positivo del caveman» solo vale bajo el segundo criterio.

**4. No convierto nada en certificado.** El Nivel II, tal como está definido, no se debería usar por separado: lo engaña un cruce (WS) y una oscilación periódica lo hace fallar con geometrías gruesas. El Nivel I es el más fiable de los dos, pero también depende de homogeneidad y es débil en el límite (ρ ≈ 1).

## 6. Propuesta para el Consejo: P1-D.3 (no ejecutada)

Hay dos defectos concretos del Nivel II y una cota del Nivel I que se pueden atacar con un experimento pequeño con grafos nuevos y predicción congelada:

1. **Estimar D_B con BFS desde nodos muestreados** (en vez de todas las distancias), para alcanzar N = 2·10⁴–5·10⁴, donde el RGG3 sí converge. Con N = 4096 no hay datos suficientes para separar convergencia de cruce.
2. **Sustituir la meseta por un criterio de convergencia suavizado a dos pasos:** que las diferencias de D_B2 decrezcan y que la última sea ≤ 0.05. Se probaría contra el cruce WS y contra la oscilación periódica con tamaños de clique nuevos (K5, K6, K10, K12).
3. **Controles de geometría inhomogénea:** RGG en una esfera S², RGG en una caja con borde y, si es barato, un espacio hiperbólico. La pregunta es qué hace el Nivel I con una geometría real que no es invariante por traslación.
4. **Más controles de decoración:** retículos de cliques 2D y 3D con lado ≥ 12, y cruces WS con β = 0.001, 0.002 y 0.005.

Sería hipótesis generada con los mismos datos de P1-D.2; por eso necesita prerregistro propio y grafos nuevos. **No lo he ejecutado.**

## 7. Escala de afirmaciones

| Nivel | Afirmación |
|---|---|
| 1 | Cálculos exactos de los conteos de retículos (Fase 2); los generadores están verificados por tests |
| Verificado en el dominio (N = 4096, toros planos) | 13/13 geometrías homogéneas pasan I y II; 23/23 L ∪ NL no pasan ambos; el Nivel II aislado es engañado por WS β = 0.003; la oscilación de D_B aparece con decoraciones periódicas |
| 2 (explicación) | El Nivel I mide invariancia de traslación estadística; la oscilación de D_B es la huella del período de la decoración |
| 4 (extrapolación, sin demostrar) | Que «geometría ⇒ Nivel I» valga para espacios inhomogéneos; que D_B2 convergente caracterice geometrías de dimensión ≥ 3 sin N muy grande |
