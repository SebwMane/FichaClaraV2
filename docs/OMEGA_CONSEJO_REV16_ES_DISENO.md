# Consejo Ω, Revisión 16: diseño de la prueba de emergencia y selección (ES) con dinámicas de Ω (2026-10-19)

- **Mandato:** R15-3 firmada. Diseñar la prueba de H_E (emergencia efectiva) y H_sel (selección de fase) con una dinámica de Ω; **primero la ruta entrópica de Ω-E y luego la acción de Benincasa-Dowker de Ω-C**.
- **Estado:** DISEÑO. Nada preregistrado ni ejecutado. Las decisiones de §6 son del usuario.
- **Requisitos (Rev. 15, enmienda 5):**
  - un ensamble o dinámica de Ω;
  - un mapa de grano grueso justificado por el modelo;
  - una ventana de escalas;
  - observables prospectivos;
  - controles sellados de familias nuevas;
  - y, para H_sel, **medir el peso o la probabilidad de la fase bajo la dinámica**.

## 1. Hallazgo previo: Ω-E no tiene dinámica

En el repositorio, Ω-E es un **criterio de reconocimiento** (la entropía de entrelazamiento de Sorkin-Yazdi con el vacío de Sorkin-Johnston, y su «ley de área»), no una medida sobre órdenes. Para probar H_sel con Ω-E hay que **definir** la dinámica entrópica. Esa definición cambia la hipótesis y la firma el usuario (R16-1).

**Propuesta del Consejo (nivel 2):** la medida entrópica de Ω-E sobre órdenes de N elementos

  μ_λ(C) ∝ exp(λ · S_SJ(C)),

donde S_SJ(C) es la entropía de entrelazamiento de Sorkin-Yazdi, con doble truncamiento, de la región intrínseca de E1′ (`interval_balanced`, |R| ≈ N/4), exactamente el instrumento de Ω-E.
- **λ = 0** es la medida de conteo puro (la «entropía de Boltzmann» de los microestados, OMEGA_MASTER, línea 24).
- **λ > 0** favorece los órdenes con más entrelazamiento; **λ < 0**, los que tienen menos.
- **Justificación:** es la única forma de convertir el objeto de Ω-E en una dinámica sin agregar ingredientes externos. La ruta de Jacobson pide que el sustrato aporte una entropía, y aquí esa entropía pesa directamente a los microestados.

**Advertencias declaradas antes de nada:**
- **Nivel 3 (Kleitman y Rothschild):** con λ = 0 en órdenes generales, casi todos los órdenes son de tipo KR. Se espera que la fase geométrica pese casi 0. Ese es el **problema de la entropía** de los conjuntos causales, y ES-E lo pone a prueba con la entropía de Ω-E como contrapeso.
- **El instrumento es de 2D** (campo sin masa en 2D, G_R = C/2 y umbral √n/4π). Aplicarlo a órdenes generales es una convención, como ya declaró E1.
- **Imponer el espacio de 2-órdenes impone la dimensión:** un 2-orden uniforme ya es un sprinkling 2D. El método del proyecto prohíbe imponer la dimensión (OMEGA_MASTER, línea 34). Por eso ES-E en 2-órdenes es **calibración del instrumento**, no evidencia de H_sel; la evidencia debe venir de órdenes generales.

## 2. ES-E: diseño con la dinámica entrópica de Ω-E

### 2.1 Espacios de estados

| Espacio | Movimientos | Papel |
|---|---|---|
| E-a: 2-órdenes de N elementos | intercambio en la segunda permutación (`mcmc.py`) | **calibración**: con λ = 0 la muestra es un sprinkling 2D y la «fase geométrica» debe pesar casi 1; prueba que el pipeline mide bien |
| E-b: órdenes parciales generales | quitar un link o agregar una relación que deje el conjunto transitivo, con etiquetado natural fijo (generador de movimientos de `mcmc_general.py`, con el peso nuevo) | **evidencia**: aquí se mide H_sel |

### 2.2 Muestreo

- **Metropolis** con aceptación min(1, exp(λ ΔS_SJ)). Costo medido (piloto declarado, semilla 820001, sin evidencia): de 1 a 22 ms por evaluación de S_SJ con N de 64 a 256. Con N = 128, una cadena de 10⁵ pasos dura unos 11 minutos. **Cuello de botella declarado:** S_SJ no tiene actualización incremental, así que cada propuesta la recalcula completa. Con N = 128, un barrido de N(N−1)/2 propuestas tarda unos 54 s, y 10⁵ propuestas son apenas unos 12 barridos. La mezcla, no el costo por evaluación, decide la N alcanzable. ES-0 lo mide.
- **Cadenas por (λ, N):** 4 independientes desde inicios distintos (orden vacío, cadena total, sprinkling 2D y KR). Diagnóstico de mezcla: R̂ < 1.05 en S_SJ y en los observables, tiempo de autocorrelación integrado y tamaño efectivo de muestra ≥ 200.
- **Lección de Ω-C:** si las cadenas no mezclan (la histéresis de L-3c-E″ o la fase vidriosa del sándwich), el resultado es **INDETERMINADO**, no un veredicto. Se preregistra Wang-Landau en S_SJ como método alternativo.
- **Rejilla:** λ ∈ {−4, −2, −1, 0, 1, 2, 4} (unidades de 1/S_SJ típica, que se fijan en la calibración E-a) y N ∈ {32, 64, 128} (256 si el costo lo permite).

### 2.3 Fase geométrica, H_sel y H_E

- **Fase geométrica (definición prospectiva, congelada antes de muestrear):** un orden está en la fase si sus observables intrínsecos (los de R14-B v2 §4 y los de R14-C v2 que sobrevivan a la verificación) caen en la banda de equivalencia de sprinklings de **alguna** dimensión d ∈ {2, 3, 4} con N igual, con bandas calibradas por validación cruzada. Es una variable que la dinámica no conoce (no se usa en el peso).
- **H_sel (lo que se mide):** P_λ,N(fase geométrica) bajo μ_λ, con su intervalo por cadenas. La pregunta es cómo escala con N: si crece o se mantiene, la fase es seleccionada; si cae a 0, no lo es (KR gana).
- **H_E (dentro de las muestras):** a los órdenes muestreados en la fase se les aplica el mapa Φ_q de R14-B v2 para ver si sus observables convergen. A los que están fuera, para ver si alguno entra en banda al hacer grano grueso. Esto separa «la dinámica produce geometría microscópica» (H_L) de «produce algo que se ve geométrico a gran escala» (H_E).
- **Controles sellados:** un guardián nuevo define familias de órdenes no geométricos y geométricos que **no** salen de la dinámica, para verificar que la definición de fase no se deja engañar (falsos positivos) antes de medir P. Se abre después de congelar la definición de fase.

### 2.4 Predicciones congeladas (nivel 2 o 3)

- **E-a con λ = 0:** P(fase, d = 2) cercano a 1, por construcción (calibración).
- **E-b con λ = 0:** P(fase) cercano a 0 y cayendo con N (Kleitman y Rothschild).
- **E-b con λ ≠ 0: sin predicción.** Es la pregunta abierta de ES-E: ¿existe algún λ en el que la entropía de entrelazamiento de Ω-E compense la entropía de conteo de KR?

## 3. ES-C: diseño con la acción de Benincasa-Dowker de Ω-C (después de ES-E)

- **Medida:**
  - **rama euclídea:** μ_β(C) ∝ exp(−β S_BD(C)), con β = b/N (R3-5);
  - **rama lorentziana:** el peso e^{iS} no es una probabilidad. Para H_sel se informa el «peso» |Σ e^{iS}| por fase como en Ω-C, con su problema de signo declarado, o se deja para después.
- **Infraestructura existente:** `mcmc.py` (2-órdenes, 2D), `mcmc_general.py` (órdenes generales, ΔS₄ incremental), `wanglandau.py`.
- **Lo que ya se sabe de Ω-C (nivel 1):**
  - en 2-órdenes, el mínimo de S es un sándwich vidrioso;
  - P4 (euclídea 4D) sin ventana;
  - la rama lorentziana es incompatible en 2D y 4D con los criterios de entonces.
  - ES-C **no repite** esas pruebas: mide P(fase geométrica) con la definición nueva y la escala con N, y aplica H_E a las muestras.
- **Mismo protocolo de §2** (cadenas, diagnóstico, fase, H_E, controles sellados), con su rejilla de b y N y sus entropías propias.

## 4. Fases de trabajo

1. **ES-0, factibilidad (piloto declarado, sin evidencia):** mezcla de las cadenas en E-a y E-b con N = 32 y 64, y el costo real. Decide si la rejilla es viable.
2. **ES-E1, exploratoria (sin voto):** E-a completa (calibración) y E-b con la rejilla; definición de fase congelada antes de empezar.
3. **ES-E2, confirmatoria (si ES-E1 muestra algo):** preregistro nuevo, semillas nuevas y controles sellados nuevos.
4. **ES-C:** lo mismo con la acción de Benincasa-Dowker.

## 5. Dependencias

- La definición de «fase geométrica» usa los observables de R14-B v2 y R14-C v2. **R14-C v2 se congela cuando termine la verificación R15-2** (fórmulas de cadenas, cadena más larga y Benincasa-Dowker).
- La rama lorentziana de ES-C hereda el problema de signo; no se resuelve aquí.

## 6. Decisiones del usuario (pendientes)

| # | Pregunta |
|---|---|
| R16-1 | ¿Adoptar como «dinámica entrópica de Ω-E» la medida μ_λ ∝ exp(λ·S_SJ) de §1? Cambia la hipótesis: Ω-E deja de ser solo un criterio y pasa a ser un peso. |
| R16-2 | ¿Aceptar que los 2-órdenes (E-a) sirvan solo como calibración y que la evidencia de H_sel venga de órdenes generales (E-b)? |
| R16-3 | ¿Ejecutar ES-0 (factibilidad: mezcla y costo, piloto declarado) antes de preregistrar ES-E1? |
| R16-4 | Para ES-C, ¿solo la rama euclídea (β = b/N), o también el «peso» lorentziano con su problema de signo declarado? |

## 7. Firmas del usuario y verificación del acta de la sesión 3 (2026-10-19)

**Firmadas R16-1 a R16-4**, con los matices del Consejo (`docs/consejo_rev16/ACTA_SESION3.md`), y orden de ejecutar ES-0.

| # | Decisión firmada |
|---|---|
| R16-1 | μ_λ ∝ exp(λ·S_SJ) se adopta como **ensamble exploratorio nuevo** del proyecto, no como dinámica física derivada de Ω-E. Extender la ley de área de Sorkin y Yazdi (sprinklings 1+1D, truncamiento espectral) a órdenes generales es hipótesis del proyecto. |
| R16-2 | E-a (2-órdenes) solo calibra; E-b (órdenes generales) es la evidencia. **Medida de etiquetado declarada:** E-b muestrea órdenes con etiquetado natural fijo, así que con λ = 0 cada orden sin etiquetas pesa según su número de extensiones lineales, no uniforme. El control KR de ES-E1 se genera con esa misma medida. |
| R16-3 | ES-0 antes de preregistrar ES-E1, solo costo y mezcla, con un muestreador nuevo que use S_SJ (`mcmc_general.py` está escrito para la acción BD 4D). |
| R16-4 | ES-C empieza por la rama euclídea. La lorentziana queda como estudio posterior independiente. |

**Verificación del orquestador:**
- **«R15-2 sigue pendiente» está desactualizado.** R15-2 cerró en el commit a7d1fd9 (`docs/consejo_rev15/VERIFICACION_R15_2.md`), y R14-C v2 quedó preparada con las fórmulas verificadas. El Consejo deliberó con una copia anterior.
- **El cálculo de barridos del Consejo coincide** con §2.2: N(N−1)/2 = 8128 propuestas por barrido en N = 128, y 10⁵ propuestas son 12,3 barridos.
- **Corrección a §1 de este documento (error mío):** escribí que la región de E1′ tiene |R| ≈ N/4. E1′ congeló **N/8 ± 5%** (`OMEGA_E1PRIMA_PREREGISTRO.md` §2), porque KR no tiene intervalos de N/4. ES-0 usa N/8, que es el instrumento real. Consecuencia: con N ∈ {32, 64, 128}, la región tiene unos 4, 8 y 16 elementos. ES-0 mide si S_SJ tiene sentido con regiones tan chicas.
