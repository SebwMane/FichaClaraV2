# Ω-P3 — L-P3-0b: resultados y estado del programa

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-p3`.
- **Prerregistro:** `docs/OMEGA_P3_L0.md` §6.
- **Datos:** `results/p3_l0b/` (56 grafos: 32 del panel y 24 finales C0 informativos).
- **Reparto:** código de Sonnet, revisado antes de correr; verificación de referencias por Haiku; análisis y decisión de Opus.

## 1. Veredicto oficial: KD1 FALLA → la familia D muere; P3 queda sin familia viable

| Criterio | Resultado | Datos (mediana de κ de Ollivier; f_neg = fracción con κ < −0.1) |
|---|---|---|
| **KD1** (¿los destinos tienen signo?) | **FALLA, solo por el árbol** | Caveman +0.50 ✔. ER −0.47 y RR −0.50 ✔. T³ y cuadrado 0.000 ✔. **Árbol −0.167, −0.167 y 0.000** ✘ (exigía ≤ −0.2). Su signo es mixto: f_neg 0.51, f_pos 0.37 |
| **KD2** (¿ve el pegado incoherente?) | **PASA con margen amplio** | RGG3 f_neg = 0.35. Retazos-2³ 0.67–0.68 y retazos-3³ 0.73–0.79, umbral 0.40. Medianas: −0.22 y −0.32 frente a −0.02 |
| **KD3** (informativo) | RGG3 **≈ plano** (−0.015); RGG2 +0.09 | No favorece cristales |

**Transparencia.** El agente que escribió el código observó en pruebas previas a la corrida oficial que la mediana de κ de un árbol aleatorio es ≈ 0, porque las aristas de hoja y los tramos de camino son planos, y lo avisó antes de ejecutar. **No se enmendó KD1**, porque el fallo es sustantivo, no técnico: la planitud local no excluye árboles ni cadenas, que es justo un destino que la familia D debía evitar.

## 2. Contraste con las predicciones congeladas

| # | Predicción | Resultado |
|---|---|---|
| 1 | KD1 pasa | ✘ (árbol) |
| 2 | KD2 pasa | ✔ |
| 3 | RGG con κ ≈ +0.2…+0.4 (por los triángulos) | ✘: RGG3 −0.015, RGG2 +0.09. La geometría aleatoria es casi plana en media a escala 1, con una distribución ancha (f_neg 0.35, f_pos 0.25) |
| 4 | Finales C0 positivos | ✘: mediana +0.03, f_neg 0.25, f_pos 0.33. **Casi planos, como el RGG3** |
| 5 | WS mayormente plano o positivo, con una cola negativa pequeña | ✔: mediana +0.17, f_neg 0.01 |

## 3. Lo que cambia respecto a Ω-D

La curvatura de Ollivier a escala 1 es el **primer observable local** que separa a la vez cuatro de los cinco destinos degenerados:

| Destino | Curvatura de Ollivier |
|---|---|
| Clique (0D) | **+0.5** |
| Expansor | **−0.5** |
| Retazos | **−0.2 a −0.3**, con 70 % de cola negativa |
| Geometría (RGG3, retículos) | **≈ 0** |
| Árbol | **Mixto**: no se separa |
| Atajos al 1 % | No se distinguen del RGG3 (−0.03 a −0.06) |
| **Finales C0 (no geométricos)** | **+0.03: no se distinguen del RGG3** |

La difusión (L-ΩD-0) veía los retazos como «más coherentes» que la geometría; la curvatura los ve claramente incoherentes.

**Por qué el árbol queda mezclado (nivel 1–2).** En una arista de camino (los dos extremos de grado 2) la curvatura perezosa es 0, igual que en el ciclo; las aristas de hoja también dan 0, como encontró el agente. Solo las ramificaciones dan curvatura negativa. Un árbol aleatorio es mayoritariamente caminos y hojas, así que su mediana es ≈ 0.

## 4. Mapa negativo acumulado del programa

| Hipótesis | Prueba | Resultado | Destino observado o derivado |
|---|---|---|---|
| S0 (triángulos) | Bloque L | ✘ | Vacío o cliques; las geometrías no son equilibrio (0/1680 persistentes) |
| Ω-B | Bloque L | ✘ | Uniforme o cliques |
| C0 (competencia y saturación) | L1–L4, F1, Θ | ✘ geometría; ✔ localidad | 1D (U), expansor (E), d ≈ 2–2.5 no variedad (R); frágil a Θ ∝ 1/N |
| C1 de 4-ciclos | Análisis previo | ✘ | Dimensión como parámetro (T³ q = 4 frente a Q₆ q = 5) |
| Ω-D (difusión) | L-ΩD-0 | ✘ premisa | No ve los retazos a τ corto; la recompensa por persistencia favorece 1D; el α objetivo impone d |
| P3-A, B y C (crecimiento genérico) | L-P3-0 (analítico) | ✘ | Mundo pequeño (T1), 1D o mundo pequeño (T2), hiperbólico o polímero ramificado (T3) |
| P3-D (planitud de Ollivier) | L-P3-0b | ✘ premisa | La planitud no excluye árboles |

**Conjetura de trabajo reforzada (nivel 4):** ninguna regla local sin un entero impuesto selecciona robustamente una geometría extendida de dimensión finita mayor que 1. La planitud es un punto crítico, por T5 (Gauss–Bonnet, demostrado).

## 5. Lo que se rescata

1. **La curvatura de Ollivier como observable de coherencia.** Separa la incoherencia de los retazos y de los expansores mejor que la difusión. **Pero no detecta la no-geometría de C0** (mediana +0.03, igual que el RGG3): su fallo es mesoscópico y la curvatura de escala 1 es local. Sería, como mucho, un complemento del Nivel I (ρ, que sí rechaza C0), nunca un sustituto. Requiere un prerregistro de validación sobre el panel de P1-D.3.
2. **El mapa negativo** (§4) es un resultado en sí mismo: siete hipótesis descartadas con criterio prerregistrado y sus mecanismos de fallo identificados.
3. **T1–T5** (L-P3-0) delimitan el espacio de reglas de crecimiento admisibles.

## 6. Decisión del cerebro y propuesta al Consejo

1. **P3 queda sin familia viable según el prerregistro.** No se ejecuta ninguna dinámica.
2. **No rescato D modificándola** (por ejemplo, «planitud + mínimo de dos padres para crear ciclos»). Aunque la exigencia de ciclos tiene base estructural (T4), formularla ahora estaría informada por estos datos. Solo podría abrirse como hipótesis nueva, D′, con prerregistro propio, datos nuevos y ratificación explícita de «planitud como ley». Además seguiría abierta la selección de dimensión: un ciclo 1D también es plano.
3. **Propongo al Consejo elegir entre:**
   - **(a) Consolidar.** Cerrar formalmente la línea de «mecanismo local» con el mapa negativo y la conjetura de §4, y preparar un informe técnico autocontenido del programa (métodos, instrumento validado, siete negativos).
   - **(b) Validar la curvatura como instrumento.** Barato: el panel de P1-D.3 ya existe. Incorporarla a la batería si pasa.
   - **(c) Abrir D′.** Solo con ratificación de la planitud como ley y prerregistro nuevo.
   
   **Mi recomendación: (a) + (b).** El programa ha llegado al punto en que su resultado más sólido es el negativo bien delimitado, y conviene fijarlo antes de abrir más hipótesis.
