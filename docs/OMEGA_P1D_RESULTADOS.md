# Ω — P1-D: resultados de la preprueba diagnóstica de crecimiento de bolas

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-p1-diagnostico`.
- **Prerregistro:** `docs/OMEGA_P1D_PRERREGISTRO.md` (§1–§6 y enmienda P1-A1).
- **Datos:** `results/p1_diagnostic/` (N ≤ 1728, 127 grafos) y `results/p1_diagnostic_ext/` (N = 4096, 20 grafos).

## 1. Veredictos oficiales

| Panel | Ventana | Veredicto | Sensibilidad | Especificidad |
|---|---|---|---|---|
| Principal (N ≤ 1728) | N/4 (primaria) | INDETERMINADO | 15/15 | 3/4: solo 4 grafos L/NL evaluables (< 5) |
| Principal | N/2 (secundaria) | MUERE | 21/21 | 0/11: la saturación homogeneiza todo |
| **Unión principal + extensión (P1-A1)** | **N/4** | **PARCIAL** | **21/21 = 1.00** | **17/20 = 0.85** |
| Secundario C0 (N = 343) | N/4 | INDETERMINADO | 3 G evaluables | 38/39 |
| Secundario C0 | N/2 | SOBREVIVE | 6/6 | 55/59 |

**Veredicto de P1-D: PARCIAL.** La firma detecta toda la geometría evaluada y separa casi todos los controles no geométricos, pero no alcanza la especificidad ≥ 0.9.

## 2. Por familia (ventana N/4; ρ = CV(r_max)/CV(2))

| Familia | Etiqueta | Estado | ρ |
|---|---|---|---|
| RGG3 k8 / k12 / diluido (N = 729) | G | HOMOGENEIZA 9/9 | 0.65–0.75 |
| RGG3 k12 (N = 1728 / 4096) | G | HOMOGENEIZA 6/6 | 0.54–0.66 / 0.35–0.42 |
| RGG2 k12 (N = 729 / 4096) | G | HOMOGENEIZA 6/6 | 0.52–0.70 / 0.24–0.29 |
| RGG3 k16, RGG4 k12 (N = 729) | G | SIN_VENTANA | — |
| Retículos (T³, cuadrado, triangular, anillo) | G | CV≡0 | — |
| **Retazos-8 / 27 (N = 4096)** | L | **NO_HOMOGENEIZA 6/6** | 1.31–1.35 / 1.11–1.13 |
| Retazos-64 (N = 4096) | L | NO_HOMOGENEIZA 2/3 | 0.997–1.025 (en el límite) |
| RGG3 + 1 % de atajos (N = 4096) | L | NO_HOMOGENEIZA 3/3 | 1.03–1.04 (en el límite) |
| Watts–Strogatz β = 0.05 (N = 729) / β = 0.01 (N = 4096) | L | NO_HOMOGENEIZA 6/6 | 1.01–1.13 / 2.5–2.6 |
| **Caveman K8** (N = 728 / 4096) | L | **HOMOGENEIZA 2/2** | 0.085 |
| Retazos, atajos 10 %, WS β = 0.2, ER, RR (N = 729) | L / NL | SIN_VENTANA | — |
| C0 U / E / R (N = 343) | L | NO_HOMOGENEIZA 38/39 evaluables | 1.04–2.6 |

## 3. Contraste con la predicción del cerebro

Predije especificidad baja: creía que la mezcla entre bloques autopromediaría las bolas de los retazos, y el veredicto sería PARCIAL o MUERE por esa vía.

**Esa predicción falló en su mecanismo.** Los retazos **no** homogeneizan: su CV sube (por ejemplo, retazos-8: 0.27 → 0.48 → 0.64). El veredicto quedó en PARCIAL por otros motivos, el caveman y un retazo de bloques pequeños en el límite, no por el que predije.

**Explicación propuesta (nivel 2).** En una geometría, todos los nodos ven estadísticamente el mismo entorno a cualquier escala, y la ley de los grandes números concentra |B_r|. En un pegado incoherente, el acceso a las conexiones de largo alcance depende de dónde está el nodo (cerca o lejos de una frontera de bloque, de un atajo), y esa desigualdad se amplifica con r.

La firma mide, por tanto, **invariancia estadística de traslación a escala intermedia**. Es una propiedad necesaria de una geometría homogénea y no fija la dimensión.

## 4. Límites

- **Caveman (falso positivo).** Un anillo de cliques es, a escala gruesa, un retículo 1D casi simétrico, y la firma no distingue «cuasi-retículo de cliques» de geometría. Es una limitación real: detecta homogeneidad, no carácter de variedad.
- **Separación estrecha cerca de la geometría.** Con un 1 % de atajos, ρ ≈ 1.03; con bloques pequeños, ρ ≈ 1.0. Las perturbaciones débiles de una geometría quedan en la frontera.
- **Separación que crece con N.** El ρ de las RGG baja con N (0.75 → 0.6 → 0.38 para N = 729, 1728 y 4096), mientras el de los controles se queda ≥ 1.
- **Tamaño mínimo.** Con N ≤ 729, los grafos con mezcla rápida no tienen ventana. Para esos grafos, el discriminador es la velocidad de crecimiento (D_B, mundo pequeño), que ya mide el certificado.

**Observación *post hoc* (no prerregistrada, solo hipótesis).**
- En las geometrías, D_B(r) forma una **meseta**: RGG3 a N = 4096 → 2.87, 2.95, 2.99; RGG2 → 1.93, 1.96, 1.98.
- En los falsos positivos y en los casos límite no la forma: el caveman oscila entre 0.5 y 1.7; los atajos suben de 2.0 a 4.1; los retazos-8 llegan a 4.5.
- Una firma conjunta («CV decreciente + meseta de D_B») podría corregir el caveman. Solo puede afirmarse en una prueba nueva prerregistrada con controles nuevos.

## 5. Dictamen del cerebro

1. **P1 como diagnóstico: útil, pero insuficiente por sí solo.** Distingue la geometría de la localidad no geométrica en casi todos los casos, incluidos los controles diseñados para engañarla (retazos), cosa que no esperaba. No basta como certificado: falla con cuasi-retículos de cliques y es frágil cerca de la frontera.
2. **Se queda como observador permanente** (`omega/diagnostics/ball_growth.py`), asignado al punto **B** del preflight geométrico. No entra en el certificado ni en la energía, según el mandato del Consejo.
3. **La pregunta puente queda abierta**, como pidió el Consejo: ¿qué mecanismo dinámico podría producir espontáneamente la invariancia estadística de traslación a escala intermedia, sin imponerla? No propongo funcional. Señalo que es una propiedad **global** de coherencia, y que los mecanismos puramente locales probados (S0, C0) no la generan.
4. **Propuesta opcional y barata para el Consejo:** prerregistrar la firma conjunta (CV↓ + meseta de D_B) con controles nuevos (anillos de cliques de varios tamaños, retículos decorados, retazos con otras geometrías) y comparar ρ con su nulo de grados fijos al mismo N.
