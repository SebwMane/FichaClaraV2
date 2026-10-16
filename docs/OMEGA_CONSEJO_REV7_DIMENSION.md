# Consejo Ω, Revisión 7: por qué la pendiente no sigue a la dimensión en 3D y 4D (preregistro, 2026-10-12)

- **Pedido del usuario:** que el Consejo investigue y se hagan pruebas para entender el fallo de F2.
- **Antecedente:** F2 (`OMEGA_F1F2_RESULTADOS.md`): β = −1.006 en 2D, pero −0.766 en 3D (objetivo −0.667) y −0.632 en 4D (objetivo −0.5), sin convergencia.

## 1. Hipótesis del Consejo

**H1, propagador equivocado (físico teórico).** iΔ se construyó con la matriz causal C, que es el propagador sin masa en 2D. En 4D el propagador retardado sin masa vive solo sobre el cono de luz, y su análogo discreto es la **matriz de links** L (Johnston 2008, 2010; revisión de Surya 2019). En 3D no hay un propagador establecido. Si H1 es cierta, la pendiente medida en 3D y 4D corresponde a otro operador, y no hay razón para que siga −2/d.

**H2, tamaño finito o tramo de ajuste (estadístico).** El ajuste usa los rangos 2 a max(6, ⌈√N/4⌉). Si el régimen de ley de potencias todavía no se estableció, β dependerá del tramo elegido y de N.

**H3, borde de la región (escéptico).** El diamante tiene bordes nulos y esquinas. En 4D una fracción grande de los elementos está cerca del borde; los modos podrían estar dominados por él.

## 2. Pruebas preregistradas

Instrumento común: autovalores positivos de iΔ_X = (i/2)(X − Xᵀ), con X = C o X = L (L = pares relacionados sin elementos entre ellos). La escala del propagador (por ejemplo 1/(2π√6) en 4D) no cambia β y no se aplica. 5 semillas por celda.

| Prueba | Qué se mide | Entropía |
|---|---|---|
| **P-H1** | 4D en diamante, N ∈ {800, 1600, 2400, 3200}: β con C y con L (tramo estándar) | 20261060 |
| **P-H2** | Mismas corridas: β con tramos K ∈ {max(6, ⌈√N/8⌉), max(6, ⌈√N/4⌉), max(6, ⌈√N/2⌉)} | 20261060 |
| **P-H3** | 4D en bloque T³ × [0, 0.4] (sin borde espacial), N = 2400: β con C y con L | 20261061 |
| **D-3D** (descriptiva) | 3D en diamante, N ∈ {800, 1600, 2400}: β con C y con L | 20261062 |

## 3. Criterios (congelados)

- **H1 apoyada** si, en 4D con N = 3200, la mediana de β con L está a ≤ 0.10 de −0.5 **y** más cerca de −0.5 que la mediana con C.
- **H2 apoyada** si, para alguna matriz (C o L) y algún N, las medianas de β con los tres tramos difieren en más de 0.10 entre sí (el resultado depende del tramo, señal de régimen no asintótico).
- **H3 apoyada** si, con N = 2400 y la misma matriz, |β_bloque − β_diamante| > 0.10 para C o para L.
- **D-3D** no vota: en 3D no hay propagador establecido.
- Las tres hipótesis pueden quedar apoyadas a la vez; no son excluyentes.

## 4. Lo que no cambia

F1 y la firma en 2D no se tocan: en 2D la matriz causal es el propagador correcto.

## 5. Resultados (2026-10-12; datos en `results/omega_rev7/`)

| Hipótesis | Veredicto según §3 |
|---|---|
| H1, propagador de links en 4D da −0.5 | **No apoyada** |
| H2, la pendiente depende del tramo de ajuste | **Apoyada** |
| H3, la pendiente depende del borde | **Apoyada** |

**4D en diamante, medianas de β (objetivo −0.5):**

| N | C, K/8 | C, K/4 | C, K/2 | L, K/8 | L, K/4 | L, K/2 |
|---|---|---|---|---|---|---|
| 800 | −0.521 | −0.651 | −0.621 | −0.239 | −0.292 | −0.299 |
| 1600 | −0.484 | −0.673 | −0.653 | −0.269 | −0.316 | −0.295 |
| 2400 | −0.606 | −0.631 | −0.670 | −0.316 | −0.312 | −0.301 |
| 3200 | −0.655 | −0.623 | −0.668 | −0.343 | −0.317 | −0.307 |

**Bloque sin borde espacial frente a diamante (N = 2400, K/4):** con C, −0.314 frente a −0.631 (diferencia 0.32); con L, −0.274 frente a −0.312 (diferencia 0.04).

**3D (descriptivo, K/4):** con C, −0.704 / −0.739 / −0.763 para N = 800 / 1600 / 2400; con L, −0.214 / −0.217 / −0.221.

## 6. Lectura del Consejo

1. **Con la matriz causal (propagador 2D) en 3D y 4D, la pendiente no es una propiedad estable (nivel 1).** Cambia hasta 0.19 según el tramo de ajuste y se reduce a la mitad al pasar del diamante al bloque sin borde. F2 midió algo frágil, dominado por la forma de la región.
2. **Con la matriz de links (el propagador correcto en 4D), la pendiente sí es estable (nivel 1):** ≈ −0.31 en 4D casi sin depender de N, del tramo ni de la región; ≈ −0.22 en 3D.
3. **Pero ese valor no es −2/d (nivel 1).** Y la tendencia va al revés de la ley de Weyl: con links, la pendiente se hace más pronunciada al subir la dimensión (−0.22 en 3D, −0.31 en 4D), mientras −2/d se suaviza.
4. **Respuesta a la pregunta del usuario (nivel 2).** La pendiente no sigue a la dimensión por dos razones combinadas: (a) en 3D y 4D se estaba usando el propagador de 2D, que produce una pendiente frágil y dominada por el borde; y (b) la predicción −2/d viene de la ley de Weyl para operadores elípticos (espacios euclídeos), y el operador de Pauli-Jordan es lorentziano: cuenta los modos de otra forma, en particular a lo largo de las direcciones nulas. Para este operador no hay una ley de Weyl establecida, y los datos no apoyan −2/d.
5. **Consecuencia para la firma.** En 2D la firma usa el propagador correcto y no cambia. Para extenderla a 4D hay que usar links y recalibrar: la «firma 4D» candidata sería β_L ≈ −0.31, estable y casi insensible al borde. Si separa la geometría 4D de sus rivales es una pregunta nueva.

## 7. Siguientes pasos propuestos (no preregistrados)

- **F3:** firma 4D con links frente a rivales (KR, órdenes en capas, órdenes aleatorios de dimensión 4), con ventanas derivadas solo de la geometría.
- **F4:** medir β_L también en 2D y 3D con varios tamaños, para ver si existe una relación entre β_L y la dimensión distinta de −2/d.
