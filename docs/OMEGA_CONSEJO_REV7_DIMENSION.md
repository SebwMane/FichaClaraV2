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
