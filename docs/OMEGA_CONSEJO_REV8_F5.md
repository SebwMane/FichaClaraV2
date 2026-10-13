# Consejo Ω, Revisión 8: preparación de F5 (coherencia entre medidas de dimensión)

- **Fecha:** 2026-10-12. **Estado:** PREPARADA, no ejecutada. Pedido del usuario: «que el Consejo prepare esa prueba F5».
- **Motivo:** F3 mostró que la firma con links reconoce «orden continuo de dimensión d» pero no la estructura lorentziana: el orden producto de dimensión 4 la imita.

## 1. Idea de F5

Dos medidas de dimensión que suponen cosas distintas. **Myrheim-Meyer (d_MM)** usa la fracción de pares relacionados y **supone conos de luz**. La **dimensión espectral (d_spec)** se obtiene de β_L comparándola con la curva de calibración de F4, y no supone la forma del cono. En un espacio-tiempo deberían coincidir; en un orden de otra estructura, quizás no.

## 2. Cuenta previa del Consejo (nivel 1, sin correr nada)

- Orden producto de dimensión 4: dos puntos al azar están relacionados con probabilidad exacta 2/2⁴ = 1/8. Con r(d) de Myrheim-Meyer, r = 0.125 da **d_MM ≈ 3.7**.
- Su β_L (−0.271 en F3) cae entre las medianas de F4 para 3D (−0.218) y 4D (−0.310): **d_spec ≈ 3.6**.
- **Expectativa preregistrada: el brazo A probablemente NO separa al impostor**, porque ambas medidas le asignan casi la misma dimensión. El Consejo lo prepara igual, como aprobó el usuario, porque es la prueba directa de la idea, y propone dos brazos más.
- Para dimensiones 3 y 5 del orden producto: r = 1/4 da d_MM ≈ 2.9; r = 1/16 da d_MM ≈ 4.6.

## 3. Diseño

### 3.1 Calibración (solo geometría lorentziana, entropía 20261066)

Sprinklings en diamante con d ∈ {2, 3, 4, 5} y N ∈ {1600, 2400}, 10 semillas por celda. De ahí salen: (a) la curva β_L(d) por N, para convertir β_L en d_spec por interpolación lineal entre dimensiones vecinas; (b) la dispersión geométrica de cada brazo, para fijar sus tolerancias por regla.

### 3.2 Brazo A, aprobado por el usuario (vota)

- **Coherencia:** |d_MM − d_spec| ≤ τ_A, con τ_A = max(0.15, 3 × 1.4826 × MAD de (d_MM − d_spec) en la geometría 4D de calibración).
- **Clasificación «espacio-tiempo 4D»:** pasa la firma de F3 (sus ventanas, sin cambios) **y** es coherente.

### 3.3 Brazo B propuesto (Myrheim-Meyer frente a punto medio)

La dimensión de punto medio (instrumento de C0-4) mide cómo escala el tamaño de un intervalo al partirlo por la mitad. En el orden producto un intervalo es una caja, y su punto medio da exactamente d (4 en el caso del impostor), mientras d_MM da ≈ 3.7: brecha ≈ 0.3. Riesgo conocido: en C0-4 la dimensión de punto medio tuvo sesgo en 4D (≈ 4.2 con N = 1024), así que la brecha geométrica podría ser de un tamaño parecido.

### 3.4 Brazo C propuesto (perfil de intervalos frente a la referencia lorentziana)

El discriminador C0-4b (validado en N = 1024) compara la cantidad de intervalos de cada tamaño con la de un sprinkling lorentziano de la dimensión estimada. Esa distribución depende de la forma del cono, así que es el candidato más directo para detectar un cono cúbico. Ventana por la regla de §9 de Rev. 3, recalibrada en la geometría 4D de calibración.

### 3.5 Prueba (entropía 20261067; N ∈ {1600, 2400}, 5 semillas por N)

| Estructura | Papel |
|---|---|
| Diamante 4D, bloque 4D sin borde espacial, de Sitter 4D | positivas |
| Diamante 3D | positiva (control de otra dimensión; se clasifica contra 3D) |
| Orden producto de dimensión 3, 4 y 5 | negativas |
| KR, percolación p = 0.05, capas gruesas | negativas (ya rechazadas por la firma; se mantienen como control) |

**PASA** (por brazo) si todas las positivas se clasifican como espacio-tiempo de su dimensión y ninguna negativa. Costo estimado: 25 a 35 minutos.

## 4. Decisiones pendientes

| # | Pregunta |
|---|---|
| R8-1 | ¿Ejecutar F5? |
| R8-2 | Brazos B y C: ¿votan, o solo se informan junto al brazo A? El Consejo recomienda que **C vote** (es el más directo para detectar la forma del cono) y que **B se informe** (por el sesgo conocido del punto medio en 4D). |

## 5. Escala de afirmación

Cualquier resultado de F5 es nivel 1 para N ≤ 2400 y estas estructuras. Si algún brazo separa al orden producto, la combinación «firma + ese brazo» sería el primer criterio que reconoce estructura lorentziana en 4D sin competencia.
