# Ω-C en 4D: resultados del piloto P4 (rama euclídea sobre órdenes generales)

- **Fecha:** 2026-10-12. **Commit de ejecución:** `27e74d5`. **Especificación:** `OMEGA_2_PROPUESTA_CAUSAL.md` §21 con enmiendas P4-A1 y P4-A2. **Datos:** `results/p4/` (incluye el intento inválido de N = 64, conservado).

## 1. Recorrido honesto

1. La primera corrida de N = 64 reveló un **error del muestreador** (A-22): aceptaba siempre los movimientos que bajan S, lo cual es incorrecto con β < 0. Se corrigió (P4-A1) y se agregó una validación exacta con β = −0.5.
2. Esa validación nueva **falló** con χ² (p ≈ 0). El diagnóstico mostró que el muestreador era correcto y que el problema era la prueba: muestras muy autocorrelacionadas. Se cambió a medias por bloques (P4-A2, declarada capaz de crear un positivo) y pasó (peor desvío 3.15 SE).
3. Solo entonces se corrieron N = 64 y 128.

## 2. Resultado oficial: **SIN VENTANA** (N = 64 y N = 128)

Clasificador 4D validado en ambos N (positivos 96% y 100%, KR y bicapas 0%).

| β | N = 128: geometría | d_MM | fracción relacionada | altura | ⟨S₄⟩ |
|---|---|---|---|---|---|
| −0.1 | 0 | 2.54 | 0.33 | 6.4 | **+9 056** |
| −0.03 | 0 | 2.80 | 0.27 | 5.0 | −317 |
| 0 | 0 | 2.89 | 0.25 | 3.0 | −1 841 |
| +0.03 | 0 | 2.59 | 0.32 | 5.1 | −2 782 |
| +0.1 | 0 | 2.69 | 0.29 | 4.0 | **−13 505** |

En todos los β, partiendo de un sprinkling 4D, la geometría se pierde dentro de los primeros 100 barridos. Con β cercano a 0, la dinámica va a órdenes de pocas capas (lo que domina por entropía).

## 3. Lo inesperado (exploración declarada, no evidencia)

Con β = −0.1 la dinámica va a una fase de **acción muy alta**. En una exploración con N = 64: 4 capas de 23, 9, 24 y 8 elementos con relaciones dispersas (25%, 29% y 6% entre capas consecutivas), y N₁ = 381. Los pares con **exactamente un elemento entre ellos** valen +9 cada uno en la acción 4D, y esta estructura los multiplica: es la imagen especular del sándwich de 2D, que explotaba un coeficiente negativo.

**Lección (nivel 2).** La acción BD 4D tiene familias no geométricas que la llevan a valores extremos **en ambos sentidos**: capas densas hacia abajo, capas dispersas con capas delgadas hacia arriba. La geometría queda en una zona intermedia que ningún signo de β selecciona, y cerca de β = 0 gana la entropía de los órdenes en pocas capas. Por eso la rama euclídea sobre órdenes generales no tiene, en este piloto, ninguna ventana geométrica.

## 4. Límites

Piloto: solo N = 64 y 128, 100 barridos de termalización, una cadena por celda, medida de órdenes etiquetados naturalmente (cada orden pesa según sus extensiones lineales). Un estudio completo pediría más tamaños, cadenas y quizá otra medida; pero la pérdida de la geometría en todos los β, incluido β = 0, deja poco margen.
