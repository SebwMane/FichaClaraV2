# F1 y F2: resultados de las pruebas de la firma espectral (2026-10-12)

- **Especificación:** `OMEGA_F1F2_PREREGISTRO.md`. **Datos:** `results/omega_f1f2/`. Ventanas C5-A sin cambios.

## 1. F1, prueba ciega: **PASA**

| Estructura | Papel | β mediana | Rodilla mediana | ¿Compatible? | ¿Correcto? |
|---|---|---|---|---|---|
| de Sitter 2D (curva) | positiva | −1.000 | 0.48 | Sí | Sí |
| Cilindro 2D (otra topología) | positiva | −1.046 | 0.37 | Sí | Sí |
| 2-orden en 200 capas | negativa | −1.151 | 1.62 | No | Sí |
| 2-orden en 400 capas | negativa | **−1.148** | **0.70** | No | Sí |
| Orden aleatorio de dimensión 3 | negativa | −0.708 | 1.51 | No | Sí |
| Orden aleatorio de dimensión 4 | negativa | −0.563 | 3.59 | No | Sí |
| Percolación p = 0.02 | negativa | **−0.905** | **12.85** | No | Sí |
| Percolación p = 0.10 | negativa | −1.634 | 5.36 | No | Sí |

**Lectura (nivel 1 para estas 8 estructuras, N ∈ {800, 1600}).** La firma reconoció como geometría 2D un espacio-tiempo curvo y uno con espacio periódico, sin haberlos visto, y rechazó a las seis estructuras no geométricas. Es la primera prueba ciega que pasa en todo el programa.

**Márgenes que hay que mirar:**
1. **Las dos mitades del criterio fueron necesarias.** El 2-orden en 400 capas tiene la rodilla **dentro** de la ventana (0.70) y solo lo rechaza la pendiente, por 0.044. La percolación con p = 0.02 tiene la pendiente **dentro** (−0.905, en el borde) y solo la rechaza la rodilla. Ninguna de las dos mitades habría bastado sola.
2. **El imitador por capas se acerca.** Su rodilla baja con el número de capas: 2.28 (100 capas), 1.62 (200), 0.70 (400). Con más capas podría entrar en la ventana; su pendiente se mantiene en ≈ −1.15, fuera por poco.

## 2. F2, dependencia de la dimensión: **FALLA**

| d | Objetivo −2/d | β en N = 400 / 800 / 1600 / 2400 | ¿Pasa? |
|---|---|---|---|
| 2 | −1.000 | −0.994 / −0.994 / −1.005 / −1.006 | **Sí** |
| 3 | −0.667 | −0.739 / −0.714 / −0.744 / −0.766 | No: se aleja |
| 4 | −0.500 | −0.469 / −0.641 / −0.673 / −0.632 | No: a 0.13 y sin tendencia clara |

**Lectura.** En 2D la pendiente es −1 con gran estabilidad. En 3D y 4D **no** converge a −2/d en estos tamaños: la hipótesis de la ley de Weyl, tal como la formuló el Consejo, **no queda confirmada** para d ≥ 3. Puede ser tamaño finito (en 4D, N = 2400 da solo ≈ 7 elementos por lado) o puede que la relación no sea −2/d; no se puede distinguir aquí.

## 3. Un hallazgo lateral

Los órdenes aleatorios de dimensión 3 y 4, que **no** son lorentzianos (orden producto), tienen pendientes muy parecidas a las de los sprinklings lorentzianos de la misma dimensión (−0.71 frente a −0.74; −0.56 frente a −0.63). Lo que los separa es la rodilla. Sugiere que la pendiente mide algo como la dimensión del orden, y que la rodilla es la que distingue la estructura lorentziana (nivel 2).

## 4. Estado de la firma

- **En 2D:** pasó una prueba ciega (F1) con márgenes estrechos en dos casos, y su pendiente es estable en −1 (F2). Sube a **nivel 1 para 2D y estas estructuras**.
- **En 3D y 4D:** la relación con la dimensión no se confirmó; la firma no está lista para dimensiones mayores.
- **Riesgo abierto:** 2-órdenes con más de 400 capas.
