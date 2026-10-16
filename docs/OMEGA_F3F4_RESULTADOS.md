# F3 y F4: resultados (2026-10-12)

- **Especificación:** `OMEGA_F3F4_PREREGISTRO.md`. **Datos:** `results/omega_f3f4/`.

## 1. F3, firma 4D con links: **FALLA**

Ventanas calculadas por la regla, solo con la calibración geométrica: β ∈ [−0.363, −0.263]; rodilla ∈ [1.511, 6.046].

| Estructura | Papel | β_L | Rodilla | ¿Compatible? | ¿Correcto? |
|---|---|---|---|---|---|
| Diamante 4D | positiva | −0.315 | 3.05 | Sí | Sí |
| Bloque 4D sin borde espacial | positiva | **−0.263** | 3.53 | Sí, en el borde exacto de la ventana | Sí |
| de Sitter 4D (curvo) | positiva | −0.313 | 3.12 | Sí | Sí |
| KR | negativa | −0.021 | 10.9 | No | Sí |
| **Orden aleatorio de dimensión 4** | negativa | **−0.271** | **3.31** | **Sí** | **No** |
| Percolación p = 0.05 | negativa | −0.095 | 1.57 | No | Sí |
| Capas gruesas | negativa | −1.312 | 10.4 | No | Sí |

**Lectura (nivel 1).** La firma 4D reconoce las tres geometrías lorentzianas (plana, sin borde espacial y curva) y rechaza KR, la percolación y las capas gruesas. Pero **no distingue la geometría 4D del orden aleatorio de dimensión 4**, un orden producto (cada par se relaciona si sus cuatro coordenadas crecen): espacio continuo de dimensión 4, pero con conos cúbicos en vez de conos de luz. Además, el bloque sin borde quedó justo en el límite de la ventana.

## 2. F4, pendiente con links frente a la dimensión

| d | β_L con N = 400 / 800 / 1600 / 2400 | Q1, estable | β_C (comparación) |
|---|---|---|---|
| 2 | −0.122 / −0.091 / −0.047 / −0.052 | Sí | −0.995 / −1.015 / −1.003 / −1.009 |
| 3 | −0.246 / −0.225 / −0.220 / −0.218 | Sí | −0.739 / −0.705 / −0.739 / −0.763 |
| 4 | −0.240 / −0.295 / −0.309 / −0.310 | Sí | −0.507 / −0.625 / −0.680 / −0.633 |
| 5 | −0.305 / −0.307 / −0.316 / −0.338 | Sí | −0.330 / −0.510 / −0.570 / −0.626 |

**Q1 (estabilidad): sí en las cuatro dimensiones. Q2 (monotonía 3 > 4 > 5): sí.**

**Lectura (nivel 1, descriptiva).** Con links, la pendiente es estable con el tamaño y crece en magnitud con la dimensión: −0.05, −0.22, −0.31, −0.34 de 2D a 5D, con incrementos cada vez menores. Con la matriz causal, en cambio, solo 2D es estable; en 5D la pendiente deriva de −0.33 a −0.63. No se ajustó ninguna fórmula (sería post hoc).

## 3. Qué significa

1. **Lo que mide la pendiente con links (nivel 2).** Una propiedad estable del orden que depende de la dimensión de forma monótona. Es un instrumento de dimensión más confiable que el de la matriz causal fuera de 2D, aunque su relación exacta con d no se conoce.
2. **Lo que la firma no ve (nivel 1).** La diferencia entre la estructura lorentziana (conos de luz) y otras estructuras continuas de la misma dimensión. En 2D no aparece este problema porque allí el orden producto de dos coordenadas **es** el orden lorentziano (coordenadas de cono de luz). En 4D ya no lo es, y la firma no los separa.
3. **Consecuencia para el caso 3.** La firma reconoce «orden continuo de dimensión d», no «espacio-tiempo». Para reconocer espacio-tiempo hace falta otro ingrediente que detecte la estructura de los conos de luz.
