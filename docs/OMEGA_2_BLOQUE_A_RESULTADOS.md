# Ω-C: resultados del Bloque A (vía analítica de la rama lorentziana)

- **Fecha:** 2026-10-12. **Commit de ejecución:** `a6d68d1`. **Especificación:** `OMEGA_2_PROPUESTA_CAUSAL.md` §13 y §15 (detalle congelado en `4cce482`). **Datos:** `results/bloque_a/`.

## 1. Resultados oficiales

| Prueba | Resultado |
|---|---|
| **A-1** ventana universal de supresión de familias por capas, λε² mod π ∈ (π/3, 2π/3) | Derivada a orden dominante. **Nivel 2**, con respaldo externo de nivel 3 para KR (P. Carlip, S. Carlip y Surya 2022). |
| **A-2a** código (bipartito, fórmula exacta) | **PASA** (12 puntos a ≤ 3·SE) |
| **A-2b** KR con N = 8, las 65 536 configuraciones exactas | **PASA**: la diferencia con el orden dominante queda bajo la cota en los 12 puntos |
| **A-3** exponente de Var(S_ε/ε²) | **PASA**: γ = 1.41 ± 0.05 (ε = 0.25) y γ = 1.56 ± 0.13 (ε = 0.1); ambos con γ + 2·SE < 2 |
| **A-3** supuesto gaussiano en el régimen medible | **PASA** en los 18 puntos utilizables |

## 2. Lectura

**Lo que se sostiene (nivel 2 para N ≤ 1024).** La cancelación de fase de la clase geométrica crece como N^γ con γ ≈ 1.4–1.6, más lento que la de las familias por capas (∝ N²). Junto con la ventana de A-1, eso cumple la **condición necesaria** para que, dentro de la ventana, la geometría supere a las capas en la suma lorentziana. Es la misma conclusión cualitativa de Carlip y colaboradores, ahora reproducida con la acción y los instrumentos de Ω-C.

**A-2b en detalle.** Incluso con solo 8 elementos, dentro de la ventana el promedio de fase exacto de KR es 0.0002 a 0.007, frente a 0.33 a 0.60 fuera de ella (λε² = 0.3). La cota absoluta es holgada (0.2 a 6.5), como se anticipó en §15.1: confirma la derivación pero no certifica la supresión por sí sola.

**Extrapolación (nivel 4, sin voto).** Con los ajustes de A-3, la clase geométrica superaría a la bicapa en λε² = 1.2 desde N* ≈ 96 (ε = 0.25) y N* ≈ 176 (ε = 0.1).

## 3. Límites

1. **Solo familias por capas.** Otras clases no geométricas no están cubiertas.
2. **ε fijo al crecer N** (A-19). Con una escala de no localidad fija en volumen (ε ∝ 1/N) el exponente podría ser otro.
3. **El supuesto gaussiano** solo se verificó donde la fase es medible (N ≤ 64, λε² ≤ 0.1). En la ventana, donde importa, sigue siendo un supuesto.
4. Los exponentes locales son ruidosos (para ε = 0.1 van de 0.95 a 2.26 entre tamaños consecutivos); el ajuste global es el criterio preregistrado.

## 4. Estado de la compuerta

| Rama | Estado |
|---|---|
| E | Geometría localmente estable con β = b/N (L-3c-E′ PASA); equilibrio no localizado (L-3c-E″ INDETERMINADO por histéresis fuerte). |
| L | Condición necesaria cumplida frente a las familias por capas (Bloque A PASA, nivel 2). Sin prueba frente a otras familias. |

**Ninguna rama está cerrada con éxito** en el sentido de §5: ambas sobreviven sus pruebas, cada una con un hueco declarado. R3-1 sigue pendiente.
