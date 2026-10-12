# Ω-C: resultados de L-ε (rama lorentziana con ε = K/N)

- **Fecha:** 2026-10-12. **Commit:** `985acb7` (árbol limpio, 9.5 s). **Especificación:** `OMEGA_2_PROPUESTA_CAUSAL.md` §18. **Datos:** `results/l_eps/summary.json`.

## 1. Resultado oficial: **INCOMPATIBLE (estructural)**

| N | D = sd(S_ε/ε²), K = 4 | D, K = 8 |
|---|---|---|
| 64 | 27 | 30 |
| 128 | 58 | 61 |
| 256 | 130 | 130 |
| 512 | 359 | 376 |
| 1024 | 904 | 875 |
| **η (D ∝ N^η)** | **1.27 ± 0.04** | **1.23 ± 0.05** |

Para que la geometría fuera coherente dentro de la ventana haría falta D ≤ 0.955 y que no creciera. D es ≈ 900 en N = 1024 y crece.

## 2. Lectura

- **Nivel 1 (medición):** con la no localidad fija en volumen, la dispersión de la acción geométrica, medida en las unidades que importan dentro de la ventana, crece como N^1.25.
- **Nivel 1 (corolario §18.1):** fuera de la ventana, con λ fijo, las capas dejan de cancelarse.
- **Conclusión preregistrada:** en 2D, con la acción BD suavizada, **ningún escalamiento de (λ, ε) cumple a la vez** «cancelar las capas» y «mantener coherente la geometría». La rama lorentziana de Ω-C **queda cerrada en d = 2**.
- **Observación (nivel 4):** con ε ∝ 1/N la situación empeora respecto de ε fijo. La cancelación gaussiana de la geometría crecería como D² ∝ N^2.5, más rápido incluso que la de las capas (∝ N²).

## 3. Contexto

2D es un caso especial: allí la acción de Einstein-Hilbert es topológica (la media de la acción BD es 1 para cualquier diamante plano), así que no hay una «geometría clásica» que la fase seleccione. Los argumentos de supresión de la literatura (Loomis-Carlip y sucesores) se formulan sobre todo en dimensión 4. El cierre de la rama L vale **para d = 2**; no dice nada sobre d = 4.
