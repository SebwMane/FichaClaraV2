# Ω — Auditoría (b): ¿qué información independiente añade la curvatura de Ollivier a la batería P1-D.3? — Prerregistro

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-coherencia`.
- **Mandato del Consejo:** integrar Ollivier como diagnóstico **complementario**, no como certificado, ni criterio de dimensión, ni mecanismo. La pregunta no es «¿identifica la geometría?», cosa que ya sabemos que no hace (K ≈ 0 en árboles y en C0). Es: **¿qué información independiente añade K a la batería P1-D.3?**
- **Estado:** prerregistrado antes de calcular K sobre este panel.

## 1. Datos

- Los 77 grafos del panel P1-D.3, regenerados con sus mismas claves y generadores (`tools/p1d3_panel.py::specs()` y `build()`), con N ≈ 2·10⁴ y 5·10⁴.
- Se reutilizan sus medidas de batería de `results/p1d3/graphs.jsonl`: ρ, estado del Nivel I, D_conv, estado del Nivel II, D_s y categoría.
- **Curvatura:** κ de Ollivier perezosa (idleness 0.5), coste en saltos, en 600 aristas muestreadas por grafo. Las distancias entre vecinos de x y de y (≤ 3) se calculan con BFS local, versión dispersa equivalente a `omega.curvature.ollivier.ollivier_edge`.
- **Resumen por grafo:** mediana de κ, f_neg (κ < −0.1) y f_pos (κ > 0.1).

## 2. Preguntas y criterios

**KA1 (no redundancia).** Correlación de Spearman entre la mediana de κ y ln ρ, sobre los grafos con ρ definido.

| Resultado | Condición |
|---|---|
| **No redundante** | \|corr\| < 0.7 |

**KA2 (valor añadido sobre los casos difíciles conocidos de la batería).**

Banda «K-geométrica»: mediana de κ ∈ [mín, máx] de las medianas de G-hom ± 0.02, **y** f_neg ≤ (máx f_neg de G-hom) + 0.05.

| Caso | Por qué es difícil para la batería | ¿Qué sería que K lo corrige? |
|---|---|---|
| RGG3 + 0.1 % de atajos | Pasa el Nivel I; cae en CRUCE en el Nivel II | Queda fuera de la banda |
| ER k12 | Pasa el Nivel I | Queda fuera de la banda |
| WS β = 0.001 y 0.002 | Pasan el Nivel II (falso D ≈ 3) | Quedan fuera de la banda |
| RGG3 en caja y RGG2 en caja | Falsos negativos del Nivel I | Quedan **dentro** de la banda |

- **Por familia:** «K corrige» si ≥ 2/3 de sus grafos quedan del lado correcto.
- **Control de no rotura:** fracción de G-hom dentro de la banda (debe ser 1, por construcción de la banda) y fracción de R (cliques) y L/NL restantes fuera de la banda.

**Veredicto:**

| Veredicto | Condición |
|---|---|
| **K APORTA** | KA1 no redundante **y** K corrige ≥ 1 familia difícil |
| **K REDUNDANTE** | \|corr\| ≥ 0.7 **y** no corrige ninguna |
| **K MARGINAL** | cualquier otro caso |

## 3. Predicciones del cerebro (congeladas)

1. KA1: no redundante (\|corr\| ≈ 0.3–0.6). K mide la escala 1 y ρ la escala intermedia.
2. KA2:
   - K **corrige las cajas**: un RGG con borde tiene la misma curvatura local que sin borde;
   - K **corrige ER** (κ ≈ −0.47);
   - **no** corrige WS β = 0.001 / 0.002: el anillo k12 está en G-hom con κ ≈ +0.25, y WS hereda su curvatura positiva;
   - **no** corrige los atajos al 0.1 %: el 99.9 % de las aristas son de RGG3.
3. Veredicto previsto: **K APORTA**, de forma limitada a la escala local.
