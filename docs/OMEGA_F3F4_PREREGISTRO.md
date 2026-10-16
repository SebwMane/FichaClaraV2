# F3 y F4: preregistro congelado (Consejo Rev. 7 §7, 2026-10-12)

- **Pedido del usuario:** F3 y F4 en paralelo.
- **Instrumento común:** autovalores positivos de iΔ_L = (i/2)(L − Lᵀ), con L la matriz de links (propagador sin masa en 4D, Johnston); β por ajuste en rangos 2 a max(6, ⌈√N/4⌉) y rodilla R-rodilla-2 (`knee_from_spectrum`) sobre ese espectro.

## F3: firma 4D con links (vota)

**Calibración (solo geometría 4D, entropía 20261063):** 10 sprinklings en diamante 4D (5 con N = 1600, 5 con N = 2400). Regla congelada: ventana de β = mediana ± max(0.05, 6·MAD); ventana de rodilla = [0.5, 2] × mediana.

**Prueba (entropía 20261064, 5 semillas con N = 1600 y 5 con N = 2400 por estructura):**

| Estructura | Papel |
|---|---|
| Diamante 4D (semillas nuevas) | positiva |
| Bloque 4D T³ × [0, 0.4] (sin borde espacial) | positiva |
| de Sitter 4D (coordenadas conformes, η ∈ [−1.5, −0.5], densidad ∝ 1/η⁴) | positiva |
| KR | negativa |
| Orden aleatorio de dimensión 4 (orden producto) | negativa |
| Percolación p = 0.05 | negativa |
| Capas gruesas (8 capas) | negativa |

**Compatible:** β mediano y rodilla mediana dentro de las ventanas. **PASA** si las tres positivas son compatibles y las cuatro negativas no. **FALLA** en otro caso. Se informa cada estructura.

## F4: pendiente con links frente a la dimensión (descriptiva, con dos preguntas que votan)

Sprinklings en diamante, d ∈ {2, 3, 4, 5}, N ∈ {400, 800, 1600, 2400}, 5 semillas (entropía 20261065). Se registran β_L y, como comparación, β_C.

- **Q1, estabilidad:** para cada d, |mediana β_L(2400) − mediana β_L(800)| ≤ 0.05. Vota por dimensión.
- **Q2, monotonía:** las medianas de β_L en N = 2400 son estrictamente decrecientes de d = 3 a d = 5 (pendiente más pronunciada al subir d, como sugirió Rev. 7). Vota.
- **No se ajusta ninguna fórmula con voto:** cualquier ley propuesta después de ver los datos sería post hoc; solo se informan los valores.
- Nota: en 2D el propagador correcto es la matriz causal, no los links; β_L en 2D se mide solo como referencia.
