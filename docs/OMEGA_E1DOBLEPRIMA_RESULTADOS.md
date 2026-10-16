# E1″: resultados de la réplica con criterio simétrico

- **Fecha:** 2026-10-12. **Commit de ejecución:** `f1b05d0` (230 corridas). **Especificación:** `OMEGA_E1DOBLEPRIMA_PREREGISTRO.md` (Consejo Rev. 5, decisión C5: EJECUTAR). **Datos:** `results/omega_e1pp/`.

## 1. Resultado oficial: **DISCRIMINA PARCIAL (solo huella)**

### Criterio principal (entropía, simétrico): **NO DISCRIMINA**

| Familia | a ± SE | Diferencia con la geometría | Límite de distinción | ¿Coincide con la geometría? |
|---|---|---|---|---|
| Geometría (150 corridas) | 0.380 ± 0.042 | | | control válido: |a − 1/3| = 0.047 ≤ 0.084 |
| **KR** | **0.407 ± 0.009** | **0.027** | 0.086 | **Sí** |
| Percolación | −5.96 ± 0.29 | 6.34 | 0.59 | No |
| Capas gruesas | 77.1 ± 2.1 (el lineal ajusta mejor) | 76.7 | 4.20 | No |

**Con un criterio que trata igual los errores de ambas, KR es indistinguible de la geometría.** Era la expectativa declarada en Rev. 5 §5.

### Criterio secundario (huella espectral, nivel 2): **DISCRIMINA**

| Familia | β mediana | Rodilla mediana | ¿Compatible? |
|---|---|---|---|
| Geometría | −1.005 | 0.45 | Sí |
| KR | −0.024 | 9.75 | No |
| Percolación | −1.647 | 11.8 | No |
| Capas gruesas | −1.132 | 9.37 | No |

## 2. Lectura

1. **El DISCRIMINA de E1′ dependía de la asimetría del criterio (nivel 1).** Con la tolerancia simétrica, el coeficiente de KR (0.407) cae dentro del margen de la geometría (0.380). La ley logarítmica de la entropía truncada, con esta precisión, **no** separa la geometría de KR.
2. **Lo que sí se repite (nivel 1).** La geometría da un coeficiente estable entre réplicas (0.375 en E1′, 0.380 aquí) y compatible con 1/3; KR da un logaritmo muy limpio con coeficiente ≈ 0.41 a 0.43. La diferencia existe, pero es menor que la incertidumbre de la geometría.
3. **La huella espectral separa a las cuatro familias** con márgenes enormes: KR tiene un espectro casi plano y las demás tienen rodillas veinte veces más altas. Pero sus límites se fijaron conociendo esos valores (C5-A los derivó solo de la geometría, aunque no fue a ciegas), así que su resultado es **nivel 2**: confirmación con semillas nuevas, no una prueba ciega.
4. **Una debilidad del control.** En la geometría, la ventaja del ajuste logarítmico sobre el lineal es muy pequeña (R² 0.358 frente a 0.353). El control fue válido según la regla, pero la evidencia de forma logarítmica en la geometría es débil.

## 3. Consecuencia para el caso 3 de Ω-0.2

La **entropía de entrelazamiento truncada**, tal como se mide aquí, no basta para reconocer geometría: un orden de capas sin geometría la imita. Lo que sí distingue a la geometría es su **espectro** (pendiente −1 y quiebre cerca de 0.45). El caso 3 no queda cerrado, pero su forma cambia: el criterio prometedor ya no es «S ∝ área», sino una propiedad espectral del sustrato, que todavía no tiene una conexión demostrada con la ruta de Jacobson.

## 4. Pendiente

**R5-4 (usuario):** si el resultado de E1′ se rebaja a «dependiente del criterio» en el informe final. Este resultado cumple la condición prevista para esa decisión.
