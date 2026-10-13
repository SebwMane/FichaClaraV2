# Consejo Ω, Revisión 5: réplica de E1′ con criterio simétrico (E1″)

- **Fecha:** 2026-10-12.
- **Estado:** PROPUESTA para el Consejo, **no congelada ni ejecutada**. Pedido del usuario: «prepara esa réplica para el Consejo».
- **Antecedente:** `OMEGA_E1PRIMA_RESULTADOS.md` (E1′: DISCRIMINA, con la asimetría señalada en §3.1).

## 1. Motivo

E1′ separó la geometría de KR, pero el criterio da más tolerancia a la familia más ruidosa: la geometría (SE = 0.059) recibió un margen de 0.168 y KR (SE = 0.009) uno de 0.068. La réplica debe comprobar si la separación sobrevive a un criterio **simétrico**.

## 2. Transparencia: lo que ya se sabe

El coeficiente de KR ya fue observado (≈ 0.42 a 0.43 en tres umbrales). Cualquier ventana fija elegida ahora alrededor de 1/3 podría ajustarse, aun sin intención, para incluir o excluir ese valor. Por eso el criterio principal **no usa ninguna ventana elegida por el orquestador**: compara a KR directamente con la geometría medida en la misma réplica.

## 3. Diseño propuesto de E1″

- **Instrumento y región:** idénticos a E1′ (SSEE, doble truncamiento con c = 1, región N/8 con control de forma).
- **Tamaños:** N_L ∈ {800, 1200, 1600, 2000, 2400}.
- **Familias y semillas:** geometría con **30 semillas por tamaño** (150 corridas, entropía 20261054); KR con 10 semillas por tamaño (50 corridas, entropía 20261055). Percolación y capas gruesas no se repiten salvo decisión R5-3: en E1′ quedaron lejos de cualquier criterio.
- **Costo estimado:** 40 a 50 minutos en tramos.

## 4. Criterio principal (simétrico, vota)

Sea a_g ± SE_g el coeficiente de la geometría y a_K ± SE_K el de KR, ajustados igual que en E1′.

1. **Validez del control:** forma logarítmica (R² log > lineal) y |a_g − 1/3| ≤ 2·SE_g.
2. **Distinción:** |a_K − a_g| > 2·√(SE_g² + SE_K²).

- **DISCRIMINA:** se cumplen 1 y 2.
- **NO DISCRIMINA:** se cumple 1 y no 2 (KR es indistinguible de la geometría con esta precisión).
- **INDETERMINADO:** no se cumple 1.

## 5. Análisis de potencia y expectativa (declarados antes de correr)

- Con 30 semillas por tamaño se espera SE_g ≈ 0.059 · √(50/150) ≈ 0.034; con SE_K ≈ 0.009, el umbral de distinción queda en ≈ 0.07.
- Si los valores de E1′ se repiten (a_g ≈ 0.375, a_K ≈ 0.43), la diferencia (≈ 0.055) **no supera** ese umbral: el resultado más probable es **NO DISCRIMINA**.
- Si el valor verdadero de la geometría está más cerca de 1/3 (V′ dio 0.529 y E1′ 0.375: la dispersión entre réplicas es grande), la diferencia podría superar el umbral.
- **Expectativa honesta:** incierta, inclinada a NO DISCRIMINA. Si eso ocurre, la conclusión correcta será que el DISCRIMINA de E1′ dependía de la asimetría del criterio.

## 6. Criterio secundario propuesto: huella espectral (E-S)

E1′ mostró que la pendiente espectral separa a las familias con claridad (geometría β = −1.00 y rodilla 0.45; KR β = −0.02; percolación −1.66; capas gruesas −1.13, con rodillas de 9 a 12). Propuesta: **geometría-compatible** si β ∈ [−1.10, −0.90] y la rodilla mediana ∈ [0.2, 1.0].

**Advertencia:** estos límites se proponen después de ver los datos de todas las familias; no sería una prueba ciega. Sirve como confirmación exploratoria con semillas nuevas, y su resultado debe leerse en nivel 2 aunque se registre como «pasa» o «no pasa».

## 7. Decisiones pendientes del usuario

| # | Pregunta |
|---|---|
| R5-1 | ¿Ejecutar E1″ tal como se propone? |
| R5-2 | ¿La huella espectral (§6) vota como criterio secundario, o solo se informa? |
| R5-3 | ¿Repetir también percolación y capas gruesas (unos 25 minutos más)? |
| R5-4 | Si E1″ da NO DISCRIMINA, ¿se rebaja el resultado de E1′ a «dependiente del criterio» en el informe final? |

## 8. Escala de afirmación

Un DISCRIMINA en E1″ sería nivel 1 para 2D, N ≤ 2400 y KR, con criterio simétrico, y daría solidez real al caso 3. Un NO DISCRIMINA sería nivel 1 sobre la insuficiencia de la entropía truncada para separar la geometría de KR con esta precisión.
