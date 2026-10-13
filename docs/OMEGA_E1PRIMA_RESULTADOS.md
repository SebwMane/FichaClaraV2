# E1′: resultados (caso 3 de Ω-0.2, ruta entrópica)

- **Fecha:** 2026-10-12. **Commit de ejecución:** `495d5af` (300 corridas, 6 familias × 5 tamaños × 10 semillas, entropía 20261053). **Especificación:** `OMEGA_E1PRIMA_PREREGISTRO.md`, congelado en el paso 4, con V′ aprobada y R4-1 firmada. **Datos:** `results/omega_e1p/e1_*.json` y `e1_summary.json`.

> **R5-4 (usuario, 2026-10-12): resultado rebajado a «DISCRIMINA, dependiente del criterio».** La réplica simétrica E1″ (`OMEGA_E1DOBLEPRIMA_RESULTADOS.md`) mostró que KR es indistinguible de la geometría cuando los errores se tratan igual. El veredicto original se conserva abajo sin cambios.

## 1. Resultado oficial: **DISCRIMINA**

| Familia | ¿Cumple el criterio (c = 1)? | Coeficiente a ± SE | R² log / lineal |
|---|---|---|---|
| **Geometría (control)** | **Sí** | 0.375 ± 0.059 (distancia a 1/3: 0.042; límite 0.168) | 0.46 / 0.44 |
| KR | No | 0.430 ± 0.009 (distancia 0.097; límite 0.068) | 0.98 / 0.96 |
| Percolación p = 0.05 | No | −5.97 ± 0.26 (la entropía baja con el tamaño) | 0.92 / 0.88 |
| Capas gruesas | No | 75.0 ± 2.0 (ley de volumen: el lineal ajusta mejor) | 0.97 / 0.99 |
| Bicapa | No aplicable | sin regiones extensas | |
| Sándwich | No aplicable | sin regiones extensas | |

## 2. Diagnósticos (sin voto)

| Familia | Barrido c = 0.5 / 1 / 2 (a) | Pendiente espectral β | Rodilla | Complementariedad |
|---|---|---|---|---|
| Geometría | 1.02 / 0.375 / 0.415 | −1.00 | 0.45 | 23% |
| KR | 0.418 / 0.430 / 0.421 | −0.02 | 9.74 | sin complemento |
| Percolación | −12.1 / −5.97 / −1.96 | −1.66 | 11.7 | sin complemento |
| Capas gruesas | 130 / 75.0 / 13.5 | −1.13 | 9.37 | sin complemento |

## 3. Lectura

**Nivel 1 (2D, N ≤ 2400, estas seis familias).** Con regiones centradas de tamaño N/8 y el criterio congelado, solo la geometría muestra una entropía truncada logarítmica con coeficiente compatible con 1/3. Ninguna de las cuatro familias que vencieron a la geometría en las competencias de Ω-C cumple el criterio. Por primera vez en el programa, un criterio distingue geometría sin hacerla competir.

**Lo que no hay que sobreinterpretar:**

1. **KR tiene una ley logarítmica genuina y robusta**, con coeficiente ≈ 0.42 a 0.43 para los tres umbrales (no es un artefacto del umbral, que era la sospecha de E1). Lo que la separa es el **valor** del coeficiente, no la forma. Queda excluida porque su ajuste es tan preciso (SE = 0.009) que la tolerancia 2·SE + 0.05 resulta estrecha; la geometría, más ruidosa, recibe más tolerancia. Esa asimetría es propia del criterio elegido y debe tenerse en cuenta (nivel 2).
2. **El criterio no dice nada sobre bicapa y sándwich**: quedan fuera por no tener regiones donde medir, no por su entropía.
3. **En la geometría el coeficiente depende del umbral** (1.02, 0.375, 0.415), mientras en KR no. Con regiones de N/8 el umbral c = 0.5 no reproduce el comportamiento esperado.
4. No prueba la ruta de Jacobson completa: muestra que la ley de área con coeficiente 1/3 es, entre estos rivales, exclusiva de la geometría.

## 4. Hallazgo lateral: la huella espectral

La pendiente del espectro separa a las familias con mucha más claridad que la entropía: la geometría tiene β = −1.00 con rodilla en 0.45, mientras KR tiene un espectro casi plano (β = −0.02) y las demás tienen rodillas en 9 a 12. Como descriptor fue preregistrado sin voto, así que aquí solo se registra (nivel 1 descriptivo); podría convertirse en un criterio propio en una prueba futura.

## 5. Siguientes pasos posibles (no preregistrados)

1. Repetir E1′ con otra semilla y un criterio simétrico frente al ruido (por ejemplo, una ventana fija alrededor de 1/3 calibrada en V′) para ver si la separación de KR depende de la asimetría señalada en §3.1.
2. Llevar la huella espectral (β y rodilla) a criterio preregistrado.
3. Etapa E2: justificar el truncamiento desde el propio sustrato.
4. Repetir en 4D.
