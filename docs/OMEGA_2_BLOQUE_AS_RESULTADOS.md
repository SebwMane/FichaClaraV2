# Ω-C: resultados del Bloque A-S (familia sándwich en la rama lorentziana)

- **Fecha:** 2026-10-12. **Especificación:** `OMEGA_2_PROPUESTA_CAUSAL.md` §17 (preregistro `56784b7`). **Commit de ejecución:** posterior a la auditoría de `1af6e5e`. **Datos:** `results/bloque_as/`.

## 1. Resultados oficiales

| Prueba | Resultado |
|---|---|
| **AS-2** fórmula exacta de tipos frente a enumeración completa (N = 6 y 8, tres ε, dos λ) | **PASA**: acciones y sumas de fase idénticas |
| **AS-3** cota rigurosa U(N) | Calculada para N = 64 a 1024 |
| **AS-4** ¿U < ln|Z_geo| en la ventana para ambos ε en N = 1024? | **NO CONCLUYENTE**: en ningún punto |

| N = 1024 | U (cota del sándwich) | ln|Z_geo| (supuesto gaussiano) |
|---|---|---|
| ε = 0.25, λε² = 1.2 | 2 125 | −23 793 |
| ε = 0.25, λε² = π/2 | 2 128 | −45 105 |
| ε = 0.1, λε² = 1.9 | 2 125 | −118 775 |

Referencia: ln(número de configuraciones del sándwich) ≈ 2 135.

## 2. Lectura más allá del veredicto (declarada; no lo cambia)

**El sándwich casi no se cancela (nivel 1 hasta N = 24, nivel 2 hasta N = 1024).** El cálculo exacto con 80 dígitos queda entre 0.35 y 3.0 unidades por debajo de U en todos los casos con N ≤ 24, así que la cota es ajustada. Y U mismo queda solo ≈ 9 unidades bajo el número de configuraciones en N = 1024. La familia sándwich suma sus fases casi en coherencia.

**La geometría sí se cancela, y mucho (nivel 4).** Con ε fijo, la dispersión de S_ε/ε² crece como N^1.4–1.6 (A-3), y en la ventana de supresión de las capas eso da ln|Z_geo| ≪ 0.

**Consecuencia.** Con la formulación actual (ε fijo y λε² en la ventana que suprime a las capas), **la rama lorentziana no puede seleccionar geometría**: una familia pequeña pero coherente como el sándwich la supera por un margen enorme. Es la forma fuerte de la tensión que §11.5 había anticipado. El PASA de A-3 sigue siendo cierto en lo que medía (la geometría supera a las familias por capas grandes), pero se ve que esa comparación no bastaba.

## 3. Qué quedaría para la rama L

La alternativa física que A-19 dejó sin probar: **escala de no localidad fija en volumen (ε ∝ 1/N)**, que es la que usa la literatura para el límite continuo. Con ella la dispersión de la acción geométrica podría quedar acotada, y las ventanas de supresión cambian. Elegir cómo escala ε con N es un cambio de hipótesis y requiere firma del usuario (R3-6).

## 4. Estado de la compuerta

| Rama | Estado |
|---|---|
| E | Estabilidad local con β = b/N (PASA); equilibrio no resuelto (fase sándwich inaccesible a movimientos locales). |
| L | Con ε fijo, **adversa**: el sándwich no se suprime frente a la geometría. Viva solo si se reformula el escalamiento de ε (R3-6). |
