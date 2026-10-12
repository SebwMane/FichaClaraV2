# Ω-C: resultados de L-3c-E″ (histéresis de la rama euclídea)

- **Fecha:** 2026-10-12. **Commit de ejecución:** `3240c7a` (siete tramos, un commit). **Especificación:** `OMEGA_2_PROPUESTA_CAUSAL.md` §14. **Datos:** `results/l3ce_hyst/`.

## 1. Resultado oficial: **INDETERMINADO** (piso en N = 128)

| N | b_c↑ (desde geometría, §12) | b_c↓ (desde bicapa) | brecha (pasos) |
|---|---|---|---|
| 32 | 4.00 | 4.00 | 0 |
| 64 | 4.76 | 2.83 | 3 |
| 128 | 6.73 | ≤ 2.00 (piso) | ≥ 7 |

## 2. Lo que se observa (nivel 1, descriptivo, sin cambiar el veredicto)

- **La histéresis existe y crece rápido con N.** En N = 32 ambos inicios coinciden; en N = 64 se separan 3 pasos; en N = 128, al menos 7.
- **Desde la bicapa las cadenas se congelan.** En N = 64 con b ≥ 2.83 y en N = 128 con b ≥ 2.38 la aceptación es ≈ 0 y la acción se queda junto a la de la bicapa (−1984 y −8064). La dinámica local no logra salir de ese estado.
- Por tanto ninguno de los dos b_c mide la transición de equilibrio: b_c↑ es un límite de **estabilidad de la geometría** y b_c↓ un límite de **estabilidad de las capas**. La transición real está en algún punto entre ambos.

## 3. Consecuencia para la compuerta

- El PASA de L-3c-E′ **se mantiene con su significado literal**: la geometría es localmente estable hasta un b_c que no cae con N. **No** demuestra que la geometría sea la fase de equilibrio en esa ventana.
- La pregunta de equilibrio queda **abierta por método**: Metropolis local no basta en una transición con metaestabilidad tan fuerte.

## 4. Orientación (nivel 4, no es resultado)

Un conteo grosero (entropía de las permutaciones v, ≈ N ln N − N, frente a βS de la bicapa ≈ −bN/2) ubica la transición de equilibrio cerca de b* ≈ 2(ln N − 1): 4.9, 6.3 y 7.7 para N = 32, 64, 128. Está cerca de b_c↑ y crece como ln N, lo que sugeriría que la ventana geométrica sería de equilibrio y no solo metaestable. Es un argumento de orden de magnitud que ignora todas las familias intermedias.

## 5. Siguiente prueba posible (no preregistrada)

Un método de equilibrio para localizar b*: templado paralelo (cadenas a varios b que intercambian estados) o integración termodinámica de la diferencia de energía libre entre ambas fases.
