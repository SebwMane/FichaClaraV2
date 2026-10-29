# R18-1: muestreo multicanónico en la altura. Resultados (2026-10-19)

- **Preregistro:** `OMEGA_R18_WL_ALTURA_PREREGISTRO.md` (commit 290001b), con la aclaración §7 y la previsión declarada antes de correr.
- **Código:** commit 6aface1, árbol limpio, 15 partes, 14 min de pared.
- **Datos:** `results/r18_wl/`. **Figura:** `docs/figuras/r18_wl_altura.png`.
- **Tipo:** piloto técnico. **No es evidencia.**

## 1. Veredicto: **NO CONVERGIDO** con N = 32, 64 y 128 (se cumplió la previsión)

| N | Fase A (Wang-Landau) | R̂ (relaciones / h / S₄) | Viajes por cadena | Concordancia | Validación con ES-C0 | Veredicto |
|---|---|---|---|---|---|---|
| 32 | 1 halving, ln f = 0.5 | 3.1 / 25 / 1.8 | 0, 0, 0, 0 | falla | **falla**: 0.379 frente a 0.348 (z = 10) | no convergido |
| 64 | 0 halvings, ln f = 1 | 8.9 / ∞ / 2.7 | 0, 0, 0, 0 | falla | no aplica | no convergido |
| 128 | 0 halvings, ln f = 1 | 31 / 11 / 2.7 | 0, 0, 0, 0 | falla | no aplica | no convergido |

## 2. Lo que muestra la figura

- **Fase A:** W(h) aprendió diferencias enormes entre alturas vecinas (de 10⁴ con N = 32 a 10⁶ con N = 64 y 128). Wang-Landau sumó ln f = 1 en cada propuesta mientras la cadena se quedaba en una sola altura, y nunca vio el histograma plano.
- **Fase B:** con ese W congelado, cada cadena quedó **clavada** en su altura de partida. Ninguna hizo un solo viaje entre sectores. Al reponderar, todo el peso cae en una altura (6 con N = 32, 3 con N = 64 y 128), así que el promedio reponderado no representa la medida de conteo. Por eso falla la validación con N = 32.
- **El sesgo empeoró la mezcla respecto de ES-C0:** con N = 32 y b = 0, ES-C0 sí mezclaba sin sesgo. Un W mal aprendido crea barreras nuevas en vez de quitarlas.

## 3. Lectura (nivel 2)

- La barrera entre las mesetas de altura 2 y 3 o más es **cinética**. Para bajar de altura hay que romper a la vez todas las cadenas de 3 elementos, y los movimientos de una relación no lo logran en estos tiempos. Ni la prueba corta con una actualización más suave de W (que no se adoptó) volvió a la altura 2 con N = 32.
- **Wang-Landau sobre la altura no basta con este movimiento local.** Aplanar la altura exige poder cambiarla, y el movimiento local casi no la cambia hacia abajo.
- Esto no prueba que el espacio esté desconectado: el movimiento local es irreducible. Lo que muestra es que la mezcla práctica con N de 32 a 128 no se logra ni con Metropolis simple (ES-0, ES-C0) ni con sesgo en la altura (R18-1).

## 4. Decisión que el preregistro devuelve al usuario

- **(ii) Rediseñar el movimiento colectivo:** pasados y futuros generales, con probabilidades de propuesta calculables. Es investigación de algoritmos con riesgo alto.
- **(iii) Cerrar la afirmación de mezcla en este rango:** registrar como resultado de factibilidad que, con los movimientos disponibles, la medida de órdenes con etiquetado natural no se puede muestrear de forma fiable con N ≥ 64. Así H_sel queda indeterminada en Ω-E y en Ω-C.
