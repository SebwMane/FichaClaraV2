# Ω-C: resultados de L-3c-E′ (rama euclídea con β = b/N)

- **Fecha:** 2026-10-12.
- **Commit de ejecución:** `399b136` (siete tramos, un único commit limpio).
- **Especificación congelada:** `OMEGA_2_PROPUESTA_CAUSAL.md` §12 y enmienda L3cE′-A1.
- **Datos:** `results/l3ce_prime/`.

## 1. Resultado oficial: **PASA (ventana estable con β = b/N)**

Clasificador validado en los tres tamaños (positivos 97–100%, capas 0%). Ningún tamaño tocó el piso de la malla.

| b | N = 32 | N = 64 | N = 128 |
|---|---|---|---|
| 2.00 | 0.95 | 1.00 | 1.00 |
| 2.83 | 0.85 | 1.00 | 1.00 |
| 3.36 | 0.70 | 0.90 | 1.00 |
| 4.00 | **0.30** ← b_c | 0.75 | 0.95 |
| 4.76 | 0.00 | **0.40** ← b_c | 0.90 |
| 5.66 | 0.00 | 0.05 | 0.70 |
| 6.73 | 0.00 | 0.00 | **0.30** ← b_c |
| 8.00 | 0.00 | 0.00 | 0.00 |

(Fracción de muestras clasificadas como geométricas; filas intermedias en `results/l3ce_prime/`.)

b_c = 4.00, 4.76 y 6.73: no cae al duplicar N; sube uno y dos pasos.

## 2. Lectura

**Nivel 1 (N ≤ 128):** con el acoplamiento β = b/N y partiendo de geometría, la fase geométrica se mantiene hasta un b_c que no disminuye con N. La rama euclídea **sigue viva** bajo R3-5.

**Por qué no se puede leer más (límites que pesan):**
1. **Histéresis probable.** Todas las cadenas parten de geometría y usan 1000 barridos de termalización a todo N. Un sistema más grande puede tardar más en abandonar la geometría, así que b_c(128) puede estar sobrestimado más que b_c(32). El aumento de b_c podría ser parcialmente un artefacto; el PASA (b_c no cae) es más robusto que la tendencia creciente.
2. **Potencia declarada.** La prueba no distingue β_c ∝ 1/N de β_c ∝ ln N / N.
3. Una cadena por celda; tres tamaños.

## 3. Siguiente prueba recomendada (por preregistrar)

**L-3c-E″ (histéresis):** repetir la malla partiendo de un 2-orden bicapa y medir b_c «de bajada». Si ambos b_c (desde geometría y desde capas) quedan dentro de un paso a cada N, la transición está bien localizada y el PASA se fortalece; si se separan y la separación crece con N, el PASA se rebaja a «estabilidad metaestable».

## 4. Estado de la compuerta

| Rama | Estado |
|---|---|
| E | **Viva** bajo R3-5, pendiente de la prueba de histéresis. |
| L | Bloqueada por muestreo; Bloque A (vía analítica) preregistrado, sin ejecutar. |

R3-1 sigue pendiente: ninguna rama está cerrada con éxito.
