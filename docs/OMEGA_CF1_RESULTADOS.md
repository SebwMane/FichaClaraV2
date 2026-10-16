# Ω — CF-1: resultados (crecimiento por compleción flag causal)

- Prerregistro: `docs/OMEGA_CF1_PRERREGISTRO.md` (commit 1e06862), con la adenda de ejecución §8 (bac572e), escrita antes de la corrida completa.
- Código: `tools/cf1_growth.py`; 22/22 tests pasan. Implementación de un agente Sonnet; revisión del cerebro.
- Datos: `results/cf1/{runs.jsonl, summary.json}` (51 ejecuciones más V0). Límite de 1200 s de crecimiento por ejecución.
- **Incidencia:** el contenedor se reinició con 38/51 ejecuciones hechas. El cerebro relanzó la herramienta, que es reanudable por ejecución, con los mismos parámetros, y completó las 51. No se cambió nada del diseño.

## 1. Veredicto formal

**CF-1: NO INTERPRETABLE.** Fallan dos calibraciones congeladas (§5):

| Calibración | Resultado | Causa |
|---|---|---|
| K1 (EL W6-inválida) | ✔ | Las 3 semillas abortan por coste antes de N₂ (cuentan como no W6-válidas, según la adenda §8) |
| **K2** (V4: δ ≈ 1/c) | **✘** | c = 2: δ₁₂ = 0 (tope de radio); c = 3: δ₁₂ = 0.87; c = 4: aborta |
| **K3** (V5 no W6-válida) | **✘** | V5 es W6-válida (grado saturado en 4.0), aunque RC-3 la excluye (X1 + X3) |

Por la regla congelada no se lee ningún resultado positivo. Tampoco lo había: **ninguna ejecución de V1 llegó a N₂**, así que W6 no es evaluable en ninguna (§2).

**Diagnóstico de K2 y K3 (nivel 1, verificado en `runs.jsonl`): las dos calibraciones estaban mal diseñadas por el cerebro.**
- En V4 y V5 **no hubo ni una extensión** (n_extend = 0): la compleción nunca terminó.
- Con la cota (V4) y con solo pares (V5), la estructura es un **tubo que crece sin fin por pura compleción**, igual que el tubo periódico ya demostrado en CF-0.1 para C₂ con w = 3.
- Un tubo es 1D. RC-3 lo excluye bien (X1 + X3) y su grado satura, por eso es W6-válido.
- Supuse que «cota c ⇒ producto de rango c» y que «solo pares ⇒ proliferación». Las dos suposiciones eran falsas. Los controles no medían lo que tenían que medir. Error del diseño, no de los instrumentos.

## 2. Lo que se observó (descriptivo, no decide; nivel 1 solo para estas instancias)

**V0 pasa.** B₇ (128 elementos) y B₈ (256) cierran con los tres planificadores, y el relabelado da estructuras isomorfas para w = 3…6. La coherencia causal de CF-0.1b se extiende a w = 7, 8.

**V1, V1-R y V3 (crecimiento): inflación de rango, la misma para cualquier semilla.**

| n | r_loc (mediana de la subida en el interior) | Subida máxima |
|---|---|---|
| 500 | 3 | 6–8 |
| 2000 | 4 | 8–10 |
| ≈ 8000 | 5 | 9–12 |

- Esta trayectoria es **idéntica para w₀ = 2, 3 y 4**, y con los planificadores LOCAL y RANK.
- Con unos 10⁴ elementos, el grado medio de Hasse es ≈ 9.9. Una rejilla de rango 2 tendría 4.
- Más del 99 % de los elementos vienen de completaciones: en V1 con w₀ = 2, 57 extensiones dieron 10 695 completaciones.
- La fracción de esquinas ambiguas es 0.56–0.65.
- r_loc sube ≈ 1 cada vez que n se multiplica por 4. Es compatible con un **rango ∝ log n**, la firma de una estructura tipo hipercubo (mundo pequeño).
- Todas las ejecuciones abortan por coste entre 7 000 y 11 000 elementos.
- **V3 (extensión libre) se comporta igual que V1.** El criterio relacional «subida < bajada» no frena nada.

**V1-A (asíncrono):** otro degenerado. r_loc = 1 pero la subida máxima llega a 19–24: concentradores.

**V2 (dos semillas):** r_loc = 4–5 tras la fusión. No se distingue de la inflación de V1, que alcanza los mismos valores, así que **CF-T5 (suma de rangos) no es contrastable** con este diseño.

**V7 (semilla aleatoria):** subida máxima de hasta 93. Explosión.

**D1 (olvido de la semilla):** la semilla se olvida (la trayectoria no depende de w₀), pero **hacia la inflación, no hacia un valor finito común**. No es la señal que pide B y no activa nada.

## 3. Lectura del cerebro: el dilema de la identidad de dirección (nivel 2)

Con CF-0, CF-0.1 y CF-1, las reglas de compleción sin etiquetas tienen exactamente tres comportamientos (nivel 1 para las reglas e instancias ensayadas):

| Régimen | Qué hace | Geometría |
|---|---|---|
| Compleción flag desde una semilla, sin extensión | Cierra en B_w y se detiene | Finita: no hay crecimiento |
| Compleción sin cierre (pares, o cota) | Tubo periódico que crece | **1D** (excluido por RC-3) |
| Compleción flag con extensión (ER o EL) | Cada cubierta nueva se completa con **todas** las demás cubiertas de su elemento | **Inflación** de rango ∝ log n (excluido; aborta por coste) |

**Por qué (nivel 2).** La regla flag, que hace falta para la coherencia (CF-0.1), trata cada paso nuevo como una **dirección independiente** de todas las existentes. Para crecer sin inflar, una extensión tendría que identificarse como «continuación de una dirección ya existente». Pero en una estructura sin etiquetas **no hay información local que diga a qué dirección pertenece un paso**: la identidad de dirección es justo la información que aportaría un alfabeto de direcciones etiquetadas, es decir, una presentación de ℤ^w. Eso es S2 FAIL en E6.

**Dilema de la identidad de dirección:**
- Sin etiquetas, la compleción coherente o no crece (cubo), o crece en 1D (tubo), o infla.
- Con etiquetas de dirección, crece como ℤ^w, pero el alfabeto codifica w.

Encaja con todo el programa: el rango de conmutación es la d (ARQ-0, R3, W-0), y la confluencia sin etiquetas no tiene forma de fijarlo durante el crecimiento.

**Alcance.** Es un resultado sobre **esta familia** de reglas (C_flag, C₂, cota; ER, EL; planificadores causales y asíncrono) y estos tamaños (≤ 2.5·10⁴), no un teorema. La firma log n es una extrapolación de nivel 4 más allá de 10⁴.

## 4. Contraste de predicciones del cerebro

| Predicción | Resultado |
|---|---|
| V0: B₇ y B₈, relabelado isomorfo | Acertada |
| V1: CF1-neg por esquinas o inflación; aborta antes de N₂ (adenda) | **Acertada** en la sustancia. El desenlace formal es NO INTERPRETABLE |
| «Más probable A+ para w = 2» (§6) | **Fallida**: w = 2 infla igual |
| V1-A peor que V1 | Acertada (otro degenerado) |
| V2: suma de rangos | No contrastable |
| V3: inflación | Acertada |
| **V4: δ ≈ 1/c** | **Fallida**: no hay producto; es un tubo de pura compleción (error de diseño) |
| **V5: proliferación** | **Fallida**: es un tubo 1D (error de diseño; CF-0.1 ya lo mostraba para C₂ con w = 3 y no lo usé) |
| V7: explosión | Acertada |
| D1: no habrá olvido | **Fallida en la forma**: hay olvido, pero hacia la inflación |

**Error metodológico reconocido (lección 33).** Diseñé los controles V4 y V5 a partir de lo que *creía* que harían, sin comprobarlo con resultados que ya tenía. CF-0.1 ya mostraba que C₂ con w = 3 da un tubo infinito. Un control positivo debe verificarse en pequeño **antes** de congelarlo como calibración.

## 5. Decisión del cerebro y propuesta al Consejo

1. **CF-1 se registra como NO INTERPRETABLE en lo formal y NEGATIVO en lo descriptivo** para el problema A. Ninguna regla de compleción sin etiquetas ensayada produjo una geometría extensa de dimensión finita mayor que 1.
2. **No propongo CF-2 con controles corregidos.**
   - El valor de la información es bajo: la inflación aparece igual con todas las semillas, planificadores y criterios de extensión, y aborta antes de N₂.
   - Unos controles corregidos solo convertirían el desenlace formal en CF1-neg.
   - Si el Consejo quiere el desenlace formal limpio, el coste es una fase con controles verificados en pequeño y la misma regla. Lo dejo a su decisión.
3. **Propuesta:**
   - **(a)** ratificar el dilema de la identidad de dirección (nivel 2) y añadirlo al mapa negativo;
   - **(b)** registrar que el problema A queda **sin candidato** también por la vía de la confluencia;
   - **(c)** cualquier propuesta futura para A que use direcciones etiquetadas debe justificar sus etiquetas por un principio independiente (criterio 1 de reapertura de B, aplicado aquí a A), porque las etiquetas son S2;
   - **(d)** decidir si el programa entra en pausa con A y B documentados, que es mi recomendación.
4. **B no se toca.** Nada de CF-1 lo afecta.
