# Ω — L-A-R1-1 (dinámica SQ): resultados

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-r1`.
- **Prerregistro:** `docs/OMEGA_R1_1_PRERREGISTRO.md`.
- **Datos:** `results/r1_1/{runs.jsonl, pairs.jsonl, summary.json}`. Son 1212 corridas y 606 pares: 270 de la etapa A y 336 de la etapa B, que se activó por regla. ≈ 4 h de CPU; 18/18 pruebas.
- **Código:**
  - `omega/dynamics/sq.py` y `tools/r1_1.py`, escritos por un agente Sonnet y revisados por el cerebro;
  - las órdenes SYNC/ASYNC/RAND, los disparadores B-s, B-t y B-st y el reenganche R2/R3 se ajustan a §1;
  - la clasificación es X0 seguido de RC-3 congelado;
  - s(u,v) se calcula por A³ y se contrasta con el cálculo por conjuntos en las pruebas.
- **Nota de trazabilidad:** el commit WIP `e7b3902` lo hizo el cerebro por el stop hook, a mitad de la corrida. El estado final añade una guarda de lista vacía en `_dominant`, que no cambia ningún resultado.

## 1. Veredicto: **R1-1 NEGATIVA**. Nulo Ω6 = **0/6** pares (B, R) positivos

| (B, R) | Clase dominante (etapa A, 36 pares amorfos) | Resultado G6 | Etapa B |
|---|---|---|---|
| B-s / R2 | FRAGMENTADO 0.75 (resto: NO-EXCLUIDO 9, todos de I1) | MULTIESTABLE | MULTIESTABLE (FRAG. 0.69–0.75) |
| B-s / R3 | FRAGMENTADO 0.69 (NO-EXCLUIDO 9 de I1; X2 2) | MULTIESTABLE | Ídem |
| B-st / R2 | FRAGMENTADO 0.75 (NO-EXCLUIDO 9 de I1) | MULTIESTABLE | Ídem |
| B-st / R3 | FRAGMENTADO 0.75 (NO-EXCLUIDO 9 de I1) | MULTIESTABLE | Ídem |
| B-t / R2 | FRAGMENTADO 0.94 | Atractor FRAGMENTADO | — |
| B-t / R3 | FRAGMENTADO 0.94 | Atractor FRAGMENTADO | — |

- Ningún par (B, R) tiene un atractor NO-EXCLUIDO, ni en la etapa A ni en la B.
- **Control I5 (Z³):**
  - fijo (absorbido en el barrido 0, aristas finales = iniciales) bajo B-s y B-st, en todas las órdenes, semillas y tamaños;
  - destruido bajo B-t (FRAGMENTADO; X2 en 2 casos);
  - es la sanidad del código prevista.

## 2. El residuo NO-EXCLUIDO de I1: inspección (G7) y clase identificada

Todos los NO-EXCLUIDO amorfos vienen de **I1** (inicio denso, grado ≈ 20) bajo B-s y B-st. Inspección del cerebro sobre los datos, con B-s/R2 y semilla 0:

| | N₀ = 5000 | 8N₀ = 40 000 |
|---|---|---|
| E/N final | 7.8–7.9 | 1.95–2.1 (aún decreciendo en el barrido 100) |
| R | 2.92 | 7.4–7.5 |
| Anillos / coherencia | SIN_VENTANA / SIN_VENTANA (r_w = 2) | CONEXO / **INCOHERENTE (γ = 0.11, f_low = 1.0)** |
| Censo | s ≈ 1.0, t ≈ 0.05 (aleatorio: t ≈ k²/N) | s ≈ 0.05 |

**Clase:** **grafo aleatorio con densidad dependiente de N**, un degenerado conocido de tipo expansor. No es un residuo nuevo.

**Mecanismo (nivel 2).**
- En G(N, k/N) el número esperado de 4-ciclos por arista es ≈ k³/N:
  - 1.6 a N₀, así que ~80 % de las aristas sobreviven a B-s;
  - 0.2 a 8N₀, así que el borrado es masivo.
- El «punto casi fijo» de I1 es un **efecto de tamaño finito**: un grafo denso aleatorio tiene cuadrados por azar mientras k³ ≳ N.
- Para N → ∞ con densidad inicial fija, colapsa. 8N₀ ya está en otro régimen.

### W5 — punto ciego nuevo de RC-3 (no previsto)

- RC-3 supone que el par (N, 8N) pertenece al **mismo régimen** del proceso.
- Aquí, el grado medio baja de 15.7 a 3.9 entre los dos tamaños. Eso produce δ = 0.45, que aparenta crecimiento «2D», solo por el **cambio de densidad**, no por escalado.
- X4 no puede activarse porque en N la coherencia se abstiene (SIN_VENTANA), aunque en 8N el grafo es fuertemente incoherente (γ = 0.11).
- **No se cambia RC-3.** Se registra como punto ciego W5 en el mapa de exclusiones. Se propone al Consejo, sin aplicarla a este resultado, una regla de uso futura: «el par solo es válido si E/N difiere < 25 % entre los dos tamaños; si no, se clasifica como CAMBIO-DE-RÉGIMEN».
- El veredicto de R1-1 **no depende** de W5: el criterio G6 ya da MULTIESTABLE, y la inspección G7 identifica la clase.

## 3. Dinámica observada (trayectorias típicas, I2, N₀)

| Orden | Comportamiento |
|---|---|
| SYNC | Rotación sin absorción: E/N ≈ 1.0 y una componente gigante que oscila (0.28–0.65). Es el ciclo «borrar todo → reenganchar al azar» previsto |
| ASYNC | Bosque fragmentado (gigante ≈ 0.005) con motivos cerrados pequeños |
| RAND | **Absorción** (barrido 22–65) con la fracción s → 0.7–1.0, pero con componentes diminutas (gigante ≈ 0.005–0.01) |

- **Lectura:** la nucleación de motivos ricos en cuadrados **sí ocurre** en las órdenes secuenciales (RAND/ASYNC). Lo que falla es la **coalescencia**: los núcleos absorben y se congelan aislados.
- **Corrección a la predicción** del cerebro, que atribuía el fallo a la nucleación (∝ 1/N): el cuello de botella es **crecer y unir dominios**, no formarlos.

## 4. Contraste con las predicciones (§9)

| Predicción | Resultado |
|---|---|
| 1. B-s / B-st desde I1–I4: rotación; FRAGMENTADO o X1 | ✔ en I2–I4 (FRAGMENTADO, con absorción en RAND). ✘ en I1: NO-EXCLUIDO por W5 (la clase es grafo aleatorio, previsto como X1) |
| 2. B-t destruye I5 y fragmenta | ✔ |
| 3. I5 fijo bajo B-s / B-st | ✔ |
| 4. Ω6 = 0/6 | ✔ |
| 5. R1-1 NEGATIVA | ✔. **Mecanismo de fallo corregido:** coalescencia, no nucleación |

## 5. Dictamen del cerebro

1. **R1-1 NEGATIVA:** la familia SQ (18 variantes) no produce una fase estable no degenerada.
   - Destinos: fragmentación (bosques o núcleos aislados), rotación sin absorción y, desde inicios densos, un grafo aleatorio cuya supervivencia es un efecto de tamaño finito.
2. **Alcance (correcciones del Consejo y de R1-0):**
   - SQ era el **primer candidato falsable**, no «la regla mínima correcta».
   - Lo falsado es **SQ y sus variantes**, junto con la idea de que un disparador de estabilidad local ciego a d (presencia de cuadrados) baste por sí solo.
   - «Ninguna regla de estabilidad local funciona» es una **extrapolación de nivel 4**.
   - Por la rama prerregistrada en R1-0 §6, R1 queda **cerrada para reglas de estabilidad local de este tipo**.
3. **Lo que este resultado aporta al problema A (nivel 2):**
   - Con una regla de estabilidad ciega a d, los dominios ricos en cuadrados **se forman, pero no crecen ni se unen**.
   - El ingrediente ausente es un mecanismo de **coalescencia o crecimiento coherente** de dominios.
   - Las clases conocidas que crecen (P3, crecimiento local) terminan en mundo pequeño, 1D, hiperbólico o polímero ramificado.
   - El candidato que queda en el mapa es **R3** (orden o estructura adicional). Por ejemplo, un orden de agregación que haga crecer un dominio por su frontera, con la auditoría E6 obligatoria (que el orden no fije d).
4. **Siguiente paso propuesto al Consejo:**
   - **(a)** ratificar el cierre de R1 con este alcance;
   - **(b)** decidir sobre W5 (regla de uso CAMBIO-DE-RÉGIMEN para el futuro);
   - **(c)** autorizar **R3-0**, un análisis sin dinámica, como R1-0: qué estructuras de orden pueden hacer crecer y unir dominios sin codificar d. El punto de partida es el diagnóstico de coalescencia.
