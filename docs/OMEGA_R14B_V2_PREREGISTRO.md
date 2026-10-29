# R14-B v2: respuesta al grano grueso de órdenes dados (H_L frente a H_E). Preregistro exploratorio (2026-10-19)

- **Mandato:** R15-1 firmada (redactar v2 con las enmiendas 1 a 4 de la Rev. 15, conservar v1 como histórica y volver a decidir después si se ejecuta). **Estado:** PREPARADO, **no ejecutado**. Ejecutarlo requiere una nueva orden explícita del usuario.
- **Reemplaza a** `OMEGA_R14B_PREREGISTRO.md` (v1, histórica, sin cambios).
- **Sin voto.** No reabre F6.

## 0. Qué cambió respecto de v1 y por qué

| Enmienda (Rev. 15) | Cambio en v2 |
|---|---|
| 1. Acotar la afirmación | El título y el alcance pasan a «respuesta de órdenes dados a un mapa de grano grueso». **Esta prueba no dice nada sobre si la dinámica de Ω produce esos órdenes** (eso es H_sel, prueba ES, Rev. 16). |
| 2. Referencia y equivalencia | Cada geometría curva tiene su propia referencia de sprinkling (no la plana). El margen δ deja de ser «0.5 desviaciones estándar a ojo»: se calibra (§3). Se define qué es converger en una ventana discreta de q. |
| 3. Jerarquía de datos | Las cuatro escalas q de una misma realización se tratan como **medidas emparejadas**. La incertidumbre se da por realización y por familia. Las estructuras ya vistas no sirven para confirmar. |

## 1. Pregunta

Dado un orden C de una familia ya vista, ¿cómo cambian sus observables bajo la dilución aleatoria Φ_q (conservar cada elemento con probabilidad q)? En particular:
- ¿entran en las bandas de una geometría al hacer el grano más grueso (comportamiento compatible con H_E)?
- ¿o ya están en banda a la escala de los elementos (H_L)?

## 2. Mapa, escalas y qué se mantiene fijo

- **Φ_q:** dilución aleatoria independiente; orden inducido. Lleva un sprinkling de densidad ρ a otro de densidad qρ.
- **Escalas:** q ∈ {1, 1/2, 1/4, 1/8}, desde N₀ = 8000 (ventana de N efectivo entre 8000 y 1000).
- **Fijo entre escalas:** la región (la misma estructura completa) y los observables, todos **intensivos o normalizados por su referencia del mismo N efectivo**.
- **Convergencia en la ventana discreta (definición):** un observable «converge a banda» si está dentro de la banda en q = 1/4 y en q = 1/8, y su distancia a la banda no crece de q = 1/2 a q = 1/8. «En banda a escala microscópica» significa dentro en q = 1.

## 3. Referencias y bandas de equivalencia

- **Referencia por geometría objetivo:** 60 sprinklings independientes de la **misma** geometría (plana, de Sitter, de Sitter fuerte, anti de Sitter, FRW polvo, D2, Kasner C0 y D0) en cada N efectivo. La referencia plana no se usa para geometrías curvas.
- **Bandas conjuntas:** para el vector de observables O se usa la distancia de Mahalanobis respecto de la referencia, con su matriz de covarianza estimada en la referencia, lo que tiene en cuenta la dependencia entre observables. Banda = percentil 95 de esa distancia en la referencia, **calibrado por validación cruzada**: 40 sprinklings fijan la banda y 20 distintos miden su cobertura real, que se informa.
- **Margen de equivalencia:** se fija desde la cobertura (la banda que contiene el 95% de los 20 de validación); no hay δ elegido a ojo.

## 4. Observables (intensivos)

1. Fracción de pares relacionados y dimensión de Myrheim-Meyer.
2. r por sub-intervalo con m ∈ {10, 20, 40} (regla m ≤ 0.1·m_max; referencia exacta 1/10 en 4D plano y la referencia de su geometría en las curvas).
3. Perfil de intervalos normalizado (N_k/N_pares, k = 0..8).
4. Fracción de pares con igual pasado y futuro (detector de gemelos).

## 5. Estructuras (todas vistas; exploratorias)

- **Geometrías:** plana, de Sitter, de Sitter fuerte, anti de Sitter, FRW polvo, D2, C0 y D0.
- **Defectos:** C5, D4 y D5.
- **Conos:** D6, D7, D8, C4 y octaédrico.
- **Producto 4.**
- **Realizaciones:** 10 por estructura, cada una diluida a los cuatro q con generadores independientes.

## 6. Informe (sin clasificar)

- Por estructura: la trayectoria emparejada de la distancia de Mahalanobis frente a q, con su intervalo por bootstrap de realizaciones.
- La etiqueta descriptiva (en banda microscópica, converge a banda, no converge, indeterminada), **sin voto**.
- **Predicciones congeladas** (las de v1 §1): las geometrías quedan en banda en todos los q; D5 converge a banda (gemelos que sobreviven con probabilidad q²); los conos no convergen.

## 7. Convenciones

- **EXPERIMENT_ID 323; entropía 20261088** (reasignadas de v1, que nunca corrió).
- Costo: unas 3 a 4 horas, incluidas las referencias por geometría.
- Herramienta: agente de código con handoff y auditoría; árbol limpio.

## 8. Orden de ejecución y aclaración A1 (2026-10-19, antes de mirar datos)

- **Orden del usuario:** «Ejecuta R14-B v2 y luego R14-C v2».
- **Generadores:** los mismos de F6 y F7, citados en `tools/r14b.py` (tabla de estructuras). Las estructuras que estaban selladas en F6 (C0, C4, C5, D0, D2 y D4 a D8) se regeneran con la entropía de este experimento. Son vistas y sus semillas no coinciden con las de F6.
- **Vector de Mahalanobis (p = 14):**
  - fracción de relaciones;
  - r₁₀, r₂₀ y r₄₀ (media por sub-intervalo);
  - N_k/N_pares para k = 0 a 8;
  - √(pares gemelos).

  d_MM se informa, pero queda fuera del vector: es función monótona de la fracción de relaciones y volvería singular la covarianza. Cada componente se estandariza con las 40 referencias de ajuste, y la covarianza se encoge con Ledoit-Wolf porque p = 14 está cerca de n = 40. Para los gemelos se usa un piso de desviación de 0.5.
- **Banda:** radio = máx(percentil 95 de las 40 de ajuste, 19.ª de las 20 de validación), que es la lectura literal de §3 («la banda que contiene el 95% de los 20»). Se informan las dos coberturas y el margen. Con datos gaussianos sintéticos, esta regla cubre cerca del 90% de los puntos nuevos (19/21); se declara así. La distancia se informa como D/radio, de modo que estar en banda significa D/radio ≤ 1.
- **Sub-intervalos:** el instrumento de F6 sin cambios. Se toman hasta 300 intervalos con exactamente m elementos interiores, con m ≤ 0.1·m_max y al menos 30 intervalos. Si un m no se puede medir en un N efectivo (sobre todo m = 40 con N ≈ 1000), esa componente se quita y la banda se recalcula sin ella. Se registra cuáles se quitaron.
- **Referencias por q:** N ~ Binomial(8000, q), que es exactamente la distribución de Φ_q aplicada a un sprinkling de 8000 elementos.
- **Objetivo principal:**
  - cada geometría se compara con su propia referencia;
  - defectos, conos y producto se comparan con la referencia plana;
  - además, toda estructura se compara con las 8 referencias.
- **Etiquetas:** por q se toma la mediana de D/radio sobre las 10 realizaciones, con bootstrap emparejado de 1000 réplicas.
  - **Estado en un q:** «dentro» si el extremo superior del IC 95% es ≤ 1; «fuera» si el inferior es > 1; si no, ambiguo.
  - **«Crece»:** el IC de mediana(1/8) − mediana(1/2) queda entero por encima de 0.
  - Las etiquetas de §2 y §6 se aplican sobre esos estados.
