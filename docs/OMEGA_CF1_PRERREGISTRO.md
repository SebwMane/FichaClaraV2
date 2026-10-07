# Ω — CF-1: crecimiento y coalescencia por compleción flag causal (prerregistro; solo problema A)

- **Fecha:** 2026-10-07.
- **Rama:** `claude/omega-cf1`. Base congelada: `claude/omega-cf01-congelado`.
- **Estado:** PRERREGISTRO, escrito antes de cualquier código de crecimiento.
- **Mandato del Consejo:** «CF-1 AUTORIZABLE», con una condición: debe ser un experimento del problema A **extremadamente hostil**, no una demostración de que Ω funciona.
- **Pregunta del Consejo, adoptada literalmente:**

> ¿Puede una regla local de completitud, aplicada respetando una precedencia causal local, hacer crecer y fusionar estructuras sin introducir explícitamente una dimensión, una aridad máxima o un planificador global, produciendo una geometría macroscópicamente coherente y estable?

**B sigue cerrado.** Nada de CF-1 lo reabre. Solo un olvido completo de la semilla (D1) llevaría a revisar los criterios de reapertura, y nunca a reabrir B automáticamente.

---

## 0. Acta del Consejo sobre CF-0.1/0.1b y auditoría

### 0.1 Ratificado

| Punto | Registro |
|---|---|
| CF-0.1 válido; W-T4 confirmado para las reglas ensayadas; «ternas ⇒ cualquier rango» refutado; flag funciona bajo orden causal (w = 2–6); paralelismo dentro de un nivel compatible, entre niveles no; planificador global innecesario | Adoptado |
| «Causalidad como ley física»: no demostrado. «Emergencia de dimensión»: no. B: no se reabre | Adoptado |
| Lección 31 reformulada por el Consejo: «El orden causal no se añade como una dimensión geométrica; actúa como una condición de precedencia para que una regla local de completitud pueda converger» | **Sustituye** al texto anterior de la lección 31 |
| Tres capas: causalidad → regla de completitud → estructura relacional; y solo después, estructura → geometría efectiva | Adoptado como esquema de lectura |
| Pruebas exigidas: olvido de la semilla (S₂, S₃, S₄), isomorfismo y relabelado, representación, w = 7, 8, fronteras, varias semillas | Incorporadas (§3) |

### 0.2 Correcciones y respuestas del cerebro

**Ω-1. «¿De dónde viene el orden causal?»** La pregunta es justa, pero tiene una respuesta precisa en este diseño.
- En CF, el orden **es** la estructura que crece: los elementos de Ω y su relación de cubrimiento. No hay un reloj externo ni coordenadas temporales.
- Lo que sí se introduce es la **semilla** (una raíz) y el **planificador**, que lee el cono pasado.
- Por tanto CF-1 no genera causalidad: usa la relación de la propia estructura como precedencia. Se declara así, como pide el Consejo.
- Además, en CF-0.1b el planificador usaba el **rango desde la raíz**, que es una cantidad global. En CF-1 se sustituye por una regla de **cono pasado** (CAUSAL-LOCAL, §2) y se comprueba que las dos son equivalentes.

**Ω-2. Instrumentos.** El Consejo pide medir D_eff, D_s, k̄, σ_k, C, L y G. Por las lecciones del programa (sobre todo la 18: el instrumento no es el mecanismo), **solo deciden los instrumentos validados**: RC-3 con W5 y W6. El resto se reporta sin peso decisorio, con la dimensión espectral marcada como no calibrada en este programa.

**Ω-3. La frontera es la incógnita principal (CF-T3).** Un elemento de esquina al que le falta más de una dirección no sabe cuál completa. El Consejo la llama «el enemigo silencioso». CF-1 la mide de forma explícita (§4, observable de esquinas).

---

## 1. Objetos

- **Estructura:** orden parcial J dado por las relaciones de cubrimiento (`up` y `down`).
- **Geometría medida:** el **grafo de Hasse no dirigido del propio J** (espacio-tiempo), **no** su retículo de cortes. El retículo de cortes es un espacio de configuraciones (lección 29).
- **Semilla S_w:** una raíz con w cubiertas.
- **Semillas hostiles:**
  - **R₂₀:** un orden aleatorio de 2 dimensiones de 20 elementos, que tiene varios minimales.
  - **D₂₊₂:** dos raíces independientes S₂ ⊔ S₂.

## 2. Regla (sin constantes de aridad ni de dimensión)

**C_flag.** Para un elemento z y un subconjunto S de sus cubiertas superiores con |S| ≥ 2, la instancia (z, S) es aplicable si todos los subsupremos J(S ∖ {s}) existen y J(S) no existe. Aplicarla crea un único elemento que cubre esos |S| subsupremos. No hay cota de aridad. Es la regla de `tools/cf01_completion.py` con k = flag.

**Planificadores:**
- **CAUSAL-LOCAL (principal).** Una instancia en z es *elegible* si no hay ninguna instancia pendiente en un elemento estrictamente anterior a z (cono pasado). Entre las elegibles se elige al azar de forma uniforme.
- **CAUSAL-RANK (control de equivalencia).** El de CF-0.1b: menor rango de z primero, después mayor aridad, después identificador.

**Extensión relacional ER.** Cuando no queda ninguna instancia pendiente, se elige al azar uniforme un elemento p con **subida(p) < bajada(p)** y se le añade una cubierta superior nueva. Sin constantes.

**Variantes del diseño:**

| Id | Cambio | Para qué |
|---|---|---|
| ASYNC | En cada paso, con probabilidad 1/2 se extiende aunque haya completaciones pendientes; si no, se completa | Robustez frente al orden (E6-D D3). El valor 1/2 es una variante de robustez, no un parámetro del mecanismo |
| EL | Extensión en cualquier elemento | Control negativo de inflación (CF-T3) |
| CAP(c) | Ninguna operación puede llevar la subida de un elemento por encima de c, con c ∈ {2, 3, 4}, semilla S₆ | Control codificado de tipo Ω3: d debería seguir a c (calibración de E6-D D2) |
| PAIR | C₂ en lugar de C_flag, semilla S₃ | Control: proliferación |

**RNG:** `rng_from_key((MASTER_CF1, id_variante, w, semilla, k))` con `MASTER_CF1 = 20261025`. k = 0 para el planificador, 1 para la extensión, 2 para la semilla aleatoria y 3 para ASYNC.

## 3. Diseño

| Bloque | Contenido | Semillas |
|---|---|---|
| **V0, cierre** | C_flag causal desde S₇ y S₈. ¿Da B₇ y B₈? Prueba de relabelado: permutar los identificadores de los elementos de S_w (w = 3…6) y comprobar el isomorfismo del resultado | — |
| **V1, principal** | ER + CAUSAL-LOCAL, semillas S₂, S₃, S₄ | 3 por w |
| V1-R | Igual que V1, con CAUSAL-RANK | 3 por w |
| V1-A | Igual que V1, con ASYNC | 3 por w |
| **V2, varias semillas** | D₂₊₂ y S₂ ⊔ S₃. Se crece cada componente hasta N₁/2; después **un** elemento nuevo cubre a un maximal de cada componente (elegidos al azar); se sigue con V1 | 3 |
| V3 | EL desde S₂ | 3 |
| V4 | CAP(c) desde S₆ | 3 por c |
| V5 | PAIR desde S₃ | 3 |
| **V7, semilla hostil** | ER + CAUSAL-LOCAL desde R₂₀ | 3 |

**Tamaños:** instantáneas de la **misma estructura creciente** con N₁ = 10⁴, N₂ = 8·10⁴ y N₃ = 6.4·10⁵ elementos. Es el dominio de W6.
- Si una variante explota en coste, se registra el tamaño alcanzado y se declara «no evaluable a N₃». No se cambia el diseño.
- La ejecución empieza con un humo (N ≤ 2·10³) y una **estimación de coste**, que se reporta antes de la corrida completa.

## 4. Medidas

**Deciden:**
- RC-3 (X1–X4) en los pares (N₁, N₂) y (N₂, N₃).
- W5 en ambos pares.
- **W6:** grado medio de la componente gigante saturado (|k₃ − k₂| / k₂ ≤ 0.05).

**Se reportan sin decidir:**
- δ₁₂ y δ₂₃;
- r_loc (mediana de la subida en el interior, es decir, elementos con subida = bajada);
- el producto r_loc · δ;
- elementos por nivel de rango;
- σ_k;
- **fracción de esquinas ambiguas** (elementos de frontera a los que les falta más de una dirección según la comparación de subida y bajada);
- recuento de supremos ambiguos (como en CF-0.1);
- dimensión espectral (no calibrada).

**Tiempo de ejecución** por bloque.

## 5. Criterios (congelados)

**A-positivo** para (variante, w):
- en al menos 2 de 3 semillas la estructura es **W6-válida**;
- y no está excluida por RC-3 en ninguno de los dos pares.

**Robustez frente al orden** (requisito para leer V1 como resultado):
- V1, V1-R y V1-A tienen la misma clase W6 y RC-3;
- y |δ̄₁₂(V1) − δ̄₁₂(V1-x)| ≤ 0.04 para x = R, A, con medias sobre las semillas.

**Trazabilidad de d:** se declara «**inicio**» si en V1 δ̄₁₂ está en 1/w ± 0.06 para cada w ∈ {2, 3, 4}.

**Calibración (deben cumplirse para que CF-1 sea interpretable):**
- **K1.** V3 (EL) es W6-inválida en al menos 2 de 3 semillas: W6 detecta la inflación en esta arquitectura.
- **K2.** V4: δ̄₁₂ sigue a c, con |δ̄₁₂ − 1/c| ≤ 0.08 para cada c. E6-D detecta el DIAL de la cota.
- **K3.** V5 (PAIR) no es W6-válida o no alcanza N₁ por proliferación.

Si K1, K2 o K3 fallan, **CF-1 es NO INTERPRETABLE**: se documenta y no se lee ningún resultado positivo.

**Desenlaces:**

| Desenlace | Condición | Lectura |
|---|---|---|
| **CF1-A+** | V1 A-positivo para al menos un w, robusto frente al orden, K1–K3 cumplidas | Primera regla de Ω sin constantes que produce geometría no degenerada desde una semilla. d trazado al **inicio** |
| **CF1-A±** | A-positivo solo con un planificador | Coherencia dependiente del planificador. No cuenta como A+ |
| **CF1-neg** | V1 no es A-positivo | Resultado negativo, con el mecanismo de fallo documentado: frontera o esquinas, inflación, proliferación o anisotropía |
| **D1 (olvido)** | En V1, la d trazada difiere de w₀ para algún w₀ y converge a un valor común | **No reabre B.** Activa E6-D completo y el análisis contra los criterios de reapertura |

**V2 y V7 no deciden A.** Responden a CF-T5: ¿suma de rangos (r_loc en el futuro común ≈ w₁ + w₂), fragmentación o explosión?

## 6. Predicciones del cerebro (congeladas; no deciden)

| Bloque | Predicción |
|---|---|
| V0 | B₇ y B₈ cierran; el relabelado da estructuras isomorfas |
| V1 | **Incierta.** Más probable para w = 2: A-positivo con δ ≈ 0.5. Para w = 3, 4: riesgo alto de defectos de esquina, con grado creciente y W6 inválida. Si las esquinas no dañan, A+ con d = w₀ (**inicio**) |
| V1-R | Igual que V1 |
| V1-A | Puede romper la coherencia: si se extiende antes de cerrar, la extensión deja de ser causal. Predicción: peor que V1 para w ≥ 3 |
| V2 | r_loc ≈ 4 en el futuro común (suma de rangos): la mezcla no es W6-válida |
| V3 | Inflación: grado creciente, W6 inválida (K1) |
| V4 | δ ≈ 1/c (K2) |
| V5 | No alcanza N₁ o no es W6-válida (K3) |
| V7 | Varios minimales actúan como varias semillas: suma de rangos y explosión. Degenerado |
| Desenlace global | Lo más probable es CF1-A+ solo para w = 2 o CF1-neg por esquinas. **D1: no habrá olvido** |

## 7. Ejecución

- Agente Sonnet: `tools/cf1_growth.py` (reutiliza C_flag de `tools/cf01_completion.py` y la medida de `tools/rc3.py`), tests y `results/cf1/`.
- Antes de la corrida completa: humo y estimación de coste, reportados al cerebro.
- Revisión del cerebro y resultados en `docs/OMEGA_CF1_RESULTADOS.md`.
