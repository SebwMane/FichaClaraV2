# Ω-D — Puerta L-ΩD-0: resultados y reorientación

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-d`.
- **Prerregistro:** `docs/OMEGA_D_L0_PRERREGISTRO.md`.
- **Datos:** `results/omega_d_l0/`. 54 grafos: 30 del panel y 24 finales C0 informativos. Duración: 2 min.
- **Reparto:** auditoría de lectura con Haiku; código de Sonnet, revisado antes de correr; análisis y decisión de Opus.

## 1. Veredicto oficial: K1 FALLA → Ω-D (definición D1) muere en su premisa

| Criterio | Resultado | Datos |
|---|---|---|
| **K1** (¿la difusión distingue el pegado incoherente?) | **FALLA** | f_nd: RGG3 = 0.44. Retazos-2³ = 0.36–0.40 y retazos-3³ = 0.31–0.34, es decir, **menos** aristas «no difusivas» que la geometría |
| **K2** (¿ley de escala libre de dimensión?) | **FALLA** | Retículos ✔ (T³ s = +0.01, cuadrado −0.03). Pero RGG3 −0.40, RGG2 −0.70 a −0.78, anillo k12 −0.87. ER y RR −0.03 a −0.06, **casi nulas**: los expansores parecen «difusivos» en τ ≤ 8 |
| **K3** (¿cliques degeneradas?) | **PASA** | Q₈: K10 = 0.0002, caveman = 0.0014, RGG3 = 1.03 |

## 2. Contraste con las predicciones

| # | Predicción | Resultado |
|---|---|---|
| 1 | K1 pasa | ✘ |
| 2 | K2 pasa en G; Q del retículo ≈ d/2 | **Mitad.** La fórmula analítica se confirma: T³ Q₈ = 1.59 frente a 1.5, cuadrado 1.007 frente a 1.0. Pero la ley de escala falla en geometrías aleatorias y en el anillo, y los expansores **no** caen en τ ≤ 8 |
| 3 | Cliques degeneradas | ✔ |
| 4 | Retículo de cliques mixto | ✔ (f_nd = 0.90, Q₈ = 0.11) |
| 5 | Finales C0 entre RGG3 y retazos | Distinto: s = −0.77, f_nd = 0.65; más no difusivo que ambos |

**Por qué falló (nivel 2).**
1. **Transitorio perezoso.** El operador perezoso conserva la mitad de la masa en el nodo, y eso produce un transitorio no gaussiano cuya duración crece con el alcance de los vecinos (k = 12). Las geometrías con aristas de longitud variable (RGG, anillo k12) no alcanzan el régimen difusivo en τ ≤ 8; los retículos de vecinos inmediatos sí.
2. **Las aristas entre bloques no son anómalas a τ corto.** A esas escalas se ven como aristas geométricas largas (Q mayor, Q₈ = 1.3–2.1), no como no difusivas: la mezcla entre bloques ocurre después.

## 3. Observación *post hoc* (no prerregistrada; no rescata la puerta)

La firma predicha aparece, pero **más tarde**, con τ = 16–32:

| Familia | Q₈ | Q₃₂ | Q₃₂/Q₈ |
|---|---|---|---|
| Geometrías (RGG3, RGG2, T³, cuadrado, anillo) | 0.42–1.59 | 0.41–1.24 | **0.76–0.99** |
| Retazos-2³ / 3³ | 1.41 / 1.81 | 0.21 / 0.04 | **0.15 / 0.02** |
| ER, RR | 2.55–2.61 | 0.00 | **0.00** |
| RGG3 + 1 % de atajos | 1.05 | 0.79 | 0.75 (no distinguible) |
| WS β = 0.01 | 0.46 | 0.46 | 0.98 (no distinguible) |
| Finales C0 (N = 729) | 0.67 | 0.40 | 0.59 |
| Cliques, caveman | ≈ 0 | ≈ 0 | degenerados |

**Lectura.** La persistencia de la distancia de difusión a τ largo separa la geometría de los retazos y de los expansores. **No** separa los atajos escasos.

Pero lo que mide es, en esencia, la **mezcla lenta** («mundo grande»). El certificado (F3) y el instrumento P1-D.3 ya cubren esa propiedad.

## 4. Por qué no rescato Ω-D con esa observación (decisión del cerebro)

1. **Metodología.** Cambiar la ventana de τ después de ver los datos es un camino bifurcado. La regla prerregistrada dice «K1 falla → P3».
2. **Argumento de fondo (nivel 2, sin depender de los datos).** Toda acción que premie la persistencia de la difusión premia la mezcla lenta. En un espacio de dimensión d, el tiempo de mezcla escala como N^{2/d}: es **máximo en d = 1**. Una Ω-D de ese tipo empujaría hacia **cadenas**, como hizo C0 desde el inicio U.
   - Para frenar ese sesgo haría falta un término que favorezca la expansión, como la competencia de C0.
   - Entonces la dimensión saldría de un **equilibrio entre dos exponentes ajustado por parámetros**: es la situación que prohíbe el requisito R4 (meseta C1-D).
   - Lo mismo pasa con el diseño «α objetivo» (L-ΩD-T2, confirmado aquí: Q ≈ d/2 en retículos).
3. **Conclusión.** La difusión **contiene** información mesoscópica útil como **observable**, pero no da por sí sola un principio que **seleccione** la dimensión. Ω-D se cierra como mecanismo; su observable queda como diagnóstico secundario, junto a D_s.

## 5. Lo que se rescata

| Hallazgo | Nivel | Uso |
|---|---|---|
| Q ≈ d/2 en retículos (T³ 1.59, cuadrado 1.007) | 1–2 (fórmula gaussiana confirmada) | Cualquier valor objetivo de solapamiento impone la dimensión. Advertencia para diseños futuros |
| Cliques degeneradas en difusión (Q → 0) | Verificado | Definición limpia de «degeneración» para futuros términos anti-clique |
| Persistencia Q₃₂/Q₈ | *Post hoc* | Diagnóstico candidato de mezcla lenta; necesita prerregistro con datos nuevos si se quiere usar |
| El transitorio perezoso enmascara la difusión en τ corto | Verificado | Cualquier uso futuro del núcleo de calor debe usar τ ≫ (alcance local)² |

## 6. Reorientación del programa

### 6.1 El patrón común (conjetura de trabajo, nivel 4)

Todos los mecanismos probados caen en uno de cuatro destinos:

| Destino | Ejemplos |
|---|---|
| **0D** (cliques, vacío) | S0, Ω-B; recompensas monótonas en la difusión (T1) |
| **1D** (cadenas) | C0 desde U; recompensas de persistencia |
| **∞-D** (expansor o mundo pequeño) | C0 desde E; retazos |
| **d ajustada por parámetros** | Objetivo α*; q* de 4-ciclos; c* de C0 |

**Conjetura:** ninguna regla local (energía o crecimiento), sin un entero impuesto, selecciona de forma robusta una dimensión finita mayor que 1.

### 6.2 Precedente externo (nivel 3; referencias por verificar antes de citarlas formalmente)

En triangulaciones dinámicas euclídeas, los ensambles solo producen fases degeneradas: «arrugada», de dimensión efectivamente infinita, y de «polímeros ramificados», de dimensión fractal baja. Una fase extendida de dimensión intermedia apareció **solo** al imponer una estructura causal (triangulaciones dinámicas causales). Es exactamente nuestro patrón de destinos.

### 6.3 Propuesta para el Consejo: P3 con orden de construcción como estructura causal mínima

- El Consejo aplazó la causalidad. Los datos de Ω, más ese precedente, sugieren que **el orden podría ser el ingrediente que selecciona la dimensión**.
- P3 (crecimiento) proporciona ese orden de forma natural: es la pregunta 4 del borrador de P3.
- Se tomaría **solo como orden combinatorio**, sin interpretación física de tiempo, horizontes ni otros objetos, y sin mecánica cuántica.

**Riesgos conocidos de P3 que su prerregistro debe enfrentar:**

| Regla de crecimiento | Resultado esperado | Nivel |
|---|---|---|
| Ancla uniforme con adhesión local | Mundo pequeño (diámetro logarítmico) | 3 |
| Pegado de símplices en la frontera | Geometría hiperbólica emergente (también mundo pequeño), según la literatura de redes con geometría emergente | 3 |
| Pegado restringido a variedad | La dimensión la fija la del símplice: impuesta | — |

La pregunta de diseño concreta es: ¿puede un orden de construcción, por ejemplo por capas o rebanadas, más una regla local sin entero impuesto, evitar los cuatro destinos?

### 6.4 Siguiente paso concreto

1. El Consejo ratifica (o no) que P3 incluya el orden de construcción como estructura causal mínima.
2. Con esa respuesta y las seis preguntas del borrador, se prerregistra **L-P3-0**: análisis de la regla frente a los cuatro destinos, con las referencias de 6.2–6.3 verificadas. Sin dinámica.
3. No se ejecuta nada más hasta entonces.
