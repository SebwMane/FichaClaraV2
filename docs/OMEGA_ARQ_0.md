# Ω — L-ARQ-0: síntesis de arquitectura (sesión conceptual, sin dinámica nueva)

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-arq`.
- **Base congelada:** `claude/omega-ldim1-congelado`.
- **Mandato del Consejo:**
  - L-DIM-1 cerrada (SELECCIÓN-NO-VIABLE), D5 refutado como mecanismo estático y L-DIM-2 no autorizada;
  - **no construir otro funcional**;
  - hacer primero una síntesis de arquitectura: «¿en qué nivel tendría que aparecer la dimensionalidad para que no sea un parámetro disfrazado?».
- **Contenido:**
  - ratificación con correcciones;
  - un resultado nuevo de nivel 1 (el lema del dial), verificado sobre datos ya cerrados;
  - la jerarquía con la evidencia acumulada;
  - un análisis de los candidatos para «?»;
  - criterios de entrada para cualquier propuesta futura.

## 0. Ratificación y correcciones al dictamen

| Punto del dictamen | Contraste | Decisión |
|---|---|---|
| L-DIM-1 cerrada; D5-B refutado como mecanismo estático; L-DIM-2 no autorizada | Coincide con la regla de muerte prerregistrada | Ratificado |
| d* = f(p) asintóticamente, f(θ, N, p) a tamaño finito | Correcto. §1 lo convierte en un teorema general (lema del dial) | Ratificado y ampliado |
| Nulo B: «L1–L3 solo preguntan por propiedades internas, así que permutar las etiquetas no debería afectarlas» | **Impreciso.** L1–L3 **sí** dependen de las etiquetas: d* es el argmin sobre los **grupos que definen las etiquetas**, y L2 compara estratos. El nulo B mide con qué frecuencia una partición **aleatoria** de las mismas instancias en grupos «d» produce un mínimo interior invariante: 8–9 %. Es una medida de **potencia** del criterio | Se adopta el **principio de no identificación**, en esta lectura: una propiedad estable que sobrevive a una asignación aleatoria de la etiqueta dimensional no es evidencia de selección. Operativo: toda prueba de selección futura exige un nulo de etiquetas permutadas ≤ 5 % |
| NC-5: predicción fallida, dirección compatible, magnitud no confirmada | Correcto; ya estaba así en los resultados | Sin cambios |
| Pólya / Nash-Williams asimétrico; la transitoriedad no es ley general de Ω | Correcto (corrección ya registrada) | Sin cambios |
| Separar el problema A (exclusión de degenerados) del B (selección de d) | Correcto, y §1 muestra que B no puede resolverse con pesos escalares | Ratificado (§4) |
| «No concluir que Ω no puede producir dimensiones» | Correcto. Lo falsado es la vía variacional escalar sobre bolas y cáscaras | Ratificado |
| Frase oficial de cierre de D5 | — | Adoptada literalmente en el mapa negativo |

## 1. Lema del dial (L-ARQ-T1, nivel 1)

**Enunciado.**
- Sea un conjunto finito de candidatos (dimensiones) con valores (A_d, B_d) y un funcional F_θ(d) = A_d + θ B_d, con θ > 0.
- (a) Si θ₁ < θ₂, d₁ minimiza F_{θ₁} y d₂ minimiza F_{θ₂}, entonces B_{d₂} ≤ B_{d₁}.
- (b) Los únicos d que son argmin para algún θ son los **vértices de la envolvente convexa inferior** de los puntos {(B_d, A_d)}. Cada vértice se selecciona en un **intervalo** de θ, cuyos extremos son las pendientes de las aristas adyacentes.

**Demostración.**
- Se suman F_{θ₁}(d₁) ≤ F_{θ₁}(d₂) y F_{θ₂}(d₂) ≤ F_{θ₂}(d₁).
- Resulta (θ₂ − θ₁)(B_{d₂} − B_{d₁}) ≤ 0, que es (a).
- (b) es la dualidad estándar entre argmin de funcionales lineales y vértices de la envolvente inferior.
- No se usa ninguna monotonía de A ni de B en d.

**Generalización** (comparativa estática monótona, nivel 3, Topkis 1978, véase §5):
- Si F(d, θ) tiene diferencias crecientes en (d, θ), el argmin es monótono en θ.
- Con varios pesos, F = A + Σ θ_i B_i, los seleccionables son los vértices de la envolvente inferior en ℝ^{k+1}, y el «dial» pasa a ser una partición del espacio de parámetros.

**Consecuencias:**
1. **Toda competencia escalar entre presiones es un dial.** El entero seleccionado lo fija el peso. «Un único θ para todas las d» no lo impide: solo elige un punto del dial.
2. **Lo que puede seleccionarse es una propiedad de los valores (A_d, B_d), no de la estructura**, y esos valores dependen de N y de la discretización. Por eso la escalera se desplaza con N (L7) y las d sin vértice no se seleccionan nunca.
3. Una selección **no de dial** requiere una de estas dos cosas:
   - **(α)** que el peso esté fijado por un principio independiente, documentado antes de mirar resultados;
   - **(β)** un mecanismo **no escalarizable**, es decir, una restricción cuyo conjunto factible contenga una sola clase dimensional, sin un balance entre costes.
4. T1 de L-DIM-1 es el caso asintótico. Allí, (B_d, A_d) ∝ (r^{−p(d−1)}, d/r), y la envolvente la fija p.

**Verificación sobre datos cerrados** (`tools/arq_dial_check.py` → `results/arq/dial_check.json`). Es la comprobación de un teorema, no un experimento nuevo.

| p | Vértices (seleccionables) | Ventanas de θ | Escalera predicha = observada (25/25 θ) |
|---|---|---|---|
| 1/2 | 1, 2, 3 | 1: [0, 0.32] · 2: [0.32, 40] · 3: [40, ∞) | ✔ |
| 1 | 1, 2, 3 | 1: [0, 0.47] · 2: [0.47, 254] · 3: [254, ∞) | ✔ |
| 2 | 1, 2, 3 | 1: [0, 1.09] · 2: [1.09, 1.9·10⁴] · 3: [1.9·10⁴, ∞) | ✔ |

- d = 4 y d = 5 **no son vértices**. A N ≈ 2·10⁴, su ψ medio no es menor que el de d = 3, porque r* es corto. Por eso **nunca** podían seleccionarse, sea cual sea θ.
- La «meseta de d = 2» es solo la arista más larga de la envolvente: 2.1, 2.7 y 4.2 décadas.
- **La meseta (L3) mide la geometría de la envolvente, no una selección.** Es la explicación de nivel 1 del principio de no identificación.

## 2. La jerarquía con la evidencia acumulada

| Nivel | Qué fija Ω | Instrumento | Excluye | No implica el nivel siguiente (contraejemplo) |
|---|---|---|---|---|
| Conectividad | — | Componente gigante | Fragmentos | ER, RR (sin localidad) |
| Localidad | C0 la produce (253/675) | `classify_c0`, H_null | Expansores, ER | C0: complejos de cliques locales sin geometría |
| Coherencia (Følner, homogeneización) | Ninguna dinámica de Ω la produjo | ρ (Nivel I), γ (B-bis) | Retazos, mundo pequeño, cruces lentos | Árboles críticos y cadenas (coherentes) |
| Estructura topológica (un extremo a toda escala) | Ninguna dinámica la produjo | B-ter (anillos), K | Árboles, cactus, cliques, expansores | 2-árbol (falla el instrumento); Z², Z³, Z⁴… comparten la clase |
| **?** | — | — | — | — |
| Dimensión | Medible (D_B2, Nivel II) | P1-D.3 | — | D ≈ 3 fabricable (WS β = 0.001) |
| Geometría tipo variedad | — | Certificado Ω-1.1 | — | — |

**Lectura.**
- Ω tiene **instrumentos** en todos los niveles, pero **mecanismos** solo en los dos primeros.
- El hueco «?» no es solo «qué selecciona d». Está **detrás** de dos niveles que ninguna dinámica de Ω ha alcanzado todavía: coherencia y un extremo. Abrir B antes que A es saltarse niveles otra vez, lo mismo que L-DIM-2 respecto de L-DIM-1.

## 3. Candidatos para «?» (dónde podría vivir la dimensión sin ser un parámetro)

Se evalúan con los criterios E1–E4 del §4.

| Locus | Idea | Estado | Por qué |
|---|---|---|---|
| **(i) Variacional sobre geometría ya formada** | Un balance de costes elige d | **Falsado** (L-DIM-1 + lema del dial) | Todo balance escalar es un dial (§1) |
| **(ii) Algebraico: número de direcciones independientes que conmutan** | En espacios homogéneos de crecimiento polinómico, d es un **conteo**: el grado de crecimiento Σ k·rango(G_k/G_{k+1}) (Bass–Guivarc'h); en Z^d, el número de generadores independientes que conmutan | **Abierto; candidato principal a «?»** | Un conteo es entero por naturaleza y **no puede moverse de forma continua**: no es un dial (tipo β). Pero desplaza la pregunta a «¿por qué ese número de direcciones?». Si el alfabeto de movimientos de la regla local lo prescribe, es **codificación por alfabeto** |
| **(iii) Estructura adicional (orden o causalidad)** | Una relación extra (por ejemplo, un orden parcial) excluye fases degeneradas | **Paralelo externo**, pendiente de verificar (§5) | En triangulaciones dinámicas euclídeas aparecen exactamente las dos degeneraciones de Ω: una fase **arrugada** (≈ clique, dimensión de Hausdorff muy grande) y otra de **polímero ramificado** (≈ árbol, d_H = 2). La versión causal obtiene extensión **imponiendo** una foliación, y además el bloque (el símplice) fija la dimensión: es codificación por alfabeto, como (ii). En conjuntos causales dominan entrópicamente los órdenes de Kleitman–Rothschild (no variedad), análogos al fondo entrópico de la Fase 2. Se toma solo el contenido combinatorio-estadístico; no se adopta mecánica cuántica (prohibición vigente) |
| **(iv) Definida por un proceso** (dimensión espectral) | d como exponente de retorno de un paseo | **Es una medida, no un selector** | Depende de la escala por construcción; ya se usa como observador (D_s) |

**Síntesis (nivel 2).**
- La literatura externa y el mapa negativo de Ω coinciden: los ensembles relacionales sin estructura adicional caen en **compactificación** (clique / arrugado) o en **ramificación** (árbol / polímero ramificado).
- Cuando se obtiene extensión, la dimensión entra por el **alfabeto** (el bloque o los generadores) o por una **estructura impuesta** (foliación).
- El candidato que no es un dial, (ii), tiene la forma «la dimensión es el **rango** de un sistema de relaciones conmutativas». La pregunta científica abierta se puede escribir así:

> **¿Puede un sistema de relaciones locales generar relaciones conmutativas (cuadrados coherentes) cuyo rango no esté prescrito por el alfabeto de la regla?**

Esto enlaza tres resultados previos de Ω:
- L-DIM-0: la conmutatividad local es lo que mata el crecimiento exponencial.
- C1: los 4-ciclos saturados hacen de la dimensión un parámetro, porque allí el **número** de cuadrados se prescribía.
- CH-T6: homogeneidad + coherencia ⇒ dimensión entera.

## 4. Criterios de entrada para cualquier propuesta futura de dimensionalidad

| Id | Criterio | Operativo |
|---|---|---|
| **E1** | Existencia | Produce una clase dimensional entera medible (Nivel II convergente) |
| **E2** | **Prueba del dial** (obligatoria) | Se barre cada parámetro continuo del mecanismo durante ≥ 4 décadas. Si aparecen ≥ 2 enteros distintos, es un **dial** y se rechaza, salvo que el parámetro esté fijado por un principio independiente documentado **antes** de mirar resultados (α) |
| **E3** | Exclusión de degenerados sin d (problema A) | Clique, árboles, 2-árbol, expansores y retazos excluidos por mecanismo, con verificación por B-ter / γ / K |
| **E4** | Escala | El mismo resultado con N, 2N, 4N y 8N |
| **E5** | No identificación | Nulo de etiquetas permutadas ≤ 5 % con potencia suficiente (≥ 5 instancias por d y estrato) |
| **E6** | Alfabeto | Si la regla tiene un alfabeto finito de tipos de movimiento o de bloque, d no puede coincidir con una función de su tamaño, salvo con una justificación independiente (codificación por alfabeto) |

**Orden lógico.** A (E3) antes que B (E1, E2, E5, E6). Ninguna propuesta de selección dimensional se ejecuta mientras no exista una dinámica que alcance los niveles «coherencia» y «un extremo».

## 5. Referencias externas

Verificación bibliográfica de un agente Haiku, revisada por el cerebro.
- El agente **no aportó URL**, así que la verificación es bibliográfica débil.
- Los enunciados los confirma el cerebro con la literatura estándar.
- Todas se usan como **nivel 3**, y solo por su contenido combinatorio y estadístico.

| Referencia | Uso en este documento | Estado |
|---|---|---|
| Topkis, *Oper. Res.* **26**, 305–321 (1978) | Generalización del lema del dial (diferencias crecientes ⇒ argmin monótono) | ✔ |
| Bass, *Proc. London Math. Soc.* (3) **25**, 603–614 (1972); Guivarc'h, *Bull. SMF* **101**, 333–379 (1973) | Grado de crecimiento = Σ k·rango(G_k/G_{k+1}); en Z^d vale d | ✔ |
| Ambjørn–Jurkiewicz, *Phys. Lett. B* **278**, 42 (1992) | Triangulaciones dinámicas euclídeas 4D: fase arrugada y fase de polímero ramificado | ✔ |
| Agishtein–Migdal, *Mod. Phys. Lett. A* **7**, 1039 (1992); véase también *Nucl. Phys. B* **385** (1992) | Ídem. Ellos interpretaron inicialmente la transición como continua | Débil: el agente dudó de la revista. Confianza moderada |
| Bialas–Burda–Krzywicki–Petersson, *Nucl. Phys. B* **472**, 293 (1996) | La transición es de primer orden: no hay límite continuo entre ambas fases | ✔ |
| Ambjørn–Jurkiewicz–Loll, *PRL* **93**, 131301 (2004) y *PRL* **95**, 171301 (2005) | Con foliación causal impuesta aparece un universo extendido; dimensión espectral ≈ 2 a escala corta y ≈ 4 a escala larga | ✔ |
| Ambjørn–Durhuus–Jonsson, *Quantum Geometry*, CUP (1997) | Polímeros ramificados genéricos: d_H = 2 | ✔ |
| Kleitman–Rothschild, *Trans. AMS* **205**, 205–220 (1975) | Casi todos los órdenes parciales tienen tres capas (dominio entrópico no variedad) | ✔ |
| Myrheim (CERN TH-2538, 1978); Meyer (tesis MIT, 1988) | Estimador de dimensión por fracción de orden | Fuentes primarias no consultadas; citado de forma consistente |

**Lectura de nivel 3 para Ω.**
- Un programa externo independiente, con 30 años de simulaciones, encontró **las mismas dos degeneraciones** que el mapa negativo de Ω: la compactificación y la ramificación.
- Ese programa solo las evitó **añadiendo estructura** (la foliación causal), con la dimensión fijada por el bloque.
- No es una prueba de que Ω deba hacer lo mismo. Sí es evidencia externa de que **ningún ensemble relacional sin estructura adicional conocido** ha producido geometría extendida de forma espontánea, y es coherente con la conjetura de nivel 4 del §1 del mapa negativo.

## 6. Decisión del cerebro

1. **No hay ni dinámica ni funcional nuevos en esta fase.** Lo prohíbe el lema del dial: cualquier funcional escalar fallaría E2 por construcción.
2. **Se cierra la línea L-DIM** con el lema del dial como resultado principal. Junto con T1 y T2, completa el cierre de D5-B.
3. **Propuesta al Consejo (sin ejecutar):** abrir **L-ARQ-1**, una sesión conceptual sobre el locus (ii): conmutatividad emergente con rango no prescrito. Antes de cualquier código respondería a tres preguntas:
   - (q1) ¿Existe una medida intrínseca del rango conmutativo, sin coordenadas, distinta del grado de crecimiento? Por ejemplo, contando cuadrados independientes por vértice frente a la razón de crecimiento.
   - (q2) ¿Qué familias de reglas de reescritura local producen cuadrados **sin** prescribir su número (contraste con C1)?
   - (q3) ¿Cómo se evita la codificación por alfabeto (E6)?
4. La alternativa equivalente en el orden lógico es volver primero al **problema A**: una dinámica que alcance coherencia y un extremo. Ninguna de las dinámicas de Ω lo ha logrado. Es la decisión que el Consejo debe tomar entre (ii)-B y A.
