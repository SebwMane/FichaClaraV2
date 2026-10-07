# Ω — CF-0: confluencia como mecanismo de coalescencia (análisis sin dinámica y propuesta formal al Consejo)

- **Fecha:** 2026-10-07.
- **Rama:** `claude/omega-cf0`. Base congelada: `claude/omega-cierre-b-congelado`.
- **Origen:** pregunta del usuario («¿qué opinas de confluencia/conmutatividad como mecanismo de coalescencia?») y encargo de plantearla formalmente, auditar los resultados y añadir lo que merezca explorarse.
- **Estado:** B cerrado y ninguna dinámica nueva autorizada. Este documento **no ejecuta nada**: analiza y propone. Todo lo que pide ejecución (CF-1) queda pendiente del Consejo.

---

## 1. Auditoría de lo que el programa ya sabe sobre la conmutatividad

| Registro | Lo que dice | Lectura en esta auditoría |
|---|---|---|
| L-DIM-0 (`OMEGA_LDIM_0.md`) | La conmutatividad local («dos pasos en direcciones distintas se cierran en un cuadrado») es lo que mata la expansión exponencial. No basta con tener ciclos: el cactus y el 2-árbol tienen muchos | Confirmado y vigente (nivel 3 + 2) |
| ARQ-0 (`OMEGA_ARQ_0.md` §3, locus ii) | Candidato principal al «?»: el rango de relaciones conmutativas. Pregunta abierta: «¿puede un sistema de relaciones locales generar cuadrados coherentes cuyo rango no esté prescrito por el alfabeto?». Propuso **L-ARQ-1** con la pregunta (q1) sobre una medida intrínseca del rango | **Hallazgo A1: el programa ha vuelto al punto que ARQ-0 ya señalaba, y L-ARQ-1 nunca se ejecutó.** R1 y R3 lo rodearon. La confluencia es la versión *dinámica* de ese locus |
| C1 (Fase 2) | Los 4-ciclos saturados convertían la dimensión en parámetro, porque se prescribía el **número** de cuadrados | Prescribir cuadrados es codificar. Completarlos sin cuota es otra cosa (§2) |
| R1-1 / SQ | La regla de **estabilidad** (borrar aristas sin 4-ciclo) nucleó dominios con cuadrados, pero no coalesció; Z³ quedó fijo | **Hallazgo A2: SQ no era una regla de confluencia.** Conservaba cuadrados existentes, pero no completaba cuñas abiertas. La coalescencia que faltó en R1 es exactamente la que aporta la compleción |
| E6-S v1.1 | P0-LAT(2), «sin triángulos y c₄ = deg(deg − 2)/2», que es la firma local de conmutación completa, sale **CIEGO** (descansa en Z^d para todo d) | La firma de conmutatividad no codifica d. Coincide con la lección 23 |
| R3-T1 | La coalescencia sin defectos es la unión de down-sets: retículo distributivo, con d = anchura asintótica | Correcto en lo matemático. Pero véase el hallazgo A3 |
| T-Mat (R3-0 §2) | Los cortes consistentes de una computación distribuida forman un retículo distributivo (Mattern) | **Hallazgo A3: la geometría que midió R3 es un espacio de configuraciones, no un espacio físico** (§3.1) |
| RC-3 | Heisenberg: predije X4 y salió INTERMEDIO, no excluido | Es el único objeto del panel donde el crecimiento (4) supera al número de direcciones locales (2). Véase §3.4 |
| Cierre B | d se reduce al número de direcciones independientes que conmutan | La confluencia **no** puede reabrir B. Su papel posible está en el problema A (§4) |

---

## 2. Definiciones

Se trabaja sobre un **orden causal J** que crece; no hay grafos con coordenadas. La geometría candidata es el **grafo de Hasse del propio J** (espacio-tiempo), no el de su retículo de cortes. La razón está en A3.

- **Compleción de diamante, CD (aridad 2).** Si un elemento z tiene dos cubiertas superiores distintas x e y que no tienen ninguna cubierta superior común, se añade un elemento nuevo u que cubre exactamente a x y a y. Se aplica hasta el cierre.
- **Extensión relacional, ER.** Se elige uniformemente un elemento p con **grado de subida < grado de bajada** y se le añade una cubierta superior nueva; después se aplica CD hasta el cierre. Es un criterio relacional, sin constantes, del mismo tipo que P0-LAT(2).
- **Extensión libre, EL.** Como ER, pero p se elige entre todos los elementos.
- **Semilla de anchura w₀.** Una raíz con w₀ cubiertas superiores.

---

## 3. Resultados del análisis

### 3.1 CF-T1: unir historias independientes **suma** dimensiones (nivel 1)

Para órdenes disjuntos, L(J₁ ⊔ J₂) ≅ L(J₁) × L(J₂). Es inmediato: un down-set de la unión es un par de down-sets.

**Consecuencia para R3 (corrección interpretativa, no matemática).**
- En la geometría de cortes, cada proceso independiente añade una dimensión, no un lugar.
- Es el espacio de **configuraciones** de un sistema distribuido (Mattern): su dimensión cuenta procesos.
- «d = w» en R3 era, en parte, la afirmación trivial de que la dimensión del espacio de estados es el número de procesos.
- Esto explica por qué B se redujo a un recuento con tanta limpieza. **Se registra como lección 29.**

### 3.2 CF-T2: los pares no bastan; hace falta la condición del cubo (nivel 1 + 3)

**Comprobación exacta.** `tools/cf0_cd_check.py` → `results/cf0/cd_check.json`, sin azar. Se aplica CD de aridad 2 con elementos **nuevos** desde una raíz de anchura w:

| w | CD fresca | Identificación ingenua | Retículo booleano |
|---|---|---|---|
| 2 | 4 elementos, cierra (el diamante) | 4, cierra | 4 |
| 3 | **no cierra** (> 2000 elementos) | 5: colapsa a un solo techo | 8 |
| 4 | **no cierra** (> 2000; 1743 maximales) | 6: colapsa | 16 |

**Lectura (nivel 1 para estas reglas).** Con w ≥ 3, la conmutación por pares completada localmente **no genera el producto**.
- Si cada cuña se cierra con un elemento nuevo, cada par crea su propio supremo, y esos supremos forman cuñas nuevas sin fin. Es el defecto de esquina del §3.3, que aparece ya en la propia semilla.
- Si se reutilizan elementos de forma ingenua, la estructura colapsa.
- Para obtener el producto hace falta una regla de **coherencia de orden superior**: cuando tres cuadrados comparten una esquina por pares, se rellena un único cubo.
  - En los complejos cúbicos CAT(0) es la condición de enlace de Gromov (enlaces de bandera).
  - En reescritura es la condición del cubo, o coherencia de ternas críticas.
  - Nivel 3, **pendiente de verificación bibliográfica** antes de cualquier prerregistro.
- Con esa regla, la compleción genera el orden producto de w₀ cadenas. Su Hasse es una caja de ℤ^{w₀}, con δ ≈ 1/w₀ como J1 en R3-0b.
- Con w₀ = 2, el producto de dos cadenas es el **diamante causal de Minkowski 1+1** en coordenadas de cono de luz (nivel 3, estándar).
  - La función de rango hace de tiempo, y las anticadenas de rango fijo son las rebanadas espaciales.
  - **La confluencia coherente produce directamente estructura causal con tiempo.** Ninguna dinámica anterior de Ω lo había hecho.

**CF-T6, escalera de aridad (nivel 2; advertencia de numerología).** *[Enmendado: véase §7. El texto original se conserva como historia.]*
- La coherencia por pares basta para rango ≤ 2; el rango ≥ 3 necesita la condición de ternas.
- Si las ternas bastan para todo rango superior (como sugieren la condición del cubo y la teoría de ternas críticas, que es nivel 3 sin verificar), la escalera de aridad **se detiene en 3**.
- **Esto no selecciona d = 3.** Con la regla de ternas son accesibles todos los rangos ≥ 1, y d sigue siendo w₀.
- Lo que sí muestra es que **la aridad de la ley de coherencia acota el rango coherente alcanzable**: con solo pares, rango ≤ 2.
- Es otra instancia de W-T4: la aridad de la ley elegida determina qué dimensiones son posibles. Se registra con advertencia explícita, porque un 3 que aparece en la aridad es justo el tipo de resultado que el programa prohíbe leer como selección.

### 3.3 CF-T3: el dilema de la extensión (nivel 2)

En un producto de w cadenas, un elemento del interior tiene w cubiertas superiores. Extenderlo le añade una cubierta w + 1, y CD la combina con las demás: **nace una dirección nueva** que se propaga a todo su cono futuro.

Por tanto, con información local:

| Regla de extensión | Resultado previsto | Clase |
|---|---|---|
| EL (en cualquier sitio) | Cada extensión interior suma una dirección: **inflación de rango** y explosión | Degenerado; W6.3 debería detectarlo (grado creciente) |
| Con cota de grado k | El rango satura en ⌊k/2⌋ | **Codificación Ω3** (k = 2d); E6 lo detecta como P0-DEG, que es DIAL |
| ER (subida < bajada) | Solo se extiende la frontera, sin constantes. El rango queda en w₀ | **CONSERVADOR**: d es un fósil de la semilla |
| Saber «extender solo si hay menos de w cubiertas» | Conserva el rango, pero nombra w | S1 FAIL |

**ER es, hasta donde llega este análisis, la primera regla de Ω sin constantes para la que se predice geometría no degenerada.** Su d procede del inicio, de acuerdo con el cierre de B.

**Riesgo no resuelto (nivel 4).** En una esquina de la frontera falta más de una dirección. La cubierta nueva no sabe cuál de ellas es, y CD puede identificarla de forma inconsistente en esquinas distintas. Eso daría **defectos de grano dentro de una sola semilla**. El lema de Newman (confluencia local + terminación ⇒ confluencia global) no aplica, porque el crecimiento no termina. **Es la incógnita principal que CF-1 resolvería.**

### 3.4 CF-T4: la conmutación parcial no da un término medio, salvo la nilpotencia (nivel 3)

- **Grupos de Artin rectangulares.** Si solo conmutan algunos pares de generadores:
  - si el grafo de conmutación es completo, el grupo es ℤ^n;
  - si no, contiene un grupo libre de rango 2 y crece exponencialmente.
- **Nilpotencia** (conmutar salvo un elemento central):
  - el grado de crecimiento es Σ k·rango(γ_k/γ_{k+1}) (Bass–Guivarc'h);
  - Heisenberg: 2 direcciones locales y grado 4.
  - Es el único caso conocido en que **la dimensión de crecimiento supera al número de direcciones locales**.
  - Sigue fijado por la presentación elegida (ley), así que no reabre B.
- **Respuesta parcial a (q1) de ARQ-0** (nivel 2): el rango local r_loc (grado de subida en el interior) y el grado de crecimiento 1/δ coinciden en los productos abelianos, y difieren por los pesos de la serie central en los nilpotentes.
  - El cociente r_loc · δ ≈ 1 sería una **prueba de coherencia de rango**.
  - Si r_loc · δ < 1, señala nilpotencia o defectos. Se propone como observable exploratorio, no decisorio.

### 3.5 CF-T5: la coalescencia de varias semillas suma rangos (nivel 2)

- Sean dos semillas sin pasado común, de anchuras w₁ y w₂.
- Si un elemento nuevo cubre a uno de cada una, CD combina las direcciones de ambas en su cono futuro, de manera que el futuro común es localmente un producto de rango w₁ + w₂ (versión dinámica de CF-T1).
- **Consecuencia:**
  - con nucleación múltiple, el rango crece con el número de semillas que se unen, es decir, de forma extensiva en N: explosión;
  - una d finita estable exige **un único origen causal común** (monogénesis) y ninguna ramificación interior posterior.
- **Lectura (nivel 4, solo como hipótesis):** en un universo confluente, la dimensión sería un fósil de la ramificación inicial.

---

## 4. Auditoría S0–S2 de E6 (escrita, antes de cualquier código)

| Elemento | Valor | Riesgo |
|---|---|---|
| Constantes enteras | Aridad de la coherencia: pares (CD) y ternas (cubo) | **Constantes discretas (W-T4, CF-T6).** «Solo pares» limita el rango a ≤ 2. Se barren pares frente a pares + ternas |
| Criterio de extensión | Relacional (subida < bajada), sin constante | Ninguno, por la misma lógica que P0-LAT(2), que salió CIEGO |
| Cantidades conservadas | Ninguna | — |
| Condición inicial | Anchura de la semilla w₀ | **Fuente declarada de d** (trazabilidad: inicio) |
| S2 (pares inversos que conmutan, productos de cadenas) | CD **genera** productos de cadenas; no los presupone | Si la semilla ya fuera un producto, sería S2 FAIL. La semilla es solo una raíz con w₀ cubiertas |
| E6-S | No aplicable: los órdenes no están en el dominio de v1.1 | Solo S0–S2 y E6-D. La prueba D1 (olvido) es la decisiva y **se predice CONSERVADOR** |

---

## 5. Propuesta formal al Consejo

### 5.1 Qué se pide

1. **Ratificar CF-0** como análisis: hallazgos A1–A3 y CF-T1 a CF-T6 con sus niveles.
2. **Registrar la corrección interpretativa de R3 (A3).** La geometría de cortes es un espacio de configuraciones, y su dimensión cuenta procesos.
3. **Verificar la bibliografía** de la condición del cubo y de la coherencia de ternas (Haiku), antes de prerregistrar.
4. **Autorizar CF-1**, la primera dinámica desde el cierre de B, **solo para el problema A**. Prerregistro pendiente; diseño en §5.2.

### 5.2 Diseño mínimo de CF-1 (a prerregistrar si se autoriza)

| Variante | Qué prueba | Predicción del cerebro |
|---|---|---|
| V1: ER + condición del cubo, una semilla, w₀ = 2, 3, 4 | ¿Coalescencia sin defectos desde un origen? | **Incierta.** O bien W6-válida con δ ≈ 1/w₀ y r_loc · δ ≈ 1 (CONSERVADOR), o bien defectos de esquina con grado creciente (W6-inválida) |
| V2: ER, dos semillas unidas por un evento | ¿Suma de rangos (CF-T5)? | Rango w₁ + w₂ en el futuro común |
| V3: EL | ¿Inflación (CF-T3)? | Grado creciente, W6-inválida, X1 a tamaño grande |
| V4: ER con cota de grado k | ¿Codificación Ω3? | d = ⌊k/2⌋, DIAL en k |
| V5: ER solo con pares (sin cubo) | Escalera de aridad (CF-T6) | w₀ = 2: producto; w₀ ≥ 3: proliferación sin cierre, degenerado |

**Medidas:**
- RC-3 + W5 + W6, con tres tamaños de la misma estructura creciente (en un crecimiento es natural).
- r_loc · δ como observable exploratorio.
- E6-D D1 (olvido de w₀).
- Trazabilidad de d.

**Calibración dinámica de E6-D (exigida por su diseño):** V4 es el control codificado conocido.

### 5.3 Qué puede y qué no puede salir de CF-1

- **Puede dar:** la primera estructura no degenerada y W6-válida producida por una regla de Ω sin constantes. Sería progreso en A, con d trazado al inicio.
- **Puede dar:** la confirmación de nivel 1 de CF-T3 y CF-T5. «Una d finita exige un origen común» pasaría de argumento a hecho, para estas reglas.
- **Puede dar:** la localización del obstáculo real de la coalescencia, si V1 falla por defectos de esquina.
- **No puede dar:** progreso en B. Se predice CONSERVADOR, en línea con el cierre. Si saliera algo distinto (d ≠ w₀ estable y que olvida el inicio), no reabriría B automáticamente: se analizaría primero contra los criterios de reapertura.

### 5.4 Coste estimado

- Análogo a R3-0b: enumeración y crecimiento de órdenes de 10⁴–6.4·10⁵ elementos, más RC-3.
- Del orden de 1–2 h de cálculo y una sesión de implementación.

---

## 6. Lo que merece explorarse (hallazgos que añado a la propuesta)

1. **Configuración frente a espacio (A3, CF-T1).** La distinción entre la geometría de cortes (dimensión = número de procesos) y la geometría del propio orden causal (espacio-tiempo) no estaba en el programa. Cambia la lectura de R3 y sugiere que cualquier propuesta futura declare cuál de las dos mide.
2. **Monogénesis (CF-T5).** Si se confirma, «d finita ⇒ origen causal común» es un enunciado estructural nuevo y comprobable. No dice nada sobre qué d, pero sí sobre qué tipo de historia permite una d finita.
3. **La escalera de aridad (CF-T6).** La coherencia local por pares no sostiene rango ≥ 3: hace falta una regla de ternas. Es un hecho estructural verificable sobre qué leyes permiten qué dimensiones, y hay que leerlo con la advertencia de numerología del §3.2.
4. **Coherencia de rango r_loc · δ (CF-T4).** Responde en parte a la pregunta (q1) que ARQ-0 dejó abierta. Además separa productos abelianos, estructuras nilpotentes y estructuras con defectos con un solo número intrínseco.
5. **Heisenberg como caso límite.** Es la única vía conocida por la que la dimensión de crecimiento supera al recuento de direcciones locales. No reabre B (fija la presentación), pero es el ejemplo más cercano a «dimensión que no se ve localmente», y merece una nota en el mapa de exclusiones.

---

## 7. Acta del Consejo y enmienda

**Veredicto del Consejo:** CF-0 RATIFICADO CON ENMIENDA; CF-1 autorizable tras la enmienda y el prerregistro, pero **no ejecutar todavía**. Se ordena una microfase CF-0.1 de formalización, sin dinámica nueva.

### 7.1 Ratificado

| Punto | Registro |
|---|---|
| Corrección de R3: espacio de configuraciones ≠ espacio físico | Adoptada (lección 29) |
| SQ no era confluente, con la precisión del Consejo: «SQ no era una prueba contra la confluencia; era una prueba contra una regla concreta que preserva cuadrados sin implementar una condición de confluencia» | Adoptada literalmente |
| La comprobación CD como contraejemplo constructivo | Adoptada |
| «Estructura causal con tiempo» rebajado a «una compleción cúbica coherente puede generar una estructura causal de tipo producto en los casos construidos; admite una interpretación causal, no es tiempo físico» | **Sustituye** a la frase del §3.2 |
| Monogénesis queda en nivel 4 hasta CF-1 | Adoptado |
| r_loc · δ es un observable, no un certificado (Heisenberg ≈ 1/2) | Adoptado |
| B sigue congelado. Solo habría señal si inicios w₀ = 1, 2, 3, 4… convergieran a una misma d, olvidando el inicio, sin codificación en la regla | Adoptado (coincide con E6-D, prueba D1) |
| Red Team: desempates, orden de actualización, definición de interior, simultaneidad | Adoptado como obligación del prerregistro de CF-1. Una sola secuencia de actualización no basta |

### 7.2 CF-T6 enmendado (texto del Consejo, adoptado)

> **CF-T6, escalera de coherencia de orden superior.** La compleción por pares genera correctamente el producto de dos cadenas, pero no basta para construir productos de tres o más cadenas mediante la regla CD ingenua. Para w ≥ 3 se requiere coherencia de orden superior. La literatura sobre complejos cúbicos relaciona esta coherencia con la condición de enlace *flag* de Gromov, y en concurrencia y reescritura aparecen condiciones de cubo relacionadas. Sin embargo, no está demostrado todavía que una condición ternaria particular sea suficiente para todos los w. La suficiencia de ternas, o la necesidad de aridades superiores, queda como cuestión matemática de CF-0.1.

### 7.3 Auditoría del acta: el Consejo tiene razón, y el problema es más fuerte de lo que dice

**Argumento del cerebro (nivel 2).** La condición *flag* exige rellenar cliques de **todo** tamaño. Un enlace con los cuatro triángulos de un tetraedro rellenos y el interior vacío satisface «toda terna se rellena», pero no es *flag*. En términos generativos, con w = 4:
- La regla de ternas, aplicada en cada esquina x_i, crea un techo de 3-cubo para cada esquina.
- Sin una identificación de aridad 4, los techos de esquinas distintas son elementos **distintos**.
- Sus pares comparten cubiertas inferiores sin cubierta común, y vuelve la proliferación.

**Predicción:** la coherencia de aridad acotada k cierra **exactamente hasta w = k**. Es el tercer escenario del Consejo: «la ley necesita una aridad que crece con w». Por W-T4, eso convierte a k en un **dial del rango máximo**. La única regla sin constante sería la compleción *flag* completa (aridad no acotada).

**Riesgo adicional (Red Team).** La comprobación previa ya mostró que, si los pares se completan *antes* que los cubos, la proliferación aparece aunque exista la regla de cubos. **El orden de aplicación puede decidir el resultado.** CF-0.1 lo mide.

---

## 8. CF-0.1: formalización de la completitud (prerregistro; cálculo exacto, sin azar y sin dinámica de crecimiento)

**Reglas.** Todas actúan desde una semilla: una raíz con w cubiertas. Cada elemento nuevo se representa por su conjunto de cubiertas inferiores.

- **C₂:** si z tiene dos cubiertas x, y sin cubierta superior común, se crea un elemento nuevo que cubre x e y.
- **C_k (k ≥ 3):** C₂, más la regla de aridad m para 3 ≤ m ≤ k. Si z tiene m cubiertas S cuyos m sub-supremos de tamaño m − 1 existen y no tienen cubierta común, se crea un único elemento que los cubre a todos (la regla «7 de 8 → el octavo», generalizada).
- **C_flag:** C_k sin cota de aridad.

**Órdenes de aplicación** (exigencia del Red Team):
- **ALTA:** en cada ronda se aplica primero la mayor aridad disponible.
- **BAJA:** primero los pares.
- **RONDA:** se calculan todas las aplicaciones posibles sobre el estado actual y se aplican juntas.
- Dentro de cada orden, los empates se resuelven por identificador y, como prueba de robustez, por 3 permutaciones de identificadores generadas con `rng_from_key((MASTER_CF, w, k, perm))`, con `MASTER_CF = 20261024`.

**Barrido:** w = 2…6, k = 2…5 y flag, con los 3 órdenes. Límite de 5000 elementos.

**Medidas:**
- si el proceso termina y con cuántos elementos;
- número de maximales;
- **isomorfismo con el retículo booleano B_w** (2^w elementos, un techo, función de rango con los coeficientes binomiales);
- para los que no terminan, el recuento de elementos por nivel de rango, y un **certificado de periodicidad**: si dos niveles consecutivos tienen la misma configuración local salvo isomorfismo, la no terminación queda demostrada por inducción (nivel 1).

**Predicciones del cerebro (congeladas):**

| Regla | Predicción |
|---|---|
| C₂ | Cierra en B₂ con w = 2; no termina con w ≥ 3 (certificado de periodicidad con w = 3) |
| C_k con ALTA y RONDA | Cierra en B_w si y solo si w ≤ k |
| C_flag con ALTA y RONDA | Cierra en B_w para todo w = 2…6 |
| Con BAJA | Falla ya con w = 3 para todo k, incluido flag: **dependencia del orden** |
| Permutaciones de desempate | No cambian el veredicto dentro de un mismo orden |

**Criterio de lectura (congelado):**

| Resultado | Lectura |
|---|---|
| **CF01-a** | Ninguna C_k acotada cierra para todo w y C_flag sí. La coherencia exige aridad no acotada: la única ley sin constante es flag (k ↔ rango máximo, W-T4). CF-1 debe usar C_flag |
| **CF01-b** | C₃ cierra para todo w = 2…6. Las ternas bastan: el Consejo tenía la duda y mi predicción falla |
| **CF01-c** | Ni siquiera C_flag cierra para algún w. La compleción local generativa no basta: CF-1 no procede en esta forma |
| **Orden** | Si el veredicto depende del orden, CF-1 debe incluir el orden como variable de E6-D (D3) y el resultado se lee como «coherencia condicionada al orden» |
