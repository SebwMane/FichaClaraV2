# Ω — L-A-R1-0: espacio de reglas, auditoría de codificación y análisis previo (sin dinámica)

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-r1`.
- **Base congelada:** `claude/omega-rc3-congelado`.
- **Mandato:** sesión extraordinaria del Consejo tras RC-3. Se adopta R1 (reescritura local, fuera del equilibrio) como **última clase abierta que merece una prueba**, en tres etapas: R1-0 (análisis de reglas, sin dinámica), R1-1 (dinámica mínima) y R1-2 (prueba de fase).
- **Este documento es R1-0.**
- **Única ejecución:** un censo estático exacto (`tools/r1_wedge_census.py` → `results/r1_0/wedge_census.json`) que verifica las afirmaciones del §3. No hay dinámica.

## 0. Evaluación del acta del Consejo

### 0.1 Ratificado (con motivo)

| Punto | Motivo |
|---|---|
| RC-3 es un instrumento A− parcial: ni certificado ni RC-4 inmediato | Coincide con el veredicto y con la regla «no tocar tras ver resultados» |
| Generador ≠ instrumento: R1 no puede optimizar contra RC-3 | Evita el aprendizaje del juez (circularidad). Se añade como regla G1 (§5) |
| Generación → clasificación → medición de d, nunca al revés | Es la lección de L-DIM-1 |
| Estabilidad de fase (N, semillas, orden de actualización, condiciones iniciales, olvido) en lugar de una foto | Se adopta íntegra como protocolo de R1-2 |
| «Residuo no clasificado», en lugar de «geometría» | Coincide con O2 de RC-3 |
| Clasificación de las conservaciones C0–C4 | Se adopta, con las correcciones Ω2–Ω3 |
| Orden de actualización ≠ causalidad física; robustez síncrona / asíncrona / aleatoria | Se adopta como prueba D de R1-2 |

### 0.2 Omisiones y correcciones (justificadas)

**Ω1. Una reescritura reversible es un ensemble de equilibrio, es decir, una clase ya cerrada** (nivel 1).
- Una cadena de Markov reversible respecto de π tiene π como estacionaria.
- Con propuestas simétricas y aceptación dependiente de diferencias locales, π ∝ e^{−H}, con H local.
- Eso es exactamente la clase «energía o entropía local de equilibrio», cerrada por C0-T1, la Fase 2 y el paralelo DT.
- **Por tanto, R1 solo tiene contenido nuevo si rompe el balance detallado:**
  - **(a)** dinámica **irreversible con estados absorbentes**, en la que los atractores son puntos fijos;
  - **(b)** o un estado estacionario **forzado** que no sea de Gibbs.
- El acta habla de «reescritura + conservación» sin esta condición. Sin ella, R1 sería una repetición de S0/C0 con otro nombre.

**Ω2. C1 se reduce a C2 cuando N es constante** (nivel 1).
- El rango cíclico es β₁ = E − N + c.
- Con N y c fijos, conservar β₁ es conservar E.
- La conservación «topológica» de ciclos es una conservación de densidad disfrazada.

**Ω3. La densidad no determina d, pero puede codificarla a través de los puntos fijos de la regla** (nivel 2).
- En L-DIM-1 (NC-1) hubo geometrías de d = 2, 3 y 4 con el mismo grado. Así que fijar k̄ **restringe, no codifica**.
- **Pero** si los puntos fijos de una regla son rígidos (tipo retículo), el grado queda ligado a d. En Z^d, k = 2d. Entonces, con E conservado, d = k̄/2 queda **fijado por la condición inicial**.
- El caso límite son los complejos cúbicos CAT(0) cuyos enlaces son esferas octaédricas: localmente, Z^d.
- **El criterio de auditoría no es «¿se conserva E?», sino «¿determina el valor conservado la d de los puntos fijos?».** Se responde analizando los puntos fijos y, operativamente, con la prueba del dial (E2) sobre el valor conservado.
- La caracterización «complejo cúbico CAT(0) simplemente conexo con todos los enlaces octaédricos ⇒ cubulación estándar de ℝ^d» es una afirmación de nivel 3 **pendiente de verificación bibliográfica**. Ω3 no depende de ella: basta el ejemplo de Z^d con k = 2d.

**Ω4. El grano grueso (renormalización) no supera el límite de resolución.**
- Agrupar bloques de radio b equivale a mirar escalas mayores. Lo que no resuelve un grafo de tamaño N tampoco lo resuelve su versión de grano grueso (O3, W2).
- Es una medida útil de **estabilidad de clase**, pero es otra forma de la prueba de escalado, no una prueba independiente.
- Además, cualquier regla de bloqueo introduce una escala b. Como **medición** es admisible; como mecanismo, no.

**Ω5. Una «región abierta de parámetros» debe combinarse con la prueba del dial.**
- Si, dentro de la región, la d efectiva cambia con el parámetro, la «fase» es un dial (L-ARQ-T1).
- **Condición de fase:** misma clase **y** misma d_eff (dentro de la tolerancia de Nivel II) en toda la región.

**Ω6. Hace falta un nulo del espacio de reglas** (comparaciones múltiples).
- Si se prueban M reglas y una «funciona», hay que informar la fracción de reglas del espacio enumerado, **fijado antes**, que deja residuo no excluido.
- Es el análogo del nulo A de L-DIM-1. Un éxito aislado en un espacio grande no es específico.

**Ω7. W_ij ~ U(0,1) sobre el grafo completo es un inicio de clique.**
- La prueba de olvido (E del acta) exige inicios **estructuralmente distintos**: completo con pesos, ER disperso, árbol, anillo y retículo.
- Además, con una conservación, los inicios deben compartir el valor conservado. Eso acopla Ω7 con Ω3: si se conserva E, el olvido nunca puede abarcar densidades distintas.

## 1. Espacio de reglas considerado

- Grafos simples sin pesos, con **N fijo** (C0).
- Una regla es un par de **disparadores locales**:
  - uno de borrado de aristas;
  - otro de creación de aristas, entre vértices a distancia ≤ 3.
- Ambos dependen solo de invariantes combinatorios del entorno de radio ≤ 2:
  - t(u,v) = número de vecinos comunes;
  - s(u,v) = número de 4-ciclos que pasan por la arista;
  - el grado.
- Los umbrales son **enteros** (no hay diales continuos).
- La actualización es irreversible: hasta la absorción, o con un forzamiento declarado.

## 2. Auditoría de invariantes (conservación)

| Clase | Ejemplo | Dictamen |
|---|---|---|
| C0: N constante | Todas las reglas del §1 | **Permitida** |
| C1: ciclos | β₁ constante | **Se reduce a C2** (Ω2) |
| C2: E constante | Deslizamiento uv → uw | **Rechazada en la primera prueba.** Riesgo Ω3: con puntos fijos rígidos, d = k̄/2 queda heredado del inicio, y el olvido queda limitado (Ω7) |
| C3: grados | Intercambio doble | **Rechazada.** Sin sesgo, la medida máxima de entropía con grados dados es un grafo aleatorio localmente arbóreo: un expansor (nivel 3, configuración / RR). Con sesgo, Ω1 o Ω3 |
| C4: con parámetro ajustable | — | Rechazada (acta) |

**Decisión sobre la pregunta abierta del Consejo (§20 del acta): qué se conserva.**
- **Solo N.** La densidad, el grado y los ciclos son **salidas** de la dinámica, no restricciones.
- Es la lectura operativa de «invariante global emergente»: el invariante, si existe, debe **aparecer** en el atractor, no imponerse.
- **Motivo:** es la única opción con la que el resultado no puede heredar d por la vía Ω3.

## 3. Análisis previo de puntos fijos (verificado con el censo estático)

### Censo exacto

| Grafo | t en aristas | t en cuñas (pares a distancia 2) |
|---|---|---|
| Z² / Z³ / Z⁴ | 0 (100 %) | 1: 1/d · 2: (d−1)/d (0.50 / 0.33, 0.67 / 0.25, 0.75) |
| Triangular | 2 (100 %) | 1: 0.5 · 2: 0.5 |
| RGG3 k12 | Distribución ancha (0–9+) | Distribución ancha |
| Árbol, ER k4, RR k6 | 0 (≥ 99.4 %) | 1 (≥ 99.7 %) |
| 2-árbol | 1–8 | 1: 0.92 · 2: 0.08 |

### R1-T1 (nivel 1 sobre el censo): las reglas de umbral sobre t no pueden olvidar inicios arbóreos

- Para que todas las Z^d (d ≥ 2) sean puntos fijos, que es lo mínimo para un conjunto de puntos fijos ciego a d:
  - el disparador de creación sobre cuñas debe excluir t ∈ {1, 2};
  - el de borrado sobre aristas debe excluir t = 0.
- Entonces los árboles, ER y RR, con aristas t = 0 y cuñas t = 1 salvo una fracción ≤ 0.6 %, también son (casi) puntos fijos.
- **Ninguna regla de umbral puro sobre t borra la memoria de un inicio localmente arbóreo** sin destruir al mismo tiempo algún Z^d.
- Con información de radio 1 (t), geometría y árbol son indistinguibles localmente. Es la versión dinámica de O4.

### R1-T2 (nivel 1): toda regla que apunte a una **fracción** de cuñas cerradas en cuadrados codifica d

- En Z^d esa fracción es (d−1)/d (censo).
- Un objetivo «fracción = f» selecciona d = 1/(1−f): es codificación, como el objetivo de 4-ciclos C1 de la Fase 2.
- **Solo son admisibles condiciones de presencia** (≥ 1), no de valor ni de fracción.

### R1-T3 (nivel 2): el disparador mínimo ciego a d que separa geometría de árbol es la presencia de cuadrados por arista

| Grafo | s(u,v) ≥ 1 |
|---|---|
| Z^d, d ≥ 2 | Todas las aristas |
| Árbol | Ninguna |
| ER / RR dispersos | Fracción → 0 |

- La regla «**borrar las aristas que no estén en ningún 4-ciclo**» desestabiliza los inicios arbóreos y deja fijas todas las Z^d. Es ciega a d porque usa presencia, no número.
- **Pero su conjunto estable también contiene degenerados:**
  - tubos y escaleras 1D (C × P₂: toda arista en un cuadrado);
  - productos de árboles (no amenables, δ ≈ 0);
  - C × C₅.
- Coincide con el locus (ii) de L-ARQ-0. La condición de cuadrados es necesaria para la conmutatividad, pero admite rango 1 (tubos) y ramificación (productos de árboles).
- **El juez A− (X1, X3) es el que debe decidir.** Se prevé que no basta.

### Restricciones sobre el disparador de creación (derivadas)

- **(i)** No debe activarse en ningún Z^d, porque si no Z^d deja de ser fijo.
  - Ejemplo: completar cuadrados sobre caminos u–v–w–x crearía atajos de distancia 3 en Z^d (es bipartito). Queda **excluido**.
- **(ii)** Debe evitar el vacío: si solo se borra, el estado absorbente desde ER es casi vacío.
- **Disparador mínimo que cumple (i) y (ii) y es ciego a d:** «un vértice de grado ≤ 1 se une a un vértice uniforme a distancia 2; si está aislado, a un vértice uniforme del grafo».
  - Ningún Z^d tiene vértices de grado ≤ 1.
  - La unión de un vértice aislado es **no local**: es un forzamiento declarado.

## 4. Familia candidata para R1-1 y predicción congelada

**Regla SQ** (sin parámetros continuos; N fijo; actualización por arista o vértice al azar):
1. Borrar toda arista con s(u,v) = 0.
2. Todo vértice de grado ≤ 1 se reengancha según §3(ii).

Se detiene por absorción, o tras T barridos (T se fija en R1-1).

| Aspecto | Contenido |
|---|---|
| Espacio de reglas para el nulo Ω6 | Variantes discretas fijadas de antemano: umbral de presencia s ≥ 1 frente a t ≥ 1 frente a ambos; reenganche a distancia 2 frente a 3; tres órdenes de actualización |
| Auditoría | C0 (solo N). Sin d, k, ρ ni L. Solo presencia (R1-T2). Z^d de cualquier d es punto fijo (por diseño, ciego a d) |
| **Predicción (nivel 2)** | Desde inicios ER o árbol, el borrado masivo (casi ninguna arista está en un 4-ciclo) fragmenta el grafo. El reenganche recrea un grafo aleatorio disperso; se forma un régimen estacionario de **bosque aleatorio con rotación**. La nucleación de racimos ricos en cuadrados exige que dos reenganches cierren un cuadrado, lo que es improbable (∝ 1/N). **Destino previsto: degenerado** (X1/X2) |

- La predicción es negativa. **Su valor es decisivo:** SQ es la regla mínima que satisface todas las restricciones derivadas.
- **Si falla,** el fallo quedará atribuido a la cinética (nucleación), no a la regla de estabilidad. Eso señalaría que R1 necesita un mecanismo de **nucleación o crecimiento** (P3), cuyos destinos ya están cerrados, o estructura adicional (R3), tal como anticipó el auditor del Consejo.

## 5. Reglas de protocolo para R1-1 y R1-2 (heredadas del acta y de este análisis)

| Id | Regla |
|---|---|
| **G1** | La dinámica no ve RC-3, d_eff, δ ni ninguna salida de un instrumento |
| **G2** | Espacio de reglas enumerado y congelado; se informa la tasa de éxito sobre todo el espacio (Ω6) |
| **G3** | Inicios: completo con pesos (W ~ U(0,1) umbralizado), ER k4, árbol de Prüfer, anillo y Z³. Prueba de olvido: clase final independiente del inicio |
| **G4** | N, 2N, 4N y 8N (N = 5·10³); semillas 0–4; tres órdenes de actualización |
| **G5** | Clasificación: primero A− (RC-3 con sus reglas de uso), después residuo; d_eff solo al final, y nunca como criterio |
| **G6** | Atractor = misma clase en ≥ 80 % de semillas × inicios × órdenes; si no, multiestabilidad, que se informa |
| **G7** | Un resultado no clasificable se informa como **residuo no clasificado** y vuelve a análisis conceptual |

## 6. Decisión del cerebro (pasos lógicos)

1. **R1-0 queda cerrado** con este documento: Ω1–Ω7, auditoría C0–C4, R1-T1..T3, familia SQ.
2. **Siguiente paso: prerregistro completo de R1-1** (los once puntos del auditor) con la regla SQ y su espacio de variantes. Después se ejecuta y se juzga con A−. No se introduce ninguna variante después de ver resultados.
3. **Rama de decisión prerregistrada:**
   - **SQ degenerada en todas las variantes:** R1 queda cerrada para reglas locales de estabilidad. El siguiente candidato es **R3** (orden o estructura adicional), con su auditoría de codificación (E6).
   - **Residuo no excluido estable (G6):** R1-2 (prueba de fase completa: región abierta con d constante, Ω5) y solo después la medición de d.
