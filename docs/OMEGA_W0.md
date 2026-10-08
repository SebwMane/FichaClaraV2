# Ω — W-0: ¿puede autoseleccionarse la anchura w? (análisis sin dinámica)

- **Fecha:** 2026-10-07.
- **Rama:** `claude/omega-w0`. Base congelada: `claude/omega-r3-0-congelado`.
- **Mandato del Consejo.** R3 reducido queda ratificado, W6 adoptado y W-0 autorizado como análisis conceptual. No hay dinámica nueva.
- **Pregunta.** ¿Puede emerger w sin estar codificado como parámetro, condición inicial o cociente de tasas?
- **Categorías:**
  - **C1:** w fijado explícitamente (codificación).
  - **C2:** w fijado por cociente de tasas o parámetros (dial).
  - **C3:** w seleccionado por una inestabilidad o un principio estructural sin parámetro que lo controle (candidato).

---

## 0. Acta del Consejo y auditoría

### 0.1 Ratificado

| Decisión | Registro |
|---|---|
| **R3 reducido**, con la formulación del Consejo: «R3 demuestra que la coalescencia sin defectos puede producir una dimensión efectiva igual a la anchura del orden. No demuestra un mecanismo que seleccione espontáneamente esa anchura» | Adoptada literalmente (mapa negativo, fila 13) |
| **Lema falsado:** «que un juez basado en un único par de escalas pueda detectar de forma fiable el crecimiento intermedio mediante la explosión de cortes» | Adoptado literalmente (mapa negativo) |
| **W6 adoptado** como condición de uso de RC-3, con su definición congelada antes del próximo candidato | Formalizado en `docs/OMEGA_W6_PRERREGISTRO.md`, con calibración y potencia |
| No cerrar el programa de selección. No hay teorema de imposibilidad: «los mecanismos examinados producen d como parámetro, recuento, condición inicial o dial» | Adoptado. W-0 **clasifica**; no pretende demostrar imposibilidad |
| W-0 no busca «la regla que produzca 3» | Adoptado. Además, todo resultado de W-0 se formula para un d genérico (prohibición vigente) |
| Sin nueva dinámica por ahora | Adoptado |

### 0.2 Omisiones (justificadas)

**Ω-1. C3 necesita un criterio operativo.** Sin él, «sin parámetro que lo controle» es una afirmación no comprobable. Criterio propuesto, que une criterios ya existentes:
- un w está en C3 si es **invariante** frente a todo parámetro continuo (E2, barrido de ≥ 4 décadas);
- y frente a toda constante entera de la regla (E6, DIAL), **incluida la elección de la ley o propiedad** que define el mecanismo.

Esto último es lo que el Consejo no menciona. Es decisivo en §3.

**Ω-2. Los umbrales estructurales son de un solo lado.** Las propiedades sin parámetro que conocemos con un entero crítico (recurrencia, intersección de trayectorias) valen para d ≤ d_c y fallan por encima.
- Una propiedad así solo selecciona d_c si existe además una **presión monótona hacia w mayor**, sin acotar, que se detenga exactamente donde la propiedad falla.
- Si la presión se equilibra con otro término, vuelve el lema del dial.

**Ω-3. Localidad.** Esas propiedades son leyes 0-1 asintóticas y globales. Una regla local en tiempo finito solo ve un cruce suave. Esto se trata en §3.3.

**Ω-4. El criterio de grado saturado necesita su potencia** (lección 25). Se calcula en la calibración de W6, no se supone.

---

## 1. Planteamiento

- Un **proceso de hilos** es un orden causal J que crece por tres operaciones: avance de un hilo, bifurcación (un hilo da dos) y fusión (dos hilos se unen; el elemento nuevo tiene dos predecesores).
- La anchura asintótica w* es el número de hilos que avanzan independientemente sin límite (R3-T1).
- Por R3-T1, la geometría del espacio de cortes tiene dimensión w*.
- La pregunta de W-0 es de qué depende w*.

## 2. C1 y C2: lo que ya está cerrado

**W-T1 (nivel 1). En bifurcación–fusión con tasas, el número de hilos lo fija un cociente de tasas.**

*Prueba.* Sea n el número de hilos activos, con bifurcación a tasa β por hilo y fusión a tasa γ por par.
- Es una cadena de nacimiento y muerte con tasa de nacimiento βn y de muerte γ n(n − 1)/2.
- Equilibrio detallado: π(n+1)/π(n) = 2β / (γ(n + 1)).
- La moda está en n ≈ 2β/γ. La distribución estacionaria es de Poisson con media 2β/γ condicionada a n ≥ 1.
- El número típico es una función continua de β/γ. Al barrer β/γ cuatro décadas recorre todos los enteros: es un **dial** (E2).

**Agravante.** Una distribución de Poisson no selecciona un entero, sino una distribución ancha. No hay un w*, sino un n que fluctúa.

**W-T2 (nivel 1; lema del dial L-ARQ-T1).** Toda selección de w como argmin de A(w) + θB(w) solo puede elegir vértices de la envolvente convexa inferior de (B_w, A_w), y el elegido cambia monótonamente con θ.
- Una meseta en θ es una arista de la envolvente, no una selección.
- Es C2, salvo que θ esté fijado por un principio independiente documentado antes de mirar resultados (cláusula α de E2).

**W-T3 (nivel 1). Fijar w con la condición inicial y conservarlo es C1.** Si las operaciones conservan el número de hilos (sin bifurcación ni fusión netas), w* = w₀. Es la clase CONSERVADOR de E6-D (prueba D1).

**Conclusión §2.** En procesos de hilos con tasas, conservación u optimización de un funcional con peso, w* es C1 o C2. Queda solo la vía de §3.

---

## 3. C3: umbrales estructurales sin parámetro

### 3.1 Los hechos (nivel 3; verificación en §6)

Hay propiedades intrínsecas de un grafo, definibles con paseos aleatorios **sobre el propio grafo** (sin coordenadas), que cambian a un entero crítico sin ningún parámetro continuo:

| Propiedad (en grafos de crecimiento polinómico de grado d, o en ℤ^d) | Vale si y solo si | Fuente |
|---|---|---|
| R1: el paseo aleatorio simple es recurrente | d ≤ 2 | Pólya (1921); Varopoulos para grupos de crecimiento polinómico |
| R2: dos paseos independientes **colisionan** (están en el mismo sitio al mismo tiempo) infinitas veces | d ≤ 2 (en ℤ^d) | Folklore (diferencia de paseos); Krishnapur y Peres (2004) para grafos |
| R3: las **trayectorias** de dos paseos independientes se cortan infinitas veces | d ≤ 4 | Erdős y Taylor (1960); Lawler |
| R4: las trayectorias de k paseos independientes tienen intersección común infinita | d < 2k/(k − 1) en el caso browniano; versión de red por verificar | Dvoretzky, Erdős, Kakutani y Taylor (1950–57) |

### 3.2 W-T4: en C3 la elección de la ley es un dial discreto (nivel 2)

**Argumento.** Supongamos un mecanismo «presión monótona hacia w mayor, detenida cuando falla la propiedad P».
- Selecciona w = d_c(P). Pero d_c depende de **qué propiedad** se elija:
  - recurrencia o colisión: 2;
  - intersección de 2 trayectorias: 4;
  - k trayectorias: decreciente en k, de 3 (k = 2, browniano) a 2 (k → ∞).
- La identidad de la propiedad, en particular su aridad k, es una constante entera del mecanismo.
- Por la regla DIAL de E6 (§3.3 del prerregistro E6), al barrer k cambia el conjunto seleccionado. **El mecanismo es DIAL**, salvo que k esté fijado por un principio independiente previo (cláusula α).

**El único resquicio identificado.** Un principio de **aridad mínima**: «la interacción relevante es la de pares, k = 2, porque es la mínima no trivial».
- Si ese principio es independiente o es una elección **no se puede demostrar dentro de Ω**: es una decisión de principio.
- Además, incluso con k = 2 queda otra elección discreta: colisión (d_c = 2) o intersección de trayectorias (d_c = 4, o 3 en el continuo). Son dos lecturas de «dos hilos se encuentran» que dan enteros distintos.

### 3.3 W-T5: obstáculo de localidad y tamaño finito (nivel 2)

- R1–R4 son leyes 0-1 sobre el comportamiento **asintótico** de paseos infinitos.
- A tamaño N finito, la probabilidad del evento es una función suave de d con correcciones lentas en d_c. Por ejemplo, en d = d_c para R3 la probabilidad de corte decae como 1/log, que es nivel 3 por verificar.
- Una regla local, que solo ve una vecindad y un tiempo finito, ve un cruce suave. El punto donde se detiene la presión depende de la escala temporal del mecanismo, es decir, de un parámetro. **A tamaño finito, C3 se degrada a C2.**
- Solo un mecanismo con información global o de tiempo infinito realizaría el umbral exacto. Eso choca con la localidad, que es un requisito de Ω desde C0.

### 3.4 W-T6: un umbral sin presión no selecciona (nivel 1)

- Una propiedad P que vale para todo d ≤ d_c es compatible con cualquier d ≤ d_c, incluida la clase 1D que RC-3 excluye.
- Sin presión hacia w mayor, el mecanismo cae en el extremo degenerado: E6-S lo clasificaría como MONÓTONO.
- Con presión, hay que justificar que no se equilibra con otro término (W-T2) y que no explota (crecimiento intermedio, W6).

---

## 4. Veredicto de W-0

**W-0: C3 NO VACÍA EN PRINCIPIO, REDUCIDA A UNA ELECCIÓN DE LEY.**

1. **C1 y C2 (nivel 1).** Todo proceso de hilos gobernado por tasas, conservación o un funcional con peso fija w* por parámetro, condición inicial o dial (W-T1 a W-T3).
2. **C3 existe como estructura (nivel 3).** Hay propiedades intrínsecas de grafos, sin coordenadas ni parámetros continuos, con entero crítico (§3.1).
3. **Pero C3 traslada la selección a la elección de la propiedad (nivel 2).** El entero depende de la ley (recurrencia, colisión, intersección) y de su aridad. Por la propia regla DIAL de E6 eso es un dial discreto, salvo que un principio independiente fije la ley (W-T4).
4. **A tamaño finito, C3 se degrada a C2 por localidad (nivel 2, W-T5).**
5. **Formulación del resultado estructurado, para el Consejo:**

   > Los mecanismos examinados producen d como parámetro, recuento, condición inicial o dial. La única vía sin parámetro continuo, los umbrales estructurales de paseos sobre el propio grafo, selecciona un entero crítico que depende de qué ley se elija, y a tamaño finito necesita una escala temporal que actúa de nuevo como dial. Que exista un principio que fije la ley es una cuestión de principio, no demostrable dentro de Ω.

## 5. Propuesta al Consejo

- **(a)** Ratificar el veredicto de W-0 con su alcance: clasificación, no imposibilidad.
- **(b)** Decidir si se abre **W-1**, una prueba estática y la última posible antes de cerrar.
  - **Pregunta:** ¿son detectables R1–R4 a los tamaños accesibles (10⁴–10⁶) con paseos intrínsecos sobre las geometrías del panel RC-3/E6 de d = 1…6?
  - **Diseño:** simétrico en d, con todas las leyes medidas y ninguna privilegiada. Contraste de potencia previo: margen entre d_c y d_c + 1 frente a la dispersión entre semillas.
  - **Si no son detectables:** C3 es inaccesible para un mecanismo Ω local a tamaños realizables, y el programa de selección se cierra con el resultado estructurado del §4.5.
  - **Si lo son:** la elección de ley queda como única pregunta abierta, que es de principio y corresponde al Consejo.
- **(c)** W6 se calibra en paralelo (`docs/OMEGA_W6_PRERREGISTRO.md`). No depende de W-0.

**Recomendación del cerebro:** (a) y (b). W-1 es barata, estática y decide si C3 es siquiera operativa. Sin ella, el cierre dependería de W-T5, que es solo de nivel 2.

## 6. Verificación bibliográfica

Agente Haiku con búsqueda web; revisión del cerebro.

| Id | Resultado | Decisión del cerebro |
|---|---|---|
| R1 | Pólya (1921): VERIFICADO. Varopoulos: VERIFICADO **solo para grupos finitamente generados** (grafos de Cayley); en grafos generales hace falta además una condición isoperimétrica | **Corregido.** En grafos generales de crecimiento polinómico solo vale una dirección sin hipótesis extra: un volumen ≲ r² con grado acotado implica recurrencia (criterio tipo Nash-Williams; recuerdo del cerebro). La transitoriedad para d > 2 exige isoperimetría. Consecuencia para W-1: el umbral R1 es fiable en redes y vértice-transitivos, no en cualquier grafo de crecimiento d |
| R2 | En ℤ^d: VERIFICADO. Krishnapur y Peres, *Recurrent graphs where two independent random walks collide finitely often*, Electron. Commun. Probab. 9 (2004) 72–81, arXiv:math/0406487 | Se añade un matiz importante: en grafos generales, ser recurrente **no** implica colisiones infinitas (el peine es un contraejemplo). R2 depende de algo más que la dimensión de crecimiento, y por tanto sería un mal selector incluso en principio |
| R3 | Erdős y Taylor (1960) y Lawler, *Intersections of Random Walks* (1991): VERIFICADO (d ≤ 4) | Se mantiene. Cita correcta: Acta Math. Acad. Sci. Hungar. 11 (1960). La URL que dio el agente (Comm. Math. Phys. 86) **no corresponde** y se descarta. El decaimiento en d = 4 es logarítmico. La forma exacta que dio el agente, (π²/8)/log n, **no** se adopta: no hay fuente primaria y depende de la variante (uno o dos lados). W-T5 solo usa que el decaimiento es logarítmico |
| R4 | DEK(T): VERIFICADO para el caso browniano (d < 2k/(k − 1)). Versión de red: el informe se contradice («también vale» y, en la cita, «se anula si D > 2k/(k − 1)») | Lectura del cerebro, coherente con R3 (k = 2, d ≤ 4): en red, **d ≤ 2k/(k − 1)**, incluido el valor crítico; en el continuo es estricto. Queda como nivel 3 **no verificado con fuente primaria**. La diferencia red/continuo en k = 2 (4 frente a 3) es otra elección discreta, y refuerza W-T4 |

**Efecto:** el veredicto se mantiene y W-T4 sale reforzado. R2 depende de la geometría fina (peine) y R1 solo es un umbral limpio en grafos homogéneos. Las leyes de C3 no son funciones solo de la dimensión de crecimiento: también dependen de la homogeneidad del grafo. W-1, si se aprueba, debe incluir geometrías no homogéneas (peine, retazos) en el panel.

---

## 7. Enmienda del Consejo (tras su revisión; aplicada en W-1)

- La formulación del §4.5 queda **sustituida** por: «Entre las familias examinadas, las rutas sin parámetro continuo no han producido todavía un selector de dimensión; las candidatas basadas en leyes de paseo trasladan la selección a una elección discreta de ley, salvo que exista un principio independiente que determine esa ley».
- W-0 queda aprobada como **reducción, no como imposibilidad**.
- W-T4 se lee con la condición explícita del Consejo: una elección de ley es un dial solo si hay varias leyes admisibles, ninguna fijada por un principio independiente previo, y las elecciones dan resultados distintos.

---

## 8. Adenda tras la sesión del Consejo sobre «R3 reducido» (caso IV y dilema de independencia)

### 8.1 Auditoría del acta

**Desfase.** El acta vota «W-0 aprobado» como paso pendiente. En el repositorio, W-0 ya se ejecutó y se aprobó como reducción (§7), y W-1 también se ejecutó: W1-D (`OMEGA_W1_RESULTADOS.md`). Se toma el acta como revisión de fondo, no como orden de repetir W-0.

**Correspondencia de sus casos con W-0:**

| Caso del Consejo | Categoría de W-0 | Estado |
|---|---|---|
| I, dial | C2 (W-T1, W-T2) | Cerrado, nivel 1 |
| II, condición inicial | C1 / CONSERVADOR (W-T3) | Cerrado, nivel 1 |
| III, ley elegida | C3 → elección de ley (W-T4), con W-1 descriptivo en N₂ (corte 2 frente a 4) | Reducido |
| IV, atractor genuino | **No tratado de forma explícita en W-0** | Es la omisión que se cubre aquí |

**Lo genuinamente nuevo del acta es su §13, la «última trampa».** Si w* = 3 sale de una regla concreta, entonces «por qué 3» se responde con «porque la regla tiene su atractor ahí», y el problema solo se ha desplazado. Formalmente:
- Toda regla sin parámetros es un punto de un espacio discreto de reglas, y su atractor w*(regla) es una función sobre ese espacio.
- Elegir la regla equivale a elegir w*, que es W-T4 generalizado.
- **El caso IV solo escapa si w* es universal sobre una clase amplia de reglas** (como los exponentes críticos son universales frente a los detalles microscópicos), no si es el atractor de una regla.

Se adopta esa definición estricta del caso IV.

### 8.2 W-T7: el dilema de independencia (nivel 2)

Para que una dinámica local regule el número de hilos asintóticamente independientes, un hilo tiene que recibir información de los demás. Hay tres casos, exhaustivos según el desfase entre hilos:

1. **Sin comunicación.**
   - Ningún hilo puede saber cuántos hay, así que w queda fijado por el inicio (C1) o crece con la bifurcación sin cota.
   - En el segundo caso hay explosión o crecimiento intermedio, como en J3/J4 (R3-0b, nivel 1).
2. **Comunicación con desfase acotado.**
   - Los hilos quedan acoplados y el espacio de cortes es un tubo: w_eff = 1 (R3-T1, corolario 2; J2 da X3 en 3/3, nivel 1).
   - La regulación destruye justo la independencia que define w.
3. **Comunicación con desfase no acotado** (por ejemplo, difusivo).
   - El vector de desfases de w hilos es un paseo en ℤ^{w−1}. Los encuentros son eventos de recurrencia (Pólya, nivel 3).
   - El encuentro de un par es recurrente en una dimensión: un regulador de tipo «fusionar al contactar» lleva a w → 1.
   - Un evento de sincronización simultánea de m hilos es recurrente si y solo si m − 1 ≤ 2.
   - Cualquier regulador basado en encuentros hereda un umbral que **depende de la aridad del evento elegido**. Eso es W-T4, una elección de ley.
   - Un regulador temporal («bifurcar si aislado durante τ») introduce una escala, es decir, C2.

**Mezclas.** Si unos pares tienen desfase acotado y otros no, los hilos se agrupan en racimos acoplados. w_eff es el número de racimos con desfase mutuo no acotado, y el análisis se repite al nivel de racimo.

**Conclusión de W-T7 (nivel 2, no teorema).** En procesos de hilos con información local, toda vía de regulación cae en I (dial o escala), II (inicio) o III (ley o aridad del evento), o colapsa a w = 1. El caso IV requeriría universalidad de w* sobre clases de reglas que escapen al dilema. No conozco ninguna, y el dilema indica dónde tendría que estar: en un mecanismo de información **no local**, que choca con la localidad exigida a Ω desde C0.

### 8.3 Por qué no se ejecuta un «W-2» de búsqueda de universalidad

Un experimento que enumere una gramática de reglas sin parámetros y mida la distribución de w*:
- Tendría que elegir la gramática. Por la lección 23, cualquier universalidad sería relativa a esa gramática, igual que «codificar d» es relativo al panel.
- El dilema predice que cada regla cae en I, II o III, o colapsa.
- Y el Consejo pidió no inventar mecanismos indefinidamente.

**Decisión del cerebro:** no se ejecuta. Si en el futuro aparece un **principio independiente** que proponga una clase de reglas concreta, por ejemplo un principio de aridad mínima o uno no local con justificación propia, esa clase se prerregistraría como fase nueva.
