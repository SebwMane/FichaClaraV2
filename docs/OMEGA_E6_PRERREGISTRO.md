# Ω — E6: criterio de no codificación dimensional (prerregistro)

Estado: **PRERREGISTRO**, escrito antes de cualquier código de esta fase.
Rama de trabajo: `claude/omega-e6`. Se congelará como `claude/omega-e6-congelado`.

Origen:
- El Consejo autorizó R3-0 y aprobó W5 como regla de validez futura.
- Encargó diseñar E6: decidir si un orden de agregación contiene una codificación oculta de la dimensión antes de que aparezca la geometría.
- Este documento recoge la auditoría de la propuesta del Consejo (§1), el criterio de decisión (§2–§5) y la calibración que se ejecuta ahora (§6).

Prohibiciones vigentes (sin cambios):
- Nada de coordenadas ni embebidos dentro de mecanismos Ω. Los generadores del panel sí pueden usarlas.
- d = 3 o k = 6 nunca es un objetivo.
- Los umbrales no se cambian después de ver resultados.
- RNG solo con `rng_from_key`.

---

## 1. Auditoría de la propuesta del Consejo

La propuesta es acertada en lo esencial:
- no basta con que «la regla no mencione d»;
- hace falta un control positivo (P0);
- el juez debe poder decir «indeterminado»;
- lo que se busca es información dimensional *previa* a la dinámica.

Se adoptan esas cuatro ideas. Hay, sin embargo, problemas y omisiones que cambian el diseño. Cada uno se justifica a continuación.

| # | Punto de la propuesta | Problema | Decisión del cerebro |
|---|---|---|---|
| A1 | Las firmas se formulan como «R3 → 3D» | Contradice la prohibición: 3D no es objetivo. Un criterio asimétrico en d es en sí un sesgo | E6 se formula **simétrico en d**. Lo prohibido es *seleccionar cualquier d* desde la regla, no «producir 3D» |
| A2 | E6-B / E6-III: clasificador sobre descripciones de reglas, con AUC ≈ 0.5 | (i) Necesita una población de reglas con d emergente conocida, y las etiquetas exigen ejecutar dinámica y certificar d, para lo que no existe instrumento A+. Es circular. (ii) Sirve para una familia de reglas, no para un candidato único. (iii) «El clasificador no encuentra información» es ausencia de evidencia, con potencia desconocida | Se sustituye por una versión operativa y aplicable a una sola regla: **aplicar la regla como función local sobre un panel de geometrías de d conocida** y medir si su veredicto local distingue d (§3, E6-S). Mide exactamente I(D; R) antes de la dinámica, sin población de reglas ni etiquetas dinámicas |
| A3 | N6 (destruir correlaciones): «si N6 pierde 3D, eso favorece la emergencia» | **No discrimina.** Destruir correlaciones destruye por igual una codificación de red (los pares de interfaces son correlaciones) y una emergencia. Una regla codificada también «pierde 3D» bajo N6 | Regla general: **una familia nula solo tiene peso decisorio si su resultado predicho difiere entre un control codificado (P0+) y uno ciego (P0−)**, y esa diferencia se verifica en calibración. N6 queda sin peso decisorio mientras no se demuestre lo contrario |
| A4 | E6-F / N5: MDL, L(R \| D = 3) | La complejidad de Kolmogorov no es computable. Cualquier aproximación depende de un lenguaje de referencia elegido a mano, y L(R \| D) no está definida sin él | Sin peso decisorio. Se conserva como pregunta de auditoría escrita: «¿se describe la regla más corto en términos de d direcciones?» |
| A5 | Firmas con «D_eff = 3.0 ± 0.1» y I(D; R_final) | **No existe instrumento A+.** RC-3 solo excluye (A−) y tiene sesgo dimensional conocido (W1). La d de salida no se puede medir con esa precisión hoy | La rama dinámica (§5) solo compara clases RC-3 y la distribución de δ (≈ 1/d en redes; estadístico no certificante) entre la regla y sus nulos. Ninguna afirmación de certificación hasta que exista A+ |
| A6 | Iguales que se escriben como «≈ R3» | La no-significación no prueba igualdad. Además hay 8 familias × varios estadísticos | Igualdad solo por **equivalencia (TOST) con margen prerregistrado**. Corrección de Holm sobre las familias con peso decisorio |
| A7 | E6-D (selector oculto q = 1, 2, 3) | Es el criterio E2 (dial) y el lema del dial L-ARQ-T1, ya establecidos | Se unifica: **toda constante entera o umbral de la regla se barre** sobre una rejilla prerregistrada. Si el conjunto seleccionado se mueve con la constante, el resultado es DIAL (§3.3) |
| A8 | No se distingue «favorece redes regulares» de «selecciona d» | Una regla puede fijar la *reticularidad* sin fijar d. Ejemplo: «cuadrados por vértice = grado·(grado − 2)/2» se cumple en Z^d para **todo** d y no nombra ninguna constante dimensional | E6 no debe castigar la reticularidad como si fuera codificación de d. La calibración incluye este caso (P0-LAT), con predicción CIEGO |
| A9 | No hay panel de referencia | Sin geometrías ricas en triángulos en cada d, un criterio de «arista en triángulo» parecería seleccionar d = 2, 3 (artefacto). Es la lección 17 otra vez: un estadístico local se lee contra lo que el panel contenga | El panel cubre, **en cada d = 1…6**, una red bipartita (Z^d), una geometría rica en triángulos (RGG_d) y una red con diagonales (§3.1) |
| A10 | No aparecen las cantidades conservadas | Ω3 de R1-0: conservar E codifica d por k = 2d en puntos fijos rígidos | Las cantidades conservadas entran en la auditoría del alfabeto (S0) y su valor se barre como cualquier constante |
| A11 | **Omisión principal: falta la prueba de olvido** | Una regla estáticamente ciega puede simplemente *conservar* la d con la que empieza. En R1-1, Z³ quedó fija bajo B-s. Eso no es generar dimensión, es heredarla. Ninguna de N1–N8 lo detecta, porque todas cambian la regla y no el inicio | Prueba obligatoria D1 (§5): **arrancar desde geometrías de d distinta**. Una regla generadora debe llevar inicios con d₀ = 2, 3, 4 al mismo estado final. Si el estado final recuerda d₀, la clase es CONSERVADORA, no emergente |
| A12 | Reglas con etiquetas (interfaces, estados, orientaciones) | El panel sin etiquetas no las evalúa | Auditoría estructural S2 (§3.4): un alfabeto con pares inversos que conmutan es una presentación de Z^r (Bass–Guivarc'h: crecimiento de rango r) y es FAIL directo. Para órdenes, el producto de r cadenas tiene dimensión de orden r (Dushnik–Miller) y entra en el panel de R3 |
| A13 | «Si R3 necesita estructura 2D/3D explícita para coalescer, cuenta como negativo» | Correcto | Adoptado como regla: E6 FAIL ⇒ el resultado de R3 es **negativo**, no un «éxito con salvedad» |

---

## 2. Arquitectura de E6

E6 tiene dos brazos y una calibración obligatoria.

| Brazo | Cuándo | Pregunta | Peso |
|---|---|---|---|
| **E6-S** (estático) | Antes de cualquier dinámica | ¿El veredicto local de la regla distingue d sobre geometrías conocidas? | Puerta: SELECTOR o DIAL ⇒ FAIL y la regla no se ejecuta |
| **E6-D** (dinámico) | Solo si E6-S pasa y la dinámica produce estructuras no excluidas por RC-3 | ¿El estado final depende del inicio, de una constante, del orden de actualización o de la representación? | Clasifica: EMERGENTE-CANDIDATO / CODIFICADO / CONSERVADOR / INDETERMINADO |
| **Calibración** | Antes de usar cualquiera de los dos brazos | ¿El brazo detecta codificaciones conocidas y deja pasar reglas ciegas conocidas? | Si falla, el brazo es INVÁLIDO y no puede usarse |

E6-S pasado es condición **necesaria, no suficiente**.
- Una regla estáticamente ciega que, tras la dinámica, produce una d estable, olvida el inicio y no tiene dial es justamente el caso interesante: la información dimensional aparece en la dinámica y no estaba en la regla.
- Es la formulación operativa de la frase del Consejo «I(D; R) ≈ 0 e I(D; R_final) > 0».

---

## 3. E6-S: brazo estático

### 3.1 Panel de referencia (congelado)

Para cada d = 1…6, todas las geometrías periódicas (sin borde):

| d | Geometrías |
|---|---|
| 1 | ciclo C_10000 (Z¹); circulante C_10000(±1, ±2, ±3) (rico en triángulos); RGG₁ anillo k = 10 (N = 10000) |
| 2 | Z² toro 100²; triangular toro 90²; panal toro 100 × 100; RGG₂ k = 8 y k = 16 (N = 8000); Z² + 20 % diagonales (lado 90) |
| 3 | Z³ toro 20³; FCC (L = 12); BCC (L = 10); diamante periódico (≈ 8000 nodos); RGG₃ k = 8 y 16; Z³ + 20 % diagonales (lado 20) |
| 4 | Z⁴ toro 10⁴; RGG₄ k = 8 y 16; Z⁴ + 20 % diagonales (lado 10) |
| 5 | Z⁵ toro 6⁵; RGG₅ k = 8 y 16; Z⁵ + 20 % diagonales (lado 6) |
| 6 | Z⁶ toro 5⁶; RGG₆ k = 8 y 16; Z⁶ + 20 % diagonales (lado 5) |

Notas:
- Lado ≥ 5 en todos los toros, para que la periodicidad no cree 4-ciclos espurios.
- Los grafos RGG y con diagonales usan 2 semillas cada uno: `rng_from_key((MASTER_E6, id_geometría, semilla))`, con `MASTER_E6 = 20261020`.

Resolución declarada:
- Con d = 1…6, una selección de d₀ ∈ {2, 3, 4, 5} es detectable como SELECTOR.
- Una preferencia por d = 1 o d = 6 solo aparece como MONÓTONA (efecto de borde del panel).

### 3.2 Estadístico

Una regla Ω expone un **veredicto local** `reposo / acción` sobre un vértice o una arista, calculado solo con la estructura del grafo (sin coordenadas). Para cada geometría g:
- a(g) es la fracción de sitios en reposo.
- Reglas de vértice: todos los vértices. Reglas de arista: 3000 aristas muestreadas con `rng_from_key((MASTER_E6, id_geometría, semilla, 7))`.
- Para geometrías con semillas, a(g) es la media entre semillas.

Por dimensión: **m(d) = max_{g de dimensión d} a(g)**. Responde a «¿existe alguna geometría de dimensión d en la que la regla descansa?».

También se reportan ā(d), la media, y la tabla completa; esas columnas no deciden.

### 3.3 Clases (umbrales congelados)

- R = {d : m(d) ≥ 0.9} (dimensiones de reposo).
- Q = {d : m(d) ≤ 0.5} (dimensiones de acción).

Las clases se asignan en este orden:
1. **CIEGO:** max m − min m ≤ 0.10.
2. **SELECTOR:** existe d ∈ R y d′ < d < d″ con d′, d″ ∈ Q (la regla descansa en un d acotado por ambos lados).
3. **MONÓTONO:** R no vacío, R es un intervalo que contiene d = 1 o d = 6, y toda d fuera de R está en Q.
4. **INDETERMINADO:** cualquier otro caso.

**DIAL** (barrido de constantes, criterio A7). Si la regla tiene una constante entera o umbral c:
- Se evalúa sobre su rejilla prerregistrada.
- La familia es DIAL si los conjuntos R no vacíos toman **≥ 2 valores distintos** y al menos uno de esos valores de c da SELECTOR.

### 3.4 Auditoría estructural (S0–S2), por escrito, antes de ejecutar

- **S0, inventario.** Constantes enteras y umbrales, cantidades conservadas, tamaño del alfabeto, aridades y orden de actualización. Cada constante y cada cantidad conservada necesita su rejilla de barrido (DIAL).
- **S1, explícita.** Si aparece d o una función declarada de d (2d, d + 1, 3^d − 1, …), es FAIL.
- **S2, estructural.** Si el alfabeto de etiquetas contiene pares inversos {x, x⁻¹} con relaciones de conmutación, presenta Z^r, y es FAIL.
  - Para órdenes, si el reposo exige productos de cadenas o r órdenes lineales que se intersectan, es FAIL.
  - Si la regla tiene etiquetas, E6-S además se ejecuta con etiquetados del panel generados por (i) asignación uniforme y (ii) permutación de las etiquetas (N1). Esa extensión la define el prerregistro de cada candidato.

### 3.5 Decisión de la puerta

| Clase E6-S | Decisión |
|---|---|
| SELECTOR o DIAL | **FAIL**: codificación estática; la regla no se ejecuta |
| CIEGO | PASA |
| MONÓTONO | PASA CON AVISO: la preferencia empuja hacia un extremo degenerado (clique o árbol / 1D). La rama dinámica debe mostrar que RC-3 no lo excluye |
| INDETERMINADO | No pasa sin un análisis escrito del perfil m(d) y una decisión del Consejo |

---

## 4. Calibración de E6-S (se ejecuta en esta fase)

E6-S solo se puede usar si detecta codificaciones conocidas y deja pasar reglas ciegas conocidas. Todas las reglas de calibración salen de resultados previos de Ω.

| Regla | Tipo | Veredicto «reposo» | Rejilla | Predicción |
|---|---|---|---|---|
| C-TRIV | vértice | siempre | — | CIEGO |
| C-SQ | arista | la arista está en un 4-ciclo (criterio de reposo de SQ, B-s; R1-T3) | — | CIEGO (circulante en d = 1; Z^d para d ≥ 2) |
| C-TRI | arista | la arista está en un triángulo (B-t) | — | CIEGO (hay geometrías ricas en triángulos en todo d) |
| P0-LAT(q) | vértice | sin triángulos en v y c₄(v) = deg(deg − q)/2 | q ∈ {1, 2, 3} | q = 2: CIEGO (reposo en Z^d para todo d). q = 1, 3: reposo casi nulo. Familia: **no** DIAL |
| P0-DEG(k) | vértice | deg = k (Ω3: k = 2d en puntos fijos rígidos) | k ∈ {2, 3, 4, 6, 8, 10, 12} | k = 3 → SELECTOR {2}; 4 → SELECTOR {2, 3}; 8 → SELECTOR {3, 4}; 10 → SELECTOR {5}; 12 → SELECTOR {3, 6}; 2 y 6 → MONÓTONO. Familia: **DIAL** |
| P0-SQV(c) | vértice | c₄(v) = c | c ∈ {4, 12, 24, 40} = 2d(d − 1) para d = 2…5 | SELECTOR en {2}, {3}, {4}, {5} respectivamente. Familia: **DIAL** |
| P0-SQF(f) | vértice | \|f₄(v) − f\| ≤ 0.01, con f₄ la fracción de pares de vecinos que cierran un 4-ciclo sin pasar por v (R1-T2) | f ∈ {2/3, 4/5, 6/7, 8/9} = 2(d − 1)/(2d − 1) | SELECTOR en d = 2…5. Familia: **DIAL** |

Notas sobre las reglas:
- c₄(v) = Σ_{a < b ∈ N(v)} (\|N(a) ∩ N(b)\| − 1) es el número de 4-ciclos que pasan por v.
- En Z^d (lado ≥ 5): deg = 2d, c₄ = 2d(d − 1), f₄ = 2(d − 1)/(2d − 1).

**Condiciones de validez (congeladas):**
- **V1.** P0-DEG y P0-SQV dan DIAL.
- **V2.** P0-SQV(12) da SELECTOR con 3 ∈ R, y P0-SQV(24) da SELECTOR con 4 ∈ R.
- **V3.** C-TRIV, C-SQ, C-TRI y P0-LAT(2) no dan SELECTOR, y la familia P0-LAT no da DIAL.

Veredicto:
- **E6-S VÁLIDO** si se cumplen V1, V2 y V3.
- Si no, **E6-S INVÁLIDO**: no se ajustan umbrales; se rediseña como una fase nueva.
- Las predicciones por regla fuera de V1–V3 (P0-DEG por valor de k, P0-SQF) son predicciones del cerebro y no deciden. Se registran para contrastarlas.

---

## 5. E6-D: brazo dinámico (diseño congelado; se ejecutará con el primer candidato R3)

Solo se aplica a una regla que pasó E6-S y cuya dinámica deja estructuras **no excluidas** por RC-3, con pares válidos según W5.

| Prueba | Qué cambia | Qué conserva | Estadístico | Peso |
|---|---|---|---|---|
| **D1 olvido** (A11) | Inicio: Z², Z³, Z⁴ y además los inicios I1–I4 de R1-1 | Regla | Clase RC-3 y δ finales | **Decisorio.** Si el δ final depende de d₀ (TOST con margen 0.04 fallido y diferencia significativa) ⇒ CONSERVADOR |
| **D2 dial** (E2, A7) | Cada constante y cada cantidad conservada sobre ≥ 4 valores | Todo lo demás | δ final | **Decisorio.** δ monótono con la constante y ≥ 2 clases ⇒ CODIFICADO |
| **D3 orden** (N8) | SYNC / ASYNC / RAND / orden invertido | Regla | δ, clase | **Decisorio.** Si el resultado solo aparece en un orden ⇒ INDETERMINADO, y se localiza la causa |
| **D4 representación** (N1, N2) | Permutación de etiquetas; renumeración de vértices; segunda implementación independiente | Estructura | δ, clase | **Decisorio.** Si no hay equivalencia (TOST, margen 0.04) ⇒ error de implementación o codificación por etiquetas |
| **D5 tamaño** (E4, W5) | N, 8N (y 2N, 4N si procede) | Todo | δ, clase, E/N | **Decisorio.** Solo cuentan pares con E/N dentro del 25 % |
| D6 nulos de inventario (N3, N4, N6) | Reasignación de aridades, frecuencias, correlaciones | Inventario | δ, clase | Peso **solo si** la calibración dinámica muestra que discriminan P0+ de P0− (A3) |
| **Calibración dinámica** | Dinámica codificada conocida: recableado hacia grado k = 2d₀ conservando N, con d₀ ∈ {2, 3, 4} | — | ¿δ sigue a k? | Si D1–D2 no la marcan CODIFICADO, **E6-D es INVÁLIDO** |

Clases (orden de aplicación):
1. **CODIFICADO:** D2 positivo o E6-S FAIL.
2. **CONSERVADOR:** D1 positivo.
3. **INDETERMINADO:** D3 o D4 fallan, o la calibración dinámica no se ejecutó.
4. **EMERGENTE-CANDIDATO:** todo lo anterior negativo, D5 válido y E6-S CIEGO o MONÓTONO.

Alcance de la afirmación:
- EMERGENTE-CANDIDATO es una afirmación de **nivel 1** sobre independencia del inicio, de las constantes, del orden y de la representación.
- **No certifica** ninguna d. Eso exige un instrumento A+ que no existe.
- Holm se aplica sobre D1–D5.

---

## 6. Ejecución de esta fase

1. Código: `tools/e6_static.py` (panel §3.1, reglas §4, clases §3.3), `tests/test_e6_static.py` y resultados en `results/e6/`.
   - Implementación por un agente Sonnet; revisión del cerebro.
   - Los tests verifican, en Z^d pequeño, los valores analíticos deg = 2d, c₄ = 2d(d − 1) y f₄ = 2(d − 1)/(2d − 1).
   - También verifican las clases sobre perfiles sintéticos de m(d).
2. Ejecución de la calibración y veredicto VÁLIDO / INVÁLIDO según V1–V3.
3. Resultados en `docs/OMEGA_E6_RESULTADOS.md`; actualización del mapa negativo y de `OMEGA_ARQ_0.md` §4 (E6 sustituido por esta definición).
4. Solo entonces empieza **R3-0**, que debe pasar sus candidatos de orden por S0–S2 y, si son evaluables, por E6-S.

## 7. Predicciones del cerebro (registradas antes de ejecutar)

- E6-S será VÁLIDO.
- Riesgo principal: que algún P0-SQV coincida en una geometría ajena (FCC, BCC o triangular tienen muchos 4-ciclos) y convierta un SELECTOR en INDETERMINADO o MONÓTONO. Si ocurre en c = 12 o c = 24, V2 falla y E6-S es inválido. Ese riesgo se acepta.
- C-SQ dará CIEGO. Con ello queda registrado que el fracaso de SQ en R1-1 **no** se debió a codificación de d.
