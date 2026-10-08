# Ω — L-CICLOS-0: redundancia cíclica multiescala (sin dinámica) y validación fuera de muestra de la regla combinada — Prerregistro

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-ciclos`.
- **Base congelada:** `claude/omega-coherencia-congelado`.
- **Mandato del Consejo:**
  1. congelar la fase de coherencia con tres conclusiones;
  2. **validar fuera de muestra** la regla combinada antes de cualquier dinámica;
  3. abrir **L-CICLOS-0**, conceptual y sin dinámica: ¿puede definirse una redundancia cíclica multiescala que excluya árboles, cliques y expansores sin imponer dimensión, escala ni densidad objetivo?
- **Estado:** este documento se commitea antes de escribir código y antes de generar ningún grafo del panel nuevo.

## 0. Contraste del dictamen y cierre de la fase de coherencia

| Afirmación del Consejo | Contraste | Decisión |
|---|---|---|
| «Coherencia + curvatura local ⇏ geometría» | **Más fuerte de lo demostrado.** Lo probado es que cada una *por separado* es insuficiente. La *combinación* no se ha refutado: en los paneles corridos, cada clase de falso positivo la atrapa al menos un diagnóstico. Es una combinación **no validada**, que no es lo mismo que una combinación refutada | Se corrige: «C ⇏ G; K ⇏ G; C ∧ K sin validar» |
| «K ayuda contra estructuras ramificadas» | Matiz: el árbol (f_neg ≈ 0.5) se midió en L-P3-0b con N = 1000, no en el panel P1-D.3. Que caiga fuera de la banda K (f_neg_max = 0.428, de P1-D.3) es una comparación **entre paneles**, *post hoc* | Se adopta con esa etiqueta |
| «El problema no es reconocer la geometría; eso ya lo conseguimos» | **Sobreafirmación.** El instrumento es VÁLIDO solo en su dominio (cerrado y casi homogéneo) y la regla combinada no está validada fuera de muestra | Por eso la validación (§3) va antes de cualquier dinámica |
| Reformulación de mi frase sobre los ciclos | Correcta. La mía era demasiado fuerte | Se adopta: «una geometría extendida parece requerir una organización no degenerada de ciclos y ramificación a través de escalas» |
| Hipótesis A (riqueza de ciclos) frente a B (los ciclos como síntoma) | Bien planteadas. L-CICLOS-0 solo puede evaluar A como **observable**; no distingue A de B | Se registra; B queda abierta |
| La clique como adversario principal | Correcto: «muchos ciclos ⇒ geometría» lo destruye la clique | Incluida; además, el **árbol de cliques** (ciclos abundantes y ramificación) |
| No convertir C en energía; no abrir «coherencia + ciclos» como dinámica | Coincide con CH-T5 y con la regla del preflight | Se adopta |

**Fase de coherencia congelada** (`claude/omega-coherencia-congelado`) con tres conclusiones:
1. **C_multiescala es útil pero insuficiente.** Corrige los bordes y los cruces lentos; los árboles críticos son amenables y pasan.
2. **K es complementaria pero insuficiente.** No redundante con ρ; ciega a los atajos escasos y a C0.
3. **La combinación debe validarse fuera de muestra antes de diseñar dinámica** (§3).

## 1. Formalización: redundancia cíclica multiescala como «un solo extremo a toda escala»

**Idea.** Contar ciclos no sirve (la clique tiene muchísimos), y rellenar ciclos a escala fija tampoco (el cactus de triángulos tiene todos sus ciclos rellenos por triángulos). Lo que distingue una geometría extendida de dimensión ≥ 2 es que, a **toda escala r**, existen ciclos que **rodean** cualquier región. Es decir, el anillo alrededor de un nodo es conexo.

**Definición.** Para un nodo u y una escala r, el anillo es

    A_r(u) = {x : r ≤ d(u, x) ≤ 2r}

y se toma el subgrafo inducido por A_r(u). Sean c_r(u) su número de componentes y f_r(u) la fracción de nodos del anillo en la componente mayor. Las clases son:

| Clase | Condición | Ejemplos |
|---|---|---|
| **CONEXO** (un extremo a esa escala) | El anillo es conexo: hay ciclos que rodean B_{r−1}(u) | Geometría de dimensión ≥ 2 |
| **DOS_EXTREMOS** | Dos componentes de tamaño comparable | Geometría 1D, incluidos los tubos gruesos |
| **RAMIFICADO** | Muchas componentes que crecen con r | Árbol, cactus, árbol de cliques |
| **SIN_VENTANA** | No existe 2r dentro de la ventana | Clique, expansor (la bola se satura antes) |

**Lo que no contiene:** ni dimensión, ni escala fija (se exige a todas las r de la ventana), ni densidad objetivo, ni coordenadas.

**Anclaje teórico (nivel 3, en verificación).** Es la versión de escala finita del **número de extremos** de un grafo (Freudenthal, Hopf, Stallings):
- los grupos finitamente generados tienen 0, 1, 2 o ∞ extremos;
- Z tiene 2;
- Z^d con d ≥ 2 tiene 1;
- los árboles y los grupos libres tienen ∞.

La relación con el relleno de ciclos (función de Dehn: cuadrática en Z^d con d ≥ 2, lineal en espacios hiperbólicos) queda como motivación, no como criterio.

## 2. Análisis previo

| Id | Enunciado | Nivel |
|---|---|---|
| **CIC-T1** | En Z^d con d ≥ 2, todo anillo r ≤ \|x\|_∞ ≤ 2r es conexo (se rodea por caminos de coordenadas). En Z, el anillo tiene 2 componentes | 1 |
| **CIC-T2** | En un árbol con ramificación ≥ 3 el número de componentes crece exponencialmente con r. En un árbol crítico (uniforme) el número de ramas que cruzan de r a 2r es una variable aleatoria de orden 1–pocas, a menudo > 1 | 1 / 2 |
| **CIC-T3** | **Cactus y árbol de cliques:** los ciclos son locales y la estructura a gran escala es arbórea. Los anillos se rompen en ramas: RAMIFICADO, aunque todos sus ciclos estén rellenos a escala fija. **Esto es lo que la redundancia multiescala añade frente al relleno a escala fija** (`omega.topology.betti.short_cycle_betti1`, que daría 0 ciclos sin rellenar en un cactus) | 1–2 |
| **CIC-T4** | **Clique y expansores:** su diámetro es ≤ 2 o ~log N, y la bola de radio 2r supera N/4 antes de que exista anillo. Quedan SIN_VENTANA, excluidos por **ausencia de escalas**, no por un umbral de densidad | 1–2 |
| **CIC-T5 (límite: no selecciona dimensión)** | Toda dimensión ≥ 2 da CONEXO, y también el grupo de Heisenberg y el plano hiperbólico, que es una variedad con un extremo. **Separa ≥ 2D de 1D y de los árboles, no 2 de 3.** El plano hiperbólico queda excluido por la coherencia (no amenable), no por esta observable. Por eso son complementarias | 2 |
| **CIC-T6 (riesgo)** | El retazo no es su objetivo: el reemparejado aleatorio conecta los anillos. Probablemente CONEXO. Lo atrapan γ y K | 2 |

**Respuesta a la pregunta del Consejo:** sí, puede definirse sin dimensión, escala ni densidad objetivo, siempre que sea válida empíricamente (§4). **No** selecciona dimensión (CIC-T5) ni es mecanismo.

## 3. Prerregistro A: L-CIC-0b (observable de anillos)

**Observable** (nuevo módulo `omega/diagnostics/annulus.py`):
- 200 centros u muestreados de la componente gigante; un BFS desde cada uno.
- r_w = mayor r con mediana de |B_r| ≤ N/4.
- Escalas evaluadas: r = 2 … ⌊r_w/2⌋, con un máximo de 12 escalas tomadas espaciadas si hay más. Si no hay ninguna, `SIN_VENTANA`.
- Para cada (u, r): c_r(u) y f_r(u) en el subgrafo inducido por A_r(u). Anillos con menos de 10 nodos se omiten.
- **Por grafo:**
  - f̃ = mediana de f_r(u) sobre todos los pares (u, r);
  - c̃ = mediana de c_r(u);
  - f2̃ = mediana de la fracción de nodos del anillo en las **dos** componentes mayores.

| Estado | Condición |
|---|---|
| `SIN_VENTANA` | ninguna escala evaluable |
| `CONEXO` | f̃ ≥ 0.95 |
| `DOS_EXTREMOS` | no CONEXO, c̃ = 2 y f2̃ ≥ 0.95 |
| `RAMIFICADO` | cualquier otro caso |

**Criterios de L-CIC-0b:**
- **Sensibilidad (V):** fracción de geometrías de dimensión ≥ 2 que salen CONEXO.
- **Exclusión (X):** fracción de grafos X que **no** salen CONEXO.

| Veredicto | Condición |
|---|---|
| **REDUNDANCIA-VÁLIDA** | sensibilidad ≥ 0.9 y exclusión ≥ 0.9 |
| **REDUNDANCIA-INVÁLIDA** | alguna < 0.7 |
| **REDUNDANCIA-PARCIAL** | cualquier otro caso |

Además se informa:
- **Independencia de la densidad:** RGG3 con k = 8, 12 y 20 deben dar el mismo estado.
- **Independencia de la dimensión:** las familias V de d = 2 y d = 3 deben dar el mismo estado.

## 4. Prerregistro B: validación fuera de muestra de la regla combinada (diagnósticos existentes, congelados)

**Regla RC-1** (fijada aquí, antes de ver el panel nuevo). Un grafo es **GEOMÉTRICO-RC1** si cumple las cuatro condiciones:

| # | Condición | Origen |
|---|---|---|
| (i) | **Localidad A:** fracción de aristas en algún triángulo o 4-ciclo ≥ 0.5 | Definición de `classify_c0` |
| (ii) | **Coherencia** γ: estado `COHERENTE` | L-COH-0b |
| (iii) | **Curvatura K:** mediana de κ ∈ [−0.0566, 0.2700] y f_neg ≤ 0.4283 | Banda de la auditoría (b), sin reajustar |
| (iv) | **Nivel II:** estado `CONVERGE` de D_B2 | P1-D.3 |

**Medidas:**
- **Sensibilidad:** fracción de V que es GEOMÉTRICO-RC1.
- **Especificidad:** fracción de X ∪ L (no geométricos) que **no** lo es.
- El veredicto usa los mismos umbrales: VÁLIDA / PARCIAL / INVÁLIDA.

**Fuera de muestra significa:** las familias marcadas como **nuevas** en el panel no se usaron en el diseño de ningún diagnóstico. Las familias repetidas (RGG, retículos, cajas) se marcan como «repetidas» y se informan aparte.

## 5. Panel nuevo (N ≈ 2·10⁴, disperso; clave (20261013, familia, semilla); semillas 0–2 en las aleatorias, 1 en las deterministas)

| Grupo | Familias (★ = nueva) |
|---|---|
| **V** (geometrías d ≥ 2) | RGG2 k12, RGG3 k12, ★RGG3 k8, ★RGG3 k20, RGG en S², RGG3 en caja, RGG2 en caja, toro cuadrado 141², T³ 27³ |
| **X** (árboles, cliques, expansores, ciclos abundantes no geométricos) | Árbol uniforme (Prüfer); árbol binario subdividido; cactus de triángulos; ★**árbol de cliques K6** (cada clique nueva comparte un nodo con una clique existente elegida al azar); ★**2-árbol aleatorio** (cada nodo nuevo se une a los dos extremos de una arista existente al azar); ★**grafo aleatorio con triángulos** (N/3 triángulos aleatorios más aristas aleatorias hasta k ≈ 12); ER k12; RR k12; ★**unión de cliques K12 enlazadas al azar** (cada clique con 2 aristas a cliques aleatorias) |
| **U** (1D: deben dar DOS_EXTREMOS, no CONEXO) | Anillo k12; caveman K8; C₁₀₀ × RR₂₀₀ (tubo); WS β = 0.001 |
| **E** (informativos) | Heisenberg H₃(ℤ₂₇) (se espera CONEXO; no euclídeo); retazos 2³ (se espera CONEXO; no es objetivo de esta observable) |

Para la regla RC-1, «no geométricos» = X ∪ U ∪ retazos. Heisenberg no entra en el veredicto: es un caso de estatus geométrico ambiguo.

## 6. Predicciones del cerebro (congeladas)

**L-CIC-0b:**
1. V: CONEXO en todos, incluidas las cajas, la esfera y las tres densidades.
2. X:
   - árbol subdividido, cactus y árbol de cliques: RAMIFICADO;
   - árbol uniforme: RAMIFICADO (confianza moderada; ver CIC-T2);
   - 2-árbol, grafo con triángulos, ER, RR y cliques enlazadas: SIN_VENTANA.
3. U: DOS_EXTREMOS en anillo, caveman y tubo; WS β = 0.001 DOS_EXTREMOS o CONEXO (los atajos pueden conectar los dos arcos).
4. E: Heisenberg CONEXO; retazos CONEXO.
5. Veredicto previsto: **REDUNDANCIA-VÁLIDA**.

**RC-1:**
- especificidad ≈ 1.0, porque cada adversario nuevo cae al menos por una condición;
- sensibilidad 0.8–0.95, porque T³ ya era INTERMEDIO en γ (0.69) y RGG3 está cerca del umbral (0.72).

Veredicto previsto: **PARCIAL o VÁLIDA**, limitado por la sensibilidad.

## 7. Lo que no se hace

- Ni dinámica ni energía.
- Ningún diagnóstico entra en la acción.
- No se reajusta ningún umbral existente (γ, K, Nivel II).
- La observable de anillos **no** se añade a RC-1 en esta fase: se valida por separado.
