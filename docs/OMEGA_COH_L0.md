# Ω — L-COHERENCIA-0: formalización de la coherencia multiescala, análisis previo y prerregistro de la prueba contra falsos positivos

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-coherencia`.
- **Base congelada:** `claude/omega-p3-congelado`.
- **Mandato del Consejo:**
  - (a) consolidar el mapa negativo;
  - (b) integrar Ollivier como diagnóstico complementario, con una auditoría sobre el panel P1-D.3;
  - D′ no aprobada;
  - nueva fase **L-COHERENCIA-0**, conceptual y sin dinámica: ¿puede definirse la coherencia topológica multiescala sin introducir geometría por la puerta trasera?
- **Estado:** este documento se commitea antes de escribir el código y antes de evaluar ningún grafo con la nueva observable.

## 0. Contraste del dictamen con lo que sabemos

| Punto | Contraste | Decisión |
|---|---|---|
| Redacción de la conjetura | Correcto. Lo que demuestran los datos es más limitado: **«las familias de reglas locales que hemos definido y analizado no han producido un mecanismo que seleccione espontáneamente una geometría extendida de dimensión finita»**. La versión fuerte (imposibilidad) sigue siendo una conjetura de nivel 4 | Se congela esta redacción |
| Tabla histórica: «C1 / variantes — controles y atajos» | **Inexacto.** C1 fue la funcional de 4-ciclos, rechazada analíticamente porque convierte la dimensión en parámetro (T³ con q = 4 frente a Q₆ con q = 5). Los atajos y los controles son adversarios del instrumento, no un mecanismo | Corregido en el mapa (a) |
| «P3-D: árbol y C0 pasan» | Precisión: el árbol tiene curvatura **mixta** (mediana −0.17/0; 51 % negativa, 37 % positiva) y C0 es casi plano (+0.03). Ambos quedan en la banda «plana» | Se adopta: K ≈ 0 es condición insuficiente |
| «No hay dimensión objetivo» | Coincide con P1-D.2 y P1-D.3: la clase queda abierta. Se añade aquí un resultado que lo matiza (CH-T6) | Se adopta |
| «Ingrediente faltante: coherencia de vecindades; B_r(u) ≈ B_r(v) para u, v próximos, sin igualdad» | Es formalizable (§1). **Pero el análisis previo la debilita como mecanismo:** como objetivo favorece 1D (CH-T5), y su valor de prefactor depende de d (CH-T1). Sirve como **restricción de clase**, no como selector de dimensión | §2 |
| «Quizá los pares no bastan; relaciones T_ijk» | Coherente con C0-T1/T2 (las acciones en (k, c) no seleccionan). Pero la observable de §1 sí es computable sobre pares, porque usa bolas de radio r | Se registra como pregunta abierta |
| Predicción: «la dinámica debería curar un retazo» | Coincide con R6 (decisión post-P1-D.2) | Pendiente de un mecanismo |
| La física estadística duda de que la causalidad sea el ingrediente | Correcto: tras L-P3-0 no hay razón independiente para el orden | Se adopta |

## 1. Formalización (D1 de la coherencia)

Sea G un grafo no dirigido (componente gigante), con distancia de saltos d y bola B_r(x) = {y : d(x, y) ≤ r}.

**Solapamiento de cartas a lo largo de una arista e = (u, v):**

    J_r(e) = |B_r(u) ∩ B_r(v)| / |B_r(u) ∪ B_r(v)|,     I_r(e) = 1 − J_r(e)   (incompatibilidad).

**Ley de coherencia multiescala.** Se define el exponente de decaimiento por arista

    γ_e = − d ln I_r(e) / d ln r,   en la ventana intermedia r ∈ [2, r_w],

donde r_w es el mayor r con mediana de |B_r| ≤ N/4.

- Una estructura es **coherente** si sus vecindades vecinas se vuelven compatibles al crecer la escala con la ley de superficie/volumen, I_r ∝ 1/r, es decir γ ≈ 1.
- Es **degenerada** si I_r = 0 (bolas idénticas: nodos gemelos, cliques).
- Es **incoherente** si I_r no decae (γ ≈ 0).

**Qué no contiene la definición:** no usa coordenadas, ni distancias geométricas, ni un número de vecinos objetivo. La escala r no se fija: la ley se mide sobre toda la ventana.

## 2. Análisis previo (antes de medir nada con esta observable)

| Id | Enunciado | Nivel |
|---|---|---|
| **CH-T1** | En Z^d, B_r(u) y B_r(v) para vecinos difieren en una «cáscara», de modo que I_r(e) ≈ c_d / r con c_d ∝ d: **γ = 1 para todo d**. La ley no depende de d, pero el **prefactor** sí: un valor objetivo de I_r a escala fija impondría d. Es el mismo patrón que L-ΩD-T2 | 1–2 (superficie/volumen) |
| **CH-T2** | En grafos **no amenables** (árboles regulares, expansores, mundo pequeño a gran escala), \|∂B_r\|/\|B_r\| está acotado inferiormente (constante de Cheeger), así que I_r no tiende a 0: **γ → 0** | 2–3 (Følner y amenabilidad; referencias en verificación) |
| **CH-T3** | En una clique y entre nodos gemelos, B_r(u) = B_r(v) para todo r ≥ 1: I_r = 0, caso **degenerado**. Caveman: la mayoría de las aristas internas son entre gemelos | 1 |
| **CH-T4** | En un retazo, las aristas entre bloques tienen I_r grande hasta la escala de bloque. Más allá, el grafo de bloques es un expansor (reemparejado aleatorio) y el crecimiento deja de ser polinómico, lo que da γ → 0 también en las aristas internas a gran escala | 2 |
| **CH-T5 (advertencia central)** | **Como objetivo, la coherencia favorece 1D.** A presupuesto fijo, I_r ≈ c_d/r es mínimo para d = 1 (la cadena tiene la menor incompatibilidad posible no degenerada). Una dinámica que maximice la coherencia colapsaría a cadenas. La coherencia es una **restricción de clase** (separa amenable/polinómico de no amenable y degenerado); **no selecciona la dimensión** | 2 |
| **CH-T6 (resultado que matiza «no hay dimensión objetivo»)** | En grafos **transitivos por vértices** con crecimiento polinómico (los coherentes homogéneos), el grado de crecimiento es **un entero** (Gromov; Trofimov; fórmula de Bass–Guivarc'h). Homogeneidad más coherencia implican dimensión entera, sin fijar cuál. Pero ese entero **no** tiene por qué ser euclídeo: el grupo de Heisenberg discreto tiene crecimiento 4 con dimensión topológica 3 | 3 (en verificación) |

**Respuestas a M1–M7 del Consejo:**

| Pregunta | Respuesta |
|---|---|
| **M1** (¿1D satisface la coherencia?) | **Sí, perfectamente** (CH-T5). La coherencia no puede ser el único principio. Hace falta otro que impida el colapso a 1D sin fijar d: la pregunta abierta |
| **M2** (¿la clique?) | Solo como **degenerada** (I = 0). Se excluye con la condición de no degeneración I_1 > 0 |
| **M3** (¿expansor?) | **No** (CH-T2) |
| **M4** (¿retazo?) | **No** a escala mayor que el bloque (CH-T4). Es el test central y se mide en §3 |
| **M5** (¿escala r?) | No se fija. La ley se mide en toda la ventana; la ventana la determina N |
| **M6** (¿depende de D?) | No explícitamente. El prefactor sí (CH-T1): queda prohibido usar valores objetivo de I |
| **M7** (¿W → 0 o W → 1?) | W → 0 da grafo desconectado o sin ventana; W → 1 da clique degenerada. Ambos excluidos |

**Conclusión previa.** La coherencia multiescala está bien definida, no introduce geometría por la puerta trasera y separa de forma limpia las clases amenable, no amenable y degenerada. **No resuelve la selección de dimensión** (CH-T5). Su valor potencial es doble:
1. **instrumento:** podría corregir el falso negativo con borde de P1-D.3, porque mide compatibilidad entre vecinos, no homogeneidad global;
2. **restricción de diseño** para cualquier mecanismo futuro, combinada con un principio aún desconocido contra el colapso a 1D.

## 3. Prerregistro de L-COH-0b: prueba contra falsos positivos (medición sobre grafos fijos, sin dinámica)

### 3.1 Observable (nuevo módulo `omega/diagnostics/coherence.py`)

- 400 aristas muestreadas uniformemente de la componente gigante, con un RNG derivado de la clave.
- Para cada arista, BFS desde u y desde v sobre el grafo disperso.
- Para r = 1 … r_cap (el primer r con mediana de |B_r(u)| ≥ N/2, máximo 200): |B_r(u)|, |B_r(v)|, intersección y unión, e I_r(e).
- **Ventana:** r_w es el mayor r con mediana de |B_r(u)| ≤ N/4.
- γ_e es la pendiente por mínimos cuadrados de −ln I_r(e) frente a ln r en r ∈ [2, r_w]. Exige r_w ≥ 4; si no, `SIN_VENTANA`.
- **Arista degenerada:** I_r(e) = 0 para algún r ∈ [1, r_w].
- **Por grafo:** mediana de γ (aristas no degeneradas); f_deg; f_low = fracción de aristas no degeneradas con γ_e < 0.5; mediana de I_1.

**Estado por grafo (umbrales fijados a priori por CH-T1/T2: γ = 1 ideal, γ = 0 no amenable):**

| Estado | Condición |
|---|---|
| `SIN_VENTANA` | r_w < 4 |
| `DEGENERADO` | f_deg > 0.5 |
| `COHERENTE` | mediana de γ ≥ 0.7 y f_low ≤ 0.10 y f_deg ≤ 0.10 |
| `INCOHERENTE` | mediana de γ < 0.5 o f_low > 0.25 |
| `INTERMEDIO` | cualquier otro caso |

### 3.2 Panel de falsos positivos (N ≈ 2·10⁴, disperso; clave (20261012, familia, semilla); semillas 0–2 en las aleatorias, 1 en las deterministas)

| Grupo | Familias | Función |
|---|---|---|
| **V** (verdaderos: geometrías) | RGG2 k12, RGG3 k12, RGG en S² k12, T³ 27³, toro cuadrado 141², **RGG3 en caja**, **RGG2 en caja**, RGG3 con gradiente | Sensibilidad; incluye las dos familias con borde que P1-D.3 rechazaba |
| **F** (falsos positivos diseñados: no geométricos) | Árbol aleatorio (Prüfer); **árbol binario subdividido** (cada arista de un árbol binario completo de profundidad 11, subdividida en 8 tramos; N ≈ 1.8·10⁴); **cactus de triángulos** (árbol de triángulos de Husimi: cada nodo en 2 triángulos); ER k12; RR k12; retazos 2³, 4³ y 8³; WS β = 0.01 y 0.005; RGG3 con 1 % de atajos; **C_100 × RR(200, k = 4)** (producto de ciclo y expansor) | Especificidad. Cada uno ataca un punto: ciclos locales con crecimiento exponencial (cactus), cruce de escala (árbol subdividido, WS), pegado incoherente (retazos), producto parcialmente coherente |
| **E** (casos de estatus especial; se informan, no deciden) | Anillo k12 (1D coherente); caveman K8 (degenerado); retículo de cliques 3D K6 L15; **Heisenberg discreto H₃(ℤ₂₇)** (Cayley con x^{±1}, y^{±1}; N = 19 683; amenable, crecimiento 4, no euclídeo); WS β = 0.001 y RGG3 con 0.1 % de atajos (cruces lentos: límite conocido) | Alcance: lo que la coherencia **no** puede distinguir por construcción |

### 3.3 Criterios

| Medida | Definición |
|---|---|
| **Sensibilidad** | Fracción de grafos V con estado `COHERENTE` |
| **Especificidad** | Fracción de grafos F **no** `COHERENTE` |

| Veredicto | Condición |
|---|---|
| **COHERENCIA-VÁLIDA** | sensibilidad ≥ 0.9 y especificidad ≥ 0.9 |
| **COHERENCIA-INVÁLIDA** | sensibilidad < 0.7 o especificidad < 0.7 |
| **COHERENCIA-PARCIAL** | cualquier otro caso |

Además se informa por separado:
- si las dos familias con borde salen `COHERENTE`, lo que corregiría el falso negativo de P1-D.3;
- el estado de cada caso E.

### 3.4 Predicciones del cerebro (congeladas)

1. **V:** `COHERENTE`, con γ ≈ 0.8–1.0, **incluidas las cajas**: el borde no rompe la compatibilidad entre vecinos.
2. **F:** árboles, cactus, ER, RR y el producto, `INCOHERENTE` (γ < 0.3). Retazos 4³ y 8³ `INCOHERENTE`; retazos 2³ quizá `INTERMEDIO` (bloques grandes). WS β = 0.01 y 0.005, `INCOHERENTE`. Atajos al 1 %, `INTERMEDIO` o `INCOHERENTE`.
3. **E:**
   - anillo `COHERENTE` (γ ≈ 1);
   - caveman `DEGENERADO`;
   - **Heisenberg `COHERENTE`** (la coherencia no implica geometría euclídea);
   - WS β = 0.001 y atajos al 0.1 %, `COHERENTE` en la ventana (límite conocido).
4. **Veredicto previsto:** COHERENCIA-VÁLIDA, con el alcance «separa amenable/polinómico de no amenable; no distingue 1D, nilpotente ni cruces lentos».

## 4. Lo que no se hace

- No hay mecanismo, acción ni dinámica.
- La coherencia **no** entra en ninguna energía: CH-T5 muestra que colapsaría a 1D.
- No se modifica el certificado. La observable, si es válida, se propondría como diagnóstico adicional del preflight (punto B-bis).
