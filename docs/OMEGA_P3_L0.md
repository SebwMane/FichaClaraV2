# Ω-P3 — L-P3-0: análisis previo de la regla de crecimiento (sin dinámica) y prerregistro de L-P3-0b

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-p3`.
- **Base congelada:** `claude/omega-d-congelado`.
- **Mandato del Consejo:** P3 avanza a L-P3-0, un **análisis de la regla, sin dinámica**. Hay que responder P3-1…P3-6 y la matriz de destinos. **Veto:** si P3 necesita introducir una propiedad equivalente a «vecindad 3D», queda rechazado antes de simular. El precedente de triangulaciones dinámicas **no** se usa para diseñar la regla.

## 0. Contraste del dictamen con lo que sabemos

| Punto del Consejo | Contraste | Decisión |
|---|---|---|
| Cierre de Ω-D como «resultado negativo limpio» | Coincide con L-ΩD-0. Matiz: lo que falló es la **premisa generativa**; la difusión sigue siendo un observable útil (Q ≈ d/2 en retículos; persistencia *post hoc*) | Se adopta |
| «Observable geométrico ≠ mecanismo geométrico» | Es la lección acumulada: D_s, ρ, Q y la curvatura **miden**; ninguno **genera** por sí solo | Principio de diseño |
| Quinto destino: **árbol / fractal** | Correcto y necesario. Un árbol aleatorio era ya control NO_LOCAL en C0 (sin ciclos). Los polímeros ramificados son el árbol crítico | Matriz de 7 clases |
| Geometría: «preguntamos si D_eff ≈ 3 y D_s ≈ 3» | **Contradice una decisión anterior del Consejo** (P1-D.2: «no fijar D = 3; la clase dimensional queda abierta»). Además, P1-D.3 mostró que D ≈ 3 se fabrica con WS | Se mantiene la clase abierta (1, 2, 3+); el juez es la batería P1-D.3 |
| «Si la regla conoce la vecindad geométrica que queremos, la geometría no emerge» | Correcto. Coincide con el veto | Se aplica a cada familia (§3) |
| Orden parcial mejor que total | Correcto. El crecimiento da un orden total de inserción, pero el orden con contenido es la **dependencia**: v ≻ u si v se enganchó a u (clausura transitiva), lo que forma un DAG | Se adopta (§2) |
| «No usar el precedente de triangulaciones para diseñar» | Correcto. Las referencias se están verificando aparte y solo se compararán al final | Se adopta |
| Lenguaje: «orden de precedencia combinatoria», no «causalidad física» | Correcto | Se adopta |

## 1. Pregunta de L-P3-0

¿La regla combinatoria propuesta contiene, antes de evolucionar, un mecanismo suficiente para evitar los destinos degenerados conocidos? Esos destinos son siete:

1. vacío / clique (0D);
2. cadena (1D);
3. **árbol / fractal**;
4. expansor;
5. mundo pequeño;
6. dimensión impuesta por parámetros;
7. geometría (el único destino aceptable).

Todo ello sin coordenadas, sin distancias geométricas y sin dimensión objetivo.

## 2. Definiciones comunes (respuestas P3-1, P3-2, P3-5, P3-6)

| Pregunta | Respuesta |
|---|---|
| **P3-1 (orden)** | Orden **parcial** de dependencia: v ≻ u si u es padre de v (clausura transitiva). La inserción solo es un índice de construcción, sin interpretación física |
| **P3-2 (información)** | El nodo nuevo solo conoce el **vecindario de radio ≤ 2 de sus candidatos a padre** y su posición en el orden (si son maximales). El planificador elige *dónde* ocurre el siguiente evento (uniforme sobre el conjunto activo); eso equivale a actualizaciones paralelas asíncronas y **no** es información que use la regla |
| **P3-5 (crecimiento)** | Crece \|V\| (un nodo por evento) y \|E\| (las aristas a los padres). Variantes con aristas entre nodos antiguos o con eliminación se declaran aparte, porque cambian la clase del modelo |
| **P3-6 (juez)** | La batería P1-D.3 (Niveles I y II por bolas; certificado como Nivel III) con su dominio declarado. Si la regla usa curvatura (familia D), el informe del certificado se da **con y sin** su componente de curvatura (`manifold_proxy_ok` la incluye), para que mecanismo y juez no coincidan |

## 3. Resultados analíticos

| Id | Enunciado | Nivel |
|---|---|---|
| **T1 (árbol recursivo de expansión)** | Si cada nodo nuevo se engancha al menos a un ancla elegida entre Θ(N) nodos sin correlación con la estructura (uniforme, preferencial…), las aristas de ancla forman un árbol recursivo aleatorio de expansión. Por tanto diámetro O(log N) con alta probabilidad: **mundo pequeño**, lo cual es incompatible con D finita | 1, dado el teorema externo de altura ≈ e·ln N del árbol recursivo (nivel 3, en verificación) |
| **T2 (balance de capacidad en la frontera)** | Con capacidad máxima de grado κ_c y m enlaces por nodo nuevo, la capacidad libre de la frontera cambia en κ_c − 2m por evento. Tres casos: **lineal**, con frontera Θ(N) que reproduce T1 (mundo pequeño); **constante**, con frontera finita que da un tubo **1D**; **negativa**, y el crecimiento se extingue. Una frontera ∝ N^{(d−1)/d} exige una ley de superficie, es decir, imponer d | 2 (campo medio) |
| **T3 (crecimiento por capas = ramificación)** | Si las capas por orden crecen como L_{t+1} = μ_eff L_t, hay tres casos: μ > 1 da crecimiento **exponencial** (hiperbólico, expansor); μ < 1 da extinción; μ = 1 crítico da un **árbol crítico, de tipo polímero ramificado** (fractal). Un crecimiento polinómico L_t ∝ t^{d−1} exige μ_eff = 1 + (d−1)/t: ajuste fino o realimentación | 2, más un hecho externo sobre árboles críticos (nivel 3, en verificación) |
| **T4 (ciclos)** | Con un padre por nodo el resultado es **exactamente un árbol** (rango cíclico 0). Con m padres, el rango cíclico es ≈ (m − 1)N. Si los padres son siempre mutuamente adyacentes, se pegan cliques K_{m+1} (símplices), con riesgo de 0D o de dimensión fijada por m | 1 |
| **T5 (planitud = ajuste fino)** | Gauss–Bonnet combinatorio en superficies trianguladas: Σ_v (6 − deg v) = 6χ, de modo que la curvatura media nula ⇔ grado medio 6. **El k = 6 que el Consejo prohibió es exactamente la condición de planitud 2D.** Una regla local genérica da un grado medio ≠ 6, es decir, curvatura positiva (cierre, esfera finita, cliques) o negativa (hiperbólica, exponencial, como T3 con μ > 1). La geometría euclídea es el **punto crítico** entre ambas | 1 para Gauss–Bonnet (en verificación); 2 para la inferencia sobre el crecimiento |

**Síntesis (nivel 2).** Una regla de crecimiento local genérica cae en mundo pequeño (T1), en tubo 1D (T2), en crecimiento hiperbólico o en cierre y clique (T3, T5), o en árbol o polímero ramificado (T3 crítico, T4 con m = 1). La única forma de evitarlo **sin un entero impuesto** es un **principio de selección de la planitud**: una realimentación local que lleve la curvatura a cero sin fijar qué dimensión es plana. Es exactamente la pregunta 6 del Consejo («¿qué selecciona la escala?») reformulada como «¿qué selecciona la planitud?».

## 4. Familias de reglas evaluadas

| Familia | Regla | Destino según §3 | ¿Usa información prohibida? | Veredicto L-P3-0 |
|---|---|---|---|---|
| **A** | Ancla uniforme; el nodo nuevo se enlaza al ancla y a m vecinos del ancla | Mundo pequeño (T1) | No | **MUERE** |
| **B** | Frontera por saturación (capacidad κ_c, m enlaces) | 1D o mundo pequeño (T2) | No | **MUERE** (salvo κ_c = 2m: 1D) |
| **C** | Capas por orden con coalescencia local (padres: un maximal y maximales a distancia ≤ 2, cada uno con probabilidad p) | Hiperbólico, extinción o polímero ramificado según μ_eff(p) (T3); planitud solo con p ajustado | No | **MUERE como genérica**. La dimensión, si apareciera, dependería de p, lo que violaría la meseta C1-D |
| **D** | **Crecimiento homeostático de curvatura:** como C, pero el nodo nuevo elige sus padres entre los candidatos locales de modo que la **curvatura de Ollivier** de las aristas afectadas se acerque a 0 | No queda excluida por T1–T5: es precisamente un principio de selección de planitud. La planitud de Ollivier no distingue dimensiones (Z^d tiene κ = 0 para todo d, y el ciclo también) | Usa curvatura local de radio 2: información local permitida. **No** usa D, coordenadas ni k objetivo | **SOBREVIVE CONDICIONADA** |

**Condiciones de la familia D (riesgos declarados):**
1. **«Planitud como ley» es una elección de principio**, análoga a las ecuaciones de vacío (Ricci = 0). No es una propiedad 3D, así que el veto no se activa, pero **el Consejo debe ratificarla** como admisible.
2. **La dimensión sigue abierta**, porque planitud no implica d. Puede salir 1D, ya que el ciclo es plano. Hay que medir la meseta en parámetros (C1-D).
3. **A escala 1, la curvatura de Ollivier favorece estructuras sin triángulos.** Los retículos tienen κ = 0, pero una geometría desordenada con triángulos (RGG) tiene κ > 0 local. D produciría estructuras cristalinas, no RGG. Eso es aceptable (son geometrías), pero hay que medirlo.
4. **Solapamiento juez–mecanismo.** El certificado usa curvatura, de ahí la regla P3-6.

**Respuesta a P3-3 (ciclos):** en las familias C y D, los ciclos los crea la coalescencia (m ≥ 2 padres) sin exigir un número fijo de triángulos.

**Respuesta a P3-4 (número de enlaces):** en D no hay m fijo. El número de padres lo decide la realimentación de curvatura. Solo hay una tolerancia de curvatura, que es continua; la robustez frente a ella es la prueba de meseta.

## 5. Veredicto de L-P3-0

| Familia | Veredicto |
|---|---|
| A, B, C | **Rechazadas analíticamente**. Caen en destinos degenerados sin necesidad de simular. Es la «victoria metodológica» que pedía el Consejo |
| **D** | **Única familia viable**, condicionada a que el Consejo ratifique la planitud como ley local |

Antes de definir su regla exacta (L-P3-1) se comprueba su **premisa**: que la curvatura local distinga los destinos degenerados y el pegado incoherente. Esa prueba es una medición sobre grafos fijos, sin dinámica, igual que L-ΩD-0.

## 6. Prerregistro de L-P3-0b: ¿ve la curvatura de Ollivier los destinos y los retazos?

- **Panel** (N ≈ 1000; clave PCG64 (20261011, familia, semilla); semillas 0–2 en las familias aleatorias):
  - mismas familias que L-ΩD-0, salvo la unión de cliques, sustituida por caveman K8;
  - más el **árbol aleatorio** (N = 1000), por el destino árbol;
  - informativo: los 24 finales C0 de N = 729.
- **Cálculo:** `omega.curvature.ollivier.ollivier_edge` con perezosidad 0.5 y coste de saltos sobre la componente gigante, en hasta 1500 aristas muestreadas por grafo.
  - Por grafo: mediana de κ, f_neg = fracción con κ < −0.1 y f_pos = fracción con κ > 0.1.

| Criterio | Pregunta | Condición |
|---|---|---|
| **KD1** | ¿Los destinos tienen signo? | Caveman K8: mediana de κ ≥ +0.2. ER, RR y árbol: mediana ≤ −0.2. T³ y toro cuadrado: \|mediana\| ≤ 0.02 |
| **KD2** | ¿Ve el pegado incoherente? | En las 3 semillas de retazos-2³ y de retazos-3³, f_neg ≥ f_neg(mediana RGG3) + 0.05 |
| **KD3** (informativo, condiciona el diseño) | ¿Es plana a escala 1 una geometría desordenada? | Mediana de κ de RGG3 y RGG2. Si es ≥ +0.1, la familia D a escala 1 favorece cristales; se anota para L-P3-1 |

| Resultado | Consecuencia |
|---|---|
| KD1 y KD2 pasan | La premisa de D es viable. Se propone L-P3-1 (regla exacta y análisis de sus puntos fijos) al Consejo, junto con la ratificación de «planitud como ley» |
| **KD2 falla** | La curvatura local, como la difusión, no ve la incoherencia a escala 1. D queda **debilitada**: necesitaría curvatura gruesa. El Consejo decide |
| **KD1 falla** | **D muere**, y P3 queda sin familia viable |

**Predicciones del cerebro (congeladas):**
1. KD1 pasa.
2. KD2 pasa: las aristas entre bloques se comportan como aristas de árbol, con κ negativo.
3. KD3: RGG3 y RGG2 tienen κ positivo, ≈ +0.2…+0.4, por los triángulos.
4. Los finales C0 tienen κ positivo, por ser complejos de cliques.
5. WS β = 0.01: la mayoría de aristas planas o positivas, más una cola negativa pequeña de atajos.

## 7. Lo que no se hace

- No se simula ninguna regla de crecimiento.
- No se define todavía la regla exacta de D: eso es L-P3-1, tras la ratificación.
- No se usan las referencias externas para diseñar.
