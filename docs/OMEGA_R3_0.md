# Ω — R3-0: estructuras de orden, crecimiento y coalescencia (análisis sin dinámica)

- **Fecha:** 2026-10-07.
- **Rama:** `claude/omega-r3`. Base congelada: `claude/omega-e6-congelado`.
- **Mandato:** el Consejo autorizó R3-0 con esta condición: «determinar qué estructuras de orden permiten crecimiento y coalescencia sin especificar ninguna dimensión». Un orden que necesite estructura 2D/3D explícita para coalescer cuenta como negativo.
- **Única ejecución:** un cálculo estático (R3-0b, §5) que verifica las afirmaciones de nivel 1 del §3 con el juez RC-3. No hay dinámica.

---

## 0. Acta del Consejo sobre E6 v1.1 y su auditoría

### 0.1 Ratificado

| Decisión del Consejo | Estado |
|---|---|
| E6 v1: INVÁLIDA y congelada; no cuenta como evidencia | Adoptada (`claude/omega-e6-v1-congelado`) |
| E6 v1.1: válida como instrumento operacional, **con validación débil** | Adoptada |
| Etiqueta obligatoria para cualquier pase: «E6-compatible dentro del dominio de calibración; neutralidad dimensional no demostrada universalmente» | Adoptada literalmente |
| Prohibido cambiar umbrales, añadir controles tras ver resultados de R3, leer «pasa E6» como «dimensión emergente» o usar E6 para certificar 3D | Adoptado |
| E6 queda congelado; lo que falta es validación **fuera de muestra** | Adoptado |
| Cada candidato de orden declara antes de ejecutarse las propiedades que podrían seleccionar d (estados, aridad, interfaces, profundidad, niveles, anchura, composición, predecesores, dependencias) | Adoptado como **S0-ORD** (§4) |

### 0.2 Omisiones y correcciones (justificadas)

**Ω-a. El panel de órdenes no puede validar fuera de muestra a E6-S v1.1.**
- E6-S v1.1 está definido para reglas locales sobre **grafos sin etiquetas**: sus cantidades son grado, c₄, |S₂| y triángulos.
- Un orden parcial es otro dominio: relación dirigida y transitiva, con cantidades propias como anchura, altura y número de predecesores.
- Aplicar E6-S a órdenes exige un brazo nuevo, **E6-S-ORD**. Su calibración sería su *primera* validación, no una validación fuera de muestra de v1.1.

Decisión:
- (1) La deuda de v1.1 se salda con un **panel de grafos reservado**, prerregistrado a ciegas, antes del primer uso de E6-S sobre un candidato de tipo grafo. No se ejecuta antes porque hoy no hay ningún candidato de ese tipo, y el Consejo pidió no seguir desarrollando E6.
- (2) E6-S-ORD solo se diseña si algún orden sobrevive a R3-0. En ese caso su calibración se prerregistra **a ciegas**, lo que la haría una validación fuerte desde el principio.

**Ω-b. La lista S0-ORD del Consejo es correcta, pero omite cuál de sus elementos es decisivo.**
- Según §3 (R3-T1), en la única clase de orden que coalesce sin defectos, la dimensión del espacio generado **es** la anchura asintótica del orden generador.
- La anchura no es un elemento más de la lista: es el portador exacto de d en esa clase.

**Ω-c. La esperanza final del Consejo choca con un hecho estructural.** La esperanza era «un orden que pasa E6, crece, coalesce y luego produce una dimensión estable».
- Coalescer sin defectos exige confluencia (distributividad, §3).
- La confluencia fija d como el número de hilos independientes.
- Por tanto, para que el orden pase S2 (no fijar ese número), el número tiene que **autoseleccionarse**. Ese es el problema B, reformulado. R3 no esquiva B; lo hace exacto.

---

## 1. Qué se entiende por «estructura de orden» y por «coalescencia»

- **Orden generador J:** un conjunto parcialmente ordenado (causal: «x precede a y») que crece añadiendo elementos.
- **Espacio de cortes L(J):** el conjunto de *down-sets* de J (conjuntos cerrados hacia el pasado: los «presentes» posibles). Dos cortes son adyacentes si difieren en un elemento. La geometría candidata es el grafo de cubrimientos (Hasse) de L(J), sin coordenadas: es la estructura de «estados globales consistentes» de un proceso causal.
- **Coalescencia sin defectos:** dos dominios crecidos por separado se unen en un único dominio. En L(J) es la unión de down-sets, que siempre es un down-set. Unir sin reglas adicionales es exactamente la operación de supremo de un **retículo distributivo**.
- **Coalescencia con defectos:** cualquier unión que no sea la unión de conjuntos, es decir, supremos en retículos no distributivos o pegado arbitrario de dominios.
- **Relación con lo ya sabido:**
  - El pegado arbitrario de dominios ricos en cuadrados es la familia de los *retazos*, ya excluida por RC-3.
  - SQ (R1-1) formó dominios que no se unían.
  - En la Hasse de un retículo distributivo, cada «cuña» que sale de un corte hacia dos cortes superiores se cierra en un cuadrado (propiedad del diamante). Es la **conmutación local** de L-DIM-0 y de la hipótesis (ii) de ARQ-0 (rango de relaciones conmutativas).

---

## 2. Hechos externos usados (nivel 3; verificación bibliográfica en §7)

| Id | Enunciado | Fuente |
|---|---|---|
| T-Bir | Todo retículo distributivo finito es isomorfo al retículo de down-sets de su orden de elementos irreducibles por supremo | Birkhoff (1937) |
| T-Dil | La dimensión de orden del retículo de down-sets L(P) es igual a la anchura de P | Dilworth (1950) |
| T-BFR | Un retículo finito es plano si y solo si su dimensión de orden es ≤ 2 | Baker, Fishburn y Roberts (1971) |
| T-Whi | El retículo libre con 3 generadores es infinito; el modular libre con 4 generadores es infinito | Whitman (1941); Dedekind (1900) |
| T-KR | Casi todos los órdenes parciales de n elementos tienen 3 niveles | Kleitman y Rothschild (1975) |
| T-Post | El orden de percolación transitiva con p fijo tiene «postes» (elementos comparables con todos) con densidad positiva | Bollobás y Brightwell (1997); véase §7 sobre la atribución |
| T-CSG | El crecimiento secuencial clásico con covarianza discreta y causalidad de Bell está parametrizado por una sucesión t_n ≥ 0; la percolación transitiva es el caso t_n = tⁿ. Que algún modelo CSG produzca órdenes de tipo variedad es un **problema abierto** | Rideout y Sorkin (2000); revisión de Surya (2019) |
| T-Mat | Los estados globales consistentes de una computación distribuida forman un retículo distributivo | Mattern (1989) |

Todos se usan solo por su contenido combinatorio. No se usa mecánica cuántica: los modelos CSG citados son estocásticos clásicos.

---

## 3. Resultados

### R3-T1. En la coalescencia sin defectos, d es la anchura asintótica del generador

Nivel 3 (T-Bir, T-Dil) más un corolario de nivel 1.

- **Enunciado.** Sea J de anchura w. Por el teorema de Dilworth, J se parte en w cadenas C₁…C_w. Cada corte D ∈ L(J) queda determinado por el vector (|D ∩ C₁|, …, |D ∩ C_w|) ∈ ℕ^w, y dos cortes adyacentes difieren en +1 en una sola coordenada.
- **Corolario 1.** La Hasse de L(J) es un subgrafo inducido por cubrimientos de la rejilla ℤ^w. Su dimensión de crecimiento es ≤ w.
- **Corolario 2.** La igualdad no está garantizada. Lo que fija la dimensión es el número de cadenas que pueden **avanzar independientemente sin límite**:
  - Si J son w cadenas disjuntas, L(J) es la caja [0, n]^w, de dimensión w.
  - Si las cadenas se acoplan con desfase acotado ℓ (x_{i,k} < x_{j,k+ℓ} para todo i ≠ j), un corte cumple |c_i − c_j| ≤ ℓ. L(J) es un tubo alrededor de la diagonal, de dimensión 1.

Verificación de los corolarios: inmediata a partir de la definición de down-set.
- Si x_{j,k+ℓ} ∈ D, entonces x_{i,k} ∈ D para todo i.
- Por tanto c_i ≥ c_j − ℓ.

**Lectura.**
- La coalescencia sin defectos resuelve el problema A para cualquier anchura fija: la caja [0, n]^w es una geometría que RC-3 no excluye para w = 2, 3.
- Pero mueve d a un recuento: el número de hilos asintóticamente independientes del proceso causal.

### R3-T2. Las tres maneras naturales de fijar la anchura codifican, colapsan o explotan

Nivel 1 para las clases (a) y (b); nivel 3 o 4 para (c) y (d).

| Generador J | Qué fija la anchura asintótica | Resultado para L(J) | Clase E6 |
|---|---|---|---|
| (a) w cadenas independientes (proceso con w hilos) | Un parámetro w | Dimensión w | **S2 FAIL**: es el producto de w cadenas, que presenta ℤ^w (prerregistro E6 §3.4) |
| (b) Hilos acoplados por mensajes con desfase acotado | El acoplamiento | Dimensión 1 (tubo) | Pasa S2, pero da un degenerado 1D (excluido por X3) |
| (c) Orden aleatorio de anchura creciente (orden aleatorio de 2 dimensiones, percolación transitiva con p pequeño, órdenes típicos) | Nada: la anchura crece con n | Número de cortes superpolinómico en la altura: «explosión», tipo mundo pequeño (predicción, §5) | Pasa S2; degenerado |
| (d) Percolación transitiva con p fijo | Postes de densidad positiva (T-Post) | Asintóticamente una cadena de bloques finitos, de dimensión 1 | Pasa S2; degenerado 1D, posiblemente tras un régimen transitorio de explosión (riesgo W5) |

Además, el **orden 2-dimensional aleatorio** (intersección de 2 órdenes lineales) introduce el 2 en su construcción. Como generador candidato sería S2 FAIL; aquí solo se usa como objeto de panel.

### R3-T3. La coalescencia con defectos lleva a crecimiento exponencial o a codificación de 2D

Nivel 3.

- **Retículos no distributivos generados libremente:** son infinitos ya con 3 generadores (T-Whi). Los modulares típicos, como los retículos de subespacios, crecen exponencialmente. Clase hiperbólica o de mundo pequeño, excluida por X1.
- **Imponer planaridad para frenar ese crecimiento:** por T-BFR, eso equivale a dimensión de orden ≤ 2. Es una codificación de 2D (S2 FAIL).
- **Pegado arbitrario de dominios:** retazos, ya excluidos (RC-3).
- Esto enlaza con L-DIM-0: lo que mata el crecimiento exponencial es la conmutación local. En órdenes, eso es la distributividad.

### R3-T4. El crecimiento causal estocástico estándar ya está cerrado

Nivel 3.
- Los órdenes típicos son de 3 niveles (T-KR): diámetro acotado, X1.
- La percolación transitiva tiene postes (T-Post), así que es 1D asintóticamente.
- En la familia CSG general no se conoce ningún modelo de tipo variedad (T-CSG). Además, su espacio de parámetros (t_n) es infinito-dimensional, de modo que cualquier éxito tendría que pasar la prueba del dial E2 sobre una sucesión entera.

---

## 4. Auditoría S0-ORD de las clases de orden (adoptada del Consejo, ampliada)

| Propiedad (lista del Consejo) | Cómo puede codificar d | Clase (a) | Clase (b) | Clase (c) | Clase (d) | Retículos no distributivos |
|---|---|---|---|---|---|---|
| **Anchura (número de hilos)** | Por R3-T1 **es** d en la clase confluente | **fija = d** | acotada → 1 | crece | postes → 1 | — |
| Profundidad / niveles / altura | Con altura acotada → diámetro acotado (T-KR) | libre | libre | — | libre | — |
| Número de predecesores (aridad del supremo) | Un supremo binario o k-ario libre explota (T-Whi) | 1 | ≤ w | aleatorio | aleatorio | k |
| Reglas de composición (supremo) | Solo la unión (distributiva) evita defectos | unión | unión | unión | unión | supremo libre |
| Número de estados / interfaces | w interfaces con pares que conmutan = ℤ^w (S2) | w | — | — | — | — |
| Dependencias (estructura causal) | Acoplamiento acotado → 1D | ninguna | acotada | aleatoria | aleatoria | — |
| **Veredicto S0–S2** | | **FAIL (S2)** | pasa, pero 1D | pasa, pero explota | pasa, pero 1D | explota o FAIL (planaridad) |

---

## 5. R3-0b: verificación estática con el juez (prerregistrado aquí, antes del código)

**Objetivo.** Comprobar con RC-3, el instrumento congelado, que la dimensión de L(J) sigue a la anchura asintótica, y que los generadores sin anchura fijada dan degenerados. No se usa E6-S: es un brazo de grafos sin etiquetas y L(J) es aquí un objeto del panel, no una regla.

**Construcción:**
- Para cada J se enumeran todos los down-sets por búsqueda en anchura desde ∅ (añadiendo elementos minimales del complemento).
- La geometría es el grafo no dirigido de cubrimientos de L(J).
- Tamaños: se elige el menor n (número de elementos de J) tal que |L(J)| ≥ 10⁴ (tamaño N) y ≥ 8·10⁴ (tamaño 8N). El par se mide con `tools/rc3.py::measure`, y se calculan δ, X1–X4 y la validez W5 (E/N dentro del 25 %).
- Si |L(J)| crece tan deprisa que el par se sale de [N, 1.5 N] en cualquiera de los dos tamaños, se registra y se mide igual.
- RNG: `rng_from_key((MASTER_R3, id_familia, semilla, índice_de_tamaño))`, con `MASTER_R3 = 20261021` y semillas 0–2 para las familias aleatorias.

| Id | J | Parámetros |
|---|---|---|
| J1(w) | w cadenas disjuntas de longitud n | w = 1, 2, 3, 4 |
| J2(w) | w cadenas con x_{i,k} < x_{j,k+2} (desfase ℓ = 2) | w = 2, 3, 4 |
| J3 | Orden 2-dimensional aleatorio: n puntos uniformes en [0, 1]², orden producto (objeto de panel) | — |
| J4(p) | Percolación transitiva: i < j con probabilidad p para i < j, y cierre transitivo | p = 0.02, 0.1, 0.5 |

**Predicciones del cerebro y criterio de corroboración (congelados):**

| Id | Predicción |
|---|---|
| J1(1) | Camino: X3 |
| J1(2), J1(3) | No excluidos; δ dentro de 1/w ± 0.08 |
| J1(4) | X4-marginal (sesgo W1 de RC-3 para d ≥ 4) |
| J2(w) | X3 (δ > 0.75 o DOS_EXTREMOS) para todo w |
| J3 | X1 (δ < 0.125) o abstención con señales de explosión (|L| superpolinómico en n) |
| J4(0.5) | X3 |
| J4(0.1) | X3 |
| J4(0.02) | Incierto: puede parecer explosivo a este tamaño. Si sale no excluido, se inspecciona como residuo con W5 |

R3-T1/T2 quedan **corroborados** (nivel 1, para estas instancias) si se cumplen las tres condiciones:
- (i) J1(2) y J1(3) no están excluidos, con δ dentro de 1/w ± 0.08;
- (ii) todos los J2 están excluidos por X3;
- (iii) ningún par válido según W5 de J3 o J4 queda no excluido con δ ∈ [0.2, 0.6].

Si (iii) falla, hay un **residuo**: un generador sin anchura fijada que parece dar una dimensión finita mayor que 1. Se inspecciona antes de concluir nada. Sería el único resultado de R3-0 que abriría una vía.

---

## 6. Veredicto previsto y propuesta (condicionados a R3-0b)

Si R3-0b corrobora:

1. **R3 queda REDUCIDO, no cerrado.** Las estructuras de orden permiten crecimiento y coalescencia sin especificar d solo de una manera, la confluente (distributiva). En ella d es la anchura asintótica. Las tres formas naturales de fijar esa anchura la codifican (S2 FAIL), la colapsan a 1D o la hacen explotar.
2. **Síntesis con lo ya sabido (nivel 2).** En las tres vías del programa d es un **recuento de direcciones independientes que conmutan**:
   - grupos: rango de Bass–Guivarc'h;
   - redes rígidas: k = 2d (Ω3);
   - órdenes: anchura de Dilworth.

   El «?» de la jerarquía de ARQ-0 queda identificado: **rango de conmutación**. El problema B pasa a ser exactamente: ¿qué dinámica ciega a d autoselecciona el número de direcciones independientes?
3. **Advertencia del lema del dial.** Un número estacionario de hilos fijado por el equilibrio entre bifurcación y fusión es un cociente de tasas, es decir, un dial (E2), salvo que algún principio independiente fije las tasas *antes* de mirar resultados. Esto es nivel 2 y se registra como el obstáculo principal de cualquier R3-1.
4. **Propuesta al Consejo:**
   - (a) ratificar la reducción;
   - (b) decidir si se abre **W-0**: análisis sin dinámica de la autoselección de anchura en procesos de hilos que se bifurcan y se fusionan, con la prueba del dial obligatoria desde el diseño;
   - (c) alternativamente, cerrar el programa de selección con el resultado negativo estructurado: A es resoluble dado un recuento; B no tiene mecanismo conocido que no sea un dial.

## 7. Verificación bibliográfica

Agente Haiku con búsqueda web; revisión del cerebro.

| Id | Resultado del agente | Decisión del cerebro |
|---|---|---|
| T-Bir | VERIFICADO. Birkhoff, *Rings of sets*, Duke Math. J. 3 (1937) 443–454, https://doi.org/10.1215/S0012-7094-37-03409-9 | Se mantiene |
| T-Dil | VERIFICADO. Dilworth, Ann. of Math. 51 (1950) 161–166 (sin URL) | Se mantiene |
| T-BFR | VERIFICADO. Baker, Fishburn y Roberts, *Partial orders of dimension 2*, Networks 2 (1971) 11–28, https://doi.org/10.1002/net.3230020103 | Se mantiene |
| T-Whi | VERIFICADO (Whitman, Ann. of Math. 42 (1941) 325–329, https://doi.org/10.2307/1968774). Matiz de Dedekind: el modular libre con **3** generadores es **finito** (28 elementos, 30 con ⊥ y ⊤); con 4 o más es infinito | Lo escrito («con 4 generadores es infinito») era correcto. Se añade el matiz: la frontera 3/4 refuerza que la finitud modular es una excepción de pocos generadores |
| T-KR | VERIFICADO. Tres niveles con proporciones ≈ 1:2:1 | Se mantiene |
| T-Post | El agente afirma que Alon, Bollobás, Brightwell y Janson, *Linear extensions of a random partial order*, Ann. Appl. Probab. 4 (1994) 108–123, trata de extensiones lineales y no de postes. Propone Bollobás y Brightwell, *The structure of random graph orders*, SIAM J. Discrete Math. 10 (1997) 318–335 | **Atribución corregida** a Bollobás y Brightwell (1997) como fuente principal. Mi recuerdo es que el artículo de 1994 usa los postes en su análisis, pero no lo puedo confirmar; queda como atribución secundaria no verificada. El enunciado (postes con densidad positiva para p fijo) se mantiene como nivel 3 |
| T-CSG | Rideout y Sorkin, Phys. Rev. D 61 (2000) 024002, https://arxiv.org/abs/gr-qc/9904062: VERIFICADO. «t_n = tⁿ para la percolación transitiva»: no encontrado textualmente. «Ningún modelo CSG conocido es de tipo variedad»: Surya (Living Rev. Relativ. 22 (2019) 5, https://arxiv.org/abs/1903.11544) lo presenta como problema abierto, no como imposibilidad | **Enunciado debilitado** en §2 a «problema abierto». La correspondencia tⁿ (con t = p/(1 − p)) queda como recuerdo del cerebro, no verificado textualmente. Ninguna conclusión de R3-0 depende de ella: R3-T4 solo usa que el espacio de parámetros es una sucesión infinita (E2) |
| T-Mat | VERIFICADO. Mattern, *Virtual time and global states of distributed systems* (1989) (sin URL) | Se mantiene |

Efecto en las conclusiones: ninguno. R3-T1 descansa en T-Bir y T-Dil, que están verificados, más la prueba directa del §3. R3-T4 queda con un enunciado más débil (problema abierto), suficiente para el veredicto.
