# Ω — L-DIM-0: principios de selección de dimensionalidad (sesión conceptual; sin código, sin dinámica, sin parámetros ajustados)

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-ldim`.
- **Base congelada:** `claude/omega-ciclos-congelado`.
- **Mandato del Consejo:**
  - cerrar la instrumentación (no hay RC-2 ni más filtros);
  - no diseñar dinámica todavía;
  - responder D1–D6: ¿qué principio puede seleccionar espontáneamente una dimensión finita sin imponer ninguna dimensión ni favorecer a priori un entero?
- **Niveles de afirmación:**
  1. demostrado;
  2. explicación;
  3. teorema externo (verificación en §8);
  4. extrapolación o conjetura.

## 0. Contraste del dictamen

| Afirmación | Contraste | Decisión |
|---|---|---|
| Cerrar la instrumentación; RC-1 como rechazador conservador | Coincide | Instrumentación **cerrada** |
| El 2-árbol como clase nueva: «topología no degenerada sin geometría macroscópica» | Correcto | Se añade al mapa negativo |
| Tres niveles: conectividad → organización local/mesoscópica → organización macroscópica | Útil, pero **choca con los nombres de los Niveles I/II/III del instrumento** | Se adopta como **estratos** S1/S2/S3 para no confundir |
| «La pregunta correcta: ¿qué permite una jerarquía macroscópica de escalas extensas?»; crecimiento polinómico frente a exponencial | Bien orientada, pero **incompleta: el crecimiento subexponencial no basta**. El árbol uniforme (polímero ramificado) crece como r² y no es geometría. Las cadenas crecen como r¹. Es la fase de polímeros ramificados de las triangulaciones euclídeas | Se reformula en §1 como tres presiones, no una |
| «Patrón común: ninguna regla controla la expansión a gran escala» | **Parcialmente falso.** Algunos fallos sí son de expansión (P3-A/T1, WS, expansores, 2-árbol), pero otros **controlan la expansión y aun así fallan**: C0 desde U (cadenas 1D, crecimiento polinómico), C0-R (D_L 2.7–4.1, finito, no variedad) y el árbol uniforme (r²) | El control de la expansión es **necesario, no suficiente** |
| D es emergente, no parámetro; no seleccionar 3 | Coincide | Se mantiene |
| Formular «una dimensión finita sin favorecer un entero» antes de «¿por qué 3?» | Correcto | Estructura de D1–D6 |

## 1. Reformulación: tres presiones, no una

Los destinos degenerados observados se agrupan en tres fallos independientes:

| Fallo | Destinos | Lo que falta |
|---|---|---|
| **F-exp** (expansión) | Expansores, mundo pequeño, 2-árbol, P3-A | Control de la expansión: crecimiento subexponencial |
| **F-ram** (ramificación) | Árboles (también críticos), cactus, árbol de cliques | Un extremo a toda escala: ciclos que rodean |
| **F-baja** (colapso dimensional) | Cadenas 1D, tubos, C0-U, cliques (0D) | Algo que impida reducir d hasta 1 o hasta 0 |

Una geometría extendida de dimensión finita ≥ 2 exige evitar **las tres a la vez**. El dictamen del Consejo apunta solo a F-exp.

## 2. Respuestas D1–D6

### D1. ¿Puede una dinámica local producir crecimiento polinómico |B_r| ~ r^D?

**Sí, pero hay que ver qué lo produce.**
- **(Nivel 3, Gromov.)** En estructuras homogéneas, el crecimiento polinómico equivale a una estructura de grupo virtualmente nilpotente. El caso canónico, Z^d, se obtiene de un grupo libre (crecimiento exponencial) añadiendo **relaciones de conmutación**: cuadrados xy = yx entre generadores.
- **(Nivel 2.)** El ingrediente que mata la expansión exponencial es la **conmutatividad local**: «dos pasos en direcciones distintas se cierran en un cuadrado». **No basta con tener ciclos.** El cactus y el 2-árbol tienen muchos, pero no son ciclos de conmutación.

**Advertencia (nivel 1–2).** En un retículo conmutativo homogéneo, grado = 2d. Seleccionar d equivale a seleccionar el **número de direcciones independientes que conmutan**. Fijar el grado (k = 6) fija d = 3, que es la prohibición original del Consejo. Una dinámica tendría que **descubrir** cuántas direcciones independientes hay, no recibirlo.

### D2. ¿Qué evita el crecimiento exponencial?

**Curvatura no negativa (nivel 3).** En geometría riemanniana, Ric ≥ 0 implica |B_r| ≤ C r^n (Bishop–Gromov). En grafos existen análogos discretos (curvatura-dimensión de Bakry–Émery; §8). La curvatura negativa sistemática corresponde a la ramificación, la hiperbolicidad y la expansión.

**Pero la curvatura positiva colapsa (nivel 3).** Con la curvatura de Ollivier, κ ≥ κ₀ > 0 en todas partes implica **diámetro ≤ 2/κ₀**, un Bonnet–Myers discreto: se pierde la escala y la estructura tiende a clique o esfera finita.

**Síntesis (nivel 2).** El control de la expansión sin colapso exige **curvatura gruesa ≈ 0**, es decir, planitud. Es exactamente el punto crítico de T5 (Gauss–Bonnet).

Lo que añade este análisis frente a L-P3-0b: allí se midió la **mediana** de κ ≈ 0, y los árboles pasaban. Lo relevante es la **ausencia de curvatura negativa sistemática a toda escala**, porque los puntos de ramificación son negativos. Esta lectura está informada por los datos de L-P3-0b y se marca como tal.

### D3. ¿Qué evita el colapso a D = 1?

En 1D, un separador de tamaño acotado parte el espacio en dos extremos (anillo de dos componentes). En d ≥ 2 los separadores deben crecer como r^{d−1} y el anillo es conexo (CIC-T1).

Una presión anti-1D natural es la **robustez**: penalizar los separadores pequeños o favorecer caminos alternativos (redundancia relacional, hipótesis B del Consejo). **Pero la robustez máxima la tienen los expansores** (constante de Cheeger máxima). La presión anti-1D empuja hacia F-exp.

**Resultado estructural (nivel 2):** las presiones contra F-exp y contra F-baja son **opuestas**. Toda propuesta debe equilibrarlas.

### D4. ¿Qué evita el colapso a cliques?

Ya está resuelto en el programa:
- curvatura positiva → diámetro acotado (Bonnet–Myers);
- saturación y presupuesto de grado al estilo de C0 → no hay densificación.

La clique es la versión 0D de F-baja.

### D5. ¿Puede el entero D ser consecuencia de estabilidad, sin imponerlo?

Es la respuesta central de esta sesión: **sí es posible en principio, por un mecanismo de dos pasos.**

1. **La integralidad la da la homogeneidad (nivel 3).** En estructuras transitivas (homogéneas) coherentes (crecimiento polinómico), el grado de crecimiento es **necesariamente un entero** (Gromov; Trofimov; Bass–Guivarc'h; CH-T6). No hace falta imponer «entero»: lo impone la homogeneidad.
   - Corolario (nivel 2): las estructuras **inhomogéneas** pueden tener dimensión fraccionaria. Es lo que se observó en C0-R (D_s ≈ 2.2–2.5) y en los polímeros ramificados.
   - Por tanto, el Nivel I del instrumento (ρ, homogeneización) **no es un detalle**: es la condición que hace posible una dimensión entera.
2. **La selección entre enteros la daría un equilibrio de presiones, con mesetas (nivel 2–4).**
   - Sea F(d) un coste por nodo que combina una presión anti-expansión, creciente con d (por ejemplo, superficie/volumen, grado o redundancia que hay que mantener), y una presión anti-baja, decreciente con d (fragilidad de los separadores ~ r^{−(d−1)}).
   - En el continuo, el mínimo d* = f(B/A) variaría **continuamente** con los parámetros. Es el patrón de «dimensión ajustada» que el Consejo prohíbe.
   - **Pero si las estructuras homogéneas realizables solo tienen d entero** (paso 1), el mínimo es arg min_{d ∈ ℕ} F(d). Ese mínimo es **constante a trozos** en B/A: **mesetas**.
   - En cada intervalo abierto de parámetros se selecciona el mismo entero. **Eso satisface el criterio de meseta C1-D**: la dimensión no está ajustada, porque es robusta en regiones abiertas.

   **Implicación de diseño (nivel 4):** un mecanismo candidato necesitaría
   - (i) generar **homogeneidad** estadística;
   - (ii) **planitud** (curvatura gruesa ≈ 0, D2);
   - (iii) un **equilibrio** entre expansión y robustez cuya forma no contenga d.

   Entonces la dimensión aparecería como **la meseta en la que cae el sistema**.

**Riesgo honesto.**
- Con dinámicas desordenadas (no transitivas), la integralidad no está garantizada. Un mecanismo podría producir d fraccionario estable, como C0, y entonces el paso 1 falla.
- **Qué meseta** se selecciona depende de los parámetros. Esto traslada la pregunta «¿por qué 3?» a «¿por qué los parámetros caen en la meseta 3?».

### D6. ¿Hay una razón física para D* = 3?

**Desde Ω actual: no determinable.** Ω no tiene campos, materia ni ondas, y las razones conocidas para d = 3 dependen de ellos (nivel 3, §8):

| Argumento | Qué dice |
|---|---|
| Estabilidad orbital y atómica (Ehrenfest; Tangherlini) | Bajo un potencial ∝ r^{−(d−2)}, solo en d = 3 hay órbitas y estados ligados estables |
| Nudos | Solo en ℝ³ existen nudos no triviales de círculos |
| Huygens | La propagación nítida de ondas solo vale en dimensiones impares ≥ 3 |

**Posición del cerebro.** Estas razones son **selección a posteriori** (dado el contenido físico), no un mecanismo de emergencia. Dentro de Ω, la única vía interna sería D5. Si un funcional produce mesetas, se puede **calcular qué anchura tiene la meseta d = 3** en el espacio de parámetros. Si fuera la más ancha o la única estable frente a fluctuaciones, sería un resultado. Hoy es pregunta abierta, y las tres respuestas del Consejo (sí / no / no determinable) siguen vigentes.

## 3. Lo que esta sesión establece

1. **El control de la expansión es necesario pero no suficiente.** Hay tres presiones (F-exp, F-ram, F-baja), y las de F-exp y F-baja se oponen.
2. **La conmutatividad local** (cuadrados de conmutación), no la mera abundancia de ciclos, es lo que convierte la expansión exponencial en polinómica en estructuras homogéneas (nivel 3 + 2).
3. **Planitud gruesa = control de la expansión sin colapso** (D2). Encaja con T5 y con K.
4. **Una dimensión entera robusta puede emerger sin imponerla:** la homogeneidad da integralidad (Gromov–Trofimov–Bass–Guivarc'h) y un equilibrio continuo de presiones entre enteros da mesetas. Es el primer esquema del programa que **satisface por construcción el criterio de meseta C1-D** en lugar de violarlo. Nivel 4: es un esquema, no un mecanismo.
5. **«¿Por qué 3?» no es determinable desde Ω actual.** Las razones físicas conocidas son a posteriori.

## 4. Propuesta: L-DIM-1 (sin dinámica, analítico y numérico sobre familias fijas)

**Objetivo:** comprobar si el esquema D5 puede tener contenido. Para cualquier funcional candidato que el Consejo proponga, se evaluaría su **paisaje dimensional** sobre la familia {T^d_N : d = 1, …, 5} con N fijo, junto a sus competidores degenerados (clique, expansor, árbol, 2-árbol, retazos). Las preguntas son:
- ¿el mínimo es interior (2 ≤ d ≤ 4)?
- ¿el entero óptimo es constante en regiones abiertas de parámetros (mesetas)?
- ¿los competidores degenerados quedan por encima?

Es la generalización del método C0-T1/L-2 al eje dimensional. **Prohibición explícita:** la forma del funcional no puede contener d, ni k = 2d, ni un valor objetivo de solapamiento o de 4-ciclos (por L-ΩD-T2, CH-T1 y la Fase 2 §1.3).

**Lo que hace falta del Consejo antes de L-DIM-1:**
1. Ratificar o rechazar el esquema D5 (homogeneidad, planitud y equilibrio expansión–robustez) como marco.
2. Proponer, sin mirar datos, **qué forma** podría tener la presión de robustez que no favorezca expansores. Es el punto más débil del esquema (D3).

## 5. Lo que no se hace

Ni código, ni dinámica, ni nuevos filtros. La instrumentación queda cerrada.
