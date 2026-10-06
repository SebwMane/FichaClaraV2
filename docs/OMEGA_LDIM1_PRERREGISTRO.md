# Ω — L-DIM-1: prueba de invariancia de la selección dimensional (sin dinámica) — Prerregistro

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-ldim1`.
- **Base congelada:** `claude/omega-ldim-congelado`.
- **Mandato del Consejo:**
  - L-DIM-0 aceptada como sesión conceptual;
  - D5 queda como conjetura de nivel 4, no validada;
  - L-DIM-1 autorizada, sin dinámica, como **test de invariancia**: «¿existe una regla que pueda seleccionar d sin saber qué significa d?».
  - **Veto:** si la meseta solo aparece tras introducir, aunque sea indirectamente, grado, densidad, escala o algo equivalente a d, es un fracaso, aunque salga d = 3.
- **Estado:** este documento se commitea antes de escribir el código y antes de evaluar ningún funcional sobre ningún grafo.

## 0. Contraste del dictamen

| Punto | Contraste | Decisión |
|---|---|---|
| D2: «la curvatura no negativa evita el crecimiento exponencial» es demasiado fuerte | Correcto. Mi propia revisión de referencias mostró que hace falta CDE′(0, n), una hipótesis que Ω no tiene | Se rebaja: «las restricciones de curvatura no negativa **pueden** controlar el crecimiento **bajo hipótesis adicionales**» |
| D4: «ya resuelto» | Correcto: reconocer una degeneración ≠ impedirla dinámicamente | Se corrige |
| D3 como cuello de botella; la robustez máxima lleva a expansores | Coincide con L-DIM-0. Se formaliza en L-DIM-T2 | — |
| D5: separar la parte A (homogeneidad → entero) de la parte B (competencia → selección) | Correcto. La parte B es la que L-DIM-1 intenta falsar | — |
| Un mínimo interior sería un resultado no trivial, porque los extremos son fáciles | Correcto, y el análisis §2 lo **cuantifica**: el extremo compacto (clique) es el ganador natural de la competencia estática | — |
| Diseño ciego, θ único para todas las d, discretizaciones múltiples, prueba de no codificación, nulo | Se adoptan y se operacionalizan (§3–§5) | — |
| Panel con T¹…T⁵ y RGG; la propuesta de «pares con X(G_d) ≈ X(G_d′)» | La prueba de permutación se aplica sobre las variables que usa F (§4, NC-4) | — |

## 1. La familia de candidatos (forma fijada desde el principio de D5, antes de mirar grafos)

El esquema D5 dice: coste de expansión, que crece con d, más coste de fragilidad, que decrece con d. Su traducción mínima al vocabulario de bolas, sin coordenadas, sin d y sin grado, es:

- r* = max(1, ⌊r_max(W)/2⌋), escala **intrínseca** de cada grafo (mitad de su ventana intermedia; ninguna escala impuesta);
- m(r) = media sobre fuentes muestreadas de |B_r|, y |S_r| = m(r) − m(r−1) con m(0) = 1;
- **expansión:** φ = |S_{r*+1}| / m(r*), es decir, frontera/volumen. Si la bola ya saturó y m(r*+1) no existe, φ = 0 (no queda frontera);
- **fragilidad:** ψ = 1 / |S_{r*}|, el inverso del tamaño del separador que rodea la bola;
- **candidato:** F_{p,θ}(G) = φ + θ · ψ^p, con p ∈ {1/2, 1, 2}.

El exponente p es la única libertad de forma; se prerregistran los tres valores. θ recorre 25 valores log-espaciados en [10⁻³, 10³].

## 2. Análisis previo (antes de evaluar)

**L-DIM-T1 (el exponente codifica la dimensión; nivel 1 bajo escalado).**
- Supuestos: |B_r| ≈ c r^d y |S_r| ≈ c′ r^{d−1}.
- Entonces r · F = d + θ c′^{−p} r^{1−p(d−1)} + O(1/r).
- Cuando r → ∞:
  - las dimensiones con p(d−1) < 1 tienen coste divergente;
  - con p(d−1) = 1 el coste es d + θ′;
  - con p(d−1) > 1 el coste tiende a d.
- Por tanto el ganador asintótico es

      d*(p) = el menor entero con p(d−1) > 1  →  p = 1/2: d* = 4;  p = 1: d* = 3 (si θ′ > 1; si no, 2);  p = 2: d* = 2.

**La dimensión seleccionada la fija el exponente de la forma, no la estructura.** Ninguna fórmula contiene d, pero el cociente de exponentes la codifica: es **codificación estructural**.

**Matiz.** La elección p = 1 es la «lineal» (aditiva, como resistencias en paralelo a través de la cáscara). Su umbral, Σ r^{1−d} < ∞ ⇔ d ≥ 3, coincide con el de **transitoriedad del paseo aleatorio** (Pólya, nivel 3). Si existiera una razón **independiente** para que la robustez relevante de Ω fuera de tipo lineal o de escape, la selección de 3 tendría contenido. Sin esa razón, es una elección de forma.

**L-DIM-T2 (compactificación; nivel 1).**
- En una clique, o en cualquier estructura que sature antes de r* + 1, φ = 0 y ψ = 1/(N − 1).
- Por tanto F_{p,θ}(clique) ≈ θ (N − 1)^{−p}, **menor que el de cualquier geometría extendida** para todo θ dentro de la rejilla, ya que φ > 0 en geometrías.
- **La competencia estática expansión–fragilidad la gana el extremo compacto.** Un mínimo interior exige además una exigencia de **extensión**, que sería una escala o tamaño impuesto y, por tanto, vetada.

**Consecuencia prevista:** todo candidato de esta familia fallará el criterio de degenerados (L4) frente a la clique. Esto es lo que el Consejo anticipó al decir que los extremos son fáciles. **Se prerregistra como predicción**, y la prueba numérica sirve para comprobar que el análisis no se equivoca a tamaño finito.

## 3. Protocolo ciego

1. **Generación.** Cada grafo recibe un **ID opaco** (hash). El código de **características** solo recibe la adyacencia, sin nombre, d, grado ni familia, y escribe `features.jsonl` (ID → m(r), r*, φ, ψ).
2. **La verdad** (ID → familia, d, discretización, N) se escribe en un archivo aparte, `truth.jsonl`, y solo la lee el **evaluador**.
3. **El evaluador** une ambos archivos después de calcular F, con el mismo θ para todas las d.

## 4. Prueba de no codificación (NC)

| Prueba | Pregunta | Fallo |
|---|---|---|
| **NC-1 (desacoplo del grado)** | El panel tiene, para d = 2, 3, 4, discretizaciones con grado medio ≈ 8 y ≈ 16 (RGG_d) además del retículo (k = 2d). ¿Es d* el mismo en cada estrato y entre estratos? | Si d* cambia con el estrato, F lee el grado o la densidad |
| **NC-2 (equivalencia de proxy)** | Spearman entre F y el grado medio sobre las geometrías, para cada (p, θ) de la meseta | \|ρ_s\| ≥ 0.9: F es un lector de grado |
| **NC-3 (ceguera)** | ¿Se respeta el §3? | Uso de cualquier etiqueta antes de calcular F |
| **NC-4 (permutación sobre las variables de F)** | Para cada par (G_d, G_d′) con d ≠ d′, ¿hay pares con (φ, ψ) a menos del 10 % relativo? | Si los hay, (φ, ψ) no contienen información suficiente para seleccionar d, y la preferencia de F no es dimensional. Se informa el número de pares |
| **NC-5 (codificación por exponente)** | ¿d*(p) sigue a L-DIM-T1 (4, 3, 2 para p = 1/2, 1, 2)? | Si sí, la selección la fija la forma: **veto estructural**, salvo que el Consejo ratifique p por una razón independiente |

## 5. Protocolo nulo

- **Nulo A (funcionales aleatorios).** M = 500 funcionales F_ω = a·z_φ + b·z_ψ, con (a, b) uniformes en la circunferencia unidad y z = log de la característica estandarizado sobre el panel. Se mide la **fracción de nulos que cumple L1–L4**. Si un candidato pasa y la tasa nula es > 5 %, el éxito es **no específico**.
- **Nulo B (etiquetas permutadas).** Se permutan al azar las etiquetas d entre las geometrías (200 veces) y se recalculan L1–L3. La tasa de aprobación tiene que ser ≈ 0; si no, el criterio está sesgado.
- **Nulo C (panel solo de degenerados).** Si sin geometrías también aparece un «mínimo interior» (con la etiqueta d heredada), el procedimiento fabrica mínimos.

## 6. Panel (disperso; clave (20261014, familia, semilla); semillas 0–2 en las familias aleatorias)

| Grupo | Familias |
|---|---|
| **Geometrías por d** | d = 1: anillo k2 (N = 20 000) y anillo k8 · d = 2: toro cuadrado 141², RGG2 k8, RGG2 k16 · d = 3: T³ 27³, RGG3 k8, RGG3 k16 · d = 4: T⁴ 12⁴ (20 736), RGG4 k8, RGG4 k16 · d = 5: T⁵ 7⁵ (16 807) |
| **Escala** (L7) | Retículos de 4N aprox.: anillo 80 000, cuadrado 283², T³ 43³, T⁴ 17⁴; RGG3 k16 con N = 8·10⁴ |
| **Degenerados** | Clique K200; árbol uniforme; árbol binario subdividido; cactus; 2-árbol; RR k12; ER k12; caveman K8; WS β = 0.01; retazos 2³ |

## 7. Criterios de decisión (por candidato p)

| Id | Criterio | Condición |
|---|---|---|
| L1 | Existencia | Hay θ con d* ∈ {2, 3, 4} (d* = argmin_d de la media de F sobre las discretizaciones de esa d) |
| L2 / L5 / L6 | Invariancia | En esos θ, d* calculado por separado en cada estrato (retículo, RGG k8, RGG k16) coincide (NC-1) |
| L3 | Meseta | El conjunto de θ con el mismo d* interior en todos los estratos abarca ≥ 1 década |
| L7 | Escala | El mismo d* con los retículos de N y de ≈ 4N, en los θ de la meseta |
| L4 | Degenerados | En los θ de la meseta, **todos** los degenerados tienen F > F(d*) |
| L8 | No codificación | NC-2 sin bandera **y** NC-5 sin veto |
| L9 | Ceguera | NC-3 cumplida |
| — | Nulos | Tasa nula A ≤ 5 %; nulos B y C ≈ 0 |

**Veredicto:**
- **SELECCIÓN-VIABLE:** algún p pasa L1–L9 y los nulos.
- **SELECCIÓN-NO-VIABLE:** ningún p pasa. Se informa qué criterio falla para cada p.
- **Regla de muerte para D5:** si todos los p fallan, D5 queda **refutado como mecanismo estático**. Además, si NC-5 confirma L-DIM-T1, la selección es codificación por exponente aunque algún p pasara L1–L7.

## 8. Predicciones del cerebro (congeladas)

1. **L4 falla para los tres p** por la clique (L-DIM-T2).
2. Ignorando degenerados, d*(p) a tamaño finito tiende a 4, 3 y 2 para p = 1/2, 1, 2 (L-DIM-T1), con corrimientos por tamaño finito en d = 4 y 5 (ventanas cortas).
3. NC-5 confirma la codificación por exponente.
4. Nulo A: una fracción sustancial de funcionales aleatorios produce algún «mínimo interior», porque con dos características monótonas en d es fácil.
5. **Veredicto previsto: SELECCIÓN-NO-VIABLE; D5 refutado como mecanismo estático.** Además se identifica qué tendría que aportar un mecanismo no estático:
   - una exigencia de extensión que no sea escala impuesta;
   - una justificación independiente del exponente.

## 9. Lo que no se hace

- Ni dinámica ni ajuste de θ o de p tras ver los datos.
- Ningún uso de d, grado o familia en el cálculo de F.
- No se reabre la instrumentación.
