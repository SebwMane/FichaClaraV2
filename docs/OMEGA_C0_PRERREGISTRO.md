# Ω-C0 — Hipótesis de competencia relacional: prerregistro

- **Fecha:** 2026-10-05.
- **Estado:** PRERREGISTRADO. Este documento se commitea antes de escribir el código de simulación de C0 y antes de cualquier corrida.
- **Rama de trabajo:** `claude/omega-c0-competencia`. Sale de `ec448e4`.
- **Línea anterior congelada:**
  - Ω-1.1 más el bloque L quedan congelados en la rama `claude/omega-1.1-congelado`, en el commit `ec448e4`.
  - La rama `claude/hopeful-galileo-88u1l1` también quedó en `ec448e4`.
  - Nada de `omega/` anterior a este documento se modifica. C0 vive en `omega/c0/` y en `tools/c0_*.py`.
- **Mandato:** decisión del Consejo post-bloque L.
  - Cerrar S0/Ω-B como mecanismo principal.
  - Crear Ω-C0.
  - Aplicar la puerta C0-L1 → C0-L2 → C0-L3 → C0-L4 antes de cualquier simulación grande.
  - Aplicar el filtro red-team R1–R8.
  - No usar N grande.
  - Prohibido introducir coordenadas, distancia euclídea, retícula, D=3, k=6 o RGG como objetivo, embedding externo o mecánica cuántica.
  - D_eff, D_s, etc. se miden **solo** si C0-L4 encuentra una fase intermedia no trivial y local, y entonces solo como descripción.

## 0. Funcional Ω-C0 (fijada)

Estado: W simétrica, diagonal 0, pesos en [0,1].

- Fuerza (grado ponderado): k_i = Σ_j W_ij.
- Codegrado ponderado: C = W² (número ponderado de vecinos comunes).

```
S_C0[W] = − Σ_{i<j} W_ij · ψ(C_ij)  +  κ Σ_i k_i²
ψ(c)    = a + b·log(1+c) − λ·c
```

| Término | Ingrediente del Consejo | Papel |
|---|---|---|
| a | Conectividad | Incentivo básico de enlace |
| b·log(1+c) | Cierre saturado | f' > 0 y f'' < 0: rendimiento decreciente del vecino común |
| −λ·c | Redundancia | Un enlace entre nodos ya muy relacionados aporta poco o resta. Σ_{i<j} W_ij C_ij = ½ tr W³, así que es una penalización lineal de triángulos |
| κ Σ k_i² | Competencia y concentración | Penaliza grados extremos. **No fija ninguna media**: la media emerge |

- **Gauge:** b = 1, porque la escala global de S es irrelevante a Θ = 0.
- **Parámetros libres:** a, λ, κ. Las escalas derivadas, que emergen y no se imponen, son:
  - c* = argmax ψ = max(0, 1/λ − 1);
  - ψ* = ψ(c*);
  - k* = ψ*/(4κ).
- **Gradiente** (por arista, para W simétrica): G_ij = −ψ(C_ij) − (W M + M W)_ij + 2κ(k_i + k_j), con M = W ∘ ψ'(C).
- **Dinámica:**
  - descenso de gradiente proyectado en la caja [0,1] (sin restricción de masa; la densidad emerge);
  - paso adaptativo monótono: S no puede subir; si sube, dt se divide entre 2; si el paso se acepta, dt × 1.1;
  - convergencia: max|ΔW| < 1e-10 durante 50 pasos;
  - máximo de 20 000 pasos;
  - RNG solo con PCG64 explícito.

## 1. Análisis previo (cerebro de pruebas)

Lo escribo antes de construir nada, siguiendo el principio del red-team: intentar destruir la hipótesis antes de construirla.

### C0-T1 / C0-T2. Cota inferior exacta y estados que la alcanzan

**Nivel 1, demostración elemental.**

Para todo W en [0,1], con M = Σ_{i<j} W_ij:

```
S_C0[W] ≥ −ψ*·M + κ(2M)²/N ≥ −N ψ*² / (16 κ)
```

La primera desigualdad usa ψ ≤ ψ*; la segunda, Cauchy–Schwarz sobre los k_i y minimizar en M.

**Igualdad si y solo si** las dos condiciones siguientes se cumplen a la vez:
- todas las aristas con W_ij > 0 tienen codegrado ponderado exactamente c*;
- todas las fuerzas valen exactamente k*.

Consecuencias:

- **(T1, binario)** Los mínimos globales binarios son exactamente los grafos k*-regulares en los que cada arista está en c* triángulos (grafos «edge-regular» (N, k*, c*)).
  - Las uniones de cliques K_{k*+1} son mínimos globales **solo si** k* = c* + 1.
  - Fuera de esa línea, las cliques binarias pierden estrictamente. La dicotomía clique/vacío del bloque L se rompe en el sentido literal.
- **(T2, degeneración no geométrica)** Con pesos, la cota la alcanzan (o casi, salvo divisibilidad de N) estados densos y no locales:
  - c* = 0: la bipartita completa K_{N/2,N/2} con peso uniforme 2k*/N la alcanza exactamente. También la alcanza cualquier grafo k*-regular sin triángulos: T³ si k* = 6, pero igual un expansor bipartito aleatorio.
  - c* > 0: uniones de «cliques difusas» (bloques de tamaño s con peso k*/(s−1)) con (s−2)k*²/(s−1)² = c*.
  - **Predicción:** C0 no tiene preferencia energética por la geometría. En el mejor caso (c* = 0, k* = 6), T³ **empata** con K_{N/2,N/2} y con los expansores.
- **(T3, observación general)** Toda acción que sea suma de funciones de bolas de radio r es idéntica:
  - para todas las retículas periódicas suficientemente grandes;
  - para uniones disjuntas de toros pequeños;
  - para cualquier grafo con la misma distribución de bolas.

  Por eso una acción **puramente local** no puede preferir una geometría grande y conexa frente a sus versiones fragmentadas. En el mejor caso selecciona una *estructura local*.

**Verificación numérica rápida** (N = 216; solo matemática, no es dato de C0-L*):
- el gradiente pasa el chequeo por diferencias finitas;
- en c* = 0 y k* = 6, S(T³) = S(K_{108,108} óptima) = cota = −324.

### Predicciones del cerebro (registradas antes de correr)

| Prueba | Predicción |
|---|---|
| C0-L1 | Las cliques binarias **no** son extremales salvo cerca de k* = c* + 1 (celdas con c* ≥ k* − 1). C0 sobrevive a L1. |
| C0-L2 | **NO ROTA** (o EMPATE). La dominancia trivial se desplaza de «clique/vacío» a «bloques densos difusos y bipartitos». Ninguna referencia dispersa gana estrictamente. |
| C0-L3 | T³ es KKT en una región con c* = 0 y k* ≈ 6. Los vértices binarios con gradiente estrictamente firmado son mínimos locales estrictos, así que el paisaje es vidrioso. RGG3: probablemente no es KKT, porque su codegrado es heterogéneo. |
| C0-L4 | **MUERTE-A o MUERTE-B.** Desde condiciones genéricas, cabe esperar: bipartito denso o disperso sin triángulos y no local con c* = 0; cliques difusas o solapadas con c* > 0. |

**Lectura si se cumple:** el problema no es solo la recompensa lineal por triángulos. Ninguna funcional de estadísticas de radio 1 (k, c) distingue la geometría de los grafos algebraicos o aleatorios con las mismas estadísticas locales. La siguiente hipótesis tendría que codificar estructura de segundo orden o un mecanismo no energético (dinámico o entrópico). Eso queda para la propuesta C1, que **no** se ejecuta sin ratificación.

## 2. Malla (fijada)

| Eje | Valores |
|---|---|
| c* (fija λ) | {0, 1, 2, 4, 8} → λ = {2, 1/2, 1/3, 1/5, 1/9}. Con c* = 0 se usa λ = 2: todo triángulo resta |
| k* (fija κ = ψ*/(4k*)) | {4, 6, 8, 12, 16} |
| a | {0.25, 1, 4} |

- Total: **75 celdas**, con N = 216.
- k* = 6 y c* = 0 se incluyen como **celda favorable a T³ declarada**. Es la prueba adversarial al nulo, no un objetivo.

## 3. Referencias (fijadas)

Todas usan N = 216. Las semillas de grafos aleatorios son 0–4, con PCG64.

| Grupo | Referencias |
|---|---|
| **Triviales** | vacío (S = 0); uniforme; unión de cliques binarias K_s (todo divisor s de 216); unión de cliques difusas (s divisor, peso uniforme); bipartita completa K_{108,108} |
| **No triviales dispersas** | T³ 6³; RGG3 con k̄ ∈ {6, 8, 12, 16}; k-regular aleatorio con k ∈ {4, 6, 8, 12, 16}; ER con k̄ ∈ {4, 6, 8, 12, 16}; árbol aleatorio; anillo con k = 6; anillo de Watts–Strogatz con β = 0.1; retícula triangular 2D en toro 12×18 |

- Cada referencia se evalúa **con la amplitud óptima** (W → t·A, con t ∈ (0,1] en una rejilla de 2000 valores más un refinamiento de Brent).
- La cota LB = −Nψ*²/(16κ) se reporta en cada celda.

## 4. Pruebas y criterios (congelados)

### C0-L1 — ¿siguen siendo extremales las cliques?

- Por celda: ¿la mejor unión de cliques **binarias** es la mejor referencia (mínimo S entre todas)?
- **MUERTE-L1** (no se simula nada más) si eso ocurre en ≥ 50% de las celdas.

### C0-L2 — ¿se rompe la dominancia trivial?

- Por celda, «rota» si la mejor referencia no trivial cumple S < (mejor trivial) − 0.01·|LB|.
- «Empate» si queda a no más de 0.01·|LB| de la mejor trivial sin mejorarla.
- **Global:**
  - **ROTA** si ≥ 10% de las celdas están rotas.
  - **EMPATE** si no se cumple ROTA y ≥ 10% de las celdas están rotas o empatadas.
  - **NO ROTA** en otro caso.
- Es informativa: no detiene la puerta. Si sale NO ROTA, cualquier tercera fase de L4 será dinámica, no energética.

### C0-L3 — KKT de T³ y RGG3

- Con gradiente G, las condiciones KKT en la caja son:
  - G ≤ tol en aristas con W = 1;
  - G ≥ −tol en aristas con W = 0;
  - |G| ≤ tol en el interior;
  - tol = 1e-9·max(1, |G|max).
- Se evalúan sobre el binario con amplitud 1 y sobre la amplitud óptima.
- Se reportan las celdas donde T³ es KKT y donde lo es alguna de las RGG3 (k̄ ∈ {6, 8, 12, 16}, semillas 0–4).
- Si se cumplen con signo estricto, es un mínimo local estricto.
- **INDICIO NEGATIVO** si ni T³ ni RGG3 son KKT en ninguna celda.

### C0-L4 — ¿aparece un tercer atractor?

**Inicios genéricos.** No geométricos, 3 por celda, con 3 semillas cada uno:
- **U:** uniforme k*/(N−1)·(1 + 0.1ξ), ξ ~ N(0,1), recortado a [0,1].
- **E:** ER binario con k̄ = k*.
- **R:** pesos iid U(0,1)·2k*/(N−1).

Son 75 × 3 × 3 = **675 corridas**. Además, como inicio geométrico **solo descriptivo** (no vota): T³ y RGG3 (k̄ = 12) más ruido iid de amplitud 1e-2, una semilla cada uno (150 corridas).

**Clasificación del estado final** (orden de decisión). El soporte fuerte es A = {W > 0.1·max W}; la gigante es la componente mayor de A.

1. **VACÍO:** max W ≤ 1e-6.
2. **DENSO-TRIVIAL:** se cumple cualquiera de estas:
   - `classify_structure` (bloque L) ∈ {uniforme, clique_única, multi_clique, cliques_solapadas};
   - la densidad de A sobre sus nodos no aislados es ≥ 0.25.
3. **FRAGMENTADO:** gigante < 0.5·N.
4. **DISPERSO-LOCAL:** se cumplen las dos condiciones siguientes. Si no, **DISPERSO-NO-LOCAL**.
   - **H_null ≥ 1.15:** distancia media de saltos en la gigante, dividida entre la media sobre 5 recableados que preservan grados (10·|E| intercambios dobles, PCG64; solo cuentan los pares conectados).
   - **Ciclos cortos ≥ 0.5:** fracción de aristas de la gigante que están en un triángulo o en un 4-ciclo.

**Validación preregistrada del criterio de localidad** (antes de cualquier dinámica):

| Grupo | Grafos | Debe dar |
|---|---|---|
| Positivos | T³ 6³; RGG3 k̄ = 12; retícula triangular 12×18; anillo k = 6 | LOCAL |
| Negativos | ER k̄ = 6; 6-regular aleatorio; K_{108,108}; unión de K8; árbol aleatorio | ≠ LOCAL |

Si la validación falla, se registra una enmienda **antes** de la dinámica (precedente L-A5).

**Decisión global** (solo inicios genéricos):

| Resultado | Condición | Consecuencia |
|---|---|---|
| **MUERTE-A** (condición del Consejo) | Todas las corridas son VACÍO o DENSO-TRIVIAL | C0 se descarta |
| **MUERTE-B** | Hay corridas FRAGMENTADO o DISPERSO, pero ninguna celda tiene DISPERSO-LOCAL en ≥ 2/3 semillas con el mismo inicio | La tercera fase existe pero no es local; se documenta. C0 se descarta como mecanismo de P1 |
| **CONTINÚA** | Alguna celda tiene DISPERSO-LOCAL en ≥ 2/3 semillas con el mismo inicio | Pasa al filtro red-team (abajo) |

**Filtro red-team** (se aplica solo a las celdas de CONTINÚA):
- **R3:** sin hubs (k_max/k̄ ≤ 3 en A).
- **R6:** misma celda con N = 343 y N = 125, ≥ 2/3 semillas DISPERSO-LOCAL en ambos.
- **R7:** equivarianza.
- **R8:** DISPERSO-LOCAL desde ≥ 2 de los 3 inicios genéricos.

Si lo supera, es **CANDIDATO-C0**. Solo entonces se aplica la batería geométrica descriptiva:
- D_eff (saltos);
- D_s;
- homogeneidad;
- isotropía;
- H_null.

Los umbrales del certificado no se modifican.

### Red-team R1–R8 sobre toda la malla (informe obligatorio)

| Test | Qué se reporta |
|---|---|
| R1 | Región de vacío: celdas donde algún inicio genérico termina en VACÍO |
| R2 | Región de cliques: celdas con clases de cliques |
| R3 | Hubs: k_max/k̄ de A en los finales no vacíos |
| R4 | Atajos: distribución de H_null en los finales dispersos |
| R5 | Árbol y D: D_eff del árbol aleatorio y de los controles. D ≈ 3 no basta: `RULE_D3_NEVER_SUFFICIENT` |
| R6 | Tamaño: solo para candidatos (arriba) |
| R7 | Reetiquetado. Se ejecuta siempre: 10 corridas al azar con la permutación aplicada al inicio; el final permutado debe coincidir con max\|Δ\| ≤ 1e-8 |
| R8 | Inicialización: la tabla de clase por inicio |

## 5. Lo que este prerregistro NO cubre

- Θ > 0, sea en C0 o en S0/Ω-B (O-05 queda como rama secundaria).
- N grande.
- Funcionales con estructura de segundo orden (C1).

Cualquiera de ellas requiere un documento propio.

## 6. Enmiendas registradas antes de las corridas oficiales

Todas se registraron después de las corridas *smoke* de código (3 celdas de paisaje y 10 dinámicas de 2000 pasos) y antes de cualquier corrida oficial de C0-L1…L4.

### C0-A1 — R7 (reetiquetado)

- **Problema:** la equivarianza exacta está garantizada por la matemática, porque solo hay operaciones matriciales. Sin embargo, en los valles planos degenerados que predice C0-T2, el orden de suma de BLAS da diferencias numéricas de entre 1e-8 y 4e-5 en las corridas convergidas.
- **Nuevo criterio de R7 (código):** max|Δ| ≤ 1e-8 a horizonte corto (60 pasos), sobre las 10 corridas sorteadas.
- **Nuevo criterio de R7 (resultado):** en las corridas finales que convergen, misma clase y |ΔS|/|S| ≤ 1e-6. max|ΔW| se reporta, pero no vota.

### C0-A2 — Precisión sobre C0-L3

- La predicción «mínimo local estricto» de §1 era incorrecta en la celda favorable (c* = 0, k* = 6, a = 1).
- Ahí T³ es un mínimo global (alcanza la cota), con G = 0 exactamente en sus aristas: −ψ* + 2κ·2k* = 0.
- Es KKT **no estricto**: hay direcciones planas, que son la degeneración de C0-T2.
- El criterio de L3 (KKT sí/no) no cambia. Solo se corrige la lectura.

### C0-A3 — Estado `stalled` (operativo)

- Una corrida termina `stalled` cuando dt < 1e-14 porque el ruido de redondeo hace que cualquier propuesta suba S.
- Se reporta `kkt_residual` (residuo del gradiente proyectado) y la corrida se clasifica igual que las demás.
- La regla de aceptación S' ≤ S no cambia.

### C0-A4 — Operativo

- Cada proceso usa un solo hilo BLAS (`OMP_NUM_THREADS` = 1).
- No afecta a ningún criterio.

### C0-A5 — Protocolo de R6 y diagnóstico D-1

Se registró después de ver el resultado de C0-L4 (CONTINÚA) y antes de correr R6 y D-1.

**R6 (sin cambio de criterio; se fija lo que §4 no detallaba)**
- Se aplica a las celdas candidatas que superan R3 y R8. Son 30 celdas.
- Por celda se repiten los inicios que votaron LOCAL (≥ 2/3 semillas), con las semillas 0, 1 y 2 y N ∈ {125, 343}.
- Mismos (c*, k*, a), max_steps, tolerancias e inicios definidos de forma intensiva (k*/(N−1), etc.).
- La celda supera R6 si **algún** inicio que votó obtiene DISPERSO-LOCAL en ≥ 2/3 semillas **en ambos** tamaños.
- La semilla RNG incluye N.

**D-1 (diagnóstico añadido; NO vota en la decisión preregistrada)**
- 157 de los 253 finales DISPERSO-LOCAL de inicios genéricos terminaron en `max_steps`, es decir, sin converger.
- Pregunta: ¿esos estados son estacionarios o son transitorios de engrosamiento (coarsening) hacia DENSO-TRIVIAL?
- Método: se continúa cada uno desde su W final hasta 100 000 pasos en total y se reclasifica.
- Se informa la fracción que sigue siendo DISPERSO-LOCAL, la que converge y la evolución de S/LB.
- Su lectura entra en el informe con el nivel de afirmación que corresponda.

### C0-A6 — Detalle de la batería descriptiva y escalado con el tamaño

Se registró antes de terminar R6 y antes de aplicar la batería a cualquier estado de C0. La herramienta es `tools/c0_battery.py`.

**Batería (solo describe; no hay umbrales nuevos):** sobre el soporte fuerte A y su componente gigante se miden:
- D_eff en saltos;
- D_s con el paseo lazy sobre A binario;
- homogeneidad;
- isotropía con k = round(D_eff);
- salto medio y diámetro;
- descriptores de cliques: clique máxima, cliques maximales ≥ 4 por nodo y fracción de nodos que están en alguna.

**Añadido: dimensión de escalado D_L.** Se define por L(N) ∝ N^{1/D_L}, con el salto medio en N = 125, 216 y 343 para cada par (celda, inicio). Con N ≤ 343, D_eff suele quedar en `no_window` incluso para T³ y RGG3, así que D_L es el estimador más estable a este tamaño.

**Calibración en controles** (solo grafos de referencia; ningún estado de C0):

| Control | D_L |
|---|---|
| T³ | 2.86 |
| RGG3 k12 | 3.09 |
| anillo k6 | 1.03 |
| caveman K8 | 1.08 |
| árbol aleatorio | 1.37 |

**Regla de lectura** (descriptiva, no veredicto):
- D_L < 1.5 → la fase es esencialmente unidimensional (cadena o árbol).
- 1.5–2.5 → de tipo 2D.
- \> 2.5 → de dimensión mayor.

D_L ≈ 3 no basta por sí solo (`RULE_D3_NEVER_SUFFICIENT`).

### C0-A7 — Control nulo de la batería (D-2, conservador)

Se registró **después** de ver la batería. Solo puede rebajar afirmaciones; no puede promover ninguna.

**Motivación**
- Una familia de finales (inicios R y E) muestra D_L entre 2.5 y 6 y D_s ≈ 3 con N = 343.
- D_eff no tiene ventana.
- Con N ≤ 343, el crecimiento logarítmico de un expansor también produce un D_L finito y grande.

**Control D-2:** para cada final candidato con N = 125, 216 y 343:
- se recablea su soporte fuerte A conservando los grados (10·|E| intercambios, PCG64);
- se recalculan el salto medio, D_L y D_s.

**Lectura:**
- Si el nulo reproduce D_L y D_s dentro de ±0.5, esas cifras no indican geometría. Se deben al grado y al tamaño.
- Solo una separación clara respecto del nulo (D_L(final) < D_L(nulo) − 0.5) se informa como «estructura de escala no trivial».

## 7. C0-F1 — Escalado y certificado de la subfamilia «R-3D» (prerregistro)

Este apartado se registra antes de cualquier corrida con N > 343. Lo autoriza la condición de continuación del Consejo: «solo si aparece una fase intermedia no trivial se ejecuta la dinámica completa y se mide D_eff, D_s, …».

**Objeto.** Las 10 celdas candidatas donde el inicio R da, con N = 343, D_s ≥ 2.5 y D_L entre 2.4 y 3.7, separado del nulo:

| c* | k* | Celdas |
|---|---|---|
| 1 | 6 | 18, 19, 20 |
| 2 | 8 | 36, 37, 38 |
| 4 | 12 | 55, 56 |
| 8 | 16 | 73, 74 |

**Corridas**
- Inicio R (mismo generador que en R6, con la semilla dependiente de N), semillas 0, 1 y 2.
- N = 512 en las 10 celdas.
- N = 729 en 4 celdas (19, 37, 55 y 73; a = 1, una por c*).
- max_steps = 40 000; el resto de la dinámica es idéntico.
- Total: 30 + 12 = 42 corridas.

**Medidas**
- Clase C0 (§4).
- Batería C0-A6.
- Control nulo C0-A7.
- **Certificado Ω-1.1 completo**, con umbrales sin modificar (`collect_run_evidence` y `assess_run`; mismo camino que el informe de L-3b).
- Referencia del certificado con el mismo N: RGG3 k12 (3 semillas) y T³ (8³ = 512 y 9³ = 729).

**D_L.** Se ajusta con N ∈ {216, 343, 512} para todas las celdas y con N ∈ {216, 343, 512, 729} para las 4 celdas con 729. Se usa la mediana de las semillas.

**Decisión de F1** (por celda; global = la mejor celda):

| Veredicto | Condición |
|---|---|
| **F1-POSITIVO** | En N máximo, ≥ 2/3 semillas DISPERSO-LOCAL **y** certificado sin códigos (veredicto candidato) |
| **F1-INDETERMINADO** | Hay DISPERSO-LOCAL, 2.5 ≤ D_L ≤ 3.5 y separación del nulo ≥ 0.5, pero el certificado falla. Todos sus códigos, salvo F10 (no convergencia), están también en la RGG3 de referencia con el mismo N: el tamaño no basta para discriminar |
| **F1-NEGATIVO** | Cualquier otro caso: D_L fuera de [2.5, 3.5]; códigos que la RGG3 de referencia no tiene; o pérdida de DISPERSO-LOCAL |

Un F1-POSITIVO **no** declara geometría emergente. Abre el paquete de confirmación (N-1, reproducibilidad, Θ > 0), que requiere ratificación del Consejo.
