# CONSEJO Ω — DECISIÓN DE SIGUIENTE FASE

> **Revisión 2 (2026-10-05).** Ver «REVISIÓN 2» más abajo: ratificación del usuario, afirmaciones rebajadas en 4 niveles, nueva prueba L-3 y criterios congelados del bloque L.

- **Fecha:** 2026-10-05.
- **Rama:** `claude/hopeful-galileo-88u1l1`.
- **Base auditada:** `539770f`.
- **Naturaleza:** documento de decisión. **No contiene código ni autoriza implementar nada** hasta que el usuario ratifique las resoluciones (§XII).
- **Principio rector:** «No queremos demostrar que Ω es correcta; queremos descubrir qué puede producir Ω y dónde los experimentos pueden demostrar que estamos equivocados.»

---

## REVISIÓN 2 (2026-10-05) — incorpora la crítica del usuario y el análisis Opus

**Esta revisión prevalece sobre el texto original donde haya conflicto.** Las secciones sustituidas llevan la marca «⟶ Rev. 2». El texto original se conserva para la trazabilidad.

### R2.0 Ratificación del usuario (2026-10-05)

| # | Resolución | Estado |
|---|---|---|
| D-1 | D2/R1/R3/R4/R8: APROBAR con el cambio propuesto (R1 modificada; R4 = veto necesario, nunca suficiente) | **RATIFICADA** |
| D-2 | s06: APROBAR (cadena completa sobre el commit congelado final + RP-1; la compuerta no se relaja) | **RATIFICADA** |
| D-3 | O-04b reducido: APROBAR | **RATIFICADA** (con la enmienda de clases de R2.4) |
| D-4 | O-06 → O-05: APROBAR CON MODIFICACIÓN. **L-1, L-2 y L-3 se ejecutan antes de comprometer el coste de O-05.** O-06 no modifica criterios de O-05 | **RATIFICADA** |
| D-5 | Ω-1.2 por ρ ≤ w_min: RECHAZAR; regla de lectura preregistrada | **RATIFICADA** |

**Instrucción global del usuario:** no ejecutar C0–C12 hasta corregir el documento:
- incorporar L-3;
- rebajar las dos afirmaciones teóricas;
- separar resultado demostrado → explicación propuesta → teorema externo aplicable → extrapolación no demostrada.

Esta revisión cumple esa instrucción. A partir de ella se ejecuta **solo** el bloque L (R2.5), en el orden de R2.3.

### R2.1 Análisis Opus de las objeciones (resumen)

El análisis Opus (solo lectura, con verificación numérica en el scratchpad) da la razón al usuario en las dos objeciones.

**Afirmación sobre cliques.** El documento original mezclaba tres cosas:
1. **Una cota de vecindario** para grafos k-regulares binarios: T ≤ N·C(k,2)/3, con igualdad si y solo si el grafo es la unión disjunta de K_{k+1}. Es correcta, pero **no es Kruskal–Katona**.
2. **Kruskal–Katona/Lovász** (a m fijo, binario): el extremal es **una sola clique** (colex), no una unión de cliques.
3. **Una extrapolación no demostrada** a grados libres, pesos continuos y mínimos locales.

Además:
- el RGG3 no es regular, así que la clase de comparación «mismos grados» está prácticamente vacía;
- a m fijo con grados libres, Σ C(k_i,2) es convexa;
- con pesos, S_dens = Σw² no es constante a Σw fija;
- confundir mínimo global con mínimo local (o con metaestabilidad) es exactamente lo que L-3 debe cerrar.

**Afirmación Θ > 0.** Chatterjee–Diaconis (2013) trata ERGM **binarios**. Ω tiene pesos continuos con medida de Lebesgue y energía cuadrática. La traslación es una **analogía**, no una prueba. Además, la cercanía en métrica de corte a un grafón constante no controla estructura con o(N²) aristas, y el teorema no dice nada sobre N finito ni sobre equilibración.

### R2.2 Escala de afirmación (sustituye «Hallazgo central» y «Hallazgo para Θ > 0»)

#### Afirmación «triángulos → cliques»

| Nivel | Contenido |
|---|---|
| **1. Demostrado** (prueba corta; L-1 y L-3a lo verifican en código) | (i) S_smooth = 2Σk² − 12T en binario y 2Σ_i k_i s_i − 2 tr W³ en ponderado. (ii) A secuencia de grados binaria fija, S = const − (α + 12η)T. (iii) Cota de vecindario k-regular con su caso de igualdad. (iv) **Lemas KKT de S0** (gradiente G = −αW² + 2βW + 2γ(u_i+u_j), verificado contra `gradient.py`): **P-KKT** (γ = 0): todo cero con codegrado > 0 viola KKT, así que el soporte de cualquier punto KKT, con pesos continuos, es una unión disjunta de cliques. **P-KKT-reg**: en estados regulares γ desaparece del gradiente, así que **ningún γ estabiliza T³**. **P-KKT-γ**: un cero entre nodos con u_i + u_j ≤ 0 exige codegrado 0. **Lema de escala**: toda arista de peso 1 exige c_ij ≥ (N−2)/α̂ (más el término γ), por lo que las cliques KKT son extensivas. (v) **Ω-B** (multiplicador λ' = ν/dt del código): condición de separación de codegrados. Con λ' > 0, los ceros admiten codegrado hasta λ'/α, así que **en Ω-B el soporte KKT NO es necesariamente una unión disjunta de cliques**. Los halos de O-04 son puntos KKT exactos **solo con γ > 0**. |
| **2. Explicación propuesta** | F0/F1 en S0 y las cliques (disjuntas, solapadas, con halo) en Ω-B ocurren porque la recompensa por triángulo exige un codegrado extensivo ≥ (N−2)/α̂, que ningún grafo geométrico de grado acotado tiene. |
| **3. Teorema externo aplicable** | Kruskal–Katona/Lovász: entre grafos binarios con m fijo, T es máximo en la clique colex. Es el mínimo global binario de Ω-B con γ = η = 0. |
| **4. Extrapolación NO demostrada** | «Toda geometría tiene más acción» con pesos, grados libres o restricciones arbitrarias. «Ninguna dinámica de la familia selecciona geometría», en particular en Ω-B con λ' > 0 y γ̂ alto, donde el piloto Opus encontró aristas RGG que se autosostienen localmente (R2.6). |

#### Afirmación «Θ > 0 → sin estructura»

| Nivel | Contenido |
|---|---|
| **1. Demostrado** | Todos los términos escalan en régimen denso (αT/Θ, βΣw²/Θ y γS_deg/Θ son O(N²)). Un término ≥ 0 que vale 0 en un maximizador constante no lo desplaza, **condicionado** a que el resultado con γ = 0 se cumpla. |
| **2. Explicación propuesta** | MF-1: aristas casi i.i.d. con una inclinación autoconsistente. |
| **3. Teorema externo (por analogía)** | La teoría de ERGM densos (Chatterjee–Diaconis 2013) da resultados de degeneración que hacen **plausible** que la rama Θ > 0 tampoco produzca geometría. **La extensión al modelo continuo de Ω requiere demostración o validación computacional.** |
| **4. Extrapolación NO demostrada** | Pesos continuos con Lebesgue, N finito (64–100), equilibración real y ausencia de estructura con o(N²) aristas. |

**Frase del hallazgo central (reemplaza la anterior y la conclusión de §IX):**

> «La evidencia actual indica que S0 y Ω-B no producen una fase geométrica en el dominio estudiado. Existe además una hipótesis fuerte de que la forma actual del funcional favorece estados densos y altamente triangulados frente a estructuras geométricas. Esta hipótesis debe someterse a pruebas extremales explícitas (L-1 → L-2 → L-3) antes de concluir que la familia completa de dinámicas es incapaz de producir geometría.»

### R2.3 Orden revisado (sustituye §VI en su orden de ejecución)

```text
L-1  identidades y cotas exactas del funcional (INCORPORAR, tests)            [minutos]
  ↓
L-2  paisaje extremal: acción de referencias concretas, incluidas
     geométricas ponderadas optimizadas a Σw fija (PROBAR)                     [< 1 h]
  ↓
L-3  ¿los estados que la dinámica selecciona son incompatibles con la geometría 3D?
     L-3a analítica: lemas KKT (INCORPORAR, tests)                           [minutos]
     L-3b dinámica: estabilidad desde arranques geométricos (PROBAR)           [~2 h, 4 CPU]
  ↓
┌──────────────────────────────┬─────────────────────────────────────────────┐
│ SÍ (incompatible)            │ NO / METAESTABLE / NO-decorado               │
↓                              ↓                                             
limitar la familia;            investigar el mecanismo; O-05 tiene valor     
O-05 solo como prueba de       real, pero SE ENMIENDA ANTES (arranques       
MF-1 (analogía de CD);         geométricos y celdas L-3b) → nuevo Consejo    
nueva hipótesis dinámica                                                     
```

- **En paralelo con L-2 y L-3 (independiente):** O4B-1, con las clases enmendadas (R2.4).
- **Después de L-3 y de congelar el código:** la cadena de la compuerta (RP-1 → s06 → O-06), y O-05 según el resultado de L-3.
- **Nada del bloque L modifica umbrales del certificado ni criterios de O-05.**

### R2.4 Clasificación estructural compartida (congelada; la usan O4B-1 y L-3b)

Para un estado W (N×N, simétrico, en [0,1]):
1. **Vacío:** max W ≤ 1e-6.
2. **Uniforme:** cv de las entradas fuera de la diagonal < `uniform_cv_max` (1e-6).
3. **Grafo fuerte** A_s = {W > 1/2}. **Estrato débil** = {1e-6 < W ≤ 1/2}; si no está vacío se marca **+halo** y se informa si su soporte es a su vez una unión de cliques.
4. **Componentes de A_s de tamaño ≥ 2:**
   - **clique** si su densidad interna es 1;
   - **cliques solapadas** si no es una clique y todas sus aristas pertenecen a cliques maximales de tamaño ≥ 4, con a lo sumo 10 de ellas;
   - **otro** en cualquier otro caso.
5. **Clase del estado:**
   - `vacío`;
   - `uniforme`;
   - `clique_única` (una componente clique);
   - `multi_clique` (≥ 2 componentes, todas clique);
   - `cliques_solapadas` (≥ 1 componente solapada y ninguna «otro»);
   - `otro` (alguna componente «otro», o A_s vacío con estrato débil no uniforme).
6. **Diagnósticos sin voto:**
   - κ_inj = t_inj(K3)/t_inj(K2)^{3/2} por componente (densidades inyectivas; vale 1 en una clique);
   - número de P3 abiertos en el soporte {W > 1e-6};
   - residuos KKT.

Esta definición se congela antes de ver ningún resultado de O4B-1 ni de L-3b. **Enmienda respecto de la ficha O4B-1 original:** se añade la clase «cliques solapadas», porque el piloto Opus (exploratorio, R2.6) mostró que existe en Ω-B. Sin ella, O4B-1 podría fallar por un motivo no geométrico.

### R2.5 Fichas revisadas (sustituyen L-1, L-2 y O4B-1 de §VII; L-3a y L-3b son nuevas)

#### L-1 — Identidades y cotas exactas · INCORPORAR (tests analíticos)
- **Prueba SOLO estas afirmaciones:**
  1. S_smooth, comparada con la definición directa Σ_ij W_ij Σ_k (W_ik − W_jk)² (tolerancia relativa 1e-10), en binario (= 2Σk² − 12T) y en ponderado (= 2Σk_i s_i − 2 tr W³).
  2. Bajo double-edge swaps que preservan grados, ΔS = −(α + 12η)ΔT exactamente.
  3. Cota de vecindario T ≤ N·C(k,2)/3 para grafos k-regulares, con igualdad si y solo si es la unión disjunta de K_{k+1}. Se verifica exhaustivamente sobre el atlas de grafos de ≤ 7 nodos.
  4. Lovász: si m = C(x,2), entonces T ≤ C(x,3), sobre todo el atlas.
  5. tr W³ ≤ (ΣW²)^{3/2} ≤ (ΣW)^{3/2} en [0,1] (200 matrices U(0,1) con PCG64 y semilla 20261005, más casos extremos).
- **No prueba nada** sobre grados libres, pesos óptimos ni mínimos locales.
- **Éxito:** 100%. **Fracaso:** cualquier excepción. En ese caso se detiene el bloque L y se vuelve al Consejo.

#### L-2 — Paisaje extremal sobre referencias concretas · PROBAR
- **Pregunta:** a Σw fija (la restricción de Ω-B, la más favorable a estructura), ¿alguna referencia geométrica **evaluada** tiene menor acción que la mejor referencia no geométrica evaluada?
- **N = 216.** Densidades ρ ∈ {6/215 (la de T³ 6³), 12/215 (RGG3 k12), 0.1}.
- **Malla:** α̂ ∈ {0.5, 1, 1.5, 2, 3}; γ̂ ∈ {0, 1, 10, 100}; η̂ = η(N−2)/β ∈ {0, 0.1, 1}; μ = 0.
- **Geométricas:**
  - T³ 6³ binaria;
  - RGG3 en el toro, binaria (k12; 5 semillas);
  - **RGG3 ponderada optimizada:** w_ij = min(1, a·φ(d_ij/r)) con φ ∈ {escalón, lineal 1 − x, suave (1 − x²)}₊ y r en una malla de 12 valores (0.5–3 veces el radio k12). a se ajusta por bisección para cumplir Σw = W0 y se toma el mínimo de S sobre (φ, r);
  - **T³ decorada 3³ ⊗ K8:** 27 bloques K8; los bloques de sitios vecinos se unen por emparejamiento vértice a vértice; grado 13.
- **No geométricas:**
  - clique colex con la misma masa;
  - unión de cliques K_{k+1} (con k el de la referencia geométrica comparada);
  - caveman conectado;
  - ER G(n,m);
  - uniforme w = ρ.
- **Igualación de masa:** las binarias con m ≠ W0 se igualan por la arista parcial colex (o se descartan si no pueden igualarse; se informa).
- **Observable:** ΔS = [min S(geométricas) − min S(no geométricas)] / (β·C(N,2)), por punto y semilla.
- **Éxito («la hipótesis energética se sostiene para las referencias evaluadas»):** ΔS > 0 en el 100% de los puntos, en ≥ 3/5 semillas.
- **Fracaso:** ΔS < 0 en algún punto, en ≥ 3/5 semillas. Esos puntos se marcan como prioritarios para L-3b; no cambian ningún criterio.
- **Alcance:** compara representantes, no minimiza globalmente. Su conclusión se redacta como «ninguna de las referencias geométricas evaluadas…», nunca como «ninguna geometría…».

#### L-3a — Lemas KKT · INCORPORAR (tests analíticos)
1. La función de residuos KKT (S0 y Ω-B, con λ' estimado) coincide con las condiciones del gradiente del código.
2. **P-KKT:** S0 con γ = 0, N = 30, 20 semillas y α̂ ∈ {0.5, 1.5, 2.5, 4}, evolucionado hasta convergencia. El soporte {W > 1e-6} tiene 0 P3 abiertos y residuo KKT ≤ 1e-7 relativo.
3. **P-KKT-reg:** T³ 6³ y 8³ binarias no son KKT en S0 para ningún (α̂, γ̂) de la malla de L-3b (violación > 0).
4. **RGG3 k12** (N = 216 y 512, 3 semillas) binaria no es KKT en S0 para ningún punto de la malla.
5. **Lema de escala:** en todos los finales convergidos del punto 2, toda arista con W = 1 cumple α c_ij ≥ 2β + 2γ(u_i+u_j) − tolerancia.
6. **Ω-B:** el estado «K22 + halo» (N = 100, f = 0.5) es KKT con γ̂ = 10 y **no** lo es con γ̂ = 0. T³ binaria con su propia ρ no es KKT en Ω-B para α > 0.

- **Éxito:** 100%. **Fracaso:** cualquier excepción. Se detiene y se vuelve al Consejo.

#### L-3b — Estabilidad de arranques geométricos · PROBAR (criterios congelados aquí)
- **Pregunta:** ¿los estados que S0 y Ω-B seleccionan **cuando arrancan cerca de una geometría** siguen siendo geométricos?
- **Inputs (N = 216):**
  - T³ 6³;
  - RGG3 k12 en el toro (3 grafos; la semilla s usa el grafo s);
  - T³ decorada 3³ ⊗ K8;
  - control negativo ER G(n,m) con la misma m que el RGG3;
  - control de clases: unión de 27 K8.
- **S0:** α̂ ∈ {0.5, 1, 1.5, 1.9, 2.1, 2.5, 3} × γ̂ ∈ {0, 1, 10, 100} (28 celdas).
- **Ω-B:** ρ = la densidad del input; α̂ = f·α̂_c(γ̂, ρ, N), con α̂_c = `uniform_state_threshold`; f ∈ {0.25, 0.5, 1, 1.5, 3} × γ̂ ∈ {0, 1, 10, 100} (20 celdas).
- **Ruido:**
  - W0 = clip(A + εξ) con ξ gaussiano simétrico (PCG64, `master_entropy = 20261005`), proyectado en Ω-B;
  - ε = 1e-2 (3 semillas) y ε = 1e-3 (3 semillas);
  - ε = 0, solo como diagnóstico del subespacio simétrico, sin voto.
- **Dinámica:** `evolve` / `evolve_fixed_density` con `DynamicsConfig` por defecto salvo `max_steps = 40000`.
- **Prueba de estabilidad:** todo final clasificado «geométrico-persistente» se re-perturba con ε = 1e-2 y se re-evoluciona.
- **Observables:**
  - clase R2.4;
  - R_orig (fracción de Σw que queda en las aristas del input);
  - ρ_S = Spearman(W_ij, −d0_ij) sobre los pares con distancia de saltos del input d0 ≤ 4;
  - fracción de la componente gigante de {W > 0.1·max W};
  - P3 abiertos, residuos KKT y λ';
  - estado de la corrida y deriva de R_orig en los últimos 10% de pasos si termina en MAX_STEPS.
  - El certificado completo **solo** se calcula sobre los finales geométrico-persistentes (como confirmación de nivel 3; no se exige para la decisión).
- **«Geométrico-persistente» (congelado):** R_orig ≥ 0.5 **y** ρ_S ≥ 0.5 **y** componente gigante ≥ 0.5N **y** clase ∉ {clique_única, multi_clique, cliques_solapadas, vacío, uniforme}.
- **Validación del clasificador (antes de las corridas):**
  - los inputs T³, RGG3 y decorada **sin evolucionar** deben clasificarse como geométrico-persistentes;
  - ER y la unión de K8 no.
  - Si falla, se vuelve al Consejo sin ajustar.
- **Resultado por celda e input:** persistente si ≥ 2/3 semillas con ε = 1e-2 son geométrico-persistentes **y** estables tras la re-perturbación.
- **Decisión global (congelada):**
  - **SÍ (incompatible):** ninguna celda persistente con ε = 1e-2, para ningún input. Además, todo final convergido con residuo KKT ≤ 1e-7 relativo (si no, se informa como «no KKT», sin voto).
  - **NO:** alguna celda persistente con un input **no decorado** (T³ o RGG3) → la geometría es un mínimo local; O-05 gana valor real y se enmienda antes de ejecutarse.
  - **NO-decorado:** solo persiste la T³ decorada → se reabre el control diferido de la retícula decorada (es geometría a gran escala con cliques locales; la pregunta pasa a ser si el certificado la reconoce).
  - **METAESTABLE:** alguna celda persiste con ε = 1e-3 pero no con 1e-2, o termina en MAX_STEPS con |deriva de R_orig| > 0.01 → se informa la vida media; si ocurre con γ̂ alto, O-05 solo tiene valor si **antes** se enmienda para incluir arranques geométricos (desde U(0,1) nunca visitaría esa cuenca).
- **Coste:** S0 ≈ 980 corridas de segundos. Ω-B ≈ 700 corridas de ~10–60 s, en 4 procesos. Total ≈ 1.5–3 h. N = 512 solo para celdas no-SÍ (nuevo Consejo).

#### O4B-1 — enmendada
- Igual que en §VII, pero con la clasificación de R2.4 (incluye `cliques_solapadas` y +halo) y residuos KKT.
- **Éxito:** 100% de los finales Ω-B de O-04 en {vacío, uniforme, clique_única, multi_clique, cliques_solapadas} (± halo) **y** residuo KKT ≤ 1e-7 relativo en los convergidos.
- **Fracaso:** algún `otro`, o un convergido que no sea KKT → explicación debilitada; se reabre el mapa de cuenca.
- **κ:** solo informe.
- **Input:** los `.npz` de `runs/omega11/o04_ablation/passports`, verificados por sha256 frente a `array_digests` (si alguno no coincide, se excluye y se informa).
- **Los finales S0 de O-04 se clasifican también**, como informe.

### R2.6 Declaración del piloto exploratorio (transparencia M§44)

Antes de congelar L-3b, el agente Opus ejecutó **pilotos exploratorios** en el scratchpad (no versionados, sin pasaporte):
- S0 desde T³ 6³ y RGG3 N = 512 → F0 en todas las celdas probadas.
- Ω-B desde T³ → uniforme o K37 colex.
- Ω-B desde RGG3 N = 512 con f = 0.5 y γ̂ = 10 → 45% del peso en aristas originales, repartido en cliques localizadas y cliques solapadas, con 254 componentes y sin geometría conexa.

La malla y los umbrales de L-3b **no se ajustaron para cambiar esos resultados**. La malla incluye γ̂ alto y f ≤ 1 **porque** el piloto señaló esa región como la única no excluida por escala. Esto es una elección de dónde mirar, declarada; no de qué cuenta como éxito. Los resultados del piloto no cuentan como evidencia de L-3b.

### R2.7 Estado de los compromisos

| Ítem | Estado tras la Rev. 2 |
|---|---|
| L-1, L-3a, L-2, L-3b, O4B-1 | **APPROVED (ejecución inmediata)**, en el orden de R2.3 |
| RP-1 + s06 + O-06 | APPROVED; se ejecutan **después** del bloque L, sobre el commit congelado final |
| O-05 (+ MF-1, PG-1, piloto) | **ESPERA el resultado de L-3** (D-4 modificada). Si L-3 = SÍ, se ejecuta como prueba de la analogía; si no, se enmienda antes |
| N-1 | APPROVED, después del bloque L |
| Resto de §V | Sin cambios |


### R2.8 Enmiendas operativas del bloque L (registradas ANTES de cualquier resultado oficial)

| ID | Ítem | Cambio | Motivo | ¿Puede crear un positivo? |
|---|---|---|---|---|
| L-A1 | L-2, igualación de masa | Una referencia geométrica binaria cuya masa propia difiera de W0 en más de 10% (`MATCH_REL_TOL`) se **descarta** en ese ρ y se informa. | Operacionaliza «no pueden igualarse»: igualar T³ (648 aristas) a ρ=0.1 (2322) añadiría 1674 aristas colex y ya no sería geometría. | No: solo reduce el conjunto geométrico. |
| L-A2 | L-2, diagnóstico | Se añade `delta_S_local` **sin voto**: igual que ΔS, pero solo con las geométricas binarias y las ponderadas con r ≤ 1.5·r12. | En el smoke, la RGG ponderada con r grande degenera en casi uniforme, y ΔS > 0 podía deberse solo a eso (ΔS ~1e-5). La decisión congelada no cambia. | No (no vota). |
| L-A3 | L-3a, punto 6 | El test K22 + halo se evalúa con ρ = 0.05 (el halo h ≈ 0.0055 de O-04) y con ρ = 0.1. Se registra un matiz: **con el mismo α̂ = 0.5·α̂_c(10) y γ̂ = 0, el halo es KKT con ρ = 0.1 pero no con ρ = 0.05.** | Error de transcripción en el encargo (ρ = 0.1 da h = 0.088). | — |

**Corrección del nivel 1 de R2.2:** «los halos de O-04 son puntos KKT exactos **solo con γ > 0**» se mantiene para las celdas de O-04 (ρ = 0.05, h ≈ 0.0055). **No es una afirmación general:** con ρ = 0.1 y el mismo α̂ el halo es KKT sin γ.

**Hallazgo de L-3a (punto 2):** se cumple, pero es casi vacuo. Desde U(0,1) con N = 30 y γ = 0, todos los finales son el vacío o K_N. Se añadió un test con 3 K10 disjuntas, que es KKT si y solo si α̂ ≥ (N−2)/(m−2) = 3.5. Así P-KKT y el lema de escala tienen contenido no trivial.

---

## 0. AUDITORÍA DE ESTADO

### 0.1 Hechos verificados en la repo (no en la documentación)

| Ítem | Valor verificado | Cómo |
|---|---|---|
| Rama | `claude/hopeful-galileo-88u1l1` | `git branch` |
| HEAD | `539770f` (27 commits) | `git log` |
| Árbol | limpio | `git status` |
| Versión | `1.1.0` (`omega/__init__.py`, `pyproject.toml`) | lectura |
| Baseline Ω-1.0 | `BASELINE_COMMIT = 9dfbea7…` (`omega/baseline.py`), manifest sha256 en `tests/expected_baseline/`; tag no publicado (el proxy bloquea tags) | lectura |
| Commit de ejecución de los pasos 1–5 | `167efca` (`code_commit` en cada `summary.json`, `git_dirty=false`) | lectura |
| Código entre `167efca` y HEAD | **ningún cambio** en `omega/`, `tools/`, `tests/` ni `pyproject.toml`; solo `docs/` y `results/` | `git diff --stat 167efca HEAD -- omega tools tests pyproject.toml` (vacío) |
| Tests | 423 funciones `def test` en `tests/` (más los `test_*` de `omega/experiments/`); suite rápida re-ejecutada en esta auditoría (resultado en §0.4) | `pytest -m "not slow"` |
| Experimentos v11 | p00, O-00…O-11 y s06 en `omega/experiments/v11/`. O-07 y O-09 no están en `STEP_ORDER` (comparten la compuerta de O-08) | lectura |
| Ejecutados (full) | p01, p02, O-00, O-01, O-02, O-03 y O-04, todos `complete=true` | `results/omega11_steps_1_5/*/summary.json` |
| No ejecutados | s06, O-06, O-05, O-07…O-11 | ausencia de `summary.json` |
| O-04b | **no existe** en código; es solo una propuesta en `results/omega11_steps_1_5/ANALYSIS.md` §7 | grep |
| Pasaportes | v1.1 JSON (`passports_json.tar.gz`). Los `.npz` están en `runs/omega11/`, ignorado por git y presente en este contenedor | `ls` |
| Puerta (gate) | `STEP_ORDER = p01, p02, o00…o04, s06, o06, o05, o08, o10, o11`. `require_prerequisites` exige que **cada** summary previo tenga `code_commit == HEAD` | `gate.py:136-158` |
| Funcional | `action` = −αT + βS_dens + γS_deg + ηS_smooth − μS_mass. η = μ = 0 en **todos** los experimentos (valores por defecto de `FunctionalParams`; `point_params` no los fija) | `functional.py`, `settings.py:69-72`, `scan.py:106` |
| Escala de Θ | Θ = Θ̂·β (`langevin.temperature`): escala densa, con la energía por arista de orden β | `langevin.py:35` |
| `omega/causality` | stub vacío («reservado, sin código») | lectura |

### 0.2 Contradicciones detectadas (documentación ↔ código ↔ resultados)

| # | Contradicción | Consecuencia | Resolución propuesta |
|---|---|---|---|
| **K-1** | La puerta exige `code_commit == HEAD` para todos los pasos previos. HEAD (`539770f`) ≠ `167efca` por commits de documentación. | **Hoy no se puede ejecutar s06 desde HEAD**, aunque el código sea idéntico. | §II-D2: re-ejecutar la cadena completa sobre el commit congelado final (prueba de reproducibilidad RP-1). No se relaja la puerta. |
| **K-2** | La instrucción y el resumen tratan la isotropía MDS como «solo diagnóstico», pero en el código `isotropy_ok` es **uno de los 14 campos obligatorios** del certificado (`taxonomy.py:139,221`). | Funciona como **veto** (condición necesaria), no como diagnóstico puro. | §II-D1/R4: se mantiene como veto necesario y nunca suficiente. Degradarlo relajaría el criterio después de ver resultados y queda prohibido. |
| **K-3** | O-05 preregistra que **no puede existir región geométrica con N ≤ 100** (R1), y su malla full es N ∈ {64, 100}. | Su variable de decisión (`dominant_geometric_region`) es **vacua por construcción**: solo puede responder «no». O-05 tal como está no puede detectar P1. | §Decisión B: añadir una pantalla pre-geométrica preregistrada (PG-1) y una predicción de campo medio (MF-1). |
| **K-4** | O-06 usa en Ω-B ρ = 0.1 = w_min (`OMEGA_B_RHO`). | Es la zona del artefacto de redondeo documentado en `ANALYSIS.md` §5. | §II-D5: regla de lectura preregistrada, sin nueva versión. |
| **K-5** | «Decisiones 1–5» no aparece con ese rótulo en ningún archivo. | Riesgo de inventar contenido. | §II: se reconstruyen solo desde la lista numerada de 5 decisiones del último informe y de `ANALYSIS.md` §7. |
| **K-6** | R2, R5, R6 y R7 (`DESIGN.md` §6) se implementaron tal como se propusieron, pero no consta aprobación explícita del usuario (solo D1 y la ejecución constan en `AMENDMENTS.md`). | Decisiones implícitas. | §II-D1: ratificación explícita. |
| **K-7** | O-02 cumplió sus expectativas **vacuamente** (ningún estado Ω tuvo ventana de escala). | La suite de distancias no se ha ejercido sobre un estado Ω no trivial. | Se registra. No se repite hasta que exista un estado no trivial. |

### 0.3 Las cuatro categorías

**CONFIRMADO** (computacionalmente, en el dominio estudiado):
1. **Umbrales homogéneos de S0:** α̂ = 1 (estabilidad local de W = 1), 3/2 (H(1) < H(0)) y 2 (frontera de cuenca desde U(0,1)). Confirmados por tests analíticos y por O-01: F0 en el 100% con α̂ ≤ 1.9 y F1 en el 100% con α̂ ≥ 2.1, para todo γ̂ (N = 200, 10 réplicas).
2. **γ es inerte** en la trayectoria homogénea: NO_DEGREE ≡ FULL (O-04).
3. **T es el único término que genera estructura:** T = 0 da F0 en el 100% (N = 100 y 200); TRIANGLES_ONLY da F1.
4. **β vacía la red;** NO_DENSITY da F1 incluso con α̂ = 0.5.
5. **Ω-B determinista** produce cliques en la frontera (puntos KKT), reproducibles bit a bit, o el estado uniforme. Nunca ventana dimensional.
6. **El certificado no da falsos positivos sobre los nulos probados:** 0/65 nulos pasan; d_volume y manifold_proxy_ok tienen 0 falsos positivos. WS p = 0.2 (D_s = 3.22, D_W = 2.92) queda rechazado por F3.
7. **Control positivo:** RGG3 k12 pasa 5/5 (N = 800) y, con A-17, 3/3 a N = 3000.
8. **Reproducibilidad:** gate, pasaportes, config_hash y golden de Ω-1.0 funcionan; las corridas son deterministas por semilla.

**REFUTADO** (en el dominio estudiado):
1. **P1 para S0 determinista (Θ = 0)** con N = 100–200, todo (α̂, γ̂) de la malla y desde U(0,1).
2. **P1 para Ω-B determinista** con ρ ∈ {0.05, 0.1, …}, γ̂ ∈ {0, 1, 10} y N = 100.
3. **La predicción preregistrada «uniforme» de Ω-B** (γ̂ = 10, f = 0.5): falló en 10/10 semillas. Se obtuvieron cliques.
4. **D ≈ 3 como evidencia suficiente:** WS y árboles lo refutan (regla `RULE_D3_NEVER_SUFFICIENT`).

**ABIERTO:**
1. **S0 con Θ > 0** (O-05). No ejecutado.
2. **La explicación Kruskal–Katona de Ω-B.** Es una hipótesis: no se ha verificado sobre los estados finales.
3. **Tamaño finito** (O-06, N ≥ 800 para la dinámica). No ejecutado.
4. **Términos η (S_smooth) y μ (S_mass) de M§15.** Existen en el código, pero nunca se han activado.
5. **Especificidad del certificado** frente a nulos con igual grado **y** igual clustering (no probado).
6. **Falsos negativos del certificado:** RGG2 k10 da F9 y RGG3 k8 falla la prueba de anillos (0.25).

**NO PERMITIDO TODAVÍA:**
1. Interpretar cualquier resultado como nivel ≥ 4 («geometría continua efectiva»). Ningún estado Ω ha alcanzado siquiera el nivel 2 robusto.
2. Modificar S0, los umbrales del certificado o D_eff (Ω-1.0 congelado) por los resultados 1–5. El criterio F0/F1 como resultado válido está preregistrado por el usuario.
3. Ejecutar O-08, O-10 y O-11 sobre estados triviales: serían vacuos, como O-02 (K-7).
4. Introducir mecánica cuántica.
5. Escalar la dinámica a N ≥ 1500, o escribir un backend disperso, sin una señal previa que lo justifique.
6. Crear Ω-1.2 por redondeo o por ρ ≤ w_min.

### 0.4 Suite rápida en esta auditoría
`pytest -m "not slow"` sobre `539770f`: **598 passed, 28 deselected (slow), 0 fallos, 317 s**. El código coincide con el estado documentado.

---

## Deliberación del Consejo (resumen adversarial)

### Hallazgo central (teórico de grafos y físico teórico) ⟶ Rev. 2 (R2.2): sobreafirmado; reemplazado

**Proposición P-Ω1, el paisaje favorece cliques** (argumento; se verifica en L-1).

Toda la familia de M§13–M§15 con α, η ≥ 0 premia los triángulos:
- −αT premia los triángulos directamente.
- Para grafos binarios, S_smooth = Σ_ij W_ij Σ_k (W_ik − W_jk)² = 2Σ_i k_i² − 12T. Con η > 0 añade −12ηT y penaliza las cuñas abiertas.
- S_dens y S_mass solo dependen del número de aristas m.
- S_deg solo depende de la secuencia de grados.

**Consecuencia.** Entre grafos binarios con igual N, m y secuencia k-regular:
- todos los términos salvo T coinciden;
- T ≤ N·C(k,2)/3, con igualdad si y solo si el grafo es la unión disjunta de K_{k+1}.

**Comparación con referencias geométricas.**
- El RGG3 tiene clustering asintótico 15/32 ≈ 0.47; por tanto T_RGG3 ≈ 0.47 T_max.
- La retícula cúbica T³ (k = 6) tiene T = 0.

Por eso, **para cualquier α > 0 o η > 0, una partición en cliques tiene acción estrictamente menor que cualquier grafo geométrico con los mismos N, m y grados.**

**Caso con peso fijo.** El análogo ponderado (Kruskal–Katona para grafones, t(K3) ≤ t(K2)^{3/2}) explica Ω-B. Además:
- los puntos estacionarios KKT de Ω-B con multiplicador λ son estacionarios de S0 − λ·S_mass;
- así que la rama μ comparte los puntos estacionarios de Ω-B. No es una ruta nueva a Θ = 0.

**Lectura del físico.** No es que la dinámica «no encuentre» el mínimo geométrico: **la funcional no tiene un mínimo geométrico.** Ninguna dinámica de descenso, determinista o con otra inicialización, puede seleccionar geometría a Θ = 0 dentro de esta familia y a densidad y grados comparables.

### Hallazgo para Θ > 0 (físico teórico y experto en Monte Carlo) ⟶ Rev. 2 (R2.2): rebajado a analogía

**Teorema de partida (Chatterjee–Diaconis 2013, Ann. Stat., Thm 6.1).** Para modelos exponenciales arista–triángulo en escala densa con coeficiente de triángulos ≥ 0, el maximizador del problema variacional es un grafón **constante**. El grafo es asintóticamente indistinguible de un Erdős–Rényi; hay una transición de primer orden entre ramas de densidad baja y alta.

**Por qué aplica a S0.**
- La escala de Θ en S0 es la densa (§0.1).
- Añadir γ·S_deg, que es ≥ 0 y vale 0 en los grafones constantes, no puede desplazar el maximizador.
- La extensión a pesos continuos en [0,1] es análoga (reducción puntual), **pero no está verificada formalmente aquí**.

**Predicción a priori para O-05 (η = μ = 0).**
- Estados tipo ER ponderado (F3/F4 en grafos densos, F0/F1 en los extremos).
- Metaestabilidad y no-equilibrio a Θ̂ bajo cerca de α̂ ∈ [1.5, 2].
- **Ninguna estructura geométrica.**

**Objeción del revisor metodológico (aceptada).** El teorema es asintótico. A N = 64–100, la comprobación empírica es barata en relación con su valor de cierre, y es la única forma legítima de cerrar la rama Θ > 0 de S0. O-05 se mantiene, pero como **prueba del teorema**, no como búsqueda de geometría (K-3).

### Objeciones por rol (registradas)

- **Geómetra:**
  - Ninguna medida dimensional se activó en 1290 corridas.
  - Cualquier positivo futuro de una funcional que premie triángulos estará «decorado» con cliques locales.
  - Hay que saber si el certificado reconoce una retícula 3D decorada con cliques (diferido, condicional).
- **Estadístico:**
  - O-05 tiene 240 celdas por motor. Toda bandera nueva debe replicarse en N = 64 **y** N = 100 y contrastarse con un nulo por permutación de pesos.
  - Las ventanas y umbrales de D_eff y del certificado ya están congelados (R3). Prohibido elegir ventanas a posteriori.
- **Monte Carlo:**
  - Metropolis es la referencia exacta; Langevin tiene sesgo O(dt) (A-2) y solo cuenta con confirmación de Metropolis (A-12).
  - A Θ̂ ≤ 0.03 se esperan cadenas no equilibradas cerca de la transición: se informan como `CHAIN_NOT_EQUILIBRATED`, nunca como fase.
  - El tempering paralelo se difiere.
- **HPC:**
  - La dinámica densa es O(N³) por paso. N ≤ 100 para O-05 es trivial; N = 800 cuesta horas por corrida.
  - No hay justificación para un backend disperso mientras ninguna rama produzca una señal.
- **Red-team:**
  - Una cadena de cliques (caveman conectado) es lo que esta familia tendería a producir si se forzara la conexión. Es casi 1D. Debe estar entre las referencias.
  - Falta un nulo con igual grado **e igual clustering** que RGG3: sin él no sabemos si el certificado premia la geometría o solo el clustering local.
- **Revisor metodológico:**
  - Toda modificación de O-05 ocurre **antes** de ejecutarlo y se motiva en teoría y en resultados de otros experimentos, no de O-05. Se registra como enmienda.
  - Las nuevas observables de O-05 **no votan en el certificado**: solo deciden si se escala a N grande (nuevo preregistro).

---

## I. Estado actual

| Área | Estado | Evidencia | Confianza |
|---|---|---|---|
| Umbrales homogéneos S0 (1, 3/2, 2) | CONFIRMADO | tests analíticos y O-01 (100%) | Alta |
| S0, Θ = 0 → F0/F1 | CONFIRMADO / P1 REFUTADA (N = 100–200) | O-01, O-04 | Alta |
| T como único motor de estructura | CONFIRMADO | O-04 (T = 0 → F0 al 100%) | Alta |
| γ inerte (trayectoria homogénea) | CONFIRMADO | O-04 (NO_DEGREE ≡ FULL) | Alta |
| Ω-B determinista → cliques KKT | CONFIRMADO / P1 REFUTADA (N = 100) | O-04, bit a bit | Alta |
| Explicación Kruskal–Katona de Ω-B | HIPÓTESIS | argumento P-Ω1; sin test | Media |
| Paisaje de la familia favorece cliques (P-Ω1) | HIPÓTESIS FUERTE | argumento combinatorio; sin test | Media-alta |
| S0, Θ > 0 | ABIERTO; prior teórico nulo (Chatterjee–Diaconis) | sin corridas | Media (prior) |
| Términos η, μ | ABIERTO; P-Ω1 predice que no ayudan | sin corridas | Media (prior) |
| Especificidad del certificado (nulos clásicos) | CONFIRMADO | O-03: 0/65 | Alta |
| Especificidad frente a igual grado y clustering | ABIERTO | — | — |
| Sensibilidad del certificado (falsos negativos) | LIMITADA | RGG2 k10 F9, RGG3 k8 anillos | Media |
| Tamaño finito / N ≥ 800 dinámico | ABIERTO | — | — |
| Reproducibilidad | CONFIRMADO (determinista por semilla) | pasaportes y golden | Alta |
| Puerta de ejecución | BLOQUEA s06 en HEAD (K-1) | `gate.py:155` | Alta |

---

## II. Decisiones 1–5

**Procedencia.** Ningún archivo usa el rótulo «Decisión 3/4/5». Las cinco decisiones son las de la lista numerada del último informe al usuario, que coincide con `results/omega11_steps_1_5/ANALYSIS.md` §7, «Decisiones del usuario»:
1. D2;
2. autorizar s06;
3. O-04b;
4. O-05 frente a O-06;
5. tratamiento de ρ ≤ w_min en Ω-1.2.

La «D1» de `AUDIT.md` (anillos) es otra numeración. Está **cerrada**: el usuario la aprobó y se aplicó como A-17. No se reabre.

### Decisión 1 — D2: confirmar R1, R3, R4 y R8 (`DESIGN.md` §6)

- **Pregunta:** ¿se aceptan el plan de tamaños (R1), el rigor de los umbrales (R3), la isotropía MDS (R4) y la expectativa honesta (R8)?
- **Evidencia:** RGG3 no tiene ventana con N ≤ 500; los umbrales están congelados (`config_hash`) y son estrictos (falsos negativos documentados); `isotropy_ok` es campo obligatorio (K-2); los resultados 1–5 coinciden con R8.
- **R1 — APROBAR CON MODIFICACIÓN:**
  - Dinámica hasta N = 800 en esta fase (O-05 con N ≤ 100).
  - N ∈ {1500, 3000} solo para controles RGG3 en O-06.
  - El backend disperso queda **DIFERIDO**. Se reabre si una celda supera la pantalla PG-1 (§VII) y exige O-05b con N ≥ 800.
- **R3 — APROBAR:**
  - Se mantienen todos los umbrales de `settings11.py`.
  - La sensibilidad de A-17 ([2, r_hi−1] frente a [2, r_hi−2]) queda como diagnóstico no decisorio.
  - Relajar cualquier umbral exigiría Ω-1.2 y repetir todos los controles (O-00, O-03).
  - Se aceptan explícitamente los falsos negativos conocidos (el sesgo conservador es preferible).
- **R4 — APROBAR CON ACLARACIÓN (K-2):**
  - La isotropía MDS es **condición necesaria (veto), nunca suficiente**.
  - Su inmersión no alimenta ninguna otra medida.
  - Degradarla a diagnóstico puro relajaría el certificado después de ver resultados: **RECHAZADO**.
- **R8 — APROBAR:** el texto de expectativa honesta se preregistra en §IX.
- **Ratificación de R2, R5, R6 y R7 (K-6):**
  - APROBAR R2 (ζ_R vota; D_res solo se informa), R5 (≥/> heredado y documentado) y R6 (regla de deriva de Ω-1.1) tal como están implementados.
  - R7 (Ω-B estocástica): **DIFERIDO**. Exige definir una medida sobre la variedad restringida. Se reabre si O-04b-lite encuentra estados «otros».
- **Pendiente para cerrar:** firma del usuario (§XII).

### Decisión 2 — Autorizar s06

- **Pregunta:** ¿se ejecuta s06 (diagrama de fases denso de S0, Θ = 0)?
- **Evidencia:**
  - Su resultado está preregistrado: F0 al 100% con α̂ ≤ 1.9, F1 al 100% con α̂ ≥ 2.1, 0 candidatos y χ ≤ 1e-12.
  - Valor informativo bajo; cierra el orden de P§22.
  - Es requisito de la puerta para O-06 y O-05.
  - K-1: hoy no puede correr desde HEAD.
- **Recomendación: APROBAR CON MODIFICACIÓN.**
  - s06 se ejecuta **sin cambios** dentro de la cadena completa p00 → O-04 → s06 → O-06 → O-05, sobre el commit congelado final de esta fase.
  - Re-ejecutar p00–O-04 en ese commit es a la vez la prueba de reproducibilidad **RP-1**: los sha256 de los arreglos deben coincidir bit a bit con los pasaportes congelados de `167efca`.
  - Coste ≈ 1.6 h + 35 min.
  - No se relaja la puerta.
- **Pendiente:** ninguno, salvo la ratificación.

### Decisión 3 — O-04b (Kruskal–Katona y cuenca de Ω-B)

- **Pregunta:** ¿se implementa O-04b (verificar la cota y mapear la cuenca ρ × γ̂ × f × c, ~3–6 h)?
- **Evidencia:** Ω-B ya está refutada para P1 en el dominio; la explicación no está probada.
- **Recomendación: APROBAR CON MODIFICACIÓN.**
  - Se aprueba **O-04b-lite** (PROBAR/DIAGNÓSTICO): reanálisis de los estados finales **ya existentes** de O-04, sin dinámica nueva, ~30 min.
  - Fichas en §VII (L-1 y O4B-1).
  - El mapa completo de cuenca queda **DIFERIDO**, porque no cambia la estrategia: Ω-B ya no es una ruta hacia P1.
- **Se reabre si:** O4B-1 encuentra estados finales que no son ni uniformes ni uniones de cliques («otros»), o si L-1 refuta P-Ω1.

### Decisión 4 — Orden y presupuesto de O-05 frente a O-06

- **Pregunta:** ¿qué va primero y con qué presupuesto?
- **Evidencia:**
  - `STEP_ORDER` fija O-06 antes que O-05.
  - O-06 no modifica ningún criterio de O-05: los umbrales del certificado están congelados y `size_robust` solo se usa para escalar.
  - O-05 es vacuo en su variable de decisión (K-3).
- **Recomendación: APROBAR CON MODIFICACIÓN.**
  - Se conserva el orden de la puerta: O-06 sin cambios (~1–3 h), y después O-05.
  - O-05 mantiene **su malla preregistrada intacta**.
  - Se le añaden, como enmiendas previas a su ejecución, la pantalla **PG-1** (DIAGNÓSTICO que decide solo el escalado) y la predicción **MF-1** (expectativa analítica).
  - Antes de lanzarlo se hace un **piloto de tiempos** (una celda Metropolis con N = 100, no contabilizada). Si estima más de 24 h, se vuelve al Consejo; no hay recortes silenciosos.
- **Pendiente:** ratificación; validación de PG-1 sobre controles **antes** de congelar (§VII).

### Decisión 5 — Tratamiento de ρ ≤ w_min (¿Ω-1.2?)

- **Pregunta:** ¿se crea Ω-1.2 para corregir la lectura de las celdas Ω-B con ρ ≤ w_min?
- **Evidencia:**
  - El artefacto no cambia ninguna conclusión (`ANALYSIS.md` §5).
  - Corregirlo cambiaría la regla de umbral (semántica de observable) solo para leer mejor un caso degenerado.
- **Recomendación: RECHAZAR como versión.** En su lugar, una **regla de lectura preregistrada** (DIAGNÓSTICO, solo documentación):
  - toda celda con ρ ≤ w_min **y** cv(W) < `uniform_cv_max` se reporta como «uniforme por cv; primario binario no interpretable»;
  - se aplica también a O-06 Ω-B (ρ = 0.1, K-4).
- **Se reabre si:** alguna decisión futura dependiera de una de esas celdas.

---

## III. Qué incorporamos ahora (APROBADO; requiere ratificación)

| ID | Qué | Etiqueta | Por qué no es nueva versión |
|---|---|---|---|
| INC-1 | Test analítico de P-Ω1 (L-1) en `tests/analytical/` | INCORPORAR | Es un test; no cambia la semántica. |
| INC-2 | Regla de lectura ρ ≤ w_min (Decisión 5) en `AMENDMENTS.md` | INCORPORAR (documentación) | Es solo lectura. |
| INC-3 | Ratificación de R1–R8 (Decisión 1) en `AMENDMENTS.md` | INCORPORAR (documentación) | — |
| INC-4 | Comprobación de reproducibilidad RP-1 (comparar sha256 frente a los pasaportes congelados) en `tools/` | INCORPORAR | Herramienta de auditoría. |
| INC-5 | MF-1 (predicción de campo medio) como expectativa analítica de O-05 | INCORPORAR | Es una predicción; no vota en el certificado. |

**Versión:** ninguna de estas cosas cambia la dinámica, el estado, los estimadores ni los criterios del certificado. **Se queda en Ω-1.1** (el número 1.1.0 no cambia).

## IV. Qué probamos ahora (por prioridad científica)

1. **L-2 — paisaje energético** (PROBAR). Decide si la familia S puede siquiera favorecer geometría. Es el test más decisivo y el más barato.
2. **O-05 + PG-1 + MF-1** (PROBAR). Única prueba empírica de la rama Θ > 0 de S0.
3. **O4B-1 — O-04b-lite** (PROBAR/DIAGNÓSTICO). Cierre explicativo de Ω-B.
4. **N-1 — nulo con igual grado y clustering** (PROBAR). Especificidad del certificado frente a positivos futuros.
5. **Cadena de la puerta:** RP-1 + s06 + O-06 sin cambios (obligatoria por `STEP_ORDER`).

## V. Qué dejamos fuera (y por qué)

| Ítem | Resolución | Razón | Qué lo reabre |
|---|---|---|---|
| Mapa completo de cuenca de Ω-B | DEFERRED | No cambia la estrategia | O4B-1 encuentra «otros» o L-1 refuta P-Ω1 |
| Ω-B estocástica (R7) | DEFERRED | Requiere medida en la variedad | Ídem |
| Dinámica con η > 0 (M§15) | DEFERRED | P-Ω1 predice que η premia triángulos | L-2 muestra un orden geométrico < cliques con η > 0 |
| Dinámica con μ ≠ 0 | DEFERRED | Comparte los puntos estacionarios KKT con Ω-B | Ídem |
| O-05b (N ≥ 800 a Θ > 0) | DEFERRED | Sin señal que escalar | Una celda supera PG-1 con primario ∉ {F0, F1, F2}, replicada en N = 64 y 100 |
| Backend disperso | DEFERRED | Coste sin señal | O-05b aprobado |
| Tempering paralelo / MALA | DEFERRED | Infraestructura sin señal | > 50% de las celdas Metropolis con α̂ ∈ [1.5, 2.5] y Θ̂ ≥ 0.1 sin equilibrar |
| O-08, O-10, O-11 (y O-07, O-09) | DEFERRED | Vacuos sobre F0/F1 (K-7) | Primer estado que pase la pantalla PG-1 o un campo dimensional |
| Retícula 3D decorada con cliques (control positivo) | DEFERRED (condicional) | Solo importa si hay escalado | O-05b aprobado |
| Nuevos controles «D ≈ 3 artificial» | REJECTED | Redundantes: WS p = 0.2 y árboles ya cubren el caso (0/65) | Un falso positivo nuevo |
| Degradar `isotropy_ok` a diagnóstico | REJECTED | Relajación post hoc | — |
| Ω-1.2 por ρ ≤ w_min | REJECTED | No semántico | Decisión 5 |
| Nueva hipótesis dinámica (fuera de M§13–15) | **DEFERRED, necesita al usuario** | No se inventa funcional. Cualquier propuesta debe pasar antes L-2 (§VI) | Documento del usuario con una funcional nueva |

## VI. Dónde bifurca Ω ⟶ Rev. 2 (R2.3): orden L-1 → L-2 → L-3

```text
Ω  (familia M§13–15: −αT + βS_dens + γS_deg + ηS_smooth − μS_mass)
│
├── [B0] Paisaje energético a Θ=0 — L-1 (teorema) + L-2 (tabla) ── BIFURCACIÓN PRINCIPAL
│     ├── P-Ω1 se confirma → TODA la familia con α,η ≥ 0 favorece cliques a (N,m,k) fijos
│     │     ├── S0, Θ=0 ........ CERRADA (O-01/O-04 + teorema)
│     │     ├── Ω-B / μ ......... CERRADA para P1; O4B-1 cierra la explicación
│     │     └── η > 0 ........... CERRADA sin dinámica (L-2 con η)
│     └── P-Ω1 se refuta (alguna referencia geométrica gana) → reabrir η/μ como PROBAR
│
├── [B1] S0, Θ>0, escala densa — O-05 + MF-1 + PG-1   (prior: grafón constante, ER)
│     ├── MF-1 se cumple y PG-1 es negativa en todo → CERRADA (empírica a N≤100 + teorema)
│     └── PG-1 positiva y replicada, primario ∉{F0,F1,F2} → O-05b (N≥800, nuevo Consejo)
│
├── [B2] Certificado (instrumento) — N-1 (igual clustering), O-06 (RGG3 por tamaño)
│     └── Si N-1 pasa algún nulo → el certificado no es específico; ningún positivo es creíble hasta repararlo
│
└── [B3] Nueva hipótesis dinámica (fuera de la familia) — SOLO con documento del usuario
      └── Puerta de entrada obligatoria: L-2 aplicado a la funcional propuesta
          (una referencia geométrica debe tener menor acción que cliques/ER/uniforme a igual restricción)
```

**Corrección del Consejo al árbol propuesto.** La bifurcación relevante no es Θ = 0 frente a Θ > 0, sino **si la funcional tiene geometría en su paisaje (B0)**.
- Si no la tiene, a Θ = 0 ninguna inicialización ni dinámica la encuentra.
- A Θ > 0 solo podría aparecer por entropía, y en escala densa el teorema de Chatterjee–Diaconis apunta a ER.

Por eso **B0 va primero**: cuesta minutos y condiciona el valor de todo lo demás.

## VII. Nuevos tests (fichas)

### L-1 — Proposición P-Ω1 (paisaje favorece cliques), test analítico · INCORPORAR ⟶ Rev. 2 (R2.5)
- **Hipótesis:** a N, m y secuencia k-regular fijos, para todo α > 0, η ≥ 0, β, γ, μ ≥ 0, la unión disjunta de K_{k+1} minimiza S. Además, con pesos en [0,1]: t(K3) ≤ t(K2)^{3/2}.
- **Qué intenta demostrar:** que la familia es estructuralmente buscadora de cliques.
- **Qué intenta falsar:** la afirmación implícita de M§13–15 de que S puede favorecer una geometría.
- **Input:** grafos pequeños exhaustivos y construidos.
- **Parámetros:**
  - N ∈ {6, 7, 8}: enumeración exhaustiva de grafos k-regulares para k ∈ {2, 3};
  - N = 216 (6³) y 729 (9³): T³ (k = 6) frente a la unión de K_7;
  - malla (α, η) ∈ {0.1, 1, 10}².
- **Seeds:** ninguna (determinista). Para los pesos en [0,1]: 1000 matrices U(0,1) con PCG64, semilla fija 20261005.
- **Algoritmo:**
  - evaluar `action` con la implementación existente;
  - comprobar la identidad S_smooth = 2Σk² − 12T en binario;
  - comprobar la desigualdad de grafones sobre las matrices aleatorias.
- **Observables:** S, T, S_smooth y t(K3)/t(K2)^{3/2}.
- **Criterio de éxito:** el mínimo exhaustivo es la unión de cliques en el 100% de los casos con k + 1 | N; la identidad se cumple a 1e-12; la razón es ≤ 1 + 1e-12 en el 100%.
- **Criterio de fracaso:** cualquier excepción. P-Ω1 queda refutada y se reabre la rama η.
- **Controles:** K_N y el grafo vacío (casos de igualdad conocidos).
- **Coste estimado:** < 1 min.
- **Dependencias:** ninguna.
- **¿Puede modificar otra decisión?** Sí: B0, η y μ (§V).

### L-2 — Paisaje energético sobre referencias · PROBAR ⟶ Rev. 2 (R2.5)
- **Hipótesis:** en las mismas condiciones (N, m y grado medio coincidentes), ninguna referencia geométrica (RGG3 k12, T³, RGG2 k10) tiene acción menor que la mejor no geométrica (unión de cliques, caveman conectado, ER, uniforme ponderado con la misma Σw), para ningún punto de la malla.
- **Qué intenta demostrar:** que la conclusión de P-Ω1 se extiende a grafos no regulares y reales.
- **Qué intenta falsar:** «existe una región de parámetros de la familia donde la geometría es energéticamente preferida».
- **Input:**
  - referencias de `omega/controls` y `omega/experiments/reference_graphs`;
  - unión de cliques y caveman construidos con el mismo m (resto repartido en la última clique);
  - uniforme ponderado w = m / C(N,2).
- **Parámetros:**
  - N ∈ {216, 800};
  - α̂ ∈ {0.5, 1, 1.5, 2, 3}, γ̂ ∈ {0, 1, 10}, η̂ = η(N−2)/β ∈ {0, 0.1, 1}, μ = 0;
  - la acción se compara a m fijo (μ no altera el orden).
- **Seeds:** 5 por referencia aleatoria (PCG64, `master_entropy = 20261005`).
- **Algoritmo:** generar las referencias, igualar m (exactamente, cuando hay pesos), evaluar `action` y ordenar.
- **Observables:** ΔS = S(mejor geométrica) − S(mejor no geométrica), por arista y en unidades de β; también T, S_deg y S_smooth.
- **Criterio de éxito (P-Ω1 se extiende):** ΔS > 0 en el 100% de los puntos y semillas.
- **Criterio de fracaso:** ΔS < 0 en algún punto, en ≥ 3/5 semillas. Se reabre esa región como PROBAR (dinámica con esos parámetros), con nuevo Consejo.
- **Controles:** N y m deben coincidir exactamente (aserción); la identidad de S_smooth (L-1).
- **Coste estimado:** < 30 min (dominado por RGG3 con N = 800).
- **Dependencias:** L-1.
- **¿Puede modificar otra decisión?** Sí: B0, η, μ y la puerta de entrada de B3.
- **Nota:** solo es condición necesaria a Θ = 0. A Θ > 0 el ΔS por arista se informa para compararlo con diferencias de entropía (no decide).

### O4B-1 — O-04b-lite: cierre explicativo de Ω-B · PROBAR/DIAGNÓSTICO ⟶ Rev. 2 (R2.4, R2.5)
- **Hipótesis:** todos los estados finales de Ω-B en O-04 son (a) uniformes (cv < `uniform_cv_max`) o (b) uniones de cliques KKT con pesos de frontera, cada una saturando Kruskal–Katona ponderado (κ_c = t(K3)/t(K2)^{3/2} ≈ 1 en su subgrafo inducido).
- **Qué intenta demostrar:** que Ω-B falla precisamente por el mecanismo de Kruskal–Katona.
- **Qué intenta falsar:** la explicación misma.
- **Input:** los `.npz` de O-04 existentes (`runs/omega11/o04_ablation`, verificados por sha256 frente a los pasaportes congelados). Si faltan, se reconstruyen con `reconstruct`.
- **Parámetros:** ninguno nuevo. **Seeds:** las de O-04. **N:** 100.
- **Algoritmo:** por estado final:
  1. clasificar según cv y binarización en 1/2 (no en w_min) en `uniforme`, `clique única`, `multi-clique` u `otro`;
  2. calcular κ_c = t(K3)/t(K2)^{3/2} **por componente**, sobre su subgrafo inducido (κ_c = 1 es la saturación; el κ global de r cliques iguales vale r^{-1/2} y no se usa para decidir);
  3. obtener las componentes de la parte de peso > 1/2 y la densidad interna de cada una.
- **Observables:** clase, κ_c por componente (y κ global, solo informe), número y tamaño de las cliques, densidad interna y fracción de peso en halos < 0.1.
- **Criterio de éxito:**
  - 100% de los estados en {uniforme, clique única, multi-clique};
  - toda componente de una mono- o multi-clique con densidad interna ≥ 0.99 y κ_c ≥ 0.9.
- **Criterio de fracaso:** algún estado `otro`, o κ_c < 0.9 en alguna componente de un estado no uniforme. La explicación KK queda debilitada y se reabre el mapa de cuenca.
- **Controles:** K_n ⊕ vacío sintético (κ_c = 1) y ER binario con p = 0.5 (κ ≈ p^{3/2} ≈ 0.35 < 0.9).
- **Coste estimado:** < 30 min.
- **Dependencias:** ninguna.
- **¿Puede modificar otra decisión?** Sí: mapa de cuenca y R7.

### MF-1 — Predicción de campo medio para S0, Θ > 0 · INCORPORAR (expectativa) / DIAGNÓSTICO (informe)
- **Hipótesis (Chatterjee–Diaconis, extensión ponderada):**
  - en equilibrio, las aristas son aproximadamente i.i.d. con densidad p_q(w) ∝ exp(−(w² − 2α̂q²w)/Θ̂) en [0,1];
  - q = E_q[w] minimiza F(q) = −Θ̂·log Z(q) + (4α̂/3)·q³ (energía por arista en unidades de β);
  - comprobaciones: con Θ̂ → 0, F → 1 − 2α̂/3 en q = 1 (reproduce el umbral 3/2); con α̂ = 0, F se reduce a la validación Gibbs ya existente;
  - γ̂ no entra a orden dominante.
- **Qué intenta demostrar:** que S0 a Θ > 0 es un grafo aleatorio de grafón constante (sin estructura).
- **Qué intenta falsar:** «la temperatura induce en S0 una estructura que no es ER».
- **Input:** celdas de O-05 (Metropolis) equilibradas.
- **Parámetros:** la malla de O-05, sin cambios.
- **Seeds:** las de O-05. **N:** 64 y 100.
- **Algoritmo:**
  - cuadratura de Z(q) en [0,1] (scipy, tol 1e-12);
  - todas las raíces de q = E_q[w] y su F;
  - coexistencia si hay dos mínimos locales.
- **Observables:** q* (mínimo global), el conjunto de mínimos locales y ⟨m⟩ medido con su SE.
- **Criterio de éxito:** |⟨m⟩ − q| ≤ 3·SE + 0.02 para algún mínimo local q, en ≥ 90% de las celdas equilibradas, y en particular para q* fuera de las celdas de coexistencia.
- **Criterio de fracaso:** < 90%, o desviaciones que **crecen** de N = 64 a N = 100. El campo medio no describe S0 y se refuerza el interés de PG-1.
- **Controles:** α̂ = 0 (validación Gibbs exacta existente).
- **Coste estimado:** segundos (analítico) más el coste de O-05.
- **Dependencias:** O-05.
- **¿Puede modificar otra decisión?** Sí, el cierre de B1. **No** toca el certificado.

### PG-1 — Pantalla pre-geométrica (no-ER) en O-05 · DIAGNÓSTICO (decide solo el escalado)
- **Hipótesis:** a Θ > 0, S0 no produce estructura distinguible de su nulo de pesos permutados.
- **Qué intenta demostrar:** si existe algo que escalar a N ≥ 800.
- **Qué intenta falsar:** «Θ > 0 induce en S0 correlaciones locales y una brecha espectral pequeña (prerrequisitos de geometría)».
- **Input:** estados de O-05 (finales y muestras post-burn-in separadas ≥ 2τ_int).
- **Parámetros (congelados ahora):**
  - K = 50 permutaciones de los pesos de la triangular superior;
  - z_τ > 5 y g < 0.5;
  - fracción de muestras ≥ 0.8;
  - replicación en N = 64 **y** N = 100.
- **Seeds:** las permutaciones usan un flujo PCG64 derivado de la clave de la celda (`seed_key`, flujo nuevo 4000).
- **N:** 64 y 100.
- **Algoritmo:** por muestra:
  - τ = t(K3)/⟨w⟩³, con z_τ respecto del nulo;
  - g = λ2(L_norm ponderado) observado / mediana de λ2 del nulo.
  - La muestra se marca si z_τ > 5 **y** g < 0.5. La celda se marca si ≥ 80% de sus muestras lo están, en ambos N.
- **Observables:** z_τ, g y la fracción de muestras marcadas.
- **Criterio de éxito (B1 cerrada):** ninguna celda equilibrada marcada.
- **Criterio de fracaso (B1 abierta):** alguna celda marcada con primario modal ∉ {F0, F1, F2}. Se propone O-05b (N ≥ 800) a un nuevo Consejo. **No** es evidencia de geometría (nivel ≤ 1).
- **Controles (deben pasar ANTES de congelar):**
  - ER ponderado i.i.d.: no se marca en ≥ 19/20 semillas;
  - RGG3 k12 con N = 100: se marca en ≥ 19/20;
  - unión de cliques: se marca (y el taxonomista la clasifica F2).
  - **Si un control falla, PG-1 vuelve al Consejo; los umbrales no se reajustan en silencio.**
- **Coste estimado:** < 5% sobre O-05.
- **Dependencias:** O-05.
- **¿Puede modificar otra decisión?** Sí: O-05b, el backend disperso y la retícula decorada. **Nunca** el certificado.

### N-1 — Nulo con igual grado e igual clustering que RGG3 · PROBAR
- **Hipótesis:** el certificado reconoce la geometría y no solo el par (grados, clustering).
- **Qué intenta demostrar:** la especificidad del certificado.
- **Qué intenta falsar:** «un positivo del certificado implica nivel 3 (propiedades compatibles con geometría)».
- **Input:** RGG3 k12, N = 800, 5 semillas.
- **Algoritmo:**
  - double-edge swaps que preservan grados, con aceptación Metropolis sobre |C − C_RGG3|, hasta |ΔC| < 0.01 y ≥ 10·m swaps aceptados;
  - después, `collect_run_evidence` y `assess_run` sin cambios.
- **Seeds:** PCG64, `master_entropy = 20261005`, flujo 4100.
- **Observables:** el veredicto, los 14 campos y los códigos F; D_vol, D_s y D_W.
- **Criterio de éxito:** 5/5 nulos NOT_CANDIDATE (algún campo de corrida falla).
- **Criterio de fracaso:** ≥ 1/5 pasa todos los campos de corrida. El certificado no es específico; se congela cualquier positivo futuro y se vuelve al Consejo.
- **Controles:** el RGG3 original debe pasar (como en O-03); el rewiring sin restricción de clustering debe fallar (como en O-03).
- **Coste estimado:** ~1–2 h.
- **Dependencias:** ninguna.
- **¿Puede modificar otra decisión?** Sí, la credibilidad de todo positivo futuro.

### RP-1 — Reproducibilidad bit a bit de los pasos 1–5 · INCORPORAR
- **Hipótesis:** re-ejecutar p00–O-04 en el commit congelado final reproduce exactamente los arreglos de `167efca`.
- **Qué intenta falsar:** la afirmación de reproducibilidad (M§§ de pasaporte).
- **Input:** los pasaportes congelados (`passports_json.tar.gz`). **Parámetros:** idénticos (`run_steps_1_5.py`).
- **Algoritmo:** comparar los sha256 de cada `.npz` y los `config_hash` por paso.
- **Éxito:** 100% de coincidencia. **Fracaso:** cualquier diferencia, lo que obliga a detener la cadena y auditar.
- **Coste estimado:** ~1.6 h (obligatorio de todos modos por K-1).
- **¿Puede modificar otra decisión?** Sí: si falla, se suspende toda la fase.

## VIII. Criterios congelados

- **Inmutables en esta fase:**
  - todos los umbrales de `omega/config/settings11.py` (`CertificateThresholds`, `LocalStructureConfig`, `TopologyConfig`, `CurvatureConfig`, `WeylConfig`, `DistanceSuiteConfig`, equilibrio R̂ ≤ 1.05 / ESS ≥ 100 / |z| ≤ 3);
  - los 14 campos del certificado y la precedencia F10 > F0 > F1 > F2 > F3 > F4 > F5 > F7 > F6 > F8 > F9;
  - `RULE_D3_NEVER_SUFFICIENT`;
  - D_eff de Ω-1.0;
  - A-17 (anillos);
  - las mallas de s06, O-06 y O-05;
  - `master_entropy = 20240901` para la cadena de la puerta.
- **Se congelan con este documento:** los umbrales de L-1, L-2, O4B-1, MF-1 (3·SE + 0.02; 90%), PG-1 (K = 50; z_τ > 5; g < 0.5; 80%; replicación en dos N) y N-1 (|ΔC| < 0.01; 5/5).
- **Prohibido:**
  - elegir ventanas de ajuste, N o celdas después de ver resultados;
  - re-etiquetar celdas no equilibradas;
  - usar Langevin sin confirmación de Metropolis (A-12);
  - tratar PG-1 o MF-1 como voto de candidatura.
- **Separación de niveles:**
  - PG-1 positiva = nivel 1 («red compleja no-ER»);
  - campo dimensional = nivel 2;
  - certificado completo = nivel 3;
  - **los niveles 4–6 no son alcanzables con ningún experimento de esta fase.**

## IX. Condiciones de falsación (qué nos obliga a abandonar cada rama) ⟶ Rev. 2: la conclusión preregistrada de R8 se sustituye por la frase de R2.2

| Rama | Se abandona si | Se reabre si |
|---|---|---|
| S0, Θ = 0 | Ya abandonada (O-01/O-04 + umbrales) | Nunca dentro de M§13–15; solo con una funcional nueva |
| Ω-B / μ | O4B-1 tiene éxito (ya refutada para P1) | O4B-1 encuentra «otros» |
| η > 0 | L-1 y L-2 tienen éxito con η > 0 | L-2 encuentra ΔS < 0 |
| S0, Θ > 0 | MF-1 tiene éxito **y** PG-1 no marca ninguna celda | Una celda marcada, replicada, con primario ∉ {F0, F1, F2} |
| Certificado como instrumento | N-1 pasa algún nulo, o RP-1 falla | Reparación aprobada por el Consejo |
| Toda la familia M§13–15 | B0 y B1 se cierran a la vez | Ninguna dentro de la familia |

**R8 — expectativa honesta (preregistrada):**
- **¿Qué esperamos si P1 es falsa en esta familia?**
  - L-1 y L-2: cliques en todos los puntos;
  - O4B-1: solo uniformes y cliques con κ_c ≈ 1;
  - O-05: MF-1 se cumple, PG-1 negativa, F0/F1/F3/F4 y cadenas no equilibradas a Θ̂ bajo;
  - s06 y O-06: F0/F1;
  - N-1: 5/5 NOT_CANDIDATE.
  - **Es el resultado más probable.**
- **¿Qué observación inesperada cambiaría la dirección?** Por orden de fuerza:
  1. ΔS < 0 en L-2 (la geometría es preferida en algún punto), lo que reabre η/μ;
  2. una celda PG-1 replicada y no clique a Θ > 0, lo que lleva a O-05b;
  3. desviaciones de MF-1 que crecen con N;
  4. un estado «otro» en O4B-1.
- **Si ninguna ocurre, la conclusión preregistrada es:** «La familia M§13–15 no puede producir una fase geométrica: a Θ = 0 su paisaje favorece cliques, y a Θ > 0 (escala densa) su equilibrio es de grafón constante. **No implementar nuevas dinámicas dentro de esta familia; el siguiente intento legítimo requiere una nueva hipótesis dinámica que pase L-2 antes de simularse.**»

## X. Coste computacional (orden por valor científico / coste)

| # | Experimento | Valor | Coste | Valor/coste |
|---|---|---|---|---|
| 1 | L-1 | Muy alto (decide B0) | < 1 min | ★★★★★ |
| 2 | L-2 | Muy alto (B0, η, μ, puerta de B3) | < 30 min | ★★★★★ |
| 3 | O4B-1 | Medio (cierre explicativo) | < 30 min | ★★★★ |
| 4 | MF-1 (analítico) | Alto (expectativa de O-05) | segundos | ★★★★ |
| 5 | N-1 | Medio (credibilidad de positivos futuros) | 1–2 h | ★★★ |
| 6 | RP-1 + s06 | Medio (reproducibilidad) / bajo (s06) | ~2.2 h | ★★ (obligatorio por la puerta) |
| 7 | O-06 | Bajo-medio (calibración de RGG3) | 1–3 h | ★★ (obligatorio por la puerta) |
| 8 | O-05 + PG-1 | Medio (prior nulo; única prueba empírica de B1) | 5–20 h (piloto previo) | ★★ |

**Regla de parada (economía).** Los ítems 1–4 se ejecutan primero, fuera de la puerta. Si L-1 o L-2 refutan P-Ω1, el Consejo se reúne **antes** de gastar las 8–25 h de los ítems 6–8, porque la prioridad cambiaría hacia η.

## XI. Plan de implementación (solo tras la ratificación; commits pequeños) ⟶ Rev. 2 (R2.3, R2.7): el bloque L va primero

| Commit | Contenido | Etiqueta | Validación antes del push |
|---|---|---|---|
| C0 | Este documento | — | — |
| C1 | `AMENDMENTS.md`: ratificación R1–R8 (D-1), regla ρ ≤ w_min (D-5), registro A-18 (PG-1) y A-19 (MF-1) como enmiendas previas a O-05 | INCORPORAR (docs) | — |
| C2 | `tests/analytical/test_landscape_cliques.py` (L-1) | INCORPORAR | pytest, mypy --strict, test de arquitectura (RNG PCG64) |
| C3 | `tools/landscape.py` (L-2) + ejecución + `results/landscape/` | PROBAR | árbol limpio al ejecutar; `summary.json` preregistrado → final |
| C4 | `tools/o04b_lite.py` (O4B-1) + ejecución + `results/o04b_lite/` | PROBAR | sha256 de los `.npz` frente a los pasaportes |
| — | **Punto de control del Consejo** si L-1, L-2 u O4B-1 fallan | — | — |
| C5 | `omega/statistics/mean_field.py` (MF-1) + tests (límites Θ̂ → 0 y α̂ = 0) | INCORPORAR | pytest y mypy |
| C6 | `omega/statistics/pregeometric.py` (PG-1) + tests de controles (ER, RGG3 y cliques, 20 semillas) | DIAGNÓSTICO | **los controles deben pasar; si no, se vuelve al Consejo** |
| C7 | Integración de MF-1 y PG-1 en el informe de O-05 (sin tocar `assess_run` ni el certificado) + smoke | DIAGNÓSTICO | suite rápida completa |
| C8 | `omega/controls/clustered_rewiring.py` + tests + `tools/n1_null.py` (N-1) | PROBAR | pytest y mypy |
| C9 | Ejecución de N-1 + `results/n1_null/` | PROBAR | — |
| C10 | **Congelación:** config_hash de s06, O-06 y O-05 en `AMENDMENTS.md`; `tools/run_chain.py` (p00 → O-04 → s06 → O-06 → piloto → O-05) | INCORPORAR | suite completa; árbol limpio |
| — | **Ejecución de la cadena sobre C10 (sin commits durante la corrida, K-1)**; RP-1 al terminar O-04 | — | — |
| C11 | Resultados de la cadena + informe de RP-1 | — | — |
| C12 | Análisis de la fase y cierre (o no) de B0, B1 y B2 | — | — |

**Versión:** todo lo anterior queda en **Ω-1.1**. Ninguna pieza cambia la dinámica, la definición de estado, los estimadores ni los criterios del certificado. Ningún experimento anterior deja de ser comparable.

---

## XII. Resoluciones que requieren firma del usuario ⟶ ratificadas (R2.0)

| # | Resolución del Consejo | Estado |
|---|---|---|
| D-1 | R1 APPROVED WITH MODIFICATION; R3 APPROVED; R4 APPROVED (veto necesario, K-2); R8 APPROVED; R2, R5 y R6 APPROVED; R7 DEFERRED | propuesta |
| D-2 | s06 APPROVED WITH MODIFICATION (cadena completa en el commit congelado + RP-1) | propuesta |
| D-3 | O-04b APPROVED WITH MODIFICATION (O4B-1 lite); mapa de cuenca DEFERRED | propuesta |
| D-4 | O-06 → O-05 en el orden de la puerta; O-05 APPROVED WITH MODIFICATION (+PG-1, +MF-1, piloto, malla intacta) | propuesta |
| D-5 | Ω-1.2 REJECTED; regla de lectura preregistrada | propuesta |
| A | O-04b → como D-3 | propuesta |
| B | O-05 → como D-4; definición operativa de «fase geométrica» sin cambios (celda equilibrada con Metropolis y ≥ 80% de estados que pasan todos los campos); PG-1 solo decide el escalado | propuesta |
| C | O-06 sin cambios; solo su parte RGG3 importa para un eventual O-05b; no modifica ningún criterio de O-05 | propuesta |
| D | L-1 y L-2 APPROVED; N-1 APPROVED; «D ≈ 3 artificial» REJECTED; retícula decorada DEFERRED | propuesta |
| E | Árbol §VI con B0 (paisaje) como bifurcación principal | propuesta |

Nada se implementa hasta recibir la ratificación (total o por partes).

## XIII. Nota de cierre de la auditoría
Suite rápida: 598 passed, 28 deselected, exit 0 (317 s), sobre `539770f`. No se escribió código en esta fase.
