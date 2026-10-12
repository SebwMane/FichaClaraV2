# Consejo Ω, Revisión 3: hipótesis candidata Ω-C (orden causal y acción discreta)

- **Fecha:** 2026-10-11.
- **Rama:** `claude/hopeful-galileo-88u1l1`, base `ec448e4`.
- **Naturaleza:** documento de decisión. Separa lo que se ejecuta ahora (bloque C0, validación de instrumentos, sin dinámica) de lo que **requiere ratificación del usuario** (la hipótesis dinámica Ω-C y su compuerta).
- **Versión del simulador:** Ω-1.1 sigue intacto. Nada de este documento toca S0, el certificado, D_eff ni ningún umbral de `settings11.py`.
- **Principio rector:** «No queremos demostrar que Ω es correcta; queremos descubrir qué puede producir Ω y dónde los experimentos pueden demostrar que estamos equivocados.»

## 1. Por qué existe este documento

El bloque L terminó en la rama «SÍ» (`docs/OMEGA_BLOQUE_L_RESULTADOS.md` §3). La consecuencia preregistrada fue:

> «No implementar nuevas dinámicas dentro de esta familia; el siguiente intento legítimo requiere una nueva hipótesis dinámica que pase L-2 antes de simularse.»

La hipótesis nueva debía venir de un documento del usuario. En la conversación del 2026-10-11 el usuario y el orquestador revisaron la literatura de espacio-tiempo emergente y el usuario aprobó seguir la recomendación de trabajar con **orden causal y una acción discreta**. Este documento la formaliza como candidata. No la da por aprobada.

## 2. Diagnóstico del bloque L que motiva el cambio

Nivel 2 (explicación propuesta, ya con soporte computacional): la recompensa por triángulo de S0 exige codegrado extensivo, por lo que apaga o colapsa en cliques cualquier geometría de grado acotado. Dos lecciones independientes de la literatura apuntan en la misma dirección:

1. **Un funcional que solo premia densidad local de relaciones no tiene por qué seleccionar variedades.** En conjuntos causales, el número de órdenes de N elementos crece como e^{N²/4} (Kleitman y Rothschild) y casi todos son de tres capas, sin geometría. La selección tiene que venir de una acción construida para que su límite continuo sea Einstein-Hilbert, no de un premio genérico.
2. **Existe una acción con esa propiedad.** La acción de Benincasa-Dowker (BD) se calcula solo contando intervalos de orden y, en promedio sobre sprinklings, tiende a la acción de Einstein-Hilbert más términos de borde (Benincasa y Dowker 2010; Dowker y Glaser 2013). En 2D, su media sobre una región causalmente convexa de Minkowski dividida con líneas nulas es la suma alternada de vértices, aristas y caras del teselado (Benincasa, Dowker y Schmitzer 2011), es decir, 1 para un intervalo causal.

## 3. Hipótesis candidata Ω-C (REQUIERE RATIFICACIÓN)

- **Estado:** Ω-C = (V, ≺), un orden parcial finito sobre V. No hay coordenadas, pesos ni métrica. El orden no es geometría: no define distancias ni dimensión por sí mismo (compatible con M§2 y M§6).
- **Acción:** S_BD en la dimensión de prueba (primero d = 2, por coste y porque allí existen resultados exactos).
- **Pregunta de bifurcación:** ¿la acción BD distingue órdenes con aproximación de variedad de los órdenes no geométricos dominantes por entropía, y en qué sentido?
- **Advertencia preregistrada (lección de L-2):** en el bloque L el criterio fue «menor acción gana». Para Ω-C ese criterio **no se traslada directamente**. La suma cuántica usa el peso e^{iS} (lorentziano), no e^{-S}. Un orden con acción muy negativa no gana por eso; puede quedar suprimido por interferencia. La compuerta L-2c debe definirse para el peso que se vaya a usar, y esa elección es del usuario (§6).

## 4. Bloque C0: validación de instrumentos causales (se ejecuta ahora)

Equivalente a M§22 (Experimento 0): si los instrumentos no recuperan lo conocido, nada de lo que midan sobre Ω-C es confiable. **C0 no simula ninguna dinámica ni dice nada sobre Ω.**

Universo de entradas:
- **Sprinklings:** N puntos uniformes en un diamante causal de Minkowski de dimensión d ∈ {2, 3, 4}. La relación se calcula una vez desde las coordenadas y después **las coordenadas se descartan**: todo estimador recibe solo la matriz de orden. Las coordenadas solo sirven como referencia en C0-2.
- **Controles no geométricos:** órdenes de Kleitman-Rothschild (KR, tres capas N/4, N/2, N/4, relaciones entre capas consecutivas con p = 1/2, cierre transitivo) y percolación transitiva (aristas i<j con probabilidad p, cierre transitivo).

### C0-1 Dimensión de Myrheim-Meyer · PROBAR
- Estimador: fracción de pares relacionados r = R / C(N,2), invertida con r(d) = Γ(d+1)Γ(d/2) / (2Γ(3d/2)).
- Malla: d ∈ {2, 3, 4}, N ∈ {256, 1024}, 20 semillas por celda.
- **Éxito:** en N = 1024, |media(d_MM) − d| ≤ 0.10 para los tres d. Justificación: la esperanza de r es exactamente r(d) para N fijo (sin sesgo de tamaño), así que la tolerancia solo cubre varianza.
- **Fracaso:** cualquier d fuera de tolerancia. Se audita el estimador antes de seguir.

### C0-2 Número = volumen · PROBAR
- Para 200 pares relacionados (p, q) por semilla, elegidos al azar, se compara el número de elementos estrictamente dentro de I(p, q), contado **solo con el orden**, con el valor esperado (N − 2)·V(I)/V(diamante), donde V(I) = c_d τ^d se calcula con coordenadas.
- Malla: d ∈ {2, 3, 4}, N = 1024, 20 semillas.
- **Éxito:** la razón agregada Σ conteos / Σ esperados está en [0.98, 1.02] para los tres d.
- Es la versión discreta de ∫ d^dx √(−g) ≈ N/ρ.

### C0-3 Acción BD 2D sobre el diamante plano · PROBAR
- S = N − 2N₀ + 4N₁ − 2N₂, donde N_k es el número de pares relacionados con exactamente k elementos estrictamente entre ellos (N₀ son los links).
- Malla: d = 2, N ∈ {64, 256, 1024}, 200 semillas por N.
- **Valor de literatura:** 1 (suma alternada V − E + F del intervalo sin dividir).
- **Éxito:** en N = 1024, |media(S) − 1| ≤ 3·SE.
- **Fracaso:** se audita primero la fórmula y la definición de N_k. No es una falsación de nada físico.

### C0-4 Discriminación de órdenes no geométricos · PROBAR
- Segundo estimador independiente: dimensión de punto medio. En el intervalo más grande I(a, b) se busca m que maximiza min(|I(a,m)|, |I(m,b)|) y d_mid = log₂(|I(a,b)| / min).
- Bandera «compatible con variedad»: |d_MM − d_mid| ≤ τ.
- **Éxito:** sprinklings d = 2, 3, 4 con N = 1024 marcados en ≥ 90% de las semillas, **y** KR y percolación transitiva (N = 1024) no marcados en ≥ 90%.
- **τ queda PENDIENTE de un piloto declarado** (el sesgo de tamaño finito de d_mid no tiene valor cerrado conocido). El piloto usa otra entropía y no cuenta como evidencia. τ se congela en una enmienda registrada antes de la corrida oficial.

### C0-D Acción BD de los controles · DIAGNÓSTICO (sin voto)
- Se registra S_BD de KR y de la percolación con N = 1024. Sirve para diseñar L-2c y para comprobar la advertencia de §3.

### Semillas y reproducibilidad
- `PCG64(SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, celda, semilla)))`, como en `omega.config.seeds`.
- Entropía oficial: 20261011. Entropía del piloto: 20261012 (distinta, declarada).
- `summary.json` en fase «preregistered» antes de calcular y «final» después. Modo full exige árbol git limpio.

## 5. Compuerta de Ω-C (definida aquí, ejecución sujeta a ratificación)

- **L-2c (paisaje):** comparar S_BD de órdenes con aproximación de variedad frente a KR y percolación a N fijo, bajo el peso que el usuario elija. Sin dinámica.
- **L-3c (estabilidad):** partir de un sprinkling 2D y aplicar una dinámica Metropolis sobre órdenes 2D con el peso elegido; medir si la bandera de C0-4 persiste.
- Solo si L-2c y L-3c no excluyen la geometría se discute simular desde órdenes aleatorios.

## 6. Resoluciones que requieren firma del usuario

| # | Resolución | Estado |
|---|---|---|
| R3-1 | Adoptar Ω-C como hipótesis candidata (no como hipótesis aceptada) | propuesta |
| R3-2 | Peso de la compuerta: (a) euclídeo e^{-βS} como en la literatura de MCMC 2D, (b) lorentziano e^{iS} vía medidas de interferencia, o (c) ambos con criterios separados | propuesta |
| R3-3 | Mantener d = 2 hasta que L-2c y L-3c estén cerrados | propuesta |
| R3-4 | C0 se ejecuta sin esperar R3-1 a R3-3, porque solo valida instrumentos | **aprobada por el usuario en chat (2026-10-11)** |

## 7. Escala de afirmación aplicable a C0

- Un éxito de C0 es **nivel 1 sobre instrumentos**: «los estimadores recuperan valores conocidos en sprinklings y detectan los controles». No dice nada sobre si Ω-C produce geometría.
- Un fracaso de C0 bloquea toda la rama Ω-C hasta reparar el instrumento.

## 8. Enmiendas operativas del bloque C0 (registradas ANTES de la corrida oficial)

Origen: auditoría del orquestador (`docs/handoff/AUDITORIA_C0.md`), el smoke (N ≤ 128, entropía 20261013) y el piloto declarado de C0-4 (N = 1024, 10 semillas, entropía 20261012). Ningún resultado oficial existía al registrarlas.

| ID | Ítem | Cambio | Motivo | ¿Puede crear un positivo? |
|---|---|---|---|---|
| C-A1 | C0-4, control de percolación | p = 0.05 queda congelado. | El handoff no fijaba p; lo eligió el agente de código y solo se había visto en el smoke con N = 128. | Sí en principio; se declara. |
| C-A2 | C0-3 | Se añade un segundo criterio con voto: las medias de N₀, N₁ y N₂ en N = 1024 deben estar a ≤ 3·SE de su esperanza exacta para 2-órdenes aleatorios (`omega/causality/exact.py`). | La auditoría midió sd(S) ≈ 322 en N = 1024, de modo que 3·SE ≈ 68 y C0-3 casi no tiene potencia. La esperanza exacta (integral independiente del código) da E[S] = 1 para todo N, lo que confirma la fórmula y el valor de literatura. | No: es más estricto. |
| C-A3 | C0-4, valores indefinidos | Si d_mid o d_MM es nan, la bandera es «no compatible con variedad». | En KR el punto medio de cualquier intervalo máximo está a un link de un extremo, por lo que min(|I(a,m)|, |I(m,b)|) = 0 y d_mid no existe. En un intervalo con aproximación de variedad siempre existe un punto medio con mitades de tamaño extensivo, así que la ausencia es una señal legítima, no un artefacto. | Sí: hace que KR no quede marcado. Se declara; vale en ambos sentidos (un sprinkling con nan también cuenta en contra). |
| C-A4 | C0-4, τ | τ = 0.75. Regla: 1.1 × el mayor |d_MM − d_mid| de los sprinklings del piloto (0.677, en d = 4), redondeado hacia arriba a 0.05. | Fijar τ era el único objetivo del piloto. | Ver nota. |

**Nota preregistrada sobre C0-4 (expectativa honesta).** En el piloto, la percolación con p = 0.05 dio d_MM ≈ 1.21 y d_mid ≈ 1.14 en las 10 semillas, con gaps ≤ 0.10. Quedará marcada como «compatible con variedad» para cualquier τ ≥ 0.10, así que **se espera que C0-4 FRACASE en la corrida oficial**. El criterio no se modifica para evitarlo. La lectura preregistrada de ese fracaso es: «la coherencia entre dos estimadores de dimensión no basta para distinguir un orden no geométrico de uno con aproximación de variedad; la percolación transitiva densa se comporta como un objeto coherente de dimensión ≈ 1.2». Es el mismo aviso de M§11 (un D_eff por sí solo no demuestra geometría) aplicado a órdenes. El rediseño (C0-4b) se preregistra en una ronda nueva con entropía nueva.

## 9. C0-4b: discriminador por perfil de abundancias · PROBAR (preregistrado 2026-10-12, sin ratificación necesaria)

**Por qué.** C0-4 fracasó porque dos estimadores de dimensión pueden coincidir en un orden no geométrico. C0-4b exige algo más fuerte: que el orden tenga el **perfil completo de tamaños de intervalo** que tendría un sprinkling de su propia dimensión estimada. Es una prueba de forma, no de un número.

**Declaración de transparencia.** Este diseño se hizo después de ver que la percolación con p = 0.05 engañó a C0-4. Para no ajustar el criterio a ese control: (a) δ se calibra **solo con sprinklings**, nunca con controles; (b) se añaden dos controles que nadie ha visto todavía, percolación con p = 0.01 y un orden de capas gruesas; (c) todas las constantes quedan fijadas aquí antes de escribir código.

**Definición (congelada).**
1. d̂ = d_MM del orden. Si es nan, la bandera es «no compatible».
2. d_ref = d̂ redondeado al entero más cercano y recortado a {1, 2, 3, 4}.
3. Perfil de referencia E_ref[N_k] para k = 0..8 al mismo N:
   - d_ref = 1: cadena exacta, N_k = N − 1 − k.
   - d_ref = 2: valor exacto de `omega/causality/exact.py`.
   - d_ref = 3, 4: media de 40 sprinklings de referencia (entropía 20261015).
4. Distancia D = max_{k ≤ 8} |ln((N_k + 1) / (E_ref[N_k] + 1))|.
5. **Calibración de δ (solo positivos):** para d = 2, 3, 4 se generan 40 sprinklings de calibración (entropía 20261016) y δ_d = 1.5 × percentil 95 de su D. Para d_ref = 1 se usa δ_1 = δ_2.
6. Bandera «compatible con variedad» ⇔ d̂ finito **y** D ≤ δ_{d_ref}.

**Corrida oficial** (entropía 20261017, N = 1024, 20 semillas por celda):
- Positivos: sprinklings d = 2, 3, 4.
- Controles: KR; percolación p = 0.05; percolación p = 0.01; capas gruesas (8 capas de N/8, relaciones entre capas consecutivas con p = 0.5, cierre transitivo).

**Éxito:** cada positivo marcado en ≥ 90% de las semillas **y** cada control marcado en ≤ 10%.
**Fracaso:** cualquier celda fuera. Se vuelve al Consejo; no se recalibra δ ni se cambia K.
**Escala de afirmación:** un éxito es nivel 1 sobre instrumentos («el clasificador separa estos positivos de estos controles a N = 1024»). No garantiza que separe cualquier orden no geométrico.

## 10. Ratificación R3-2 y preregistro de L-2c (2026-10-12)

### 10.1 Estado de las resoluciones

| # | Estado |
|---|---|
| R3-1 | **Pendiente.** Ω-C sigue siendo candidata hasta que la compuerta (L-2c, L-3c) dé su resultado. |
| R3-2 | **RATIFICADA por el usuario (2026-10-12): opción (c), ambos pesos con criterios separados.** Cada rama tiene su propio veredicto; el resultado de una no puede compensar el de la otra. |
| R3-3 | **No firmada.** Se usa d = 2 porque es el valor por defecto de §3 y el único con resultados exactos; el usuario puede cambiarlo. |

La ratificación de R3-2 autoriza la compuerta definida en §5 (L-2c ahora; L-3c después de leer L-2c). No autoriza simular Ω-C desde órdenes aleatorios fuera de la compuerta.

### 10.2 Hechos estructurales usados (nivel 1, prueba corta)

- **Paridad.** S = N − 2N₀ + 4N₁ − 2N₂ implica S ≡ N (mod 2) para todo orden. Por eso e^{iλS} no se cancela en λ = π para ninguna clase.
- **Bipartito completo aleatorio** (n/2 + n/2 elementos, cada par entre capas relacionado con p = 1/2): todas las relaciones son links, S = N − 2·Bin(n²/4, 1/2) y |E e^{iλS}| = |cos λ|^{n²/4} exactamente.
- **KR de tres capas.** Igual que el bipartito, salvo los pares capa 1–capa 3 con exactamente 1 o 2 vecinos comunes, cuya probabilidad es ≤ C·N^4·(3/4)^{N/2}, despreciable en N = 1024. Con ~2^{N²/4} órdenes en la clase, su peso total es ≈ |2 cos λ|^{N²/4}: **suprimido si y solo si |cos λ| < 1/2, es decir λ mod π ∈ (π/3, 2π/3)**. Nivel 1 para el bipartito; nivel 2 (con cota explícita) para KR.

### 10.3 L-2c-E (rama euclídea, e^{−βS} restringido a 2-órdenes) · PROBAR

- **Pregunta:** en el límite β → ∞ (el análogo exacto de L-2: «gana la menor acción»), ¿es mínima la acción de los 2-órdenes geométricos?
- **Geométricos:** 2-órdenes aleatorios marcados por C0-4b (N = 1024, 20 semillas, entropía 20261018).
- **No geométricos (2-órdenes construidos):** k antichains de tamaño N/k apiladas y totalmente relacionadas, k ∈ {2, 4, 16, 128}, etiquetas permutadas. Validación previa: ninguna debe quedar marcada por C0-4b; si alguna queda marcada, se excluye y se informa.
- **Estadístico:** ΔS_E = min S(no geométricas) − media S(geométricas).
- **Decisión:** ΔS_E > 0 → «geometría preferida a β → ∞»; ΔS_E < 0 → «geometría no preferida a β → ∞».
- **Expectativa honesta preregistrada:** «no preferida». La bicapa da S = N − N²/2 ≈ −5·10⁵ frente a S ≈ 1. Coincide con la fase «cristalina» de la literatura de MCMC 2D.
- **Lectura preregistrada (distinta de L-2):** un «no preferida» **no cierra** la rama euclídea, porque en β = 0 el 2-orden típico es geométrico (C0-4b: 100%). La pregunta pasa a ser si existe una ventana 0 < β < β_c con fase geométrica, y eso solo lo responde L-3c-E (Metropolis sobre 2-órdenes). Un «preferida» cerraría la compuerta E con éxito sin necesidad de L-3c-E.

### 10.4 L-2c-L (rama lorentziana, e^{iλS} sobre todos los órdenes) · PROBAR

- **L-i, instrumento (vota):** el código reproduce |cos λ|^{n²/4} para el bipartito con n ∈ {8, 12}, λ ∈ {0.2, 0.5, 1.0, π/2}, 20 000 muestras: |φ̂ − fórmula| ≤ 3·SE en todos los puntos. Si falla, se audita el código.
- **L-ii, ventana de supresión de KR:** analítica, §10.2. Ventana W = (π/3, 2π/3).
- **L-iii, coherencia de la clase geométrica (vota):** para 2-órdenes aleatorios con N ∈ {16, 32, 64, 128, 256}, 20 000 muestras por N (entropía 20261019), se estima |φ₂(λ)| = |media de e^{iλS}| en λ ∈ {1.2, π/2, 1.9} (dentro de W). Solo se usan puntos con |φ̂₂| > 3/√M. Se ajusta −ln|φ₂(λ)| ≈ c·N^b por mínimos cuadrados en log-log.
  - **PASA:** en al menos uno de los tres λ hay ≥ 3 puntos utilizables y el ajuste da b < 2 con su extremo superior (b + 2·SE) < 2. Lectura: dentro de W, la cancelación de fase de la clase geométrica crece más lento que la de KR (∝ N²), así que la clase geométrica **domina a KR** asintóticamente.
  - **NO PASA:** b + 2·SE ≥ 2 en todos los λ utilizables.
  - **INDETERMINADO:** ningún λ con 3 puntos utilizables.
- **Declaración de transparencia:** los valores de S de C0-3 (N = 64, 256, 1024; 200 semillas) existen y habrían permitido estimar φ₂; **no se usaron**. Antes de este preregistro solo se comprobó la paridad de S (§10.2).
- **Escala de afirmación:** L-2c-L compara la clase geométrica solo con KR y con la bicapa (que es un bipartito completo, ver §10.2). Un PASA es nivel 2 para N ≤ 256 y **nivel 4 (extrapolación) para N grande**. No dice nada sobre otras clases no geométricas.

### 10.5 Enmiendas operativas de L-2c (registradas ANTES de la corrida oficial)

| ID | Ítem | Cambio | Motivo | ¿Puede crear un positivo? |
|---|---|---|---|---|
| L2c-A1 | L-iii | Un punto con φ̂ ≥ 1 − 10⁻¹² no es utilizable. | −ln φ̂ = 0 no admite el ajuste log-log. Con la dispersión observada de S (57 valores distintos ya en N = 16) es prácticamente imposible. | En teoría sí; en la práctica es inerte. |
| L2c-A2 | L-i | Entropía 20261020. | §10 no la fijaba; la eligió el agente de código. | No: L-i solo valida el instrumento. |
| L2c-A3 | L-2c-E | Un ΔS_E indefinido se lee como INDETERMINADO. | El código lo leía como «no preferida». | No: elimina una lectura falsa. |

## 11. Ratificación R3-4, sondeo de la acción suavizada y preregistro de C1 y L-3c-E (2026-10-12)

### 11.1 Resolución

| # | Estado |
|---|---|
| R3-4 | **RATIFICADA por el usuario (2026-10-12):** la rama lorentziana de Ω-C usa la acción BD suavizada con escala de no localidad ε. Es un parámetro nuevo de la hipótesis. La rama euclídea conserva la acción sin suavizar, tal como se preregistró L-2c-E, para no cambiar la definición de la rama a mitad de la compuerta. |

### 11.2 Acción suavizada (definición congelada)

S_ε = ε · ( N − 2ε · Σ_n N_n f(n, ε) ),  f(n, ε) = (1−ε)^n [1 − 2εn/(1−ε) + ε² n(n−1) / (2(1−ε)²)].

Normalizada para que ε → 1 dé exactamente S = N − 2N₀ + 4N₁ − 2N₂ (la literatura usa un factor 2 global distinto). Para 2-órdenes aleatorios, con x = ab y M = N − 2:

E[S_ε] = ε ( N − 2ε · C(N,2) · ½ ∫∫ 4(1−a)(1−b) [ (1−εx)^M − 2εMx(1−εx)^{M−1} + ε² M(M−1) x² (1−εx)^{M−2} / 2 ] da db ).

### 11.3 Sondeo declarado (no es evidencia)

Entropía 31337, 150 2-órdenes por celda, N ∈ {64, 128, 256}, ε ∈ {0.5, 0.25, 0.1}. Se registró la media y la dispersión de S_ε y de S_ε/ε². Hallazgo que motiva §11.5: sd(S_ε) baja mucho con ε (0.90 en N = 256, ε = 0.1), pero sd(S_ε/ε²) queda entre 29 y 90, igual que la acción sin suavizar.

### 11.4 C1: validación del instrumento suavizado · PROBAR

- **C1-1:** S_ε con ε = 1 − 10⁻⁹ coincide con S en ≤ 10⁻³ (test unitario).
- **C1-2:** para N ∈ {64, 256} y ε ∈ {0.5, 0.25, 0.1}, 200 2-órdenes (entropía 20261021): |media(S_ε) − E[S_ε]| ≤ 3·SE en las seis celdas.
- **Fracaso:** se audita fórmula y código; la rama L no usa S_ε hasta que pase.

### 11.5 Rama L: lectura estructural (nivel 2; no se ejecuta prueba nueva)

Con S_ε, los links de KR aportan −2ε² cada uno, así que KR queda suprimido sí y solo sí λε² mod π ∈ (π/3, 2π/3). En esa misma ventana la fase de la clase geométrica es λS_ε = (λε²)(S_ε/ε²), con sd(S_ε/ε²) ≫ 1, así que también se cancela. Comparar ambas cancelaciones exige medir promedios de fase exponencialmente pequeños: el **problema de signo**. El muestreo Monte Carlo no puede decidir la rama L. Decidirla exige un método analítico (como el de Loomis y Carlip para KR) o una continuación analítica, que es otra hipótesis. **Estado de la rama L: bloqueada por método, no por resultado.**

### 11.6 L-3c-E: estabilidad de la geometría bajo e^{−βS} en 2-órdenes · PROBAR

- **Dinámica:** Metropolis sobre 2-órdenes con u fijo y v una permutación; movimiento = transponer dos valores de v elegidos al azar; aceptación min(1, e^{−βΔS}) con la acción sin suavizar.
- **Inicio:** un 2-orden aleatorio (geométrico). Es una prueba de estabilidad, como L-3b.
- **Malla:** N ∈ {32, 64}; β = 0.025·(√2)^k, k = 0..12 (0.025 a 1.6); una cadena por celda; 1000 barridos de termalización (1 barrido = N movimientos) y 1000 de medición, una muestra cada 50 barridos (20 muestras). Entropía 20261022.
- **Clasificador:** bandera C0-4b recalibrada a cada N con la regla de §9 (δ = 1.5 × p95 de 40 2-órdenes aleatorios, entropía 20261023; referencia exacta d = 2; d_ref = 1 usa δ₂).
- **β_c(N):** el menor β de la malla con fracción marcada < 0.5. Si no existe, β_c > 1.6.
- **Decisión:**
  - **Ventana estable:** β_c(64) ≥ β_c(32).
  - **Ventana que se angosta:** β_c(64) < β_c(32).
- **Expectativa honesta (nivel 4):** un argumento de conteo (bicapa con S ≈ −N²/2 frente a ~2N ln N de entropía de 2-órdenes) sugiere β_c ∝ ln N / N, es decir, que se angosta.
- **Lectura preregistrada de «se angosta»:** con β fijo, la geometría no sobrevive al límite N → ∞. Mantener la rama E exigiría un acoplamiento que escale con N, lo cual es un cambio de hipótesis que requiere firma del usuario. **No cierra la rama por sí solo.**
- **Escala:** nivel 1 para N ∈ {32, 64}; cualquier afirmación sobre N grande es nivel 4.

### 11.7 Enmienda operativa de L-3c-E (registrada ANTES de la corrida oficial)

| ID | Ítem | Cambio | Motivo | ¿Puede crear un positivo? |
|---|---|---|---|---|
| L3cE-A1 | Clasificador a N = 32 y 64 | Antes de leer β_c, el clasificador recalibrado debe marcar ≥ 90% de 100 2-órdenes aleatorios y ≤ 10% de 20 2-órdenes en capas para cada k ∈ {2, 4, 8} (entropía 20261024). Si falla en algún N, L-3c-E es INDETERMINADO. | §11.6 recalibraba δ a N pequeño pero no exigía comprobar que siguiera separando geometría de capas (el δ₂ de N = 1024 no vale a N chico, ver A-13). En el scratchpad de auditoría (entropía 5150) sí separa: 98–99% frente a 0%. | No: solo puede bloquear. |

## 12. Ratificación R3-5 y preregistro de L-3c-E′ (2026-10-12)

| # | Estado |
|---|---|
| R3-5 | **RATIFICADA por el usuario (2026-10-12):** en la rama euclídea el acoplamiento es β = b/N con b fijo. Es un supuesto nuevo de la hipótesis que queda por justificar físicamente. |

**Pregunta.** Con β = b/N, ¿la ventana geométrica deja de angostarse al crecer N?

- **Dinámica, movimiento, inicio y clasificador:** idénticos a §11.6, con la enmienda L3cE-A1 (validación del clasificador a cada N, entropía 20261027).
- **Malla:** N ∈ {32, 64, 128}; b = 2.0 · 2^{k/4}, k = 0..8 (2.0 a 8.0, pasos de ≈ 19%); β = b/N; 1000 barridos de termalización y 1000 de medición, una muestra cada 50. Entropías: cadena 20261025, calibración 20261026, validación 20261027.
- **b_c(N):** el menor b de la malla con fracción marcada < 0.5; si no existe, b_c > 8.
- **Decisión:**
  - **PASA (ventana estable con β = b/N):** b_c(64) ≥ b_c(32)·2^{−1/4} y b_c(128) ≥ b_c(64)·2^{−1/4}. Es decir, b_c no cae más de un paso al duplicar N.
  - **FALLA:** b_c cae más de un paso en alguna duplicación.
  - **INDETERMINADO:** el clasificador no se valida en algún N, o b_c > 8 en algún N.
- **Declaraciones.** (a) b_c(32) y b_c(64) se pueden anticipar por L-3c-E (β_c·N ≈ 4.5, con resolución de un paso de √2); la información nueva es N = 128 y la malla más fina. (b) **Potencia limitada:** un escalamiento β_c ∝ ln N / N daría b_c creciendo ≈ 17% por duplicación, menos de un paso; la prueba no distingue entre «1/N» y «ln N / N». Solo detecta un angostamiento más rápido que 1/N. (c) La parte con N = 128 se ejecuta en cinco tramos por el límite de tiempo del entorno; `final` verifica un único commit.

## 13. Bloque A: vía analítica para la rama lorentziana (preregistrado, sin ejecutar)

**Método de referencia.** Loomis y Carlip (2017) y sus extensiones (Mathur, Singh y Surya 2020; P. Carlip, S. Carlip y Surya 2022): en vez de muestrear, se cuenta exactamente cuántos órdenes de una familia tienen cada valor de acción y se suma la fase de forma cerrada. En órdenes por capas la acción depende casi solo del número de links, y la suma se factoriza par por par.

**A-1 (derivación, nivel 1 con cota de error).** Para cualquier familia por capas (bicapa, KR, k capas) con la acción suavizada S_ε, los pares entre capas no consecutivas tienen n ≳ N/k elementos entre ellos y f(n, ε) ≤ C n² (1−ε)^n, despreciable a ε fijo. Sumando sobre todos los grafos de relaciones entre capas consecutivas, cada par aporta (1 + e^{−2iλε²}), de modo que el peso de la familia es ≈ |2 cos(λε²)|^{P}, con P el número de pares entre capas consecutivas (∝ N²). **Supresión universal de las familias por capas ⇔ λε² mod π ∈ (π/3, 2π/3).** Hay que escribir la cota explícita del término despreciado.

**A-2 (verificación numérica, vota).** Para bicapa y KR con N ≤ 16 y λε² pequeño (donde la fase es medible), el promedio de fase de S_ε coincide con la fórmula de A-1 más su corrección exacta a N chico, a ≤ 3·SE.

**A-3 (clase geométrica, vota; nivel 2 con supuesto declarado).** La clase geométrica no admite conteo exacto. Se usa la aproximación de cumulantes −ln|φ| ≈ (λε²)² Var(S_ε/ε²)/2, que es un **supuesto** (nivel 4 hasta verificarse). Se mide Var(S_ε/ε²) ∝ N^γ para N ∈ {64, …, 1024} y ε ∈ {0.25, 0.1}. Además se verifica la aproximación donde |φ| sí es medible (λε² pequeño, N chico).
- **Condición necesaria para que la rama L viva:** γ + 2·SE < 2. Con la entropía de los 2-órdenes (≈ 2N ln N) y la cancelación de las capas (∝ N²), eso implicaría que la geometría supera a las capas dentro de la ventana.
- **No es suficiente:** no cubre familias no geométricas que no sean por capas.
