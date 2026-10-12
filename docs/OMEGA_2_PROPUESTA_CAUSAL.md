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
