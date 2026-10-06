# Ω — Fase 2: de la localidad a la geometría — Prerregistro

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-fase2-theta`.
- **Base congelada:** `claude/omega-c0-congelado` (`a25b280`), con Ω-C0 completo; `claude/omega-1.1-congelado` (`ec448e4`).
- **Estado:** este documento se commitea **antes** de escribir el código de la fase. Nada de lo que sigue usa resultados de Θ > 0.

## 0. Decisiones del Consejo (sesión posterior a C0-F1)

| # | Decisión | Estado en esta fase |
|---|---|---|
| D1 | C0 es un **resultado positivo parcial, no una teoría exitosa**. Logra la localidad (A); no logra variedad (B) ni dimensión (C) | Registrado |
| D2 | Congelar C0 antes de modificar nada | Hecho: `claude/omega-c0-congelado` |
| D3 | **No implementar C1** sin análisis previo del paisaje (C1-A, B, C) | Solo se registra el protocolo (§3); no hay dinámica C1 |
| D4 | **Θ > 0 sobre C0** como experimento barato, sin esperar geometría | Se ejecuta: protocolo C0-Θ (§2) |
| D5 | No declarar «dimensión 2.5»: la fase es local y no variedad, con dimensión efectiva baja y no certificada | Lenguaje obligatorio en resultados |
| D6 | No añadir agujeros negros ni otros objetos físicos hasta tener geometría y causalidad | Registrado |
| D7 | O-05 sigue como rama secundaria; no hay N grande ni barridos paramétricos | Registrado |

## 1. Análisis del cerebro (antes de correr)

### 1.1 Reformulación: tres propiedades separables

- **A. Localidad:** las aristas se concentran en vecindarios. Se mide con H_null y ciclos cortos.
- **B. Variedad:** las bolas crecen de forma homogénea e isótropa, el certificado no da F9 y la dimensión es estable en escala.
- **C. Dimensión:** la de B vale ≈ 3 y no está impuesta.

C0 demostró experimentalmente que A no implica B (C0-F1).

### 1.2 Por qué una acción local de baja orden no puede dar B: generalización de C0-T1/T2 (nivel 2)

C0-T1 muestra que, cuando la acción solo depende de estadísticas locales por arista (k, c), los estados de igualdad de la cota son los grafos **homogéneos por arista**: todas las aristas tienen los mismos valores (k*, c*).

Esperamos que esto valga para cualquier acción construida igual (penalización convexa de la fuerza más una recompensa por arista que solo depende de conteos locales):
- sus mínimos globales son estructuras **cristalinas**, homogéneas en esos conteos;
- una geometría desordenada tipo RGG3 nunca es mínimo global, como mucho metaestable;
- la degeneración no desaparece mientras existan grafos no geométricos con los mismos conteos.

En C0 los hay en abundancia: K_{m,m}, cliques difusas y complejos de cliques.

Es nivel 2 porque cada nueva acción necesita su propia demostración (C1-A).

### 1.3 El riesgo central de C1: convertir la dimensión en un parámetro (cálculo exacto, nivel 1)

Sea q = número de 4-ciclos por arista. En la familia de retículos hipercúbicos Z^d (y sus toros, con lado ≥ 5) cada arista tiene

    (k, c, q) = (2d, 0, 2(d − 1)).

El hipercubo Q_k tiene (k, 0, k − 1).

Una recompensa saturante en q con óptimo q*, más el término de grado de C0 con óptimo k*, tiene estos mínimos homogéneos:
- **Z^d** si (k*, q*) = (2d, 2(d − 1));
- **Q_k** si q* = k* − 1;
- otros grafos con esos conteos en el resto de casos.

Por tanto, C1 en su forma ingenua **elige la dimensión a través de los parámetros**: d = q*/2 + 1 en la familia reticular. Eso es imponer D, no hacerlo emerger, y el Consejo lo prohibió (k = 6, D = 3 como objetivo).

Existe un resultado externo de **rigidez local** (nivel 3, *referencia a verificar antes de citarla formalmente*): un grafo cuyas bolas de radio r ≥ r₀ coinciden con las de Z^d es un cociente de Z^d. Significa que la información de segundo orden **sí puede fijar** una variedad, pero solo si los parámetros ya codifican la de Z^d.

**Criterio nuevo, C1-D (meseta), propuesto por el cerebro para que el Consejo lo ratifique:** una dimensión cuenta como emergente solo si se selecciona sobre una **región abierta** del espacio de parámetros, en la que D es constante, y no en un punto ajustado a los conteos de un retículo concreto.

### 1.4 Predicción para Θ > 0 (nivel 2, cálculo de campo medio)

Para un par lejano (i, j) que no es arista (C_ij ≈ 0, sin vecinos comunes) en un estado con fuerzas ≈ k*, el gradiente es

    G_bg = ψ* − a = b·ln(1 + c*) − λ c*  ≥ 0     (= 0 si c* = 0).

Con Langevin a temperatura Θ, cada uno de esos ~N pares por nodo adquiere un peso medio ≈ Θ / G_bg. El fondo aporta así una fuerza ≈ N Θ / G_bg, que iguala la fuerza objetivo k* cuando

    Θ_N := G_bg · k* / N.

**Predicciones registradas:**
- **P-Θ1.** Con Θ ≤ 0.1 Θ_N la fase local sobrevive.
- **P-Θ2.** Con Θ ≥ 3 Θ_N la fase se destruye: el fondo domina la fuerza y el soporte fuerte deja de ser DISPERSO-LOCAL.
- **P-Θ3.** La frontera escala como 1/N. Expresadas en Θ/Θ_N(N), las curvas de N = 216 y N = 343 coinciden con una tolerancia de un punto de la rejilla.

Consecuencia si P-Θ1–3 se cumplen: a Θ fija, la fase local de C0 desaparece cuando N → ∞. Sería una fase de Θ → 0 o de tamaño finito: la entropía de los ~N²/2 pares de fondo gana. Esto cuantifica la previsión de `OMEGA_C0_RESULTADOS.md` §3.

| Celda | c* | k* | a | G_bg | Θ_N (N = 216) | Θ_N (N = 343) |
|---|---|---|---|---|---|---|
| 19 | 1 | 6 | 1 | 0.1931 | 0.00536 | 0.00338 |
| 37 | 2 | 8 | 1 | 0.4319 | 0.01600 | 0.01007 |
| 55 | 4 | 12 | 1 | 0.8094 | 0.04497 | 0.02832 |
| 73 | 8 | 16 | 1 | 1.3083 | 0.09691 | 0.06103 |

## 2. Protocolo C0-Θ (se ejecuta)

### 2.1 Integrador

Langevin con reflexión sobre el vector triangular u (una variable por par), igual que el motor de Ω-1.1 pero con la funcional C0:

    u ← B(u − dt·G(u) + sqrt(2 Θ dt)·ξ),   ξ ~ N(0, 1),   B = reflexión en [0, 1].

- G es el gradiente por par de S_C0 (`grad_c0`, la misma derivada que usa `evolve_c0`).
- La densidad estacionaria es ∝ exp(−S_C0 / Θ).
- Θ es una temperatura estadística, no tiempo físico.
- **dt = 0.02**, fijo.
  - **Regla de respaldo:** si en el brazo Θ = 0 alguna corrida deja de ser DISPERSO-LOCAL, o S/LB empeora (sube) más de un 1 % relativo entre el inicio y el final, dt se reduce a la mitad y se repite toda la tabla. Se permiten como máximo dos reducciones, y se informan.
- Un solo hilo BLAS.
- RNG: PCG64 con clave (20261006, bloque, N, celda, inicio, semilla, índice de Θ).

### 2.2 Diseño

- **Celdas:** 19, 37, 55, 73. Es la diagonal R-3D de C0, todas CANDIDATO-C0, con c* = 1, 2, 4 y 8.
- **Rejilla:** Θ / Θ_N(N) ∈ {0, 0.1, 0.3, 1, 3}.

**Θ-S (supervivencia):** parte de los estados finales oficiales de C0 a Θ = 0.
- N = 216: finales de C0-L4 (`runs/c0/…`) con inicios U y R, semillas 0–2.
- N = 343: finales de R6 (`runs/c0_redteam/out/R6_c{celda}_R_n343_s{s}.npz`) con inicio R, semillas 0–2.
- 20 000 pasos.
- Total: 4 × (2 × 3 × 5 + 1 × 3 × 5) = **180 corridas**.

**Θ-F (formación):** parte de un inicio R genérico nuevo (`build_input("R")`).
- N = 216, Θ/Θ_N ∈ {0.3, 1, 3}, semillas 0–2, 20 000 pasos.
- Total: **36 corridas**.
- La comparación a Θ = 0 es el resultado oficial de C0-L4 (no se vuelve a correr).

### 2.3 Observables (por corrida, al final y cada 1000 pasos)

- **clase:** `classify_c0(W)`, sin cambios.
- **f_bg:** fracción de la fuerza total en pares con W < 0.1·max W (el «gas» de fondo).
- **J:** índice de Jaccard entre el soporte fuerte final y el inicial. Solo en Θ-S.
- **m_rel:** número de aristas del soporte fuerte / (N k*/2).
- **S/LB:** al final y en serie temporal.
- **Estacionariedad (informe):** diferencia relativa entre la media de S en el último cuarto y en el penúltimo cuarto.

### 2.4 Criterios

Para cada combinación (celda, N, inicio, Θ):
- **SOBREVIVE:** ≥ 2/3 de las semillas son DISPERSO-LOCAL con f_bg ≤ 0.5, y además, en Θ-S, J ≥ 0.5.
- **DESTRUIDA:** ≥ 2/3 de las semillas no son DISPERSO-LOCAL o tienen f_bg > 0.5.
- **MIXTA:** cualquier otro caso.

**Frontera.** Θ_c / Θ_N es el mayor punto de la rejilla con SOBREVIVE, siempre que todos los puntos menores también sobrevivan; si alguno menor no sobrevive, Θ_c / Θ_N = 0. Se calcula por celda, N e inicio.

**Veredicto global de Θ-S**, sobre las 4 celdas con inicio R a N = 216 y 343:

| Veredicto | Condición | Lectura |
|---|---|---|
| **FRÁGIL-1/N** | En las 4 celdas Θ_c/Θ_N está en {0.1, 0.3, 1} y difiere entre N = 216 y 343 en ≤ 1 punto de rejilla | P-Θ1–3 confirmadas: fase de Θ → 0 |
| **ROBUSTA** | SOBREVIVE en Θ = 3 Θ_N en ≥ 3 de las 4 celdas a los dos N | La entropía no la destruye en esta escala; la heurística §1.4 queda refutada |
| **INESTABLE** | En ≥ 3 de las 4 celdas Θ_c/Θ_N = 0 a los dos N | La fase solo existe a Θ = 0 |
| **NO CONCLUYENTE** | Cualquier otro caso | — |

**Θ-F, por celda:**
- **FORMA** si en Θ = 0.3 Θ_N la combinación SOBREVIVE (≥ 2/3 semillas DISPERSO-LOCAL con f_bg ≤ 0.5);
- **NO FORMA** en caso contrario.

Lo que **no** se hace: batería geométrica ni certificado. La pregunta es solo de estabilidad (D4); la batería se reserva para el caso ROBUSTA.

## 3. Protocolo C1 de análisis previo (registrado, **no ejecutado**)

Se ejecutará solo cuando el Consejo elija la forma funcional de C1. Las pruebas son de paisaje, sin dinámica.

| Prueba | Pregunta | Método | Muerte |
|---|---|---|---|
| **C1-A** | ¿Qué estructuras favorece? | Cota tipo C0-T1 y sus estados de igualdad; S/LB sobre T³, Q_6, K_{m,m}, cliques, ER, RGG3, toro triangular y los finales C0-R de N = 729 (`runs/c0_f1/out/`) | Si el mejor no geométrico empata o gana a T³ en toda la rejilla |
| **C1-B** | ¿Hay mínimos geométricos estables? | KKT y **hessiana** en T³ (no basta un punto crítico, lección de C0-L3) | Si T³ no es KKT o tiene direcciones descendentes en toda la rejilla |
| **C1-C** | ¿Distingue la geometría de la mera localidad? | S_C1(T³ y RGG3) frente a S_C1(finales C0 de N = 729, que son locales y no variedad), por arista y a igual densidad | Si S(G_local) ≤ S(G_geom) |
| **C1-D** | ¿La dimensión es emergente o impuesta? (§1.3) | Si C1-A/B pasan, se mapea la región de parámetros en la que el mínimo homogéneo es cada Z^d o Q_k | Si D = 3 solo aparece en el punto ajustado (k*, q*) = (6, 4) y no en una región abierta |

Solo si C1-A…D pasan, C1 entra en la puerta dinámica de C0 (L4, red-team, F1 con N = 729).

**Información de apoyo (descriptiva, no es C1):** la herramienta `tools/c1_local_stats.py` mide la distribución por arista de (k, c, q) en las referencias y en los finales C0, para saber qué conteos tendría que distinguir cualquier C1.

## 4. Exclusiones

- Ninguna funcional nueva se simula en esta fase.
- No se cambian umbrales del certificado ni código congelado.
- Ni mecánica cuántica, ni coordenadas, ni embedding, ni D = 3 o k = 6 como objetivo.

## 5. Enmiendas

### F2-A1: paso dt (registrada tras el smoke y antes de la corrida oficial; no se ha visto ningún dato de Θ > 0 oficial)

**Qué mostró el smoke.** El brazo de control Θ = 0 con dt = 0.02, aplicado al final C0 de la celda 37 (inicio R, semilla 0), empeora S/LB de 0.9997 a 0.945, un 5.5 %, y la cifra queda estable en el tiempo.

**Diagnóstico.** El déficit es lineal en dt:

| dt | Déficit de S/LB |
|---|---|
| 0.02 | 5.5 % |
| 0.01 | 2.7 % |
| 0.005 | 1.35 % |
| 0.0025 | 0.68 % |

Es el sesgo O(dt) de la reflexión en la pared W = 1. Las aristas saturadas rebotan hasta 1 − dt·|G|. No es una inestabilidad.

**Problema con la regla original.** La regla de §2.1 solo permitía dos reducciones, hasta dt = 0.005, y con ese paso el control sigue sin cumplir el criterio del 1 %.

**Enmienda:**
- dt = 0.0025 (tres reducciones);
- n_steps = 40 000, es decir, tiempo de integración 100 frente a 400 en el diseño original;
- el resto del protocolo no cambia.

**Limitaciones que se informarán:**
- «sobrevive» significa que sobrevive durante un tiempo 100. La observable de estacionariedad indica si sigue derivando;
- el sesgo de pared (≈ 0.7 % en S) afecta igual a todos los brazos.

Observación del smoke (cuatro corridas de 400 pasos con dt = 0.02, **no oficiales**): con Θ = Θ_N, la celda 37 pasó a DENSO-TRIVIAL en un tiempo de 8. Se registra solo para transparencia y no altera los criterios.
