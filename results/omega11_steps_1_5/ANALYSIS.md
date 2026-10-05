# Ω-1.1 — Análisis de los pasos 1–5 (P§22)

Analista: cerebro orquestador (Opus). Guardado por la sesión principal.

- **Base de ejecución:** commit `167efca`, árbol limpio.
- **Datos:** `summary.json` de esta carpeta y los `.npz` de `runs/omega11/`, que solo se leyeron.
- **Diagnósticos propios:** los marcados *exploratorios* **no estaban preregistrados**.

## 0. Resumen
- **Integridad.** Ningún paso relajó la compuerta y todo corrió sobre el mismo commit. Los `config_hash` de entrada coinciden con la tabla de enmiendas. Los de cada `summary.json` son de la configuración derivada de cada experimento, como estaba previsto.
- **Resultado global.** De las **1290 corridas dinámicas** (O-01 + O-04), **ninguna pasa los campos de corrida**. Ninguna corrida de dinámica Ω activó jamás un campo dimensional (d_volume, d_spectral, d_weyl). Códigos primarios: F0 = 560, F1 = 530, F2 = 200.
- **Única predicción fallida.** Ω-B con γ̂=10, factor 0.5 y ρ ∈ {0.05, 0.1}. El estado es reproducible bit a bit, es un punto KKT y forma una **clique (o varias cliques)**, en algunas semillas más un halo de peso bajo. No es geometría: es justamente la estructura que advertía P§3.

## 1. Resultados por paso
| Paso | Pregunta | Resultado | Preregistro |
|---|---|---|---|
| p01/p02 | Tests analíticos (w*, α_c, tres umbrales) y golden de Ω-1.0 | rc = 0 | cumplido |
| O-00 | ¿Reproducen los detectores las tablas de calibración? | 31/31 | cumplido |
| O-01 | S0 con la taxonomía v1.1 | 18/18 puntos; 0 corridas pasan | cumplido |
| O-02 | ¿Discrimina la suite de distancias? | 4/4 (una cumplida vacuamente, §3) | cumplido |
| O-03 | Falsos positivos en nulos y controles positivos | 15/15; 0/65 nulos pasan; RGG3 5/5 | cumplido |
| O-04 | Ablación y Ω-B | 49/51 | **2 celdas fallidas** |

**O-00, valores medidos:**
- **D_Weyl:**
  - anillos: 1.01;
  - toros 2D: 2.13, 2.08 y 2.07;
  - toros 3D: 2.88, 2.98, 3.08 y 3.15.
- **RGG3 con N=800 y k=12:** D_vol 3.04–3.08, D_s 2.62–2.64 y D_W 2.85–2.92. Pasa 3/3.
- **Controles que deben fallar, y su primario:**
  - WS p=0.05: F3 (D_vol 4.28, D_s 1.90, D_W 1.64);
  - ER k12: F3, con F4 en el conjunto;
  - árbol: F3;
  - slab: F4;
  - K_100: F1;
  - vacío: F0.
- **Cautela:** los toros 2D y 3D también pasan los campos *de corrida*. La clase 3 solo se exige a nivel de punto: «pasa la corrida» ≠ «3D».

## 2. Baseline S0 (O-01)
Malla: N=200, 10 semillas, α̂ ∈ {0.5, 1, 1.5, 1.9, 2.1, 3} × γ̂ ∈ {0, 1, 10}.

| α̂ | Primario (cualquier γ̂) | ⟨W⟩ final | Ω-1.0 | Estado |
|---|---|---|---|---|
| 0.5–1.9 | F0 en 100% | 1e-11 a 2e-9 | A ×10 | CONVERGED |
| 2.1, 3 | F1 en 100% | 1.0 exacto | E ×10 | CONVERGED |

- **Se confirman los teoremas del baseline:** la frontera de cuenca está en α̂=2 (w* = 1/α̂ frente a ⟨W0⟩ = 1/2), y γ̂ no cambia el resultado.
- **Campos que llegan a pasar:** solo `connected` y `locality`, y únicamente en F1, donde K_N es trivialmente local. Ningún campo dimensional pasa nunca.
- **Cautela de lectura:** en las corridas triviales, F2–F9 aparecen porque la evidencia ausente cuenta como False. Solo el código primario es interpretable ahí.

**O-03, nulos** (N=800, 13 familias, 65 corridas; ninguna pasa). Tasa de verdadero por campo:

| Campo | Tasa | Familias que lo activan |
|---|---|---|
| nontrivial, connected | 1.00 | todas (por construcción) |
| locality | 0.15 | ER k30, WS p0.01 |
| **d_volume** | **0.00** | — |
| d_spectral | 0.06 | WS p0.01 |
| d_weyl | 0.15 | WS p0.01, WS p0.05 |
| metric_robust | 0.08 | WS p0.05 |
| isotropy_ok | 0.14 | WS p0.01, WS p0.05 |
| homogeneity_ok | 0.15 | WS p0.05, árbol |
| topology_stable | 0.23 | ER k30, WS p0.01, árbol |
| **manifold_proxy_ok** | **0.00** | — |

- **Discriminadores fuertes:** d_volume y manifold_proxy_ok, con 0 falsos positivos.
- **Red-team:** WS p=0.2 da D_s = 3.22 y D_W = 2.92. Lo atrapan la localidad (F3) y la falta de ventana en D_vol, lo que confirma la regla «D≈3 nunca basta».
- **Controles positivos:**
  - RGG3 k12: pasa 5/5 (D* = 2.91, clase 3).
  - **RGG2 k10: falla 5/5 con F9** (D* = 1.86). La conectividad de anillos da 0.64, 0.36 y 0.20 en r = 2, 3 y 4.
  - El proxy de manifold está calibrado para k≈12 en 3D y **da un falso negativo en una geometría 2D genuina con k=10**. No afecta a P1, que es 3D, pero documenta el riesgo R3: estados geométricos con grado bajo podrían rechazarse.

## 3. O-02: sensibilidad métrica
| Familia | Fracción sensible | Motivo real |
|---|---|---|
| Final binario de Ω-1.0 (α̂=1.5, γ̂=1) | 1.0 | gigante trivial (F0) |
| Final del integrador sigmoide | 1.0 | MAX_STEPS, ⟨W⟩=0.017, ninguna arista > w_min |
| U(0,1) estático, w_min 0.01–0.5 | 1.0 | sin ventana (denso) |
| Langevin, Θ̂ = 0.1 y 0.3 | 1.0 | sin ventana (⟨W⟩ 0.20/0.38, tipo ER denso) |
| Toro 9³ con pesos aleatorios | **0.0** | dispersión 0.30 |
| Toro 28² con pesos aleatorios | **0.0** | dispersión 0.04–0.07 |
| RGG euclidiano k12, N=800 | **1.0** | dispersión 0.65–0.69 (HOP 2.86, INV 2.15–2.22) |

- **La suite sí discrimina en los controles.**
- **En las familias Ω, «sensible» significa «sin ventana».** Ningún estado de la dinámica Ω tuvo ventana de escala, así que **la suite nunca se ejerció sobre un estado Ω no trivial**.
- La expectativa `max_finite_spread == 0` se cumplió **vacuamente**, porque no hubo ninguna dispersión finita.

## 4. O-04: ablación (N=100; FULL y NO_TRIANGLES también con N=200)
| Funcional | Resultado | Lectura |
|---|---|---|
| FULL | F0 con α̂ ≤ 1.5; F1 con α̂ ≥ 2.5 (todo γ̂, N = 100 y 200) | umbral α̂=2 |
| NO_TRIANGLES (T=0) | **F0 en 100%** (N = 100 y 200) | cumplido |
| TRIANGLES_ONLY | **F1 en 100%** | cumplido |
| NO_DENSITY | F1 en 100%, incluso con α̂=0.5 | sin β nada frena a T |
| NO_DEGREE | idéntico a FULL | γ es inerte |

**Respuesta a P§3/M§33: no hay ninguna propiedad emergente que atribuir.**
- T es el único motor de estructura, y lo que produce es densificación: F1 en S0 y cliques F2 en Ω-B.
- β produce vaciamiento (F0).
- γ no cambia ningún código.
- Con T=0 no hay estructura.

La preocupación red-team de P§3 queda **confirmada empíricamente**: −αT empuja hacia clustering → hiperdensidad o cliques, nunca hacia una estructura tipo retícula.

## 5. Ω-B: la predicción fallida y los artefactos
### 5.1 Tabla (N=100, 10 semillas; α̂ = f·α̂_c, con α̂_c = (1+γ̂)(N−2)/(ρ(N−4)))
| ρ | γ̂ | f=0.5 (preregistro: F1) | f=1.5 | f=3 |
|---|---|---|---|---|
| 0.05 | 0 / 1 | primario F0, con F1 en el conjunto (uniforme, ver §5.4) | F2, clique de 23 | F2 |
| 0.05 | 10 | **F2 10/10 (fallida)** | F2 | F2 |
| 0.1 | 0 / 1 | F1 (uniforme) | F2, clique de 32 | F2 |
| 0.1 | 10 | **F2 10/10 (fallida)** | F2 | F2 |
| 0.2 | 0 / 1 / 10 | F1 (uniforme) | F2, clique de 45 | F2 |

### 5.2 Caracterización del estado fallido
- **Reproducibilidad:** se reconstruyeron 3 corridas (002145, 002146 y 002244); w0 y w_final coinciden bit a bit.
- **ρ=0.05 (α̂=112):**
  - 7/10: clique de 22 a peso 1, más un nodo unido a ella con peso 0.75 y 77 nodos aislados. Es el análogo ponderado del grafo colex, el mismo estado que aparece con γ̂=0 por encima del umbral.
  - 1/10: la misma clique más un halo uniforme de 0.0055 sobre 3003 aristas.
  - 2/10: dos cliques disjuntas.
- **ρ=0.1 (α̂=56):**
  - 4/10: clique de 32.
  - 4/10: clique de 30–31 más un halo de 0.013–0.025.
  - 2/10: varias cliques.
- **Rangos corregidos** (el README congelado decía «5–9%» y «23–30%»): las aristas por encima de w_min son el 4.7–5.2% con ρ=0.05 y el 8.8–24.4% con ρ=0.1. La componente gigante abarca el 18–23% y el 28–32%, respectivamente.
- **Son puntos KKT en la frontera:**
  - aristas libres: |g+ν| ≤ 1e-7;
  - aristas en 0: g+ν ≥ 0;
  - aristas en 1: g+ν ≤ 0.
- **La fragmentación no es un artefacto del umbral:** el soporte (W > 0) también está fragmentado.

### 5.3 Por qué falló la predicción lineal (diagnósticos *exploratorios*)
1. **El uniforme sigue siendo localmente estable.** Desde W = ρ(1+εξ) vuelve a uniforme en 12/12 casos. El criterio lineal es correcto; lo que falló fue suponer que U(0,1) proyectado cae en su cuenca.
2. **El uniforme es metaestable, no global.** H(uniforme) frente a H(clique):
   - ρ=0.05, γ̂=10: −34 frente a −2744;
   - ρ=0.1, γ̂=10: −136 frente a −3030;
   - incluso ρ=0.05, γ̂=0, f=0.5: +8.2 frente a −104.5.
3. **La cuenca depende de la amplitud inicial c y del α̂ absoluto:**
   - ρ=0.05, γ̂=10: c=0.2 lleva al uniforme; c ≥ 0.4, a la clique.
   - ρ=0.1, γ̂=10: c ≤ 0.6 lleva al uniforme; c=1, a la clique.
   - Con γ̂=10, la frontera en f está entre 0.25 y 0.5.
4. **Hipótesis de mecanismo:**
   - γ solo endurece los modos de grado.
   - Los ~N²/2 modos de ciclo tienen una rigidez ≈ 2.2β, independiente de γ̂.
   - Mantener f fijo exige α̂ ∝ (1+γ̂), así que el término cúbico crece ~11 veces sin que crezca esa rigidez.
5. **Kruskal–Katona (forma de Lovász).** Con m aristas binarias fijas, T ≤ (√2/3)·m^{3/2}, con igualdad asintótica en la clique. Como S_dens = m es constante, **el mínimo global binario de Ω-B con γ=0 es la clique o el grafo colex**. Con γ>0 solo cabe repartir la masa en cliques más regulares. **Ω-B a densidad fija favorece estructuralmente las cliques.** Es riguroso para estados binarios y heurístico para los ponderados.

### 5.4 Artefactos de medición (no son hallazgos)
- **ρ=0.05, γ̂ ∈ {0, 1}, f=0.5: primario F0.**
  - El estado real es el uniforme exacto: todas las aristas valen 0.05, con cv ~1e-8.
  - Con w_min = 0.1 > ρ, el grafo binario queda vacío. Se dispara F0, que por precedencia gana a F1.
  - La expectativa se dio por cumplida porque juzga «código ∈ conjunto».
  - El impacto en las conclusiones es nulo, pero en una tabla de primarios se lee «vacío» donde hay un K_N uniforme.
- **ρ = 0.1 = w_min (celdas F1).**
  - Hay redondeos de ±1e-9 alrededor de 0.1, y el grafo binario resulta aleatorio, fabricado por el error numérico.
  - Las medidas binarias de esas corridas se calculan sobre ruido. El primario F1 es correcto por la regla de cv.
- **Halos de 0.0055–0.025:** son invisibles a las medidas binarias. Aquí no cambian el código.
- **Tratamiento:** solo documentación. Las celdas Ω-B con ρ ≤ w_min se marcan «primario no interpretable; uniforme por cv». Una corrección futura sería una enmienda preregistrada de Ω-1.2 y es decisión del usuario.

## 6. Estado de P1
- **Códigos dominantes:**
  - F0: S0 con α̂ < 2, T=0 y NO_DEGREE con α̂ < 2;
  - F1: S0 con α̂ > 2, TRIANGLES_ONLY, NO_DENSITY y el uniforme de Ω-B;
  - F2: cliques de Ω-B;
  - en los nulos: F3 y F4.
- **P1: sin evidencia a favor** en todo el dominio explorado. El fracaso está explicado analíticamente: umbrales homogéneos, γ inerte, T=0 ⇒ F0 y Kruskal–Katona en Ω-B.
- **Alcance de la refutación:** P1 queda refutada para la dinámica de gradiente S0 y para Ω-B determinista en este rango. No queda refutada globalmente: faltan Θ>0, N grande y otras inicializaciones.
- **Ninguna región candidata para s06.** Su preregistro ya anticipa F0 en 100% con α̂ ≤ 1.9, F1 en 100% con α̂ ≥ 2.1, 0 candidatos y χ_max ≤ 1e-12. Su valor informativo es bajo: confirma O-01 en una malla más densa y cierra el orden de P§22.

## 7. Siguientes pasos
| # | Paso | Coste | Comentario |
|---|---|---|---|
| 1 | s06 sin cambios | ~35 min | Bloqueado por D2 (confirmar R1, R3, R4 y R8). |
| 2 | O-04b, experimento nuevo para Ω-B | ~3–6 h | (a) test analítico de Kruskal–Katona; (b) mapa de cuenca ρ × γ̂ × f × c con N = 100 y 200. Predicción: F1 o F2, nunca ventana dimensional. Si se confirma, Ω-B queda cerrada como vía hacia P1 bajo S0. Decisión del usuario. |
| 3 | O-05, ensembles con Θ>0 | ~5–15 h | La única vía abierta para S0 (R8). Conviene un piloto de tiempos. Predicción: F0/F1 con Θ̂ bajo; ER denso sin ventana (F3/F4) con Θ̂ intermedio. La región dominante debe confirmarse con Metropolis (A-12). |
| 4 | O-06, tamaño finito | ~1–3 h | Sobre F0/F1 no aporta geometría. Su valor es calibrar `size_robust` con RGG3 en N = 800, 1500 y 3000 (R1). |
| 5 | Documentación | — | Artefacto ρ ≤ w_min, expectativa vacía de O-02, F5/F9 por evidencia ausente, falso negativo de RGG2 y rangos del README. |

**Decisiones del usuario:**
- D2;
- autorizar s06;
- O-04b y el teorema de Ω-B;
- orden y presupuesto de O-05 frente a O-06;
- tratamiento de ρ ≤ w_min en Ω-1.2.

Según el criterio del usuario, nada de esto justifica modificar la funcional S0: F0, F1 y F2 son el resultado científico del baseline.
