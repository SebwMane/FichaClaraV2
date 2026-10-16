# Ω-C0 — Resultados de la puerta de competencia relacional

- **Fecha:** 2026-10-05.
- **Rama:** `claude/omega-c0-competencia`.
- **Especificación:** `docs/OMEGA_C0_PRERREGISTRO.md`, con §0–§7 y las enmiendas C0-A1…A7. Cada enmienda está fechada respecto de los datos que se habían visto al registrarla.
- **Línea anterior:** Ω-1.1 y el bloque L están congelados en `claude/omega-1.1-congelado` (`ec448e4`) y no se modificaron.
- **Código:**
  - módulos y tests: `omega/c0/` y `tests/test_c0_*.py` (39 tests junto con los de arquitectura);
  - herramientas: `tools/c0_landscape.py`, `c0_dynamics.py`, `c0_redteam.py`, `c0_battery.py`, `c0_null_battery.py` y `c0_f1.py`.

## 1. Resultados oficiales

| Prueba | Decisión | Datos |
|---|---|---|
| **C0-L1** ¿siguen siendo extremales las cliques? | **SOBREVIVE** | La clique binaria es la mejor referencia en 5/75 celdas, todas con c* = 8 y k* ≤ 8, la región que predice T1 |
| **C0-L2** ¿se rompe la dominancia trivial en el paisaje? | **EMPATE** | 0 celdas rotas, 12 empatadas. T³ (c* = 0) y el toro triangular (c* = 1, 2) alcanzan la cota, pero también la alcanzan K_{108,108} y las cliques difusas |
| **C0-L3** KKT | **KKT PRESENTE** (no estricto) | T³ es KKT en 6 celdas (c* = 0, k* ∈ {4, 6}), con direcciones planas. RGG3 no es KKT en ninguna celda |
| **C0-L4** ¿aparece un tercer atractor? | **CONTINÚA** | De 675 corridas genéricas: VACÍO 0, DENSO-TRIVIAL 178, FRAGMENTADO 30, **DISPERSO-LOCAL 253**, DISPERSO-NO-LOCAL 214. 43 celdas votan LOCAL |
| R3 hubs | Pasa | k_max/k̄ ≤ 2.26 en todos los finales locales |
| R8 inicialización | 30/43 celdas | LOCAL desde ≥ 2 de los 3 inicios genéricos |
| R7 reetiquetado (C0-A1) | Pasa | Código: 10/10 con ≤ 1e-8 a 60 pasos. Resultado: misma clase en 10/10; |ΔS|/S ≤ 2e-9 en las convergidas |
| **R6 tamaño** (N = 125 y 343) | **24/30 celdas** | Con N = 343, 208/213 corridas son LOCAL. Las 6 celdas que fallan, fallan en N = 125, no en 343: la localidad **crece** con N |
| D-1 convergencia (100 000 pasos) | Persisten | 156/157 finales LOCAL no convergidos siguen LOCAL; 64 convergen. No son transitorios hacia cliques |

**Resultado de la puerta:** las 24 celdas son **CANDIDATO-C0**.

### Cómo se distribuye la fase

| c* | Corridas genéricas DISPERSO-LOCAL |
|---|---|
| 0 | 1/135 (el régimen c* = 0 da bipartito denso o disperso no local) |
| 1 | 60/135 |
| 2 | 73/135 |
| 4 | 71/135 |
| 8 | 48/135 |

## 2. Qué es la fase: batería descriptiva (C0-A6) y control nulo (C0-A7)

Todas las cifras son descriptivas. No se usa ningún umbral nuevo de geometría.

| Familia | Dónde | Estructura | Escala |
|---|---|---|---|
| **«Cadenas»** (inicio U) | casi todas las celdas candidatas | cadenas de cliques solapadas (clique máxima 10–30) | D_L ≈ 0.8–1.4; D_s ≈ 0.9: esencialmente **1D** |
| **«Expansor decorado»** (inicio E) | c* ≥ 1 | cliques pequeñas (6–10) con ~2 cliques por nodo | D_L ≈ 3.7–6, cerca del nulo (≈ 5). En 4/13 no se separa del nulo: crecimiento ~log N |
| **«R-3D»** (inicio R) | celdas 18–20, 36–38, 55–56, 73–74; diagonal (c*, k*) = (1, 6), (2, 8), (4, 12), (8, 16) | k̄ ≈ 12–18, clique máxima 6–11, muchas cliques solapadas pequeñas | **D_L ≈ 2.5–3.7, D_s(N = 343) ≈ 2.6–3.3**, separado del nulo (D_L nulo ≈ 4.6–5.1) |

**Control nulo:** 36/40 combinaciones (celda, inicio) se separan de su recableado con grados fijos (D_L < D_L nulo − 0.5). La fase no es un expansor disfrazado, salvo parte de la familia E.

**Certificado Ω-1.1** con N = 343 (sin cambios de umbral; informe, no voto):

| Grafo | Códigos |
|---|---|
| T³ 7³ | ninguno (pasa) |
| RGG3 k12 | F4, F5, F9 |
| 2 finales R-3D de la celda 37 | F10, F4, F5, F9 |

Con N = 343, el certificado **no distingue** la subfamilia R-3D de un RGG3 auténtico. El código añadido, F10, se debe a la falta de convergencia. Hace falta más N: esa es la prueba **C0-F1** (prerregistro §7).

## 3. Escala de afirmación

**Demostrado (nivel 1)**
- C0-T1 y C0-T2: la cota S ≥ −Nψ*²/(16κ) y sus estados de igualdad, con verificación numérica.
- Consecuencia: una acción que solo depende de (k, c) **no puede preferir energéticamente la geometría**. T³ y el toro triangular son mínimos globales **degenerados** con estados densos no locales.

**Verificado en el dominio probado (N ≤ 343, Θ = 0)**
- C0 rompe la dicotomía vacío/clique del bloque L.
- Desde inicios genéricos aparece una fase dispersa, local (H_null ≥ 1.15, ciclos cortos), sin hubs y persistente en tiempo (D-1) y en tamaño (R6).
- Se separa del nulo de grados.

**Explicación propuesta (nivel 2)**
- La fase es **dinámica, no energética**. Por C0-T2, los estados locales no tienen menor acción que los triviales: son mínimos locales y valles planos que el descenso alcanza desde condiciones genéricas.
- El tipo de estructura (1D, expansor o R-3D) depende de la **condición inicial**. Eso indica un paisaje vidrioso con muchas cuencas, no una fase única seleccionada por la funcional.

**Abierto, NO demostrado (nivel 4)**
- Que la subfamilia R-3D sea geométrica. D_L ≈ 3 y D_s ≈ 3 no bastan (`RULE_D3_NEVER_SUFFICIENT`). Con N = 343 tampoco el RGG3 real pasa el certificado.
- Θ > 0. La degeneración de C0-T2 hace prever que la entropía favorezca los estados densos o expansores, pero no está probado.

## 4. Contraste con las predicciones del cerebro (registradas antes de correr)

| Prueba | Predicción | Resultado | ¿Acertó? |
|---|---|---|---|
| L1 | sobrevive | sobrevive | ✔ |
| L2 | NO ROTA / EMPATE | EMPATE | ✔ |
| L3 | T³ KKT estricto con c* = 0 | KKT pero no estricto | ✘ en lo de «estricto» (C0-A2) |
| L4 | MUERTE-A/B | **CONTINÚA** | ✘ |

El fallo de L4 es el resultado más informativo de C0: la dinámica encuentra estructuras locales que el paisaje no premia.

## 5. C0-F1: escalado y certificado (N = 512 y 729) — **F1-NEGATIVO**

| Fase | Resultado |
|---|---|
| F1 oficial | Global INDETERMINADO: 4 celdas INDETERMINADO, que solo llegaron a N = 512, donde la RGG3 de referencia tampoco pasa; 6 NEGATIVO |
| F1b (C0-A9) | Las 4 celdas restantes con N = 729 son NEGATIVO |
| **Global** | **F1-NEGATIVO, 10/10 celdas** |

**Referencias del certificado con el mismo N:**

| Grafo | N = 512 | N = 729 |
|---|---|---|
| T³ (8³ y 9³) | pasa | pasa |
| RGG3 k12 | F4, F5, F9 | **pasa** (D_eff ≈ 3.03, D_s ≈ 2.6) |

**Subfamilia R-3D con N = 729 (8 celdas × 3 semillas):**
- Las 24 corridas siguen DISPERSO-LOCAL.
- El certificado da 0/24: siempre F4 (artefacto de dimensión), F9 (no-variedad) y casi siempre F5 (dependencia de la métrica). La RGG3 auténtica, con el mismo N, no tiene esos códigos.
- Dimensiones:
  - D_s = 1.77–2.76 (mediana 2.20);
  - D_eff = 2.22–2.62 (mediana 2.47; 15/24 con ventana);
  - D_L = 2.7–4.1. Las celdas con c* alto derivan hacia el valor del nulo (≈ 4.8).

**Lectura:** la «R-3D» de N = 343 era un efecto de tamaño. Al crecer N, la fase sigue siendo local y no trivial, pero tiene dimensión baja e inhomogénea (≈ 2–2.5) y no es una variedad. Es un complejo de cliques pequeñas solapadas, no una geometría 3D.

## 6. Conclusión de Ω-C0

**Lo que C0 consiguió (verificado):**
- rompe la dicotomía vacío/clique del bloque L;
- produce una fase intermedia local, persistente en tiempo y tamaño, sin hubs y separada del nulo de grados.

Es el primer «tercer atractor» del proyecto.

**Lo que C0 no consiguió:**
- **Paisaje:** por C0-T1 y C0-T2, la geometría nunca gana energéticamente; solo empata con estados densos.
- **Dinámica:** la fase que produce es un complejo de cliques solapadas cuya estructura depende de la condición inicial:

  | Inicio | Estructura |
  |---|---|
  | U | cadenas 1D |
  | E | expansor decorado |
  | R | dimensión ≈ 2–2.5, no variedad |

  En el primer tamaño que discrimina, ninguna supera el certificado.

**Decisión según el árbol preregistrado:**
- C0 superó la puerta L1–L4 y el filtro red-team (24 CANDIDATO-C0), y falló F1.
- **C0 queda limitada:** produce localidad, no geometría.
- La rama no se reabre con búsqueda paramétrica; esa decisión del Consejo sobre S0/Ω-B vale igual aquí.

**Hipótesis de trabajo para el Consejo (nivel 2):** la localidad puede surgir de la competencia y la saturación; la dimensión y la homogeneidad no. Dos datos lo apoyan:
- la degeneración de C0-T2: una acción que depende solo de (k, c) no distingue estructuras con la misma estadística local;
- la dependencia de la condición inicial: el paisaje es vidrioso y no hay selección.

Para seleccionar una dimensión haría falta un ingrediente que distinga la estructura a escala intermedia: relaciones de segundo orden, como sugirió el especialista en emergencia, o un mecanismo entrópico o dinámico.

## 7. Pendiente para el Consejo (no ejecutado)

1. **C1** — estructura de segundo orden: término sobre 4-ciclos o codegrado de no-aristas, con saturación. Antes de simular tiene que pasar la misma puerta (L1–L4, red-team y F1 a N = 729). Advertencia analítica previa: una recompensa cóncava de 2-caminos favorece los expansores e hipercubos frente a T³; la forma funcional debe someterse primero a un análisis tipo C0-T1.
2. **Θ > 0 sobre C0:** ¿la entropía estabiliza o destruye la fase local? Es una pregunta legítima y barata a N ≤ 343.
3. **O-05** (S0/Ω-B, Θ > 0): sigue como rama secundaria.
