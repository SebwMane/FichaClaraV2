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

## 5. Siguiente paso en curso

**C0-F1** (§7 del prerregistro):
- N = 512 en las 10 celdas R-3D y N = 729 en 4 de ellas;
- 40 000 pasos;
- certificado Ω-1.1 completo frente a RGG3 y T³ con el mismo N;
- D_L con 3–4 tamaños.

Decisión preregistrada: F1-POSITIVO, INDETERMINADO o NEGATIVO. Un POSITIVO no declara geometría; abre el paquete de confirmación, que requiere ratificación del Consejo.

## 6. Pendiente para el Consejo

- Ratificar la lectura «fase dinámica, no energética» y decidir si es un problema (falta de selección) o un rasgo (las cuencas dinámicas *son* el mecanismo).
- Si F1 sale POSITIVO o INDETERMINADO:
  - N-1 con nulos con clustering;
  - Θ > 0 sobre C0, para comprobar si la entropía destruye la fase;
  - reproducibilidad.
- Propuesta C1 (no ejecutada): añadir estructura de segundo orden para romper la degeneración de C0-T2.
