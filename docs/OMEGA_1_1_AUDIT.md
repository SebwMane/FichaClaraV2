# Ω-1.1 — Auditoría de la implementación (OMEGA_1_1_AUDIT)

Auditor: cerebro orquestador (Opus). Auditado en HEAD `4b036d7` (rama `claude/hopeful-galileo-88u1l1`). No modificó el repo; las pruebas se ejecutaron en el scratchpad. Guardado por la sesión principal.

## 1. Veredicto global
**Todavía no se pueden ejecutar los pasos full. Bastan parches pequeños y localizados.**

Lo que está bien:
- Ω-1.0 sigue intacto.
- El certificado resiste los ataques adversariales.
- Los controles de evidencia se comportan como se diseñó.
- El pasaporte es reconstruible.

Lo que bloquea:
- **B1 y B2:** varios experimentos calculan la curvatura de Ollivier completa sobre estados F1 densos. Eso añade unas 30 h de CPU inútil.
- **B3 y B4:** el control positivo R1 con N=3000 falla por construcción.
- **B5 y B6:** la semántica de los códigos por punto es errónea, y una expectativa de s06 es falsa.
- **B7 y B8:** la compuerta no cubre los pasos 1–2 ni comprueba que el commit coincida, y solo la réplica 0 recibe pasaporte.

### Verificaciones
- **(a) Ω-1.0 intacto.** El diff contra 9dfbea7 solo cambia la versión. Hay 55 sha256 en el manifiesto y el golden pasa. Pasan 570 tests y mypy queda limpio en 87 archivos.
- **(b) Certificado adversarial.** Se partió de 10 corridas sintéticas que lo cumplen todo, con D_vol = D_s = D_W = 3, y se falsificó un solo campo cada vez (22 casos: connected, densidad, uniforme, locality, small-world, plateau, d_s, d_weyl, distance_sensitive, ζ, ζ=NaN, isotropía, cv, cv=NaN, topología, anillos, κ̄=NaN, b1, no convergido, clase None, D*=NaN, G=NaN).
  - Los 22 casos dan NOT_CANDIDATE.
  - También dan NOT_CANDIDATE, a nivel de punto: tener 1 o 2 tamaños, un nulo ausente, un nulo que pasa, un nulo con el mismo D*, deriva en N, o una sola semilla.
  - Con 8/10 corridas limpias hay candidato; con 7/10, no.
  - Un NaN nunca se convierte en True.
- **(c) Compuerta.** s06 en modo full lanza PrerequisiteError antes de escribir nada, y un resumen smoke no cuenta. Pero la compuerta no cubre los pasos 1–2 ni el commit, la config ni un árbol limpio (ver B7).
- **(d) Controles medidos por la ruta real:**

  | Control | Resultado |
  |---|---|
  | RGG3 N=800 k12, 3 semillas | PASS (D*=2.89 / 2.93 / 2.87, ζ=0.30, iso p10 0.52–0.59, κ̄≈−0.02; unos 10 s por corrida) |
  | WS p=0.05 | F3 (también F4 y F9) |
  | ER k12 | F3 (también F4, F5 y F9) |
  | Árbol | F3 |
  | K_50 | F1 |
  | Vacío | F0 |
  | RGG2 | F9 |
  | RGG3 N=1500 | PASS |
  | RGG3 N=3000 | F5 + F9 (ver B3 y B4) |

- **(e) Revisión del código.**
  - No hay RNG global.
  - No se muta ninguna entrada: probado con `writeable=False`.
  - Θ = Θ̂β es coherente.
  - dt = dt_safety/L.
  - Las ventanas son inclusivas y correctas.
  - Los umbrales ≥ y > coinciden con §1.11–1.12.
- **(f) Pasaporte.** Tiene los 14 campos de P§17, y `reconstruct` reproduce los digests (ver B8 y B15 sobre su cobertura).

## 2. Desviaciones reportadas: decisiones
| # | Desviación | Decisión |
|---|---|---|
| 1 | β1^(4)(T³ 4³)=0 y κ=1/6 | **Aceptar.** El error estaba en el diseño: con lado 4, los ciclos que dan la vuelta al toro miden 4 y quedan rellenados. Los tests con lado ≥6 son correctos. Hay que corregir la tabla del diseño. |
| 2 | Langevin tiene sesgo O(dt); ¿usar MALA? | **Aceptar sin MALA en Ω-1.1.** Metropolis es el único muestreador de Gibbs de referencia (medido: −0.25 y +0.51 SE). Langevin es la SDE del panel a dt finito y se informa con su sesgo (+21.8 y +33.6 SE; con dt/2 baja a la mitad; con Richardson queda en ~1–2 SE). MALA pasa a Ω-1.2. Condiciones: B16 y B17. |
| 3 | Configuration model con reparación; el greedy no da N/2 exacto; make_null WS p=0.1 y RGG T² | **Aceptar** con dos condiciones: informar n'/n por nivel, y que `null_separated` rechace RANDOM_GEOMETRIC como nulo (es un control positivo; prioridad baja). |
| 4 | seed_robust ≥2; size_robust exige TODOS los tamaños; salvaguarda F4 | seed_robust ≥2 y la salvaguarda F4: **aceptar**. size_robust con "todos los tamaños": **corregir** (B9). |
| 5 | E3 con diferencias hacia adelante; tests analíticos lentos | **Aceptar.** Es exacto, porque con una sola arista (W²)_ab=0. Los pasos 1–2 tardan 66 s. |
| 6 | Fuente única de resistance_consistent; isotropía p10 de RGG3 0.54–0.60 | **Aceptar.** |
| 7 | ER k12 con primario F3; sigmoide acotada; O-04 solo con N=100; s06 con un solo tamaño; memoización | F3 por precedencia, sigmoide acotada y memoización: **aceptar**. O-04: **corregir**, añadiendo N=200 para FULL y NO_TRIANGLES. s06 con un solo tamaño: **aceptar** junto con B5. |
| 8a | `evidence_cfg` (25 aristas de Ollivier si ρ_bin > 0.1) | **Corregir el criterio y aplicarlo de forma uniforme** (B1 y B2). Con 25 aristas, la SE de la fracción de la cola (~0.034) es mayor que el umbral (0.03), y el rango 0.1 < ρ < 0.5 contiene estados no triviales. |
| 8b | Correcciones al preregistro en O-05 y O-08 (también O-06 "N≥64" y O-02) | **No violan M§44**: no tocan el modelo ni los umbrales de decisión, y ocurrieron en smoke antes de cualquier corrida full. Deben registrarse como **enmiendas** (ver `docs/OMEGA_1_1_AMENDMENTS.md`). `GIBBS_SE_FLOOR=0.005` sí es un ensanchamiento no autorizado (B17). |

## 3. Bugs
- **B1 · Alta (coste).** O-00, O-01, O-03, O-04 y s06 no usan `evidence_cfg`.
  - Medido: K_200 tarda 357 s (10 s con `evidence_cfg`); U(0,1) con N=200 tarda 363 s.
  - Archivos: `test_s06_phase_diagram_v11.py:99,118`, `test_o01_baseline.py:64`, `test_o04_ablation.py:86,118`, `test_o03_nulls.py:120`, `test_o00_validation_v11.py:121`.
  - Parche: `c_ev = evidence_cfg(cfg_caso, w)`, que se pasa a `collect_run_evidence` y a `assessment_row`.
- **B2 · Media.** `finite_size.py:47-52`: el recorte solo debe aplicarse cuando F1 está garantizado. Parche:
  ```python
  v = w[iu]
  dense = bool(np.mean(v > cfg.base.graph.w_min) >= cfg.certificate.dense_rho
               or float(v.mean()) >= cfg.certificate.dense_meanw)
  ```
  Retirar `DENSE_RHO` de `__all__`.
- **B3 · Alta.** `settings11.py:130`: `resistance_max_nodes: int = 2000` → 4000. Con N=3000, RESISTANCE se omite y aparece F5. Es un tope de coste, no un umbral (la pinv de N=3000 tarda ~15 s y ocupa 72 MB; ζ=0.24).
- **B4 · Alta (requiere la decisión D1 del usuario).**
  - Problema: la conectividad de anillos en r = r_hi se fragmenta por saturación del toro. Con RGG3 N=3000 da 0.80–0.875 (F9), mientras que T³ 15³ da 1.0.
  - Parche en `local_structure.py:192-199`: `r_top = r_hi - 1 if r_hi > cfg.annulus_r_min else r_hi`.
  - Con ese parche, las 3 semillas de N=3000 dan ≥0.906. WS, árbol y anillo siguen en 0.
- **B5 · Media.** `certificate.py:261-276` añade F6/F7/F8 aunque el resultado modal sea un fallo: un punto 100% F9 reporta primario F7.
  - Parche: añadir F6/F7/F8 solo si el modal es PASS.
  - Si `count/n < frac`, el código es F7; si el modal es un fallo, se añade ese código.
  - El veredicto no cambia.
- **B6 · Media.** En s06, el peso de F0 es ~1e-10, lo que da χ~4e-21 y un λ* espurio. Parche: `CHI_ZERO_TOL = 1e-12` en `test_s06_phase_diagram_v11.py:148,160`. Registrarlo como enmienda.
- **B7 · Media.** En `gate.py`, `STEP_ORDER` no incluye los pasos 1–2, y un resumen viejo, de otro commit o de un árbol sucio satisface la compuerta. Parche:
  - (i) Añadir `"p01_analytical"` y `"p02_golden"` al inicio de `STEP_ORDER`, con un nuevo `omega/experiments/v11/test_p00_prereq.py`. Su `run()` ejecuta `pytest.main([-q, -p no:cacheprovider, -m "slow or not slow", tests/analytical, tests/test_baseline_golden.py])` y escribe ambos resúmenes con `complete = (rc == 0)`.
  - (ii) `write_summary` guarda `code_commit`, `git_dirty` y `config_hash`.
  - (iii) `require_prerequisites` exige `code_commit` == HEAD y que el árbol no esté sucio.
  - (iv) Un `run(mode="full")` con el árbol sucio se rechaza.
- **B8 · Media (P§17).** Solo se escribe pasaporte para `rep == 0`, en 10 sitios. Parche: quitar `if rep == 0:` en las corridas dinámicas de o01, o04, o06 y s06.
- **B9 · Media.** `size_robust` (`certificate.py:100-140`) debe usar el sufijo contiguo de tamaños limpios (fracción ≥ seed_fraction) que incluye al mayor, con longitud ≥ min_sizes. Registrar `size_suffix_min` en la evidencia.
- **B10 · Baja-media.** En `null_separated` (`certificate.py:174`): `need = max(c.null_sigma * sigma, c.dim_tol)`.
- **B11 · Baja.** En `ensemble.py:144-172`, una cadena congelada da `equilibrated=True`. Parche: si alguna serie primaria tiene varianza 0, `rhat_ok = False`.
- **B12 · Baja.** En `ollivier.py:116-118`, sin aristas se devuelve mean y tail 0; deben ser NaN (ok=False se mantiene).
- **B13 · Baja.** En `provenance.py:132-133`, las fórmulas deben ser `"d=1/(W+ε)"` y `"d=1-ln(max(W,log_floor))"`.
- **B14 · Baja (enmienda recomendada, pendiente).** `homogeneity` descarta en silencio los nodos con S_i(r)=0. Propuesta: exigir n_valid ≥ 0.9·n_gigante. Requiere un campo nuevo de configuración y una versión nueva del diseño.
- **B15 · Baja.** En O-04, guardar `ablation` en `results`. A futuro, un `replay(passport) -> w_final` por corrida.
- **B16 · Baja.** En O-05 (`test_o05_ensemble.py:293`), una región dominante requiere que la celda Metropolis homóloga también lo sea. Si no, se marca `langevin_only_unconfirmed`.
- **B17 · Baja.** Eliminar `GIBBS_SE_FLOOR = 0.005` (`test_o05_ensemble.py:75,304`).
- **B18 · Rendimiento.** `fiedler_length` recalcula un espectro que `weyl_dimension` ya calculó.

## 4. Orden y coste (4 núcleos)
| Paso | Contenido | Tiempo con B1 | Tiempo sin B1 |
|---|---|---|---|
| 0 | Parches, suite rápida (~6 min) y commit limpio | — | — |
| 1–2 | Analítica y golden | 66 s | — |
| — | O-00 full | ~10 min | — |
| — | O-01 full | ~15 min | ~6 h |
| 3 | O-02 | 20–40 min | — |
| 4 | O-03 | 30–45 min | — |
| 5 | O-04 | 30–40 min | — |
| **1–5** | **Total** | **≈2–2.5 h** | **≈13 h** |
| 6 | s06 | ≈1 h | ~20 h |

- **O-06:** medir primero una corrida Ω-B con N=800 antes de lanzar la malla completa.
- **O-05:** ≈18–20 h.

## 5. Decisiones del usuario
- **D1:** la regla de anillos de B4, es decir, evaluar r ∈ [2, r_hi−1] en lugar de [2, r_hi].
  - Relaja el criterio solo en el radio contaminado por la saturación del toro.
  - La alternativa es mantener la regla actual y aceptar que el certificado 3D es inalcanzable con tamaños XL.
- **D2:** confirmar R1, R3, R4 y R8 del §6 del diseño antes de s06.
