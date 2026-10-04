# Ω-1.1 — Dictamen del panel (fuente de requisitos de Ω-1.1)

Transcripción condensada y fiel del dictamen entregado por el usuario tras los resultados de Ω-1.0 (numeración original). Es la fuente de verdad de los requisitos de Ω-1.1, junto con `docs/OMEGA_MASTER.md`.

**Regla rectora:** no se reemplaza el simulador. Ω-1.0 (WP1–WP4, commit `9dfbea7`, tag `OMEGA_1_0_BASELINE`) se congela como control histórico y Ω-1.1 se construye alrededor. Así se podrá decir: «el modelo original produjo X; la modificación Y produjo Z», y no «modificamos el modelo hasta que funcionó».

**Regla escrita en el código:** D≈3 jamás será suficiente por sí solo para declarar geometría emergente.

## 1. Dictamen sobre Ω-1.0
- 233 pruebas rápidas pasan.
- mypy strict limpio.
- Gradiente validado.
- Experimento 0 valida la medición de dimensión.
- La red aleatoria produce D_eff espurio.
- La dinámica colapsa a vacío o a grafo completo, y esto no es solo un problema numérico.

Subfamilia homogénea W_ij = w:
- T = C(N,3)·w³, S_dens = C(N,2)·w², S_deg = 0.
- H(w) = −α·C(N,3)·w³ + β·C(N,2)·w².
- Punto estacionario interior w* = 2β/(α(N−2)): es un MÁXIMO, no un mínimo.
- W=1 queda favorecido en la frontera cuando α(N−2) > 2β. Con β=1: α_c ≈ 0.0204 (N=100) y ≈ 0.0101 (N=200).

Por tanto, el colapso de la malla con α/β ≥ 0.1 al grafo completo era esperable matemáticamente. γ no puede arreglarlo en la trayectoria homogénea, porque ahí S_deg = 0 y su gradiente desaparece. El papel de γ como estabilizador geométrico es limitado.

> Nota de la sesión (verificada numéricamente): sobre la variedad homogénea hay tres umbrales distintos, todos en α̂ = α(N−2)/(2β):
> - α̂ = 1: W=1 pasa a ser localmente estable (el umbral del panel).
> - α̂ = 3/2: H(1) < H(0), es decir, el grafo completo es globalmente preferido.
> - α̂ = 2: frontera de cuenca desde U(0,1) con ⟨W⟩ = 1/2, porque w* = 1/α̂. Es el umbral medido en Ω-1.0.
>
> Con N=200 y β=1 dan α = 0.0101, 0.0152 y 0.0202. Los tests analíticos deben distinguirlos.

## 2. Ω-MATH: el resultado se congela como teorema del baseline
- `tests/analytical/`: `test_empty_state.py`, `test_complete_state.py`, `test_uniform_manifold.py` y `test_parameter_scaling.py`.
- W=0 es punto estacionario. W=1 es estado frontera cuya estabilidad depende de α(N−2) − 2β.
- `expected_baseline/` contiene `empty_phase`, `complete_phase` y `uniform_barrier`.
- Ω-1.0 ya tiene un comportamiento matemáticamente caracterizable. No se corrige en silencio.

## 3. Ω-REDTEAM: sesgo del término de triángulos
−αT recompensa los triángulos. Triangulación, clustering y 3D no son equivalentes:
- la red cúbica 3D tiene cuadrados pero ningún triángulo;
- una clique tiene una cantidad enorme de triángulos.

Maximizar T favorece la secuencia clustering → hiperdensidad. T es un mecanismo dinámico experimental, no la hipótesis de que la geometría sea triangular, así que la ablación de T es crucial.

## 4. Ω-NET: distinguir geometría de small-world
Suite de distancias `DistanceMode`: HOP, WEIGHTED_INVERSE, WEIGHTED_LOG, RESISTANCE. Cada simulación produce D_hop, D_inv, D_log y D_resistance. Se busca una señal geométrica que sobreviva a cambios razonables de métrica. Si D≈3 aparece solo con ℓ = 1/W, el resultado se marca `distance_sensitive` y no `geometric_phase`.

## 5. El umbral como filtración
G_θ = (V, E_θ) para θ_1…θ_m. Se estudian D(θ), β_i(θ), G(θ) y C(θ). El umbral pasa a ser variable experimental. Módulos `omega/topology/`: `filtration.py`, `betti.py`, `persistence.py`. Referencia: huellas topológicas por números de Betti bajo coarse-graining (EPJC 2026, doi 10.1140/epjc/s10052-026-15322-x).

## 6. Ω-GEOMETRY: tercer estimador independiente
Se conservan D_H (volumen) y D_s. Se añade D_Weyl: con los autovalores del Laplaciano 0 = λ0 ≤ λ1 ≤ …, la ley de Weyl da N(Λ) ∝ Λ^{D/2}, y D_Weyl = 2·d ln N(Λ) / d ln Λ. D_eff, D_s y D_Weyl son mediciones parcialmente independientes (cf. CDT, PRD 100, 026014).

## 7. Certificado 3D Ω
Ningún código devuelve «D = 3.04 → 3D candidate». Existe `GeometryCertificate` con los campos connected, nontrivial, locality, d_volume, d_spectral, d_weyl, metric_robust, seed_robust, size_robust, null_separated, isotropy_ok, homogeneity_ok, topology_stable y manifold_proxy_ok. Solo si se cumplen todos: GEOMETRIC_CANDIDATE. Nunca `D3 = True`.

## 8. Ω-STAT: separar mínimos de fases
dW/dτ = −∇H dice qué mínimos encuentra el algoritmo; no demuestra fases termodinámicas. Por eso hay dos motores:
- `dynamics/` (gradiente, evolución);
- `statistics/` (`metropolis.py`, `langevin.py`, `ensemble.py`), con dW = −∇H dτ + √(2Θ) dB, donde Θ es temperatura estadística y no tiempo.

Pregunta nueva: ¿existe una región estadísticamente dominante con propiedades geométricas?

## 9. Ω-STAT: rama de densidad controlada
- **Ω-A, densidad espontánea:** H0 exacto. ¿Qué densidad genera Ω por sí sola?
- **Ω-B, geometría condicionada:** se fija ρ = 2/(N(N−1))·Σ_{i<j} W_ij. Dada una estructura no trivial, ¿puede la dinámica organizarla geométricamente?

Separa «emergió geometría» de «aumentó la cantidad de conexiones».

## 10. Ω-REDTEAM: controles nulos obligatorios
`omega/controls/`: `erdos_renyi.py`, `configuration_model.py`, `degree_preserving_rewire.py`, `small_world.py`, `random_geometric.py`. El RGG es un control positivo: tiene geometría inicial y comprueba que los detectores la identifican.

## 11. Isotropía y homogeneidad
D_i(r) por nodo, con μ_D(r) y σ_D(r). Un candidato debe cumplir σ_D/μ_D → 0 en el régimen macroscópico, no solo ⟨D⟩ ≈ 3. Isotropía: comparar el crecimiento local de vecindades entre nodos y regiones. Un patrón «A→3, B→7, C→1» no es homogéneo.

## 12. Ω-TOPO (WP5): manifoldicidad
- Se usa `MANIFOLD_PROXY`, no «manifold».
- Se miden β0, β1 y β2 a distintas escalas y χ = β0 − β1 + β2 − ….
- Se añade curvatura discreta. Opción pertinente: Quantum Ricci Curvature (PRD 97, 046008).
- Topología, dimensión y curvatura forman el detector.

## 13. Ω-CURVATURE: curvatura sin coordenadas
Rama paralela W → d → curvatura discreta, y solo después g_μν. Evita que un embedding 3D fabrique la geometría.

## 14. Ω-COARSE
Ω_N → Ω_{N/2} → Ω_{N/4} … con una regla fijada ANTES del resultado. Se miden D, β_i, C, R y homogeneidad. Se busca que micro, meso y macro den la misma clase geométrica.

## 15. Ω-FINITE-SIZE
N = 64, 100, 150, 200, 300, 500, 800, … según coste. Se grafican D(N), C(N), G(N)/N, ρ(N) y ξ(N). Pregunta: ¿sobrevive la fase cuando N→∞?

## 16. Ω-STAT: transición real
Con una variable de control λ: d⟨O⟩/dλ, susceptibilidad χ_O = ⟨O²⟩ − ⟨O⟩², histéresis y distribución del parámetro de orden.

## 17. Ω-COMPUTE: pasaporte ampliado
model_version, code_commit, git_branch, python_version, dependency_lock_hash, config_hash, initialization_distribution, distance_definition, threshold_rule, dimension_algorithm, coarse_graining_algorithm, null_model, random_seed y termination_reason. Cada corrida se identifica como «Ω-EXP-004821» y se reconstruye exactamente.

## 18. Ω-QG: prohibición
Nada de Ŵ ni |Ω⟩ todavía. Se mantiene un modelo clásico/estadístico de pregeometría.

## 19. Ω-GR: preguntas separadas
- Q1: ¿puede emerger una geometría ESPACIAL no trivial sin geometría inicial?
- Q2: ¿métrica efectiva?
- Q3: ¿causalidad?
- Q4: ¿GR?

## 20. Jerarquía experimental
| Nivel | Contenido |
|---|---|
| Ω-0 | Validación matemática |
| Ω-1 | Baseline S0 |
| Ω-2 | Falsación de distancia |
| Ω-3 | Modelos nulos |
| Ω-4 | Ablaciones |
| Ω-5 | Ensamble estadístico |
| Ω-6 | Finite-size scaling |
| Ω-7 | Dimensionalidad múltiple |
| Ω-8 | Topología |
| Ω-9 | Manifold proxy |
| Ω-10 | Curvatura |
| Ω-11 | Coarse-graining |
| Ω-12 | Geometría efectiva |

Más adelante: Ω-13 causalidad, Ω-14 Lorentzianidad, Ω-15 límite GR, Ω-16 fenomenología.

## 21. Qué no tocar
WP1–WP4 no se cambian en silencio: son `OMEGA_1_0_BASELINE`. `OMEGA_1_1` se construye encima.

## 22. Orden exacto antes de los `slow`
1. Tests analíticos (w*, α_c).
2. Congelar el baseline.
3. Distance suite.
4. Null suite.
5. Ablación completa, especialmente T=0.
6. Solo entonces, Experimento 3.
7. Finite-size.
8. Ensamble estadístico.
9. Topología y curvatura.
10. Coarse-graining.

## 23. Proposición P1 (a demostrar o destruir)
Existe una región robusta del espacio de parámetros donde una dinámica relacional sin coordenadas ni geometría iniciales produce una estructura no trivial, conectada, local, aproximadamente homogénea e isotrópica, cuya dimensionalidad efectiva converge hacia una clase tridimensional y cuya señal sobrevive a cambios razonables de métrica, tamaño, semilla, algoritmo y modelo nulo.

## 24. Taxonomía de fracaso
| Código | Significado |
|---|---|
| Ω-F0 | trivial-empty |
| Ω-F1 | trivial-complete |
| Ω-F2 | fragmented |
| Ω-F3 | small-world |
| Ω-F4 | dimension-artifact |
| Ω-F5 | metric-dependent |
| Ω-F6 | finite-size artifact |
| Ω-F7 | seed-dependent |
| Ω-F8 | null-reproducible |
| Ω-F9 | non-manifold |
| Ω-F10 | no stable phase |

Solo entonces: Ω-CANDIDATE. El sistema nunca convierte un número interesante en una «victoria».

## 25. Contexto bibliográfico
- Quantum Graphity (PRD 77, 104029): grafos dinámicos sin geometría de fondo.
- Modelos estadísticos de grafos con estructuras tipo variedad, donde la valencia afecta fuertemente la dimensionalidad (PRD 87, 084011).

La prueba de originalidad es: ¿qué produce específicamente esta dinámica Ω que no sea una propiedad genérica de redes dinámicas conocidas?

## 26. Meta de Ω-1.1
Ω-1.1 no intenta demostrar la teoría. Responde: ¿existe realmente una clase geométrica emergente en nuestro modelo, o solamente estamos midiendo artefactos de redes?
