# ANALYSIS.md — Verificación, riesgos, decisiones de diseño y contratos del simulador Ω (v1)

Fuente de verdad: `docs/OMEGA_MASTER.md` (**M§n**). `docs/SPEC_EXTRACT.md` (**S§n**) se usó como apoyo y se corrige donde discrepa (§2). Las cifras numéricas proceden de scripts de prueba (Python 3.11, numpy/scipy/networkx, semillas `SeedSequence` explícitas) ejecutados fuera del repo. Autor del análisis: agente Opus; guardado por la sesión principal.

## 0. Resumen ejecutivo
1. **El gradiente de M§17 es EXACTO** con la convención "una variable por arista no ordenada a<b". Al derivar k̄ respecto a W aparece un término −(4/N)Σ_i(k_i−k̄), que es idénticamente 0. Diferencias finitas: error relativo máximo ≈ 2·10⁻⁹.
2. **Hallazgo crítico (preliminar):** con η=μ=0 y W₀~U(0,1), la dinámica es biestable y termina en **W≡0 (A)** o en **W≡J−I (grafo completo, E)**. El umbral es α_c ≈ 4β/(N−2), es decir α̂ := α(N−2)/(2β) ≈ 2. Con N=100–300, toda la malla de M§25 con α/β ≥ 0.1 cae en E y α=0 cae en A; γ solo cambia la velocidad. Es el riesgo de M§20: debe registrarse como resultado negativo. Barrer en (α̂, γ̂) (D-15). Preliminar (N=100 y 200, 3 semillas): confirmar en el Experimento 2.
3. **Hallazgo crítico (medición):** en la línea base aleatoria, un ajuste log-log ingenuo de N(r) da **D ≈ 2.8–4.2** según W_min. Es un artefacto, porque las distancias abarcan solo un rango ~2×. Hace falta una ventana mínima de escala y un estado "no definido".
4. El estimador "bola" (ln N vs ln r) está muy sesgado a N≈300 (toro 3D 7³: **1.98**). El estimador "cáscara" (D = 1 + dlnS/dlnr) da **2.84**. Para D_s: paseo *lazy* 3.12; núcleo de calor con Laplaciano normalizado, sobreimpulso de 3.7–4.3. Por defecto: cáscara y lazy.

## 1. Verificación matemática (M§13, M§17)

### 1.1 Convención
Las variables son las M=N(N−1)/2 aristas w_ab con a<b, de modo que W_ab=W_ba=w_ab y W_aa=0. Las sumas de M§13 recorren i<j e i<j<k. El gradiente se representa como una matriz simétrica G con diagonal 0, donde G_ab=∂S0/∂w_ab.

### 1.2 Derivación
- **T** = Σ_{i<j<k} W_ij W_jk W_ki = tr(W³)/6 (verificado contra el triple bucle: 1e-14). ∂T/∂w_ab = Σ_{k≠a,b} W_ak W_kb = (W²)_ab, porque W_aa=W_bb=0. Exacto.
- **S_dens**: ∂/∂w_ab = 2W_ab. Exacto.
- **S_deg (derivada completa)**: ∂k_i/∂w_ab = δ_ia+δ_ib y k̄ = (2/N)Σ_{i<j}W_ij, luego ∂k̄/∂w_ab = 2/N.
  ∂S_deg/∂w_ab = Σ_i 2(k_i−k̄)(δ_ia+δ_ib−2/N) = 2[(k_a−k̄)+(k_b−k̄)] − (4/N)Σ_i(k_i−k̄) = 2[(k_a−k̄)+(k_b−k̄)], porque Σ(k_i−k̄)≡0.
  **EXACTO.** Congelar k̄ da el mismo resultado (verificado: 1.7e-9). Se dejaría de cumplir con la variante ((k_i−k̄)/k̄)² de S§7.14, que NO se adopta.
- **S_smooth (η, apagado en v1)**: C_ij = s_i+s_j−2(W²)_ij con s_i=(W²)_ii, de donde Σ_{i,j}W_ij C_ij = 2Σ_i k_i s_i − 2tr(W³), con la suma sobre pares ordenados como en M§15. ∂/∂w_ab = η[2(s_a+s_b) + 4W_ab(k_a+k_b) − 12(W²)_ab]. Verificado (1.7e-9).
- **S_mass** = −μΣ_{i<j}W_ij, con derivada −μ.
- **Sigmoide**: ∂S0/∂θ_ab = G_ab·W_ab(1−W_ab).

### 1.3 Convención i<j frente a pares ordenados
Implementar Σ_{i≠j}W² (=2S_dens) o tr W³ (=6T) cambia los coeficientes efectivos por factores 2 y 6, lo que desplaza el mapa de fases y el Δτ efectivo. Derivar respecto a una sola entrada W_ab rompe la simetría y da un valor sin sentido (0.463 frente a −0.312 de la fórmula de arista). **Regla:** funcional literal (i<j, i<j<k) y gradiente por arista.

### 1.4 Verificación numérica
N=12, (α,β,γ,η)=(0.7,1,0.3,0.2), diferencias centrales con h=1e-6 sobre las 66 aristas.

| Término | Error absoluto máx. | Error relativo |
|---|---|---|
| T | 5.1e-9 | 1.4e-9 |
| S_dens | 1.5e-9 | 7.5e-10 |
| S_deg | 1.8e-9 | 3.7e-10 |
| S0 | 3.4e-9 | 1.6e-9 |
| S0 + ηS_smooth | 1.3e-8 | 1.6e-9 |

### 1.5 Fórmula que el código DEBE usar
```python
k=W.sum(1); u=k-k.mean(); W2=W@W
G=-alpha*W2+2*beta*W+2*gamma*(u[:,None]+u[None,:])
G+=eta*(2*(s[:,None]+s[None,:])+4*W*(k[:,None]+k[None,:])-12*W2)  # s=diag(W2); solo si eta!=0
G+=-mu                                                              # solo si mu!=0
np.fill_diagonal(G,0.0)
```
Actualización: W ← clip(W−ΔτG, 0, 1) seguida de diag=0. Es descenso proyectado sobre [0,1]^M.

### 1.6 Hessiano y Δτ
H = 2βI + 2γBᵀ(I−11ᵀ/N)B − αH_T, con B la matriz de incidencia y (H_T)_{ab,ac}=W_bc.
- λ_max(Bᵀ(I−11ᵀ/N)B) = N−2; verificado para N=8, 12 y 20.
- ‖H_T‖ ≤ 2(N−2).
- **L = 2β + 2(γ+|α|)(N−2)** (η=0). La estabilidad exige Δτ < 2/L; por defecto Δτ = 0.5/L.

### 1.7 Equivariancia
S0(PWPᵀ)=S0(W) y G(PWPᵀ)=P G Pᵀ (prueba obligatoria, M§32).

## 2. Discrepancias de SPEC_EXTRACT (corregidas)
- **S2.6:** el maestro dice arista "enorme" en *longitud* (peso pequeño), no "de gran peso".
- **S2.20:** M§29 pide probar distribuciones distintas como robustez. U(0,1) es la condición por defecto, no una prohibición.
- **Exp 0:** el umbral "fallo > 0.5" es inventado (ver §6.3).
- **Exp 1:** "D_eff≈∞" y "fallo si D≈3" no están en el maestro y son numéricamente falsos con un estimador ingenuo (2.8–4.2). Lo correcto es "sin ventana". La afirmación ER sí es correcta: A ~ G(N, 1−W_min).
- **Exp 2:** "W fuera de [0,1]" es imposible con clip; el fallo real son los valores no finitos. Calcular observables en cada paso es caro (D-4).
- **4.1 y 4.4:** los umbrales 0.3 y 0.1 son inventados. Bajo permutación la igualdad debe ser exacta (≤1e-9).
- **4.3:** la histéresis es evidencia posible, no un criterio de éxito.
- **4.5:** M§33 lista cuatro combinaciones, todas con −αT. "Solo β" y "solo γ" son controles extra; "solo β" es trivial (mínimo único en 0).
- **§5:** los "parámetros típicos", G∈[0.1,0.5] y k̄≈N−1 no están en M§19. "A: α grande" contradice la numérica (α grande lleva a E).
- **F con las 12 propiedades por corrida:** mezcla M§19 con M§46. Se separa en F-candidata (por corrida) y F-confirmada (por punto).
- **G:** M§19-C usa G/N, así que G es un tamaño; se guardan tamaño y fracción.
- **`update_weight` mutante y `save_weights` en network/:** violan la pureza; la persistencia pasa a io/.
- **Simetrización (X+Xᵀ)/2:** produce una distribución triangular, no U(0,1). Se usa la triangular superior espejada.
- **Huecos omitidos:** μ/W0; Δτ con η>0; el nombre `io`; pytest recolecta `experiments/test_*`.

## 3. Riesgos científicos

**R1. Fase trivial y biestabilidad.**
- Con α=0, S0 es estrictamente convexa y ≥0=S0(0): W=0 es el mínimo global único para todo γ≥0.
- Con α>0, W=0 es mínimo local (Hessiano en 0 ≻ 0).
- J−I es punto KKT atractor si α̂ ≥ 1.
- El estado uniforme w* = 1/α̂ es un punto silla que separa ambos destinos. Con ⟨W₀⟩=1/2, el umbral es α̂_c ≈ 2.
- Numérica (N=100 y 200; γ̂∈{0,10,100}; 3 semillas): 100% W≡0 con α̂≤1.98, 100% W≡1 con α̂≥2.02; transición dependiente de la semilla en [1.99, 2.01]. No hubo estados mixtos.
- T~N³ y S_dens~N² no escalan igual, así que α crudo no es comparable entre N. **Mitigación:** parámetros reducidos (D-15), registrar el resultado negativo, y abrir la rama M§20 de forma preregistrada.

**R2. Umbral W_min.**
- G, k̄_bin y C_bin dependen arbitrariamente de W_min.
- Las distancias apenas: una arista con W≲0.5 nunca es geodésica (dos saltos fuertes ≈2 < 1/W); solo el 52% de las aristas son geodésicas con W_min=0.
- Con clip, el estado final está en {0,1}.
- **Mitigación:** observables ponderados (k, C_W) y barrido obligatorio W_min∈{0.01,0.05,0.1,0.2,0.5}; una etiqueta solo se acepta si es estable en ≥4/5.

**R3. Tamaño finito en D_eff.**
- Un toro 3D con N≈343 solo da r≤3 útil, menos de media década.
- Bola: 0.93 / 1.63 / 1.98. Cáscara: 1.00 / 2.00 / 2.84. Cáscara con 15³: 2.92. Toro abierto 7³: 2.40.
- Con N≤300, D_eff solo sirve para cribar; H2 exige N≥1000–3000 y la curva D_eff(r).

**R4. Tamaño finito en D_s.**
- El núcleo de calor con L_norm en redes bipartitas sobreimpulsa (pendiente local en 3D de 3.7–4.3).
- Restar 1/N hace que la pendiente diverja cerca de 1/λ₂.
- Lazy (q=1/2), con t≥4 y P̄≥5/N, da 0.99 / 2.03 / 3.12 a N≈300.
- La ventana debe cortarse por saturación.

**R5. Sesgo de clip.**
- Crea masas en 0 y 1. Los estados finales son puntos KKT, no estacionarios.
- La convergencia se mide con el gradiente proyectado; se registran las fracciones en 0 y en 1.
- La sigmoide sirve como prueba de robustez del algoritmo, no como sustituto.

**R6. Circularidad de la métrica 1/(W+ε).**
- La longitud mínima es ≈1. En un grafo casi completo, N_i(r) ≈ 1+(N−1)P(W≥1/r): D_eff mide la distribución de pesos y se puede "obtener" cualquier D.
- Explica el D≈3 espurio de la línea base.
- **Mitigación:** razón de escala ≥3 y saturación ≤0.2·N_gc; nulos de pesos barajados y recableado con grados preservados; D_eff también en saltos, contrastado con D_s.

**R7. Red densa o completa.**
- Con ρ≥0.5 el diámetro en saltos es ≤2 y no hay escalas: D no está definida y E se decide por topología.
- W≡1 produce un N(r) en escalón: el estimador debe devolver `no_window` sin fallar.
- Grafo vacío: G=1/N y estado `insufficient_component`.

**R8. Estadística.** Varias semillas; probabilidad de fase con IC; ralentización crítica cerca de α̂_c (el estado MAX_STEPS da etiqueta U).

**R9. Ajuste retrospectivo (M§44).** Congelar umbrales y `schema_version` antes del Experimento 3; cualquier cambio implica nueva versión y re-ejecución.

## 4. Decisiones de diseño (todas en un único `OmegaConfig` frozen y guardadas en el pasaporte)

| # | Hueco | Valor/regla por defecto | Justificación | Campo |
|---|---|---|---|---|
| D-0 | Simetrización | U(0,1) en la triangular superior, espejada; diag 0 | Conserva la marginal | init.symmetrization |
| D-1 | ε | 1e-9 | Solo evita dividir entre 0; error relativo <1e-7 si W>0.01 | graph.epsilon |
| D-2 | W_min | 0.1 + barrido {0.01,0.05,0.1,0.2,0.5} | R2 | graph.w_min, w_min_sensitivity |
| D-3 | Δτ | auto: 0.5/L con L=2β+2(γ+\|α\|)(N−2); fijo obligatorio si η≠0 | §1.6 | dynamics.dt_mode, dt_safety, dt_fixed |
| D-4 | Convergencia | ‖ΔW‖∞<1e-10 durante 10 pasos; max_steps=200000; estado ∈ {CONVERGED, MAX_STEPS, NONFINITE}; τ_max=pasos·Δτ. Escalares en cada paso; instantáneas en {0,1,2,4,…} y final; geometría solo en instantáneas | R5 | dynamics.tol_step, patience, max_steps, snapshot_schedule |
| D-5 | Semillas | Exp1: 10; Exp3: 10 por punto; robustez: 30 | Coste frente a error | seeds.replicates |
| D-6 | Rango de r | Enteros si las distancias lo son (tol. 1e-6); si no, 40 log-espaciados en [r_lo=mediana del vecino más cercano, d_max]. r_hi: N̄(r)≤0.2·N_gc. Se exigen ≥3 puntos y r_hi/r_lo≥3; si no, no_window | Prototipo | dimension.* |
| D-7 | Rango de σ | Tiempos enteros log en [1, 10·N_gc] (200 puntos); ventana t≥4 y P̄≥5/N_gc; razón ≥3 | R4 | spectral.* |
| D-8 | Malla | α̂∈{0,0.25,…,4} + [1.5,2.5] con paso 0.02; γ̂∈{0,0.1,0.3,1,3,10,30,100}; la malla cruda sigue disponible | R1 | scan.* |
| D-9 | N | Exp1–3: N=200; tamaños {100,200,300} y luego {500,1000} | M§27 | init.n |
| D-10 | C | C_W = tr(W³)/Σ_i(k_i²−Σ_j W_ij²) ∈ [0,1] y C_bin (media local sobre A); se informa C_bin/ρ | Forma cerrada, coherente con T | — |
| D-11 | L | Media de d sobre pares de la componente gigante; L en saltos; diámetro aparte | Evita infinitos | — |
| D-12 | P(σ) | Lazy discreto sobre la componente gigante de W⊙A: P=(1−q)I+qD⁻¹(W⊙A), q=1/2. Retorno al nodo exacto, media uniforme: P̄(t)=(1/N)Σ_kμ_kᵗ (eigh de la forma simétrica). Alternativas: calor con L_norm, binario | R4; forma cerrada | spectral.method, laziness, graph |
| D-13 | Sigmoide | Logística; θ₀=logit(clip(W₀,1e-6,1−1e-6)); ∂/∂θ=G⊙W(1−W) | M§16 | dynamics.integrator, sigmoid_clip |
| D-14 | Normalización de S_deg | Ninguna | Fidelidad; gradiente exacto | — |
| D-15 | Normalización de S0 | Funcional literal; experimento en α̂=α(N−2)/(2β) y γ̂=γ(N−2)/β con β=1; se guardan ambos | Reparametrización, no cambio de modelo; w*=1/α̂ | scan.scaling ∈ {reduced, raw} |
| D-16 | Redes de referencia | Anillo y toros 2D/3D periódicos con W=1; abiertas solo informe; RGG en toro; nulos: U(0,1), completo, ER denso, 3-regular | §6.3 | experiments.reference_graphs |
| D-17 | G≈1 | ≥0.9 (C, D), ≥0.95 (F) | Margen | phases.g_connected, g_geometric |
| D-18 | cv bajo | σ_k/k̄ ≤ 0.2 (U(0,1) con N=200 da ≈0.04) | — | phases.cv_homogeneous |
| D-19 | Estabilidad de escala | Desviación estándar de las pendientes locales ≤0.3 ⇒ plateau | R3 | dimension.plateau_tol |
| D-20 | Variaciones pequeñas | ±10% en α̂ y γ̂ | — | phases.param_perturbation |
| D-21 | Permutación | Diferencia absoluta ≤1e-9 | Determinista | — |
| D-22 | Agregación | MCO log-log en la ventana (pendiente, error estándar, R²) + curva local + estimador bola como validación cruzada | M§29 | dimension.estimator |
| D-23/24 | MDS y métrica | Fuera de v1 | M§35 | — |
| D-25 | Gradiente | Analítico; diferencias finitas solo en pruebas (h=1e-6, tol. rel. 1e-6) | §1.4 | — |
| D-26 | N grande | Denso hasta ~2000 | M§27 | — |
| D-28 | Semillas | SeedSequence(entropy, spawn_key=(exp_id, point, rep)) → Generator(PCG64); prohibido np.random global | Independiente del orden | seeds.* |
| D-30 | M§20 | μ=0, sin restricción en v1 | M§20 | functional.mu |

### 4.1 Fases por corrida (preregistradas)

| Fase | Regla |
|---|---|
| E | ρ≥0.5 o ⟨W⟩≥0.5 |
| A | G_frac<0.1 y k̄_bin<1 |
| B | ≥2 componentes con ≥10% de N y G_frac<0.9 |
| F-cand | G_frac≥0.95, cv≤0.2, D_eff ok con plateau, D_s ok, \|D_s−D_eff\|≤0.5, no E |
| D | G_frac≥0.9, C_bin≥0.3 y C_bin/ρ≥3 |
| C | G_frac≥0.9 y cv≤0.2 |

- Precedencia: E>A>B>F>D>C>U. Se devuelven también todas las banderas. Una corrida no convergida recibe U.
- **F-confirmada** (por punto) exige: ≥80% de semillas F-cand; desviación estándar de D_eff ≤0.2; \|ΔD\|≤0.3 entre N∈{100,200,300} sin deriva monótona; estable ante ±10% de parámetros; permutación exacta; estable en ≥4/5 valores de W_min; distinta del nulo con pesos barajados.
- Ablación y coarse-graining se informan aparte.

## 5. Contratos de API

**Reglas globales:**
- Paquete `omega/`. Una función = una tarea. Funciones puras que no mutan sus entradas; sin estado global.
- RNG explícito como `np.random.Generator`. La E/S solo está permitida en io/ y experiments/.
- `mypy --strict`. Las funciones públicas validan y lanzan ValueError/TypeError. Todo es float64.
- `omega/types.py`: `FloatArray=npt.NDArray[np.float64]`, `BoolArray`, `IntArray`. Una "WeightMatrix" es (N,N) con N≥2, finita, simétrica (atol 1e-12), diag 0 y valores en [0,1].

**5.1 config/settings.py** (frozen, slots, `__post_init__` valida)
```python
FunctionalParams(alpha:float, beta:float=1.0, gamma:float=0.0, eta:float=0.0, mu:float=0.0)  # beta>0, gamma>=0, eta>=0
InitConfig(n:int, distribution:Literal["uniform"]="uniform", symmetrization:Literal["upper_mirror"]="upper_mirror")
DynamicsConfig(integrator:Literal["clip","sigmoid"]="clip", dt_mode:Literal["auto","fixed"]="auto", dt_safety:float=0.5,
  dt_fixed:float|None=None, max_steps:int=200_000, tol_step:float=1e-10, patience:int=10,
  snapshot_schedule:Literal["log2","none"]="log2", sigmoid_clip:float=1e-6)
GraphConfig(w_min:float=0.1, epsilon:float=1e-9, w_min_sensitivity:tuple[float,...]=(0.01,0.05,0.1,0.2,0.5))
DimensionConfig(estimator:Literal["shell","ball"]="shell", n_radii:int=40, saturation:float=0.2, min_scale_ratio:float=3.0,
  min_points:int=3, plateau_tol:float=0.3, integer_tol:float=1e-6)
SpectralConfig(method:Literal["lazy_walk","heat_normalized"]="lazy_walk", laziness:float=0.5, t_min:float=4.0,
  saturation_factor:float=5.0, min_scale_ratio:float=3.0, n_times:int=200,
  graph:Literal["thresholded_weighted","weighted","binary"]="thresholded_weighted")
PhaseThresholds(g_dispersed=0.1, kbin_dispersed=1.0, large_component_frac=0.1, g_connected=0.9, g_geometric=0.95,
  cv_homogeneous=0.2, rho_hyperdense=0.5, meanw_hyperdense=0.5, c_clustered=0.3, c_ratio_clustered=3.0, ds_deff_tol=0.5,
  f_seed_fraction=0.8, f_seed_std=0.2, f_size_tol=0.3, param_perturbation=0.1)
SeedConfig(master_entropy:int, experiment_id:int, replicates:int=10)
ScanConfig(alpha_hat:tuple[float,...], gamma_hat:tuple[float,...], scaling:Literal["reduced","raw"]="reduced")
OmegaConfig(init, functional, dynamics, graph, dimension, spectral, phases, seeds, scan:ScanConfig|None=None, schema_version:str="1.0")
  # valida: eta!=0 => dt_mode=="fixed"
```

**5.2 config/convert.py**
```python
config_to_dict(cfg)->dict[str,Any]
config_from_dict(d)->OmegaConfig             # rechaza claves desconocidas y schema distinto
reduced_to_raw(alpha_hat, gamma_hat, n, beta=1.0)->FunctionalParams
raw_to_reduced(p, n)->tuple[float,float]
```

**5.3 config/seeds.py**
```python
SeedKey(entropy:int, spawn_key:tuple[int,...])
seed_key(cfg:SeedConfig, point_index:int, replicate:int)->SeedKey
make_rng(key)->np.random.Generator            # PCG64(SeedSequence(entropy, spawn_key))
```

**5.4 types.py** (dataclasses frozen, sin lógica)
- `RunStatus{CONVERGED, MAX_STEPS, NONFINITE}`; `EstimateStatus = Literal["ok","no_window","insufficient_component"]`; `PhaseLabel{A..F, U}`.
- `DimensionEstimate(value, stderr, status, window, scales, profile, local_slopes, plateau, method)`
- `TopologyObservables(n, mean_strength, std_strength, cv_strength, mean_binary_degree, binary_density, mean_weight, giant_size, giant_fraction, n_components, n_large_components, clustering_weighted, clustering_binary, frac_at_zero, frac_at_one)`
- `GeometryObservables(path_length, path_length_hops, diameter, d_eff, d_eff_ball, d_eff_hops, d_s)`
- `Observables(topology, geometry)`
- `Trajectory(w_final, status, steps, dt, tau, scalars:Mapping[str,FloatArray], snapshot_steps, snapshots[(n_snap, M), triangular superior])`
- `PhaseAssessment(label, flags)`
- `RunResult(seed, params, alpha_hat, gamma_hat, w0, trajectory, observables, assessment)`

**5.5 network/**
- **weights.py**
  ```python
  validate_weight_matrix(w,*,atol=1e-12)->None
  upper_triangle(w)->FloatArray
  from_upper_triangle(v,n)->FloatArray
  clip_unit(w)->FloatArray
  permute(w, perm:IntArray)->FloatArray
  ```
- **initialization.py**
  ```python
  random_uniform_weights(n:int, rng:Generator)->FloatArray
  ```
- **topology.py**
  ```python
  adjacency(w, w_min)->BoolArray
  strength(w)->FloatArray
  binary_degree(a)->IntArray
  component_labels(a)->tuple[int,IntArray]
  component_sizes(labels)->IntArray
  giant_component_nodes(labels)->IntArray      # desempate determinista
  submatrix(w, nodes)->FloatArray
  weighted_clustering(w)->float
  binary_clustering(a)->float
  topology_observables(w, w_min, large_component_frac)->TopologyObservables
  ```

**5.6 dynamics/**
- **functional.py**
  ```python
  triangles(w)->float                           # tr(W³)/6
  density(w)->float
  degree_irregularity(w)->float
  smoothness(w)->float
  mass(w)->float
  action(w, p:FunctionalParams)->float
  ```
- **gradient.py** (derivadas por arista; simétricas con diag 0)
  ```python
  grad_triangles(w); grad_density(w); grad_degree_irregularity(w); grad_smoothness(w); grad_mass(w)
  grad_action(w, p)->FloatArray
  lipschitz_bound(n, p)->float                  # lanza si eta!=0
  theta_gradient(g, w)->FloatArray
  ```
- **evolution.py**
  ```python
  step_size(n, p, cfg)->float
  projected_step(w, g, dt)->FloatArray
  sigmoid(theta); logit(w, delta)
  sigmoid_step(theta, p, dt)->FloatArray
  snapshot_steps(max_steps, schedule)->IntArray
  step_scalars(w, w_prev, p)->dict[str,float]
  is_converged(max_dw_recent, tol, patience)->bool
  evolve(w0, p, cfg:DynamicsConfig)->Trajectory    # no calcula geometría
  ```

**5.7 geometry/**
- **distances.py**
  ```python
  edge_length_matrix(w, a, epsilon)->scipy.sparse.csr_array
  distance_matrix(w, w_min, epsilon)->FloatArray     # Dijkstra; inf entre componentes
  hop_distance_matrix(a)->FloatArray
  ball(d, i, r)->IntArray
  ball_counts(d, radii)->IntArray(n,K)
  path_length(d)->float
  diameter(d)->float
  ```
- **dimension.py**
  ```python
  radius_grid(d, cfg)->FloatArray
  mean_profile(counts)->FloatArray
  scaling_window(scales, profile, upper_limit, min_ratio, min_points)->tuple[int,int]|None
  local_slopes(x, y)->FloatArray
  fit_loglog(x, y)->tuple[slope, stderr, r2]
  shell_profile(scales, profile, integer)->tuple[FloatArray, FloatArray]
  effective_dimension(d, cfg)->DimensionEstimate
  ```
- **spectral.py**
  ```python
  walk_operator_spectrum(w, laziness)->FloatArray
  heat_spectrum(w)->FloatArray
  time_grid(n, cfg)->FloatArray
  mean_return_probability(spectrum, times, method)->FloatArray
  spectral_dimension(w, cfg)->DimensionEstimate
  ```
- **observables.py**
  ```python
  geometry_observables(w, graph, dim, spec)->GeometryObservables   # restringido a la componente gigante
  ```

**5.8 phases/**
- **classification.py**
  ```python
  is_hyperdense / is_dispersed / is_fragmented / is_geometric_candidate / is_clustered / is_connected_homogeneous(o, t)->bool
  phase_flags(o, t)->dict[str,bool]
  classify_run(o, status, t)->PhaseAssessment
  ```
- **scan.py**
  ```python
  parameter_grid(scan)->tuple[tuple[int,float,float],...]
  point_params(alpha_hat, gamma_hat, n, scan, beta)->FunctionalParams
  simulate(cfg, p, key)->RunResult
  scan_phase_diagram(cfg)->tuple[RunResult,...]
  phase_probabilities(results)->dict
  ```
- **stability.py**
  ```python
  seed_summary(results, key)->tuple[mean, std, n]
  permutation_invariance(w, perm, cfg)->dict[str,float]
  shuffled_weight_null(w, w_min, rng)->FloatArray
  degree_preserving_null(w, w_min, n_swaps, rng)->FloatArray
  hysteresis_sweep(w_start, alpha_hat_path, gamma_hat, cfg)->tuple[RunResult,...]   # continuación
  hysteresis_gap(fwd, bwd, key)->FloatArray
  confirm_geometric(point, by_size, neighbors, t)->PhaseAssessment
  ```

**5.9 statistics/summary.py**
```python
bootstrap_ci(values, rng, n_boot=2000, level=0.95)
proportion_ci(k, n, level=0.95)                  # Wilson
```
`coarse_graining/`, `curvature/` y `causality/` llevan solo un `__init__.py` con docstring.

**5.10 io/passport.py** (`omega.io` siempre con imports absolutos)
```python
array_digest(a)->str                             # sha256
software_versions()->dict
build_passport(cfg, result, versions, git_commit)->dict
save_run(result, cfg, directory:Path, run_id:str, git_commit)->tuple[Path,Path]
load_run(json_path)->tuple[dict, dict[str,FloatArray]]   # verifica digests
```
- **JSON:** schema_version, run_id, created_utc, git_commit, versions, config completa, seed{entropy, spawn_key, "PCG64"}, params_raw, params_reduced, N, dt, tau_max, steps, status, omega0{distribution, symmetrization, sha256}, observables (k̄, σ_k, C_W, C_bin, L, G_size, G_frac, D_eff y D_s con error, ventana y estado), phase{label, flags}, curvature:null, npz_file, npz_sha256.
- **NPZ comprimido:** w0_upper, w_final_upper, snapshot_steps, snapshots_upper, scalars_*, deff_scales, deff_profile, ds_times, ds_profile.

**5.11 experiments/**
- `reference_graphs.py`: periodic_lattice, open_lattice, random_geometric_torus, complete_graph, random_regular. Solo para validación; está prohibido importarlo desde dynamics/ y network/.
- `test_00_validation.py`: pytest de §6.3.
- `test_01` a `test_05`: `run(cfg, out_dir)->Path` + prueba de humo; las corridas completas van con `@pytest.mark.slow`.
- `pyproject.toml`: testpaths=["tests","omega/experiments"], addopts="-m 'not slow'".

## 6. Plan de pruebas

### 6.1 Invariantes transversales
Simetría, diag 0, valores finitos y en [0,1]; no mutación de entradas; determinismo bit a bit por SeedKey; invarianza (escalares) y equivariancia (matrices) por permutación.

### 6.2 Por módulo
- **config:** ida y vuelta; valores inválidos rechazados; reducido↔crudo es la identidad (1e-12); semillas distintas producen flujos distintos.
- **initialization:** prueba KS de la marginal U(0,1), p>0.01 con N=200.
- **topology:** grafo vacío G=1/N; completo con G=ρ=C_bin=C_W=1; dos cliques dan n_large=2; C_W∈[0,1].
- **functional:** T por bucle = traza (1e-12); smoothness por bucle = forma cerrada; en el grafo completo T=C(N,3) y S_deg=0.
- **gradient:** frente a diferencias finitas (N=12, tol. rel. 1e-6); equivariancia; λ_max(H_deg)=2γ(N−2).
- **evolution:** S0 no creciente en 200 pasos con Δτ auto; α=0 lleva a W→0; α̂=3 lleva a J−I; el punto uniforme 1/α̂ es estacionario; un NaN produce NONFINITE; longitudes del registro correctas.
- **distances:** valores analíticos en una cadena; desigualdad triangular; inf entre componentes; ball_counts monótono.
- **dimension:** con N=c·r^D (D∈{1, 2.5, 3}) recupera D (1e-6); un escalón da no_window.
- **spectral:** espectro en [1−2q, 1]; P̄ decreciente hacia 1/N; grafo completo da no_window.
- **classification:** Observables sintéticos por fase; precedencia; MAX_STEPS da U.
- **scan/stability:** reproducibilidad; permutación ≤1e-9; los nulos conservan pesos o grados.
- **io:** ida y vuelta idéntica; digest alterado lanza excepción; campos de M§49 presentes.
- **Arquitectura (prueba estática):** ningún uso de np.random global ni importación de reference_graphs fuera de su ámbito.

### 6.3 Aceptación del Experimento 0
D_eff con estimador cáscara sobre saltos (W=1); D_s lazy; ventanas automáticas. Valor del prototipo entre paréntesis.

| Red periódica | D_eff ± tol | D_s ± tol |
|---|---|---|
| anillo 300 | 1±0.05 (1.00) | 1±0.10 (0.99) |
| anillo 1000 | 1±0.05 (1.00) | 1±0.10 (1.00) |
| 17² | 2±0.10 (2.00) | 2±0.15 (2.03) |
| 32² | 2±0.10 (2.00) | 2±0.15 (2.01) |
| 7³ | 3±0.30 (2.84) | 3±0.25 (3.12) |
| 10³ | 3±0.25 (2.87) | 3±0.25 (3.12) |
| 15³ | 3±0.15 (2.92) | 3±0.25 (3.09) |

- **RGG en toro** (métrica euclidiana, N=1000, 4 semillas): 2D 2±0.15 (≈1.96); 3D 3±0.3 (≈2.85).
- **Controles negativos** (deben dar no_window): U(0,1) con N∈{100,300} y W_min∈{0,0.5,0.9} (12/12 en el prototipo) y grafo completo. 3-regular: no_window o plateau=False (pendiente de verificar).
- **Solo informe:** redes abiertas (7³ ≈ 2.40) y estimador bola.
- **Invariancias:** permutación ≤1e-9; independencia de W_min en redes con W=1.
- **Regla de M§22:** si falla una fila obligatoria, no se ejecutan los Experimentos 1–5.

## 7. Paquetes de trabajo (archivos exclusivos; los contratos de §5 quedan congelados)

| WP | Archivos | Depende de |
|---|---|---|
| WP1 Fundamentos | pyproject.toml, omega/__init__.py, omega/types.py, omega/config/{__init__,settings,convert,seeds}.py, omega/io/{__init__,passport}.py, omega/network/{__init__,initialization,weights,topology}.py, omega/statistics/{__init__,summary}.py, omega/{coarse_graining,curvature,causality}/__init__.py, tests/{conftest,test_config,test_seeds,test_io,test_network,test_statistics,test_architecture}.py | — |
| WP2 Dinámica | omega/dynamics/{__init__,functional,gradient,evolution}.py, tests/{test_functional,test_gradient,test_evolution}.py | WP1 (types, config) |
| WP3 Geometría + Exp 0 | omega/geometry/{__init__,distances,dimension,spectral,observables}.py, omega/experiments/{__init__,reference_graphs,test_00_validation}.py, tests/{test_distances,test_dimension,test_spectral}.py | WP1 (types, config, topology) |
| WP4 Fases + experimentos | omega/phases/{__init__,classification,scan,stability}.py, omega/experiments/test_0{1..5}_*.py, tests/{test_classification,test_scan,test_stability}.py | WP1, WP2, WP3 |

**Orden:**
1. WP1 se fusiona primero.
2. WP2 y WP3 en paralelo, programando contra los contratos congelados.
3. WP3 debe pasar el Experimento 0 antes de que WP4 ejecute nada pesado.
4. WP4 integra; sus pruebas de clasificación con datos sintéticos no esperan a nadie.

Ningún archivo aparece en dos paquetes.
