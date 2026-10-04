# Ω-1.1 — Diseño científico y de ingeniería (OMEGA_1_1_DESIGN)

Autor: cerebro orquestador (Opus). Fuentes: `docs/OMEGA_1_1_PANEL.md` (P§n), `docs/OMEGA_MASTER.md` (M§n), `docs/ANALYSIS.md` (A§n) y el código real.

**Estado verificado del repo:** HEAD 6a8ea1b. Entre 9dfbea7 y HEAD el código es idéntico (`git diff 9dfbea7 HEAD -- omega tests pyproject.toml` sale vacío). Pasan 233 tests rápidos. `mypy --strict omega` está limpio; `mypy --strict tests` ya tenía 2 errores heredados, que no se tocan.

**Baseline:** commit `9dfbea7d6be28a8f17fbf336579c23628ac5669d`. El tag local `OMEGA_1_0_BASELINE` es un objeto anotado (be4b729) que apunta a ese commit.

**Prototipos:** se ejecutaron 22 scripts en el scratchpad. Todas las cifras de este documento salen de ellos.

---

## 0. Resumen ejecutivo y hallazgos de prototipos

1. **Los tres umbrales homogéneos y γ están verificados.**
   - Con N=200 y β=1: α = 0.010101 (α̂=1), 0.015152 (α̂=3/2) y 0.020202 (α̂=2).
   - H(1)−H(0) = βC(N,2)(1−2α̂/3). w*=1/α̂ es estacionario con |G| ≤ 5e-15.
   - **γ no altera la trayectoria homogénea:** con un dt común estable, la diferencia máxima en instantáneas es ≤1.4e-14 para γ̂ ∈ {0,1,100}.
   - **Advertencia:** si se usa el dt automático de γ=0 con γ̂=100, la trayectoria diverge. Es inestabilidad numérica (dt·2γ(N−2) > 2), no física. Los tests deben usar un dt estable para el γ máximo.
   - Desde U(0,1), J−I es inaccesible para α̂∈(1,2) aunque sea localmente estable (y globalmente preferido si α̂>3/2).
2. **RESISTANCE no puede dar D por crecimiento de bolas en d≥2. Es matemática, no un bug:** la resistencia efectiva está acotada en grafos transitorios.
   - Medido: anillo D_res=0.99; toros 2D y 3D dan `no_window`, con rango R∈[0.33, 0.48] en 3D.
   - Decisión: D_res se informa, pero no vota. La resistencia entra mediante el exponente ζ_R=ln(R̄(r_hi)/R̄(1))/ln r_hi.
   - Valores de ζ_R: 1D≈0.92–0.97; 2D 0.37–0.48; 3D 0.21–0.30; Watts–Strogatz p=0.05: 0.51.
3. **Hallazgos red-team con D≈3 espurio:**
   - Watts–Strogatz (N=800, k=6, p=0.05) da D_hop=3.17 ("ok", sin plateau). Con p=0.2 da 3.92.
   - Un árbol 3-ario de altura 5 da D_hop=2.83 *con plateau*.
   - Los atrapan: D_s/D_Weyl (WS p=0.05: 1.85/1.55), la localidad, la cola de curvatura negativa y la conectividad de anillos.
4. **D_Weyl** usa el Laplaciano combinatorio sobre el grafo binario, la escalera "upper" (cuenta de autovalores no nulos ≤ nivel) y una ventana de cuenta [10, 0.2N].

   | Grafo | anillo 300/1000 | 17² | 28² | 32² | 7³ | 9³ | 10³ | 15³ | RGG3 (N=800, k=12) | RGG2 (N=800, k=10) |
   |---|---|---|---|---|---|---|---|---|---|---|
   | D_Weyl | 1.01 | 2.13 | 2.08 | 2.07 | 2.88 | 2.98 | 3.08 | 3.15 | 2.96 | 1.83–1.90 |

   - El normalizado es peor en grafos irregulares (RGG3: 2.61). La pendiente local no sirve como plateau (std 1–8 por degeneración), así que la calidad se mide con R² ≥ 0.95 (retículas 0.959–0.994; RGG 0.998).
   - Nulos densos U(0,1) dan D_W de 13 a 121 (nunca ≈3). K_n da `no_window`.
5. **Suite de distancias.** El estimador "shell" de Ω-1.0 (malla log) es bueno para pesos aleatorios en toros: 15³ con W∈U(0.3,1) da INV 2.89, LOG 2.83, HOP 2.92.
   - Falla LOG en un RGG con pesos euclidianos porque quedan cáscaras vacías (perfil escalonado). Un reintento preregistrado con n_radii=12 lo resuelve (2.49).
   - Una malla lineal adaptativa se probó y se rechazó: sesgo sistemático de +0.6 (15³: 3.49).
   - **En estados binarios (W∈{0,1}, que es el caso de clip en Ω-1.0), HOP=INV=LOG exactamente:** la suite solo discrimina en estados no binarios (Langevin, sigmoide, Ω-B, líneas base).
6. **Curvatura de Ollivier** (idleness 0.5, métrica de saltos, W1 por LP HiGHS, unos 2.5 ms por arista):
   - Toros hipercúbicos 2D/3D: exactamente 0.
   - S² frente a T² (RGG, N=300, k=20, 3 semillas): 0.104 vs 0.092 (≈4σ). Con N=400 la diferencia es 0.005 (≈2σ).
   - Colas: RGG3 k12 media −0.015, frac(κ<−0.5)=0.005; WS p=0.05: 0.043; ER k12 media −0.43; RGG recableado −0.43.
   - **QRC (Klitgaard–Loll) no es usable con δ≤3** (artefactos de retícula: en T², d̄/δ=1.88/1.68/1.65). Se implementa solo como informe.
   - **Forman se descarta:** en T³ da 4−d_u−d_v=−8, no 0, sin 2-celdas.
7. **Topología:**
   - El complejo de cliques de una retícula cúbica no tiene triángulos: en T³ 9³ da β1=1459.
   - **β1 con relleno de ciclos de longitud ≤4 (triángulos+cuadrados, sobre Z/2)** da exactamente β1(T³)=3, β1(T²)=2 y anillo 1.
   - La densidad β1/E separa RGG3 k12 (0.016) de WS p=0.05 (0.052) y ER k12 (0.33). Coste <0.3 s con N=800.
   - El complejo de cliques de un RGG a N≤800 es ruidoso (T³: (1,25,13) con k=20), así que no se exigen valores concretos de β.
8. **Coarse-graining** (emparejamiento de arista pesada + agregación max):
   - Conserva D en 1D (1.00 en 4 niveles) y en 2D (2.00→2.10→2.17).
   - En 3D con N≤1000 la ventana desaparece tras un nivel. Es una limitación de tamaño finito.
9. **Ω-B (densidad fija):**
   - El estado uniforme W=ρ es estacionario (multiplicador de Lagrange).
   - Es linealmente inestable si y solo si **α̂ρ(N−4)/(N−2) > 1+γ̂**. Deducido de los autovalores del grafo de líneas de K_N: N−4 en modos de grado y −2 en modos de ciclo.
   - Verificado con N=80 y ρ=0.1: α̂=8 queda uniforme; α̂=15 da un único clique de 26 nodos (C(26,2)=325≈ρM) más nodos aislados; con γ̂=10 sigue uniforme hasta α̂=40.
10. **Tamaño finito:**
    - El control positivo RGG3 k12 no tiene ventana 3D con N∈{300, 500}. Sí la tiene con N=800 (3.06), 1500 (3.11) y 3000 (3.16).
    - **Con N≤800 es imposible certificar `size_robust` en 3D** (riesgo R1).

---

## 1. Decisiones científicas (valor por defecto y justificación)

Todos los valores se congelan en `omega/config/settings11.py` y entran en `config_hash` antes de ejecutar el paso 6 (Exp 3).

### 1.1 DistanceMode (P§4)

| Modo | Longitud de arista sobre A=(W>w_min) | Notas |
|---|---|---|
| HOP | 1 | Es `hop_distance_matrix` de Ω-1.0 (D_eff_hops). |
| WEIGHTED_INVERSE | 1/(W+ε), con ε=1e-9 | **Llama a `distance_matrix` de Ω-1.0.** Su D es bit a bit el `d_eff` de Ω-1.0 cuando este es "ok". |
| WEIGHTED_LOG | ℓ = 1 − ln(max(W, 1e-12)) | W=1 da ℓ=1 (límite de saltos, sin colapso); W→0 da ∞; monótona. Se rechaza ℓ=−ln W porque daría ℓ=0 con W=1 y colapsaría las retículas de validación a un punto. |
| RESISTANCE | R_ij = L⁺_ii + L⁺_jj − 2L⁺_ij, con L el Laplaciano de W⊙A sobre la componente gigante | pinv hermítica, coste O(n³) (0.5 s con n=800); se rechaza si n>2000. |

- **Ámbito:** la suite opera sobre la componente gigante de `threshold_adjacency(w, w_min)` (umbral estricto >), igual que `geometry_observables`.
- **Estimador común:** `effective_dimension` de Ω-1.0 (shell, malla log, ventana D-6). Si un modo no entero (INV, LOG, RES) da `no_window`, se reintenta una sola vez con `n_radii=12` y se marca `fallback_used`.
- **`distance_sensitive`** es True salvo que HOP, INV y LOG estén todos "ok" y su dispersión máxima por pares sea ≤ `metric_tol`=0.5.
  - Calibración: binarios 0; 9³ pesado 0.21; 15³ pesado 0.09; 28² pesado 0.10.
  - El RGG3 con pesos euclidianos a N=800 da 0.6 y es distance-sensitive. Se documenta como sensibilidad real de esa ponderación a ese N.
- **`resistance_consistent(ζ, clase)`:**
  - clase 1: ζ ≥ 0.8;
  - clase 2: 0.35 ≤ ζ ≤ 0.65;
  - clase ≥3: ζ < 0.35;
  - clase None: False.
- **`metric_robust`** = ¬distance_sensitive ∧ resistance_consistent.

### 1.2 Filtración θ, Betti, χ y persistencia (P§5, P§12)

- G_θ = (V, {W > θ}), con umbral estricto como en las distancias.
- Malla por defecto: θ ∈ (0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01), descendente.
- Por cada θ se miden β0 (componentes, incluidos aislados), fracción gigante, k̄, aristas, β1^(4) y densidad b1.
- **β1^(4):** dimensión del espacio de ciclos sobre Z/2 módulo el subespacio generado por triángulos y 4-ciclos:
  - β1^(4) = E − rank∂1 − rank(caras);
  - presupuesto: aborta con `over_budget` si #triángulos + #C4 > 500 000;
  - fórmulas cerradas para contar: T = tr(A³)/6 y C4 = (tr(A⁴) − 2Σk_i² + 2m)/8 (comprobado: K4 da 3).
- **Complejo de cliques** (Z/2, truncado en tetraedros, de modo que β0, β1 y β2 son exactos para el complejo):
  - presupuesto: aborta con `over_budget` si #triángulos + #tetraedros > 300 000;
  - se informan χ_β = β0 − β1 + β2 y χ_cuentas = V − E + T − K4;
  - solo informe: no se exigen valores concretos (ruidoso a N≤800 y sesgado por retículas sin triángulos, P§3).
- **Rango sobre GF(2):** reducción de columnas con máscaras de bits en enteros de Python. Hay que convertir a `int` antes de desplazar bits (`np.int64` desborda).
- **Persistencia simple:**
  - H0 exacta por Kruskal sobre aristas W>0 ordenadas por (−W, i, j). Barras (nacimiento=1.0, muerte=W_e); las esenciales son las componentes con W>0.
  - Curvas de Betti β0(θ) y β1^(4)(θ).
  - No hay código de barras H1 (coste y complejidad no justificados a este N).

### 1.3 D_Weyl (P§6)

- Grafo: componente gigante de A estricta, binaria. Laplaciano combinatorio L = D − A.
- Se descarta exactamente un autovalor nulo. Los niveles se agrupan con tolerancia relativa 1e-9 y N(λ_u) = #{λ_j>0 : λ_j ≤ λ_u}.
- Ventana: niveles con cuenta ∈ [10, 0.2·n_gc] y al menos 4 niveles; si no, `no_window`. Con n<10 el estado es `insufficient_component`.
- D_W = 2·pendiente MCO de ln N frente a ln λ; stderr = 2·se.
- `plateau := R² ≥ 0.95`. La pendiente local no se usa como plateau (std 1–8 por degeneración, medido).
- Validación (Ω-0), valor medido ± tolerancia:

  | Grafo | anillo 300 | anillo 1000 | 17² | 28² | 32² | 7³ | 9³ | 10³ | 15³ (slow) | RGG3 N=800 k12, 3 semillas | RGG2 N=800 k10 |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | D_Weyl | 1.01 ±0.05 | 1.01 ±0.05 | 2.13 ±0.2 | 2.08 ±0.2 | 2.07 ±0.2 | 2.88 ±0.3 | 2.98 ±0.25 | 3.08 ±0.25 | 3.15 ±0.25 | 2.96 ±0.25 | 1.83–1.90 (2±0.25) |

  Además: K_n da `no_window`; U(0,1) con N=300 y w_min ∈ {0.1, 0.5, 0.9} da D_W > 5.
- También se exporta `fiedler_length = λ2^{-1/2}`, la longitud de correlación ξ para tamaño finito.

### 1.4 Homogeneidad e isotropía (P§11)

- **Homogeneidad:**
  - D_i por nodo = 1 + pendiente de ln S_i(r) frente a ln r, sobre los radios enteros de la ventana ajustada de HOP (la misma que da la suite).
  - Si algún S_i(r)=0, D_i es NaN.
  - `homogeneity_ok := σ_D/μ_D ≤ 0.2` con ≥ 10 nodos válidos. Medido: toros 0, RGG3 0.14.
  - A nivel de punto se exige además que la cv no crezca con N (ver `size_robust`).
- **Isotropía:** proxy MDS local.
  - Por cada fuente i (64 fuentes con rng) se toma la bola B_i(R) en saltos, con R = r_hi(HOP) + 1. Se aplica MDS clásico (G = −½JD²J) y se calcula ratio = λ_k/λ_1, con k = clase de dimensión.
  - `isotropy_ok := mediana ≥ 0.5 ∧ percentil 10 ≥ 0.3`. Con clase 1 es trivialmente True.
  - Medido (k=3): toro 1.00; RGG3 k12 0.70/0.58; RGG3 k8 0.63/0.53; slab 32×5×5 0.27 (falla); slab 20×8×5 0.57; WS 0.60–0.77.
  - Es un diagnóstico sobre una métrica ya medida; no alimenta la dinámica.
  - **Rechazado:** el proxy de "conos geodésicos" daba RGG3 0.60, peor que el slab anisótropo (0.15).

### 1.5 Localidad, small-world y manifold proxy

- **Localidad:** fracción de aristas (gigante, A estricta) con desvío ≤ 3 saltos en G−e. Fórmula exacta:
  - (A²)_uv > 0 ∨ (A³)_uv − k_u − k_v + 1 > 0.
  - `locality := fracción ≥ 0.97`. Medido: toros 1.0; RGG3 k12 0.998; RGG3 k8 0.979; WS p=0.05 0.949; WS p=0.2 0.82; ER k12 0.86.
- **Small-world** (bandera que forma parte de locality): small_world := (C_bin/ρ ≥ 3) ∧ (frac(κ<−0.5) > 0.03), donde C_bin y ρ vienen de los TopologyObservables de Ω-1.0. El campo final es `locality_field := locality ∧ ¬small_world`.
- **Conectividad de anillos:**
  - Para r ∈ [2, r_hi(HOP)], A_i(r) = {j : r ≤ d_ij ≤ r+1}; subgrafo inducido conexo; 64 fuentes.
  - `annulus_ok := min_r fracción ≥ 0.9`. Si r_hi < 2, es False ("escala insuficiente").
  - Medido: 9³ 1.0; RGG3 k12 0.97 (r=2); RGG3 k8 0.25; WS p=0.05 0.00; árbol 0.00; anillo 0.00.
- **Curvatura OK:** `curvature_ok := κ̄ ≥ −0.1 ∧ frac(κ<−0.5) ≤ 0.03`.
- **Manifold proxy:** `manifold_proxy_ok := annulus_ok ∧ curvature_ok ∧ b1_density(w_min) ≤ 0.03`. En el código se llama `MANIFOLD_PROXY`, nunca "manifold".

### 1.6 Curvatura (P§12–13)

- **Ollivier-Ricci** con medida m_x = 0.5·δ_x + 0.5·W_xy/Σ_y W_xy sobre los vecinos en A.
  - Coste de suelo: distancia del modo `CurvatureConfig.distance_mode` (HOP por defecto) sobre la componente gigante.
  - W1 por `scipy.optimize.linprog(method="highs")`, problema de transporte. Hasta 1000 aristas muestreadas con rng, sin reemplazo.
  - Salida: media, desviación, SE, fracción de cola y curvatura escalar por nodo.
  - Validación: toros hipercúbicos κ=0 en cada arista (≤1e-9); ciclo C_n (n≥6) κ=0; K_n κ = 1 − |α − (1−α)/(n−1)| (≤1e-9); S² > T² (slow).
- **QRC:** d̄(S_δ(x), S_δ(y))/δ con δ ∈ {1,2,3} y hasta 200 pares. Solo informe; no entra en el certificado.

### 1.7 Coarse-graining preregistrado (P§14)

La regla `heavy_edge_matching/max/v1` se congela antes de Ω-11.

1. Aristas con W > w_min ordenadas por (−W, u), con u~U(0,1) del rng (desempate aleatorio: los estados binarios tienen empates masivos y un desempate por índice violaría M§32).
2. Emparejamiento greedy maximal; los nodos no emparejados quedan como bloques unitarios.
3. W'_AB = max_{i∈A, j∈B} W_ij. En binarios equivale a la contracción (grafo cociente), no diluye frente al umbral y se queda en [0,1]. "mean" es la alternativa informada.
4. Bloques canónicos: tuplas ordenadas, ordenadas por su mínimo.
5. Niveles hasta n' < 50 o 4 niveles.
6. 3 réplicas de rng. La métrica es la consistencia de clase (P§14): todos los niveles medibles con la misma `dimension_class` y sin nuevo código de fallo. Los niveles con `no_window` cuentan como "indeterminado", no como fallo.

### 1.8 Motores estadísticos (P§8, P§16)

- **Temperatura:** Θ = Θ̂·β (escala por arista del término βw²). Malla Θ̂ ∈ {0.01, 0.03, 0.1, 0.3, 1}.
- **Langevin:**
  - Sobre el vector triangular u: u ← B(u − dt·g + √(2Θdt)·ξ), ξ ~ N(0,1).
  - dt = 0.1/L(n,p) (L de Ω-1.0; exige η=0).
  - **Frontera por defecto: reflexión** B(x) = y si y≤1, 2−y si no, con y = x mod 2. Así la densidad estacionaria es ∝ exp(−H/Θ) en [0,1]^M.
  - "project" (clip) queda solo como alternativa documentada: crea átomos en 0 y 1 y no es Gibbs. Hay un test que lo demuestra.
- **Metropolis:**
  - Una propuesta por arista en orden rng.permutation(M): w' = B(w + s·U(−1,1)) (simétrica con reflexión).
  - ΔH exacto en O(N) (fórmulas en §3, WP-D).
  - El paso s se adapta solo durante el burn-in (cada 10 sweeps: ×1.1 si aceptación > 0.5, /1.1 si < 0.3, acotado a [1e-4, 1]). Después queda fijo, para no romper el balance detallado.
- **Ensemble:**
  - 4 cadenas desde U(0,1) independientes; burn-in = 50%.
  - τ_int de Sokal con ventana automática c=5 (τ = ½ + Σρ_t; se elige el menor M ≥ c·τ(M)). ESS = n/(2τ).
  - R̂ dividido (split-R̂) y Geweke z (primer 10% frente a último 50%).
  - `equilibrated := R̂ ≤ 1.05 (acción y peso medio) ∧ ESS_total ≥ 100 ∧ |z| ≤ 3` en cada cadena.
  - Probabilidad de fase: `assess_run` sobre el estado final más hasta 5 estados post-burn-in separados por ≥2τ_int, con IC de Wilson.
- **Transición** (λ = α̂ o Θ̂):
  - m = peso medio; χ_m = M·Var(m); Binder U4 = 1 − ⟨δm⁴⟩/(3⟨δm²⟩²).
  - Coeficiente de bimodalidad de Sarle (b > 5/9 indica bimodal); d⟨m⟩/dλ por diferencias centrales; λ* = argmax χ.
  - Histéresis con `hysteresis_sweep` de Ω-1.0.
  - En corridas de gradiente, la varianza entre semillas se etiqueta `seed_variance` (no es térmica).

### 1.9 Rama Ω-B: densidad fija (P§9)

- Proyección exacta sobre C = {x ∈ [0,1]^M : Σx = W0}: x = clip(v − ν, 0, 1).
  - ν se obtiene ordenando los 2M puntos de quiebre {v_i, v_i − 1}, bisección sobre el índice e interpolación lineal en el segmento.
  - Precisión: |Σx − W0| ≤ 1e-9·max(1, W0). Si W0 ∉ [0, M] se lanza ValueError.
- W0 = ρ·N(N−1)/2. Si ρ es None se usa ρ = peso medio de W0.
- **El gradiente se extrae con `g[np.triu_indices(n,1)]`, no con `upper_triangle`**, que valida [0,1] y rechaza gradientes.
- Solo integrador clip. Predicción analítica: inestabilidad si y solo si α̂ρ(N−4)/(N−2) > 1+γ̂.

### 1.10 Controles nulos (P§10)

- **ER:** G(n,m) con el mismo m que A del candidato; pesos binarios o barajados.
- **Configuration model:** networkx con semilla del rng y "erased" (sin lazos ni multiaristas).
- **Recableado con grados preservados:** envuelve `degree_preserving_null` de Ω-1.0 con n_swaps = 10·E.
- **Watts–Strogatz:** implementación propia con rng; k = par más cercano a k̄_bin; p ∈ {0.01, 0.05, 0.1, 0.2, 1}.
- **RGG (control positivo):** en el toro T^d (binario o euclidiano) y en la esfera S^d. El ángulo de casquete se obtiene por bisección (1e-12): S² usa (1−cos r)/2 = k̄/(n−1); S³ usa (r − sin r cos r)/π.
- **Árbol balanceado:** nulo red-team (D_hop 2.83).
- Batería por defecto contra cada candidato: {DEGREE_PRESERVING_REWIRE, CONFIGURATION_MODEL, ERDOS_RENYI, SHUFFLED_WEIGHTS}. WS y RGG se usan en Ω-3.

### 1.11 Taxonomía Ω-F0..Ω-F10 y precedencia (P§24)

**Banderas por corrida** (todas se evalúan; evidencia ausente o no finita cuenta como False):

| Bandera | Regla | Código si falla |
|---|---|---|
| converged | status == CONVERGED | F10 |
| empty | giant_fraction < 0.1 ∧ mean_binary_degree < 1 (regla A de Ω-1.0) | F0 si True |
| dense_or_uniform | binary_density ≥ 0.5 ∨ mean_weight ≥ 0.5 (regla E) ∨ (cv de pesos ≤ 1e-6 ∧ media > 0) | F1 si True |
| connected | giant_fraction ≥ 0.95 | F2 |
| locality (con small-world incluido) | §1.5 | F3 |
| d_volume / d_spectral / d_weyl | §1.12 | F4 |
| metric_robust | §1.1 | F5 |
| isotropy_ok / homogeneity_ok / topology_stable / manifold_proxy_ok | §1.4, §1.5, §1.12 | F9 |

- nontrivial = ¬empty ∧ ¬dense_or_uniform. El estado uniforme es F1 porque el grafo relacional es K_N uniforme; cubre el uniforme de Ω-B.
- **Por punto:**
  - F10 si la fracción CONVERGED < 0.8 o el ensemble no está equilibrado.
  - F7 si ninguna salida de corrida (código primario o PASS) alcanza 0.8; si el modal es un fallo, ese código.
  - F6: `size_robust` falso. F7: `seed_robust` falso. F8: `null_separated` falso.
- **Precedencia del código primario:** F10 > F0 > F1 > F2 > F3 > F4 > F5 > F7 > F6 > F8 > F9.
  - Se informa siempre el conjunto completo de códigos.
  - Racional: primero las trivialidades; la tendencia en tamaño solo se evalúa sobre una señal estable entre semillas; la separación de nulos, sobre una señal robusta; la manifoldicidad (lo más exigente) al final.
- Ω-CANDIDATE solo con conjunto vacío de códigos y certificado completo.

### 1.12 GeometryCertificate (P§7): campos y evidencia

D* = mediana de {D_vol = d_eff de Ω-1.0 (INV shell), D_s de Ω-1.0, D_Weyl}, solo si los tres están "ok"; si no, D* = NaN. `dimension_class` = round(D*) si |D* − round| ≤ 0.35; si no, None.

**Campos de nivel de corrida.** A nivel de punto, cada uno es True si es True en ≥ 80% de las semillas.

| Campo | Regla |
|---|---|
| connected | giant_fraction ≥ 0.95 |
| nontrivial | ¬empty ∧ ¬dense_or_uniform |
| locality | detour3 ≥ 0.97 ∧ ¬small_world |
| d_volume | D_vol ok ∧ plateau ∧ \|D_vol − D*\| ≤ 0.35 |
| d_spectral | D_s ok ∧ plateau ∧ \|D_s − D*\| ≤ 0.35 |
| d_weyl | D_W ok ∧ R² ≥ 0.95 ∧ \|D_W − D*\| ≤ 0.35 |
| metric_robust | §1.1 |
| isotropy_ok | §1.4 |
| homogeneity_ok | §1.4 |
| topology_stable | para θ ∈ `graph.w_min_sensitivity`, s(θ) = (giant_fraction(θ) ≥ 0.95 ∧ b1_density(θ) ≤ 0.03 ∧ status ok); True si s(θ) se cumple en ≥ 80% de los θ |
| manifold_proxy_ok | §1.5 |

Calibración de d_*: RGG3 da 3.06/2.60/2.91 con D*=2.91 y pasa; WS p=0.05 da 3.17/1.85/1.55 con D*=1.85 y falla d_volume.

**Campos solo de nivel de punto:**
- `seed_robust`: ≥ 80% de las corridas sin código ∧ std(D*) ≤ 0.2 ∧ ≥ 80% con la misma `dimension_class`.
- `size_robust`:
  - ≥ 3 tamaños con ≥ 80% de corridas sin código;
  - rango de la media de D* ≤ 0.3;
  - |dD*/d ln N| (MCO) ≤ 0.1. Calibración RGG3: 3.06/3.11/3.16, pendiente 0.075. **Sustituye a la regla de deriva monótona de Ω-1.0**, que rechazaría al control positivo;
  - cv de homogeneidad en el N mayor ≤ cv en el N menor + 0.02.
- `null_separated`: para cada nulo de la batería, fracción de réplicas nulas sin código ≤ 20% ∧ (si D*_null es finito) |media D*_cand − media D*_null| ≥ 3σ_null, o ≥ 0.35 si σ = 0.

**Veredicto:** GEOMETRIC_CANDIDATE ("Ω-CANDIDATE") si y solo si los 14 campos son True ∧ `dimension_class` no es None ∧ no hay F10. Si falta evidencia (tamaños o nulos sin ejecutar), el campo es False y el veredicto es NOT_CANDIDATE con `missing` en la evidencia.

**Regla de código:** la constante `RULE_D3_NEVER_SUFFICIENT` y un test exhaustivo sobre las 2^14 combinaciones. Ningún identificador d3/D3/is_3d/three_d en `omega/`. P1 (clase 3) se informa como `dimension_class == 3`, nunca como un booleano.

### 1.13 Pasaporte ampliado (P§17)

Módulo nuevo `omega/io/provenance.py`. El pasaporte v1.0 queda intacto y se embebe en el v1.1.

| Campo | Definición |
|---|---|
| model_version | "Ω-1.1" |
| baseline | {name: "OMEGA_1_0_BASELINE", commit: 9dfbea7d…} |
| code_commit, git_branch | `git rev-parse HEAD` y `--abbrev-ref HEAD` por subprocess; "unknown" si falla |
| python_version | versión de Python |
| dependency_lock_hash | sha256 de "\n".join(sorted(f"{n}=={v}")) para python, numpy, scipy y networkx; se guarda también el dict |
| config_hash | sha256(canonical_json(config11_to_dict(cfg))), con canonical_json = json.dumps(sort_keys=True, separators=(",",":"), ensure_ascii=False, allow_nan=False); incluye el OmegaConfig base |
| initialization_distribution | p. ej. "uniform/upper_mirror" |
| distance_definition | modos y fórmulas |
| threshold_rule | literal: "distances/curvature/filtration: W>w_min (estricto); topology_observables Ω-1.0: W>=w_min" |
| dimension_algorithm | "d_vol: shell log-grid Ω-1.0; d_s: lazy q=0.5; d_weyl: combinatorial binary staircase [10,0.2N] R2>=0.95" |
| coarse_graining_algorithm | "heavy_edge_matching/max/v1" o null |
| null_model | modelo nulo o null |
| random_seed | {entropy, spawn_key, "PCG64"} |
| termination_reason | ∈ {CONVERGED, MAX_STEPS, NONFINITE, STATIC, CHAIN_COMPLETE, CHAIN_NOT_EQUILIBRATED} |
| experiment_id | "Ω-EXP-%06d" (ver abajo) |
| reconstruct | {entrypoint "módulo:función", config, seed, w0_sha256} |
| omega10 | `build_passport` de Ω-1.0 si es una corrida de gradiente |

- **experiment_id:** la etiqueta JSON es "Ω-EXP-%06d"; el nombre de archivo, "OMEGA-EXP-%06d" (el `_RUN_ID` de Ω-1.0 es solo ASCII).
- La asignación es atómica: `os.open(registry/ids/NNNNNN.claim, O_CREAT|O_EXCL)` con siguiente = max + 1 y reintento ante colisión.
- **Reconstrucción:** `reconstruct(passport) → (cfg, SeedKey)`; volver a ejecutar el entrypoint reproduce los digests de los arreglos.

### 1.14 Tamaño finito, susceptibilidad y transición (P§15–16)

- Tamaños: N ∈ (64, 100, 150, 200, 300, 500, 800). Tamaños XL (1500, 3000) solo `@slow`, y solo para puntos que pasan a nivel de corrida con N=800 (ver R1).
- Observables: D*, cada D, C_bin, G/N, ρ binaria, ξ_F = λ2^{-1/2} y L_hop.
  - Ajustes: pendiente log-log y ajuste en ln N.
  - ξ_F ~ N^{1/D} en estructuras geométricas y aproximadamente constante en expansores; es el discriminante F3/F6.
- Susceptibilidad y transición: §1.8. Para transiciones en ensemble se informan χ_max(N) y λ*(N).

---

## 2. Tests analíticos: `tests/analytical/` (P§2, paso 1 de P§22)

**Tolerancias:**
- Identidades algebraicas: |a−b| ≤ 1e-12·max(1, |a|, |b|).
- Trayectorias: ≤ 1e-12 absoluto.
- Estructuras exactas: igualdad bit a bit cuando se indica.
- Por defecto, β=1 y parámetros vía `reduced_to_raw`.

### `test_empty_state.py`
- **E1.** `grad_action(0, p)` es exactamente 0.0 en todas las entradas, para N ∈ {5, 12, 40}, (α̂, γ̂) ∈ {0, 1, 3}², η ∈ {0, 0.2} y μ=0.
- **E2.** `action(0, p) == 0.0`.
- **E3.** Hessiano en 0 por diferencias centrales del gradiente (N=8, M=28, h=1e-6). Coincide con 2βI + 2γBᵀ(I − 11ᵀ/N)B (rel 1e-6) y λ_min = 2β (rel 1e-6). T no contribuye en 0.
- **E4.** Con α=0: S0(W) ≥ 0 en 200 matrices W aleatorias (N=12, U(0,1) y versiones dispersas), con γ ∈ {0, 1}, y S0 > 0 si W ≠ 0.
- **E5.** W0 = 0.05·U(0,1), N=40, α̂ ∈ {0.5, 1, 1.5, 3}, γ̂ ∈ {0, 1}: CONVERGED, max w_final ≤ 1e-6 y etiqueta A de Ω-1.0. W=0 es localmente estable para todo α̂.
- **E6.** W0 = (0.9/α̂)(J−I) con α̂ ∈ {1.2, 2, 3}: converge a 0.

### `test_complete_state.py`
- **C1.** T(J−I) = C(N,3), S_dens = C(N,2) y S_deg = 0, para N ∈ {5, 40, 200}.
- **C2.** `grad_action(J−I)` = 2β(1−α̂)(J−I) entrada a entrada, para γ̂ ∈ {0, 1, 100} (independiente de γ).
- **C3 (KKT, umbral α̂=1).** `projected_step(J−I, G, dt)` == J−I bit a bit para α̂ ∈ {1.0, 1.1, 3}. Para α̂=0.9 todas las aristas bajan exactamente dt·2β·0.1.
- **C4 (estabilidad local).** N=40:
  - W0 = 0.99(J−I): α̂=1.05 lleva exactamente a J−I con CONVERGED; α̂=0.95 lleva a 0 (≤ 1e-6).
  - W0 = (J−I) − 0.01·U(0,1) (simetrizado, no homogéneo) con α̂=1.2 lleva a J−I.
- **C5 (umbral α̂=3/2).** action(J−I) − action(0) = βC(N,2)(1 − 2α̂/3) para α̂ ∈ {1.4, 1.5, 1.6}. Los signos son +, 0 (≤ 1e-9·C(N,2)) y −.

### `test_uniform_manifold.py`
- **U1.** action(w(J−I)) = −αC(N,3)w³ + βC(N,2)w² para γ ∈ {0, 1, 100} y w ∈ {0.1, 0.37, 0.8}.
- **U2.** grad = 2βw(1 − α̂w)(J−I), y max|grad_degree_irregularity| ≤ 1e-12.
- **U3.** w* = 1/α̂ es estacionario para α̂ ∈ {1.25, 2, 3}: max|G| ≤ 1e-12·2β. Para α̂ ≤ 1 no hay estacionario en (0,1).
- **U4.** w* es un MÁXIMO de H: H(w* ± 0.01) < H(w*) y H''(w*) = −2βC(N,2) (diferencias finitas, rel 1e-6).
- **U5 (tres umbrales distintos).**
  - α(α̂=1, 1.5, 2; N=200) redondeado a 4 decimales = 0.0101, 0.0152, 0.0202; α(α̂=1; N=100) = 0.0204.
  - α̂=1: estabilidad local de J−I (cubierta por C3).
  - α̂=3/2: H(1) = H(0).
  - α̂=2 con w0 = ½ homogéneo (N=200): G ≡ 0 (≤ 1e-12); α̂=1.98 lleva a 0; α̂=2.02 lleva a J−I.
  - Desde U(0,1) (N=100 rápido; N=200 `@slow`; semillas `seed_key(SeedConfig(20261004, 950, 3), 0, r)`): α̂ ∈ {1.2, 1.6, 1.9} da A en 3/3 y α̂=2.1 da E en 3/3. Se afirma explícitamente que J−I es localmente estable (y globalmente preferido con α̂ > 3/2) pero inaccesible desde U(0,1) para α̂ < 2.
- **U6 (γ no cambia la trayectoria homogénea).**
  - N=50; pares (α̂, w0) ∈ {(0.5, 0.9), (1.2, 0.8), (1.2, 0.9), (2, 0.3), (2, 0.49), (2, 0.51)}.
  - `DynamicsConfig(dt_mode="fixed", dt_fixed=0.5/L(γ̂=100), max_steps=20000)`, con γ̂ ∈ {0, 1, 100}.
  - Mismos status y número de pasos; instantáneas ≤ 1e-12 (medido ≤ 1.4e-14).
  - Coincide con la recursión escalar u_{t+1} = clip(u_t − dt·2βu_t(1 − α̂u_t), 0, 1), ≤ 1e-12 en cada paso de instantánea.
- **U7.** L(γ̂) = 2β(1 + γ̂ + 2α̂) exactamente. El dt automático depende de γ̂, y esa es la única vía por la que γ entra en la trayectoria homogénea.

### `test_parameter_scaling.py`
- **S1.** `lipschitz_bound(n, reduced_to_raw(α̂, γ̂, n))` = 2β(1 + γ̂ + 2|α̂|) para N ∈ {10, 64, 200, 800}. No depende de N.
- **S2.** Las trayectorias homogéneas son invariantes en N en parámetros reducidos:
  - N ∈ {20, 50, 120}, α̂ = 2.5, γ̂ ∈ {0, 3}, w0 ∈ {0.35, 0.45}, dt automático;
  - mismos pasos y status; escalares y w_final iguales ≤ 1e-12.
- **S3.** Al escalar (α, β, γ) por λ ∈ {0.5, 2, 10}: la acción y el gradiente se multiplican por λ, y `evolve` con dt automático da trayectorias idénticas (≤ 1e-12; N=30, misma W0 U(0,1), max_steps=300).
- **S4.** α_c(N) = 2βα̂_c/(N−2) para α̂_c ∈ {1, 1.5, 2} y N ∈ {64, 100, 200, 800}.
- **S5.** En la malla cruda de M§25 (α/β ≥ 0.1), α̂ = 0.05(N−2) ≥ 2 para N ≥ 42: se afirma α̂ > 2 en todos los puntos con α > 0 y N ∈ {100, 200}. Explica el colapso E de la malla cruda de Ω-1.0.

**Golden baseline** (P§2 y paso 2; lo hace WP-0, antes que nada).
- `tools/freeze_baseline.py` se ejecuta sobre el árbol actual (= 9dfbea7).
- Escribe `tests/expected_baseline/` con floats en repr de 17 dígitos y `"generated_from_commit"`.

| Fixture | Contenido |
|---|---|
| empty_phase.json | `default_config(40, master_entropy=20261004, experiment_id=900, replicates=3)` con `reduced_to_raw(1.0, 1.0, 40)` y `simulate` para rep 0..2. Se guardan: status, steps, dt, tau, label, flags, frac_at_zero/one, mean_strength, `scalars["action"]` en índices [0, 1, 2, −1], `scalars["max_dw"][-1]`, digest de w0 (exacto) y w_final_upper como lista (N=40 da 780 valores) |
| complete_phase.json | igual con α̂=3.0 y γ̂=1.0 |
| uniform_barrier.json | N=40, γ̂=0, `DynamicsConfig()` por defecto. Arranques homogéneos w0 ∈ {w*−0.01, min(w*+0.01, 1)} para α̂ ∈ {1.0, 1.5, 2.0}: steps, status, w final, action en [0, 1, −1]. Más U(0,1) en rep 0 con α̂ ∈ {1.9, 2.1}: etiqueta A/E |
| geometry_reference.json | `scan.observe` sobre `periodic_lattice((7,7,7))` (n=343) y sobre U(0,1) con N=60, SeedKey(20261004, (901, 0, 0)): todos los campos de Topology/GeometryObservables (value, stderr, status, window, plateau) |
| passport_keys.json | conjunto de claves de `build_passport` (nivel superior y `observables`) |
| omega10_manifest.json | sha256 de cada archivo de Ω-1.0 (`omega/**/*.py`, `tests/*.py`, `omega/experiments/test_0*.py`, `pyproject.toml`). Excepciones permitidas: `omega/__init__.py` y `pyproject.toml` (solo versión), `omega/curvature/__init__.py` y `omega/coarse_graining/__init__.py` (solo docstring) |

`tests/test_baseline_golden.py` recalcula cada fixture:
- enteros, cadenas y booleanos exactos;
- floats y arreglos ≤ 1e-12 absoluto;
- digest de w0 exacto (PCG64 es portable);
- manifiesto: sha256 iguales salvo las excepciones;
- duración < 20 s.

---

## 3. Contratos de API

**Reglas globales** (como A§5):
- `mypy --strict`, float64, funciones puras sin mutación ni estado global.
- RNG explícito (`np.random.Generator`); el test de arquitectura existente ya escanea todo `omega/` y `tests/`.
- E/S solo en `omega/io/` y `omega/experiments/`; imports absolutos en `io`.
- ValueError/TypeError en entradas inválidas. Basenames de tests únicos (no hay `__init__` en `tests/`).
- **Ningún archivo de Ω-1.0 se modifica** (lo garantiza el manifiesto). Ω-1.1 solo añade archivos.

### 3.0 WP-0: tipos compartidos (los crea la sesión principal antes de lanzar agentes)

```python
# omega/baseline.py
BASELINE_NAME: Final = "OMEGA_1_0_BASELINE"
BASELINE_COMMIT: Final = "9dfbea7d6be28a8f17fbf336579c23628ac5669d"
BASELINE_SCHEMA_VERSION: Final = "1.0"
MODEL_VERSION: Final = "Ω-1.1"
```

```python
# omega/config/settings11.py  (frozen, slots, __post_init__ valida, mismo estilo que settings.py)
class DistanceMode(Enum): HOP="hop"; WEIGHTED_INVERSE="weighted_inverse"; WEIGHTED_LOG="weighted_log"; RESISTANCE="resistance"
GEODESIC_MODES: Final[tuple[DistanceMode, ...]] = (DistanceMode.HOP, DistanceMode.WEIGHTED_INVERSE, DistanceMode.WEIGHTED_LOG)
class NullModel(Enum): ERDOS_RENYI="erdos_renyi"; CONFIGURATION_MODEL="configuration_model"; DEGREE_PRESERVING_REWIRE="degree_preserving_rewire"; SHUFFLED_WEIGHTS="shuffled_weights"; SMALL_WORLD="small_world"; RANDOM_GEOMETRIC="random_geometric"
class Engine(Enum): GRADIENT="gradient"; LANGEVIN="langevin"; METROPOLIS="metropolis"; FIXED_DENSITY="fixed_density"

DistanceSuiteConfig(modes: tuple[DistanceMode,...] = (HOP, WEIGHTED_INVERSE, WEIGHTED_LOG, RESISTANCE), metric_tol: float = 0.5,
    log_floor: float = 1e-12, fallback_n_radii: int = 12, resistance_max_nodes: int = 2000,
    zeta_1d_min: float = 0.8, zeta_2d_min: float = 0.35, zeta_2d_max: float = 0.65, zeta_3d_max: float = 0.35)
WeylConfig(laplacian: Literal["combinatorial","normalized"] = "combinatorial", graph: Literal["binary","thresholded_weighted"] = "binary",
    count_min: int = 10, count_max_frac: float = 0.2, min_levels: int = 4, degeneracy_rtol: float = 1e-9, r2_min: float = 0.95, min_nodes: int = 10)
TopologyConfig(thetas: tuple[float,...] = (0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05,0.01), clique_max_dim: int = 2,
    max_simplices: int = 300_000, short_cycle_length: int = 4, max_faces: int = 500_000, b1_density_max: float = 0.03, stable_fraction: float = 0.8)
CurvatureConfig(idleness: float = 0.5, distance_mode: DistanceMode = DistanceMode.HOP, max_edges: int = 1000, mean_min: float = -0.1,
    tail_cut: float = -0.5, tail_frac_max: float = 0.03, qrc_deltas: tuple[int,...] = (1,2,3), qrc_max_pairs: int = 200)
LocalStructureConfig(homogeneity_cv_max: float = 0.2, homogeneity_min_nodes: int = 10, isotropy_sources: int = 64, isotropy_radius_offset: int = 1,
    isotropy_median_min: float = 0.5, isotropy_p10_min: float = 0.3, locality_min: float = 0.97, annulus_sources: int = 64,
    annulus_min: float = 0.9, annulus_r_min: int = 2, small_world_clustering_ratio: float = 3.0)
CertificateThresholds(g_connected: float = 0.95, dim_tol: float = 0.35, class_tol: float = 0.35, seed_fraction: float = 0.8, seed_std: float = 0.2,
    size_range_tol: float = 0.3, size_slope_tol: float = 0.1, min_sizes: int = 3, homogeneity_size_slack: float = 0.02,
    null_pass_max: float = 0.2, null_sigma: float = 3.0, uniform_cv_max: float = 1e-6, empty_g: float = 0.1, empty_kbin: float = 1.0,
    dense_rho: float = 0.5, dense_meanw: float = 0.5)
CoarseGrainConfig(rule: Literal["heavy_edge_matching"] = "heavy_edge_matching", aggregation: Literal["max","mean"] = "max",
    n_min: int = 50, max_levels: int = 4, replicates: int = 3)
LangevinConfig(theta_hat: float, n_steps: int = 20_000, dt_safety: float = 0.1, boundary: Literal["reflect","project"] = "reflect",
    thin: int = 10, burn_in_fraction: float = 0.5)
MetropolisConfig(theta_hat: float, n_sweeps: int = 2_000, step_init: float = 0.2, acceptance_low: float = 0.3, acceptance_high: float = 0.5,
    adapt_every: int = 10, thin: int = 1, burn_in_fraction: float = 0.5)
EnsembleConfig(n_chains: int = 4, autocorr_c: float = 5.0, rhat_max: float = 1.05, ess_min: float = 100.0, geweke_z_max: float = 3.0, phase_samples: int = 5)
FixedDensityConfig(rho: float | None = None, projection_tol: float = 1e-9)
FiniteSizeConfig(sizes: tuple[int,...] = (64,100,150,200,300,500,800), xl_sizes: tuple[int,...] = (1500, 3000), replicates: int = 10)
Omega11Config(base: OmegaConfig, distance: DistanceSuiteConfig = ..., weyl: WeylConfig = ..., topology: TopologyConfig = ...,
    curvature: CurvatureConfig = ..., local: LocalStructureConfig = ..., certificate: CertificateThresholds = ..., coarse: CoarseGrainConfig = ...,
    ensemble: EnsembleConfig = ..., finite_size: FiniteSizeConfig = ..., engine: Engine = Engine.GRADIENT,
    langevin: LangevinConfig | None = None, metropolis: MetropolisConfig | None = None, fixed_density: FixedDensityConfig | None = None,
    null_models: tuple[NullModel,...] = (DEGREE_PRESERVING_REWIRE, CONFIGURATION_MODEL, ERDOS_RENYI, SHUFFLED_WEIGHTS),
    schema_version: str = "1.1")
  # valida: base.schema_version == "1.0"; el engine exige su config (LANGEVIN => langevin no None, etc.);
  # Langevin/Metropolis exigen base.functional.eta == 0.
  # Los valores por defecto se escriben con field(default_factory=...).
```

```python
# omega/contracts.py  (dataclasses frozen+slots; sin lógica salvo validación y propiedades triviales)
class FailureCode(Enum):
    F0="Ω-F0"; F1="Ω-F1"; F2="Ω-F2"; F3="Ω-F3"; F4="Ω-F4"; F5="Ω-F5"; F6="Ω-F6"; F7="Ω-F7"; F8="Ω-F8"; F9="Ω-F9"; F10="Ω-F10"
FAILURE_MEANING: Final[Mapping[FailureCode, str]]  # trivial-empty, trivial-complete, fragmented, small-world, dimension-artifact,
                                                   # metric-dependent, finite-size artifact, seed-dependent, null-reproducible, non-manifold, no stable phase
FAILURE_PRECEDENCE: Final = (F10, F0, F1, F2, F3, F4, F5, F7, F6, F8, F9)
class Verdict(Enum): GEOMETRIC_CANDIDATE="Ω-CANDIDATE"; NOT_CANDIDATE="NOT_CANDIDATE"
CERTIFICATE_FIELDS: Final = ("connected","nontrivial","locality","d_volume","d_spectral","d_weyl","metric_robust","seed_robust",
                             "size_robust","null_separated","isotropy_ok","homogeneity_ok","topology_stable","manifold_proxy_ok")
FIELD_CODE: Final[Mapping[str, FailureCode]]  # connected→F2, nontrivial→F1 (F0 lo resuelven las banderas), locality→F3, d_*→F4,
                                              # metric_robust→F5, size_robust→F6, seed_robust→F7, null_separated→F8, resto→F9
DistanceSuiteResult(estimates: Mapping[DistanceMode, DimensionEstimate], fallback_used: Mapping[DistanceMode, bool],
    resistance_exponent: float, metric_spread: float, distance_sensitive: bool, reason: str)
HomogeneityReport(mu: float, sigma: float, cv: float, n_valid: int, ok: bool)
IsotropyReport(median_ratio: float, p10_ratio: float, radius: int, k: int, n_sources: int, ok: bool)
LocalityReport(detour_fraction: float, ok: bool)
AnnulusReport(fractions: Mapping[int, float], ok: bool)
BettiResult(betti: tuple[int, ...], counts: tuple[int, ...], euler_betti: int, euler_counts: int,
    status: Literal["ok","over_budget"], field: str = "Z2")
ShortCycleBetti(b0: int, b1: int, n_edges: int, n_faces: int, b1_density: float, max_length: int, status: Literal["ok","over_budget"])
FiltrationLevel(theta: float, n_edges: int, n_components: int, giant_fraction: float, mean_degree: float)
BettiCurves(thetas: FloatArray, beta0: IntArray, beta1_short: IntArray, b1_density: FloatArray, giant_fraction: FloatArray,
    status: tuple[str, ...])
PersistenceH0(births: FloatArray, deaths: FloatArray, n_essential: int)
TopologySummary(curves: BettiCurves, h0: PersistenceH0, at_w_min: ShortCycleBetti, clique_at_w_min: BettiResult,
    stable_flags: tuple[bool, ...], stable: bool)
CurvatureSummary(mean: float, std: float, se: float, tail_fraction: float, n_edges: int, sampled: bool, edge_values: FloatArray, ok: bool)
QRCProfile(deltas: IntArray, ratio_mean: FloatArray, ratio_se: FloatArray, n_pairs: IntArray)
RunEvidence(status: RunStatus, topology: TopologyObservables, geometry: GeometryObservables, suite: DistanceSuiteResult,
    d_weyl: DimensionEstimate, homogeneity: HomogeneityReport, isotropy: IsotropyReport, locality: LocalityReport,
    annulus: AnnulusReport, topo: TopologySummary, curvature: CurvatureSummary, qrc: QRCProfile | None,
    weight_mean: float, weight_cv: float, consensus_dimension: float, dimension_class: int | None,
    clustering_ratio: float, fiedler_length: float, n: int)
RunAssessment11(codes: tuple[FailureCode, ...], primary: FailureCode | None, flags: Mapping[str, bool], passes: bool)
GeometryCertificate(connected: bool, nontrivial: bool, locality: bool, d_volume: bool, d_spectral: bool, d_weyl: bool,
    metric_robust: bool, seed_robust: bool, size_robust: bool, null_separated: bool, isotropy_ok: bool, homogeneity_ok: bool,
    topology_stable: bool, manifold_proxy_ok: bool, consensus_dimension: float, dimension_class: int | None, n_runs: int,
    evidence: Mapping[str, float | int | str | bool | None])
  # __post_init__: type(x) is bool para los 14 campos (rechaza np.bool_)
  # @property satisfied_all -> bool ; failed_fields -> tuple[str, ...]
  # @property verdict -> Verdict   (GEOMETRIC_CANDIDATE sii satisfied_all and dimension_class is not None)
PointVerdict(certificate: GeometryCertificate, codes: tuple[FailureCode, ...], primary: FailureCode | None, verdict: Verdict,
    run_outcome_fractions: Mapping[str, float])
ChainResult(engine: Engine, theta: float, samples: Mapping[str, FloatArray], sample_steps: IntArray, w_final: FloatArray,
    states: tuple[FloatArray, ...], acceptance: float, step_size: float, n_steps: int)
EnsembleSummary(means: Mapping[str,float], ses: Mapping[str,float], rhat: Mapping[str,float], tau_int: Mapping[str,float],
    ess: Mapping[str,float], geweke_ok: bool, equilibrated: bool)
TransitionSummary(lambdas: FloatArray, mean: FloatArray, derivative: FloatArray, susceptibility: FloatArray,
    binder: FloatArray, bimodality: FloatArray, peak_lambda: float)
CoarseLevel(level: int, w: FloatArray, blocks: tuple[tuple[int, ...], ...], n: int)
CoarseConsistency(classes: tuple[int | None, ...], codes: tuple[FailureCode | None, ...], measured_levels: int,
    undetermined_levels: int, consistent: bool)
```

**Otros archivos de WP-0:**
- `__init__.py` con docstring para `omega/topology/`, `omega/controls/`, `omega/certificate/` y `omega/experiments/v11/`.
- `tools/freeze_baseline.py`, `tests/expected_baseline/*.json` y `tests/test_baseline_golden.py`.
- `__version__ = "1.1.0"` en `omega/__init__.py` y `version` en `pyproject.toml` (solo metadatos; los fixtures excluyen `versions`).
- `docs/OMEGA_1_1_DESIGN.md`.

### 3.1 WP-A: analítica, ablación y Ω-B

```python
# omega/dynamics/fixed_density.py
def density_to_total(rho: float, n: int) -> float
def project_capped_simplex(v: FloatArray, total: float) -> FloatArray
def project_fixed_density(w: FloatArray, total: float) -> FloatArray
def fixed_density_step(w: FloatArray, g: FloatArray, dt: float, total: float) -> FloatArray
def evolve_fixed_density(w0: FloatArray, p: FunctionalParams, cfg: DynamicsConfig, fd: FixedDensityConfig) -> Trajectory
    # proyecta W0 primero; dt como step_size de Ω-1.0; mismos criterios de parada; escalares de step_scalars + "multiplier"; solo "clip"
def uniform_state_threshold(gamma_hat: float, rho: float, n: int) -> float   # α̂_c = (1+γ̂)(N−2)/(ρ(N−4)); n>=5
def uniform_state_unstable(alpha_hat: float, gamma_hat: float, rho: float, n: int) -> bool
```

```python
# omega/phases/ablation.py
class AblationSpec(Enum): FULL="full"; NO_TRIANGLES="no_triangles"; NO_DENSITY="no_density"; NO_DEGREE="no_degree"; TRIANGLES_ONLY="triangles_only"
@dataclass(frozen=True, slots=True)
class AblatedParams: alpha: float; beta: float; gamma: float; spec: AblationSpec   # beta>=0 permitido (FunctionalParams exige >0)
def effective_coefficients(ap: AblatedParams) -> tuple[float, float, float]
def ablated_action(w: FloatArray, ap: AblatedParams) -> float          # usa triangles/density/degree_irregularity de Ω-1.0
def ablated_gradient(w: FloatArray, ap: AblatedParams) -> FloatArray   # usa grad_* de Ω-1.0
def ablated_lipschitz(n: int, ap: AblatedParams) -> float              # 2β_e + 2(γ_e+|α_e|)(N−2); ValueError si es 0
def evolve_ablated(w0: FloatArray, ap: AblatedParams, cfg: DynamicsConfig) -> Trajectory   # bucle clip; FULL ≡ evolve() ≤1e-12
def ablation_set(alpha_hat: float, gamma_hat: float, n: int, beta: float = 1.0) -> tuple[AblatedParams, ...]
```

**Tests de WP-A:**
- `tests/test_fixed_density.py`:
  - proyección: Σ = total ≤ 1e-9, rango [0,1] y optimalidad (‖x−v‖ ≤ ‖y−v‖ para 200 puntos factibles aleatorios);
  - W = ρ(J−I) es punto fijo (≤ 1e-12);
  - umbral: N=40, ρ=0.1, α̂_c = 10.556. Con 0.9·α̂_c, una perturbación de 1e-3 vuelve (max|w−ρ| < 1e-6); con 1.3·α̂_c se aleja (> 0.1); con γ̂=1 el umbral se duplica;
  - con γ̂=0 por encima del umbral: un único clique (componente con densidad interna ≥ 0.99 y el resto aislado).
- `tests/test_ablation.py`:
  - FULL ≡ `evolve` (≤ 1e-12);
  - TRIANGLES_ONLY desde U(0,1) lleva a J−I (no decrece ninguna arista: ∂/∂w = −α(W²) ≤ 0);
  - NO_TRIANGLES lleva a 0 (E4);
  - NO_TRIANGLES con densidad fija lleva al uniforme ρ (estrictamente convexo en el conjunto afín).
- Más los 4 archivos de `tests/analytical/` (§2).

### 3.2 WP-B: geometría

```python
# omega/geometry/distance_suite.py
def edge_length_matrix_mode(w: FloatArray, a: BoolArray, mode: DistanceMode, epsilon: float, log_floor: float) -> csr_array
def resistance_matrix(w: FloatArray, a: BoolArray, max_nodes: int = 2000) -> FloatArray  # inf entre componentes
def mode_distance_matrix(w: FloatArray, w_min: float, mode: DistanceMode, graph: GraphConfig, suite: DistanceSuiteConfig) -> FloatArray
def mode_dimension(d: FloatArray, dim: DimensionConfig, suite: DistanceSuiteConfig, integer_metric: bool) -> tuple[DimensionEstimate, bool]
def resistance_exponent(r: FloatArray, hop: FloatArray, r_hi: float) -> float   # NaN si r_hi<2
def metric_spread(estimates: Mapping[DistanceMode, DimensionEstimate]) -> float  # sobre GEODESIC_MODES; inf si alguno no ok
def distance_suite(w: FloatArray, graph: GraphConfig, dim: DimensionConfig, suite: DistanceSuiteConfig) -> DistanceSuiteResult
```

`distance_suite` restringe a la gigante exactamente como `geometry_observables` (`threshold_adjacency`, `component_labels`, `giant_component_nodes`, `submatrix`). Invariantes con test:
- estimates[INV] == `geometry_observables(...).d_eff` (valor, ventana y estado) cuando este es "ok";
- estimates[HOP] == `d_eff_hops`.

```python
# omega/geometry/weyl.py
def laplacian_spectrum(w: FloatArray, kind: Literal["combinatorial","normalized"]) -> FloatArray
def weyl_staircase(eigs: FloatArray, rtol: float) -> tuple[FloatArray, FloatArray]   # (niveles, cuentas) sin el modo cero
def weyl_dimension(w: FloatArray, w_min: float, cfg: WeylConfig) -> DimensionEstimate  # plateau := R²>=r2_min; method "weyl_<lap>_<graph>"
def fiedler_length(w: FloatArray, w_min: float) -> float   # λ2^{-1/2}, Laplaciano combinatorio binario sobre la gigante
```

```python
# omega/geometry/local_structure.py   (a, d_hop: YA restringidos a la gigante)
def node_dimensions(d_hop: FloatArray, radii: FloatArray) -> FloatArray
def homogeneity(d_hop: FloatArray, hop_est: DimensionEstimate, cfg: LocalStructureConfig) -> HomogeneityReport
def ball_mds_ratio(d_hop: FloatArray, center: int, radius: int, k: int) -> float
def isotropy(d_hop: FloatArray, hop_est: DimensionEstimate, k: int | None, cfg: LocalStructureConfig, rng: np.random.Generator) -> IsotropyReport
def edge_detour_fraction(a: BoolArray) -> float
def locality(a: BoolArray, cfg: LocalStructureConfig) -> LocalityReport
def annulus_connectivity(a: BoolArray, d_hop: FloatArray, hop_est: DimensionEstimate, cfg: LocalStructureConfig, rng: np.random.Generator) -> AnnulusReport
```

**Tests de WP-B:**
- `tests/test_distance_suite.py`:
  - en retículas con W=1, HOP=INV=LOG exactamente;
  - RESISTANCE: anillo de 300 da 1±0.05; 17² y 7³ dan `no_window`;
  - ζ: anillo > 0.8; 20² ∈ [0.35, 0.65]; 7³ < 0.35;
  - en 9³ con pesos U(0.3,1), la dispersión ≤ 0.5 y `distance_sensitive` es False;
  - un W sintético donde solo INV tiene ventana da True;
  - igualdad bit a bit con Ω-1.0;
  - la resistencia satisface la desigualdad triangular.
- `tests/test_weyl.py`: tabla rápida (anillo 300, 17², 28², 7³, 9³ con tolerancias de §1.3); K_n da `no_window`; U(0,1) N=100 da > 5; la escalera es monótona y sin modo cero.
- `tests/test_local_structure.py`:
  - detour: el anillo da 0 (C_n con n ≥ 5); 9³, K_n y la retícula triangular dan 1; la fórmula coincide con BFS en G−e en 5 grafos aleatorios pequeños;
  - homogeneidad: 9³ da cv = 0;
  - isotropía: 9³ ≥ 0.99; slab 32×5×5 da ok=False;
  - anillos: 9³ ok; árbol 3-ario ok=False.

### 3.3 WP-C: topología y curvatura

```python
# omega/topology/filtration.py
def threshold_graph(w: FloatArray, theta: float) -> BoolArray            # W > θ
def descending_edges(w: FloatArray) -> tuple[IntArray, IntArray, FloatArray]
def filtration_levels(w: FloatArray, thetas: Sequence[float]) -> tuple[FiltrationLevel, ...]
# omega/topology/betti.py
def gf2_rank(columns: Sequence[int]) -> int
def count_triangles(a: BoolArray) -> int                                  # tr(A³)/6
def count_four_cycles(a: BoolArray) -> int                                # (tr A⁴ − 2Σk² + 2m)/8
def clique_complex_betti(a: BoolArray, max_dim: int = 2, max_simplices: int = 300_000) -> BettiResult
def short_cycle_betti1(a: BoolArray, max_length: int = 4, max_faces: int = 500_000) -> ShortCycleBetti
# omega/topology/persistence.py
def h0_persistence(w: FloatArray) -> PersistenceH0
def betti_curves(w: FloatArray, cfg: TopologyConfig) -> BettiCurves
def topology_stable(w: FloatArray, thetas: Sequence[float], g_connected: float, cfg: TopologyConfig) -> tuple[bool, tuple[bool, ...]]
def topology_summary(w: FloatArray, w_min: float, sensitivity: Sequence[float], g_connected: float, cfg: TopologyConfig) -> TopologySummary
# omega/curvature/ollivier.py
def wasserstein1(mu: FloatArray, nu: FloatArray, cost: FloatArray) -> float
def neighbor_measure(w: FloatArray, a: BoolArray, x: int, idleness: float) -> tuple[IntArray, FloatArray]
def ollivier_edge(w: FloatArray, a: BoolArray, d: FloatArray, x: int, y: int, idleness: float) -> float
def ollivier_curvature(w: FloatArray, w_min: float, graph: GraphConfig, suite: DistanceSuiteConfig, cfg: CurvatureConfig,
                       rng: np.random.Generator) -> CurvatureSummary
# omega/curvature/qrc.py
def average_sphere_distance(d: FloatArray, x: int, y: int, delta: int) -> float
def quantum_ricci_profile(d: FloatArray, deltas: Sequence[int], max_pairs: int, rng: np.random.Generator) -> QRCProfile
```

**Tests de WP-C:**
- Betti del complejo de cliques: K4 da (1,0,0); C5 da (1,1,0); el grafo del octaedro da (1,0,1) con χ=2; dos triángulos disjuntos dan β0=2.
- β1^(4): T² 6×6 da 2; T³ 4³ da 3; anillo da 1; K_n da 0.
- `count_four_cycles` coincide con la enumeración.
- `over_budget` se activa con K_60 y presupuesto bajo.
- H0: número de barras = N−1 en un grafo conexo; las muertes son las del árbol generador máximo (frente a networkx).
- Curvas: monótonas (β0 no crece al bajar θ).
- Ollivier: T² 6×6 y T³ 4³ dan κ ≡ 0 (≤ 1e-9); C_8 da 0; K_5 y K_10 cumplen la forma cerrada; `wasserstein1` coincide con `linear_sum_assignment` en medidas uniformes del mismo tamaño.
- `@slow`: S² > T² (RGG N=300, k=20, 5 semillas por geometría; diferencia > 2·SE combinado de semillas).
- QRC: en el anillo, d̄/δ es exacto a mano en un caso pequeño.

### 3.4 WP-D: motores estadísticos

```python
# omega/statistics/langevin.py
def reflect_unit(x: FloatArray) -> FloatArray
def temperature(theta_hat: float, beta: float) -> float
def langevin_step(u: FloatArray, n: int, p: FunctionalParams, dt: float, theta: float, rng: np.random.Generator,
                  boundary: Literal["reflect","project"]) -> FloatArray
def chain_observables(w: FloatArray, p: FunctionalParams, w_min: float) -> dict[str, float]
    # action, mean_weight, binary_density (W>w_min), triangle_density = T/C(N,3), strength_cv
def run_langevin(w0: FloatArray, p: FunctionalParams, cfg: LangevinConfig, w_min: float, rng: np.random.Generator,
                 n_states: int = 0) -> ChainResult
# omega/statistics/metropolis.py
def delta_action_edge(w: FloatArray, a: int, b: int, new_value: float, p: FunctionalParams, k: FloatArray) -> float
    # δ=w'-w: −αδ(W²)_ab + β(w'²−w²) + γ[(2δk_a+δ²)+(2δk_b+δ²) − ((K+2δ)²−K²)/N] − μδ ; K=Σk ; exige η=0
def metropolis_sweep(w: FloatArray, p: FunctionalParams, theta: float, step: float, rng: np.random.Generator) -> tuple[FloatArray, int]
def run_metropolis(w0: FloatArray, p: FunctionalParams, cfg: MetropolisConfig, w_min: float, rng: np.random.Generator,
                   n_states: int = 0) -> ChainResult
# omega/statistics/ensemble.py
def autocorrelation(x: FloatArray) -> FloatArray                          # FFT
def integrated_autocorr_time(x: FloatArray, c: float = 5.0) -> float
def effective_sample_size(x: FloatArray, c: float = 5.0) -> float
def split_rhat(chains: Sequence[FloatArray]) -> float
def geweke_z(x: FloatArray, first: float = 0.1, last: float = 0.5, c: float = 5.0) -> float
def post_burn_in(x: FloatArray, fraction: float) -> FloatArray
def ensemble_summary(chains: Sequence[ChainResult], keys: Sequence[str], cfg: EnsembleConfig, burn_in_fraction: float) -> EnsembleSummary
# omega/statistics/transition.py
def fluctuation(x: FloatArray) -> float
def susceptibility(m: FloatArray, n_edges: int) -> float
def binder_cumulant(m: FloatArray) -> float
def bimodality_coefficient(m: FloatArray) -> float
def order_parameter_histogram(m: FloatArray, bins: int = 40) -> tuple[FloatArray, FloatArray]
def transition_summary(samples: Mapping[float, FloatArray], n_edges: int) -> TransitionSummary
```

**Tests de WP-D:**
- `delta_action_edge` coincide con action(w') − action(w) (rel 1e-10) en 200 pruebas.
- Gibbs de aristas independientes (α=γ=0, N=8, M=28): la media de w coincide con la cuadratura de ∫w·e^{−βw²/Θ} en [0,1] (3·SE vía ESS) para Θ̂ ∈ {0.1, 1}, con Langevin (reflect) y con Metropolis.
- "project" produce átomos (frac(w==0) > 0) y "reflect" no.
- Caso interactivo pequeño (N=6, α̂=1, γ̂=1, Θ̂=0.3): Langevin y Metropolis concuerdan (3·SE).
- τ_int de AR(1) con φ=0.9 y n=2e5: ≈ 9.5 ± 10%.
- R̂ de cadenas iid ≤ 1.01; con medias distintas > 1.1.
- Binder de una gaussiana ≈ 0 (±0.02); bimodalidad: bimodal > 5/9 > unimodal.
- Determinismo por semilla.

### 3.5 WP-E: controles, regla de coarse-graining y procedencia

```python
# omega/controls/erdos_renyi.py
def erdos_renyi_gnp(n: int, p: float, rng: np.random.Generator) -> FloatArray
def erdos_renyi_gnm(n: int, m: int, rng: np.random.Generator) -> FloatArray
def matched_erdos_renyi(w: FloatArray, w_min: float, rng: np.random.Generator, weights: Literal["binary","shuffled"] = "binary") -> FloatArray
# omega/controls/configuration_model.py
def configuration_model(degrees: IntArray, rng: np.random.Generator) -> FloatArray
def matched_configuration_model(w: FloatArray, w_min: float, rng: np.random.Generator) -> FloatArray
# omega/controls/degree_preserving_rewire.py
def degree_preserving_rewire(w: FloatArray, w_min: float, rng: np.random.Generator, swaps_per_edge: int = 10) -> FloatArray
# omega/controls/small_world.py
def ring_lattice(n: int, k: int) -> FloatArray
def watts_strogatz(n: int, k: int, p: float, rng: np.random.Generator) -> FloatArray
def matched_small_world(w: FloatArray, w_min: float, p: float, rng: np.random.Generator) -> FloatArray
# omega/controls/random_geometric.py
def rgg_torus(n: int, dim: int, mean_degree: float, rng: np.random.Generator, weights: Literal["binary","euclidean"] = "binary") -> FloatArray
def cap_angle(n: int, dim: int, mean_degree: float) -> float
def rgg_sphere(n: int, dim: int, mean_degree: float, rng: np.random.Generator) -> FloatArray
def balanced_tree(branching: int, height: int) -> FloatArray
# omega/controls/battery.py
def null_seed_key(base: SeedKey, kind: NullModel) -> SeedKey      # spawn_key + (1000 + índice del enum,)
def make_null(w: FloatArray, w_min: float, kind: NullModel, rng: np.random.Generator) -> FloatArray
def null_battery(w: FloatArray, w_min: float, kinds: Sequence[NullModel], base: SeedKey) -> dict[NullModel, FloatArray]
# omega/coarse_graining/rule.py
def heavy_edge_matching(w: FloatArray, w_min: float, rng: np.random.Generator) -> tuple[tuple[int, ...], ...]
def aggregate(w: FloatArray, blocks: Sequence[Sequence[int]], how: Literal["max","mean"]) -> FloatArray
def coarse_grain(w: FloatArray, w_min: float, cfg: CoarseGrainConfig, rng: np.random.Generator, level: int = 1) -> CoarseLevel
def coarse_grain_hierarchy(w: FloatArray, w_min: float, cfg: CoarseGrainConfig, rng: np.random.Generator) -> tuple[CoarseLevel, ...]
# omega/config/convert11.py
def config11_to_dict(cfg: Omega11Config) -> dict[str, Any]     # enums -> .value; usa config_to_dict de Ω-1.0 para base
def config11_from_dict(d: Mapping[str, Any]) -> Omega11Config   # rechaza claves desconocidas; schema "1.1"
def canonical_json(obj: Any) -> str
def config_hash(cfg: Omega11Config) -> str
# omega/io/provenance.py
TERMINATION_REASONS: Final[frozenset[str]]
def dependency_lock(packages: Sequence[str] = ("numpy","scipy","networkx")) -> dict[str, str]   # incluye "python"
def dependency_lock_hash(lock: Mapping[str, str]) -> str
def git_info(repo: Path) -> tuple[str, str]
def format_experiment_id(number: int) -> str      # "Ω-EXP-004821"
def experiment_file_stem(number: int) -> str      # "OMEGA-EXP-004821"
def parse_experiment_id(label: str) -> int
def allocate_experiment_id(registry_dir: Path) -> int
def build_passport_v11(cfg: Omega11Config, *, experiment_number: int, seed: SeedKey, termination_reason: str,
    initialization_distribution: str, null_model: NullModel | None, coarse_graining: bool, results: Mapping[str, Any],
    code_commit: str, git_branch: str, entrypoint: str, w0_sha256: str | None,
    omega10_passport: Mapping[str, Any] | None) -> dict[str, Any]
def save_passport_v11(passport: Mapping[str, Any], arrays: Mapping[str, FloatArray], directory: Path) -> tuple[Path, Path]
def load_passport_v11(json_path: Path) -> tuple[dict[str, Any], dict[str, FloatArray]]   # verifica sha256 (reusa array_digest de Ω-1.0)
def reconstruct(passport: Mapping[str, Any]) -> tuple[Omega11Config, SeedKey]
```

**Tests de WP-E:**
- `tests/test_controls.py`:
  - G(n,m) tiene exactamente m aristas;
  - el recableado conserva grados y el multiconjunto de pesos;
  - CM: grados ≤ los objetivo y ≥ 95% de los stubs;
  - WS con p=0 es la retícula en anillo exacta;
  - `rgg_sphere` tiene grado medio dentro de ±10%;
  - determinismo por semilla;
  - `null_seed_key` es distinto por tipo.
- `tests/test_coarse_rule.py`:
  - los bloques forman una partición de tamaño ≤ 2;
  - el emparejamiento es maximal (ninguna arista entre singletons);
  - max sobre binario = `networkx.quotient_graph`;
  - el anillo de 800 mantiene D=1.00±0.02 en todos los niveles;
  - equivariancia exacta bajo permutación con pesos distintos (sin empates).
- `tests/test_convert11.py`: ida y vuelta; hash estable ante el orden de claves; hash distinto si cambia un umbral.
- `tests/test_provenance.py`:
  - estén presentes los 14 campos de P§17;
  - asignación concurrente (8 hilos) da IDs únicos;
  - formato y parseo;
  - un digest alterado lanza excepción;
  - `reconstruct` y una nueva ejecución dan digests iguales (simulación Ω-1.0 con N=12).

### 3.6 WP-F: certificado y taxonomía (lógica pura)

```python
# omega/certificate/taxonomy.py
def consensus_dimension(estimates: Sequence[DimensionEstimate], class_tol: float) -> tuple[float, int | None]
def resistance_consistent(zeta: float, dimension_class: int | None, suite: DistanceSuiteConfig) -> bool
def run_flags(ev: RunEvidence, cfg: Omega11Config) -> dict[str, bool]
    # claves: converged, empty, dense_or_uniform, nontrivial, connected, small_world, locality, d_volume, d_spectral, d_weyl,
    #         metric_robust, isotropy_ok, homogeneity_ok, topology_stable, manifold_proxy_ok
def run_failure_codes(flags: Mapping[str, bool]) -> tuple[FailureCode, ...]
def primary_code(codes: Iterable[FailureCode]) -> FailureCode | None
def assess_run(ev: RunEvidence, cfg: Omega11Config) -> RunAssessment11
# omega/certificate/certificate.py
RULE_D3_NEVER_SUFFICIENT: Final = "D≈3 jamás basta por sí sola: GEOMETRIC_CANDIDATE exige los 14 campos del certificado"
def seed_robust(runs: Sequence[RunEvidence], cfg: Omega11Config) -> tuple[bool, dict[str, float]]
def size_robust(by_size: Mapping[int, Sequence[RunEvidence]], cfg: Omega11Config) -> tuple[bool, dict[str, float]]
def null_separated(cand: Sequence[RunEvidence], nulls: Mapping[NullModel, Sequence[RunEvidence]], cfg: Omega11Config) -> tuple[bool, dict[str, float]]
def build_certificate(runs: Sequence[RunEvidence], by_size: Mapping[int, Sequence[RunEvidence]],
                      nulls: Mapping[NullModel, Sequence[RunEvidence]], cfg: Omega11Config) -> GeometryCertificate
def point_verdict(runs: Sequence[RunEvidence], by_size: Mapping[int, Sequence[RunEvidence]],
                  nulls: Mapping[NullModel, Sequence[RunEvidence]], cfg: Omega11Config,
                  ensemble: EnsembleSummary | None = None) -> PointVerdict
```

**Tests de WP-F:**
- `tests/factories_v11.py` (sin prefijo `test_`; se importa como `tests.factories_v11`): `make_evidence(**overrides) -> RunEvidence` con un estado "pasa todo" (D_vol = D_s = D_W = 3.0, etc.).
- `tests/test_taxonomy.py`:
  - cada bandera aislada produce su código;
  - precedencia;
  - un no-convergido da F10;
  - el uniforme da F1;
  - clase None hace fallar metric_robust.
- `tests/test_certificate.py`:
  - las 2^14 combinaciones: el veredicto es candidato solo si todo es True;
  - **D=3 exacto en los tres estimadores con locality False da NOT_CANDIDATE y primario F3**;
  - evidencia ausente da False;
  - `size_robust` acepta la serie RGG (3.06, 3.11, 3.16) y rechaza (2.6, 3.0, 3.4);
  - `null_separated` rechaza un nulo que pasa;
  - `np.bool_` es rechazado.
- `tests/test_architecture_v11.py`:
  - regex `(?i)\b(d3|is_?3d|three_?d)\b` sin coincidencias en identificadores AST de `omega/`;
  - `GeometryCertificate` solo se instancia en `omega/certificate/` (AST).

### 3.7 WP-G (oleada 2): integración

```python
# omega/certificate/evidence.py
def evidence_rng(key: SeedKey) -> np.random.Generator     # make_rng(SeedKey(e, spawn_key + (2000,)))
def weight_stats(w: FloatArray) -> tuple[float, float]     # (media, cv) del triángulo superior
def collect_run_evidence(w: FloatArray, status: RunStatus, cfg: Omega11Config, rng: np.random.Generator, *, with_qrc: bool = False) -> RunEvidence
def evidence_summary(ev: RunEvidence) -> dict[str, Any]
# omega/coarse_graining/analysis.py
def hierarchy_evidence(w: FloatArray, cfg: Omega11Config, rng: np.random.Generator) -> tuple[RunEvidence, ...]
def class_consistency(levels: Sequence[RunEvidence], cfg: Omega11Config) -> CoarseConsistency
# omega/phases/finite_size.py
def size_config(cfg: Omega11Config, n: int) -> Omega11Config
def fss_observables(ev: RunEvidence) -> dict[str, float]     # D*, D_vol, D_s, D_W, C_bin, G/N, rho, xi_F, L_hop
def fss_fit(sizes: Sequence[int], values: Sequence[float]) -> dict[str, float]
# omega/experiments/v11/gate.py
STEP_ORDER: Final = ("o00_validation","o01_baseline","o02_distance","o03_nulls","o04_ablation","s06_phase_diagram",
                     "o06_finite_size","o05_ensemble","o08_topology","o10_curvature","o11_coarse")
class PrerequisiteError(RuntimeError): ...
def write_summary(out_root: Path, name: str, data: Mapping[str, Any], *, mode: Literal["smoke","full"]) -> Path
def require_prerequisites(out_root: Path, step: str) -> None   # exige out_root/<prev>/summary.json con mode=="full" y complete==True
```

---

## 4. Experimentos Ω-0…Ω-11 (`omega/experiments/v11/`)

**Convenciones:**
- Cada archivo expone `run(cfg: Omega11Config, out_root: Path, *, mode: Literal["smoke","full"]) -> Path` y escribe `out_root/<nombre>/summary.json` más pasaportes v1.1.
- `test_smoke` dura ≤ 30 s con N pequeño.
- `@pytest.mark.slow test_full_run` llama primero a `require_prerequisites` (salvo O-00).
- Raíz de salida: `os.environ.get("OMEGA_RUNS_DIR", "runs/omega11")`, ignorada por git.

| Archivo | Nivel | Qué mide / afirma | Smoke | Full (`@slow`) |
|---|---|---|---|---|
| `test_o00_validation_v11.py` | Ω-0 | Tablas de §1.3 (Weyl), §1.1 (suite y ζ), §1.6 (κ=0, K_n, S²>T²), §1.2 (β), §1.4–1.5 (calibraciones), y expectativas a nivel de corrida: RGG3 N=800 k12 binario pasa todos los campos de corrida; WS p=0.05 da primario F3; ER k12 da F4; árbol 3-ario no pasa; K_n da F1; vacío da F0; slab 32×5×5 no pasa | anillo 100, T² 10×10, T³ 6³, K_10, octaedro, árbol pequeño | tablas completas, 3 semillas RGG |
| `test_o01_baseline.py` | Ω-1 | S0 de Ω-1.0 con taxonomía v1.1. Con N=200, α̂ ∈ {0.5, 1, 1.5, 1.9, 2.1, 3}, γ̂ ∈ {0, 1, 10} y 10 semillas: α̂ ≤ 1.9 da F0 en 100%; α̂ ≥ 2.1 da F1 en 100% (registro negativo) | N=32, 2 semillas, α̂ ∈ {1.5, 3} | malla indicada |
| `test_o02_distance.py` | Ω-2 | Fracción `distance_sensitive` y dispersión por familia: finales Ω-1.0 (binarios: 0 por construcción, se informa), finales del integrador sigmoide, U(0,1) estático con w_min de sensibilidad, estados Langevin con Θ̂ ∈ {0.1, 0.3}, toros con pesos aleatorios (no sensibles), RGG euclidiano (sensible con N=800, documentado) | N=40 | N=200–800 |
| `test_o03_nulls.py` | Ω-3 | Batería: ER (k̄ 6/12/30), CM y recableado (de RGG3), WS (p ∈ {0.01, 0.05, 0.1, 0.2, 1}), RGG2/3 y árbol. Se registran códigos y la tasa de falsos positivos de cada campo; RGG3 debe pasar a nivel de corrida | N=100 | N=800, 5 semillas |
| `test_o04_ablation.py` | Ω-4 | Las cinco AblationSpec × α̂ ∈ {0.5, 1.5, 2.5, 4} × γ̂ ∈ {0, 1, 10}. **T=0 da F0 siempre (analítico)**; TRIANGLES_ONLY da F1. Ω-B: ρ ∈ {0.05, 0.1, 0.2} × α̂ ∈ {0.5, 1.5, 3}·α̂_c(γ̂) × γ̂ ∈ {0, 1, 10}: se afirma uniforme (F1) bajo el umbral y clique más aislados (F2) con γ̂=0 por encima | N=24 | N=100/200, 10 semillas |
| `test_s06_phase_diagram_v11.py` | paso 6 (Exp 3) | Malla reducida D-8 con N=200 y 10 semillas: distribución de códigos, certificado por punto y λ* de seed_variance. **Gate: O-00..O-04 en modo full** | N=24, 2×2 puntos (sin gate) | completa |
| `test_o06_finite_size.py` | Ω-6 (paso 7) | N ∈ (64…800) en α̂ ∈ {1.9, 2.1, 3} × γ̂ ∈ {0, 10}, puntos Ω-B y puntos que pasen; RGG3 k12 en N ∈ (800, 1500, 3000) como control del mecanismo `size_robust`. Curvas D(N), C(N), G/N, ρ(N), ξ_F(N) y L_hop(N) | N ∈ (24, 32, 48) | N completo; XL solo para candidatos |
| `test_o05_ensemble.py` | Ω-5 (paso 8) | Langevin y Metropolis, 4 cadenas: N ∈ {64, 100}, Θ̂ ∈ {0.01, 0.03, 0.1, 0.3, 1}, α̂ ∈ {0, 0.5, 1, 1.5, 2, 2.5, 3, 4}, γ̂ ∈ {0, 1, 10}. R̂, τ_int, distribución de códigos, ⟨m⟩, χ, U4, bimodalidad, λ*(Θ̂) e histéresis. Pregunta: ¿hay una región dominante con propiedades geométricas? | N=16, 200 sweeps | completa |
| `test_o07_dimensions.py` | Ω-7 | D_vol, D_hop, D_log, D_ball, D_s y D_Weyl con curvas por escala; acuerdo por pares en estados Ω y controles | N=40 | N=800 |
| `test_o08_topology.py` | Ω-8 (paso 9) | Curvas de Betti(θ), barras H0, β1^(4) y complejo de cliques (con χ) en estados no binarios y controles; `topology_stable` | N=40 | N=800 |
| `test_o09_manifold.py` | Ω-9 | Homogeneidad, isotropía, localidad, anillos y `MANIFOLD_PROXY` en estados Ω y controles | N=60 | N=800 |
| `test_o10_curvature.py` | Ω-10 (paso 9) | Ollivier (media y cola) y QRC solo informe; S² vs T² repetido como control | 50 aristas | completo |
| `test_o11_coarse.py` | Ω-11 (paso 10) | Jerarquía (3 réplicas) en anillo, 2D, 3D, RGG y estados Ω; `CoarseConsistency` | anillo 200, 16² | N=800 y XL |

**Orden de P§22, codificado en el gate:**
1. analítica (pytest);
2. golden;
3. O-02;
4. O-03;
5. O-04;
6. s06 (bloqueado sin los pasos 1–5 en modo full);
7. O-06;
8. O-05;
9. O-08 y O-10;
10. O-11.

O-00 y O-01 son prerrequisitos de O-02. O-07 y O-09 no tienen orden propio: se ejecutan con 9.

---

## 5. Paquetes de trabajo (archivos exclusivos)

**WP-0** lo hace la sesión principal, en secuencia y antes de lanzar nada: §3.0 completo, golden y manifiesto.
- Primero se generan los fixtures sobre el árbol intacto; después se sube la versión.
- Sus pruebas deben pasar (`pytest -q`: 233 + golden) y `mypy --strict omega` debe estar limpio.
- Desde ese momento, `settings11.py`, `contracts.py` y `baseline.py` quedan CONGELADOS. Un cambio exige decisión del orquestador y una nueva versión del documento.

**Oleada 1** (6 agentes Sonnet en paralelo; ningún archivo compartido):

| WP | Archivos exclusivos | Depende de |
|---|---|---|
| A Analítica, ablación, Ω-B | `tests/analytical/{test_empty_state,test_complete_state,test_uniform_manifold,test_parameter_scaling}.py`, `omega/dynamics/fixed_density.py`, `omega/phases/ablation.py`, `tests/{test_fixed_density,test_ablation}.py` | WP-0, Ω-1.0 |
| B Geometría | `omega/geometry/{distance_suite,weyl,local_structure}.py`, `tests/{test_distance_suite,test_weyl,test_local_structure}.py` | WP-0 |
| C Topología y curvatura | `omega/topology/{filtration,betti,persistence}.py`, `omega/curvature/{ollivier,qrc}.py`, `tests/{test_filtration,test_betti,test_persistence,test_ollivier,test_qrc}.py` | WP-0. Las distancias de Ollivier se calculan internamente (HOP con `hop_distance_matrix`, INV con `distance_matrix` de Ω-1.0). Para LOG y RESISTANCE se usa un import diferido de `omega.geometry.distance_suite.mode_distance_matrix`, solo si `distance_mode` ≠ HOP/INV; los tests de WP-C usan HOP |
| D Motores estadísticos | `omega/statistics/{langevin,metropolis,ensemble,transition}.py`, `tests/{test_langevin,test_metropolis,test_ensemble,test_transition}.py` | WP-0 |
| E Controles, coarse y procedencia | `omega/controls/{erdos_renyi,configuration_model,degree_preserving_rewire,small_world,random_geometric,battery}.py`, `omega/coarse_graining/rule.py`, `omega/config/convert11.py`, `omega/io/provenance.py`, `tests/{test_controls,test_coarse_rule,test_convert11,test_provenance}.py` | WP-0 |
| F Certificado y taxonomía | `omega/certificate/{taxonomy,certificate}.py`, `tests/{factories_v11,test_taxonomy,test_certificate,test_architecture_v11}.py` | WP-0 (solo contracts y settings11) |

**Oleada 2** (tras fusionar A–F; puede partirse en G1 y G2):
- **G1:** `omega/certificate/evidence.py`, `omega/experiments/v11/{gate,test_o00_validation_v11,test_o01_baseline,test_o02_distance,test_o03_nulls,test_o04_ablation,test_s06_phase_diagram_v11}.py`, `tests/test_evidence.py`.
- **G2:** `omega/coarse_graining/analysis.py`, `omega/phases/finite_size.py`, `omega/experiments/v11/{test_o05_ensemble,test_o06_finite_size,test_o07_dimensions,test_o08_topology,test_o09_manifold,test_o10_curvature,test_o11_coarse}.py`, `tests/{test_finite_size,test_coarse_analysis}.py`.
- G2 programa contra la firma congelada de `collect_run_evidence`; G1 se fusiona primero.
- `tests/test_evidence.py` exige:
  - que el RGG3 binario (N=800) dé `assess_run(...).passes`;
  - que WS p=0.05 dé primario F3;
  - igualdad de `evidence.geometry` con `scan.observe` de Ω-1.0;
  - determinismo por SeedKey;
  - un tiempo de evidencia con N=800 < 20 s.

**Criterios de aceptación** (cada WP):
1. `python -m pytest -q <sus tests>` pasa.
2. `python -m pytest -q` completo (rápido) pasa, incluidos `test_baseline_golden` y los 233 de Ω-1.0.
3. `python -m mypy --strict omega` está limpio.
4. `python -m mypy --strict --explicit-package-bases <sus archivos de test nuevos>` está limpio.
5. `git diff --name-only` ⊆ su lista de archivos.
6. Ningún uso de `np.random` global; sin E/S fuera de `io/` y `experiments/`.
7. Los smoke tests nuevos de su paquete suman ≤ 90 s.

Las tolerancias numéricas y los umbrales son exactamente los de este documento: ningún agente los ensancha. Si un valor de validación no se reproduce, el agente lo informa y no lo ajusta.

---

## 6. Riesgos y decisiones abiertas (requieren decisión del usuario)

- **R1. Tamaño finito frente a coste. Es la decisión principal.**
  - En 3D, la ventana de escala solo existe con N ≥ 800 (RGG3 k12: no_window con N=300 y 500). Con N ≤ 800, `size_robust`, la coherencia de coarse-graining y la señal de curvatura de esfera (≈1%) no pueden certificarse.
  - Propuesta: mantener N=64…800 (P§15) para tendencias y añadir N ∈ {1500, 3000} solo para puntos que pasen a nivel de corrida con N=800. Coste: dinámica densa O(N³) de ~1 s/paso con N=3000, es decir, horas por corrida.
  - Alternativa: backend disperso en Ω-1.2. Sin esta decisión, el resultado esperable de Ω-1.1 en 3D es como mucho NOT_CANDIDATE con F6.
- **R2. RESISTANCE.** D_resistance por volumen está indefinida en d ≥ 2 (es un teorema), así que no puede "sobrevivir" como pide literalmente P§4. Se propone informarla y que vote mediante ζ_R (clase-consistencia). Necesita aprobación.
- **R3. Umbrales calibrados en controles.**
  - Los umbrales se calibraron sobre RGG3 k12, WS, ER, árboles y slabs a N ≤ 800: locality 0.97, cola de curvatura 0.03, b1/E 0.03, isotropía 0.5/0.3, anillos 0.9, tolerancia dimensional 0.35 y métrica 0.5.
  - Son estrictos: el RGG3 con k=8 (más disperso) falla la prueba de anillos (0.25) y el RGG euclidiano es distance-sensitive con N=800.
  - Se congelan antes del paso 6. Pedir confirmación de que se acepta este rigor.
- **R4. Isotropía por MDS local.** Usa una inmersión local solo como diagnóstico de una métrica ya medida, nunca como entrada (P§13 previene sobre inmersiones). Necesita aprobación; si no, isotropy_ok queda sin proxy válido y el certificado nunca podría pasar.
- **R5. Inconsistencia heredada de umbral en Ω-1.0.** `topology.adjacency` usa W ≥ w_min y `threshold_adjacency` usa W > w_min. No se corrige (Ω-1.0 congelado): se documenta en `threshold_rule` y Ω-1.1 usa > en todo lo nuevo. ¿Se acepta?
- **R6. Versión y regla de size_robust.**
  - La subida a 1.1.0 cambia el metadato `versions.omega` de los pasaportes v1.0 que se generen desde ahora (no su contenido científico).
  - La regla de deriva de Ω-1.1 (pendiente ≤ 0.1 por e-fold) es distinta de la de Ω-1.0 (deriva monótona), que rechazaría al control positivo; `confirm_geometric` de Ω-1.0 no se toca.
- **R7. Ω-B estocástica.** Langevin o Metropolis sobre {Σw = W0} no están en el alcance (exigiría medida en la variedad restringida). ¿Basta Ω-B determinista para Ω-1.1?
- **R8. Expectativa honesta.**
  - La dinámica S0 con clip produce estados binarios F0/F1 en todo el plano (α̂, γ̂) (analítico más Ω-1.0).
  - Con T=0 no hay estructura (F0 o F1 uniforme).
  - En Ω-B, con γ̂=0, la predicción es un clique (F2).
  - La única vía abierta a una región no trivial es Θ > 0 (Ω-5) o Ω-B con γ̂ > 0 cerca del umbral α̂ρ ≈ 1+γ̂. Se recomienda preregistrar esta expectativa en el diario de resultados antes del paso 6.

**Prototipos** (fuera del repo): `/tmp/claude-0/-home-user-FichaClaraV2/84579788-2356-5174-b7ed-bc9de278c22a/scratchpad/p1.py` … `p22.py`.