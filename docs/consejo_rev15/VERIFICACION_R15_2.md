# Informe R15-2: verificación de citas y fórmulas para el preregistro R14-C

Verificador de fuentes, 19/10/2026. No se ha tocado el repositorio. Scripts y salidas en
`/tmp/claude-0/-home-user-FichaClaraV2/90c7523c-b9d2-5eee-afba-260b3848fae0/scratchpad/verif_r152/`
(`chains.py`, `chains_out.txt`; descargas en `dl/`).

Estados usados:
- **V-P**: verificado en la fuente primaria (texto completo del artículo).
- **V-A**: verificado solo en el abstract o en metadatos (arXiv, Crossref).
- **V-S**: verificado solo en fuente secundaria (revisión que lo cita).
- **NE**: no encontrado.
- **OTRA**: la fuente dice otra cosa, o matiza de forma relevante lo atribuido.

## A. Citas

| # | Cita | Existencia y datos bibliográficos | Lo que se le atribuye | Estado |
|---|---|---|---|---|
| 1 | Müller, arXiv:2503.01719 | Olaf Müller, «On the Hauptvermutung of Causal Set Theory». v1 3/3/2025, v2 29/12/2025, v3 24/8/2026, **v4 3/9/2026** (coincide con «versión del 3 de septiembre de 2026»). | Teorema 3: «la Hauptvermutung finita es falsa para d⁻»; Teorema 6: «la Hauptvermutung uniforme a escala de Planck es verdadera para d×»; Teorema 2: «la Hauptvermutung numerable es verdadera». La introducción advierte que estos resultados «muy probablemente aún no corresponden exactamente» a lo que la teoría espera. | **V-P** (HTML de v4). Matiz: la versión «verdadera» es la *uniforme a escala de Planck* con la distancia d×; conviene nombrarla así. |
| 2 | Eichhorn, Mack, Le, Wagner (2026), arXiv:2605.27514 | Astrid Eichhorn, Harald Mack, Kim Tuyen Le, Fabian Wagner, «Charting causal set configuration space with graph observables», gr-qc, 26/5/2026. | Abstract: nueve clases; la distribución de grado de links, los autovalores del Laplaciano del Hasse simetrizado y la abundancia de intervalos las distinguen. | **V-P con matiz (OTRA parcial)**. El cuerpo matiza: (a) solo el espectro del Laplaciano separa las nueve clases; la distribución de grado y las abundancias separan bien retículos, cuasicristales, órdenes por capas y órdenes de conectividad fija, pero detectan cambios de topología o de frontera solo cuando son fuertes (las abundancias son «pobres» para eso); (b) las clases manifold-like son sprinklings en **d = 2** (caja 1+1), n = 2048, 10⁴ causets por clase; la dependencia con la dimensión solo se mira en el apéndice B para la clase trivial; (c) las comparaciones son sobre todo de cada clase contra la clase nula (topología trivial), no por pares completos. Extrapolar a Ω en d = 4 no está respaldado por el artículo. |
| 3 | Boguñá y Krioukov | Marián Boguñá, Dmitri Krioukov, «Local d'Alembertian for causal sets», arXiv:2506.18745 (23/6/2025), **Phys. Rev. D 113, 024046 (2026)**, doi 10.1103/np89-qzjp (Crossref: publicado 22/1/2026). | Operador d'Alembertiano local construido con distancias temporales y espaciales intrínsecas del causet; converge al continuo para campos arbitrarios en Minkowski; los d'Alembertianos no locales previos fallan para algunos campos. | **V-A**. Falta el número de artículo en el acta: añadir «113, 024046». |
| 4 | Carlip (2024) | Steven Carlip, «Causal sets and an emerging continuum», **Gen. Relativ. Gravit. 56, 95 (2024)**, doi 10.1007/s10714-024-03281-1; arXiv:2405.14059 (v3 19/7/2024). Es una **revisión**. | Abstract: casi todos los causets no manifold-like están muy suprimidos en la integral de caminos gravitatoria; «esto no demuestra del todo la emergencia de un continuo». | **V-A**. Coincide con lo atribuido en el acta. Matiz: la supresión la prueban trabajos citados (Loomis y Carlip 2017; Carlip, Carlip y Surya 2023), no esta revisión. |
| 5 | Glaser y Surya (2013) | Lisa Glaser, Sumati Surya, «Towards a definition of locality in a manifoldlike causal set», **Phys. Rev. D 88, 124026 (2013)**, arXiv:1309.3403. | Localidad definida por la abundancia de intervalos de m elementos; valores esperados analíticos en intervalos de Alexandrov de Minkowski d-dimensional; evidencia de regiones locales como condición *necesaria* de manifoldlikeness; nuevo estimador de dimensión. | **V-A**. Nótese: condición necesaria, no suficiente. |
| 6a | Bombelli, Lee, Meyer, Sorkin (1987) | «Space-time as a causal set», **Phys. Rev. Lett. 59, 521-524 (1987)**, doi 10.1103/PhysRevLett.59.521. | Propuesta fundacional. | **V-A** (Crossref y APS). |
| 6b | Surya (2019) | «The causal set approach to quantum gravity», **Living Rev. Relativ. 22, 5 (2019)**, doi 10.1007/s41114-019-0023-1, arXiv:1903.11544. | Revisión. | **V-P** (texto en ar5iv; usado para las fórmulas B1 a B3). |

Citas auxiliares comprobadas al verificar las fórmulas:
- Benincasa y Dowker, «Scalar Curvature of a Causal Set», Phys. Rev. Lett. 104, 181301 (2010), arXiv:1001.2725. **V-P** (texto del PDF v4).
- Brightwell y Gregory, «Structure of random discrete spacetime», Phys. Rev. Lett. 66, 260-263 (1991). **V-A** (solo metadatos; texto no leído).
- Bollobás y Winkler, «The longest chain among random points in Euclidean space», Proc. AMS 103, 347-353 (1988). **V-A**.

## B. Fórmulas

### B1. Número esperado de k-cadenas ⟨C_k⟩ en un diamante de Minkowski (Myrheim 1978; Meyer 1988)

Fuente: Surya (2019), ec. 19, atribuida a Meyer (1988). Proceso de Poisson de densidad ρ en un intervalo de Alexandrov de volumen V:

  ⟨C_k⟩ = χ_k (ρV)^k,  χ_k = (1/k) · (Γ(d+1)/2)^(k−1) · Γ(d/2) Γ(d) / [ Γ(kd/2) Γ((k+1)d/2) ].

Con **n fijo** (n puntos uniformes, que es lo que hace `sprinkle_diamond`), (ρV)^k se sustituye por el factorial descendente:

  ⟨C_k⟩_n = n(n−1)···(n−k+1) · χ_k.

(Derivación: C_k suma sobre k-subconjuntos ordenados por el tiempo; E[N(N−1)···(N−k+1)] = (ρV)^k para Poisson. No está escrita así en la fuente; la confirma el Monte Carlo de la sección C, que distingue las dos formas.)

Para d = 4, χ_k = 12^(k−1) · 6 / [ k (2k−1)! (2k+1)! ]:

| k | χ_k (d = 4) | ⟨C_k⟩, n fijo | ⟨C_k⟩, Poisson |
|---|---|---|---|
| 2 | 1/20 | n(n−1)/20 | (ρV)²/20 |
| 3 | 1/2100 | n(n−1)(n−2)/2100 | (ρV)³/2100 |
| 4 | 1/705600 | n(n−1)(n−2)(n−3)/705600 | (ρV)⁴/705600 |
| 5 | 1/582120000 | n(n−1)···(n−4)/582120000 | (ρV)⁵/582120000 |

Relación con la fracción de orden (Myrheim 1978, Surya ec. 14 y 18): r = 2R/[n(n−1)] con R = C_2; ⟨r⟩ = Γ(d+1)Γ(d/2)/[2Γ(3d/2)]. Para d = 4: 24·1/(2·120) = **1/10**. **Correcto**: el 10 % de los pares (no ordenados) está relacionado. Cuidado con la normalización: Meyer usa f₀ = ⟨R⟩/⟨n⟩² = 1/20 en d = 4, que es la mitad de r. El repo (`mm_fraction`, `ordering_fraction`) usa r y es coherente con 1/10.

Estado: **V-S** (Surya 2019, revisión; Myrheim 1978 es un preprint del CERN y Meyer 1988 una tesis del MIT, no consultados) + **confirmado numéricamente** (sección C).

### B2. Cadena más larga

Fuente: Surya (2019), atribuido a Brightwell y Gregory (1991); existencia del límite en el cubo por Bollobás y Winkler (1988).

  lim_{ρ→∞} ⟨L⟩ / (ρV)^(1/d) = m_d,

con cotas 2^(1−1/d)/Γ(1+1/d) ≤ m_d ≤ 2^(1−1/d) · e · Γ(1+d)^(1/d) / d.

Valores numéricos de las cotas (calculados): d = 2: [1.596, 2.718]; d = 3: [1.778, 2.614]; d = 4: **[1.856, 2.530]**. El «1.77 ≤ m_d ≤ 2.62» de Surya es la cota uniforme para d ≥ 3. m_2 = 2 es exacto (el diamante 2D es un cuadrado en coordenadas de cono de luz y L es la subsucesión creciente más larga; resultado de Logan-Shepp y Vershik-Kerov 1977, no consultado aquí). Para d > 2 **solo hay cotas**; no se conoce m_d en forma cerrada. Con n fijo, ρV = n.

Aviso práctico: con n = 100 y 200 el Monte Carlo da ⟨L⟩/n^(1/4) = 1.52 y 1.61, **por debajo de la cota inferior asintótica 1.856**. Las correcciones de tamaño finito son grandes; para m del orden de los sub-intervalos de R14-C hay que usar una referencia Monte Carlo del diamante plano, no la constante asintótica.

Estado: **V-S** (Brightwell y Gregory y Bollobás y Winkler solo en metadatos y abstract).

### B3. Acción y operador de Benincasa-Dowker en 4D

Fuente primaria: Benincasa y Dowker (2010), ec. 2, 3, 8, 9, 14; Surya (2019) ec. 27 a 37.

- Operador (densidad ρ = ℓ⁻⁴): Bφ(x) = (4/(√6 ℓ²)) [ −φ(x) + ( Σ_{L1} − 9 Σ_{L2} + 16 Σ_{L3} − 8 Σ_{L4} ) φ(y) ], con L_i = {y ≺ x : n(x,y) = i − 1} y n(x,y) = |I(x,y)| − 2 (elementos estrictamente entre). Coeficientes **(1, −9, 16, −8)** y prefactor **4/√6** verificados.
- Límite: ⟨Bφ⟩ → □φ en Minkowski; en curvado → □φ − ½Rφ. BD: «aplicado al campo constante φ = −2» da la curvatura escalar. Por tanto, con s(x) = 1 − N₀(x) + 9N₁(x) − 16N₂(x) + 8N₃(x) (N_k(x) = número de y ≺ x con k elementos entre): **R(x) ≈ (8/√6) ℓ⁻² s(x)**, es decir ⟨s(x)⟩ → (√6/8) ℓ² R(x).
- Acción (ε = 1): S/ℏ = (4/√6)(ℓ/ℓ_p)² [ N − N₁ + 9N₂ − 16N₃ + 8N₄ ] en la notación de BD, donde N_i = número de intervalos inclusivos de (i+1) elementos (N₁ = links), y ℓ_p² = 8πGℏ. Comprobación de normalización: (4/√6)(ℓ²/ℓ_p²) · ρ∫(√6/8)ℓ²R dV = (1/(2ℓ_p²))∫R√−g, que es Einstein-Hilbert. **Hasta términos de frontera**: en un diamante plano ⟨S⟩ ≠ 0. BD escriben la acción de su ec. 14 «salvo factores de orden uno»; el prefactor 4/√6 y (ℓ/ℓ_p)² aparece explícito en Surya ec. 34 y 35.
- Versión suavizada: f(n,ε) = (1−ε)^n [ 1 − 9εn/(1−ε) + 8ε² n(n−1)/(1−ε)² − 4ε³ n(n−1)(n−2)/(3(1−ε)³) ], ε = (ℓ/ℓ_k)⁴, B_k φ = (4/(√6 ℓ_k²)) [ −φ + ε Σ_{y≺x} f(n(x,y),ε) φ(y) ]. f(n,1) da (1, −9, 16, −8, 0, ...). La acción de la familia es Σ_x B_k(−1) por ℏℓ²(ℓ²/ℓ_p²), luego S_ε/ℏ = (4/√6)(ℓ/ℓ_p)² **√ε** [ N − ε Σ_n f(n,ε) N_n ].

Comparación con el repo:
- `omega/causality/estimators.py::bd_action_4d`: N − N0 + 9N1 − 16N2 + 8N3 con N_k = pares con k elementos entre (`interval_abundances`), es decir N_k(repo) = N_{k+1}(BD). **Correcto**, con el prefactor (4/√6)(ℓ/ℓ_p)² absorbido, como dice el docstring.
- `smearing_weight_4d`: coincide término a término con la ec. 9 de BD, con n = elementos estrictamente entre. **Correcto.**
- `bd_action_4d_smeared`: devuelve N − ε Σ f. **Falta el factor √ε** de la acción de BD. Es irrelevante si ε es fijo, pero cambia la escala relativa si se comparan acciones o β entre distintos ε. (La versión 2D, `bd_action_2d_smeared`, sí incluye su factor ε, que es el análogo correcto en 2D, donde ℓ²/ℓ_k² = ε.) Es una inconsistencia de normalización, no un error de signo ni de coeficientes.
- `omega/causality/torus.py::retarded_action_4d`: s(x) = 1 − Σ_{y≺x} coef[n] con coef = (1, −9, 16, −8), o 1 − εΣf. Es s(x) = (√6 ℓ²/4) B_k(−1)(x). **Correcto**; para pasar a curvatura, R ≈ (8/√6) ℓ_k⁻² s(x) (con ε = 1, ℓ_k = ℓ).

Estado: **V-P** (BD 2010, texto completo; la ec. 4 del valor medio no se leyó limpia en el PDF y no se usa aquí).

### B4. Observables de Eichhorn et al. (definiciones exactas, sección III)

- Matriz de links L_ij = 1 si e_i ≺* e_j (relación sin intermedios; diagrama de Hasse). Ec. 13: d_out,j = Σ_i L_ji, d_in,j = Σ_i L_ij.
- Grado no dirigido, ec. 14: **d_j = ½ Σ_i (L_ji + L_ij)**. Distribución, ec. 15: P_j = (1/n) Σ_i δ_{d_i, j}, Σ_j P_j = 1.
- Laplaciano, ec. 17 a 19: L_S = L + Lᵀ; Δ = 𝟙 − 𝒟^(−1/2) L_S 𝒟^(−1/2) (Laplaciano simétricamente normalizado), con 𝒟 diagonal de grados y 𝒟^(−1/2) como pseudoinversa de Moore-Penrose (vértices aislados ↦ 0). Se usan todos los autovalores ordenados y el par reducido (λ₂, λ_n); los extremos se calculan con Lanczos reiniciado implícitamente.
- **Ambigüedad a resolver antes de implementar**: la ec. 17 define el grado «a partir de la matriz de adyacencia A_S» (suma de filas, d_in + d_out), y la ec. 19 remite a la ec. 14, que lleva un ½. Con el ½, Δ = 𝟙 − 2·(normalizado estándar) y el espectro pasaría de [0, 2] a [−1, 3]; el texto afirma que λ_max = 2 para órdenes por capas (bipartitos), que solo cuadra con la convención estándar (sin ½). Recomendación: usar 𝒟 = diag(suma de filas de L_S) y anotarlo, o mirar el código citado (QuantumGrav, CausalSetZoology, CausalSets.jl). Para la distribución de grado, el ½ solo reescala el eje (y da semienteros).
- Abundancia de intervalos, ec. 16: S_m = [ n / Σ_ij C_ij ] · Σ_{pares} δ_{|I[e_i,e_j]|, m}, con I inclusivo (un link tiene |I| = 2), solo pares relacionados. Con esa normalización Σ_m S_m = n; la distribución de probabilidad en m/n es S_m/n. La relación con el repo: S_m/n = N_{m−2}(repo) / R.
- Tamaños: n = 2048 en todas las clases; manifold-like en d = 2 con frontera de caja.

Estado: **V-P** (HTML del preprint), con la ambigüedad del ½ señalada.

## C. Comprobación numérica de ⟨C_k⟩ (d = 4)

`chains.py`: SeedSequence(830001) con `spawn(2)`, `sprinkle_diamond` y `causal_matrix` del repo, cadenas por programación dinámica (v₁ = 1, v_k = Cᵀ v_{k−1}, C_k = Σ v_k). 4000 repeticiones con n = 100 y 2000 con n = 200.

| n | k | teoría n fijo | media MC | error estándar | z |
|---|---|---|---|---|---|
| 100 | 2 | 495.00 | 492.87 | 1.56 | −1.36 |
| 100 | 3 | 462.00 | 458.57 | 3.69 | −0.93 |
| 100 | 4 | 133.38 | 132.78 | 2.23 | −0.27 |
| 100 | 5 | 15.52 | 15.78 | 0.60 | +0.43 |
| 200 | 2 | 1990.00 | 1994.46 | 6.29 | +0.71 |
| 200 | 3 | 3752.57 | 3770.52 | 29.92 | +0.60 |
| 200 | 4 | 2200.17 | 2216.56 | 34.23 | +0.48 |
| 200 | 5 | 522.71 | 527.77 | 15.50 | +0.33 |

- Todas las |z| < 1.4: la fórmula de n fijo se confirma. La forma de Poisson con ρV = n queda descartada para n fijo pequeño (n = 100, k = 2: 500 frente a 492.9 ± 1.6, z ≈ −4.6; k = 5: 17.18 frente a 15.78 ± 0.60, z ≈ −2.3).
- Fracción de orden r: 0.09957 ± 0.00032 (n = 100) y 0.10022 ± 0.00032 (n = 200); teoría 1/10.
- Cadena más larga: ⟨L⟩ = 4.81 (n = 100) y 6.05 (n = 200); ⟨L⟩/n^(1/4) = 1.52 y 1.61, por debajo de la cota asintótica 1.856 (ver B2).

## D. Lo que NO se pudo verificar

1. Myrheim (1978, preprint CERN) y Meyer (1988, tesis MIT): no consultados; la fórmula de ⟨C_k⟩ se verifica a través de Surya (2019) y del Monte Carlo.
2. Brightwell y Gregory (1991): texto del PRL no leído; las cotas de m_d vienen de Surya (2019). No se ha comprobado si la cota y el límite en el diamante para d > 2 son de Brightwell y Gregory o de Bollobás y Brightwell («Box-spaces and random partial orders»); una fuente secundaria atribuye la asintótica a Bollobás y Brightwell.
3. m_2 = 2 (Logan-Shepp, Vershik-Kerov 1977): conocido, no consultado en fuente primaria.
4. Boguñá y Krioukov, Carlip, Glaser y Surya, BLMS: solo abstract y metadatos; el contenido técnico (definición del d'Alembertiano local, fórmulas de abundancias de Glaser y Surya) no se ha leído.
5. Eichhorn et al.: la convención efectiva del ½ en el grado (ec. 14 frente a ec. 17) no se ha comprobado en su código.
6. La ec. 4 de BD (valor medio con e^(−ξ)(1 − 9ξ + 8ξ² − 4ξ³/3)) no se leyó limpia en el PDF.
7. Carlip: la supresión la demuestran trabajos citados por la revisión, que no se han leído.

Fuentes: arxiv.org/abs/2503.01719 (y html v4); arxiv.org/abs/2605.27514 (y html); arxiv.org/abs/2506.18745; arxiv.org/abs/2405.14059; arxiv.org/abs/1309.3403; arxiv.org/abs/1903.11544 (ar5iv); arxiv.org/pdf/1001.2725; api.crossref.org para los DOI 10.1007/s10714-024-03281-1, 10.1007/s41114-019-0023-1, 10.1103/np89-qzjp, 10.1103/PhysRevD.88.124026, 10.1103/PhysRevLett.59.521, 10.1103/PhysRevLett.66.260, 10.1103/PhysRevLett.104.181301; ams.org/proc/1988-103-02/S0002-9939-1988-0943043-6.
