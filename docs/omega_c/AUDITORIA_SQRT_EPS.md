# Ω-C: auditoría de √ε y de la acción BD incremental (Rev. 17 §8, pasos 1 y 2)

- **Fecha:** 2026-10-19. **Mandato:** `OMEGA_CONSEJO_REV17_ES0.md` §8, pasos 1 y 2. Antecedente: `consejo_rev15/VERIFICACION_R15_2.md` B3 y `OMEGA_CONSEJO_REV15_SESION2.md` (hallazgo lateral).
- **Alcance:** solo lectura del código y de los documentos. No se modificó ningún módulo ni documento existente. No se hizo commit.
- **Archivos nuevos:** este informe y `tests/test_bd_incremental.py`. Chequeos numéricos baratos en el scratchpad de la sesión (`sqrt_eps/check.py`), resumidos en §1.5 y §3.
- **Resultado corto:** el factor √ε existe en la literatura y falta en el código 4D. No cambia ningún veredicto registrado: en todas las herramientas el factor se absorbe en el acoplamiento o se cancela en la normalización. La acción incremental de `mcmc_general` coincide exactamente con la completa.

Marcas: **V-P** = verificado en la fuente primaria (texto completo y, donde el texto extraído era dudoso, imagen de la página). **V-T** = deducido de fórmulas verificadas. **NV** = no verificado.

## 1. Fórmula verificada con fuentes

### 1.1 Fuentes leídas

| Clave | Referencia | Lectura |
|---|---|---|
| BD10 | Benincasa y Dowker, «The Scalar Curvature of a Causal Set», arXiv:1001.2725v4 | texto completo |
| DG13 | Dowker y Glaser, «Causal set d'Alembertians for various dimensions», arXiv:1305.2588v2 | texto completo; páginas 8 y 14 vistas como imagen |
| BBD16 | Belenchia, Benincasa y Dowker, arXiv:1510.04656v2 | ec. 1 a 3 |
| S19 | Surya, «The causal set approach to quantum gravity», Living Rev. Rel., arXiv:1903.11544 | §4.5 (ec. 27 a 37) y §6 (ec. 90) |

### 1.2 Operador 4D sin suavizar (V-P)

- BD10 ec. 2 y 3; BBD16 ec. 3: Bφ(x) = (4/(√6 ℓ²)) [ −φ(x) + ( Σ_{L1} − 9 Σ_{L2} + 16 Σ_{L3} − 8 Σ_{L4} ) φ(y) ].
- Capas: L_i = {y ≺ x : n(x, y) = i − 1}, con n(x, y) = |I(x, y)| − 2. Es decir, n cuenta los elementos **estrictamente entre** x e y (BD10 tras ec. 3; DG13 p. 1; BBD16 antes de ec. 1).
- DG13 ec. 24 y Tabla 2: la misma forma general, B^(d) = (1/ℓ²)(α_d φ + β_d Σ_i C_i Σ_{L_i} φ). Para d = 4: α₄ = −2/Γ(3/2) · (π/6)^{1/2} = −4/√6 y β₄ = 4/√6 (V-T, cálculo directo). C^(4) = (1, −9, 16, −8).

### 1.3 Operador 4D suavizado y definición de ε (V-P)

- BD10 ec. 8: B_k φ(x) = (4/(√6 ℓ_k²)) [ −φ(x) + ε Σ_{y≺x} f(n(x, y), ε) φ(y) ].
- BD10, junto a ec. 8: **ε = (ℓ/ℓ_k)⁴**, con ℓ_k ≥ ℓ la escala de no localidad.
- BD10 ec. 9: f(n, ε) = (1 − ε)^n [ 1 − 9εn/(1 − ε) + 8ε² n(n − 1)/(1 − ε)² − 4ε³ n(n − 1)(n − 2)/(3(1 − ε)³) ].
- BD10 tras ec. 9: «B_k reduce a B cuando ε = 1».
- DG13 ec. 25 (general d): B_ε^(d) φ = (ε^{2/d}/ℓ²) ( α_d φ + β_d ε Σ_{y≺x} f_d(n, ε) φ(y) ), con **ε := (ℓ/ξ)^d** (ξ = escala de no localidad). Como ε^{2/d}/ℓ² = 1/ξ², coincide con BD10 ec. 8 en d = 4.
- DG13 ec. 26: f_d(n, ε) = (1 − ε)^n Σ_i C_i^(d) binom(n, i − 1) (ε/(1 − ε))^{i−1}. Para d = 4 da exactamente BD10 ec. 9 (16·n(n−1)/2 = 8n(n−1); 8·n(n−1)(n−2)/6 = (4/3)n(n−1)(n−2)). V-T.
- S19 ec. 36 y 37: misma f, con ε = ρ_κ/ρ_c = (ℓ_c/ℓ_κ)^d (glosario de S19). **Atención:** S19 ec. 36 escribe B_κ φ = (4/√6)[−φ + ε Σ f φ] **sin el factor 1/ℓ_κ²** y con n = |I(e, e′)| (inclusivo). Las dos cosas difieren de BD10 y DG13. La propuesta del repo (§19) cita a S19 para f₄; es la explicación más probable de la omisión de √ε (NV como causa; la diferencia de textos sí es V-P).
- DG13 nota 2 (p. 8): avisa de un error de convención de n en la fórmula del operador no local de BD10. Con n = elementos estrictamente entre, el límite ε → 1 da los coeficientes correctos (ver §1.5), que es lo que fija la convención. A qué versión de BD10 se refiere la nota: NV.

### 1.4 Acción 4D suavizada (V-P para la receta y para DG13 ec. 42)

- BD10, párrafo antes de ec. 13: la familia de acciones se obtiene sumando B_k(−1) sobre los elementos, por ℏℓ² (unidades), por ℓ²/ℓ_p² (número de orden uno en 4D), con ℓ_p² = 8πGℏ.
- Aplicado a BD10 ec. 8 (V-T):
  S_ε/ℏ = (ℓ/ℓ_p)² ℓ² Σ_x B_k(−1)(x) = (4/√6)(ℓ/ℓ_p)² (ℓ/ℓ_k)² [ N − ε Σ_{y≺x} f(n, ε) ],
  y (ℓ/ℓ_k)² = **√ε**.
- DG13 ec. 42 (vista en la imagen de la p. 14): S_ε^(d)/ℏ = ζ_d [ ε^{2/d} N + ε^{(d+2)/d} (β_d/α_d) Σ_x Σ_{y≺x} f_d(n(x, y), ε) ], con ζ_d = −α_d (ℓ/ℓ_p)^{d−2} (DG13 tras ec. 41).
- Para d = 4: ζ₄ = (4/√6)(ℓ/ℓ_p)², β₄/α₄ = −1, ε^{2/4} = √ε, ε^{6/4} = ε^{3/2}. Por tanto (V-T):

  **S_ε^(4)/ℏ = (4/√6)(ℓ/ℓ_p)² √ε [ N − ε Σ_{pares y≺x} f₄(n(x, y), ε) ].**

  Equivalentemente, el prefactor es (4/√6) ℓ⁴/(ℓ_p² ℓ_k²).
- Sin suavizar (ε = 1): BD10 ec. 14, S/ℏ = N − N₁ + 9N₂ − 16N₃ + 8N₄ «salvo factores de orden uno», con N_i = intervalos inclusivos de i + 1 elementos. DG13 ec. 41 y S19 ec. 34 dan el prefactor explícito (4/√6)(ℓ/ℓ_p)². En notación del repo (N_k = pares con k elementos entre): N − N₀ + 9N₁ − 16N₂ + 8N₃.
- **Signo:** positivo, S = +Σ_x (algo proporcional a B(−1)), y el aporte de un link es −1 (ε = 1) o −ε^{3/2} (en unidades del prefactor). V-T.
- **Unidades:** ℏ y ℓ_p entran solo en el prefactor global (4/√6)(ℓ/ℓ_p)². La parte dependiente de ε es √ε, adimensional. V-P (BD10 receta; DG13 ec. 41 y 42).
- **Peso estadístico:** S19 ec. 90 usa e^{iS/ℏ} y dice que la continuación β → iβ permite MCMC; eso da e^{−βS} con el mismo signo de S. V-P en S19 solo; las fuentes originales (Surya 2012, Glaser et al. 2018) NV.

### 1.5 Límite ε → 1 (V-T y numérico)

- f₄(n, 1) vale (1, −9, 16, −8, 0, 0, …) para n = 0, 1, 2, 3, 4, …, porque (1 − ε)^n (ε/(1 − ε))^{i−1} → 1 solo si n = i − 1.
- Numérico: f₄(n, 1 − 10⁻⁶) = (1, −9, 16, −7.9999, 0, 0).
- Con √1 = 1, S_ε → (4/√6)(ℓ/ℓ_p)² [N − N₀ + 9N₁ − 16N₂ + 8N₃]: **prefactor 4/√6 y coeficientes (1, −9, 16, −8)**. La fórmula corregida es continua en ε = 1. Numérico (sprinkling 4D, N = 200): √ε·S_repo(ε = 1 − 10⁻⁷) = 281.9998 frente a bd_action_4d = 282.
- Numérico: DG13 ec. 42 sin ζ₄ coincide con √ε · `bd_action_4d_smeared` a 6 decimales para ε ∈ {0.5, 0.25, 0.05}.

### 1.6 Caso 2D (V-P y V-T)

- BD10 ec. 10 (media en 2D): prefactor 2/ℓ_k², pesos e^{−ξ}(1 − 2ξ + ξ²/2), es decir C^(2) = (1, −2, 1). V-P.
- DG13 Tabla 2: α₂ = −2, β₂ = 4. DG13 ec. 25 con d = 2: B_ε = (ε/ℓ²)(−2φ + 4ε Σ f₂ φ) = (2/ℓ_k²)(−φ + 2ε Σ f₂ φ), con **ε = (ℓ/ℓ_k)²**. V-T.
- DG13 ec. 42 con d = 2: ζ₂ = 2 (ℓ/ℓ_p)⁰ = 2, ε^{2/2} = ε, ε^{4/2} = ε², β₂/α₂ = −2:
  **S_ε^(2)/ℏ = 2ε [ N − 2ε Σ f₂(n, ε) ]**, con f₂(n, ε) = (1 − ε)^n [1 − 2εn/(1 − ε) + ε² n(n − 1)/(2(1 − ε)²)]. V-T.
- Sin suavizar: BD10 ec. 13, N − 2N₁ + 4N₂ − 2N₃ (notación BD), y DG13 ec. 41 con ζ₂ = 2.
- En 2D el factor análogo a √ε es ε¹ = (ℓ/ℓ_k)². El repo lo incluye (ver §2).

## 2. Diferencia con el código

| Función | Código actual | Fórmula de la fuente (prefactor global quitado igual que en ε = 1) | Veredicto |
|---|---|---|---|
| `estimators.bd_action_4d` | N − N₀ + 9N₁ − 16N₂ + 8N₃ | igual (sin (4/√6)(ℓ/ℓ_p)²) | correcto |
| `estimators.smearing_weight_4d` | f₄ de BD10 ec. 9 con n = elementos entre | BD10 ec. 9 = DG13 ec. 26 | correcto |
| `estimators.bd_action_4d_smeared` | N − ε Σ f₄ | **√ε** (N − ε Σ f₄) | **falta √ε**; nada más difiere |
| `torus.retarded_action_4d` | s(x) = 1 − ε Σ_{y≺x} f₄ (por elemento) | s(x) = (√6 ℓ_k²/4) B_k(−1)(x); acción = √ε Σ_x s(x) | correcto como densidad local; su suma **no** es la acción hasta multiplicar por √ε |
| `mcmc_general.action4` | N + Σ COEF4[k] sobre relaciones, COEF4 = (−1, 9, −16, 8, 0) | igual a `bd_action_4d` | correcto (probado en §4) |
| `estimators.bd_action_2d` | N − 2N₀ + 4N₁ − 2N₂ | DG13 ec. 41 / 2 | correcto |
| `estimators.bd_action_2d_smeared` | ε (N − 2ε Σ f₂) | DG13 ec. 42 (d = 2) / 2 | **correcto**: lleva el factor ε, análogo 2D de √ε |
| `tools/c4_l4.py::s4_from_between` | copia de N − ε Σ f₄ | √ε (N − ε Σ f₄) | falta √ε (misma omisión) |
| `reference4d.expected_action_4d` | N − ε C(N,2)·E[f₄] | √ε por lo mismo | falta √ε, pero es coherente con el estimador que valida |

Detalle de lo que **no** difiere:
- Coeficientes de f₄: idénticos a BD10 ec. 9 (comprobado término a término y con el límite ε → 1).
- Argumento de f₄: el repo usa `between_counts` (elementos estrictamente entre), que es la convención de BD10 y DG13. La forma inclusiva de S19 ec. 36 sería incorrecta.
- Signo: el mismo que la fuente (un link aporta −ε en el repo, −ε^{3/2} con √ε).
- Término N: el mismo.
- La rama ε = 1 de `bd_action_4d_smeared` delega en `bd_action_4d`: con √1 = 1 coincide con la fórmula corregida. El repo no mezcla convenciones en ε = 1.

Conclusión de §2: el repo implementa S_repo(ε) = S_ε/√ε. Es una normalización distinta por cada ε, no un error de signo ni de coeficientes. Con ε fijo equivale a reescalar el acoplamiento: β_fuente = β_repo/√ε (o λ_fuente = λ_repo/√ε en la rama lorentziana).

## 3. Impacto por herramienta

Búsqueda: `grep` de `bd_action_4d_smeared`, `smearing_weight_4d`, `retarded_action_4d`, `bd_action_2d_smeared` en `tools/`, `omega/` y `tests/`. Solo dos herramientas usan la acción 4D suavizada: `c4_l4.py` y `l4b.py`. Las otras tres de la lista usan la 2D (correcta). `estimators.bd_action_4d_smeared` no lo usa ninguna herramienta; solo `tests/test_causality.py` (línea 229).

### 3.1 `tools/c4_l4.py` (C4 y L4; `OMEGA_2_C4_L4_RESULTADOS.md`)

Código:
- `s4_from_between`: `return float(n - eps * float(np.sum(smearing_weight_4d(vals.astype(np.float64), eps))))` (ε < 1); rama ε = 1 entera sin suavizar.
- C4: `ref, ref_se = expected_action_4d(n, e, v)` y `"ok": abs(mean - ref) <= tol`.
- L4: `"D4": float(np.std(v, ddof=1)) / sched[k]`, esquemas fijos (1, 0.5, 0.25, 0.1, 0.05) y ε = K/N con K ∈ {16, 64}.
- Veredicto: `D4_2048 <= D_MAX` (3/(2π)) y `eta + 2*se <= 0`.

Análisis:
- **C4.** Estimador y referencia usan la misma convención. Con √ε, media, SE, referencia y su SE se multiplican todos por √ε, y la tolerancia también. Las seis celdas siguen igual. Los números de la tabla de §1 del documento están en unidades del repo (por ejemplo, N = 128, ε = 0.25: 68.8 pasaría a 34.4, y la referencia 70.0 a 35.0). **Reescalado, sin efecto en el veredicto.**
- **L4.** Aquí sí se comparan esquemas de ε distintos, y los esquemas K/N cambian ε con N. Pero D₄ no es sd(S) a secas. La propuesta §19.2 lo define como la dispersión de fase en unidades del aporte de un link: «Per link, S₄,ε aporta −ε … la fase de la geometría es (λε)·(S₄,ε/ε)». Con la acción corregida, un link aporta −ε^{3/2}, la ventana es λ′ε^{3/2} ∈ W₄ y la fase de la geometría es (λ′ε^{3/2})·(S_ε/ε^{3/2}). Como S_ε = √ε·S_repo, se tiene sd(S_ε)/ε^{3/2} = sd(S_repo)/ε: **D₄ es exactamente el mismo** para cada ε y cada N. Los ajustes η, el umbral 3/(2π) y el veredicto «INCOMPATIBLE (esquemas probados)» no cambian.
- Chequeo numérico (40 sprinklings 4D, N = 128): D₄ repo = D₄ corregido = 51.357 (ε = 0.25), 56.379 (ε = 0.05), 43.162 (ε = 16/N). Solo una mezcla incoherente (acción con √ε pero normalizada por ε) cambiaría: 25.7, 12.6 y 15.3, y en los esquemas K/N la pendiente η bajaría 0.5. Ninguna herramienta hace esa mezcla.
- Afirmaciones del documento que dependen de esto: la tabla de L4 (D₄ y η) y la lectura «la dispersión … está entre 600 y 2 800 en N = 2048». Siguen válidas tal como están, porque D₄ está en unidades del aporte por link.

**Veredicto c4_l4: reescalado de λ por 1/√ε solamente (λ_fuente = λ_repo/√ε). C4 y L4 sin efecto.**

### 3.2 `tools/l4b.py` (C4b y L4b; `OMEGA_2_L4B_RESULTADOS.md`)

Código:
- `samples[e].append(float(retarded_action_4d(c, e).sum()))` y `ref, ref_se = expected_action_4d(n, e, v)` (C4b).
- `"D_tot": float(np.std(tot[k], ddof=1)) / sch[k]`, `"D_int": float(np.std(intr[k], ddof=1)) / sch[k]` (L4b).

Análisis: igual que L4. `retarded_action_4d(c, e).sum()` es S_repo (falta √ε), y la referencia también. D_int y D_tot se dividen por ε, el aporte por link en las mismas unidades. Con la acción corregida ambos factores llevan √ε y se cancelan. η_int, η_tot, el diagnóstico «η_int compatible con 0.5» y el veredicto «INCOMPATIBLE» no cambian. Los valores de C4b (39.05 frente a 40.17; 358.06 frente a 358.17) se multiplicarían por √ε juntos.

**Veredicto l4b: reescalado de λ por 1/√ε solamente. Sin efecto en C4b ni en L4b.**

### 3.3 `tools/bloque_a.py` (`OMEGA_2_BLOQUE_A_RESULTADOS.md`)

Solo 2D: `bd_action_2d_smeared(c, eps)` y `smeared_from_between`: `eps * (n - 2 * eps * sum f)`, con `lam = x / eps**2` (aporte por link −2ε²). Es la fórmula correcta de DG13 ec. 42 con d = 2 (salvo ζ₂ = 2, constante). **Sin efecto.**

### 3.4 `tools/bloque_as.py` (`OMEGA_2_BLOQUE_AS_RESULTADOS.md`)

Solo 2D: `action_from_types`: `eps * (n - 2 * eps * (links + f1 * n1 + f2 * n2))`, validada contra `bd_action_2d_smeared`; fase con `lam = x / eps**2`. Fórmula 2D correcta. **Sin efecto.**

### 3.5 `tools/c1_l3ce.py` (C1 y L-3c-E)

C1 compara `bd_action_2d_smeared(c, eps)` con `expected_smeared_action_2order(n, eps)`: 2D, misma convención en ambos lados. L-3c-E usa `omega/causality/mcmc.py`, que muestrea con `bd_action_2d` sin suavizar. **Sin efecto.**

### 3.6 Otras piezas

- `tools/p4.py` y `OMEGA_2_P4_RESULTADOS.md`: usan `mcmc_general` con la acción 4D sin suavizar (ε = 1, √1 = 1). **Sin efecto.**
- L-ε (`OMEGA_2_PROPUESTA_CAUSAL.md` §18): es 2D, con factor ε correcto. **Sin efecto.** Esto resuelve la duda de `OMEGA_CONSEJO_REV15_SESION2.md` («puede afectar a las comparaciones entre esquemas de ε de L4 y L-ε»): no afecta, porque L-ε no usa la acción 4D y L4 normaliza por el aporte de un link.
- `tests/test_causality.py:229`: `retarded_action_4d(c, 0.25).sum() == bd_action_4d_smeared(c, 0.25)`. Coherente hoy (las dos sin √ε). Si se aplica la propuesta de §5 hay que actualizarla.
- `OMEGA_R14C_V2_PREREGISTRO.md` fila 8: usa la versión sin suavizar. Sin efecto.

### 3.7 Resumen

| Herramienta | Acción | ¿Compara ε distintos? | Veredicto |
|---|---|---|---|
| c4_l4.py | 4D suavizada sin √ε | sí, pero D₄ normalizado por aporte de link | reescalado de λ por 1/√ε; veredictos intactos |
| l4b.py | 4D suavizada sin √ε | sí, igual que L4 | reescalado de λ por 1/√ε; veredictos intactos |
| bloque_a.py | 2D suavizada (correcta) | sí | sin efecto |
| bloque_as.py | 2D suavizada (correcta) | sí | sin efecto |
| c1_l3ce.py | 2D suavizada (correcta) y 2D sin suavizar | no mezcla | sin efecto |
| p4.py | 4D sin suavizar | no | sin efecto |

Riesgo futuro: la rama euclídea 4D suavizada (Ω-C) sí compararía β entre ε distintos, por ejemplo una malla (β, ε) o umbrales β_c(ε). Ahí el factor importa: β_c(ε) del repo sería √ε veces el de la fuente. Por eso conviene congelar la fórmula antes del preregistro (§5).

## 4. Incremental frente a completa

### 4.1 Qué se probó

`tests/test_bd_incremental.py` (57 pruebas, unos 4 s):

- **Movimiento a movimiento, con el código real de `sweep`.** Un generador de prueba fuerza que la primera propuesta sea el par (i, j) y se acepte; con β = NaN todas las demás se rechazan (−β·ΔS ≥ 0 y u < e^{NaN} son falsos). Se comprueba:
  - aceptado ⇔ el movimiento es válido según una definición independiente (quitar si es link; agregar si el cierre transitivo no agrega nada más);
  - el orden resultante es el esperado y sigue siendo un orden parcial;
  - **st.s − S(antes) = bd_action_4d(después) − bd_action_4d(antes), exacto en enteros**;
  - la matriz B mantenida es igual a C·C recalculada **en todos los pares** (también los no relacionados, donde vale 0).
- **Exhaustivo para N ≤ 5:** todos los órdenes etiquetados naturalmente (2, 7, 40 y 357 para N = 2, 3, 4, 5; el 40 coincide con P4-0) y todos los pares i < j: 3 833 propuestas, 2 714 válidas.
- **Aleatorio para N ∈ {8, 12, 20, 30}:** sprinklings 2D y 4D en diamante, 2-órdenes aleatorios, percolación transitiva (p = 0.1 y 0.3), todos reetiquetados por una extensión lineal, y órdenes transitivos generados en forma triangular. 24 órdenes; en cada uno hasta 12 links, 12 pares agregables y 6 inválidos.
- **Tres formas de la acción completa:** `mcmc_general.action4(c, B)` = `estimators.bd_action_4d(c)` = `torus.retarded_action_4d(c, 1).sum()` en los 24 órdenes. COEF4 = (−1, 9, −16, 8, 0).
- **Cadenas reales:** `sweep` con β ∈ {−0.5, −0.1, 0, 0.05, 0.5} (incluye negativos, por A-22), 4 barridos desde un sprinkling 4D (N = 16), un 2-orden (N = 12) y la anticadena (N = 10). Tras cada barrido: orden parcial triangular superior, st.s == bd_action_4d(st.c) exacto y B == C·C.
- **Control de sensibilidad (scratchpad):** al cambiar en `_c` el coeficiente de k = 3 de 8 a 7, la prueba exhaustiva falla en N = 5, como debe (k = 3 exige 5 elementos).

### 4.2 Resultado

- Todas las pruebas pasan. `ruff check` y `ruff format --check` limpios en el archivo nuevo.
- **La actualización incremental de `mcmc_general` es exacta** para la acción 4D sin suavizar, en todos los movimientos de agregar y quitar probados, y mantiene B correcta en todos los pares.
- `mypy` del entorno no tiene stubs de numpy ni pytest; da los mismos tipos de error en los tests existentes. No se usó como criterio.

### 4.3 Lo que no existe

**No hay muestreador incremental de la acción suavizada.** `mcmc_general` solo implementa COEF4 (ε = 1). `bd_action_4d_smeared` y `retarded_action_4d` son solo recálculos completos, O(N³). Una rama euclídea 4D suavizada necesitará un ΔS_ε nuevo: con la fórmula de §5, para agregar o quitar (i, j),
ΔS_ε = √ε · (−ε) · [ ±f₄(0, ε) + Σ_{s ≻ j} (f₄(B_is ± 1) − f₄(B_is)) + Σ_{p ≺ i} (f₄(B_pj ± 1) − f₄(B_pj)) ],
con los signos de `sweep` (+ al agregar). Es una propuesta, no verificada ni implementada; ya no es entera, así que su prueba debe usar tolerancia relativa y además comparar con `bd_action_4d_smeared` corregida. El arnés de `tests/test_bd_incremental.py` sirve tal cual.

## 5. Propuesta de fórmula y código congelados (NO aplicada)

### 5.1 Fórmula

Para ε ∈ (0, 1], adimensional, quitando el mismo prefactor que en `bd_action_4d`, es decir (4/√6)(ℓ/ℓ_p)²:

**S_ε^(4)(C) = √ε · [ N − ε Σ_{x ≻ y} f₄(n(x, y), ε) ]**, n(x, y) = número de elementos estrictamente entre,

f₄(n, ε) = (1 − ε)^n [ 1 − 9εn/(1 − ε) + 8ε² n(n − 1)/(1 − ε)² − 4ε³ n(n − 1)(n − 2)/(3(1 − ε)³) ], ε = (ℓ/ℓ_k)⁴.

- Fuentes: BD10 ec. 8 y 9 más la receta antes de ec. 13; DG13 ec. 25, 26, 41 y 42 con d = 4.
- S_1 = N − N₀ + 9N₁ − 16N₂ + 8N₃ = `bd_action_4d`. Continua en ε → 1.
- Peso euclídeo: e^{−β S_ε}, con β = (4/√6)(ℓ/ℓ_p)² en unidades de ℏ (positivo en la continuación de S19 ec. 90). Relación con lo ya medido: β_nuevo = β_viejo/√ε a ε fijo.
- 2D, sin cambios: S_ε^(2) = ε [N − 2ε Σ f₂] (DG13 ec. 42 / ζ₂).

### 5.2 Código propuesto

```python
def bd_action_4d_smeared(c: BoolMatrix, eps: float) -> float:
    """S4_eps = sqrt(eps) (N - eps sum_{pares} f4(n, eps)) (BD 2010 ec. 8; Dowker-Glaser 2013 ec. 42, d = 4).

    Prefactor (4/sqrt(6)) (l/l_p)^2 absorbido como en bd_action_4d; eps = (l/l_k)^4; eps = 1 da bd_action_4d.
    """
    if eps == 1.0:
        return float(bd_action_4d(c))
    if not 0 < eps < 1:
        raise ValueError("eps debe estar en (0, 1]")
    b = between_counts(c)
    vals = b[b >= 0].astype(np.float64)
    return float(math.sqrt(eps) * (c.shape[0] - eps * float(np.sum(smearing_weight_4d(vals, eps)))))
```

Cambios asociados (todos NO aplicados):
1. `torus.retarded_action_4d`: dejarla como densidad local s(x) y documentar en el docstring que la acción es √ε · Σ_x s(x) y que R(x) ≈ (8/√6) ℓ_k⁻² s(x).
2. `tests/test_causality.py:229`: pasar a `math.sqrt(0.25) * retarded_action_4d(c, 0.25).sum() == approx(bd_action_4d_smeared(c, 0.25))`.
3. Pruebas nuevas congeladas: (a) identidad con DG13 ec. 42 (d = 4) en un sprinkling; (b) continuidad en ε = 1 − 10⁻⁷ frente a `bd_action_4d`; (c) escalado S_ε/S_repo = √ε.
4. `tools/c4_l4.py`, `tools/l4b.py` y `reference4d.expected_action_4d`: no tocarlos. Son herramientas de resultados ya registrados y, por §3, sus veredictos no dependen del factor. Añadir solo una nota en sus documentos de resultados: «unidades del repo sin √ε; D₄ y D_int son invariantes».
5. Cualquier muestreador euclídeo 4D suavizado de Ω-C debe usar la función corregida y declarar el acoplamiento como β de la fuente (no β·√ε).

Estado: **propuesta para firma**. Nada de §5 está implementado.
