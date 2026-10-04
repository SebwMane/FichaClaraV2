# SPEC_EXTRACT.md — Extracción estructurada de OMEGA_MASTER.md

Documento de referencia técnica fiel. Extrae únicamente lo definido en OMEGA_MASTER.md sin adiciones. Secciones numeradas para trazabilidad.

---

## 1. Definiciones y fórmulas exactas

### Representación fundamental

**Ω = (V, E, W)**
- V: grados de libertad fundamentales (nodos)
- E: relaciones entre grados de libertad
- W: matriz de pesos asociados a relaciones

**Propiedades de la matriz W:**
- Simétrica: W_ij = W_ji
- Rango de valores: 0 ≤ W_ij ≤ 1
- Diagonal nula: W_ii = 0
- Inicialización: W_ij ~ U(0,1) con simetrización

### Funcional dinámico (S0)

**Formulación base:**
```
S0[W] = −α T + β S_dens + γ S_deg
```

**Componentes:**

**T (Triángulos):**
```
T = Σ_{i<j<k} W_ij W_jk W_ki
```
- Término −αT favorece relaciones cerradas
- Factor α: parámetro de acoplamiento

**S_dens (Densidad):**
```
S_dens = Σ_{i<j} W_ij²
```
- Término +β penaliza pesos grandes
- Factor β: parámetro de acoplamiento (referencia de escala: β = 1)

**S_deg (Regularidad de conectividad):**
```
k_i = Σ_j W_ij                     (grado de nodo i)
k̄ = (1/N) Σ_i k_i                 (grado promedio)
S_deg = Σ_i (k_i − k̄)²
```
- Término +γ favorece homogeneidad de conectividad
- Factor γ: parámetro de acoplamiento
- Usa k̄ como cantidad emergente, NO fija k0

**Término adicional opcional (inicialmente desactivado):**
```
S_smooth = η Σ_ij W_ij C_ij
C_ij = Σ_k (W_ik − W_jk)²
```
- η = 0 en primera versión
- Podría favorecer geometría suave

### Distancia y conectividad

**Longitud de arista:**
```
l_ij = 1 / (W_ij + ε)
```
- Relación fuerte (W_ij grande) → separación pequeña
- ε: constante de regularización (valor no especificado en documento)

**Matriz de adyacencia con umbral:**
```
A_ij = 1  si W_ij > W_min
A_ij = 0  si W_ij ≤ W_min
```
- W_min: umbral de conectividad (valor numérico no especificado)
- Usado para evitar considerar aristas triviales

**Distancia global:**
```
d(i,j) = min sobre caminos γ:i→j de Σ_{(a,b)∈γ} l_ab
```
- Camino más corto en la red ponderada
- Calcula distancias sobre relaciones existentes (A_ij = 1)

### Dimensionalidad efectiva

**Bola (vecindario):**
```
B_i(r) = {j : d(i,j) ≤ r}
N_i(r) = |B_i(r)|    (número de nodos dentro de distancia r desde i)
```

**Dimensión efectiva (dependiente de escala):**
```
D_eff(r) = d ln N(r) / d ln r
```
- Relación esperada en espacio D-dimensional: N(r) ~ r^D
- Ejemplos:
  - N ~ r     ⟹ D_eff ≈ 1 (unidimensional)
  - N ~ r²    ⟹ D_eff ≈ 2 (bidimensional)
  - N ~ r³    ⟹ D_eff ≈ 3 (tridimensional)
- NO asumir D=3: descubrir qué valor emerge
- D_eff puede variar con la escala (D_micro ≠ D_macro posible)

### Dimensión espectral

**Probabilidad de retorno por difusión:**
```
P(σ) = probabilidad de retorno en tiempo abstracto σ
```
- σ: parámetro de tiempo de difusión (unidad no especificada)

**Dimensión espectral:**
```
D_s(σ) = −2 d ln P(σ) / d ln σ
```
- Obtenida por difusión sobre la red
- Puede revelar estructura que D_eff no detecta
- Criterio de robustez: D_eff ≈ 3 Y D_s ≈ 3 en rangos compatibles
- Nota: D_eff ≈ 3 solo NO demuestra geometría tridimensional

### Dinámica

**Ecuación de evolución:**
```
dW_ij/dτ = −∂S0/∂W_ij
```
- τ: parámetro abstracto (NO es tiempo físico)

**Actualización numérica (primera implementación):**
```
W_ij ← clip[W_ij − Δτ ∂S0/∂W_ij, 0, 1]
```
- Δτ: paso de evolución
- clip[·, 0, 1]: mantiene 0 ≤ W_ij ≤ 1

**Parametrización alternativa (posterior):**
```
W_ij = σ(θ_ij)
```
- σ: función sigmoide (definición exacta no dada)
- θ_ij: parámetros subyacentes
- Ventaja: mantiene automáticamente 0 < W_ij < 1

### Gradiente de S0

**Derivada de T respecto a arista (a,b):**
```
∂T/∂W_ab = Σ_{k≠a,b} W_ak W_kb
```

**Derivada de S_dens respecto a arista (a,b):**
```
∂S_dens/∂W_ab = 2 W_ab
```

**Derivada de S_deg respecto a arista (a,b):**
```
∂S_deg/∂W_ab = 2[(k_a − k̄) + (k_b − k̄)]
```

**Gradiente total de S0 respecto a arista (a,b):**
```
∂S0/∂W_ab = −α Σ_{k≠a,b} W_ak W_kb + 2β W_ab + 2γ[(k_a − k̄) + (k_b − k̄)]
```

### Espacio de parámetros

**Parámetros relevantes (normalización β = 1):**
- α/β: razón triángulos/densidad
- γ/β: razón regularidad/densidad
- η/β: razón suavidad/densidad (si se introduce S_smooth)

**Mapa de fases 2D:** (α/β, γ/β) → fase
**Mapa de fases 3D (con suavidad):** (α/β, γ/β, η/β) → fase

### Geometría efectiva

**Secuencia de emergencia:**
```
Ω → W → d → D → x → g_{μν}
```
- Ω: estructura pregeométrica
- W: matriz de pesos
- d: distancias globales
- D: dimensionalidad
- x: coordenadas reconstruidas
- g_{μν}: métrica efectiva

**Reconstrucción de métrica (posterior a geometría):**
```
d_ij² ≈ g_{μν} Δx^μ Δx^ν
```

**Curvatura de Riemann (posterior a métrica):**
```
R_{μνρσ} emergente de: δW → δd → δg → R
```

---

## 2. Restricciones metodológicas prohibidas ("NO hacer")

1. NO usar red cúbica, lattice 2D, lattice hexagonal como condición inicial
2. NO asignar coordenadas espaciales x, y, z a los nodos
3. NO introducir distancia convencional euclidiana como entrada
4. NO fijar conectividad objetivo k0 = 6 (ni ningún k0 predeterminado)
5. NO imponer geometría 3D preexistente en la condición inicial
6. NO considerar automáticamente todos los pares de nodos conectados con arista de gran peso
7. NO asumir D = 3; descubrir qué dimensionalidad emerge
8. NO omitir registro de la evolución temporal completa (guardar solo el estado final es insuficiente)
9. NO clasificar fases usando una sola variable observable
10. NO introducir suavidad (S_smooth con η > 0) desde el inicio
11. NO fijar la masa total Σ_{i<j} W_ij antes de explorar dinámicas
12. NO ajustar retroactivamente el modelo después de observar resultados
13. NO omitir semillas de aleatoriedad múltiples para cada condición
14. NO confundir la red abstracta Ω con un espacio físico incrustado
15. NO ignorar posibles transiciones de fase histéresis (ida ≠ vuelta en parámetros)
16. NO usar un único algoritmo para estimar dimensión sin validación cruzada
17. NO despreciar fase trivial (W_ij → 0): es información sobre dinámicas
18. NO omitir controles (p.ej., pesos puramente aleatorios sin dinámica)
19. NO calcular distancias sobre aristas triviales; usar umbral W_min
20. NO inicializar con distribuciones distintas a U(0,1) sin justificación explícita

---

## 3. Observables a medir por simulación y pasaporte de reproducibilidad

### Observables principales (O)

**Conectividad:**
- k̄: conectividad media (promedio de Σ_j W_ij sobre todos los nodos)
- σ_k: dispersión de conectividad (desviación estándar de k_i)

**Topología local:**
- C: coeficiente de clustering (definición exacta no dada en documento)

**Escala global:**
- L: longitud característica de caminos (definición exacta no dada)
- G: tamaño de componente gigante (fracción de nodos en la componente conexa principal)

**Geometría:**
- D_eff: dimensión efectiva (valor final o promedio en rango de escalas)
- D_s: dimensión espectral (valor final o promedio en rango de tiempo)

**Otros:**
- Curvatura: observable a calcular posterior a reconstrucción de métrica

### Pasaporte de reproducibilidad (P)

Cada simulación debe guardar:
```
P = {
  N,              parámetro de tamaño del sistema (número de nodos)
  α,              parámetro de triángulos
  β,              parámetro de densidad
  γ,              parámetro de regularidad
  η,              parámetro de suavidad (si aplicable)
  Ω0,             estado inicial (matriz W_0 o especificación de distribución)
  τ_max,          tiempo final de evolución
  Δτ,             paso de tiempo
  D_eff,          dimensión efectiva medida
  D_s,            dimensión espectral medida
  k̄,              conectividad media final
  σ_k,            dispersión de conectividad final
  C,              clustering final
  G,              tamaño de componente gigante
  curvatura,      observable de curvatura (si calculado)
  seed            semilla de aleatoriedad
}
```

---

## 4. Experimentos 0–3 y pruebas de robustez

### Experimento 0: Validación

**Propósito:** Verificar que el algoritmo de medición recupera dimensiones conocidas.

**Entrada:**
- Redes construidas manualmente con dimensionalidad conocida:
  - Red 1D (cadena, árbol linear)
  - Red 2D (lattice cuadrado u otro)
  - Red 3D (cubic lattice u otro)

**Procedimiento:**
1. Construir cada red con N nodos
2. Calcular N_i(r) para rango de r relevantes
3. Estimar D_eff(r)
4. Calcular D_s mediante difusión
5. Comparar con valor esperado

**Salida esperada:**
- D_eff ≈ 1 para red 1D
- D_eff ≈ 2 para red 2D
- D_eff ≈ 3 para red 3D
- D_s debe ser consistente con D_eff

**Criterio de éxito:** Recuperación correcta de D para todas las redes de control.

**Criterio de fallo:** Desviaciones > 0.5 en D_eff respecto a valor esperado; indicaría problemas en algoritmo de medición.

---

### Experimento 1: Baseline aleatorio

**Propósito:** Establecer línea base de propiedades de redes aleatorias sin dinámicas.

**Entrada:**
- W_ij ~ U(0,1) sin evolución dinámica
- N = 100–300 (según disponibilidad)
- Múltiples semillas (seed_1, seed_2, ..., seed_n; n ≥ 10 sugerido)

**Procedimiento:**
1. Generar matriz W con distribución uniforme
2. Medir: k̄, σ_k, C, L, G, D_eff, D_s
3. Calcular matriz de distancias d(i,j) sobre umbral W_min

**Salida esperada:**
- Valores de conectividad y distribución característica de redes Erdős–Rényi
- D_eff ≈ ∞ o muy alto (sin estructura geométrica)
- D_s ≈ ∞ o muy alto

**Criterio de éxito:** Baseline replicable y sin estructura geométrica emergente.

**Criterio de fallo:** Aparición espontánea de D_eff ≈ 3 sin dinámica (indicaría error en algoritmo).

---

### Experimento 2: Dinámica

**Propósito:** Estudiar evolución de estructura bajo dinámica S0.

**Entrada:**
- Pesos aleatorios iniciales: W_ij ~ U(0,1)
- Parámetros fijos a explorar: (α/β, γ/β) = valores seleccionados (p.ej. α/β = 0.1, γ/β = 0.5)
- β = 1 (normalización)
- η = 0 (sin suavidad en primera versión)
- Evolución: τ = 0 a τ_max con paso Δτ
- N = 100–300 inicialmente
- Semillas múltiples

**Procedimiento:**
1. Inicializar W ~ U(0,1)
2. Evolucionar según dW_ij/dτ = −∂S0/∂W_ij
3. En cada paso τ:
   - Registrar W(τ)
   - Calcular observables O(τ) = (k̄, σ_k, C, L, G, D_eff, D_s)
4. Continuar hasta convergencia o τ_max

**Salida esperada:**
- Evolución temporal de estructura
- Posible emergencia de fases distintas según parámetros
- Posible formación de componente gigante

**Criterio de éxito:**
- Evolución suave y reproducible entre semillas
- Convergencia a estado final estable
- Diferencias respecto a baseline aleatorio

**Criterio de fallo:**
- Divergencia numérica (W_ij fuera de [0,1])
- Comportamiento caótico no reproducible entre semillas
- Colapso a fase trivial W_ij → 0 para todos los parámetros

---

### Experimento 3: Barrido de parámetros

**Propósito:** Mapear región del espacio de parámetros y sus fases asociadas.

**Entrada:**
- Malla de parámetros: α/β ∈ {0, 0.1, 0.2, 0.3, ...} (resolución adaptable)
- γ/β ∈ {0, 0.1, 0.2, 0.3, ...} (resolución adaptable)
- Para cada punto (α/β, γ/β): múltiples semillas

**Procedimiento:**
1. Para cada punto (α, γ) de la malla:
   - Ejecutar Experimento 2 con esos parámetros
   - Registrar estado final y observables
2. Compilar tabla: (α/β, γ/β) → fase observada

**Salida esperada:**
- Mapa de fases 2D mostrando regiones de distintas fases
- Posibles transiciones entre fases

**Criterio de éxito:**
- Mapa reproducible y sin discontinuidades artificiales
- Identificación clara de regiones de fase

**Criterio de fallo:**
- Mapa ruidoso o no reproducible entre semillas
- Imposibilidad de clasificar fases de forma consistente

---

### Pruebas de robustez

#### 4.1 Robustez a semillas

**Procedimiento:**
- Para región de parámetros interesante (especialmente si D_eff ≈ 3):
  - Ejecutar 10–100 semillas independientes
  - Calcular media y desviación de D_eff, D_s, G, k̄

**Criterio de éxito:**
- Baja varianza entre semillas (σ_D_eff < 0.3 sugerido)
- Misma fase emergente independientemente de condición inicial
- Interpretación: fase es atractor robusto

**Criterio de fallo:**
- Alta varianza o fases distintas según semilla
- Indicaría resultado accidental, no estructural

#### 4.2 Robustez a tamaño del sistema

**Procedimiento:**
- Replicar simulación interesante con N = 100, 200, 300, 500, 1000
- Medir D_eff(r) como función de escala para cada N

**Criterio de éxito:**
- D_eff estable o ligeramente decreciente con N
- Convergencia aparente a valor límite

**Criterio de fallo:**
- D_eff divergente o fuertemente dependiente de N
- Indicaría efecto de tamaño finito relevante

#### 4.3 Histéresis y transiciones de fase

**Procedimiento:**
- Fijar γ/β = constante
- Variar α/β lentamente: α/β = 0 → 1 → 0 (ida y vuelta)
- En cada α: registrar observables (D_eff, G, k̄, etc.)

**Criterio de éxito:**
- Ida ≠ vuelta (histéresis presente) → transición de fase verdadera
- Ida = vuelta → equilibrio, no histéresis

**Criterio de fallo:**
- Ruido numérico que oculta histéresis real

#### 4.4 Invariancia ante renombramiento (permutación)

**Procedimiento:**
- Simular hasta convergencia con semilla inicial S0
- Obtener matriz final W
- Generar permutación P (reordenamiento de índices)
- Calcular W' = P W P^T
- Medir observables de W y W'

**Criterio de éxito:**
- D_eff(W) = D_eff(W') (igualdad o diferencia < 0.1)
- G(W) ≈ G(W')
- k̄(W) = k̄(W')
- Indicaría que geometría es topológica, no artefacto de indexación

**Criterio de fallo:**
- D_eff(W) ≠ D_eff(W')
- Indicaría que "dimensionalidad" es sensible a orden de índices (error grave)

#### 4.5 Ablación

**Procedimiento:**
- Replicar simulación con cada término de S0 aislado:
  - Solo −αT
  - Solo +βS_dens
  - Solo +γS_deg
  - −αT + βS_dens
  - −αT + γS_deg
  - −αT + βS_dens + γS_deg (completo)
- Usar mismo punto (α, γ) y semilla

**Criterio de éxito:**
- Cada término contribuye a aspecto distinto de estructura
- Término completo produce estructura más robusta que cualquier término solo

**Criterio de fallo:**
- Un solo término domina; otros no contribuyen
- Indicaría sobrecorrección o parámetros desbalanceados

#### 4.6 Coarse-graining (renormalización)

**Procedimiento:**
- Obtener estructura final Ω de simulación interesante
- Agrupar nodos en "supernodos" progresivamente:
  - Ω → Ω' (nivel 1)
  - Ω' → Ω'' (nivel 2)
- Calcular D_eff en cada nivel

**Criterio de éxito:**
- D_eff_micro ≈ D_eff_macro (o dependencia gradual suave)
- Indicaría estructura autoafín o fractal suave

**Criterio de fallo:**
- D_eff cambia abruptamente con coarse-graining
- Indicaría estructura frágil o artefacto de escala

---

## 5. Fases A–F con criterios cualitativos

### Fase A: Dispersa

**Descripción:** Red con baja conectividad, componentes pequeñas.

**Criterios cualitativos:**
- k̄ → 0 (conectividad muy baja)
- Componentes pequeñas y aisladas
- G << 1 (componente gigante ausente o trivial)
- Topología fragmentada

**Parámetros típicos:** α/β grande, γ/β bajo

---

### Fase B: Fragmentada

**Descripción:** Múltiples componentes grandes but distintas.

**Criterios cualitativos:**
- Varias componentes conexas grandes
- G ~ 0.1–0.5 (no hay componente que domine)
- k̄ intermedio
- No hay una estructura única dominante

**Parámetros típicos:** α/β intermedio, γ/β bajo–intermedio

---

### Fase C: Conectada y homogénea

**Descripción:** Una componente gigante con conectividad uniforme.

**Criterios cualitativos:**
- G ≈ 1 (casi todos los nodos en componente principal)
- Baja dispersión relativa de grados: σ_k/k̄ bajo
- k̄ moderado
- Topología regular

**Parámetros típicos:** α/β bajo–intermedio, γ/β alto

---

### Fase D: Altamente agrupada

**Descripción:** Estructura con clustering elevado y posible modularidad.

**Criterios cualitativos:**
- Coeficiente de clustering C elevado
- Posibles subcomunidades o módulos
- G ≈ 1
- Longitud de caminos L posiblemente elevada

**Parámetros típicos:** α/β alto (favorece triángulos), γ/β intermedio

---

### Fase E: Hiperdensa

**Descripción:** Red sobreconectada con pesos W_ij cercanos a 1.

**Criterios cualitativos:**
- k̄ muy alto (cercano a k̄_max = N–1)
- Matriz W aproximadamente llena
- G ≈ 1
- Geometría poco diferenciada

**Parámetros típicos:** α/β bajo, γ/β bajo, o término de masa activo

---

### Fase F: Candidata geométrica

**Descripción:** Estructura con propiedades sugiriendo geometría emergente (posible D ≈ 3).

**Criterios cualitativos:**
1. Componente gigante: G ≈ 1
2. Estructura local homogénea: σ_k/k̄ bajo
3. Escalas geométricas: N_i(r) muestra crecimiento monotónico suave
4. D_eff estable: D_eff(r) ≈ constante en rango de escalas (p.ej. r ∈ [r_min, r_max])
5. Dimensión espectral compatible: D_s ≈ D_eff en rango de tiempos compatibles
6. Robustez: propiedades persisten ante:
   - Variaciones pequeñas de parámetros α, γ
   - Múltiples semillas iniciales
   - Cambios en N (tamaño del sistema)
7. Reproducibilidad: mismo comportamiento con algoritmos distintos de medición
8. Ablación: propiedades sobreviven a omisión de términos secundarios
9. Transición de fase: evidencia de comportamiento crítico
10. Compatible con coarse-graining: microestructura complicada pero macroscopía suave
11. Geometría efectiva: distancias d(i,j) pueden incrustarse en espacio continuo
12. Métrica emergente: posible reconstrucción de tensor métrico g_{μν}

**Parámetros típicos:** rango específico de (α/β, γ/β) a identificar por Experimento 3

**Criterio de éxito (fase F):** Todas las propiedades 1–12 presentes simultáneamente.

**Criterio de fallo:** Ausencia de propiedades críticas o falta de robustez a semillas/N.

---

## 6. Arquitectura de carpetas propuesta (módulo por módulo)

Basado en sección 48 de OMEGA_MASTER.md.

### Estructura completa

```
OMEGA/
├── network/
├── dynamics/
├── geometry/
├── phases/
├── coarse_graining/
├── curvature/
├── causality/
├── statistics/
└── experiments/
```

### Responsabilidades por módulo y archivo

#### network/

**Archivo: initialization.py**
- Responsabilidad: Generar matriz W inicial
- Funciones esperadas:
  - `generate_random_weights(N, seed)` → W ~ U(0,1) simétrica
  - Aplicar simetrización tras generación aleatoria
  - Asegurar W_ii = 0

**Archivo: weights.py**
- Responsabilidad: Manipulación de matriz W (lectura, escritura, actualización)
- Funciones esperadas:
  - `clip_weights(W)` → mantener 0 ≤ W_ij ≤ 1
  - `update_weight(W, i, j, delta)` → W_ij ← W_ij + delta
  - `save_weights(W, path)` → serializar
  - `load_weights(path)` → deserializar

**Archivo: topology.py**
- Responsabilidad: Análisis topológico de red
- Funciones esperadas:
  - `adjacency_matrix(W, threshold)` → A_ij = 1 si W_ij > threshold else 0
  - `degree(W)` → vector k
  - `mean_degree(W)` → k̄
  - `degree_variance(W)` → σ_k
  - `giant_component(W)` → G (fracción o tamaño)
  - `connected_components(W)` → lista de componentes

#### dynamics/

**Archivo: functional.py**
- Responsabilidad: Definir y evaluar funcional S0
- Funciones esperadas:
  - `triangles(W)` → T = Σ_{i<j<k} W_ij W_jk W_ki
  - `density(W)` → S_dens = Σ_{i<j} W_ij²
  - `regularity(W)` → S_deg = Σ_i (k_i − k̄)²
  - `smoothness(W)` → S_smooth (si η > 0)
  - `functional_S0(W, alpha, beta, gamma, eta)` → S0[W]
  - `control_term(W, mu)` → −μ Σ_{i<j} W_ij (si usado)

**Archivo: gradient.py**
- Responsabilidad: Calcular gradiente de S0 respecto a W
- Funciones esperadas:
  - `grad_triangles(W)` → ∂T/∂W_{ij} para todas (i,j)
  - `grad_density(W)` → ∂S_dens/∂W_{ij}
  - `grad_regularity(W)` → ∂S_deg/∂W_{ij}
  - `grad_smoothness(W)` → ∂S_smooth/∂W_{ij}
  - `grad_functional(W, alpha, beta, gamma, eta)` → ∂S0/∂W_{ij}

**Archivo: evolution.py**
- Responsabilidad: Integrador temporal de dinámica
- Funciones esperadas:
  - `evolve_step(W, grad_S0, dt)` → W ← clip[W − dt · grad_S0, 0, 1]
  - `evolve_sigmoidal(theta, grad_theta, dt)` → θ ← θ − dt · grad_theta (si parametrización sigmoide)
  - `sigmoid(theta)` → W = σ(θ) (si parametrización sigmoide)
  - `evolve_full(W0, alpha, beta, gamma, tau_max, dt, convergence_criterion)` → (W_f, W_history)

#### geometry/

**Archivo: distances.py**
- Responsabilidad: Calcular distancias en red
- Funciones esperadas:
  - `edge_lengths(W, epsilon)` → l_ij = 1/(W_ij + ε) para A_ij = 1
  - `shortest_paths(W, threshold)` → d(i,j) para todos los pares (Dijkstra u otro)
  - `distance_matrix(W, threshold)` → matriz D con d(i,j)
  - `neighborhood(D, i, r)` → B_i(r) = {j : d(i,j) ≤ r}

**Archivo: dimension.py**
- Responsabilidad: Estimar D_eff
- Funciones esperadas:
  - `neighborhood_sizes(D, r_values)` → N_i(r) para todos los i y r
  - `effective_dimension(N_values, r_values)` → D_eff = d ln N / d ln r (regresión lineal u otro)
  - `dimension_by_scale(D_matrix, scale_range)` → D_eff(r) como función de r
  - `validate_dimension_1D()`, `validate_dimension_2D()`, `validate_dimension_3D()` para Experimento 0

**Archivo: spectral.py**
- Responsabilidad: Estimar D_s mediante difusión
- Funciones esperadas:
  - `random_walk_matrix(A)` → matriz de transición P para random walk sobre A
  - `return_probability(P, time_values)` → P(σ) para σ en time_values
  - `spectral_dimension(P_values, time_values)` → D_s = −2 d ln P / d ln σ
  - `diffusion_profile(A, steps)` → evolución de distribución de random walker

#### phases/

**Archivo: scan.py**
- Responsabilidad: Ejecutar barrido sistemático de parámetros (Experimento 3)
- Funciones esperadas:
  - `parameter_grid(alpha_min, alpha_max, gamma_min, gamma_max, resolution)` → malla
  - `scan_phase_diagram(grid, N, tau_max, seeds_per_point)` → tabla resultados
  - `run_parameter_point(alpha, beta, gamma, N, seeds)` → lista de observables finales

**Archivo: classification.py**
- Responsabilidad: Clasificar fases observadas según criterios cualitativos
- Funciones esperadas:
  - `classify_phase(observables)` → etiqueta fase (A, B, C, D, E, F)
  - `criteria_phase_A(obs)` → True/False si cumple fase A
  - `criteria_phase_F(obs)` → True/False si cumple fase F (incluyendo todas 12 propiedades)

**Archivo: stability.py**
- Responsabilidad: Analizar estabilidad y robustez de fases
- Funciones esperadas:
  - `stability_across_seeds(alpha, beta, gamma, N, num_seeds)` → (media_obs, std_obs)
  - `stability_across_sizes(alpha, beta, gamma, N_values, seed)` → observables(N)
  - `hysteresis_test(beta, gamma, alpha_path_forward, alpha_path_backward)` → detectar histéresis

#### coarse_graining/

**Archivo: (genérico)**
- Responsabilidad: Técnicas de renormalización y análisis multiescala
- Funciones esperadas:
  - `coarse_grain_level1(W)` → W' (nivel más grueso)
  - `dimension_across_scales(W)` → D_eff(escala_0), D_eff(escala_1), ...
  - `autoaffinity_test(D_values)` → True/False si estructura es autoafín

#### curvature/

**Archivo: (genérico)**
- Responsabilidad: Estimar curvatura de Riemann (posterior a métrica)
- Funciones esperadas:
  - `metric_reconstruction(d_matrix, coordinates)` → g_{μν}
  - `riemann_tensor(g)` → R_{μνρσ}
  - `scalar_curvature(R)` → R (escalar)

#### causality/

**Archivo: (genérico)**
- Responsabilidad: Estructura causal (posterior a geometría espacial)
- Funciones esperadas:
  - Ampliación futura: no implementar en versión 1

#### statistics/

**Archivo: (genérico)**
- Responsabilidad: Análisis estadístico de resultados
- Funciones esperadas:
  - `compute_observable_statistics(results)` → mean, std, correlaciones
  - `plot_phase_diagram()` → visualización mapa de fases
  - `reproducibility_passport(sim_result)` → P (estructura completa)

#### experiments/

**Archivo: test_00_validation.py**
- Responsabilidad: Experimento 0
- Contenido:
  - Construir redes 1D, 2D, 3D de referencia
  - Verificar D_eff ≈ 1, 2, 3 respectivamente
  - Verificar D_s consistente

**Archivo: test_01_random.py**
- Responsabilidad: Experimento 1
- Contenido:
  - Generar W ~ U(0,1) sin dinámica
  - Medir observables
  - Compilar línea base

**Archivo: test_02_dynamics.py**
- Responsabilidad: Experimento 2
- Contenido:
  - Evolucionar desde pesos aleatorios
  - Registrar trazas temporales
  - Guardar evolución completa

**Archivo: test_03_phase_diagram.py**
- Responsabilidad: Experimento 3
- Contenido:
  - Ejecutar barrido de parámetros
  - Compilar mapa de fases
  - Clasificar regiones

**Archivo: test_04_dimension.py**
- Responsabilidad: Prueba específica de robustez dimensional
- Contenido:
  - Múltiples semillas
  - Múltiples tamaños N
  - Calcular media y varianza de D_eff

**Archivo: test_05_robustness.py**
- Responsabilidad: Pruebas de robustez (semillas, tamaño, histéresis, permutación, ablación, coarse-graining)
- Contenido:
  - Subrutinas para cada prueba
  - Compilar reporte de robustez

---

## 7. Ambigüedades y huecos sin definir en OMEGA_MASTER.md

### Definiciones numéricas no especificadas

1. **Valor de ε en longitud de arista:** l_ij = 1/(W_ij + ε)
   - ¿Cuál es ε exacto? (p.ej. 1e−6, 1e−10, pequeño pero no cero)

2. **Valor de umbral W_min:** A_ij = 1 si W_ij > W_min, else 0
   - ¿Cuál es W_min? (p.ej. 0.01, 0.05, 0.1, o dinámico)

3. **Paso de tiempo Δτ:** W_ij ← clip[W_ij − Δτ ∂S0/∂W_ij, 0, 1]
   - ¿Valor específico o adaptativo?
   - ¿Cómo seleccionar para convergencia sin divergencia?

4. **Criterio de convergencia:** ¿Cuándo detener evolución?
   - ¿||dW/dτ|| < ε_conv?
   - ¿Tiempo máximo τ_max fijo?
   - ¿Umbral de change_rate en observables?

5. **Número de semillas por punto:** Experimentos 1–3
   - Documento sugiere "múltiples semillas" but no especifica número mínimo
   - ¿10, 100, 1000 semillas?

### Rangos no especificados

6. **Rango de r para D_eff:** B_i(r) = {j : d(i,j) ≤ r}
   - ¿r mínimo? (p.ej. r_min = 1 o r_min dinámica)
   - ¿r máximo? (p.ej. r_max = diámetro red)
   - ¿Resolución de r? (número de valores entre r_min y r_max)

7. **Rango de σ para D_s:** Difusión temporal
   - ¿σ mínimo (primeros pasos)?
   - ¿σ máximo (tiempo de mezcla)?
   - ¿Resolución temporal?

8. **Resolución de malla de parámetros:** Experimento 3
   - Documento dice "resolución adaptable al coste"
   - ¿0.01, 0.05, 0.1 en (α/β, γ/β)?

9. **Rango de N:** Tamaño del sistema
   - Documento menciona N = 100–300, luego 500–1000, luego 10³–10⁴+
   - ¿Cuál es N para experimentos base?

### Definiciones cualitativas sin formalización

10. **Definición exacta de C (coeficiente de clustering):**
    - ¿C = Σ_i c_i / N donde c_i = (triángulos en i) / (posibles triángulos)?
    - ¿Cómo se define en red ponderada?

11. **Definición exacta de L (longitud característica de caminos):**
    - ¿L = promedio de d(i,j) sobre todos los pares?
    - ¿Diámetro = max d(i,j)?

12. **Definición exacta de P(σ) (probabilidad de retorno):**
    - ¿Es probabilidad de walker regresar a punto inicial exacto?
    - ¿O probabilidad de estar en vecindario original?
    - ¿Cómo se define en red ponderada vs. no ponderada?

### Parametrizaciones no detalladas

13. **Función sigmoide σ(θ) en parametrización alternativa:**
    - ¿σ(θ) = 1 / (1 + exp(−θ))?
    - ¿Rango de θ ∈ ℝ?
    - ¿Cómo se inicializa θ desde W_0?

14. **Normalización de pesos en término S_deg:**
    - ¿Se normaliza k_i por k_max = N−1?
    - ¿Se usa (k_i − k̄)² o ((k_i − k̄)/k̄)²?

15. **Normalización de términos en S0:**
    - Cada término (T, S_dens, S_deg, S_smooth) tiene dimensión distinta
    - ¿Hay factores de normalización implícitos?
    - ¿Se normalizan por N o por máximo posible?

### Construcción de estructuras de referencia (Experimento 0)

16. **Construcción exacta de red 1D, 2D, 3D:**
    - Para 1D: ¿cadena líneal o árbol?
    - Para 2D: ¿lattice cuadrado o triangular?
    - Para 3D: ¿cubic lattice simple o estructura más compleja?
    - ¿Condiciones de contorno periódicas o abiertas?

### Criterios de clasificación de fases

17. **Umbral numérico para "G ≈ 1":** Fase F criterio 1
    - ¿G > 0.9, > 0.95, > 0.99?

18. **Umbral numérico para "σ_k/k̄ bajo":** Fase F criterio 2
    - ¿σ_k/k̄ < 0.1, < 0.2, < 0.5?

19. **Definición de "rango de escalas compatible":** Fase F criterio 4
    - ¿Qué significa "rango de escalas"? (p.ej. r ∈ [2, N/2])
    - ¿Qué significa "estable"? (variación < 0.1, 0.3, 0.5)

20. **Definición de "pequeñas variaciones de parámetros":** Fase F criterio 6
    - ¿Δα/α < 0.1, 0.2, 0.5?
    - ¿Δγ/γ similar?

### Métricas y comparaciones

21. **¿Cómo se define "diferencia < 0.1" en invariancia ante permutación?** (Prueba 4.4)
    - ¿Diferencia absoluta |D_eff(W) − D_eff(W')|?
    - ¿Relativa (|D_eff(W) − D_eff(W')|) / D_eff(W)?

22. **¿Cómo se agrega información de múltiples escalas?** D_eff(r)
    - ¿Promediar D_eff(r) sobre rango r ∈ [r_min, r_max]?
    - ¿Usar regresión lineal en log-log?

### Posterior a estructura geométrica

23. **Algoritmo de reconstrucción de coordenadas:**
    - "Multidimensional scaling u otros" (sección 35)
    - ¿MDS euclidiano o generalizado?
    - ¿Qué espacio objetivo? (ℝ³ asumido pero no forzado)

24. **Procedimiento de ajuste de métrica:**
    - "d_ij² ≈ g_{μν} Δx^μ Δx^ν; ajustar métrica local g_eff" (sección 36)
    - ¿Ajuste global o local por punto?
    - ¿Cuál es criterio de "buen ajuste"?

### Implementación numérica

25. **Precisión numérica de gradiente:**
    - ¿Gradiente analítico o numérico (diferencias finitas)?
    - Si numérico: ¿paso de diferenciación?

26. **Manejo de matrices grandes (N > 10⁴):**
    - Documento menciona O(N²) almacenamiento
    - ¿Usar matriz densa o sparse?
    - ¿Algoritmo especial para N grande?

---

**Fin de SPEC_EXTRACT.md**

Nota: Este documento lista únicamente ambigüedades presentes en OMEGA_MASTER.md. La sección 7 NO resuelve estas ambigüedades; solo las identifica como huecos a aclarar en futuros refinamientos.
