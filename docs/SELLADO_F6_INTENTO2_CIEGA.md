# Ω, F6 intento 2: estructuras a ciegas (sellado)

- **Autor:** guardián de reserva nuevo (subagente que no participa en el diseño ni vio datos del intento 2). Fecha: 2026-10-18.
- **Uso:** este archivo se abre solo después de congelar el preregistro del intento 2 (Rev. 12 §3). Antes de eso solo se publica su SHA-256.
- **Referencias leídas:** Rev. 12 §3 (requisitos), resultados del intento 1, `SELLADO_F6_CONFIRMACION_CIEGA.md` (C0 a C7, ya vistas: ninguna se repite ni se varía trivialmente aquí), `OMEGA_F6_EXPLORACION.md` (estructuras vistas), Rev. 11 §8.1 (zona gris ε ≤ 0.20), `omega/causality/kasner.py`, `tools/f6_confirmacion.py`.
- **Convenciones (iguales a las del sellado del intento 1):** dimensión 1 + 3; T = `DIAMOND_HALF_HEIGHT` = 0.5; columna 0 = t; Δ = x_j − x_i. «Diamante estándar» = `sprinkle_diamond(n, 4, rng)`, región {|t| + ‖x‖₂ ≤ T}, volumen V = π(2T)⁴/24 = π/24. C[i, j] = True ⇔ i precede estrictamente a j.
- **Conos (D6 a D8):** gauge F (convexo, positivamente homogéneo, F(v) > 0 si v ≠ 0). Relación: j ≻ i ⇔ Δt > 0 y Δt ≥ F(Δx). Región = diamante causal del propio cono entre a = (−T, 0) y b = (T, 0): {t + T ≥ F(x) y T − t ≥ F(−x)} (si F es par, es {|t| ≤ T − F(x)}). Puntos uniformes por rechazo en la caja indicada en cada estructura, como `sprinkle_cone` (lotes de max(64, 32·faltantes) candidatos uniformes en la caja; se aceptan en orden los que caen en la región; se toman los primeros n). Transitividad: F es subaditiva por ser convexa y homogénea.
- **Deformación:** ε = max_{‖u‖₂=1} |F(u) − 1| (definición de E-R11-1, Rev. 11 §8.1). Se da la fórmula cerrada y el valor numérico comprobado sobre 2·10⁶ direcciones al azar más los 26 ejes y diagonales. Se informa además, sin efecto en el papel, ε_esc = (max F − min F)/(max F + min F) (mejor reescalado uniforme).
- **Papeles:** positiva (debe clasificarse espacio-tiempo 4D), negativa (debe rechazarse), zona gris (ε ≤ 0.20, sin voto). **En este conjunto no hay zona gris.**

## Entropías y semillas

- **20261082:** generación de puntos (y todo lo que el texto llame «generador principal»).
  `gen = np.random.Generator(np.random.PCG64(np.random.SeedSequence(20261082, spawn_key=(322, k, N, s))))`.
- **20261083:** azar auxiliar (elección de pares agregados en D4; elección de originales y permutación en D5).
  `aux = np.random.Generator(np.random.PCG64(np.random.SeedSequence(20261083, spawn_key=(322, k, N, s))))`.
- k = índice de la tabla (0..8), N = tamaño de la estructura, s = realización (0, 1, …).
- Las estructuras que no mencionan `aux` no lo consumen. En D4 y D5, `gen` se consume completo (puntos) antes de usar `aux`.
- No se usa ninguna otra entropía para estas estructuras. N, número de realizaciones, rejilla de m y regla de voto son los del preregistro congelado; estas definiciones valen para cualquier N.
- **Piloto declarado (no es evidencia):** para fijar la magnitud de D4 se contaron solo relaciones (ningún Δr, intercepto ni K) en dos realizaciones con semillas `default_rng(999/998)`, ajenas a las entropías de arriba.

## Lista (orden = índice k)

| k | Nombre | Papel | d verdadera | ε |
|---|---|---|---|---|
| 0 | D0 Kasner de vacío u = 2 permutado, p = (6/7, −2/7, 3/7), t ∈ [0.70, 1.35] | positiva | 4 | — (no conformemente plano) |
| 1 | D1 Minkowski con cizalla galileana (cono redondo inclinado), b = (0.25, −0.15, 0.10) | positiva | 4 | crudo 0.4456; afín 0 |
| 2 | D2 FRW con factor de escala exponencial a(η) = e^η (no potencial) | positiva | 4 | 0 |
| 3 | D3 FRW de materia rígida a ∝ η^{1/2}, η ∈ [0.5, 1.5] (R < 0 variable) | positiva | 4 | 0 |
| 4 | D4 Minkowski con relaciones agregadas a escala de discreción + cierre transitivo | negativa | — | — (impostor discreto) |
| 5 | D5 Minkowski con 10% de elementos gemelos (duplicados) | negativa | — | — (impostor discreto) |
| 6 | D6 cono «esfera con una tapa plana» (no reversible), s = 0.80 | negativa | — | 0.2500 |
| 7 | D7 cono «huso» (esfera + dos puntas cónicas en ±z), L = 1.45 | negativa | — | 0.3103 |
| 8 | D8 cono «suma de normas» F = 0.68‖x‖₂ + 0.32‖x‖₁ | negativa | — | 0.2342 |

Resumen por papel: **4 positivas, 5 negativas, 0 zona gris.** Ninguna positiva tiene ε en (0, 0.20].

## D0. Kasner de vacío, otros exponentes y otra ventana (positiva, no conformemente plana)

- **Métrica:** ds² = −dt² + Σᵢ t^{2pᵢ} dxᵢ², **p = (6/7, −2/7, 3/7)** (Lifshitz-Khalatnikov u = 2: (−2, 3, 6)/7, permutado para que el eje contraído sea x₂, no x₁ como en C0).
  - Σp = (6 − 2 + 3)/7 = 1; Σp² = (36 + 4 + 9)/49 = 1. **Vacío exacto** (Ricci = 0, R = 0; Weyl ≠ 0 porque p no es (1, 0, 0) ni permutación).
- **Región:** J⁺(a) ∩ J⁻(b), a = (0.70, 0, 0, 0), b = (1.35, 0, 0, 0). **t_lo = 0.70, t_hi = 1.35.**
- **Curvatura:** Kretschmann K = −16 p₁p₂p₃ / t⁴ = (576/343)/t⁴ ≈ 1.6793/t⁴ (comprobado con la suma de componentes ortonormales 4Σ[pᵢ(pᵢ−1)]² + 4Σ_{i<j}(pᵢpⱼ)²). Longitud de curvatura ℓ = K^{−1/4} ≈ 0.8785·t (0.61 a 1.19 en la ventana). Altura propia τ = 0.65; τ/ℓ(t = 1.025) ≈ 0.72. Densidad en coordenadas ∝ t: variación ×1.93 en todo el intervalo (≤ ×2, S5).
- **Generador (solo `gen`):** `p = validate_kasner_exponents((6/7, -2/7, 3/7))`; `x = sprinkle_kasner_interval(n, p, 0.70, 1.35, gen, k=800)`; `C = kasner_causal_matrix(x, p, k=800)`.
- **Por qué es positiva:** solución de vacío de Einstein con curvatura de Weyl; la relación causal no es la de Minkowski. Distinta de C0 en exponentes (u = 2 frente a u = 3/2), en el eje contraído y en la ventana temporal (curvatura algo más fuerte en el borde inferior).

## D1. Minkowski con cizalla galileana (positiva, invariancia lineal que mezcla t y x)

- **Puntos (solo `gen`):** `x = sprinkle_diamond(n, 4, gen)` (diamante estándar de coordenadas).
- **Relación:** j ≻ i ⇔ Δt > 0 y ‖Δx − b·Δt‖₂ ≤ Δt, con **b = (0.25, −0.15, 0.10)** (|b| ≈ 0.3082 < 1). Implementación equivalente: x' = x − b·t (por fila) y `causal_matrix` de Minkowski sobre (t, x').
- **Por qué es positiva:** (t, x) ↦ (t, x − bt) es lineal con jacobiano 1: es Minkowski exacto, muestreado uniformemente en volumen propio. La región, vista en las coordenadas de Minkowski, es un diamante cizallado (un convexo que no es un diamante causal). El cono en coordenadas es una bola redonda desplazada (no reversible: F(v) ≠ F(−v)).
- **ε crudo:** F(u) ∈ [1/(1+|b|), 1/(1−|b|)] = [0.7644, 1.4456] ⇒ **ε = 0.4456** (fuera de la zona gris); ε_esc = 0.308; ε afín (admitiendo transformaciones lineales de espacio-tiempo) = 0. Distinta de C1 (allí la transformación era solo espacial y reversible).

## D2. FRW con factor de escala exponencial (positiva, curvatura que cambia en el tiempo, no potencial)

- **Métrica:** ds² = a(η)²(−dη² + dx²), **a(η) = e^{η}**, con η = t ∈ [−0.5, 0.5] (coordenada temporal del diamante estándar, sin desplazamiento).
- **Muestreo (solo `gen`), patrón de `frw_power_points`:** repetir: `cand = sprinkle_diamond(4n, 4, gen)`; aceptar las filas con `gen.random(len(cand)) < exp(4·(cand[:,0] − 0.5))`; acumular; tomar las primeras n aceptadas. Densidad en coordenadas ∝ a⁴ = e^{4η} (√−g = a⁴).
- **Relación:** Minkowski (`causal_matrix`), por invariancia conforme.
- **Curvatura:** R = 6a''/a³ = 6e^{−2η}: de 6e ≈ 16.3 (η = −0.5) a 6/e ≈ 2.2 (η = +0.5), siempre positiva y decreciente (factor e² ≈ 7.4 en la ventana). Ricci ≠ 0 y Weyl = 0. Es el fluido perfecto con w = −1/3 (ρ ∝ a^{−2}; a ∝ tiempo cósmico). Variación total de densidad e⁴ ≈ 55 (menor que FRW polvo vista, 3⁸).
- **Por qué es positiva:** espacio-tiempo lorentziano suave, conformemente plano, sprinkling uniforme en volumen propio. No es a ∝ η ni a ∝ η² (ni ninguna potencia de η): pone a prueba el riesgo declarado de FRW en el colapso con un R(η) de otra forma funcional.

## D3. FRW de materia rígida (positiva, R < 0 que cambia en el tiempo)

- **Métrica:** ds² = a(η)²(−dη² + dx²), **a ∝ η^{1/2}**, con η = t + 1 ∈ [0.5, 1.5].
- **Muestreo (solo `gen`):** igual que D2 con probabilidad de aceptación `((cand[:,0] + 1)/1.5)**2` (densidad ∝ a⁴ = η²). Equivale a `frw_power_points(n, 4, 2.0, gen)` de `tools/f7.py`.
- **Relación:** Minkowski (`causal_matrix`).
- **Curvatura:** R = 6a''/a³ = −(3/2)·η^{−3}: de −12 (η = 0.5) a −0.44 (η = 1.5). Fluido rígido w = 1 (ρ ∝ a^{−6}). Variación total de densidad ×9.
- **Por qué es positiva:** espacio-tiempo suave conformemente plano. Potencia 1/2, distinta de 1 (radiación) y 2 (polvo); curvatura **negativa** variable, el caso que ninguna FRW vista cubre.

## D4. Minkowski con relaciones agregadas a escala de discreción y cierre transitivo (negativa, defecto discreto)

- **Construcción por realización:**
  1. `x = sprinkle_diamond(n, 4, gen)`; C = `causal_matrix(x)`.
  2. Escala de discreción: ℓ_N = (V/n)^{1/4} con V = π/24 (ℓ ≈ 0.0756 para n = 4000; 0.0636 para n = 8000). δ² = (1.0·ℓ_N)².
  3. **Candidatos:** pares (i, j) con t_j − t_i > 0, C[i, j] = False y ‖x_j − x_i‖₂² − (t_j − t_i)² ≤ δ² (separados como espaciales, pero a menos de una longitud de discreción del cono de luz; condición invariante de Lorentz). Listados en orden lexicográfico de (i, j) (el de `np.nonzero` sobre la matriz booleana de candidatos, por filas).
  4. **Con `aux`:** `u = aux.random(número_de_candidatos)` en ese orden; se agregan los candidatos con u < **q = 0.10**. S = C ∪ {pares agregados}.
  5. **C'' = cierre transitivo de S**: iterar S ← S ∪ (S∘S) (producto booleano, como `transitive_step`) hasta que no cambie. La prueba usa C''.
- **Por qué es orden parcial:** todo par de S cumple t_j > t_i; la composición de pares con Δt > 0 tiene Δt > 0, así que todo par de C'' tiene t_j > t_i ⇒ irreflexiva y antisimétrica (acíclica); transitiva por construcción.
- **Verificación obligatoria por realización:** `is_partial_order(C'')` y `transitive_step(C'') == C''`; informar número de candidatos, pares agregados, relaciones de C, de S y de C'' y la fracción nueva (|C''| − |C|)/|C|, y links de C y de C''.
- **Fracción afectada (piloto declarado, solo conteos):** unos 25 candidatos por elemento con n = 4000 (35 con n = 8000); pares agregados ≈ 1.3% de las relaciones de C; con el cierre, **relaciones nuevas ≈ 6% de |C|** (n = 4000). Se informa el valor exacto de cada realización.
- **Por qué es negativa:** no es el sprinkling de ninguna variedad lorentziana (un cono «borroso» a la escala de los elementos, que además arrastra relaciones por el cierre). Como la escala del defecto es ℓ_N, en un intervalo de m elementos el defecto depende de m y no de m/N, como C5, pero con el signo opuesto (agrega relaciones: imita R > 0, no R < 0) y sin borrar nada.

## D5. Minkowski con 10% de elementos gemelos (negativa, defecto discreto por duplicación)

- **Construcción por realización (n total):**
  1. n_dup = round(0.10·n), n_base = n − n_dup.
  2. `x = sprinkle_diamond(n_base, 4, gen)`; C = `causal_matrix(x)` (n_base × n_base).
  3. **Con `aux`, en este orden:** `orig = aux.choice(n_base, size=n_dup, replace=False)`; luego `perm = aux.permutation(n)`.
  4. φ: {0..n−1} → {0..n_base−1}, φ(e) = e para e < n_base y φ(n_base + q) = orig[q]. Relación previa C₀[u, v] = C[φ(u), φ(v)] (un gemelo tiene exactamente el pasado y el futuro de su original y no se relaciona con él).
  5. **C' = C₀[perm][:, perm]** (reordenamiento aleatorio de etiquetas, para que los gemelos no queden al final). La prueba usa C'.
- **Por qué es orden parcial:** u ≺ v ⇔ φ(u) ≺ φ(v) en C: irreflexiva (φ(u) ⊀ φ(u)), antisimétrica y transitiva heredadas de C. La permutación no cambia el orden.
- **Verificación obligatoria:** `is_partial_order(C')` y `transitive_step(C') == C'`; informar n_dup y el número de pares de elementos con el mismo pasado y el mismo futuro (debe ser n_dup).
- **Fracción afectada:** 10% de los elementos son gemelos (n_dup = 400 con n = 4000; 800 con n = 8000); ningún original se duplica dos veces.
- **Por qué es negativa:** en un sprinkling de una variedad lorentziana (distinguidora) dos elementos distintos nunca tienen el mismo pasado y el mismo futuro; aquí un décimo de los elementos lo tiene. El defecto es microscópico y su efecto en un intervalo depende del número de gemelos dentro, es decir de m.

## D6. Cono «esfera con una tapa plana» (negativa, no reversible, cerca de la frontera por fuera)

- **Gauge:** F(x) = max(‖x‖₂, (x·n)/s), con **n = (2, −1, 2)/3** (unitario) y **s = 0.80**. Bola unidad = bola euclídea cortada por un único semiespacio {x·n ≤ 0.80}: una sola cara plana, no simétrica por x ↦ −x.
- **Región:** {t + T ≥ F(x) y T − t ≥ F(−x)}, con F(−x) = max(‖x‖₂, −(x·n)/s). Caja de rechazo [−T, T]⁴ (F ≥ ‖x‖₂ ≥ ‖x‖_∞).
- **Relación:** Δt > 0 y Δt ≥ F(Δx) (con Δx = x_j − x_i, sin valor absoluto).
- **ε = 1/s − 1 = 0.2500** exacto (máximo 1.25 en la dirección n; mínimo 1 en todo el hemisferio opuesto). Numérico: 0.25000. ε_esc = 0.111.
- **Por qué es negativa:** la bola no es un elipsoide ni un elipsoide desplazado (tiene una cara plana con arista circular); ε > 0.20. Ninguna estructura vista rompe la simetría x ↦ −x con una forma no elipsoidal.

## D7. Cono «huso» (negativa, puntas cónicas)

- **Bola unidad:** envolvente convexa de la bola euclídea unidad y de los puntos (0, 0, ±L), **L = 1.45**.
- **Gauge (fórmula cerrada):** para x ≠ 0, r = ‖x‖₂, cos θ = |x₃|/r, θ ∈ [0, π/2], φ* = arccos(1/L):
  - si θ ≥ φ*: F(x) = r;
  - si θ < φ*: F(x) = r·(cos θ + sin θ·√(L² − 1))/L.
  (Deducida de la función soporte h(ν) = max(1, L|ν₃|); convexidad comprobada numéricamente en 2·10⁶ pares, violación máxima 4·10⁻¹¹, redondeo.)
- **Región:** {|t| ≤ T − F(x)} (F par). Caja de rechazo **[−T, T] × [−T, T] × [−T, T] × [−L·T, L·T]** (columnas t, x, y, z; F ≥ |x|, F ≥ |y|, F ≥ |z|/L).
- **ε = 1 − 1/L = 9/29 ≈ 0.3103** (mínimo 1/L en el eje z; máximo 1). Numérico: 0.31034. ε_esc = 0.184.
- **Por qué es negativa:** la bola tiene dos vértices cónicos (no es suave ni elipsoidal), ε > 0.20. Distinta del cono cilíndrico (allí hay una arista circular, aquí puntas), de los conos L^p y de los poliedros vistos.

## D8. Cono «suma de normas» (negativa, la más cercana a la frontera)

- **Gauge:** **F(x) = 0.68·‖x‖₂ + 0.32·‖x‖₁** (convexo como suma de normas; pliegues sobre los planos coordenados).
- **Región:** {|t| ≤ T − F(x)} (F par). Caja de rechazo [−T, T]⁴ (F ≥ ‖x‖₂).
- **ε = 0.32·(√3 − 1) ≈ 0.2342** (mínimo 1 en los ejes; máximo 0.68 + 0.32·√3 ≈ 1.2342 en las diagonales). Numérico: 0.23423. ε_esc = 0.105.
- **Por qué es negativa:** ε entre 0.22 y 0.35, por fuera de la zona gris; no pertenece a la familia L^p ni a las normas anidadas ni a los poliedros vistos.

## Regla de lectura (recordatorio, no criterio)

El criterio (intercepto y colapso), los umbrales, N, las realizaciones y la rejilla son los del preregistro congelado. Para el veredicto: D0, D1, D2 y D3 deben aceptarse como espacio-tiempo 4D; D4, D5, D6, D7 y D8 deben rechazarse. No hay estructura en zona gris. Si una estructura queda «no medible» según la regla congelada, cuenta como fallo de su papel salvo que el preregistro diga otra cosa antes de abrir este archivo. Si falla una verificación obligatoria de D4 o D5, la realización se rehace con s siguiente solo si el preregistro lo permite; si no, se informa como error de implementación antes de cualquier veredicto.

## Cobertura de los requisitos de Rev. 12 §3

- Positiva no conformemente plana distinta de C0: **D0** (y D1 como positiva con ε crudo fuera de la zona gris, por invariancia lineal no vista).
- Positiva con curvatura que cambia en el tiempo y no es a ∝ η ni a ∝ η²: **D2** (exponencial, R > 0) y **D3** (potencia 1/2, R < 0).
- Defecto discreto distinto del borrado de links: **D4** (relaciones agregadas + cierre) y **D5** (duplicación).
- Negativas tipo cono cerca de la frontera por fuera, ε ∈ [0.22, 0.35], distintas de las vistas: **D6** (0.250), **D7** (0.310), **D8** (0.234).
- Ninguna positiva en la zona gris; ninguna estructura en zona gris.
- Entropías 20261082 (puntos) y 20261083 (auxiliar), spawn_key (322, k, N, s).
