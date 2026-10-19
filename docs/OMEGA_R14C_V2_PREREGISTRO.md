# R14-C v2: observables causales multivariables. Preregistro exploratorio (2026-10-19)

- **Mandato:** R15-1 firmada (v2 con la enmienda 4 de la Rev. 15) y R15-2 firmada (verificación previa de citas y fórmulas: `docs/consejo_rev15/VERIFICACION_R15_2.md`). **Estado:** PREPARADO, **no ejecutado**; ejecutarlo requiere una nueva orden explícita del usuario.
- **Reemplaza a** `OMEGA_R14C_PREREGISTRO.md` (v1, histórica, sin cambios). **Sin voto.**

## 0. Cambios respecto de v1

| Enmienda 4 (Rev. 15) | Cambio en v2 |
|---|---|
| Verificar las fórmulas pendientes | Hecho en R15-2. §1 da las fórmulas verificadas y la fuente de cada una. La cadena más larga **no** tiene fórmula cerrada en 4D: se usa referencia por Monte Carlo. |
| Fijar normalización y resumen por observable | §1, columna «normalización y resumen». |
| Multiplicidad | §3: se informan todas las comparaciones, con corrección de Holm para cualquier frase del tipo «separa». |
| Agregar grado de links y espectro del Laplaciano (Eichhorn et al. 2026) | §1, observables 6 y 7, con la definición del preprint y la ambigüedad resuelta antes de mirar datos. **Matiz verificado:** su evidencia es en 2D (sprinklings con frontera de caja, n = 2048) y solo el espectro separa las nueve clases. No se extrapola a 4D: aquí se explora. |
| No elegir después una combinación ganadora | §3. |

## 1. Observables (intrínsecos), por sub-intervalo con m ∈ {20, 40, 80, 160} (regla m ≤ 0.1·m_max), salvo los globales

| # | Observable | Referencia y fórmula (estado de verificación) | Normalización y resumen |
|---|---|---|---|
| 1 | Distribución de r por intervalo | ⟨r⟩ = Γ(d+1)Γ(d/2)/[2Γ(3d/2)] = **1/10 en 4D** (Surya 2019; confirmada por Monte Carlo: 0.0996 ± 0.0003 y 0.1002 ± 0.0003) | Media, desviación estándar y distancia de energía a la distribución de referencia (Monte Carlo del diamante plano con m puntos) |
| 2 | Tamaños de intervalo completos | Abundancias normalizadas S_m/n = N_{m−2}/R (definición de Eichhorn et al., ec. 16; equivale al perfil del repo) | Distancia de energía a la referencia por Monte Carlo |
| 3 | Abundancia de k-cadenas, k = 2..5 | ⟨C_k⟩ = n(n−1)···(n−k+1)·χ_k, con χ_2 = 1/20, χ_3 = 1/2100, χ_4 = 1/705600 y χ_5 = 1/582120000 en 4D. **Factorial descendente con n fijo, no n^k** (Meyer, vía Surya 2019 ec. 19; confirmada por Monte Carlo con |z| ≤ 1.4) | Cociente C_k/⟨C_k⟩; la dimensión que lo hace igual a 1 es un estimador de d |
| 4 | Cadena más larga L(m) | **Sin forma cerrada en 4D** (solo cotas: 1.856 ≤ m₄ ≤ 2.530 asintóticamente; para m ≤ 200 el valor real queda por debajo de la cota inferior). Referencia: Monte Carlo del diamante plano | L/L_ref(m) |
| 5 | Dimensión de punto medio | Instrumento C0-4 del repo | Diferencia con la dimensión de r |
| 6 | Grado de links (global) | d_j = ½Σ_i(L_ji + L_ij); P_j = (1/n)Σ_i δ_{d_i, j} (Eichhorn et al., ec. 14-15) | Distancia de energía de la distribución de grados a la referencia por Monte Carlo del mismo N |
| 7 | Espectro del Laplaciano del Hasse simetrizado (global) | L_S = L + Lᵀ; Δ = 𝟙 − 𝒟^(−1/2) L_S 𝒟^(−1/2). **Ambigüedad resuelta antes de mirar datos:** 𝒟 = diag(sumas de filas de L_S) (convención estándar, espectro en [0, 2], coherente con λ_max = 2 en órdenes bipartitos según el propio texto). Vértices aislados con pseudoinversa | λ₂, λ_max e histograma del espectro; distancia de energía a la referencia |
| 8 | Curvatura (informativa) | Pendiente A de F6 y la acción local de Benincasa-Dowker s(x) = 1 − N₀ + 9N₁ − 16N₂ + 8N₃ (BD 2010, verificada; ⟨s⟩ → (√6/8)ℓ²R). **Nota:** la versión suavizada 4D del repo no lleva el factor √ε; aquí no se usa la versión suavizada | Media por sub-intervalo frente a la referencia |

- **Calibración de referencias:** 200 diamantes planos 4D por m (sub-intervalos) y 60 sprinklings por N y por geometría (globales), con entropía propia.
- **Geometrías curvas:** su referencia es la propia geometría (R14-B v2 §3), no la plana.

## 2. Estructuras y unidad

- **Las de v1 §3** (todas vistas): geometrías conformemente planas, Kasner, lineales, conos suaves, poliédricos, defectos, producto y percolación.
- **Realizaciones y tamaño:** 10 por estructura y N ∈ {4000, 8000}.
- **Unidad:** la realización dentro de la estructura, y la **familia** para afirmaciones generales (R14-A §7). Los sub-intervalos de una realización no son réplicas.

## 3. Multiplicidad y salvaguardas

- **Se informan todas las combinaciones** (observable × m × N × familia) en una sola tabla, sin omitir las que no separan.
- **Cualquier afirmación «el observable X separa la familia F de las geometrías»** usa el tamaño de efecto por familia con IC por bootstrap de realizaciones y la corrección de Holm sobre todas las comparaciones de la tabla.
- **No se elige una combinación ganadora para reclasificar estas estructuras.** Si sale un candidato, se congela y se prueba con familias nuevas y un guardián nuevo.

## 4. Convenciones

- **EXPERIMENT_ID 324; entropía 20261089** (de v1, que nunca corrió).
- **Costo estimado:** de 4 a 5 horas (el espectro del Laplaciano global con N = 8000 es lo más caro: diagonalizar una matriz de 8000 × 8000 lleva unos minutos por realización; se mide en el humo).
- **Herramienta:** agente de código con handoff y auditoría; árbol limpio; tandas de menos de 2 horas.
