# Ω, F6 en 4D: estructuras a ciegas para la confirmación (sellado)

- **Autor:** guardián de reserva (subagente que no participa en el diseño). Fecha: 2026-10-14.
- **Uso:** este archivo se abre solo después de congelar el criterio de F6 (Rev. 11 §6.2). Antes de eso solo se publica su SHA-256.
- **Referencias:** Rev. 11 §1, §2, §5, §6, §8.1 (ε₀ = 0.20); `consejo_rev11/sobreajuste.md` §4; `OMEGA_F6_EXPLORACION.md` (estructuras vistas: ninguna se repite aquí); handoff Kasner.
- **Convenciones:** dimensión 1 + 3; T = `DIAMOND_HALF_HEIGHT` = 0.5; columna 0 = t; Δ = y − x. «Diamante estándar» = `sprinkle_diamond(n, 4, rng)`, región {|t| + ‖x‖₂ ≤ T}. Para conos: región {|t| ≤ T − ‖x‖_cono}, puntos uniformes por rechazo en [−T, T]⁴ (todas las normas usadas cumplen ‖x‖ ≥ ‖x‖_∞, así que la caja contiene la región), relación y ≻ x ⇔ Δt > 0 y Δt ≥ ‖Δx‖_cono.
- **Deformación:** ε = max_{‖u‖₂=1} |‖u‖_cono − 1| (definición de E-R11-1). Valores exactos cuando hay fórmula cerrada; si no, numéricos sobre 2·10⁶ direcciones al azar más los 26 ejes y diagonales.
- **Papeles:** positiva (debe clasificarse como espacio-tiempo de dimensión 4), negativa (debe rechazarse), zona gris (ε ≤ 0.20: se informa, **sin voto**).

## Entropías y semillas

- **20261076:** generación de puntos de las estructuras de prueba. `SeedSequence(20261076, spawn_key=(EXPERIMENT_ID, k, N, s))`, con k = índice de la estructura en la tabla de abajo (0..7), N el tamaño y s la realización.
- **20261077:** aleatoriedad auxiliar por realización (modos de la densidad de C2, elección de links borrados de C5), misma spawn_key que arriba. Se consume ANTES de generar el resto de la realización solo en lo indicado.
- **20261078:** referencia plana exacta por m y calibración de la confirmación (diamantes planos 4D).
- **20261079:** distribución nula, remuestreo (bootstrap) y cualquier otro azar del análisis confirmatorio.
- Ninguna de las cuatro se usó antes ni puede usarse en otra parte. N, número de realizaciones, rejilla de m y regla de voto son los que fije el preregistro congelado; estas estructuras se aplican a cualquier N.

## Lista (orden = índice k)

| k | Nombre | Papel | d verdadera | ε |
|---|---|---|---|---|
| 0 | C0 Kasner de vacío, p = (−6/19, 10/19, 15/19), t ∈ [1.0, 1.9] | positiva | 4 | — (no conformemente plano) |
| 1 | C1 cono elipsoidal rotado en región diamante de coordenadas | positiva | 4 | ε crudo 0.4286; ε afín = 0 |
| 2 | C2 Minkowski conformemente reescalado, densidad suave aleatoria | positiva | 4 | 0 |
| 3 | C3 cono de sección cilíndrica | negativa | — | 0.2929 |
| 4 | C4 cono «esfera truncada por cubo», s = 0.79 | negativa | — | 0.2658 |
| 5 | C5 Minkowski con 13% de links borrados | negativa | — | — (impostor discreto) |
| 6 | C6 cono con p distinto por eje, p = (1, 8, 32) | negativa | — | 0.3144 |
| 7 | C7 cono anidado ‖(‖(x,y)‖₂, z)‖₄ | zona gris, sin voto | — | 0.1591 |

Resumen por papel: **3 positivas, 4 negativas, 1 zona gris.** Ninguna positiva tiene ε en (0, 0.20].

## C0. Kasner de vacío (positiva, no conformemente plana)

- **Métrica:** ds² = −dt² + Σᵢ t^{2pᵢ} dxᵢ², p = (p₁, p₂, p₃) = (−6/19, 10/19, 15/19) (familia de Lifshitz-Khalatnikov con u = 3/2).
  - Σp = (−6 + 10 + 15)/19 = 1; Σp² = (36 + 100 + 225)/361 = 361/361 = 1. Vacío exacto.
- **Región:** intervalo causal J⁺(a) ∩ J⁻(b), a = (t_lo, 0, 0, 0), b = (t_hi, 0, 0, 0), **t_lo = 1.0, t_hi = 1.9**.
- **Curvatura:** Kretschmann K = −16 p₁p₂p₃ / t⁴ = (14400/6859)/t⁴ ≈ 2.099/t⁴; longitud de curvatura ℓ = K^{−1/4} ≈ 0.832·t (0.83 a 1.58 en la ventana). Altura propia τ = 0.9; τ/ℓ(t = 1.45) ≈ 0.75 (moderada: de Sitter estándar ya visto tiene τ/ℓ ≈ 2.4). Densidad en coordenadas ∝ t: variación ×1.9 en todo el intervalo (dentro del límite S5 de ×2).
- **Generador:** `omega.causality.kasner.sprinkle_kasner_interval(n, p, 1.0, 1.9, rng, k=800)` y `kasner_causal_matrix(coords, p, k=800)` (K = 800 direcciones simétricas). Validar p con `validate_kasner_exponents`.
- **Por qué es positiva:** solución de vacío de Einstein, R = 0 y Ricci = 0, Weyl ≠ 0; la relación causal NO es la de Minkowski (prueba lo que la curvatura conformemente plana no prueba). Predicción H_loc: sin término τ² en vacío; exponente 4/d = 1.

## C1. Cono elipsoidal rotado (positiva, invariancia lineal)

- **Relación:** y ≻ x ⇔ Δt > 0 y Δt ≥ ‖A Δx‖₂, con A = diag(1/c₁, 1/c₂, 1/c₃)·Rᵀ, **c = (0.70, 1.00, 1.45)** (velocidad de la luz por eje propio) y **R = R_z(0.6)·R_y(0.9)·R_x(0.3)** (ángulos en radianes, rotaciones activas estándar).
- **Región:** diamante estándar de coordenadas {|t| + ‖x‖₂ ≤ 0.5} (`sprinkle_diamond(n, 4, rng)`), uniforme en coordenadas = uniforme en volumen propio (jacobiano constante).
- **Por qué es positiva:** con x' = A x es Minkowski exacto; la región, vista en x', es un convexo de Minkowski que NO es un diamante causal (los intervalos del borde quedan recortados, como en la caja). ε crudo de la norma espacial = max|1/cᵢ − 1| = 3/7 ≈ 0.4286 (fuera de la zona gris); ε afín-invariante (distancia al elipsoide más cercano) = 0.

## C2. Minkowski conformemente reescalado con densidad suave aleatoria (positiva)

- **Densidad en coordenadas:** f(x) = exp(φ(x)), φ(x) = Σ_{j=1}^{5} 0.07·cos(2π kⱼ·x + θⱼ), x = (t, x, y, z).
  - Por realización, con el generador de la entropía 20261077: kⱼ ∈ ℝ⁴ con cada componente uniforme en [−1.5, 1.5] y θⱼ uniforme en [0, 2π), sacados en el orden j = 1..5, primero las 4 componentes de kⱼ y luego θⱼ.
  - Rango: f_max/f_min ≤ e^{0.7} ≈ 2.01 en toda la región.
- **Muestreo (entropía 20261076):** candidatos uniformes del diamante estándar, aceptados con probabilidad exp(φ(x) − 0.35); se toman los primeros n aceptados.
- **Relación:** Minkowski (`causal_matrix`).
- **Por qué es positiva:** es el sprinkling uniforme en volumen propio de la métrica g = e^{φ/2} η (Ω⁴ = e^{φ}), conformemente plana y suave, con curvatura de escala dada por los modos.

## C3. Cono de sección cilíndrica (negativa, cerca de la frontera por fuera)

- **Norma:** ‖x‖ = max(√(x² + y²), |z|).
- **ε = 1 − 1/√2 ≈ 0.2929** (mínimo en las direcciones a 45° entre el plano xy y el eje z; máximo 1).
- Región y relación según las convenciones de conos.

## C4. Cono «esfera truncada por cubo» (negativa, cerca de la frontera por fuera)

- **Norma:** ‖x‖ = max(‖x‖₂, ‖x‖_∞ / s), **s = 0.79** (bola unidad = bola euclídea ∩ cubo de semilado 0.79: caras planas y casquetes redondos).
- **ε = 1/s − 1 = 21/79 ≈ 0.2658** (máximo en los ejes; mínimo 0 en direcciones donde manda la esfera).
- Región y relación según las convenciones de conos.

## C5. Minkowski con 13% de links borrados (negativa, impostor adversarial)

- **Construcción por realización:**
  1. n puntos en el diamante estándar (entropía 20261076) y C = matriz causal de Minkowski.
  2. Links: L = {(i, j) : C[i, j] y no existe k con C[i, k] y C[k, j]} (calculado desde C, sin coordenadas).
  3. Con el generador de la entropía 20261077: elegir **exactamente round(0.13·|L|)** pares de L, uniformes sin reemplazo (por ejemplo, `rng.choice(|L|, size, replace=False)` sobre L listado en orden lexicográfico de (i, j)).
  4. C' = C con esos pares puestos en False, todos a la vez. La prueba usa C'.
- **Por qué sigue siendo orden parcial:** C' ⊂ C hereda irreflexividad y antisimetría. Transitividad: si a ≺' b ≺' c, entonces a ≺ b ≺ c en C, así que (a, c) ∈ C y tiene un intermedio (b), luego (a, c) no era link y no fue borrado: (a, c) ∈ C'. Por la misma razón el cierre transitivo de C' es C' mismo: los pares borrados no son implicados por ninguna cadena (no tienen intermedios), así que el cierre no restaura nada.
- **Verificación obligatoria:** `is_partial_order(C')` y cierre(C') == C' en cada realización; informar |L|, pares borrados y la fracción de links de C' (puede aparecer algún link nuevo cuyo único intermedio dependía de un link borrado).
- **Por qué es negativa:** no es el sprinkling de ninguna variedad lorentziana; imita el hueco en k = 0 de de Sitter fuerte (falla 2 de F7) y no se «diluye» con la escala como la curvatura.

## C6. Cono con p distinto por eje (negativa)

- **Bola unidad:** B = {(x, y, z) : |x|¹ + |y|⁸ + |z|³² ≤ 1}, **p = (1, 8, 32)**, convexa porque cada pᵢ ≥ 1.
- **Norma (gauge de Minkowski):** ‖x‖ = el único λ > 0 con Σᵢ |xᵢ/λ|^{pᵢ} = 1 (bisección en λ; tolerancia relativa 1e-12).
- **ε ≈ 0.3144** (numérico; mín 0.6856, máx 1.0088), fuera de la zona gris y lejos de los conos ya vistos.
- Región y relación según las convenciones de conos.

## C7. Cono anidado casi redondo (zona gris, sin voto)

- **Norma:** ‖x‖ = ((x² + y²)² + z⁴)^{1/4}.
- **ε = 1 − 2^{−1/4} ≈ 0.1591** (≤ 0.20: zona gris). Se informa; no vota.

## Regla de lectura (recordatorio, no criterio)

El criterio, el umbral, N y la rejilla son los del preregistro congelado. Para el veredicto: C0, C1 y C2 deben aceptarse con d* = 4; C3, C4, C5 y C6 deben rechazarse; C7 no vota. Si una estructura queda «no medible» según la regla congelada, cuenta como fallo de su papel salvo que el preregistro diga otra cosa antes de abrir este archivo.
