# Ω — P1-D.3: calibración del instrumento geométrico a N grande — Prerregistro

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-p1d3`.
- **Base congelada:** `claude/omega-p1d2-congelado`.
- **Estado:** este documento se commitea **antes** de escribir el código del panel y de generar ningún grafo.

## 0. Qué pregunta responde (y cuál no)

**Pregunta:** ¿qué rasgos medibles distinguen sistemáticamente una geometría de un grafo que solo parece geométrico, a escalas donde los efectos de tamaño finito ya no dominan (N ≈ 2·10⁴ y 5·10⁴)?

Es una prueba del **instrumento**, no de Ω. No se simula ninguna dinámica Ω. Ninguna salida de este experimento puede modificar el certificado Ω-1.1 ni sus umbrales, ni convertirse en término de energía.

**Lo que no se hace (y por qué):**
- **No se re-ejecuta F1 ni se escala C0.** F1 y F1b ya se completaron (F1-NEGATIVO 10/10; 0/24 certificados con N = 729). La dimensión de los estados C0 baja con N: D_s pasa de 2.6–3.3 con N = 343 a ≈ 2.2 con N = 729. En P1-D y P1-D.2, C0 no homogeneiza en ningún grafo evaluable. Escalarlo sería un intento de rescate, que el Consejo ya desaconsejó tras la Fase 2.
- **No se repite Θ > 0.** C0-Θ ya dio FRÁGIL-1/N.
- **No hay funcional nueva.**

## 1. Reglas generales

- Todo grafo se construye **disperso** (`scipy.sparse`); nunca matrices densas N × N.
- **RNG:** PCG64 con clave (20261009, familia, semilla), mediante `omega.c0.references.rng_from_key`.
- Se trabaja sobre la componente gigante de la adyacencia binaria simétrica.
- Un solo hilo BLAS.

## 2. Observables (nuevo módulo `omega/diagnostics/sampled_growth.py`)

### 2.1 Bolas por BFS muestreado

- Se eligen S = 400 nodos fuente uniformes sin reemplazo de la componente gigante (con su propio RNG derivado de la clave).
- Distancias en saltos desde cada fuente con `scipy.sparse.csgraph.shortest_path(..., unweighted=True, indices=fuentes)`.
- |B_i(r)| = #{j : d(i, j) ≤ r} para cada fuente i, r = 1 … r_cap, donde r_cap es el primer r con mediana de |B_r| ≥ N_g/2 (máximo 200).
- m(r) = media sobre las fuentes; CV(r) = desviación típica / media; mediana(r).
- **Ventana:** W = {r ≥ 2 : mediana(r) ≤ N_g/4}.

### 2.2 Nivel I: coherencia multiescala (definición de P1-D sin cambios)

- ρ = CV(max W)/CV(2).
- **Estados:** `CV_CERO` (CV(2) < 1e-9), `SIN_VENTANA` (|W| < 2), `HOMOGENEIZA` (ρ < 1), `NO_HOMOGENEIZA`.
- **Pasa:** HOMOGENEIZA o CV_CERO.

### 2.3 Nivel II: estabilidad dimensional (nuevo criterio de convergencia, sustituye a la meseta de P1-D.2)

- **D_B2(r)** = ln(m(r+2)/m(r)) / ln((r+2)/r), para r ≥ 2 con r + 2 ≤ max W. Es la derivada logarítmica suavizada a dos pasos, que elimina la oscilación de período 2 vista en P1-D.2.
- Sea d_1, …, d_n la sucesión de D_B2 y Δ_i = d_{i+1} − d_i.
- **Estados:**
  - `SIN_VENTANA_D`: n < 4;
  - `CONVERGE`: |Δ_{n−1}| ≤ 0.05 **y** |Δ_{n−1}| ≤ |Δ_{n−2}| + 0.02 **y** |Δ_{n−2}| ≤ |Δ_{n−3}| + 0.02. Las diferencias son pequeñas al final y no crecen;
  - `CRUCE`: Δ_{n−1} > 0.05, es decir, D sigue subiendo;
  - `NO_CONVERGE`: cualquier otro caso.
- **D_conv** = d_n, el último valor.
- **Clase dimensional:** round(D_conv) si |D_conv − round(D_conv)| ≤ 0.25 y round(D_conv) ≥ 1; en otro caso, `NO_ENTERA`.
- **Pasa:** CONVERGE.

**Procedencia declarada.** Los umbrales 0.05 y 0.02 se eligieron después de ver las trazas de D_B2 de P1-D.2:
- geometrías con N = 4096: diferencias finales de 0.00 a 0.03 y decrecientes;
- WS β = 0.003: ≈ 0.2 sostenido.

Por eso aquí se aplican a grafos **nuevos** y con N mayor. El veredicto se informa también con umbral final 0.03 y 0.08 como sensibilidad.

### 2.4 Dimensión espectral (descriptiva; no entra en ningún criterio)

- Paseo perezoso P = ½(I + D⁻¹A).
- Para 64 fuentes muestreadas se propaga p_t = P^t e_i con productos dispersos para t = 1 … T, donde T = 4·(max W)², con un máximo de 4000.
- Se promedia la probabilidad de retorno p_t(i).
- D_s(t) = −2 · d ln p / d ln t, por diferencias en una rejilla logarítmica de t (cada factor 1.25).
- Se informan la serie y su mediana en la mitad superior de t.

### 2.5 Categoría de informe (pedida por el especialista en geometría)

| Categoría | Condición |
|---|---|
| `NO_GEOMETRICO` | No pasa I o no pasa II |
| `GEOMETRIA_GRUESA(d)` | Pasa I y II, con clase d |
| `GEOMETRIA_GRUESA(no entera)` | Pasa I y II, clase NO_ENTERA |
| `GEOMETRIA_VARIEDAD(d)` | Gruesa con clase d ≥ 3 **y** su análogo N ≈ 729 pasa el certificado Ω-1.1 (§5) |
| `VARIEDAD_NO_EVALUABLE` | Gruesa con clase 1 o 2: el certificado no está calibrado para esas clases, porque rechaza también RGG2 y el anillo (P1-D.2 §5) |

## 3. Panel (N ≈ 2·10⁴; semillas 0–2 en familias aleatorias, 1 en deterministas)

| Clase | Familia | Construcción |
|---|---|---|
| **G-hom** (geometrías homogéneas) | `RGG2_k12`, `RGG3_k12`, `RGG3_k8` | Toro [0,1)^d periódico, N = 20 000, radio con V_d r^d (N − 1) = k; vecinos con `cKDTree(boxsize=1)` |
| G-hom | `RGG_S2_k12` | N = 20 000 puntos uniformes en la esfera unidad; arista si la distancia de cuerda < r, con r fijado para que el grado medio esperado sea 12 (casquete: 2π(1 − cos θ)·N/(4π) = 12, cuerda r = 2 sin(θ/2)) |
| G-hom | `anillo_k12` (N = 20 000), `cuadrado_141` (141² = 19 881), `T3_27` (27³ = 19 683) | Retículos |
| **G-inh** (geometrías reales inhomogéneas) | `RGG3_caja_k12` | RGG3 en el cubo [0,1)³ **sin** periodicidad, N = 20 000, mismo radio que el toro |
| G-inh | `RGG2_caja_k12` | Igual en 2D |
| G-inh | `RGG3_gradiente` | Toro 3D; la coordenada x se muestrea con densidad ∝ 1 + 2x (normalizada); radio fijo con grado medio global ≈ 12 |
| **R** (geometría gruesa, localmente no variedad) | `caveman_K{4,5,6,8,10,12,16}` | Anillo de cliques (misma regla que `connected_caveman`), N = 20 000 (los sobrantes van a la última clique) |
| R | `cliques3D_K6_L15` (N = 20 250), `cliques3D_K8_L13` (N = 17 576) | Toro 3D de sitios con clique K_m; la dirección j enlaza el miembro j mod m con el miembro (j xor 1) mod m del vecino |
| R | `cliques2D_K6_L58` (N = 20 184) | Igual con 4 direcciones en un toro 2D |
| **L** (locales no geométricas) | `retazos_b2`, `retazos_b4`, `retazos_b8` | RGG3 k12 toroidal; bloques 2³, 4³, 8³; los extremos de las aristas entre bloques distintos se reemparejan al azar (P1-D) |
| L | `RGG3_atajos_0.1`, `RGG3_atajos_1` | RGG3 k12; se reconecta un extremo del 0.1 % / 1 % de las aristas a un nodo al azar |
| L | `WS_b{0.001, 0.002, 0.005, 0.01}` | Anillo k12; se reconecta un extremo con probabilidad β |
| **NL** | `ER_k12`, `RR_k12` | Aleatorios |

**Subpanel de escala, N ≈ 5·10⁴, semillas 0 y 1:** `RGG3_k12`, `RGG2_k12`, `retazos_b4`, `WS_b0.002`, `caveman_K8` y `T3_37` (37³ = 50 653). Sirve para ver la tendencia con N de D_conv y de ρ.

## 4. Predicciones del cerebro (congeladas)

1. **G-hom.** Pasan I y II. D_conv está a ±0.25 de la dimensión real, incluida la esfera S² (≈ 2). En el RGG3, D_conv con N = 5·10⁴ está más cerca de 3 que con N = 2·10⁴.
2. **G-inh: el Nivel I falla en geometrías reales inhomogéneas** (ρ ≥ 1) en la caja 3D y en el gradiente. Es el límite que dejó señalado P1-D.2: el Nivel I mide invariancia de traslación, no geometría. El Nivel II sí converge en ellas.
3. **R.** Pasan I y II con D_conv ≈ 1 (caveman), ≈ 3 (retículos de cliques 3D) y ≈ 2 (2D) → `GEOMETRIA_GRUESA`. El análogo 3D falla el certificado con F9 y no llega a `GEOMETRIA_VARIEDAD`.
4. **L.** Ninguna familia pasa I y II. Los más cercanos son atajos al 0.1 % y WS β = 0.001; espero CRUCE (D sigue subiendo) o NO_HOMOGENEIZA.
5. **NL.** SIN_VENTANA.
6. **D_s** (descriptiva) concuerda con D_conv en G-hom con un error ≤ 0.4. En R da valores distintos de D_conv, porque las cliques atrapan al paseante.

## 5. Nivel III (certificado Ω-1.1, sin cambios)

- Para cada familia R o G-inh con categoría `GEOMETRIA_GRUESA(d ≥ 3)` en ≥ 2/3 de sus grafos se construye un análogo con N ≈ 729 y se corre `certificate_report` con 3 claves:
  - `cliques3D_K6_L5` (N = 750);
  - `cliques3D_K8_L4` (N = 512);
  - `RGG3_caja` (N = 729);
  - `RGG3_gradiente` (N = 729).
- Referencia: `RGG3_k12` toroidal con N = 729.

## 6. Veredicto del instrumento

| Medida | Definición |
|---|---|
| Sensibilidad | Fracción de G-hom que pasa I y II con la clase correcta (±0.25) |
| Especificidad | Fracción de L ∪ NL que **no** pasa I y II |

| Veredicto | Condición |
|---|---|
| **INSTRUMENTO-VÁLIDO** | Sensibilidad ≥ 0.9, especificidad ≥ 0.9 y el Nivel III rechaza (código F9 o F4) todo análogo R de clase 3 |
| **INSTRUMENTO-INVÁLIDO** | Sensibilidad < 0.7 o especificidad < 0.7 |
| **INSTRUMENTO-PARCIAL** | Cualquier otro caso |

El veredicto se calcula con el umbral final de Nivel II = 0.05 y se informa también con 0.03 y 0.08. G-inh y R no entran en el veredicto: definen el **alcance** del instrumento y se informan por categoría.

## 7. Uso posterior (decisión ya tomada)

Si el instrumento es VÁLIDO o PARCIAL, se convierte en la batería obligatoria (puntos B y D del preflight A–H) para evaluar **cualquier** dinámica futura (funcional nueva o P3).
- Siempre informará las tres categorías: gruesa, variedad y clase.
- Nunca entrará en la acción.
- No se reabre C0.
