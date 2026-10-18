# Ω — E6-S v1.1: panel cerrado por transporte (prerregistro)

Estado: **PRERREGISTRO**, escrito después de ver el fallo de v1 (`docs/OMEGA_E6_RESULTADOS.md`) y antes de cualquier código de v1.1.

Este diseño **no es ciego** respecto del fallo de v1. Para compensarlo:
- se añaden controles frescos, que nunca se han evaluado;
- los codificadores antiguos tienen que seguir detectándose sobre el panel ampliado;
- un pase de v1.1 se registra como validación **más débil que una ciega**.

## 1. Qué cambia y qué no

**Sin cambios:**
- Umbrales: R ≥ 0.9, Q ≤ 0.5, CIEGO si max − min ≤ 0.10.
- El orden de clases y la definición de DIAL.
- Las reglas y rejillas de v1.
- El panel base §3.1 y las semillas (`MASTER_E6 = 20261020`).

**Cambio único: cierre del panel por transporte.** Se aplican tres operaciones genéricas que conservan la dimensión de crecimiento a cada Z^d (d = 1…6), y las geometrías resultantes se añaden a su d.

| Operación | Definición | Grado | Lados (d = 1…6) |
|---|---|---|---|
| **Truncación T(Z^d)** | Cada vértice se sustituye por un ciclo de 2d nodos, uno por arista incidente, en orden cíclico (+e₁, …, +e_d, −e₁, …, −e_d). Cada arista de Z^d une los nodos correspondientes de sus dos extremos | 3 | d = 1 no se aplica (gadget degenerado); d = 2: 50; 3: 12; 4: 6; 5: 5; 6: 5 |
| **Grafo de líneas L(Z^d)** | Nodos = aristas de Z^d; adyacentes si comparten extremo | 4d − 2 | 10000; 100; 20; 10; 6; 5 |
| **Producto Z^d □ K₂** | Dos copias unidas por peldaños | 2d + 1 | 5000; 70; 16; 8; 6; 5 |

Muestreo:
- Si una geometría tiene N > 20000, las reglas de vértice se evalúan sobre 3000 vértices muestreados con `rng_from_key((MASTER_E6, id, 0, 8))`.
- Si N ≤ 20000, sobre todos los vértices.
- Las reglas de arista mantienen el muestreo de v1.

## 2. Controles frescos (nunca evaluados antes de este documento)

| Regla | Tipo | Reposo | Rejilla | Predicción |
|---|---|---|---|---|
| C-5CIC | arista | la arista está en un 5-ciclo | — | CIEGO (circulante en d = 1; RGG en todo d) |
| C-PAR | vértice | grado par | — | CIEGO (Z^d en todo d) |
| P0-S2(c) | vértice | \|S₂(v)\| = c (esfera de radio 2). En Z^d, \|S₂\| = 2d² | c ∈ {8, 18, 32, 50} (d = 2…5) | SELECTOR en {2}, {3}, {4}, {5}; familia **DIAL** |

## 3. Condiciones de validez de v1.1 (congeladas)

- **V1′.** P0-DEG, P0-SQV y P0-S2 son DIAL.
- **V2′.** Clases y conjuntos de reposo requeridos:
  - P0-SQV(12) SELECTOR con 3 ∈ R;
  - P0-SQV(24) SELECTOR con 4 ∈ R;
  - P0-S2(18) SELECTOR con 3 ∈ R;
  - P0-S2(32) SELECTOR con 4 ∈ R.
- **V3′.** C-TRIV, C-SQ, C-TRI, P0-LAT(2), C-5CIC y C-PAR no son SELECTOR, y la familia P0-LAT no es DIAL.

Veredicto: **E6-S v1.1 VÁLIDO** si se cumplen V1′–V3′; si no, INVÁLIDO. No hay v1.2 sin decisión del Consejo.

## 4. Predicciones del cerebro (no deciden)

| Regla | Predicción |
|---|---|
| P0-LAT(3) | MONÓTONO, R = {2, …, 6} (panal y T(Z^d) para d ≥ 3) |
| P0-DEG(3) | MONÓTONO, R = {2, …, 6} |
| Resto de reglas de v1 | Misma clase que en v1 |
| Mayor riesgo | Que L(Z^d) o Z^d □ K₂ coincidan con algún valor codificador (c₄ o \|S₂\|) en otra d y rompan V2′ |

Comprobación previa de colisiones en Z^d □ K₂:
- c₄ = 2d(d − 1) + 2d = 2d² no coincide con ningún c de P0-SQV.
- \|S₂\| = 2d² + 2d no coincide con ningún c de P0-S2.
- L(Z^d) no se ha calculado de antemano.
