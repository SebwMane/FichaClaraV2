# R17-3: aceleración exacta del S_SJ original (E1′). Resultados técnicos (2026-10-19)

- **Mandato:** R17-3 firmada (Rev. 17 §6). Se acelera el **S original de E1′**, no el candidato de R17-1. Desarrollo técnico, sin evidencia.
- **Código:** `omega/causality/sj_fast.py`, `tests/test_sj_fast.py` y `tools/r17_3_bench.py`. **Datos:** `results/r17_3/`, con semillas SeedSequence(20261102, spawn_key = (327, ...)).

## 1. Exactitud: cumple el criterio

- **Comparación estado por estado:** 20 040 estados de cadenas E-a y E-b con N = 32, 64, 128, 256 y 512, incluidas las propuestas rechazadas.
  - |ΔS| máximo frente al cálculo completo: **8.3e-12** (el límite era 1e-9).
  - Diferencias en |R|: **0**.
- **Verificación independiente del orquestador:** con 24 órdenes de N = 48, 100 y 200 (2-órdenes y KR), |ΔS| máximo de 4.8e-13.
- **Pruebas:** 48 pasan (14 nuevas). Además, las cadenas rápidas con las mismas semillas reproducen exactamente los estados de los muestreadores de ES-0.
- **Modos cerca del umbral:** si un autovalor queda a menos de 1e-8 de un umbral de truncamiento, se usa el cálculo original. Nunca hizo falta: el margen mínimo observado fue 4.9e-7.

## 2. Cómo se logra (álgebra exacta, sin aproximaciones)

- iΔ = (i/2)K, con K = C − Cᵀ real, y (iΔ)² = KᵀK = M, que tiene entradas enteras.
- Se diagonaliza la matriz **real** M en lugar de la compleja iΔ. De ahí salen |μ| y los vectores que construyen exactamente iΔ truncado y W.
- Solo se forman los bloques de la región, no las matrices N × N completas.
- M y los conteos de intervalos se actualizan en O(N²) por movimiento.
- **No se adoptaron** las actualizaciones de bajo rango de la diagonalización: con estos N cuestan lo mismo o más que diagonalizar de nuevo, y acumulan error.

## 3. Ganancia

| N | ms por evaluación, completo → rápido | Aceleración |
|---|---|---|
| 64 | 1.1 → 0.8 | 1.4 veces |
| 128 | 4.6 → 2.1 | 2.1 veces |
| 256 | 22 → 7.3 | 3.1 veces |
| 512 | 190 → 34 | 3.5 a 5.7 veces |

**Segundos por barrido:**

| Espacio | N = 64 | N = 128 | N = 256 |
|---|---|---|---|
| E-a (2-órdenes) | 2.2 → 1.6 | 37 → 17 | 730 → 240 |
| E-b (órdenes generales) | 0.8 → 0.55 | 7.2 → 2.5 | unos 86 → 28 (estimado) |

## 4. Límite

Lo que sigue siendo O(N³) es la diagonalización de M, que es el piso de cualquier método que diagonalice la matriz completa. La aceleración es de factor 2 a 5, no de órdenes de magnitud.

Con N de miles, donde el instrumento de E1′ tiene sentido, el costo por propuesta sigue en segundos, y las cadenas necesitan millones de propuestas. **R17-3 no cierra por sí sola la brecha entre validez del instrumento y mezcla** (Rev. 17 §1, tema 3). Sí permite trabajar con N = 128 a 256 en desarrollo técnico.
