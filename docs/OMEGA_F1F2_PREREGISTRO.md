# F1 y F2: pruebas de la firma espectral (preregistro congelado, 2026-10-12)

- **Pedido del usuario:** ejecutar R5-4 y, en paralelo, someter la firma a F1 y F2 (Consejo Rev. 6 §6).
- **Instrumento:** autovalores positivos de iΔ = (i/2)(C − Cᵀ); pendiente β y rodilla con `knee_from_spectrum` (R-rodilla-2), sin cambios.

## F1: prueba ciega de discriminación (vota)

Ninguna de estas estructuras fue medida espectralmente antes:

| Estructura | Papel | Descripción |
|---|---|---|
| de Sitter 2D | positiva | sprinkling en coordenadas conformes (η, x), diamante centrado en η = −1 con semialtura 0.5, densidad ∝ 1/η² (el orden causal es el plano; la geometría es curva) |
| Cilindro 2D | positiva | espacio periódico de largo 1, tiempo [0, 0.4], relación por imagen mínima |
| 2-orden en 200 capas | negativa | `layered_2order` |
| 2-orden en 400 capas | negativa | ídem |
| Orden aleatorio de dimensión 3 | negativa | i ≺ j si sus 3 coordenadas uniformes crecen (orden producto, no lorentziano) |
| Orden aleatorio de dimensión 4 | negativa | ídem con 4 coordenadas |
| Percolación p = 0.02 | negativa | |
| Percolación p = 0.1 | negativa | |

- N ∈ {800, 1600}, 5 semillas por N (10 corridas por estructura), entropía 20261058.
- **Compatible** (ventanas C5-A, sin cambios): β mediano ∈ [−1.104, −0.904] y rodilla mediana ∈ [0.227, 0.906].
- **PASA:** las dos positivas compatibles y las seis negativas no compatibles. **FALLA:** cualquier otro caso. Se informa cada estructura.
- **Alcance:** las ventanas son 2D; F1 prueba si la firma reconoce «geometría bidimensional», incluida una curva y una con otra topología.

## F2: dependencia de la dimensión (vota)

- Sprinklings en diamantes de Minkowski, d ∈ {2, 3, 4}, N ∈ {400, 800, 1600, 2400}, 5 semillas, entropía 20261059.
- Hipótesis (ley de Weyl): β → −2/d.
- **PASA:** para cada d, |β(2400) + 2/d| < |β(400) + 2/d| (se acerca al crecer N) **y** |β(2400) + 2/d| ≤ 0.10.
- **FALLA:** cualquier d que no cumpla.
- **Expectativa honesta:** 2D probablemente pasa; 4D puede fallar la tolerancia (la exploración dio −0.64 en N = 800 frente a −0.5).
- Nota: el tramo de ajuste de la pendiente (rangos 2 a max(6, ⌈√N/4⌉)) es corto en N = 400; se declara.
