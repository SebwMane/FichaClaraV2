# Ω — Mapa final de exclusiones (juez A−, tras RC-3)

- **Fecha:** 2026-10-06.
- **Uso:** catálogo de clases degeneradas conocidas, con la prueba que excluye cada una, la evidencia y sus puntos ciegos.
- **Es de mundo abierto (O2):** una clase no catalogada puede pasar.

| Clase degenerada | Firma | Prueba RC-3 | Validada en | Punto ciego / margen |
|---|---|---|---|---|
| Clique / compacta | R acotado (δ = 0) | X1 | K200 / K1600 | — |
| Expansor (RR, ER) | Diámetro logarítmico | X1 (+ X2, X4) | RR k3, ER k4 | — |
| k-árbol, 2-árbol | Diámetro logarítmico | X1 | 4-árbol, 2-árbol | δ hasta 0.11 (margen 0.015) |
| Triangulación apilada (apoloniana) | Diámetro logarítmico | X1 | Apoloniana | **s0 con δ = 0.130: no excluida** (W3) |
| Preferencial (BA) | Diámetro logarítmico | X1 | BA m3 | — |
| Árbol (uniforme, binario) | Varios extremos | X2 | Prüfer, binario subdividido | — |
| Cactus | Varios extremos + incoherencia | X2, X4 | Cactus | — |
| Árbol × fibra | Incoherencia | X4 | Árbol binario × C₁₀ | X4-marginal |
| **Árbol de bloques grandes** | Varios extremos solo a escala ≫ bloque | — | Árbol de cubos 5³ | **No excluido a N ≤ 8·10⁴** (W2) |
| Mundo pequeño con cruce lento | Incoherencia | X4 | RGG3 + 0.1 %, WS 2D β 0.01 | **X4-marginal** (γ 0.43 / 0.30 frente a geometrías ≥ 0.56) |
| Retazos | Incoherencia | X4 | — | **SIN_VENTANA en N = 10⁴: no excluidos** |
| 1D (anillos, escaleras, tubos, cilindros) | Dos extremos / δ ≈ 1 | X3 | C × C₅, RGG1, escalera, tubo C × RR | — |

**Resolución del juez:**
- d_max = 7 por diseño (X1);
- comprobado sin exclusiones falsas en d = 2…6 con N / 8N = 10⁴ / 8·10⁴.

**Lo que A− no hace:**
- certificar geometría (eso es A+);
- seleccionar ni estimar d, aunque δ ≈ 1/d sea un subproducto;
- excluir clases fuera del catálogo.
