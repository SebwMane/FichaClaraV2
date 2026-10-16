# Ω-C: resultados de C1 (acción suavizada) y L-3c-E (estabilidad euclídea)

- **Fecha:** 2026-10-12.
- **Commit de ejecución:** `157fe07` (las cuatro partes en el mismo commit limpio).
- **Especificación congelada:** `OMEGA_2_PROPUESTA_CAUSAL.md` §11 (preregistro `d3570d1`, enmienda L3cE-A1 en `157fe07`).
- **Datos:** `results/c1_l3ce/`.

## 1. C1: instrumento de la acción suavizada · ÉXITO

| N | ε | media S_ε ± SE | exacto |
|---|---|---|---|
| 64 | 0.5 | 0.46 ± 0.53 | 1.000 |
| 64 | 0.25 | 0.86 ± 0.13 | 1.000 |
| 64 | 0.1 | 0.984 ± 0.022 | 0.999 |
| 256 | 0.5 | −0.16 ± 1.23 | 1.000 |
| 256 | 0.25 | 0.83 ± 0.33 | 1.000 |
| 256 | 0.1 | 1.04 ± 0.07 | 1.000 |

Las seis celdas a ≤ 3·SE de su valor exacto. La media de la acción suavizada sobre el diamante plano también es ≈ 1, como la sin suavizar, con dispersión mucho menor. Nivel 1 sobre el instrumento.

## 2. L-3c-E: estabilidad de la geometría bajo e^{−βS} · «VENTANA QUE SE ANGOSTA»

Validación previa del clasificador (L3cE-A1): N = 32 marca 98% de 2-órdenes aleatorios y 0% de capas; N = 64, 99% y 0%. Válido.

| β | N = 32: fracción geométrica | N = 64: fracción geométrica |
|---|---|---|
| 0.025 | 1.00 | 1.00 |
| 0.035 | 1.00 | 1.00 |
| 0.050 | 1.00 | 0.65 |
| 0.071 | 1.00 | **0.20** ← β_c(64) |
| 0.100 | 0.70 | 0.00 |
| 0.141 | **0.25** ← β_c(32) | 0.00 |
| ≥ 0.200 | 0.00 | 0.00 |

**Decisión preregistrada: β_c(64) = 0.071 < β_c(32) = 0.141, «ventana que se angosta».** Era la expectativa (nivel 4) del argumento de conteo.

**Observación (nivel 1, dos puntos, resolución de un paso de √2 en la malla):** β_c · N ≈ 4.5 en ambos tamaños. La ventana geométrica se comporta como si β_c ∝ 1/N. El argumento de conteo predecía ln N / N; con solo dos N no se distinguen.

**Lectura preregistrada.** Con β fijo, la geometría no sobrevive al límite N → ∞ en la rama euclídea con la acción sin suavizar. **La rama no se cierra**: queda viva solo si el acoplamiento escala con el tamaño, β = b/N con b fijo. Eso cambia la hipótesis y requiere firma (R3-5).

## 3. Límites honestos

- A β ≥ 0.2 la aceptación cae a ≤ 3% y la acción media no es monótona: las cadenas se congelan sin llegar al mínimo (la bicapa daría −1984 en N = 64). La fase no geométrica está mal muestreada; solo se usó para decidir «no geométrico».
- Partir siempre de geometría mide estabilidad. Si hay histéresis, β_c está sobrestimado igual en ambos N (A-16).
- Una sola cadena por celda y dos tamaños. La relación β_c ∝ 1/N es una observación, no una ley.

## 4. Estado de la compuerta

| Rama | Estado |
|---|---|
| E | Viva bajo condición: requiere decidir R3-5 (acoplamiento β = b/N). |
| L | Bloqueada por método (problema de signo, §11.5). Requiere un enfoque analítico. |

R3-1 (adoptar Ω-C) sigue pendiente.
