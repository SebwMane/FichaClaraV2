# F6 confirmatoria: preregistro congelado (2026-10-14)

- **Firma:** el usuario aprobó y firmó la propuesta de `OMEGA_F6_EXPLORACION_RESULTADOS.md` §7 («Queda aprobada y firmada»). Marca de compuerta: `R12-1 | **Firmada**`.
- **Intento 1 de 2** del presupuesto de la línea del perfil local (Rev. 11 §6). Todos los criterios con voto anteriores en 4D fallaron (F3, F5 brazos A y C, F5′ y F7).
- **Este documento se congela antes de abrir `CONFIRMACION_CIEGA.md`.** Nada de lo que sigue puede cambiar después de abrirlo.

| # | Decisión |
|---|---|
| R12-1 | **Firmada**: confirmación de F6 con el intercepto de la rama 1 como criterio único. |

## 1. Instrumento (idéntico al de la exploración, `tools/f6_exploracion.py`)

- Por realización: matriz causal, B, m_max; para m ∈ {10, 20, 40, 80, 160} con m ≤ 0.1·m_max, hasta 300 intervalos con exactamente m interiores elegidos al azar; si hay menos de 30, la celda (realización, m) no es medible.
- **r** = fracción de pares relacionados entre los m interiores (sin puntas). Δr = media de r − r_plano(m), con r_plano(m) de una referencia exacta nueva: 4000 diamantes planos 4D de m puntos más las puntas.
- **Una «estructura medida»** = 10 realizaciones con N = 4000 y 10 con N = 8000. Una celda (m, N) entra en el ajuste si al menos la mitad de sus realizaciones son medibles. Su valor es la media de Δr sobre las realizaciones medibles.
- **Ajuste:** Δr = b + A·x, con x = √(m/N), por mínimos cuadrados sobre las celdas. Si hay menos de dos valores distintos de x, la estructura es **«no medible»** y **no** se clasifica como espacio-tiempo (regla de R10-4).
- Sin dilución ni N = 16000 en la confirmación.

## 2. Calibración (antes de abrir el sello; entropía 20261080)

- **15 grupos de calibración**, cada uno una «estructura medida» completa:
  - 10 diamantes planos 4D;
  - 3 de Sitter 4D (η ∈ [−1.5, −0.5]);
  - 2 FRW de polvo 4D (a ∝ η²).
- En total, 300 realizaciones.
- **b_plano** = mediana de b en los 10 grupos planos.
- **τ_b = 1.5 × max sobre los 15 grupos de |b − b_plano|.** Con 15 grupos el percentil 95 es casi el máximo, así que se usa el máximo y se declara.

## 3. Clasificación y veredicto

- **«Espacio-tiempo 4D»** si |b − b_plano| ≤ τ_b.
- **PASA** si todas las positivas quedan clasificadas como espacio-tiempo y **ninguna** negativa con voto. Las estructuras en zona gris (ε ≤ 0.20, E-R11-1 y Rev. 11 §8.1) se informan sin voto. **FALLA** en cualquier otro caso.

## 4. Estructuras de prueba (entropía de las ya vistas: 20261081; las selladas usan la asignación de su archivo, 20261076 a 20261079)

**Selladas** (`CONFIRMACION_CIEGA.md`, SHA-256 `b8439d85fcf2ea3c6ea76b92cce8caced98483d85dac64bcb7ca43f3fbc91c79`): 3 positivas, 4 negativas y 1 en zona gris, con el papel que les asignó el guardián. Se verifica el hash antes de usarlas; si no coincide, se anula la confirmación.

**Ya vistas, con semillas nuevas:**

| Estructura | Papel |
|---|---|
| diamante 4D plano, de Sitter, de Sitter fuerte, anti de Sitter, FRW radiación, FRW polvo | positivas |
| cono cúbico, cono octaédrico, cono p = 6, orden producto 4, percolación 0.05 | negativas |
| conos p = 1.5 y p = 3 | zona gris (sin voto) |

## 5. Informado sin voto

- **Pendiente A con su signo, como predicción ciega:** + en de Sitter (más grande en el fuerte), − en anti de Sitter, ≈ 0 en FRW de radiación. En Kasner de vacío (R = 0 y R_ab = 0) la teoría predice que el término de orden τ² desaparece: **A compatible con 0**.
- **Rama 2:** D_loc por m (forma y máximo) y max-T con la nula de 30 diamantes planos por N (entropía 20261080, otra rama del spawn_key).
- **Intervalos de confianza** de b por bootstrap de realizaciones completas, y el margen |b − b_plano|/τ_b por estructura.

## 6. Convenciones

- **EXPERIMENT_ID 321.**
- **spawn_key:** (321, familia, índice, N, realización), con familias 0 = referencia, 1 = nula, 2 = calibración, 3 = vistas y 4 = bootstrap. Las selladas usan la asignación de su archivo.
- La corrida oficial se hace sobre un árbol git limpio, después de la auditoría de la herramienta.
