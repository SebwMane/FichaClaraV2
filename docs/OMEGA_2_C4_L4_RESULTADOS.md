# Ω-C en 4D: resultados de C4 (instrumento) y L4 (rama lorentziana)

- **Fecha:** 2026-10-12. **Commit:** `22de481` (árbol limpio). **Especificación:** `OMEGA_2_PROPUESTA_CAUSAL.md` §19. **Datos:** `results/c4_l4/`.

## 1. C4: acción 4D · ÉXITO

| N | ε | media ± SE | referencia independiente |
|---|---|---|---|
| 128 | 1 | 41.3 ± 13.8 | 45.2 ± 0.9 |
| 128 | 0.25 | 68.8 ± 1.0 | 70.0 ± 0.2 |
| 128 | 0.05 | 102.09 ± 0.22 | 102.10 ± 0.05 |
| 512 | 1 | 231 ± 58 | 83 ± 12 |
| 512 | 0.25 | 174.5 ± 4.4 | 180.5 ± 3.5 |
| 512 | 0.05 | 298.71 ± 0.68 | 298.54 ± 0.77 |

Las seis celdas dentro de 3σ. Con ε pequeño la coincidencia es muy precisa; con ε = 1 la prueba tiene poca potencia por las fluctuaciones de la acción sin suavizar. La referencia usa solo coordenadas (método verificado contra el exacto 2D en la auditoría). La media crece con N, como corresponde a los términos de borde en 4D.

## 2. L4: compatibilidad de la rama lorentziana · **INCOMPATIBLE (esquemas probados)**

D₄ = sd(S₄,ε)/ε; se necesitaba D₄ ≤ 0.477 sin crecer.

| esquema | N = 128 | 512 | 2048 | η |
|---|---|---|---|---|
| ε = 1 | 204 | 760 | 2 766 | 0.97 ± 0.03 |
| ε = 0.25 | 51 | 208 | 639 | 0.92 ± 0.04 |
| ε = 0.05 | 54 | 191 | 650 | 0.93 ± 0.04 |
| ε = 16/N | 42 | 226 | 1 488 | 1.31 ± 0.03 |
| ε = 64/N | 79 | 201 | 695 | 0.81 ± 0.05 |

## 3. Lectura

- **Nivel 1:** en diamantes causales 4D, la dispersión de la acción geométrica en las unidades relevantes para la ventana de supresión de las capas está entre 600 y 2 800 en N = 2048 y crece aproximadamente como N^0.8 a N^1.3, en los siete esquemas.
- **Nivel 2:** igual que en 2D, en ningún esquema probado la geometría puede ser coherente donde las capas se cancelan.

## 4. Límite importante (no cambia el veredicto)

Todas las regiones medidas son **diamantes causales, con borde**. La literatura asocia a la acción BD términos de borde importantes y fluctuaciones ligadas a ellos. No se puede descartar que buena parte de la dispersión venga del borde y que una región sin borde espacial (por ejemplo, un toro espacial) se comporte distinto. Ese es el siguiente experimento natural: **L4b**, separar la dispersión de volumen de la de borde.

## 5. Estado

| Rama | 2D | 4D |
|---|---|---|
| Lorentziana | Cerrada (L-ε) | Incompatible en diamantes; pendiente L4b (sin borde) |
| Euclídea | Estable localmente; equilibrio no resuelto | No iniciada (requiere Metropolis sobre órdenes generales) |

> **Nota de unidades (2026-10-19, auditoría de √ε, `docs/omega_c/AUDITORIA_SQRT_EPS.md`):** la acción 4D suavizada de esta corrida está en unidades del repositorio de entonces, sin el factor √ε (S_ε/√ε). D₄ y D_int son invariantes ante ese factor, porque numerador y denominador cambian igual, así que los veredictos no cambian. Desde el commit de la corrección, `bd_action_4d_smeared` incluye √ε.
