# F7: resultados (distancia de forma del perfil, 2026-10-14)

- **Especificación:** `OMEGA_CONSEJO_REV10_F6F7.md` §8 (R10-1 a R10-4 firmadas; preregistro en 67495fb). **Herramienta:** `tools/f7.py` (agente de código, auditoría §22). **Corrida oficial:** commit c1cae17, árbol limpio. **Datos:** `results/omega_f7/`.

## 1. Veredicto oficial: **FALLA**, por tres estructuras

- **de Sitter 4D fuerte (positiva, a ciegas): rechazado** (razón D_forma/δ_forma 1.55).
- **Cono p = 1.5 (negativa, a ciegas): aceptado** (0.58).
- **Cono p = 3 (negativa, a ciegas, vota por R10-2): aceptado** (0.41).

Por R10-3, **F6 queda pospuesta** y vuelve al Consejo.

## 2. Detalle (ET = clasificado como espacio-tiempo con el criterio que vota)

| Estructura | Papel | d* | Firma | D_forma/δ_forma | ET | ¿Correcto? | D/δ (F5′, informativo) |
|---|---|---|---|---|---|---|---|
| Diamante 4D | positiva | 4 | sí | 0.35 | sí | sí | 0.26 |
| Diamante 3D | positiva | 3 | sí | 0.43 | sí | sí | 0.28 |
| de Sitter 4D | positiva | 4 | sí | 0.72 | sí | sí | 0.97 |
| de Sitter 3D | positiva | 3 | sí | 0.52 | sí | sí | 0.62 |
| **de Sitter 4D fuerte** | positiva | 4 | sí | **1.55** | **no** | **no** | 2.03 |
| FRW de radiación 4D | positiva | 4 | sí | 0.49 | sí | sí | 0.55 |
| Cono cúbico (p = ∞) | negativa | 4 | sí | 3.05 | no | sí | 1.13 |
| Cono octaédrico (p = 1) | negativa | 4 | sí | 2.65 | no | sí | 0.80 |
| **Cono p = 1.5** | negativa | 4 | sí | **0.58** | **sí** | **no** | 0.38 |
| **Cono p = 3** | negativa | 4 | sí | **0.41** | **sí** | **no** | 0.34 |
| Cono p = 6 | informativa | 4 | sí | 1.49 | no | (sin voto) | 0.67 |
| Órdenes producto 3, 4, 5 | negativas | 3, 4, 4 | sí | 5.5, 6.3, 3.4 | no | sí | 3.1, 2.1, 3.3 |
| KR, percolación, capas | negativas | | | ≥ 13 | no | sí | |

Umbrales δ_forma: 0.080 y 0.052 en 4D (N = 1600 y 2400); 0.047 y 0.038 en 3D.

## 3. Lectura

**Lo que sí hizo la distancia de forma (nivel 1).** Con semillas nuevas, rechazó con margen a los conos de caras planas (cúbico 3.05, octaédrico 2.65), que la distancia de F5′ no separaba (aquí, de nuevo, el octaédrico pasa con D: 0.80). Reconoció a ciegas un espacio-tiempo curvo nuevo (FRW de radiación, 0.49).

**Falla 1: resolución (nivel 1).** Los conos p = 1.5 y p = 3 son indistinguibles del redondo con N ≤ 2400; p = 6 sí se separa (1.49). En el perfil medio (N = 2400, sin el desplazamiento parejo), el hueco de links es −0.18 y −0.20 en los conos p = 1 y ∞, −0.09 en p = 6 y solo −0.03 en p = 1.5 y p = 3, igual al ruido. El perfil detecta conos con aristas marcadas, no conos «algo deformados».

**Falla 2: confunde curvatura con forma del cono (nivel 1).** de Sitter fuerte no sube el perfil de forma pareja: lo **inclina**, con menos links (−0.09) y más intervalos grandes (+0.05). Ese hueco en k = 0 es **mayor** que el de los conos p = 1.5 y 3. Así que ningún umbral sobre D_forma puede aceptar a de Sitter fuerte y rechazar a esos conos a la vez: el perfil de intervalos, a este tamaño, no separa «curvatura» de «cono no redondo». Ya en F5′ de Sitter normal tenía D_forma mayor que el plano; con el doble de curvatura el efecto se duplica.

**Interpretación propuesta (nivel 2).** Los intervalos pequeños de un espacio-tiempo curvo son casi planos, pero los que pesan en N_k con N ≤ 2400 abarcan buena parte de la región, donde la curvatura sí se nota. Un criterio que no confunda curvatura con forma tendría que medir en intervalos chicos comparados con el radio de curvatura, que es justamente lo que propone el perfil local de F6. Esto es una hipótesis, no un resultado.

## 4. Lo que deja F7

- **Estado del programa:** la combinación firma espectral + perfil reconoce espacio-tiempo 4D frente a órdenes producto, KR, percolación, capas y conos con aristas marcadas (con D_forma), y acepta curvatura moderada. **No** distingue conos ligeramente deformados ni tolera curvatura fuerte.
- **F6 pospuesta** (R10-3). El Consejo debería reconsiderarla: el perfil local en sub-intervalos podría atacar también la falla 2, porque dentro de un intervalo chico la curvatura pesa menos. A cambio pierde resolución (m = 20 elementos), lo que agravaría la falla 1.
- **Pregunta abierta para el usuario:** si «distinguir un cono p = 3 de uno redondo» es una exigencia razonable para Ω, o si basta con reconocer espacio-tiempo frente a estructuras que no tienen conos suaves.
