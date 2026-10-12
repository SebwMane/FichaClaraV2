# Ω-C: resultados de L-2c (compuerta de paisaje, ramas E y L con criterios separados)

- **Fecha:** 2026-10-12.
- **Commit de ejecución:** `bce9d95` (árbol limpio, 54 s).
- **Especificación congelada:** `OMEGA_2_PROPUESTA_CAUSAL.md` §10 (preregistro en `bf39667`, enmiendas L2c-A1 a A3 en `bce9d95`, todas antes de la corrida).
- **Datos:** `results/l2c/summary.json` y `runs.jsonl`.

## 1. Rama euclídea (e^{−βS}, restringida a 2-órdenes)

| Orden (N = 1024) | S_BD | Marcado por C0-4b |
|---|---|---|
| 2-órdenes aleatorios (20 semillas) | media ≈ 1 (es geometría) | 100% |
| bicapa (k = 2) | −523 264 | no |
| k = 4 | −392 192 | no |
| k = 16 | −121 856 | no |
| k = 128 | −15 232 | no |

**Decisión: «geometría NO preferida a β → ∞»** (ΔS_E = −523 277). Era la expectativa preregistrada.

**Lectura preregistrada:** esto **no cierra** la rama E. A β = 0 el 2-orden típico es geométrico; a β → ∞ gana un orden en capas. La pregunta que queda, y que decide la rama, es si existe una ventana 0 < β < β_c con fase geométrica estable. Eso lo responde **L-3c-E** (Metropolis sobre 2-órdenes), ya autorizada por la ratificación de R3-2.

Escala de afirmación: nivel 1 (los valores de S de las capas son exactos: S = N − 2(k−1)(N/k)²).

## 2. Rama lorentziana (e^{iλS}, todos los órdenes)

| Prueba | Resultado |
|---|---|
| **L-i** instrumento (bipartito, fórmula exacta |cos λ|^{n²/4}) | **PASA**: 8 puntos, máxima desviación 1.37·SE |
| **L-ii** ventana de supresión de KR | **λ mod π ∈ (π/3, 2π/3)** (analítica, nivel 1 para el bipartito, nivel 2 para KR) |
| **L-iii** coherencia de la clase geométrica | **INDETERMINADO**: ningún punto utilizable |

**Qué pasó en L-iii.** El promedio de fase de los 2-órdenes geométricos quedó bajo el umbral de resolución (3/√M = 0.021) **ya en N = 16**, en los tres λ. Con 20 000 muestras no hay nada que ajustar. La dispersión de la acción explica por qué: sd(S) = 16.8, 30.8, 53.2, 88.4 y 142.0 para N = 16 a 256 (crece como ≈ N^0.77). Con una dispersión así, e^{iλS} gira tantas vueltas entre muestras que se cancela casi por completo.

**Lo que esto significa (nivel 2, explicación propuesta).** La acción BD **sin suavizar** es demasiado ruidosa para la rama lorentziana: hace cancelarse no solo a KR sino también a la propia geometría. No es un hallazgo nuevo: la literatura introduce una escala de no localidad (acción «suavizada») justamente para controlar esas fluctuaciones. Ahora lo verificamos con datos propios.

**Extrapolación NO demostrada (nivel 4), solo como orientación.** Si la cancelación de la clase geométrica creciera como N^{1.5} (lo que sugiere una aproximación gaussiana con sd ∝ N^0.77), seguiría siendo menor que la de KR (∝ N²) dentro de la ventana. Esto no es un resultado: L-iii no pudo medirlo.

## 3. Estado de la compuerta

| Rama | Estado | Siguiente paso |
|---|---|---|
| E | Viva, no decidida | L-3c-E: Metropolis sobre 2-órdenes con e^{−βS}, barrido en β, bandera C0-4b |
| L | Bloqueada por el instrumento | Requiere decidir si Ω-C usa la acción suavizada (nueva resolución R3-4) |

R3-1 (adoptar Ω-C) sigue pendiente: ninguna rama cerró todavía.
