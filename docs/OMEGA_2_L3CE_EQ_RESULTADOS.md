# Ω-C: resultados de L-3c-EQ (equilibrio de la rama euclídea por Wang-Landau)

- **Fecha:** 2026-10-12. **Especificación:** `OMEGA_2_PROPUESTA_CAUSAL.md` §16, enmiendas EQ-A1 y EQ-A2.
- **Commits de ejecución:** validación fallida en `0ef0f1b`; validación y WL de N = 32 en `0282d15`. Registro de tramos: `results/l3ce_eq/chunks.jsonl`.

## 1. Resultado oficial: **INDETERMINADO (cobertura no cumplida en N = 32)**

| Paso | Resultado |
|---|---|
| V, primer intento (5·10⁵ pasos) | FALLA (0.19 > 0.15), por ruido estadístico (EQ-A2) |
| V, segundo intento (2·10⁶ pasos, entropía nueva) | PASA (0.039) |
| WL N = 32, 6·10⁶ pasos | Mínimo visitado S = −504; sándwich conocido S = −538. **La cobertura falla.** ln f nunca bajó de 1: el histograma no llegó a ser plano. |
| WL N = 32, segunda mitad (3·10⁶ pasos) | El mínimo visitado **no se movió** de −504 |
| N = 64 | **No ejecutado (declarado):** con la cobertura de N = 32 fallida, el veredicto preregistrado es INDETERMINADO cualquiera sea el resultado de N = 64, y su coste habría sido ≈ 30 minutos. |

## 2. Lo que se aprende (nivel 1, descriptivo)

1. **La fase de baja energía no es la bicapa sino un sándwich** (EQ-A1): dos capas grandes separadas por una capa central de 2 elementos, con S = −N²/2 − N + 6. El término −2N₂ de la acción premia a los pares con exactamente 2 elementos entre ellos.
2. **Esa fase es cinéticamente inaccesible para movimientos locales.** Ni Metropolis (L-3c-E″) ni Wang-Landau con transposiciones logran explorarla: el paseo queda trabado unas 34 unidades de acción por encima del sándwich en N = 32. Es la firma de un paisaje vidrioso: llegar al orden de baja energía requiere reorganizaciones colectivas que un intercambio de dos elementos no produce.
3. Con los métodos disponibles en este entorno, **la pregunta de equilibrio de la rama E no se puede responder.** El PASA de L-3c-E′ sigue valiendo solo como estabilidad local de la geometría.

## 3. Nuevo hueco detectado en la rama L (consecuencia lateral)

La derivación A-1 del Bloque A supone que en las familias por capas los pares no consecutivos tienen muchos elementos entre ellos, de modo que f(n, ε) es despreciable. **En el sándwich la capa central tiene 2 elementos**, así que los pares extremo-extremo tienen n ∈ {0, 1, 2} y f(n, ε) no es despreciable. La ventana universal de A-1 **no cubre a la familia sándwich**. Hay que extender A-1 y A-2 a esa familia.

## 4. Qué haría falta para la rama E

- Movimientos no locales: por ejemplo, mover un bloque entero de elementos de una capa a otra, o un algoritmo de clústeres adaptado a 2-órdenes.
- O métodos de población (poblaciones con reponderación, «population annealing») que toleran mejor los paisajes vidriosos.
Ambos son trabajo de desarrollo, no una corrida.
