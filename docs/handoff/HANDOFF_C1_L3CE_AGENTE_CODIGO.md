# Handoff C1 y L-3c-E: orquestador → agente de código

- **Especificación que manda:** `docs/OMEGA_2_PROPUESTA_CAUSAL.md` §11.

## Qué construir

1. `omega/causality/estimators.py`: `bd_action_2d_smeared(C, eps)` según §11.2.
2. `omega/causality/exact.py`: `expected_smeared_action_2order(n, eps)` según §11.2.
3. `omega/causality/mcmc.py`: cadena Metropolis sobre 2-órdenes de §11.6. Estado = permutación v (u = identidad). Funciones puras, RNG inyectado. Devuelve muestras (matrices o métricas) y tasa de aceptación.
4. `tools/c1_l3ce.py`: modos smoke y full; C1-2 y L-3c-E con los veredictos de §11; summary preregistrado → final; árbol limpio en full.
5. Tests: C1-1, detalle del movimiento (v sigue siendo permutación), β = 0 acepta todo.

## Reglas

- Las de los handoffs anteriores. Entropías exactamente las de §11.
