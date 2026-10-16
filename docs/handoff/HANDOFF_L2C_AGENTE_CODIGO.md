# Handoff L-2c: orquestador → agente de código

- **Especificación que manda:** `docs/OMEGA_2_PROPUESTA_CAUSAL.md` §10.

## Qué construir

1. `omega/causality/orders.py`:
   - `layered_2order(n, k, rng)`: k antichains apiladas como 2-orden (bloques de rangos consecutivos en u y en v; dentro de cada bloque u creciente y v decreciente); etiquetas permutadas.
   - `random_bipartite(n, p, rng)`: n//2 abajo, el resto arriba, cada par entre capas relacionado con prob. p (sin cierre: ya es transitivo).
2. `omega/causality/phase.py`: `phase_average(actions, lam) -> (abs, se)` con se por bootstrap analítico: SE de |media| aproximado por sqrt(var(cos λS)+var(sin λS))/sqrt(M); `bipartite_phase_exact(n, lam) = |cos lam|^{(n//2)*(n - n//2)}`.
3. `tools/l2c_landscape.py`: calcula L-2c-E, L-i y L-iii con las entropías de §10; escribe `results/l2c/` con summary preregistrado → final y veredictos; full exige árbol limpio.
4. Tests.

## Reglas

- Las del handoff C0. No ajustar constantes; si algo de §10 es ambiguo, se elige la lectura literal y se anota en el commit.
- Para L-iii, el filtro de puntos utilizables es exactamente |φ̂| > 3/√M.
