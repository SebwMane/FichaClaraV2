# Handoff C0: orquestador → agente de código

- **Emisor:** orquestador (sesión 2026-10-11).
- **Receptor:** agente de código (rol Sonnet; en esta sesión lo ejecuta la misma instancia con separación de roles, ver nota final).
- **Especificación que manda:** `docs/OMEGA_2_PROPUESTA_CAUSAL.md` §4. Si este handoff y el documento difieren, manda el documento.

## Qué construir

1. `omega/causality/sprinkling.py`
   - `sprinkle_diamond(n, d, rng) -> NDArray[float64]` de forma (n, d), columna 0 = tiempo. Diamante {|x| ≤ T − |t|} con T = 1/2 (altura propia τ = 1). Muestreo por rechazo en la caja.
   - `interval_volume(tau, d) -> float` con c_d = π^{(d−1)/2} / (2^{d−1} d Γ((d+1)/2)).
   - `causal_matrix(coords) -> NDArray[bool_]`, C[i, j] = True si i ≺ j estricto (Δt > 0 y Δt² ≥ |Δx|²).
2. `omega/causality/orders.py`
   - `random_2order(n, rng)` (intersección de dos órdenes totales; equivale al sprinkling 2D).
   - `kr_order(n, rng, p=0.5)` y `transitive_percolation(n, p, rng)`.
   - `transitive_closure(adj)` y `is_partial_order(C)` (irreflexivo, antisimétrico, transitivo).
3. `omega/causality/estimators.py`
   - `ordering_fraction`, `mm_dimension` (inversión monótona de r(d) en d ∈ [1, 12]; fuera de rango devuelve nan).
   - `between_counts(C)`: matriz de elementos estrictamente entre i y j (C·C restringido a pares relacionados).
   - `bd_action_2d(C)` con S = N − 2N₀ + 4N₁ − 2N₂.
   - `midpoint_dimension(C)` según §4 C0-4.
4. `tools/c0_validation.py`: modos `smoke`, `pilot`, `full`; escribe `results/c0_validation/` (o `results/c0_pilot/`) con `summary.json` preregistrado → final y `runs.jsonl`. Modo full exige árbol limpio.
5. `tests/test_causality.py`.

## Reglas no negociables (del repo)

- RNG solo vía `np.random.Generator(PCG64(...))`; nada de `np.random.rand` ni `random` (lo verifica `tests/test_architecture.py`).
- `mypy --strict` limpio en los archivos nuevos.
- Ningún estimador recibe coordenadas, solo la matriz de orden.
- No tocar nada fuera de `omega/causality/`, `tools/c0_validation.py`, `tests/test_causality.py` y `results/c0_*`.
- No ajustar tolerancias. La única tolerancia abierta es τ de C0-4 y la fija el orquestador tras el piloto.

## Entregables al orquestador

- Código y tests pasando.
- Una nota breve de dudas o desviaciones en el mensaje del commit.

## Nota de transparencia

En esta sesión no había forma de lanzar un modelo Sonnet independiente. El rol de código lo ejecuta la misma instancia, ciñéndose solo a este handoff, y el orquestador audita después como si fuera código ajeno (`docs/handoff/AUDITORIA_C0.md`).
