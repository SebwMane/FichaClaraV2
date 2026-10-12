# Handoff C0-4b: orquestador → agente de código

- **Especificación que manda:** `docs/OMEGA_2_PROPUESTA_CAUSAL.md` §9. Si este handoff difiere, manda el documento.
- **Base:** commit con §9 ya congelado.

## Qué construir

1. En `omega/causality/orders.py`: `thick_layer_order(n, layers, p, rng)` (capas de tamaño n // layers, el resto a la última; relaciones entre capas consecutivas con prob. p; cierre transitivo; índices permutados).
2. `omega/causality/profile.py`:
   - `K_MAX = 8`.
   - `chain_abundances(n, k_max)` exacta.
   - `reference_abundances(d, n, k_max, rngs)`: d = 1 cadena, d = 2 exacta, d = 3, 4 media sobre los generadores dados.
   - `profile_distance(nk, ref)` = max |ln((nk+1)/(ref+1))|.
   - `d_ref(d_hat)`: redondeo recortado a {1..4}; None si nan.
3. `tools/c0_4b_validation.py` con modos `smoke` y `full`; summary preregistrado → final con veredicto calculado; full exige árbol limpio. Guardar δ_d calibrados y todas las D por semilla.
4. Tests en `tests/test_causality.py`.

## Reglas

- Mismas del handoff C0: RNG solo PCG64, mypy --strict, sin coordenadas en estimadores, no tocar Ω-1.1.
- Las entropías son las de §9; no inventar otras.
- No ajustar ninguna constante. Dudas, en el mensaje del commit.
