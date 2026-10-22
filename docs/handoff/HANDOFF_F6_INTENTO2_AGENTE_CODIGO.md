# Handoff F6 intento 2 para el agente de código (2026-10-18)

**Especificación:** `docs/OMEGA_F6_INTENTO2_PREREGISTRO.md` (§2 a §8) y `docs/SELLADO_F6_INTENTO2_CIEGA.md` (D0 a D8).

## Entregar

1. **`tools/f6_intento2.py`.**
   - Importa el instrumento y los generadores de `tools/f6_confirmacion.py` (incluidos C0 a C7, que ahora son «vistas») y de `tools/f6_exploracion.py`. No los copies.
   - **Cambios respecto del intento 1:**
     - Δr = r − 1/10 (referencia exacta; sin parte `ref` para r);
     - 30 realizaciones por N;
     - la condición de colapso K;
     - EXPERIMENT_ID 322 y las entropías y familias de §6;
     - una parte `refperfil` (entropía 20261086, familia 5) solo para el perfil de la rama 2.
   - **Partes:** `refperfil`, `nula_<N>`, `cal_<g>`, `vista_<nombre>` (las 21 de §5: 13 geometrías y conos ya vistos más C0 a C7), `sellada_<k>` (D0 a D8, spawn_key (322, k, N, s)) y `final`. Todas admiten `_<N>` y `_<N>_r<a>-<b>` (rango de realizaciones) para partir el trabajo. `--listar` debe dar una partición en la que **cada parte tarde menos de 8 minutos**; mide antes.
   - **Compuerta:**
     - el preregistro debe contener `R13-1 | **Firmada`;
     - el árbol git debe estar limpio fuera de `results/`;
     - en `sellada_*` y en `final`, el SHA-256 del sello debe ser `b64503ffb0d5f5772a0ab8c8dcf80e689371d75ec6f93a6d93f9912c930866cf` (calculado sobre `docs/SELLADO_F6_INTENTO2_CIEGA.md`).
   - **`final`:**
     - b y A con bootstrap (familia 4) y K con sus 4 dⱼ/seⱼ;
     - b_plano, τ_b y τ_K según §4;
     - la clasificación con las dos condiciones y el veredicto de §5;
     - lo informativo y las verificaciones de D4 y D5;
     - el commit de cada parte usada (todas deben tener el mismo).
     - Salida en `results/omega_f6_intento2/`, con `summary.json` y tablas CSV.
2. **Generadores D0 a D8,** al pie de la letra del sello. Anota cualquier ambigüedad.
   - D0: Kasner con k = 800.
   - D1: cizalla x' = x − b·t.
   - D2 y D3: rechazo con su densidad.
   - D4: candidatos, aux.random, q = 0.10 y cierre transitivo iterado hasta que no cambie.
   - D5: gemelos con φ y la permutación.
   - D6: gauge no par; la región usa F(x) y F(−x).
   - D7: gauge cerrado del huso y su caja de rechazo.
   - D8: 0.68‖x‖₂ + 0.32‖x‖₁.
   - Si una verificación de D4 o D5 falla, la parte debe detenerse con un error (B2).
3. **Tests rápidos** en `tests/test_f6_intento2.py`:
   - D4 y D5 son órdenes parciales y sus conteos son coherentes (D5: exactamente n_dup pares de gemelos);
   - D6, D7 y D8 dan órdenes parciales con los puntos en la región, y el gauge de D7 coincide con la definición por función soporte en unas direcciones;
   - D1 coincide con Minkowski tras la cizalla;
   - K de un caso fabricado.
4. **Verificación:** pytest de los archivos tocados; mypy sin errores nuevos.

## No hacer

No correr partes oficiales ni escribir en `results/omega_f6_intento2/`, no usar entropías oficiales fuera de las constantes (para medir usa `SeedSequence(840001)` en el scratchpad), no tocar otros archivos, no hacer commit.
