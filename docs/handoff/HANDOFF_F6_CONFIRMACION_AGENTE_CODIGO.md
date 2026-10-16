# Handoff F6 confirmatoria para el agente de código (2026-10-14)

**Especificación única:** `docs/OMEGA_F6_CONFIRMACION_PREREGISTRO.md` (§1 a §7) y `docs/SELLADO_F6_CONFIRMACION_CIEGA.md` (construcciones C0 a C7). Implementa solo eso.

## Entregar

1. **`tools/f6_confirmacion.py`.** Reutiliza el instrumento de `tools/f6_exploracion.py`: `causal_blocks`, `between_blocks`, `measure_order`, `interval_obs`, `local_ms`, `fit_poly`, la lógica de `_cells` y `fit_structure`, y los generadores de las estructuras ya vistas. Importa lo que puedas en lugar de copiar.
   - Partes:
     - `ref` (4000 diamantes por m, entropía 20261078);
     - `nula_<N>` (30 diamantes planos por N, entropía 20261080, familia 1);
     - `cal_<grupo>` (15 grupos, cada uno con 10 realizaciones en N = 4000 y 10 en N = 8000, entropía 20261080, familia 2; los grupos 0 a 9 son diamantes planos, 10 a 12 de Sitter y 13 y 14 FRW polvo);
     - `vista_<nombre>` (entropía 20261081, familia 3);
     - `sellada_<k>` (entropías 20261076 y 20261077 con el spawn_key del archivo sellado: (321, k, N, s));
     - `final`.
   - Si una parte no cabe en unos 9 minutos, pártela por N o por realización y documéntalo.
   - Compuerta: el preregistro debe contener `R12-1 | **Firmada` y el árbol git debe estar limpio fuera de `results/`.
   - Constantes: EXPERIMENT_ID 321; m ∈ {10, 20, 40, 80, 160}; N ∈ {4000, 8000}; 10 realizaciones por N; sin dilución.
   - Salida en `results/omega_f6_confirmacion/`.
2. **`final`:**
   - b y A por estructura y grupo, con el IC por bootstrap (1000 réplicas de realizaciones completas; entropía 20261080, familia 4);
   - b_plano = mediana de los 10 grupos planos;
   - τ_b = 1.5 × max sobre los 15 grupos de |b − b_plano|;
   - clasificación |b − b_plano| ≤ τ_b; «no medible» si hay menos de 2 valores distintos de x;
   - veredicto PASA/FALLA según §3, con la zona gris sin voto;
   - lo informativo de §5: signo de A, rama 2 con max-T frente a la nula, y el margen |b − b_plano|/τ_b;
   - para C5, sus verificaciones obligatorias.
   - Escribe `summary.json` y tablas CSV.
3. **Generadores sellados** (sigue el archivo al pie de la letra; si algo es ambiguo, anótalo):
   - C0: `omega.causality.kasner` con k = 800.
   - C1: relación Δt > 0 y Δt ≥ ‖AΔx‖₂ en el diamante estándar.
   - C2: densidad exp(φ) con los modos sacados del generador 20261077 en el orden indicado; aceptación exp(φ − 0.35); relación de Minkowski.
   - C3, C4, C6 y C7: conos con su norma, región {|t| ≤ T − ‖x‖} por rechazo en [−T, T]⁴. **Para C6 no hace falta bisección:** Δt ≥ ‖Δx‖ equivale a Σᵢ |Δxᵢ/Δt|^{pᵢ} ≤ 1 (y lo mismo para la región con T − |t|). Úsalo por bloques; cuidado con 0/0 y con potencias altas.
   - C5: links por B == 0 y C; borrar exactamente round(0.13·|L|) con el generador 20261077, eligiendo con `rng.choice` sobre L en orden lexicográfico; verificar el orden parcial y que el cierre de C' sea C'. Informa |L|, los borrados y la fracción de links nueva.
4. **Tests rápidos** en `tests/test_f6_confirmacion.py`:
   - cada cono sellado da un orden parcial y sus puntos están en la región;
   - la relación de C6 por suma de potencias coincide con la de la norma por bisección en unos pocos pares;
   - C1 coincide con Minkowski tras x' = Ax;
   - C5 sigue siendo orden parcial y su cierre no cambia;
   - la densidad de C2 está en [e^{−0.35}, e^{0.35}] relativa.
5. **Verificación:** pytest de los archivos tocados; mypy sin errores nuevos. Mide los tiempos (Kasner con N = 8000 y k = 800 es el más caro).

## No hacer

No correr partes oficiales ni escribir en `results/omega_f6_confirmacion/`, no usar las entropías oficiales fuera de las constantes, no cambiar criterios, no tocar otros archivos, no hacer commit.
