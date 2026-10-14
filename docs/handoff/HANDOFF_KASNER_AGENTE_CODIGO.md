# Handoff: generador Kasner (R11-4) para el agente de código (2026-10-14)

**Motivo:** toda la curvatura probada hasta hoy es conformemente plana: solo cambia la densidad, la relación causal es la de Minkowski. R11-4 pide una positiva **no** conformemente plana: Kasner de vacío.

## Física (verificar, no suponer)

- Métrica: ds² = −dt² + Σᵢ t^{2pᵢ} dxᵢ², con Σpᵢ = 1 y Σpᵢ² = 1 (vacío). El elemento de volumen es √−g = t.
- **Relación causal:** el espacio-tiempo es homogéneo en x, así que el conjunto alcanzable desde (t₁, x₁) hasta el tiempo t₂ es la suma de Minkowski (integral) de los elipsoides de velocidades {v : Σ aᵢ(t)² vᵢ² ≤ 1}, con aᵢ = t^{pᵢ}. Es convexo y su función soporte es h(u; t₁, t₂) = ∫_{t₁}^{t₂} √(Σ uᵢ²/aᵢ(t)²) dt.
- Entonces y está en el futuro causal de x si y solo si t₂ > t₁ y u·Δx ≤ h(u; t₁, t₂) para toda dirección unitaria u.
- **Implementación propuesta:** un conjunto fijo y simétrico de K direcciones (por ejemplo K ≥ 400, casi uniformes en la esfera, simétricas por ±u) y tablas acumuladas H_u(t) en una grilla fina de t, de modo que h = H_u(t₂) − H_u(t₁).
- Como la función soporte es aditiva en el tiempo, la relación aproximada «u·Δx ≤ h para las K direcciones» **es transitiva** con cualquier K. Debe verificarse con un test.

## Entregar

1. **`omega/causality/kasner.py`** (tipado estricto):
   - `kasner_causal_matrix(coords, p, K, ...)`.
   - `sprinkle_kasner_interval(n, p, t_lo, t_hi, rng, ...)`: n puntos uniformes en volumen propio (densidad ∝ t en coordenadas) dentro del intervalo causal J⁺(a) ∩ J⁻(b), con a = (t_lo, 0) y b = (t_hi, 0). Por rechazo: caja de coordenadas que contenga el intervalo, aceptación ∝ t, y pertenencia con la misma relación.
   - `kasner_exponents(...)` o validación de (p₁, p₂, p₃) que cumpla las dos sumas con tolerancia.
   - No elijas los exponentes ni los tiempos de la prueba: los decide otra persona a ciegas. Solo exige la condición de vacío.
2. **Tests en `tests/test_kasner.py`** (rápidos):
   - irreflexividad, antisimetría y transitividad sobre unos cientos de puntos (`is_partial_order`);
   - **límite plano:** p = (1, 0, 0) cumple las dos condiciones de vacío y es el Kasner plano trivial (región de Milne). La relación debe coincidir con la de Minkowski tras T = t·cosh x₁, X = t·sinh x₁ (y, z iguales), salvo una fracción ínfima de pares en el borde por la discretización en K direcciones; acotar esa fracción;
   - los puntos muestreados están dentro del intervalo y su densidad en t es ∝ t (prueba estadística simple);
   - la convergencia en K: con K = 200 y K = 800 las relaciones difieren en menos de una fracción pequeña de los pares.
3. **Rendimiento:** medir el tiempo con N = 2400 y N = 4000. Debe ser viable por bloques.
4. **Verificación:** `python -m pytest -q tests/test_kasner.py` y `python -m mypy omega/causality/kasner.py tests/test_kasner.py` limpios.

## No hacer

No tocar otros archivos (otro agente trabaja en `tools/f6_exploracion.py` y `tests/test_f6_exploracion.py`), no correr ninguna prueba oficial, no hacer commit.
