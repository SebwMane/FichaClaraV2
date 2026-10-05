# Ω — Resultados del bloque L (Consejo, Revisión 2)

- **Fecha:** 2026-10-05.
- **Rama:** `claude/hopeful-galileo-88u1l1`.
- **Especificación congelada:** `docs/OMEGA_CONSEJO_FASE_SIGUIENTE.md`, «REVISIÓN 2» (R2.4, R2.5 y enmiendas R2.8 L-A1…L-A5).
- **Esta fase no modificó nada de lo siguiente:**
  - umbrales del certificado;
  - S0;
  - D_eff;
  - criterios de O-05.
- **Versión:** sigue siendo Ω-1.1.

## 1. Resultados oficiales

| Prueba | Commit de ejecución | Decisión congelada | Datos |
|---|---|---|---|
| **L-1** identidades y cotas | suite (`tests/analytical/test_l1_identities.py`) | **100%** | Identidades de S_smooth (binaria y ponderada), ΔS = −(α+12η)ΔT bajo swaps que preservan grados, cota de vecindario k-regular con igualdad ⇔ ⊔K_{k+1} (atlas ≤ 7 nodos), Lovász, tr W³ ≤ (ΣW²)^{3/2} ≤ (ΣW)^{3/2} |
| **L-3a** lemas KKT | suite (`tests/analytical/test_l3a_kkt.py`) | **100%** | P-KKT, P-KKT-reg, lema de escala (con test no vacuo: 3×K10 es KKT ⇔ α̂ ≥ 3.5), T³ y RGG3 nunca son KKT en S0, T³ no es KKT en Ω-B, K22 + halo. Matiz L-A3 |
| **L-2** paisaje, malla absoluta α̂ ≤ 3 | `42a2c43` | **ÉXITO** | ΔS > 0 en 180/180 puntos × 5/5 semillas. **Gana el uniforme en el 100%**; la malla no alcanza el régimen de cliques |
| **L-2b** paisaje, escala Ω-B α̂ = f·α̂_c (declarado tras L-2, L-A4) | `3cba45b` | **ÉXITO** | ΔS > 0 en 180/180 × 5/5. Mejor no geométrica: clique colex 810/900, uniforme 60, unión de K13 30. Diagnóstico local (r ≤ 1.5·r12): ΔS_local ≥ 3.8e-3 en el 100% |
| **O4B-1** estados finales de O-04 | `a963b77` | **ÉXITO** | 1110/1110 sha256 verificados. Ω-B (270): uniforme 70, clique única 198, multi-clique 2, «otro» 0. Todos los convergidos son KKT. S0 (informe): vacío 420, uniforme 420 |
| **L-3b** estabilidad de arranques geométricos (N = 216, 1680 corridas) | `008cfc6` | **SÍ (incompatible)** | Persistentes con ε = 1e-2: 0. Metaestables con ε = 1e-3: 0. Deriva en MAX_STEPS: 0. Convergidos no KKT: 0/1659. Alarma de controles: no. Validación L-A5 del clasificador: OK |

### Detalle de L-3b (inputs geométricos: T³ 6³, RGG3 k12 y T³ decorada 3³⊗K8)

- **S0 (588 corridas):** todas terminan en `vacío`.
  - Los pesos residuales (~1e-7, parada por `tol_step`) conservan el patrón del input (R = J = H = 1).
  - Es un **apagado en amplitud** del estado geométrico, no una reorganización. Lo predice el lema de escala: c_ij < (N−2)/α̂ en toda la malla.
  - La clase `vacío` lo excluye por diseño.
- **Ω-B (420 corridas):**

  | Clase | Corridas |
  |---|---|
  | cliques solapadas | 158 |
  | multi-clique | 115 |
  | uniforme | 85 |
  | clique única | 45 |
  | otro | 17 |

  - En todos los no vacíos: J ≤ 0.43, H ≤ 0.59 y R_orig ≤ 0.60. **Ningún final no vacío tiene J ≥ 0.5.**
  - Los «otro» tienen J ≤ 0.08 y H ≈ 0.22–0.35: estados colapsados.
- **Controles:**
  - La unión de K8 nunca es persistente.
  - El ER tampoco: la dinámica no conserva estructura dispersa arbitraria.

## 2. Escala de afirmación tras el bloque L

**Demostrado o verificado (nivel 1)**
- Identidades y cotas de L-1.
- Lemas KKT de L-3a: el soporte de un punto KKT de S0 con γ = 0 es una unión disjunta de cliques; ningún γ estabiliza un estado regular como T³; las aristas de peso 1 exigen codegrado ≥ (N−2)/α̂.
- En el dominio probado:
  - ninguna referencia geométrica evaluada tiene menor acción que la mejor no geométrica, ni en la escala absoluta ni en la de Ω-B (L-2, L-2b);
  - ningún arranque geométrico persiste bajo S0 ni Ω-B (L-3b);
  - todos los finales de Ω-B de O-04 son uniformes o uniones de cliques KKT (O4B-1).

**Explicación propuesta (nivel 2), ahora con soporte computacional directo**
- La recompensa por triángulo exige codegrado extensivo. Por eso S0 apaga cualquier geometría de grado acotado (F0) y Ω-B la reorganiza en cliques (disjuntas, solapadas o con halo).

**Teorema externo (nivel 3)**
- Kruskal–Katona/Lovász, para el mínimo global binario de Ω-B.

**Extrapolación NO demostrada (nivel 4)**
- Pesos óptimos arbitrarios, N ≫ 216, otras restricciones y Θ > 0.
- La analogía de Chatterjee–Diaconis para Θ > 0 sigue siendo solo una analogía.

## 3. Decisión según el árbol preregistrado (R2.3)

L-1, L-2/L-2b y L-3 coinciden: **rama «SÍ»**.

> «Bajo las restricciones probadas (S0 y Ω-B, N = 216, α̂ ∈ [0.5, 3] en S0 y α̂ = f·α̂_c con f ∈ [0.25, 3] en Ω-B, γ̂ ∈ {0, 1, 10, 100}, η = μ = 0 en la dinámica), los estados que la dinámica selecciona son sistemáticamente incompatibles con la geometría 3D buscada. Partiendo de una geometría, S0 la apaga y Ω-B la convierte en cliques.»

**Consecuencia preregistrada:** limitar la familia y buscar una nueva hipótesis dinámica. O-05 pasa a ser una prueba de la analogía de campo medio (MF-1), ya no la vía principal hacia P1.

**Sigue abierto (no cerrado por el bloque L):**
1. Θ > 0. L-3 es determinista; la cuestión entrópica solo la prueba O-05.
2. η > 0 en la dinámica. L-1 y L-2/L-2b lo cubren en el paisaje, no en la dinámica.
3. N grande.
4. Funcionales fuera de M§13–15.

## 4. Lo que queda aprobado y no ejecutado (requiere la decisión del usuario sobre prioridad)

| Ítem | Coste | Valor tras el bloque L |
|---|---|---|
| RP-1 + s06 + O-06 (cadena de la compuerta) | ~4–5 h | Reproducibilidad y cierre del orden P§22. Bajo valor científico nuevo |
| O-05 + MF-1 + PG-1 (piloto previo) | 5–20 h + código | Prueba de la analogía Θ > 0. Prior nulo reforzado por el bloque L, pero es la única prueba empírica de la entropía |
| N-1 (nulo con clustering) | 1–2 h | Solo importa ante un positivo futuro |
| Nueva hipótesis dinámica | — | **Requiere un documento del usuario.** Puerta de entrada: debe superar L-2 y L-3b antes de simularse |
