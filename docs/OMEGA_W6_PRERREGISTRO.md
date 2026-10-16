# Ω — W6: condición de uso de RC-3 a tres escalas (prerregistro y calibración)

Estado: **PRERREGISTRO**, escrito antes de cualquier código o medida de esta fase.
Rama: `claude/omega-w0`. Se congelará como `claude/omega-w6-congelado`.

Mandato: el Consejo adoptó W6 como **condición de uso** de RC-3, no como teoría. El Auditor exige formalizar «grado saturado» (ventana y tolerancia) antes del próximo candidato. Por la lección 25, el umbral lleva un cálculo de potencia.

RC-3 no se modifica: sus umbrales X1–X4 siguen congelados. W6 decide cuándo un resultado «no excluido» de RC-3 puede leerse como evidencia de dimensión finita.

## 1. Definición

Un candidato G es **W6-válido** si se cumplen las cuatro condiciones:

- **W6.1, misma estructura.** Tres tamaños N₁ < N₂ < N₃ con N_{i+1}/N_i = 8 (10⁴, 8·10⁴, 6.4·10⁵), en una de dos formas:
  - de una **misma estructura creciente**: la de N_i es subestructura de la de N_{i+1};
  - o de una **familia determinista** (sin azar).
  - Las realizaciones independientes no valen.
- **W6.2, δ estable.** |δ₂₃ − δ₁₂| ≤ τ_δ, con δ₁₂ = δ(N₁, N₂) y δ₂₃ = δ(N₂, N₃) calculados como en RC-3.
- **W6.3, grado saturado.** |k₃ − k₂| / k₂ ≤ τ_k, con k_i el grado medio de la componente gigante.
- **W6.4, W5 en ambos pares.** E/N cambia menos del 25 % entre tamaños consecutivos.

Cuando haya más de tres tamaños disponibles, W6.2 y W6.3 se aplican a cada par de pares consecutivos. Además, δ no puede descender en todos los pasos (prohibida la tendencia sistemática δ → 0).

## 2. Calibración de τ_δ y τ_k (regla fijada antes de medir)

**Controles de dimensión finita** (deben ser W6-válidos):

| Control | Forma creciente |
|---|---|
| Cajas J1(2), J1(3), J1(4) (`tools/r3_0b.py`) | Anidadas por longitud de cadena |
| Toros Z², Z³, Z⁴ | Deterministas |
| Triangular abierta | Determinista, caja |
| FCC | Determinista |
| RGG₂ y RGG₃ de densidad fija en caja abierta, grado medio 8 | Anidadas: se generan los puntos de N₃ y se toman los contenidos en la subcaja de volumen N_i/N₃ con la misma esquina |

**Controles de crecimiento intermedio** (deben ser W6-inválidos): J3, semillas 0–2, misma estructura creciente. Ya están medidos en `results/r3_0b_inspect/` (I2) y se reutilizan sin volver a medir.

**Regla de umbrales:**
- τ_δ = max(0.02, 1.5 · max |δ₂₃ − δ₁₂| sobre los controles finitos).
- τ_k = max(0.05, 1.5 · max |k₃ − k₂|/k₂ sobre los controles finitos).

Los controles finitos son W6-válidos por construcción del umbral. La prueba real es la potencia.

**Potencia (condición de validez de W6):**
- **P1.** Las 3 semillas de J3 son W6-inválidas con esos umbrales.
- **P2.** Ningún umbral queda fijado por un único control atípico: si quitar el control que fija τ cambia el veredicto de algún J3, se registra como «potencia frágil».
- **P3.** τ_δ ≤ 0.10. Si fuera mayor, W6.2 no tendría resolución útil frente a δ ≈ 1/d, porque 1/3 − 1/4 = 0.083, y se registra W6.2 como no informativa.

**Veredicto:**
- **W6 VÁLIDA** si se cumple P1.
- **W6 INVÁLIDA** si no: hacen falta más tamaños o más semillas; no se ajusta nada.
- P2 y P3 condicionan la lectura, no el veredicto.

**Predicción del cerebro:**
- W6.3 (grado) aportará casi toda la potencia, y W6.2 (δ) poca.
- En J3 el grado sube ≈ 17 % por paso (I2), y en las cajas abiertas ≈ 3–5 % por efecto de borde.
- Esperado: τ_k ≈ 0.07, τ_δ ≈ 0.03–0.06, y P1 cumplida por W6.3.

## 3. Ejecución

- Agente Sonnet: `tools/w6_calib.py`, tests y `results/w6/`. Medida con `tools/rc3.py::measure`, sin cambios.
- RNG: `rng_from_key((MASTER_W6, id, semilla, i))`, con `MASTER_W6 = 20261022`.
- Revisión del cerebro. Resultados en §4 de este mismo documento.

---

## 4. Resultados (revisados por el cerebro)

- Código: `tools/w6_calib.py`; 7/7 tests pasan.
- Datos: `results/w6/`.
- Tiempo de ejecución: 16.5 min. Las filas de J3 se reutilizan de I2, sin volver a medir.

| Familia | δ₁₂ | δ₂₃ | \|Δδ\| | k₁ → k₂ → k₃ | \|Δk\|/k₂ | W6 |
|---|---|---|---|---|---|---|
| J1(2) | .505 | .048* | .458 | 3.96 → 3.99 → 4.00 | .002 | válido |
| J1(3) | .339 | .336 | .003 | 5.73 → 5.86 → 5.93 | .012 | válido |
| J1(4) | .256 | .256 | .000 | 7.20 → 7.53 → 7.71 | .025 | válido |
| Z² | .503 | .168* | .335 | 4 | 0 | válido |
| Z³ | .341 | .337 | .003 | 6 | 0 | válido |
| Z⁴ | .259 | .257 | .002 | 8 | 0 | válido |
| triangular | .511 | .136* | .375 | 5.92 → 5.99 | .003 | válido |
| FCC | .345 | .339 | .006 | 12 | 0 | válido |
| RGG₂ | .491 | .219* | .273 | 7.80 → 7.99 | .008 | válido |
| RGG₃ | .328 | .318 | .010 | 7.63 → 7.88 | .019 | válido |
| **J3 s0** | .073 | .153 | .079 | 9.41 → 11.36 → 13.26 | **.167** | **inválido** |
| **J3 s1** | .136 | .082 | .055 | 9.90 → 11.71 → 13.86 | **.184** | **inválido** |
| **J3 s2** | .102 | .092 | .010 | 10.08 → 12.44 → 14.15 | **.138** | **inválido** |

\* R = 200 = r_cap en N₃: en las familias 2D, el radio de media masa supera el tope de RC-3 (R_MAX = 200, punto ciego W4 ya declarado). Lo comprobé en `runs.jsonl`: R = 200.0 en las cuatro filas 2D a 6.4·10⁵.

**Umbrales (regla congelada):**
- τ_δ = 1.5 · 0.458 = **0.686**, fijado por J1(2).
- τ_k = **0.05** (el suelo; 1.5 · 0.025 = 0.037).

**Veredicto: W6 VÁLIDA.**
- **P1** se cumple: las 3 semillas de J3 son inválidas, todas por W6.3.
- **P2** se cumple: sin J1(2), τ_δ = 0.56 y ningún veredicto cambia.
- **P3** falla: τ_δ = 0.69 > 0.10, así que **W6.2 queda registrada como no informativa**.

### 4.1 Lectura

1. **Toda la potencia viene del grado** (W6.3), como predije.
   - En J3 el grado sube un 14–18 % por paso.
   - En los controles finitos sube como mucho un 2.5 % (efecto de borde de la caja 4D).
   - El margen es de un factor ≈ 5.5 sobre el control más cercano.
2. **W6.2 no es informativa por una razón conocida, no por ruido.**
   - En d = 2, 6.4·10⁵ nodos superan el tope de radio de RC-3 (W4).
   - En los controles bien resueltos (d = 3, 4), |Δδ| ≤ 0.010.
   - Excluir los controles 2D daría τ_δ ≈ 0.02, pero **es un cambio de la regla después de ver datos y no se aplica**.
   - Propuesta para el futuro, que requiere decisión del Consejo y una fase nueva: W6.2 solo se evalúa en pares con R < r_cap en todos los tamaños.
   - Aun con τ_δ = 0.02, J3 s2 (|Δδ| = 0.010) pasaría W6.2. **La estabilidad de δ no detecta por sí sola el crecimiento intermedio.**
3. **Alcance.** La potencia de W6.3 está demostrada frente a **una** familia de crecimiento intermedio (J3, nivel 1). Que el grado creciente acompañe siempre al crecimiento intermedio es nivel 4.
4. **Contraste de predicciones:**

   | Predicción | Resultado |
   |---|---|
   | W6.3 aporta la potencia | Acertada |
   | τ_k ≈ 0.07 | Fallida: 0.05, el suelo |
   | τ_δ ≈ 0.03–0.06 | **Fallida**: 0.69. No anticipé que el tope de radio en 2D entrara en la calibración |
   | P1 por W6.3 | Acertada |

**Regla de uso vigente para RC-3, desde hoy y no retroactiva.** Una estructura no excluida por RC-3 solo cuenta como evidencia de dimensión finita si es W6-válida:
- tres tamaños de la misma estructura creciente;
- W5 en ambos pares;
- grado saturado con |Δk|/k₂ ≤ 0.05.

W6.2 se calcula y se reporta, pero no decide.
