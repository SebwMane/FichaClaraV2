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
