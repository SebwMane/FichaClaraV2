# Ω — RC-3 (juez de exclusión A−): resultados

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-rc3`.
- **Prerregistro:** `docs/OMEGA_RC3_PRERREGISTRO.md`.
- **Datos:** `results/rc3/{pairs.jsonl, summary.json}` (57 pares N / 8N; 559 s; 9/9 pruebas).
- **Código:** `tools/rc3.py`, escrito por un agente Sonnet y revisado por el cerebro.
  - La interpolación de R, el cálculo de δ y X1–X4 se ajustan a §1.
  - Las subclaves de RNG son las prerregistradas.
  - No se activó ninguna guarda.

## 1. Cifras

| Conjunto | FE (exclusión falsa) | EP (poder de exclusión) | Por d (excluidas / n) |
|---|---|---|---|
| **Familias nuevas** | **0/10 = 0.00** | **17/18 = 0.944** | d2 0/1 · d3 0/2 · d4 0/2 · d5 0/3 · d6 0/1 |
| Familias repetidas | 0/16 | 10/13 = 0.769 | todas 0 |
| Todas | 0/26 | 27/31 = 0.871 | todas 0 |

**No excluidas:**
- nueva: árbol de cubos 5³ (instancia única);
- repetidas: apoloniana s0 (δ = 0.130, justo por encima del corte 0.125) y retazos b3 s0–s1 (SIN_VENTANA en N, así que X4 no puede activarse).

## 2. Veredicto oficial: **RC3-PARCIAL** (con lectura literal = VÁLIDA, declarada)

- **La ambigüedad.** La regla por familia («todas sus instancias excluidas salvo como máximo una») se escribió pensando en las familias con semillas.
  - En lectura literal, una familia determinista de **una sola instancia** que no se excluye la cumple.
  - Exactamente eso ocurrió con el árbol de cubos, la única familia nueva no excluida.
- **Por qué la lectura conservadora.** Elegir la lectura literal **después de ver** que salva el veredicto sería una lectura favorable a posteriori. La regla del cerebro es resolver toda ambigüedad descubierta tras los datos **en contra de la hipótesis**.
- **Con la lectura intencional,** la familia árbol de cubos incumple la regla por familia. Así, el veredicto es **PARCIAL**: FE ≤ 0.10 y EP ≥ 0.75.
- Se registran ambas lecturas.

## 3. Lo que sí queda demostrado: **cero exclusiones falsas en d = 2…6**

- Es el resultado central, y supera lo que pedía el Consejo:
  - ninguna geometría auténtica es excluida, en 26 pares con d = 2, 3, 4, 5 y 6;
  - incluye discretizaciones distintas (cuadrada, triangular, panal, FCC, BCC, diagonales y RGG), bordes abiertos y Heisenberg;
  - el juez de RC-2 rechazaba ya d = 4.
- δ reproduce 1/d en todas las geometrías:

| Dimensión | 1/d | δ medido |
|---|---|---|
| d = 2 | 0.50 | 0.49–0.52 |
| d = 3 | 0.333 | 0.31–0.34 |
| d = 4 | 0.25 | 0.24–0.26 |
| d = 5 | 0.20 | 0.20–0.21 |
| d = 6 | 0.167 | 0.17 |

- Las degeneradas de diámetro logarítmico dan δ entre 0.00 y 0.11.
- Esta precisión es un subproducto, no un criterio: **RC-3 no usa d**, aunque δ la estime.

## 4. Advertencias halladas en la revisión del cerebro (no las informó el agente)

**W1. X4 (coherencia) está sesgada en d alta a N = 10⁴, y el veredicto depende del margen.**
- **Todas** las geometrías de d ≥ 4 son **INCOHERENTES en N**: RGG4 k8, RGG4 k16 s0, T⁵, RGG5, T⁶ y Heisenberg, con γ ≈ 0.47–0.51 y f_low 0.45–1.0.
- Se salvan solo porque en 8N pasan a INTERMEDIO (γ 0.56–0.60, f_low ≤ 0.10), y X4 exige INCOHERENTE en **los dos** tamaños.
- Las familias de atajos y WS se excluyen **solo** por X4, con γ(8N) = 0.43 y 0.30.
- La separación en 8N existe, pero es estrecha: 0.56 frente a 0.43.
- Es la lección 17 otra vez: una ventana corta hace parecer incoherentes a las d altas.
- **No se cambia ningún umbral.** Sí se fija una regla de uso (§6).

**W2. Árboles de bloques grandes: escapan a tamaño finito.**
- El árbol de cubos 5³ tiene δ = 0.30 (parece 3D) y está CONEXO en los dos tamaños.
- Su f̃ de anillos baja de 1.000 a 0.955, justo sobre el umbral 0.95, y γ baja de 0.85 a 0.69.
- Las tendencias apuntan a la ramificación, pero con 80–640 bloques todavía no llega.
- Un árbol de bloques solo es degenerado a escalas mucho mayores que el bloque. Es un **límite de resolución**, análogo a O3, que el prerregistro no incluyó.

**W3. Margen estrecho de X1.**
- La apoloniana tiene δ = 0.092 / 0.130 según la semilla. La teoría del logaritmo puro daba 0.098, pero las constantes aditivas en R suben δ a tamaño finito.

**W4. El límite R_MAX = 200 del observador congelado satura R** en la escalera, el cilindro C × C₅ y el árbol de Prüfer en 8N.
- δ deja de ser un exponente allí.
- La decisión no cambia, porque X2 y X3 también se activan.
- Además, el RGG1 se fragmenta, con δ = −9.67, y se excluye por X3.

## 5. Contraste con las predicciones (§4)

| Predicción | Resultado |
|---|---|
| FE: 0–1 en d = 5–6 | ✔ 0 |
| Heisenberg excluido por X4 | ✘ (INC → INTERMEDIO) |
| Las de diámetro logarítmico, por X1 | ✔ salvo la apoloniana s0 (repetida, δ = 0.130) |
| Árboles por X2 | ✔ Prüfer, binario y cactus; ✘ árbol de cubos (no excluido); árbol × C₁₀ por X4 |
| Retazos por X4 | ✘ (SIN_VENTANA en N) |
| 1D por X3 | ✔ |
| **RGG3 + 0.1 % de atajos no excluido** | ✘: **sí se excluye**, por X4 |
| WS por X1 | Excluido, pero por X4 (δ 0.15–0.17) |
| Veredicto PARCIAL | ✔ con la lectura conservadora (✘ con la literal) |

## 6. Dictamen del cerebro

1. **RC-3 se adopta como juez A− provisional**, oficialmente PARCIAL. Es el primer juez del programa con **FE = 0 en d = 2…6**, es decir, sin sesgo dimensional demostrable en el rango validado.
2. **Regla de uso**, que no altera umbrales:
   - toda aplicación a una dinámica informa **R, δ, γ y f̃ en los dos tamaños**;
   - una exclusión que dependa **solo de X4** se marca como **«X4-marginal»** y no cuenta como exclusión firme;
   - todo «NO-EXCLUIDO» se inspecciona (O2).
3. **Catálogo de puntos ciegos declarados:**
   - árboles de bloques con bloques grandes respecto de N (W2);
   - triangulaciones apiladas cerca del corte de δ (W3);
   - retazos sin ventana en N.
4. **Mapa final de exclusiones:** `docs/OMEGA_MAPA_EXCLUSIONES.md`.
