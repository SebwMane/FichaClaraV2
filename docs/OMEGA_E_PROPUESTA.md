# Ω-E: ruta entrópica (caso 3 de Ω-0.2) · preregistro de E0 y E1 (2026-10-12)

- **Autorizado por el usuario** (2026-10-12): «preregistra y empieza a probarlo».
- **Pregunta de fondo.** ¿La ley de área de la entropía funciona como criterio para reconocer geometría sin competencia entre estructuras? Ruta de Jacobson: S ∝ área es la condición que el sustrato debe aportar (necesaria, no suficiente).
- **Método:** el de siempre (preregistro, código, auditoría, corrida limpia, resultados). Escala de afirmación en 4 niveles.

## 1. Instrumento (Sorkin y Yazdi 2018, arXiv:1611.10281)

- Orden con matriz causal C (C_xy = 1 si x ≺ y). Campo escalar sin masa en 2D: G_R = C/2, Δ = (C − Cᵀ)/2, iΔ hermítica.
- Vacío de Sorkin-Johnston en la región grande: W = parte positiva de iΔ.
- Entropía de una subregión R: λ = autovalores de (iΔ_R)⁻¹ W_R sobre la imagen de iΔ_R; S = Σ λ ln|λ|.
- **Sin truncar:** se excluye solo el núcleo numérico de iΔ_R (|ν| ≤ 10⁻⁸ · max|ν|).
- **Truncado (doble):** en la región grande se eliminan de iΔ y de W los modos con |μ| < √N_L/(4π); se restringe a R; en iΔ_R se eliminan los modos con |ν| < √n_R/(4π) y se proyecta W_R sobre los restantes.

## 2. E0: reproducir el resultado publicado · PROBAR (vota)

- Diamante 2D (el de C0), subdiamante centrado de lado mitad (|u|, |v| ≤ 1/4 en coordenadas de cono de luz), N_L ∈ {200, 400, 800, 1200, 1600}, 5 semillas por N (entropía 20261047).
- **E0a (sin truncar):** ajuste lineal S = a n + b con a ∈ [0.25, 0.40] (publicado 0.32) **y** R² lineal > R² logarítmico.
- **E0b (truncado):** ajuste S = a ln(√n/4π) + b con a ∈ [0.23, 0.43] (publicado 0.346 ± 0.028; continuo 1/3) **y** R² logarítmico > R² lineal.
- **Fracaso:** se audita el instrumento; E1 no se ejecuta.

## 3. E1: ¿el criterio separa geometría de las familias que la vencieron? · PROBAR (vota)

**Regla de región intrínseca (congelada; solo usa el orden).** R = intervalo I[a, b] (con extremos) de tamaño más cercano a N/4; empate: el par (a, b) de menores índices. La región grande es todo el orden.

**Familias** (N_L ∈ {200, 400, 800, 1200, 1600}, 5 semillas, entropía 20261048):
- control geométrico: sprinkling 2D (misma regla intrínseca, no la del subdiamante centrado);
- KR (capas N/4, N/2, N/4, p = 1/2);
- percolación transitiva p = 0.05;
- capas gruesas (8 capas, p = 1/2);
- bicapa aleatoria y sándwich (capas (N−2)/2, 2, (N−2)/2, p = 1/2).

**Bandera «ley de área»:** con el doble truncamiento, R² logarítmico > R² lineal **y** a ∈ [0.23, 0.43], ajustando S frente a n = |R| en las 25 corridas de la familia.

**Familias sin regiones extensas.** Si en N_L = 1600 el intervalo elegido tiene menos de 20 elementos en todas las semillas, la familia se registra como «no aplicable: sin regiones extensas» y cuenta como **sin** ley de área (el criterio exige una región donde medir).

**Decisión:**
- **DISCRIMINA:** el control geométrico tiene ley de área y ninguna familia no geométrica la tiene.
- **NO DISCRIMINA:** alguna familia no geométrica tiene ley de área.
- **INDETERMINADO:** el control geométrico no tiene ley de área con la regla intrínseca.

**Declaraciones.** (a) G_R = C/2 y el umbral √n/(4π) vienen de un campo sin masa en 2D; aplicarlos a órdenes no geométricos es una convención congelada aquí, no un resultado. (b) El umbral se fijó antes de calcular cualquier familia; no se ajustará. (c) Un DISCRIMINA sería nivel 1 para N ≤ 1600 en 2D, frente a estas familias; no prueba la ruta de Jacobson completa.
