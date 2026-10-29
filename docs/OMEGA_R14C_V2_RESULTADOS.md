# R14-C v2: observables causales multivariables. Resultados (2026-10-20)

- **Preregistro:** `OMEGA_R14C_V2_PREREGISTRO.md`, con la aclaración A1 (§5).
- **Código:** commit a0ae61f, árbol limpio. Son 2004 unidades en 3 tandas, todas con el mismo commit.
- **Datos:** `results/omega_r14c/` (`tabla_r14c.csv` con 784 filas y `summary.json`). **Figura:** `docs/figuras/r14c_v2.png`.
- **Tipo:** exploratorio, **sin voto**. Todas las estructuras son vistas. Lo que aparezca aquí es hipótesis, no criterio.

## 1. Calibración

- **Ninguna geometría sale significativa contra su propia referencia** (conformes, Kasner y lineales, con N = 4000 y 8000). La tabla está bien calibrada.
- **Chequeo de la fórmula de cadenas contra los diamantes de Monte Carlo:** |z| ≤ 1.9 en 15 de 16 celdas. La excepción es k = 5 con m = 80 (z = −3.0), donde C₅ es muy asimétrica (su mediana es 0). Ya se había anticipado al escribir el código.

## 2. Qué familias se separan de la plana (Holm < 0.05 sobre las 784 filas)

| Familia | Observables que separan (de 49, con N = 4000 / 8000) | Lectura |
|---|---|---|
| Percolación | 49 / 49 | separación total |
| Producto 4 | 41 / 43 | casi total |
| Poliédricos (cúbico y octaédrico) | 37 / 41 | casi total, también en la pendiente de curvatura A |
| Conos suaves (9 estructuras) | 1 / 8 | se separan solo con N = 8000, con efectos moderados (\|Δ\| de 0.5 a 2.7) |
| **Defectos (C5, D4 y D5)** | **0 / 0** | efectos enormes pero heterogéneos; no sobreviven a Holm (ver §3) |

- **Conos suaves con N = 8000:** separan el máximo del espectro del Laplaciano (Δ = −2.0), las abundancias con m = 20, el cociente de la cadena más larga con m = 20 y 160, las 3 y 4 cadenas con m = 160, y r con m = 160. **En todas, la geometría más cercana es de Sitter o de Sitter fuerte**, salvo la de densidad suave C2. La curvatura vuelve a imitar a los conos, como pasó con el perfil en F7 (pregunta P2 de v1).
- **Lo nuevo frente a F6:** el **máximo del espectro del Laplaciano** (el observable de Eichhorn et al.) es el único que separa a los conos suaves con N = 4000 y con N = 8000. Es el candidato más fuerte, pero queda solo como hipótesis.

## 3. El caso de los defectos: una debilidad de diseño declarada

- La familia de defectos tiene efectos medios grandes: por ejemplo, Δ = 53 en el espectro del Laplaciano, Δ = 13 en r con m = 160 y Δ = 23 en el grado de links. Pero sus 3 estructuras son muy distintas entre sí (C5 borra links, D4 agrega relaciones, D5 agrega gemelos).
- Con solo 3 estructuras, el bootstrap jerárquico por estructuras da intervalos muy anchos (por ejemplo, de 5.5 a 114 en el espectro del Laplaciano), y ningún p-valor sobrevive a Holm sobre 784 pruebas.
- **No significa que los defectos sean invisibles:** R14-B mostró C5, D4 y D5 muy fuera de banda con q = 1. Significa que el promedio por familia con 3 miembros heterogéneos no es una unidad útil. La regla de R14-A §7 (la familia como unidad para afirmaciones generales) es correcta para generalizar, pero cuando la familia mezcla mecanismos distintos pierde toda la potencia.

## 4. Respuestas exploratorias a las preguntas de v1 §2

- **P1** (¿qué separa a D6, D7 y D8?): a nivel de familia, el máximo del espectro del Laplaciano y los observables de sub-intervalos grandes (m = 160) con N = 8000. Por estructura no se votó.
- **P2** (¿la curvatura imita?): **sí.** Para los conos suaves, la geometría más cercana es de Sitter en la mayoría de los observables que los separan de la plana.
- **P3** (¿se detectan los gemelos sin un detector a medida?): a nivel de familia no se puede afirmar, por el problema de §3. Por estructura habría que mirarlo con una prueba nueva.
- **P4** (¿una combinación mínima?): no se elige ninguna ganadora (v2 §3). El candidato que surge para preregistrar es el **máximo del espectro del Laplaciano más un observable de sub-intervalos grandes**, contrastado contra referencias curvas y no solo contra la plana.

## 5. Qué implica

1. **Para conos «redondos casi en todas partes»,** la información existe con N = 8000, pero es débil y se confunde con curvatura positiva. Cualquier criterio futuro debe compararse con la familia de geometrías curvas, no solo con la plana.
2. **Para afirmaciones por familia,** las familias deben agrupar un mismo mecanismo. Los defectos deberían separarse en «quitar», «agregar» y «duplicar».
3. **Confirmar algo de esto** exige familias nuevas, un guardián nuevo y un preregistro con el candidato congelado (v2 §3).
