# Ω — L-A-R1-1: dinámica mínima SQ — Prerregistro

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-r1`.
- **Base:** `claude/omega-r1-0-congelado` (análisis R1-0).
- **Mandato del Consejo:** R1-1 aprobada, formulada así:

  > «SQ es un experimento de falsación sobre si una regla local, sin conservar grado, densidad, ciclos ni dimensión, puede generar espontáneamente una estructura estable con riqueza cíclica. No es todavía una hipótesis de emergencia de 3D.»

- **Correcciones del Consejo incorporadas:**
  - SQ es el **primer candidato falsable**, no «la regla mínima correcta»;
  - se comprueban componentes, extremos y comportamiento a escala creciente;
  - el éxito es «una fase estable que sobrevive a los adversarios», no «aparecieron cuadrados».
- **Estado:** este documento se commitea antes de escribir el código de la dinámica.

## 1. Espacio de reglas (punto 1 del auditor), cerrado: 3 × 2 × 3 = 18 variantes

Grafo simple sin pesos con **N fijo**. Un barrido consta de una **fase de borrado** seguida de una **fase de reenganche**.

**Disparador de borrado (B):**
- **B-s:** se borra la arista uv si s(u,v) = 0, es decir, si no está en ningún 4-ciclo. Se usa s(u,v) = (A³)_{uv} − deg u − deg v + 1.
- **B-t:** se borra si t(u,v) = 0 (no está en ningún triángulo). Es la variante de contraste: no deja fijo a Z^d.
- **B-st:** se borra si s(u,v) = 0 **y** t(u,v) = 0 (no está en ningún ciclo corto de longitud 3 o 4).

**Reenganche (R):** todo vértice con grado 0 o 1 añade **una** arista.
- **R2:** si el grado es 1, a un vértice uniforme entre los que están a distancia exactamente 2.
- **R3:** si el grado es 1, a uno uniforme entre los que están a distancia exactamente 3.
- **Si no existe candidato, o si el grado es 0:** a un vértice uniforme del grafo, distinto de sí mismo y no vecino. Es un **forzamiento no local declarado**.

**Orden de actualización (O):**
- **SYNC:** se evalúan todos los disparadores de borrado sobre el grafo actual y se aplican juntos. Después se evalúan todos los grados y se aplican los reenganches. Sus destinos se sortean sobre el grafo posterior al borrado, en orden aleatorio, sin reevaluar.
- **ASYNC:** se recorre una permutación aleatoria de las aristas actuales, reevaluando el disparador sobre el grafo vigente. Después, una permutación aleatoria de los vértices, reevaluando el grado.
- **RAND:** E + N eventos por barrido. Cada evento elige una arista uniforme (con probabilidad E/(E+N)), o un vértice uniforme, y aplica su regla sobre el grafo vigente.

## 2. Invariantes (punto 2)

- Solo N.
- E, el grado, los ciclos y las componentes son **salidas**.
- Ninguna regla contiene d, k, ρ ni L, ni ninguna salida de instrumentos (G1).

## 3. Condiciones iniciales (punto 3), de clases estructuralmente distintas (G3)

| Id | Inicio | Tipo |
|---|---|---|
| I1 | «Completo con pesos»: W_ij ~ U(0,1), se conservan las aristas con W_ij > 1 − 20/(N−1) (grado medio ≈ 20) | Amorfo denso |
| I2 | ER con grado medio 4 | Amorfo disperso |
| I3 | Árbol de Prüfer uniforme | Amorfo arbóreo |
| I4 | Anillo C_N | Amorfo 1D |
| I5 | Toro Z³ con el lado más cercano a N^{1/3} | **Control, no amorfo**: comprueba que Z³ es punto fijo con B-s y B-st. **No cuenta** como emergencia ni para el criterio de atractor |

## 4. Tamaños, semillas y pasos (puntos 4–6)

| Etapa | Contenido |
|---|---|
| **A** (obligatoria) | N₀ = 5000 y 8N₀ = 40 000, que es el par que exige RC-3 (I5: lados 17 y 34). Semillas 0–2. Las 18 variantes × 5 inicios ⇒ 540 corridas |
| **B** (condicional, fijada ya) | Para todo par (B, R) con algún NO-EXCLUIDO en la etapa A, se añaden 2N₀ y 4N₀ y las semillas 3–4 |

- **Desviación respecto de G4 de R1-0** (semillas 0–4 y los cuatro tamaños): el motivo es el coste, sin numba (~1800 corridas). Solo amplía, nunca recorta: es conservadora.
- **Pasos:** T = 100 barridos como máximo. Se detiene antes por **absorción**: ningún disparador activo y ningún vértice de grado ≤ 1.
- La clave de RNG es (20261018, id de variante, id de inicio, semilla, índice de tamaño), solo con `rng_from_key`.

## 5. Métricas (punto 7)

**Durante la dinámica** (cada 5 barridos; solo se registran y **la dinámica no las ve**, G1):
- E/N;
- fracción de aristas con s ≥ 1;
- fracción de aristas con t ≥ 1;
- número de componentes;
- fracción de la componente gigante;
- fracción de vértices de grado ≤ 1;
- número de borrados y reenganches del barrido.

**Al final:** el grafo final de cada corrida.

**No se mide d_eff durante la dinámica.** δ solo se calcula dentro de RC-3, después, y nunca como criterio de dimensión (G5).

## 6. Clasificación ciega (puntos 8 y 10), por par (N₀, 8N₀) de la misma variante, inicio y semilla

1. **X0 FRAGMENTADO:** la componente gigante es < 0.5 N en cualquiera de los dos tamaños. Clase vacío o bosque; RC-3 no se aplica.
2. Si no, **RC-3** congelado (`tools/rc3.py`: X1–X4, con la regla de uso de exclusión X4-marginal). Se informan R, δ, γ y f̃ en los dos tamaños.
3. Resultado por par:
   - FRAGMENTADO;
   - EXCLUIDO(Xi);
   - EXCLUIDO-X4-marginal (no cuenta como exclusión firme);
   - NO-EXCLUIDO.
4. **Todo NO-EXCLUIDO, y todo X4-marginal, se inspecciona** con el censo t/s y la distribución de grados. Se informa como **residuo no clasificado** si no corresponde a ninguna clase conocida (G7).

## 7. Criterio de atractor (punto 9, G6) y nulo del espacio de reglas (Ω6)

- **Clase dominante de un par (B, R):** la clase más frecuente sobre 3 órdenes × 4 inicios amorfos (I1–I4) × semillas (36 pares en la etapa A).
- **Atractor:** clase dominante con frecuencia ≥ 0.8. Si no la hay, **MULTIESTABLE**, y se informa la distribución.
- **Resultado positivo de R1-1:** algún par (B, R) cuyo atractor sea NO-EXCLUIDO (y no marginal) en la etapa A, y que lo mantenga en la etapa B (2N₀, 4N₀ y semillas 3–4).
- **Nulo Ω6:** fracción de los 6 pares (B, R) con resultado positivo. Un positivo aislado se informa con esa tasa.

## 8. Resultados ambiguos (punto 11)

- MULTIESTABLE, X4-marginal, una dinámica sin absorción con oscilación, o una clase desconocida ⇒ **residuo no clasificado**.
  - No se interpreta geométricamente.
  - Vuelve a análisis conceptual.
- Ninguna variante, parámetro ni umbral se añade o cambia después de ver resultados.

## 9. Predicciones del cerebro (congeladas)

1. **B-s y B-st desde I1–I4.**
   - En el primer barrido se borran casi todas las aristas: ER, árbol y anillo no tienen 4-ciclos; el inicio denso I1 sí tiene muchos.
   - El reenganche recrea un grafo aleatorio disperso (≈ N aristas, sobre todo no locales), que vuelve a borrarse.
   - Resultado: **rotación sin absorción**. Clase final **FRAGMENTADO** o **X1** (grafo aleatorio arbóreo), según la fase del ciclo en que acabe.
   - En I1 (denso), el borrado no es total, pero el residuo es aleatorio: **X1**.
2. **B-t:** lo mismo, y además destruye I5 (Z³ no tiene triángulos).
3. **I5 (control):** fijo bajo B-s y B-st, que es la sanidad del código. Destruido bajo B-t.
4. **Nulo Ω6:** 0/6 pares positivos.
5. **Veredicto previsto: R1-1 NEGATIVA (SQ degenerada en todas las variantes).**
   - El fallo será de **nucleación**: la frecuencia de cierre de cuadrados por reenganches es ∝ 1/N.
   - Según la rama prerregistrada en R1-0 §6, el siguiente candidato sería R3, con auditoría E6.
