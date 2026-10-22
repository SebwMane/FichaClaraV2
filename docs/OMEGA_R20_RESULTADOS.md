# R20: repetición corregida de R19-(b), alcance A. Resultados (2026-10-21)

- **Preregistro:** `OMEGA_R20_REPETICION_PREREGISTRO.md`, versión 2, firmada.
- **Sellado:** `docs/SELLADO_R20.md` (sha256 f6dc1a12…), con el recibo en `docs/sellado_r20/`.
- **Cadena de custodia:** hash c338b9c, ejecutor 3cda6bd y apertura c2e227d. La corrida usó c2e227d con el árbol limpio: 1830 unidades (1200 referencias y 630 realizaciones de estructuras), en 2 tandas.
- **Datos y logs:** `results/omega_r20/` (incluye `entorno.json` y `logs/`). **Figura:** `docs/figuras/r20.png`.

## 1. Veredicto oficial: **PASA** (14 estructuras votantes)

| Condición | Resultado |
|---|---|
| Cobertura de validación (≥ 87/100 en cada q) | 92, 92, 100 y 95 de 100: **cumple** |
| Calibración nula, L₁/₈ dentro de ±0.262 | plana 0.10 (−0.08 a 0.17) y de Sitter 0.03 (−0.05 a 0.21): **cumple** |
| Población (≥ 2 familias por tipo, con ≥ 2 identificables) | gemelos y tripletes; conos poliédricos y suaves: **cumple** |
| «Se desvanece» (6 votantes): IC por debajo de −0.262 | las 6 cumplen; L₁/₈ de −0.80 a −1.41 |
| «Autosemejante» (8 votantes): IC dentro de ±0.262 | las 8 cumplen; los IC van de −0.15 a 0.15 |

## 2. Datos por estructura (L₁/₈, mediana e IC 95%)

| Estructura | Estado | L₁/₈ |
|---|---|---|
| G1b, G1c y G1d (gemelos, 12%, 17% y 23%) | votan | −0.80, −1.01 y −1.41 (IC entre −1.51 y −0.74) |
| G2b, G2c y G2d (tripletes, 12%, 17% y 23%) | votan | −0.91, −1.11 y −1.32 (IC entre −1.40 y −0.75) |
| P1a a P1d (conos poliédricos) | votan | de −0.02 a 0.01 (IC entre −0.07 y 0.09) |
| S2a a S2d (conos suaves) | votan | de −0.07 a 0.08 (IC entre −0.15 y 0.15) |
| G1a (gemelos al 7%) | **no identificable**; además, **no evaluable** con q = 1/8 (r₂₀ no medible en el control de una realización) | |
| G2a (tripletes al 7%) | **no identificable** | −0.54 (−0.69 a −0.41); se informa, no vota |
| D5, D7 y D8 (vistas) | no votan | D5: −0.63; D7: −0.04; D8: 0.06 |

En la figura, los conos tienen las curvas del observado y del control superpuestas en los 4 q, y las copias tienen el observado muy por debajo de su control desde q = 1/2.

## 3. Qué autoriza este PASA, y qué no

- **Autoriza (nivel 1, sobre el instrumento):** en las estructuras sintéticas identificables del universo sellado, Φ_q aleja las copias de su control a densidad fija, y mantiene compatibles con su control los conos autosemejantes, con la métrica fijada (vector de 14 componentes, Mahalanobis con Ledoit-Wolf y banda conforme). **El control reinyectado separa un defecto que se desvanece de uno que solo deja de verse.** Era el problema de método que había mostrado R14-B.
- **No autoriza:**
  - que Ω produzca estas estructuras;
  - nada sobre H_E, H_L ni H_sel;
  - extender el resultado a conos que el instrumento no ve.
- **Límites declarados:**
  1. **Selección por visibilidad a nivel de familia.** La primera familia de conos suaves del guardián tuvo 0 de 4 estructuras visibles y fue reemplazada una vez (registrado en el recibo). Para conos suaves, la población es la de los visibles; los descartados no son identificables por este instrumento.
  2. **Clasificación del guardián.** La familia «conos suaves» incluye prismas de sección curva (Reuleaux, semicircular, lenticular), que tienen caras planas en una dirección.
  3. **Universo sintético.** 16 estructuras de 4 mecanismos, con N₀ = 8000 y 30 realizaciones por estructura. No hay garantía de potencia para otros mecanismos.
  4. **FALLA** no tenía control familiar del error, pero no se aplicó. **PASA** se obtuvo por intersección-unión: cumplieron las 14 votantes.

## 4. Relación con R19-(b)

R19-(b) sigue INDETERMINADO, y este resultado no lo modifica. R20 corrige sus dos causas (la cobertura y las estructuras invisibles) y además el error de semillas de R19 (errata en `OMEGA_R19_GEMELOS_RESULTADOS.md` §5). El patrón de nivel 2 que había mostrado R19 se confirma aquí de forma preregistrada.
