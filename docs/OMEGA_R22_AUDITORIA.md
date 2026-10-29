# R22: gate de auditoría de Φ_q después de R20 (fase 3 del dictamen de la sesión 11), 2026-10-26

- **Mandato:** dictamen del Consejo, sesión 11 (`OMEGA_CONSEJO_REV22_R20.md`), fase 3; el usuario pidió adoptar los pasos prudentes.
- **Tipo:** auditoría descriptiva (nivel 2). **No cambia el veredicto de R20 ni sus datos.** No usa ninguna semilla confirmatoria de R20 (20261113, 20261115 ni 20261116) para generar datos.
- **Código y datos:**
  - `tools/r22_auditoria.py` → `results/omega_r22/auditoria.json` (solo lee `results/omega_r20`);
  - `tools/r22_seleccion.py` → `results/omega_r22/seleccion/` (regenera la selección del guardián con su esquema publicado, entropía 20261114: 400 referencias y 200 realizaciones, unos 60 min de pared con 4 núcleos);
  - `results/omega_r22/seleccion/sensibilidad_banda.json` (remuestreo de bandas, `default_rng(20261117)`; próxima entropía libre: 20261118).

## 1. Dimensión y máscara

| q | 1 | 1/2 | 1/4 | 1/8 |
|---|---|---|---|---|
| Componentes en la banda del ejecutor | 14 | 14 | 14 | **13** |
| r₄₀ no medible (ajuste / calibración / validación, de 100) | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | **37 / 49 / 44** |
| Máscara si se fija solo con el conjunto de ajuste | igual | igual | igual | igual |

- **14 era un inventario nominal, no un requisito matemático.** La regla «componentes finitas en todas las referencias» venía de R14-B y R19, y estaba escrita en el ejecutor antes de correr.
- **Por qué r₄₀ sale con q = 1/8:** con N ≈ 1000, en el 43% de los diamantes planos m = 40 supera el límite de localidad 0.1·m_max, así que el observable no se mide. No es un faltante ocasional; es estructural a esa escala.
- **Faltantes en estructuras:** con q = 1/8, ningún brazo perdió una componente de la máscara, salvo uno: r₂₀ en una realización del control de G1a. Por eso G1a quedó no evaluable.

## 2. Faltantes y comparabilidad entre q

Con un vector común de 13 componentes en los 4 q:
- la cobertura de validación queda en 91, 92, 99 y 95 de 100 (todas ≥ 87);
- L₁/₈ no cambia, porque es la misma máscara;
- las medianas de L con q = 1, 1/2 y 1/4 cambian como mucho 0.03 en las 20 estructuras.

**Conclusión:** la máscara variable no afectó nada material en R20. Para estudios nuevos conviene fijar el vector común desde el principio, porque permite comparar D entre q sin esa salvedad.

## 3. Garantía bajo el selector real

- La máscara es una función **simétrica** de las 300 referencias de cada q (no depende de qué referencia esté en qué conjunto). Dada la máscara, calibración y validación siguen siendo intercambiables, así que **el diagnóstico de cobertura (87/100) es válido tal como se aplicó.**
- En los datos realizados, la máscara con solo el conjunto de ajuste es idéntica, así que la elección no dependió de calibración ni de validación.
- **Empates:** las distancias son continuas; los empates tienen probabilidad cero y la regla por índice nunca decidió nada.
- **Lo que no garantiza:** la cobertura vale para puntos planos con todas las componentes de la máscara medibles. Para un orden no plano no hay garantía de cobertura; la banda solo define qué significa «en banda».

## 4. Reproducción de la selección del guardián (hallazgo principal)

**4.1 La reproducción es exacta.** Con su esquema de semillas y su métrica (banda propia y D por par), recuperé sus radios al último decimal (4.823196459419363 con q = 1 y 4.439727163734706 con q = 1/8) y las 80 medianas del recibo con una diferencia máxima de 5·10⁻⁵ (redondeo del recibo), incluidas las 4 S1 descartadas.

**4.2 Con la métrica del ejecutor cambian dos estructuras.** Mismas realizaciones, banda conforme de R20:

| Estructura | Requisito | Guardián | Ejecutor | ¿Identificable? (G / E) |
|---|---|---|---|---|
| S2b, prisma de Reuleaux | D_obs/radio con q = 1/8 ≥ 1.25 | 1.38 | **1.17** | sí / **no** |
| S2d, prisma lenticular | D_obs/radio con q = 1/8 ≥ 1.25 | 1.31 | **1.21** | sí / **no** |
| las otras 18 | | | | sin cambios |

Las S1 siguen sin ser identificables con las dos métricas (su D con q = 1/8 está entre 0.65 y 0.90).

**4.3 La causa es el ruido de la banda, no la máscara ni el D por par.** Las dos métricas usan 13 componentes con q = 1/8 y ningún par perdió componentes, así que el D por par coincide con el D de la máscara. Lo que cambia es la muestra plana con la que se ajusta la banda. Para medirlo, rehíce la banda 300 veces con 100 de ajuste y 100 de calibración sorteadas de las 500 referencias planas disponibles por q (200 del guardián y 300 del ejecutor), sobre las mismas realizaciones de selección:

| Estructura | P(identificable) con banda remuestreada |
|---|---|
| G1b a G1d, G2c, G2d, P1a, P1c, P1d, S2a y S2c | 1.00 |
| P1b | 0.99 |
| G2b | 0.85 |
| **S2b** | **0.67** |
| **S2d** | **0.39** |
| G1a, G2a y las 4 S1 | 0.00 |

- Votarían las 14 en el 38% de las bandas; la mediana es 13 y el percentil 5 es 11.
- **El mínimo de población se cumple con cualquier banda:** S2a y S2c pasan siempre, así que S2 nunca baja de 2 identificables.
- **Las S1 no se rescatan con ninguna banda.** El reemplazo S1 → S2 no fue un efecto de la métrica.

**4.4 Qué significa para R20.**
- **El veredicto no cambia.** PASA exigía que todas las votantes cumplieran, y las 14 cumplieron. Cualquier subconjunto de votantes también cumple. Con la métrica del ejecutor votarían 12 y el veredicto sería el mismo.
- **El alcance sí se precisa, sin rescatar nada con resultados ya vistos:** la afirmación robusta a la banda vale para **12 estructuras** (G1b a G1d, G2b a G2d, P1a a P1d, S2a y S2c). S2b y S2d votaron con la métrica sellada y quedaron equivalentes, pero su identificabilidad depende de la muestra plana. Sobre conos suaves, el resultado firme es para un cono circular y un prisma semicircular.
- **Es una limitación del diseño, no un error de ejecución.** El guardián y el ejecutor usaron dos bandas independientes, como pedía la separación de agentes. Los umbrales de visibilidad (1.25 y 1.5) no tenían margen para la incertidumbre de la banda.

## 5. El borde cercano a 0.8

| q | Mediana de D/radio de la plana (validación) | Referencia χ² (mediana / cuantil 96/101) |
|---|---|---|
| 1 | 0.71 | 0.75 |
| 1/2 | 0.75 | 0.75 |
| 1/4 | 0.70 | 0.75 |
| 1/8 | 0.72 | 0.74 |

- **El «piso» es la distribución nula.** Un punto plano típico tiene D/radio cerca de 0.72; el borde de banda (1) es su cuantil 96/101. No hay saturación técnica: los D son continuos, no se amontonan en un valor y su dispersión entre realizaciones es comparable a la de la plana.
- **Las copias con q = 1/8 no llegan al piso:**

  | Votante | D_obs mediana | En banda | Mann-Whitney contra la plana |
  |---|---|---|---|
  | G1b, G1c y G1d | 0.79, 0.87 y 0.84 | 83%, 83% y 77% | p = 0.021, < 0.001 y 0.008 |
  | G2b, G2c y G2d | 0.79, 0.85 y 0.91 | 87%, 77% y 63% | p = 0.037, < 0.001 y < 0.001 |

  Quedan cerca de la distribución plana, pero todavía se distinguen de ella. **Corrige mi lectura 1 de la sesión 11**, que decía «indistinguible del plano».
- **L queda acotado por el piso:** si D_obs fuera exactamente un punto plano típico, L sería ln(0.72/D_ctrl). Para las 6 votantes, eso da valores entre 0.06 y 0.23 más negativos que el L medido. El L observado está cerca de ese límite, así que con q = 1/8 depende sobre todo de D_ctrl.
- **La no monotonía de L entre q = 1/4 y 1/8** (G1b, G1c) viene del control: D_ctrl baja de 2.57 a 1.71 en G1b porque el control también pierde resolución al bajar N, mientras D_obs ya está cerca del piso.

## 6. Salida del gate

| Opción del dictamen | Juicio |
|---|---|
| (i) Garantía formal para el procedimiento actual | **Sí para la cobertura** de puntos planos (§3). **No para la selección por visibilidad:** depende de la muestra plana cerca de los umbrales (§4). |
| (ii) Nueva calibración y validación con instrumento congelado | **Recomendada antes de cualquier estudio confirmatorio nuevo que invoque Φ_q.** |
| (iii) Solo exploratorio | Suficiente para el piloto (a), que es exploratorio por diseño. |

**Recomendación:** (i) para la cobertura y el alcance de R20 tal como se precisó en §4.4; (ii) como requisito de entrada de cualquier confirmatoria futura. Para (ii), el instrumento congelado debería tener:
1. **una sola banda** compartida por la selección y la confirmación, ajustada antes y publicada con su hash;
2. **el vector común** de 13 componentes en todos los q;
3. **más referencias de ajuste**, o umbrales de visibilidad con un margen para la variación de la banda: por ejemplo, exigir que la estructura pase en al menos el 95% de las bandas remuestreadas, como en §4.3;
4. **un estadístico de exceso sobre la nula** si se quiere medir velocidad (§5).
