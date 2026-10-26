# Ω: síntesis del programa después de R20 (2026-10-26)

- **Mandato:** fase 2 del dictamen del Consejo, sesión 11 (`docs/OMEGA_CONSEJO_REV22_R20.md`). **Aprobada por el usuario el 2026-10-26 (R22-3).**
- **Alcance:** la fase causal del programa, desde la Revisión 14 (2026-10-18), cuando se separaron H_rec, H_E, H_L y H_sel. Lo anterior está resumido en `OMEGA_ESTADO.md` (fase de grafos, problema B cerrado), `OMEGA_MAPA_NEGATIVO.md` y `OMEGA_2_SINTESIS_2D.md` (Ω-C en dos dimensiones).
- **Regla de lectura:** cada afirmación lleva su escalón (§1) y su documento. Lo que no está en la tabla de hechos (§2) no está establecido.

## 1. Escala de afirmación

Dos ejes independientes. El primero dice **sobre qué** se afirma; el segundo, **con qué fuerza** (niveles del proyecto: 1 medido, 2 explicación, 3 teorema externo verificado, 4 extrapolación).

| Escalón | Sobre qué | Ejemplo | Qué hace falta para subir al siguiente |
|---|---|---|---|
| **E0** Observación del código | El código calcula lo que dice | balance detallado por enumeración; S_ε con √ε corregido | un instrumento que use ese código y se valide contra casos conocidos |
| **E1** Evidencia del instrumento | Un instrumento distingue A de B en un universo declarado | R20: Φ_q con control reinyectado separa copias de conos visibles | aplicar el instrumento a órdenes cuyo origen no se conoce de antemano |
| **E2** Hipótesis física | Una propiedad de los órdenes (H_E, H_L) | «un orden no geométrico se vuelve geométrico al diluirlo» | una dinámica de Ω que produzca esos órdenes con peso medible |
| **E3** Generación o selección por Ω | La dinámica de Ω asigna peso a esas fases (H_sel) | ninguno | |

**Estado actual:** el programa tiene resultados en E0 y E1. **No tiene ningún resultado confirmatorio en E2 ni en E3.**

## 2. Hechos establecidos

| # | Hecho | Escalón y nivel | Universo y límites | Fuente |
|---|---|---|---|---|
| H1 | Φ_q (dilución) más el control reinyectado separa defectos que se desvanecen (copias) de defectos autosemejantes (conos), con la métrica ejecutada | E1, nivel 1, **confirmatorio** (PASA) | 14 estructuras sintéticas votantes de 4 familias visibles (G1, G2, P1, S2); N₀ = 8000; q = 1/8 con 13 componentes; conclusión limitada a S2 para los conos suaves; **robusta a la banda para 12** (S2b y S2d tienen identificabilidad frágil) | `OMEGA_R20_RESULTADOS.md` y su adenda A1 |
| H2 | Diluir confunde convergencia con pérdida de resolución: al bajar N, las bandas se ensanchan y un defecto fijo termina en banda | E1, nivel 1, exploratorio | estructuras vistas, N₀ = 8000 | `OMEGA_R14B_V2_RESULTADOS.md` |
| H3 | Percolación, producto y poliedros se separan de la plana con observables causales; los conos suaves, débilmente y confundidos con de Sitter | E1, nivel 1, exploratorio | N = 4000 y 8000; Holm sobre 784 filas | `OMEGA_R14C_V2_RESULTADOS.md` |
| H4 | El instrumento E1′ (entropía SJ) es de 2D: se estabiliza solo con N ≥ 512; KR imita el coeficiente con 2 modos; en 3D y 4D no da ley de área | E1, nivel 1 | sprinklings y KR, N de 64 a 2048 | `OMEGA_A_E1P_RESULTADOS.md` |
| H5 | El S_max (intervalo máximo) no sirve como instrumento: puntúa más alto la trampa de altura 2 que la geometría | E1, nivel 1 | N ≤ 64 | `OMEGA_R17_1_RESULTADOS.md` |
| H6 | La aceleración exacta de E1′ reproduce el cálculo completo (error ≤ 8.3e-12) con ganancias de 1.4 a 5.7 veces | E0 | N de 32 a 512 | `OMEGA_R17_3_RESULTADOS.md` |
| H7 | La acción suavizada 4D lleva el factor √ε (Dowker y Glaser 2013, ec. 42); los resultados anteriores eran S_ε/√ε y sus cocientes no cambian | E0 | todo el repositorio | `omega_c/AUDITORIA_SQRT_EPS.md` |
| H8 | Con el movimiento local, las cadenas de Metropolis sobre órdenes generales no mezclan con N ≥ 64, ni siquiera con la medida de conteo (b = 0); Wang-Landau en la altura no lo arregla; la reinserción R no es aplicable | E1 sobre el muestreador, nivel 1 | N de 32 a 128; Ω-E y Ω-C | `OMEGA_ES0_RESULTADOS.md`, `OMEGA_ESC0_RESULTADOS.md`, `OMEGA_R18_WL_ALTURA_RESULTADOS.md`, `OMEGA_CONSEJO_REV18_CIERRE_MEZCLA.md` |

## 3. Lo indeterminado

| Asunto | Estado | Por qué |
|---|---|---|
| R19-(b) | INDETERMINADO, no se reinterpreta | falló la cobertura con q = 1/8 y dos conos suaves quedaron ambiguos; R20 es otro experimento |
| Conos que el instrumento no ve | **no identificables por este instrumento** | la familia S1 de R20 (lente, tricilindro, cápsula, almohada) tuvo 0 de 4 visibles; la A2 de R19 quedó ambigua |
| Copias al 7% (G1a, G2a) | no identificables; G1a además no evaluable | su control con q = 1/8 no queda lo bastante fuera de banda |
| H_E (emergencia efectiva) | **sin evidencia confirmatoria** | R20 valida el instrumento; todavía no se aplicó a órdenes de origen desconocido |
| H_L (compatibilidad microscópica) | sin evidencia | hipótesis separada desde la Rev. 14 |
| H_sel (selección dinámica) | indeterminada; afirmación de mezcla cerrada con N de 64 a 128 | ninguna cadena mezcla en ese rango (H8); (ii), rediseñar el movimiento, queda en reserva |

## 4. Límites por escala, dimensión y mezcla

| Línea | Escala donde vale | Dimensión | Mezcla | Lo que no se puede afirmar |
|---|---|---|---|---|
| Ω-E (entropía SJ) | N ≥ 512 para el coeficiente | 2D | no con N ≥ 64 | que Ω-E seleccione geometría; que E1′ reconozca 4D |
| Ω-C (acción BD) | N = 32 mezcla | 2D (síntesis 2D) y 4D (acción) | no con N ≥ 64 | cualquier peso de Ω-C sobre fases geométricas |
| Reconocimiento (R14-C) | N = 4000 y 8000 | 4D | no aplica (órdenes dados) | separar conos suaves de curvatura |
| Instrumento Φ_q (R20) | N₀ = 8000, q hasta 1/8 (N ≈ 1000) | 4D, referencia diamante de Minkowski | no aplica | otras familias, otros N₀, otras referencias; la velocidad de decaimiento (§6) |

## 5. R20: correcciones y advertencias

**Correcciones registradas** (no cambian datos ni veredicto):
1. El riesgo de cobertura de la versión 1 (4.5%) usaba un modelo equivocado. Con el modelo correcto (Beta-Binomial(100, 96, 5)), el gate de la versión 2 (87/100) falla por azar con 0.0476 como mucho (Rev. 21).
2. El error de semillas de R19 (`{"k": 100 + j, **e}`) se corrigió en R20 con una clave aparte (`OMEGA_R19_GEMELOS_RESULTADOS.md` §5).
3. «Vector de 14 componentes» describe el inventario nominal. Con q = 1/8, el estadístico que votó usó 13: r₄₀ no es medible en el 43% de las referencias planas de ese tamaño (adenda A1).
4. La lectura 1 del informe de apertura de la sesión 11 («las copias diluidas se vuelven indistinguibles del plano») estaba exagerada. Quedan mayormente en banda, pero su distribución todavía se distingue de la plana (adenda A1 y `OMEGA_R22_AUDITORIA.md`).

**Advertencias:**
- **De implementación:** el guardián seleccionó con su propia banda y con D por par; el ejecutor votó con otra banda y con la máscara de las 300 referencias. La reproducción exacta de su selección muestra que, con la banda del ejecutor, S2b y S2d no habrían sido identificables. La causa es el ruido de la muestra plana de la banda. El veredicto no cambia, pero la afirmación robusta vale para 12 estructuras (`OMEGA_R22_AUDITORIA.md` §4).
- **De población:** la familia de conos suaves que votó (S2) reemplazó a una familia invisible (S1), y S2 incluye prismas de sección curva. La conclusión sobre conos suaves vale para S2, no para todos los conos suaves.

## 6. Interpretaciones vigentes (nivel 2, no son hechos)

- El control reinyectado funciona porque mantiene la fracción de defecto mientras el observado la pierde con q² (gemelos) o q³ (tripletes completos).
- Con q = 1/8, D_obs de las copias queda cerca de la distribución nula de la plana (mediana de D/radio cerca de 0.72 a 0.74). Ese piso es estadístico, no técnico: un punto casi plano tiene la distancia de un punto plano típico. Por eso L, con q = 1/8, depende sobre todo de cuán visible sigue el control, y no mide la velocidad de decaimiento.
- Con q = 1/2 los gemelos se separan más rápido que los tripletes a igual fracción. Es compatible con q² frente a q³, sin predicción previa.
- R14-C sugiere el máximo del espectro del Laplaciano como el observable más prometedor para conos suaves; sigue siendo hipótesis.

## 7. Hipótesis abiertas y qué las pondría a prueba

| Hipótesis | Prueba mínima que la tocaría | Requisito previo |
|---|---|---|
| H_E sobre órdenes no geométricos | opción (a): ¿el descendiente Φ_q(X) de un orden no geométrico se acerca a un blanco geométrico independiente más que un control fresco de la misma familia? | gate de la fase 3 cerrado; blanco, familias y margen preregistrados |
| Visibilidad de conos suaves | opción (b): un estadístico nuevo (por ejemplo, λ_max del Laplaciano) calibrado contra curvaturas | piloto aparte |
| Velocidad de decaimiento | opción (c): pendiente de D o de log D en una ventana fuera del piso nulo | piloto de rango dinámico |
| H_sel | un muestreador que mezcle con N ≥ 64 y una medida definida | condiciones de reingreso de la Rev. 18 |

## 8. Decisiones pendientes del usuario

R22-1 a R22-4 quedaron firmadas (`OMEGA_CONSEJO_REV22_R20.md` §6). Pendiente: la firma del preregistro del piloto (a) (`OMEGA_R23_PILOTO_A_PREREGISTRO.md`).
