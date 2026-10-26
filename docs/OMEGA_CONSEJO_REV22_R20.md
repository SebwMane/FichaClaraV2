# Consejo Ω, Revisión 22: dictamen de la sesión 11 sobre R20, auditoría y pasos adoptados (2026-10-26)

- **Fuente:** dictamen «Secuencia lógica recomendada», pegado por el usuario. Responde al informe de apertura `docs/consejo_rev22/INFORME_APERTURA.md`.
- **Mandato del usuario:** «con esta nueva sesión del consejo, auditala, adopta los pasos que consideres prudentes, puedes refinarlo si lo estimas conveniente».
- **Criterio del orquestador para «prudente»:** adopto lo que no cambia la hipótesis, no corre ningún experimento confirmatorio y no reinterpreta un resultado histórico (aclaraciones, auditorías sobre datos guardados o con semillas ya publicadas, y documentos). Lo que fija alcance o abre un estudio nuevo queda para la firma del usuario (§5).

## 1. Dictamen (resumen fiel)

**Dirección:** se conserva la recomendación de la sesión 11, (d) síntesis y después (a) exploración, dividida en fases con dependencias. La síntesis no espera datos nuevos. Ningún uso confirmatorio nuevo de Φ_q antes de resolver la auditoría de especificación y máscara.

| Fase | Contenido |
|---|---|
| 1 | **R22-1:** conservar el PASA computacional para las 14 votantes del universo final, «bajo la métrica y el procedimiento ejecutados», no como evidencia de emergencia física. **R22-2:** registrar S1 → S2 como desviación material y limitar la afirmación sobre conos a las 4 S2; si el objetivo era el universo original, esa rama es INDETERMINADA. **Correcciones del informe:** 14 observables nominales frente a 13 usados con q = 1/8; componentes por q; comparar máscaras del guardián y del ejecutor; separar el gate empírico de la garantía Beta-Binomial. |
| 2 | Síntesis del programa (opción d), en paralelo con la auditoría, separando hechos, interpretaciones, hipótesis y decisiones, con una escala de afirmación de cuatro escalones. |
| 3 | **Gate antes de reutilizar Φ_q:** (1) reconciliar la dimensión; (2) revisar faltantes y comparabilidad entre q; (3) resolver la garantía bajo el selector real; (4) reproducir la selección del guardián con la métrica final; (5) diagnosticar el borde cercano a 0.8. Salida: (i) garantía formal, (ii) nueva calibración con instrumento congelado o (iii) herramienta solo exploratoria. |
| 4 | Preparar (a) como piloto exploratorio: distancia absoluta a un blanco geométrico independiente, padre no geométrico X_N, descendiente Y_q = Φ_q(X_N), control fresco Z_q del mismo tamaño; unidad independiente = el padre; márgenes de R20 no trasladables. |
| 5 | Ramas secundarias: (c) piloto de rango dinámico; (b) si la prioridad es ampliar identificabilidad; H_sel después de resolver la mezcla. |

**Decisiones sugeridas al usuario:** confirmar el PASA en ese alcance; aceptar S1 → S2 como desviación restringida a S2; aprobar (d) y exigir el gate antes de invocar Φ_q como validado; autorizar preparar (a), con blanco 3+1 solo si esa es la dimensión objetivo. Explícitamente **no son firmas ni autorizan una corrida.**

## 2. Auditoría del dictamen

Verifiqué cada afirmación técnica contra los datos de R20 (`tools/r22_auditoria.py`, `tools/r22_seleccion.py`; resultados en `docs/OMEGA_R22_AUDITORIA.md`).

| Punto del dictamen | Verificación | Juicio |
|---|---|---|
| 14 nominales y 13 usados con q = 1/8, sin r₄₀ | **Correcto.** Las máscaras cambian solo con q = 1/8. r₄₀ falta en el 43% de las referencias planas de ese tamaño (límite de localidad), así que no es un caso de borde. | adoptado (adenda A1) |
| Máscaras del guardián y del ejecutor | Las dos excluyen r₄₀ con q = 1/8; difieren en la banda y en el D por par. Reproduje la selección con las dos métricas sobre las mismas realizaciones (§4 de la auditoría). | adoptado |
| La garantía bajo el selector real | La máscara es una función simétrica de las 300 referencias y coincide con la que da el conjunto de ajuste solo, así que calibración y validación siguen siendo intercambiables dado la máscara. El gate vale. | adoptado, con la conclusión de la auditoría §3 |
| Borde cercano a 0.8 | Es la distribución nula de D/radio de un punto plano, no una saturación. Además mostró que **mi lectura 1 de la sesión 11 estaba exagerada** (las copias diluidas todavía se distinguen de la plana). | adoptado, con corrección propia |
| Separar gate empírico y riesgo Beta-Binomial | Correcto; la Rev. 21 ya lo decía, pero el informe de apertura no lo separaba. | adoptado (adenda A1) |
| Unidad independiente = el padre en (a) | Correcto. En R20 la unidad fue la realización y el bootstrap ya remuestreaba los 4 q juntos; (a) debe hacer lo mismo por padre. | adoptado para el preregistro de (a) |
| Márgenes de R20 no trasladables | Correcto: ±ln 1.3 se fijó con el ruido del borde de la banda de R19 y quedó holgado en R20. | adoptado |

**No encontré errores en el dictamen.** Tengo cinco refinamientos (§3).

## 3. Refinamientos del orquestador

1. **El blanco 3+1 ya existe en el instrumento.** Las referencias de R20 son sprinklings nuevos de un diamante causal de Minkowski 4D (`r14b.GEN["plana"]`, de `diamante_4d`), y las bases de copias también son 4D. Elegir 3+1 como blanco de (a) no agrega un supuesto nuevo al instrumento: hace explícito uno que ya estaba. Igual es una decisión del usuario (R22-4), porque fija la dimensión física objetivo.
2. **«Entrar en banda» es un criterio débil de proximidad.** En R20, las copias con q = 1/8 quedaron entre el 63% y el 87% de las veces en banda, y aun así su distribución se distingue de la plana. Para (a) propongo definir la proximidad como **distribución de D compatible con la nula plana del mismo N** (prueba de dos muestras con un margen preregistrado), no como pertenencia a la banda.
3. **El piso nulo limita lo que se puede medir.** Una distancia de Mahalanobis a una referencia no baja de la distribución de un punto plano típico. Para (a) y (c) conviene usar un estadístico de **exceso sobre la nula** (por ejemplo, el rango de D del descendiente entre las D planas del mismo N, o D² menos su media nula), que no tiene techo.
4. **Comparar D entre q exige un vector común.** La auditoría muestra que pasar a un vector común de 13 componentes en los 4 q no cambia nada material en R20, así que (a) puede fijarlo desde el principio sin costo.
5. **Selección por visibilidad en (a).** Para órdenes no geométricos la visibilidad con q = 1 es casi segura (R14-C: percolación separa en 49 de 49 observables). El riesgo es el opuesto: familias que **nunca** se acercan, porque nada en Φ_q las transforma en 4D (por ejemplo, KR diluido sigue en 3 capas). El preregistro de (a) debe declarar la predicción por familia antes de correr, para que «no se acerca» sea un resultado y no una sorpresa.

## 4. Pasos adoptados (sin firma, por el criterio del encabezado)

| Paso | Qué se hizo | Dónde |
|---|---|---|
| Fase 1, correcciones | Adenda A1 al informe de R20 (componentes por q, comparación entre q, máscara del guardián, gate frente a garantía, formulación del PASA) y nota de corrección en el informe de apertura | `OMEGA_R20_RESULTADOS.md` §5; `consejo_rev22/INFORME_APERTURA.md` |
| Fase 2, síntesis | Borrador con escala de cuatro escalones, hechos, indeterminados, límites, correcciones, interpretaciones e hipótesis | `OMEGA_SINTESIS_R22.md` |
| Fase 3, gate | Los cinco puntos, sobre datos guardados de R20 y con la reproducción de la selección del guardián con sus semillas publicadas (entropía 20261114). No se usó ninguna semilla confirmatoria. | `OMEGA_R22_AUDITORIA.md`; `results/omega_r22/` |

**No adoptado:** el preregistro de (a). Depende de R22-3 y R22-4, y de los refinamientos 1 a 5.

## 5. Decisiones para el usuario

| # | Decisión | Recomendación del orquestador |
|---|---|---|
| R22-1 | Confirmar el PASA de R20 como PASA computacional para las 14 votantes, bajo la métrica y el procedimiento ejecutados | **sí** |
| R22-2 | Registrar S1 → S2 como desviación material; la afirmación sobre conos suaves vale solo para S2 | **sí** |
| R22-3 | Aprobar la síntesis (d) y la salida del gate de la fase 3 | **sí**, con la salida que recomienda la auditoría (§6 de `OMEGA_R22_AUDITORIA.md`) |
| R22-4 | Autorizar el preregistro de (a) como piloto exploratorio, con blanco geométrico en un diamante de Minkowski 3+1 independiente | **sí**, si 3+1 es la dimensión física objetivo del proyecto |
