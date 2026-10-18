# Consejo Ω, Revisión 14: dictamen sobre emergencia efectiva tras F6 (2026-10-18)

- **Origen:** dictamen de un Consejo externo de tres miembros convocado por el usuario fuera de esta sesión (según el texto, en otra herramienta con modelos GPT-6 Luna en esfuerzo alto). El usuario lo pegó en la conversación; se registra **tal cual** en §2. El expediente individual, las rondas de contraste y el acta quedaron fuera del repositorio, en carpetas de trabajo del usuario.
- **Estado:** RECOMENDACIÓN. **No se ejecuta nada** hasta que el usuario decida sobre R14-A a R14-C.

## 1. Verificación del orquestador contra el repositorio

| Afirmación del dictamen | Verificación |
|---|---|
| F6 intento 2: 13/13 positivas, 10/14 negativas, FALLA; D5 a D8 aceptadas | **Correcta.** Las 14 negativas con voto: 5 vistas (cúbico, octaédrico, p = 6, producto 4, percolación), C3 a C6 y D4 a D8. |
| τ_b ≈ 0.00267, τ_K ≈ 11.45, 15 grupos, N = 4000 y 8000 | **Correcta** (`results/omega_f6_intento2/summary.json`). |
| El cierre está en `bff0241`, rama `claude/hopeful-galileo-d7ffri`; la predeterminada `claude/hopeful-galileo-88u1l1` (716659a) no lo contiene; el PR #2 se fusionó antes del cierre | **Correcta.** La predeterminada contenía hasta `a410bc0`. Resuelto en esta sesión con un merge (eb757ae) que conserva los SHA citados por los datos; el cierre llegará a la predeterminada con un PR nuevo. |
| Citas: Bombelli et al. 1987, Surya 2019, Glaser y Surya 2013, Carlip 2024, Boguñá y Krioukov 2026 | **No verificadas por el orquestador.** Glaser y Surya (abundancias de intervalos como prueba de localidad) ya aparecía, sin verificar, en la opinión del asesor de curvatura de Rev. 11. Se registran como citas del dictamen, no como respaldo propio. |

**Coincidencias con el programa:** la separación entre reconocimiento (H_rec), emergencia (H_E) y selección (H_sel) es la de Rev. 11 §2. La recomendación de no abrir un tercer F6 respeta R13-4.

## 2. Texto del dictamen (registrado tal cual)

## Decisión

**Como programa de trabajo, Ω debe priorizar investigar si una estructura causal discreta produce una geometría lorentziana efectiva.** Debe conservar la hipótesis de geometría microscópica exacta como una prueba independiente. El consejo no afirma que Ω ya genere geometría efectiva ni que haya resuelto qué es fundamental en la naturaleza.

La distinción importa porque un orden causal discreto puede ser la estructura microscópica exacta sin que exista una métrica suave asignada a cada evento. En causal set theory, el objeto fundamental es un orden parcial localmente finito; el orden y el conteo de elementos pueden aproximar la causalidad y el volumen de un espacio-tiempo continuo. La geometría del continuo se formula como una aproximación dependiente de escala. La literatura mantiene abierta la cuestión de qué conjuntos admiten esa aproximación y de su unicidad. [Bombelli et al. (1987)](https://doi.org/10.1103/PhysRevLett.59.521), [Surya (2019)](https://doi.org/10.1007/s41114-019-0023-1)

Por tanto, la pregunta elegida es **qué propiedades geométricas emergen, a qué escalas, con qué precisión y mediante qué dinámica**, no si cada microestado ya es una variedad lorentziana continua. La compatibilidad microscópica H_L sigue siendo valiosa, pero debe definirse con un nulo verificable —por ejemplo, compatibilidad con un proceso de sprinkling fiel a densidad fijada— en vez de usar “exacta” sin una definición discreta.

## Qué resolvió el consejo

Tres miembros independientes, configurados con GPT-6 Luna y esfuerzo alto, estudiaron el problema desde geometría causal, emergencia/dinámica e inferencia/estadística. Después contrastaron entre sí sus dictámenes. Los tres coinciden en que la evidencia actual no decide la ontología, que H_rec, H_E y H_sel deben separarse y que se necesita un protocolo preregistrado antes de otra confirmación. Dos priorizan H_E explícitamente como agenda; el tercero acepta esa prioridad solo como decisión metodológica, no ontológica.

## F6, intento 2: evidencia y alcance

El resultado oficial está en el commit [`bff0241`](https://github.com/SebwMane/FichaClaraV2/blob/bff0241/docs/OMEGA_F6_INTENTO2_RESULTADOS.md), en la rama `claude/hopeful-galileo-d7ffri`. El `summary.json` enlazado al mismo commit registra `FALLA`, 13/13 positivos bien clasificados y 10/14 negativos bien clasificados. Cuatro negativos —D5 (gemelos), D6 (tapa plana), D7 (huso) y D8 (suma de normas)— fueron aceptados como espacio-tiempo; sus razones normalizadas para b y K quedaron bajo el umbral. Resultado agregado de este catálogo: **23/27**, o 85,2 %, pero no es una estimación de precisión poblacional porque el catálogo fue diseñado y no representa una muestra aleatoria de posibles estructuras. [Resumen de corrida](https://github.com/SebwMane/FichaClaraV2/blob/bff0241/results/omega_f6_intento2/summary.json), [prerregistro](https://github.com/SebwMane/FichaClaraV2/blob/bff0241/docs/OMEGA_F6_INTENTO2_PREREGISTRO.md)

F6 usó 15 grupos de calibración, τ_b ≈ 0,00267 y τ_K ≈ 11,45, con tamaños N=4000 y N=8000. El test de colapso sí separó C5 (links borrados) y D4 (relaciones agregadas), pero no D5–D8. El resultado demuestra que **b y K no bastan para la discriminación pretendida en esta batería**. No demuestra que toda estructura causal con esas estadísticas sea lorentziana, ni que la geometría efectiva no pueda aparecer en otra escala.

D5 permanece negativa según la etiqueta prerregistrada. Su clasificación no debe reescribirse después del resultado. La cuestión distinta —si un microestado con twins puede conservar una geometría efectiva al coarse-grain— queda abierta y debe probarse con un diseño nuevo.

Los conteos son descriptivos de este conjunto de casos, no tasas universales. Los 15 grupos de calibración y las réplicas no equivalen a una población amplia de familias de estructuras. Las celdas de intervalos se solapan; la próxima inferencia debe estimar incertidumbre por estructura/familia generadora, no tratar subintervalos correlacionados como réplicas independientes.

## Por qué la agenda efectiva es la más útil ahora

La teoría de conjuntos causales estudia explícitamente la transición entre una estructura parcial y localmente finita y una descripción continua. El resultado de unicidad de la aproximación continua —la Hauptvermutung— sigue siendo una conjetura, no una garantía para cualquier poset finito. [Surya (2019)](https://doi.org/10.1007/s41114-019-0023-1)

Los perfiles de abundancia de intervalos se han propuesto como pruebas de localidad en causal sets manifoldlike. Son más informativos que una única fracción de pares, pero no son certificados universales de geometría. [Glaser y Surya (2013)](https://doi.org/10.1103/PhysRevD.88.124026)

La literatura reciente ilustra tanto el progreso como el límite: análisis de 2024 muestran supresión dinámica de amplias clases de causal sets estratificados no manifoldlike, pero dejan sin resolver qué estructuras restantes sobreviven y no prueban una emergencia general del continuo. [Carlip (2024)](https://doi.org/10.1007/s10714-024-03281-1). En 2026, un operador d’Alembertiano local definido en causal sets demostró convergencia a su contraparte continua para campos en espacio-tiempo de Minkowski; es una recuperación concreta de un observable, no una derivación completa de la geometría o de GR en Ω. [Boguñá y Krioukov (2026)](https://arxiv.org/abs/2506.18745)

## Resolución R14

El consejo recomienda R14-A–C como auditoría analítica y fase exploratoria, sin voto confirmatorio:

1. **R14-A — inferencia y trazabilidad.** Mantener intactos los resultados de F6; revisar calibración, covarianza entre escalas, dependencia entre intervalos y definición de la unidad independiente. La calibración por máximo de 15 grupos no fija por sí sola una tasa de error para otras familias.
2. **R14-B — H_L frente a H_E.** Especificar dos modelos que hagan predicciones distintas. Para H_L, fijar geometría de referencia, densidad y compatibilidad del orden microscópico con su distribución nula. Para H_E, fijar mapa de coarse-graining y escalas antes de observar datos, y probar convergencia/estabilidad de orden, volumen, dimensión, perfil de intervalos y observables de curvatura.
3. **R14-C — observables causales multivariables.** Explorar distribuciones completas de tamaños de intervalos, abundancias de cadenas y estimadores independientes de dimensión/curvatura. Usar familias generadoras independientes, bandas de equivalencia preregistradas y controles sellados nuevos. D5–D8 son ejemplos ya vistos; no deben servir como test confirmatorio nuevo.
4. **H_sel — dinámica.** Una etapa posterior debe ejecutar la dinámica o medida de Ω y estimar cuánto peso asigna a fases que satisfacen H_E. La recuperación geométrica de estructuras introducidas a mano no responde esta pregunta.

**No se recomienda un tercer F6.** El intento 2 quedó registrado como último en R13-4. Cualquier prueba confirmatoria nueva debe nacer de R14, con protocolo nuevo y controles ciegos no usados para diseñar el reconocedor.

## Trazabilidad del repositorio

El repositorio se clonó inicialmente desde la rama por defecto `claude/hopeful-galileo-88u1l1`, cuyo `HEAD` remoto consultado fue `716659a`; no contiene el cierre F6 del 18 de octubre. El informe, código y resultados del intento 2 están en `bff0241`, rama `claude/hopeful-galileo-d7ffri`. El PR #2 aparece fusionado el 16 de octubre, antes del commit del cierre. Conviene integrar el cierre con un PR separado conservando el SHA de la corrida, el sello y sus datos; no se cambió el repositorio para este dictamen. [PR #2](https://github.com/SebwMane/FichaClaraV2/pull/2)

## Límites de esta revisión

El consejo analizó los documentos, el prerregistro, el resumen JSON, el código y resultados versionados en el commit indicado. No volvió a ejecutar las 158 partes ni realizó una reproducción numérica independiente. Por tanto, los conteos se reportan como verificados en los artefactos publicados, no como una réplica nueva. El consejo tampoco demuestra que una dinámica de Ω produzca causal sets manifoldlike. Esa es la hipótesis H_sel que queda para una fase posterior.

## 3. Decisiones pendientes del usuario

| # | Pregunta |
|---|---|
| R14-1 | ¿Adoptar la emergencia efectiva (H_E) como prioridad del programa, con H_L como hipótesis separada? |
| R14-2 | ¿Ejecutar R14-A (auditoría de inferencia de F6 sobre datos ya publicados, sin corridas nuevas)? |
| R14-3 | ¿Preparar R14-B y R14-C como preregistros exploratorios (sin voto)? |
| R14-4 | ¿Abrir un PR nuevo desde `claude/hopeful-galileo-d7ffri` hacia la rama predeterminada para integrar el cierre de F6 y este dictamen? |

## 4. Firmas (2026-10-18)

| # | Decisión |
|---|---|
| R14-1 | **Firmada:** la emergencia efectiva (H_E) es la prioridad situacional del programa; la compatibilidad microscópica (H_L) queda como hipótesis separada. |
| R14-2 | **Firmada y ejecutada:** auditoría R14-A en `OMEGA_R14A_AUDITORIA.md` (sin corridas nuevas; figura `docs/figuras/r14a_auditoria.png`). |
| R14-3 | **Firmada:** preregistros exploratorios sin voto, preparados y no ejecutados: `OMEGA_R14B_PREREGISTRO.md` y `OMEGA_R14C_PREREGISTRO.md`. |
| R14-4 | **Firmada:** PR nuevo desde `claude/hopeful-galileo-d7ffri` hacia la rama predeterminada. |
