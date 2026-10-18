# Consejo Ω, Revisión 15: acta de la sesión 2 (2026-10-18)

- **Origen:** segundo Consejo externo convocado por el usuario fuera de esta sesión (según el acta: tres dictámenes independientes con GPT-6 Luna en esfuerzo alto, crítica cruzada y resolución de su orquestador). El usuario pegó el acta en la conversación; se registra **tal cual** en §2. Los dictámenes individuales están fuera del repositorio (`work/consejo_omega_sesion_2/` del usuario).
- **Estado:** RECOMENDACIÓN. Las decisiones de §3 son del usuario.

## 1. Verificación del orquestador

| Afirmación del acta | Verificación |
|---|---|
| Punta de `claude/hopeful-galileo-d7ffri` = 853ac92; predeterminada = 716659a, sin el cierre de F6 ni la Rev. 14 | **Correcta** al momento del registro. El cierre llegará a la predeterminada si se fusiona el PR [SebwMane/FichaClaraV2#3](https://github.com/SebwMane/FichaClaraV2/pull/3). |
| F6 intento 2: FALLA, 13/13 y 10/14; R14-A con D6 aceptada en el 99.7% de los remuestreos | **Correcta** (`OMEGA_F6_INTENTO2_RESULTADOS.md`, `OMEGA_R14A_AUDITORIA.md`). |
| R14-B y R14-C: exploratorios, sin voto, no ejecutados; ninguno muestrea una dinámica de Ω | **Correcta.** R14-B diluye estructuras dadas; no prueba H_sel. |
| Citas nuevas: Müller (arXiv:2503.01719, versión de septiembre de 2026), Eichhorn, Mack, Le y Wagner (arXiv:2605.27514), Boguñá y Krioukov (PRD, doi 10.1103/np89-qzjp) | **No verificadas por el orquestador.** Se registran como citas del acta. Antes de usar los observables de Eichhorn et al. en R14-C hay que leer el preprint. |

## 2. Acta (registrada tal cual)

## Resolución

**El Consejo reafirma H_E —investigar la emergencia de una geometría lorentziana efectiva— como prioridad metodológica del programa Ω. H_L, la compatibilidad microscópica, debe mantenerse como hipótesis separada y operacionalizable.** Esta decisión prioriza una línea de investigación; no demuestra que Ω ya genere geometría lorentziana, ni determina cuál es la ontología de la naturaleza.

La microestructura de Ω puede ser exacta en sus propios términos —por ejemplo, como un orden discreto— sin que cada elemento tenga asignada una métrica suave. Para hacer contrastable H_L hay que decir qué propiedad discreta exige, frente a qué geometría de referencia y con qué tolerancia. Deben conservarse separadas tres preguntas: H_rec (qué reconoce un instrumento en estructuras dadas), H_E (qué geometría se obtiene al coarse-grain una estructura) y H_sel (si la dinámica y medida de Ω generan o favorecen esa fase).

**Decisión prescriptiva sobre el trabajo inmediato:** no ejecutar todavía R14-B ni R14-C en su versión actual. Enmendar primero sus especificaciones y preservar los prerregistros actuales como versiones históricas. Si se ejecutan luego, sus resultados seguirán siendo exploratorios y no confirmarán por sí solos que Ω produzca geometría. No reabrir F6 ni iniciar un tercer intento.

## Análisis prescriptivo (extracto: enmiendas antes de decidir sobre R14-B/C)

1. **Acotar la afirmación de R14-B.** La thinning de estructuras ya suministradas estudia cómo responde un orden a un mapa de coarse-graining; no demuestra que la dinámica de Ω produzca esos órdenes. Etiquetar el resultado como exploración de compatibilidad/coarse-graining, salvo que se especifique y use un generador de Ω.
2. **Definir referencia y equivalencia.** Justificar el margen δ=0,5 desviaciones estándar, especificar cobertura conjunta y dependencia entre observables, y usar referencias de sprinkling apropiadas para cada geometría/curvatura objetivo. La referencia plana por sí sola no sustenta compatibilidad frente a toda geometría curva. Precisar qué significa converger en la ventana discreta q y qué se mantiene fijo entre escalas.
3. **Respetar la jerarquía de datos.** Tratar q de una misma realización como medidas emparejadas; cuantificar incertidumbre por realización y familia generadora. Subintervalos correlacionados no cuentan como réplicas independientes. Separar estructuras ya vistas de familias y controles nuevos destinados a confirmación.
4. **Cerrar R14-C antes de correrlo.** Verificar las fórmulas que el propio prerregistro deja pendientes, fijar normalización y resumen para cada observable y definir cómo se informará la multiplicidad. Añadir como candidatos exploratorios la distribución de grado de links y el espectro del Laplaciano del Hasse simetrizado, motivados por el preprint de 2026; especificarlos y congelarlos antes de inspeccionar salidas. La abundancia de intervalos ya se solapa con el perfil de tamaños propuesto. Conservar resultados completos y no seleccionar después una combinación ganadora para clasificar estas mismas estructuras.
5. **Preparar una prueba auténtica de H_E/H_sel.** Para sostener que Ω produce geometría, se necesita un ensamble o dinámica de Ω, un mapa de coarse-graining justificado por el modelo, una ventana de escalas y observables prospectivos, y controles sellados de familias nuevas. H_sel requiere medir el peso/probabilidad de la fase bajo la dinámica, no solo reconocer o adelgazar estructuras diseñadas.

Tras completar esas enmiendas, corresponde una nueva decisión explícita sobre ejecutar B/C como exploraciones sin voto. Una confirmación posterior requiere familias no vistas y un protocolo congelado antes de abrir sus resultados. Un resultado negativo debe limitarse al mapa, escala y ensamble probados.

## Literatura citada por el acta (resumen del acta, sin verificar aquí)

- Müller (arXiv:2503.01719, versión del 3 de septiembre de 2026): dos versiones finitas de la Hauptvermutung, una falsa y otra verdadera, más una versión numerable; el autor advierte que podrían no coincidir con la conjetura física. «La Hauptvermutung sigue abierta» debe precisar la versión.
- Eichhorn, Mack, Le y Wagner (2026, arXiv:2605.27514): en nueve clases de causal sets, la distribución de grado de links, el espectro del Laplaciano del diagrama de Hasse simetrizado y la abundancia de intervalos distinguen las clases.
- Boguñá y Krioukov (2026): convergencia de un d'Alembertiano local intrínseco en Minkowski.
- Carlip (2024): supresión de ciertas clases no manifold-like en sumas de caminos causales, sin demostrar la emergencia del continuo.

## Límites (del acta)

La sesión cotejó documentos, resúmenes, prerregistros y código versionados; no reejecutó F6 ni R14-A, y no corrió R14-B/C. Las fuentes de conjuntos causales son antecedentes especializados; no sustituyen una formulación ni dinámica completa de Ω.

## 3. Decisiones pendientes del usuario

| # | Pregunta |
|---|---|
| R15-1 | ¿Redactar R14-B v2 y R14-C v2 con las enmiendas 1 a 4, conservando las versiones actuales como históricas, y volver a decidir después si se ejecutan? |
| R15-2 | ¿Verificar antes las citas nuevas (Müller; Eichhorn et al.) y las fórmulas pendientes de R14-C con un agente aparte? |
| R15-3 | ¿Abrir el diseño de la prueba de H_E y H_sel con una dinámica de Ω (enmienda 5)? Requiere elegir qué ensamble o dinámica de Ω se usa (por ejemplo, la acción de Benincasa-Dowker de Ω-C o la ruta entrópica de Ω-E). |
