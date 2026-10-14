# Contexto para retomar el proyecto Ω en un chat nuevo (Claude Code)

Fecha de corte: 2026-10-12. Autor del proyecto: Sebastían Ignacio Ramos Bobadilla (Sebw; el nombre se escribe «Sebastían», con tilde en la í).

## 1. Qué estamos haciendo

Proyecto personal y experimental: la **Hipótesis del Espacio Independiente (modelo Ω)**. Idea central de Sebastían (Ω-0.2): el espacio-tiempo no es fundamental; emerge de un sustrato de relaciones que puede tener fases distintas. El trabajo consiste en convertir esa idea en hipótesis verificables y probarlas con simulaciones, con una metodología de preregistro y auditoría.

Recorrido hasta hoy, en una línea por etapa:

1. **Ω-1.1 y bloque L** (antes de este chat): el funcional S0 quedó cerrado por incompatible.
2. **Ω-C** (orden causal más acción discreta de Benincasa-Dowker): probado en 2D y 4D con peso euclídeo y lorentziano. **Ninguna rama seleccionó geometría.** Hallazgos: el mínimo de la acción es un «sándwich»; fluctuaciones de la acción correlacionadas a gran distancia; extremos no geométricos en ambos sentidos en 4D.
3. **Ω-E** (caso 3 de Ω-0.2, ruta entrópica de Jacobson): la entropía con ley de área **no** separa la geometría del orden KR con un criterio simétrico. Apareció en cambio una **firma espectral** (pendiente y rodilla del espectro del operador de Pauli-Jordan).
4. **Firma espectral:** en 2D pasó una prueba ciega (F1). En 4D, con el propagador correcto (links), reconoce «orden continuo de dimensión d», pero no la forma de los conos: el orden producto la imita (F3). El **perfil de intervalos** (brazo C de F5) sí rechazó a todos los impostores, y es el candidato para reconocer espacio-tiempo lorentziano en 4D.

## 2. Dónde está todo

- **Repositorio:** GitHub `SebwMane/FichaClaraV2`, rama `claude/hopeful-galileo-88u1l1`.
- **Base:** commit `ec448e4` (fin del bloque L). Todo el trabajo de este chat son los commits posteriores, entregados en `omega_c0.bundle` y en la carpeta `patches/`.
- **Aplicar el paquete:** en la rama, `git pull --ff-only /ruta/omega_c0.bundle HEAD` (o `git am patches/*.patch`) y luego `git push origin claude/hopeful-galileo-88u1l1`. Si el push ya se hizo, verificar con `git log` que el último commit sea el de este documento.
- **Verificación hecha aquí:** suite completa 711 tests sin fallos y `mypy --strict` limpio en 104 archivos (antes de las pruebas F1 a F5; los módulos posteriores tienen sus tests en `tests/test_causality.py`).

## 3. Estado actual y lo pendiente

**Actualización 2026-10-14:** F5′ ejecutada (R9-1 y R9-2 firmadas) y **FALLA** (`docs/OMEGA_F5PRIMA_RESULTADOS.md`): las cuatro positivas quedan bien, incluido de Sitter 3D y 4D a ciegas, pero el cono octaédrico 4D pasa como espacio-tiempo. Rama de trabajo actual: `claude/hopeful-galileo-d7ffri`. Próximo libre: EXPERIMENT_ID 319, entropía 20261072.

Otras líneas abiertas, sin preparar: **F6** (comprometida por R9-2) (perfil medido solo en el interior, para regiones con bordes y topología distintos); el riesgo del imitador de 2-órdenes con más de 400 capas en 2D; la etapa E2 (justificar el truncamiento desde el sustrato); y la pregunta que Sebastían no ha respondido: qué parte de Ω-0.2 quiere conservar sobre todo.

## 4. Metodología (no negociable)

- **Roles separados:** el orquestador preregistra y escribe el handoff; el agente de código implementa solo eso; el orquestador audita como código ajeno y documenta. En el chat no había subagentes, así que la misma instancia hizo los roles y lo declaró. En Claude Code se pueden usar subagentes reales.
- **Preregistro:** criterios, semillas y umbrales congelados en un commit antes de correr. Cambios posteriores solo como **enmiendas fechadas**, indicando si pueden crear un resultado positivo.
- **Pilotos declarados:** semillas distintas de las oficiales; no cuentan como evidencia.
- **Corrida oficial** sobre árbol git limpio; los resultados negativos se registran con la misma visibilidad que los positivos; los intentos fallidos se conservan, no se borran.
- **Escala de afirmación:** 1 demostrado o medido; 2 explicación propuesta; 3 resultado externo; 4 extrapolación.
- **El Consejo Ω** (físico teórico, estadístico, metodólogo, escéptico y presidencia) delibera cuando el usuario delega una decisión; sus documentos son `docs/OMEGA_CONSEJO_REV*.md`.
- **Las decisiones que cambian la hipótesis las firma el usuario** (R3-x, R4-x, …). Nada de dinámica ni de prueba nueva sin su aprobación.

## 5. Mapa de documentos (en `docs/`)

**Puntos de entrada:** `handoff/HANDOFF_SIGUIENTE.md` (bitácora de estado, la más reciente arriba), `OMEGA_2_INFORME_FINAL.md` (informe del programa Ω-C con anexo de Ω-E; versión editable en el Claude Doc «Informe final del programa Ω-C»), `OMEGA_2_SINTESIS_2D.md` (síntesis de Ω-C con anexos 4D).

**Ω-C (hipótesis, enmiendas y resultados):** `OMEGA_2_PROPUESTA_CAUSAL.md` (§1 a §21, todo el preregistro de Ω-C), y los resultados `OMEGA_2_C0_RESULTADOS.md`, `OMEGA_2_L2C_RESULTADOS.md`, `OMEGA_2_C1_L3CE_RESULTADOS.md`, `OMEGA_2_L3CE_PRIME_RESULTADOS.md`, `OMEGA_2_L3CE_HYST_RESULTADOS.md`, `OMEGA_2_L3CE_EQ_RESULTADOS.md`, `OMEGA_2_BLOQUE_A_RESULTADOS.md`, `OMEGA_2_BLOQUE_AS_RESULTADOS.md`, `OMEGA_2_L_EPS_RESULTADOS.md`, `OMEGA_2_C4_L4_RESULTADOS.md`, `OMEGA_2_L4B_RESULTADOS.md`, `OMEGA_2_P4_RESULTADOS.md`.

**Ω-E y firma espectral:** `OMEGA_E_PROPUESTA.md` y `OMEGA_E_RESULTADOS.md` (E0, E1); `OMEGA_CONSEJO_REV4_E1PRIMA.md` (diseño de E1′, pilotos, rodilla); `OMEGA_E1PRIMA_PREREGISTRO.md` y `OMEGA_E1PRIMA_RESULTADOS.md` (E1′, rebajado a «dependiente del criterio»); `OMEGA_CONSEJO_REV5_REPLICA.md`, `OMEGA_E1DOBLEPRIMA_PREREGISTRO.md` y `OMEGA_E1DOBLEPRIMA_RESULTADOS.md` (réplica simétrica); `OMEGA_CONSEJO_REV6_FIRMA.md` (qué es la firma); `OMEGA_F1F2_PREREGISTRO.md` y `OMEGA_F1F2_RESULTADOS.md`; `OMEGA_CONSEJO_REV7_DIMENSION.md` (matriz causal frente a links); `OMEGA_F3F4_PREREGISTRO.md` y `OMEGA_F3F4_RESULTADOS.md`; `OMEGA_CONSEJO_REV8_F5.md` y `OMEGA_F5_RESULTADOS.md`; `OMEGA_CONSEJO_REV9_F5PRIMA.md` (pendiente).

**Auditorías y handoffs:** `handoff/AUDITORIA_C0.md` (todas las auditorías, §1 a §20) y los `handoff/HANDOFF_*_AGENTE_CODIGO.md`.

## 6. Resumen de veredictos

| Prueba | Veredicto | Documento |
|---|---|---|
| C0, C0-4b, C1, C4, C4b (instrumentos) | Éxito (C0-4 retirado) | `OMEGA_2_C0_RESULTADOS.md`, `OMEGA_2_C1_L3CE_RESULTADOS.md`, `OMEGA_2_C4_L4_RESULTADOS.md`, `OMEGA_2_L4B_RESULTADOS.md` |
| Rama euclídea 2D | Estable localmente con β = b/N; equilibrio sin resolver | `OMEGA_2_L3CE_*` |
| Rama lorentziana 2D y 4D | Cerrada en 2D; incompatible en 4D | `OMEGA_2_L_EPS_RESULTADOS.md`, `OMEGA_2_C4_L4_RESULTADOS.md`, `OMEGA_2_L4B_RESULTADOS.md` |
| P4 (euclídea 4D) | Sin ventana | `OMEGA_2_P4_RESULTADOS.md` |
| E0 / E1 / E1′ / E1″ | Reproducido / indeterminado / dependiente del criterio / discrimina solo por huella | `OMEGA_E_*`, `OMEGA_E1PRIMA_*`, `OMEGA_E1DOBLEPRIMA_*` |
| F1 (firma 2D ciega) | **Pasa** (márgenes estrechos en capas 400 y percolación 0.02) | `OMEGA_F1F2_RESULTADOS.md` |
| F2 (β → −2/d) | Falla | ídem |
| Rev. 7 (H1 a H3) | Con matriz causal la pendiente es frágil fuera de 2D; con links es estable pero no −2/d | `OMEGA_CONSEJO_REV7_DIMENSION.md` |
| F3 (firma 4D con links) | Falla: el orden producto 4D la imita | `OMEGA_F3F4_RESULTADOS.md` |
| F4 (β_L frente a d) | Estable y monótona: −0.05, −0.22, −0.31, −0.34 de 2D a 5D | ídem |
| F5′ (perfil con d* espectral) | Falla: positivas bien (incluido de Sitter a ciegas), pero el cono octaédrico pasa | `OMEGA_F5PRIMA_RESULTADOS.md` |
| F5 (brazos A, C, B) | Fallan; C rechazó a todos los impostores y falló por diseño en el bloque y por umbral justo en de Sitter | `OMEGA_F5_RESULTADOS.md` |

## 7. Código

**Módulos (`omega/causality/`):** `sprinkling.py` (diamantes de Minkowski, matriz causal), `orders.py` (2-órdenes, KR, percolación, capas, bipartito), `estimators.py` (Myrheim-Meyer, punto medio, abundancias de intervalos, acción BD 2D y 4D, suavizadas), `exact.py` (esperanzas exactas en 2-órdenes), `profile.py` (discriminador C0-4b), `phase.py` (promedios de fase), `mcmc.py` (Metropolis en 2-órdenes), `mcmc_general.py` (Metropolis en órdenes generales, ΔS₄ incremental), `wanglandau.py`, `torus.py` (bloque T³ × [0, T]), `reference4d.py` (referencia por pares), `entropy.py` (vacío de Sorkin-Johnston, SSEE, regiones, rodilla).

**Herramientas (`tools/`):** una por prueba (`c0_validation.py` … `f5.py`), ejecutadas por partes (`--part`) para respetar el límite de tiempo de cada comando, y con `final` para calcular el veredicto. Resultados en `results/<prueba>/` (están versionados).

**Convenciones:** RNG solo con `Generator(PCG64(SeedSequence(entropía, spawn_key=(EXPERIMENT_ID, …))))`. EXPERIMENT_ID usados: 300 a 317 (el próximo libre es **319**). Entropías usadas: 20261011 a 20261068 (la próxima libre es **20261072**). Nunca reutilizar una entropía oficial.

## 8. Trampas técnicas aprendidas

- **Encadenar comandos sin detenerse ante errores** dejó pasar commits con tests o tipos rotos. Usar `bash -c` con `set -euo pipefail`; `| tail` oculta el código de salida sin `pipefail`.
- **El propagador depende de la dimensión:** matriz causal C en 2D; matriz de links L en 4D. En 3D no hay uno establecido.
- **Myrheim-Meyer y el perfil de intervalos suponen que la región es un intervalo causal (un diamante).** Fuera de eso dan dimensiones equivocadas.
- **Los criterios con tolerancia 2·SE + constante favorecen a la familia más ruidosa** (lo que hizo que E1′ pareciera discriminar).
- **Validar un muestreador en todo el rango de parámetros** (el error de β < 0 en P4 pasó porque solo se validó β ≥ 0).
- **Con pocas calibraciones, el percentil 95 es casi el máximo**: umbrales frágiles.

## 9. Preferencias de Sebastían para este proyecto

- Español, voz honesta y directa, sin lenguaje inflado; **no usar la raya (—)** ni reemplazarla con un exceso de dos puntos.
- Explicar con ejemplos simples (no es del área), y con imágenes cuando ayuda.
- Cuando algo requiere su opinión o una decisión que cambia la hipótesis, **preguntarle** en vez de suponerla.
- Presentar los resultados negativos con la misma visibilidad que los positivos.

## 10. Cómo empezar el chat nuevo

Sugerencia de primer mensaje en Claude Code, desde la carpeta del repositorio:

> «Lee `docs/CONTEXTO_NUEVO_CHAT.md` y `docs/handoff/HANDOFF_SIGUIENTE.md`. Actúa como orquestador con la metodología descrita: preregistro, subagente de código, auditoría y handoff. Estamos en F5′ (`docs/OMEGA_CONSEJO_REV9_F5PRIMA.md`); antes de ejecutar necesito decidir R9-1 y R9-2.»
