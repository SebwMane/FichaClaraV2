# Handoff para la siguiente sesión (orquestador o agente)

- **Estado al cierre:** 2026-10-12, rama `claude/hopeful-galileo-88u1l1`, HEAD tras los resultados de C0-4b.
- **Actualización 2026-10-12 (15):** diseño de E1′ entregado al Consejo (`OMEGA_CONSEJO_REV4_E1PRIMA.md`): regla de región que maximiza min(|pasado(a)|, |futuro(b)|), criterio con complementariedad, compuerta V′ y rama secundaria por rodilla. **No ejecutado; no congelado.** Decisiones R4-1 a R4-4 pendientes del usuario.
- **Actualización 2026-10-12 (14), Ω-E (ruta entrópica):** preregistro `OMEGA_E_PROPUESTA.md`. E0a PASA; E0b FALLA y E0b′ PASA (a = 0.347 ± 0.092, enmienda E-A1). E1 **INDETERMINADO**: la regla de región intrínseca produce intervalos descentrados y el control geométrico no muestra ley de área; KR casi la muestra (a = 0.436). Resultados en `OMEGA_E_RESULTADOS.md`. Siguiente posible: E1′ con regla de región que controle la forma, validada antes en sprinklings.
- **Actualización 2026-10-12 (13):** informe final escrito (`OMEGA_2_INFORME_FINAL.md`, con enlace al Claude Doc editable). Leer primero ese informe. Ningún trabajo nuevo autorizado hasta que el usuario elija camino.
- **Actualización 2026-10-12 (12), FIN DEL PROGRAMA DE PRUEBAS Ω-C:** piloto P4 SIN VENTANA en N = 64 y 128 (`OMEGA_2_P4_RESULTADOS.md`). Error A-22 del muestreador detectado y corregido; validación P4-A1b falló por autocorrelación y se rehízo con medias por bloques (P4-A2). Hallazgo inesperado: fase de acción alta con β < 0 (capas dispersas que explotan +9N₁). Todas las ramas de Ω-C quedaron probadas. Siguiente paso recomendado: informe final del programa y decisión del usuario sobre continuar Ω-C o volver a Ω-0.2.
- **Actualización 2026-10-12 (11):** L4b INCOMPATIBLE: sin borde espacial la dispersión del interior crece como N^0.84 a N^1.48 (no √N); el borde no era la causa (`OMEGA_2_L4B_RESULTADOS.md`; anexo 4D en `OMEGA_2_SINTESIS_2D.md`). Rama L incompatible en 4D en todo lo probado. Queda: rama euclídea 4D (no iniciada).
- **Actualización 2026-10-12 (10), FASE 4D:** R3-3 ratificada (d = 4). C4 (acción 4D y suavizada) ÉXITO. L4 INCOMPATIBLE en diamantes en los 7 esquemas de ε (`OMEGA_2_C4_L4_RESULTADOS.md`). Siguiente natural: **L4b**, repetir en una región sin borde espacial (toro) para separar volumen y borde. La rama euclídea 4D no se ha iniciado.
- **Actualización 2026-10-12 (9), CIERRE DEL CICLO 2D:** R3-6 ratificada. L-ε INCOMPATIBLE (D ∝ N^1.25): rama L cerrada en d = 2. Leer primero `OMEGA_2_SINTESIS_2D.md`. Decisiones abiertas del usuario: R3-3 (pasar a d = 4) o desarrollar movimientos no locales para la rama E en 2D.
- **Actualización 2026-10-12 (8):** Bloque A-S (sándwich): AS-2 PASA; AS-4 NO CONCLUYENTE. Lectura: el sándwich casi no se cancela (exacto ≈ cota ≈ ln del conteo), mientras la geometría con ε fijo se cancela fuertemente, así que con ε fijo la rama L no selecciona geometría (`OMEGA_2_BLOQUE_AS_RESULTADOS.md`). Pendiente del usuario: **R3-6, ¿ε ∝ 1/N (no localidad fija en volumen) en la rama L?**
- **Actualización 2026-10-12 (7):** Templado paralelo descartado por factibilidad (≈ 100 a 400 réplicas); se usó Wang-Landau (L-3c-EQ, §16). Hallazgo: el mínimo de S en 2-órdenes es un **sándwich** (capas a, 2, a), no la bicapa (EQ-A1). V falló una vez por presupuesto y pasó al repetirse (EQ-A2). L-3c-EQ **INDETERMINADO**: WL trabado en S = −504 sin llegar al sándwich (−538) en N = 32; N = 64 no ejecutado (declarado). Nuevo hueco en la rama L: A-1 no cubre la familia sándwich. Pendientes: (a) extender A-1/A-2 al sándwich (rama L, ejecutable ya); (b) movimientos no locales o population annealing (rama E, desarrollo).
- **Actualización 2026-10-12 (6):** Bloque A ejecutado: A-2a, A-2b y A-3 PASAN (γ = 1.41 y 1.56 < 2), A-1 en nivel 2 (`OMEGA_2_BLOQUE_A_RESULTADOS.md`). Huecos abiertos: rama E necesita un método de equilibrio (templado paralelo o integración termodinámica); rama L necesita cubrir familias no geométricas que no sean por capas y probar ε ∝ 1/N. Ningún encargo nuevo preregistrado.
- **Actualización 2026-10-12 (5):** L-3c-E″ (histéresis) INDETERMINADO por piso en N = 128; brecha 0, 3, ≥ 7 pasos (`OMEGA_2_L3CE_HYST_RESULTADOS.md`). El PASA de L-3c-E′ vale solo como estabilidad local. Equilibrio abierto: requiere templado paralelo o integración termodinámica (no preregistrado).
- **Actualización 2026-10-12 (4):** R3-5 ratificada (β = b/N). L-3c-E′ PASA: b_c = 4.00, 4.76, 6.73 para N = 32, 64, 128 (`OMEGA_2_L3CE_PRIME_RESULTADOS.md`), con riesgo de histéresis. Bloque A (vía analítica de la rama L) preregistrado en §13, sin ejecutar. **Próximos encargos autorizados, en este orden:** (1) preregistrar y correr L-3c-E″ (histéresis, partir de bicapa); (2) ejecutar Bloque A (A-1 derivación con cota, A-2 verificación, A-3 exponente de Var(S_ε/ε²)).
- **Actualización 2026-10-12 (3):** R3-4 ratificada. C1 (acción suavizada) ÉXITO. Rama L declarada bloqueada por el problema de signo (§11.5). L-3c-E: «ventana que se angosta», β_c·N ≈ 4.5 en N = 32 y 64 (`OMEGA_2_C1_L3CE_RESULTADOS.md`). Pendiente del usuario: **R3-5, ¿acoplamiento β = b/N en la rama E?** Si se firma, el siguiente encargo es L-3c-E′: barrido en b con N ∈ {32, 64, 128} y criterio de b_c estable, preregistrado antes.
- **Actualización 2026-10-12 (2):** R3-2 ratificada (ambos pesos, criterios separados). L-2c ejecutada (`OMEGA_2_L2C_RESULTADOS.md`): rama E viva, requiere L-3c-E; rama L bloqueada porque la acción sin suavizar cancela también a la geometría (L-iii INDETERMINADO). Próximo encargo autorizado: **L-3c-E** (falta preregistrarlo). Pendiente del usuario: **R3-4, ¿usar la acción suavizada en la rama L?**
- **Actualización 2026-10-12:** C0-4b ejecutado con ÉXITO (`OMEGA_2_C0_RESULTADOS.md` §5). La bandera de C0-4b es el clasificador habilitado; la de C0-4 queda retirada. El bloque C0 está **cerrado**: no queda trabajo de instrumentos pendiente que no requiera ratificación.
- **Leer primero, en este orden:** `docs/OMEGA_MASTER.md` (intuición y reglas), `docs/OMEGA_BLOQUE_L_RESULTADOS.md` (por qué S0 quedó cerrada), `docs/OMEGA_2_PROPUESTA_CAUSAL.md` (hipótesis Ω-C y enmiendas), `docs/OMEGA_2_C0_RESULTADOS.md` (resultados oficiales), `docs/handoff/AUDITORIA_C0.md`.

## Metodología vigente (no negociable)

1. **Roles separados.** El orquestador escribe el handoff con criterios congelados; el agente de código implementa solo eso; el orquestador audita como código ajeno, corrige y documenta en `docs/handoff/AUDITORIA_*.md`.
2. **Preregistro.** Criterios en el documento de decisión antes de cualquier corrida oficial. Los cambios posteriores van como enmiendas fechadas, con la columna «¿puede crear un positivo?».
3. **Pilotos declarados.** Entropía distinta de la oficial; no cuentan como evidencia; se declara para qué se usaron.
4. **Commits pequeños y ordenados:** documento → código → auditoría → corrida oficial (árbol limpio) → resultados.
5. **Escala de afirmación en 4 niveles:** demostrado, explicación propuesta, teorema externo, extrapolación no demostrada.
6. **Los resultados negativos se registran** y no se rescatan cambiando criterios.
7. **Nada de dinámica nueva sin ratificación del usuario.**

## Pendiente de decisión del usuario

| ID | Pregunta |
|---|---|
| R3-1 | ¿Adoptar Ω-C como hipótesis candidata? |
| R3-2 | **RATIFICADA: ambos, criterios separados.** Antes: peso de la compuerta: euclídeo e^{−βS}, lorentziano e^{iS}, o ambos. **Ver `OMEGA_2_C0_RESULTADOS.md` §3: con el euclídeo y el espacio de todos los órdenes, KR gana por mucho.** |
| R3-3 | **RATIFICADA: d = 4.** |
| R3-4 | **RATIFICADA.** Antes: ¿Ω-C pasa a usar la acción BD suavizada (escala de no localidad) para la rama lorentziana? Es un cambio de la hipótesis, no del instrumento |
| R3-6 | **RATIFICADA (resultado: INCOMPATIBLE).** Antes: ¿la rama L usa ε ∝ 1/N (escala de no localidad fija en volumen) en lugar de ε fijo? Cambia la hipótesis |
| R3-5 | **RATIFICADA.** Antes: ¿la rama E adopta el acoplamiento β = b/N (b fijo)? Cambia la hipótesis; se sugiere por β_c·N ≈ 4.5 en N = 32 y 64 |

## Trabajo de instrumentos

Cerrado. C0-4b pasó (ver resultados §5). Mejora opcional, no bloqueante: repetir C0-4b con N = 512 y 2048 para medir el margen en 4D.

## Trabajo que requiere ratificación

- **L-2c** (paisaje de la acción con el peso elegido en R3-2) y **L-3c** (estabilidad de un sprinkling 2D bajo Metropolis en 2-órdenes). Definidas en `OMEGA_2_PROPUESTA_CAUSAL.md` §5.

## Notas técnicas

- `omega/causality/` no depende de Ω-1.1. Ningún archivo de Ω-1.1 cambió.
- `tools/c0_validation.py --mode full` exige árbol limpio y tarda ~45 s.
- `omega/causality/exact.py` da E[N_k] y E[S_BD] exactos para 2-órdenes; úsese como prueba de regresión.
- Suite completa verificada el 2026-10-12 sobre `63b1e8a`, por grupos (el contenedor corta procesos largos): 700 passed, 0 failed (`tests/analytical` + `tests/test_[a-f]*` 398; `test_[g-o]*` 109; `test_[p-z]*` 131; `omega/experiments` 62). Ningún archivo de Ω-1.1 cambió.
