# Handoff para la siguiente sesión (orquestador o agente)

- **Estado al cierre:** 2026-10-12, rama `claude/hopeful-galileo-88u1l1`, HEAD tras los resultados de C0-4b.
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
| R3-3 | ¿Mantener d = 2 hasta cerrar la compuerta? (se usa d = 2 por defecto, sin firma) |
| R3-4 | **RATIFICADA.** Antes: ¿Ω-C pasa a usar la acción BD suavizada (escala de no localidad) para la rama lorentziana? Es un cambio de la hipótesis, no del instrumento |
| R3-5 | ¿La rama E adopta el acoplamiento β = b/N (b fijo)? Cambia la hipótesis; se sugiere por β_c·N ≈ 4.5 en N = 32 y 64 |

## Trabajo de instrumentos

Cerrado. C0-4b pasó (ver resultados §5). Mejora opcional, no bloqueante: repetir C0-4b con N = 512 y 2048 para medir el margen en 4D.

## Trabajo que requiere ratificación

- **L-2c** (paisaje de la acción con el peso elegido en R3-2) y **L-3c** (estabilidad de un sprinkling 2D bajo Metropolis en 2-órdenes). Definidas en `OMEGA_2_PROPUESTA_CAUSAL.md` §5.

## Notas técnicas

- `omega/causality/` no depende de Ω-1.1. Ningún archivo de Ω-1.1 cambió.
- `tools/c0_validation.py --mode full` exige árbol limpio y tarda ~45 s.
- `omega/causality/exact.py` da E[N_k] y E[S_BD] exactos para 2-órdenes; úsese como prueba de regresión.
- Suite completa verificada el 2026-10-12 sobre `63b1e8a`, por grupos (el contenedor corta procesos largos): 700 passed, 0 failed (`tests/analytical` + `tests/test_[a-f]*` 398; `test_[g-o]*` 109; `test_[p-z]*` 131; `omega/experiments` 62). Ningún archivo de Ω-1.1 cambió.
