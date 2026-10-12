# Handoff para la siguiente sesión (orquestador o agente)

- **Estado al cierre:** 2026-10-12, rama `claude/hopeful-galileo-88u1l1`, HEAD tras los resultados de C0-4b.
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
| R3-2 | Peso de la compuerta: euclídeo e^{−βS}, lorentziano e^{iS}, o ambos. **Ver `OMEGA_2_C0_RESULTADOS.md` §3: con el euclídeo y el espacio de todos los órdenes, KR gana por mucho.** |
| R3-3 | ¿Mantener d = 2 hasta cerrar la compuerta? |

## Trabajo de instrumentos

Cerrado. C0-4b pasó (ver resultados §5). Mejora opcional, no bloqueante: repetir C0-4b con N = 512 y 2048 para medir el margen en 4D.

## Trabajo que requiere ratificación

- **L-2c** (paisaje de la acción con el peso elegido en R3-2) y **L-3c** (estabilidad de un sprinkling 2D bajo Metropolis en 2-órdenes). Definidas en `OMEGA_2_PROPUESTA_CAUSAL.md` §5.

## Notas técnicas

- `omega/causality/` no depende de Ω-1.1. Ningún archivo de Ω-1.1 cambió.
- `tools/c0_validation.py --mode full` exige árbol limpio y tarda ~45 s.
- `omega/causality/exact.py` da E[N_k] y E[S_BD] exactos para 2-órdenes; úsese como prueba de regresión.
- Suite completa verificada el 2026-10-12 sobre `63b1e8a`, por grupos (el contenedor corta procesos largos): 700 passed, 0 failed (`tests/analytical` + `tests/test_[a-f]*` 398; `test_[g-o]*` 109; `test_[p-z]*` 131; `omega/experiments` 62). Ningún archivo de Ω-1.1 cambió.
