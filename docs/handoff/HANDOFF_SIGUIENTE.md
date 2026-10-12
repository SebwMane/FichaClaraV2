# Handoff para la siguiente sesión (orquestador o agente)

- **Estado al cierre:** 2026-10-11, rama `claude/hopeful-galileo-88u1l1`.
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

## Trabajo listo para encargar (no requiere ratificación, solo instrumentos)

**C0-4b: discriminador de órdenes no geométricos con validez de constructo.** Propuesta del orquestador para que el Consejo la congele antes de programarla:
- Añadir un tercer estimador de naturaleza distinta, por ejemplo el escalamiento de la cadena más larga entre los extremos del intervalo máximo, L ≈ m_d · N^{1/d}, medido en dos tamaños.
- Exigir además que la distribución de abundancias de intervalos N_k / E[N_k] se mantenga dentro de una banda frente a la de un sprinkling de la dimensión estimada (la percolación densa tiene un perfil de N_k muy distinto, ver S_BD = +4413).
- Controles negativos obligatorios: KR, percolación con dos valores de p, y un orden «capa gruesa» intermedio.
- Validación antes de correr, al estilo L-A5: si el clasificador no separa sprinklings sin evolucionar de los controles, se vuelve al Consejo.

## Trabajo que requiere ratificación

- **L-2c** (paisaje de la acción con el peso elegido en R3-2) y **L-3c** (estabilidad de un sprinkling 2D bajo Metropolis en 2-órdenes). Definidas en `OMEGA_2_PROPUESTA_CAUSAL.md` §5.

## Notas técnicas

- `omega/causality/` no depende de Ω-1.1. Ningún archivo de Ω-1.1 cambió.
- `tools/c0_validation.py --mode full` exige árbol limpio y tarda ~45 s.
- `omega/causality/exact.py` da E[N_k] y E[S_BD] exactos para 2-órdenes; úsese como prueba de regresión.
- No se ejecutó la suite completa del repo (~5 min); sí los tests de arquitectura y `tests/test_causality.py`.
