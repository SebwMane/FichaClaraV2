# Handoff para la siguiente sesión (orquestador o agente)

- **Estado al cierre:** 2026-10-12, rama `claude/hopeful-galileo-88u1l1`, HEAD tras los resultados de C0-4b.
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
| R3-3 | ¿Mantener d = 2 hasta cerrar la compuerta? (se usa d = 2 por defecto, sin firma) |
| R3-4 | **RATIFICADA.** Antes: ¿Ω-C pasa a usar la acción BD suavizada (escala de no localidad) para la rama lorentziana? Es un cambio de la hipótesis, no del instrumento |
| R3-6 | ¿La rama L usa ε ∝ 1/N (escala de no localidad fija en volumen) en lugar de ε fijo? Cambia la hipótesis |
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
