# R18-1: piloto técnico de muestreo multicanónico en la altura. Preregistro (2026-10-19)

- **Mandato:** dictamen del Consejo Ω, sesión 7. Siete de ocho perfiles apoyaron un piloto por altura bajo condiciones; el octavo pidió cerrar antes ES-C0b-R y definir objetivo y soporte. Firma del usuario: «Apruebo R18-1».
- **Tipo:** piloto **técnico y no oficial**. **No es evidencia** de selección física, y sus resultados no se interpretan como tal.
- **ES-C0b con el movimiento R queda retirada del plan.** El código se conserva como desarrollo verificado (commit c61dd79).

## 0. Resumen del dictamen de la sesión 7 (texto pegado por el usuario)

- **Los ceros de aplicabilidad son de R, no del movimiento local** de agregar y quitar relaciones.
- Las mesetas de ES-C0 respaldan una dificultad práctica con N de 64 a 128, pero **no prueban que el espacio esté desconectado** ni se extrapolan a otros tamaños.
- **La altura es una coordenada razonable para probar, pero por sí sola no identifica geometría** (Brightwell y Gregory; Kleitman y Rothschild).
- **Wang-Landau:** la adaptación debe quedar separada de la producción. Se congela W(h), se muestrea con el sesgo fijo y se repondera por e^{W(h)} (Wang y Landau 2001). Eso no compensa la falta de soporte ni las pocas muestras efectivas.
- **Antes de correr,** preregistrar el objetivo, los sectores, el presupuesto y los criterios de parada. Exigir viajes repetidos entre sectores además de recorrer alturas: un histograma plano no basta.
- **Si falla:** rediseñar R (opción ii) o cerrar la afirmación de mezcla en ese rango (opción iii).

## 1. Objetivo y soporte

- **Medida objetivo:** la de conteo sobre órdenes de N elementos con etiquetado natural fijo (β = 0). Es una **prueba del algoritmo**, no una medida física de Ω.
- **Soporte:** todos esos órdenes, sin exclusiones.
- **Movimiento:** solo el local de ES-C0 (quitar un link, o agregar una relación que deje el conjunto transitivo), ya verificado de forma exacta.
- **Sesgo:** se muestrea π_W(C) ∝ exp(−W(h(C))), donde h es la altura (la cadena más larga).
  - W(h) se aprende en la ventana h ∈ [2, h_max(N)].
  - Fuera de la ventana, W toma el valor del borde más cercano. Así **no se excluye ningún estado** del soporte.
  - Al reponderar, cada muestra pesa e^{W(h)}.
- **h_max(N):** 6 con N = 32, 7 con N = 64 y 8 con N = 128. Cubren con margen las alturas de las mesetas de ES-C0 (2 a 6, sin contar la cadena total congelada).

## 2. Sectores y viajes

- **Sector bajo:** h ≤ 2. **Sector alto:** h ≥ 5.
- **Un viaje completo** es ir del sector bajo al alto y volver al bajo. En ES-C0 las mesetas se repartían entre altura 2 y altura 3 a 5, y ninguna cadena las cruzaba.

## 3. Fases y presupuesto

| Fase | Qué hace | Presupuesto |
|---|---|---|
| A, adaptación | Wang-Landau con ln f desde 1, partido por 2 cada vez que el histograma de h en la ventana queda plano (cada barra ≥ 80% de la media), hasta ln f < 1e-4 o hasta agotar el presupuesto. Una sola cadena por N. | 3000 barridos con N = 32 y 64; 1500 con N = 128 |
| B, producción | W(h) **congelado**. Cuatro cadenas desde los cuatro inicios de ES-C0 (orden vacío, cadena total, sprinkling 4D y KR natural). | 2000 barridos con N = 32 y 64; 1000 con N = 128 |

- La fase A no se usa para estimar nada.
- Si la fase A agota el presupuesto sin llegar a ln f < 1e-4, se congela el W alcanzado y se declara.

## 4. Criterios (todos deben cumplirse para «mezcla», por N)

1. **R̂ dividido < 1.05 y ESS ≥ 200** en la fracción de relaciones, la altura y S₄ (la acción sin suavizar, solo como observable), en la mitad final de las 4 cadenas sesgadas.
2. **Viajes:** al menos 10 viajes completos por cadena en la producción.
3. **ESS de reponderación** ≥ 200 por N: (Σw)²/Σw², sumando las 4 cadenas.
4. **Concordancia entre cadenas:** la media reponderada de la fracción de relaciones de cada cadena queda dentro de 2 errores (por bloques) de la media conjunta.
5. **Validación externa con N = 32:** la media reponderada de la fracción de relaciones coincide, dentro de 2 errores combinados, con la de ES-C0 en N = 32, b = 0 y ε = 1. Es la única celda de órdenes generales que mezcló.

**Si se cumplen todos en algún N,** el método queda validado técnicamente en ese N. **Si fallan,** se informa «no convergido» y la decisión entre (ii) y (iii) vuelve al usuario.

## 5. Implementación y verificación previa

- **Altura incremental exacta.**
  - Al agregar i < j: h′ = max(h, prof(i) + 1 + coprof(j)).
  - Al quitar un link: se recalcula solo si ese link estaba en alguna cadena de longitud h.
  - Se verifica contra `mcmc_general.height` en cada paso de cadenas de prueba.
- **Balance detallado del núcleo sesgado:** enumeración exacta con N ≤ 5 y un W arbitrario.
- **Reponderación:** con N ≤ 5, el promedio reponderado de una cadena larga debe reproducir el promedio exacto de la medida de conteo.

## 6. Convenciones

- **EXPERIMENT_ID 331.** Entropía 20261108 para la fase A y 20261109 para la fase B.
- **Próximos libres:** EXPERIMENT_ID 332 y entropía 20261110.
- **Herramienta:** agente de código con auditoría; corrida sobre árbol limpio, aunque es no oficial, para que sea reproducible.

## 7. Aclaración antes de la corrida (2026-10-19)

- **Convención de Wang-Landau:** W[h] += ln f después de **cada propuesta** (aceptada, rechazada o inválida), solo dentro de la ventana. Es la convención estándar de Wang y Landau (2001). La planilla de §3 (ln f desde 1, partido por 2 con histograma plano) no se cambia.
- **Hallazgo de la prueba corta (no oficial), declarado antes de correr:** la planilla no alcanzó la planitud en ninguna prueba corta. Tampoco lo hizo una variante más suave (ln f por barrido), que no se adopta. Con N = 32, la cadena nunca volvió a h ≤ 2 después de los primeros barridos. La previsión es que la fase A congele un W no convergido y que el veredicto sea «no convergido». Se corre igual, porque es el resultado preregistrado; no se cambia la planilla a la vista de la prueba corta.
