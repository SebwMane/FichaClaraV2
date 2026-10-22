# Consejo Ω, Revisión 18: cierre de la sesión 7 y decisión final entre (ii) y (iii) (2026-10-19)

- **Mandato:** el usuario delegó en el orquestador la decisión final («Que el orquestador termine la sesión del consejo y delibere entre (II) y (III) de manera final»).
- **Insumos:**
  - acta de la sesión 7 (`docs/consejo_rev18/ACTA_SESION7_NOTA.md`);
  - resultados de ES-0, ES-C0, ES-C0b (prueba corta) y R18-1.

## 1. Lo que dicen los datos (resumen)

| Prueba | Método | N | Resultado |
|---|---|---|---|
| ES-0 | Metropolis local, Ω-E | 32 y 64 | N = 64 no mezcla, ni con λ = 0 |
| ES-C0 | Metropolis local, Ω-C (ε = 1 y 0.25) | 32, 64 y 128 | solo mezcla N = 32; con N ≥ 64, mesetas de capas que no se cruzan, incluso con b = 0 |
| ES-C0b (prueba corta) | reinserción R con pasados principales | 64 y 128 | R inaplicable (0 de 40 365 y 0 de 163 198 propuestas) |
| R18-1 | Wang-Landau en la altura más movimiento local | 32, 64 y 128 | no convergido; 0 viajes entre sectores; falla la validación con ES-C0 |

## 2. Cómo se aplican las reglas del acta

- **Regla de parada:** R18-1 falló en cobertura, cruces, R̂ y validación. Corresponde **cerrar según (iii) la afirmación de mezcla en el rango estudiado**.
- **Regla de reserva de (ii):** se cumple su condición de entrada, porque los sectores no se comunicaron mediante cambios de altura: el sesgo existió, pero la cadena no podía bajar de altura con movimientos de una relación. Eso habilita (ii), pero no lo obliga.

## 3. Resolución final del orquestador

**Se adopta (iii), y (ii) queda en reserva con condiciones de reingreso explícitas.**

1. **Cierre (iii), resultado de factibilidad (nivel 1, sobre los instrumentos):**
   - Con los movimientos disponibles (local, reinserción principal, sesgo en la altura), la medida de órdenes con etiquetado natural **no se pudo muestrear de forma fiable con N de 64 a 128**.
   - Esto vale con λ = 0 (Ω-E), con b = 0 y con b ≠ 0 (Ω-C, ε = 1 y 0.25).
   - **Consecuencia:** H_sel queda **indeterminada** en Ω-E y en Ω-C en este rango. **No refutada.** No se extrapola a otros N, otras medidas ni otros movimientos.
2. **Por qué no (ii) ahora:**
   - Un movimiento colectivo para pasados y futuros arbitrarios exige probabilidades de propuesta inversas calculables. Contar pasados equivale a contar anticadenas, que en general es intratable, y no hay un diseño concreto que cumpla eso sin un desarrollo de algoritmos propio.
   - Tres intentos seguidos encontraron la misma barrera. Un cuarto sin un diseño nuevo concreto sería insistir sin hipótesis.
   - Por la orden del usuario («supervivencia de la metodología sin sobreajuste»), es mejor cerrar honestamente que seguir ajustando el muestreador hasta que algo pase.
3. **Condiciones de reingreso de (ii)** (cualquier propuesta futura debe traerlas antes de escribir código):
   - un movimiento colectivo con probabilidad inversa **calculable y verificada** por enumeración exacta;
   - aplicabilidad medida en los órdenes en capas de ES-C0 (no solo desde los inicios), con un umbral fijado antes;
   - preregistro propio, con presupuesto y los mismos criterios de R18-1 §4.
4. **Lo que sigue en pie y no depende del muestreo:**
   - Ω-E / E1′ como herramienta descriptiva en 2D;
   - la corrección de √ε;
   - las verificaciones exactas de los muestreadores;
   - y las pruebas sobre órdenes **dados**, que no requieren muestrear la medida: R14-B v2 (grano grueso, H_E) y R14-C v2 (observables multivariables), preparadas y sin ejecutar.

## 4. Estado del programa tras esta revisión

- **H_sel:** indeterminada en el rango computable. ES-E1 y ES-C1 quedan **cerradas** con esta definición. No se reabren sin un muestreador que pase R18-1 §4.
- **H_E y H_L:** sus pruebas sobre órdenes dados siguen disponibles (R14-B v2 y R14-C v2), a la espera de la orden del usuario.
