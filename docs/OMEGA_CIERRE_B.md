# Ω — Cierre del programa de selección dimensional (problema B)

- **Fecha:** 2026-10-07.
- **Rama:** `claude/omega-cierre-b`; congelada como `claude/omega-cierre-b-congelado`.
- **Decisión:** la toma el cerebro por delegación expresa del usuario. Aplica la condición de cierre que fijó el propio Consejo.

## 1. La condición del Consejo y por qué se cumple

El Consejo escribió: «Si W-0 demuestra que el número de hilos siempre termina siendo parámetro, condición inicial, conservación o elección discreta de ley, el Consejo podrá decir algo mucho más fuerte […] y ahí sí, si no aparece un principio independiente que rompa esa conclusión, yo recomendaría cerrar el programa de selección dimensional de Ω».

| Requisito | Dónde se cumple | Nivel |
|---|---|---|
| Tasas ⇒ dial | W-T1 (Poisson de media 2β/γ); W-T2 (lema del dial) | 1 |
| Conservación o inicio ⇒ condición inicial | W-T3; E6-D, clase CONSERVADOR; SQ conserva Z³ (R1-1) | 1 |
| Umbrales sin parámetro ⇒ elección de ley | W-T4; W-1 descriptivo en N₂ (corte en 2 o en 4 según la ley) | 2 (más dato de nivel 1 en N₂) |
| Atractor genuino (caso IV) | W-T7, dilema de independencia: sin comunicación → inicio o explosión; desfase acotado → w = 1; desfase no acotado → aridad del evento (ley) o escala (dial) | 2 |
| Resolución empírica de la ruta restante | W-1: W1-D en el dominio prerregistrado; no universal (peine-2) | 1 |
| Principio independiente que rompa la conclusión | Ninguno propuesto | — |

**«Demuestra» se lee con el estándar del programa.** Los pasos de nivel 1 están demostrados. Los de nivel 2 son argumentos completos dentro de su planteamiento, pero no teoremas sobre toda dinámica posible. El Consejo ya aceptó que el cierre **no** es un teorema de imposibilidad.

## 2. Enunciado de cierre

Une las formulaciones aprobadas por el Consejo (W-0 §7, W-1 §5 y la sesión sobre R3):

> **En las familias estructurales examinadas, la dimensionalidad de Ω queda reducida a un grado de libertad (el número de direcciones independientes que conmutan: anchura causal, rango) cuya selección no puede emerger sin introducir una fuente equivalente de información: parámetro, estado inicial, conservación o ley elegida.**
>
> Las rutas sin parámetro continuo (umbrales de paseos sobre el propio grafo):
> - no son resolubles de forma estable a ≤ 10⁶ nodos en el dominio prerregistrado;
> - cuando se ven, su corte depende de la ley;
> - no son universales.
>
> Un atractor genuino exigiría universalidad sobre clases de reglas, y el dilema de independencia indica que, con información local, toda regulación del número de hilos cae en dial, inicio o ley, o colapsa a una sola dirección.
>
> **Medir dimensión no es generarla.**
>
> Esto **no es un teorema de imposibilidad.** Identifica dónde falla la ambición explicativa de Ω: no en producir estructuras, sino en explicar por qué una de ellas tiene un número concreto de direcciones independientes.

## 3. Cadena de evidencia

| Fase | Resultado | Lección |
|---|---|---|
| L-DIM-1 | La selección estática por un funcional es un dial (d* = f(θ, N)) | 17, dial |
| L-ARQ-0 | Lema del dial: solo se eligen vértices de la envolvente; la meseta es una arista | — |
| R1-0/R1-1 | La regla local ciega a d (SQ) nuclea pero no coalesce; conserva Z³ | 22 |
| E6 | Codificar d es relativo al panel; puerta estática calibrada | 23 |
| R3-0 | La coalescencia sin defectos fija d = anchura asintótica (Birkhoff–Dilworth) | 24, 25 |
| W6 | Regla de uso de RC-3 a tres escalas (grado saturado) | — |
| W-0 (+ W-T7) | C1 y C2 cerradas; C3 → elección de ley; caso IV → dilema de independencia | 26 |
| W-1 | W1-D; corte dependiente de la ley; recurrencia ≠ dimensión | 27 |

## 4. Lo que el cierre no dice

- No dice que Ω no pueda producir geometrías: produce redes de cualquier d dado un recuento.
- No dice que ninguna dinámica posible seleccione d. El cierre es sobre las familias examinadas, con un argumento de nivel 2 para el caso IV.
- No dice nada sobre d = 3 en particular. Ningún resultado del programa privilegia ni descarta un valor concreto.

## 5. Lo que sigue en pie

- **Problema A (exclusión de degenerados):** instrumentos validados en sus dominios.
  - RC-3, juez A−, con los puntos ciegos W1–W6 declarados y las reglas de uso W5 y W6.
  - E6-S v1.1, con validación débil; su deuda es un panel de grafos reservado y ciego.
  - El mapa de exclusiones.
- **Metodología:** 27 lecciones, el mapa negativo y la disciplina de prerregistro, potencia y no retroactividad.
- **Pregunta abierta, formulada con precisión:** ¿existe un principio pregeométrico, independiente del resultado, que fije la ley de interacción (su aridad, su carácter local o no local) y, a partir de ella, el número de direcciones independientes?

## 6. Criterios de reapertura (congelados)

El programa B solo se reabre si se presenta, **antes** de mirar ningún resultado dimensional, al menos una de estas tres cosas:

1. **Un principio independiente que fije una ley o una aridad**, justificado por razones ajenas a la dimensión que produce. Debe pasar S0–S2 de E6 y la cláusula α de E2.
2. **Una clase de reglas propuesta por un principio** (no elegida por su resultado) en la que se conjeture universalidad de w*. Se prerregistraría con nulos, calibración de E6-D y W6.
3. **Un mecanismo de información no local** con justificación propia, que escape al dilema de independencia. Exigiría revisar explícitamente el requisito de localidad de C0.

No reabren el programa:
- una regla nueva que dé un atractor concreto (la «última trampa» del Consejo);
- repetir W-1 con más tamaño (valor de la información nulo para esta decisión; W-1 §5).
