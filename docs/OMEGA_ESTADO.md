# Ω — Estado del programa (2026-10-07)

Documento de entrada para cualquier trabajo futuro. Resume y enlaza; no sustituye a los documentos de cada fase.

## 1. Estado en una tabla

| Problema | Estado | Documento |
|---|---|---|
| **B: selección dimensional** | **CERRADO** (reducción estructural negativa; no es un teorema de imposibilidad). Reapertura solo por los criterios congelados | `OMEGA_CIERRE_B.md` |
| **A: exclusión de degenerados** | **ABIERTO, sin candidato.** Los instrumentos están listos, pero ninguna dinámica de Ω ha producido una estructura no excluida y W6-válida | `OMEGA_MAPA_EXCLUSIONES.md` |
| Dinámica nueva | No autorizada por ahora (Consejo) | — |

## 2. Instrumentos vigentes y su dominio

| Instrumento | Uso | Dominio y límites |
|---|---|---|
| RC-3 (juez A−) | Excluye degenerados (X1–X4) | Puntos ciegos W1–W6 declarados; no certifica ni estima d |
| W5 | Validez de un par de tamaños | E/N con cambio < 25 % |
| W6 | Validez de una lectura de dimensión finita | Tres tamaños, misma estructura, grado saturado ≤ 5 %. Validada solo frente a la familia J3 |
| E6-S v1.1 | Puerta estática de no codificación | Grafos sin etiquetas; validación débil. Deuda: panel de grafos reservado y ciego antes del primer candidato de tipo grafo |
| E6-D | Brazo dinámico (olvido, dial, orden, representación, tamaño) | Diseño congelado; su calibración dinámica va con el primer candidato |
| Preflight A–H + E1–E6 | Requisitos de entrada | `OMEGA_PREFLIGHT_GEOMETRICO.md` |
| Trazabilidad de d | Toda estructura W6-válida declara la fuente de su d | `OMEGA_CIERRE_B.md` §7 |

## 3. Qué tendría que traer una propuesta futura

**Para A** (un candidato que produzca estructura):
1. Prerregistro antes del código, con cálculo de potencia (lección 25).
2. Preflight y auditoría escrita E6 S0–S2. Si es un grafo, se paga antes la deuda del panel reservado.
3. RC-3 con W5 y W6, más la trazabilidad de d.
4. E6-D completo si la fuente de d es «desconocida».

**Para reabrir B:** un principio independiente que fije la ley o su aridad, una clase de reglas derivada de un principio, o un mecanismo no local justificado por razones ajenas a la dimensión que produce. En los tres casos, presentado antes de mirar resultados dimensionales.

## 4. La pregunta abierta, formulada con precisión

> ¿Existe un principio pregeométrico, independiente del resultado, que fije la ley de interacción (su aridad, su carácter local o no local) y, a partir de ella, el número de direcciones independientes que conmutan?

Medir d no es explicar d. Producir una geometría no es explicar por qué tiene esa dimensión.

## 5. Ramas congeladas recientes

`claude/omega-e6-congelado` · `claude/omega-r3-0-congelado` · `claude/omega-w6-congelado` · `claude/omega-w0-congelado` · `claude/omega-w1-congelado` · `claude/omega-cierre-b-congelado`

Las 28 lecciones y las correcciones (incluidas las predicciones fallidas del cerebro) están en `OMEGA_MAPA_NEGATIVO.md`.
