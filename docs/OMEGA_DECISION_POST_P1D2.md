# Ω — Decisión del cerebro sobre la propuesta del Consejo (posterior a P1-D.2)

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-p1d3`.
- **Contexto:** el Consejo (11 especialistas) propuso la siguiente fase de simulaciones. Antes de decidir, contrasté sus premisas con los registros del repositorio. Un agente de solo lectura hizo la auditoría y verifiqué a mano el dato que salió contradictorio.

## 1. Auditoría de las premisas del Consejo

| Premisa del Consejo | Estado real en los registros | Fuente |
|---|---|---|
| «Terminar F1 con N = 512/729» | **Ya terminado.** F1 + F1b: F1-NEGATIVO en 10/10 celdas. Con N = 729, **0/24** pasan el certificado: 24/24 con F4 y F9, 19/24 también con F5 | `results/c0_f1/`, `OMEGA_C0_RESULTADOS.md` §5 |
| «Algunos estados C0 muestran escalas 2.5–3.7» | Cierto con N = 343 (D_L 2.5–3.7, D_s 2.6–3.3). **Con N = 729 bajan**: D_s 1.77–2.76 (mediana 2.20), D_eff 2.22–2.62 | `OMEGA_C0_RESULTADOS.md` §2, §5 |
| «157/253 sin converger: continuar hasta 100 000 pasos» | **Ya hecho (D-1):** 156/157 siguen locales a 100 000 pasos | `results/c0_redteam/` |
| «30 celdas pasan R3 + R8» | Sí; después R6 dejó **24 CANDIDATO-C0**, que luego fallaron F1 | ídem |
| «Θ > 0 en segunda prioridad» | **Ya ejecutado:** FRÁGIL-1/N (la frontera Θ_c ≈ 0.1·Θ_N con Θ_N ∝ 1/N) | `OMEGA_FASE2_RESULTADOS.md` |
| «Red-team de los candidatos C0 (recableado, grados)» | **Ya hecho:** control nulo 36/40 separado; retazos y atajos en P1-D | `results/c0_battery/null_control.json` |
| «¿C0 se homogeneiza a escala mayor?» | **No:** ρ = 1.64–2.53 (NO_HOMOGENEIZA en 18/18 evaluables, N = 729) | `results/p1d2/` |

**Estadística pedida por el especialista en inferencia** (ya disponible, sin simular nada):

| Condición | P(local) / P(certificado) |
|---|---|
| P(local \| inicio genérico, N = 216) | 253/675 (U 95/225, E 54/225, R 104/225) |
| P(persiste \| local sin converger) | 156/157 a 100 000 pasos |
| P(local), R6 | N = 125: 110/213; N = 343: 208/213 |
| P(local), N = 729 | 24/24 |
| **P(certificado), N = 729** | **0/24** |

*Nota de auditoría:* el agente lector informó «24/24 pasan el certificado con N = 729». Lo verifiqué a mano y es falso: 0/24 pasan. Lo corrijo aquí.

## 2. Análisis crítico de la propuesta del Consejo

**En qué tiene razón el Consejo** y lo adopto:
1. **No diseñar todavía una Ω-2.0.**
2. **No usar Θ como salvavidas.**
3. **Separar las preguntas A, B y C** (¿hay estructura?, ¿geometría?, ¿3D?) y las categorías **gruesa, variedad y clase**.
4. **N grande es obligatorio para el instrumento.** P1-D.2 mostró que con N = 4096 un cruce WS imita una dimensión de 2.7 y que el RGG3 todavía no ha convergido (D_B2 = 2.56 → 2.78 → 2.90).
5. **Controles inhomogéneos (S², caja con borde).** Es el límite que P1-D.2 dejó abierto.
6. **No modificar el certificado tras ver los resultados.**

**En qué se apoya en premisas ya superadas:**
- Su árbol de decisión empieza con «F1: ¿sobrevive el escalado? → NO: cerrar esa rama». **Esa puerta ya está respondida, y la respuesta es NO.** Lo mismo pasa con D-1, R6, Θ > 0 y el control nulo.
- Volver a correr C0 a N creciente (2·10³ → 2·10⁴) costaría O(N³) por paso: ~10 h por corrida con N = 2000 y miles de horas con N = 2·10⁴. Además respondería una pregunta que los datos ya contestaron: la dimensión de C0 **baja** con N, y C0 no homogeneiza. Sería exactamente el rescate que el propio Consejo rechazó tras la Fase 2.

**Lo que el Consejo no dijo y considero vital:**
- **El cuello de botella del programa es ahora el instrumento, no la dinámica.** Cada fase (L, C0, P1-D, P1-D.2) acabó revelando un defecto de la medida:
  - D ≈ 3 a N pequeño era un efecto de tamaño;
  - la meseta de D_B la engaña un cruce;
  - el Nivel I falla con la periodicidad y probablemente con la inhomogeneidad;
  - el certificado solo está calibrado para la clase 3.

  Sin un instrumento validado a N grande, cualquier funcional nueva volverá a producir ambigüedades. **P1-D.3 es por tanto la máxima prioridad, pero como calibración del instrumento, no como rescate de C0.**
- **La pregunta del Consejo es la correcta y la hago mía:** ¿qué mecanismo mantiene la localidad a escala pequeña y produce comportamiento dimensional a escala grande, sin colapsar en cliques ni volverse expansor? Los datos ya dicen algo sobre ella:
  - las acciones locales en (k, c) no pueden seleccionarlo (C0-T1/T2);
  - los conteos hasta segundo orden no distinguen C0 de RGG3;
  - la diferencia es mesoscópica (P1-D).

  La próxima funcional, o P3, debe diseñarse **contra** el control de retazos (geometría local idéntica, pegado incoherente). Ese es el adversario que encarna la pregunta.

## 3. Decisión

| # | Acción | Decisión | Estado |
|---|---|---|---|
| 1 | Terminar F1 | Ya está terminado (NEGATIVO) | — |
| 2 | Continuar los estados C0 hasta convergencia y llevarlos a N creciente | **No.** D-1 y F1 ya lo respondieron; sería un rescate | Cerrado |
| 3 | **P1-D.3**: instrumento a N ≈ 2·10⁴–5·10⁴, Nivel II por convergencia, S², caja, gradiente, retículos de cliques, WS, retazos, atajos y categorías gruesa / variedad / clase | **Sí, máxima prioridad** | Prerregistrado (`OMEGA_P1D3_PRERREGISTRO.md`); en ejecución |
| 4 | Certificado sin cambios después de ver los datos | Sí; Nivel III solo con análogos de N ≈ 729 | Incluido en P1-D.3 §5 |
| 5 | Θ > 0 adicional | No | — |
| 6 | Nueva funcional o P3 | **Después** de P1-D.3, en sesión del Consejo, con el instrumento validado y diseñada contra el control de retazos | Pendiente |

**Reparto de trabajo** (pedido por el usuario):
- **Opus:** análisis, decisión y prerregistro.
- **Haiku:** auditoría de solo lectura de los registros.
- **Sonnet:** implementación del código de P1-D.3, que reviso antes de ejecutar la corrida oficial.
