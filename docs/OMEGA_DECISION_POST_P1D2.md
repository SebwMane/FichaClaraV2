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

## 4. Resultado de P1-D.3 y decisión sobre cómo seguir

**P1-D.3: INSTRUMENTO-VÁLIDO** (sensibilidad 15/15, especificidad 33/33, con los tres umbrales; el Nivel III rechaza los retículos de cliques 3D). Detalle en `OMEGA_P1D3_RESULTADOS.md`.

**Alcance:** espacios cerrados homogéneos o casi homogéneos. **No reconoce geometrías con borde**, que dan falso negativo. «D ≈ 3 que converge» puede fabricarse con un anillo y 0.1 % de atajos, y solo el Nivel I y la discordancia con D_s lo desenmascaran.

### Decisión del cerebro

1. **Se cierra la fase de instrumentación.** La batería (Niveles I y II y certificado como Nivel III, con su alcance declarado) es ahora obligatoria para toda dinámica futura. No habrá más rondas P1 salvo que una dinámica produzca un caso fuera de alcance (por ejemplo, con borde).
2. **No se simula nada más hasta que el Consejo ratifique un mecanismo.** Es lo que pidió el propio Consejo («solo después de ese informe, convocar al Consejo para diseñar la siguiente funcional»). El cuello de botella ya no es la medida: es la **falta de un mecanismo de coherencia mesoscópica**.
3. **Requisitos de diseño** para el siguiente mecanismo, derivados de todo lo medido:

   | # | Requisito | Origen |
   |---|---|---|
   | R1 | Debe producir **localidad** sin cliques ni vacío | C0 lo consiguió; S0/Ω-B no |
   | R2 | Debe poder distinguir, energética o dinámicamente, una geometría de los **retazos**, que tienen la misma estadística local y un pegado incoherente | C0-T1/T2; conteos (k, c, q) de C0 ≈ RGG3 |
   | R3 | La «cercanía» debe salir de la propia red, no de coordenadas | Pregunta 3 del Consejo para P3 |
   | R4 | La escala no puede imponer la dimensión: criterio de meseta C1-D en parámetros | Fase 2 §1.3 |
   | R5 | A Θ > 0 debe tener un coste por relación que no se anule con N, o una ligadura | C0-Θ FRÁGIL-1/N |
   | R6 | Se evalúa con la batería P1-D.3 y debe superar el control de retazos **como entrada**: un mecanismo correcto debería *curar* un retazo, no conservarlo | Propuesta nueva |

4. **Candidatos que propongo al Consejo** (no ejecutados; ninguno tiene código):
   - **Ω-D, cercanía por difusión autoconsistente (recomendado primero).** Las relaciones se reajustan según el núcleo de calor de la propia red, K_τ = exp(−τL). Dos nodos se acercan si sus distribuciones de difusión a tiempo τ se solapan. Por qué lo recomiendo:
     - responde a R3 (la cercanía nace de la red);
     - es mesoscópico por construcción (R2), porque el núcleo de calor ve la escala √τ, no solo vecinos comunes;
     - permite el test R6 de forma directa: ¿un retazo, al iterar, se cura o se congela?
     - Riesgos que hay que analizar antes de simular: τ podría fijar la escala o la dimensión (R4), puede colapsar a cliques (como C0) o a expansores, y su punto fijo trivial es el grafo completo.
     - Primer paso: análisis de puntos fijos (tipo C0-L1..L3) sobre T³, RGG3, retazos, caveman y expansores.
   - **P3, crecimiento contiguo.** Sigue como alternativa y su borrador ya incluye las seis preguntas del Consejo. Ataca R2 por construcción, porque el crecimiento contiguo no puede pegar de forma incoherente. Pero arrastra la cuestión de la causalidad, que el Consejo aplazó.

5. **Siguiente acción concreta:** sesión del Consejo para elegir entre Ω-D y P3, o rechazar ambos. Con la elección se prerregistra el análisis de puntos fijos o de la regla, sin dinámica. Solo después se simula con N ≤ 343 y la batería.
