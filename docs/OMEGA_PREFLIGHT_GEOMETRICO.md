# Ω — Preflight geométrico A–H (estándar aprobado por el Consejo tras la Fase 2)

Toda funcional o dinámica nueva debe pasar este preflight **antes** de cualquier campaña grande de simulación.

La pregunta no es «¿produce algo parecido a 3D?». La pregunta es: **¿contiene un mecanismo capaz de distinguir una variedad geométrica de una red local no geométrica?**

| Punto | Propiedad | Medida operativa (herramienta existente) | Estado |
|---|---|---|---|
| **A** | Localidad | `classify_c0` (H_null ≥ 1.15, ciclos cortos), salto medio frente a nulo de grados | Implementado |
| **B** | Homogeneización intermedia (Nivel I) | `omega/diagnostics/ball_growth.py`: ρ = CV(r_max)/CV(2) en la ventana N/4 | Implementado. P1-D PARCIAL (sens 1.0, spec 0.85); P1-D.2: 13/13 geometrías homogéneas pasan. Mide invariancia de traslación estadística; no es suficiente por sí solo |
| **C** | Isotropía | `omega.geometry.local_structure.isotropy` (ratio mediano) | Implementado (certificado Ω-1.1) |
| **D** | Consistencia multiescala (Nivel II) | D_B(r), D_eff, D_s y D_L por escalado con N (`c0_battery`) | D_B implementado. **La meseta por rango relativo (P1-D.2) no es fiable**: la engaña un cruce (WS β = 0.003) y una oscilación periódica la hace fallar con geometrías gruesas. Criterio de convergencia suavizado **pendiente de prerregistro (P1-D.3)** |
| **E** | Ser variedad | Certificado Ω-1.1, código F9 | Implementado |
| **F** | Estabilidad | Atractor desde inicios genéricos (L4, R8), persistencia en tiempo (D-1), KKT/hessiana | Implementado (arnés C0) |
| **G** | Robustez | Semillas, N (R6/F1 hasta N = 729), reetiquetado (R7), parámetros vecinos (C1-D, meseta), Θ razonable (C0-Θ) | Implementado salvo C1-D |
| **H** | Control no geométrico con las mismas estadísticas locales | `tools/c1_local_stats.py` (k, c, q por arista); **retazos** (`tools/p1_diagnostic.py::patchwork`), con geometría local idéntica y pegado incoherente | Implementado |

## Reglas

1. Hay que pasar **A–H en conjunto**. Ningún punto aislado, tampoco B, cuenta como evidencia de geometría (`RULE_D3_NEVER_SUFFICIENT` sigue vigente).
2. **H es obligatorio.** Si la funcional no pierde energéticamente o dinámicamente frente a un control con las mismas estadísticas locales (retazos, finales C0), no puede seleccionar geometría.
3. Los observadores de A–H son **diagnósticos**. No pueden convertirse en términos de la acción sin una decisión explícita del Consejo y un análisis previo que descarte la ingeniería *post hoc*.
4. **Orden:** análisis del paisaje (tipo C0-T1 / C1-A..D), preflight A–H con N ≤ 343, y solo después escalado a N = 729 o más con el certificado.

## Los tres niveles (Consejo, tras P1-D.2)

| Nivel | Pregunta | Herramienta | Estado |
|---|---|---|---|
| I | Coherencia multiescala: ¿los vecindarios se homogeneizan? | ρ (punto B) | Mide invariancia de traslación; validado solo con toros planos |
| II | Dimensión estable | D_B (punto D) | Definición actual poco fiable (§ arriba) |
| III | Geometría tipo variedad | Certificado Ω-1.1 (E, C) | Calibrado solo para la clase 3; rechaza RGG2 y anillo 1D |

I → II → III: ninguno sustituye al siguiente. Hay que informar siempre en dos categorías:
- **Geometría gruesa:** cuasi-isométrica a un espacio homogéneo (caveman ≈ 1D, retículo de cliques ≈ 3D).
- **Geometría tipo variedad:** vecindarios de variedad (lo que exige F9).

Además, **la dimensión no se fija**: ningún nivel exige D = 3. D_B ≈ 3 no es requisito.

## Estado del instrumento tras P1-D.3 (VÁLIDO, con alcance)

| Componente | Herramienta |
|---|---|
| Observador a N grande | `omega/diagnostics/sampled_growth.py` (BFS muestreado, D_B2, D_s) |
| Nivel II | Criterio de **convergencia** de D_B2, que sustituye a la meseta de P1-D.2 |

**Resultado:** sensibilidad 15/15 y especificidad 33/33 con N ≈ 2·10⁴. La separación crece con N.

**Alcance:** espacios cerrados homogéneos o casi homogéneos.
- **Falso negativo con borde:** el RGG en caja sale NO_GEOMÉTRICO.
- El Nivel II aislado puede dar D ≈ 3 convergente en un anillo con 0.1 % de atajos. Nunca debe usarse solo.

**Informe obligatorio:** NO_GEOMÉTRICO / GEOMETRÍA_GRUESA(d) / GEOMETRÍA_VARIEDAD(d ≥ 3) / VARIEDAD_NO_EVALUABLE (d ≤ 2).

## Diagnósticos complementarios (tras L-COH-0b y la auditoría (b))

| Diagnóstico | Herramienta | Atrapa | No ve | Estado |
|---|---|---|---|---|
| **B-bis: coherencia multiescala γ** | `omega/diagnostics/coherence.py` | Retazos, cruces lentos (WS β = 0.001, atajos al 0.1 %), cactus; **corrige el falso negativo con borde** | Árboles críticos (amenables), tubos 1D gruesos | Diagnóstico, PARCIAL (0.95 / 0.81) |
| **K: curvatura de Ollivier (escala 1)** | `omega/curvature/ollivier_sparse.py` | ER y RR, retazos, cliques, árboles (f_neg); corrige los bordes | WS y atajos escasos; C0 | Diagnóstico, APORTA (no redundante con ρ) |

**Reglas:**
- Ninguno es certificado, criterio de dimensión ni término de energía.
- La coherencia como objetivo colapsaría a 1D (CH-T5).
- Una regla de decisión combinada queda **pendiente de prerregistro** con familias nuevas.

## B-ter: un extremo a toda escala (tras L-CIC-0b)

| Herramienta | Atrapa | No ve | Estado |
|---|---|---|---|
| `omega/diagnostics/annulus.py` (conexidad del anillo r ≤ d ≤ 2r) | Árboles (incluido el uniforme), cactus, árboles de cliques; distingue 1D (DOS_EXTREMOS); cliques y expansores sin ventana | 2-árbol y estructuras de mundo pequeño trianguladas (una sola escala); retazos; tubos a escala menor que su fibra | Diagnóstico, PARCIAL (1.00 / 0.87). No selecciona dimensión |

**Regla combinada RC-1** (localidad + γ + K + Nivel II), validada fuera de muestra: **INVÁLIDA** (sensibilidad 0.696, especificidad 0.97). Sirve como **filtro de rechazo conservador**, no como certificado de aceptación.

## Propuestas de dimensionalidad (tras L-DIM-1 y L-ARQ-0)

Toda propuesta que pretenda **seleccionar** una dimensión debe cumplir, además de A–H, los criterios **E1–E6** de `docs/OMEGA_ARQ_0.md` §4.

- **E2, prueba del dial:** se barre cada parámetro continuo durante ≥ 4 décadas. Si aparecen ≥ 2 enteros, la propuesta es un dial y se rechaza, salvo que el parámetro esté fijado de antemano por un principio independiente.
- **E5, no identificación:** el nulo de etiquetas permutadas debe quedar en ≤ 5 %.
- **E6, no codificación:** auditoría escrita S0–S2, puerta estática E6-S v1.1 (SELECTOR o DIAL ⇒ FAIL conservador, revisable por transporte) y, tras la dinámica, E6-D (olvido del inicio, dial, orden, representación y tamaño). Véase `docs/OMEGA_E6_PRERREGISTRO.md` y `docs/OMEGA_E6_RESULTADOS.md`.
- El **problema A** (exclusión de degenerados) va antes que el **B** (selección).

## Sesgo dimensional de los diagnósticos (tras L-A-0)

**RC-2 INVÁLIDA** (sensibilidad 0.632, especificidad 0.864, panel fuera de muestra).

- **La ventana de κ de RC-1/RC-2 y cualquier exigencia de «≥ s escalas» a N fijo penalizan las d altas.** No deben usarse como jueces ciegos a d.
- Todo informe de B-ter, γ o K debe declarar su resolución dimensional d_max(N).
- Toda validación debe incluir d = 2, 3 y 4, con varias discretizaciones de cada una.

## Juez A− (RC-3, tras su validación)

`tools/rc3.py`: pares N / 8N del mismo proceso, con exclusiones X1–X4 (δ, anillos, coherencia) y sin estadísticas locales.

- **Resultado:** PARCIAL; FE = 0 con d = 2…6.
- **Uso obligatorio:**
  - informar R, δ, γ y f̃ en los dos tamaños;
  - una exclusión que dependa solo de X4 es «X4-marginal»;
  - todo NO-EXCLUIDO se inspecciona.
- **Puntos ciegos:** `docs/OMEGA_MAPA_EXCLUSIONES.md`.

## Regla de trazabilidad de d (tras el cierre de B)

Toda estructura que supere RC-3 + W6 debe declarar de qué fuente procede su d (parámetro o escala, inicio, conservación, ley, o «desconocida») antes de cualquier lectura dimensional. «Desconocida» obliga a E6-D completo. Un éxito del problema A no cuenta como progreso del problema B (`docs/OMEGA_CIERRE_B.md` §7).
