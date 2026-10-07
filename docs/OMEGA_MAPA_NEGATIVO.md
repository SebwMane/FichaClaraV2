# Ω — Mapa negativo consolidado (documento técnico interno)

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-coherencia`.
- **Mandato:** aprobado por el Consejo como (a). Recoge hipótesis, prerregistros, resultados, criterios, fallos, correcciones, lo aprendido y lo que queda abierto.
- **Verificación:** todas las cifras proceden de los documentos citados. La compilación inicial la hizo un agente Haiku y el cerebro la revisó contra las fuentes; las cifras del bloque L se comprobaron línea a línea en la rama congelada.

## 1. Redacción congelada del estado del programa

> Las familias de reglas locales que hemos definido y analizado hasta ahora **no han producido** un mecanismo que seleccione espontáneamente una geometría extendida de dimensión finita.

La versión fuerte («ninguna regla local sin un entero impuesto puede hacerlo») es una **conjetura de nivel 4**, no una conclusión.

## 2. Método común

- **Prerregistro antes del código:** criterios, semillas, N y predicciones.
- **Niveles de afirmación:**
  1. demostrado;
  2. explicación;
  3. teorema externo verificado;
  4. extrapolación.
- **Enmiendas:** fechadas respecto a los datos vistos.
- **Ramas:** una congelada por fase.
- **Prohibiciones:** coordenadas, D = 3 o k = 6 como objetivo, mecánica cuántica y cambiar umbrales del certificado después de ver resultados.
- **Juez geométrico:** batería P1-D.3 (Niveles I y II por bolas) y certificado Ω-1.1 (Nivel III, calibrado para la clase 3). Diagnósticos complementarios: coherencia γ (B-bis) y curvatura de Ollivier.

## 3. Tabla de hipótesis

| # | Hipótesis / mecanismo | Rama congelada | Documento | Veredicto | Cifras clave | Mecanismo de fallo |
|---|---|---|---|---|---|---|
| 1 | **S0** (triángulos y densidad, Ω-1.1) | `claude/omega-1.1-congelado` | `OMEGA_BLOQUE_L_RESULTADOS.md` | ✘ geometría | L-2: el uniforme gana en 180/180 puntos × 5/5 semillas. L-3a: T³/RGG3 nunca KKT en S0. L-3b: 0 persistentes en 1680 corridas; 0/1659 convergidos no KKT; J ≤ 0.43 | Vacío (apagado de amplitud) o soporte de cliques (lema L-3a) |
| 2 | **Ω-B** (variante de escala) | ídem | ídem | ✘ | L-2b: mejor no geométrica, clique 810/900. O4B-1: Ω-B termina uniforme 70 / clique 198 / multi 2 | Uniforme o cliques |
| 3 | **C0** (competencia más saturación) | `claude/omega-c0-congelado` | `OMEGA_C0_RESULTADOS.md` | ✔ localidad, ✘ geometría | L4: 253/675 locales, 0 vacíos. R6: 24 candidatas. F1: 10/10 NEGATIVO; con N = 729, 0/24 certificados; D_s 1.77–2.76 | Complejo de cliques no variedad; según el inicio: 1D (U), expansor (E), d ≈ 2–2.5 (R) |
| 4 | **C0 a Θ > 0** | `claude/omega-fase2-congelado` | `OMEGA_FASE2_RESULTADOS.md` | ✘ (FRÁGIL-1/N) | Θ_c/Θ_N = 0.1 con N = 216 y 343; f_bg 0.15–0.18 con 0.1 Θ_N; destruida con 1 Θ_N | La entropía de ~N pares de fondo gana; Θ_N ∝ 1/N |
| 5 | **C1** (4-ciclos saturados) | ídem | `OMEGA_FASE2_PRERREGISTRO.md` §1.3; `RESULTADOS` §2.1 | ✘ (analítico) | T³ q = 4, Q₆ q = 5. (k, c, q): C0 ≈ RGG3 (q 34–48 frente a 51) | La dimensión pasa a ser parámetro; no separa C0 del RGG3 |
| 6 | **Ω-D** (cercanía por difusión) | `claude/omega-d-congelado` | `OMEGA_D_L0_RESULTADOS.md` | ✘ premisa | K1: f_nd RGG3 0.44, retazos 0.31–0.40. K2: RGG3 s = −0.40. Retículos Q ≈ d/2 (1.59; 1.007) | Transitorio perezoso; la persistencia favorece 1D; el α objetivo impone d |
| 7 | **P3-A** (ancla uniforme) | `claude/omega-p3-congelado` | `OMEGA_P3_L0.md` | ✘ (analítico, T1) | Diámetro O(log N) | Árbol recursivo de expansión → mundo pequeño |
| 8 | **P3-B** (frontera por saturación) | ídem | ídem | ✘ (T2) | Capacidad κ_c − 2m | 1D o mundo pequeño |
| 9 | **P3-C** (capas con coalescencia) | ídem | ídem | ✘ (T3) | μ_eff > 1 / < 1 / = 1 | Hiperbólico, extinción o polímero ramificado; d dependiente de p |
| 10 | **P3-D** (planitud de Ollivier) | ídem | `OMEGA_P3_L0B_RESULTADOS.md` | ✘ premisa | Árbol: mediana −0.167 / 0.000. Retazos: f_neg 0.67–0.79. C0: +0.03 | K ≈ 0 no excluye árboles ni C0 |
| 12 | **R1-1 / SQ** (reescritura irreversible: borrar las aristas sin 4-ciclo y reenganchar los vértices de grado ≤ 1; solo se conserva N) | `claude/omega-r1-1-congelado` | `OMEGA_R1_1_RESULTADOS.md` | ✘ NEGATIVA (Ω6 0/6) | 606 pares: FRAGMENTADO 0.69–0.94; MULTIESTABLE; residuo I1 = grafo aleatorio de tamaño finito (k³/N) | Los dominios con cuadrados nuclean pero no coalescen |
| 13 | **R3-0** (estructuras de orden: espacio de cortes de un orden causal; coalescencia por unión) | `claude/omega-r3-0-congelado` | `OMEGA_R3_0.md` | ✘ REDUCIDO | J1(w): δ = 0.489 / 0.337 / 0.254; J2: X3 ×3; J3: residuo NO CLASIFICADO (δ −0.06…0.24; ln \|L\| ~ √n); J4(0.5): X3 | La coalescencia sin defectos fija d = anchura asintótica (Dilworth); ninguna construcción natural la fija sin codificarla |
| 14 | **W-1** (umbrales de paseo perezoso sobre el propio grafo como selector: retorno, colisión, cruce de 2 y 3 trayectorias) | `claude/omega-w1-congelado` | `OMEGA_W1_RESULTADOS.md` | ✘ W1-D (por insuficiencia en N₁) | N₂: L1 y L2 cortan en 2, L3 en 4, L4 sin corte; peine-2 (d_g = 3) CRECE en L1 y L2 | Medir no es generar; el corte depende de la ley; la recurrencia ve la isoperimetría, no la dimensión |
| 11 | **D5-B** (selección estática expansión + θ·fragilidad^p) | `claude/omega-ldim1-congelado` | `OMEGA_LDIM1_RESULTADOS.md` | ✘ SELECCIÓN-NO-VIABLE | p = 1/2, 1, 2: meseta d* = 2, L4 falla siempre (clique; árboles); nulo B 0.08–0.09; NC-2 \|ρ\| ≤ 0.47 | Compactificación (T2); baja expansión degenerada; d* = f(θ, N), y de p asintóticamente (T1) |

## 4. Instrumentación (resultados positivos del programa)

| Fase | Veredicto | Lo que fijó |
|---|---|---|
| P1-D | PARCIAL (21/21, 17/20) | ρ mide la invariancia de traslación a escala intermedia |
| P1-D.2 | SEPARA-PARCIAL | La meseta de D_B no es fiable; geometría gruesa ≠ variedad |
| **P1-D.3** | **VÁLIDO** (15/15, 33/33, N ≈ 2·10⁴–5·10⁴) | Batería obligatoria; dominio: espacios cerrados casi homogéneos; D ≈ 3 fabricable con WS |
| L-COH-0b | PARCIAL (0.95, 0.81) | γ corrige los bordes y los cruces lentos; no excluye árboles críticos |
| L-CIC-0b | PARCIAL (1.00, 0.87) | Un extremo a toda escala (conexidad de anillos): excluye árboles (incluido el uniforme), cliques y expansores sin imponer d ni densidad; falla con el 2-árbol (sin escalas) |
| RC-1 (fuera de muestra) | INVÁLIDA (sens 0.696, spec 0.97) | La combinación congelada es un rechazador conservador; pierde RGG3 k8, T³ y cajas; acepta el anillo 1D |
| Auditoría (b) | K APORTA | Curvatura no redundante; corrige ER y los bordes; ciega a los atajos |

## 5. Correcciones registradas (errores propios y ajenos)

| Origen | Error | Corrección |
|---|---|---|
| Cerebro, W-1 | Predijo W1-C, L3 no resuelta, L2 no universal por el peine-1 | W1-D por insuficiencia de Z⁵ en N₁ (T_max = 144, a 16 pasos del mínimo); L3 corta en 4; el peine-1 crece (finitud solo asintótica) y la no universalidad apareció en el peine-2 |
| Cerebro, W6 | Predijo τ_δ ≈ 0.03–0.06 | τ_δ = 0.69: el tope de radio de RC-3 en 2D (W4) entró en la calibración; W6.2 no es informativa |
| Cerebro, R3-0b | Predijo J1(4) X4-marginal, J3 explosivo visible (X1 o abstención) y J4(0.1) X3 | J1(4) no excluido; J3 residuo con δ ≈ 0.1–0.24 (crecimiento intermedio, W6); J4(0.1) X1 + X4 (régimen transitorio) |
| Cerebro, inspección J3 | Umbral de descenso de δ de 0.02 con 3 semillas, sin cálculo de potencia; media de δ en [0.10, 0.25] | Descenso medio −0.005; media 0.087: NO CLASIFICADO. Prueba sin potencia (lección 25) |
| Cerebro, E6-S v1 | Predijo VÁLIDO y que P0-LAT(3) casi no descansaría; P0-DEG(3) predicho como SELECTOR genuino | El panal hizo SELECTOR a P0-LAT(3): v1 INVÁLIDO. Corregido en v1.1 con panel cerrado; «grado = 3» no codifica d |
| Cerebro, E6-S v1.1 | Predijo P0-DEG(3) MONÓTONO y las demás clases de v1 sin cambios | P0-DEG(3) CIEGO (escalera cúbica en d = 1); cambian 3 reglas fuera de V1′–V3′ |
| Cerebro, C0-L4 | Predijo MUERTE; salió CONTINÚA | Fase local real, pero no geométrica |
| Cerebro, P1-D | Predijo que los retazos homogeneizarían | No homogeneizan |
| Cerebro, P1-D.2 | Predijo que el retículo de cliques pasaría el Nivel II | La oscilación de período 2 lo impidió |
| Cerebro, L-ΩD-0 | Predijo K1 ✔ | Falló |
| Cerebro, L-P3-0b | Predijo KD1 ✔ y RGG con κ positivo | El árbol tiene signo mixto; el RGG3 es plano |
| Cerebro, L-COH-0b | Predijo que los árboles aleatorios serían incoherentes | Los árboles críticos son amenables. **CH-T2 restringido a ramificación acotada inferiormente** |
| Cerebro, L-COH-0b | C × RR diseñado con un expansor demasiado pequeño | Resultó un tubo 1D (cuasi-isométrico a un ciclo) |
| Cerebro, L-CIC-0b | Predijo el 2-árbol SIN_VENTANA, el tubo DOS_EXTREMOS y RC-1 PARCIAL o VÁLIDA | 2-árbol CONEXO (una sola escala); tubo CONEXO (r < diámetro de la fibra); RC-1 INVÁLIDA |
| Cerebro, L-DIM-1 | Predijo d*(p) → 4, 3, 2, NC-5 confirmada y nulo A sustancial | Meseta d = 2 para los tres p a N ≈ 2·10⁴; NC-5 no confirmada; nulo A 0.008 (alto solo en L1: 0.37) |
| Cerebro, prerregistro L-DIM-1 | «Σ 1/\|S_r\| < ∞ ⇔ transitoriedad» | Solo en retículos; en general Nash-Williams da solo recurrencia (Thomassen 1992) |
| Cerebro, L-A-0 | Predijo RC-2 PARCIAL y la aceptación de RGG4 k16 y de T³ con diagonales | INVÁLIDA; ambas rechazadas por los sesgos dimensionales de ventana y curvatura |
| Cerebro, L-CIC-0b | «B-ter es independiente de la dimensión» | Solo se comprobó con d = 2 y 3; a N fijo, el número de escalas decrece con d |
| Cerebro, RC-3 | Predijo que RGG3 + 0.1 % de atajos no se excluiría y que Heisenberg y los retazos sí | Atajos excluidos (X4-marginal); Heisenberg y retazos no excluidos; árbol de cubos no excluido (no previsto) |
| Agente Sonnet, RC-3 | Veredicto «VÁLIDA» por lectura literal de la regla por familia | El cerebro adopta la lectura conservadora: PARCIAL |
| Cerebro, R1-1 | Atribuyó el fallo previsto a la nucleación (∝ 1/N) | Hay nucleación en RAND/ASYNC; el fallo es de coalescencia |
| Cerebro, R1-1 | Predijo X1 para el inicio denso I1 | NO-EXCLUIDO por el punto ciego W5 de RC-3 (cambio de régimen entre tamaños); la clase es grafo aleatorio |
| Consejo | «C + K ⇏ geometría» | No demostrado: solo cada uno por separado; la combinación estaba sin validar |
| Consejo | «Ya reconocemos la geometría» | RC-1 fuera de muestra: rechazo fiable (0.97), aceptación con pérdidas (0.70) |
| Auditor Haiku | «24/24 pasan el certificado con N = 729» | Falso: 0/24 (verificado a mano) |
| Auditor Haiku | Etiquetó L-1/L-2/L-3a como «éxito» | Son análisis que demuestran destinos degenerados |
| Consejo | «C0 tiende a vacío/cliques» | C0 rompió la dicotomía; su fallo fue no ser geometría |
| Consejo | «Terminar F1», «D-1 pendiente», «Θ pendiente» | Ya estaban hechos |
| Consejo | «D_eff ≈ 3 y D_s ≈ 3» como objetivo | Contradice la decisión de clase abierta |
| Consejo | Tabla histórica: «C1 — controles y atajos» | C1 fue la funcional de 4-ciclos |

## 6. Lo aprendido (restricciones demostradas o verificadas)

1. **Acciones locales en (k, c):** no pueden preferir energéticamente la geometría (C0-T1/T2, nivel 1).
2. **Localidad ⇏ geometría** (C0).
3. **Conteos hasta segundo orden ⇏ geometría:** C0 ≈ RGG3 en (k, c, q).
4. **D ≈ 3 ⇏ espacio 3D** (WS β = 0.001; R-3D de C0 con N = 343).
5. **Difusión sola ⇏ cercanía geométrica**, en τ cortos (L-ΩD-0).
6. **Crecimiento local genérico** → mundo pequeño, 1D, hiperbólico o polímero ramificado (T1–T3).
7. **La planitud es un punto crítico** (Gauss–Bonnet, T5, nivel 1); K ≈ 0 ⇏ geometría.
8. **Coherencia de Følner ⇏ geometría:** los árboles críticos y las cadenas la cumplen. Y homogeneidad más coherencia ⇒ dimensión entera, pero no necesariamente euclídea (Bass–Guivarc'h, Trofimov; CH-T6).
9. **Valores objetivo de solapamiento, de incompatibilidad o de 4-ciclos imponen la dimensión** (L-ΩD-T2, CH-T1, Fase 2 §1.3).
11. **Clase nueva: «topología no degenerada sin geometría macroscópica»** (2-árbol aleatorio). Tiene ciclos abundantes, redundancia y estructura local rica, pero diámetro logarítmico y ninguna jerarquía de escalas extensas.
12. **Control de la expansión necesario pero no suficiente** (L-DIM-0). Hay tres presiones (F-exp, F-ram, F-baja), y las de F-exp y F-baja se oponen.
13. **Competencia estática ⇏ selección de dimensión** (L-DIM-1). Los extremos degenerados (clique, árboles) ganan cualquier funcional expansión–fragilidad. El ganador entre geometrías lo fijan θ y N, y p asintóticamente (L-DIM-T1). Dentro de funcionales de bolas con leyes de potencias, la invariancia de escala y la no codificación por exponente son incompatibles (nivel 2).
14. **Lema del dial** (L-ARQ-T1, nivel 1): toda competencia escalar A + θB entre candidatos selecciona solo vértices de la envolvente convexa inferior de {(B_d, A_d)}, y el entero lo fija θ. Una meseta mide una arista de la envolvente, no una selección (verificado 75/75 con los datos de L-DIM-1).
15. **Principio de no identificación:** una propiedad estable que sobrevive a una asignación aleatoria de la etiqueta dimensional (nulo B de L-DIM-1: 8–9 %) no es evidencia de selección.
16. **Paralelo externo (nivel 3):** las triangulaciones dinámicas euclídeas presentan las mismas degeneraciones (fase arrugada ≈ clique; polímero ramificado ≈ árbol). La extensión solo aparece con estructura impuesta (foliación causal) y con la dimensión fijada por el bloque.
17. **El juez de no degeneración a N fijo no es ciego a d** (L-A-0, RC-2 INVÁLIDA: sensibilidad 0.632, especificidad 0.864). El número de escalas disponibles decrece con d (r_w ≈ (N/4c)^{1/d}), y la curvatura de Ollivier a escala 1 decrece con d y depende de la discretización (RGG2 +0.1, RGG3 −0.03, RGG4 −0.11; T³ con diagonales −0.10). RC-1 y RC-2 tienen un sesgo hacia d baja. Todo juez futuro debe declarar su resolución d_max(N) y validarse con d = 2..4.
18. **Instrumento ≠ mecanismo** (Consejo). Un instrumento detecta G; un mecanismo explica Ω → G.
19. **Un juez de exclusión de escalado (N → 8N) puede ser ciego a d en el rango validado** (RC-3: cero exclusiones falsas con d = 2…6). Lo que falla con N fijo son las estadísticas locales y de ventana corta, no la exclusión como tal. Límite: ningún juez a N finito es ciego a d para toda d (O3), y los árboles de bloques grandes solo se excluyen a N ≫ bloque^d (W2).
20. **Las reglas de reescritura reversibles son ensembles de equilibrio (clase cerrada).** R1 solo aporta algo si rompe el balance detallado (R1-0, Ω1).
21. **Con información de radio 1 (vecinos comunes t), geometría y árbol son indistinguibles localmente.** Ninguna regla de umbral sobre t con todas las Z^d como puntos fijos borra la memoria de un inicio arbóreo (R1-T1, censo exacto). Toda regla que apunte a la fracción de cuñas cuadradas ((d−1)/d en Z^d) codifica d (R1-T2).
22. **Nuclear no es crecer** (R1-1). Una regla de estabilidad ciega a d (presencia de 4-ciclos) forma núcleos ricos en cuadrados que se absorben aislados. El cuello de botella del problema A es la **coalescencia** de dominios, no su formación.
23. **«Codificar d» es relativo a un panel de referencia** (E6). Un motivo local presente en una sola dimensión del panel hace parecer selectora a cualquier regla que descanse en él; E6-S v1 falló por eso (panal: grado 3 sin ciclos cortos, solo en d = 2). Con el panel cerrado por truncación, grafo de líneas y producto con K₂, «grado = 3» resulta CIEGO. CIEGO es estable al ampliar el panel; SELECTOR no lo es. Por eso un FAIL de E6-S es conservador y revisable, y un PASS es firme dentro de su alcance.
24. **En la coalescencia sin defectos, d es un recuento** (R3-0). El espacio de cortes de un orden causal vive en ℤ^w, con w la anchura (Birkhoff–Dilworth). Su dimensión es el número de hilos que avanzan independientemente sin límite: hilos libres dan w (codificado), hilos acoplados dan 1D, y una anchura creciente da crecimiento intermedio o explosivo. En grupos, redes rígidas y órdenes, d es el rango de direcciones que conmutan.
25. **Un umbral sin cálculo de potencia no es una prueba** (R3-0b). El descenso de δ predicho (≈ 0.03) era menor que el ruido entre pares (≈ 0.05–0.08). La inspección quedó NO CLASIFICADO por falta de potencia, no por evidencia contraria. Toda regla de decisión futura con umbral de efecto lleva su cálculo de potencia en el prerregistro.
26. **Los umbrales estructurales trasladan la selección a la elección de ley** (W-0). Las propiedades de paseos sobre el propio grafo tienen enteros críticos sin parámetro continuo, pero el entero depende de qué ley se elija (recurrencia, colisión, intersección, aridad, red o continuo). Por la regla DIAL de E6, esa elección es un dial discreto, salvo un principio independiente que no se puede demostrar dentro de Ω.
27. **La recurrencia no es dimensión** (W-1). Un peine de dimensión de crecimiento 3 es recurrente y las colisiones crecen en él como en d ≤ 2. Las leyes de paseo miden la conductancia o la isoperimetría. Además, sus umbrales en d alto solo son accesibles cerca de 10⁶ nodos, y una sola celda insuficiente puede decidir una prueba: el tamaño mínimo por celda debe comprobarse con margen en el prerregistro.
28. **Un atractor concreto solo desplaza la pregunta** (cierre B). «w* = n porque la regla tiene su atractor en n» elige n al elegir la regla. Solo la universalidad de w* sobre una clase de reglas propuesta por un principio cuenta como emergencia. Con información local, regular el número de hilos cae en dial, inicio o ley, o colapsa a una sola dirección (dilema de independencia, W-T7).
29. **El espacio de cortes es un espacio de configuraciones** (CF-0). L(J₁ ⊔ J₂) ≅ L(J₁) × L(J₂): unir historias independientes suma dimensiones. La geometría de R3 cuenta procesos, no lugares; por eso B se redujo a un recuento con tanta limpieza. Toda propuesta debe declarar si mide el espacio de configuraciones o el propio orden causal.
30. **La conmutación por pares no basta para rango ≥ 3** (CF-0, `results/cf0/cd_check.json`). Completar cada cuña con un elemento nuevo no cierra para w ≥ 3, e identificar elementos de forma ingenua colapsa la estructura. Un producto coherente necesita una condición de ternas (cubo). La aridad de la ley de coherencia acota el rango alcanzable (W-T4). Advertencia: esto no selecciona ningún d.
10. **Un extremo a toda escala separa la dimensión ≥ 2 de árboles, 1D, cliques y expansores, pero no 2 de 3** (CIC-T5). Las propiedades estructurales estudiadas son todas ciegas a la dimensión por construcción: **la selección de dimensión es la frontera abierta**.

## 7. Lo que queda abierto

- **Cierre oficial de D5 (frase del Consejo):** «L-DIM-1 no demuestra que Ω sea incapaz de generar una dimensión. Demuestra que la dimensionalidad no puede darse por emergente simplemente porque un funcional estático produzca una meseta o un mínimo interior. En la familia ensayada, la selección queda controlada por parámetros del funcional, escala finita y estructuras degeneradas. Por tanto, cualquier mecanismo futuro de dimensionalidad deberá explicar simultáneamente por qué existe una dimensión, por qué no es un parámetro oculto y por qué las estructuras degeneradas quedan excluidas sin introducir la dimensión que se pretende explicar.»
- **L-ARQ-0 (síntesis):** la jerarquía conectividad → localidad → coherencia → un extremo → ? → dimensión → geometría. Ω solo tiene mecanismos en los dos primeros niveles. El candidato no dial para «?» es el rango de un sistema de relaciones conmutativas (Bass–Guivarc'h). Criterios de entrada E1–E6 (`OMEGA_ARQ_0.md` §4). Decisión pendiente del Consejo: L-ARQ-1 (conmutatividad con rango no prescrito) o volver al problema A.
- **E6 (no codificación):** v1 INVÁLIDA; v1.1 aprobada con reserva (validación débil, etiqueta obligatoria). Deuda: panel de grafos reservado y ciego antes del primer candidato de tipo grafo.
- **CF-0 (confluencia, propuesta al Consejo):** análisis sin dinámica. Hallazgos: L-ARQ-1 nunca se ejecutó (A1); SQ no era confluente (A2); R3 midió un espacio de configuraciones (A3). CF-T1 a CF-T6: suma de rangos, condición del cubo, dilema de la extensión (inflación, codificación Ω3 o CONSERVADOR), conmutación parcial (Z^n o exponencial; Heisenberg), monogénesis y escalera de aridad. Se pide autorizar CF-1 **solo para el problema A**.
- **PROBLEMA B CERRADO** (`OMEGA_CIERRE_B.md`, 2026-10-07). En las familias examinadas, d se reduce al número de direcciones independientes que conmutan, cuya selección exige parámetro, inicio, conservación o ley elegida. No es un teorema de imposibilidad. Criterios de reapertura congelados (§6): principio independiente que fije la ley, clase de reglas propuesta por un principio, o mecanismo no local justificado. El problema A sigue con sus instrumentos validados.
- **W-1:** W1-D literal (insuficiencia en N₁). En N₂, el corte depende de la ley (2 frente a 4) y no es universal (peine-2). Propuesta del cerebro: **cerrar el programa de selección (problema B)** con el resultado negativo estructurado (`OMEGA_W1_RESULTADOS.md` §5); no repetir W-1 con más tamaño (valor de la información nulo para esa decisión). Pendiente del Consejo.
- **W6:** VÁLIDA (regla de uso de RC-3: tres tamaños, misma estructura, grado saturado ≤ 5 %). La estabilidad de δ no decide.
- **W-0 (autoselección de anchura):** C1 y C2 cerradas (Poisson 2β/γ, lema del dial, conservación); C3 no vacía en principio (umbrales de paseos: recurrencia d ≤ 2, trayectorias d ≤ 4), pero reducida a una elección de ley (dial discreto) y degradada a C2 a tamaño finito por localidad. Aprobada por el Consejo como reducción, no como imposibilidad. W-1 autorizado (estático; prerregistro `OMEGA_W1_PRERREGISTRO.md`).
- **R3-0 (órdenes):** REDUCIDO. d = anchura asintótica en la coalescencia sin defectos; R3-T2(c) falsado en forma instrumental (residuo J3, punto ciego W6). Pendiente del Consejo: ratificar, decidir W6 y elegir entre W-0 (autoselección de anchura) y el cierre del programa de selección.
- **R1-1 (SQ):** NEGATIVA. Pendiente del Consejo: cierre de R1 con alcance limitado a SQ y similares, regla W5 y autorización de R3-0 (orden que haga crecer y unir dominios sin codificar d).
- **R1-0 (análisis de reglas):** se conserva solo N; la densidad y los ciclos son salidas. Familia mínima admisible SQ: borrar las aristas sin 4-ciclo y reenganchar los vértices de grado ≤ 1. Predicción congelada: degenerada (bosque con rotación). Siguiente paso: prerregistro de R1-1.
- **RC-3 (juez A−):** PARCIAL (lectura conservadora), FE = 0 con d = 2…6, EP 0.944 en las familias nuevas. Se adopta como juez provisional con puntos ciegos declarados (`OMEGA_MAPA_EXCLUSIONES.md`). Siguiente paso en el orden lógico: dinámica A (pendiente del Consejo).
- **L-A-0 (problema A):** RC-2 INVÁLIDA. El juez de A no es ciego a d. Pendiente del Consejo: A1 (juez relativo, RC-3), A2 (solo rechazo) o A3 (evaluación por tendencias con N). Ninguna dinámica antes.
- **L-DIM-1 (cerrada):** D5-B refutado como mecanismo estático. L-DIM-2 no procede. Queda una pregunta: ¿existe una razón independiente para un exponente (por ejemplo, la transitoriedad)?
- **L-DIM-0 (sesión conceptual):** la integralidad de d puede venir de la homogeneidad (Gromov, Trofimov, Bass–Guivarc'h). Un equilibrio continuo entre la presión de expansión y la de robustez sobre los enteros produciría **mesetas** (criterio C1-D). Es un esquema de nivel 4, pendiente de ratificar. «¿Por qué 3?» no es determinable desde Ω actual.

- **Pregunta central del Consejo:** ¿qué ingrediente mínimo falta para que una dinámica sin geometría inicial estabilice una estructura geométrica de dimensión finita?
- **Candidatos que el programa no ha discriminado todavía:** conservación, competencia entre términos no locales, restricciones topológicas (riqueza de ciclos en todas las escalas), defectos, fluctuaciones, reglas de reescritura y relaciones de orden superior (T_ijk). La hipótesis «coherencia + riqueza de ciclos» (§5 de los resultados de L-COH-0b) está informada por los datos y requeriría el tratamiento de D′.
- **Instrumento:** una regla de decisión combinada (localidad + ρ + γ + K + Nivel II), pendiente de prerregistro y validación con familias nuevas.
