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

## 4. Instrumentación (resultados positivos del programa)

| Fase | Veredicto | Lo que fijó |
|---|---|---|
| P1-D | PARCIAL (21/21, 17/20) | ρ mide la invariancia de traslación a escala intermedia |
| P1-D.2 | SEPARA-PARCIAL | La meseta de D_B no es fiable; geometría gruesa ≠ variedad |
| **P1-D.3** | **VÁLIDO** (15/15, 33/33, N ≈ 2·10⁴–5·10⁴) | Batería obligatoria; dominio: espacios cerrados casi homogéneos; D ≈ 3 fabricable con WS |
| L-COH-0b | PARCIAL (0.95, 0.81) | γ corrige los bordes y los cruces lentos; no excluye árboles críticos |
| Auditoría (b) | K APORTA | Curvatura no redundante; corrige ER y los bordes; ciega a los atajos |

## 5. Correcciones registradas (errores propios y ajenos)

| Origen | Error | Corrección |
|---|---|---|
| Cerebro, C0-L4 | Predijo MUERTE; salió CONTINÚA | Fase local real, pero no geométrica |
| Cerebro, P1-D | Predijo que los retazos homogeneizarían | No homogeneizan |
| Cerebro, P1-D.2 | Predijo que el retículo de cliques pasaría el Nivel II | La oscilación de período 2 lo impidió |
| Cerebro, L-ΩD-0 | Predijo K1 ✔ | Falló |
| Cerebro, L-P3-0b | Predijo KD1 ✔ y RGG con κ positivo | El árbol tiene signo mixto; el RGG3 es plano |
| Cerebro, L-COH-0b | Predijo que los árboles aleatorios serían incoherentes | Los árboles críticos son amenables. **CH-T2 restringido a ramificación acotada inferiormente** |
| Cerebro, L-COH-0b | C × RR diseñado con un expansor demasiado pequeño | Resultó un tubo 1D (cuasi-isométrico a un ciclo) |
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

## 7. Lo que queda abierto

- **Pregunta central del Consejo:** ¿qué ingrediente mínimo falta para que una dinámica sin geometría inicial estabilice una estructura geométrica de dimensión finita?
- **Candidatos que el programa no ha discriminado todavía:** conservación, competencia entre términos no locales, restricciones topológicas (riqueza de ciclos en todas las escalas), defectos, fluctuaciones, reglas de reescritura y relaciones de orden superior (T_ijk). La hipótesis «coherencia + riqueza de ciclos» (§5 de los resultados de L-COH-0b) está informada por los datos y requeriría el tratamiento de D′.
- **Instrumento:** una regla de decisión combinada (localidad + ρ + γ + K + Nivel II), pendiente de prerregistro y validación con familias nuevas.
