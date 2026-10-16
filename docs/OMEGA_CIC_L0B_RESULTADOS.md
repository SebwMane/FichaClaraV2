# Ω — L-CICLOS-0b y validación fuera de muestra de RC-1: resultados

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-ciclos`.
- **Prerregistro:** `docs/OMEGA_CIC_L0.md`.
- **Datos:** `results/cic_l0b/` (58 grafos con N ≈ 2·10⁴; 142 s).
- **Reparto:**
  - código de un agente Sonnet, revisado por el cerebro antes de correr (anillo = {r ≤ d ≤ 2r}, componentes en el subgrafo inducido, ventana y escalas según §3);
  - referencias de un agente Haiku, revisadas por el cerebro: Freudenthal y Hopf (0, 1, 2, ∞ extremos), Stallings, Z^d con un extremo para d ≥ 2, grupos de superficie con un extremo. El agente **no aportó citas exactas**: verificación débil, contenido confirmado por el cerebro con la literatura estándar (Freudenthal, *Math. Z.* 1931; Hopf, *Comment. Math. Helv.* 1944; Stallings, *Ann. Math.* 1968).

## 1. L-CIC-0b (redundancia cíclica = un extremo a toda escala): **REDUNDANCIA-PARCIAL**

| Medida | Resultado |
|---|---|
| Sensibilidad (V → CONEXO) | **23/23 = 1.00**: RGG2, RGG3 con k = 8, 12 y 20, esfera, cajas, cuadrado y T³ |
| Exclusión (X → no CONEXO) | **20/23 = 0.87**. Única familia que falla: el **2-árbol aleatorio** (3/3 CONEXO) |
| Independencia de la densidad | ✔ (k = 8, 12 y 20, todas CONEXO) |
| Independencia de la dimensión | ✔ (d = 2 y d = 3, CONEXO) |
| Familias nuevas / repetidas | Nuevas: exclusión 0.75 (el 2-árbol). Repetidas: VÁLIDA (1.00 / 1.00) |

| Familia | Estado | Predicho |
|---|---|---|
| **Árbol uniforme** (Prüfer) | **RAMIFICADO** (f̃ = 0.46–0.49) | ✔ **Es el ángulo ciego de la coherencia, y aquí queda cubierto** |
| Árbol binario subdividido, cactus, árbol de cliques K6 | RAMIFICADO | ✔ |
| Cliques K12 enlazadas | RAMIFICADO | Predije SIN_VENTANA; igualmente excluido |
| Grafo aleatorio con triángulos, ER, RR | SIN_VENTANA | ✔ (excluidos por ausencia de escalas) |
| **2-árbol aleatorio** | **CONEXO** (f̃ = 0.998) | ✘ (predije SIN_VENTANA) |
| Anillo k12, caveman | DOS_EXTREMOS | ✔ |
| WS β = 0.001 | RAMIFICADO (no CONEXO) | ✔ |
| **Tubo C₁₀₀ × RR₂₀₀** | **CONEXO** | ✘ (predije DOS_EXTREMOS) |
| Heisenberg | CONEXO | ✔ |
| Retazos 2³ | CONEXO | ✔ (no es objetivo de esta observable; CIC-T6) |

**Por qué falla con el 2-árbol (nivel 2).** El 2-árbol aleatorio es una triangulación arbórea con diámetro ~log N: su ventana solo admite **una escala** (r_w = 5, de modo que r = 2). A esa escala la triangulación rellena el anillo. La prueba «multiescala» degenera en una prueba local, que es precisamente el defecto CIC-T3 que pretendía evitar. Las estructuras de mundo pequeño nunca dan ventana suficiente, crezca o no N.

**Por qué falla con el tubo (nivel 2).** Con r ≤ r_w/2 ≈ 8, del orden del diámetro de la fibra (≈ 5), el tubo todavía no muestra sus dos extremos; haría falta r ≫ el diámetro de la fibra.

## 2. RC-1 (regla combinada congelada, fuera de muestra): **RC1-INVÁLIDA**

| Medida | Resultado |
|---|---|
| Sensibilidad (V → GEOMÉTRICO-RC1) | **16/23 = 0.696**, por debajo del umbral 0.7 |
| Especificidad | **33/34 = 0.97** |
| Familias nuevas | sensibilidad 0.50, especificidad 1.00 |
| Familias repetidas | sensibilidad 0.76, especificidad 0.95 |

**Qué condición falla en cada geometría rechazada:**

| Familia | Condición que falla |
|---|---|
| RGG3 k8 (3/3) | (ii) coherencia INTERMEDIO |
| T³ | (ii) coherencia INTERMEDIO (γ = 0.69) |
| RGG3 en caja (3/3) | (iv) Nivel II no converge: el borde curva D_B2, ya visto en P1-D.3 |

**Único falso positivo:** el **anillo k12**, que pasa las cuatro condiciones. RC-1 no tiene ninguna condición que excluya 1D.

**Lectura.** RC-1 es un **rechazador conservador**: casi nunca acepta una no geometría (0.97), pero pierde ~30 % de las geometrías auténticas. La frase del Consejo «ya reconocemos la geometría» queda matizada: **rechazamos la no geometría de forma fiable; la aceptación tiene pérdidas.**

## 3. Contraste con las predicciones

| Predicción | Resultado |
|---|---|
| L-CIC-0b VÁLIDA | ✘ PARCIAL (2-árbol) |
| Árbol uniforme RAMIFICADO (confianza moderada) | ✔ |
| Tubo DOS_EXTREMOS | ✘ |
| RC-1: especificidad ≈ 1.0 | ✔ 0.97 |
| RC-1: sensibilidad 0.8–0.95 | ✘ 0.696 |
| RC-1: veredicto PARCIAL o VÁLIDA | ✘ INVÁLIDA |

## 4. Observaciones *post hoc* (solo hipótesis; no se reajusta nada)

- Los tres casos que RC-1 rechaza indebidamente (RGG3 k8, T³, caja) son **CONEXO** en la observable de anillos.
- Su único falso positivo (anillo k12) es **DOS_EXTREMOS**.
- El fallo de los anillos (2-árbol) lo atrapan la coherencia y el Nivel II.

Una regla «anillos + coherencia + K», sin los umbrales estrictos de γ y del Nivel II, separaría este panel. Pero se construiría mirando estos datos, así que **no se adopta**. Requeriría el prerregistro de una RC-2 y un panel nuevo.

## 5. Dictamen del cerebro

1. **Redundancia cíclica multiescala (un extremo a toda escala):**
   - Es el **primer observable** que excluye árboles (incluido el uniforme), cliques y expansores **sin imponer dimensión ni densidad**.
   - Responde **sí, parcialmente**, a la pregunta de L-CICLOS-0.
   - Su límite estructural son las estructuras de mundo pequeño trianguladas, que no tienen escalas suficientes.
   - Se incorpora como diagnóstico **B-ter**. **No** selecciona dimensión (CIC-T5) y **no** es mecanismo.
2. **Hipótesis A frente a B del Consejo:**
   - El resultado es compatible con ambas. La conexidad de anillos es una forma de **redundancia relacional**: hay caminos alternativos que rodean cualquier región.
   - No distingue si los ciclos son causa o síntoma.
3. **El instrumento:** su combinación congelada es conservadora. Se puede usar como **filtro de rechazo**, no como certificado de aceptación.
4. **Lo que este bloque deja claro para el problema central:**
   - Las propiedades estructurales estudiadas definen con bastante precisión la clase **«geometría gruesa de dimensión ≥ 2»**: localidad, coherencia (amenable o polinómica), un extremo a toda escala y curvatura no degenerada.
   - **Todas son, por diseño, independientes de la dimensión.** Para estructuras homogéneas, coherencia más un extremo implican dimensión entera ≥ 2 (CH-T6 y Stallings), **sin decir cuál**.
   - La **selección de dimensión** sigue sin tocar. Es la pregunta que queda.

## 6. Propuesta al Consejo

- **(a) Cerrar la instrumentación aquí.** El mapa de diagnósticos está completo para las clases degeneradas conocidas. Una RC-2 solo si una dinámica futura lo exige.
- **(b) Abrir una sesión conceptual**, sin código, sobre **la selección de dimensión**: ¿qué principio podría preferir un entero concreto (2, 3, …) dentro de la clase «geometría gruesa de dimensión ≥ 2», sin imponerlo? Todos los ingredientes probados hasta ahora son ciegos a la dimensión por construcción, y eso no es casual: es la frontera real del problema.
