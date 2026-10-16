# Ω — L-COH-0b y auditoría (b): resultados

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-coherencia`.
- **Prerregistros:** `docs/OMEGA_COH_L0.md` §3 y `docs/OMEGA_K_AUDIT_PRERREGISTRO.md`.
- **Datos:**
  - `results/coh_l0b/`: 62 grafos con N ≈ 2·10⁴.
  - `results/k_audit/`: los 77 grafos de P1-D.3, regenerados idénticos (0 discrepancias).
- **Reparto:**
  - código de dos agentes Sonnet, revisado por el cerebro antes de las corridas: la intersección y unión de bolas, la pendiente γ, y la curvatura dispersa, verificada contra la densa hasta 1e-9;
  - referencias y compilación de cifras por dos agentes Haiku, **revisadas por el cerebro**: las cifras del bloque L se comprobaron contra el documento congelado (180/180, 810/900, 0/1659 y J ≤ 0.43 coinciden); las etiquetas de «éxito» que el agente puso a L-1/L-2/L-3a se corrigieron, porque son análisis que muestran que S0/Ω-B favorecen cliques o el uniforme.

## 1. L-COH-0b: veredicto oficial COHERENCIA-PARCIAL

| Medida | Resultado |
|---|---|
| Sensibilidad (V → COHERENTE) | **19/20 = 0.95**. Única falla: T³ 27³, INTERMEDIO con γ = 0.691, justo bajo el umbral 0.7 |
| Especificidad (F → no COHERENTE) | **26/32 = 0.81**. Falsos positivos: **árbol de Prüfer (3/3)** y **C₁₀₀ × RR₂₀₀ (3/3)** |
| Cajas (falso negativo de P1-D.3) | **RGG3 en caja y RGG2 en caja: COHERENTE 6/6**. La coherencia corrige el falso negativo del Nivel I |

**Por familia (mediana de γ):**

| Grupo | Familia | Mediana de γ | Estado |
|---|---|---|---|
| V | RGG2 | 0.85–0.87 | COHERENTE |
| V | RGG3 | 0.72–0.73 | COHERENTE |
| V | Esfera S² | 0.90–0.91 | COHERENTE |
| V | Cajas | 0.84–0.94 | COHERENTE |
| V | Gradiente | 0.74–0.75 | COHERENTE |
| V | Cuadrado | 0.87 | COHERENTE |
| V | T³ | 0.69 | INTERMEDIO |
| F | Retazos 2³ y 4³ | 0.13–0.22 | INCOHERENTE |
| F | Retazos 8³, ER, RR | — | SIN_VENTANA (crecimiento demasiado rápido) |
| F | WS β = 0.01 / 0.005 | 0.15–0.25 | INCOHERENTE |
| F | Atajos al 1 % | 0.22 | INCOHERENTE |
| F | Cactus de triángulos | 0.14 | INCOHERENTE |
| F | Árbol binario subdividido | 0.43 | INCOHERENTE |
| F | **Árbol de Prüfer** | 0.81–0.91 | **COHERENTE** |
| F | **C₁₀₀ × RR₂₀₀** | 1.38 | **COHERENTE** |
| E | Anillo | 0.98 | COHERENTE |
| E | Caveman | 1.01 (f_deg 0.39) | INTERMEDIO |
| E | Retículo de cliques 3D | 0.71 | INTERMEDIO |
| E | Heisenberg | 0.49 | INCOHERENTE |
| E | **WS β = 0.001 y atajos al 0.1 %** | 0.50–0.52 | **INCOHERENTE** |

## 2. Contraste con las predicciones congeladas

| Predicción | Resultado | Explicación |
|---|---|---|
| V coherente, incluidas las cajas | ✔ (T³ en el límite) | El borde no rompe la compatibilidad entre vecinos |
| **Árboles aleatorios incoherentes** | ✘ **(error analítico mío)** | CH-T2 vale para árboles **con ramificación acotada inferiormente** (árbol binario: ✔ incoherente). El árbol uniforme es crítico, un polímero ramificado con \|B_r\| ~ r²: tiene **crecimiento polinómico y es amenable**, así que la ley de Følner lo da por coherente. **La coherencia no excluye el destino árbol/fractal** |
| **C × RR incoherente** | ✘ **(error de diseño mío)** | Con un factor expansor de solo 200 nodos (diámetro ≈ 5), el producto es **cuasi-isométrico al ciclo**: un «tubo» 1D grueso, y por tanto una geometría gruesa legítima, no un falso positivo de fondo |
| Heisenberg coherente | ✘ | Ventana corta (r_w = 10). A escalas pequeñas el grupo crece de forma casi libre antes de que las relaciones actúen; el régimen asintótico de grado 4 no se alcanza con N ≈ 2·10⁴. No es concluyente |
| Cruces lentos (WS β = 0.001, atajos al 0.1 %) coherentes, como límite conocido | ✘ **a favor**: son INCOHERENTES | La coherencia atrapa los cruces lentos que engañaban al Nivel II |
| Caveman degenerado | Parcial (INTERMEDIO, f_deg 0.39) | — |

## 3. Auditoría (b): K APORTA (resultado oficial)

| Pregunta | Resultado |
|---|---|
| KA1 (no redundancia) | Spearman(mediana de κ, ln ρ) = **−0.17** (p = 0.16; n = 72) → **no redundante** |
| KA2 | **K corrige ER** (κ ≈ −0.80), **RGG3 en caja** (κ ≈ −0.01) y **RGG2 en caja** (κ ≈ +0.10): 3 familias. **No** corrige WS β = 0.001 / 0.002 (κ ≈ +0.17, como el anillo) ni los atajos al 0.1 % (κ ≈ −0.03). Exactamente lo predicho |
| Control | G-hom dentro de la banda 20/20. R (cliques) fuera 10/11. Del resto L/NL, fuera 14/25: todos los retazos y RR fuera; WS y atajos al 1 % dentro |

## 4. El espacio diagnóstico combinado (lo que pedía el Consejo)

**Qué atrapa cada observable:**

| Caso | ρ (Nivel I) | D_B2 (Nivel II) | **Coherencia γ** | **Curvatura K** | Localidad (A: ciclos cortos) |
|---|---|---|---|---|---|
| Retazos | ✔ | ✔ | ✔ | ✔ | — |
| Expansores (ER, RR) | ✘ (ER pasa) | ✔ (sin ventana) | ✔ (sin ventana) | ✔ | ✔ |
| WS / atajos escasos | ✔ | ✘ (WS β = 0.001 pasa) | ✔ | ✘ | — |
| Geometría con borde | ✘ (falso negativo) | ✔ | ✔ | ✔ | — |
| Árbol aleatorio | no medido | no medido | ✘ (coherente) | ✔ (f_neg 0.51 > banda 0.43) | ✔ (sin ciclos → NO_LOCAL) |
| Cliques / caveman | ✘ (ρ ≪ 1) | — | parcial (degenerado) | ✔ (+0.5) | — |
| C0 | ✔ | sin ventana | sin ventana (N = 729) | ✘ (+0.03) | ✔ local |

**Lectura.** Ninguna observable aislada separa todos los falsos positivos, pero **cada clase de falso positivo es atrapada por al menos una**. Los ángulos ciegos son disjuntos:
- la coherencia no ve árboles, pero K y la localidad sí;
- K no ve atajos, pero la coherencia y ρ sí;
- ρ no ve los bordes, pero la coherencia y K sí.

Esta complementariedad es una observación ***post hoc*** sobre los paneles ya corridos. Para usarla como regla de decisión combinada hace falta un prerregistro nuevo con grafos nuevos.

## 5. Dictamen del cerebro

1. **Coherencia multiescala:** bien definida, sin geometría por la puerta trasera, PARCIAL como discriminador.
   - **Corrige el falso negativo con borde** y **atrapa los cruces lentos**, dos defectos conocidos de la batería.
   - **No excluye los polímeros ramificados** (árboles críticos), que son amenables.
   - **Se incorpora como diagnóstico B-bis**, nunca como energía: CH-T5 muestra que, como objetivo, colapsaría a 1D. Y ahora sabemos además que no distingue árboles críticos.
2. **Curvatura de Ollivier:** se incorpora como **diagnóstico complementario de escala 1**: no redundante con ρ, corrige ER y los bordes. **No** es certificado ni criterio de dimensión.
3. **Para el Consejo, sobre el «ingrediente faltante»:**
   - La coherencia de vecindades (la propuesta del Consejo) **es satisfecha por los árboles críticos y por las cadenas**. No puede ser por sí sola el principio de selección geométrica.
   - Lo que separa la geometría de esos dos destinos en nuestro instrumento es **la presencia de ciclos cortos en todas las escalas** (localidad, punto A) y **la ausencia de ramificación negativa** (K).
   - Eso sugiere que el ingrediente debe combinar coherencia (Følner) con **riqueza de ciclos** (complejo de ciclos que «rellena» la estructura).
   - Lo dejo como **hipótesis de trabajo (nivel 4)**, no como diseño: informada por estos datos, requeriría el mismo tratamiento que D′.
4. **Propuesta:** prerregistrar una **regla de decisión combinada** (localidad + ρ + γ + K + Nivel II) y validarla sobre un panel nuevo con familias no usadas hasta ahora: árboles críticos de distinto tipo, polímeros con ciclos, geometrías hiperbólicas con ciclos, productos con factores grandes y bordes de distinta forma. Si pasa, sustituiría al Nivel I aislado en el preflight.
