# Ω — Fase 2: resultados y dictamen del cerebro

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-fase2-theta`.
- **Prerregistro:** `docs/OMEGA_FASE2_PRERREGISTRO.md` (§0–§4 y enmiendas F2-A1 y F2-A2).
- **Base congelada:** `claude/omega-c0-congelado`.

## 1. C0-Θ: ¿sobrevive a temperatura la fase local de C0? — **FRÁGIL-1/N**

216 corridas oficiales con dt = 0.0025 y 40 000 pasos, más 16 de sensibilidad con dt = 0.00125 y 80 000 pasos.

### 1.1 Supervivencia (Θ-S), desde los finales oficiales de C0

Cada casilla es el veredicto prerregistrado (≥ 2/3 semillas). Inicio R; entre paréntesis, inicio U cuando difiere. En la parte superior, la mediana de f_bg en cada columna:

| f_bg (mediana) | 0.02 | 0.16 | 0.36 | 0.27 | 0.27 |
|---|---|---|---|---|---|
| **Celda (c\*, k\*)** | **Θ = 0** | **0.1 Θ_N** | **0.3 Θ_N** | **1 Θ_N** | **3 Θ_N** |
| 19 (1, 6), N = 216 / 343 | sobrevive | sobrevive | MIXTA: local, J ≈ 0.33 | DESTRUIDA: denso | DESTRUIDA |
| 37 (2, 8), N = 216 / 343 | sobrevive | sobrevive | MIXTA: local, J ≈ 0.33 | DESTRUIDA | DESTRUIDA |
| 55 (4, 12), N = 216 | sobrevive | sobrevive (U: MIXTA) | DESTRUIDA: no local | DESTRUIDA | DESTRUIDA |
| 55, N = 343 | sobrevive | sobrevive | MIXTA: local, J ≈ 0.32 | DESTRUIDA | DESTRUIDA |
| 73 (8, 16), N = 216 | sobrevive | sobrevive (U: MIXTA) | DESTRUIDA: no local | DESTRUIDA | DESTRUIDA |
| 73, N = 343 | sobrevive | sobrevive | MIXTA: local, J ≈ 0.32 | DESTRUIDA | DESTRUIDA |

- **Frontera:** Θ_c/Θ_N = 0.1 en las 4 celdas, en N = 216 y 343, con inicio R. Con inicio U vale 0.1 en las celdas 19 y 37, y 0 en la 55 y la 73.
- **Veredicto: FRÁGIL-1/N.**

### 1.2 Formación (Θ-F), desde un inicio genérico R con N = 216

| c\* | Celdas | Resultado con Θ = 0.3 Θ_N |
|---|---|---|
| 1, 2 | 19, 37 | **FORMA**: la fase local también aparece a temperatura |
| 4, 8 | 55, 73 | **NO FORMA**: los estados quedan dispersos y no locales |

### 1.3 Contraste con las predicciones registradas

| Predicción | Resultado | ¿Acertó? |
|---|---|---|
| P-Θ1: sobrevive con Θ ≤ 0.1 Θ_N | Sobrevive en 12/12 combinaciones con inicio R | ✔ |
| P-Θ2: se destruye con Θ ≥ 3 Θ_N | Se destruye **ya con 1 Θ_N**, en 100 % de las corridas (DENSO-TRIVIAL) | ✔, con una frontera más baja de lo previsto |
| P-Θ3: la frontera escala como 1/N | **Colapso casi exacto.** Con el mismo Θ/Θ_N, f_bg coincide entre N = 216 y 343 y entre celdas (0.15–0.18 con 0.1 Θ_N; 0.35–0.38 con 0.3 Θ_N), y Θ_c/Θ_N es idéntica en los dos N | ✔ |

El cálculo de campo medio de §1.4 del prerregistro acierta la **escala** Θ_N = G_bg·k*/N. Subestima el tamaño del fondo: con 0.1 Θ_N predecía f_bg ≈ 0.09 y se mide 0.16.

### 1.4 Validez del integrador (F2-A1 y F2-A2)

- **Control a Θ = 0:** conserva la clase y J ≥ 0.97 en todas las corridas.
  - S/LB empeora ≤ 0.9 % con N = 216, pero 1.2–2.4 % con N = 343, por encima del criterio del 1 %. Es un sesgo de pared O(dt).
- **Sensibilidad con dt/2 (N = 343):**
  - en las 16 corridas, la clase y el lado de J coinciden con la corrida oficial, y f_bg cambia ≤ 0.07;
  - el sesgo del control baja a 0.2–1.1 %.
  - Por F2-A2, el veredicto se mantiene.
- **Estacionariedad:** la deriva de S entre los dos últimos cuartos es ≤ 0.2 % en los brazos con Θ ≤ 0.3 Θ_N, y ≤ 2 % en los destruidos.

### 1.5 Lectura

**Verificado (N ∈ {216, 343}, tiempo 100):**
- la fase local de C0 existe solo para Θ ≲ 0.1–0.3 Θ_N, con Θ_N ∝ 1/N;
- por encima, el sistema pasa a estados densos de alta energía; con 3 Θ_N y c* ≥ 4, S llega a ser positiva. La entropía domina.

**Explicación (nivel 2):** cada nodo tiene ~N pares de fondo cuyo coste por par, G_bg = ψ* − a, es intensivo. La entropía del fondo crece con N y su coste no, de modo que la temperatura que tolera la fase cae como 1/N.

**Extrapolación (nivel 4, fuertemente apoyada por el colapso):** a Θ > 0 fija, la fase local de C0 desaparece cuando N → ∞. **Es un fenómeno de Θ → 0, no una fase termodinámica.**

**Matiz:** con 0.1–0.3 Θ_N la clase local persiste, pero la configuración concreta se reorganiza (J ≈ 0.3–0.6). A temperatura baja la fase se comporta como un **líquido local**, no como una estructura congelada.

## 2. Información de apoyo para C1 (descriptiva, no son pruebas de C1)

### 2.1 Conteos locales (k, c, q) por arista — `results/c1_local_stats/`

| Grafo | k | c (vecinos comunes) | q (4-ciclos por arista) | CV de q |
|---|---|---|---|---|
| T³ 9³ | 6 | 0 | 4 | 0 |
| Q₆ (hipercubo) | 6 | 0 | 5 | 0 |
| RGG3 k12, N = 729 | 12.2 | 5.8 | 51 | 0.77 |
| Finales C0-R, N = 729 (celdas 18–38) | 12.5–13.8 | 4.6–5.1 | 34–48 | 0.36–0.53 |
| ER k12 | 12 | 0.18 | 2.4 | 0.76 |

**Lectura:** hasta segundo orden, la fase C0 (no variedad, no certificada) y una geometría aleatoria auténtica (RGG3, certificada con N = 729) tienen **las mismas estadísticas locales**, y C0 es incluso más homogénea. Ninguna recompensa sobre (k, c, q) puede preferir el RGG3 sin preferir también la fase C0.

En la familia reticular, los conteos solo seleccionan dimensión por ajuste de parámetros: T³ tiene q = 4 y Q₆ tiene q = 5, con el mismo k = 6 (prerregistro §1.3).

### 2.2 Crecimiento de bolas en escala intermedia — `results/c1_ball_growth/`

CV entre nodos de |B_r|, el número de nodos a distancia ≤ r:

| Grafo | r = 2 | r = 3 | r = 4 | Tendencia |
|---|---|---|---|---|
| RGG3 k12 (3 semillas) | 0.19–0.25 | 0.13–0.19 | 0.09–0.15 | **baja**: se autopromedia |
| C0, celdas 18–38 (18 grafos) | 0.12–0.14 | 0.22–0.32 | 0.22–0.37 | **sube** en 18/18 |
| C0, celdas 55–74 (9 grafos) | 0.14–0.20 | 0.23–0.29 | 0.16–0.26 | sube y luego se estanca por encima del RGG3 |

**Lectura (descriptiva):** la diferencia entre localidad y geometría **no está en los conteos de radio ≤ 2, sino en la escala intermedia**.
- Una geometría aleatoria es irregular de cerca y se vuelve homogénea al mirar más lejos (autopromediado).
- La fase C0 es regular de cerca y se vuelve heterogénea más lejos.

Es la firma de un pegado incoherente de vecindarios locales. Coincide con los códigos F9 (no variedad) y F5 (dependencia de la métrica) del certificado.

## 3. Dictamen del cerebro

### 3.1 Sobre lo que dicen los datos

1. **C0 queda cerrada como mecanismo**, tanto a Θ = 0 (F1-NEGATIVO) como a Θ > 0 (FRÁGIL-1/N). Su aportación permanente es haber demostrado que **localidad ≠ geometría** y que la dicotomía vacío/clique no era inevitable.
2. **El problema tiene ahora tres capas medidas:**

   | Capa | Estado |
   |---|---|
   | A. Localidad | C0 la consigue, pero solo a Θ → 0 |
   | B. Variedad | Falla la homogeneidad mesoscópica: CV de |B_r| creciente |
   | C. Dimensión | No es accesible mientras B falle |

3. **Las tres lecciones que debe respetar cualquier propuesta futura:**
   - **L-mesoescala:** el ingrediente que falta actúa en radio ≥ 3, no en conteos locales. Una C1 de 4-ciclos saturados repetiría el problema de C0 (§2.1).
   - **L-entropía:** en un modelo relacional sin espacio, cada nodo tiene ~N relaciones posibles. Sin un coste por relación que crezca con N, o sin una ligadura en el número de relaciones, la entropía destruye la localidad a cualquier Θ fija (§1). En la literatura de grafos aleatorios exponenciales dispersos, el coste por arista se escala como log N. Es una referencia externa (nivel 3) que hay que verificar antes de citarla.
   - **L-meseta (C1-D):** una dimensión solo cuenta como emergente si se selecciona en una región abierta de parámetros, no en un punto ajustado.

### 3.2 Decisiones que tomo (dentro de lo que el Consejo ya autorizó)

| Decisión | Motivo |
|---|---|
| Congelar Fase 2 / C0-Θ cuando el Consejo lo revise | Mismo procedimiento que con C0 |
| **No implementar C1 en su forma de 4-ciclos saturados** | §2.1: no separa C0 del RGG3; §1.3 del prerregistro: convierte D en parámetro |
| No correr C1-A..D todavía | Requieren elegir la forma funcional, que es decisión del Consejo |
| O-05 sigue secundaria | Sin cambios |

### 3.3 Propuestas para que el Consejo elija (no ejecutadas)

**P1. Ingrediente mesoscópico de homogeneidad (recomendada).** Un término que penalice la **varianza entre nodos** del tamaño de bola |B_r| para r ≥ 3, normalizada por su media, sin fijar su valor.
- *Por qué:* es exactamente lo que separa al RGG3 de C0 en §2.2, y no impone exponente ni dimensión. Ataca la capa B sin tocar la C.
- *Riesgos que hay que analizar antes de simular (puerta tipo C1-A..D):*
  - |B_r| no es diferenciable en los pesos (habría que suavizarla, por ejemplo con un núcleo de calor de tiempo t);
  - hay homogéneos triviales: cliques, retículos de cualquier dimensión y expansores regulares también tienen CV 0. B y C quedan separados y la dimensión sigue abierta, lo cual es honesto pero no resuelve C;
  - el coste es O(N³) por evaluación.
- *Preprueba barata:* calcular la versión suavizada sobre T³, RGG3, Q_d, expansores y finales C0 antes de escribir ninguna dinámica.

**P2. Coste extensivo por relación (L-entropía).** Repetir C0-Θ con un potencial químico μ(N) por par, o con ligadura de número de aristas, para comprobar si la fase local pasa a ser termodinámica. Es barato y diagnóstico, pero no da geometría: solo resuelve la capa A a Θ > 0.

**P3. Selección dinámica o histórica.** C0 mostró que el resultado depende de la condición inicial (paisaje vidrioso). Una alternativa a «la energía selecciona la geometría» es que la selecciona el **proceso**: crecimiento secuencial o reglas locales de adición. Es un cambio de paradigma y requiere un prerregistro conceptual propio. Hay que vigilar que no se cuele el orden causal, aplazado por el Consejo.

**Orden sugerido:** P2, por ser barato y cerrar la capa A, junto con la preprueba de P1, que es descriptiva. Después, decisión del Consejo sobre P1 completa o P3.

## 4. Escala de afirmaciones

| Nivel | Afirmación |
|---|---|
| 1 (demostrado o exacto) | C0-T1/T2 (fase anterior); conteos (k, c, q) de retículos e hipercubos (§2.1) |
| Verificado en dominio (N ≤ 343 o 729) | C0-Θ FRÁGIL-1/N con colapso en Θ/Θ_N; formación con c* ≤ 2; igualdad de conteos locales C0 ≈ RGG3; CV de |B_r| creciente en C0 y decreciente en RGG3 |
| 2 (explicación) | Entropía del fondo ~N frente a coste intensivo; pegado incoherente de vecindarios |
| 3 (externo, por verificar) | Rigidez local de Z^d; escalado log N en grafos aleatorios exponenciales dispersos |
| 4 (extrapolación) | Desaparición de la fase C0 con N → ∞ a Θ fija; que P1 baste para la capa B |

## 5. Archivos

| Qué | Dónde |
|---|---|
| Código | `omega/c0/langevin.py`, `tests/test_c0_langevin.py` |
| Herramientas | `tools/c0_theta.py`, `tools/c1_local_stats.py`, `tools/c1_ball_growth.py` |
| Resultados | `results/c0_theta/`, `results/c0_theta_dt/`, `results/c1_local_stats/`, `results/c1_ball_growth/` |
| Pesos finales (no versionados) | `runs/c0_theta/out/` |
