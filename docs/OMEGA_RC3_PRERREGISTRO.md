# Ω — RC-3: juez de exclusión A− (unilateral, con abstención, de escalado) — Prerregistro

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-rc3`.
- **Base congelada:** `claude/omega-a0-congelado`.
- **Mandato del Consejo:**
  - RC-2 INVÁLIDA por sesgo dimensional demostrado;
  - RC-3 autorizado como **diseño nuevo**: no certificar geometría, sino **excluir degenerados** sea cual sea la dimensión del superviviente;
  - prioridad A3 (escalado) sobre A1 (nulo emparejado);
  - separación A− (exclusión) / A+ (certificación);
  - todavía sin dinámica.
- **Estado:** commiteado antes de escribir el código y antes de generar ningún grafo.

## 0. Evaluación del dictamen

### 0.1 Ratificado

| Punto | Motivo |
|---|---|
| RC-2 no se repara *post hoc* | Regla prerregistrada §II.3 de L-A-0 |
| Lección 17, en la formulación fuerte del Consejo | «La ausencia de d en la fórmula no demuestra independencia de d; hay que demostrarla frente a cambios de dimensión, escala y discretización» |
| La curvatura sigue siendo un instrumento **descriptivo** | Lo invalidado es la ventana absoluta como certificado. Las auditorías anteriores de K siguen vigentes |
| A− antes que A+; un juez de exclusión puede abstenerse | Es la reformulación correcta: las propiedades negativas (expansión excesiva, ramificación, varios extremos, incoherencia) no necesitan caracterizar la geometría y, por tanto, no necesitan d |
| Instrumento ≠ mecanismo | Lección 18 (§5) |
| Tabla de clases de mecanismos | Se adopta con dos precisiones (§5) |

### 0.2 Omisiones del dictamen y su justificación (lo que el Consejo no dijo y cambia el diseño)

**O1. Un juez unilateral con abstención necesita métricas propias, o es trivial.**
- Un juez que se abstiene siempre cumple «no excluye geometrías» sin hacer nada.
- Hacen falta dos medidas separadas:
  - **(a) tasa de exclusión falsa sobre geometrías, estratificada por d.** Es la prueba de ceguera a d;
  - **(b) poder de exclusión por familia degenerada.**
- La sensibilidad y la especificidad de un clasificador binario no sirven aquí (§3).

**O2. A− solo excluye clases catalogadas.**
- «No excluido» significa «no pertenece a ninguna clase degenerada **conocida**», no «no degenerado».
- Una dinámica puede acabar en una clase degenerada nueva y pasar A−. El 2-árbol ya fue un ejemplo: una clase nueva descubierta por accidente.
- Consecuencias:
  - A− es **necesario, no suficiente**;
  - el residuo de toda dinámica futura debe inspeccionarse y el catálogo ampliarse;
  - es un juez de **mundo abierto**.

**O3. Ningún juez a N finito puede ser ciego a d para toda d** (nivel 2, por escalado).
- Una geometría de dimensión d con N vértices tiene radio R ~ N^{1/d}.
- Cuando N^{1/d} es del orden de unos pocos saltos, la estructura es **indistinguible** de una de diámetro logarítmico.
- La ceguera a d solo puede pedirse **hasta una resolución declarada d_max(N)**. Por debajo de esa resolución, el juez debe **abstenerse**, no excluir.
- Respuesta a la pregunta del Consejo («¿cómo medir escalas resolubles sin d?»): con cantidades **intrínsecas del propio grafo**, sin d:
  - la ventana r_w;
  - el exponente de escalado δ entre N y 8N (§1).
- Esto sustituye el requisito «≥ 2 escalas», que excluía, por una **abstención**.

**O4. Toda estadística local a grado fijo está sesgada en d alta. Por eso el nulo emparejado (A1) no basta.**
- En una geometría discreta de grado k fijo, al crecer d la fracción de vecinos comunes de una arista decae. El volumen de la intersección de dos bolas a distancia comparable a su radio cae exponencialmente con d.
- Los vecindarios tienden así a los de un grafo aleatorio disperso, casi arbóreo.
  - Es la tendencia observada: κ de RGG2/3/4 = +0.10 / −0.03 / −0.11, con la de un expansor en −0.8.
  - La localidad por ciclos cortos tiende a cero.
- Un nulo emparejado corrige la **discretización**, pero no esta convergencia: en d alta, la geometría y su nulo se parecen localmente **de verdad**.
- **Decisión:** RC-3 **no contiene ninguna estadística local**. Fuera κ, fuera la localidad por triángulos o 4-ciclos, y fuera el nulo emparejado local.
  - Esto se desvía de la lista del Consejo, que pedía «nulos emparejados».
  - El «nulo emparejado» que se usa es interno: **el mismo proceso generador a dos tamaños**, N y 8N.

**O5. Varias clases degeneradas del Consejo tienen una sola firma de escalado.**
- Clique, expansor, mundo pequeño, 2-árbol / k-árbol, red apoloniana y BA tienen todas diámetro O(log N) o O(1). Su exponente de escalado δ ≈ 0.
- Una sola prueba cubre las seis clases.
- Reduce el problema de comparaciones múltiples y el número de umbrales.

## 1. Especificación de RC-3 (congelada)

**Entrada.** Un **par** (G_N, G_{8N}) del mismo proceso generador con los mismos parámetros. N ≈ 10⁴; las redes deterministas usan el lado más cercano.
- Para una dinámica futura: el mismo proceso a dos tamaños (como E4).

**Estadísticas.** Son todas de gran escala, intrínsecas y sin d.

1. **R(G), radio de media masa interpolado.**
   - m(r) = media de |B_r| del observador congelado `sampled_ball_profile` (400 fuentes), con m(0) = 1.
   - R es el r en que m(r) cruza n_giant/2, interpolado linealmente entre r − 1 y r.
   - Si no cruza dentro del perfil devuelto, R = r_cap.
2. **δ = ln(R(G_{8N}) / R(G_N)) / ln(n_giant(8N) / n_giant(N)).** Su valor es 1/d en una geometría y ≈ 0 con diámetro logarítmico o acotado.
3. **Anillos** (`annulus_profile`/`annulus_status`) en los dos tamaños.
4. **Coherencia** (`edge_coherence_profile`/`coherence_status`) en los dos tamaños.

**Exclusiones** (cada una es la firma de una clase degenerada del catálogo):

| Id | Clase | Condición |
|---|---|---|
| X1 | Diámetro logarítmico o acotado: clique, expansor, mundo pequeño, k-árbol, apoloniana, BA | δ < 0.125 |
| X2 | Ramificación persistente: árboles, cactus, árboles de bloques | `annulus_status` = RAMIFICADO en **los dos** tamaños |
| X3 | 1D | DOS_EXTREMOS en los dos tamaños, **o** δ > 0.75 |
| X4 | Pegado incoherente (retazos) | `coherence_status` = INCOHERENTE en los dos tamaños |

- **Veredicto por par:** EXCLUIDO(Xi…) si alguna Xi se activa; en otro caso, **NO-EXCLUIDO** (residuo, que incluye las abstenciones).
- Se informan también la resolución (r_w y número de escalas a 8N) y los valores R y δ.
- **Ausentes por diseño:** localidad, κ, Nivel II, d y cualquier umbral calibrado con una geometría concreta.

**Justificación a priori de los umbrales** (teoría, no datos):
- **X1:** con diámetro logarítmico, δ ≈ ln(ln 8N₀ / ln N₀)/ln 8 ≈ 0.098 para N₀ = 10⁴. Una geometría tiene δ = 1/d.
  - δ_c = 1/8 implica una **resolución declarada d_max = 7**. Por encima, el juez puede excluir geometrías auténticas: es el límite O3.
- **X3:** δ > 0.75 ⇔ 1/δ < 1.33. No afecta a ninguna d ≥ 2 (δ ≤ 0.5).
- **X2 y X4:** exigen persistencia en los dos tamaños, para que un artefacto de ventana corta no excluya.

## 2. Panel de validación (pares N ≈ 10⁴ y 8N)

- Clave: (20261016, familia, semilla, índice de tamaño).
- Semillas 0–1 en las familias aleatorias y 0 en las deterministas.
- \* = familia **nueva**: nunca se usó para fijar RC-1, RC-2 ni RC-3. **El veredicto se calcula solo sobre las familias nuevas.** Las repetidas se informan aparte.

| Grupo | Familias |
|---|---|
| **V, d = 2** | Toro cuadrado 100² / 283² · RGG2 k8 · triangular en caja abierta* · panal hexagonal |
| **V, d = 3** | T³ 22³ / 44³ · RGG3 k8 · red FCC* (k = 12) · red BCC* (k = 8) · RGG3 k12 en caja · T³ + 20 % de diagonales |
| **V, d = 4** | T⁴ 10⁴ / 17⁴ · RGG4 k8 · T⁴ + 20 % de diagonales* · RGG4 k16 |
| **V, d = 5** | T⁵ 6⁵ / 11⁵* · RGG5 k16* |
| **V, d = 6** | T⁶ 4⁶ / 7⁶* |
| **V-ext** (no degenerada, no euclídea) | Heisenberg 22 / 44 |
| **X** | Clique K200 / K1600* · 4-árbol aleatorio* · apoloniana · BA m = 3* · RR k3* · ER k4* (gigante) · árbol de Prüfer · árbol binario subdividido · árbol de cubos 5³* · árbol binario × C₁₀* · RGG3 k12 + 0.1 % de atajos* · cuadrado WS β = 0.01* (reconexión de un extremo) · retazos b3 · cactus · 2-árbol |
| **U** | Cilindro C × C₅* · RGG1 (anillo aleatorio) k10* · escalera C × P₄ · tubo C × RR |

## 3. Métricas y veredicto (sobre las familias nuevas)

| Medida | Definición |
|---|---|
| **FE** (exclusión falsa) | Fracción de instancias V excluidas, global y por d |
| **EP** (poder de exclusión) | Fracción de instancias X + U excluidas, global y por familia |

| Veredicto | Condición |
|---|---|
| **RC3-VÁLIDA** | FE global ≤ 0.05, **ninguna** d con más de una instancia excluida, EP global ≥ 0.90 y cada familia X/U nueva con todas sus instancias excluidas salvo como máximo una |
| **RC3-PARCIAL** | FE global ≤ 0.10 y EP global ≥ 0.75 |
| **RC3-INVÁLIDA** | En otro caso |

- Heisenberg (V-ext) se informa aparte y no entra en FE.
- **Regla del Consejo:** RC-3 no se toca después de ver resultados. Si falla, se acepta que **todavía no hay juez**, y el problema A vuelve al Consejo.

## 4. Predicciones del cerebro (congeladas)

1. **FE.**
   - El riesgo está en d ≥ 5: T⁵ y RGG5 tienen R de unos pocos saltos a N = 10⁴, y la interpolación es ruidosa. Predicción: 0–1 exclusiones falsas en d = 5–6.
   - Heisenberg excluido por X4 (era INCOHERENTE en L-CIC-0b).
2. **EP.**
   - Clique, 4-árbol, apoloniana, BA, RR3, ER4 y 2-árbol, por X1.
   - Prüfer, árbol de cubos, árbol × C₁₀ y cactus, por X2.
   - Retazos, por X4.
   - Las 1D, por X3.
3. **Fallo previsto:** RGG3 k12 + 0.1 % de atajos **no se excluye**. Su δ queda intermedio, en 0.15–0.25, porque el cruce lento no llega a régimen logarítmico entre N y 8N.
4. Cuadrado WS β = 0.01 excluido por X1, aunque con confianza baja.
5. **Veredicto previsto: RC3-PARCIAL**, porque la familia de atajos incumple la regla por familia.

## 5. Registro conceptual (se añade al mapa negativo con los resultados)

**Lección 18 (instrumento ≠ mecanismo, Consejo).** Un instrumento detecta G; un mecanismo explica Ω → G. Coherencia, curvatura, anillos y ciclos son instrumentos.

**Tabla de clases.** Se adopta la del Consejo con dos precisiones:
- «densidad/regularización» se refiere a S0 de Ω-1.1 y a Ω-B;
- «crecimiento por curvatura D» es P3-D.

Quedan **abiertas**:
- dinámica que produzca coherencia;
- dinámica que produzca un extremo;
- orden o causalidad (externa);
- rango algebraico (conceptual).
