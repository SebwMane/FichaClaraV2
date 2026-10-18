# Ω — L-DIM-1: resultados

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-ldim1`.
- **Prerregistro:** `docs/OMEGA_LDIM1_PRERREGISTRO.md` (§1–§9 + enmienda A1, commiteada antes de escribir el código).
- **Datos:** `results/ldim1/{features.jsonl, truth.jsonl, summary.json}`. Son 54 grafos (24 geometrías a N, 8 a ≈ 4N, 22 degenerados); extracción en 49 s; 5/5 pruebas.
- **Reparto:**
  - El código lo escribió un agente Sonnet. El cerebro lo revisó antes de correrlo (extractor ciego `extract(adj, id, rng)`, evaluador fiel a §7 + A1).
  - Las referencias las verificó un agente Haiku, con revisión del cerebro (§5).

## 1. Veredicto: **SELECCIÓN-NO-VIABLE**. D5 queda **refutado como mecanismo estático** (regla de muerte §7)

| p | Meseta (L3) | L1 | L2 | L3 | L4 | L7 | L8 | L9 | Nulo A | Nulo B | Nulo C | Falla |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1/2 | d* = 2, θ ∈ [0.56, 10] (1.25 déc.) | ✔ | ✔ | ✔ | ✘ | ✘ | ✔ | ✔ | 0.008 ✔ | 0.09 ✘ | 0.04 ✔ | L4, L7, nulo B |
| 1 | d* = 2, θ ∈ [1, 56] (1.75 déc.) | ✔ | ✔ | ✔ | ✘ | ✔ | ✔ | ✔ | 0.008 ✔ | 0.08 ✘ | 0.035 ✔ | L4, nulo B |
| 2 | d* = 2, θ ∈ [10, 1000] (2 déc.) | ✔ | ✔ | ✔ | ✘ | ✔ | ✔ | ✔ | 0.008 ✔ | 0.08 ✘ | 0.03 ✔ | L4, nulo B |

### d* combinado sobre la rejilla θ = 10⁻³…10³ (25 puntos)

| p | Dígitos (uno por punto de θ) |
|---|---|
| 1/2 | `1111111111122222222333333` |
| 1 | `1111111111122222222222333` |
| 2 | `1111111111111222222222222` |

### No codificación

| Prueba | Resultado |
|---|---|
| NC-2 | \|ρ_s(F, grado)\| ≤ 0.47 en toda la rejilla. F **no es un lector de grado** |
| NC-4 | 0 pares (φ, ψ) indistinguibles entre dimensiones distintas |
| NC-5 | Modas d*(p) = (2, 2, 2): **no confirmada** según su definición operativa |

## 2. Lectura

1. **L4 falla siempre, y por lo previsto.**
   - En **todos** los θ de las tres mesetas, la **clique K200** tiene F menor que la geometría ganadora: F = θ/199^p y φ = 0, que es L-DIM-T2 confirmada a tamaño finito.
   - En la mitad inferior de cada meseta también la socavan **los dos árboles** (Prüfer y binario subdividido). Tienen frontera relativa φ ≈ 0.03–0.05, menor que la de un retículo 2D.
   - La competencia la ganan los **extremos degenerados por ambos lados**: el compacto (clique) y el de baja expansión (árboles). Es decir, F-exp no distingue «poco crecimiento» de «geometría de baja dimensión».
2. **La dimensión la selecciona θ, no la estructura** (nivel 2, a partir de la tabla de d*).
   - A N fijo, d*(θ) es una **escalera monótona** 1 → 2 → 3 (→ 4). Cada entero tiene su ventana de θ, así que θ funciona como un **dial de dimensión**.
   - El requisito «un único θ para todas las d» se cumple formalmente, pero no impide esto: elegir θ es elegir d.
   - La meseta de d = 2 es la más ancha por un efecto de tamaño finito. Con N ≈ 2·10⁴, ψ apenas distingue d = 3, 4 y 5 (ψ ≈ 0.0050–0.0059), porque r* es pequeño (7, 4, 3).
3. **La escalera se desplaza con N** (L7).
   - A ≈ 4N, el estrato retículo pasa a d mayor en θ menores. Con p = 1/2 llega a d = 4 en los dos últimos θ, lo que hace fallar L7.
   - Coincide con L-DIM-T1: todo decae con r*, y θ_eff ∝ θ r*^{1−p(d−1)}. **La selección depende de (θ, N) a tamaño finito, y solo de p asintóticamente.** En ningún régimen es una propiedad estructural independiente de la parametrización.
4. **El nulo B falla (8–9 % > 5 %).**
   - Permutando las etiquetas d dentro de los estratos se pasan L1–L3 en ~8 % de los casos.
   - Con 9–15 instancias por estrato, L1–L3 tienen **poca potencia**: una meseta invariante no es, por sí sola, evidencia de selección.
   - Es una debilidad del criterio, no del candidato. Se registra para cualquier repetición.
5. **El nulo A es bajo (0.8 %) solo por L4.**
   - Funcionales aleatorios sobre (log φ, log ψ): el 37 % da un d* interior y el 11 % un d* interior invariante entre estratos.
   - **Producir un mínimo interior invariante es barato.** Lo caro es batir a los degenerados.

## 3. Contraste con las predicciones congeladas (§8)

| # | Predicción | Resultado |
|---|---|---|
| 1 | L4 falla para los tres p por la clique | ✔ (y además por los árboles) |
| 2 | d*(p) → 4, 3, 2 | ✘ a tamaño finito: la meseta es d = 2 para todos los p. ✔ en dirección: p menor desplaza la escalera hacia d mayor, y con p = 1/2 aparece d = 4 a 4N |
| 3 | NC-5 confirma la codificación por exponente | ✘ según la definición operativa (modas 2, 2, 2). La codificación observada es por **θ y N**, no (todavía) por p |
| 4 | Fracción sustancial de nulos A con «mínimo interior» | ✔ si se mira L1 (37 %); ✘ según la métrica prerregistrada L1 + L2 + L4 (0.8 %) |
| 5 | SELECCIÓN-NO-VIABLE; D5 refutado como mecanismo estático | ✔ |

## 4. Respuesta a la pregunta del Consejo

> «¿Existe una regla que pueda seleccionar d sin saber qué significa d?»

Para la familia mínima que traduce D5 (expansión + θ · fragilidad^p, con escala intrínseca), la respuesta es **no**.

- Ningún p pasa.
- Las tres causas son independientes:
  - **(i) compactificación:** clique (T2, nivel 1, confirmada);
  - **(ii) baja expansión degenerada:** árboles;
  - **(iii) d* lo fijan θ y N** (y p asintóticamente, T1).
- La **no codificación por grado sí se cumple** (NC-2 0.47, NC-4 0). El fallo no se debe a que F lea grado o densidad: es **estructural**.

## 5. Verificación de referencias (agente Haiku, revisada por el cerebro)

| Referencia | Estado |
|---|---|
| Pólya, *Math. Ann.* **84**, 149–160 (1921) | Recurrencia en Z^d para d ≤ 2, transitoriedad para d ≥ 3. ✔ |
| Nash-Williams, *Proc. Camb. Phil. Soc.* **55**, 181–194 (1959) | Σ 1/\|Π_n\| = ∞ ⇒ recurrencia. ✔ |
| Doyle–Snell, *Random Walks and Electric Networks*, Carus 22 (MAA, 1984) | ✔ |

**Corrección al prerregistro (§2, «Matiz»).**
- Que Σ 1/|S_r| < ∞ ⇔ d ≥ 3 sea el umbral de transitoriedad es cierto **en retículos** (Pólya, nivel 3).
- En grafos generales Nash-Williams solo da la dirección de recurrencia. La transitoriedad exige una desigualdad isoperimétrica: Thomassen, *Ann. Probab.* **20**, 1592–1600 (1992), con Σ f(k)^{−2} < ∞; y Varopoulos, *J. Funct. Anal.* **63**, 215–239 (1985).
- El agente no accedió a los textos completos, así que la verificación es bibliográfica; los enunciados los confirma el cerebro con la literatura estándar.

## 6. Dictamen del cerebro y propuesta al Consejo

1. **L-DIM-1 cerrada: D5-B (competencia estática → selección) refutada** para la familia mínima. La parte A de D5 (homogeneidad → integralidad) no se ha tocado y sigue siendo de nivel 3 externo.
2. **Lo que debería aportar cualquier mecanismo futuro, según los resultados:**
   - **(a) Extensión sin escala impuesta.** Tiene que excluir la clique y los árboles por razones **no energéticas**, por ejemplo topológicas, como el extremo único a toda escala de B-ter. Pero B-ter es ciego a d, así que la exclusión y la selección quedan separadas.
   - **(b) Invariancia de escala del criterio.** Un funcional cuyo ganador no dependa de N exige que los dos términos escalen igual en r, y entonces el cociente de exponentes (T1) vuelve a fijar d. **Dentro de funcionales de bolas con leyes de potencias, la selección invariante de escala y la no codificación por exponente son incompatibles** (nivel 2; T1 muestra que el ganador asintótico lo fija p).
3. **No se recomienda L-DIM-2** («¿por qué 3?»), porque su condición, que L-DIM-1 pase, no se cumple. Una línea razonable sería una sesión conceptual sobre (b): ¿hay alguna razón **independiente**, física o de información, para un exponente concreto? El ejemplo es la transitoriedad (p = 1 ↔ d ≥ 3 en retículos). Sin esa razón, cualquier selección será codificación.
