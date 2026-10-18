# Ω — W-1: ¿son resolubles los umbrales de paseo y son universales? (prerregistro, estático)

- **Fecha:** 2026-10-07.
- **Rama:** `claude/omega-w1`. Base congelada: `claude/omega-w0-congelado`.
- **Estado:** PRERREGISTRO, escrito antes de cualquier código o medida.
- **Mandato del Consejo.** W-1 aprobado: estático, simétrico en d, panel homogéneo y no homogéneo, leyes definidas de antemano, potencia previa, y separación explícita entre *detectar* y *generar* dimensión. Sin dinámica nueva.

---

## 0. Acta del Consejo y auditoría

### 0.1 Ratificado

| Decisión | Registro |
|---|---|
| W6 VÁLIDA **en su dominio calibrado**: «instrumento válido para la familia de residuos calibrada, no detector universal de todo crecimiento intermedio» | Adoptado literalmente. Se corrige el mapa de exclusiones |
| W-0 aprobada **como reducción**. «W-0 demuestra imposibilidad»: rechazado | Adoptado |
| Reformulación de W-0: «Entre las familias examinadas, las rutas sin parámetro continuo no han producido todavía un selector de dimensión; las candidatas basadas en leyes de paseo trasladan la selección a una elección discreta de ley, salvo que exista un principio independiente que determine esa ley» | **Sustituye** a la frase del §4.5 de W-0. Es más débil y más defendible |
| «Elección de ley = dial» solo si hay varias leyes admisibles, ninguna fijada por un principio independiente previo, y las elecciones dan resultados distintos | Adoptado. Coincide con la cláusula α de E2, ahora escrita explícitamente para leyes |
| Diagnóstico ≠ mecanismo | Adoptado como columna obligatoria (§5) |
| Cuatro desenlaces W1-A…D | Adoptados, con la corrección Ω-1 |

### 0.2 Omisiones y correcciones

**Ω-1. W1-A es inalcanzable por diseño.**
- W-1 es estático: aplica leyes de paseo a geometrías dadas. Por construcción mide (geometría → ley → d) y nunca genera (ley → dinámica → geometría).
- Para afirmar W1-A, «hay evidencia de que selecciona, no solo mide», haría falta una dinámica, y el Consejo la ha prohibido por ahora.
- Por tanto el máximo alcanzable en W-1 es **W1-B** (detector robusto). Y si varias leyes resolubles dan umbrales distintos, **W1-C**.
- Así queda escrito para que nadie lea un W1-B como si fuera W1-A.

**Ω-2. Hay un conflicto de tamaño que el Consejo intuye pero no cuantifica.**
- Un paseo en un toro de N nodos solo ve geometría «infinita» mientras su desplazamiento sea menor que el lado, es decir, durante T ≲ R².
- Con N = 10⁶ y d = 6, el lado es 10, así que T ≲ 25 pasos.
- Los umbrales de dimensión alta (d_c = 3 o 4) se juegan exactamente donde el tiempo disponible es menor. Es la forma cuantitativa de W-T5.
- W-1 usa un tiempo máximo **intrínseco**: T_max = ⌊R²⌋, con R el radio de media masa de RC-3, que no usa coordenadas. Se declara de antemano que las dimensiones altas pueden salir INSUFICIENTES.

**Ω-3. En la dimensión crítica, la señal es logarítmica.** El exponente de crecimiento en una ventana finita es pequeño pero positivo, y cae entre «crece» y «acotado». La regla de resolubilidad (§3) admite AMBIGUA en el valor crítico, sin exigir una separación limpia ahí.

**Ω-4. Paridad.** En grafos bipartitos, dos paseos simples que empiezan a distancia impar nunca coinciden. Por eso se usan paseos **perezosos** (quedarse con probabilidad 1/2) en todas las leyes.

---

## 1. Leyes (congeladas; ninguna se añade ni se cambia después)

Todas usan paseos perezosos sobre el propio grafo, sin coordenadas, desde un vértice inicial uniforme de la componente gigante.

| Ley | Observable C(T) en tiempo T | Umbral teórico en ℤ^d (nivel 3) |
|---|---|---|
| **L1 Retorno** | Número de visitas al origen en los pasos 1…T | crece si d ≤ 2 |
| **L2 Colisión** | Número de instantes t ≤ T en que dos paseos independientes desde el mismo origen están en el mismo vértice | crece si d ≤ 2 |
| **L3 Cruce de 2 trayectorias** | \|rango(W₁) ∩ rango(W₂)\| hasta T, con ambos paseos desde el mismo origen | crece si d ≤ 4 |
| **L4 Cruce de 3 trayectorias** | \|rango(W₁) ∩ rango(W₂) ∩ rango(W₃)\| hasta T | crece si d ≤ 3 (en red; nivel 3 sin fuente primaria) |

**Exponente:**
- a = pendiente por mínimos cuadrados de ln(1 + C̄(T)) frente a ln T en T ∈ {T_max/16, T_max/4, T_max}.
- C̄ es la media sobre M = 256 orígenes, con paseos independientes por origen.
- Error típico SE(a) por bootstrap sobre orígenes, con 500 remuestras.

**Clase por (geometría, ley):**
- CRECE si a ≥ 0.10;
- ACOTADA si a ≤ 0.05;
- AMBIGUA en otro caso;
- **INSUFICIENTE** si SE(a) > 0.03 o T_max/16 < 10.

## 2. Panel

**Capa 1, homogéneo.** Dos tamaños, N₁ ≈ 1.25·10⁵ y N₂ ≈ 10⁶:

| d | Geometrías |
|---|---|
| 1 | Toro Z¹ |
| 2 | Toro Z²; triangular; RGG₂ k = 8 (2 semillas) |
| 3 | Toro Z³; FCC; BCC; RGG₃ k = 8 (2 semillas) |
| 4 | Toro Z⁴ |
| 5 | Toro Z⁵ |
| 6 | Toro Z⁶ (lado ≥ 5) |

**Capa 2, no homogéneo.** Solo N₂. La dimensión de referencia d_g es la dimensión de crecimiento por construcción:

| Geometría | d_g |
|---|---|
| peine-1 (ciclo de L con dientes-camino de longitud L) | 2 |
| peine-2 (toro Z² de lado s con dientes-camino de longitud s) | 3 |
| cilindro C_n × C₁₀ | 1 |
| Z³ diluido (cada arista se conserva con p = 0.6; componente gigante; 2 semillas) | 3 |
| retazos `patchwork(n, 3)` | sin d_g: solo se reporta |

**RNG:** `rng_from_key((MASTER_W1, id, semilla, tamaño, k))` con `MASTER_W1 = 20261023`. k = 0 para la geometría, k = 1 para los orígenes, k = 2 para los paseos.

## 3. Criterios (congelados)

**Resolubilidad de una ley** en la capa 1, con la clase en N₂:
- La ley es **RESOLUBLE** si existe un corte d* tal que:
  - toda geometría homogénea con d < d* es CRECE;
  - toda geometría con d > d* es ACOTADA;
  - las de d = d* pueden tener cualquier clase;
  - ninguna geometría con d ≠ d* es INSUFICIENTE.
- d̂_c es el intervalo de los d* que cumplen la condición.
- **Estabilidad en N:** la misma condición debe cumplirse con las clases en N₁ y dar un intervalo que se solape. Si no, la ley es INESTABLE-N.
- Si alguna geometría necesaria es INSUFICIENTE, la ley es **NO RESUELTA**. Eso es distinto de no resoluble.

**Universalidad** (capa 2) de una ley resoluble:
- Para cada geometría no homogénea con d_g, su clase debe coincidir con la que la ley asigna a d_g en la capa 1.
- Si d_g ∈ d̂_c (valor crítico), cualquier clase se acepta.
- Una sola discordancia determinada (CRECE frente a ACOTADA) ⇒ **NO UNIVERSAL**.

**Desenlace de W-1** (en este orden):

| Desenlace | Condición |
|---|---|
| **W1-D** | Ninguna ley RESOLUBLE y estable en N, por insuficiencia o por ausencia de corte |
| **W1-C** | ≥ 2 leyes RESOLUBLES y estables con intervalos d̂_c disjuntos: la selección depende de la ley |
| **W1-B** | Al menos una ley RESOLUBLE, estable y UNIVERSAL, y todas las resolubles comparten d̂_c: detector robusto, no selector |
| **W1-B⁻** | Hay leyes resolubles, pero ninguna es universal: detector dependiente de la homogeneidad |
| W1-A | Inalcanzable en W-1 (Ω-1) |

Además se registran las subclases relevantes, por ejemplo «W1-C con L2 no universal».

## 4. Potencia (prerregistrada)

- El efecto mínimo de interés es una separación de clases: ≥ 0.05 entre el umbral CRECE (0.10) y el de ACOTADA (0.05).
- La resolución exigida es SE(a) ≤ 0.03 por celda. Con eso, una diferencia de 0.05 entre la frontera y el valor real se detecta con ≈ 1.7 SE.
- Es una potencia **modesta**, y queda declarada así.
- Las celdas que no alcancen SE ≤ 0.03 con M = 256 se declaran INSUFICIENTES. No se aumenta M después de ver datos.

## 5. Tabla obligatoria por ley (petición del Consejo)

| Pregunta | Cómo se responde |
|---|---|
| ¿Detecta dimensión? | RESOLUBLE, y con qué d̂_c |
| ¿Es robusta a la heterogeneidad? | UNIVERSAL o NO UNIVERSAL (capa 2) |
| ¿Es robusta a N? | Estable o INESTABLE-N |
| ¿Depende del número de caminantes? | Comparando L1 (1), L2 y L3 (2) y L4 (3) |
| ¿Depende de la representación red/continuo? | No se mide (W-1 solo usa redes). Se responde con nivel 3: R4 da 3 en el continuo y 4 en red para k = 2 |
| ¿Depende de la ley elegida? | Comparando los d̂_c |
| ¿Puede generar dimensión o solo medirla? | **Solo medirla** (Ω-1), para todas las leyes |

## 6. Predicciones del cerebro (no deciden)

- **Valores críticos.** Cada ley tendrá su d crítico AMBIGUO (señal logarítmica): L1 y L2 en d = 2, L3 en d = 4, L4 en d = 3.
- **Dimensiones altas.** En d ≥ 4 el T_max intrínseco será pequeño. Riesgo alto de INSUFICIENTE en Z⁵ y Z⁶ a N₂. Eso dejaría L3 NO RESUELTA.
- **Leyes resolubles.** L1 y L2 serán RESOLUBLES con d̂_c ∋ 2, y L4 probablemente con d̂_c ∋ 3.
- **Heterogeneidad.** L2 será **NO UNIVERSAL**: en el peine-1 se predice ACOTADA aunque d_g = 2 (Krishnapur y Peres). L1 será universal en el peine-1 (el peine es recurrente).
- **Desenlace más probable:** **W1-C** (L1/L2 frente a L4), con L2 no universal y L3 no resuelta.
- **Alternativa:** W1-D por insuficiencia en d alto.

## 7. Ejecución

- Agente Sonnet: `tools/w1_walks.py`, tests y `results/w1/`.
- Paseos vectorizados sobre CSR. Rangos por trayectoria con `np.unique`.
- Revisión del cerebro y resultados en `docs/OMEGA_W1_RESULTADOS.md`.
