# Ω — L-A-0, Parte II: validación fuera de muestra de RC-2. Resultados

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-a0`.
- **Prerregistro:** `docs/OMEGA_A0_PRERREGISTRO.md`.
- **Datos:** `results/a0_rc2/` (42 grafos con N ≈ 2·10⁴; 112 s; 6/6 pruebas).
- **Código:** `tools/a0_rc2.py`, escrito por un agente Sonnet y revisado por el cerebro.
  - `rc2_conditions` se ajusta a §II.1.
  - Los instrumentos se importan sin cambios y con las mismas subclaves de RNG.
- **Desviación declarada (viabilidad, no decisión).** Una guarda de memoria en la curvatura declara no evaluable una arista hub–hub (producto de soportes 10⁶). Solo afecta al 3-árbol s0, que ya rechazan (ii) y (iii), así que no cambia el veredicto.

## 1. Veredicto: **RC2-INVÁLIDA** (sensibilidad 12/19 = 0.632, especificidad 19/22 = 0.864)

| Familia | Grupo | Aceptadas | Anillos / n.º de escalas | Coherencia | κ mediana / f_neg | Rechazada por |
|---|---|---|---|---|---|---|
| Toro triangular 141² | V | 1/1 | CONEXO / 12 | C | 0.00 / 0.00 | — |
| **Panal hexagonal** | V | **0/1** | CONEXO / 12 | C | −0.33 / 1.00 | (i), (iv) |
| T⁴ 12⁴ | V | 1/1 | CONEXO / 3 | I | 0.00 / 0.00 | — |
| **RGG4 k16** | V | **0/3** | CONEXO / **1** (r_w = 5) | I | −0.11 / 0.52–0.56 | **(ii), (iv)** |
| RGG2 k20 | V | 3/3 | CONEXO / 8 | C | 0.08–0.10 / 0.14–0.18 | — |
| RGG3 k12, −15 % de aristas | V | 3/3 | CONEXO / 4 | C | −0.056…−0.029 / 0.39–0.42 | — (en el límite) |
| **T³ + 20 % de diagonales** | V | **0/3** | CONEXO / 4 | C | −0.10 / 0.45–0.52 | **(iv)** |
| Cubo abierto 27³ | V | 1/1 | CONEXO / 8 | C | 0.00 / 0.00 | — |
| RGG2 en cilindro | V | 3/3 | CONEXO / 12 | C | 0.07–0.10 / 0.16–0.22 | — |
| 3-árbol | X | 0/3 | SIN_VENTANA | SV | — | (ii), (iii), (iv) |
| Red apoloniana | X | 0/3 | CONEXO / 1 | INC | ≈ 0 / 0.37–0.42 | (ii), (iii) |
| Árbol de rejillas | X | 0/1 | RAMIFICADO / 12 | C | 0.00 / 0.01 | (ii) |
| BA m = 6 | X | 0/3 | SIN_VENTANA | SV | −0.67 / 1.00 | (ii), (iii), (iv) |
| Cuadrado + 1 % de atajos | X | 0/3 | CONEXO / 7 | INC | 0.00 / 0.08 | (iii) |
| **T³ + 0.1 % de atajos** | X | **3/3 (FP)** | CONEXO / 5 | I | 0.00 / 0.01–0.02 | — |
| Retazos b3 | X | 0/3 | CONEXO / 1 | INC | −0.15 / 0.58 | (ii), (iii), (iv) |
| Escalera, cilindro C×C₁₀, anillo k6 | U | 0/3 | DOS_EXTREMOS | C | — | (ii) |
| Heisenberg (E) | — | 0/1 | CONEXO / 4 | INC | −0.50 / 1.00 | (i), (iii), (iv) (informativo) |

**Predicciones (§II.4).**
- **Se cumplen:** 1–9, salvo la sensibilidad.
  - Panal rechazado por (i).
  - Árbol de rejillas RAMIFICADO.
  - 3-árbol y apoloniana rechazados.
  - BA SIN_VENTANA.
  - Cuadrado + 1 % rechazado.
  - T³ + 0.1 % falso positivo.
  - Las 1D, DOS_EXTREMOS.
  - Cubo abierto y cilindro RGG2 aceptados.
- **Fallan:**
  - La 10: el veredicto previsto era PARCIAL y salió INVÁLIDA.
  - El rechazo de **RGG4 k16** y de **T³ + diagonales** no se previó.

## 2. El hallazgo importante: **el juez de A no es ciego a d**

Los dos fallos no previstos no son ruido. Revelan dos **sesgos dimensionales** en instrumentos que se suponían ciegos a d.

### 2.1 Sesgo de ventana (nivel 2, argumento de escalado)

- Con N fijo, la ventana intermedia es r_w ≈ (N / 4c)^{1/d}. El número de escalas de anillo disponibles, entre 2 y r_w/2, **decrece con d**.
  - RGG4 k16 con N = 2·10⁴: r_w = 5, una sola escala.
  - T⁴ 12⁴: 3 escalas.
  - RGG3: 4–5 escalas.
  - RGG2: 8–12 escalas.
- La condición «≥ 2 escalas», añadida para excluir el 2-árbol, equivale a **una cota superior d ≤ d_max(N, c)**.
- **Toda condición «multiescala» a N fijo penaliza las d altas.** Para juzgar una d con s escalas hace falta N ≳ 4c(2(s+1))^d: el **coste de resolución crece exponencialmente con d**.

### 2.2 Sesgo de curvatura (nivel 2, dato + mecanismo)

- La curvatura de Ollivier perezosa de los RGG con k fijo **decrece con d**:
  - RGG2 k12–20: +0.08–0.12;
  - RGG3 k12: −0.03;
  - RGG4 k16: −0.11.
  - Los datos de RGG2/RGG3 vienen de L-CIC-0b; los de RGG4, de aquí.
- El mecanismo: al crecer d, los vecindarios de dos extremos de una arista se solapan menos. La «planitud» que mide K a escala 1 depende de la dimensión y de la discretización: el T³ con diagonales (k = 7.2) da −0.10, y el T³ puro da 0.
- La ventana de κ de RC-1, calibrada con datos de d ≤ 3, **es un filtro de dimensión y de discretización**.

### 2.3 Consecuencia

- **RC-1 y RC-2 tienen un sesgo hacia d baja.** Juzgada con ellas a N ≈ 2·10⁴, una dinámica futura que produjera geometría de d = 4 sería rechazada, y una de d = 2 o 3 aceptada.
- Usar este juez en el problema A **introduciría por la puerta de atrás la selección de dimensión** (problema B). Es exactamente lo que el orden lógico prohíbe.
- Esto también matiza resultados anteriores: la «independencia de la dimensión» de B-ter (L-CIC-0b) solo se comprobó con d = 2 y 3.

### 2.4 El falso positivo

- T³ + 0.1 % de atajos (≈ 60 atajos) pasa las cuatro condiciones.
- La distancia típica a un extremo de atajo es ≈ (N/120)^{1/3} ≈ 5.5 saltos, dentro de la ventana (r_w = 11). Aun así, ningún instrumento a escala intermedia lo detecta con este N.
- Es el mismo defecto de cruce lento ya visto con WS β = 0.001.

## 3. Dictamen del cerebro

1. **RC-2 queda INVÁLIDA**, como juez de A y como diagnóstico ciego a d.
   - Por regla prerregistrada (§II.3), **el problema A vuelve al Consejo antes de cualquier dinámica**.
   - No se reajusta ningún umbral.
2. **Restricción nueva para el programa** (mapa negativo, lección 17): un juez de no degeneración a N fijo no es ciego a d.
   - El número de escalas disponibles decrece con d.
   - La curvatura local depende de d y de la discretización.
   - Cualquier juez futuro de A debe cumplir dos cosas:
     - **(a)** declarar su **resolución dimensional** d_max(N);
     - **(b)** validarse con **d = 2, 3, 4 y varias discretizaciones** de cada una.
3. **Propuesta al Consejo, para elegir una:**
   - **(A1) Juez relativo a la resolución.** Las condiciones se expresan en escalas relativas: número de escalas ≥ s, con N creciendo con d en la validación. La curvatura se sustituye por un criterio **relativo a un nulo con el mismo grado y la misma discretización**, por ejemplo κ(G) frente a κ(G reconectado localmente) en lugar de una ventana absoluta. Requiere un prerregistro nuevo (RC-3) y un panel con d = 2..5.
   - **(A2) Juez solo de rechazo, sin ventana de curvatura.** Solo B-ter y coherencia, con resolución declarada. Acepta perder especificidad (expansores con triángulos, cruces lentos) a cambio de no sesgar d.
   - **(A3) Pausa de instrumentación.** Si cada nuevo juez descubre un sesgo nuevo, el problema puede ser que el reconocimiento de no degeneración a tamaño finito sea intrínsecamente dependiente de la escala. En ese caso, la búsqueda de dinámicas debería evaluarse por **tendencias con N** (N, 2N, 4N y 8N, como E4), no por umbrales a un N fijo.

   **Recomendación del cerebro: A3 + A1.** Un juez de A solo es ciego a d si es **de escalado** (compara N con kN) y **relativo** (compara con un nulo emparejado), nunca de umbral absoluto.
