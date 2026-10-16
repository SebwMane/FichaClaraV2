# Ω — L-A-0: problema A (exclusión de degenerados, ciega a d). Paisaje de mecanismos y juez RC-2 — Prerregistro

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-a0`.
- **Base congelada:** `claude/omega-arq-congelado`.
- **Mandato:** el usuario aprueba el orden lógico de L-ARQ-0 §6.4: **el problema A antes que el B**.
- **Estado:** este documento se commitea antes de escribir el código de RC-2 y antes de generar ningún grafo del panel nuevo.

## Parte I — Enunciado y paisaje (análisis, sin código)

### I.1 Problema A

Se busca un mecanismo local, sin coordenadas, sin d, sin escala impuesta y sin alfabeto que fije d (E6), cuyo resultado esté en la clase:

> **G≥2 (geometría gruesa no degenerada):** local, coherente, con un extremo a toda escala y curvatura no degenerada, para alguna d ≥ 2 **no especificada**.

Por CH-T6 y Stallings, en estructuras homogéneas eso implica una dimensión entera ≥ 2, **sin decir cuál**. Así que el problema A no prejuzga el B.

### I.2 Clases de mecanismos ya cerradas (evidencia del mapa negativo y paralelos externos)

| Clase | Destino | Fuente |
|---|---|---|
| Energía de equilibrio sobre conteos locales (k, c, q) | No puede preferir geometría | C0-T1/T2 (nivel 1); C1 |
| Entropía (Θ > 0) sobre esos mismos ensembles | Gana el fondo entrópico | Fase 2. Paralelo: órdenes de Kleitman–Rothschild (nivel 3) |
| Ensembles de equilibrio de variedades (condición de enlace local) | Arrugado o polímero ramificado; transición de primer orden | Paralelo DT euclídeo (nivel 3, L-ARQ-0 §5) |
| Crecimiento con anclaje local (uniforme, por frontera, por capas) | Mundo pequeño, 1D, hiperbólico o polímero ramificado | P3 T1–T3 (nivel 1) |
| Cercanía por difusión | Premisa falsa | Ω-D |
| Planitud de Ollivier como objetivo | No excluye árboles ni C0 | P3-D |
| Competencia escalar entre presiones | Dial | L-ARQ-T1 (nivel 1) |

### I.3 Clases todavía no discriminadas

Ninguna se ejecuta en esta fase.

- **R1. Reescritura local fuera del equilibrio con una cantidad conservada.** No es energía ni crecimiento.
- **R2. Restricciones (tipo β) que no son de equilibrio.** Ejemplo: una condición de enlace local ciega a d («el enlace de cada vértice es una esfera combinatoria de alguna dimensión»). Por conexidad, la dimensión del enlace es constante, así que la condición es ciega a d y excluye cliques y expansores. **Pero** el paralelo DT indica que sola no basta: los polímeros ramificados de DT son variedades. Predicción de nivel 2: R2 sola falla con los árboles de bloques.
- **R3. Estructura adicional** (orden u foliación). El paralelo CDT funciona, pero impone la foliación y la dimensión del bloque (E6). Haría falta una versión ciega a d.

### I.4 Requisito previo

Cualquier búsqueda en R1–R3 necesita un **juez de A validado**:
- RC-1 resultó INVÁLIDA (L-CIC-0b).
- La regla «anillos + coherencia + K» se observó *post hoc* sobre el panel de L-CIC-0b, así que solo puede validarse fuera de muestra.

Esa es la única parte ejecutable de esta fase (Parte II).

## Parte II — RC-2: juez de A, validación fuera de muestra

### II.1 Regla congelada (antes de ver el panel nuevo)

Instrumentos congelados, con los umbrales de L-CIC-0b, las mismas claves de RNG por instrumento y el mismo número de muestras:

| Condición | Contenido | Herramienta |
|---|---|---|
| (i) Localidad | Fracción de aristas en un triángulo o un 4-ciclo ≥ 0.5 | `cic_l0b.locality_fraction` |
| (ii) Un extremo multiescala | `annulus_status` = CONEXO **y** ≥ 2 escalas en `annulus["scales"]` | `annulus.py` |
| (iii) Coherencia no rota | `coherence_status` ∈ {COHERENTE, INTERMEDIO} | `coherence.py` |
| (iv) Curvatura no degenerada | Mediana de κ en [−0.0566, 0.2700] y f_neg ≤ 0.4283 | `cic_l0b.curvature_summary` |

- **A-NO-DEGENERADO** si se cumplen (i)–(iv).
- **No hay ninguna condición de dimensión ni de Nivel II.** RC-2 es ciega a d por construcción.

**Origen *post hoc* declarado.**
- (ii) añade «≥ 2 escalas» por el fallo del 2-árbol.
- (iii) relaja COHERENTE para admitir INTERMEDIO, por T³ y RGG3 k8.
- (iv) se mantiene como en RC-1.
- Sobre el panel antiguo, RC-2 separaría por construcción. **Por eso solo cuenta el panel nuevo.**

### II.2 Panel nuevo

- Ninguna familia se usó para fijar RC-1 ni RC-2.
- N ≈ 2·10⁴. Clave (20261015, familia, semilla). Semillas 0–2 en las familias aleatorias y 0 en las deterministas.

| Grupo | Familias |
|---|---|
| **V** (no degeneradas, d ≥ 2) | Toro triangular 2D 141² (k = 6) · panal hexagonal 2D toroidal (k = 3, sin triángulos ni 4-ciclos) · T⁴ 12⁴ · RGG4 k16 · RGG2 k20 · RGG3 k12 con el 15 % de aristas eliminadas al azar (componente gigante) · T³ 27³ con un 20 % de diagonales de cara añadidas al azar · cubo abierto 27³ (sin periodicidad) · RGG2 k12 en cilindro (periódico solo en x, con razón de aspecto 1) |
| **X** (degeneradas) | 3-árbol aleatorio · red apoloniana (triangulación apilada) · árbol de rejillas (árbol binario de parches 10×10 pegados por una arista) · Barabási–Albert m = 6 · cuadrado 141² + 1 % de atajos · T³ + 0.1 % de atajos · retazos con blocks = 3 |
| **U** (1D) | Escalera C₅₀₀₀ × P₄ · cilindro C₂₀₀₀ × C₁₀ · anillo k6 |
| **E** (informativo, no entra en el veredicto) | Heisenberg (no local por diseño) |

### II.3 Criterios

- Sensibilidad = fracción de instancias V → A-NO-DEGENERADO.
- Especificidad = fracción de instancias X + U → rechazadas.

| Veredicto | Condición |
|---|---|
| **RC2-VÁLIDA** | Sensibilidad ≥ 0.8 y especificidad ≥ 0.95 |
| **RC2-PARCIAL** | Sensibilidad ≥ 0.7 y especificidad ≥ 0.9 |
| **RC2-INVÁLIDA** | En otro caso |

Se informa qué condición rechaza cada instancia. **Ningún umbral se toca después de ver resultados.** Si RC-2 falla, la Parte I se queda sin juez y el problema A vuelve al Consejo antes de cualquier dinámica.

### II.4 Predicciones del cerebro (congeladas)

1. **Panal hexagonal:** falso negativo por (i). La localidad por triángulos o 4-ciclos es ciega a las geometrías de cintura 6. Es un límite conocido del instrumento.
2. **Árbol de rejillas:** rechazado por (ii) RAMIFICADO. Las escalas llegan más allá del tamaño de parche.
3. **3-árbol y red apoloniana:** rechazados (≤ 1 escala o SIN_VENTANA).
4. **BA:** SIN_VENTANA.
5. **Cuadrado + 1 % de atajos:** rechazado (coherencia o anillos).
6. **T³ + 0.1 % de atajos:** riesgo de falso positivo (cruce lento).
7. **Escalera y cilindro:** DOS_EXTREMOS.
8. **Anillo k6:** DOS_EXTREMOS.
9. **Cubo abierto y cilindro RGG2:** aceptados (RC-2 no usa Nivel II).
10. **Veredicto previsto: RC2-PARCIAL**, con sensibilidad ≈ 0.85–0.9 (pierde el panal) y especificidad ≈ 0.93–1.0.
