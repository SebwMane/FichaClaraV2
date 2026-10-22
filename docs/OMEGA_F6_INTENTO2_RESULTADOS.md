# F6, intento 2 (último): resultados (2026-10-18)

- **Especificación:** `OMEGA_F6_INTENTO2_PREREGISTRO.md`: congelado en 9b75b21, hash del sello en c0709bd, sello abierto después; aclaraciones B1 a B3. **Herramienta:** `tools/f6_intento2.py` (auditoría §26). **Corrida oficial:** commit 1619162, árbol limpio, 158 partes y `final`.
- **Incidentes de la corrida:** se hizo en 4 tandas. Un reinicio del contenedor cortó una parte en curso, que no había guardado nada; se repitió completa sobre el mismo commit. Antes de reanudar se verificaron las 98 partes guardadas: todas legibles y con el mismo commit.
- **Datos:** `results/omega_f6_intento2/`. **Figura:** `docs/figuras/f6_intento2.png`.
- **Calibración:** b_plano = +0.00001; τ_b = 0.0027; τ_K = 11.5. Los 15 grupos son válidos.

## 1. Veredicto oficial: **FALLA**

**Cuatro negativas selladas quedaron clasificadas como espacio-tiempo:** D5 (10% de elementos gemelos) y los tres conos nuevos cerca de la frontera, D6 (tapa plana), D7 (huso) y D8 (suma de normas). Todas las positivas, incluidas las cuatro selladas, quedaron bien.

| Estructura | Papel | \|b − b_plano\|/τ_b | K/τ_K | Clasificación | ¿Correcta? |
|---|---|---|---|---|---|
| 6 geometrías vistas (plano, de Sitter, de Sitter fuerte, anti de Sitter, FRW radiación, FRW polvo) | positivas | 0.05 a 0.84 | 0.15 a 0.47 | espacio-tiempo | sí (6/6) |
| C0 Kasner, C1 cono elíptico, C2 densidad ondulada | positivas (vistas) | 0.05 a 0.09 | 0.17 a 0.79 | espacio-tiempo | sí (3/3) |
| **D0 Kasner u = 2, D1 cizalla, D2 FRW exponencial, D3 FRW rígida** | positivas selladas | 0.01 a 0.16 | 0.11 a 0.46 | espacio-tiempo | **sí (4/4)** |
| Conos cúbico, octaédrico, p = 6; producto 4; percolación | negativas (vistas) | 2.8 a 134 | | rechazadas | sí (5/5) |
| C3, C4, C6 (conos) | negativas (vistas) | 2.0 a 4.3 | | rechazadas | sí |
| **C5 links borrados** | negativa (vista) | 0.56 | **3.04** | **rechazada por el colapso** | sí |
| **D4 relaciones agregadas** | negativa sellada | 49 | 109 | rechazada | sí |
| **D5 gemelos 10%** | negativa sellada | 0.83 | 0.35 | espacio-tiempo | **NO** |
| **D6 tapa plana (ε = 0.25)** | negativa sellada | 0.19 | 0.15 | espacio-tiempo | **NO** |
| **D7 huso (ε = 0.31)** | negativa sellada | 0.66 | 0.57 | espacio-tiempo | **NO** |
| **D8 suma de normas (ε = 0.23)** | negativa sellada | 0.71 | 0.16 | espacio-tiempo | **NO** |

Las verificaciones obligatorias de D4 y D5 pasaron en las 60 realizaciones de cada una.

## 2. Lectura (nivel 1 lo medido; nivel 2 lo marcado)

- **La prueba de colapso hizo lo que debía:** atrapó a C5 (K/τ_K = 3.0) y a D4, sin rechazar ninguna positiva. Esto incluye las dos FRW nuevas: el riesgo declarado de FRW no se materializó.
- **El instrumento reconoce bien los espacio-tiempos:** 13 de 13 positivas, entre ellas dos Kasner de vacío, una cizalla que mezcla tiempo y espacio y tres FRW con curvatura variable, dos de ellas a ciegas.
- **Lo que falló es la resolución frente a conos «suaves en partes».** D6 tiene una sola cara plana, D7 dos puntas y D8 pliegues suaves. Sus interceptos (0.0005 a 0.0019) quedan dentro del ruido de la calibración (τ_b = 0.0027), mientras C4, con ε parecido (0.27), dio 0.0054. **La medida ε (desviación máxima del cono redondo) no predice qué tan visible es un cono** (nivel 2): importa cuánta superficie del cono se aparta de la esfera, no el máximo. Una cara plana o dos puntas cambian poco el volumen de los intervalos.
- **D5 (gemelos):** 10% de elementos duplicados cambia r en apenas −0.002 (cada par de gemelos es un par no relacionado más, diluido entre muchos pares) y el defecto escala igual a toda densidad, así que el colapso no lo ve. El intercepto queda justo bajo el umbral (0.83).

## 3. Cierre de la línea

- **Este era el último intento** (R13-4). La línea del perfil local (F6 con intercepto y colapso) **vuelve al Consejo para una revisión completa**; no se prueba otra variante sin ella.
- **El conjunto de reserva** (`RESERVA.md`) **no se abre:** solo se abría tras un PASA.
- **Lo que queda establecido (nivel 1, en este alcance: regiones con forma de intervalo, 4D, N ≤ 8000):**
  1. el instrumento acepta todos los espacio-tiempos probados (planos, curvos, no conformemente planos, con curvatura variable, con transformaciones lineales);
  2. rechaza órdenes producto, KR, percolación, capas, conos con aristas marcadas y los defectos discretos que cambian la escala (C5, D4);
  3. no distingue conos con deformaciones «parciales» de hasta ε ≈ 0.31, ni una duplicación del 10% de los elementos.
- **Preguntas para la revisión del Consejo:**
  - ¿El límite de resolución se mueve con N (más densidad, menos ruido en b)?
  - ¿Hay un observable mejor que r para la forma del cono (por ejemplo, la distribución de tamaños de intervalos dentro de cada diamante chico)?
  - ¿Los gemelos son un impostor relevante para Ω, o una «degeneración» aceptable del sustrato (como preguntaba OMEGA_MASTER §4 para los links borrados)?
