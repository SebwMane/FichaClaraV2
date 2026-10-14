# F6 confirmatoria: resultados (intento 1 de 2, 2026-10-14)

- **Especificación:** `OMEGA_F6_CONFIRMACION_PREREGISTRO.md` (congelado en f7ce889; sello abierto después, hash verificado; aclaraciones A1 y A2). **Herramienta:** `tools/f6_confirmacion.py` (auditoría §25). **Corrida oficial:** commit 891755f, árbol limpio, 52 partes. La corrida se cortó una vez por el límite de 2 horas de un comando; la parte en curso no había guardado nada y se repitió completa sobre el mismo commit. **Datos:** `results/omega_f6_confirmacion/`. **Figura:** `docs/figuras/f6_confirmacion.png`.
- **Calibración:** b_plano = +0.0004; τ_b = 0.0039 (15 grupos válidos, ninguno no medible).

## 1. Veredicto oficial: **FALLA**, por una sola estructura

**C5, Minkowski con 13% de links borrados (negativa sellada), quedó clasificada como espacio-tiempo** (|b − b_plano|/τ_b = 0.002). Las otras 19 estructuras con voto quedaron bien clasificadas.

| Estructura | Papel | |b − b_plano|/τ_b | Clasificación | ¿Correcta? |
|---|---|---|---|---|
| Diamante plano, de Sitter, de Sitter fuerte, anti de Sitter, FRW radiación, FRW polvo | positivas (vistas) | 0.00 a 0.23 | espacio-tiempo | sí (6/6) |
| **C0 Kasner de vacío (no conformemente plano)** | positiva sellada | 0.07 | espacio-tiempo | **sí** |
| **C1 cono elíptico rotado** | positiva sellada | 0.16 | espacio-tiempo | **sí** (el riesgo de borde anotado no se materializó) |
| **C2 Minkowski con densidad ondulada** | positiva sellada | 0.02 | espacio-tiempo | **sí** |
| Conos cúbico, octaédrico, p = 6 | negativas (vistas) | 4.85, 4.37, 2.28 | rechazados | sí |
| Orden producto 4, percolación | negativas (vistas) | 6.2, 92 | rechazados | sí |
| **C3 cono cilíndrico (ε = 0.29)** | negativa sellada | 2.79 | rechazado | **sí** |
| **C4 esfera ∩ cubo (ε = 0.27)** | negativa sellada | 1.43 | rechazado | **sí** (el margen más estrecho) |
| **C5 links borrados 13%** | negativa sellada | **0.002** | **espacio-tiempo** | **NO** |
| **C6 p por eje (ε = 0.31)** | negativa sellada | 2.78 | rechazado | **sí** |
| Conos p = 1.5, p = 3, C7 | zona gris | 0.21, 0.78, 0.58 | espacio-tiempo | sin voto |

## 2. Predicciones informadas (sin voto)

- **Signo de la pendiente A:**
  - de Sitter fuerte: +0.033, IC [+0.025, +0.042]. Cumple.
  - anti de Sitter: −0.042, IC [−0.054, −0.030]. Cumple.
  - de Sitter normal: +0.009, pero su IC incluye 0. Solo cumple en el signo.
  - A del fuerte es mayor que A del normal: cumple.
- **Valores que debían ser ≈ 0:**
  - FRW de radiación: −0.002, IC [−0.011, +0.008]. Cumple; esta vez sin la pendiente negativa al límite que había dado la exploración.
  - **Kasner de vacío:** −0.003, IC [−0.016, +0.010]. Cumple la predicción de H_loc para vacío (sin término τ²). Es la primera prueba de la línea en un espacio-tiempo que no es conformemente plano.
- **C5:** las verificaciones obligatorias pasaron en las 20 realizaciones (orden parcial y cierre igual a C′); se borraron unos 22 100 links con N = 4000 y unos 66 500 con N = 8000.

## 3. Por qué pasó C5 (nivel 1 lo medido; nivel 2 la explicación; observado después del veredicto)

Δr de C5 por m (N = 4000 / N = 8000):

| m | 10 | 20 | 40 | 80 | 160 |
|---|---|---|---|---|---|
| N = 4000 | +0.000 | −0.002 | −0.007 | −0.008 | −0.009 |
| N = 8000 | −0.000 | −0.002 | −0.008 | −0.008 | −0.008 |

- **El defecto desaparece en los diamantes más chicos** (m = 10), porque un intervalo chico rara vez contiene un link borrado entre sus interiores. Por eso el ajuste lineal en √(m/N) le asigna intercepto ≈ 0 y una pendiente grande (A = −0.057, mayor en módulo que la de anti de Sitter): **el instrumento lo lee como curvatura negativa fuerte.**
- **Pero no depende de la densidad:** con el mismo m, N = 4000 y N = 8000 dan lo mismo. La curvatura sí depende (de Sitter fuerte y anti de Sitter se acercan a 0 al subir N, figura B). Es la firma que separó curvatura de forma en la exploración (`OMEGA_F6_EXPLORACION_RESULTADOS.md` §2), que **no votaba** en esta confirmación.
- **Pregunta de Ω que esto vuelve concreta** (OMEGA_MASTER §4 y asesor de ramas, ronda 2): un defecto microscópico que se borra en las escalas chicas pero no se «diluye» con la densidad, ¿es un impostor o un microestado de la misma geometría? Para el reconocedor, hoy es un impostor; para Ω es una pregunta abierta.

## 4. Lo que deja F6

- **Primer resultado de la línea con estructuras a ciegas y curvatura no conformemente plana:** el intercepto separó 9 de 9 positivas (incluidas tres selladas) de 9 de 10 negativas (incluidas tres selladas cerca de la frontera de la zona gris). La única fuga fue un impostor diseñado para imitar a la curvatura, y lo logró.
- **Presupuesto: queda un intento** (Rev. 11 §6). Un segundo intento debería volver al Consejo con una propuesta que use la prueba de la densidad (la desviación de la curvatura escala con m/N; un defecto o una forma, solo con m). Esa propuesta nace de ver C5: tendrá que preregistrarse antes, probarse con semillas nuevas y con estructuras selladas nuevas, porque C0 a C7 ya están vistas.
- **El conjunto de reserva** (`RESERVA.md`) **sigue sellado**: solo se abre tras un PASA.
- **Bordes y topología** (la pregunta original de F6) siguen abiertos y requieren N ≳ 30 000.
