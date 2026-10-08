# Ω — P1-D: preprueba diagnóstica del crecimiento de bolas — Prerregistro

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-p1-diagnostico`.
- **Base congelada:** `claude/omega-fase2-congelado`.
- **Mandato del Consejo:** P1 solo como diagnóstico, **no** como energía. Pregunta: ¿la firma «CV(|B_r|) decreciente a escala intermedia» distingue de forma sistemática la geometría de la localidad no geométrica? Si no la distingue, P1 muere aquí, barato.
- **Estado:** este documento se commitea antes de escribir el código y antes de generar ningún grafo del panel.

## 1. Origen de la hipótesis y conjunto de descubrimiento

La firma se observó en `results/c1_ball_growth/` (Fase 2 §2.2) sobre:
- 3 RGG3 k12 con N = 729, que se autopromedian (el CV baja);
- 27 finales C0-R con N = 729, en los que el CV sube.

Ese es el **conjunto de descubrimiento**. **No se usa para decidir**; solo se informa.

## 2. Observable (`omega/diagnostics/ball_growth.py`, observador permanente)

Sobre la componente gigante del soporte binario (los grafos de C0 usan el soporte fuerte, W > 0.1·max W), con N_g su número de nodos:

- |B_i(r)| = #{j : d(i, j) ≤ r}, para r = 1 … r_cap, donde r_cap es el primer r con mediana de |B_r| ≥ N_g/2 (como máximo 15);
- m(r) = media de |B_r|; CV(r) = desviación típica / media entre nodos;
- D_B(r) = ln(m(r+1)/m(r)) / ln((r+1)/r), derivada logarítmica discreta. Es descriptiva;
- **ventana intermedia:** W = {r ≥ 2 : mediana de |B_r| ≤ N_g/4}. Excluye r = 1, que es el grado, y la saturación por tamaño finito.

**Estadístico primario:** ρ = CV(r_max) / CV(2), con r_max = max W.

**Estados:**
- **CV≡0:** CV(2) < 1e-9 (homogéneo por simetría, como los retículos). Se informa aparte.
- **SIN_VENTANA:** |W| < 2.
- **HOMOGENEIZA:** ρ < 1. Es la firma de geometría candidata.
- **NO_HOMOGENEIZA:** ρ ≥ 1.

**Secundario (informe):**
- ρ con cota N_g/2 en lugar de N_g/4;
- fracción de pasos de W en los que CV baja;
- rango de D_B en W.

## 3. Panel de validación (todo nuevo)

- RNG: PCG64 con clave (20261007, familia, semilla).
- N = 729 salvo indicación; semillas 0–2 en las familias aleatorias.

| Etiqueta | Familia | Descripción |
|---|---|---|
| **G** | RGG3 k8, k12, k16 | Toro 3D, binario |
| G | RGG3 k12, N = 1728 | Ventana más amplia |
| G | RGG2 k12 | Toro 2D |
| G | RGG4 k12 | Toro 4D |
| G | RGG3 diluido | RGG3 k16 con el 25 % de aristas borradas al azar (k ≈ 12) |
| G (CV≡0 esperado, aparte) | Retículos | T³ 9³; toro cuadrado 27×27; toro triangular 27×27; anillo k12 (N = 729) |
| **L** | Retazos-27 | RGG3 k12, toro dividido en 3×3×3 bloques. Los extremos de las aristas que cruzan bloques se reemparejan al azar entre sí (preserva grados y la geometría interna de cada bloque). Geometría local, pegado global incoherente |
| L | Retazos-8 | Igual, con 2×2×2 bloques |
| L | Retazos-64, N = 1728 | 4×4×4 bloques |
| L | RGG3 + atajos | RGG3 k12 con el 10 % de las aristas reconectadas: un extremo se cambia por un nodo al azar |
| L | Watts–Strogatz | Anillo k12 con β = 0.05 y β = 0.2 |
| L | Caveman | Caveman conectado K8, N = 728 |
| **NL** | ER k12, RR k12 | Controles no locales |

**Panel secundario C0 (N = 343, fuera del conjunto de descubrimiento):**
- L: todos los finales DISPERSO-LOCAL de R6 con N = 343, semilla 0, inicios U, E y R (familias de cadenas, expansor decorado y R-3D);
- G de referencia al mismo N: RGG3 k12 y RGG2 k12 con N = 343, semillas 0–2.

## 4. Criterio (sobre el panel principal; se excluyen CV≡0 y SIN_VENTANA, que se informan)

- **Sensibilidad** = fracción de grafos G con HOMOGENEIZA.
- **Especificidad** = fracción de grafos L ∪ NL con NO_HOMOGENEIZA.

| Veredicto | Condición |
|---|---|
| **P1-D SOBREVIVE** | sensibilidad ≥ 0.9 **y** especificidad ≥ 0.9 |
| **P1-D MUERE** | sensibilidad < 0.7 **o** especificidad < 0.7 |
| **P1-D PARCIAL** | cualquier otro caso. Se informa por familia qué detecta y qué no |

El panel secundario C0 se evalúa igual y se informa por separado. Sirve para saber si la firma detecta al menos todas las familias C0 y no solo la R.

**Si algún grupo (G o L ∪ NL) queda con menos de 5 grafos evaluables**, el veredicto es INDETERMINADO.

## 5. Predicción del cerebro (registrada antes de correr)

- **Sensibilidad alta (≥ 0.9).** Las RGG en cualquier dimensión se autopromedian.
- **Especificidad baja.** Predigo que **la firma no es específica de la geometría**. El autopromediado es una consecuencia genérica de la mezcla: cuando una bola cubre muchas regiones independientes, su tamaño se concentra (ley de los grandes números). Por eso los retazos, los atajos y Watts–Strogatz homogeneizarán. Es posible que algunos queden SIN_VENTANA porque crecen rápido.
- **Predicción de veredicto: PARCIAL o MUERE por especificidad.** El CV creciente sería un detector específico de la heterogeneidad tipo C0 (complejos de cliques), no un certificado de geometría.

## 6. Lo que no se hace

- Ninguna funcional ni dinámica nueva. El observador no entra en el certificado.
- No se cambia ningún umbral existente.

## 7. Enmiendas

### P1-A1: panel de extensión con N = 4096 para resolver el INDETERMINADO (registrada tras ver el resultado oficial)

**Resultado oficial (sin cambios):**
- ventana N/4: principal INDETERMINADO. La sensibilidad es 15/15, pero solo hay 4 grafos L/NL evaluables; retazos, atajos, ER, RR y WS β = 0.2 quedan SIN_VENTANA a N = 729 porque crecen demasiado rápido;
- ventana N/2: MUERE, con especificidad 0/11.

**Motivo.** Los controles más duros (geometría local pegada de forma incoherente) no tienen ventana intermedia con N = 729. Con N = 4096 y bloques grandes, sí la tienen.

**Panel de extensión** (N = 4096; semillas 0–2 salvo indicación; clave PCG64 (20261007, 30 + i, s)):

| Etiqueta | Familia |
|---|---|
| G | RGG3 k12, RGG2 k12 |
| L | Retazos-8 (bloques de ~512 nodos), retazos-27 (~152), retazos-64 (~64) |
| L | RGG3 con 1 % de atajos |
| L | Watts–Strogatz anillo k12, β = 0.01 |
| L | Caveman K8 (N = 4096, una semilla) |
| NL | ER k12 (una semilla) |

**Criterio:** el mismo de §4, ventana N/4, aplicado a la **unión** del panel principal y la extensión. El veredicto de la unión sustituye al INDETERMINADO. Si sigue habiendo menos de 5 grafos evaluables en algún grupo, P1-D queda **INDETERMINADO definitivo** y se informa así.

**Predicción del cerebro:** los retazos-8 y retazos-27 HOMOGENEIZAN, porque la mezcla entre bloques autopromedia las bolas. La especificidad cae por debajo de 0.9 y el veredicto será **PARCIAL o MUERE**.
