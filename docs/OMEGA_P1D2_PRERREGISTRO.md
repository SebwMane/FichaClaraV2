# Ω — P1-D.2: coherencia multiescala y estabilidad dimensional — Prerregistro

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-p1d2`.
- **Base congelada:** `claude/omega-p1d-congelado` (P1-D PARCIAL).
- **Mandato del Consejo:** experimento diagnóstico pequeño, con controles nuevos y predicción congelada.
  - **No** es un certificado y **no** se llama «firma geométrica».
  - **No** fija D = 3: puede aparecer D = 1, 2, 3, 4…
  - **No** modifica C0 ni entra en la energía.
- **Estado:** este documento se commitea antes de escribir el código del panel y antes de generar ningún grafo.

## 0. Los tres niveles (Consejo)

| Nivel | Pregunta | Herramienta | En este experimento |
|---|---|---|---|
| **I** | ¿Los vecindarios se homogeneizan a escala intermedia? | ρ = CV(r_max)/CV(2) (P1-D) | Se mide |
| **II** | ¿Existe una dimensión efectiva estable? | meseta de D_B(r) | Se mide |
| **III** | ¿Se comporta como una geometría local (variedad)? | Certificado Ω-1.1 y sus banderas | Solo como comprobación, para las familias que pasen I y II (§5) |

I → II → III: ninguno sustituye al siguiente.

## 1. Definición de Nivel II (nueva; es lo único que se añade)

Sobre el perfil de P1-D (componente gigante, ventana N/4, W = {r ≥ 2 : mediana |B_r| ≤ N_g/4}):

- D_B(r) = ln(m(r+1)/m(r)) / ln((r+1)/r), definido para r ∈ [2, max W − 1] (se exige que r y r+1 estén en W).
- **Meseta:** hay una racha de ≥ 3 valores consecutivos de D_B con (máx − mín)/media ≤ δ.
  - Valor de la meseta: D_plat = media de D_B en la racha más larga (si empatan, la de mayor r).
- **Estados:** `PLATEAU`; `SIN_PLATEAU` (hay ≥ 3 valores de D_B pero ninguna racha cumple); `SIN_VENTANA_D` (menos de 3 valores de D_B).
- **δ primario = 0.15** (relativo). Sensibilidad: δ = 0.10 y δ = 0.25.

**Declaración de procedencia.** El umbral δ se eligió **después** de ver las trazas de D_B de P1-D:
- RGG3 a N = 4096 → rango relativo ≈ 0.04;
- caveman, con oscilaciones de 0.5 a 1.7 → rango relativo > 1.

Por eso es *post hoc*. Se informa el veredicto con los tres δ y no se afirma nada que dependa del valor elegido.

**Nivel I en este experimento:** pasa si el estado es `HOMOGENEIZA` o `CV_CERO`. CV_CERO son los retículos, perfectamente homogéneos por simetría.

**Pasa conjunto:** I y II.

## 2. Panel (todo con N = 4096; clave PCG64 (20261008, familia, semilla))

| Clase | Familia | Semillas |
|---|---|---|
| **G** (geometrías) | RGG2 k12, RGG3 k12, RGG3 k8 | 0–2 |
| G | Retículos: anillo k12 (1D), toro cuadrado 64×64 (2D), toro triangular 64×64, T³ 16³ (3D) | 1 |
| **R** (regulares, geometría gruesa pero localmente no variedad) | Caveman K4, K8, K16 (anillos de cliques) | 1 |
| R | Retículo de cliques 3D: toro 8³ de sitios, cada sitio es una K8; la dirección j (0…5) enlaza el miembro j de un sitio con el miembro j^1 del vecino en esa dirección | 1 |
| **L** (locales no geométricas) | Retazos-8, retazos-27, retazos-64 (RGG3 k12 con reemparejado de aristas entre bloques; P1-D) | 0–2 |
| L | RGG3 k12 con 1 % de atajos | 0–2 |
| L | Watts–Strogatz k12 con β = 0.003, 0.01 y 0.03 | 0–2 |
| **NL** (aleatorias) | ER k12, RR k12 | 1 |
| *Informativo* | Finales C0-R, F1, N = 729 (27 grafos de `runs/c0_f1/out/`) | 0–2 |

Total del panel de decisión: 40 grafos.

**Aviso sobre la clase R.** Un anillo de cliques es cuasi-isométrico a Z, y un retículo de cliques 3D es cuasi-isométrico a Z³. A escala gruesa son geometrías; localmente no son variedades, porque los vecindarios son cliques. El Consejo los listó como «regulares no geométricas»; aquí se mantienen **aparte** de las L y NL precisamente porque su estatus geométrico depende de la definición que se adopte (cuasi-isometría o vecindarios de variedad). Esa ambigüedad es lo que el Nivel III debe resolver.

## 3. Predicciones del cerebro (congeladas antes de generar ningún grafo)

1. **G.** RGG2 y RGG3 pasan I y II en ≥ 2/3 de las semillas. Los retículos son CV_CERO con meseta cerca de 1, 2 y 3 (tolerancia 0.25).
2. **R, caveman (K4, K8, K16).** Pasan el Nivel I y **fallan el Nivel II**, por la oscilación de período 2 de D_B, que persiste con el tamaño de clique.
3. **R, retículo de cliques 3D.** Aquí me juego el diagnóstico. Predigo que **sí pasa I y II** (con una meseta cercana a 3), porque es cuasi-isométrico a Z³.
   - **Nivel III esperado:** F9 por las cliques locales.
   - **Consecuencia:** I y II solos no separan una geometría gruesa de una geometría con vecindarios de variedad.
4. **L.** Retazos, atajos y WS no pasan I (ρ ≥ 1), como en P1-D, y por tanto no pasan el conjunto.
5. **NL.** ER y RR no tienen ventana (SIN_VENTANA).
6. **Nivel II aislado.** Añade poco sobre el Nivel I. Solo separa la familia caveman.

**Predicción global:** el panel de decisión da **SEPARA-PARCIAL**. Con una geometría gruesa de la clase R pasando I y II (el retículo de cliques), el Nivel III hace el trabajo restante.

## 4. Criterio

Sobre el panel de decisión de 40 grafos:

- **Sensibilidad** = fracción de G que pasa el conjunto (I y II).
- **Especificidad** = fracción de L ∪ NL que **no** pasa el conjunto.
  - Un grafo sin ventana (SIN_VENTANA) cuenta como «no pasa»: no puede reproducir la firma. Se informa aparte cuántos grafos no son evaluables.

| Veredicto | Condición |
|---|---|
| **SEPARA** | sensibilidad ≥ 0.9 **y** especificidad ≥ 0.9 **y** ninguna familia R pasa el conjunto |
| **SEPARA-PARCIAL** | sensibilidad ≥ 0.9 y especificidad ≥ 0.9, pero alguna familia R pasa el conjunto |
| **NO SEPARA** | sensibilidad < 0.7 o especificidad < 0.7 |
| **INDETERMINADO** | cualquier otro caso, o menos de 5 grafos evaluables en L ∪ NL |

El veredicto se calcula con δ = 0.15 y se informa también para δ = 0.10 y 0.25.

## 5. Nivel III (comprobación, solo para familias R que pasen el conjunto)

- Para cada familia R que pase I y II, se construye su análogo con N ≈ 729 y se corre `certificate_report` con 3 claves.
  - Caveman K8 con N = 728.
  - Retículo de cliques 3D con K6 sobre un toro 5³, N = 750.
- Referencias al mismo N: RGG3 k12 (729), RGG2 k12 (729) y anillo k12 (729).

**Limitación declarada (medida antes de prerregistrar, sobre una semilla):** el certificado, tal como está calibrado, valida el RGG3 a N = 729 pero **también rechaza geometrías reales de otra dimensión**:

| Grafo (N ≈ 729) | Códigos |
|---|---|
| RGG3 k12 | ninguno (pasa) |
| RGG2 k12 | F9 (`manifold_proxy_ok`) |
| Anillo k12 | F5, F9 |
| Caveman K8 | F3, F4, F9 |

El Nivel III solo discrimina a partir de la clase 3. En D = 1 y D = 2 los umbrales del certificado rechazan incluso geometrías auténticas, de modo que su fallo ahí **no** es evidencia contra la geometría. No se modifican sus umbrales.

## 6. Lo que no se hace

- Ni certificado nuevo ni término de energía.
- Ni requisito D_B ≈ 3.
- Ni cambios en el certificado Ω-1.1 ni en código congelado.
- Los finales C0 de N = 729 son solo informativos: su ventana a N = 729 no admite el Nivel II.
