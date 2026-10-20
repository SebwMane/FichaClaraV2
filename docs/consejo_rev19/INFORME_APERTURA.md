# Consejo Ω, sesión 8: informe de apertura (R14-B v2 y R14-C v2)

- **Fecha:** 2026-10-20. **Preparado por:** el orquestador, a pedido del usuario («quiero que el consejo revise estos resultados antes de decidir cuál preparar»).
- **Mandato del Consejo:** revisar R14-B v2 y R14-C v2 y recomendar cuál de las tres pruebas confirmatorias candidatas preparar (o ninguna). La decisión es del usuario.
- **Fuentes en el repositorio:**
  - informes: `docs/OMEGA_R14B_V2_RESULTADOS.md` y `docs/OMEGA_R14C_V2_RESULTADOS.md`;
  - preregistros: `docs/OMEGA_R14B_V2_PREREGISTRO.md` y `docs/OMEGA_R14C_V2_PREREGISTRO.md` (con sus aclaraciones A1);
  - datos: `results/omega_r14b/` y `results/omega_r14c/`;
  - figuras: `docs/figuras/r14b_v2.png` y `docs/figuras/r14c_v2.png`;
  - tabla descriptiva por estructura (calculada después de ver los datos, sin voto): `docs/consejo_rev19/descriptivo_por_estructura_r14c.json`.
- **Commits de ejecución:** R14-B en 73e4925 y R14-C en a0ae61f, ambos con el árbol limpio y un solo commit por corrida. Informes en c6d79b6 y e0a12fc.

## 0. Contexto del programa

- **Prioridad firmada (R14-1):** emergencia efectiva (H_E); compatibilidad microscópica (H_L) como hipótesis separada.
- **H_sel (selección dinámica):** indeterminada, no refutada. La afirmación de mezcla con N de 64 a 128 quedó cerrada (Rev. 18) tras ES-0, ES-C0, ES-C0b y R18-1.
- **Por qué R14-B y R14-C:** prueban órdenes **dados**, sin muestrear ninguna dinámica. Son exploratorias y **sin voto**, y todas sus estructuras ya se habían visto en F6 y F7.

## 1. R14-B v2: respuesta a la dilución aleatoria Φ_q

**Pregunta:** si se conserva cada elemento con probabilidad q, ¿los observables de un orden entran en la banda de una geometría (lo que H_E esperaría), o ya estaban en banda desde el principio (H_L)?

**Diseño:**
- N₀ = 8000 elementos y q ∈ {1, 1/2, 1/4, 1/8}, es decir, N efectivo de 8000, 4000, 2000 y 1000.
- 10 realizaciones por estructura.
- Distancia de Mahalanobis con 14 observables.
- Referencias propias para cada una de las 8 geometrías (60 sprinklings por N efectivo).
- Banda calibrada con 40 referencias de ajuste y 20 de validación.

### 1.1 Datos (mediana de D/radio; la banda está en 1)

| Estructura | q = 1 | 1/2 | 1/4 | 1/8 | Etiqueta | Caída de q = 1 a 1/8 |
|---|---|---|---|---|---|---|
| plana (propia) | 0.91 | 0.75 | 0.58 | 0.73 | converge a banda | |
| de Sitter (propia) | 0.62 | 0.69 | 0.85 | 1.03 | en banda microscópica | |
| de Sitter fuerte (propia) | 0.82 | 0.89 | 0.74 | 0.81 | converge a banda | |
| anti de Sitter (propia) | 0.76 | 0.62 | 0.85 | 0.61 | en banda microscópica | |
| FRW polvo (propia) | 0.86 | 0.68 | 0.77 | 0.73 | en banda microscópica | |
| D2 FRW exponencial (propia) | 0.82 | 0.67 | 0.83 | 0.82 | en banda microscópica | |
| C0 Kasner (propia) | 0.84 | 0.76 | 0.69 | 0.74 | en banda microscópica | |
| D0 Kasner (propia) | 0.65 | 0.69 | 0.68 | 0.88 | en banda microscópica | |
| C5, 13% de links borrados | 4.33 | 2.12 | 0.94 | 0.95 | indeterminada | 4.6× |
| D4, 10% de relaciones agregadas | 30.2 | 10.3 | 4.29 | 2.27 | no converge | 13× |
| **D5, 10% de gemelos** | **5.40** | **1.71** | **0.66** | **0.74** | **converge a banda** | **7.3×** |
| D6, cono con tapa plana | 1.05 | 0.64 | 0.63 | 0.80 | converge a banda | 1.3× |
| D7, cono en huso | 1.91 | 1.02 | 0.65 | 0.74 | converge a banda | 2.6× |
| D8, cono suma de normas | 1.99 | 1.09 | 0.66 | 0.75 | converge a banda | 2.7× |
| C4, cono esfera y cubo | 5.22 | 2.12 | 0.85 | 0.89 | indeterminada | 5.9× |
| cono octaédrico | 13.9 | 5.08 | 2.15 | 1.48 | no converge | 9.4× |
| producto 4 | 24.1 | 10.3 | 5.12 | 3.65 | no converge | 6.6× |

Las 8 geometrías se comparan con su propia referencia; las demás estructuras, con la plana.

### 1.2 Predicciones congeladas

| Predicción | Resultado |
|---|---|
| Las geometrías quedan en banda en todo q | Sí en lo esencial: ninguna queda fuera de banda en ningún q. En 4 de 8 algún q queda ambiguo, porque el IC toca el borde. |
| D5 converge a banda | Se cumple. |
| Los conos no convergen | **Falla en 4 de 5.** Solo el octaédrico no converge. |

### 1.3 Hallazgo de método

- **Diluir confunde convergencia con pérdida de resolución.** Al pasar de N = 8000 a 1000, el ruido de la referencia crece como 1/√N. Un defecto de tamaño fijo, medido en desviaciones estándar, se achica unas √8 ≈ 2.8 veces solo por eso.
- **D7 y D8 caen 2.6 y 2.7 veces,** justo lo que daría el ruido. Encajan con un defecto fijo que deja de verse, no con un orden que se vuelve geométrico.
- **D5 cae 7.3 veces,** bastante más que el ruido, y hay un mecanismo que lo explica: un par de gemelos sobrevive a la dilución con probabilidad q².
- **Cautela:** D/radio combina 14 componentes, así que la escala con √N es aproximada.
- **Calibración de bandas:** con p = 14 y n = 40, el percentil 95 de las referencias de ajuste cubre solo entre el 50% y el 95% de las de validación. La regla declarada (la 19.ª de las 20 de validación) corrige eso, pero el borde queda ruidoso.

## 2. R14-C v2: observables causales multivariables

**Pregunta:** ¿qué observables intrínsecos separan cada familia de impostores de las geometrías, con N = 4000 y 8000?

**Diseño:**
- 28 estructuras en 8 familias, 10 realizaciones por estructura y N.
- 49 observables por N:
  - 11 por sub-intervalo, con m ∈ {20, 40, 80, 160}: distribución de r, abundancias, 2 a 5-cadenas, cadena más larga, punto medio y acción BD;
  - pendiente de curvatura;
  - grado de links;
  - espectro del Laplaciano del diagrama de Hasse.
- Δ de Glass frente a 60 referencias de la geometría objetivo.
- Holm sobre las 784 filas de la tabla.

### 2.1 Datos por familia (Holm < 0.05)

| Familia | Observables que separan de la plana (de 49), N = 4000 / 8000 |
|---|---|
| Conformes, Kasner y lineales (frente a su propia referencia) | 0 / 0 (calibración correcta) |
| Percolación | 49 / 49 |
| Producto 4 | 41 / 43 |
| Poliédricos | 37 / 41 |
| Conos suaves (9 estructuras) | 1 / 8 |
| Defectos (C5, D4 y D5) | 0 / 0 |

- **Conos suaves con N = 8000:** separan λ_max del Laplaciano (Δ = −2.0), las abundancias con m = 20, L/L_ref con m = 20 y 160, C₃ y C₄ con m = 160, y r con m = 160. Con N = 4000 solo separa λ_max.
- **Defectos:** los efectos medios son enormes (Δ de 13 a 53), pero los intervalos son muy anchos (por ejemplo, de 5.5 a 114), porque la familia junta tres mecanismos distintos con solo 3 estructuras.

### 2.2 Tabla descriptiva por estructura con N = 8000 (calculada después de ver los datos, sin voto)

z = (media de la estructura − media de la referencia)/desviación estándar de la referencia, frente a la plana. La última columna es la geometría que más se le parece solo en λ_max.

| Estructura | z(λ_max) | z(r, m = 160) | z(r, m = 20) | Geometría más cercana en λ_max (\|z\|) |
|---|---|---|---|---|
| C3, cono cilíndrico | −3.51 | 4.91 | 2.98 | de Sitter fuerte (0.40) |
| C4, cono esfera y cubo | −2.79 | 3.19 | 0.90 | D2 (0.67) |
| C6, cono p por eje | −3.59 | 5.41 | 2.78 | de Sitter fuerte (0.29) |
| cono p = 6 | −3.67 | 5.76 | 2.59 | de Sitter fuerte (0.18) |
| cono cúbico | −2.39 | 8.21 | 4.84 | de Sitter (0.44) |
| cono octaédrico | −2.97 | 6.92 | 4.52 | D2 (0.71) |
| D6, D7 y D8 | −0.57 a −1.06 | 0.71 a 1.10 | 0.37 a 0.79 | FRW y D2 (0.09 a 0.22) |
| C5, links borrados | 2.32 | −6.72 | −1.12 | D2 (0.64) |
| D4, relaciones agregadas | −4.68 | 18.2 | 26.5 | de Sitter fuerte (1.15) |
| **D5, gemelos** | **0.49** | **−0.01** | **0.02** | D2 (0.17) |
| de Sitter fuerte (geometría) | −4.03 | 4.27 | 0.24 | de Sitter fuerte (0.28) |
| de Sitter (geometría) | −1.67 | 2.54 | 0.53 | de Sitter (0.30) |

**Lo que muestra esta tabla (nivel 2):**
1. **λ_max solo no sirve como criterio.** La firma de los conos en λ_max (de −2.4 a −3.7) es la misma que la de de Sitter fuerte (−4.0). Para cada impostor hay una geometría a menos de 1.2 desviaciones en λ_max. La curvatura positiva imita a los conos, como en F7.
2. **D6, D7 y D8 son casi invisibles** para todos estos observables (|z| ≤ 1.1), con N = 8000.
3. **D5 (gemelos) es invisible para todos los observables genéricos** (|z| ≤ 0.5). R14-B solo lo detectó porque su vector incluye el detector de gemelos.

## 3. Las tres candidatas, con la evidencia a favor y en contra

| Candidata | A favor | En contra | Costo estimado |
|---|---|---|---|
| **(a)** λ_max del Laplaciano más sub-intervalos grandes, con referencias curvas | Es el único observable que separa los conos suaves en los dos N (a nivel de familia) | Por estructura, de Sitter y de Sitter fuerte reproducen la misma firma (§2.2). Solo, no clasifica. Haría falta un criterio multivariable contra toda la familia curva, y R14-A ya mostró que esa ruta se sobreajusta con facilidad. | 6 h o más (espectro con N = 8000) |
| **(b)** Convergencia de los gemelos bajo dilución, con control de ruido | Responde a la prioridad firmada (H_E frente a H_L). Hay un mecanismo físico explícito (supervivencia q²). Su caída de 7.3× supera la de ruido (2.8×). Ataca directamente el hallazgo de método de R14-B. | D5 es un defecto construido por el proyecto: es una prueba del **mapa** de grano grueso, no de Ω. Hace falta diseñar el control. | 1 a 2 h |
| **(c)** Defectos reorganizados por mecanismo | Corrige una debilidad de diseño evidente de R14-C | Es una regla de diseño más que una pregunta: se puede adoptar como norma sin una prueba propia | sin costo propio |

## 4. Recomendación del orquestador (para deliberar; nada se adopta sin firma)

1. **Preparar (b) primero,** con un control de ruido concreto: el **defecto reinyectado**.
   - Para cada realización de D5, además de diluir el orden con gemelos, se diluye solo la base y después se agregan gemelos en la misma proporción. Así se tiene una estructura con el defecto a densidad **fija** en cada q.
   - **Predicción congelada:** si la caída de D5 es convergencia real, D/radio de D5 diluido cae claramente más rápido que el del control reinyectado. Si es solo pérdida de resolución, caen igual.
   - **El mismo control para los conos:** un cono es autosemejante, así que su control es la misma estructura generada directamente con el N diluido. Eso permite medir si D7 y D8 «convergen» más allá de la resolución.
   - **Unidades:** la realización como unidad, familias nuevas para confirmar (por ejemplo, gemelos al 5% y al 20% y tripletes), un guardián que no haya visto estos datos y semillas nuevas.
2. **No preparar (a) en su forma actual.** λ_max solo queda descartado como criterio por §2.2. Si el Consejo quiere seguir esa ruta, que la redefina como pregunta: «¿algún observable separa los conos de **toda** la familia de Sitter?». Debería preregistrarse con las geometrías curvas como nulo.
3. **Adoptar (c) como regla de diseño** para cualquier prueba futura: las familias se arman por mecanismo (quitar, agregar, duplicar, deformar el cono, curvar), y el análisis por estructura se informa siempre junto al de la familia.
4. **Mantener todo lo de R14-B y R14-C como exploratorio.** Nada de esto confirma H_E ni H_L.

## 5. Preguntas para el Consejo

| # | Pregunta |
|---|---|
| R19-1 | ¿Se prepara (b) con el control de defecto reinyectado como siguiente prueba confirmatoria? |
| R19-2 | ¿Se retira λ_max solo como candidato, y (a) se reformula, si sigue, con las geometrías curvas como nulo? |
| R19-3 | ¿Se adopta (c), familias por mecanismo e informe por estructura, como regla de diseño? |
| R19-4 | ¿Qué familias nuevas y qué guardián para la confirmatoria? |

## 6. Límites declarados

- Todas las estructuras son vistas; las dos pruebas son exploratorias y sin voto.
- La tabla de §2.2 y la comparación de caídas de §1.3 se calcularon **después** de ver los datos: sirven para generar hipótesis, no para confirmarlas.
- Con N ≤ 8000 no se resuelven los conos «redondos casi en todas partes» (D6 a D8), con ningún observable probado.
- El piso de error de Wald en R14-C se agregó después de una prueba corta no oficial. Solo vuelve más conservadoras las pruebas, y está declarado en la aclaración A1.

## 7. Firmas del usuario (2026-10-20)

- **R19-1, R19-2 y R19-3 aprobadas.**
  - Se prepara (b) con el control de defecto reinyectado: `docs/OMEGA_R19_GEMELOS_PREREGISTRO.md`.
  - λ_max solo queda retirado como candidato; (a), si sigue, se reformula con las geometrías curvas como nulo.
  - Las familias se arman por mecanismo y se informa por estructura (regla de diseño).
- **R19-4** (familias nuevas y guardián): pendiente. Hay una propuesta en el preregistro de (b), §6.
