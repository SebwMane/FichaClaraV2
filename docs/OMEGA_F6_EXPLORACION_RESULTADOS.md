# F6 exploratoria: resultados y comparación de ramas (2026-10-14)

- **Estado:** EXPLORATORIA (R11-3). **Sin voto. Nada de esto es evidencia**; sirve para diseñar la confirmación. Protocolo: `OMEGA_F6_EXPLORACION.md` (congelado en e6de7be). Herramienta: `tools/f6_exploracion.py` (auditoría §24). Corrida: commit b0ba414, árbol limpio, 39 partes, unos 70 minutos. Datos: `results/omega_f6_exploracion/`. Figura: `docs/figuras/f6_exploracion.png`.
- Niveles: 1 medido (aquí, descriptivo y exploratorio), 2 explicación propuesta.

## 1. Rama 1: intercepto al achicar los diamantes (Δr = b + A·√(m/N))

| Estructura | Papel | b (IC 95%) | A (IC 95%) | Lectura |
|---|---|---|---|---|
| Diamante plano | geometría | −0.0017 [−0.0027, −0.0007] | +0.008 [−0.000, +0.016] | b casi 0, con un sesgo pequeño (§5) |
| de Sitter | geometría, R > 0 | −0.0007 [−0.0023, +0.0008] | **+0.017** [+0.006, +0.030] | b ≈ 0; A positiva |
| de Sitter fuerte | geometría, R > 0 | −0.0005 [−0.0016, +0.0007] | **+0.033** [+0.023, +0.044] | b ≈ 0; A casi el doble |
| anti de Sitter | geometría, R < 0 | +0.0004 [−0.0006, +0.0013] | **−0.049** [−0.056, −0.042] | b ≈ 0; **A negativa** |
| FRW radiación | geometría, R = 0 | +0.0003 [−0.0012, +0.0018] | −0.012 [−0.024, +0.000] | b ≈ 0; A compatible con 0 al límite |
| FRW polvo | geometría, R > 0 suave | −0.0012 [−0.0026, −0.0000] | −0.002 [−0.012, +0.009] | b ≈ 0; A ≈ 0 (efecto por debajo del ruido, como se predijo) |
| Cono p = 1.5 | zona gris | +0.0013 [−0.0001, +0.0028] | | casi como la geometría |
| Cono p = 3 | zona gris | +0.0016 [+0.0002, +0.0031] | | casi como la geometría |
| Cono p = 6 | impostor | **+0.0062** [+0.0051, +0.0074] | | separado de 0 |
| Cono octaédrico | impostor | **+0.0184** [+0.0172, +0.0194] | | muy separado |
| Cono cúbico | impostor | **+0.0180** [+0.0161, +0.0197] | | muy separado |
| Orden producto 4 | impostor | **+0.0248** [+0.0236, +0.0260] | ≈ 0 | separado, sin pendiente |
| Percolación | impostor | **+0.358** | | lejísimos |

**Lectura (nivel 1, exploratorio):**
- Las seis geometrías tienen intercepto entre −0.0017 y +0.0004.
- Los impostores con aristas quedan entre +0.006 y +0.36.
- Los conos de la zona gris quedan en +0.0013 y +0.0016, pegados a la geometría, lo que es coherente con E-R11-1.
- **El signo de la pendiente sigue al de la curvatura:** positivo en de Sitter (y se duplica al duplicar la curvatura), negativo en anti de Sitter, y compatible con cero en FRW de radiación (R = 0) y en FRW de polvo (curvatura suave).
- Las cinco predicciones del asesor de curvatura de Rev. 11 §2 se cumplieron en la exploración. Al ser exploratorio, esto no es confirmación.

## 2. Descubrimiento: la densidad separa curvatura de forma (nivel 1)

Con m = 160 fijo, al pasar de N = 4000 a N = 8000 (doble densidad, diamantes físicamente más chicos):

| | N = 4000 | N = 8000 | Diluida a 4000 desde 8000 |
|---|---|---|---|
| de Sitter fuerte | +0.0066 | **+0.0042** (se achica) | +0.0063 |
| anti de Sitter | −0.0088 | **−0.0067** (se achica) | −0.0082 |
| Cono octaédrico | +0.0051 | **+0.0076** (no se achica) | +0.0052 |
| Cono cúbico | +0.0070 | **+0.0086** (no se achica) | +0.0055 |
| Orden producto 4 | +0.0253 | +0.0244 (igual) | +0.0241 |

- **La curvatura se achica al subir la densidad; la forma del cono no** (incluso sube un poco). Es la separación predicha (Rev. 11 §2), vista por primera vez.
- **La dilución colapsa:** diluir 8000 a 4000 reproduce el valor de N = 4000 en todas las estructuras (diferencias de 1σ o menos en casi todas las celdas; tres de unas 35 llegan a unos 2σ: anti de Sitter con m = 10, cono cúbico con m = 160 y FRW de radiación con m = 160, lo esperable por azar).
- **Interpretación propuesta (nivel 2):** en un cono L^p los intervalos no son todos equivalentes (no hay boosts), así que la mezcla de orientaciones de los intervalos con m elementos cambia con su tamaño relativo. Por eso la desviación del cono depende de m y de N, y crece hacia su intercepto cuando el diamante se achica. En Minkowski todos los intervalos son equivalentes por boosts y ese efecto no existe.

## 3. Rama 2: perfil local (max-T frente a la nula, N = 8000)

| | max-T (mediana) | p95 de la nula |
|---|---|---|
| Geometrías (plano, de Sitter, de Sitter fuerte, anti de Sitter, FRW) | 0.88 a 1.01 | 1.19 |
| Conos p = 1.5 y p = 3 (zona gris) | 0.86 a 0.90 | |
| Cono p = 6 | 1.00 | |
| Cono octaédrico | 1.26 | |
| Cono cúbico | 1.37 | |
| Orden producto 4 | 2.29 | |
| Percolación | 34 | |

**Lectura (nivel 1):**
- **Lo bueno:** el perfil local **no confunde curvatura con forma**: todas las geometrías quedan bajo el ruido, incluido de Sitter fuerte, que hizo fallar a F7. Medir en diamantes chicos resolvió esa falla, como se esperaba.
- **Lo malo:** perdió resolución. Los conos con aristas quedan apenas sobre el ruido (1.26 y 1.37 frente a 1.19) y el cono p = 6 no se distingue.

## 4. Comparación de las dos ramas

| | Rama 1 (intercepto de r) | Rama 2 (perfil local) |
|---|---|---|
| Separa geometría de conos con aristas | **sí, con margen amplio** (≈ 10 veces el sesgo del plano) | apenas |
| Detecta el cono p = 6 (ε = 0.31) | **sí** (b = +0.006) | no |
| Robusta a la curvatura | sí: la curvatura va a la pendiente, no al intercepto | sí |
| Mide la curvatura (signo y tamaño) | **sí** (A ∝ curvatura, con signo) | no |
| Órdenes producto y percolación | sí | sí |
| Perillas | una (el modelo lineal) | varias (k ≤ m/2, forma, max-T) |

**Conclusión exploratoria (nivel 1 descriptivo, sin voto):** la rama 1 domina. Hace lo que hacía la rama 2 y además separa conos con más margen y mide la curvatura. La rama 2 queda como control: un orden que pase la rama 1 pero tenga un perfil local anómalo merece atención.

## 5. Sesgos y límites que aparecieron (los que el usuario pidió vigilar)

- **S1 y S2 confirmados, la pregunta original de F6 no se puede responder con estos N.** El bloque solo es medible con m = 10 (N = 8000) y m = 10, 20 (N = 16000): su ajuste tiene 3 puntos y un intercepto de +0.001 ± 0.016, inútil. La caja solo es medible con m = 10 en N = 8000.
- **S3 visible: el borde imita a un impostor.** En la caja (49% de intervalos recortados) Δr = +0.007 con m = 10, del tamaño del intercepto del cono p = 6. Una región con borde podría rechazarse por su borde, no por su forma.
- **Sesgo pequeño del plano:** el diamante plano tiene b = −0.0017 (el IC excluye 0), y la nula tiene un Δr medio de unos −0.0005. Es pequeño frente a los impostores, pero debe entrar en la calibración del umbral y no suponerse cero.
- **S4 (solapamiento):** los errores se calcularon entre realizaciones, como se pidió. Con 10 realizaciones algunos IC son estrechos; la confirmación debe usar más.
- **FRW de radiación:** A = −0.012 con un IC que toca 0. Si se repite con signo negativo, la teoría (R = 0) no lo explica a primer orden (habría que mirar R_ab T^a T^b, que no es cero en FRW).

## 6. Ramas y preguntas nuevas que dejan los datos

1. **«Curvímetro» intrínseco:** A da el signo y el tamaño de la curvatura sin coordenadas. Pregunta: ¿A es proporcional a R conocida en una familia (de Sitter con varios radios)?
2. **La dependencia en N de los conos:** ¿hacia dónde converge su desviación al achicar? ¿El intercepto mide la «anisotropía» del cono (relación con ε)?
3. **Bordes y topología (pregunta original de F6):** requiere N ≳ 30 000 o regiones de otra proporción. Queda abierta y separada.
4. **El sesgo del plano:** ¿viene de la selección de intervalos con m exacto o de la referencia? Es barato de investigar.

## 7. Propuesta del Consejo para la confirmación (requiere la firma del usuario)

- **Criterio único que vota:** el intercepto b de la rama 1, con el modelo Δr = b + A·√(m/N), m ∈ {10, 20, 40, 80, 160} (regla m ≤ 0.1·m_max) y N ∈ {4000, 8000}. Clasificación «espacio-tiempo 4D» si |b − b_plano| ≤ τ_b. La rama 2 y la pendiente A se informan sin voto (A con su signo como predicción ciega: de Sitter, anti de Sitter).
- **Umbral:** τ_b por regla, con b calculado en grupos de calibración (diamantes planos y curvos con entropía nueva, al menos 150 realizaciones) y fijado antes de abrir el conjunto sellado.
- **Alcance:** regiones con forma de intervalo. Los bordes y la topología quedan para una prueba aparte con N mayor (§6.3).
- **Estructuras:** el conjunto sellado `CONFIRMACION_CIEGA.md` (verificando su hash) más las geometrías y los impostores ya vistos con semillas nuevas.
- **Presupuesto:** intento 1 de 2 (Rev. 11 §6).
