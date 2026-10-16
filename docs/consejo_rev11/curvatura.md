# Asesor de curvatura local: F6 y la exploración de tamaños

Niveles: 1 medido (aquí, piloto declarado), 2 propuesta, 3 externo, 4 extrapolación.

## 1. Qué predice la teoría

**Externo (3).** En un intervalo causal pequeño de tiempo propio τ, el volumen y los conteos de cadenas e intervalos son los de Minkowski más correcciones de orden τ² que dependen solo de R y de R_ab T^a T^b, con T el eje del diamante (Myrheim 1978, Meyer 1988, Gibbons y Solodukhin 2007, Khetrapal y Surya 2013, Roy, Sinha y Surya 2013). Benincasa y Dowker (2010) extraen R de conteos locales. Glaser y Surya (2013) proponen justamente usar abundancias de intervalos en sub-intervalos pequeños como prueba de «parecerse a variedad»: F6 tiene antecedente directo. Los coeficientes exactos deben tomarse de esos trabajos, no de mi memoria.

**Propuesta (2).** A densidad ρ, un sub-intervalo con m elementos tiene τ^d ∝ m/ρ, así que la desviación por curvatura es ∝ R·(m/ρ)^{2/d}: en 4D, ∝ (m/N)^{1/2} con la región fija. Es lineal en la curvatura, así que cambia de signo entre de Sitter (R > 0) y anti de Sitter (R < 0). El tensor de Weyl no entra a ese orden. En de Sitter R_ab ∝ g_ab y la corrección no depende del boost del sub-intervalo; en FRW sí. Un cono L^p es plano, invariante por traslación y autosemejante: al escalar o al diluir, sus sub-intervalos de m elementos siguen igual. Su desviación depende de m, pero **no de la densidad**.

**Cómo separar (2).** Hay dos palancas:
- **Recorrer m con un observable intensivo.** La fracción de pares relacionados r (Myrheim-Meyer) del interior de un sub-intervalo tiene esperanza independiente de m en cualquier forma fija. Entonces: curvatura da r(m) − r_plano = A·(m/ρ)^{2/d} → 0; un cono da una constante b ≠ 0. La prueba es el **intercepto** b.
- **Diluir a m fijo.** Quedarse con una fracción q de los elementos al azar es intrínseco (no usa coordenadas) y equivale a otra densidad. La curvatura crece como q^{−2/d}; el cono no cambia. Sirve también para observables no intensivos (links por elemento), cuya sensibilidad crece con m.

Supuestos: muestreo uniforme en volumen propio (Poisson o binomial), región causalmente convexa, impostor autosemejante (un orden con escala propia o un cono no redondo y además curvo no se separa así), y m bastante chico para que domine el término τ².

## 2. Hipótesis completa

**En simple.** Una estructura es espacio-tiempo de dimensión d si cada diamante chico que se recorta de ella se parece a un diamante de Minkowski de dimensión d, y lo que no se parece se achica al achicar el diamante, a un ritmo fijo: en 4D, a la mitad cuando el número de elementos del diamante baja a la cuarta parte. Lo que se achica así es curvatura. Lo que no se achica es forma del cono, y eso la descarta.

**En técnico (2).** H_loc(d): para observables intrínsecos O de sub-intervalos con m elementos interiores, tras diluir con fracción q,
E[O | m, q] = O_plano,d(m) + (m/q)^{2/d}·[α_O R̄ + β_O ⟨R_ab T^a T^b⟩] + O((m/q)^{4/d}),
con O_plano,d(m) exacto (m puntos uniformes en un diamante plano más sus puntas) y (i) intercepto nulo dentro de una tolerancia δ declarada; (ii) exponente 2/d con la misma d que da Myrheim-Meyer; (iii) los datos dependen solo de m/q (colapso) en los observables intensivos; (iv) todos los observables se explican con las mismas dos cantidades de curvatura.

**La refuta:** intercepto distinto de cero más allá de δ; desviación que no baja al diluir; exponente incompatible con 2/d; d del exponente distinta de d de Myrheim-Meyer; observables que piden más de dos parámetros de curvatura. Como a Sebastían no le importa p = 3 frente a redondo, δ debe fijarse explícitamente para aceptar «casi redondo» (2); eso hay que escribirlo, no esconderlo.

## 3. Predicciones para congelar antes de ver datos

- **Exponente (2):** 2/d; 0.5 en 4D. Piloto (1): r en de Sitter fuerte da 0.60 [0.54, 0.66] y en anti de Sitter fuerte 0.42 [0.34, 0.49]; la media, 0.51. Lo explico (2) con el término siguiente (∝ R², igual signo en ambos), que empina uno y aplana otro. Recomiendo congelar el modelo b + A·x^{1/2} + B·x y probar solo b, no ajustar el exponente libre.
- **Signos (1 en piloto, 2 en general), para los observables locales:** R > 0 (de Sitter): r sube y los links por elemento suben (Myrheim-Meyer local baja). R < 0 (anti de Sitter): ambos bajan. FRW de radiación tiene R = 0: sin desplazamiento isótropo a primer orden (piloto: nada detectable). FRW de polvo (a ∝ η²): R > 0 pero, en unidades de la región, unas 10 a 25 veces menor que en de Sitter fuerte; predicción: subida de r por debajo del ruido con N ≤ 4000 (piloto: −0.004 ± 0.008 en m = 40). Estos signos son locales; el de F7 (perfil global, de Sitter fuerte con menos links) mide otra cosa y no se contradice.
- **Dilución (2):** a m fijo, duplicar la densidad multiplica la desviación por 2^{−1/2} ≈ 0.71 en 4D. Piloto (1), links por elemento en de Sitter fuerte: exponente en N 0.60, compatible.
- **Magnitud (1):** con N = 4000 y m = 20, de Sitter fuerte se desvía +1.6% en r y +0.2% en links por elemento. Es unas 14 veces menos curvatura que en la región entera ((m/N)^{1/2} ≈ 0.07). Esto es lo que podría arreglar la falla 2 de F7 (2).
- **Lo que no sé (4):** cuánto se desvían los conos L^p o poliédricos a m = 20. No lo medí (no se permite piloto con impostores). Si su desviación local es del orden del ruido, F6 rechazará curvatura con holgura y no verá conos: se cambia una falla por otra.

## 4. Riesgos

- **Tamaños chicos (1):** con m = 10 la curvatura fuerte ya es invisible (links +0.01%) y el error de r por realización es de 1% a 3%. Con m ≤ 10 nada se distingue.
- **Errores mal calculados (1), el más grave:** los sub-intervalos de una misma realización comparten elementos. Con error calculado entre intervalos, el plano mostró «tendencias» de 2σ a 3σ que no existen. Con 8 realizaciones, r del plano dio un intercepto de −1.75% a 2.7σ; con 40 realizaciones, 0.1005 ± 0.0009 frente a 0.1006 de referencia: era ruido. El error debe calcularse entre realizaciones independientes, con unas 30 o más.
- **Tamaños grandes (2):** con m/N por encima de un 2% los sub-intervalos se concentran en el centro de la región, crece el término (m/N)^{4/d} y en regiones no convexas (la caja con bordes) quedan recortados. La caja debe compararse con el bloque a los mismos m/N.
- **Densidad no uniforme en coordenadas (2):** los observables son intrínsecos, así que no es problema si el muestreo es uniforme en volumen propio. Pero R·ρ^{−2/d} varía dentro de FRW y de Sitter: se mide un promedio, y la dilución debe ser al azar sobre elementos, nunca por coordenadas.
- **Boost (2):** a m fijo hay sub-intervalos muy alargados; en FRW su corrección crece con γ². En el piloto la correlación entre alargamiento y links fue nula en todas (1), pero con más curvatura habría que informarla.
- **Sobreajuste (2):** la rejilla de m y q, el observable (r interior), el modelo con exponente fijo y δ deben congelarse antes; explorar tamaños después de ver datos es justo lo que convierte F6 en ajuste.

## 5. Qué corrí (piloto declarado, no es evidencia)

Scripts en `scratchpad/asesores/`: `piloto_curv.py`, `piloto2.py`, `control_plano.py`, `control2.py`, `analiza*.py`. Semillas SeedSequence(900101, spawn_key=...). Solo positivas: diamante plano 4D, de Sitter (η ∈ [−1.5, −0.5]), de Sitter fuerte (η ∈ [−1.2, −0.2]), anti de Sitter de Poincaré (z ∈ [0.2, 1.2], densidad ∝ z^{−4}, añadido por mí como positiva con R < 0), FRW de radiación y de polvo. N ∈ {1000, 2000, 4000, 8000}, m ∈ {5, 10, 20, 40, 80}, hasta 1500 sub-intervalos por corrida; piloto 2 con 8 realizaciones por caso; referencia con 20000 a 40000 diamantes planos. Ningún impostor. Nada en el repositorio.

**Mi evaluación:** F6 reformulada como «intercepto bajo dilución y recorrido de m» tiene base teórica y el piloto confirma la parte de curvatura. Su punto débil es la resolución ante conos a m chico, que solo una prueba preregistrada puede medir. La decisión es de Sebastían.

## Ronda 2

**Acuerdos.** Calcular los errores entre corridas, no entre intervalos. Separar una fase exploratoria de una confirmatoria de un solo disparo, con todo congelado y un solo criterio que vote (sobreajuste, ramas). Escalar N con m. Agregar impostores nuevos.

**Discrepancias.** T = max_m razón_m prueba cada m por separado y desperdicia lo que distingue curvatura de forma, que es cómo cambia la desviación con el tamaño (2). Prefiero que vote el intercepto bajo dilución, con max-T solo informado. Hay que elegir uno antes de correr.

**Objeción de tamaño finito.** Tiene razón. En el bloque plano, que los intervalos no sean locales no sesga nada, porque la referencia es exacta sin curvatura (2). En FRW bloque y de Sitter fuerte, una densidad que varía ×9 a ×12 dentro del intervalo significa que R·τ² es del orden de 1. Estamos fuera del régimen perturbativo y mi modelo b + A·x^{1/2} + B·x no vale ahí (2). Mi +1.6% es un promedio sobre intervalos de ese tipo (1). Remedio: curvatura más suave o N ≳ 800 m, y declarar no medible cualquier caso con variación de densidad por intervalo mayor que ×2 (2).

**Objeción de sobreajuste.** Es correcta y es la más importante. En geometrías conformemente planas el orden es exactamente el de Minkowski y solo cambia la densidad, así que el piloto no puso a prueba ninguna deformación de los conos (1). H_loc no cambia de forma, pero hace una predicción nueva y fuerte: en vacío no conformemente plano (Kasner, campo débil fuera de la masa) se anulan R y R_ab, el término τ² desaparece y la desviación debe ir como (m/ρ)^{4/d}, con exponente 1 en 4D (2, apoyado en 3). Si el exponente sale 1/2, o no hay desaparición, H_loc queda refutada. Hay que congelar eso antes de medir.

**Anti de Sitter.** Fue un desvío del encargo. La lista era cerrada y no consulté. Su exponente medido (0.42) inspiró el término B·x: el modelo salió en parte de datos ya vistos. Debe tratarse como exploratorio y sin voto. Si el Consejo lo quiere en la prueba oficial, que entre como positiva declarada desde ya, con semillas nuevas y sin usarse para calibrar umbrales.
