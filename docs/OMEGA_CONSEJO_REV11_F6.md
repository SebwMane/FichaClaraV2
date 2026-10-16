# Consejo Ω, Revisión 11: F6 reconsiderada con asesores nuevos (2026-10-14)

- **Estado:** deliberación cerrada. **Nada preregistrado ni ejecutado como prueba oficial.** Las decisiones las toma Sebastían (§7).
- **Pedido del usuario:** «que el Consejo reconsidere F6, llamar a nuevos asesores para discernir la propuesta; la decisión la tomaré yo». Le interesa supervisar los sesgos por diferencias de tamaño entre diamantes, la exploración de tamaños, la hipótesis completa articulada y las ramas nuevas. Orden: **supervivencia de la metodología sin sobreajuste**.
- **Asesores nuevos (subagentes independientes, dos rondas):** tamaño finito y selección; curvatura local; inferencia y sobreajuste; explorador de ramas. Los textos completos, con su ronda 2, están en `docs/consejo_rev11/`. Los pilotos declarados (solo geometrías positivas, entropías 900001 y 900101, fuera de la serie oficial) están en `docs/consejo_rev11/pilotos/`.
- **Figura:** `docs/figuras/rev11_f6_tamanos.png` (A: qué propone F6; B: F7 medido; C: hipótesis H_loc, no son datos).

## 1. Decisión del usuario registrada: enmienda E-R11-1

> **E-R11-1 (2026-10-14), del usuario.** Para Ω **no** se exige distinguir un cono casi redondo de uno redondo. Motivo del usuario: hay que explorar y analizar los datos en conjunto, dejando que surjan ramas nuevas.
> - **Se tomó después de ver F7**, y coincide con la resolución observada del instrumento: **no cuenta como éxito independiente.**
> - **¿Puede crear un positivo? Sí** (quita negativas que el criterio fallaba). Vale **solo hacia adelante**, con semillas nuevas.
> - **No reabre F7**, que sigue FALLA. Además, F7 seguiría fallando aunque se aplicara: de Sitter fuerte (1.55) no depende de esta decisión.
> - **Falta fijar el tamaño de la zona gris (R11-1).** Una medida propuesta es la deformación ε = max |‖u‖_p/‖u‖₂ − 1| sobre las direcciones: p = 3 tiene ε ≈ 0.17, p = 1.5 ≈ 0.20, p = 6 ≈ 0.31, el cubo ≈ 0.42 y el octaedro ≈ 0.73.

## 2. La hipótesis completa, articulada

**En simple.** Debajo del espacio-tiempo habría una red de relaciones «antes / después» que puede organizarse en fases, como el agua en hielo o vapor. Ω tiene tres partes:

1. **Reconocer (H_rec):** dada una red, saber si está en la «fase espacio-tiempo». Es lo que prueban F1 a F7, y aquí F6.
2. **Seleccionar (H_sel):** por qué la red estaría en esa fase. Sigue abierta: Ω-C no seleccionó geometría y, contando órdenes al azar, casi todos son de tipo KR.
3. **Desviarse (H_fen):** qué se vería distinto de la relatividad. Lejana.

**La versión local de H_rec (H_loc, del asesor de curvatura, nivel 2).** Una red es espacio-tiempo de dimensión d si cada diamante chico recortado de ella se parece a un diamante de Minkowski de dimensión d. Lo que no se parece debe **achicarse al achicar el diamante**, a un ritmo fijo: en 4D, a la mitad cuando el número de puntos baja a la cuarta parte. Lo que se achica así es **curvatura**; lo que no se achica (el **intercepto**) es **forma del cono**.

En técnico: para un observable intensivo O de sub-intervalos con m elementos (por ejemplo la fracción de pares relacionados r del interior),

  E[O | m] = O_plano,d(m) + A·(m/ρ)^{2/d} + …

con O_plano,d exacto (m puntos uniformes en un diamante plano) y A ∝ R (nivel 3, Myrheim, Meyer, Benincasa-Dowker, Roy-Sinha-Surya; Glaser y Surya ya propusieron usar abundancias en sub-intervalos).

**La refuta:**
- un intercepto distinto de cero más allá de la zona gris;
- una desviación que no baja al achicar;
- un exponente incompatible con 2/d, o con una d distinta de la de Myrheim-Meyer;
- en vacío no conformemente plano (Kasner, campo débil), un exponente distinto de 4/d (= 1 en 4D), porque allí el término τ² desaparece.

**Predicciones congelables antes de medir (nivel 2):**
- el signo de A es opuesto en de Sitter (R > 0: r sube) y anti de Sitter (R < 0: r baja);
- FRW de radiación (R = 0) no tiene término τ²;
- a m fijo, duplicar la densidad multiplica la desviación por curvatura por 0.71 en 4D.

## 3. Sesgos por tamaño de los diamantes (lo que el usuario pidió supervisar)

| # | Sesgo | Evidencia | Qué lo controla |
|---|---|---|---|
| S1 | **En el bloque, m = 20 no es local**: son casi los intervalos más altos que caben (centros en el 10% central del tiempo) | piloto, nivel 1 | m ≤ 0.1·m_max medido desde el orden; o N ≳ 800·m |
| S2 | **El número de diamantes depende de la forma de la región**, no de que sea espacio-tiempo. Con N = 4000 y m = 40: 3066 en el diamante, 34 en el bloque, 4 en FRW bloque | 1 | calibrar con el mismo número de intervalos usados (n_sel); la regla «menos de 50 = no medible» castiga primero a positivas curvas |
| S3 | **Recorte por el borde de la caja**: 49% a 67% de los intervalos recortados; +2.1% links y +3.8% relacionados, del mismo signo que de Sitter, y crece con m | 1 | estratificar interior/borde; comparar caja y bloque a igual m/m_max |
| S4 | **Solapamiento**: cada elemento aparece 3 a 17 veces; con errores «ingenuos» el plano mostró tendencias falsas de 2σ a 3σ | 1 (dos pilotos, independientes) | errores **solo entre realizaciones**, 30 o más |
| S5 | **Densidad variable dentro de un intervalo** en FRW y de Sitter fuerte (×9 a ×12): la curvatura no es pequeña, la corrección τ² no vale | 1 | curvatura más suave, o declarar no medible si la variación es mayor que ×2 |
| S6 | **Bins degenerados con m chico** (k ≥ m siempre 0 o 1): diluyen D_forma | 1 | usar k ≤ m/2, o r del interior sin puntas |
| S7 | **La distribución nula cambia con m**, no monótona | 1 | umbral por (m, N, n_sel); nunca un umbral único |
| S8 | **Dilución q**: diluir agranda los intervalos y empeora S1, S3 y S4 | 2 | comparar a igual m/m_max; remuestrear realizaciones completas |
| S9 | **Elegir el m que mejor funcione** después de mirar | metodológico | rejilla y m votante congelados antes |
| S10 | **Asimetría**: los pilotos solo miraron positivas; un umbral ajustado para aceptarlas sin saber cómo se ven los conos a m chico | metodológico | explorar también impostores ya vistos, con entropías de piloto, antes de congelar |

**Advertencia transversal del asesor de sobreajuste (nivel 1, cuenta):** en la línea F1 a F7 hay 7 bifurcaciones decididas tras ver datos y unas 8 perillas sin barrido. Los seis criterios con voto en 4D fallaron. **El riesgo de sobreajuste para la próxima prueba es alto** si no se separan exploración y confirmación.

## 4. Mapa de la interacción entre asesores

**Consenso (los cuatro):**
- errores entre realizaciones;
- un solo criterio con voto, todo congelado antes;
- F7 sigue FALLA;
- E-R11-1 vale solo hacia adelante;
- escalar N con m o restringir m;
- hacen falta impostores y positivas nuevos.

**Tensiones abiertas (las decide el usuario):**

| Tensión | Posición A | Posición B |
|---|---|---|
| **Qué vota** | **Intercepto** bajo recorrido de m (curvatura): usa lo que distingue curvatura de forma, que es cómo cambia con el tamaño | **Perfil local** con k ≤ m/2 y max-T sobre 3 valores de m (sobreajuste): menos perillas; el modelo de tres parámetros con pocos x deja el intercepto libre de absorber cualquier cosa |
| **Dilución q** | Eje de prueba (curvatura) | Solo control exploratorio (sobreajuste; tamaño finito: empeora la localidad) |
| **N** | Escalar a N ≳ 16 000 para que m = 20 sea local en el bloque (tamaño finito; costo no medido) | Restringir el bloque a m ≤ 10 y declararlo (con m ≤ 10, curvatura: «no se distingue nada») |
| **Anti de Sitter** | El asesor de curvatura lo añadió fuera de la lista y lo reconoce como desvío; su dato inspiró el término B·x del modelo | Queda «visto»: exploratorio, sin voto, o positiva declarada con semillas nuevas y sin calibrar con él |

**Correcciones mutuas (ejemplos de la ronda 2):**
- El explorador corrigió su rama R-A: la desviación de un cono sí depende de m, pero no de la densidad.
- La asesora de tamaño finito retiró su sesgo de lugar en de Sitter: es un efecto de la carta, que no sesga con muestreo uniforme en volumen propio. En FRW sí sesga.
- El asesor de curvatura aceptó que su modelo no vale en curvatura fuerte (S5).

**Lo que nadie sabe (nivel 4, el punto más importante):** cuánto se desvían los conos con aristas en diamantes chicos. Si su intercepto local es del orden del ruido, F6 cambia la falla de curvatura por una falla de resolución.

## 5. Ramas nuevas que surgieron

| Rama | Pregunta | Qué la haría fallar | Costo | Comparte instrumento con F6 |
|---|---|---|---|---|
| R-A+B (fundidas) | Intercepto (forma) y pendiente (curvatura) del perfil local frente al tamaño | intercepto ≠ 0 en geometría; pendiente no monótona en la curvatura | bajo sobre F6 | sí |
| R-C | Mapas locales de dimensión y curvatura (dominios «hielo en agua») | el mapa de FRW no sigue la ley conocida | bajo, pero ruido sin N grande | sí |
| R-D | Heterogeneidad entre sub-intervalos como discriminador | no atrapa conos ni órdenes producto (homogéneos) | casi nulo | sí |
| R-E | Clase ampliada «conos suaves» (estructuras de tipo Finsler); Lorentz exacto como pregunta aparte | | una decisión (E-R11-1) | no |
| **Signo dS/AdS** | Predicción ciega del signo de la pendiente | signo equivocado | bajo | sí |
| **Positiva no conformemente plana** (Kasner, campo débil) | ¿Exponente 4/d en vacío? Hasta hoy toda la curvatura probada solo cambiaba densidad | exponente 1/2 o sin desvanecimiento | **desarrollo** (relación causal por geodésicas) | sí |
| **Minkowski con links borrados** | Impostor que imita el hueco de links de la curvatura y se diluye; ¿defecto microscópico o microestado de la misma geometría? (OMEGA_MASTER §4) | | bajo | sí |
| R-F | Usar el reconocedor como juez de dinámicas (crecimiento secuencial, acciones, pesos entrópicos): conecta con H_sel | ninguna dinámica de un catálogo congelado produce espacio-tiempo | alto (semanas) | usa R como juez |

## 6. Plan que sale de la deliberación (propuesta, no decisión)

1. **Fase exploratoria declarada (sin voto).** Entropías de piloto; puede mirar positivas e **impostores ya vistos**. Rejilla de m declarada antes de empezar, que no se amplía. Mide el costo real de N = 8000 a 16 000. Responde lo que nadie sabe (§4).
2. **Congelar** un único criterio, la rejilla, la regla de N, n_sel, ε₀, el umbral con ≥ 150 calibraciones y el control de comparaciones múltiples.
3. **Confirmación de un solo disparo** con entropías nuevas e impostores y positivas nuevos a ciegas, definidos por un subagente que no participa en el diseño. **Presupuesto: dos confirmaciones**; luego vuelve al Consejo.
4. **Conjunto de reserva:** se publica solo su hash SHA-256 y se abre una vez, tras un PASA, con N mayor.

## 7. Decisiones del usuario (pendientes)

| # | Pregunta |
|---|---|
| R11-1 | Zona gris de E-R11-1: ¿incluye solo conos tan redondos como p = 3 (ε ≤ 0.17), también p = 1.5 (ε ≤ 0.20), o se fija con otro criterio? |
| R11-2 | ¿Qué vota en F6: el intercepto bajo recorrido de tamaños, o el perfil local con max-T? |
| R11-3 | ¿Autorizar la fase exploratoria (§6.1), incluyendo impostores ya vistos, antes de congelar nada? |
| R11-4 | ¿Incluir en la confirmación una positiva no conformemente plana (requiere desarrollo) y el conjunto de reserva con hash? |

## 8. Firmas y decisión del Consejo (2026-10-14)

| # | Decisión |
|---|---|
| R11-1 | **Delegada al Consejo** por el usuario. Decisión en §8.1. |
| R11-2 | **Firmada:** primero se explora la rama del **intercepto al achicar los diamantes**; según su resultado, se prueba el **perfil local**; se comparan ambos resultados y lo que descubran. |
| R11-3 | **Firmada:** fase exploratoria sin voto, que puede mirar impostores ya vistos, antes de congelar nada. Protocolo: `OMEGA_F6_EXPLORACION.md`. |
| R11-4 | **Firmada:** la confirmación incluirá una positiva no conformemente plana (Kasner de vacío, requiere desarrollo) y habrá conjunto de reserva sellado con hash. |

### 8.1 Decisión del Consejo sobre R11-1: zona gris ε₀ = 0.20 (incluye p = 1.5 y p = 3)

- **Argumento (físico y estadístico).** El usuario fijó un criterio geométrico: «casi redondo». La medida ε no distingue p menores ni mayores que 2; mide cuánto se aparta el cono del redondo. p = 1.5 (0.20) y p = 3 (0.17) se apartan casi lo mismo. Fijar ε₀ = 0.17 haría depender la frontera de una diferencia de 0.03, menor que la resolución demostrada del instrumento (F7): sería diseñar una prueba que se sabe que fallará por una distinción que el usuario declaró irrelevante.
- **Objeción (escéptico y sobreajuste), aceptada como condición.** Es la segunda reducción de alcance después de ver F7. Se registra así y cuenta en el presupuesto de intentos. Para que la frontera no quede «a medida»:
  1. ε₀ = 0.20 se fija **hoy, antes de cualquier dato de la exploración**, y no se mueve.
  2. Los conos con ε ≤ 0.20 se informan como **«zona gris»**: no votan ni como positivas ni como negativas.
  3. La confirmación debe incluir negativas a ciegas **cerca de la frontera por fuera** (ε entre 0.25 y 0.35, por ejemplo un cono de sección cilíndrica, ε ≈ 0.29), y ninguna positiva puede caer en la zona gris.
- **¿Puede crear un positivo?** Sí. Vale solo hacia adelante. F7 sigue FALLA.
