# F6 reconsiderada: opinión de la asesora de tamaño finito y selección

## Piloto declarado (no es evidencia)

`scratchpad/asesores/piloto_tamano.py` (salida `.json`), `SeedSequence(900001, …)`, N = 4000. Solo positivas: diamante 4D y bloque T³ × [0, 0.45] (6 corridas cada uno); caja con bordes, de Sitter, de Sitter fuerte y bloque FRW a ∝ η (4 cada uno). m ∈ {5, 10, 20, 40, 80, 160}, hasta 300 sub-intervalos por corrida, referencia con 3000 diamantes de m puntos más las puntas.

## 1. Sesgos de selección por tamaño

**Cuántos hay depende de la forma de la región, no de que sea espacio-tiempo (1).** Mediana de intervalos con exactamente m interiores:

| m | diamante | bloque | caja | FRW bloque |
|---|---|---|---|---|
| 20 | 6572 | 904 | 572 | 382 |
| 40 | 3066 | 34 | 14 | 4 |
| 80 | 1214 | 0 | 0 | 0 |

En el bloque m_max ≈ 0.012 N para cualquier toro con duración < lado/2 (2). Con N = 4000 y m = 20, los centros caen en el 10% central del tiempo (cuantiles 0.45 a 0.55; en el diamante 0.38 a 0.63) (1). **En el bloque, m = 20 no es local: son casi los intervalos más altos que caben** (1).

**Exactitud «Poisson condicionada».** Dado (x, y) con m interiores, esos puntos son uniformes en su diamante (3). Se rompe en dos casos:
- **Caja:** el 49% (m = 5) a 67% (m = 20) de los intervalos elegidos tienen el diamante geométrico recortado por el borde (1). A m = 20: links +2.1% en recortados frente a −0.6% en enteros; fracción de relacionados +3.8% (1). Pequeño, sistemático, crece con m.
- **Densidad variable:** en coordenadas conformes la densidad dentro de un intervalo cambia ×9 (FRW bloque, m = 20) y hasta ×12 (de Sitter fuerte, m = 20) (1). Incluso m = 5 abarca ×3.6 en de Sitter (1): en 4D los intervalos con pocos elementos suelen ser largos y casi nulos (2). En de Sitter, por simetría, solo importa el volumen propio (2); en FRW no.

**Sesgo de lugar (1).** Con m fijo, las zonas más densas aportan más intervalos: en de Sitter fuerte la mediana temporal está en 0.64, no en 0.5. El promedio pesa más la parte tardía.

**Señal de curvatura frente a m (1).** En de Sitter fuerte, el desplazamiento de links v₀ sube con m: +0.007 (20), +0.020 (40), +0.027 (80), +0.033 (160); la fracción de relacionados, de +1.9% a +5.1%. En el diamante plano todo queda en ±0.01. **La curvatura se reduce al bajar m, pero lentamente**: en 4D el tamaño lineal va como (m/N)^{1/4}, así que m = 20 con N = 4000 ocupa ≈ 27% del alto de la región (2).

**Solapamiento (1).** Fracción de pares de sub-intervalos que comparten algún elemento: 0.12 (m = 20), 0.25 (40), 0.53 (80), 0.79 (160). Cada elemento aparece en promedio 3.3, 5.5, 9.7 y 17.5 veces. Con m ≥ 80 la dispersión entre corridas de v₀ es 1.5 a 3 veces el error ingenuo (1, con pocas corridas). **Los 300 intervalos no son 300 datos independientes.**

## 2. Ruido y potencia según m

- m pequeño: con m < 8 los bins k > m son cero y k = m vale siempre 1 (las puntas). D_forma promedia esos ceros y se achica (1: nula 0.029 con m = 5, 0.050 con m = 20).
- m intermedio: el ruido de D_forma lo pone el bin más escaso (k = 8 con m = 20: dispersión 0.05 entre corridas, frente a 0.004 en links) (1). El máximo sobre k elige el bin más ruidoso.
- m grande: pocos intervalos fuera del diamante, más solapamiento, más curvatura.
- **La distribución nula de D_forma cambia con m de forma no monótona** (1). Un umbral común a varios m es incorrecto.

## 3. Cómo explorar tamaños

1. Grilla fijada antes de mirar: m ∈ {10, 20, 40, 80}, y 160 en el diamante. **Escalar N con m** para que m ≤ 0.1 m_max en el bloque (N ≳ 800 m; m = 20 pide N ≈ 16 000; costo O(N³) sin verificar, 4). Si no es viable, declarar que el bloque solo se mide con m ≤ 10.
2. Calibrar por (m, N) y por número de intervalos usados: submuestrear la calibración plana al mismo n_sel que tenga cada estructura, para que una positiva con 34 intervalos no se rechace solo por ruido.
3. Controlar comparaciones múltiples con un solo estadístico: T = max_m (razón_m), con su umbral calibrado como máximo sobre las mismas corridas planas. Alternativa: un m primario vota y el resto se informa.
4. Usar el promedio de distancias estudentizado por bin (dividir v_k por su error de calibración) o limitar k a k ≤ m/2, para que el bin escaso no domine.
5. Errores siempre entre corridas, nunca el error ingenuo dentro de una corrida.
6. Reportar por estructura y m: intervalos (mediana, mínimo), cuantiles del tiempo central, fracción recortada, rango de densidad interna, multiplicidad, v_k con error entre corridas, fracción de relacionados, d*_loc, D_forma y razón.
7. Rama nueva (2, no probada): los conos L^p son autosemejantes, su defecto de links no debería bajar con m; el de la curvatura sí. **La pendiente de v₀ frente a m** podría separarlos. Requiere preregistro.

## 4. Riesgos a vigilar

1. Llamar «local» a un m que en el bloque es casi el intervalo máximo.
2. La regla «menos de 50 = no medible» depende de la forma: castiga primero a las positivas curvas y al bloque (FRW bloque tiene 4 intervalos con m = 40).
3. Recortes en la caja: sesgo que crece con m y puede leerse como «forma».
4. Densidad variable: la referencia «exacta» deja de serlo en FRW.
5. Ponderación temporal hacia zonas densas.
6. Error subestimado por solapamiento, que hace parecer significativas diferencias de ruido.
7. Umbral único para varios m, o elegir el m que mejor resulte después de ver (sobreajuste directo).
8. Comparar estructuras con distinto n_sel bajo el mismo umbral.
9. Bins vacíos o deterministas (k ≥ m) que diluyen D_forma.
10. Con N fijo, m chico gana localidad pero pierde potencia frente a conos poliédricos; no está medido cuánto (4).

**Evaluación.** F6 tiene sentido para curvatura y bordes si se escala N con m, se calibra por (m, N, n_sel) y se controla la familia de m. Con N = 4000 y m = 20 fijo, bloque y FRW quedan casi sin intervalos locales y cualquier resultado sería difícil de leer. Decide Sebastían.

## Ronda 2

**(1) De acuerdo.** Errores entre realizaciones independientes (30 o más): el asesor de curvatura vio lo mismo que yo, el solapamiento fabrica «tendencias» de 2σ a 3σ en el plano (1). También con max-T sobre m con ≥ 150 calibraciones y con un solo criterio que vote. Usar r del interior sin las puntas evita mis bins degenerados, porque su esperanza no depende de m en el plano (3). Corrijo un punto mío: en de Sitter, el «sesgo de lugar» y el rango de densidad son efectos de la carta conforme y no sesgan con muestreo uniforme en volumen propio (2). En FRW sí sesgan, porque R cambia con el tiempo.

**(2) Discrepo.** Treinta realizaciones no bastan si el bloque y la caja siguen con N = 4000: la escalera en m solo llega a ≈ 0.012 N. Además, el exponente 0.60 de de Sitter fuerte puede incluir centrado y borde, no solo el término R² (2).

**(3) Sesgos de tamaño ajenos.**
- **Dilución q:** a m fijo, diluir agranda los intervalos. En el bloque, m_max baja a ≈ 0.012 qN (2). Con q = 0.5 y N = 4000, m = 20 casi desaparece, los centros se concentran más, sube la fracción recortada de la caja y sube el solapamiento (2). Una desviación que «crece al diluir» en la caja puede ser recorte, que tiene el mismo signo que de Sitter (+3.8% en relacionados, 1). Propongo comparar a igual m/m_max, con m_max medido desde el orden (el mayor intervalo), y excluir m > 0.1 m_max por una regla fijada antes.
- **q anidados:** los niveles de q de una realización están correlacionados. Los errores del ajuste b + A x^{1/2} + B x deben salir de remuestrear realizaciones completas.
- **Max-T:** se calibra con diamantes, que tienen 300 intervalos en cada m. El bloque tiene 34 o 4 en los m altos. Sin calibrar a igual n_sel, el m más ralo domina el máximo y rechaza positivas por ruido.

**(4) Pregunta al asesor de curvatura:** con N = 4000, ¿qué pares (m, q) cumplen m ≤ 0.1 m_max en el bloque y la caja? Si no queda ninguno, ¿el intercepto en regiones sin diamante grande exige N ≳ 16 000?
