# Ramas y preguntas nuevas tras F5′ y F7 (asesor explorador, 2026-10-14)

Niveles: 1 medido, 2 propuesta, 3 externo, 4 extrapolación.

## 1. La hipótesis completa

**En simple.** Debajo del espacio y el tiempo habría algo más básico: una red de relaciones «antes / después». Esa red puede organizarse de varias maneras (fases), como el agua puede ser hielo o vapor. En una de esas maneras, mirada de lejos, la red es indistinguible de un espacio-tiempo. Hoy se puede comprobar una parte: dada una red, ¿sabemos decir si está en la «fase espacio-tiempo»? Sigue abierta la otra: por qué la red estaría en esa fase y no en otra.

**En técnico.** Ω-0.2 se separa en tres afirmaciones:

- **H_rec (reconocer).** Existe un criterio intrínseco R, que solo usa el orden, tal que R = 1 si y solo si el orden se parece a un sprinkling de una variedad de dimensión d con conos suaves. Estado (nivel 1, intervalos, N ≤ 2400): firma espectral más perfil acepta diamantes 3D y 4D, de Sitter moderado y FRW de radiación a ciegas; rechaza órdenes producto, KR, percolación, capas y, con D_forma, conos con aristas. Falla con conos p = 1.5 y 3 (los acepta) y con de Sitter fuerte (lo rechaza).
- **H_sel (seleccionar).** Alguna dinámica o medida del sustrato tiene una fase robusta con R = 1. Estado: Ω-C no seleccionó geometría (nivel 1); en Ω-E la ley de área no es exclusiva de la geometría (nivel 1). Además, contando órdenes al azar, casi todos son de tipo KR (nivel 3, Kleitman-Rothschild). Abierta.
- **H_fen (desviarse).** En esa fase existe δg ≠ 0 observable. Lejana (nivel 4).

## 2. Ramas nuevas

**R-A. Perfil local por tamaños (F6 ampliada).** Medir D_forma en sub-intervalos de m ∈ {10, 20, 40, 80}. Idea (nivel 2): un cono deformado pero plano se ve igual a toda escala (achicar un cono L^p da otro cono L^p), así que su desviación no depende de m; la curvatura sí, porque crece con el tamaño del intervalo frente al radio de curvatura (nivel 3, fórmulas de volumen de diamantes pequeños). Ejemplo: una hoja arrugada se ve plana con lupa; un cuadrado se ve cuadrado a cualquier aumento.
Pregunta: ¿la desviación de de Sitter fuerte baja al achicar m y la de los conos con aristas no? Falla si las dos curvas son paralelas. Sesgo esperado: con m chico se pierde resolución y los conos con aristas se acercan al ruido (empeora la falla 1; tolerable solo si se firma R-E). Costo: el de F6 por cuatro tamaños, unos 20 a 40 minutos.

**R-B. Curvatura como señal.** Estimar la curvatura escalar con la inclinación del perfil local frente a m, usando una fórmula externa (sin ajustar a los datos). Pregunta: ¿el estimador ordena bien de Sitter normal < fuerte y da el signo y la variación correctos en FRW? Falla si no es monótono con la curvatura verdadera. Costo: bajo si reusa los datos de R-A.

**R-C. Mapas locales de dimensión y curvatura.** Asignar d y curvatura a cada sub-intervalo y ver cómo cambian a lo largo de la estructura (por ejemplo, por altura en la cadena más larga). En FRW la curvatura cambia con η de forma conocida: es un control a ciegas natural. Para Ω importa más: un sustrato podría tener **dominios** geométricos y no geométricos, como hielo flotando en agua. Pregunta: ¿el mapa de FRW sigue la ley conocida? Costo: bajo sobre R-A y R-B.

**R-D. Heterogeneidad como discriminador.** Comparar la varianza del perfil entre sub-intervalos con la de la referencia exacta de Poisson. Pregunta: ¿percolación y capas muestran exceso de varianza? Advertencia (nivel 2): los órdenes producto y los conos L^p son homogéneos, así que esta rama no los atrapa. Costo: casi cero con los datos de R-A.

**R-E. Clase ampliada «conos suaves».** Si no se exige separar p = 3 de p = 2, lo que se reconoce es «orden continuo de dimensión d con conos suaves y estrictamente convexos» (en la literatura, estructuras de tipo Finsler; nivel 3). Consecuencia honesta: la pregunta no desaparece, se traslada. La invariancia de Lorentz está acotada con mucha precisión por observaciones (nivel 3), y el argumento de que un sprinkling respeta Lorentz (Bombelli, Henson, Sorkin; nivel 3) no vale para conos L^p. Propuesta: firmar la clase ampliada como enmienda, dejar p = 1.5 y 3 como informativos en pruebas futuras y anotar «Lorentz exacto» como pregunta abierta para una etapa posterior. Costo: una decisión. Lo que **no** cambia: F7 sigue FALLANDO también por de Sitter fuerte.

**R-F. Puente con H_sel y Ω-E.** Usar el reconocedor como examen de dinámicas candidatas: crecimiento secuencial clásico (nivel 3, Rideout y Sorkin), acciones, pesos entrópicos. Medir también la firma de Ω-E en sub-intervalos. Pregunta: ¿alguna dinámica sin geometría impuesta produce órdenes que pasen R? Falla si ninguna lo hace en un catálogo preregistrado. Costo: alto (semanas). Es la rama que más vale para Ω, pero solo después de que R sea confiable.

## 3. Orden por valor frente a costo

| Orden | Rama | Valor para Ω | Costo | Comparte instrumento con F6 |
|---|---|---|---|---|
| 1 | R-E (decidir la clase) | alto: define los negativos | nulo | no aplica, va antes |
| 2 | R-A (tamaños) | alto: ataca la falla de curvatura | bajo | **sí, es F6** |
| 3 | R-B (curvatura) | alto: vuelve el ruido en medición | bajo | **sí** |
| 4 | R-D (heterogeneidad) | medio | casi nulo | **sí** |
| 5 | R-C (mapas) | medio ahora, alto para dominios | bajo | **sí** |
| 6 | R-F (dinámicas) | el mayor a largo plazo | alto | usa R como juez |

Recomendación (nivel 2): **una sola corrida preregistrada** de F6 por tamaños que guarde el perfil de cada sub-intervalo. Un único criterio vota (por ejemplo, D_forma local en el m más chico medible); R-B, R-C y R-D se informan sin voto y, si algo aparece, se prueba después con semillas nuevas.

## 4. Riesgos de sobreajuste

1. **Mover el umbral de F7 después de verlo.** Entre de Sitter fuerte (1.55) y los conos con aristas (2.65 y 3.05) cabe un umbral que «pasaría», pero se eligió mirando los datos. No vale como evidencia.
2. **Elegir el m que mejor funcione.** Fijar los m y el que vota antes de correr.
3. **Corregir la curvatura con de Sitter y luego aceptar a de Sitter.** La fórmula debe ser externa, y la validación en una curvatura no usada.
4. **Contar sub-intervalos como independientes.** Se solapan; el error real es mayor. Calcular errores por corrida, no por sub-intervalo.
5. **Sesgo de selección por tamaño.** Elegir intervalos con exactamente m elementos favorece regiones densas o de cierta forma; en la caja con bordes, los intervalos recortados cerca del borde pesan distinto. Medir y declarar.
6. **Rescatar F7 con la clase ampliada.** R-E vale hacia adelante; F7 queda como FALLA.
7. **Muchas métricas informadas.** Con cinco medidas sin voto, alguna parecerá discriminar por azar. Tratarlas solo como generadoras de hipótesis.
8. **Ajustar una dinámica hasta que pase R (R-F).** Fijar el catálogo de dinámicas y sus parámetros antes de aplicar el juez.

En una frase: la curvatura que hoy hace fallar al reconocedor puede volverse su mejor herramienta, si se mide en varios tamaños con reglas congeladas antes de mirar.

## Ronda 2

**(1) Acuerdos.**
- Los errores se calculan entre realizaciones, con 30 o más.
- Vota un solo criterio, con la rejilla de m (y de q) congelada.
- F7 queda como FALLA.

**(2) Discrepancias y correcciones a mi texto.**
- Mi R-A suponía que la desviación de un cono no depende de m. El asesor de curvatura me corrige: depende de m, pero no de la densidad. Conviene diluir.
- Mi rango m ∈ {10, …, 80} con N = 4000 no es local (tamaño finito: hace falta N ≳ 800 m).
- R-D por sí sola vale poco: diluir también cambia a los órdenes con escala propia.

**(3) Efecto de H_loc.** Sí cambia la articulación (nivel 2). H_rec pasa a ser operativa: «lo que se desvanece al diluir con exponente 2/d es curvatura; lo que queda (intercepto) es forma».
- Hay que cuidar una cosa: la d del exponente tiene que coincidir con la de Myrheim-Meyer. Sin esa condición, cualquier defecto que se diluya se haría pasar por curvatura.
- Cambio de orden: R-A y R-B se funden en una sola rama, porque el mismo ajuste da b (vota) y A (curvatura, informada).
- R-E sigue primero, porque fija δ y ε₀.
- R-C baja: sin N grande, el mapa es ruido.
- Sigue faltando lo que nadie midió: cuánto se desvían los conos a m chico. Si es del orden del ruido, H_loc cambia una falla por otra.

**(4) Ramas nuevas de la interacción (todas nivel 2).**
- **Signo dS/AdS como predicción ciega.** Congelar el signo esperado de A antes de correr. Barata y difícil de sobreajustar.
- **Positiva no conformemente plana (Kasner de vacío).** Como R_ab = 0, H_loc predice A ≈ 0 e intercepto nulo. En campo débil, la dependencia con el boost daría una «dirección de Ricci» medible.
- **Minkowski con links borrados.** Al diluir, los links de los sobrevivientes son en su mayoría relaciones antiguas no borradas, así que el defecto se diluye como la curvatura. Lo atraparía solo el exponente, que sería distinto de 2/d. Esto abre una pregunta propia de Ω (OMEGA_MASTER §4): ¿un defecto microscópico que desaparece al hacer grano grueso es un impostor, o un microestado distinto de la misma geometría? Propongo que vote como impostor y que quede además como pregunta abierta de H_sel.
