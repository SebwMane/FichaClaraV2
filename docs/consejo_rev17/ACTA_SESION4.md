# Consejo Ω, sesión 4: resultados ES-0 y decisiones (texto pegado por el usuario, 2026-10-19)

Registrado tal cual. Los dictámenes individuales citados (`work/consejo_omega_sesion_4/`) no están en este repositorio. La verificación del orquestador está en `docs/OMEGA_CONSEJO_REV17_ES0.md` §2.

---

**Mensaje del usuario que acompañó el acta:**

El Consejo deliberó con tres especialistas independientes y una ronda de crítica cruzada. La decisión principal es pausar ES-E1 tal como está diseñado. El piloto no ofrece una estimación fiable de selección: E-b no mezcló en N=64 y quedó bajo el criterio de ESS en N=32. Esto deja H_sel indeterminada; no la refuta.

Las tres decisiones del Consejo son:

1. Sin región: retirar la regla `S_SJ = 0`; no adoptar todavía ninguna opción como definición oficial. Excluir estados condicionaría la pregunta; elegir el intervalo máximo cambia el instrumento; normalizar por tamaño define otro ensamble. El máximo intervalo puede probarse como prototipo técnico, y cualquier normalización requiere referencia independiente y una nueva especificación.
2. Muestreo: probar templado paralelo en λ después de definir el estadístico y el soporte. No garantiza mezcla ni corrige la trampa por sí solo. Wang–Landau queda como alternativa posterior.
3. Escala: N≤64 solo para desarrollo técnico, sin afirmaciones sobre selección o ley de área. Antes de reabrir ES-E1, hay que encontrar tamaños donde coincidan validez del instrumento y mezcla práctica.

También corregí un dato de estado: R15-2 ya se cerró en `a7d1fd9` y R14-C v2 está preparada, pero no ejecutada. La portada de Rev.16 conserva algunas frases anteriores a ES-0; la nueva acta las señala para corregirlas.

---

# Consejo Ω, sesión 4: resultados ES-0 y decisiones

## Acta de deliberación

**Fecha:** 19 de octubre de 2026
**Mandato:** analizar los resultados oficiales de ES-0 y decidir cómo tratar los estados sin región, qué muestreador probar en órdenes generales y qué escalas permitir antes de ES-E1.
**Método:** tres dictámenes independientes (Ω-E/entropía, ensambles/MCMC y estadística/epistemología) seguidos de crítica cruzada y resolución del orquestador. No se corrieron nuevos experimentos ni se modificó FichaClaraV2.

## Resolución ejecutiva

**Pausar ES-E1 en su definición actual.** ES-0 detectó problemas de mezcla, una trampa creada por la convención `sin región → S_SJ = 0` y un rango de regiones demasiado pequeño para interpretar el estadístico como el instrumento de ley de área validado. El piloto no estima si Ω-E selecciona geometría: en las cadenas generales de N=64 no se pudo estimar la distribución objetivo.

El Consejo resuelve las tres preguntas nuevas así:

| Decisión | Resolución del Consejo |
|---|---|
| **1. Estados sin región** | **Invalidar `S_SJ = 0` y no firmar todavía ninguna opción (a), (b) o (c) como definición oficial.** La exclusión cambia la pregunta a una medida condicionada; usar el máximo intervalo cambia E1′ y no resuelve antichains ni modos insuficientes; dividir por E[S∣\|R\|] es inestable si la media es pequeña y cambia el ensamble. Se puede prototipar el máximo intervalo como alternativa técnica, sin convertirlo en la medida de ES-E1. Toda normalización/residualización será una nueva versión de R16-1 con referencia independiente y preregistro propio. |
| **2. Mezcla en órdenes generales** | **Probar templado paralelo en λ como primer prototipo canónico, después de definir una función objetivo total y válida.** No es una garantía ni una reparación de la trampa. Mantener Wang–Landau como alternativa/contraste posterior, no como solución ya disponible para SSEE en E-b. |
| **3. Escala** | **Limitar N≤64 al desarrollo técnico, sin evidencia de H_sel, H_E ni ley de área.** Pausar el estudio inferencial hasta encontrar una intersección entre validez instrumental, tamaño de región y mezcla. Priorizar una actualización más barata que preserve S_SJ; un proxy nuevo requiere una hipótesis nueva. |

Las firmas R16 originales permanecen registradas, pero ES-E1 no debe ejecutarse con la convención que ES-0 invalidó. Una modificación del estadístico o del soporte requiere una enmienda y firma nuevas; esta acta formula la recomendación del Consejo, no la firma del usuario.

## Resultados informativos de ES-0

El informe publicado marca ES-0 como **piloto oficial de factibilidad, no evidencia de H_sel ni H_E**. Las 26 partes se registran en el commit del código `331a105`, con árbol limpio y `official = true`; el informe, resumen y figura están asociados a `4302cc8`.

### Mezcla por espacio y tamaño

| Espacio | N | Resultado reportado | Lectura del Consejo |
|---|---:|---|---|
| E-a, 2-órdenes | 32 | λ=0 y negativo mezclan; con λ>0 una cadena se queda en S≈1.5, R̂_S=1.56. | Calibración útil, pero no completamente estable en la rejilla. |
| E-a, 2-órdenes | 64 | Mezcla en las tres celdas probadas. | Sigue siendo calibración restringida a 2D, no evidencia de selección general. |
| E-b, órdenes generales | 32 | "Casi": ESS de la fracción de relaciones 146; el informe estima que faltarían alrededor de 1 100–1 500 barridos en las celdas afectadas. | No alcanza el umbral prerregistrado; no contar como mezcla demostrada. |
| E-b, órdenes generales | 64 | No mezcla en ninguno de los tres λ; R̂ de relaciones 7.7 en λ=0 y hasta 17 en λ>0. | La probabilidad de fase es indeterminada en este tamaño/presupuesto. |

El control de cordura E-a con λ=0 da fracciones de relaciones 0.504 y 0.501 en N=32 y 64. El informe también registra que la enumeración exacta en cuatro elementos reproduce la distribución objetivo y detecta una elección errónea de λ. Esto verifica casos pequeños y la implementación examinada; no establece mezcla práctica ni validez para E-b.

### Trampa por ausencia de región

En E-b N=64, una cadena iniciada en el orden vacío alcanza órdenes de altura 2 que no tienen un intervalo en la ventana objetivo E1′. La implementación les asigna S_SJ=0. El informe observa que no sale de esa zona durante 200 barridos con λ=±2/σ₀; con λ<0, el valor cero queda favorecido frente a valores positivos típicos. En varias celdas, los movimientos dentro de la zona tienen ΔS=0; el 25% de las muestras finales de esas celdas carece de región. Este es un artefacto de la convención y una dificultad de mezcla, no una fase física que ES-0 haya descubierto.

### Tamaño del instrumento y costo

La región objetivo es E1′ **N/8 ± 5%**, no N/4. ES-0 obtiene |R| promedio aproximado de 3, 7–8 y 15 en N=32, 64 y 128. En varias celdas, |R| explica entre 0.65 y 0.93 de la varianza de S_SJ; una celda N=128 tiene solo tres muestras y su η² no es precisa. El umbral de truncamiento se calibró en regiones de cientos de elementos, así que N≤128 no valida el régimen instrumental que sustenta la ley de área.

El cálculo completo tarda 39 s por barrido en E-a N=128 y 7.2 s por barrido en E-b N=128, según el informe; 1 000 barridos E-b se estiman en unas dos horas por cadena. En E-a, el costo equivalente es de unas once horas. El informe concluye que llegar a N de miles, donde el instrumento tenía su régimen, no es viable con el recálculo completo actual. Más duración de cadena no corrige un estadístico fuera de su rango de calibración.

Hay además un detalle de periodicidad en E-a, λ=0: cada propuesta se acepta y cada barrido tiene un número par de transposiciones, por lo que el muestreo al final de cada barrido conserva la paridad de la permutación inicial. Los observables medidos no mostraron efecto en ES-0, pero una siguiente cadena debe usar paso perezoso o muestrear en instantes que no fijen esa paridad para cubrir el espacio completo.

La escala de λ también requiere congelarse mejor: σ₀ varía de 0.07 a 0.29 entre celdas del piloto, de modo que ±2/σ₀ equivale a valores absolutos de aproximadamente ±7 a ±28. No comparar esas celdas como si compartieran un mismo acoplamiento absoluto sin declarar la estandarización.

## Evaluación de las tres opciones

### (a) Excluir estados sin región

Esto define una distribución sobre el subconjunto de órdenes que sí tiene región admisible. La pregunta pasa a ser "¿qué fase aparece condicionada a que exista una región medible?", no cuánto peso recibe la geometría sobre todos los órdenes generales. Puede además quitar del control precisamente estructuras no geométricas. Si se reporta, debe ser un brazo condicional separado, y hay que probar conectividad y mezcla en el soporte restringido. No es la opción primaria del Consejo.

### (b) Elegir el intervalo más grande

Puede probarse como **prototipo del selector de región** para comparar cómo varían tamaño, forma y espectro. No es equivalente a E1′. El máximo puede tener tamaño muy distinto entre órdenes; en antichains no hay par comparable, y algunos tamaños pueden no dejar modos después del truncamiento. Por tanto, requiere reglas explícitas para antichains, empates y estados sin modos. El Consejo permite evaluarlo solo como desarrollo técnico; no lo adopta todavía como definición del peso.

### (c) Dividir por E[S_SJ∣|R|]

La esperanza condicional no es única hasta fijar la medida base, el selector de región, N y el origen de los datos de referencia. Si la media condicional es cercana a cero o cambia de signo, el cociente puede volverse inestable. Dividir no garantiza eliminar toda dependencia de |R|, por ejemplo diferencias de varianza dentro de cada tamaño. Además, esa transformación cambia los cocientes de peso y define un ensamble distinto. Un residuo condicionado o un puntaje estandarizado podría analizarse como otra variante, con referencia independiente y muestras separadas de calibración y evaluación. **El Consejo no firma (c) como corrección inmediata.**

Ninguna de estas tres opciones resuelve por sí sola la combinación de soporte total, significación espectral, tamaño de región y mezcla. Hasta definir un estadístico válido para todo el soporte elegido, la decisión correcta es pausar ES-E1 y retirar la convención `sin región → 0`.

## Muestreo: recomendación para el siguiente prototipo

El Consejo prefiere templado paralelo en λ frente a Wang–Landau como primera prueba canónica porque mantiene réplicas dirigidas a los ensambles π_i(C) ∝ ν(C)·e^{λ_i S(C)}. Para dos réplicas con estados x, y, un intercambio simétrico debe aceptar con

min{1, exp[(λ_i − λ_j)(S(y) − S(x))]}.

Antes de usarlo para estimar H_sel hay que definir S y el soporte, verificar la razón de intercambio y el balance local, fijar una escalera con solapamiento, y reportar tasas de intercambio, round trips de réplicas y mezcla en S, |R|, fracción de relaciones, altura y observables estructurales. Que las réplicas intercambien con frecuencia no basta: E-b N=64 tampoco mezcla en la referencia λ=0, y PT no ayuda con barreras que no se reflejen en S.

El `wanglandau.py` existente no es una solución lista para este problema: el código citado es para la acción BD 2D en 2-órdenes, no SSEE continuo en órdenes generales. Wang–Landau puede servir después como contraste de densidad de estados, una vez definidos bins/rango y validados errores y mezcla entre estructuras con S parecido.

## Corrección de estado documental

Rev.16 §7 en `4302cc8` corrige la información anterior: **R15-2 quedó cerrada en `a7d1fd9` y R14-C v2 está preparada, no ejecutada**. La versión anterior del acta que decía R15-2 pendiente quedó desactualizada.

Hay dos restos de texto antiguo que conviene limpiar en una actualización documental: la portada de Rev.16 aún dice que nada fue ejecutado, aunque §7 y el informe ES-0 registran firmas y corrida; el §1 todavía menciona |R|≈N/4, aunque §7 corrige a N/8 ± 5%. El PR #3 sigue abierto al corte consultado.

## Secuencia prescriptiva

1. Conservar ES-0 como resultado oficial de factibilidad, sin reetiquetar sus datos como evidencia.
2. Mantener intacta la pregunta H_sel sobre órdenes generales. No excluir estados sin región por conveniencia; si se explora esa opción, llamarla probabilidad condicionada.
3. Desarrollar y comparar una regla total de región (el máximo intervalo puede ser el primer prototipo técnico) y definir antichains y validez espectral sin imputar cero.
4. Investigar normalización/residualización solo como variante de ensamble. Fijar medida de referencia y estimarla en cadenas o conjuntos independientes, con validación separada. No elegir la variante que produzca mayor geometría en los mismos datos.
5. En N≤64, comprobar balance, paridad, soporte y templado paralelo; medir mezcla estructural desde vacío/antichain, cadena, sprinkling y KR. Registrar "no convergido" como resultado posible.
6. Buscar una aceleración de S_SJ que concuerde con el cálculo completo en tamaños compartidos. Si se propone un proxy, darle hipótesis y calibración propias.
7. Solo reabrir ES-E1 cuando existan tamaños donde el instrumento sea válido y las cadenas mezclen. Congelar la definición de fase, rango de λ, controles KR alineados con la medida de etiquetado natural y semillas nuevas en un preregistro separado.

## Límites

Los especialistas revisaron los documentos, el resumen oficial y el código identificado por los artefactos del commit; no reejecutaron las 26 partes ni volvieron a calcular R̂/ESS o costos. Las conclusiones cuantitativas son las reportadas por ES-0, no una réplica independiente. La falta de mezcla hace indeterminada la probabilidad de fase; el piloto no apoya ni refuta la hipótesis física de selección.

## Registro del Consejo

Los dictámenes individuales y las críticas cruzadas están en `work/consejo_omega_sesion_4/`. Esta acta contiene la resolución colectiva y los pasos recomendados; no se ejecutó una nueva corrida ni se alteró el PR.
