# ES-C: selección con la acción de Benincasa-Dowker euclídea (Ω-C). Preregistro (2026-10-19)

- **Mandato:** R16-4 (rama euclídea primero), Rev. 17 §8 (priorizar Ω-C; auditoría de √ε primero) y firma del usuario del 2026-10-19 («firmo las dos»): corrección de √ε congelada (commit cf128f9) y preregistro del peso.
- **Etapas:**
  - **ES-C0**, diagnóstico de mezcla: se ejecuta ahora y **no es evidencia**.
  - **ES-C1**, exploratoria, mide P(fase): solo si ES-C0 muestra mezcla y después de una revisión del Consejo, con preregistro propio que completa este documento.

## 0. Lo que ya se sabe y lo que ES-C agrega

- **P4** (`OMEGA_2_P4_RESULTADOS.md`, nivel 1) ya corrió la acción 4D **sin suavizar** sobre órdenes generales, con N = 64 y 128 y β entre −0.1 y 0.1. Resultado: **SIN VENTANA**. Partiendo de un sprinkling 4D, la geometría se pierde en los primeros 100 barridos con todos los β.
- **Pero P4 corrió una sola cadena por celda, sin diagnóstico de mezcla.** No se puede separar «la dinámica no selecciona geometría» de «la cadena no llegó al equilibrio».
- **Corrección de mi recomendación anterior:** te propuse empezar con ε = 1, sin recordar que P4 ya cubrió ese caso. Repetirlo sin más no agrega nada.
- **ES-C agrega dos cosas:**
  1. **mezcla verificada:** cuatro cadenas desde inicios opuestos, con R̂ y ESS;
  2. **la acción suavizada** con ε < 1, ya con la fórmula corregida, que P4 no probó en 4D. La versión suavizada reduce las fluctuaciones de la acción (Benincasa y Dowker 2010), así que puede cambiar el paisaje.

## 1. Ensamble

- **Espacio:** órdenes de N elementos con etiquetado natural fijo (E-b, igual que P4 y ES-0). Con b = 0, cada orden sin etiquetas pesa según su número de extensiones lineales (medida declarada, R16-2).
- **Peso:** π_b(C) ∝ exp(−(b/N)·S_ε(C)), con S_ε = √ε·[N − ε·Σ f₄(n, ε)] (fórmula congelada, `bd_action_4d_smeared`); con ε = 1 es la acción sin suavizar.
- **Movimiento:** el de `mcmc_general.py` (quitar un link, o agregar una relación que deje el conjunto transitivo). La actualización incremental de S_ε generaliza la verificada en `tests/test_bd_incremental.py`: los mismos pares cambian su número de elementos intermedios, con f₄(n, ε) en lugar de los coeficientes truncados.

## 2. Rejilla de ES-C0

- **ε ∈ {1, 0.25}.**
- **b ∈ {−12, −4, 0, 4, 12}.** Con N = 128 cubre los β de P4: ±0.03 y ±0.1 equivalen a b de ±3.8 y ±12.8.
- **N ∈ {32, 64, 128}.**
- **Cadenas:** 4 por celda, desde el orden vacío, la cadena total, un sprinkling 4D del diamante con etiquetado natural y KR (`orders.kr_order`, reetiquetado de forma natural).
- **Largo:** 2000 barridos con N = 32, 1000 con N = 64 y 500 con N = 128. Se mide en cada barrido.

## 3. Observables (independientes de la acción, salvo S como diagnóstico)

- **Bandera de fase 4D:** el clasificador de P4, recalibrado en cada N con entropía propia. Marca un orden si su dimensión de Myrheim-Meyer redondeada es 4 y su perfil de intervalos queda dentro de δ₄ (1.5 × el percentil 95 de 40 sprinklings 4D).
  - Se valida antes de mirar las cadenas: al menos 90% de 100 sprinklings 4D marcados, y como máximo 10% de 20 KR y de 20 bicapas.
  - Si no valida en algún N, ese N queda **INDETERMINADO**.
- **Otros observables:** fracción de pares relacionados, dimensión de Myrheim-Meyer, altura, links/N y S_ε.

## 4. Criterio de mezcla de ES-C0 (decide si ES-C1 es posible; no la hipótesis)

- **Una celda mezcla** si, en la mitad final de las cuatro cadenas, R̂ dividido < 1.05 y ESS ≥ 200 en S_ε, en la fracción de relaciones y en la altura. Se usa el mismo código de diagnóstico que ES-0.
- **Se informa además, sin voto:** la fracción marcada por cadena y el barrido en que la cadena del sprinkling 4D deja de estar marcada.
- **Si no mezcla:** se informa cuántos barridos harían falta, o «no estimable». Las alternativas (templado paralelo en b, Wang-Landau en S) quedan para una decisión posterior.

## 5. ES-C1 (solo se esboza; requiere completar el preregistro)

- Mide P_b,N,ε(bandera) con su intervalo por cadenas en las celdas que mezclen, y su tendencia con N.
- Lleva controles sellados nuevos (familias no geométricas definidas por un guardián) para medir los falsos positivos de la bandera antes de abrir P.
- **Predicción congelada (nivel 2), a partir de P4:**
  - con ε = 1, P ≈ 0 en todos los b;
  - con ε = 0.25, sin predicción.

## 6. Convenciones

- **EXPERIMENT_ID 329.**
- **Entropías:** 20261104 para las cadenas, 20261105 para la calibración del clasificador y 20261106 para su validación.
- **Próximos libres:** EXPERIMENT_ID 330 y entropía 20261107.
- **Herramienta:** agente de código con auditoría del orquestador; pruebas exactas de la actualización incremental de S_ε (enumeración con N ≤ 5 y balance detallado) antes de correr; corrida oficial sobre árbol limpio; tandas de menos de 2 horas.

## 7. Aclaración A1 (2026-10-19, antes de correr cualquier cadena)

- **Qué pasó:** al calibrar el clasificador con las semillas oficiales, N = 32 **no valida** (marca 82 de 100 sprinklings 4D, por debajo del 90%; KR y bicapas: 0%). N = 64 y 128 validan (98% y 99%; KR y bicapas: 0%).
- **Aclaración:** §3 dice que «ese N queda INDETERMINADO». Se aplica solo a la **bandera de fase**: con N = 32, las fracciones marcadas se informan como INDETERMINADAS y no se interpretan. Las cadenas de N = 32 sí corren, porque el criterio de mezcla de §4 usa S_ε, la fracción de relaciones y la altura, que no dependen del clasificador.
- **Por qué no es sobreajuste:** se decide sin haber visto ninguna cadena, y no cambia ninguna regla de mezcla ni ninguna predicción.
