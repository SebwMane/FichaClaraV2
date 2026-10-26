# R23: piloto (a), descendientes de órdenes no geométricos frente a un blanco 3+1. Preregistro exploratorio (2026-10-26)

- **Mandato:** R22-4 firmada (`OMEGA_CONSEJO_REV22_R20.md` §6): preparar (a) como piloto exploratorio con blanco en un diamante de Minkowski 3+1. El diseño sigue la fase 4 del dictamen de la sesión 11 y los refinamientos 1 a 6 de la Rev. 22.
- **Tipo:** **exploratorio, sin voto.** Sirve para ajustar el instrumento y los márgenes. **No es evidencia confirmatoria de H_E**, cualquiera sea el resultado. No necesita guardián, porque no hay veredicto que proteger; sí preregistro antes del código, código con pruebas, auditoría y corrida sobre árbol limpio en un solo commit.
- **Estado:** **PREPARADO, no ejecutado.** Para correr hace falta la firma del usuario.

## 1. Pregunta

**¿El descendiente Φ_q(X) de un orden de una familia dada se acerca, en distancia absoluta, a un blanco geométrico independiente, más de lo que se acerca un orden nuevo de la misma familia y del mismo tamaño?**

- Si se acerca y además mejora frente al control nuevo, el patrón es **compatible con acercamiento**: el orden «se ve geométrico» a escala gruesa por algo que la dilución hace y que el tamaño solo no hace.
- Si descendiente y control se mueven igual, el cambio es un **efecto de tamaño o de familia**, no de grano grueso.
- Un L < 0 relativo, como el de R20, no basta: puede reflejar tamaño o ensamble. Por eso el estimando principal es la distancia **absoluta** al blanco.

## 2. Blanco geométrico

- **Ensemble:** sprinklings de Poisson en un diamante causal de Minkowski 3+1 (`r14b.GEN["plana"]`, el mismo generador de las referencias de R20), con N ~ Binomial(N₀, q) y N₀ = 8000.
- **Instrumento (actual, permitido para un piloto por la salida del gate de R22):**
  - **vector común de 13 componentes en los 4 q** (frac_rel, r₁₀, r₂₀, perfil N_k/N_pares con k de 0 a 8 y √gemelos; sin r₄₀);
  - media, escala con pisos y Ledoit-Wolf con las 100 referencias de **ajuste** de R20 de cada q;
  - radio conforme (rango 96 de 100) con las 100 de **calibración** de R20 de cada q.
- **Conjunto nulo nuevo e independiente:** 100 sprinklings planos nuevos por q (entropía 20261118). Nunca se usan para ajustar nada. Dan la distribución nula de D con la que se compara cada brazo.
- **Componentes no medibles en un brazo:** D se calcula con la banda restringida a las componentes finitas **en los dos brazos del par** (descendiente y control); con esas mismas componentes se recalculan el radio conforme y las D del conjunto nulo (como `r19.pair_D`). Si quedan menos de 10 componentes, el par es **no evaluable** en ese q. Se informa siempre qué componentes faltaron y con qué frecuencia.

## 3. Familias (parámetros congelados) y predicciones

| Familia | Generador (código existente) | Papel | Predicción con q = 1/8, declarada antes |
|---|---|---|---|
| plana | `r14b.GEN["plana"]` | control nulo | no visible al inicio; descendiente y control en la nula |
| D5, gemelos al 10% | `f6_intento2.d5_gemelos` | control positivo | compatible con acercamiento: AUC cerca de 0.58 (como en R20) y contraste negativo |
| C5, 13% de links borrados | `f6_confirmacion.c5_links_borrados` | **pregunta principal** | sin predicción de clase. Los pares borrados son links de la escala de N₀; diluido, el defecto por par es menor que en un C5 nuevo de N_q, así que un acercamiento es posible |
| D4, 10% de relaciones agregadas cerca del cono | `f6_intento2.d4_agregadas` | pregunta secundaria | sin predicción de clase. El defecto está a la escala ℓ de N₀; R14-B sugiere que sigue fuera de banda con q = 1/8 |
| KR (3 capas, p = 1/2) | `orders.kr_order` | no geométrico | lejos; contraste cerca de 0 (la altura ≤ 3 se conserva al diluir) |
| percolación transitiva, p = 0.05 fijo | `orders.transitive_percolation` | no geométrico | lejos; sin predicción de contraste (con p fijo, cambiar N cambia la familia efectiva) |
| producto de 4 órdenes totales | `f1f2.random_dorder(4)` | no geométrico, autosemejante | lejos; contraste igual a 0 (diluir un d-orden aleatorio da otro d-orden aleatorio) |
| capas gruesas (6 capas, p = 1/2) | `orders.thick_layer_order(n, 6, 0.5)` | no geométrico | lejos; contraste cerca de 0 |

- **No se sustituye ninguna familia después de ver datos.** Una familia no visible al inicio se informa como «no aplicable», no se reemplaza.

## 4. Brazos, unidad y muestreo

- **Padre:** X, un orden de N₀ = 8000 de la familia, generado con semillas propias.
- **Descendiente:** Y_q = Φ_q(X) (`r14b.dilute`), con q ∈ {1, 1/2, 1/4, 1/8}.
- **Control nuevo:** Z_q, un orden nuevo de la misma familia y los mismos parámetros, con **N = |Y_q|** exacto. No es otra máscara del mismo padre. Con q = 1, Z₁ es un segundo orden independiente de tamaño N₀.
- **Unidad independiente: el padre.** Los 4 q de un padre son medidas repetidas. **30 padres por familia.** El bootstrap remuestrea padres completos (los 4 q y los dos brazos juntos), con 2000 réplicas. El conjunto nulo queda fijo, así que los IC son condicionales al instrumento y al conjunto nulo.
- **Precisión esperada:** con 30 unidades, los IC de R20 tuvieron semiancho de 0.05 a 0.3 en ln D. Para el AUC frente a 100 nulos, el semiancho esperado es de 0.1 a 0.15 cerca de 0.6. Alcanza para distinguir «cerca» (≤ 0.75) de «lejos» (≥ 0.95), no para afinar dentro de «cerca».

## 5. Estadísticos

Para cada padre p, cada q y cada brazo B ∈ {Y, Z}:
- **D_B** = distancia de Mahalanobis a la referencia plana de ese q dividida por el radio (§2).
- **u_B** = fracción del conjunto nulo con D menor que D_B (empates cuentan la mitad). Es la distancia expresada como rango contra la plana, y **no tiene techo** (refinamiento 3 de la Rev. 22).

Por familia y q:
- **Proximidad:** AUC_B = media de u_B sobre los padres. Si B fuera plano, vale 0.5.
- **Contraste:** Δ = mediana sobre padres de ln(D_Y/D_Z), emparejados por padre.
- Se informan también la mediana de D_Y y de D_Z, la fracción en banda y las componentes faltantes.

## 6. Clasificación descriptiva (márgenes provisionales)

Los márgenes salen de datos previos de R20 (AUC de las copias diluidas con q = 1/8 contra la plana: entre 0.58 y 0.79; de Sitter, 0.79; controles visibles, de 0.95 a 1). **Son provisionales y el piloto sirve para evaluarlos.** Los márgenes de R20 (±ln 1.3 y su cobertura) no se trasladan como garantía; ln 1.3 se usa solo como referencia descriptiva para el contraste.

| Etiqueta | Condición (IC 95% bootstrap por padre) |
|---|---|
| Visible al inicio | IC de AUC_Y con q = 1 entero ≥ 0.95 |
| Cerca con q = 1/8 | IC de AUC_Y con q = 1/8 entero ≤ 0.75 |
| Lejos con q = 1/8 | IC de AUC_Y con q = 1/8 entero ≥ 0.95 |
| Mejora frente al control | IC de Δ con q = 1/8 entero < −ln 1.3 |
| Igual que el control | IC de Δ con q = 1/8 dentro de ±ln 1.3 |

**Lectura por familia:**
- **Compatible con acercamiento:** visible al inicio, cerca y mejora.
- **Efecto de tamaño o familia:** visible al inicio e igual que el control (se informa si queda cerca o lejos).
- **Sin acercamiento:** visible al inicio, lejos y sin mejora.
- **No aplicable:** no visible al inicio.
- **INDETERMINADO:** cualquier otra combinación, o falta de resolución (pares no evaluables en más del 20% de los padres).

Siempre se informan los cuatro q. **No se infiere velocidad de decaimiento** de estas cifras.

## 7. Lo que el piloto puede y no puede decir

- **Puede:** describir cómo Φ_q mueve familias conocidas respecto de un blanco 3+1, comparar con el control nuevo, evaluar si los márgenes provisionales separan los controles (plana y D5) y preparar una confirmatoria.
- **No puede:** confirmar H_E; decir que Ω produce estos órdenes (H_sel); extender nada a otros N₀, otras dimensiones u otros blancos (de Sitter y FRW serían blancos separados).
- **Si una familia queda «compatible con acercamiento»,** el siguiente paso sería una confirmatoria con el instrumento congelado de una sola banda (salida (ii) del gate), un guardián y familias nuevas, no reinterpretar este piloto.

## 8. Convenciones

- **EXPERIMENT_ID 334.**
- **Entropías:**
  - 20261118: conjunto nulo plano (100 por q);
  - 20261119: padres, diluciones, controles nuevos, intervalos de medición y bootstrap.
- **Próximos libres:** EXPERIMENT_ID 335 y entropía 20261120.
- **Código:** `tools/r23.py` y `tests/test_r23.py`. Reutiliza `r14b` (Φ_q, `measure`, `obs_vector`, `Band`), `r20` (radio conforme, cuantiles extendidos) y los generadores citados, sin modificarlos. Pruebas mínimas:
  1. cada generador devuelve exactamente N elementos y un orden parcial;
  2. el control nuevo tiene N = |Y_q|;
  3. u y el AUC dan 0.5 en promedio con datos nulos sintéticos;
  4. la regla de componentes por par reproduce el D de la máscara completa cuando no falta nada.
- **Costo estimado:** 400 conjuntos nulos (unos 10 min de pared) más 8 familias × 30 padres × unos 50 s, en total unos 50 a 60 min de pared con 4 núcleos, en tandas de menos de 2 horas.

## 9. Aclaración A1 (2026-10-26, después de escribir y auditar el código, antes de correr)

Ningún dato del piloto se generó todavía. Las mediciones de tiempo y de componentes usaron una entropía de prueba, no 20261119.

1. **Costo corregido.** Un padre cuesta en promedio unos 123 s, no 50: generar el control Z con q = 1 es otro orden completo de N₀, y C5 y D4 tardan en construirse (de 158 a 256 s por padre). Total estimado: **unas 2.2 a 2.3 h de pared con 4 núcleos**, en dos tandas de menos de 2 h. El código reanuda.
2. **KR y capas_6 se miden con 11 componentes.** No tienen ningún intervalo con exactamente 10 o 20 elementos interiores: su estructura en capas solo produce intervalos grandes. Pasa en el 100% de las muestras, con N = 8000 y con N ≈ 1000. Por la regla de §2, sus pares usan |S| = 11, con la banda, el radio y la nula recalculados en esas 11. **Consecuencia declarada:** la ausencia de r₁₀ y r₂₀ es, en sí misma, una señal no geométrica que el instrumento descarta, así que la distancia de esas dos familias queda subestimada. Las demás familias tienen las 13 componentes.
3. **La regla de |S| < 10 no puede activarse con este vector:** frac_rel, el perfil y √gemelos siempre son finitos, así que |S| ≥ 11. Por eso la vía «INDETERMINADO por más del 20% de no evaluables» no se activará en la práctica. Se evalúa con q = 1 y con q = 1/8, los dos q que usa la lectura.
4. **Lecturas superpuestas.** «Visible, lejos e igual que el control» cumple a la vez «efecto de tamaño o familia» y «sin acercamiento». Se informa la primera, en el orden de la tabla de §6, y la segunda queda registrada como «también cumple». Es el caso que se predijo para KR, producto_4 y capas_6.
5. **Detalles fijados por el código:** el bootstrap de cada familia usa la clave (334, 7, f); el margen de «igual» es un intervalo cerrado [−ln 1.3, ln 1.3] y el de «mejora» es estricto; AUC y Δ se calculan sobre los padres evaluables de cada q. Las 100 referencias de validación de R20 no se usan.
6. **Código:** `tools/r23.py` y `tests/test_r23.py` (11 pruebas pasan), escritos por el agente de código y auditados por el orquestador.
