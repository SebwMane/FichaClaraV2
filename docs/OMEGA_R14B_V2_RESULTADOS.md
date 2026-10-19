# R14-B v2: respuesta de órdenes dados a la dilución aleatoria. Resultados (2026-10-19)

- **Preregistro:** `OMEGA_R14B_V2_PREREGISTRO.md`, con la aclaración A1 (§8, antes de mirar datos).
- **Código:** commit 73e4925, árbol limpio. Son 1920 referencias y 170 realizaciones, todas con el mismo commit.
- **Datos:** `results/omega_r14b/` (`summary.json`, `tabla_bandas.csv`, `tabla_trayectorias.csv`). **Figura:** `docs/figuras/r14b_v2.png`.
- **Tipo:** exploratorio, **sin voto**. Todas las estructuras son vistas.
- **Incidente de ejecución:** la segunda tanda no arrancó, porque los archivos de R14-C, todavía sin subir, ensuciaban el árbol. Se apartaron sin cambios y se relanzó la misma tanda sobre el mismo commit. Ningún dato se perdió ni se mezcló.

## 1. Trayectorias (mediana de D/radio sobre 10 realizaciones; banda = 1)

| Estructura | Objetivo | q = 1 | 1/2 | 1/4 | 1/8 | Etiqueta |
|---|---|---|---|---|---|---|
| plana | plana | 0.91 | 0.75 | 0.58 | 0.73 | converge a banda |
| de Sitter | propia | 0.62 | 0.69 | 0.85 | 1.03 | en banda microscópica |
| de Sitter fuerte | propia | 0.82 | 0.89 | 0.74 | 0.81 | converge a banda |
| anti de Sitter | propia | 0.76 | 0.62 | 0.85 | 0.61 | en banda microscópica |
| FRW polvo | propia | 0.86 | 0.68 | 0.77 | 0.73 | en banda microscópica |
| D2 (FRW exponencial) | propia | 0.82 | 0.67 | 0.83 | 0.82 | en banda microscópica |
| C0 (Kasner) | propia | 0.84 | 0.76 | 0.69 | 0.74 | en banda microscópica |
| D0 (Kasner) | propia | 0.65 | 0.69 | 0.68 | 0.88 | en banda microscópica |
| C5 (13% de links borrados) | plana | 4.33 | 2.12 | 0.94 | 0.95 | indeterminada |
| D4 (10% de relaciones agregadas) | plana | 30.2 | 10.3 | 4.29 | 2.27 | no converge |
| D5 (10% de gemelos) | plana | 5.40 | 1.71 | 0.66 | 0.74 | converge a banda |
| D6 (cono con tapa plana) | plana | 1.05 | 0.64 | 0.63 | 0.80 | converge a banda |
| D7 (cono en huso) | plana | 1.91 | 1.02 | 0.65 | 0.74 | converge a banda |
| D8 (cono suma de normas) | plana | 1.99 | 1.09 | 0.66 | 0.75 | converge a banda |
| C4 (cono esfera y cubo) | plana | 5.22 | 2.12 | 0.85 | 0.89 | indeterminada |
| cono octaédrico | plana | 13.9 | 5.08 | 2.15 | 1.48 | no converge |
| producto 4 | plana | 24.1 | 10.3 | 5.12 | 3.65 | no converge |

## 2. Predicciones congeladas (v1 §1)

| Predicción | Resultado |
|---|---|
| Las geometrías quedan en banda en todo q | **Sí en lo esencial.** Ninguna geometría queda fuera de banda en ningún q. En 4 de 8, algún q queda ambiguo (el IC toca el borde): por eso plana y de Sitter fuerte se etiquetan «converge a banda» y no «en banda microscópica». |
| D5 converge a banda | **Se cumple.** |
| Los conos no convergen | **Falla en 4 de 5.** Solo el octaédrico no converge. D6, D7 y D8 convergen y C4 queda indeterminado. |

## 3. Lectura (nivel 2): la dilución mezcla convergencia con pérdida de resolución

- **El punto central:** al diluir, N efectivo baja de 8000 a 1000 y las bandas de referencia se ensanchan, porque el ruido crece como 1/√N. Un defecto de tamaño **fijo** se vuelve indistinguible solo por eso. Con 8 veces menos elementos, un efecto fijo medido en desviaciones estándar se achica unas √8 ≈ 2.8 veces.
- **Por eso «converge a banda» no prueba por sí solo H_E.** Los conos D6 a D8, que por autosemejanza no deberían cambiar al diluir, entran en banda porque la referencia deja de verlos, no porque se vuelvan geométricos. Es la misma lección de R14-A: con N chico, el instrumento no ve los conos «redondos casi en todas partes».
- **Lo que sí distingue algo:** cuánto cae D/radio de q = 1 a q = 1/8, comparado con ese factor de unas 2.8 veces:

  | Estructura | Caída de q = 1 a q = 1/8 |
  |---|---|
  | D7 | 2.6× |
  | D8 | 2.7× |
  | D6 | 1.3× (ya estaba casi en banda) |
  | C5 | 4.6× |
  | C4 | 5.9× |
  | producto 4 | 6.6× |
  | D5 | 7.3× |
  | octaédrico | 9.4× |
  | D4 | 13× |

  - D5 cae mucho más que el factor de ruido. Es compatible con la predicción: los gemelos sobreviven con probabilidad q², así que el defecto se desvanece de verdad.
  - D7 y D8 caen como el ruido. Es compatible con un defecto fijo que se deja de ver.
  - D4, el octaédrico y el producto caen más que el ruido, pero siguen fuera de banda.
  - **Cautela:** la distancia de Mahalanobis combina 14 componentes, así que la escala con √N es solo aproximada. Esta comparación es exploratoria y no se votó.
- **Calibración de las bandas:** la cobertura del percentil 95 de las 40 referencias de ajuste sobre las 20 de validación fue de 0.50 a 0.95 (muchas veces cerca de 0.75, como predecía la simulación con p = 14 y n = 40). La regla declarada (la 19.ª de las 20 de validación) corrige eso, pero con solo 20 referencias el borde de la banda es ruidoso.

## 4. Qué implica

1. **Para separar H_E de una pérdida de resolución, el mapa de grano grueso necesita un control.** Una opción es una referencia de «defecto fijo», o comparar la caída observada con la caída esperada solo por ruido, calibrada con la misma estructura sin diluir y submuestreada. Sin ese control, «converge a banda» es ambiguo.
2. **D5 (gemelos) es el único caso cuya caída supera claramente el ruido con un mecanismo que la explica.** Es el mejor candidato a «microestado compatible con H_E y no con H_L», pero es una observación exploratoria. Confirmarlo exige un preregistro nuevo con el control del punto 1.
