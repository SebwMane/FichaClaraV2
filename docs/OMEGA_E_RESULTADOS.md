# Ω-E: resultados de E0 y E1 (ruta entrópica, caso 3 de Ω-0.2)

- **Fecha:** 2026-10-12. **Commit de ejecución:** `54cf515` (todas las partes de E0b′ y E1). **Especificación:** `OMEGA_E_PROPUESTA.md` §1 a §4 con la enmienda E-A1. **Datos:** `results/omega_e/`.

## 1. E0: reproducir a Sorkin y Yazdi

| Prueba | Resultado |
|---|---|
| E0a, sin truncar | **PASA**: S = 0.315 n − 3.1 (publicado 0.32 n − 6.64); R² lineal 0.990 frente a 0.908 |
| E0b, truncado, primer intento (N_L 200 a 1600) | **FALLA**: a = 0.480 ± 0.057 |
| E0b′, truncado (N_L 800 a 2400, enmienda E-A1) | **PASA**: a = 0.347 ± 0.092 (publicado 0.346 ± 0.028), con ventaja pequeña del ajuste logarítmico (R² 0.384 frente a 0.361) |

**Lectura (nivel 1).** El instrumento reproduce la ley de volumen sin truncar y, con el doble truncamiento, el coeficiente logarítmico publicado. La reproducción del caso truncado tiene poca potencia estadística.

## 2. E1: ¿la ley de área distingue geometría? · **INDETERMINADO**

| Familia | ¿Ley de área? | Ajuste del caso truncado | Región en N_L = 2400 |
|---|---|---|---|
| Control geométrico | **No** | log a = 0.666; R² lineal 0.464 > log 0.422 | 600 elementos |
| KR | No (por poco) | log a = **0.436** (ventana hasta 0.43); R² log 0.989 > lineal 0.977 | 371 |
| Percolación p = 0.05 | No | entropía que **disminuye** con n (log a = −5.6) | 600 |
| Capas gruesas | No | ley de volumen: lineal a = 0.094, R² 0.996 | 600 |
| Bicapa | No aplicable | sin regiones extensas (intervalo máximo de 2) | 2 |
| Sándwich | No aplicable | sin regiones extensas (intervalo máximo de 4) | 4 |

(El campo `n_medio_1600` del summary corresponde en realidad a N_L = 2400, el último de la malla tras E-A1; error de nombre, no de cálculo.)

**Decisión preregistrada: INDETERMINADO**, porque el control geométrico no mostró ley de área con la regla intrínseca.

## 3. Por qué (diagnóstico declarado, no evidencia)

La regla «intervalo de tamaño más cercano a N/4» elige regiones del tamaño correcto pero **descentradas y alargadas**: en tres semillas con N_L = 1200, los rangos en coordenadas de cono de luz fueron, por ejemplo, u ∈ [−0.26, 0.45] y v ∈ [−0.11, 0.24], casi tocando el borde del diamante grande. E0 usaba un subdiamante centrado y cuadrado. El coeficiente ≈ 1/3 corresponde a esa forma, así que la ventana de E1 no aplicaba a las regiones que la regla producía. **Es un defecto de diseño del preregistro**, no del instrumento.

## 4. Observaciones que valen por sí mismas (nivel 1, descriptivas)

1. **KR casi pasa el criterio.** Su entropía truncada crece logarítmicamente, con mejor ajuste que el de la geometría y un coeficiente de 0.436, a 0.006 del borde de la ventana. Una ventana algo más ancha lo habría contado como «ley de área». En la forma probada, el «termómetro» no separa con claridad geometría de KR.
2. **Las familias se comportan distinto entre sí:** capas gruesas siguen ley de volumen incluso truncadas; la percolación da una entropía que baja con el tamaño; bicapa y sándwich no tienen regiones donde medir.
3. El umbral √n/(4π) proviene del espectro de un diamante 2D; en órdenes no geométricos no tiene justificación (declarado en §3 de la propuesta), y el caso KR muestra que puede producir apariencias de ley de área.

## 5. Siguiente paso posible (no preregistrado)

**E1′**: una regla de región intrínseca que controle la forma, no solo el tamaño (por ejemplo, el intervalo de tamaño N/4 cuyo punto medio, en el sentido de C0-4, esté más equilibrado, y con validación previa en sprinklings de que reproduce el coeficiente de E0b′). Y un criterio que no dependa solo de una ventana del coeficiente, porque KR quedó al borde.
