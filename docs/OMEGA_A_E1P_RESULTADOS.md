# A-E1′: auditoría del instrumento E1′ en su dominio y viabilidad de escala (2026-10-19)

- **Mandato:** Rev. 17 §7.2, plan firmado en §7.3. Es **descriptivo y sin pesos**: no es evidencia de H_sel ni de H_E.
- **Código:** commit 0129e1e (árbol limpio), `tools/a_e1p_auditoria.py` y `omega/causality/sj_multicut.py`. Evaluación exacta: |ΔS| ≤ 4e-12 frente al cálculo de referencia. EXPERIMENT_ID 328, entropía 20261103. 426 corridas.
- **Datos:** `results/omega_a_e1p/`. **Figura:** `docs/figuras/a_e1p_auditoria.png`.

## 1. Resultado en cuatro frases

1. **En sprinklings 2D, el instrumento solo se estabiliza con N ≥ 512.** Con el corte de E1′ (c = 1), el coeficiente de la ley de área baja de 0.98 a 0.61 y luego a 0.37 al pasar de la ventana 64 a 256, a la de 256 a 1024 y a la de 512 a 2048. Solo en esta última es compatible con 1/3.
2. **KR imita la ley de área con un solo par de modos.** En la ventana 512 a 2048 su coeficiente es 0.41 ± 0.01, indistinguible del de la geometría (diferencia +0.04 ± 0.05; confirma E1″). Pero su región conserva **siempre 2 modos**, mientras la geometría conserva de 22 a 48. El coeficiente no los separa; el número de modos los separa sin ambigüedad.
3. **En 3D y 4D, el instrumento no da ley de área logarítmica.** Sus coeficientes van de 3 a 12 y crecen con N. El instrumento es de 2D, como estaba declarado, y no reconoce la geometría de la dimensión que interesa al proyecto.
4. **La escala útil es incomputable en este entorno.** Con N = 512, 4 cadenas de 1000 barridos cuestan unas 950 horas de CPU en órdenes generales y unas 6300 en 2-órdenes. Con N = 1024, unas 19 000 y 127 000.

## 2. Detalle por familia (c = 1, medianas)

| Familia | N | ¿Hay región? | \|R\| | Modos de la región | S |
|---|---|---|---|---|---|
| Sprinkling 2D | 64 / 256 / 1024 / 2048 | 100% | 8 / 31 / 122 / 246 | 4 / 14 / 32 / 48 | 0.98 / 1.67 / 2.09 / 2.20 |
| KR de 3 capas | 64 a 2048 | 100% | 8 a 256 | **2 en todos los N** | 0.27 a 0.79 |
| Altura 2 | todos | **0%** | sin región | | |
| Cadena total | 64 a 2048 | 100% | 7 a 243 | 4 a 22 | 0.74 a 1.04 |
| Sprinkling 3D | 64 a 2048 | 100% | 8 a 250 | 4 a 71 | 1.06 a 13.7 |
| Sprinkling 4D | 64 a 2048 | 80% a 100% | 8 a 254 | 4 a 58 | 0.80 a 13.8 |

- **Coeficiente de la ley de área por ventana, sprinkling 2D:**

  | Corte | 64 a 256 | 256 a 1024 | 512 a 2048 |
  |---|---|---|---|
  | c = 0.5 | 2.77 | 2.15 | 1.20 |
  | c = 1 | 0.98 ± 0.06 | 0.61 ± 0.05 | **0.37 ± 0.05** |
  | c = 2 | 0.65 | 0.42 | 0.41 |

  El valor depende del corte y del tramo de N; el corte c = 2 se estabiliza antes.
- **La fracción de modos que conserva la región cae con N** (0.57 con N = 64 y 0.20 con N = 2048, para c = 1), porque el corte crece como √|R|.
- **La región de KR empieza en un elemento mínimo** (fracción de pasado 0): el desempate de E1′ no logra centrarla.
- **Margen mínimo al umbral:** 8.6e-05. Ningún caso necesitó el cálculo de respaldo.

## 3. Escala

- **Costo medido por evaluación (sprinkling 2D, un hilo):** 0.7, 2.2, 7.6, 43, 217 y 1481 ms con N = 64, 128, 256, 512, 1024 y 2048. Ajuste: crece como N^2.2.
- **N mínima según el criterio de §7.2:** 512, con la ventana de 512 a 2048. Es la primera donde el coeficiente con c = 1 queda a menos de 0.1 de 1/3 con error menor que 0.05 y la región conserva al menos 10 modos.
- **Costo de cadenas en esa escala** (4 cadenas × 1000 barridos):

  | N | Órdenes generales (supuesto: 15% de propuestas válidas) | 2-órdenes |
  |---|---|---|
  | 512 | 946 h de CPU | 6309 h |
  | 1024 | 18 979 h | 126 524 h |

- **Registro de factibilidad:** la validez del instrumento (N ≥ 512) y la mezcla práctica (N ≤ 64 a 128 en ES-0) **no se solapan** con los recursos de este entorno. Hay más de un orden de magnitud de diferencia, aun con la aceleración exacta de R17-3.

## 4. Lectura para el paso 3 (nivel 2; la decisión es del usuario)

- **El instrumento de Ω-E es de 2D y en su dominio no separa la geometría de KR por el coeficiente.** Un peso exp(λ·S) basado en él no tiene cómo preferir la geometría frente a KR. Para H_sel en 4D tampoco reconoce la geometría 4D.
- **Observación nueva (exploratoria, se congela solo como hipótesis):** el **número de modos de la región** separa por completo a KR (2 modos) de las geometrías (de 14 a 71 con N ≥ 256). No se usa para nada sin un preregistro nuevo con controles nuevos.
- **Opciones que quedan sobre la mesa:** (a) restringir Ω-E a reconocer sprinklings 2D, sin papel en la selección; (b) proponer un observable entrópico nuevo, con hipótesis y calibración propias, que cuente modos o use la dimensión correcta; (c) priorizar Ω-C (Benincasa-Dowker euclídea, R16-4), cuyo instrumento sí tiene versión 4D y actualización incremental exacta (`mcmc_general.py`).
