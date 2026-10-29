# ES-0: resultados del piloto de factibilidad (2026-10-19)

- **Plan:** `OMEGA_ES0_PILOTO.md` (commit 8ec2cb7, antes del código). **Código:** commit 331a105, árbol limpio; las 26 partes llevan `official = true`.
- **Datos:** `results/omega_es0/` (partes por cadena y `summary.json`). **Figura:** `docs/figuras/es0_piloto.png`.
- **Tipo:** piloto declarado. **No es evidencia** de H_sel ni de H_E. Lo que sigue decide el plan de ES-E1, no la hipótesis.

## 1. Resultado en una frase

**La rejilla de ES-E1 tal como la diseñó la Rev. 16 no es viable.** Los 2-órdenes (calibración) mezclan bien con N = 32 y 64. Los órdenes generales (evidencia) casi mezclan con N = 32, pero con N = 64 no mezclan, y aparece una trampa causada por la convención «sin región, S = 0». Además, con estos N la región tiene entre 3 y 9 elementos, y S_SJ mide en buena parte el tamaño de la región, no el entrelazamiento.

## 2. Mezcla (R̂ < 1.05 y ESS ≥ 200 en S_SJ, fracción de relaciones y altura)

| Espacio | N | λ = 0 | λ = −2/σ₀ | λ = +2/σ₀ |
|---|---|---|---|---|
| E-a (2-órdenes) | 32 | mezcla | mezcla | **no** (una cadena salta a S ≈ 1.5 y se queda; R̂ de S = 1.56) |
| E-a | 64 | mezcla | mezcla | mezcla |
| E-b (órdenes generales) | 32 | casi (R̂ 1.05, ESS 146 en relaciones; harían falta unos 1100 barridos, 0.1 h de CPU) | casi (unos 1100 barridos) | casi (unos 1500 barridos) |
| E-b | 64 | **no** (R̂ 7.7 en relaciones; no estimable) | **no** (R̂ hasta 7.3) | **no** (R̂ hasta 17) |

**Control de cordura:** en E-a con λ = 0, la fracción de relaciones es 0.504 (N = 32) y 0.501 (N = 64), como debe ser en un 2-orden uniforme. La prueba de balance detallado por enumeración exacta pasa.

## 3. Lo que muestran las cadenas (figura)

- **Trampa sin región (E-b, N = 64, λ ≠ 0):** la cadena que empieza en el orden vacío llega a órdenes de altura 2 (dos capas), que no tienen ningún intervalo de unos 8 elementos. Ahí S_SJ vale 0 por convención y todos los movimientos dentro de esa zona tienen ΔS = 0. Con λ = 0 la cadena sale de ahí hacia el barrido 80. Con λ = ±2/σ₀ no sale en 200 barridos (el 25% de las muestras finales no tiene región). Con λ < 0 la zona sin región es además favorecida, porque S = 0 está por debajo de lo típico. **La convención del piloto crea una fase artificial**, y ES-E1 no puede usarla.
- **E-b con λ = 0 y N = 64:** las cuatro cadenas se estacionan en fracciones de relaciones distintas (0.26, 0.30, 0.38 y 0.39) y no se juntan en 200 barridos. La relajación es lenta, no solo por la trampa.
- **E-b con λ = +2/σ₀ y N = 64:** las cadenas encuentran órdenes más altos (altura 6) con S_SJ entre 3 y 4.6, cada una en un nivel distinto. Es el comportamiento típico de un paisaje con varios valles. Con 200 barridos no se puede decir cuál pesa más.

## 4. El instrumento con regiones diminutas

- |R| promedio: unos 3 con N = 32, unos 7 a 8 con N = 64 y unos 15 con N = 128.
- **S_SJ depende mucho del tamaño de la región en varias celdas:** la fracción de la varianza de S explicada por |R| (η²) es 0.93 en E-a con N = 32 y λ > 0, 0.65 en E-b con N = 32 y λ = 0, y 0.89 en E-a con N = 128 (con solo 3 muestras, poco fiable). En esas celdas, el peso exp(λS) premia sobre todo el tamaño de la región. En otras celdas η² es bajo (0.01 en E-a con N = 64).
- El umbral √|R|/(4π) se calibró con regiones de cientos de elementos. Con 3 a 15 elementos, el truncamiento deja muy pocos modos. **Con N ≤ 128, S_SJ no es todavía el objeto de Ω-E** (la ley de área de E0 se vio con N de miles).

## 5. Costo

| Espacio | N | ms por evaluación de S_SJ | s por barrido | Fracción de propuestas válidas |
|---|---|---|---|---|
| E-a | 32 / 64 / 128 | 0.30 / 0.96 / 4.7 | 0.15 / 2.0 / 39 | 1 |
| E-b | 32 / 64 / 128 | 0.31 / 0.92 / 5.0 | 0.07 / 0.91 / 7.2 | 0.42 / 0.48 / 0.18 |

Con N = 128, 1000 barridos de E-b son unas 2 horas por cadena, y en E-a unas 11 horas. Con N del orden de los miles, donde el instrumento tiene sentido, el recálculo completo (O(N³) por propuesta) lo hace imposible.

## 6. Notas declaradas

- **Paridad en E-a con λ = 0:** se aceptan todas las propuestas y cada barrido tiene un número par de transposiciones, así que al final de cada barrido la cadena está siempre en permutaciones de la misma paridad que su inicio. No afecta a los observables medidos (la fracción de relaciones da 0.50), pero ES-E1 debe medir en momentos que no fuercen la paridad o usar un paso perezoso.
- **Unidad de λ:** σ₀ sale de la etapa 1 y es distinta en cada celda (de 0.07 a 0.29), así que λ va de ±7 a ±28 en unidades absolutas.
- **Figura:** después de la corrida solo se corrigió la superposición del título (`tight_layout`), sin tocar datos ni análisis.

## 7. Qué implica para ES-E1 (para decidir; nada se adopta sin firma)

1. **La convención para órdenes sin región tiene que cambiar.** Opciones: (a) excluir esos estados del espacio (el peso es 0); (b) usar como región el intervalo más grande disponible, con su |R| registrado; (c) normalizar S por su valor esperado dado |R|, lo que además ataca el problema del §4. La (c) cambia la definición del peso de R16-1 y la firma el usuario.
2. **El tamaño de región y el N chico vuelven a S_SJ casi una medida de |R|.** O se acepta trabajar con un S normalizado por |R| (opción c), o se necesita N mucho mayor, lo que exige una actualización incremental de S_SJ o un estimador más barato (problema abierto).
3. **Mezcla en E-b:** con N = 64 hace falta un método mejor que Metropolis simple: templado paralelo en λ (las réplicas con λ = 0 sacan a las demás de las trampas) o Wang-Landau en S_SJ, como ya preveía la Rev. 16 §2.2.
4. **E-a (calibración) es viable** con N ≤ 64 tal como está.
