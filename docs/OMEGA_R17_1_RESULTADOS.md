# R17-1: candidato S_SJ^max (intervalo máximo). Resultados técnicos (2026-10-19)

- **Mandato:** R17-1 firmada (Rev. 17 §6): anticadena fuera del soporte; regiones sin modos con la variante (ii).
- **Tipo:** desarrollo técnico con N ≤ 64. **No es evidencia** de H_sel ni de H_E. El candidato está separado del E1′ firmado.
- **Código:** `omega/causality/sj_max.py`, `tests/test_sj_max.py` y `tools/r17_1_describe.py`. **Descripción:** `results/r17_1/r17_1_describe.json`, con semillas SeedSequence(20261101, spawn_key = (326, familia, N, réplica)) y 50 órdenes por familia.

## 1. El muestreador es correcto

- Sin la anticadena, el espacio sigue conectado con N = 4 y 5, en los dos espacios. Con E-b son 39 y 356 estados; con E-a, 23 y 119.
- El núcleo exacto de una propuesta cumple balance detallado con el peso exp(λ·S_SJ^max), con error de 1e-14.
- Las cadenas reproducen la distribución exacta: la distancia de variación total está entre 0.011 y 0.026, con λ = ±6. Con λ = ±1 el peso casi no se aparta del uniforme, porque con N ≤ 5 el valor de S^max solo va de 0 a 0.61.
- Las propuestas que llevan a la anticadena se rechazan sin evaluar S, y las cadenas nunca la visitan.
- La regla del modo mínimo no se activó en ningún estado: ni en la enumeración ni en la descripción.

## 2. El candidato no sirve como instrumento (hallazgo principal)

| Familia | N | \|R\| máx (mediana) | S^max (media ± de) | S E1′ (media ± de) |
|---|---|---|---|---|
| Sprinkling 2D | 32 | 21 | 0.45 ± 0.22 | 0.51 ± 0.07 |
| Sprinkling 2D | 64 | 50 | 0.42 ± 0.19 | 0.98 ± 0.29 |
| Altura 2 (la trampa de ES-0) | 32 | 2 | **0.95 ± 0.12** | sin región |
| Altura 2 | 64 | 2 | **1.27 ± 0.10** | sin región |
| KR de 3 capas | 32 | 10 | 0.16 ± 0.04 | 0.35 ± 0.11 |
| KR de 3 capas | 64 | 17 | 0.21 ± 0.03 | 0.30 ± 0.05 |
| Cadena total | 32 y 64 | todo el orden | 0 | 0.67 y 0.74 |

- **Los órdenes de altura 2 obtienen el S^max más alto de todas las familias,** y el valor crece con N. Su región es un solo link de 2 elementos. Con λ > 0, el peso empujaría las cadenas hacia la misma trampa que mostró ES-0, ahora favorecida en vez de neutra.
- **En los sprinklings 2D, la región ocupa casi todo el orden** (21 de 32 y 50 de 64). Una región cercana al todo está cerca de un estado puro, así que S^max es chico y casi no se correlaciona con E1′ (r = −0.07 y 0.21). En la cadena total la región es el orden completo y S^max = 0 exactamente.
- **El desempate de E1′ nunca actúa.** El intervalo más grande siempre va de un elemento mínimo a uno máximo, así que min(|pasado(a)|, |futuro(b)|) = 0 siempre y decide el desempate por índices. Lo verifica una prueba.
- **Lectura (nivel 2):** el intervalo máximo mide cuán lejos está la región más grande de cubrir todo el orden, no el entrelazamiento de una región de tamaño controlado. Premia los órdenes fragmentados. Como regla de región para el peso exp(λS), queda **descartado como candidato** por razones de instrumento, antes de muestrear nada.

## 3. Consecuencia para R17-2

R17-2 debía hacer templado paralelo sobre el estadístico de R17-1. Sobre S^max no tiene sentido: validaría un método sobre un peso que ya se sabe que favorece la trampa. **R17-2 queda detenido** hasta que haya un candidato de región que pase esta misma descripción. Esta decisión vuelve al usuario y al Consejo.

## 4. Decisiones de implementación registradas

1. Los intervalos se definen por la relación (c[a, b]), no por el orden de las etiquetas, para que funcione con etiquetas desordenadas.
2. Si el espectro de iΔ_R es nulo (|ν| máximo < 1e-12), la suma queda vacía: S = 0, con 0 modos. No ocurrió nunca.
3. El truncamiento global no cambia (√N/(4π)), y el modo mínimo solo actúa sobre la región.
4. El soporte se aplica con barridos envoltorio nuevos. `mcmc_sj.py` no cambió, y ES-0 sigue reproducible: con un soporte que lo acepta todo, los barridos nuevos dan exactamente los mismos estados y el mismo flujo aleatorio que los originales.
