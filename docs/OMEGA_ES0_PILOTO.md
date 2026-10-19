# ES-0: piloto de factibilidad del ensamble entrópico de Ω-E (2026-10-19)

- **Mandato:** R16-3 firmada y orden de ejecutar ES-0 (Rev. 16 §7).
- **Tipo:** piloto declarado. **No es evidencia** de H_sel ni de H_E: mide costo y mezcla para decidir si la rejilla de ES-E1 es viable. Nada de lo que salga aquí se usa para fijar la definición de fase geométrica.
- **EXPERIMENT_ID 325; entropía 20261100.** spawn_key (325, espacio, N, índice de λ, cadena).

## 1. Ensamble

μ_λ(C) ∝ exp(λ·S_SJ(C)), con S_SJ el instrumento de Ω-E tal como lo usó E1″:
1. iΔ = (i/2)(C − Cᵀ), diagonalizada una vez;
2. truncamiento global: se quitan los modos con |μ| < √N/(4π) (c = 1);
3. región R = `interval_balanced(C, N // 8)` (N/8 ± 5%, E1′ §2);
4. S = `ssee` sobre R con truncamiento √|R|/(4π).

**Convención declarada para el piloto:** si C no tiene ningún intervalo del tamaño pedido, S_SJ(C) := 0, y se registra la fracción de muestras sin región. ES-E1 debe revisar esta convención con lo que muestre el piloto.

## 2. Espacios y movimientos

| Espacio | Estado | Movimiento (simétrico) | Inicios de las 4 cadenas |
|---|---|---|---|
| E-a: 2-órdenes | permutación v | transponer v_a y v_b | cadena total, anticadena, dos 2-órdenes al azar |
| E-b: órdenes generales con etiquetado natural fijo | matriz triangular superior | par i < j al azar: quitar si es link, agregar si deja el conjunto transitivo; si no se puede, la propuesta cuenta como rechazada | orden vacío, cadena total, sprinkling 2D, orden KR de 3 capas (N/4, N/2, N/4; capas 1 y 3 todas relacionadas; capas contiguas con probabilidad 1/2) |

Aceptación de Metropolis min(1, exp(λ ΔS)). S_SJ se recalcula completa en cada propuesta válida (no hay actualización incremental). Un **barrido** son N(N−1)/2 propuestas.

## 3. Rejilla del piloto

- **N = 32 y N = 64**, los dos espacios.
- **Etapa 1:** λ = 0. **Etapa 2:** λ = ±2/σ₀, con σ₀ la desviación estándar de S_SJ en la mitad final de las cadenas de la etapa 1 del mismo espacio y N (fija la unidad de λ; regla declarada).
- **Largo:** 400 barridos con N = 32 y 200 con N = 64, se mida cada barrido.
- **N = 128:** solo costo (una cadena por espacio, λ = 0, 5 barridos).

## 4. Qué se mide

- **Costo:** milisegundos por evaluación de S_SJ, segundos por barrido, fracción de propuestas válidas (E-b) y tasa de aceptación.
- **Trazas por barrido:** S_SJ, |R| (y si hay región), fracción de pares relacionados, número de links y altura (cadena más larga).
- **Mezcla,** sobre la mitad final de cada cadena: R̂ dividido entre las 4 cadenas, tiempo de autocorrelación integrado τ (ventana automática de Sokal) y tamaño efectivo de muestra ESS, para S_SJ, la fracción de relaciones y la altura.

## 5. Criterio de viabilidad (declarado antes de correr; decide el plan, no la hipótesis)

- Una celda (espacio, N, λ) **mezcla** si R̂ < 1.05 y ESS ≥ 200 en los tres observables.
- Si no mezcla, se informa cuántos barridos harían falta (200·τ, o «no estimable» si las cadenas no se juntan) y cuánto tiempo costarían, y se proponen alternativas (templado paralelo en λ, Wang-Landau en S_SJ, N menor).
- **Controles de cordura** (deben cumplirse para confiar en el código, no son resultados): en E-a con λ = 0 la fracción de relaciones debe rondar 1/2 (2-orden uniforme); el muestreador debe reproducir una distribución exacta en un caso enumerable (prueba unitaria).

## 6. Riesgos declarados

- **Regiones diminutas:** con N/8, |R| ≈ 4 en N = 32 y 8 en N = 64. El umbral √|R|/(4π) se calibró con regiones de cientos de elementos. El piloto informa la distribución de |R| y de S_SJ, y si S_SJ se reduce casi a una función de |R|.
- **Etiquetado:** en E-b con λ = 0 la medida no es uniforme sobre órdenes sin etiquetas (Rev. 16 §7, R16-2).

## 7. Herramienta

`omega/causality/mcmc_sj.py` (muestreador con peso genérico), `tools/es0_piloto.py` (ejecutor reanudable por partes), pruebas en `tests/test_mcmc_sj.py`. Agente de código con handoff y auditoría; corrida oficial sobre árbol limpio; tandas de menos de 2 horas.
