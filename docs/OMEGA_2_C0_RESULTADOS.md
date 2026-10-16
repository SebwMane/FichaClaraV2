# Ω-C: resultados del bloque C0 (validación de instrumentos causales)

- **Fecha:** 2026-10-11.
- **Commit de ejecución:** `bdae5f8` (árbol limpio, entropía 20261011, 43 s).
- **Especificación congelada:** `docs/OMEGA_2_PROPUESTA_CAUSAL.md` §4 con las enmiendas C-A1 a C-A4 de §8, registradas antes de esta corrida.
- **Datos:** `results/c0_validation/summary.json` y `runs.jsonl`.
- **Alcance:** C0 valida instrumentos. No simula dinámica y no dice nada sobre si Ω-C produce geometría.

## 1. Resultados oficiales

| Prueba | Decisión | Cifras |
|---|---|---|
| **C0-1** Myrheim-Meyer | **ÉXITO** | N = 1024: d = 2 → 2.007 ± 0.036; d = 3 → 2.987 ± 0.039; d = 4 → 3.998 ± 0.085 (tolerancia 0.10) |
| **C0-2** número = volumen | **ÉXITO** | razones 0.997 (2D), 1.002 (3D), 1.014 (4D) (banda [0.98, 1.02]) |
| **C0-3** acción BD 2D | **ÉXITO** | N = 1024: media −37.8 ± 25.4 frente a 1 (dentro de 3·SE, criterio de baja potencia) |
| **C0-3 exacto** (C-A2) | **ÉXITO** | N = 1024: N₀, N₁, N₂ = 5653.9 / 4624.2 / 4125.3 frente a 5648.9 / 4632.4 / 4127.4 exactos (todas ≤ 1.5·SE) |
| **C0-4** discriminación | **FRACASO (preregistrado como esperado)** | sprinklings marcados 100% en d = 2, 3, 4; KR marcado 0%; **percolación marcada 100%** |

## 2. Lectura

**Nivel 1 (verificado sobre instrumentos).** Contando solo relaciones de orden se recupera la dimensión de sprinklings en 2, 3 y 4 dimensiones, el volumen de intervalos causales (la versión discreta de ∫ dᵈx √(−g) ≈ N/ρ) y la acción discreta 2D con su valor esperado exacto. El diccionario «orden más número igual a geometría» funciona como instrumento en este rango.

**Lo que C0-4 enseña.** La coherencia entre dos estimadores de dimensión no basta para reconocer un orden no geométrico. La percolación transitiva densa (p = 0.05, N = 1024) se presenta como un objeto coherente de dimensión ≈ 1.2 (d_MM 1.22, d_mid 1.16). Es el aviso de M§11 trasladado a órdenes: un número de dimensión, incluso confirmado por dos vías, no prueba geometría.

## 3. Diagnóstico sin voto para la compuerta (C0-D)

| Orden (N = 1024) | S_BD 2D medio |
|---|---|
| sprinkling 2D | −36 (esperanza exacta 1) |
| KR | **−261 039** |
| percolación p = 0.05 | +4 413 |

(La acción 2D sobre sprinklings 3D y 4D no tiene significado físico y no se interpreta.)

**Consecuencia para R3-2.** Si la compuerta L-2c usara el criterio del bloque L, «menor acción gana», o un peso euclídeo e^{−βS} con β > 0, los órdenes KR ganarían por cientos de miles de unidades. Se repetiría exactamente el modo de fallo de S0: un premio que favorece estructuras densas no geométricas. La literatura de conjuntos causales evita esto de dos maneras: restringiendo el espacio a 2-órdenes (donde KR no existe) o usando el peso lorentziano e^{iS}, donde una acción enorme produce oscilaciones que se cancelan. Esa elección es del usuario y está pendiente.

## 4. Decisión según el plan

- C0-1, C0-2 y C0-3 habilitan los instrumentos de dimensión, volumen y acción 2D para la rama Ω-C.
- C0-4 bloquea el uso de la bandera «compatible con variedad» como clasificador hasta que C0-4b pase.
- La rama Ω-C sigue **sin ratificar** (R3-1 a R3-3).

## 5. C0-4b: discriminador por perfil de abundancias (corrida oficial)

- **Commit de ejecución:** `caa6a2b` (árbol limpio, 49 s). Especificación: `OMEGA_2_PROPUESTA_CAUSAL.md` §9, congelada en `9c97244` antes de escribir código.
- **Decisión: ÉXITO.**

| Celda (N = 1024, 20 semillas) | d̂ | d_ref | D mediana [mín, máx] | Marcadas |
|---|---|---|---|---|
| sprinkling 2D | 2.00 | 2 | 0.033 [0.018, 0.056] | 100% |
| sprinkling 3D | 3.00 | 3 | 0.040 [0.027, 0.073] | 100% |
| sprinkling 4D | 3.97 | 4 | 0.071 [0.030, **0.208**] | 100% |
| KR | 2.38 | 2 | 8.44 | 0% |
| percolación p = 0.05 | 1.22 | 1 | 1.87 [1.81, 1.92] | 0% |
| percolación p = 0.01 (control nuevo) | 2.70 | 3 | 1.42 [1.39, 1.45] | 0% |
| capas gruesas, 8 capas (control nuevo) | 1.41 | 1 | 6.93 | 0% |

δ calibrados solo con sprinklings: δ₂ = 0.070, δ₃ = 0.101, δ₄ = 0.215 (δ₁ = δ₂).

**Lectura (nivel 1 sobre instrumentos).** Exigir el perfil completo de tamaños de intervalo de la dimensión estimada separa los cuatro controles de los positivos por un factor ≥ 6 en D. Los dos controles nuevos, que nadie había visto al diseñar el criterio, quedan igual de lejos que los conocidos, lo que reduce la sospecha de que el criterio se ajustó a la percolación con p = 0.05.

**Límites honestos.**
- En 4D el positivo más alejado quedó en 0.208 frente a δ₄ = 0.215. El clasificador funciona en 4D, pero sin holgura; con N distinto podría no hacerlo.
- Probado solo en N = 1024 y en espacio-tiempo plano. No dice nada sobre sprinklings curvos ni sobre órdenes que imiten deliberadamente el perfil.

**Consecuencia.** La bandera de C0-4b reemplaza a la de C0-4 como clasificador habilitado para la rama Ω-C (para L-3c, si se ratifica). La de C0-4 queda retirada.
