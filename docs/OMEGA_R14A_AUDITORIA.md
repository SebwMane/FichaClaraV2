# R14-A: auditoría de inferencia y trazabilidad de F6 (2026-10-18)

- **Mandato:** Rev. 14, R14-2 firmada por el usuario. Auditoría analítica **sin corridas nuevas y sin voto**: no cambia ningún veredicto de F6.
- **Fuente:** `results/omega_f6_intento2/` (partes del commit 1619162). **Herramienta:** `tools/r14a_auditoria.py`. **Salida:** `results/omega_r14a/summary.json`.
- **Azar:** solo un bootstrap de grupos, con la entropía nueva 20261087. No se genera ninguna estructura.

## 1. Unidad independiente (A1)

Por celda (m, N), en los 10 grupos planos de calibración, se compara la varianza de r **entre realizaciones** con la varianza «ingenua» que darían 300 intervalos independientes (r_sd²/n_sel).

| N \\ m | 10 | 20 | 40 | 80 | 160 |
|---|---|---|---|---|---|
| 4000 | 1.24 | 1.38 | 1.44 | 2.32 | **4.34** |
| 8000 | 1.20 | 1.14 | 1.05 | 1.41 | 1.73 |

(Medianas de la razón entre varianzas; rango por grupo en el JSON.)

**Lectura (nivel 1):** el solapamiento de los intervalos infla la varianza hasta 4.3 veces, más cuanto mayor es m/N. Tratar los intervalos como réplicas habría subestimado los errores hasta un factor 2 en desviación estándar. **F6 usó la realización como unidad, que era lo correcto;** la auditoría lo confirma con números.

## 2. Covarianza entre escalas (A2)

Correlación entre realizaciones de Δr en distintos m del mismo N (planos, centrado por grupo):
- **N = 4000:** de 0.15 a 0.58. La mayor es entre m = 80 y m = 160.
- **N = 8000:** de 0.03 a 0.41.

**Lectura:**
- **Ajuste lineal (nivel 1):** las celdas no son independientes. Trata cada celda como un punto, pero su intervalo de confianza venía de remuestrear realizaciones completas, lo que conserva esta correlación. El IC es válido; la estimación puntual no es la más eficiente.
- **K (nivel 2):** los pares 3 y 4 comparten realizaciones de N = 8000 correlacionadas, así que K no es exactamente χ²₄.

## 3. Distribución nula de K (A3)

- Los 15 valores de calibración van de 0.75 a 7.63, con media 3.25. La media esperada con χ²₄ es 4.
- Un Kolmogorov-Smirnov frente a χ²₄ da p = 0.40.
- Con 15 grupos, el máximo esperado de χ²₄ es de unos 9.7, y τ_K = 11.45 queda cerca.

**Lectura (nivel 1):** pese a la correlación, la calibración es compatible con χ²₄. La regla 1.5 × max es algo más permisiva que un percentil 95 nominal (9.5).

## 4. Estabilidad de los umbrales (A4)

Se remuestrean los 15 grupos con reposición (4000 réplicas, regla 1.5 × max): τ_b ∈ [0.0013, 0.0030] y τ_K ∈ [6.6, 11.45] (percentiles 5 y 95).

**Frecuencia con que cambiaría la clasificación:**

| Estructura | Papel | Clasificación oficial | Cambia en |
|---|---|---|---|
| **D6 tapa plana** | negativa | aceptada (error) | **0.3%** |
| D7 huso | negativa | aceptada (error) | 27% |
| D8 suma de normas | negativa | aceptada (error) | 28% |
| D5 gemelos | negativa | aceptada (error) | 36% |
| FRW polvo | positiva | aceptada | 36% |
| C0 Kasner | positiva | aceptada | 36% |
| Conos p = 1.5, p = 3 y C7 | zona gris | | 19% a 31% |
| Resto | | | menos de 1% |

**Lectura (nivel 1):**
- **El veredicto FALLA es robusto:** D6 queda aceptada en el 99.7% de los umbrales remuestreados.
- **Los márgenes son frágiles en ambas direcciones:** con otra muestra de calibración, D5, D7 y D8 podrían haberse rechazado, y FRW polvo y Kasner C0 podrían haber fallado. **Con 15 grupos, la regla del máximo no fija una tasa de error estable**, como advirtió el dictamen.

## 5. Sensibilidad a la regla (A5; descriptivo, no rescata nada)

| Regla | Errores |
|---|---|
| 1.0 × max | FRW polvo y C0 Kasner (positivas rechazadas); D6 y D7 (negativas aceptadas) |
| 1.5 × max (oficial) | D5, D6, D7 y D8 |
| 2.0 × max | D5, D6, D7 y D8 |

**Lectura (nivel 1):** **ninguna regla de esta familia clasifica todo bien.** El intercepto de D6 (+0.0005 sobre b_plano) es menor que el de positivas como FRW polvo (−0.0022 en módulo). Las distribuciones de b de geometrías y de algunos conos se solapan. El problema no es el umbral, es el observable: **b y K no contienen la información necesaria** para esta batería. Esto confirma la conclusión del dictamen con otro argumento.

## 6. Forma del cono frente a b (A6)

| Cono | ε | b − b_plano |
|---|---|---|
| C7 anidado (zona gris) | 0.159 | +0.0028 |
| p = 3 (zona gris) | 0.167 | +0.0020 |
| p = 1.5 (zona gris) | 0.201 | +0.0016 |
| D8 suma de normas | 0.234 | +0.0019 |
| **D6 tapa plana** | 0.250 | **+0.0005** |
| C4 esfera ∩ cubo | 0.266 | +0.0054 |
| C3 cilíndrico | 0.293 | +0.0116 |
| p = 6 | 0.307 | +0.0074 |
| **D7 huso** | 0.310 | **+0.0017** |
| C6 p por eje | 0.314 | +0.0104 |
| cúbico | 0.423 | +0.0179 |
| octaédrico | 0.732 | +0.0174 |

**Lectura:**
- **Nivel 1:** ε no ordena b. Conos con ε casi igual difieren en b por un factor de 10 (D6 frente a C4; D7 frente a p = 6).
- **Nivel 2:** los conos que el instrumento no ve son **redondos casi en todas partes**, con el defecto concentrado en una región pequeña de direcciones: una tapa en D6, dos puntas en D7, pliegues delgados en D8. r promedia sobre todas las direcciones, así que un defecto localizado se diluye. Una medida adecuada sería **el promedio** de la desviación sobre la esfera, no el máximo. Hipótesis para R14-C, no un resultado.

## 7. Trazabilidad

- Todas las partes del intento 2 tienen el commit 1619162 (verificado por `final` y aquí). Los resultados están en bff0241.
- La rama predeterminada no contenía el cierre: se integró por merge (eb757ae), conservando los SHA, y se abre un PR nuevo (R14-4).
- **Unidad independiente para el futuro (recomendación):** la realización dentro de cada estructura, y además **la familia generadora** como unidad para afirmaciones generales. F6 tuvo 27 estructuras con voto que caben en solo 8 familias amplias (geometrías conformemente planas, Kasner, transformaciones lineales de Minkowski, conos de norma suave, conos poliédricos, órdenes producto, percolación y defectos discretos). Cualquier tasa de error debe darse por familia, no por estructura.

## 8. Qué entrega R14-A a R14-B y R14-C

1. Seguir usando la realización como unidad y el bootstrap de realizaciones completas.
2. Reemplazar «1.5 × max de 15 grupos» por bandas de equivalencia preregistradas con más grupos, o con un percentil con intervalo.
3. Usar observables que promedien la forma del cono de otra manera, o que la resuelvan por dirección: la distribución completa de tamaños de intervalos y las cadenas (R14-C).
4. D5 a D8 quedan como **ejemplos vistos**, útiles para explorar y no para confirmar.
