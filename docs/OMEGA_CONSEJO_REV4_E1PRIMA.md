# Consejo Ω, Revisión 4: diseño de E1′ (ruta entrópica) para consideración futura

- **Fecha:** 2026-10-12.
- **Estado:** PROPUESTA DE DISEÑO, **no congelada ni ejecutada**. Pedido del usuario: «diseña E1′ para salir de la duda, pasa la información al Consejo para ser considerado más adelante».
- **Antecedentes que se deben leer:** `OMEGA_E_PROPUESTA.md`, `OMEGA_E_RESULTADOS.md`, `OMEGA_2_INFORME_FINAL.md`.

## 1. Lo que dejó E1 y por qué hace falta E1′

1. **E1 terminó INDETERMINADO** porque la regla de región («intervalo de tamaño más cercano a N/4») elige regiones descentradas y alargadas, para las que no vale el coeficiente ≈ 1/3 que E0 validó en subdiamantes centrados.
2. **KR casi cumplió el criterio**: entropía truncada logarítmica con coeficiente 0.436, a 0.006 del borde de la ventana. Esto sugiere que el umbral √n/(4π), calibrado en diamantes 2D, puede fabricar apariencias de ley de área en órdenes no geométricos.
3. **La duda que E1′ debe resolver:** con regiones bien formadas y un criterio que no dependa solo de una ventana del coeficiente, ¿la ley de área es exclusiva de los órdenes geométricos? Si no lo es, el caso 3 de Ω-0.2 pierde su principal ventaja (un criterio de geometría sin competencia).

## 2. Diseño propuesto de E1′

### 2.1 Regla de región intrínseca con control de forma

Entre los intervalos I[a, b] (con extremos) de tamaño dentro de N/4 ± 5%, elegir el que **maximiza min(|pasado(a)|, |futuro(b)|)**, donde |pasado(a)| es el número de elementos que preceden a a.

Justificación (nivel 1, cálculo directo): en un diamante 2D de coordenadas de cono de luz [−1/2, 1/2]², un rectángulo de lados Δu, Δv con área fija ΔuΔv = 1/4 tiene pasado de su extremo inferior (1/2 − Δu/2)(1/2 − Δv/2) si está centrado; con el área fija, ese valor es máximo para Δu = Δv, y descentrarlo reduce el pasado de un extremo o el futuro del otro. El óptimo es el subdiamante centrado y cuadrado de E0, con |pasado(a)| = |futuro(b)| ≈ N/16.

### 2.2 Criterio de «ley de área genuina» (tres condiciones, todas obligatorias)

1. **Forma logarítmica:** en el caso truncado, R² del ajuste S = a ln(√n/4π) + b mayor que el del lineal.
2. **Coeficiente con su incertidumbre:** |a − 1/3| ≤ 2·SE(a) + 0.05.
3. **Complementariedad:** con el doble truncamiento, la entropía de R debe igualar la de su complemento causal R^⊥ = {x : x no está relacionado con ningún elemento de R}: media de |S(R) − S(R^⊥)| / S(R) ≤ 0.15. El truncamiento de R^⊥ usa √|R^⊥|/(4π). Sorkin y Yazdi reportan que la igualdad se cumple en la geometría solo tras truncar; es una propiedad de estado puro que una coincidencia en el coeficiente no imita.

### 2.3 Compuerta de validación V′ (antes de mirar ninguna familia no geométrica)

En 10 sprinklings 2D por tamaño, N_L ∈ {800, 1200, 1600, 2000, 2400}:
- **Forma:** |pasado(a)|/N y |futuro(b)|/N en [0.045, 0.080] (centro 1/16 = 0.0625) en ≥ 90% de las corridas.
- **Reproducción:** las tres condiciones de 2.2 se cumplen en el control geométrico.
- Si V′ falla, E1′ no se ejecuta y el diseño vuelve al Consejo.

### 2.4 Familias y decisión

- Familias: control geométrico (sprinklings 2D), KR, percolación p = 0.05, capas gruesas (8 capas), bicapa, sándwich; 10 semillas por tamaño; misma malla que V′.
- Familias sin intervalos de tamaño N/4 ± 5%: «no aplicable: sin regiones extensas» (cuenta como sin ley de área).
- **DISCRIMINA:** el control geométrico cumple las tres condiciones y ninguna familia no geométrica las cumple.
- **NO DISCRIMINA:** alguna familia no geométrica cumple las tres.
- **INDETERMINADO:** V′ falla, o el control geométrico no cumple las tres en E1′.

### 2.5 Rama secundaria (descriptiva, sin voto)

Truncamiento por la «rodilla» del espectro de iΔ en lugar del umbral fijo, como sugieren Sorkin y Yazdi. La regla para ubicar la rodilla debe fijarse antes de correr; en órdenes sin rodilla no se trunca. Sirve para ver si el caso KR depende del umbral calibrado en 2D.

## 3. Riesgos conocidos

1. **Complementariedad en el borde de la validación:** si en la geometría la igualdad S(R) = S(R^⊥) no se cumple con la tolerancia propuesta, V′ falla. Antes de congelar conviene fijar la tolerancia con un piloto declarado.
2. **El umbral √n/(4π)** sigue siendo una convención 2D aplicada a todo; la rama secundaria es la forma de vigilarlo.
3. **Potencia:** E0b′ pasó con ventaja pequeña del ajuste logarítmico (R² 0.384 frente a 0.361). Por eso E1′ usa 10 semillas.
4. **Costo:** unas 300 diagonalizaciones con N hasta 2400, del orden de 45 a 60 minutos de cómputo en este entorno, en tramos.

## 4. Decisiones pendientes del usuario

| # | Pregunta |
|---|---|
| R4-1 | ¿Ejecutar E1′ (V′ primero) cuando se retome Ω-E? |
| R4-2 | ¿Aceptar la complementariedad como condición obligatoria del criterio? |
| R4-3 | ¿Incluir la rama secundaria de truncamiento por rodilla? |
| R4-4 | Si E1′ da NO DISCRIMINA, ¿cerrar el caso 3 o reformularlo (por ejemplo, ley de área con otra definición de entropía)? |

## 5. Escala de afirmación aplicable

Un DISCRIMINA sería nivel 1 para 2D, N ≤ 2400 y estas familias: mostraría que el criterio separa geometría de los rivales conocidos, no que la ruta de Jacobson funcione completa. Un NO DISCRIMINA sería nivel 1 sobre la insuficiencia de esta medida de entropía como criterio de geometría.
