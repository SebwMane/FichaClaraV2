# Consejo Ω, Revisión 6: evaluación de la firma espectral (2026-10-12)

- **Pedido del usuario:** antes de decidir R5-4, que el Consejo evalúe la firma espectral y explique qué es lo que distingue a la geometría.
- **Insumo nuevo:** exploración declarada del Consejo (no evidencia), semilla 777, N = 800, 3 semillas por estructura.

## 1. Qué se mide

Sobre cualquier orden se construye la matriz iΔ = (i/2)(C − Cᵀ), que codifica, para cada par, si uno está en el pasado o en el futuro del otro. Sus autovalores positivos, ordenados de mayor a menor (λ₁ ≥ λ₂ ≥ …), son los **modos** que un campo puede tener sobre ese orden: el primero es el más «largo» y los siguientes son cada vez más finos, como los armónicos de una cuerda. La firma registra dos cosas: cómo decaen (la pendiente β de ln λ_k frente a ln k) y dónde ese decaimiento se quiebra (la rodilla).

## 2. Exploración nueva

| Estructura | β mediana | Rodilla | Observación |
|---|---|---|---|
| Geometría 2D | −1.02 | 0.39 a 0.50 | cientos de modos en escalera regular |
| Geometría 3D | −0.73 | 0.92 a 1.12 | decae más lento |
| Geometría 4D | −0.64 | 0.80 a 1.09 | decae aún más lento |
| 2-orden en 2 capas | sin β | sin rodilla | casi sin modos |
| 2-orden en 16 capas | sin β | sin rodilla | solo 8 modos positivos |
| 2-orden en 100 capas | −1.17 | 2.28 | pendiente parecida a 2D, rodilla 5 veces más alta |
| KR | −0.04 | 6.5 a 6.7 | un modo enorme y luego una meseta plana |
| Bicapa | −0.04 | 6.6 a 6.7 | igual que KR |
| Sándwich | sin β | sin rodilla | casi sin modos |

## 3. Evaluación por roles

**Físico teórico.** La pendiente está midiendo la **dimensión**. En un espacio continuo de dimensión d, el número de modos con autovalor mayor que λ crece como una potencia de 1/λ fijada por d (ley de Weyl), y eso da λ_k ∝ k^(−2/d): β = −1 en 2D, −0.67 en 3D y −0.5 en 4D. Los valores medidos (−1.02, −0.73, −0.64) siguen esa tendencia, con desvíos en 3D y 4D atribuibles al tamaño finito (nivel 2 para la interpretación; la fórmula exacta en el orden discreto no está demostrada). En otras palabras: la firma es una **dimensión espectral** del sustrato, emparentada con la que las triangulaciones dinámicas causales miden con caminatas aleatorias.

**Estadístico.** Lo que separa a los rivales no es solo β. KR y la bicapa tienen un espectro de **matriz aleatoria**: un modo dominante y una meseta, porque sus relaciones se sortean sin estructura. Los 2-órdenes en pocas capas casi no tienen modos. El caso de 100 capas es el más cercano a la geometría: β = −1.17 sería compatible con la ventana C5-A, y lo separa la rodilla (2.28 frente a ≈ 0.45). La rodilla carga una parte importante de la discriminación.

**Metodólogo.** La firma no ha tenido ninguna prueba ciega. Los límites de C5-A se fijaron tras ver todas las familias de E1′, y esta exploración usa 3 semillas con N = 800. Hasta una prueba preregistrada con estructuras nuevas, su nivel es 2.

**Escéptico.** Dos objeciones. (1) La ventana C5-A está calibrada en 2D: una geometría 3D o 4D **no la pasaría** (β = −0.73 y −0.64). Tal como está, el criterio reconoce «geometría bidimensional», no geometría en general. (2) El 2-orden en 100 capas muestra que una estructura no geométrica puede acercarse a la pendiente correcta; si con más capas la rodilla también se acerca, la firma podría ser imitada. Esto hay que probarlo, no suponerlo.

## 4. Qué distingue a la geometría (síntesis)

1. **Muchos modos ordenados en escalera:** cientos de modos que decaen como una potencia regular, igual que los armónicos de un espacio continuo. Los rivales tienen pocos modos o una meseta sin estructura.
2. **La pendiente codifica la dimensión** (−2/d aproximadamente): es una medida de cuántas direcciones tiene el espacio que el campo «ve».
3. **El quiebre ocurre en modos de escala fina** (rodilla baja): la estructura regular se mantiene hasta casi la escala de los elementos individuales.

## 5. Conexión con la ruta de Jacobson (nivel 4, hipótesis)

La ley de área de la entropía depende de cómo se cuentan los modos cerca de una frontera, y la ley de Weyl es justamente ese conteo. Si la firma mide que el sustrato tiene el conteo de modos de un espacio continuo de dimensión d, sería el ingrediente previo que la ley de área necesita. No está demostrado; es la pregunta que la firma abre.

## 6. Recomendaciones del Consejo

1. **F1, prueba ciega preregistrada:** estructuras nunca medidas espectralmente (2-órdenes con 200 a 400 capas, órdenes aleatorios de dimensión fija, sprinklings en espacio-tiempo curvo, percolación con otros p), con ventanas fijadas antes.
2. **F2, dependencia de la dimensión:** comprobar con varios N si β tiende a −2/d en 2D, 3D y 4D. Si se confirma, reemplazar la ventana fija por una que dependa de la dimensión estimada por Myrheim-Meyer.
3. **Sobre R5-4:** la evaluación no cambia lo que mostró E1″ sobre la entropía; el Consejo sigue recomendando rebajar E1′ a «dependiente del criterio». La firma no entra en esa decisión porque no ha sido probada a ciegas.
