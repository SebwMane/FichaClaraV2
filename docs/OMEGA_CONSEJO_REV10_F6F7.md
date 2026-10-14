# Consejo Ω, Revisión 10: preparación de F6 y F7 (2026-10-14)

- **Estado:** PREPARADAS, no ejecutadas. Pedido del usuario: «que el Consejo prepare F6 y F7».
- **Antecedentes:** `OMEGA_F5PRIMA_RESULTADOS.md` (F5′ FALLA: el cono octaédrico pasa; el perfil no distingue conos redondos de poliédricos) y R9-2 (F6 comprometida: regiones con otros bordes y topologías).

## 1. Las dos preguntas

| Prueba | Pregunta | En simple |
|---|---|---|
| **F7** | ¿Un criterio sobre la **forma** del perfil de intervalos distingue un cono redondo de uno con caras planas, sin rechazar espacio-tiempos curvos? | F5′ medía «cuánto se aleja» el perfil. F7 mide «de qué manera se aleja». |
| **F6** | ¿El perfil, medido **dentro de intervalos pequeños** de la estructura, funciona en regiones que no son un intervalo (un bloque con espacio periódico o una caja con bordes)? | En vez de mirar la estructura entera, se miran muchos «diamantes chicos» dentro de ella. |

El Consejo recomienda **F7 primero**: F6 necesita un perfil que rechace a los conos poliédricos, y eso es justamente lo que F7 pone a prueba.

## 2. F7: distancia de forma

### 2.1 Idea y cuenta previa (declarada)

Con v_k = ln((N_k + 1)/(ref_k + 1)), k = 0..8, la distancia de F5′ es D = max_k |v_k|. La **distancia de forma** es

  D_forma = max_k |v_k − media_k(v)|,

es decir, se descuenta el desplazamiento parejo del perfil (que es lo que produce la curvatura) y se mide solo lo que cambia la forma.

**Cuenta previa sobre datos ya publicados de F5′ (no es evidencia; la idea salió de esos datos):**

| | D_forma típica |
|---|---|
| Diamantes planos y FRW de calibración 4D (p95) | 0.031 a 0.046 |
| Diamante 4D / de Sitter 4D | 0.027 / 0.041 |
| Conos cúbico / octaédrico | 0.178 / 0.166 |
| Órdenes producto 3, 4, 5 | 0.24 / 0.40 / 0.91 |

Con la regla de umbral de F5′ (1.5 × p95), δ_forma(4D) ≈ 0.06 a 0.07: los conos poliédricos quedarían unas 2.5 veces por encima y de Sitter en ≈ 0.6. **Riesgos que ve el escéptico:** (a) de Sitter 4D ya tiene D_forma mayor que el plano (0.041 frente a 0.027), así que con más curvatura podría cruzar el umbral; (b) un cono «casi redondo» puede ser indistinguible con N ≤ 2400. F7 incluye una prueba a ciegas de cada riesgo.

### 2.2 Diseño

Igual que F5′ (`tools/f5prima.py`, §6 de Rev. 9) salvo lo que se indica: misma firma espectral, misma d* espectral por mediana, mismo N ∈ {1600, 2400}, mismas referencias y calibración (40 planos + 20 FRW con a ∝ η²), **entropías nuevas**.

- **Criterio que vota:** firma para d* **y** mediana de D_forma/δ_forma(d*, N) ≤ 1, con δ_forma = 1.5 × p95 de D_forma en la calibración (misma regla que F5′).
- **Informado sin voto:** la clasificación con D de F5′ (réplica de F5′ con semillas nuevas).

**Estructuras de prueba (5 semillas por N):**

| Estructura | Papel | ¿Vista antes? |
|---|---|---|
| Diamante 4D, diamante 3D | positivas | sí |
| de Sitter 4D y 3D (η ∈ [−1.5, −0.5]) | positivas | sí (F5′) |
| **de Sitter 4D con más curvatura** (η ∈ [−1.2, −0.2]; el factor de escala cambia ×6 en vez de ×3) | positiva | **no** |
| **FRW de radiación 4D** (a ∝ η, η ∈ [0.5, 1.5]) | positiva | **no** (la calibración usa a ∝ η²) |
| Conos cúbico (p = ∞) y octaédrico (p = 1) | negativas | sí (F5′) |
| **Cono L^p con p = 1.5** | negativa | **no** |
| **Cono L^p con p = 3** | negativa (ver R10-2) | **no** |
| Cono L^p con p = 6 | solo informativo (curva de resolución) | **no** |
| Órdenes producto 3, 4, 5; KR; percolación 0.05; capas gruesas | negativas | sí |

Cono L^p: puntos uniformes en {|t| ≤ T, ‖x‖_p ≤ T − |t|}; x ≺ y si Δt > 0 y Δt ≥ ‖Δx‖_p. Es un orden parcial para p ≥ 1 por la desigualdad triangular. p = 2 es Minkowski. A medida que p se acerca a 2, el cono se vuelve más redondo y más difícil de detectar.

**PASA** si las seis positivas quedan clasificadas como espacio-tiempo de su dimensión y ninguna negativa con voto. Costo estimado: 12 a 15 minutos.

## 3. F6: perfil local

### 3.1 Obstáculo y salida

El perfil y Myrheim-Meyer solo valen en intervalos. El bloque T³ × [0, T] es «ancho y bajo»: su mayor intervalo es apenas 1% a 2% del total. Esto es geométrico: para que la imagen mínima del toro sea única hace falta duración < lado/2, así que ningún bloque sin borde deja un intervalo grande.

**Salida propuesta:** medir **muchos intervalos chicos**. Si un intervalo [x, y] contiene exactamente m elementos, en un sprinkling esos m puntos están repartidos uniformemente dentro del diamante entre x e y (propiedad de Poisson). Entonces el perfil promedio de esos sub-intervalos tiene una **referencia exacta**: m puntos uniformes en un diamante más sus dos puntas. Esto vale en cualquier región causalmente convexa (diamante, bloque periódico), sin importar el borde ni la topología.

**Piloto de factibilidad (declarado; semilla no oficial 777; sin medir discriminación):** con N = 4000 hay ≈ 970 intervalos con m = 20 en el bloque (T = 0.45), ≈ 590 en una caja con borde y más de 6000 en un diamante. Con N = 2400 hay pocos (84 en el bloque). Costo de `between_counts` con N = 4000: 0.4 s.

### 3.2 Diseño

- **N = 4000; m = 20** (vota) y m = 10 (informado). Por corrida se toman hasta 300 intervalos al azar entre los que tienen exactamente m elementos interiores.
- **Dimensión local d*_loc:** fracción de pares relacionados promediada sobre los sub-intervalos, comparada con la de referencia para d ∈ {2, 3, 4, 5} y el mismo m; se elige la d más cercana.
- **Perfil local:** media de N_k (k = 0..8) sobre los sub-intervalos (cada uno con sus dos puntas, m + 2 elementos), comparada con la referencia de d*_loc. Distancia: **la que F7 valide** (D_forma si F7 PASA; D de F5′ si F7 FALLA, y en ese caso se espera que los conos poliédricos pasen, ver R10-3).
- **Umbral:** 1.5 × p95 de la distancia local en la calibración: 40 diamantes planos y 20 FRW (a ∝ η²) por d con N = 4000.
- **Si una estructura tiene menos de 50 intervalos con m interiores** en la mediana de sus corridas, es «no medible» y **no** se clasifica como espacio-tiempo (ver R10-4).
- **Firma espectral:** sin voto en F6. La pendiente depende de la forma de la región (bloque −0.26 frente a diamante −0.31, F5), así que exigirla castigaría justamente lo que F6 quiere probar.

**Estructuras (5 semillas, N = 4000):**

| Estructura | Papel | Qué prueba |
|---|---|---|
| Diamante 4D | positiva | control: el perfil local funciona donde el global ya funcionaba |
| Bloque T³ × [0, 0.45] | positiva | sin borde espacial, topología de toro |
| Caja [0, 0.45] × [0, 1]³ con bordes | positiva | borde espacial; **la caja no es causalmente convexa** y los intervalos cerca del borde quedan recortados (riesgo declarado) |
| **Bloque FRW** (a ∝ η, η ∈ [0.5, 0.95], es decir duración 0.45; espacio periódico; densidad ∝ η⁴) | positiva | curvatura + toro, **a ciegas** |
| Orden producto 4 | negativa | |
| Bloque con cono cúbico y bloque con cono octaédrico (toro, imagen mínima) | negativas | impostores con otra topología |
| KR, percolación 0.05, capas gruesas | negativas | se esperan «no medibles» |

**PASA** si las cuatro positivas quedan clasificadas como espacio-tiempo 4D y ninguna negativa. Costo estimado: 8 a 12 minutos.

## 4. Lo que diría cada resultado

- **F7 PASA:** primer criterio que separa conos redondos de poliédricos en 4D, con curvatura a ciegas. Nivel 1 para intervalos, N ≤ 2400, y los p probados.
- **F7 FALLA en p = 3:** el perfil tiene un límite de resolución; se informa dónde (con p = 6 como referencia).
- **F7 FALLA en de Sitter fuerte o FRW de radiación:** la forma del perfil también cambia con la curvatura; el criterio confundiría curvatura con forma del cono.
- **F6 PASA:** el reconocimiento no depende del borde ni de la topología (en lo probado).
- **F6 FALLA en la caja pero no en el bloque:** el problema es el borde que recorta intervalos, no la topología.

## 5. Decisiones pendientes

| # | Pregunta |
|---|---|
| R10-1 | ¿Ejecutar F7? |
| R10-2 | ¿El cono p = 3 vota (prueba exigente, puede fallar por resolución) o solo se informa junto a p = 6? El escéptico pide que vote; el metodólogo advierte que un FALLA ahí no distinguiría «criterio malo» de «N chico». **Recomendación del Consejo: que vote**, y que el resultado se lea con p = 6 al lado. |
| R10-3 | ¿Ejecutar F6 después de F7, con la distancia que F7 valide (y si F7 falla, igual correr F6 con D para responder solo la pregunta de bordes y topología)? |
| R10-4 | En F6: ¿aceptar que la firma espectral no vote y que «no medible» cuente como rechazo? |

## 6. Convenciones previstas

F7: EXPERIMENT_ID 319, entropías 20261072 (referencias), 20261073 (calibración), 20261074 (prueba). F6: EXPERIMENT_ID 320, entropías 20261075 a 20261077. Se congelan al firmar.
