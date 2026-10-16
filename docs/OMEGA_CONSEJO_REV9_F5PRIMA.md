# Consejo Ω, Revisión 9: preparación de F5′ (2026-10-12)

- **Estado:** PREPARADA, no ejecutada. Pedido del usuario: «que el Consejo prepare F5′».
- **Antecedente:** `OMEGA_F5_RESULTADOS.md`. El brazo C (perfil de intervalos) rechazó a los seis impostores, incluidos los tres órdenes producto, pero falló en el bloque 4D (dimensión asignada con Myrheim-Meyer fuera de un intervalo) y en de Sitter 4D (umbral justo, razón 1.21).

## 1. Obstáculo detectado al preparar

La corrección «medir dentro del mayor intervalo causal» no sirve para el bloque T³ × [0, 0.4]: el mayor diamante que cabe tiene altura 0.4 y volumen π·0.4⁴/24 ≈ 0.0034, es decir ≈ 0.8% del bloque (≈ 20 elementos con N = 2400). Ni el perfil de intervalos ni Myrheim-Meyer están validados fuera de regiones con forma de intervalo.

**Propuesta del Consejo:** F5′ se limita a regiones con forma de intervalo. La robustez frente a bordes y topología queda como pregunta separada (F6), que necesitaría un perfil medido solo en el interior.

## 2. Correcciones respecto de F5

1. **Dimensión asignada d* por vía espectral**, no por Myrheim-Meyer: d* = d_spec redondeada, con la curva de calibración de β_L. En Rev. 7 la pendiente con links fue casi insensible a la forma de la región (bloque −0.274 frente a diamante −0.312).
2. **Umbral del perfil con más muestras y con curvatura:** δ(d, N) = 1.5 × percentil 95 de D sobre 40 diamantes planos **más** 20 diamantes en un espacio-tiempo curvo de calibración (FRW con a(η) ∝ η², coordenadas conformes). de Sitter **no** entra en la calibración: se prueba a ciegas.
3. **Brazo A y brazo B** se informan sin voto (A falló como se predijo; B separa poco).

## 3. Diseño

**Calibración (entropías nuevas):** curva de β_L con diamantes planos d ∈ {2, 3, 4, 5}, N ∈ {1600, 2400}, 10 semillas; referencias de perfil con 20 diamantes por (d, N); umbral del perfil como en §2.2 (en 3D y 4D, las dimensiones de las positivas).

**Prueba (entropía nueva; N ∈ {1600, 2400}, 5 semillas por N):**

| Estructura | Papel | ¿Vista antes? |
|---|---|---|
| Diamante 4D plano | positiva | sí (semillas nuevas) |
| Diamante 3D plano | positiva | sí (semillas nuevas) |
| Diamante en de Sitter 4D | positiva | sí en F5; **no en la calibración** |
| Diamante en de Sitter 3D | positiva | **no** |
| Cono cúbico 4D: x ≺ y si Δt ≥ max(|Δx|, |Δy|, |Δz|) | negativa | **no** |
| Cono octaédrico 4D: x ≺ y si Δt ≥ |Δx| + |Δy| + |Δz| | negativa | **no** |
| Órdenes producto de dimensión 3, 4 y 5 | negativas | sí (en F5) |
| KR, percolación p = 0.05, capas gruesas | negativas | sí |

Los dos conos poliédricos son impostores nuevos: tienen eje temporal como un espacio-tiempo, pero sus conos tienen caras planas en vez de ser redondos. Son la prueba más directa de si el perfil detecta la **forma** del cono.

**Clasificación:** «espacio-tiempo de dimensión d*» si pasa la firma para d* (4D: ventanas de F3; 3D: regla de F3 sobre la calibración) **y** el perfil cumple D ≤ δ(d*, N), con mediana por estructura de la razón D/δ.

**PASA** si las cuatro positivas quedan clasificadas como espacio-tiempo de su dimensión y ninguna negativa. Costo estimado: 30 a 40 minutos.

## 4. Decisiones pendientes

| # | Pregunta |
|---|---|
| R9-1 | ¿Ejecutar F5′? |
| R9-2 | ¿Aceptar el alcance acotado a regiones con forma de intervalo (bloque fuera, F6 aparte)? |

## 5. Escala de afirmación

Un PASA sería nivel 1 para regiones con forma de intervalo, N ≤ 2400, 3D y 4D y estas estructuras: el primer criterio del programa que reconoce espacio-tiempo lorentziano en 4D, incluido uno curvo, frente a impostores con otros conos. No diría nada sobre regiones con otros bordes o topologías.

## 6. Firmas y congelamiento (2026-10-14)

| # | Decisión |
|---|---|
| R9-1 | **Firmada: ejecutar F5′**, en el orden referencias → calibración → prueba → final. |
| R9-2 | **Firmada: alcance acotado a regiones con forma de intervalo.** El usuario pidió primero incluir el bloque; al explicarle que el perfil no tiene un intervalo útil dentro del bloque, eligió «F5′ ahora, F6 después». **F6 (bordes y topología, con un perfil medido en el interior) queda como siguiente paso comprometido.** |

### 6.1 Cuenta previa sobre datos ya publicados de F5 (declarada; sin corridas nuevas)

Se aplicó la regla de d* espectral de §2.1 a `results/omega_f5/` (calibración y prueba de F5, ya publicadas). No es evidencia de F5′; sirvió para fijar detalles operativos.

- **Entre 4D y 5D la curva casi no separa en N = 1600** (medianas −0.313 y −0.316). Por corrida individual, d_spec es inestable cerca de 4 (un diamante 4D dio 6.6; varios diamantes 5D dan 4). **Por eso d* se toma de la mediana por estructura**, no de cada corrida: con F5, diamante 4D 3.98, de Sitter 4D 3.97, diamante 3D 2.99.
- **Bloque 4D:** d_spec ≈ 3.44, es decir d* = 3. La premisa de §2.1 («pendiente casi insensible a la forma de la región») no se sostiene para el bloque. Fuera del alcance de F5′; es un dato para F6.
- **Órdenes producto 4 y 5:** d_spec ≈ 3.59 y 3.55, cerca del corte 3.5. Cualquiera de los dos redondeos los compara con una referencia lorentziana; como son negativas, el veredicto depende solo de que el perfil los rechace.
- **Riesgo declarado:** una positiva 4D con β_L mediana más baja que −0.316 en N = 1600 recibiría d* = 5 y fallaría por el instrumento de dimensión, no por el perfil. Si ocurre, se registra así, sin cambiar la regla.

### 6.2 Detalles operativos congelados

- **EXPERIMENT_ID = 318.** Entropías: referencias 20261069, calibración 20261070, prueba 20261071. `spawn_key = (318, familia, d o índice de estructura, N, semilla)`.
- **N ∈ {1600, 2400}** en todo.
- **Referencias de perfil:** 20 diamantes planos por (d, N), d ∈ {2, 3, 4, 5}; referencia = media de N_k, k = 0..8 (`K_MAX`).
- **Calibración plana:** 40 diamantes por (d, N), d ∈ {2, 3, 4, 5}. Para todas: abundancias N_k. Para las semillas 0 a 9: además β_L, rodilla_L (links, `knee_from_spectrum`), d_MM y d_mid.
- **Calibración curva (FRW):** 20 diamantes por (d, N), d ∈ {3, 4}, en coordenadas conformes con a(η) ∝ η²: η = t + 1 ∈ [0.5, 1.5] (diamante de `sprinkle_diamond` desplazado), densidad ∝ a^d = η^{2d}, por rechazo con probabilidad (η/1.5)^{2d}. Relación causal = la de Minkowski en (η, x) (invariancia conforme). Solo abundancias.
- **Curva β_L(d) por N:** mediana de β_L de las 10 corridas planas con espectro. **d_spec** por interpolación lineal entre dimensiones vecinas, extrapolación con el tramo extremo (misma función `d_spec` de F5).
- **Dimensión asignada:** d* = mediana de d_spec por corrida sobre las 10 corridas de la estructura (5 por N), redondeada y recortada a {2, 3, 4, 5}.
- **Firma:** 4D, ventanas de F3 sin cambios (β ∈ [−0.363, −0.263], rodilla ∈ [1.511, 6.046]); 2D, 3D y 5D, ventanas con la regla de F3 (`windows`) sobre las 10 corridas planas con espectro de esa dimensión (ambos N juntos).
- **Umbral del perfil:** D = `profile_distance` a la referencia (d, N). δ(d, N) = 1.5 × percentil 95 de D sobre las 40 planas de (d, N) más, si d ∈ {3, 4}, las 20 FRW de (d, N). Para d ∈ {2, 5}, solo las 40 planas.
- **Clasificación:** «espacio-tiempo de dimensión d*» si β_L mediana y rodilla_L mediana caen en las ventanas de d* **y** la mediana sobre las 10 corridas de D/δ(d*, N) es ≤ 1.
- **Veredicto (vota solo esta clasificación):** **PASA** si las cuatro positivas quedan clasificadas como espacio-tiempo con d* igual a su dimensión verdadera y ninguna de las ocho negativas queda clasificada como espacio-tiempo. **FALLA** en otro caso.
- **Sin voto, solo informados:** brazo A (|d_MM − d_spec| ≤ τ_A, regla de F5) y brazo B (|d_MM − d_mid| ≤ τ_B, regla de F5), y el margen de cada estructura (razón D/δ y posición en las ventanas).

### 6.3 Estructuras de prueba congeladas (5 semillas por N)

| # | Nombre | Papel | d verdadera | Construcción |
|---|---|---|---|---|
| 0 | diamante_4d | positiva | 4 | `sprinkle_diamond(n, 4)` |
| 1 | diamante_3d | positiva | 3 | `sprinkle_diamond(n, 3)` |
| 2 | de_sitter_4d | positiva | 4 | η = t − 1 ∈ [−1.5, −0.5], aceptación (0.5/|η|)^4 (igual que F3 y F5) |
| 3 | de_sitter_3d | positiva | 3 | ídem con exponente 3 |
| 4 | cono_cubico_4d | negativa | | uniforme en {|t| ≤ T, max(|x|,|y|,|z|) ≤ T − |t|}; x ≺ y si Δt > 0 y Δt ≥ max(|Δx|, |Δy|, |Δz|) |
| 5 | cono_octaedrico_4d | negativa | | uniforme en {|t| ≤ T, |x|+|y|+|z| ≤ T − |t|}; x ≺ y si Δt > 0 y Δt ≥ |Δx| + |Δy| + |Δz| |
| 6 | producto_3 | negativa | | `random_dorder(3)` |
| 7 | producto_4 | negativa | | `random_dorder(4)` |
| 8 | producto_5 | negativa | | `random_dorder(5)` |
| 9 | kr | negativa | | `kr_order(n)` |
| 10 | percolacion_005 | negativa | | `transitive_percolation(n, 0.05)` |
| 11 | capas_gruesas | negativa | | `thick_layer_order(n, 8, 0.5)` |

T = 0.5 como en `sprinkle_diamond`. Los dos conos poliédricos se muestrean en su propio intervalo (la región entre el punto más bajo y el más alto según su orden), para que el perfil se compare en igualdad de condiciones.

### 6.4 Escala y lectura previa

Un PASA es nivel 1 para intervalos, 3D y 4D, N ≤ 2400 y estas doce estructuras. Un FALLA se informa por estructura, indicando si falló la dimensión asignada, la firma o el perfil.
