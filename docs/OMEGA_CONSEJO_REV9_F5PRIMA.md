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
