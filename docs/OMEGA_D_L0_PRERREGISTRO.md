# Ω-D — Cercanía por difusión: análisis previo y puerta L-ΩD-0 (premisa) — Prerregistro

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-d`.
- **Base congelada:** `claude/omega-p1d3-congelado`.
- **Mandato:** el usuario y el Consejo eligen Ω-D como candidato principal, solo como **hipótesis matemática** por ahora. Se exige D1 (definición), D2 (paisaje), D3 (estabilidad), D4 (escalas) y D5 (prueba de retazos) antes de cualquier simulación dinámica.
- **Estado:** este documento se commitea antes de escribir el código y antes de evaluar ningún grafo.

## 0. Auditoría del mensaje del Consejo

| Afirmación | Contraste con los registros | Se conserva |
|---|---|---|
| «C0 tiende a vacío/cliques» | **Incorrecto.** C0 **rompió** esa dicotomía: en C0-L4 hubo 0 VACÍO y 253/675 DISPERSO-LOCAL. Su fallo fue otro: localidad **sin geometría** (dimensión ≈ 2–2.5, no variedad, no homogénea, frágil a Θ > 0). La tendencia a vacío o cliques era de S0 y Ω-B | La lección correcta: C0 consiguió **localidad**, y Ω-D debe conservarla (requisito R1) |
| «P1-D.3 mostró límites en bordes, gradientes y cruces» | **Parcial.** El borde sí es un falso negativo. El gradiente **pasó** los Niveles I y II (GRUESA(3)); lo inestable fue el certificado (1/3 claves). Los cruces (WS) engañan al Nivel II aislado, pero el instrumento I+II los rechaza (33/33) | Dominio declarado: espacios cerrados casi homogéneos |
| «D ≈ 3 no significa espacio 3D (WS)» | **Correcto.** WS β = 0.001 da D_conv = 2.93–3.28 convergente | Sí |
| «La propiedad que falta es mesoscópica» | **Correcto.** Conteos (k, c, q) de C0 ≈ RGG3; P1-D: los retazos no homogeneizan | Sí: es la premisa que aquí se pone a prueba |
| D5, dibujo de los retazos como bloques **desconectados** | **Hay que precisarlo.** Un retazo desconectado se distingue trivialmente por conectividad (F2). El adversario correcto es el de P1-D: bloques con la **misma** geometría interna y la **misma** secuencia de grados, pegados entre sí por aristas reemparejadas al azar | Se usa el retazo de P1-D |
| D4: «dos escalas explícitas» | **Riesgo.** Una escala fija τ junto con un valor objetivo puede **imponer la dimensión** (§2, L-ΩD-T2) | Se adopta como **covarianza de escala** (§3), no como escala fija |
| D5bis: partir de un retazo y ver si la dinámica lo **cura** | Coincide con el requisito R6 de la decisión anterior | Sí, para después de L-ΩD-1..3 (requiere dinámica) |
| Árbol: Ω-D primero, P3 si fracasa, nada a N grande al principio | Coherente con el preflight A–H | Sí |
| «Ω-1.1: las geometrías ni siquiera eran puntos de equilibrio» | Se está verificando contra los documentos del bloque L (auditoría de lectura en curso); se informará en los resultados | Lección D3: comprobar estabilidad antes de simular |

## 1. D1: definición de la cercanía por difusión

Para pesos W simétricos en [0, 1], con fuerzas k_i = Σ_j W_ij y K = diag(k):

- **Operador de difusión perezoso simétrico:** M = ½ (I + K^{−1/2} W K^{−1/2}). Es semejante al paseo perezoso P = ½(I + K⁻¹W), tiene espectro en [0, 1] y no hay oscilaciones de paridad.
- **Perfil de difusión** de i a tiempo τ: u_i(τ) = M^τ e_i.
- **Solapamiento de difusión** (la «cercanía»):

      α_τ(i, j) = ⟨u_i(τ), u_j(τ)⟩ / (‖u_i(τ)‖ ‖u_j(τ)‖) = (M^{2τ})_ij / sqrt((M^{2τ})_ii (M^{2τ})_jj)  ∈ [0, 1].

  Vale 1 si las nubes de difusión coinciden y 0 si no se tocan. Es adimensional y no usa coordenadas.
- **Distancia de difusión normalizada por escala:** Q_τ(i, j) = −τ · ln α_τ(i, j).

## 2. Análisis previo (antes de evaluar ningún grafo)

**L-ΩD-T1 (clique).** Toda recompensa monótona en α favorece las cliques. En una unión de cliques, las nubes de los miembros coinciden tras pocos pasos (α → 1, Q → 0): la clique es el máximo de cualquier ψ(α) creciente. Una Ω-D con «premiar aristas de alto solapamiento» repetiría S0 (vacío/cliques). Hace falta un ingrediente contra la degeneración.

**L-ΩD-T2 (imponer la dimensión con un valor objetivo).** En el retículo Z^d, el paseo perezoso tiene varianza por coordenada τ/(2d). Para vecinos a distancia 1, la aproximación gaussiana da

    α_τ(vecinos) ≈ exp(−d/(2τ)),   es decir,   Q_τ(vecinos) ≈ d/2.

Una funcional que premie un valor objetivo α* (como el c* de C0) a escala fija τ selecciona d ≈ 2τ·ln(1/α*): **la dimensión sale de los parámetros.** Violaría el requisito R4 y el criterio de meseta C1-D. Queda prohibido el diseño «α objetivo a escala fija».

**L-ΩD-T3 (generalización de C0-T1, nivel 2).** Si la acción solo depende de α_τ sobre las aristas más un presupuesto de fuerza, los estados de igualdad son grafos «α-regulares», con todas las aristas en el mismo α. Ahí entran toros y también grafos de Cayley no geométricos. La degeneración persiste mientras la acción mire una sola escala.

## 3. La salida que propongo: covarianza de escala (diffusive scaling)

Por L-ΩD-T2, en una geometría (cualquier d) Q_τ(i, j) es **independiente de τ** mientras τ ≪ tiempo de mezcla: depende de la distancia, no de la escala de observación. En cambio:

| Estructura | Comportamiento de Q_τ al crecer τ | Interpretación |
|---|---|---|
| Expansor (ER, RR) | Cae a 0: las nubes convergen exponencialmente a la estacionaria | No difusivo |
| Clique | ≈ 0 desde el principio | Degenerado: puntos que coinciden |
| Retazo | Las aristas entre bloques empiezan con Q grande y luego Q cae, cuando la mezcla entre bloques es rápida | Incoherente |

Por tanto, una cercanía **coherente** es la que cumple la **ley de escala difusiva**, s_e = d ln Q_τ(e) / d ln τ ≈ 0, y además no es degenerada (Q > 0).

Este criterio no fija d: cualquier dimensión cumple s ≈ 0, y el valor de Q no se impone. Es una propiedad **multiescala** (D4) y **mesoscópica** (premisa central). Antes de construir ninguna acción sobre él, hay que comprobar que la premisa es cierta en grafos de referencia: eso es L-ΩD-0.

## 4. Puerta L-ΩD-0: ¿contiene la difusión la información que falta? (medición sobre grafos fijos, sin dinámica)

### 4.1 Panel (N ≈ 1000; clave PCG64 (20261010, familia, semilla); semillas 0–2 en las familias aleatorias)

| Clase | Grafos |
|---|---|
| G | RGG3 k12; RGG2 k12; T³ 10³; toro cuadrado 32² (N = 1024); anillo k12 |
| A (adversarios con la misma estadística local) | Retazos-2³ y retazos-3³ (RGG3 k12, de P1-D); RGG3 con 1 % de atajos; WS anillo k12 β = 0.01 |
| R (degenerados / regulares) | Unión de cliques K10 (100 cliques); caveman K8; retículo de cliques 3D K8, L = 5 (N = 1000) |
| NL | ER k12, RR k12 |
| Informativo | Los 24 finales C0-R de N = 729 (`runs/c0_f1/out/`), con su soporte fuerte |

### 4.2 Cálculo

- M densa sobre la componente gigante.
- M^{2τ} por cuadrados sucesivos para τ ∈ {1, 2, 4, 8, 16, 32}.
- α_τ y Q_τ sobre **todas** las aristas, y sobre una muestra de 5000 pares a distancia 2 como contraste.
- Por arista, s_e = pendiente de ln Q_τ(e) frente a ln τ en τ ∈ {2, 4, 8}, por mínimos cuadrados. Si α = 1 numéricamente (Q = 0), la arista se marca como **degenerada**.

**Medidas por grafo:**
- mediana de s_e;
- f_nd = fracción de aristas no difusivas (s_e < −0.5 o degeneradas);
- mediana de Q_τ(e) para cada τ;
- mediana y percentil 10 de α_τ(e).

### 4.3 Criterios (puerta de muerte de Ω-D)

| Criterio | Pregunta | Condición |
|---|---|---|
| **K1 (primario)** | ¿La difusión distingue el pegado incoherente? | f_nd(retazo) ≥ f_nd(RGG3) + 0.05, en las 3 semillas, para retazos-2³ **y** retazos-3³ (se compara cada semilla del retazo con la mediana de f_nd de las RGG3). Los atajos al 1 % y WS se informan, no deciden |
| **K2** | ¿Hay ley de escala difusiva independiente de la dimensión? | En G, \|mediana de s_e\| ≤ 0.15. En NL, mediana de s_e ≤ −0.5 |
| **K3** | ¿Degeneración de las cliques? | En la unión de cliques y en caveman K8, la mediana de Q_8(e) ≤ 0.1 · mediana de Q_8(e) de RGG3 |

**Decisión:**

| Resultado | Consecuencia |
|---|---|
| **K1 falla** | **Ω-D muere** en su premisa: la difusión no contiene la información del pegado incoherente. Se pasa a P3 |
| **K1 pasa y K2 falla** | La difusión distingue, pero no existe una formulación libre de dimensión con N ≈ 1000. Solo quedarían diseños «α objetivo», prohibidos por T2. Ω-D queda **en suspenso** y el Consejo decide |
| **K1 y K2 pasan** | Se define la acción candidata en L-ΩD-1 (covarianza de escala + no degeneración + presupuesto tipo C0) y se pasa a paisaje (D2), estabilidad (D3) y prueba de retazos energética (D5). Sigue sin dinámica |
| K3 | Es informativo para el diseño: confirma, o no, que hace falta el término de no degeneración |

### 4.4 Predicciones del cerebro (congeladas)

1. **K1 pasa.** Las aristas entre bloques de los retazos son no difusivas: su Q cae al cruzar τ el tiempo de mezcla entre bloques.
2. **K2 pasa en G** con N ≈ 1000 y τ ≤ 8. Q_τ en el retículo cumple la aproximación ≈ d/2 con un error ≤ 30 %: T³ ≈ 1.5; cuadrado ≈ 1.0. NL cae fuertemente.
3. **K3 pasa:** las cliques son degeneradas (Q ≈ 0).
4. **Retículo de cliques 3D** (geometría gruesa): aristas internas degeneradas y aristas entre sitios difusivas. Mezcla; se informa.
5. **Finales C0:** f_nd intermedio, entre el RGG3 y los retazos.

## 5. Lo que no se hace

- Ni dinámica ni acción.
- La acción solo se define después de K1/K2, y en otro prerregistro (L-ΩD-1).
- No se toca el certificado ni la batería P1-D.3.
- **Separación declarada:** la batería mide bolas en saltos y Ω-D usa el núcleo de difusión. Pero D_s de la batería también deriva de la difusión, así que en las fases siguientes el juez geométrico será la batería de bolas (Niveles I y II) más el certificado, **no** D_s. Así se evita que el juez y el mecanismo coincidan.
