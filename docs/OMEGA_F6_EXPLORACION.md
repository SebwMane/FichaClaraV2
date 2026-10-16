# F6, fase exploratoria: protocolo congelado antes de correr (2026-10-14)

- **Estado:** EXPLORATORIA (R11-3). **Sin voto.** Nada de lo que muestre cuenta como evidencia; sirve para diseñar la confirmación. Puede mirar positivas e impostores **ya vistos**. No incluye ninguna estructura reservada para la confirmación (Kasner, conos nuevos, links borrados, conjunto de reserva).
- **Orden firmado (R11-2):** rama 1, intercepto al achicar los diamantes; rama 2, perfil local; luego comparación.
- **EXPERIMENT_ID 320; entropía 20261075** (declarada exploratoria: nunca se usará en una confirmación).

## 1. Instrumento común

Para cada realización (una estructura con N elementos):

1. Matriz causal C y conteos «entre» B (`between_counts`).
2. **m_max:** el mayor número de elementos interiores de un intervalo de la estructura (medido del orden, sin coordenadas).
3. Para cada m de la rejilla con **m ≤ 0.1·m_max** (regla de localidad S1, fijada antes): hasta **300** intervalos [x, y] al azar entre los que tienen exactamente m interiores. Si hay menos de 30, esa (estructura, N, m) queda «no medible».
4. Por intervalo, sobre sus m interiores **sin las puntas**:
   - **r**: fracción de pares relacionados;
   - **perfil local** N_k para k = 0..⌊m/2⌋, contando los pares relacionados con k elementos entre ellos dentro del sub-orden de los m + 2 elementos (con puntas).
5. Por realización y m: media de r y media del perfil sobre sus intervalos (la unidad estadística es la realización, S4).

**Referencia exacta por m:** 4000 diamantes planos 4D de m puntos uniformes (más las puntas), con la misma medición. Da r_plano(m) y el perfil plano por m.

## 2. Rejilla y estructuras (fijadas antes)

- **m ∈ {10, 20, 40, 80, 160}** (cada uno solo si cumple la regla de localidad).
- **N ∈ {4000, 8000}** con 10 realizaciones por estructura; **N = 16000** con 5 realizaciones para el subconjunto marcado (★).
- **Dilución (control, sin papel de eje):** en N = 8000, además se mide la misma realización diluida al azar a q = 0.5 (4000 elementos), para ver si la desviación depende solo de m/(qN) (colapso).

| Estructura | Papel en la exploración | ★ |
|---|---|---|
| diamante 4D plano | geometría | ★ |
| de Sitter 4D (η ∈ [−1.5, −0.5]) | geometría curva R > 0 | |
| de Sitter 4D fuerte (η ∈ [−1.2, −0.2]) | geometría curva R > 0 | ★ |
| anti de Sitter 4D (Poincaré, z ∈ [0.2, 1.2], densidad ∝ z⁻⁴) | geometría curva R < 0 (ya vista por el asesor) | ★ |
| FRW radiación (a ∝ η) | geometría, R = 0 | |
| FRW polvo (a ∝ η²) | geometría curva R > 0 suave | |
| bloque T³ × [0, 0.45] | geometría sin borde | ★ |
| caja [0, 0.45] × [0, 1]³ | geometría con borde | |
| cono cúbico (p = ∞), cono octaédrico (p = 1) | impostores ya vistos, ε > 0.20 | ★ (octaédrico) |
| cono p = 6 | impostor ya visto, ε = 0.31 | |
| conos p = 1.5 y p = 3 | **zona gris** (ε ≤ 0.20) | |
| orden producto 4 | impostor ya visto | |
| percolación p = 0.05 | impostor ya visto | |

## 3. Rama 1: intercepto al achicar

- **Desviación:** Δr(m, N) = media de r − r_plano(m).
- **Eje:** x = √(m/N) (4D: tamaño lineal relativo del diamante chico).
- **Ajuste por estructura:** Δr = b + A·x, por mínimos cuadrados sobre todos los (m, N) medibles; errores por **remuestreo de realizaciones completas** (bootstrap de 1000). El término B·x² se informa aparte (fue sugerido por un dato ya visto: anti de Sitter).
- **Lectura esperada (H_loc, nivel 2):** geometrías con b ≈ 0 y signo de A según R (positivo en de Sitter, negativo en anti de Sitter, ≈ 0 en FRW radiación); impostores con b ≠ 0.
- **Colapso:** comparar Δr de la realización diluida (q = 0.5, N = 8000) con Δr a N = 4000 y el mismo m.

## 4. Rama 2: perfil local

- **v_k(m) = ln((N_k + 1)/(N_k^plano + 1))**, k = 0..⌊m/2⌋; **D_loc(m) = max_k |v_k − media_k(v)|** (forma, como F7) y también max_k |v_k| (como F5′).
- **Nula por m:** la misma medición sobre 30 diamantes planos (entropía exploratoria, otra rama del spawn_key) con N = 4000 y 8000.
- **Informe:** D_loc(m) por estructura y m, con la nula al lado; y el máximo sobre m (estilo max-T) como descriptor.

## 5. Qué se informa

Por estructura y N: m_max, intervalos disponibles por m, fracción recortada (caja), Δr(m) con error entre realizaciones, b y A con intervalos de bootstrap, D_loc(m) frente a la nula. Comparación de las dos ramas: qué separa cada una, dónde coinciden, dónde no, y qué sesgos de §3 de Rev. 11 aparecieron.

## 6. Lo que esta fase no puede hacer

No da veredictos. No toca Kasner, los impostores nuevos ni el conjunto de reserva. Cualquier regla que se proponga a partir de ella se congela después en un preregistro y se prueba con semillas y estructuras nuevas.
