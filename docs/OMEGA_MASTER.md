# Proyecto Ω — Documento maestro de bifurcación e investigación

> **Estado (2026-10-07):** el programa de selección dimensional (problema B) está **cerrado** con un resultado negativo estructurado; el problema A sigue abierto, sin candidato. Punto de entrada: `docs/OMEGA_ESTADO.md`.

(Transcripción fiel del documento entregado por el usuario; la notación LaTeX se pasó a texto plano. Las secciones conservan su numeración original 0–51.)

## 0. Propósito
Concentra lo esencial de la investigación conceptual y matemática sobre una hipótesis de estructura pregeométrica Ω. Sirve como punto de bifurcación para nuevas conversaciones, experimentos, modelos y simulaciones. NO es una teoría física establecida: es un programa especulativo que busca convertir una intuición filosófica en una hipótesis matemática, computacional y eventualmente falsable.

**Regla metodológica principal:** No queremos demostrar que Ω tiene razón. Queremos descubrir qué puede producir Ω y dónde el universo podría demostrar que estamos equivocados.

## 1. Intuición inicial
El mapa observable del universo no debe confundirse con la estructura fundamental. Cadena de observación: universo → fenómeno físico → interacción → señal → instrumento → datos → modelo → interpretación. La luz sería principalmente un mensajero de información. "No observable electromagnéticamente" ≠ "físicamente inexistente". Pregunta: ¿podría la geometría observable ser una descripción emergente de una estructura más profunda?

## 2. Hipótesis Ω
Ω es una estructura física pregeométrica de la que podría emerger la geometría del espacio-tiempo. NO es "un espacio dentro del espacio". Inicialmente NO se asigna: coordenadas x,y,z; distancia convencional; geometría euclidiana; métrica de Einstein; tiempo físico convencional.
Representación mínima: Ω = (V, E, Q): V grados de libertad fundamentales; E relaciones; Q estados físicos/cuánticos.
Geometría efectiva: g_eff_{μν} = F(Ω). Esquema: Ω → g_{μν} → {materia, luz, gravedad}.

## 3. Condición para que Ω sea teoría física
Una interpretación que reproduzca exactamente GR para siempre sería solo reinterpretación ontológica. Se necesita: (1) dinámica definida; (2) geometría efectiva; (3) recuperar GR en algún límite macroscópico; (4) idealmente desviaciones observables. Conceptualmente g_eff = g_GR + δg; pregunta final: δg ≠ 0 ¿dónde y cómo?

## 4. Microestados de Ω
Podrían existir muchos microestados Q_1..Q_N con aproximadamente la misma geometría macroscópica (analogía termodinámica: temperatura emerge del comportamiento colectivo). microestructura → comportamiento colectivo → geometría. Guía conceptual S ~ k_B ln(Ω_microstates), no definición: la entropía correcta debe derivarse del modelo. Pregunta: ¿relación entre microestados de Ω y geometría emergente?

## 5. Distancia como propiedad emergente
La distancia no es fundamental. Dadas relaciones R_ij, podría existir D(i,j) = F(R_ij, R_ik, R_kj, ...). relaciones → distancia → geometría. Evita introducir espacio antes de probar que la dinámica lo produce.

## 6. Primera representación concreta
Ω = (V, E, W): nodos, relaciones, pesos. Matriz W simétrica (W_ij = W_ji), 0 ≤ W_ij ≤ 1, W_ii = 0.
**Importante:** la red dibujada NO representa espacio físico. Los nodos no tienen coordenadas. Los índices 1..N no tienen significado espacial.

## 7. Condición inicial fundamental
W_ij ~ U(0,1) con simetrización. NO usar: red cúbica, 2D, lattice hexagonal, coordenadas, vecinos espaciales, estructura diseñada a mano, geometría 3D preexistente. Razón: si aparece fase 3D, no debe poder atribuirse a la condición inicial. Toda dimensionalidad macroscópica debe ser propiedad de la dinámica.

## 8. Primera noción de distancia
Longitud de arista l_ij = 1/(W_ij + ε). Relación fuerte = separación pequeña. Distancia global d(i,j) = min sobre caminos γ:i→j de Σ_{(a,b)∈γ} l_ab. Ω → W → d.
Precaución: no considerar automáticamente todos los pares conectados con arista enorme. Usar umbral: A_ij = 1 si W_ij > W_min, 0 si W_ij ≤ W_min; calcular distancias sobre las relaciones existentes.

## 9. Dimensionalidad emergente
Bola B_i(r) = {j : d(i,j) ≤ r}, N_i(r) = |B_i(r)|. En espacio D-dimensional N(r) ~ r^D. D_eff(r) = d ln N(r) / d ln r. Ejemplos: N~r → D≈1; r² → 2; r³ → 3. NO asumir D=3: descubrir qué valor emerge.

## 10. Dimensionalidad dependiente de escala
Posible D_micro ≠ D_macro (p.ej. 1→3, 2→4, o ninguna transición). Medir D_eff(r) como función de la escala, no un único número.

## 11. Dimensión espectral
Difusión sobre la red; probabilidad de retorno P(σ). D_s(σ) = −2 d ln P(σ) / d ln σ. Puede revelar lo que la dimensión basada en distancias no detecta. Una candidata sólida debería mostrar D_eff ≈ 3 y D_s ≈ 3 en rangos de escala compatibles. D_eff ≈ 3 por sí solo NO demuestra geometría tridimensional.

## 12. Pregunta central sobre fases
¿Puede Ω poseer fases colectivas distintas? Ω_A → g^(A), Ω_B → g^(B). La misma dinámica podría generar organizaciones distintas según parámetros. Pregunta del primer experimento: ¿qué universos geométricos puede generar Ω? No se intenta aún demostrar que uno sea nuestro universo.

## 13. Primera funcional dinámica
S0[W] = −α T + β S_dens + γ S_deg
- 13.1 Triángulos: T = Σ_{i<j<k} W_ij W_jk W_ki. El término −αT favorece relaciones cerradas.
- 13.2 Densidad: S_dens = Σ_{i<j} W_ij². Con +β penaliza pesos grandes.
- 13.3 Regularidad: k_i = Σ_j W_ij; k̄ = (1/N) Σ_i k_i; S_deg = Σ_i (k_i − k̄)². Favorece homogeneidad de conectividad.

## 14. Por qué NO fijar k0 = 6
Alternativa S_deg = γ Σ (k_i − k0)² con k0=6 introduciría preferencia por estructuras 3D antes de comprobar que 3D emerge. k0 NO se fija; se usa k̄ como cantidad emergente. Preguntas: ¿qué conectividad produce espontáneamente la dinámica? ¿hay relación con D_eff ≈ 3?

## 15. No introducir inicialmente suavidad
Término propuesto: C_ij = Σ_k (W_ik − W_jk)², S_smooth = η Σ_ij W_ij C_ij. Podría favorecer geometría suave, pero **η = 0** en la primera versión. Primero ver qué produce −αT+βS_dens+γS_deg sin ayudar artificialmente. Luego introducir η y comparar.

## 16. Dinámica
dW_ij/dτ = −∂S0/∂W_ij, τ parámetro abstracto (no tiempo físico). Implementación inicial: W_ij ← clip[W_ij − Δτ ∂S0/∂W_ij, 0, 1]. Posteriormente parametrización W_ij = σ(θ_ij) para mantener 0 < W_ij < 1.

## 17. Gradiente
Para arista (a,b):
- ∂T/∂W_ab = Σ_{k≠a,b} W_ak W_kb
- ∂S_dens/∂W_ab = 2 W_ab
- ∂S_deg/∂W_ab = 2[(k_a − k̄) + (k_b − k̄)]
- ∂S0/∂W_ab = −α Σ_{k≠a,b} W_ak W_kb + 2β W_ab + 2γ[(k_a − k̄) + (k_b − k̄)]
Base matemática del primer simulador.

## 18. Primer mapa de fases
Redundancia de escala: fijar β = 1. Parámetros relevantes α/β, γ/β. Mapa (α/β, γ/β) → fase. Con η: espacio 3D (α/β, γ/β, η/β).

## 19. Posibles fases (provisionales; no asumir que existen)
- A Dispersa: k̄→0 con componentes pequeñas.
- B Fragmentada: varias componentes grandes.
- C Conectada y homogénea: G/N≈1, baja dispersión relativa de grados.
- D Altamente agrupada: clustering elevado.
- E Hiperdensa: conectividad excesiva.
- F Candidata geométrica: componente gigante; estructura local homogénea; escalas con crecimiento de bolas; D_eff estable; luego D_s compatible; robustez ante semillas y parámetros; reproducible.

## 20. Problema de la fase trivial
+βS_dens y +γS_deg pueden favorecer W_ij → 0 (fase trivial). No ocultarlo: indica tendencia hacia solución dispersa. Segunda rama experimental: término de control S_mass = −μ Σ_{i<j} W_ij, o restricción Σ_{i<j} W_ij = W0. No introducirlo desde el principio.

## 21. Qué medir
O = (k̄, σ_k, C, L, G, D_eff, D_s, ...): conectividad media, dispersión de conectividad, clustering, longitud característica de caminos, tamaño de componente gigante, dimensión efectiva, dimensión espectral. No clasificar fases con una sola variable.

## 22. Experimento 0 — Validación
Construir redes conocidas 1D, 2D, 3D y recuperar D_eff ≈ 1, 2, 3. Si el algoritmo no lo logra, los resultados de Ω no son confiables.

## 23. Experimento 1 — Baseline aleatorio
W_ij ~ U(0,1) sin dinámica. Medir conectividad, clustering, componente gigante, distancias, D_eff, D_s.

## 24. Experimento 2 — Dinámica
Aplicar S0 desde pesos aleatorios. Registrar toda la evolución, no solo el estado final.

## 25. Experimento 3 — Barrido de parámetros
Malla α/β = 0, 0.1, 0.2, ... y γ/β = 0, 0.1, 0.2, ... (resolución adaptable al coste). Cada punto con múltiples semillas. Resultado: mapa parámetros → estructura/fase.

## 26. Repetición con semillas
Una semilla no basta. Para cada punto interesante seed_1..seed_n. ¿Emerge la misma fase independientemente de la condición aleatoria? Crítico si aparece D≈3.

## 27. Tamaño del sistema
Etapa 1: N = 100–300. Luego N = 500–1000. Si hay región interesante: N → 10³–10⁴+ con estructuras dispersas. Matriz densa: almacenamiento O(N²); no comenzar con N = 100 000.

## 28. Computador personal
Sí. Para N~100–300 matriz densa es razonable. Dificultad posterior: muchos nodos, semillas, barridos enormes, simulaciones largas, coarse-graining, dinámica cuántica, búsqueda estadística exhaustiva. No necesarias para el primer experimento.

## 29. Robustez de una fase
Probar: parámetros (variaciones pequeñas), semillas (muchas), tamaño (N=100,200,300,500,...), algoritmo (distintos métodos de estimar dimensión), evolución (cambiar Δτ y criterios de convergencia), condición inicial (siempre aleatoria, distintas distribuciones y semillas).

## 30. ¿La fase 3D es un atractor?
Si muchas condiciones aleatorias distintas terminan en la misma clase geométrica, la interpretación es mucho más interesante que una trayectoria única.

## 31. Histéresis y transiciones de fase
Variar un parámetro gradualmente λ1<λ2<... y luego en sentido contrario. Si ida ≠ vuelta puede haber histéresis, evidencia de transición de fase.

## 32. Invariancia ante renombramiento
Los índices no tienen significado físico. Si W → P W P⁻¹ (permutación), los observables deben permanecer invariantes: D_eff(W) = D_eff(P W P⁻¹). Detecta geometría accidental por orden de índices.

## 33. Ablación
Comparar: −αT; −αT+βS_dens; −αT+γS_deg; −αT+βS_dens+γS_deg. ¿Qué término es responsable de cada propiedad emergente? Evita ajuste retrospectivo.

## 34. Coarse-graining
Si aparece estructura interesante, estudiarla a distintas escalas: Ω_micro → Ω_meso → Ω_macro; D_micro → D_macro. Resultado especialmente interesante: microestructura complicada que tras coarse-graining se comporte como geometría suave de D≈3.

## 35. Reconstrucción de coordenadas
Solo después de descubrir estructura geométrica: d_ij → {x_i} (multidimensional scaling u otros). Las coordenadas son resultado de la geometría, no entrada de la simulación.

## 36. Métrica emergente
Con coordenadas reconstruidas: d_ij² ≈ g_{μν} Δx^μ Δx^ν; ajustar métrica local g_eff. Secuencia: Ω → W → d → D → x → g_{μν}. Objetivo central.

## 37. Curvatura emergente
δW → δd → δg → R_{μνρσ}. ¿Puede la curvatura surgir de cambios en la organización microscópica de Ω?

## 38. Materia y geometría
Hipótesis futura T_{μν} → Ω' → g'_{μν}. NO imponer G_{μν} = (8πG/c⁴) T_{μν} desde el inicio; buscar un mecanismo microscópico cuya aproximación macroscópica lo reproduzca.

## 39. Causalidad
La red inicial no contiene causalidad. Posible ampliación Ω = (V, E_rel, E_causal, Q). Meta: estructura relacional + orden causal → espacio-tiempo efectivo. Posterior a la emergencia de geometría espacial efectiva.

## 40. Fases del proyecto completo (niveles)
1 Dinámica relacional; 2 Fases topológicas; 3 Dimensionalidad emergente; 4 Geometría efectiva; 5 Métrica; 6 Curvatura; 7 Causalidad; 8 Límite relativista; 9 Predicciones nuevas.

## 41. Jerarquía de hipótesis
H0 no emerge geometría estable. H1 emergen estructuras geométricas estables. H2 aparece D_eff≈3. H3 la fase 3D es robusta ante condiciones iniciales. H4 geometría efectiva aproximadamente continua. H5 deformaciones microscópicas producen curvatura. H6 puede emerger causalidad. H7 el límite macroscópico reproduce GR. H8 existe δg ≠ 0 en algún régimen (convertiría el modelo en teoría diferenciable experimentalmente).

## 42. Posibles lugares de falsación futura
Curvatura extrema (agujeros negros); ondas gravitacionales (dispersión, amortiguamiento, polarizaciones, ringdown); propagación cosmológica (desviaciones acumulativas); vacíos (T_{μν}≈0 pero Ω≠∅; NO implica automáticamente materia/energía oscura).

## 43. Preguntas abiertas
Sobre Ω: qué es un grado de libertad; relaciones discretas o continuas; clásicas o cuánticas; qué determina los pesos; por qué existe la estructura. Geometría: ¿emerge la distancia? ¿la dimensión? ¿por qué 3? ¿dimensión no entera? ¿dependiente de escala? Fases: ¿estables? ¿transiciones? ¿histéresis? ¿atractores? ¿universalidad? Física: ¿métrica, curvatura, causalidad, estructura lorentziana, GR? Observación: ¿qué desviación produciría Ω, en qué régimen, qué la falsaría, qué observación sería incompatible?

## 44. Principio epistemológico
No estamos buscando demostrar que tenemos razón, sino el punto exacto donde el universo podría demostrar que estamos equivocados. Explícitamente: no imponer D=3; no imponer k0=6; no comenzar con red 3D; no introducir coordenadas espaciales; no ajustar el modelo tras ver el resultado; usar controles; repetir semillas; hacer ablaciones; comprobar invariancias; registrar resultados negativos; buscar fases inesperadas; intentar destruir el modelo.

## 45. Resultado débil
D_eff≈3 en una sola simulación.

## 46. Resultado fuerte (dentro del programa computacional)
1 Parte de pesos aleatorios. 2 Sin coordenadas iniciales. 3 Sin estructura 3D inicial. 4 Aparece componente gigante. 5 D_eff≈3. 6 D_s≈3. 7 Ambos estables en rango de escalas. 8 Con múltiples semillas. 9 Para distintos N. 10 Persiste ante pequeñas variaciones de parámetros. 11 No depende de un algoritmo particular. 12 Sobrevive a ablaciones. 13 Comportamiento de transición de fase. 14 Compatible con coarse-graining. 15 Admite geometría efectiva suave. 16 Eventualmente permite reconstruir una métrica.
No demostraría que el universo funciona así; sería evidencia computacional de que una geometría 3D puede emerger de una dinámica que no la introdujo.

## 47. Qué NO demostraría
Que Ω exista físicamente; que el universo esté hecho de redes; que la gravedad o GR sean emergentes; que se explique el origen del universo; que nuestro universo corresponda a cierta fase; materia oscura; energía oscura.

## 48. Arquitectura propuesta del primer simulador
```text
OMEGA/
├── network/        initialization.py, weights.py, topology.py
├── dynamics/       functional.py, gradient.py, evolution.py
├── geometry/       distances.py, dimension.py, spectral.py
├── phases/         scan.py, classification.py, stability.py
├── coarse_graining/
├── curvature/
├── causality/
├── statistics/
└── experiments/    test_00_validation.py, test_01_random.py, test_02_dynamics.py,
                    test_03_phase_diagram.py, test_04_dimension.py, test_05_robustness.py
```

## 49. Pasaporte de reproducibilidad
Cada simulación guarda P = {N, α, β, γ, η, Ω0, τ_max, Δτ, D_eff, D_s, k̄, σ_k, C, G, curvatura, seed}.

## 50. Primer objetivo práctico
Construir un simulador de Ω que empiece exclusivamente con pesos aleatorios y descubra el catálogo de fases producido por S0. Experimento completo: pesos aleatorios → dinámica → estado final → métricas topológicas → D_eff → D_s → clasificación de fase. Después: mapa de fases. Solo después: geometría emergente.

## 51. Frase de estado
Queremos descubrir qué estructuras geométricas pueden emerger espontáneamente de una dinámica relacional sin coordenadas espaciales iniciales, y determinar si alguna puede conducir posteriormente a una teoría física falsable. La pregunta inmediata no es "¿Ω genera nuestro universo?" sino "¿qué universos geométricos puede generar Ω?". La primera herramienta: el mapa de fases producido por pesos iniciales completamente aleatorios.
