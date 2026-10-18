# R14-B: compatibilidad microscópica (H_L) frente a emergencia efectiva (H_E). Preregistro exploratorio (2026-10-18)

- **Mandato:** Rev. 14, R14-3 firmada (preparar R14-B como preregistro exploratorio, **sin voto**). **Estado:** PREPARADO, **no ejecutado**. Su ejecución requiere una orden del usuario.
- **Prioridad del programa (R14-1, firmada):** emergencia efectiva (H_E) como prioridad situacional; compatibilidad microscópica (H_L) como hipótesis separada.
- **Lo que esta fase no puede hacer:** no da veredictos, no usa conjuntos sellados y no reabre F6. Cualquier regla que proponga se congela después en un preregistro confirmatorio con controles ciegos nuevos.

## 1. Las dos hipótesis, con definiciones discretas

**Nulo de referencia S(M, g, ρ):** el sprinkling de Poisson en una región con forma de intervalo de un espacio-tiempo (M, g) con densidad ρ. Para cada observable O se construye su distribución bajo S con al menos 100 sprinklings de referencia del mismo N.

- **H_L (compatibilidad microscópica):** un orden C es H_L-compatible con (M, g) si, **a la escala de los elementos**, todos sus observables caen en las bandas de equivalencia de S(M, g, ρ). Las bandas son el intervalo central del 95% de la referencia, ensanchado por un margen δ_O fijado abajo.
- **H_E (emergencia efectiva):** C cumple H_E si existe un **mapa de grano grueso** Φ_q tal que Φ_q(C) es H_L-compatible con algún (M, g) para todo q ≤ q*, y los observables de Φ_q(C) **convergen** al hacerse el grano más grueso.
- **Mapa de grano grueso (fijado antes de mirar datos):** dilución aleatoria. Se conserva cada elemento con probabilidad q, de forma independiente, y se toma el orden inducido. Es intrínseco (no usa coordenadas) y lleva un sprinkling de densidad ρ a otro de densidad qρ.
- **Escalas:** q ∈ {1, 1/2, 1/4, 1/8}, partiendo de N₀ = 8000 (de 8000 a unos 1000 elementos).

**Predicciones que distinguen H_L de H_E (nivel 2, congeladas aquí):**

| Estructura (ya vista) | H_L | H_E | Predicción |
|---|---|---|---|
| Geometrías (planas, curvas, Kasner) | sí | sí | observables en banda para todo q |
| D5 gemelos 10% | no | **posible** | un par de gemelos sobrevive con probabilidad q², así que la fracción de gemelos cae como q·10% y el defecto debería desvanecerse al diluir |
| C5 links borrados 13% y D4 relaciones agregadas | no | por determinar | sus defectos están ligados a la escala de los elementos; la exploración mide si escalan con q |
| Conos D6 a D8 y conos de norma | no | no | la forma del cono no cambia al diluir (autosemejanza); la desviación debería mantenerse con q |

**Lo que diría cada resultado (sin voto):**
- D5 converge a banda: los gemelos serían un microestado compatible con H_E, no con H_L. Así se separan las dos hipótesis.
- D6 a D8 no convergen: impostores de H_E también.
- Alguna geometría falla: el mapa de grano grueso o las bandas son inadecuados.

## 2. Observables (todos intrínsecos), en cada (estructura, q)

1. **Fracción de pares relacionados** del orden completo y su dimensión de Myrheim-Meyer d_MM.
2. **r por sub-intervalo** con m ∈ {10, 20, 40, 80}, con la regla m ≤ 0.1·m_max y la referencia exacta r = 1/10 en 4D (instrumento de F6).
3. **Perfil de intervalos** del orden completo (N_k, k = 0..8) frente a la referencia S.
4. **Pendiente A de curvatura** (ajuste Δr = b + A·√(m/N)), como observable de curvatura.
5. **Fracción de pares de elementos con igual pasado y futuro** (detector directo de gemelos; para geometrías es casi 0).

## 3. Estructuras (todas ya vistas)

- **Geometrías:** diamante plano 4D, de Sitter, de Sitter fuerte, anti de Sitter, FRW polvo, C0 y D0 (Kasner), D1 (cizalla), D2 (FRW exponencial).
- **Defectos:** C5, D4, D5.
- **Conos:** D6, D7, D8, C4 y el octaédrico.
- **Producto 4** (control).
- **Realizaciones:** 10 por estructura en N₀ = 8000; cada realización se diluye a los cuatro q con generadores independientes.
- **Unidad estadística:** la realización (R14-A §1). Las afirmaciones generales se dan por familia (R14-A §7).

## 4. Bandas y márgenes (fijados antes)

- Referencia S: 100 sprinklings del diamante plano 4D por cada N efectivo (8000, 4000, 2000 y 1000); para las geometrías curvas, la referencia es la geometría plana más la curvatura medida en A (descriptivo).
- δ_O = 0.5 × la desviación estándar de O en la referencia (margen de equivalencia, fijado aquí sin ver datos).
- Informe: la trayectoria O(q) con su banda, y si entra en banda y en qué q.

## 5. Convenciones

- **EXPERIMENT_ID 323; entropía 20261088** (exploratoria: nunca se usará para confirmar).
- spawn_key (323, familia, índice, N₀, realización, nivel de q).
- **Costo estimado:** unas 2 a 3 horas (18 estructuras, 10 realizaciones y 4 niveles, con N ≤ 8000), más las referencias.
- **Herramienta:** a cargo de un agente de código, con handoff y auditoría, como siempre; árbol limpio; tandas de menos de 2 horas.
