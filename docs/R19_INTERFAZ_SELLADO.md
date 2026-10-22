# R19-(b): interfaz del archivo sellado (para el guardián y el ejecutor)

El archivo sellado `sellado_r19.py` es un módulo de Python que define:

```python
FAMILIAS: list[dict]   # cada una: {"nombre": str, "tipo": "se_desvanece" | "autosemejante", "estructuras": list[str]}
ESTRUCTURAS: dict[str, dict]
```

Cada entrada de `ESTRUCTURAS` es de uno de dos tipos:

1. `{"tipo": "copias", "k": int, "frac": float}`
   - **Defecto de copias** sobre una base de sprinkling 4D del diamante plano (`sprinkle_diamond(n_base, 4, rng)` y relación de Minkowski).
   - Se agregan copias de elementos originales distintos, con pasado y futuro idénticos al original, al estilo de `d5_build` en `tools/f6_intento2.py`.
   - Cada original elegido recibe **k** copias. k = 1 son gemelos; k = 2, tripletes.
   - **frac** es la fracción de copias sobre el total, entre 0.05 y 0.25.
   - Con un total N: n_copias = round(frac·N) y n_base = N − n_copias.
   - El **ejecutor** implementa el generador y el control reinyectado: diluye la base y después agrega copias con la misma k y la misma frac sobre el N diluido.
2. `{"tipo": "region", "generar": callable}`
   - `generar(n: int, rng: np.random.Generator) -> np.ndarray` (bool, n×n): un sprinkling de Poisson de **n** puntos en una región fija con una relación causal fija (por ejemplo, un cono de gauge F como `f6_intento2._gauge_gen`). La densidad escala con n y la región no cambia.
   - Debe ser distinta de D6, D7, D8 y de los conos de F6 y F7: cilindro, esfera y cubo, p por eje, anidado, L^p con p = 1.5, 3 y 6, cúbico y octaédrico.
   - El ejecutor usa `generar(N0, …)` diluido como observado y `generar(Binomial(N0, q), …)` como control.

**Reglas para el guardián:**
- **Familias:** al menos 2 de cada tipo, con 2 o 3 estructuras cada una.
- **Qué no leer:** `results/omega_r14b`, `results/omega_r14c`, `docs/OMEGA_R14*_RESULTADOS.md` ni `docs/consejo_rev19/`.
- **Verificación con sus propias semillas** (ninguna entropía 2026111x): para cada estructura de copias, el control con N = 1000 queda fuera de la banda plana. Para eso usa el código de observables y bandas de `tools/r14b.py` (solo el código), con referencias planas propias.
