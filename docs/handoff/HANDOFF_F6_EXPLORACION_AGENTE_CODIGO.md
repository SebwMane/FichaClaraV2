# Handoff F6 exploratoria para el agente de código (2026-10-14)

**Especificación única:** `docs/OMEGA_F6_EXPLORACION.md` (congelada en e6de7be). Implementa solo eso; lo no especificado va al informe.

## Entregar

1. **`tools/f6_exploracion.py`**, con partes `--part`:
   - `ref` (referencia exacta por m: 4000 diamantes planos 4D de m puntos más las puntas; r_plano(m) y perfil plano N_k(m), k = 0..⌊m/2⌋);
   - `nula_<N>` (30 diamantes planos 4D completos con N elementos, medidos con el instrumento);
   - `est_<nombre>_<N>` (las realizaciones de una estructura en un N);
   - `final` (ajustes, bootstrap, tablas y `summary.json`).
   - Compuerta: el documento debe contener `R11-3 | **Firmada` y el árbol git debe estar limpio fuera de `results/`.
   - Salida en `results/omega_f6_exploracion/`.
   - Constantes: EXPERIMENT_ID 320, entropía 20261075. spawn_key `(320, familia, índice, N, realización)`, con familia 0 = referencias, 1 = nula, 2 = estructuras, 3 = dilución.
2. **Generadores:** reutiliza los de `tools/f7.py` y `tools/f5prima.py`: diamante, de Sitter (`de_sitter_eta`), FRW (`frw_power_points`, potencia d para radiación y 2d para polvo), conos (`polyhedral_cone`, `lp_cone`), `random_dorder(4)`, `transitive_percolation`. Agrega:
   - anti de Sitter de Poincaré: diamante desplazado con z = x₃ + 0.7 ∈ [0.2, 1.2], aceptación (0.2/z)⁴, relación de Minkowski. Puedes mirar `docs/consejo_rev11/pilotos/piloto_curv.py`.
   - bloque: `omega/causality/torus.py` con duración 0.45.
   - caja [0, 0.45] × [0, 1]³ con borde y relación de Minkowski. Informa la fracción de intervalos recortados: el diamante geométrico del intervalo se sale de la caja (se calcula con coordenadas solo para el informe).
3. **Rendimiento:** con N = 16000 la memoria importa (hay 15 GB). Construye C sin matrices float64 completas innecesarias (por bloques o float32) y calcula B una sola vez. Mide y reporta los tiempos.
4. **Instrumento por intervalo:** sub-orden de los m interiores más las puntas, tomado de C. r se calcula sobre los m interiores sin puntas; el perfil N_k(k ≤ ⌊m/2⌋) con puntas. m_max = máximo de B. Regla m ≤ 0.1·m_max; mínimo 30 intervalos; hasta 300 al azar.
5. **Ajustes en `final`:**
   - Δr = b + A·x con x = √(m/N), por estructura, con bootstrap de 1000 sobre realizaciones;
   - B·x² informado aparte;
   - el colapso de la dilución;
   - D_loc(m) en sus dos versiones (forma y máximo), frente a la nula por m.
6. **Tests rápidos** en un archivo nuevo `tests/test_f6_exploracion.py`:
   - en un diamante plano chico, r medio cerca de r_plano(m);
   - el sub-orden con puntas tiene m + 2 elementos y las puntas son mínimo y máximo;
   - el generador AdS da un orden parcial dentro del diamante;
   - la regla de localidad excluye m > 0.1·m_max.
7. **Verificación:** `python -m pytest -q tests/` en los archivos tocados y `python -m mypy` sin errores nuevos (hay 6 previos conocidos).

## No hacer

No ejecutar partes oficiales ni escribir en `results/omega_f6_exploracion/`, no usar la entropía 20261075 fuera de la constante, no tocar otros archivos (en especial `omega/causality/kasner.py` y `tests/test_kasner.py`, que hace otro agente), no hacer commit.
