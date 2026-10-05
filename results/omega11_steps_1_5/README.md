# Ω-1.1 — Resultados congelados de los pasos 1–5 (P§22)

- **Commit de ejecución:** `167efca` (árbol limpio; `code_commit` y `config_hash` en cada `summary.json`).
- **Script:** `tools/run_steps_1_5.py`, con `master_entropy=20240901`.
- **Contenido:** `summary.json` por paso (fase `preregistered` → `final`), pasaportes v1.1 en JSON (`passports_json.tar.gz`), registro de IDs Ω-EXP y log de tiempos.
- **Arreglos `.npz`:** no se versionan, por tamaño (~170 MB). Cada pasaporte guarda su sha256, y las corridas son deterministas por semilla, así que se pueden reconstruir con `omega.io.provenance.reconstruct`.

| Paso | Duración | Completo | Expectativas preregistradas |
|---|---|---|---|
| p01/p02 analítica + golden | 67 s | sí | cumplidas |
| O-00 validación | 77 s | sí | cumplidas |
| O-01 baseline S0 | 789 s | sí | cumplidas |
| O-02 falsación de distancia | 569 s | sí | cumplidas |
| O-03 modelos nulos | 520 s | sí | cumplidas |
| O-04 ablación + Ω-B | 3579 s | sí | **49/51 cumplidas** |

## Incidencias de ejecución
- **La compuerta de prerrequisitos no se relajó en ningún momento.** Todo se ejecutó sobre el mismo commit, sin cambiar parámetros.
- **Primer corte:** la tarea en segundo plano se detuvo por el límite de tiempo del entorno durante O-03. Se relanzó desde O-03.
- **Segundo corte:** un reinicio de la sesión mató O-04 al ~92%. Su salida parcial se archivó (no versionada) como `_aborted_o04_ablation_restart_2314` y O-04 se repitió completo.

## Predicción preregistrada fallida (O-04, Ω-B)
- **Celdas afectadas:** ρ ∈ {0.05, 0.1}, γ̂ = 10, factor = 0.5 (por debajo del umbral lineal α̂ρ(N−4)/(N−2) = 1+γ̂).
- **Esperado:** F1 (estado uniforme).
- **Obtenido:** F2 en 10/10 semillas.
- **Estado final:** CONVERGED y heterogéneo. Del 5% al 9% de las aristas supera w_min, la componente gigante abarca el 23–30% de los nodos y se conserva ⟨W⟩ = ρ.
- **Interpretación:** la predicción se basaba solo en la estabilidad lineal del estado uniforme. Desde U(0,1) proyectado, la dinámica con γ̂ alto alcanza otros puntos KKT en la frontera.
- **Tratamiento:** se registra como predicción fallida y no se modifica nada (M§44). No es geometría.
