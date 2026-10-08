# Guía local: instalar, probar y simular (Proyecto Ω)

## 1. Instalación

Requisitos: Python ≥ 3.11 y git.

```bash
git clone https://github.com/SebwMane/FichaClaraV2.git
cd FichaClaraV2
git checkout claude/omega-c0-competencia     # o la rama que quieras reproducir (ver §6)
python3 -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[dev]"                      # numpy, scipy, networkx + pytest, mypy, matplotlib
```

**Importante: un solo hilo BLAS.** Las herramientas fijan `OMP_NUM_THREADS=1` (enmienda C0-A4); con varios hilos, la dinámica era unas 30 veces más lenta. Si escribes tus propios scripts, exporta lo mismo:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
```

## 2. Pruebas

```bash
pytest                                   # suite rápida (las marcadas 'slow' se excluyen)
pytest tests/test_c0_*.py -v             # solo Ω-C0: funcional, localidad, dinámica
pytest tests/test_architecture*.py       # reglas de arquitectura (p. ej. RNG solo vía PCG64)
pytest -m slow                           # corridas completas (lentas)
mypy --strict --explicit-package-bases omega/c0
```

### Añadir una prueba propia

Crea `tests/test_<tema>.py`. Ejemplo mínimo, que comprueba la cota C0-T1 (S ≥ −Nψ*²/16κ) sobre un estado aleatorio:

```python
import numpy as np
from omega.c0.functional import action_c0, params_from_targets

def test_cota_c0_t1() -> None:
    p = params_from_targets(c_star=2, k_star=8, a=1.0)
    rng = np.random.Generator(np.random.PCG64(0))   # RNG siempre explícito con PCG64
    n = 64
    w = rng.random((n, n)) * 0.2
    w = (w + w.T) / 2
    np.fill_diagonal(w, 0.0)
    assert action_c0(w, p) >= p.lower_bound(n) - 1e-9
```

Reglas que hacen cumplir los tests de arquitectura:
- nunca `np.random.seed` ni `np.random.rand`: usa `np.random.Generator(np.random.PCG64(...))`;
- no modifiques el código congelado de Ω-1.0/1.1 ni los umbrales del certificado.

## 3. Mapa del código

| Ruta | Qué hay |
|---|---|
| `omega/c0/functional.py` | Acción S_C0, gradiente, `params_from_targets(c*, k*, a)`, cota inferior, KKT |
| `omega/c0/dynamics.py` | `evolve_c0`: descenso de gradiente proyectado en [0,1] con paso adaptativo |
| `omega/c0/locality.py` | `classify_c0` (VACIO / DENSO_TRIVIAL / FRAGMENTADO / DISPERSO_LOCAL / DISPERSO_NO_LOCAL) |
| `omega/c0/references.py` | Grafos de referencia: T³, RGG3, cliques, ER, toro triangular, anillo, árbol… |
| `omega/geometry/`, `omega/certificate/` | Dimensión, homogeneidad, isotropía, certificado Ω-1.1 |
| `tools/c0_*.py` | Scripts de cada etapa de Ω-C0 |
| `docs/OMEGA_C0_PRERREGISTRO.md` | Especificación congelada (parámetros, criterios, enmiendas) |
| `docs/OMEGA_C0_RESULTADOS.md` | Resultados e interpretación |
| `results/` | Resultados oficiales versionados (JSON y JSONL) |
| `runs/` | Pesos W finales (`.npz`). **Ignorado por git**: un clon nuevo no los trae |

## 4. Simulaciones

### 4.1 Prueba rápida (recomendado para empezar)

Los modos `--smoke` usan pocas celdas y semillas, escriben en `results/*_smoke` y no tocan los resultados oficiales:

```bash
python tools/c0_landscape.py --smoke --allow-dirty            # L1–L3, ~minutos
python tools/c0_dynamics.py  --smoke --allow-dirty --procs 2  # L4, 2 celdas × 3 inicios, 2000 pasos
```

Para escribir en otra carpeta: `--out /ruta/salida` (y `--runs-dir` para los `.npz` en dinámica).

### 4.2 Experimento a mano

```python
import numpy as np
from omega.c0.functional import params_from_targets
from omega.c0.dynamics import evolve_c0
from omega.c0.locality import classify_c0

p = params_from_targets(c_star=2, k_star=8, a=1.0)          # celda 37
rng = np.random.Generator(np.random.PCG64(1))
n = 216
w0 = rng.random((n, n)) * 2 * 8 / (n - 1)                   # inicio "R"
w0 = np.triu(w0, 1)
w0 = w0 + w0.T
ev = evolve_c0(w0, p, max_steps=5000)
print(ev["status"], ev["s_final"] / p.lower_bound(n))
print(classify_c0(ev["w"], np.random.Generator(np.random.PCG64(2))))
```

Una corrida de 20 000 pasos tarda ~1 min con N=216 y ~3 min con N=343.

### 4.3 Reproducir el pipeline oficial completo

**Orden obligatorio** (cada etapa lee la salida de la anterior):

| # | Comando | Lee | CPU aprox. |
|---|---|---|---|
| 1 | `python tools/c0_landscape.py` | — | minutos |
| 2 | `python tools/c0_dynamics.py` | — | ~5 h |
| 3 | `python tools/c0_redteam.py` | 2 (+ `.npz` de 2) | ~14 h |
| 4 | `python tools/c0_battery.py` | 2, 3 (+ `.npz`) | ~1 h |
| 5 | `python tools/c0_null_battery.py` | 4 (+ `.npz` de 3) | <1 h |
| 6 | `python tools/c0_f1.py` | 4, 5 | ~21 h |

Los tiempos son de CPU; divídelos por el número de procesos (`--procs`, por defecto 4).

Antes de lanzar una reproducción oficial, ten en cuenta dos protecciones:
- **Árbol limpio:** los scripts exigen un árbol git limpio (sin cambios sin commitear) y guardan el commit en cada resultado.
- **Resultados existentes:** si `results/<etapa>/runs.jsonl` ya existe, el script se niega a correr. Como esos archivos son los resultados oficiales, reprodúcelos en una rama o clon aparte (por ejemplo, `git checkout -b repro` y borra `results/c0_*`), nunca sobre la rama oficial.

Otras opciones:
- `--resume` retoma una corrida interrumpida;
- `--summarize-only` recalcula el resumen a partir del `runs.jsonl` existente, sin simular.

## 5. Leer resultados

| Archivo | Responde a |
|---|---|
| `results/c0_landscape/summary.json` | L1–L3: ¿las cliques o T³ son mínimos? ¿Hay KKT? |
| `results/c0_dynamics/summary.json` | L4: clases de los estados finales (675 corridas, N=216) |
| `results/c0_redteam/summary.json` | R6 (N=125, 343) y D-1 (100 000 pasos): celdas CANDIDATO-C0 |
| `results/c0_battery/summary.json`, `null_control.json` | D_eff, D_s, D_L y separación del nulo con grados fijos |
| `results/c0_f1/summary.json` | Veredicto F1 con N=512 y 729 y códigos del certificado |

```python
import json
rows = [json.loads(l) for l in open("results/c0_dynamics/runs.jsonl")]
gen = [r for r in rows if r["kind"] == "generic"]
local = [r for r in gen if r["cls"]["class"] == "DISPERSO_LOCAL"]
print(len(local), "/", len(gen))
print(local[0]["cls"])   # H_null, fracción de ciclos cortos, kmax/kmean...
```

Campos útiles de cada fila:
- `cell_idx`, `c_star`, `k_star`, `a`, `init` (U/E/R), `seed`;
- `status` (converged / max_steps / stalled);
- `S_over_LB` (acción / cota C0-T1);
- `kkt_residual`;
- `cls.class`;
- `cls.H_null`: salto medio relativo al grafo recableado con los mismos grados (mayor que 1 indica localidad).

## 6. Ramas

| Rama | Contenido |
|---|---|
| `claude/omega-1.1-congelado` | Ω-1.1 y bloque L, congelados |
| `claude/omega-c0-congelado` | Ω-C0 completo, congelado |
| `claude/omega-c0-competencia` | Rama de trabajo de Ω-C0 |
| `claude/omega-fase2-theta` | Fase 2: Θ>0 sobre C0 y preparación de C1 (`docs/OMEGA_FASE2_PRERREGISTRO.md`) |

Regla del proyecto: primero se prerregistra (documento con parámetros, semillas y criterios commiteado) y después se escribe el código y se corre.

## 7. Problemas frecuentes

| Problema | Solución |
|---|---|
| `ModuleNotFoundError: omega` | `pip install -e .` desde la raíz |
| "arbol sucio" | commitea, o usa `--smoke --allow-dirty` |
| "runs.jsonl ya existe" | `--resume`, o `--out` a otra carpeta |
| Falta un `.npz` en `runs/` | No están en git: corre la etapa anterior |
| Va lentísimo | Comprueba `OMP_NUM_THREADS=1` |
