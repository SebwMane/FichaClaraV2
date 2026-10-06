# Guía de Instalación y Desarrollo Local — Proyecto Omega

## 1. Requisitos previos

- **Python 3.11+**
- **pip** o equivalente
- **Git**

## 2. Instalación

### Clonar el repositorio

```bash
git clone https://github.com/SebwMane/FichaClaraV2.git
cd FichaClaraV2
```

### Crear un entorno virtual (recomendado)

```bash
python3 -m venv venv
source venv/bin/activate  # En Windows: venv\Scripts\activate
```

### Instalar dependencias

```bash
pip install -e .                    # Instalar el paquete en modo desarrollo
pip install -e ".[dev]"            # Instalar también dependencias de desarrollo
```

Esto instala:
- **Núcleo**: numpy ≥1.26, scipy ≥1.11, networkx ≥3.0
- **Desarrollo**: pytest ≥7.4, matplotlib ≥3.7, mypy ≥1.8

## 3. Estructura del proyecto

```
FichaClaraV2/
├── docs/
│   ├── OMEGA_C0_PRERREGISTRO.md      # Especificación preregistrada de C0
│   └── OMEGA_C0_RESULTADOS.md        # Resultados finales
├── omega/
│   ├── c0/                           # Funcional y dinámica de Ω-C0
│   │   ├── functional.py             # Acción, gradiente, parámetros
│   │   ├── dynamics.py               # Evolución por descenso de gradiente
│   │   ├── locality.py               # Clasificación y medidas locales
│   │   └── references.py             # Grafos de referencia
│   ├── config/
│   │   ├── seeds.py                  # Generación de claves RNG (SeedKey)
│   │   └── settings*.py              # Configuración de medidas geométricas
│   ├── geometry/                     # Dimensión, distancia, estructura local
│   └── experiments/v11/
│       └── gate.py                   # Herramientas del certificado Ω-1.1
├── tools/
│   ├── c0_landscape.py               # Prueba L1–L3 (puerta analítica)
│   ├── c0_dynamics.py                # Prueba L4 (dinámica a N=216)
│   ├── c0_battery.py                 # Descriptores geométricos (C0-A6)
│   ├── c0_null_battery.py            # Control nulo (C0-A7)
│   ├── c0_redteam.py                 # R6 (escalado N=125,343), D-1 (persistencia)
│   ├── c0_f1.py                      # F1 y F1b (N=512,729, certificado)
│   └── l3b_stability.py              # Certificado Ω-1.1
├── tests/
│   ├── test_c0_*.py                  # Tests unitarios (18 pruebas)
│   └── ...                           # Tests de arquitectura
├── results/                          # Salida de simulaciones (datos JSON, JSONLines)
├── runs/                             # Pesos finales (archivos .npz comprimidos)
├── DEVELOPMENT.md                    # Este archivo
└── pyproject.toml                    # Configuración de proyecto
```

## 4. Pruebas unitarias

Ejecutar todos los tests:

```bash
pytest
```

Incluye 39 tests en total (18 en `tests/test_c0_*.py` + tests de arquitectura).

Tipos de test:
- **Unitarios**: funcionalidad de módulos (`test_c0_functional.py`, `test_c0_locality.py`, etc.)
- **Integración**: fin a fin (`test_c0_integration.py`)
- **Arquitectura**: garantías del sistema (`omega/experiments/v11/`)

Ver detalles con:

```bash
pytest -v
pytest -k "c0_functional"       # Solo tests de funcional
pytest -k "locality"            # Solo tests de localidad
pytest --co                     # Listar tests sin ejecutar
```

Verificar tipos con mypy (estricto):

```bash
mypy --strict --explicit-package-bases omega/c0
```

## 5. Simulaciones de Ω-C0

Las simulaciones se corren en etapas, cada una prerregistrada en `docs/OMEGA_C0_PRERREGISTRO.md`. Deben ejecutarse **en orden** y la salida de cada etapa alimenta la siguiente.

### 5.1 Landscape (C0-L1, L2, L3)

**Qué hace**: Analiza la energía y el equilibrio KKT en 75 celdas (parrilla de parámetros).

```bash
python tools/c0_landscape.py
```

- **Salida**: `results/c0_landscape/{runs.jsonl, summary.json, log.txt}`
- **Datos**: energía, gradiente, KKT, codegrado en referencias (T³, RGG3, toro triangular)
- **Tiempo**: ~5–15 min

### 5.2 Dinámica (C0-L4)

**Qué hace**: Corre descenso de gradiente desde inicios genéricos (U, E, R) en N=216 para todas las celdas y clasifica los estados finales.

```bash
python tools/c0_dynamics.py [--resume] [--summarize-only] [--procs 4]
```

- **Salida**: `results/c0_dynamics/{runs.jsonl, summary.json, log.txt}`; pesos en `runs/c0/out/`
- **Flags**:
  - `--resume`: continuar desde donde paró
  - `--summarize-only`: solo regenerar resumen sin correr
  - `--procs 4`: número de procesos paralelos (ajustar a tu CPU)
- **Tiempo**: ~30 min (CPU-bound, paralelizable)

### 5.3 Batería descriptiva (C0-A6) + Control nulo (C0-A7)

**Qué hace**: Computa dimensión efectiva, dimensión espectral, homogeneidad, isotropía, salto promedio en estados DISPERSO-LOCAL.

```bash
python tools/c0_battery.py
```

- **Requisito previo**: `c0_dynamics.py` debe haber terminado
- **Salida**: `results/c0_battery/{summary.json, finals.jsonl}`
- **Tiempo**: ~20 min

Tabla de diagnóstico (D-1 + control nulo):

```bash
python tools/c0_null_battery.py
```

- **Salida**: `results/c0_battery/null_control.json`
- **Tiempo**: ~10 min

### 5.4 Red-team (R6, D-1)

**Qué hace**:
- **R6**: escalado a N=125 y 343 en celdas candidatas
- **D-1**: continúa corridas no convergidas a 100 000 pasos para verificar persistencia

```bash
python tools/c0_redteam.py [--resume] [--summarize-only] [--procs 4]
```

- **Requisito previo**: `c0_dynamics.py` completo
- **Salida**: `results/c0_redteam/{runs.jsonl, summary.json, log.txt}`; pesos en `runs/c0_redteam/out/`
- **Tiempo**: ~1–2 horas

### 5.5 Prueba F1 (Certificado + escalado a N=512 y 729)

**Qué hace**: Corre la subfamilia R-3D en N=512 y N=729, aplica el certificado Ω-1.1 y compara contra referencias (T³, RGG3).

```bash
python tools/c0_f1.py [--resume] [--summarize-only] [--procs 4]
```

- **Requisito previo**: `c0_redteam.py` completo
- **Salida**: `results/c0_f1/{runs.jsonl, summary.json, log.txt}`; pesos en `runs/c0_f1/out/`
- **Veredicto final**: F1-POSITIVO / F1-INDETERMINADO / F1-NEGATIVO
- **Tiempo**: ~2–3 horas

## 6. Ejecución completa (pipeline)

Para correr todas las simulaciones en orden:

```bash
set -e  # Parar si algo falla

python tools/c0_landscape.py
python tools/c0_dynamics.py
python tools/c0_battery.py
python tools/c0_null_battery.py
python tools/c0_redteam.py
python tools/c0_f1.py

echo "✓ Pipeline completo terminado"
```

**Tiempo total**: 4–6 horas (dependiendo de CPU y paralelización).

## 7. Entender los resultados

### Archivos principales

| Archivo | Contenido | Consultar para... |
|---------|-----------|-------------------|
| `results/c0_landscape/summary.json` | Energía, KKT en 75 celdas | Ver si cliques y T³ son extremales |
| `results/c0_dynamics/summary.json` | L4, clasificación de 675 estados | Contar fases DISPERSO-LOCAL vs triviales |
| `results/c0_battery/summary.json` | Dimensión D_L, D_s, homogeneidad | Saber qué familias geométricas emergen |
| `results/c0_battery/null_control.json` | Comparación vs grafo recableado | Verificar que no es solo ruido |
| `results/c0_redteam/summary.json` | 24 celdas candidatas, escalado en N | Confirmar persistencia de localidad |
| `results/c0_f1/summary.json` | Veredicto F1 y códigos de certificado | Decidir si la fase pasa como geométrica |
| `docs/OMEGA_C0_RESULTADOS.md` | Interpretación completa | Leer el análisis humano |

### Estructura de un resultado (JSONLines)

Cada línea en `runs.jsonl` es un diccionario JSON con:

```json
{
  "id": "c{cell}_U_n216_s0",          // Identificador único
  "cell_idx": 0,                       // Índice en la parrilla
  "c_star": 0,                         // Parámetro codegrado objetivo
  "k_star": 4,                         // Parámetro grado objetivo
  "a": 0.25,                           // Parámetro de offset
  "init": "U",                         // Inicio (U/E/R)
  "seed": 0,                           // Semilla RNG
  "n": 216,                            // Tamaño del grafo
  "status": "converged",               // CONVERGED / MAX_STEPS / stalled
  "steps": 1248,                       // Pasos dados
  "S_over_LB": -0.95,                  // Acción / cota inferior
  "kkt_residual": 1.2e-11,             // Residuo KKT
  "cls": {
    "class": "DISPERSO_LOCAL",         // Clasificación L4
    "H_null": 1.42,                    // Nulidad de código
    "kmax_over_kmean": 1.8             // Índice de hubs
  }
}
```

### Filtrar y explorar resultados

Con Python:

```python
import json

# Leer todas las filas
with open("results/c0_dynamics/runs.jsonl") as f:
    rows = [json.loads(line) for line in f]

# Filtrar locales
local = [r for r in rows if r["cls"]["class"] == "DISPERSO_LOCAL"]
print(f"Estados locales: {len(local)}/{len(rows)}")

# Por celda
cell18 = [r for r in rows if r["cell_idx"] == 18]
print(f"Celda 18: {len(cell18)} corridas")

# Estadísticas
dims = [r["battery"]["D_L"] for r in local if "battery" in r and r["battery"].get("D_L")]
print(f"D_L (dimensión de escalado): {sorted(dims)}")
```

## 8. Personalizar parámetros

Para modificar los parámetros de las simulaciones, edita las constantes en cada herramienta:

- `c0_landscape.py`: `MASTER`, `TOL`, `MAX_STEPS`
- `c0_dynamics.py`: `N` (tamaño), `SEEDS`, número de celdas
- `c0_redteam.py`: `SIZES` (125, 343), `D1_TOTAL` (100 000)
- `c0_f1.py`: `CELLS_512`, `CELLS_729`, `MAX_STEPS` (40 000)

Véase `docs/OMEGA_C0_PRERREGISTRO.md` para la especificación oficial de todos los parámetros.

## 9. Troubleshooting

### ImportError: omega.c0

Asegúrate de instalar en modo desarrollo:

```bash
pip install -e .
```

### Los tests fallan

Verifica el entorno:

```bash
python --version              # ≥3.11
pip show numpy scipy networkx  # Versiones correctas
```

Si fallan tests específicos de arquitectura, puede ser por el RNG o por precisión numérica. Consulta `omega/experiments/v11/`.

### Las simulaciones son muy lentas

- Reduce `--procs` si el disco está saturado
- Usa `--summarize-only` para regenerar solo estadísticas
- Verifica que los threads de BLAS estén limitados (las herramientas lo hacen automáticamente via `OMP_NUM_THREADS=1`)

### Out of memory en F1 con N=729

- Reduce `--procs` a 1–2
- Ejecuta `python tools/c0_f1.py --resume` después de un reinicio

## 10. Próximos pasos

Una vez completada la prueba C0 y revisados los resultados:

1. **Lee** `docs/OMEGA_C0_RESULTADOS.md` para la interpretación
2. **Decide**: ¿C1 (segundo orden)? ¿Θ>0 (entropía)? ¿O-05?
3. **Prerregistra** la nueva dinámica en `docs/` antes de simular
4. **Crea** los scripts correspondientes en `tools/`
5. **Ejecuta** y analiza

---

**Última actualización**: 2026-10-06  
**Rama principal**: `claude/omega-c0-competencia`  
**Estado**: C0 completo; pendiente decisión del Consejo para C1.
