# R14-C: observables causales multivariables. Preregistro exploratorio (2026-10-18)

- **Mandato:** Rev. 14, R14-3 firmada (preparar R14-C como preregistro exploratorio, **sin voto**). **Estado:** PREPARADO, **no ejecutado**. Su ejecución requiere una orden del usuario.
- **Motivo (R14-A §5 y §6):** b y K no contienen la información necesaria. Promediar r sobre todas las direcciones diluye los defectos localizados del cono (una tapa, dos puntas, pliegues). Hacen falta observables que vean **distribuciones**, no solo promedios.
- **Lo que esta fase no puede hacer:** no da veredictos. D5 a D8 y C0 a C7 son **ejemplos vistos**: sirven para explorar, no para confirmar. Una confirmación futura necesitará familias generadoras nuevas, bandas preregistradas y un guardián nuevo.

## 1. Observables candidatos (todos intrínsecos), por sub-intervalo con m ∈ {20, 40, 80, 160}

1. **Distribución de r por intervalo**, no solo su media. Un cono con una tapa plana mezcla intervalos normales con intervalos alterados según su orientación, así que la distribución debería ensancharse o volverse bimodal aunque la media apenas cambie (hipótesis, nivel 2). Medidas: varianza, asimetría, curtosis y la distancia de energía a la distribución exacta del diamante plano con m puntos.
2. **Distribución completa de tamaños de intervalos** dentro de cada sub-intervalo (todas las k, no solo k ≤ 8), normalizada y comparada con la exacta.
3. **Abundancias de cadenas** C_k (número de k-cadenas, k = 2..5) por sub-intervalo. Su valor esperado en un diamante de Minkowski de dimensión d con m puntos tiene forma cerrada (Myrheim; Meyer; nivel 3, a verificar antes de usarla). Dan estimadores de dimensión independientes de r.
4. **Cadena más larga** L(m) frente a m: en Minkowski crece como m^{1/d} con constante conocida asintóticamente (Brightwell y Gregory; nivel 3, a verificar).
5. **Dimensión de punto medio** por sub-intervalo (instrumento de C0-4), como estimador independiente.
6. **Curvatura:** pendiente A de F6 y, como contraste, el estimador de Benincasa-Dowker de la acción local por sub-intervalo (nivel 3, a verificar).

## 2. Preguntas exploratorias (sin voto)

- **P1.** ¿Qué observable separa a D6, D7 y D8 de las geometrías, a igual N? Se mide como tamaño de efecto por familia, con la realización como unidad.
- **P2.** ¿La separación de P1 se mantiene en las geometrías curvas, o la curvatura la imita, como pasó con el perfil en F7?
- **P3.** ¿Los gemelos (D5) se detectan por la distribución de r o por las cadenas, sin un detector hecho a medida?
- **P4.** ¿Qué combinación mínima de observables (2 o 3) separa todas las familias vistas? Su elección queda registrada como **hipótesis**, y la confirmación la hace otra prueba.

## 3. Estructuras y calibración (todas ya vistas)

- **Geometrías (familias):**
  - conformemente planas: diamante, de Sitter, de Sitter fuerte, anti de Sitter, FRW polvo, D2, D3 y C2;
  - Kasner: C0 y D0;
  - lineales: C1 y D1.
- **Impostores (familias):**
  - conos suaves: D6, D7, D8, C3, C4, C6 y los L^p;
  - poliédricos: cúbico y octaédrico;
  - defectos: C5, D4 y D5;
  - producto 4;
  - percolación.
- **Realizaciones y tamaño:** 10 por estructura, N ∈ {4000, 8000}, hasta 300 intervalos por celda.
- **Bandas exploratorias:** el percentil 95 por familia de geometrías. **Se informa la superposición, no se clasifica.**

## 4. Salvaguardas contra el sobreajuste

- **Se informan todos los observables y todas las comparaciones,** no solo los que separan. Con 6 observables y varias m, alguno separará por azar: lo que se vea aquí es una hipótesis.
- **Unidad:** la realización dentro de la estructura, y la familia para afirmaciones generales (R14-A §7).
- **Si de aquí sale un criterio:** se congela con bandas de equivalencia y con controles sellados de familias nuevas, definidos por un guardián que no vea estos datos.

## 5. Convenciones

- **EXPERIMENT_ID 324; entropía 20261089** (exploratoria).
- spawn_key (324, familia, índice, N, realización).
- **Costo estimado:** unas 3 a 4 horas. Las cadenas se cuentan por sub-intervalo, que es barato, porque los sub-intervalos tienen como máximo 162 elementos.
- **Antes del código:** verificar las fórmulas cerradas de cadenas, de cadena más larga y de Benincasa-Dowker con un agente aparte y anotar su fuente. Si alguna no se verifica, se usa una referencia por Monte Carlo del diamante plano.
- **Herramienta:** a cargo de un agente de código, con handoff y auditoría; árbol limpio; tandas de menos de 2 horas.
