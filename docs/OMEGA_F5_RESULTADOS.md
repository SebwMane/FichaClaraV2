# F5: resultados (coherencia entre medidas de dimensión, 2026-10-12)

- **Especificación:** `OMEGA_CONSEJO_REV8_F5.md` §3 y §6 (R8-1 y R8-2 firmadas). **Datos:** `results/omega_f5/`.
- Tolerancias calculadas por la regla: τ_A = 0.73; τ_B = 0.363.

## 1. Veredictos oficiales

| Brazo | Papel | Veredicto |
|---|---|---|
| A, Myrheim-Meyer frente a dimensión espectral | vota | **FALLA** |
| C, perfil de intervalos | vota | **FALLA** |
| B, Myrheim-Meyer frente a punto medio | informativo | FALLA |

## 2. Detalle (ET = clasificado como espacio-tiempo; ✗ = clasificación incorrecta)

| Estructura | Papel | d* | Firma | A | C (razón D/δ) | B |
|---|---|---|---|---|---|---|
| Diamante 4D | positiva | 4 | sí | ET | ET (0.34) | ET |
| Bloque 4D sin borde espacial | positiva | **5** | no | ✗ | ✗ (4.51) | ✗ |
| de Sitter 4D | positiva | 4 | sí | ET | ✗ (**1.21**) | ET |
| Diamante 3D | positiva | 3 | sí | ET | ET (0.43) | ET |
| Orden producto 3 | negativa | 3 | sí | ET ✗ | **rechazado (4.93)** | ET ✗ |
| Orden producto 4 | negativa | 4 | sí | ET ✗ | **rechazado (2.37)** | ET ✗ |
| Orden producto 5 | negativa | 5 | sí | rechazado | **rechazado (2.48)** | rechazado |
| KR | negativa | 2 | no | rechazado | rechazado | rechazado |
| Percolación p = 0.05 | negativa | 2 | sí | rechazado | rechazado (12.0) | ET ✗ |
| Capas gruesas | negativa | 2 | no | rechazado | rechazado | rechazado |

## 3. Lectura

**Brazo A (nivel 1).** Falló como anticipó la cuenta previa del Consejo: los órdenes producto de dimensión 3 y 4 tienen dimensiones de Myrheim-Meyer y espectral casi iguales, y pasan como espacio-tiempo.

**Brazo C (nivel 1 descriptivo; no cambia el veredicto).** **Rechazó a las seis negativas, incluidos los tres órdenes producto**, con razones de 2.4 a 4.9 veces el umbral, mientras los dos diamantes (3D y 4D) quedaron muy por debajo (0.34 y 0.43). Es el único instrumento del programa que distingue la forma de los conos. Falló en dos positivas, por razones distintas:

1. **Bloque 4D: un error de diseño del preregistro.** La dimensión asignada d* se tomó de Myrheim-Meyer, que supone que la región completa es un intervalo causal (un diamante). El bloque no lo es, y Myrheim-Meyer le asigna 5 en vez de 4; desde ahí la firma y el perfil se comparan con la dimensión equivocada. El mismo error afecta a los tres brazos.
2. **de Sitter 4D: el umbral quedó justo.** Razón 1.21: la curvatura cambia el perfil de intervalos más de lo que admite un umbral calibrado solo con espacio plano, y ese umbral salió del percentil 95 de apenas 10 valores (declarado como debilidad en §6.1).

**Brazo B (informativo).** Separa poco: los órdenes producto 3 y 4 y la percolación pasan.

## 4. Lo que deja F5

- **Hallazgo principal (nivel 1, descriptivo):** el perfil de intervalos detecta la estructura de los conos. Combinado con la firma espectral, es el candidato para reconocer espacio-tiempo en 4D.
- **Dos correcciones necesarias**, ambas identificadas aquí y no aplicadas a este resultado: medir dentro del mayor intervalo causal de cada estructura (para que Myrheim-Meyer y el perfil se apliquen donde sus supuestos valen), y calibrar el umbral del perfil con más muestras, incluyendo geometrías curvas.
- **Siguiente prueba posible (no preregistrada): F5′**, con esas dos correcciones y estructuras de prueba nuevas.
