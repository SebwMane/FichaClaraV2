# Ω-C en 4D: resultados de L4b (bloque con espacio toroidal)

- **Fecha:** 2026-10-12. **Commit:** `654ff70`. **Especificación:** `OMEGA_2_PROPUESTA_CAUSAL.md` §20. **Datos:** `results/l4b/`.

## 1. C4b: instrumento en el toro · ÉXITO

ε = 0.25: 39.05 ± 2.61 frente a 40.17 ± 2.05; ε = 0.05: 358.06 ± 0.66 frente a 358.17 ± 0.47 (N = 512).

## 2. L4b · **INCOMPATIBLE**

D_int = sd(acción de la mitad superior del bloque)/ε; se necesitaba ≤ 0.477 sin crecer.

| esquema | N = 256 | 512 | 1024 | 2048 | η_int | η_tot |
|---|---|---|---|---|---|---|
| ε = 1 | 120 | 325 | 811 | 1 747 | 1.29 ± 0.05 | 1.30 ± 0.05 |
| ε = 0.25 | 72 | 145 | 322 | 553 | 1.00 ± 0.05 | 1.02 ± 0.05 |
| ε = 0.05 | 170 | 268 | 475 | 970 | 0.84 ± 0.06 | 1.13 ± 0.05 |
| ε = 16/N | 142 | 389 | 972 | 3 207 | 1.48 ± 0.06 | 1.37 ± 0.08 |
| ε = 64/N | 72 | 165 | 434 | 1 270 | 1.38 ± 0.05 | 1.31 ± 0.05 |

## 3. Lectura

- **El borde espacial no era la causa (nivel 1 en lo medido).** Sin bordes laterales y midiendo solo el interior, la dispersión es del mismo orden que en los diamantes y crece igual.
- **No es una suma de términos independientes (nivel 1).** En ningún esquema η_int es compatible con 0.5; crece entre N^0.84 y N^1.48. Las contribuciones locales s(x) están correlacionadas a larga distancia.
- **Explicación propuesta (nivel 2).** La no localidad de la acción discreta: cada elemento tiene links a lo largo de todo su cono de luz pasado, así que los s(x) comparten información a grandes distancias en direcciones nulas y sus fluctuaciones no se promedian.
- **Límite declarado:** la acción «interior» sigue incluyendo pares casi nulos que llegan hasta t = 0. Eliminar eso del todo exige otra geometría (por ejemplo, una región mucho más larga en el tiempo), con coste mayor.

## 4. Conclusión para la rama lorentziana 4D

En diamantes y en bloques con espacio toroidal, para los siete esquemas de ε probados, la geometría no puede ser coherente donde las familias por capas se cancelan. La rama lorentziana queda **incompatible en 4D en todo lo probado**. Esto no refuta la teoría de conjuntos causales en general: se ajusta a un problema conocido de la literatura (las fluctuaciones de la acción y el papel de la no localidad), y lo cuantifica con datos propios.

> **Nota de unidades (2026-10-19, auditoría de √ε, `docs/omega_c/AUDITORIA_SQRT_EPS.md`):** la acción 4D suavizada de esta corrida está en unidades del repositorio de entonces, sin el factor √ε (S_ε/√ε). D₄ y D_int son invariantes ante ese factor, porque numerador y denominador cambian igual, así que los veredictos no cambian. Desde el commit de la corrección, `bd_action_4d_smeared` incluye √ε.
