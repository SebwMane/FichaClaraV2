# F5′: resultados (perfil de intervalos con dimensión espectral, 2026-10-14)

- **Especificación:** `OMEGA_CONSEJO_REV9_F5PRIMA.md` §6 (R9-1 y R9-2 firmadas; preregistro en 1d8c770). **Herramienta:** `tools/f5prima.py` (agente de código, auditoría §21). **Corrida oficial:** commit 6d48b67, árbol limpio. **Datos:** `results/omega_f5prima/`.

## 1. Veredicto oficial: **FALLA**

Las cuatro positivas quedaron bien clasificadas, pero **el cono octaédrico 4D fue aceptado como espacio-tiempo 4D** (razón 0.92 ≤ 1). Las demás negativas fueron rechazadas, el cono cúbico por muy poco (razón 1.03).

## 2. Detalle (ET = clasificado como espacio-tiempo)

| Estructura | Papel | d_spec | d* | Firma | Razón D/δ | ET | ¿Correcto? |
|---|---|---|---|---|---|---|---|
| Diamante 4D | positiva | 3.96 | 4 | sí | 0.24 | sí | sí |
| Diamante 3D | positiva | 2.99 | 3 | sí | 0.25 | sí | sí |
| de Sitter 4D (a ciegas) | positiva | 3.97 | 4 | sí | 0.89 | sí | sí |
| de Sitter 3D (a ciegas) | positiva | 2.99 | 3 | sí | 0.56 | sí | sí |
| Cono cúbico 4D | negativa | 3.73 | 4 | sí | **1.03** | no | sí (al límite) |
| **Cono octaédrico 4D** | negativa | 3.83 | 4 | sí | **0.92** | **sí** | **no** |
| Orden producto 3 | negativa | 2.83 | 3 | sí | 3.11 | no | sí |
| Orden producto 4 | negativa | 3.63 | 4 | sí | 2.13 | no | sí |
| Orden producto 5 | negativa | 3.47 | 3 | sí | 14.9 | no | sí |
| KR | negativa | 1.79 | 2 | no | 193 | no | sí |
| Percolación p = 0.05 | negativa | 2.30 | 2 | sí | 11.5 | no | sí |
| Capas gruesas | negativa | (fuera de rango) | 2 | no | 193 | no | sí |

Umbrales: δ(4, 1600) = 0.174, δ(4, 2400) = 0.160; δ(3, ·) = 0.101 y 0.071. Brazos sin voto: A y B también aceptan a los dos conos poliédricos.

## 3. Lectura

**Lo que funcionó (nivel 1, en este alcance).** Las dos correcciones de F5 hicieron lo que debían. La dimensión espectral asignó bien la dimensión a las cuatro positivas, como anticipó la cuenta previa (§6.1). Con la calibración curva, de Sitter 4D pasó (0.89, antes 1.21), y de Sitter 3D también, ambos sin haber entrado en la calibración. Los tres órdenes producto siguen claramente rechazados.

**Lo que falló (nivel 1).** El perfil de intervalos **no distingue un cono redondo de uno poliédrico** con este umbral. Las distancias crudas lo muestran (N = 2400):

| | D típica |
|---|---|
| Diamantes planos de calibración | 0.04 (p95: 0.12) |
| de Sitter 4D | 0.14 a 0.20 |
| Cono cúbico | 0.16 a 0.21 |
| Cono octaédrico | 0.15 a 0.21 |
| Orden producto 4 | 0.31 a 0.39 |

Los conos poliédricos quedan tan lejos de la referencia como de Sitter, y ambos dentro del margen que deja el factor 1.5. **Lo que el perfil rechazaba en F5 era el cono de un orden producto (un ortante, cuya sección espacial es un tetraedro), no cualquier cono no redondo.** La lectura de F5 («detecta la forma del cono») queda rebajada: detecta conos muy distintos del redondo, no conos parecidos.

**Observación posterior a ver los datos (nivel 2, sin voto).** Las desviaciones tienen formas distintas. Los conos poliédricos tienen unos 17% menos links (k = 0) y casi el mismo resto del perfil; de Sitter se desplaza parejo, entre 9% y 17% en todos los k. La distancia máxima no ve esa diferencia de forma. Un criterio que mire la forma del perfil (por ejemplo, la desviación en k = 0 relativa al resto) podría separarlos, pero eso se propuso **después** de ver estos datos: solo puede probarse en una prueba nueva, preregistrada y con semillas nuevas, y no cambia este veredicto.

**Riesgo anotado.** En de Sitter 4D la distancia crece con N (mediana 0.12 en N = 1600 y 0.17 en N = 2400), así que su aceptación podría no sostenerse con N mayor.

## 4. Lo que deja F5′

- **Ningún criterio del programa reconoce todavía espacio-tiempo lorentziano en 4D frente a todos los impostores.** La combinación firma espectral + perfil de intervalos reconoce «orden continuo de dimensión d con conos cercanos a los redondos», incluidos los curvos, y rechaza órdenes producto, KR, percolación y capas. No separa conos poliédricos.
- **F6** (bordes y topología) sigue comprometida por R9-2.
- **Posible F7** (no preparada; requiere al Consejo y la firma del usuario): un criterio de forma del perfil, con los conos poliédricos como impostores de calibración ya vistos y otros conos nuevos a ciegas.
