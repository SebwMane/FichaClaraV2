# R19-(b): convergencia bajo dilución con control reinyectado. Resultados (2026-10-20)

- **Preregistro:** `OMEGA_R19_GEMELOS_PREREGISTRO.md`, con la aclaración A1 (§9).
- **Sellado:** `docs/SELLADO_R19.md` (sha256 e3c0e40c…, verificado al abrirlo).
- **Código:** ejecutor congelado en f512b02 y archivo sellado abierto en f7f9244. La corrida usó f7f9244 con el árbol limpio: 600 unidades en 49 min.
- **Datos:** `results/omega_r19/`. **Figura:** `docs/figuras/r19.png`.

## 1. Veredicto oficial: **INDETERMINADO**, por dos razones independientes

1. **Falló la calibración de cobertura con q = 1/8:** 14 de 20 sprinklings planos nuevos quedaron en banda, y el mínimo era 16. Es el riesgo que declaré en A1 antes de correr (una probabilidad del 13% al 16% de que falle por azar). Con q = 1, 1/2 y 1/4, la cobertura fue de 16, 19 y 19 de 20.
2. **Aun con la cobertura bien, el veredicto no habría sido PASA:** las dos estructuras selladas de la familia A2 (conos suaves) quedaron **ambiguas**, porque sus IC cruzan el borde del margen ±ln 1.3:
   - A2a: L = 0.22, con IC de 0.05 a 0.39;
   - A2b: L = −0.03, con IC de −0.29 a 0.09.

No se reinterpreta nada: el resultado oficial es INDETERMINADO.

## 2. Datos (L₁/₈ = ln(D_obs/D_ctrl), mediana e IC 95%; margen ±0.262)

| Estructura | Grupo | Papel | L₁/₈ | Clase |
|---|---|---|---|---|
| plana | calibración | nula | 0.06 (−0.19 a 0.22) | equivalente |
| de Sitter | calibración | nula | 0.05 (−0.09 a 0.15) | equivalente |
| D7 | vista | autosemejante | −0.01 (−0.11 a 0.22) | equivalente |
| D8 | vista | autosemejante | −0.08 (−0.18 a 0.08) | equivalente |
| D5, gemelos al 10% | vista | se desvanece | −0.62 (−0.78 a −0.48) | separa |
| S1a, gemelos | **sellada** | se desvanece | −0.95 (−1.09 a −0.88) | separa |
| S1b, gemelos | **sellada** | se desvanece | −1.44 (−1.50 a −1.27) | separa |
| S2a, tripletes | **sellada** | se desvanece | −0.76 (−0.86 a −0.66) | separa |
| S2b, tripletes | **sellada** | se desvanece | −1.02 (−1.18 a −0.93) | separa |
| A1a, cono tetraédrico | **sellada** | autosemejante | 0.06 (−0.09 a 0.10) | equivalente |
| A1b, cono dodecaédrico | **sellada** | autosemejante | 0.05 (−0.12 a 0.19) | equivalente |
| A1c, cono rombododecaédrico | **sellada** | autosemejante | −0.03 (−0.21 a 0.14) | equivalente |
| A2a, cono suma de elipsoides | **sellada** | autosemejante | 0.22 (0.05 a 0.39) | ambiguo |
| A2b, cono esfera sesgada | **sellada** | autosemejante | −0.03 (−0.29 a 0.09) | ambiguo |

## 3. Lectura (nivel 2, sin cambiar el veredicto)

- **Las cuatro estructuras selladas de copias** (gemelos y tripletes) y D5 se separan con claridad: L va de −0.76 a −1.44, con IC enteramente por debajo de −0.262. **El control reinyectado sí distingue un defecto que se desvanece de uno que solo deja de verse.** En la figura, la curva del observado cae muy por debajo de la del control desde q = 1/2.
- **Los conos visibles** (A1a a A1c, D7 y D8) y las dos calibraciones quedan equivalentes: observado y control coinciden en todo q, como debe pasar con estructuras autosemejantes.
- **Por qué A2 quedó ambigua:** A2a y A2b son conos **casi invisibles**. Su D/radio está por debajo de 1 en todo q, incluso con q = 1 (de 0.72 a 0.94). Entonces L es el cociente de dos distancias que están dentro de la banda, y es puro ruido. El guardián verificó que los **controles de copias** quedaran fuera de banda, pero el preregistro no le exigía lo mismo a las estructuras autosemejantes. **Es un hueco de diseño mío:** una estructura autosemejante que el instrumento no ve no puede poner a prueba la equivalencia.
- **La cobertura con q = 1/8** falló por la razón que se declaró: el radio A1 cubre alrededor del 90% de los puntos nuevos.

## 4. Qué implica

1. **El método funciona donde puede ver.** Para defectos visibles, el control reinyectado separa «se desvanece» de «autosemejante» de forma clara y en el sentido predicho, con familias que no se habían visto. Es un resultado de nivel 2: no vota, porque el veredicto oficial es INDETERMINADO.
2. **Para una confirmación** hacen falta dos correcciones de diseño antes de repetir, ambas visibles en este diseño y ninguna escogida mirando qué familia pasó:
   - **(i)** exigir que **toda** estructura sellada (también las autosemejantes) quede fuera de la banda con q = 1, verificado por el guardián con sus propias semillas;
   - **(ii)** fijar la cobertura de calibración con un criterio coherente con el radio, por ejemplo 15 de 20 cuando el radio cubre el 90%, o calcular el radio con más referencias.

   Una repetición necesita guardián nuevo, familias nuevas y semillas nuevas.
