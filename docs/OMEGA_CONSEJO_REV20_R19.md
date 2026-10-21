# Consejo Ω, Revisión 20: revisión de R19-(b) (sesión 9, 2026-10-20)

- **Acta:** el usuario la pegó completa. Los dictámenes individuales están en rutas locales del usuario, fuera de este repositorio. Aquí se registran la resolución, las propuestas por perfil y la verificación del orquestador.
- **Pedido del usuario:** «Que el Consejo revise el resultado antes de repetir».

## 1. Resolución (unánime, 8 perfiles)

- **Se mantiene el veredicto oficial: INDETERMINADO.** Las reglas congeladas se aplicaron bien. Hay dos causas, cada una suficiente por sí sola:
  - la cobertura con q = 1/8 fue de 14 de 20, por debajo del mínimo de 16;
  - A2a y A2b son ambiguas, porque ya estaban dentro de la banda con q = 1, y así L no tiene señal de partida.
- **Hallazgos de nivel 2, secundarios:**
  - las 4 estructuras selladas de copias y D5 se separan;
  - A1a a A1c, D7, D8, plana y de Sitter quedan equivalentes.

  La lectura defendible es condicional: Φ_q parece distinguir, en las estructuras visibles, copias que se desvanecen de controles a densidad fija. **No prueba que Ω produzca esas estructuras ni que seleccione geometría.**
- **Recomendación:** preparar una repetición confirmatoria corregida, con preregistro, guardián, familias y semillas nuevos. **No autoriza ejecutarla.**

## 2. Recomendaciones para la repetición (resumen del acta)

1. No retocar R19-(b) después de verlo: no bajar 16/20 a 15/20 ni retirar A2.
2. **Cobertura:**
   - separar el ajuste de media y covarianza, la calibración del radio y la validación independiente;
   - fijar el tamaño, la cobertura objetivo y una regla binomial simultánea para los 4 q, que tenga en cuenta el riesgo conjunto;
   - si falla, el resultado es indeterminado, sin recalibrar con el holdout.
3. **Población identificable:**
   - el guardián aplica una regla numérica de visibilidad con q = 1, con semillas de selección distintas de las confirmatorias, y sobre un universo de candidatos fijado de antemano;
   - los excluidos se informan como «no identificables por este instrumento», no como equivalentes;
   - la conclusión queda limitada a estructuras detectables.
4. Si se quieren incluir los autosemejantes invisibles, primero hace falta otro estadístico, calibrado en un piloto independiente.
5. **Unidad:** la realización es la unidad, y los 4 q son medidas repetidas. Para generalizar a un mecanismo hacen falta varias estructuras o generadores por mecanismo. Gemelos (q²) y tripletes (q³) van separados.
6. **Condiciones por estructura antes de sellar:**
   - la estructura votante es medible con q = 1;
   - los controles reinyectados siguen visibles en la escala de decisión;
   - los autosemejantes detectables pueden cumplir la regla de equivalencia.

   Se congelan el margen, L, los intervalos, las exclusiones, la multiplicidad y el trato de los casos no evaluables.

**Propuestas por perfil:**
- **Geometría causal:** A2 como no identificable, con criterio previo y conteo de excluidos.
- **Geometría causal (2):** evitar que la selección por visibilidad sesgue la curva emparejada; usar semillas de selección separadas.
- **Dinámica emergente:** el resultado es una propiedad de Φ_q sobre órdenes dados.
- **Dinámica emergente (2):** verificar la visibilidad de estructuras y controles en las escalas que importan.
- **Ensambles:** gemelos y tripletes separados; q² y q³ orientan la predicción, pero no fijan el orden de los efectos.
- **Estadística:** no cambiar 16/20 después de ver 14/20; el ajuste y la validación son incertidumbres distintas.
- **Epistemología:** mantener separados el veredicto formal y el patrón de nivel 2.
- **Protocolo:** guardián nuevo; registrar versiones, logs y salidas. El hash acredita la integridad, no la ausencia de filtraciones.

## 3. Verificación del orquestador

- **Probabilidades binomiales del acta, recalculadas con cobertura de 0.90 y n = 20:**
  - P(menos de 16) = 0.0432;
  - P(que falle al menos una de 4 escalas) = 0.162;
  - P(14 o menos) = 0.0113.

  Coinciden con el acta. Mi aclaración A1 decía «del 13% al 16%»; el valor con cobertura de 0.90 es 16.2%.
- **Una nota que el acta no menciona:** observar 14 de 20 o menos tiene probabilidad 1.1% si la cobertura real es del 90%. Es compatible con azar (se miraron 4 escalas), pero también con que la cobertura real con q = 1/8 (N ≈ 1000, con una componente menos en el vector) sea menor que 90%. La repetición debe estimarla con más validación, como recomienda el punto 2.

## 4. Firma pendiente

| # | Decisión |
|---|---|
| R20-1 | ¿Autorizar la **preparación** del preregistro corregido, sin autorizar su ejecución? |
| R20-2 | Alcance: **(A) acotado** a estructuras autosemejantes detectables con q = 1, con filtro previo y excluidos informados como no identificables; o **(B) amplio**, que incluya los invisibles y por eso exige antes desarrollar y calibrar otro estadístico en un piloto aparte. |

## 5. Firmas del usuario (2026-10-21)

- **R20-1 y R20-2 (opción A, acotada) firmadas.**
- Preregistro preparado: `docs/OMEGA_R20_REPETICION_PREREGISTRO.md`. No ejecutado. La ejecución requiere el sellado del guardián nuevo y una orden explícita.
