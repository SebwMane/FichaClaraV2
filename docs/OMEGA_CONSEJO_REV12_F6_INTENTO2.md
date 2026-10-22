# Consejo Ω, Revisión 12: propuesta del intento 2 de F6 (2026-10-17)

- **Estado:** PROPUESTA, no preregistrada ni ejecutada. Pedido del usuario: «continúa con tu propuesta para el intento 2 de F6».
- **Antecedente:** `OMEGA_F6_CONFIRMACION_RESULTADOS.md`. El intento 1 FALLÓ solo por C5 (Minkowski con 13% de links borrados): el intercepto lo leyó como curvatura negativa.
- **Presupuesto:** es el **último intento** de la línea del perfil local (Rev. 11 §6). Si falla, la línea vuelve al Consejo y no se prueba otra variante sin una revisión completa.
- **Figura:** `docs/figuras/rev12_colapso.png`.

## 1. La idea: la prueba de colapso, sin modelo

En la rejilla de F6, el tamaño relativo √(m/N) se repite exactamente entre los dos N:

| Par | N = 4000 | N = 8000 | √(m/N) |
|---|---|---|---|
| 1 | m = 10 | m = 20 | 0.050 |
| 2 | m = 20 | m = 40 | 0.071 |
| 3 | m = 40 | m = 80 | 0.100 |
| 4 | m = 80 | m = 160 | 0.141 |

**Predicción de H_loc (Rev. 11 §2, condición 3, escrita antes de C5):** en una geometría, la desviación depende solo del tamaño relativo, así que los dos miembros de cada par deben dar el mismo Δr, **cualquiera sea la curvatura y sin ajustar ninguna función**. C5 no cumple eso: su defecto depende de m, no de m/N.

**Estadístico:**

  K = Σ_pares (d_j / se_j)², d_j = Δr(2m, 8000) − Δr(m, 4000),

con se_j calculado entre realizaciones. Una geometría da K pequeño (cerca de 4 si no hay sesgos); un defecto como C5 da K grande.

**Honestidad sobre el origen:** la condición estaba en H_loc desde Rev. 11, pero **la decisión de hacerla votar nace de ver C5**. Por eso C0 a C7 pasan a ser «vistas» y el intento 2 necesita un conjunto sellado nuevo.

## 2. Cuenta previa sobre datos ya vistos (confirmación del intento 1; declarada, no es evidencia)

1. **Ruido compartido en la referencia.** La referencia empírica r_plano(20) de la confirmación salió 0.0980 (−2.5σ; azar, no un defecto del código), y ese error se restó igual en todas las estructuras. Además fabricó un patrón falso en la prueba de colapso: los diamantes planos de la calibración llegaban a K = 25.
2. **Mejora:** en un diamante 4D de Minkowski, la fracción de pares relacionados entre m puntos uniformes vale **exactamente 1/10** para todo m (dos veces la fracción de Myrheim-Meyer f(4) = 1/20; nivel 3). Con la referencia exacta, la calibración da K entre 0.2 y 19.3.
3. **Con 10 realizaciones, C5 no se atrapa:** K = 23.4, frente a un umbral de 1.5 × 19.3 ≈ 29. Con 30 realizaciones por N, el error baja a la mitad y K de C5 subiría a unos 70; el de una geometría sin sesgo seguiría cerca de 4 (nivel 2, extrapolación simple).
4. **Riesgo declarado: FRW de polvo**, que es positiva, dio K = 14.6 con los pares 2 a 4 todos del mismo signo. Puede ser azar o un efecto real (en FRW la curvatura cambia con el tiempo y la mezcla de lugares de los intervalos podría depender de N). Con más realizaciones, si es real, crecerá y la positiva fallaría. El intento 2 lo pone a prueba; no se corrige antes.
5. **Los conos no se atrapan con el colapso** (K entre 2.9 y 10.8): para ellos sirve el intercepto. Las dos condiciones se complementan.

## 3. Diseño propuesto

- **Clasificación «espacio-tiempo 4D»:** las dos condiciones a la vez.
  1. **Intercepto:** |b − b_plano| ≤ τ_b, como en el intento 1, pero con la referencia exacta 1/10.
  2. **Colapso:** K ≤ τ_K, con τ_K = 1.5 × max de K sobre los grupos de calibración, la misma regla que τ_b.
  - «No medible» (menos de dos x distintos, o menos de dos pares medibles) = no espacio-tiempo.
- **Realizaciones:** **30 por N** (antes 10). Se mantienen la rejilla de m, N ∈ {4000, 8000}, la regla m ≤ 0.1·m_max y el máximo de 300 intervalos.
- **Calibración:** 15 grupos (10 planos, 3 de Sitter, 2 FRW polvo), cada uno con 30 + 30 realizaciones, y entropía nueva.
- **Estructuras:**
  - un **conjunto sellado nuevo**, definido por un **guardián nuevo** que no vio nada de este diseño;
  - las vistas con semillas nuevas, que incluyen ahora **C0 a C7** (C5 como negativa vista).
- **Requisitos para el guardián nuevo:**
  - al menos una positiva no conformemente plana distinta de C0;
  - al menos un **defecto discreto distinto** del borrado de links (por ejemplo, relaciones agregadas o puntos duplicados);
  - una positiva curva con curvatura que varíe en el tiempo, para probar el riesgo de FRW;
  - negativas cerca de la frontera de la zona gris;
  - ninguna positiva en la zona gris.
- **Informado sin voto:** pendiente A con su signo; rama 2.
- **Costo estimado:** unas 6 horas de cómputo, repartidas en varios comandos de menos de 2 horas.
- **El conjunto de reserva** sigue sellado; solo se abre si este intento PASA.

## 4. Qué diría cada resultado

- **PASA:** primer criterio de la línea que reconoce espacio-tiempo 4D frente a todos los impostores probados, incluidos defectos discretos y estructuras selladas. Habilita abrir la reserva con N = 8000.
- **FALLA por una positiva curva en el colapso** (por ejemplo FRW): la condición «solo depende de m/N» es demasiado estricta para curvatura que varía, y H_loc necesita revisarse.
- **FALLA por un defecto nuevo:** hay impostores que imitan a la vez la escala y la forma; la línea vuelve al Consejo con esa pregunta.

## 5. Decisiones del usuario (pendientes)

| # | Pregunta |
|---|---|
| R13-1 | ¿Agregar la prueba de colapso como segunda condición con voto, junto al intercepto? |
| R13-2 | ¿30 realizaciones por N (unas 6 h) o 20 (unas 4 h, menos potencia: K de C5 cerca de 47)? |
| R13-3 | ¿Convocar un guardián nuevo para un conjunto sellado nuevo, con los requisitos de §3? |
| R13-4 | ¿Confirmas que este es el último intento y que, si falla, la línea vuelve al Consejo sin probar otra variante? |
