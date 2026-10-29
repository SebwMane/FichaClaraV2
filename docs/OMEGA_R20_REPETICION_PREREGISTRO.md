# R20: repetición corregida de R19-(b), alcance acotado (A). Preregistro confirmatorio (2026-10-21)

- **Mandato:** el usuario firmó el 2026-10-21 «R20-1 y R20-2 opción A, prepara el preregistro». Se basa en las recomendaciones del Consejo, sesión 9 (`docs/OMEGA_CONSEJO_REV20_R19.md`).
- **Versión 2 (2026-10-21):** incorpora la enmienda de la Rev. 21 (cobertura, visibilidad, emparejamiento e inferencia), **pendiente de la firma explícita del usuario** antes de sellar. Las diferencias con la versión 1 (commit 5068497) quedan en el historial de git.
- **Estado:** **PREPARADO, no ejecutado.** Para ejecutarlo hacen falta el sellado de un **guardián nuevo** y una orden explícita del usuario.
- **R19-(b) no se toca:** sigue INDETERMINADO. Este es un experimento nuevo, con familias, guardián y semillas nuevos.
- **Alcance (opción A):** la conclusión vale solo para estructuras **identificables por el instrumento** (§3). Las que no lo son se informan como «no identificables por este instrumento», nunca como equivalentes. La conclusión es sobre cómo Φ_q transforma órdenes sintéticos. **No** dice que Ω los produzca, ni resuelve H_E, H_L ni H_sel.

## 1. Pregunta y método (los mismos de R19-(b), sin cambios)

- **Hipótesis H_conv:** con el control a densidad fija, un defecto que se desvanece bajo Φ_q se aleja de su control, y uno autosemejante no.
- **Control (versión 2, R21-3(a)): emparejado por N.** N es el tamaño total del orden observado ya diluido, Φ_q(X).
  - **Copias:** n_copias = k·round(frac·N/k). La base del control es un subconjunto uniforme, de tamaño N − n_copias, de la base de la misma realización, con el orden inducido. Encima se agregan copias nuevas con la misma k y la misma frac, y con un generador auxiliar aparte.
  - **Regiones:** el control es generar(N, rng_nuevo), con el mismo N.
  - Con q = 1, N = N₀.
- **Estadístico:** L = ln(D_obs/D_ctrl) con q = 1/8, usando el vector de 14 componentes y la distancia de Mahalanobis con Ledoit-Wolf de R14-B v2.
- **Margen de equivalencia:** ±ln 1.3.
- **Unidad:** la realización, con **30 realizaciones** por estructura. Los cuatro q son medidas repetidas de la misma realización.
- **Escala:** N₀ = 8000 y q ∈ {1, 1/2, 1/4, 1/8}.

## 2. Corrección 1: la banda y su cobertura (causa 1 del INDETERMINADO de R19)

Por cada q se generan **tres conjuntos disjuntos** de sprinklings planos, con semillas propias:

| Conjunto | Tamaño | Uso |
|---|---|---|
| Ajuste | 100 | media y covarianza con Ledoit-Wolf, y estandarización |
| Calibración | 100 | **radio conforme**: el estadístico de orden número ⌈0.95·101⌉ = 96 de las distancias de calibración. Garantiza una cobertura marginal ≥ 0.95 para un punto nuevo intercambiable. |
| Validación | 100 | solo para verificar la cobertura; **nunca** se usa para recalibrar |

**Regla de cobertura simultánea (versión 2, enmienda de la Rev. 21; la cifra de la versión 1 era errónea):** se exige que al menos **87 de 100** sprinklings de validación queden en banda **en cada uno de los 4 q**.
- **Qué garantiza el radio:** el rango 96 de 100 da una cobertura **marginal** ≥ 96/101 para un punto nuevo intercambiable, sin empates. No fija la cobertura de la banda ya calibrada, que es aleatoria y sigue Beta(96, 5). Como las 100 validaciones comparten el radio, su conteo sigue una **Beta-Binomial(100, 96, 5)**.
- **Riesgo de que el gate falle por azar:** 0.0119 por q. Que falle en algún q: 0.0468 si los q son independientes, y como mucho 0.0476 por la cota de unión.
- **Lo que la versión 1 decía mal:** daba 0.045 con un gate de 90/100. Esa cifra suponía una cobertura condicional fija de 0.95. Con el modelo correcto, el gate de 90/100 fallaba con 0.193.
- **Sensibilidad del diagnóstico:** el gate de 87/100 detecta una cobertura real de 0.85 en un q con probabilidad 0.65, y una de 0.88 con 0.31. Es un **diagnóstico** y no certifica una cobertura mínima.
- **Empates de distancia** en el corte de calibración: se rompen por el índice de la referencia, de forma determinista.

Si la regla falla en algún q, el resultado es **INDETERMINADO**, sin recalibrar.

## 3. Corrección 2: población identificable (causa 2)

- **Universo fijado de antemano por el guardián:** a lo sumo 16 estructuras candidatas, repartidas en al menos 4 familias por mecanismo:
  - **gemelos** (k = 1) y **tripletes** (k = 2), como familias separadas;
  - **conos poliédricos** y **conos suaves**, también separados.
- **Cada familia** debe tener al menos 3 estructuras candidatas.
- **Requisitos de visibilidad**, evaluados por el guardián con **semillas de selección propias** (entropía 20261114), distintas de las confirmatorias, sobre 10 realizaciones de selección y con su propia banda calibrada igual que en §2. Las medianas se calculan **por brazo**, con la banda plana del mismo q; se informan cada cantidad, la semilla, la banda y N:

| Tipo | Requisito para ser identificable |
|---|---|
| Copias («se desvanece») | mediana de D/radio del **observado con q = 1** ≥ 2.0, y mediana de D/radio del **control con q = 1/8** ≥ 1.5 |
| Región («autosemejante») | mediana de D/radio del **observado con q = 1** ≥ 2.0, y con **q = 1/8** mediana ≥ 1.25 en **los dos brazos**, observado y control (versión 2, R21-2) |

- **Qué se sella:** el guardián sella el universo completo, con el resultado de cada requisito por estructura. **Solo votan las identificables.**
- **Mínimo para votar:** al menos 2 familias de cada tipo, con al menos 2 estructuras identificables cada una. Si no se llega, el resultado es INDETERMINADO, por población insuficiente.
- **Informe obligatorio:** el número de excluidas por familia y su lista, como «no identificables por este instrumento». Las excluidas se miden igual en la corrida confirmatoria, se informan y no votan.

## 4. Reglas de decisión (congeladas)

- **Por estructura votante (versión 2, R21-4):** la mediana de L₁/₈ sobre 30 realizaciones.
  - **IC 95%:** bootstrap percentil de 2000 réplicas, que remuestrea **realizaciones completas** (observado y control emparejados, los 4 q juntos).
  - **Alcance de los IC:** son condicionales al instrumento ajustado (las referencias no se remuestrean).
  - **Valores no finitos:** se tratan como límites extendidos, sin descartarlos.
    - D_obs = 0 < D_ctrl da L = −∞.
    - D_ctrl = 0 < D_obs da L = +∞.
    - Esos valores entran en la mediana y en el bootstrap como valores ordenables.
    - Si en alguna realización D_obs = D_ctrl = 0, o si L no está definido por otro motivo (por ejemplo, un observable no medible en algún brazo), la estructura queda **no evaluable**. No hay exclusiones selectivas.
- **Calibración** (si falla, el resultado es INDETERMINADO):
  - se cumple la regla de cobertura de §2;
  - plana y de Sitter tienen el IC de L₁/₈ dentro de ±ln 1.3. Plana queda dentro de la banda, así que solo verifica el procedimiento; de Sitter es visible.
- **PASA:**
  - en todas las estructuras votantes de «se desvanece», el IC queda entero por debajo de −ln 1.3;
  - en todas las votantes «autosemejantes», el IC queda entero dentro de ±ln 1.3.
- **FALLA:** es una contradicción **observada** en el conjunto, sin afirmar control del error familiar (versión 2, R21-4): alguna votante de «se desvanece» tiene el IC entero por encima de −ln 1.3, o alguna votante «autosemejante» tiene el IC entero fuera del margen.
- **INDETERMINADO:** en cualquier otro caso.
- **Multiplicidad:** PASA exige que todas las estructuras votantes cumplan a la vez (prueba de intersección-unión). Eso controla el falso PASA sin corrección adicional, pero **no** controla los falsos FALLA. Se informa cuántas estructuras votaron.
- **Informe:** por estructura y por familia, con L_q en los cuatro q y las trayectorias de D_obs y D_ctrl. Se separan siempre el veredicto formal y cualquier patrón de nivel 2.

## 5. Predicciones (nivel 2, declaradas antes)

- **Gemelos:** L₁/₈ claramente negativo. El motivo completo sobrevive con probabilidad q².
- **Tripletes:** también negativo. El triplete completo sobrevive con q³ y cada par interno con q², así que el efecto puede ser parecido o mayor que en los gemelos. **No se predice un orden fijo** entre las dos familias (Consejo, sesión 9).
- **Conos identificables, plana y de Sitter:** L ≈ 0.

## 6. Guardián y protocolo

- **Guardián nuevo:** un agente distinto del de R19, que no lee ningún resultado de R14-B, R14-C ni R19, ni los informes del Consejo.
  - Recibe solo §1 a §4, la interfaz `docs/R19_INTERFAZ_SELLADO.md` con la regla de copias y las exclusiones de conos ya vistos, que ahora incluyen también los 5 conos sellados de R19.
  - Entrega el hash del archivo sellado, un **recibo** (versión del código, semillas de selección y log de los requisitos) y nada más.
- **Orden:**
  1. Sellar: el hash se registra en `docs/SELLADO_R20.md`.
  2. Congelar el ejecutor en un commit propio.
  3. Abrir el archivo sellado en otro commit, verificando el hash.
  4. Correr sobre el árbol limpio.
- **Registro:** se conservan los logs de la corrida y las versiones de las dependencias junto a los resultados.

## 7. Convenciones

- **EXPERIMENT_ID 333.**
- **Entropías:**
  - 20261113: estructuras visibles, diluciones, controles y bootstrap;
  - 20261114: selección del guardián (requisitos de visibilidad);
  - 20261115: puntos de las estructuras selladas en la corrida confirmatoria;
  - 20261116: referencias planas (ajuste, calibración y validación).
- **Próximos libres:** EXPERIMENT_ID 334 y entropía 20261117.
- **Costo estimado:** 300 referencias por q (unas 1 a 1.5 h de CPU) más hasta 18 estructuras × 30 realizaciones (unas 6 a 7 h de CPU). Son unas 2 a 2.5 h de pared con 4 núcleos, en tandas de menos de 2 horas.
- **Código:** se extiende `tools/r19.py` con los tres conjuntos de referencias, el radio conforme, la regla de cobertura simultánea, la lectura de los requisitos sellados y el informe de no identificables. Lleva pruebas de la cobertura conforme con datos sintéticos.

## 8. Firma de la versión 2 y orden de ejecución (2026-10-21)

- **Firma del usuario:** «Lo firmo, lanza el guardián, todo lo demás aprobado, corre sin problema».
- **Interfaz del archivo sellado:** es la misma de `docs/R19_INTERFAZ_SELLADO.md`, con un agregado: `IDENTIFICABLES: dict[str, dict]`, con `{"identificable": bool, "medidas": {...}, "semillas": ..., "N": ...}` por estructura del universo, calculado por el guardián con la entropía 20261114 y la banda de §2.
- **Separación de agentes:** el guardián nuevo y el agente del ejecutor trabajan por separado. El ejecutor no ve el archivo sellado. El guardián no lee resultados de R14-B, R14-C ni R19, ni actas o informes del Consejo; solo puede leer el código (incluido `tools/sellado_r19.py`, para no repetir los conos de R19).
