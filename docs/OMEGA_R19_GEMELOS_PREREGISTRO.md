# R19-(b): convergencia de defectos bajo dilución, con control de defecto reinyectado. Preregistro confirmatorio (2026-10-20)

- **Mandato:** firma del usuario del 2026-10-20: «Apruebo R19-1 a R19-3, prepara el preregistro de (b)». Fuente: el informe de apertura `docs/consejo_rev19/INFORME_APERTURA.md`.
- **R19-3 aplicada:** las familias se arman por mecanismo y se informa siempre por estructura además de por familia.
- **Estado:** **PREPARADO, no ejecutado.** Para ejecutarlo hacen falta:
  - la firma de R19-4 (familias nuevas y guardián, §6);
  - el sellado del guardián;
  - una orden explícita del usuario.
- **Lo que esta prueba puede y no puede decir:** prueba el **mapa de grano grueso** (la dilución Φ_q) como instrumento para separar H_E de H_L. **No prueba que la dinámica de Ω produzca esos órdenes** (eso es H_sel, indeterminada; Rev. 18). Un PASA dice que «Φ_q distingue un defecto que se desvanece de uno que solo deja de verse». No dice que «la naturaleza tenga gemelos».

## 1. Hipótesis

R14-B mostró que diluir confunde convergencia con pérdida de resolución: al bajar N, las bandas se ensanchan y cualquier defecto fijo termina pareciendo en banda.

- **H_conv (lo que se pone a prueba):** existe un control que separa las dos cosas.
  - Un defecto que **se desvanece** bajo Φ_q, como los gemelos (un par sobrevive con probabilidad q²), se aleja de su control a densidad fija a medida que q baja.
  - Un defecto que **no se desvanece** (un cono autosemejante) queda igual a su control.

## 2. El control

Para cada estructura X, cada realización y cada q se construyen dos órdenes con el mismo N efectivo:

| Tipo de defecto | Observado: Φ_q(X) | Control: defecto a densidad fija |
|---|---|---|
| **Gemelos y copias** (defecto agregado a una base de sprinkling) | se diluye el orden completo | se diluye **solo la base** y después se **reinyecta** el defecto con la misma fracción y el mismo procedimiento que al generar X |
| **Conos y geometrías** (sprinkling de una región con su propia relación causal) | se diluye el sprinkling de N₀ | un sprinkling **nuevo** de la misma región con N ~ Binomial(N₀, q) elementos |

- **El caso de los conos y las geometrías:** diluir un sprinkling da exactamente otro sprinkling de densidad q·ρ, así que observado y control tienen la misma distribución. Por eso esa fila es la **calibración nula** del método.
- **Estadístico por realización:** Lq = ln[D_obs(q)/D_ctrl(q)].
  - D es la distancia de Mahalanobis a la referencia plana de ese N efectivo, normalizada por el radio de la banda, como en R14-B v2.
  - El vector es el mismo de R14-B v2: 14 componentes, incluido el detector de gemelos.
- **Estadístico principal:** L = L_{1/8}. Se informan también L_{1/2} y L_{1/4}.

## 3. Estructuras

| Grupo | Estructuras | Papel |
|---|---|---|
| **Calibración nula, vistas** | plana; de Sitter; D7 y D8 (conos de gauge F) | se espera L ≈ 0 |
| **Defecto que se desvanece, visto** | D5 (800 gemelos sobre 7200 elementos de base) | caso de desarrollo; **no vota** |
| **Familias selladas nuevas** (las define el guardián, §6) | familias de «se desvanece» y de «autosemejante», con 2 o 3 estructuras cada una | **votan** |

- **Unidad:** la realización, con 20 realizaciones por estructura.
- **Escala:** N₀ = 8000 y q ∈ {1, 1/2, 1/4, 1/8}.

## 4. Reglas de decisión (congeladas)

- **IC 95% de la mediana de L sobre las 20 realizaciones,** por bootstrap con 2000 réplicas.
- **Margen de equivalencia:** |L| ≤ ln 1.3. La banda de R14-B, calibrada con 20 referencias, tiene un ruido de borde del orden del 20% al 30%, así que diferencias menores no son distinguibles.
- **Calibración** (debe cumplirse; si no, el resultado es **INDETERMINADO**):
  - plana y de Sitter: el IC de L queda dentro de [−ln 1.3, ln 1.3];
  - la referencia plana cubre al menos 16 de 20 sprinklings nuevos de validación con cada N efectivo.
- **PASA** si se cumplen las dos condiciones:
  1. en **cada** estructura sellada de «se desvanece», el IC de L queda entero **por debajo de −ln 1.3**;
  2. en **cada** estructura sellada «autosemejante», el IC de L queda entero **dentro de [−ln 1.3, ln 1.3]**.
- **FALLA** si alguna estructura sellada queda del lado contrario: un «se desvanece» con IC por encima de −ln 1.3, o un «autosemejante» con IC entero fuera del margen. Cualquier otro caso es **INDETERMINADO**.
- **Informe:** por estructura y por familia (R19-3), con L_q en los cuatro q y las trayectorias de D_obs y D_ctrl.

## 5. Predicciones (nivel 2), declaradas antes

- **Gemelos al 10%:** con q = 1/8, el observado conserva pocos pares completos (cerca de q² de los 800 originales). El control conserva el 10%. Se espera L de −1 a −2 (en R14-B, D5 diluido a 1/8 queda en 0.74, y un control con 10% de gemelos debería seguir fuera de banda).
- **Conos D7 y D8, plana y de Sitter:** L cerca de 0.
- **Riesgo declarado:** si el control reinyectado ya queda en banda con q = 1/8 (por ejemplo, con fracciones de defecto chicas), L se acerca a 0 aunque el defecto sí se desvanezca. Por eso las densidades selladas deben dejar el control **fuera de banda en todo q**, algo que el guardián verifica con su generador y no con estos datos (§6).

## 6. Familias nuevas y guardián (R19-4, **firma pendiente**)

**Propuesta del orquestador:**
- **Guardián:** un agente aparte que no lee R14-B, R14-C ni este informe de resultados.
  - Recibe solo §2, §3 y §4 y la lista de mecanismos permitidos.
  - Escribe las familias en un archivo sellado y entrega solo su hash.
  - El archivo se abre después de congelar el código y antes de correr.
- **Mecanismos permitidos:**
  - **Se desvanece:** copias de elementos con pasado y futuro idénticos (gemelos, tripletes, etc.), con una fracción de entre 5% y 25%, sobre una base de sprinkling 4D plano.
  - **Autosemejante:** conos de gauge F distintos de D6, D7, D8 y de los vistos en F6 y F7, con su caja de sprinkling.
- **Condición que el guardián verifica con sus propios generadores:** con N = 1000 (N efectivo de q = 1/8), el control de cada familia «se desvanece» queda fuera de la banda plana. Lo hace con referencias propias, sin datos de este proyecto.

## 7. Convenciones

- **EXPERIMENT_ID 332.**
- **Entropías:**
  - 20261110: estructuras y diluciones;
  - 20261111: generador sellado del guardián;
  - 20261112: referencias y validación.
- **Próximos libres:** EXPERIMENT_ID 333 y entropía 20261113.
- **Código:** reutiliza `tools/r14b.py` (observables, Mahalanobis con Ledoit-Wolf, banda) y agrega el control, el estadístico L y las reglas de §4.
- **Antes de correr:**
  - pruebas que verifiquen que el control conserva la fracción de defecto en cada q;
  - pruebas que verifiquen que el control de un sprinkling tiene la distribución de Φ_q.
- **Costo estimado:** de 1 a 2 horas, con árbol limpio y un solo commit.
