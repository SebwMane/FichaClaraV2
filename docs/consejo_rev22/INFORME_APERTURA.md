# Consejo Ω, sesión 11: informe de apertura sobre R20

- **Fecha:** 2026-10-21.
- **Pedido del usuario:** «prepara el informe y su apertura para el consejo para revisar los resultados y así decidimos el siguiente paso».
- **Mandato del Consejo:**
  1. auditar el PASA de R20 y su alcance;
  2. recomendar el siguiente paso del programa.

  Las dos decisiones son del usuario.
- **Fuentes en el repositorio:**
  - **Informe de resultados:** `docs/OMEGA_R20_RESULTADOS.md` (commit 7b28bf6).
  - **Preregistro:** `docs/OMEGA_R20_REPETICION_PREREGISTRO.md`, versión 2, firmada (commits 0983987 y 1000201).
  - **Sellado y recibo del guardián:** `docs/SELLADO_R20.md` y `docs/sellado_r20/` (RECIBO.txt, universo v1 descartado, identificables v1).
  - **Datos y logs:** `results/omega_r20/` (summary.json, tabla_L.csv, tabla_bandas.csv, tabla_familias.csv, entorno.json, logs/).
  - **Figura:** `docs/figuras/r20.png`.
  - **Código:** `tools/r20.py` y `tests/test_r20.py`.
  - **Antecedentes:**
    - Rev. 20 (`OMEGA_CONSEJO_REV20_R19.md`) y Rev. 21 (`OMEGA_CONSEJO_REV21_R20.md`);
    - R19-(b) y su errata (`OMEGA_R19_GEMELOS_RESULTADOS.md` §5).
- **Cadena de custodia:**

  | Paso | Commit |
  |---|---|
  | Firma de la versión 2 | 1000201 |
  | Hash registrado | c338b9c |
  | Ejecutor congelado | 3cda6bd |
  | Archivo sellado abierto (hash verificado) | c2e227d |
  | Corrida, con un solo commit y el árbol limpio | c2e227d |
  | Resultados | 7b28bf6 |

## 1. Veredicto oficial: **PASA**

| Condición (§4, versión 2) | Resultado |
|---|---|
| Cobertura de validación ≥ 87/100 en cada q (banda conforme, radio = 96.º de 100) | 92, 92, 100 y 95 de 100 |
| Calibración nula: IC de L₁/₈ dentro de ±0.262 | plana 0.10 (−0.08 a 0.17); de Sitter 0.03 (−0.05 a 0.21) |
| Población: ≥ 2 familias por tipo, con ≥ 2 identificables | gemelos 3/4, tripletes 3/4, conos poliédricos 4/4, conos suaves 4/4 |
| «Se desvanece»: IC entero por debajo de −0.262 | las 6 votantes |
| «Autosemejante»: IC dentro de ±0.262 | las 8 votantes |

PASA se decide por intersección-unión: todas las 14 votantes cumplen.

## 2. Datos

### 2.1 L_q = ln(D_obs/D_ctrl), mediana sobre 30 realizaciones (q = 1, 1/2, 1/4, 1/8)

| Estructura | Estado | q = 1 | 1/2 | 1/4 | 1/8 | D_obs / D_ctrl con q = 1/8 |
|---|---|---|---|---|---|---|
| G1b, gemelos al 12% | vota | −0.01 | −0.65 | −0.93 | −0.80 | 0.79 / 1.71 |
| G1c, gemelos al 17% | vota | 0.01 | −0.77 | −1.27 | −1.01 | 0.87 / 2.38 |
| G1d, gemelos al 23% | vota | 0.00 | −0.86 | −1.39 | −1.41 | 0.84 / 3.19 |
| G2b, tripletes al 12% | vota | 0.02 | −0.50 | −0.87 | −0.91 | 0.79 / 1.92 |
| G2c, tripletes al 17% | vota | 0.01 | −0.58 | −0.85 | −1.11 | 0.85 / 2.51 |
| G2d, tripletes al 23% | vota | −0.01 | −0.59 | −1.03 | −1.32 | 0.91 / 3.40 |
| P1a a P1d, conos poliédricos | votan | ±0.01 | ±0.01 | de −0.04 a 0.03 | de −0.02 a 0.01 | de 1.32 a 2.82 (los dos brazos) |
| S2a a S2d, conos suaves | votan | ±0.01 | ±0.01 | de 0.00 a 0.04 | de −0.07 a 0.08 | de 1.25 a 2.18 (los dos brazos) |
| G1a, gemelos al 7% | no identificable; **no evaluable** con q = 1/8 | 0.02 | −0.57 | −0.70 | sin dato | r₂₀ no medible en el control de 1 realización |
| G2a, tripletes al 7% | no identificable | 0.01 | −0.43 | −0.71 | −0.54 | 0.74 / 1.29 |
| D5, D7 y D8 (vistas) | no votan | | | | −0.63, −0.04 y 0.06 | |
| plana y de Sitter (calibración) | | | | | 0.10 y 0.03 | |

### 2.2 Lecturas de nivel 2 que el informe de resultados no destaca

1. **En las copias, el observado diluido queda dentro de la banda plana con q = 1/8:** D_obs va de 0.79 a 0.91, mientras el control va de 1.7 a 3.4. El defecto no solo se aleja de su control: **se vuelve indistinguible del plano**, que es lo que H_E esperaría de un defecto microscópico.
2. **Efecto techo en L.** D_obs choca con el piso de ruido de la banda, cerca de 0.8, así que la magnitud de L con q = 1/8 está dominada por cuán visible sigue el **control**, no por cuánto decayó el observado. Por eso en G1b y G1c, |L| es mayor con q = 1/4 que con q = 1/8. **L mide bien el signo y la separación, pero no la velocidad de decaimiento.**
3. **Gemelos frente a tripletes con q = 1/2:** a igual fracción, los gemelos se separan más rápido (−0.65, −0.77 y −0.86) que los tripletes (−0.50, −0.58 y −0.59). Es compatible con que el motivo completo sobreviva con q² para los gemelos y q³ para los tripletes, pero cada par dentro de un triplete sobrevive con q², así que un triplete diluido conserva más pares gemelos. No había predicción de orden (sesión 9); es descriptivo.

## 3. Límites que el Consejo debe auditar

1. **Selección a nivel de familia.** La primera familia de conos suaves tuvo 0 de 4 estructuras visibles. El guardián la reemplazó **una vez** (registrado). La conclusión sobre conos suaves vale solo para los visibles, y los descartados quedan como «no identificables por este instrumento».
2. **Clasificación de la familia S2.** Incluye prismas de sección curva (Reuleaux, semicircular y lenticular), con caras planas en una dirección. Los únicos conos «suaves» en todas las direcciones son de tipo circular.
3. **Distancia de selección distinta de la confirmatoria.** El guardián usó los componentes finitos por par; el ejecutor, los finitos en las 300 referencias. Declarado antes de correr; en la práctica coinciden (r₄₀ queda fuera con q = 1/8).
4. **G1a no evaluable.** Le faltó un observable en una realización del control. No votaba, así que no afecta el veredicto, pero muestra que la regla «sin exclusiones selectivas» sí se aplicó.
5. **Inferencia condicional al instrumento ajustado.** Las referencias no se remuestrearon. FALLA no tenía control del error familiar, y no ocurrió.
6. **El margen ±ln 1.3** se fijó en R19 a partir del ruido del borde de la banda. Con la banda conforme, los IC de los conos son mucho más estrechos (de ±0.05 a ±0.15). El margen fue holgado.
7. **Alcance.** Son 16 estructuras sintéticas de 4 mecanismos, con N₀ = 8000. **No dice nada sobre Ω, H_E, H_L ni H_sel**: es una validación del instrumento Φ_q más el control reinyectado.

## 4. Estado del programa

| Línea | Estado |
|---|---|
| H_sel (selección dinámica) | indeterminada; mezcla cerrada con N de 64 a 128 (Rev. 18) |
| Ω-E (entropía SJ) | herramienta descriptiva en 2D; no separa KR por el coeficiente (A-E1′) |
| Ω-C (acción BD) | √ε corregido; no mezcla con N ≥ 64 (ES-C0) |
| Reconocimiento (R14-C) | percolación, producto y poliedros se separan; conos suaves débiles y confundidos con de Sitter |
| **Instrumento de grano grueso** | **validado (R20, PASA)** en estructuras sintéticas visibles |

## 5. Opciones para el siguiente paso (recomendación del orquestador al final)

| Opción | Qué haría | A favor | En contra |
|---|---|---|---|
| **(a)** Aplicar el instrumento validado a órdenes no geométricos ya conocidos (KR, percolación, producto, órdenes en capas de ES-C0) | Preguntar si alguno se vuelve indistinguible de una geometría bajo Φ_q **más rápido que su control** (señal H_E) o si es autosemejante | Usa lo que R20 validó; responde a la prioridad firmada H_E; no necesita muestrear la dinámica | Para órdenes sin «defecto», el control es la misma familia generada con N chico. Si el orden no es autosemejante (por ejemplo, la percolación con p fijo), L ≠ 0 mide eso, no H_E. Hay que definir con cuidado qué cuenta como señal. |
| **(b)** Ampliar el instrumento a conos invisibles (opción B de R20-2) | Desarrollar y calibrar otro estadístico (por ejemplo, multivariable con λ_max del Laplaciano y sub-intervalos grandes, contra la familia curva) | Cierra el hueco de los conos «redondos casi en todas partes» | R14-C mostró que de Sitter imita esas firmas; requiere un piloto aparte antes de cualquier confirmatoria |
| **(c)** Medir la velocidad de decaimiento con un estadístico sin techo | Reemplazar L por un estadístico que no sature con el ruido (por ejemplo, la pendiente de D_obs/D_ctrl en q mientras D_obs > 1) | Aprovecharía la diferencia q² frente a q³ | Es un refinamiento del instrumento, no una pregunta nueva sobre Ω |
| **(d)** Síntesis del programa | Un documento de estado que reúna lo establecido, lo indeterminado y lo cerrado, con la escala de afirmación | Ordena el proyecto antes de abrir otra línea; útil para revisión externa | No produce datos nuevos |

**Recomendación del orquestador:** **(d) primero, después (a) como exploratoria y sin voto.**
- La línea de grano grueso llegó a su primer resultado confirmatorio, y conviene fijar su lugar en el programa antes de usarla.
- Luego, (a) con un preregistro que defina antes qué es señal H_E para órdenes sin defecto. Mi propuesta: comparar el decaimiento de D(q) de la familia frente a su propia referencia generada con N chico.
- (b) y (c) quedan como refinamientos del instrumento, a pedido.

## 6. Preguntas para el Consejo

| # | Pregunta |
|---|---|
| R22-1 | ¿Se confirma el PASA de R20 tal como está, con los límites de §3? ¿Hace falta alguna salvedad adicional en el informe de resultados? |
| R22-2 | ¿El reemplazo de la familia de conos suaves por el guardián (§3.1) es aceptable en este alcance, o debe constar como desviación del espíritu de §3? |
| R22-3 | ¿Qué siguiente paso: (a), (b), (c), (d) o una combinación? |
| R22-4 | Si es (a): ¿qué definición previa de «señal H_E» para órdenes sin defecto? |
