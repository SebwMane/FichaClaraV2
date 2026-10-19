# ES-C0b: movimiento colectivo y rejilla recalibrada. Preregistro de diagnóstico (2026-10-19)

- **Mandato:** firma del usuario del 2026-10-19 («Lo firmo, aprobado, ejecuta») sobre la opción (a) de `OMEGA_ESC0_RESULTADOS.md` §6: movimientos colectivos verificados con balance detallado exacto, primero con b = 0, más la recalibración de la rejilla de b para ε < 1.
- **Tipo:** desarrollo técnico y diagnóstico de mezcla. **No es evidencia** de H_sel.
- **Lo que no cambia de ES-C:** la medida (órdenes con etiquetado natural), el peso exp(−β·S_ε) con √ε corregido, los cuatro inicios, el clasificador de ES-C0 (validado con N = 64 y 128, entropías 20261105 y 20261106) y los criterios de mezcla de ES-C §4.

## 1. El movimiento colectivo R (reinserción de un elemento)

ES-C0 mostró que las cadenas quedan en mesetas de capas, porque mover un elemento de una capa a otra exige muchos cambios de una relación. El movimiento R reubica un elemento entero en un solo paso.

1. Se elige un elemento x al azar.
2. Se quitan todas sus relaciones. Queda C′, que sigue siendo un orden, sobre los demás elementos.
3. Se arma la familia F(x, C′) de pares (D, U):
   - D es el pasado nuevo de x: vacío, o el pasado principal de un elemento y con etiqueta menor que x (y junto con todo lo que precede a y en C′);
   - U es el futuro nuevo de x: vacío, o el futuro principal de un elemento z con etiqueta mayor que x;
   - el par es **válido** si todo elemento de D precede a todo elemento de U en C′ (así el resultado es transitivo).
4. Si el pasado y el futuro actuales de x no son un par de F(x, C′), el movimiento **se rechaza** sin cambiar nada.
5. Si lo son, se propone un par de F(x, C′) al azar y uniforme, y se acepta con min(1, exp(−β·ΔS_ε)).

**Por qué cumple balance detallado:** C′ es el mismo en la ida y en la vuelta, así que la familia y la probabilidad de proponer cada par (1/|F|) son iguales en los dos sentidos. La propuesta es simétrica.

- **Núcleo mixto:** cada propuesta es R con probabilidad 0.1 y el movimiento local de ES-C0 con probabilidad 0.9. Un barrido sigue siendo N(N−1)/2 propuestas.
- **Verificación antes de correr:**
  - enumeración exacta con N ≤ 5 de balance detallado, estacionariedad e irreducibilidad del núcleo mixto, con varios ε y β (incluido β negativo);
  - prueba de variación total en una cadena larga;
  - ΔS_ε del movimiento R igual a la acción recalculada;
  - el orden sigue siendo transitivo y con etiquetado natural después de cada R aceptado.

## 2. Etapa 1: medida de conteo (b = 0)

- **Celdas:** ε ∈ {1, 0.25} (con b = 0 el ε no cambia la medida, pero sí el observable S_ε) y N ∈ {64, 128}.
- **Cadenas:** las mismas 4 de ES-C0 (orden vacío, cadena total, sprinkling 4D y KR natural), con 1000 barridos para N = 64 y 500 para N = 128.
- **Pregunta:** ¿el núcleo mixto junta las mesetas que ES-C0 dejó separadas con b = 0? Se compara con ES-C0, que tenía el mismo largo y los mismos inicios.

## 3. Etapa 2: rejilla recalibrada (solo si la etapa 1 mezcla con algún N)

- **Unidad de β:** σ₀(ε, N) es la desviación estándar de S_ε en la mitad final de las 4 cadenas de la etapa 1 de la misma celda (la regla de ES-0). Corrige el error de ES-C0 §4.
- **Rejilla:** β = κ/σ₀, con κ ∈ {−4, −1, 1, 4}.
- **Mismas cadenas y largos** que la etapa 1, solo en los N donde la etapa 1 mezcló.

## 4. Qué se informa

- Mezcla por celda (R̂ dividido < 1.05 y ESS ≥ 200 en S_ε, en la fracción de relaciones y en la altura).
- Tasa de aceptación de R y fracción de propuestas R aplicables (pasado y futuro principales).
- Comparación directa con ES-C0 para b = 0.
- La bandera 4D (fracción por cadena y barrido de salida), solo descriptiva.

## 5. Convenciones

- **EXPERIMENT_ID 330.** Entropía de las cadenas: 20261107.
- **Clasificador:** reutiliza el de ES-C0 (sus archivos oficiales de `results/omega_esc0/`).
- **Próximos libres:** EXPERIMENT_ID 331 y entropía 20261108.
- **Herramienta:** agente de código con auditoría; árbol limpio.
