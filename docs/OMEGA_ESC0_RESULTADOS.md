# ES-C0: diagnóstico de mezcla de Ω-C euclídea (2026-10-19)

- **Preregistro:** `OMEGA_ESC_PREREGISTRO.md` (commit 4d8d45e, antes del código), con la aclaración A1.
- **Código:** commit 33aa25f, árbol limpio, 120 cadenas oficiales, 869 s de pared.
- **Datos:** `results/omega_esc0/`. **Figura:** `docs/figuras/esc0_piloto.png`.
- **Tipo:** diagnóstico. **No es evidencia** de H_sel.

## 1. Resultado en una frase

**Solo mezcla N = 32.** Con N = 64 y 128 ninguna celda mezcla. Esto incluye b = 0, que es la medida de conteo pura, sin peso de acción. El obstáculo está en el espacio de órdenes y en el movimiento, no en la acción.

## 2. Mezcla por celda

| ε | N | b = −12 | b = −4 | b = 0 | b = +4 | b = +12 |
|---|---|---|---|---|---|---|
| 1 | 32 | no | no | **sí** | no | no |
| 1 | 64 | no | no | no (R̂ de relaciones 7.3) | no | no |
| 1 | 128 | no | no | no | no | no |
| 0.25 | 32 | **sí** | **sí** | **sí** | **sí** | casi (ESS 117) |
| 0.25 | 64 | no | no | no | no | no |
| 0.25 | 128 | no | no | no | no | no |

## 3. Lo que muestran las cadenas (figura)

- **Mesetas que no se cruzan (N = 64 y 128).** Cada cadena se queda en la meseta a la que llega primero, y ninguna pasa a otra en 500 a 1000 barridos:
  - la que parte de KR queda en órdenes de **3 capas**, con fracción de relaciones 0.38 y altura 3;
  - la del orden vacío y la del sprinkling 4D quedan en órdenes de **2 o 3 capas más ralas** (0.22 a 0.31);
  - con ε = 1 y |b| grande, la cadena total queda casi congelada (0.97 a 0.99).
- **Lectura (nivel 2):** cambiar el número de capas de un orden en capas exige muchos movimientos coordinados. Quitar un link o agregar una relación de a uno no basta en estos tiempos. Es la misma relajación lenta que vio ES-0 con λ = 0 en N = 64, que usaba la misma medida base.
- **La geometría 4D se pierde de inmediato.** La cadena que parte de un sprinkling 4D deja de estar marcada en el **primer barrido**, en todas las celdas. Con N = 64 y 128 la fracción marcada es 0 en la mitad final de todas las cadenas, igual que en P4. Como no hay mezcla, esto **no** mide la probabilidad de fase: solo dice que la geometría no es un estado de partida estable bajo esta dinámica.
- **Con N = 32 la bandera es INDETERMINADA** (aclaración A1); no se interpreta.

## 4. Un problema de diseño que es mío: la escala de b con ε = 0.25

Con ε = 0.25 la acción es mucho más chica: en la mitad final, S_ε va de unos −30 a −56 con N = 64, frente a cientos o miles con ε = 1. Con β = b/N, el exponente β·S_ε queda de orden 1 o menos aun con |b| = 12. **La rejilla de b de ε = 0.25 es casi la medida de conteo** y no explora el peso. La rejilla debió calibrarse con la desviación de S_ε, como hizo ES-0 con σ₀. Queda registrado; no se corrige sobre estos datos.

## 5. Costo

0.006, 0.027 y 0.105 s por barrido con N = 32, 64 y 128. El costo no es el obstáculo: con N = 128, 10 000 barridos son unos 17 minutos por cadena. Pero las mesetas no muestran tendencia a juntarse, así que más barridos solos probablemente no bastan.

## 6. Qué implica (para decidir; nada se adopta sin firma)

1. **ES-C1 no es viable tal como está**, salvo con N = 32, donde el clasificador no valida.
2. **El problema de mezcla es de la medida base de órdenes con etiquetado natural y del movimiento local**, y afecta por igual a Ω-E y a Ω-C. Opciones:
   - (a) movimientos colectivos que muevan elementos entre capas, verificados con balance detallado exacto;
   - (b) templado paralelo en b, que no ayuda si la barrera está también en b = 0;
   - (c) Wang-Landau o muestreo multicanónico en la altura o en el número de capas, que es la coordenada lenta que muestra la figura;
   - (d) aceptar N = 32 con un clasificador recalibrado para N chico (como L3cE-A1), sabiendo que N = 32 es muy pequeño para 4D.
3. **La rejilla de b de ε < 1** debe calibrarse con la escala de S_ε antes de cualquier repetición.
