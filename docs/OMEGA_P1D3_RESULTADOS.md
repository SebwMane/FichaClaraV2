# Ω — P1-D.3: resultados de la calibración del instrumento a N grande

- **Fecha:** 2026-10-06.
- **Rama:** `claude/omega-p1d3`.
- **Prerregistro:** `docs/OMEGA_P1D3_PRERREGISTRO.md`.
- **Datos:** `results/p1d3/` (`graphs.jsonl`, `summary.json`, `level3.json`).
- **Panel:** 66 grafos con N ≈ 2·10⁴ y 11 con N ≈ 5·10⁴ (77 en total). Duración: ~20 min con 3 procesos.
- **Reparto:** código de un agente Sonnet, revisado por el cerebro antes de la corrida oficial; auditoría previa de un agente Haiku.

## 1. Veredicto oficial: INSTRUMENTO-VÁLIDO

| Umbral final de Nivel II | Veredicto | Sensibilidad (G-hom, clase correcta) | Especificidad (L ∪ NL) | Nivel III |
|---|---|---|---|---|
| 0.05 (primario) | **VÁLIDO** | 15/15 | 33/33 | Rechaza los retículos de cliques 3D |
| 0.03 | VÁLIDO | 15/15 | 33/33 | ídem |
| 0.08 | VÁLIDO | 15/15 | 33/33 | ídem |

## 2. Resultado por familia (umbral 0.05)

| Clase | Familia | I | II | D_conv | ρ | D_s | Categoría |
|---|---|---|---|---|---|---|---|
| G-hom | RGG2 k12 | 3/3 | 3/3 | 2.02–2.03 | 0.10–0.12 | 2.02–2.04 | GRUESA(2) |
| G-hom | RGG3 k12 / k8 | 6/6 | 6/6 | 3.04–3.11 | 0.21–0.27 | 2.95–3.05 | GRUESA(3) |
| G-hom | RGG en S² | 3/3 | 3/3 | 1.85–1.86 | ≈ 0.10 | 1.96–1.99 | GRUESA(2) |
| G-hom | Anillo / cuadrado / T³ 27³ | CV≡0 | 3/3 | 1.00 / 1.98 / 2.87 | — | 1.00 / 2.00 / 3.01 | GRUESA(1/2/3) |
| G-inh | **RGG3 en caja** | **0/3** | 0/3 | 2.41–2.47 (NO_CONVERGE) | 1.14–1.22 | 2.45–2.57 | **NO_GEOMÉTRICO** |
| G-inh | **RGG2 en caja** | **0/3** | 3/3 | 1.67–1.73 (no entera) | 1.06–1.24 | 1.77–1.83 | **NO_GEOMÉTRICO** |
| G-inh | RGG3 con gradiente de densidad | 3/3 | 3/3 | 2.96–2.98 | 0.65–0.69 | 3.00–3.04 | GRUESA(3); certificado: 1/3 claves pasan |
| R | Caveman K4…K16 (7 tamaños) | 7/7 | 7/7 | 1.00 | ≈ 0.007 | 1.00–1.05 | GRUESA(1) |
| R | Retículos de cliques 3D (K6, K8) | 2/2 | 2/2 | 2.88–2.91 | 0.32 | 3.14–3.21 | GRUESA(3); **certificado F3, F4, F5, F9** → no variedad |
| R | Retículo de cliques 2D (K6) | 1/1 | 1/1 | 1.99 | 0.09 | 2.03 | GRUESA(2) |
| L | Retazos de 2³, 4³ y 8³ bloques | 0/9 | 0/9 | — | 1.06–1.94 | 3.1–5.7 | NO_GEOMÉTRICO |
| L | RGG3 con 0.1 % de atajos | **3/3** | 0/3 (**CRUCE**, D sube hasta 3.7) | — | 0.79–0.86 | 2.93–2.97 | NO_GEOMÉTRICO |
| L | RGG3 con 1 % de atajos | 0/3 | 0/3 | — | 1.15–1.18 | 3.06–3.21 | NO_GEOMÉTRICO |
| L | **WS β = 0.001** | 0/3 | **3/3 (D 2.93–3.28)** | — | 7.4–16.8 | **1.28–1.33** | NO_GEOMÉTRICO |
| L | WS β = 0.002 / 0.005 / 0.01 | 0/9 | 2/9 | — | 2.6–13.5 | 1.4–1.6 | NO_GEOMÉTRICO |
| NL | ER k12 / RR k12 | 3/6 | 0/6 (sin ventana) | — | 0.97 / 1.5–1.6 | 7.7–8.0 | NO_GEOMÉTRICO |

**Subpanel con N ≈ 5·10⁴:**

| Familia | N ≈ 2·10⁴ | N ≈ 5·10⁴ | Lectura |
|---|---|---|---|
| RGG3, D_conv | 3.04–3.05 | 3.06–3.07 | Se mantiene |
| RGG2, D_conv | 2.02 | 2.02–2.03 | Se mantiene |
| T³, D_conv | 2.87 | 2.89 | Se mantiene |
| ρ de RGG3 | 0.21–0.27 | 0.16 | Baja con N |
| ρ de retazos-4³ | 1.22–1.34 | 1.61–1.68 | Sube con N |
| Caveman K8 | 1.00 | 1.00 | Igual |
| WS β = 0.002 | — | II 1/2, D ≈ 4.3–4.4 | Sigue sin ser geométrico |

**La separación entre geometría y retazos crece con N.**

## 3. Contraste con las predicciones congeladas

| # | Predicción | Resultado |
|---|---|---|
| 1 | G-hom pasa I y II con clase correcta (±0.25), S² incluida | ✔ 15/15. Pero «con N = 5·10⁴, RGG3 más cerca de 3» ✘: 3.05 → 3.065, un leve exceso estable |
| 2 | El Nivel I falla en geometrías inhomogéneas (caja 3D, gradiente); el Nivel II converge | **Mitad.** La caja 3D y la caja 2D fallan el Nivel I ✔. **El gradiente pasa el Nivel I** ✘: un gradiente suave de densidad no rompe el autopromediado. El Nivel II **no** converge en la caja 3D ✘: el borde curva D_B2 hacia abajo |
| 3 | R → GRUESA(1/2/3); los análogos 3D fallan el certificado (F9) | ✔ completo. El suavizado a dos pasos eliminó el artefacto de oscilación de P1-D.2 |
| 4 | Ninguna L pasa I y II; los atajos al 0.1 % y WS β = 0.001 como casos difíciles | ✔. Los atajos al 0.1 % pasan I y caen en CRUCE; WS β = 0.001 pasa II y cae en I |
| 5 | NL sin ventana | Parcial: ER sí tiene ventana en el Nivel I (ρ ≈ 0.97), pero no en el Nivel II |
| 6 | D_s ≈ D_conv en G-hom (≤ 0.4); distinto en R | ✔ en G-hom (≤ 0.14). En R: caveman D_s ≈ D_conv ✘; retículos 3D +0.26–0.30 ✔ parcial |

## 4. Hallazgos que importan para Ω

1. **El instrumento es válido dentro de un alcance declarado:** espacios **cerrados y homogéneos o casi homogéneos** (toros, esfera, gradientes suaves de densidad). Ahí distingue geometría gruesa de clase 1, 2 y 3 frente a retazos, atajos, mundo pequeño y aleatorios, y mejora con N.
2. **Falso negativo con borde.** Una geometría auténtica con borde (RGG en caja) sale NO_GEOMÉTRICO: el Nivel I mide invariancia de traslación y el borde la rompe. Si Ω produjera una geometría con borde, este instrumento no la reconocería. Es una **limitación conocida**, no un fallo a corregir *post hoc*.
3. **«D ≈ 3 estable» puede fabricarse con un anillo y 0.1 % de atajos.** WS β = 0.001 pasa el Nivel II con D_conv = 2.93–3.28, clase 3 en una semilla. Es la confirmación más fuerte hasta ahora de `RULE_D3_NEVER_SUFFICIENT`: una dimensión que converge no basta. Lo salvan:
   - el Nivel I (ρ = 7–17);
   - el desacuerdo con D_s ≈ 1.3, observación descriptiva no prerregistrada como criterio.
4. **Geometría gruesa ≠ variedad.**
   - Los retículos de cliques 3D son GRUESA(3) con D_s ≈ 3.2, y el certificado los rechaza (F3, F4, F5, F9). La separación entre las dos categorías funciona en la clase 3.
   - El RGG3 con gradiente pasa el certificado solo con 1/3 claves con N = 729. El certificado es inestable ante una inhomogeneidad suave.
5. **Sesgo de D_conv.** Toma el último valor de la ventana, que está sesgado por saturación y curvatura (S²: el pico es 1.99 y el final 1.85). No cambió ninguna clase, pero conviene saberlo.

## 5. Observaciones *post hoc* (solo hipótesis, no criterio)

- **Concordancia de dimensiones.** |D_s − D_conv| ≤ 0.14 en todas las geometrías homogéneas, 0.26–0.30 en retículos de cliques y 1.6–2.3 en WS. Un tercer nivel de «concordancia de dimensiones» podría cerrar el hueco de WS sin depender del Nivel I. Necesitaría prerregistro y grafos nuevos.
- **Estimador.** Sustituir D_conv = último valor por la mediana de la cola convergida.

## 6. Estado

P1-D.3 queda **cerrada**. El instrumento (observador muestreado, Niveles I y II y certificado como Nivel III, con su alcance) se adopta como **batería obligatoria** para cualquier dinámica futura, según el §7 del prerregistro. No entra en la acción. Rama congelada: `claude/omega-p1d3-congelado`.
