# E1′: preregistro congelado (paso 4, 2026-10-12)

- **Fuente de decisiones:** `OMEGA_CONSEJO_REV4_E1PRIMA.md` (§2, §6, §6.1b, §8, §9). Este documento las reúne y las cierra; desde su commit no se ajusta nada.
- **Ejecución:** la compuerta V′ (paso 5) puede correr tras la auditoría; **E1′ (paso 6) requiere la firma R4-1 del usuario.**

## 1. Instrumento

SSEE con vacío de Sorkin-Johnston (`omega/causality/entropy.py`, validado en E0a y E0b′): iΔ = (i/2)(C − Cᵀ); W = parte positiva; S = Σ λ ln|λ|. Doble truncamiento con c · √N_L/(4π) en la región grande y c · √n_R/(4π) en la restringida. **Una sola diagonalización de iΔ por corrida**, reutilizada para todos los umbrales y para el descriptor espectral.

## 2. Decisión de congelamiento: tamaño de región N/8 (declarada)

Al consolidar se detectó que con el objetivo N/4 ± 5% (Rev. 4 §2.1) **KR no tiene intervalos candidatos**: sus intervalos máximos tienen ≈ N/8 + 2 elementos (en E1, 371 de 2400). KR quedaría «no aplicable» y E1′ no podría responder la pregunta que lo motivó. **Se congela el objetivo N/8 ± 5% para todas las familias.**

En la geometría, una región de N/8 centrada y cuadrada es un subdiamante de lado s = √(1/8) ≈ 0.354, con |pasado(a)|/N = (1/2 − s/2)² ≈ 0.104. El coeficiente esperado sigue siendo ≈ 1/3; esto no fue validado en E0 (que usó N/4) y **lo verifica V′** antes de mirar ninguna otra familia.

## 3. Regla de región

Entre los intervalos I[a, b] (con extremos) con |I| ∈ [⌊0.95 · N/8⌋, ⌈1.05 · N/8⌉], el que maximiza min(|pasado(a)|, |futuro(b)|); empate por menores índices. Sin candidatos: la corrida no tiene región.

## 4. Criterio que vota (umbral c = 1)

Ajustes sobre las corridas de una familia (n_R frente a S truncada):
1. **Forma logarítmica:** R² de S = a ln(√n/4π) + b mayor que R² de S = a′n + b′.
2. **Coeficiente:** |a − 1/3| ≤ 2 · SE(a) + 0.05.

## 5. Diagnósticos sin voto

- **Complementariedad:** |S(R) − S(R^⊥)| / S(R) con c = 1; base de la geometría ≈ 30% (piloto 1).
- **Barrido del umbral:** criterio de §4 evaluado también con c = 0.5 y c = 2. Lectura: «dependiente del umbral» si una familia no geométrica pasa con c = 1 y no con c = 0.5 ni c = 2 mientras la geometría pasa con los tres; «robusto» si ambas pasan con los tres.
- **Descriptor espectral:** rodilla R-rodilla-2 y pendiente β del espectro de la región grande (rechazada para truncar, §8.4).

## 6. Compuerta V′ (paso 5): solo sprinklings 2D

- N_L ∈ {800, 1200, 1600, 2000, 2400}, 10 semillas, entropía 20261052.
- **Forma:** |pasado(a)|/N y |futuro(b)|/N en [0.085, 0.125] en ≥ 90% de las corridas.
- **Reproducción:** las dos condiciones de §4 con c = 1.
- Si falla, E1′ no se ejecuta y vuelve al Consejo.

## 7. E1′ (paso 6, requiere R4-1)

- Familias: geometría (control propio), KR, percolación p = 0.05, capas gruesas (8), bicapa, sándwich. Misma malla, 10 semillas, entropía 20261053.
- Familia sin región en ninguna semilla de N_L = 2400: «no aplicable: sin regiones extensas» (sin ley de área).
- **DISCRIMINA:** la geometría cumple §4 y ninguna familia no geométrica lo cumple.
- **NO DISCRIMINA:** alguna familia no geométrica cumple §4.
- **INDETERMINADO:** la geometría no cumple §4 en E1′.
- Riesgo aceptado por el usuario (Rev. 4 §6.1b): el criterio §4.2 puede admitir a KR.

## 8. Escala de afirmación

DISCRIMINA o NO DISCRIMINA son nivel 1 para 2D, N ≤ 2400 y estas familias. No prueban ni refutan la ruta de Jacobson completa.

## 9. Aprobación del usuario (2026-10-12)

El usuario aprobó la decisión de congelamiento de §2 (tamaño de región N/8). No se modifica ningún criterio. R4-1 sigue pendiente.
