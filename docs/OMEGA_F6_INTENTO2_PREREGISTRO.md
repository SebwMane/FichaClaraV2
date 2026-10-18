# F6, intento 2 (último): recomendaciones del Consejo y preregistro congelado (2026-10-18)

- **Pedido del usuario:** «Que el Consejo arme un grupo de recomendaciones para las preguntas y ejecuta». El usuario delegó R13-1 a R13-4 en el Consejo y autorizó la ejecución. Marca de compuerta: `R13-1 | **Firmada**`.
- **Base:** `OMEGA_CONSEJO_REV12_F6_INTENTO2.md`.
- **Este documento se congela antes de abrir el conjunto sellado nuevo** (`CONFIRMACION_CIEGA_2.md`, definido por un guardián nuevo). Nada de §2 a §6 puede cambiar después de abrirlo.

| # | Firma |
|---|---|
| R13-1 | **Firmada** (por delegación del usuario en el Consejo; ver la tabla) |

## 1. Recomendaciones del Consejo (decisiones delegadas)

| # | Pregunta | Recomendación y decisión | Argumento | Objeción registrada |
|---|---|---|---|---|
| R13-1 | ¿Prueba de colapso como segunda condición con voto? | **Sí. Firmada.** | La condición de colapso estaba en H_loc desde Rev. 11 (§2, condición 3) y no necesita modelo: cada par compara la misma escala relativa. El intercepto atrapa conos y el colapso atrapa defectos que se borran con el tamaño; se complementan. | (Sobreajuste) La decisión de que vote nace de ver C5. Por eso C0 a C7 pasan a «vistas» y hay conjunto sellado nuevo. Una condición más también aumenta el riesgo de rechazar positivas: se acepta como el precio de no confundir defectos con curvatura. |
| R13-2 | ¿30 o 20 realizaciones por N? | **30. Firmada.** | Con 20, C5 quedaría cerca de 47, con poco margen sobre un umbral que también depende de los grupos de calibración. Con 30 el error baja a la mitad y la prueba tiene potencia real; el costo (unas 7 horas) es aceptable. | (Escéptico) Más realizaciones hacen visible cualquier sesgo pequeño de una geometría; si el riesgo de FRW es real, se verá. Se acepta: es lo que el intento debe averiguar. |
| R13-3 | ¿Guardián nuevo con conjunto sellado nuevo? | **Sí. Firmada** (ya convocado). | Sin estructuras nuevas a ciegas, un PASA no valdría nada: el criterio nuevo se diseñó viendo C5. | Ninguna. |
| R13-4 | ¿Último intento? | **Sí. Firmada.** | Es el presupuesto de Rev. 11 §6. Si falla, la línea del perfil local vuelve al Consejo para una revisión completa, sin probar otra variante. | Ninguna. |

**Cambio técnico del Consejo (no vota distinto, se declara):** la referencia de r pasa a ser el valor exacto **r_plano = 1/10** para todo m (dos veces la fracción de Myrheim-Meyer f(4) = 1/20; nivel 3). Motivo: la referencia empírica del intento 1 tenía un error compartido de −0.002 en m = 20 (Rev. 12 §2). **¿Puede crear un positivo?** Solo elimina un ruido común; afecta por igual a la calibración y a la prueba.

## 2. Instrumento

- Igual que el intento 1 (`tools/f6_confirmacion.py`, que importa de `tools/f6_exploracion.py`), salvo la referencia: **Δr = media de r − 1/10**.
- m ∈ {10, 20, 40, 80, 160}, regla m ≤ 0.1·m_max, entre 30 y 300 intervalos por celda.
- Una **estructura medida** = **30 realizaciones con N = 4000 y 30 con N = 8000**. Una celda (m, N) entra si al menos la mitad de sus realizaciones son medibles; su valor es la media de Δr sobre las medibles y su error se = desviación estándar entre realizaciones / √(medibles).

## 3. Clasificación

Una estructura es **«espacio-tiempo 4D»** si cumple las dos condiciones:

1. **Intercepto:** ajuste Δr = b + A·√(m/N) por mínimos cuadrados sobre las celdas y |b − b_plano| ≤ τ_b.
2. **Colapso:** K = Σⱼ (dⱼ/seⱼ)² ≤ τ_K, con dⱼ = Δr(2m, 8000) − Δr(m, 4000) para m ∈ {10, 20, 40, 80} y seⱼ = √(se₁² + se₂²). Solo cuentan los pares con las dos celdas presentes.

**No medible** (no es espacio-tiempo): menos de dos valores distintos de x en el ajuste, o menos de dos pares en K.

## 4. Calibración y umbrales (entropía 20261084)

- **15 grupos** (10 diamantes planos, 3 de Sitter con η ∈ [−1.5, −0.5] y 2 FRW polvo), cada uno una estructura medida completa (30 + 30 realizaciones).
- **b_plano** = mediana de b en los 10 grupos planos; **τ_b** = 1.5 × max sobre los 15 grupos de |b − b_plano|.
- **τ_K** = 1.5 × max sobre los 15 grupos de K.
- Si algún grupo no es medible, el veredicto es FALLA.

## 5. Estructuras y veredicto

- **Selladas nuevas** (`CONFIRMACION_CIEGA_2.md`, guardián nuevo; entropías 20261082 para la generación y 20261083 para el azar auxiliar, spawn_key (322, k, N, s)), con el papel que les asignó el guardián.
- **Vistas con semillas nuevas** (entropía 20261085):

| Grupo | Positivas | Negativas | Zona gris (sin voto) |
|---|---|---|---|
| Geometrías y conos ya vistos | diamante, de Sitter, de Sitter fuerte, anti de Sitter, FRW radiación, FRW polvo | cono cúbico, cono octaédrico, cono p = 6, orden producto 4, percolación 0.05 | conos p = 1.5 y p = 3 |
| Selladas del intento 1, ahora vistas | C0 Kasner, C1 cono elíptico, C2 densidad ondulada | C3, C4, **C5 (links borrados)**, C6 | C7 |

- **PASA** si todas las positivas son espacio-tiempo y ninguna negativa con voto lo es. **FALLA** en cualquier otro caso.
- **Informado sin voto:** pendiente A y su signo, rama 2 (D_loc y max-T frente a una nula de 30 diamantes planos por N, entropía 20261084), y los márgenes |b − b_plano|/τ_b y K/τ_K por estructura.

## 6. Convenciones

- **EXPERIMENT_ID 322.**
- **Entropías:**
  - 20261082 y 20261083: selladas;
  - 20261084: calibración, nula y bootstrap;
  - 20261085: vistas;
  - 20261086: referencia del perfil local (solo para la rama 2, informativa).
- **spawn_key:** (322, familia, índice, N, realización), con familias 1 = nula, 2 = calibración, 3 = vistas, 4 = bootstrap y 5 = referencia del perfil. Las selladas usan (322, k, N, s).
- La corrida oficial se hace sobre un árbol limpio, después de la auditoría, en tandas de menos de 2 horas.

## 7. Hash del conjunto sellado nuevo (registrado antes de abrirlo)

`CONFIRMACION_CIEGA_2.md`, SHA-256 `b64503ffb0d5f5772a0ab8c8dcf80e689371d75ec6f93a6d93f9912c930866cf`.
- **Contenido:** 4 positivas, 5 negativas, 0 en zona gris.
- **Aviso del guardián:** para fijar la intensidad de un defecto hizo una prueba previa que solo contó relaciones (sin Δr, sin intercepto y sin K), con semillas ajenas a las entropías oficiales.
