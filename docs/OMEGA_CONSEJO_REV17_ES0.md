# Consejo Ω, Revisión 17: resultados de ES-0 y pausa de ES-E1 (2026-10-19)

- **Acta:** `docs/consejo_rev17/ACTA_SESION4.md` (sesión 4 del Consejo, registrada tal cual).
- **Estado:** **ES-E1 en pausa.** Se retira la convención «sin región, S_SJ = 0». Las firmas R16 siguen registradas. Lo que sigue en §3 son recomendaciones del Consejo **pendientes de firma**.

## 1. Resolución del Consejo (resumen)

| # | Tema | Resolución |
|---|---|---|
| 1 | Estados sin región | Invalidar S_SJ = 0. No adoptar todavía (a) excluir, (b) intervalo máximo ni (c) normalizar. (b) se puede prototipar como técnica; (c) sería una nueva versión de R16-1 con referencia independiente. |
| 2 | Mezcla en órdenes generales | Templado paralelo en λ como primer prototipo, **después** de definir un estadístico válido en todo el soporte. Wang-Landau, después y como contraste. |
| 3 | Escala | N ≤ 64 solo para desarrollo técnico, sin afirmaciones sobre H_sel, H_E ni ley de área. ES-E1 se reabre solo donde coincidan validez del instrumento y mezcla. |

**Lectura:** H_sel queda **indeterminada**, no refutada.

## 2. Verificación del orquestador

- **Las cifras citadas coinciden con `results/omega_es0/summary.json`:**
  - R̂ de S = 1.56 en E-a con N = 32 y λ > 0;
  - ESS de 146 en relaciones para E-b con N = 32 y λ = 0;
  - R̂ de 7.7 y hasta 17 en E-b con N = 64;
  - 25% de muestras finales sin región;
  - η² entre 0.65 y 0.93;
  - 39 s y 7.2 s por barrido con N = 128.
- **Wang-Landau:** confirmado. `omega/causality/wanglandau.py` trabaja con la acción BD 2D entera sobre 2-órdenes; no sirve tal cual para S_SJ continuo en órdenes generales.
- **Restos documentales de la Rev. 16, corregidos en este commit:**
  - la portada ahora indica firmas, ES-0 ejecutado y pausa;
  - §1 dice N/8 ± 5%, con remisión a §7.
- **Una precisión al acta:** dice que en «varias celdas» los movimientos tienen ΔS = 0 dentro de la zona sin región. En los datos, la trampa persistente está en una cadena (la que parte del orden vacío) en cada una de las dos celdas de E-b con N = 64 y λ ≠ 0. Esa cadena aporta el 25% de las muestras finales de su celda. No cambia la conclusión.

## 3. Pasos recomendados por el Consejo (pendientes de firma)

| # | Paso | Tipo |
|---|---|---|
| R17-1 | Prototipo técnico de una **regla total de región** (intervalo máximo como primer candidato), con reglas explícitas para anticadenas, empates y estados sin modos tras el truncamiento; sin imputar cero | desarrollo técnico, N ≤ 64, sin evidencia |
| R17-2 | **Templado paralelo en λ** sobre ese estadístico: balance local y razón de intercambio verificados, escalera con solapamiento, tasas de intercambio, viajes de ida y vuelta de las réplicas, mezcla en S, \|R\|, relaciones, altura y observables estructurales; paso perezoso para la paridad de E-a | desarrollo técnico, N ≤ 64 |
| R17-3 | **Aceleración de S_SJ** que coincida con el cálculo completo en tamaños compartidos (por ejemplo, actualización de bajo rango de la diagonalización tras agregar o quitar una relación). Si en cambio se propone un sustituto aproximado, lleva hipótesis y calibración propias | desarrollo técnico |
| R17-4 | Normalización o residualización por \|R\| solo como **variante de ensamble** (nueva versión de R16-1), con referencia independiente y validación separada | requiere firma propia y preregistro |

**Secuencia:** R17-1 → R17-2, y R17-3 en paralelo. ES-E1 se reabre con un preregistro nuevo, con semillas y controles KR nuevos alineados con la medida de etiquetado natural.
