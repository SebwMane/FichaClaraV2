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

## 4. Dictamen del Consejo sobre R17-1 a R17-4 (2026-10-19, texto pegado por el usuario)

Las tres revisiones independientes coinciden: **sí a R17-1, R17-2 y R17-3 con límites explícitos; R17-4 va a una propuesta y un preregistro separados.** El Consejo aclara que su dictamen **no registra una firma en nombre del usuario**.

- **Precisión clave:** el prototipo de intervalo máximo no cambia R16-1 por sí solo. Pero ponderar órdenes con S_SJ^max en el templado paralelo **sí crea un ensamble distinto** del firmado en R16-1. Ese piloto valida el método, no prueba H_sel. ES-E1 sigue en pausa.
- **R17-1:**
  - fijar reglas reproducibles para anticadenas, empates y casos sin modos;
  - no asignar cero por defecto;
  - definir el soporte antes de muestrear;
  - mantener el candidato separado del E1′ firmado.
- **R17-2:**
  - empezar solo después de fijar el objetivo y su soporte;
  - verificar balance detallado, intercambios y viajes completos entre valores de λ;
  - medir la mezcla también en la estructura, no solo en S;
  - el paso perezoso ayuda con la paridad, pero no demuestra mezcla.
- **R17-3:**
  - investigar en paralelo, comparando estado por estado con el cálculo completo;
  - dejar explícito si se acelera el S original o el candidato;
  - una aproximación que cambie el valor o los pesos requiere otra definición.
- **Orden:** R17-1 y R17-3 en paralelo; R17-2 después; revisión del Consejo antes de cualquier nuevo preregistro de ES-E1. N ≤ 64 sigue siendo solo desarrollo técnico.
- **Corrección factual confirmada** (§2): una cadena, la que parte del orden vacío, en cada una de las dos celdas de E-b con N = 64 y λ ≠ 0.

## 5. Especificación propuesta por el orquestador para R17-1 y R17-3 (borrador; espera firma, nada implementado)

### 5.1 R17-1: candidato S_SJ^max (separado de E1′)

| Caso | Regla propuesta |
|---|---|
| Región | Entre todos los intervalos I[a, b] con a < b (con extremos), el de mayor \|I\|. Así existe región en todo orden con al menos una relación, incluidos los de altura 2 de la trampa de ES-0 (sus links dan \|I\| = 2). |
| Empates | Igual que E1′: se maximiza min(\|pasado(a)\|, \|futuro(b)\|); si persiste el empate, los menores índices (a, b). |
| Anticadena | Es el único orden sin relaciones (un solo estado con etiquetado natural). **Soporte declarado:** todos los órdenes de N elementos **menos la anticadena**. Es una exclusión de un estado, explícita y verificable, no una imputación. La conectividad del espacio sin la anticadena se prueba (los movimientos de quitar un link nunca la necesitan como paso intermedio salvo desde órdenes con una sola relación, que se conectan entre sí pasando por órdenes de dos relaciones; se verifica por enumeración para N ≤ 5). |
| Sin modos tras el truncamiento | No se imputa cero. Se proponen dos variantes, a elegir antes de muestrear: **(i)** el soporte excluye además los órdenes sin modos, con la fracción medida en muestreo uniforme y conectividad verificada; **(ii)** el truncamiento de la región se aplica con un mínimo de un modo (el de mayor \|ν\|). Recomiendo **(i)** si la fracción sin modos es pequeña en la medida base; si no lo es, el candidato se descarta. |
| Registro | Cada evaluación guarda \|R\|, número de modos conservados y la forma del intervalo (altura y fracción de relaciones dentro de R). |

**Pruebas antes de muestrear:** enumeración exacta para N ≤ 5 (soporte, conectividad y balance detallado con el peso nuevo); comparación de S_SJ^max con S_SJ de E1′ en sprinklings 2D donde ambas regiones existen, solo descriptiva.

### 5.2 R17-3: aceleración exacta de S_SJ original (E1′)

- **Objetivo:** el **S original de E1′** (el candidato comparte la parte global y se beneficia después). Criterio de aceptación: diferencia con el cálculo completo ≤ 1e-9 en cada estado de cadenas de prueba (E-a y E-b, N = 32 a 256), incluidas las propuestas rechazadas.
- **Paso 1, perfil:** medir cuánto tiempo va a la diagonalización global, a `interval_balanced` (producto de matrices) y a `ssee`, con N de 64 a 512.
- **Paso 2, candidatos exactos:**
  - conteos de intervalos actualizados en O(N) por movimiento (ya existe en `mcmc_general.py`);
  - actualización de bajo rango de la diagonalización global: agregar o quitar una relación cambia iΔ en rango 2;
  - recálculo de `ssee` solo cuando cambian la región o los modos que la afectan.
- **Límite declarado:** el truncamiento hace discontinuo a S, porque un modo puede cruzar el umbral. Una actualización exacta debe reproducir esos saltos; si un método solo los aproxima, deja de ser R17-3 y pasa a ser un sustituto con hipótesis propia.
