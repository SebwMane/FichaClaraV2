# Ω — P3: selección por crecimiento — BORRADOR conceptual (no ejecutable)

- **Estado:** el Consejo pidió **preparar, no ejecutar**. Este borrador no es un prerregistro.
- **Antes de escribir código:** el Consejo ratifica la forma de la regla y se prerregistra en un documento propio.

## 1. Hipótesis

C0 mostró que el resultado depende de la condición inicial (paisaje vidrioso) y que los mecanismos locales de energía no generan coherencia global. P3 cambia la pregunta: la geometría no la selecciona el mínimo de una acción, sino la **historia de construcción**,

    Ω_0 → Ω_1 → Ω_2 → … ,

en la que la estructura ya construida decide dónde pueden aparecer elementos o relaciones nuevas.

## 2. Por qué podría atacar el punto B (y por qué podría no hacerlo)

**A favor.** Un crecimiento por adición local, en el que un nodo nuevo se conecta a un vecindario existente y sus vecinos, extiende la estructura de forma **contigua**. Eso evita por construcción el pegado incoherente que P1-D detecta en los retazos.

**En contra (riesgos conocidos):**
- los procesos de adhesión local suelen producir árboles o estructuras fractales (dimensión no entera, como la agregación limitada por difusión) o mundos pequeños, si hay cierre triádico con aleatoriedad;
- la dimensión podría quedar fijada por la regla, por ejemplo por el número de vecinos a los que se adhiere el nodo nuevo. Eso violaría el criterio de meseta C1-D.

## 3. Preguntas que el Consejo debe responder antes del prerregistro

### 3.1 Las seis preguntas del Consejo (tras P1-D.2)

1. ¿Qué crece? Nodos, aristas, pesos, regiones u otra cosa.
2. ¿Qué determina una nueva conexión? Una regla independiente de «conectar lo que está cerca en 3D».
3. **¿De dónde surge la noción de cercanía?** Es probablemente la más importante: «añadir conexiones cercanas» presupone la geometría que se quiere explicar.
4. ¿Qué evita las cliques? S0 y C0 mostraron una tendencia fuerte a agruparse.
5. ¿Qué hace que la estructura sea globalmente coherente? No basta con vecindarios razonables.
6. **¿Qué principio selecciona la escala?** Aun con crecimiento coherente, hay que explicar por qué aparece una región donde D_B ≈ constante. Esta pregunta forma parte obligatoria del prerregistro de P3.

Nota: P1-D.2 muestra que existen geometrías gruesas con decoración periódica (caveman, retículo de cliques) que son coherentes a escala intermedia. Un proceso de crecimiento que produzca solo eso respondería a 5 y 6 pero no sería una variedad; el prerregistro debe decir cuál de los dos objetivos persigue (gruesa o de variedad).

### 3.2 Preguntas operativas previas

1. ¿El crecimiento añade **nodos**, **relaciones** o ambos?
2. ¿La regla de adhesión usa solo información relacional (prohibido: coordenadas, embedding, D, k = 6)?
3. ¿Hay **relajación** (dinámica tipo C0) entre adiciones, o es puro crecimiento?
4. **Causalidad.** El orden de adición induce un orden parcial, y el Consejo aplazó la causalidad. ¿Se trata ese orden como mero índice de construcción, sin interpretación física?
5. ¿Qué control de **tamaño** y qué control **nulo**? Por ejemplo, el mismo proceso con adhesión a nodos elegidos al azar.

## 4. Puertas propuestas (por analogía con C0)

| Paso | Contenido |
|---|---|
| P3-L1 | Análisis de la regla: ¿qué estructuras deterministas genera? Si genera árboles o cliques, muere |
| P3-L2 | Preflight A–H a N ≤ 343, con especial atención a B (ρ), D (meseta de D_B) y H (retazos) |
| P3-L3 | Robustez frente a variaciones de la regla (meseta C1-D) |
| P3-F1 | N = 729 y 4096 con certificado |

**Muerte:** si el crecimiento solo produce árboles, cliques, mundos pequeños o una dimensión fijada por un parámetro entero de la regla.

## 5. Observación sobre el estado y la regla de formación

El Consejo señala que Ω_final podría no contener información suficiente para explicar por qué tiene geometría: dos redes con estadísticas locales similares pueden haber surgido por un proceso coherente o por ensamblaje aleatorio (el panel de retazos es ese caso). Si esto es cierto, la historia dinámica es una variable físicamente relevante. **Es una hipótesis, no un resultado**; P3 sería la prueba.
