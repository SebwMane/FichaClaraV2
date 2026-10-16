# Asesor de inferencia: sobreajuste en la línea F1 a F7 y plan para F6

## 1. Grados de libertad ya gastados

**Criterios con voto en la línea 4D:** F3, F5 (brazos A, B, C), F5′ y F7. Son 6 criterios y todos FALLAN. F1 (2D) pasa con márgenes estrechos.

**Bifurcaciones decididas después de ver datos (7):** (1) links en vez de matriz causal (Rev. 7); (2) d* espectral en vez de Myrheim-Meyer (tras F5); (3) FRW en la calibración (tras de Sitter 1.21); (4) bloque fuera del alcance (tras F5); (5) d* por mediana de estructura (tras la cuenta previa sobre F5); (6) D_forma (tras F5′); (7) ahora, p = 3 fuera de la exigencia (tras F7).

**Perillas fijadas sin barrido declarado (~8):** K_MAX = 8, el +1 del logaritmo, norma máxima, factor 1.5, percentil 95 sobre 60 calibraciones (casi el máximo), N ∈ {1600, 2400}, corte de redondeo de d*, regla de ventanas de F3.

**Diagnóstico.** El registro oficial está limpio (cada FALLA quedó como FALLA); el riesgo está en lo que viene. El criterio actual está esculpido contra dos familias de impostores ya vistas (órdenes producto y conos L^p) y dos de curvatura ya vistas (de Sitter y FRW). Toda la curvatura probada es conformemente plana: solo cambia la densidad, nunca la relación causal. Un PASA sobre estructuras ya vistas no valdría casi nada. **Riesgo de sobreajuste para la próxima prueba: alto**, salvo que se separen exploración y confirmación.

## 2. Registro de la decisión sobre p = 3

Texto propuesto:

> **Enmienda E-R11-1 (fecha), firmada por Sebastían.** Para Ω no se exige distinguir un cono casi redondo del redondo. Se fija la tolerancia geométrica ε = max sobre la esfera unidad de |‖u‖_p / ‖u‖_2 − 1|. Los conos con ε ≤ ε₀ son «zona gris»: no votan como positivas ni como negativas. ε₀ = (valor) se fija hoy y no se modifica después.
> **Se tomó después de ver F7**, y la frontera elegida coincide con la resolución observada del instrumento. Por eso no cuenta como éxito independiente.
> **¿Puede crear un positivo?** Sí: quita negativas que el criterio fallaba. Por eso vale solo para pruebas futuras con semillas nuevas.
> **No reabre F7.** F7 sigue siendo FALLA.

Lo que **no** puede hacerse con la enmienda:
- Releer F7 como PASA. Además, F7 seguiría fallando: de Sitter fuerte (1.55) es independiente de esta decisión, y p = 1.5 tiene una deformación comparable a la de p = 3 (ε ≈ 0.20 frente a 0.17; son exponentes conjugados). Si se quiere cubrir p = 1.5, hay que decirlo explícitamente, como otra decisión tomada después de ver los datos.
- Mover ε₀ después de ver la próxima prueba.
- Sacar del alcance la curvatura fuerte sin declararlo como otra reducción posterior a una falla.

## 3. Explorar tamaños m sin sobreajuste

**Sesgos propios de los sub-intervalos.** Con m fijo, el volumen propio del sub-intervalo varía con la densidad. En regiones curvas, los intervalos de m elementos tienen tamaños distintos, y el efecto de la curvatura escala como (tamaño/radio)². Los intervalos cerca del borde de la caja quedan recortados. Los sub-intervalos se solapan, así que no son 300 muestras independientes: la unidad estadística es la corrida, no el sub-intervalo. Hay que informar por estratos de posición (interior o borde) y de volumen.

**Fase exploratoria.** Usa un bloque de entropías declarado como piloto y no tiene voto. Puede mirarlo todo: una rejilla fija m ∈ {6, 10, 15, 20, 30, 40, 60}, las distancias D, D_forma y la componente k = 0, y todas las estructuras ya vistas. La rejilla se declara antes de empezar y no se amplía porque un resultado quedó «casi».

**Qué se congela antes de confirmar** (en un commit con el hash del código):
- un m principal (los demás m se informan sin voto);
- la regla de selección de sub-intervalos y el máximo por corrida;
- K_MAX relativo a m;
- la distancia;
- la regla del umbral, con al menos 150 calibraciones para que el p95 no sea el máximo;
- la regla de d*_loc;
- el mínimo para «no medible»;
- N y ε₀;
- la lista de estructuras con su papel;
- la regla de veredicto.

**Comparaciones múltiples.** Si más de un m vota, no se acepta con «cualquier m pasa», porque eso multiplica las fugas de impostores. Lo limpio es calibrar el umbral sobre el máximo de las razones D/δ a través de los m (estilo max-T), de modo que la tasa de rechazo de positivas de calibración sea 5% para la familia completa.

**Reglas de parada y presupuesto.** La confirmación es un solo disparo. Si falla, se registra, y un reintento exige entropías nuevas e impostores a ciegas nuevos. Presupuesto propuesto: dos confirmaciones para la línea del perfil local, y después vuelve al Consejo. Cualquier afirmación futura cita el número de intentos («PASA en el intento 2 de 2; seis criterios previos fallaron»).

**Reserva (holdout).** Un subagente que no participa en el diseño define ya un conjunto de estructuras y un bloque de entropías. En el repositorio se publica solo el hash SHA-256 del archivo, que se guarda fuera del repositorio con Sebastían. El conjunto se abre una sola vez, después de un PASA confirmatorio, e idealmente con N mayor (de 4000 a 8000), para comprobar el riesgo de que la distancia de de Sitter crezca con N.

## 4. Estructuras nuevas para que la confirmación valga

Las define el subagente a ciegas, con ε lejos de ε₀ salvo un par cerca de la frontera.

**Negativas:**
- Cono de sección cilíndrica, Δt ≥ max(‖(Δx,Δy)‖₂, |Δz|), con ε ≈ 0.29: suave en parte y con aristas en parte.
- Cono con p distinto por eje.
- Minkowski con 10% a 15% de links borrados, que sigue siendo un orden parcial. Es el impostor adversarial directo de la falla 2, porque imita el hueco en k = 0 de de Sitter fuerte.

**Positivas:**
- Cono elipsoidal (velocidad de la luz distinta por eje). Es Minkowski tras un cambio lineal, así que debe pasar: es una prueba de invariancia.
- Minkowski con densidad variable al azar y suave. Es conformemente plano, así que también es un espacio-tiempo.
- **Al menos uno no conformemente plano** (Kasner o un campo débil), con la relación causal calculada por geodésicas. Sin esto, «reconoce curvatura» solo significa «tolera cambios de densidad».

## 5. Sesgos que conviene vigilar al supervisar

1. **Mover el blanco:** reducir el alcance justo donde el instrumento falló. Cuenten las reducciones.
2. **Reutilizar lo visto:** impostores o positivas ya mirados presentados como evidencia nueva.
3. **Elegir la perilla ganadora:** el m, la distancia, K_MAX o el factor que mejor separa, sin corrección por comparaciones múltiples.
4. **Pseudorréplicas:** contar sub-intervalos solapados como muestras independientes.
5. **Sesgo de tamaño:** m fijo implica volumen variable, y los intervalos de borde están recortados.
6. **Parada oportunista:** «una semilla más» o «un N mayor» cuando una razón queda en 0.95 o en 1.05.
7. **Escrutinio asimétrico:** auditar los FALLA con más rigor que los PASA.
8. **Umbrales frágiles:** pocas calibraciones y un factor 1.5 elegido a ojo.

## Ronda 2

**1. Acuerdos.** Coincido con lo siguiente:
- un solo criterio con voto (ramas);
- calibrar con T = max_m razón_m y medir los errores entre corridas, no entre sub-intervalos (tamaño finito);
- tomar el exponente de la corrección, 2/d, de una teoría externa y no del ajuste (curvatura).

El hallazgo de que con m = 20 en el bloque los intervalos no son locales confirma el sesgo de tamaño que señalé.

**2. Discrepancias.** Curvatura propone dos ejes, m y dilución q, y dos observables, r y el perfil. Si los dos votan, se repite el jardín de F5 con sus tres brazos. El modelo b + A·x^{1/2} + B·x tiene tres parámetros: con pocos valores de x, el intercepto queda libre para absorber cualquier desviación. En ramas, R-B a R-D «informadas» generan relatos que después se convierten en criterios.

**3. Bifurcaciones nuevas.** Cuento unas 7: la rejilla de m, el eje q, la elección del observable, la forma del modelo, la elección entre estudentizar o usar k ≤ m/2, la regla de N según m y n_sel. Propongo recortar así:
- **Un solo observable con voto:** el perfil local con k ≤ m/2. No estudentizar, porque agrega una varianza estimada y otra perilla. r queda como informativo.
- **Un solo eje:** m, con tres valores fijados antes y max-T. q va solo a la fase exploratoria. Diluir un Poisson sigue siendo Poisson, así que q es útil como control, pero no debe votar.
- **El modelo, solo informativo,** con el exponente fijado en 1/2 y al menos 5 valores de x. No vota en la primera confirmación.
- **Una regla de N** (N ≳ 800 m o el máximo que permita el costo) y un n_sel fijo, congelados.
- **Ramas sin voto,** anotadas en un registro de intentos. Cada una que pase a votar cuenta en el presupuesto.

**4. Pilotos.** Usaron entropías fuera de la serie oficial (900001 y 900101) y estructuras positivas ya vistas, así que no contaminan la confirmación. Son la fase exploratoria y deben declararse como tal, con todo lo que se miró. Hay dos costos:
- **Anti de Sitter** queda «visto»: ya no sirve como positiva a ciegas, y su desviación inclinó la elección del modelo.
- **Asimetría:** solo se miraron positivas, así que el umbral quedaría ajustado para aceptarlas sin saber cómo se comportan los conos a m chico.

Antes de congelar conviene explorar también los impostores ya vistos (conos p = 1, ∞, 1.5 y 3, órdenes producto) con entropías de piloto. Los impostores nuevos a ciegas se reservan para la confirmación.
