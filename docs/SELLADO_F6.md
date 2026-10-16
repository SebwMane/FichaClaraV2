# Conjuntos sellados para la confirmación de F6 (R11-4, 2026-10-14)

Definidos a ciegas por un subagente independiente (guardián de reserva), que no participa en el diseño. **El orquestador no los ha leído.** El repositorio guarda solo los hashes; los archivos los guarda Sebastían.

| Archivo | SHA-256 | Positivas | Negativas | Zona gris | Entropías | Cuándo se abre |
|---|---|---|---|---|---|---|
| `CONFIRMACION_CIEGA.md` | `b8439d85fcf2ea3c6ea76b92cce8caced98483d85dac64bcb7ca43f3fbc91c79` | 3 | 4 | 1 | 20261076 a 20261079 | después de congelar el criterio de la confirmación |
| `RESERVA.md` | `179cf6acf5f3c067d09f14dda5c6201c0ca7e47c96274e03bd1b50b16aa6f562` | 3 | 4 | 1 | 20261090 a 20261099 | una sola vez, después de un PASA confirmatorio (N = 8000) |

Según el guardián, la confirmación cubre:
- una positiva Kasner de vacío;
- una positiva de invariancia lineal y una conformemente plana con densidad aleatoria suave;
- negativas con ε entre 0.25 y 0.35;
- Minkowski con 10% a 15% de links borrados;
- un cono con p distinto por eje;
- ninguna positiva en la zona gris.

Al abrir cada archivo se verifica su hash con `sha256sum` antes de usarlo. Un hash distinto invalida el conjunto.
