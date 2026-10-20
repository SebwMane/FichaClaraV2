# Sellado de R19-(b) (2026-10-20)

- **Archivo sellado:** `sellado_r19.py`, escrito por el guardián fuera del repositorio. Se copiará a `tools/sellado_r19.py` después de congelar el ejecutor.
- **sha256:** `e3c0e40c44dbecab8a97eb204a2caffe89d944b1b31d770e0dfd9ccb474dd579`
- **Contenido declarado sin abrirlo:**
  - 2 familias «se desvanece», con 2 estructuras de copias cada una;
  - 2 familias «autosemejante», con 3 y 2 estructuras de región.
- **Verificación del guardián, con sus propias semillas** (99190001, ajenas al proyecto):
  - todas las regiones son órdenes parciales con n = 500 y n = 8000;
  - todos los controles de copias quedan fuera de la banda plana con N ≈ 1000 (D/radio mínimo entre 1.39 y 2.46, con 10 realizaciones cada uno).
  - Dos elecciones iniciales no cumplían la condición; el guardián subió su fracción dentro del rango permitido antes de sellar.
- **Asunto de interfaz señalado por el guardián:** con k ≥ 2, round(frac·N) no siempre es múltiplo de k. **Regla fijada por el orquestador antes de abrir el archivo:** n_originales = round(frac·N/k) y n_copias = k·n_originales.
