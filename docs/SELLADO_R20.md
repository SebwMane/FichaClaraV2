# Sellado de R20 (2026-10-21)

- **Archivo sellado:** `sellado_r20.py`, escrito por el guardián nuevo fuera del repositorio.
- **sha256:** `f6dc1a12c41c5e7438286fe31f4c9fab0fd974343ffa17b58e4f876f70f27b2d`
- **Universo:** 16 estructuras en 4 familias, con 4 candidatas cada una.

  | Familia | Tipo | Identificables |
  |---|---|---|
  | gemelos | se desvanece | 3 de 4 |
  | tripletes | se desvanece | 3 de 4 |
  | conos poliédricos | autosemejante | 4 de 4 |
  | conos suaves | autosemejante | 4 de 4 |

  Votan 14 estructuras. **Se cumple el mínimo de §3.**
- **Cambio de universo declarado por el guardián:** su primera familia de conos suaves tuvo 0 de 4 identificables. El guardián la reemplazó **una sola vez** por otra familia de conos suaves, diseñada y evaluada una vez, sin tocar las otras tres familias. La versión descartada y su registro de umbrales quedan en el recibo, que se abre junto con el archivo sellado.
  - **Nota del orquestador:** esto es una selección por visibilidad a nivel de familia, permitida por la instrucción al guardián (registrar todo cambio). La población objetivo de la familia de conos suaves queda limitada a la que el instrumento puede ver; los conos suaves descartados son, de hecho, «no identificables por este instrumento». Se informará así.
- **Diferencia técnica declarada antes de correr:** para los umbrales, el guardián calculó D con los componentes finitos en los dos vectores de cada par (como `r19.pair_D`). El ejecutor de R20 usa los componentes finitos en las 300 referencias de cada q, y declara no evaluable un brazo al que le falte alguno. La única diferencia práctica esperada es r_40 con q = 1/8, que en los dos casos queda fuera.
- **Aviso esperado del ejecutor:** el validador heredado de R19 advierte por las familias de 4 estructuras (allí la regla era de 2 o 3). En R20 la regla es ≥ 3, así que el aviso es inocuo.
