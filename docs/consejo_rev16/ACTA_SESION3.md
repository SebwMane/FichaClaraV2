# Consejo Ω, sesión 3: dictamen sobre R16-1 a R16-4 (texto pegado por el usuario, 2026-10-19)

Registrado tal cual. El enlace al acta completa apuntaba a una ruta local no accesible desde este entorno; se conserva solo el resumen pegado. La verificación del orquestador está en `docs/OMEGA_CONSEJO_REV16_ES_DISENO.md` §7.

---

La sesión concluyó con **tres dictámenes independientes** (Ω-E, ensambles/MCMC y Ω-C) y una ronda de crítica cruzada.

El Consejo recomienda:

- **R16-1:** aceptar μ_λ ∝ e^{λ S_SJ} como un ensamble exploratorio nuevo, no como una dinámica física derivada de Ω-E. El artículo de Sorkin y Yazdi obtuvo la ley de área con un truncamiento espectral específico en sprinklings 1+1D; extender ese resultado a órdenes generales es una hipótesis del proyecto. (arXiv:1611.10281)
- **R16-2:** usar 2-órdenes para calibración y órdenes generales para evaluar selección, declarando la medida de etiquetado y alineando con ella el control KR.
- **R16-3:** hacer ES-0 antes de preregistrar ES-E1, solo como piloto de costo y mezcla. Hay que implementar primero un muestreador que use S_SJ: el `mcmc_general.py` revisado apunta a la acción BD 4D. Además, 100 000 propuestas equivalen a unas 12,3 barridas nominales en N = 128; eso no demuestra mezcla.
- **R16-4:** empezar Ω-C por la rama euclídea. Mantener la lorentziana como estudio posterior independiente: sus amplitudes complejas y su módulo no son probabilidades de fase. Trabajos previos estudian supresión lorentziana para clases específicas, sin resolver el caso general de Ω. (Loomis y Carlip, arXiv:1709.00064)

El piloto **no se ejecutó** y el repositorio **no se modificó**. R15-2 sigue reportada como pendiente y R14-C v2 depende de su cierre.
