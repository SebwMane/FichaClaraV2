"""CF-0: comprobación exacta de la compleción de diamante (CD) de aridad 2 desde una raíz de anchura w (docs/OMEGA_CF0.md §3.2).

Sin azar ni coordenadas. Dos variantes:
  * fresca: si z tiene cubiertas x, y sin cubierta superior común, se añade un elemento NUEVO que cubre exactamente a x e y;
  * identificacion_ingenua: antes de crear, se reutiliza una cubierta u de x que ya cubre a otra cubierta de z (ilustra el
    exceso de identificación).
Se cierra hasta punto fijo o hasta MAXN elementos. Referencia: el retículo booleano de w átomos tiene 2^w elementos.
Salida: results/cf0/cd_check.json.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

MAXN = 2000


def close(w: int, naive_identify: bool) -> dict:
    up: dict[int, set[int]] = {0: set()}
    down: dict[int, set[int]] = {0: set()}
    n = 1
    for _ in range(w):
        up[n], down[n] = set(), {0}
        up[0].add(n)
        n += 1
    changed = True
    while changed and n < MAXN:
        changed = False
        for z in list(up):
            for x, y in itertools.combinations(sorted(up[z]), 2):
                if up[x] & up[y]:
                    continue
                if naive_identify:
                    cand = [u for u in up[x] if any(v in up[z] and v != x for v in down[u])]
                    if cand:
                        u = min(cand)
                        up[y].add(u)
                        down[u].add(y)
                        changed = True
                        continue
                up[n], down[n] = set(), {x, y}
                up[x].add(n)
                up[y].add(n)
                n += 1
                changed = True
                if n >= MAXN:
                    break
            if n >= MAXN:
                break
    tops = [v for v in up if not up[v]]
    return {"elements": n, "closed": n < MAXN, "maximal_elements": len(tops), "boolean_size": 2**w}


def main() -> None:
    out = {f"w={w}": {"fresca": close(w, False), "identificacion_ingenua": close(w, True)} for w in (2, 3, 4)}
    p = Path(__file__).resolve().parent.parent / "results" / "cf0" / "cd_check.json"
    p.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
