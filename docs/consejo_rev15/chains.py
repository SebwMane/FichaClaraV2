"""Comprobacion de <C_k> (Meyer) en el diamante 4D con n fijo. Solo lectura del repo."""
import sys, math
from fractions import Fraction
import numpy as np
sys.path.insert(0, "/home/user/FichaClaraV2")
from omega.causality.sprinkling import sprinkle_diamond, causal_matrix

D, KMAX = 4, 5

def chi(k, d):
    lg = (k - 1) * (math.lgamma(d + 1) - math.log(2)) + math.lgamma(d / 2) + math.lgamma(d) \
         - math.log(k) - math.lgamma(k * d / 2) - math.lgamma((k + 1) * d / 2)
    return math.exp(lg)

def chi_exact_d4(k):
    # chi_k = 12^(k-1) * 6 / (k (2k-1)! (2k+1)!)
    return Fraction(12 ** (k - 1) * 6, k * math.factorial(2 * k - 1) * math.factorial(2 * k + 1))

def falling(n, k):
    return math.prod(range(n - k + 1, n + 1))

def chain_counts(c, kmax):
    # c[i,j] = i precede a j; v_k[j] = numero de k-cadenas que terminan en j
    ct = c.T.astype(np.float64)
    v = np.ones(c.shape[0]); out = [float(v.sum())]
    for _ in range(2, kmax + 1):
        v = ct @ v; out.append(float(v.sum()))
    return out

def longest(c, t):
    order = np.argsort(t); cs = c[np.ix_(order, order)]
    h = np.ones(len(t), dtype=int)
    for j in range(len(t)):
        p = np.flatnonzero(cs[:j, j])
        if p.size: h[j] = h[p].max() + 1
    return int(h.max())

print("chi_k exactos d=4:", {k: str(chi_exact_d4(k)) for k in range(1, KMAX + 1)})
for k in range(1, KMAX + 1):
    assert abs(float(chi_exact_d4(k)) - chi(k, D)) / chi(k, D) < 1e-12
ss = np.random.SeedSequence(830001)
reps = {100: 4000, 200: 2000}
children = ss.spawn(len(reps))
for (n, R), child in zip(reps.items(), children):
    rng = np.random.default_rng(child)
    data = np.empty((R, KMAX)); L = np.empty(R); frac = np.empty(R)
    for r in range(R):
        x = sprinkle_diamond(n, D, rng); c = causal_matrix(x)
        data[r] = chain_counts(c, KMAX)
        L[r] = longest(c, x[:, 0])
        frac[r] = c.sum() / (n * (n - 1) / 2)
    print(f"\nn={n}, repeticiones={R}")
    print(" k |   teoria n!/(n-k)! chi_k |  media MC  |  err.est. |  z")
    for k in range(2, KMAX + 1):
        th = falling(n, k) * chi(k, D)
        m = data[:, k - 1].mean(); se = data[:, k - 1].std(ddof=1) / math.sqrt(R)
        print(f" {k} | {th:24.4f} | {m:10.4f} | {se:9.4f} | {(m - th) / se:+.2f}")
    print(f" fraccion de orden r: media {frac.mean():.5f} +- {frac.std(ddof=1)/math.sqrt(R):.5f} (teoria 0.1)")
    print(f" cadena mas larga: media {L.mean():.3f}; L/n^(1/4) = {L.mean()/n**0.25:.3f}")
    # Poisson variant only differs by replacing falling(n,k) with (rho V)^k = n^k
    print(" (Poisson: (rhoV)^k chi_k =", [round(n ** k * chi(k, D), 3) for k in range(2, KMAX + 1)], ")")
lo = 2 ** (1 - 1 / D) / math.gamma(1 + 1 / D); hi = 2 ** (1 - 1 / D) * math.e * math.gamma(1 + D) ** (1 / D) / D
print(f"\ncotas Brightwell-Gregory para m_4: {lo:.4f} <= m_4 <= {hi:.4f}")
for d in (2, 3, 4):
    print(d, round(2 ** (1 - 1 / d) / math.gamma(1 + 1 / d), 4), round(2 ** (1 - 1 / d) * math.e * math.gamma(1 + d) ** (1 / d) / d, 4))
