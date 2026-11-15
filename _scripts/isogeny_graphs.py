#!/usr/bin/env python3
"""Compute the supersingular isogeny graphs drawn on the home page.

Writes assets/data/isogeny-graphs.json. Pure Python, no dependencies
(numpy is used only for the optional spectral statistics).

    python3 _scripts/isogeny_graphs.py            # the default list of primes
    python3 _scripts/isogeny_graphs.py 431 1009   # your own primes

Three graphs are computed for each prime p:

  "2"    the supersingular 2-isogeny graph G_p(2)
  "23"   the (2,3)-isogeny graph G_p({2,3}): same vertices, 2- and 3-isogenies
  "2b3"  the 2-isogeny graph with Borel (Gamma_0(3)) level structure:
         vertices are pairs (E, C), C < E of order 3, up to isomorphism;
         a 2-isogeny phi: E -> E' gives an edge (E, C) -> (E', phi(C)).

Method. Everything is done with short Weierstrass models y^2 = x^3 + ax + b
over F_{p^2} and Velu's formulas, so only F_{p^2} arithmetic and root finding
for polynomials of degree <= 4 are needed.

  * One supersingular curve is found by reducing a class-number-one CM curve
    at a prime that is inert or ramified in its CM field; the rest are found
    by walking the (connected) 2-isogeny graph.
  * A subgroup C of order 3 is {O, +-P}; it is determined by x(P) = x0, a root
    of the 3-division polynomial psi_3 = 3x^4 + 6ax^2 + 12bx - a^2.
    Isomorphisms act by (a, b, x0) -> (t^2 a, t^3 b, t x0), so for x0 != 0 the
    pair (a/x0^2, b/x0^3) is a complete invariant of (E, C). (x0 = 0 forces
    a = 0, i.e. j = 0, and is a single vertex.)
  * Edges are *directed* and counted by kernels, so every vertex has
    out-degree ell + 1 (3 for the level-structure graph). At vertices with
    extra automorphisms (j = 0, 1728 and their lifts) the adjacency matrix is
    not symmetric; the drawing shows max(A_ij, A_ji) undirected edges.

Checks run on every prime (the script aborts if one fails):
  * number of supersingular j-invariants agrees with the class number formula;
  * every Velu output is a root of the classical modular polynomial
    Phi_2 or Phi_3;
  * pushed-forward level structure is again 3-torsion (psi_3 vanishes);
  * A_ij * |Aut_j| = A_ji * |Aut_i| (the Brandt-matrix symmetry);
  * number of level-structure vertices agrees with an orbit count giving
    total mass 4 * (p - 1) / 24.
"""
from __future__ import annotations

import json
import random
import sys
from fractions import Fraction
from pathlib import Path

DEFAULT_PRIMES = [101, 157, 211, 307, 431, 599, 827, 1009]
# The level-structure graph has about four times as many vertices, so it is only
# drawn (and only stored) up to this prime. Keep in step with maxP in
# assets/js/isogeny-hero.js. It is still computed and checked for every prime.
LEVEL_MAX_P = 431

# j-invariants of the CM curves with class number one: discriminant -> j
CM_J = {
    -3: 0, -4: 1728, -7: -3375, -8: 8000, -11: -32768, -19: -884736,
    -43: -884736000, -67: -147197952000, -163: -262537412640768000,
}


# --------------------------------------------------------------------------
# F_{p^2} = F_p[i] / (i^2 - r), r a quadratic non-residue. Elements: (a, b).
# --------------------------------------------------------------------------
class Fp2:
    def __init__(self, p: int):
        self.p = p
        if p % 4 == 3:
            r = p - 1  # i^2 = -1
        else:
            r = next(n for n in range(2, p) if pow(n, (p - 1) // 2, p) == p - 1)
        self.r = r
        self.q = p * p
        self.zero = (0, 0)
        self.one = (1, 0)

    def el(self, n: int):
        return (n % self.p, 0)

    def add(self, x, y):
        return ((x[0] + y[0]) % self.p, (x[1] + y[1]) % self.p)

    def sub(self, x, y):
        return ((x[0] - y[0]) % self.p, (x[1] - y[1]) % self.p)

    def neg(self, x):
        return ((-x[0]) % self.p, (-x[1]) % self.p)

    def mul(self, x, y):
        p = self.p
        return ((x[0] * y[0] + self.r * x[1] * y[1]) % p,
                (x[0] * y[1] + x[1] * y[0]) % p)

    def smul(self, n: int, x):
        return ((n * x[0]) % self.p, (n * x[1]) % self.p)

    def inv(self, x):
        p = self.p
        n = (x[0] * x[0] - self.r * x[1] * x[1]) % p
        ni = pow(n, p - 2, p)
        return ((x[0] * ni) % p, (-x[1] * ni) % p)

    def div(self, x, y):
        return self.mul(x, self.inv(y))

    def pow(self, x, n: int):
        out = self.one
        while n:
            if n & 1:
                out = self.mul(out, x)
            x = self.mul(x, x)
            n >>= 1
        return out

    def frob(self, x):
        return (x[0], (-x[1]) % self.p)

    def rand(self):
        return (random.randrange(self.p), random.randrange(self.p))

    def show(self, x) -> str:
        a, b = x
        if b == 0:
            return str(a)
        if a == 0:
            return "i" if b == 1 else f"{b}i"
        return f"{a} + i" if b == 1 else f"{a} + {b}i"


# --------------------------------------------------------------------------
# Polynomials over F_{p^2}: lists of coefficients, lowest degree first.
# --------------------------------------------------------------------------
def ptrim(F, f):
    while f and f[-1] == F.zero:
        f = f[:-1]
    return f


def pmod(F, f, m):
    f = list(f)
    dm = len(m) - 1
    lead_inv = F.inv(m[-1])
    while len(f) - 1 >= dm and f:
        c = F.mul(f[-1], lead_inv)
        shift = len(f) - 1 - dm
        for k in range(dm + 1):
            f[shift + k] = F.sub(f[shift + k], F.mul(c, m[k]))
        f = ptrim(F, f)
    return f


def pmulmod(F, f, g, m):
    if not f or not g:
        return []
    out = [F.zero] * (len(f) + len(g) - 1)
    for a, fa in enumerate(f):
        if fa == F.zero:
            continue
        for b, gb in enumerate(g):
            out[a + b] = F.add(out[a + b], F.mul(fa, gb))
    return pmod(F, ptrim(F, out), m)


def ppowmod(F, f, n, m):
    out = [F.one]
    f = pmod(F, f, m)
    while n:
        if n & 1:
            out = pmulmod(F, out, f, m)
        f = pmulmod(F, f, f, m)
        n >>= 1
    return out


def psub(F, f, g):
    n = max(len(f), len(g))
    f = f + [F.zero] * (n - len(f))
    g = g + [F.zero] * (n - len(g))
    return ptrim(F, [F.sub(a, b) for a, b in zip(f, g)])


def pgcd(F, f, g):
    f, g = ptrim(F, list(f)), ptrim(F, list(g))
    while g:
        f, g = g, pmod(F, f, g)
    if not f:
        return f
    li = F.inv(f[-1])
    return [F.mul(c, li) for c in f]


def peval(F, f, x):
    out = F.zero
    for c in reversed(f):
        out = F.add(F.mul(out, x), c)
    return out


def roots(F, f):
    """Distinct roots of f in F_{p^2}."""
    f = ptrim(F, list(f))
    x = [F.zero, F.one]
    g = pgcd(F, psub(F, ppowmod(F, x, F.q, f), x), f)  # split, squarefree part
    out = []

    def split(g):
        d = len(g) - 1
        if d == 0:
            return
        if d == 1:
            out.append(F.neg(F.div(g[0], g[1])))
            return
        while True:
            delta = F.rand()
            h = ppowmod(F, [delta, F.one], (F.q - 1) // 2, g)
            h = pgcd(F, psub(F, h, [F.one]), g)
            if 0 < len(h) - 1 < d:
                break
        split(h)
        # g / h by repeated root removal is unnecessary: gcd with the cofactor
        cof = pdivexact(F, g, h)
        split(cof)

    split(g)
    return sorted(out)


def pdivexact(F, f, g):
    f = list(f)
    dg = len(g) - 1
    li = F.inv(g[-1])
    quo = [F.zero] * (len(f) - dg)
    for k in range(len(f) - dg - 1, -1, -1):
        c = F.mul(f[k + dg], li)
        quo[k] = c
        for t in range(dg + 1):
            f[k + t] = F.sub(f[k + t], F.mul(c, g[t]))
    assert not ptrim(F, f), "inexact division"
    return quo


# --------------------------------------------------------------------------
# Elliptic curves y^2 = x^3 + ax + b
# --------------------------------------------------------------------------
def j_inv(F, a, b):
    a3 = F.smul(4, F.mul(a, F.mul(a, a)))
    den = F.add(a3, F.smul(27, F.mul(b, b)))
    return F.div(F.smul(1728, a3), den)


def model(F, j):
    """A fixed model for each j, so that level structures can be compared."""
    if j == F.zero:
        return (F.zero, F.one)
    if j == F.el(1728):
        return (F.one, F.zero)
    k = F.div(j, F.sub(F.el(1728), j))
    return (F.smul(3, k), F.smul(2, k))


def two_torsion(F, a, b):
    return roots(F, [b, a, F.zero, F.one])


def three_torsion(F, a, b):
    """x-coordinates of the four subgroups of order 3."""
    return roots(F, [F.neg(F.mul(a, a)), F.smul(12, b), F.smul(6, a), F.zero, F.el(3)])


def psi3(F, a, b, x):
    return peval(F, [F.neg(F.mul(a, a)), F.smul(12, b), F.smul(6, a), F.zero, F.el(3)], x)


def velu2(F, a, b, e):
    """Quotient by <(e, 0)>. Returns (a', b', xmap)."""
    t = F.add(F.smul(3, F.mul(e, e)), a)
    w = F.mul(e, t)
    a2 = F.sub(a, F.smul(5, t))
    b2 = F.sub(b, F.smul(7, w))
    return a2, b2, (lambda x: F.add(x, F.div(t, F.sub(x, e))))


def velu3(F, a, b, x0):
    """Quotient by the subgroup {O, +-P} with x(P) = x0."""
    t = F.smul(2, F.add(F.smul(3, F.mul(x0, x0)), a))
    y2 = F.add(F.mul(x0, F.add(F.mul(x0, x0), a)), b)
    u = F.smul(4, y2)
    w = F.add(u, F.mul(x0, t))
    return F.sub(a, F.smul(5, t)), F.sub(b, F.smul(7, w))


def phi2(F, x, y):
    m, a, s = F.mul, F.add, F.smul
    x2, y2 = m(x, x), m(y, y)
    terms = [
        m(x2, x), m(y2, y), F.neg(m(x2, y2)),
        s(1488, a(m(x2, y), m(x, y2))),
        s(-162000, a(x2, y2)),
        s(40773375, m(x, y)),
        s(8748000000, a(x, y)),
        F.el(-157464000000000),
    ]
    out = F.zero
    for t in terms:
        out = a(out, t)
    return out


def phi3(F, x, y):
    m, a, s = F.mul, F.add, F.smul
    x2, y2 = m(x, x), m(y, y)
    x3, y3 = m(x2, x), m(y2, y)
    terms = [
        m(x2, x2), m(y2, y2), F.neg(m(x3, y3)),
        s(2232, a(m(x3, y2), m(x2, y3))),
        s(-1069956, a(m(x3, y), m(x, y3))),
        s(36864000, a(x3, y3)),
        s(2587918086, m(x2, y2)),
        s(8900222976000, a(m(x2, y), m(x, y2))),
        s(452984832000000, a(x2, y2)),
        s(-770845966336000000, m(x, y)),
        s(1855425871872000000000, a(x, y)),
    ]
    out = F.zero
    for t in terms:
        out = a(out, t)
    return out


def kronecker(d: int, p: int) -> int:
    d %= p
    if d == 0:
        return 0
    return 1 if pow(d, (p - 1) // 2, p) == 1 else -1


def n_supersingular(p: int) -> int:
    return p // 12 + {1: 0, 5: 1, 7: 1, 11: 2}[p % 12]


def aut_order(F, j) -> int:
    if j == F.zero:
        return 6
    if j == F.el(1728):
        return 4
    return 2


# --------------------------------------------------------------------------
# The graphs
# --------------------------------------------------------------------------
def supersingular_js(F):
    p = F.p
    seed = next(F.el(j) for d, j in CM_J.items() if kronecker(d, p) != 1)
    seen, queue = {seed}, [seed]
    while queue:
        j = queue.pop()
        a, b = model(F, j)
        tors = two_torsion(F, a, b)
        assert len(tors) == 3
        for e in tors:
            a2, b2, _ = velu2(F, a, b, e)
            j2 = j_inv(F, a2, b2)
            if j2 not in seen:
                seen.add(j2)
                queue.append(j2)
    js = sorted(seen, key=lambda x: (x[1] != 0, x))  # F_p-rational first
    assert len(js) == n_supersingular(p), (p, len(js), n_supersingular(p))
    return js


def isogeny_adjacency(F, js, ell):
    idx = {j: k for k, j in enumerate(js)}
    n = len(js)
    A = [[0] * n for _ in range(n)]
    for j in js:
        a, b = model(F, j)
        if ell == 2:
            tors = two_torsion(F, a, b)
            assert len(tors) == 3
            for e in tors:
                a2, b2, _ = velu2(F, a, b, e)
                j2 = j_inv(F, a2, b2)
                assert phi2(F, j, j2) == F.zero, "Velu/Phi_2 mismatch"
                A[idx[j]][idx[j2]] += 1
        else:
            tors = three_torsion(F, a, b)
            assert len(tors) == 4
            for x0 in tors:
                a2, b2 = velu3(F, a, b, x0)
                j2 = j_inv(F, a2, b2)
                assert phi3(F, j, j2) == F.zero, "Velu/Phi_3 mismatch"
                A[idx[j]][idx[j2]] += 1
    aut = [aut_order(F, j) for j in js]
    for s in range(n):
        for t in range(n):
            assert A[s][t] * aut[t] == A[t][s] * aut[s], "Brandt symmetry fails"
    return A


def level_invariant(F, a, b, x0):
    if x0 == F.zero:
        return ("zero",)
    xi = F.inv(x0)
    xi2 = F.mul(xi, xi)
    return (F.mul(a, xi2), F.mul(b, F.mul(xi2, xi)))


def borel3_graph(F, js):
    """2-isogeny graph on pairs (E, C), C of order 3."""
    verts = {}    # invariant -> index
    reps = []     # (a, b, x0, j, j')
    orbit = []    # number of x0 on the fixed model giving this vertex
    for j in js:
        a, b = model(F, j)
        tors = three_torsion(F, a, b)
        assert len(tors) == 4
        for x0 in tors:
            inv = level_invariant(F, a, b, x0)
            if inv in verts:
                orbit[verts[inv]] += 1
                continue
            verts[inv] = len(reps)
            a3, b3 = velu3(F, a, b, x0)
            reps.append((a, b, x0, j, j_inv(F, a3, b3)))
            orbit.append(1)
    n = len(reps)
    # |Aut(E, C)| = |Aut(E)| / orbit size; total mass must be 4 (p - 1) / 24.
    aut = [Fraction(aut_order(F, r[3]), orbit[k]) for k, r in enumerate(reps)]
    assert all(x.denominator == 1 for x in aut)
    assert sum(1 / x for x in aut) == Fraction(4 * (F.p - 1), 24), "mass formula fails"
    A = [[0] * n for _ in range(n)]
    for s, (a, b, x0, j, _) in enumerate(reps):
        for e in two_torsion(F, a, b):
            a2, b2, xmap = velu2(F, a, b, e)
            x1 = xmap(x0)
            assert psi3(F, a2, b2, x1) == F.zero, "image is not 3-torsion"
            # move to the fixed model of the target j-invariant
            inv = level_invariant(F, a2, b2, x1)
            assert inv in verts, "unknown level-structure vertex"
            A[s][verts[inv]] += 1
    for s in range(n):
        assert sum(A[s]) == 3
        for t in range(n):
            assert A[s][t] * aut[t] == A[t][s] * aut[s], "level Brandt symmetry fails"
    fp = [F.frob(r[3]) == r[3] and F.frob(r[4]) == r[4] and
          all(F.frob(c) == c for c in level_invariant(F, r[0], r[1], r[2]) if c != "zero")
          for r in reps]
    labels = [[F.show(r[3]), F.show(r[4])] for r in reps]
    return A, labels, fp


def links_from(A, ell):
    n = len(A)
    out = []
    for s in range(n):
        if A[s][s]:
            out.append([s, s, ell, A[s][s]])
        for t in range(s + 1, n):
            m = max(A[s][t], A[t][s])
            if m:
                out.append([s, t, ell, m])
    return out


def second_eigenvalue(mats):
    try:
        import numpy as np
    except ImportError:
        return None
    M = sum(np.array(A, dtype=float) for A in mats)
    ev = sorted(abs(np.linalg.eigvals(M)), reverse=True)
    return round(float(ev[1]), 4) if len(ev) > 1 else None


def compute(p: int):
    F = Fp2(p)
    js = supersingular_js(F)
    A2 = isogeny_adjacency(F, js, 2)
    A3 = isogeny_adjacency(F, js, 3)
    AB, lab_b, fp_b = borel3_graph(F, js)
    lab = [F.show(j) for j in js]
    fp = [j[1] == 0 for j in js]
    lam = {"2": second_eigenvalue([A2]), "23": second_eigenvalue([A2, A3]),
           "2b3": second_eigenvalue([AB])}
    # Ramanujan-type bounds, checked numerically rather than assumed.
    for key, bound in (("2", 2 * 2 ** 0.5), ("23", 2 * 2 ** 0.5 + 2 * 3 ** 0.5),
                       ("2b3", 2 * 2 ** 0.5)):
        if lam[key] is not None:
            assert lam[key] <= bound + 1e-6, (p, key, lam[key], bound)
    out = {
        "p": p,
        "nonresidue": F.r if F.r != p - 1 else -1,
        "n_level": len(lab_b),
        "graphs": {
            "2": {"nodes": lab, "fp": fp, "links": links_from(A2, 2), "lambda2": lam["2"]},
            "23": {"nodes": lab, "fp": fp,
                   "links": links_from(A2, 2) + links_from(A3, 3), "lambda2": lam["23"]},
            "2b3": {"nodes": lab_b, "fp": fp_b, "links": links_from(AB, 2),
                    "lambda2": lam["2b3"]},
        },
    }
    if p > LEVEL_MAX_P:
        del out["graphs"]["2b3"]
    return out


def main():
    random.seed(0)
    primes = [int(a) for a in sys.argv[1:]] or DEFAULT_PRIMES
    data = {"primes": primes, "by_prime": {}}
    for p in primes:
        res = compute(p)
        data["by_prime"][str(p)] = res
        g = res["graphs"]
        print(f"p = {p:5d}: {len(g['2']['nodes']):4d} curves, "
              f"{res['n_level']:4d} level-structure vertices, "
              f"lambda2 = {g['2']['lambda2']}, {g['23']['lambda2']}")
    out = Path(__file__).resolve().parent.parent / "assets" / "data" / "isogeny-graphs.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    print(f"wrote {out} ({out.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
