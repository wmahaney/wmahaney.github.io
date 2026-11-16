#!/usr/bin/env python3
"""Compute the supersingular isogeny graphs drawn on the home page.

Writes assets/data/isogeny/: one small file p<prime>.json per prime (the page loads
only the prime you pick) and an index.json listing them. Pure Python, no dependencies
(numpy is used only for the optional spectral statistics).

    python3 _scripts/isogeny_graphs.py            # the default list of primes
    python3 _scripts/isogeny_graphs.py 431 1009   # your own primes (any primes >= 5)

Three graphs are computed for each prime p:

  "2"    the supersingular 2-isogeny graph G_p(2)
  "23"   the (2,3)-isogeny graph G_p({2,3}): same vertices, 2- and 3-isogenies
  "2b3"  the 2-isogeny graph with Borel (Gamma_0(3)) level structure:
         vertices are pairs (E, C), C < E of order 3, up to isomorphism;
         a 2-isogeny phi: E -> E' gives an edge (E, C) -> (E', phi(C)).
  "23b5" the (2,3)-isogeny graph with Borel (Gamma_0(5)) level structure:
         vertices are pairs (E, C), C < E of order 5; 2- and 3-isogenies
         push C forward.

Every graph also carries "sigma", the permutation induced by Frobenius
(E -> E^(p), C -> C^(p)). Its fixed vertices are the F_p-rational ones (the "spine");
the page draws the folded graph (the quotient by sigma) from it.

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
  * number of F_p-rational vertices (the spine) agrees with the class number
    formula (h(-4p)/2, h(-p), or 2 h(-p) according to p mod 8);
  * Frobenius is an involutive automorphism of every graph;
  * for level 5: the pushed-forward kernel polynomial is a factor of psi_5, every
    quotient curve is supersingular, and the mass is (p - 1)/4;
  * number of level-structure vertices agrees with an orbit count giving
    total mass 4 * (p - 1) / 24.
"""
from __future__ import annotations

import json
import random
import sys
from fractions import Fraction
from pathlib import Path

DEFAULT_PRIMES = [101, 157, 211, 307, 431, 599, 827, 1009, 1499, 1999, 2503,
                  2999, 3511, 4001, 4507, 4999]
# The level-structure graph has about four times as many vertices. It is stored for
# primes up to this bound (recorded in index.json, which the page reads). It is
# computed and checked for every prime either way.
LEVEL_MAX_P = 5000
# Same for the Gamma_0(5) graph, which has six times as many vertices as G_p(2).
LEVEL5_MAX_P = 3000

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


def pmul(F, f, g):
    if not f or not g:
        return []
    out = [F.zero] * (len(f) + len(g) - 1)
    for a, fa in enumerate(f):
        if fa == F.zero:
            continue
        for b, gb in enumerate(g):
            out[a + b] = F.add(out[a + b], F.mul(fa, gb))
    return ptrim(F, out)


def pmonic(F, f):
    li = F.inv(f[-1])
    return [F.mul(c, li) for c in f]


def quadratic_factors(F, f):
    """Factor a squarefree f, all of whose irreducible factors over F_{p^2} have
    degree 1 or 2, into monic quadratics (as (s, q): x^2 - s x + q).

    Linear factors cannot be paired up here (that needs the curve), so they are
    returned separately as roots."""
    f = pmonic(F, ptrim(F, list(f)))
    x = [F.zero, F.one]
    g1 = pgcd(F, psub(F, ppowmod(F, x, F.q, f), x), f)
    lin = roots(F, g1) if len(g1) > 1 else []
    f2 = pdivexact(F, f, g1) if len(g1) > 1 else f
    quads = []

    def split(g):
        d = len(g) - 1
        if d == 0:
            return
        if d == 2:
            g = pmonic(F, g)
            quads.append((F.neg(g[1]), g[0]))
            return
        while True:
            delta = F.rand()
            h = ppowmod(F, [delta, F.one], (F.q * F.q - 1) // 2, g)
            h = pgcd(F, psub(F, h, [F.one]), g)
            if 0 < len(h) - 1 < d:
                break
        split(h)
        split(pdivexact(F, g, h))

    if len(f2) > 1:
        g2 = pgcd(F, psub(F, ppowmod(F, x, F.q * F.q, f2), x), f2)
        assert len(g2) == len(f2), "factor of degree > 2 in the 5-division polynomial"
        split(g2)
    return lin, sorted(quads)


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


def class_number(D: int) -> int:
    """h(D) for a negative discriminant D, by counting reduced forms."""
    from math import gcd
    n, a = 0, 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or gcd(gcd(a, abs(b)), c) != 1 or (a == c and b < 0):
                continue
            n += 1
        a += 1
    return n


def n_spine(p: int) -> int:
    """Number of F_p-rational supersingular j-invariants (Delfs-Galbraith)."""
    if p % 4 == 1:
        return class_number(-4 * p) // 2
    return class_number(-p) * (1 if p % 8 == 7 else 2)


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
    sigma = [verts[frob_inv(F, level_invariant(F, r[0], r[1], r[2]))] for r in reps]
    check_sigma(A, sigma)
    labels = [[F.show(r[3]), F.show(r[4])] for r in reps]
    return A, labels, sigma


# --------------------------------------------------------------------------
# Level 5. A subgroup C of order 5 is {O, +-P, +-2P}; it is determined by the
# kernel polynomial x^2 - s x + q with roots x(P), x(2P) (a quadratic over F_{p^2}
# even when x(P) is not in F_{p^2}, because Frobenius is a scalar and C is stable).
# Isomorphisms act by (a, b, s, q) -> (t^2 a, t^3 b, t s, t^2 q).
# --------------------------------------------------------------------------
def psi5(F, a, b):
    f = [b, a, F.zero, F.one]
    g = [F.sub(F.smul(-8, F.mul(b, b)), F.mul(a, F.mul(a, a))),
         F.smul(-4, F.mul(a, b)), F.smul(-5, F.mul(a, a)), F.smul(20, b),
         F.smul(5, a), F.zero, F.one]
    p3 = [F.neg(F.mul(a, a)), F.smul(12, b), F.smul(6, a), F.zero, F.el(3)]
    t1 = [F.smul(32, c) for c in pmul(F, pmul(F, f, f), g)]
    t2 = pmul(F, p3, pmul(F, p3, p3))
    return psub(F, t1, t2)


def x_double(F, a, b, x):
    num = F.add(F.sub(F.mul(F.mul(x, x), F.mul(x, x)), F.smul(2, F.mul(a, F.mul(x, x)))),
                F.add(F.smul(-8, F.mul(b, x)), F.mul(a, a)))
    den = F.smul(4, F.add(F.mul(x, F.add(F.mul(x, x), a)), b))
    return F.div(num, den)


def five_subgroups(F, a, b):
    """Kernel polynomials (s, q) of the six subgroups of order 5."""
    lin, quads = quadratic_factors(F, psi5(F, a, b))
    seen = set(lin)
    for r in lin:
        if r in seen and r is not None:
            r2 = x_double(F, a, b, r)
            assert r2 in seen and r2 != r, "5-torsion roots do not pair up"
            seen.discard(r)
            seen.discard(r2)
            quads.append((F.add(r, r2), F.mul(r, r2)))
    assert not seen
    quads = sorted(quads)
    assert len(quads) == 6, len(quads)
    return quads


class Kalg:
    """F_{p^2}[X] / (X^2 - s X + q): elements c0 + c1 X as pairs."""

    def __init__(self, F, s, q):
        self.F, self.s, self.q = F, s, q

    def mul(self, u, v):
        F = self.F
        c11 = F.mul(u[1], v[1])
        return (F.sub(F.mul(u[0], v[0]), F.mul(self.q, c11)),
                F.add(F.add(F.mul(u[0], v[1]), F.mul(u[1], v[0])), F.mul(self.s, c11)))

    def add(self, u, v):
        return (self.F.add(u[0], v[0]), self.F.add(u[1], v[1]))

    def scal(self, c, u):
        return (self.F.mul(c, u[0]), self.F.mul(c, u[1]))

    def trace(self, u):
        F = self.F
        return F.add(F.smul(2, u[0]), F.mul(u[1], self.s))

    def norm(self, u):
        F = self.F
        return F.add(F.add(F.mul(u[0], u[0]), F.mul(F.mul(u[0], u[1]), self.s)),
                     F.mul(F.mul(u[1], u[1]), self.q))

    def inv(self, u):
        F = self.F
        n = self.norm(u)
        assert n != F.zero, "kernel meets the isogeny kernel"
        conj = (F.add(u[0], F.mul(u[1], self.s)), F.neg(u[1]))
        return self.scal(F.inv(n), conj)

    def X(self):
        return (self.F.zero, self.F.one)

    def const(self, c):
        return (c, self.F.zero)

    def char_poly(self, u):
        return self.trace(u), self.norm(u)


def push2(F, a, b, e, s, q):
    """(a', b', s', q') for the image of C under the 2-isogeny with kernel (e, 0)."""
    a2, b2, _ = velu2(F, a, b, e)
    t = F.add(F.smul(3, F.mul(e, e)), a)
    K = Kalg(F, s, q)
    inv = K.inv(K.add(K.X(), K.const(F.neg(e))))
    img = K.add(K.X(), K.scal(t, inv))
    return a2, b2, *K.char_poly(img)


def push3(F, a, b, x0, s, q):
    a2, b2 = velu3(F, a, b, x0)
    t = F.smul(2, F.add(F.smul(3, F.mul(x0, x0)), a))
    u = F.smul(4, F.add(F.mul(x0, F.add(F.mul(x0, x0), a)), b))
    K = Kalg(F, s, q)
    inv = K.inv(K.add(K.X(), K.const(F.neg(x0))))
    img = K.add(K.add(K.X(), K.scal(t, inv)), K.scal(u, K.mul(inv, inv)))
    return a2, b2, *K.char_poly(img)


def kernel_vanishes(F, a, b, s, q):
    """Does psi_5 vanish at the roots of x^2 - s x + q?"""
    K = Kalg(F, s, q)
    X = K.X()
    out = (F.zero, F.zero)
    for c in reversed(psi5(F, a, b)):
        out = K.add(K.mul(out, X), K.const(c))
    return out == (F.zero, F.zero)


def velu5(F, a, b, s, q):
    p2 = F.sub(F.mul(s, s), F.smul(2, q))
    p3 = F.sub(F.mul(s, F.mul(s, s)), F.smul(3, F.mul(s, q)))
    v = F.add(F.smul(6, p2), F.smul(4, a))
    w = F.add(F.add(F.smul(10, p3), F.smul(6, F.mul(a, s))), F.smul(8, b))
    return F.sub(a, F.smul(5, v)), F.sub(b, F.smul(7, w))


def level5_invariant(F, a, b, s, q):
    if s != F.zero:
        si = F.inv(s)
        s2 = F.mul(si, si)
        return ("s", F.mul(a, s2), F.mul(b, F.mul(s2, si)), F.mul(q, s2))
    qi = F.inv(q)
    return ("q", F.mul(a, qi), F.mul(F.mul(b, b), F.mul(qi, F.mul(qi, qi))))


def frob_inv(F, inv):
    return tuple(F.frob(c) if isinstance(c, tuple) else c for c in inv)


def check_sigma(A, sigma):
    n = len(A)
    assert sorted(sigma) == list(range(n))
    assert all(sigma[sigma[v]] == v for v in range(n)), "Frobenius is not an involution"
    for s_ in range(n):
        for t_ in range(n):
            assert A[sigma[s_]][sigma[t_]] == A[s_][t_], "Frobenius is not an automorphism"


def borel5_graph(F, js):
    """(2,3)-isogeny graph on pairs (E, C), C of order 5."""
    verts, reps, orbit = {}, [], []
    for j in js:
        a, b = model(F, j)
        for (s, q) in five_subgroups(F, a, b):
            inv = level5_invariant(F, a, b, s, q)
            if inv in verts:
                orbit[verts[inv]] += 1
                continue
            verts[inv] = len(reps)
            a5, b5 = velu5(F, a, b, s, q)
            reps.append((a, b, s, q, j, j_inv(F, a5, b5)))
            orbit.append(1)
    n = len(reps)
    aut = [Fraction(aut_order(F, r[4]), orbit[k]) for k, r in enumerate(reps)]
    assert all(x.denominator == 1 for x in aut)
    assert sum(1 / x for x in aut) == Fraction(F.p - 1, 4), "mass formula fails (level 5)"
    jset = set(js)
    assert all(r[5] in jset for r in reps), "Velu quotient by C is not supersingular"
    A2 = [[0] * n for _ in range(n)]
    A3 = [[0] * n for _ in range(n)]
    for k, (a, b, s, q, j, _) in enumerate(reps):
        assert kernel_vanishes(F, a, b, s, q)
        for e in two_torsion(F, a, b):
            a2, b2, s2, q2 = push2(F, a, b, e, s, q)
            assert kernel_vanishes(F, a2, b2, s2, q2), "image is not 5-torsion"
            inv = level5_invariant(F, a2, b2, s2, q2)
            assert inv in verts, "unknown level-5 vertex (2-isogeny)"
            A2[k][verts[inv]] += 1
        for x0 in three_torsion(F, a, b):
            a3, b3, s3, q3 = push3(F, a, b, x0, s, q)
            assert kernel_vanishes(F, a3, b3, s3, q3), "image is not 5-torsion"
            inv = level5_invariant(F, a3, b3, s3, q3)
            assert inv in verts, "unknown level-5 vertex (3-isogeny)"
            A3[k][verts[inv]] += 1
    for A, deg in ((A2, 3), (A3, 4)):
        for k in range(n):
            assert sum(A[k]) == deg
            for m in range(n):
                assert A[k][m] * aut[m] == A[m][k] * aut[k], "level-5 Brandt symmetry fails"
    sigma = [verts[frob_inv(F, level5_invariant(F, r[0], r[1], r[2], r[3]))] for r in reps]
    check_sigma(A2, sigma)
    check_sigma(A3, sigma)
    labels = [[F.show(r[4]), F.show(r[5])] for r in reps]
    return A2, A3, labels, sigma


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
    AB, lab_b, sig_b = borel3_graph(F, js)
    lab = [F.show(j) for j in js]
    idx = {j: k for k, j in enumerate(js)}
    sig = [idx[F.frob(j)] for j in js]
    for A in (A2, A3):
        check_sigma(A, sig)
    fp = [sig[k] == k for k in range(len(js))]
    assert fp == [j[1] == 0 for j in js]
    if p > 3:
        assert sum(fp) == n_spine(p), "spine size disagrees with the class number formula"
    fp_b = [sig_b[k] == k for k in range(len(sig_b))]
    A25 = None
    if p != 5 and p <= LEVEL5_MAX_P:
        A25, A35, lab_5, sig_5 = borel5_graph(F, js)
        fp_5 = [sig_5[k] == k for k in range(len(sig_5))]
    lam = {"2": second_eigenvalue([A2]), "23": second_eigenvalue([A2, A3]),
           "2b3": second_eigenvalue([AB]),
           "23b5": second_eigenvalue([A25, A35]) if A25 else None}
    # Ramanujan-type bounds, checked numerically rather than assumed.
    for key, bound in (("2", 2 * 2 ** 0.5), ("23", 2 * 2 ** 0.5 + 2 * 3 ** 0.5),
                       ("2b3", 2 * 2 ** 0.5),
                       ("23b5", 2 * 2 ** 0.5 + 2 * 3 ** 0.5)):
        if lam[key] is not None:
            assert lam[key] <= bound + 1e-6, (p, key, lam[key], bound)
    out = {
        "p": p,
        "nonresidue": F.r if F.r != p - 1 else -1,
        "n_level": len(lab_b),
        "graphs": {
            "2": {"nodes": lab, "fp": fp, "links": links_from(A2, 2), "sigma": sig,
                  "lambda2": lam["2"]},
            "23": {"nodes": lab, "fp": fp,
                   "links": links_from(A2, 2) + links_from(A3, 3), "sigma": sig,
                   "lambda2": lam["23"]},
            "2b3": {"nodes": lab_b, "fp": fp_b, "links": links_from(AB, 2),
                    "sigma": sig_b, "lambda2": lam["2b3"]},
        },
    }
    if A25 is not None:
        out["n_level5"] = len(lab_5)
        out["graphs"]["23b5"] = {
            "nodes": lab_5, "fp": fp_5, "sigma": sig_5,
            "links": links_from(A25, 2) + links_from(A35, 3), "lambda2": lam["23b5"]}
    if p > LEVEL_MAX_P:
        del out["graphs"]["2b3"]
    return out


def is_prime(n: int) -> bool:
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def main():
    random.seed(0)
    primes = sorted({int(a) for a in sys.argv[1:]} or DEFAULT_PRIMES)
    bad = [p for p in primes if not is_prime(p) or p < 5]
    if bad:
        sys.exit(f"not primes >= 5: {bad}")
    out = Path(__file__).resolve().parent.parent / "assets" / "data" / "isogeny"
    out.mkdir(parents=True, exist_ok=True)
    for old in out.glob("p*.json"):      # drop files for primes no longer listed
        if old.stem[1:].isdigit():
            old.unlink()
    total = 0
    for p in primes:
        res = compute(p)
        g = res["graphs"]
        print(f"p = {p:5d}: {len(g['2']['nodes']):4d} curves, "
              f"{res['n_level']:4d} level-structure vertices, "
              f"lambda2 = {g['2']['lambda2']}, {g['23']['lambda2']}")
        f = out / f"p{p}.json"
        f.write_text(json.dumps(res, separators=(",", ":")) + "\n")
        total += f.stat().st_size
    (out / "index.json").write_text(json.dumps(
        {"primes": primes, "level_max_p": LEVEL_MAX_P, "level5_max_p": LEVEL5_MAX_P}, separators=(",", ":")) + "\n")
    print(f"wrote {len(primes)} files to {out} ({total} bytes in all)")


if __name__ == "__main__":
    main()
