#!/usr/bin/env python3
"""Compute the supersingular isogeny graphs drawn on the home page.

Writes assets/data/isogeny/: one file p<prime>-<graph>.json per prime and graph (the
page loads only the one you pick), an index.json listing them, and spectra.json. Pure Python, no dependencies
(numpy is needed for the eigenvalues).

    python3 _scripts/isogeny_graphs.py            # the default list of primes
    python3 _scripts/isogeny_graphs.py 431 1009   # your own primes (any primes >= 5)
    python3 _scripts/isogeny_graphs.py --skip-spectra   # leave spectra.json alone

Five graphs are computed for each prime p:

  "2"     the supersingular 2-isogeny graph G_p(2)
  "23"    G_p({2,3}): the same vertices, 2- and 3-isogenies
  "235"   G_p({2,3,5}): the same vertices, 2-, 3- and 5-isogenies
  "2b5"   G_p(2) with Borel (Gamma_0(5)) level structure: vertices are pairs (E, C),
          C < E of order 5, up to isomorphism; a 2-isogeny phi: E -> E' gives the
          edge (E, C) -> (E', phi(C))
  "23b5"  G_p({2,3}) with Borel level-5 structure (2- and 3-isogenies)

and, in spectra.json, the eigenvalues of the Hecke operators T_2, T_3, T_5 on
S_2(Gamma_0(p)) for a longer list of primes (the Brandt matrices of G_p(l) without
their eigenvalue l + 1), for the spectral-law animation.

Every graph also carries
  "sigma"  the permutation of the vertices induced by Frobenius (E, C) -> (E^(p), C^(p));
  "fp"     whether a vertex is fixed by Frobenius;
  "jfp"    whether its j-invariant (of E) lies in F_p (informational: without level
           structure jfp = fp; with level structure fp is the smaller set).
The F_p-subgraph (also called the spine) is the subgraph induced on the vertices fixed
by Frobenius ("fp"), which acts on the level structure as well.
  "fold"   the links of the quotient by Frobenius: vertices are
           Frobenius-orbits of vertices (named by the smaller index), edges are orbits
           of edges. This needs Frobenius on edges, computed and checked below.

Method. Everything is done with short Weierstrass models y^2 = x^3 + ax + b
over F_{p^2} and Velu's formulas, so only F_{p^2} arithmetic and root finding
for polynomials of degree <= 4 are needed.

  * One supersingular curve is found by reducing a class-number-one CM curve
    at a prime that is inert or ramified in its CM field; the rest are found
    by walking the (connected) 2-isogeny graph.
  * A subgroup of order 3 is {O, +-P}, determined by x(P), a root of the
    3-division polynomial psi_3 = 3x^4 + 6ax^2 + 12bx - a^2. A subgroup of order 5
    is determined by its kernel polynomial x^2 - s x + q (see below).
  * Edges are *directed* and counted by kernels, so every vertex has
    out-degree the sum of ell + 1 over the ells used. At vertices with
    extra automorphisms (j = 0, 1728 and their lifts) the adjacency matrix is
    not symmetric; the drawing shows max(A_ij, A_ji) undirected edges.

Checks run on every prime (the script aborts if one fails):
  * number of supersingular j-invariants agrees with the class number formula;
  * every Velu output is a root of the classical modular polynomial
    Phi_2 or Phi_3; every 5-isogeny lands on a supersingular curve;
  * A_ij * |Aut_j| = A_ji * |Aut_i| (the Brandt-matrix symmetry);
  * number of F_p-rational vertices (the spine) agrees with the class number
    formula (h(-4p)/2, h(-p), or 2 h(-p) according to p mod 8);
  * Frobenius is an involution on vertices and on edges and respects the source
    and target of every edge;
  * for level 5: the pushed-forward kernel polynomial is a factor of psi_5, every
    quotient curve is supersingular, and the mass is (p - 1)/4;
  * the second-largest eigenvalue modulus is within the Ramanujan bound
    (the sum of 2 sqrt(l) over the l used), and in spectra.json every eigenvalue
    of T_l on S_2(Gamma_0(p)) is at most 2 sqrt(l) in absolute value.
"""
from __future__ import annotations

import json
import random
import sys
from fractions import Fraction
from pathlib import Path

DEFAULT_PRIMES = [101, 157, 211, 307, 431, 599, 827, 1009, 1499, 1999, 2503,
                  2999, 3511, 4001, 4507, 4999, 6007, 8009, 10007, 15013, 20011]
# Largest p for each graph (recorded in index.json, which the page reads). The graphs
# with level structure have 6 (Gamma_0(5)) or 12 (Gamma_1(5)) times as many vertices as
# G_p(2); the full level-2 graph is a small example.
CAPS = {"2": None, "23": None, "235": None, "2b5": 5000, "23b5": 5000, "2g15": 5000, "3g2": 431}
LEVEL5_MAX_P = 5000
# Primes for the spectral-law animation (spectra.json): roughly geometric up to 20000.
SPECTRA_PRIMES = [101, 131, 173, 227, 293, 383, 499, 653, 853, 1109, 1447, 1889, 2459,
                  3203, 4177, 5443, 7079, 9221, 12007, 15641, 20369]

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


def arcs_from(A, ell):
    """Directed links [s, t, ell, A[s][t]] (for directed graphs)."""
    n = len(A)
    return [[s_, t_, ell, A[s_][t_]] for s_ in range(n) for t_ in range(n) if A[s_][t_]]


def check_sigma(A, sigma):
    """Frobenius is a permutation of the vertices of order dividing 4 (order 2, except
    with Gamma_1(5) structure) that preserves the adjacency matrix."""
    n = len(A)
    assert sorted(sigma) == list(range(n))
    assert all(sigma[sigma[sigma[sigma[v]]]] == v for v in range(n)), "Frobenius^4 is not the identity"
    for s_ in range(n):
        for t_ in range(n):
            assert A[sigma[s_]][sigma[t_]] == A[s_][t_], "Frobenius is not an automorphism"


def orbits(perm):
    """The cycles of a permutation, each starting at its smallest element."""
    seen, out = set(), []
    for v in range(len(perm)):
        if v in seen:
            continue
        cyc, w = [], v
        while w not in seen:
            seen.add(w)
            cyc.append(w)
            w = perm[w]
        out.append(cyc)
    return out


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
    return A2, A3, labels, sigma, reps


# --------------------------------------------------------------------------
# Frobenius on edges, and the folded graph.
#
# An edge out of a vertex is a kernel on that vertex's fixed model: a 2-torsion
# x-coordinate e (2-isogeny) or the x-coordinate x0 of a 3-torsion point (3-isogeny).
# Frobenius sends the edge (v, K) to (sigma(v), K^(p)). The fixed model of sigma(v) is
# the conjugate of the model of v (model() has coefficients in F_p(j)), so on the
# 2- and (2,3)-graphs K^(p) is just the conjugate kernel. With level structure the
# stored representative of sigma(v) can differ from the conjugate representative by
# an automorphism of the model, x -> c x with c = u^2 (c = -1 needs j = 1728,
# c a cube root of unity needs j = 0); the kernel is moved by the same c.
# The folded graph is the quotient by Frobenius: its vertices are the orbits of
# vertices, its directed edges the orbits of directed edges. It is drawn, like the
# others, with max(B_ij, B_ji) undirected edges between two vertices.
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# F_{p^4} = F_{p^2}[t] / (t^2 - c), c a non-square in F_{p^2}: needed for points of
# order 5, whose x-coordinates need not lie in F_{p^2}. Elements: pairs (u, v) = u + v t.
# --------------------------------------------------------------------------
class Fp4:
    def __init__(self, F):
        self.F = F
        rng = random.Random(1)
        while True:
            c = (rng.randrange(F.p), rng.randrange(F.p))
            if c != F.zero and F.pow(c, (F.q - 1) // 2) != F.one:
                break
        self.c = c
        self.zero = (F.zero, F.zero)
        self.one = (F.one, F.zero)

    def emb(self, x):
        return (x, self.F.zero)

    def add(self, x, y):
        F = self.F
        return (F.add(x[0], y[0]), F.add(x[1], y[1]))

    def sub(self, x, y):
        F = self.F
        return (F.sub(x[0], y[0]), F.sub(x[1], y[1]))

    def scal(self, k, x):                  # k in F_{p^2}
        F = self.F
        return (F.mul(k, x[0]), F.mul(k, x[1]))

    def mul(self, x, y):
        F = self.F
        return (F.add(F.mul(x[0], y[0]), F.mul(self.c, F.mul(x[1], y[1]))),
                F.add(F.mul(x[0], y[1]), F.mul(x[1], y[0])))

    def inv(self, x):
        F = self.F
        nrm = F.sub(F.mul(x[0], x[0]), F.mul(self.c, F.mul(x[1], x[1])))
        ni = F.inv(nrm)
        return (F.mul(x[0], ni), F.neg(F.mul(x[1], ni)))

    def pow(self, x, n):
        out = self.one
        while n:
            if n & 1:
                out = self.mul(out, x)
            x = self.mul(x, x)
            n >>= 1
        return out

    def frob(self, x):                      # x -> x^p
        return self.pow(x, self.F.p)

    def quadratic_roots(self, s, q):
        """The two roots of x^2 - s x + q (s, q in F_{p^2}) in F_{p^4}."""
        F = self.F
        rs = roots(F, [q, F.neg(s), F.one])
        if len(rs) == 2:
            out = [self.emb(r) for r in rs]
        else:
            assert not rs, "repeated root"
            disc = F.sub(F.mul(s, s), F.smul(4, q))
            w = roots(F, [F.neg(F.div(disc, self.c)), F.zero, F.one])
            assert len(w) == 2
            half = F.inv(F.el(2))
            out = [(F.mul(s, half), F.mul(w[0], half)), (F.mul(s, half), F.neg(F.mul(w[0], half)))]
        for x in out:
            assert self.add(self.sub(self.mul(x, x), self.scal(s, x)), self.emb(q)) == self.zero
        return out


def gamma1_invariant(F, F4, a, b, x0, j):
    """A complete invariant of (E, +-P) with x(P) = x0: isomorphisms act by
    (a, b, x) -> (u^4 a, u^6 b, u^2 x)."""
    if x0 == F4.zero:
        return ("z", j)
    xi = F4.inv(x0)
    xi2 = F4.mul(xi, xi)
    return ("x", F4.scal(a, xi2), F4.scal(b, F4.mul(xi2, xi)))


def gamma1_graph(F, js):
    """G_p(2) with Gamma_1(5) structure: vertices (E, +-P), P of order 5, up to
    isomorphism; a 2-isogeny phi gives (E, +-P) -> (E', +-phi(P)).

    This graph is directed: the dual isogeny sends +-phi(P) to +-2P, not +-P. The
    Brandt symmetry becomes A[v][w] |Aut w| = A[w][<2>v] |Aut v|, where the diamond
    operator <2> sends (E, +-P) to (E, +-2P)."""
    F4 = Fp4(F)
    verts, reps, orbit, other = {}, [], [], []
    for j in js:
        a, b = model(F, j)
        for s_, q_ in five_subgroups(F, a, b):
            xs = F4.quadratic_roots(s_, q_)          # x(P) and x(2P): two classes +-P
            for x0, x2 in (xs, xs[::-1]):
                inv = gamma1_invariant(F, F4, a, b, x0, j)
                if inv in verts:
                    orbit[verts[inv]] += 1
                    continue
                verts[inv] = len(reps)
                a5, b5 = velu5(F, a, b, s_, q_)
                reps.append((a, b, x0, j, j_inv(F, a5, b5)))
                other.append(gamma1_invariant(F, F4, a, b, x2, j))
                orbit.append(1)
    n = len(reps)
    aut = [Fraction(aut_order(F, r[3]), orbit[k]) for k, r in enumerate(reps)]
    assert all(x.denominator == 1 for x in aut)
    assert sum(1 / x for x in aut) == Fraction(F.p - 1, 2), "mass formula fails (Gamma_1(5))"
    edges = []
    for v, (a, b, x0, j, _) in enumerate(reps):
        for e in two_torsion(F, a, b):
            a2, b2, _ = velu2(F, a, b, e)
            t = F.add(F.smul(3, F.mul(e, e)), a)
            x1 = F4.add(x0, F4.scal(t, F4.inv(F4.sub(x0, F4.emb(e)))))
            inv = gamma1_invariant(F, F4, a2, b2, x1, j_inv(F, a2, b2))
            assert inv in verts, "unknown Gamma_1(5) vertex"
            edges.append((v, verts[inv], 2, e))
    A = [[0] * n for _ in range(n)]
    for s_, t_, _, _ in edges:
        A[s_][t_] += 1
    dia = [verts[o] for o in other]          # the diamond operator <2>
    for v in range(n):
        assert sum(A[v]) == 3
        assert dia[dia[v]] == v, "<2> is not an involution on +-P classes"
        for w in range(n):
            assert A[v][w] * aut[w] == A[w][dia[v]] * aut[v], "Gamma_1(5) Brandt symmetry fails"

    def conj_inv(inv):
        return ("z", F.frob(inv[1])) if inv[0] == "z" else ("x", F4.frob(inv[1]), F4.frob(inv[2]))
    sig = [verts[conj_inv(gamma1_invariant(F, F4, r[0], r[1], r[2], r[3]))] for r in reps]
    check_sigma(A, sig)
    # Frobenius squared is the diamond operator <p>: the identity if p = +-1 (mod 5),
    # and <2> if p = +-2 (mod 5).
    sq = dia if F.p % 5 in (2, 3) else list(range(n))
    assert all(sig[sig[v]] == sq[v] for v in range(n)), "Frobenius^2 is not <p>"
    cs = []                                 # automorphism scalings, as in level_sigma
    for v, r in enumerate(reps):
        r2 = reps[sig[v]]
        fx = F4.frob(r[2])
        for c in scalings(F, r2[0], r2[1]):
            if F4.scal(c, fx) == r2[2]:
                cs.append(c)
                break
        else:
            raise AssertionError("no automorphism matches the conjugate point")
    se = edge_sigma(F, edges, sig, lambda v, key: F.mul(cs[v], F.frob(key)))
    labels = [[F.show(r[3]), F.show(r[4])] for r in reps]
    return A, labels, sig, reps, edges, se


def full2_graph(F, js):
    """G_p(3) with full level-2 structure: vertices (E, (P, Q)), (P, Q) an ordered basis
    of E[2], up to isomorphism. The Legendre invariant lambda = (e3 - e1)/(e2 - e1)
    (e_i the x-coordinates of P, Q, P + Q) is a complete invariant, so the vertices are
    the supersingular lambda, the (p - 1)/2 roots of the Hasse polynomial. A 3-isogeny
    phi gives (E, (P, Q)) -> (E', (phi P, phi Q))."""
    verts, reps, orbit = {}, [], []
    for j in js:
        a, b = model(F, j)
        es = two_torsion(F, a, b)
        assert len(es) == 3
        for e1 in es:
            for e2 in es:
                if e2 == e1:
                    continue
                e3 = next(e for e in es if e not in (e1, e2))
                lam = F.div(F.sub(e3, e1), F.sub(e2, e1))
                if lam in verts:
                    orbit[verts[lam]] += 1
                    continue
                verts[lam] = len(reps)
                reps.append((a, b, e1, e2, j, lam))
                orbit.append(1)
    n = len(reps)
    assert n == (F.p - 1) // 2, "number of supersingular lambda is not (p - 1)/2"
    aut = [Fraction(aut_order(F, r[4]), orbit[k]) for k, r in enumerate(reps)]
    assert sum(1 / x for x in aut) == Fraction(F.p - 1, 4), "mass formula fails (Gamma(2))"
    edges = []
    for v, (a, b, e1, e2, j, _) in enumerate(reps):
        es = two_torsion(F, a, b)
        e3 = next(e for e in es if e not in (e1, e2))
        for x0 in three_torsion(F, a, b):
            a3, b3 = velu3(F, a, b, x0)
            t = F.smul(2, F.add(F.smul(3, F.mul(x0, x0)), a))
            u = F.smul(4, F.add(F.mul(x0, F.add(F.mul(x0, x0), a)), b))
            push = lambda x: F.add(F.add(x, F.div(t, F.sub(x, x0))),
                                   F.div(u, F.mul(F.sub(x, x0), F.sub(x, x0))))
            f1, f2, f3 = push(e1), push(e2), push(e3)
            assert set(two_torsion(F, a3, b3)) == {f1, f2, f3}, "image is not the 2-torsion"
            lam = F.div(F.sub(f3, f1), F.sub(f2, f1))
            edges.append((v, verts[lam], 3, x0))
    A = [[0] * n for _ in range(n)]
    for s_, t_, _, _ in edges:
        A[s_][t_] += 1
    for v in range(n):
        assert sum(A[v]) == 4
        for w in range(n):
            assert A[v][w] * aut[w] == A[w][v] * aut[v], "Gamma(2) Brandt symmetry fails"
    sig = [verts[F.frob(r[5])] for r in reps]
    check_sigma(A, sig)
    se = level_sigma(F, edges, reps, sig, lambda r: [(r[2], 1), (r[3], 1)], {3: A})
    labels = [[F.show(r[4]), F.show(r[5])] for r in reps]
    return A, labels, sig, reps, edges, se


def scalings(F, a, b):
    """The values c = u^2 for automorphisms (x, y) -> (u^2 x, u^3 y) of y^2 = x^3 + ax + b."""
    out = [F.one]
    if b == F.zero:
        out.append(F.neg(F.one))
    if a == F.zero:
        out += roots(F, [F.one, F.one, F.one])          # c^2 + c + 1 = 0
    return out


def edge_sigma(F, edges, sig_v, transport):
    """edges: list of (src, dst, ell, key). transport(src, key) -> key on sigma(src).
    Returns the permutation of edges induced by Frobenius, after checking it."""
    where = {(s, ell, key): k for k, (s, _, ell, key) in enumerate(edges)}
    sig_e = []
    for s, t, ell, key in edges:
        k2 = where.get((sig_v[s], ell, transport(s, key)))
        assert k2 is not None, "Frobenius image of an edge not found"
        assert edges[k2][1] == sig_v[t], "Frobenius does not respect the edge's target"
        sig_e.append(k2)
    assert all(sig_e[sig_e[sig_e[sig_e[k]]]] == k for k in range(len(edges))), \
        "Frobenius^4 is not the identity on edges"
    assert all(sig_e[sig_e[k]] == k for k in range(len(edges))) == \
        all(sig_v[sig_v[v]] == v for v in range(len(sig_v))), "Frobenius^2 differs on edges and vertices"
    return sig_e


def fold_links(edges, sig_v, sig_e, directed=False):
    """Undirected links of the quotient graph, with vertices named by the smaller
    member of each Frobenius orbit."""
    name = {}
    for cyc in orbits(sig_v):
        for v in cyc:
            name[v] = cyc[0]
    B = {}
    for cyc in orbits(sig_e):
        s, t, ell, _ = edges[cyc[0]]
        key = (name[s], name[t], ell)
        B[key] = B.get(key, 0) + 1
    if directed:
        return [[a, b, ell, m] for (a, b, ell), m in sorted(B.items())]
    out = []
    for a, b, ell in sorted({(min(a, b), max(a, b), ell) for a, b, ell in B}):
        m = B.get((a, b, ell), 0) if a == b else max(B.get((a, b, ell), 0), B.get((b, a, ell), 0))
        out.append([a, b, ell, m])
    return out


def base_edges(F, js, ell):
    idx = {j: k for k, j in enumerate(js)}
    edges = []
    for s, j in enumerate(js):
        a, b = model(F, j)
        if ell == 2:
            for e in two_torsion(F, a, b):
                a2, b2, _ = velu2(F, a, b, e)
                edges.append((s, idx[j_inv(F, a2, b2)], 2, e))
        elif ell == 3:
            for x0 in three_torsion(F, a, b):
                a2, b2 = velu3(F, a, b, x0)
                edges.append((s, idx[j_inv(F, a2, b2)], 3, x0))
        else:                  # ell = 5: the kernel is given by x^2 - s x + q
            for sq in five_subgroups(F, a, b):
                a2, b2 = velu5(F, a, b, *sq)
                j2 = j_inv(F, a2, b2)
                assert j2 in idx, "a 5-isogeny leaves the supersingular locus"
                edges.append((s, idx[j2], 5, sq))
    sig_v = [idx[F.frob(j)] for j in js]
    for s, j in enumerate(js):
        assert model(F, F.frob(j)) == tuple(F.frob(c) for c in model(F, j))
    conj = (lambda key: tuple(F.frob(c) for c in key)) if ell == 5 else F.frob
    sig_e = edge_sigma(F, edges, sig_v, lambda s, key: conj(key))
    return edges, sig_v, sig_e


def adjacency_from(edges, n, ell, aut):
    """Adjacency matrix of the ell-edges, with the Brandt-matrix symmetry checked."""
    A = [[0] * n for _ in range(n)]
    for s, t, e_ell, _ in edges:
        if e_ell == ell:
            A[s][t] += 1
    for i in range(n):
        assert sum(A[i]) == ell + 1
        for k in range(n):
            assert A[i][k] * aut[k] == A[k][i] * aut[i], "Brandt symmetry fails"
    return A


def join_edges(*parts):
    """Concatenate (edges, sigma on edges) pairs."""
    edges, sig = [], []
    for e, se in parts:
        off = len(edges)
        edges += e
        sig += [k + off for k in se]
    return edges, sig


def keep_ells(edges, se, ells):
    """The edges with ell in ells, and Frobenius on them (it preserves ell)."""
    keep = [k for k, e in enumerate(edges) if e[2] in ells]
    new = {k: i for i, k in enumerate(keep)}
    return [edges[k] for k in keep], [new[se[k]] for k in keep]


def level_sigma(F, edges, reps, sig_v, kernel_of, A_by_ell):
    """Frobenius on the edges of a level-structure graph. kernel_of(rep) lists the
    level data with their weights under x -> c x: [(x0, 1)] or [(s, 1), (q, 2)]."""
    for ell, A in A_by_ell.items():          # the edge list must match the matrices
        n = len(A)
        cnt = [[0] * n for _ in range(n)]
        for s, t, e_ell, _ in edges:
            if e_ell == ell:
                cnt[s][t] += 1
        assert cnt == A, "edge list disagrees with the adjacency matrix"
    cs = []
    for v, r in enumerate(reps):
        r2 = reps[sig_v[v]]
        assert (F.frob(r[0]), F.frob(r[1])) == (r2[0], r2[1])
        for c in scalings(F, r2[0], r2[1]):
            if all(F.mul(F.pow(c, w), F.frob(x)) == y
                   for (x, w), (y, _) in zip(kernel_of(r), kernel_of(r2))):
                cs.append(c)
                break
        else:
            raise AssertionError("no automorphism matches the conjugate level structure")
    return edge_sigma(F, edges, sig_v, lambda v, key: F.mul(cs[v], F.frob(key)))


def level5_edges(F, reps, sig_v, A2, A3):
    verts = {level5_invariant(F, r[0], r[1], r[2], r[3]): k for k, r in enumerate(reps)}
    edges = []
    for v, (a, b, s_, q_, _, _) in enumerate(reps):
        for e in two_torsion(F, a, b):
            a2, b2, s2, q2 = push2(F, a, b, e, s_, q_)
            edges.append((v, verts[level5_invariant(F, a2, b2, s2, q2)], 2, e))
        for x0 in three_torsion(F, a, b):
            a3, b3, s3, q3 = push3(F, a, b, x0, s_, q_)
            edges.append((v, verts[level5_invariant(F, a3, b3, s3, q3)], 3, x0))
    return edges, level_sigma(F, edges, reps, sig_v, lambda r: [(r[2], 1), (r[3], 2)],
                              {2: A2, 3: A3})


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


RAMANUJAN = {2: 2 * 2 ** 0.5, 3: 2 * 3 ** 0.5, 5: 2 * 5 ** 0.5}


def compute(p: int, level5: bool = True):
    F = Fp2(p)
    js = supersingular_js(F)
    n = len(js)
    lab = [F.show(j) for j in js]
    aut = [aut_order(F, j) for j in js]
    E = {ell: base_edges(F, js, ell) for ell in (2, 3, 5)}
    sig = E[2][1]
    A = {ell: adjacency_from(E[ell][0], n, ell, aut) for ell in (2, 3, 5)}
    assert A[2] == isogeny_adjacency(F, js, 2) and A[3] == isogeny_adjacency(F, js, 3)
    for M in A.values():
        check_sigma(M, sig)
    fp = [sig[k] == k for k in range(n)]
    assert fp == [j[1] == 0 for j in js]
    assert sum(fp) == n_spine(p), "number of j in F_p disagrees with the class number formula"
    graphs = {}

    def add(key, nodes, fixed, jfp, sigma, mats, edges, sig_e, directed=False):
        lam = second_eigenvalue([M for _, M in mats])
        bound = sum(RAMANUJAN[ell] for ell, _ in mats)
        if lam is not None:
            assert lam <= bound + 1e-6, (p, key, lam, bound)
        lf = arcs_from if directed else links_from
        graphs[key] = {
            "nodes": nodes, "fp": fixed, "jfp": jfp, "sigma": sigma,
            "links": [lk for ell, M in mats for lk in lf(M, ell)],
            "fold": fold_links(edges, sigma, sig_e, directed), "lambda2": lam}
        if directed:
            graphs[key]["directed"] = True

    for key, ells in (("2", (2,)), ("23", (2, 3)), ("235", (2, 3, 5))):
        edges, se = join_edges(*((E[ell][0], E[ell][2]) for ell in ells))
        add(key, lab, fp, fp, sig, [(ell, A[ell]) for ell in ells], edges, se)
    out = {"p": p, "nonresidue": F.r if F.r != p - 1 else -1, "graphs": graphs}
    jset_fp = {j for j in js if j[1] == 0}
    if level5 and p <= CAPS["2g15"]:
        A1, lab_1, sig_1, reps_1, e1, se1 = gamma1_graph(F, js)
        add("2g15", lab_1, [sig_1[k] == k for k in range(len(sig_1))],
            [r[3] in jset_fp for r in reps_1], sig_1, [(2, A1)], e1, se1, directed=True)
    if level5 and p <= CAPS["3g2"]:
        Ag, lab_g, sig_g, reps_g, eg, seg = full2_graph(F, js)
        add("3g2", lab_g, [sig_g[k] == k for k in range(len(sig_g))],
            [r[4] in jset_fp for r in reps_g], sig_g, [(3, Ag)], eg, seg)
    if level5 and p <= LEVEL5_MAX_P:
        A25, A35, lab_5, sig_5, reps_5 = borel5_graph(F, js)
        fp_5 = [sig_5[k] == k for k in range(len(sig_5))]
        jfp_5 = [r[4] in jset_fp for r in reps_5]
        e5, se5 = level5_edges(F, reps_5, sig_5, A25, A35)
        add("2b5", lab_5, fp_5, jfp_5, sig_5, [(2, A25)], *keep_ells(e5, se5, {2}))
        add("23b5", lab_5, fp_5, jfp_5, sig_5, [(2, A25), (3, A35)], e5, se5)
        out["n_level5"] = len(lab_5)
    return out


def hecke_spectra(p: int):
    """Eigenvalues of T_l on S_2(Gamma_0(p)) for l = 2, 3, 5: the eigenvalues of the
    Brandt matrix (the adjacency matrix of G_p(l)) with the eigenvalue l + 1 of the
    Eisenstein series removed. A D is symmetric (D = diag |Aut|), so
    D^(-1/2) (A D) D^(-1/2) is a symmetric matrix similar to A."""
    import numpy as np
    F = Fp2(p)
    js = supersingular_js(F)
    n = len(js)
    aut = [aut_order(F, j) for j in js]
    d = np.sqrt(np.array(aut, dtype=float))
    out = {}
    for ell in (2, 3, 5):
        A = np.array(adjacency_from(base_edges(F, js, ell)[0], n, ell, aut), dtype=float)
        S = (A * np.array(aut, dtype=float)[None, :]) / d[:, None] / d[None, :]
        assert np.allclose(S, S.T)
        ev = np.sort(np.linalg.eigvalsh(S))
        assert abs(ev[-1] - (ell + 1)) < 1e-8, "the top eigenvalue is not l + 1"
        cusp = ev[:-1]
        assert np.all(np.abs(cusp) <= 2 * ell ** 0.5 + 1e-8), "an eigenvalue exceeds 2 sqrt(l)"
        out[str(ell)] = [round(float(x), 4) for x in cusp]
    return {"p": p, "dim": n - 1, "eigenvalues": out}


def is_prime(n: int) -> bool:
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def main():
    random.seed(0)
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    primes = sorted({int(a) for a in args} or DEFAULT_PRIMES)
    bad = [p for p in primes if not is_prime(p) or p < 5]
    if bad:
        sys.exit(f"not primes >= 5: {bad}")
    out = Path(__file__).resolve().parent.parent / "assets" / "data" / "isogeny"
    out.mkdir(parents=True, exist_ok=True)
    for old in out.glob("p*.json"):      # drop old files
        old.unlink()
    total = 0
    for p in primes:
        res = compute(p)
        g = res["graphs"]
        print(f"p = {p:5d}: {len(g['2']['nodes']):4d} curves, "
              f"{res.get('n_level5', '-'):>4} level-5 vertices, "
              f"lambda2 = {g['2']['lambda2']}, {g['23']['lambda2']}, {g['235']['lambda2']}")
        total += write_prime(out, res)
    write_index(out, primes)
    print(f"wrote {len(primes)} files to {out} ({total} bytes in all)")
    if "--skip-spectra" not in sys.argv:
        spectra = []
        for p in SPECTRA_PRIMES:
            spectra.append(hecke_spectra(p))
            print(f"spectra: p = {p}, dim S_2 = {spectra[-1]['dim']}")
        (out / "spectra.json").write_text(json.dumps(spectra, separators=(",", ":")) + "\n")


def write_prime(out: Path, res) -> int:
    """One file per graph: p<prime>-<graph>.json. Returns the bytes written."""
    total = 0
    for kind, g in res["graphs"].items():
        f = out / f"p{res['p']}-{kind}.json"
        f.write_text(json.dumps({"p": res["p"], "nonresidue": res["nonresidue"], "graph": g},
                                separators=(",", ":")) + "\n")
        total += f.stat().st_size
    return total


def write_index(out: Path, primes) -> None:
    """index.json: the primes on offer, the caps, and the size of each F_p-subgraph
    (the number of vertices fixed by Frobenius)."""
    spine = {}
    for p in primes:
        spine[str(p)] = {}
        for f in out.glob(f"p{p}-*.json"):
            spine[str(p)][f.stem.split("-", 1)[1]] = sum(json.loads(f.read_text())["graph"]["fp"])
    (out / "index.json").write_text(json.dumps(
        {"primes": primes, "caps": CAPS, "spine": spine},
        separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
