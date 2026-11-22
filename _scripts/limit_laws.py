#!/usr/bin/env python3
"""Compute the data for the Limit laws page (animations/limits.qmd).

Writes assets/data/limits/:

  gamma1.json    T_2 on S_2(Gamma_1(101)), split by the even characters chi mod 101:
                 for each chi its value chi(2) and the eigenvalues of T_2 on S_2(101, chi).
  quadratic.json T_2 on S_2(N, (./N)) for primes N = 5 (mod 8), where (2/N) = -1 and the
                 eigenvalues are purely imaginary: histograms of Im(lambda)/sqrt 2.
  planar.json    T_2 on all of S_2(Gamma_1(N)) for growing primes N: the normalised
                 eigenvalues lambda/sqrt 2 as points of the complex plane (a fixed random
                 sample of 3000 of them when there are more), with statistics of all.
  diagonal.json  T_l on S_2(Gamma_0(p)) with l growing with p (l the largest prime
                 <= (log p)^2), from the Brandt matrix of G_p(l): histograms of lambda/sqrt l.
  weight.json    T_l on S_k(SL_2(Z)) for l = 2, 3, 5 and growing weight k: histograms of
                 lambda / l^((k-1)/2).

Needs Python 3, numpy and cypari2 (pip install cypari2). Run from the repository root:

    python3 _scripts/limit_laws.py              # everything (under an hour on one core)
    python3 _scripts/limit_laws.py planar       # only some of the files
    PARI_STACK=6000000000 python3 _scripts/limit_laws.py planar   # planar needs ~6 GB

Method.
  * Spaces with character (gamma1, quadratic, planar): PARI's mfinit([N, 2, chi], 1) and
    mfheckemat, embedded into C by mfembed. For the Conrey label a, PARI's chi is the
    character chareval(znstar(N, 1), a, .) and the embedding is the one with those values;
    the check lambda = chi(2) * conj(lambda) below confirms that the two agree.
    For planar.json, one Galois orbit of characters at a time (gamma1_by_orbits): T_2 over
    Q(chi) = Q[t]/Phi(t), with t replaced by each complex root of Phi; checked against the
    one-character-at-a-time method at N = 101.
  * Brandt matrices (diagonal): the supersingular j-invariants in F_{p^2}, found by a
    breadth-first search with the modular polynomial Phi_2 from one CM j-invariant, and
    B(l)_{ij} = the multiplicity of j_j as a root of Phi_l(j_i, Y) (PARI's polmodular),
    which counts the subgroups K of order l of E_i with E_i/K = E_j.
  * Level one (weight): the basis Delta^i E_4^a E_6^b (i = 1, ..., d; 4a + 6b = k - 12 i)
    of S_k(SL_2(Z)), whose q-expansions start at q^i, and T_l on q-expansions,
    a_n(T_l f) = a_{ln}(f) + l^(k-1) a_{n/l}(f). The characteristic polynomial is exact
    and its roots are found with enough digits.

Checks (the script stops if one fails):
  * dim S_2(Gamma_1(N)) = (N - 5)(N - 7)/24, the genus of X_1(N), for primes N >= 5;
  * every eigenvalue satisfies |lambda| <= 2 sqrt(l) (Eichler-Shimura, Deligne);
  * lambda = chi(l) * conj(lambda) for every eigenvalue on S_2(N, chi) (since
    T_l^* = conj(chi(l)) T_l for the Petersson product), and the eigenvalues for conj(chi)
    are the conjugates of those for chi;
  * the quadratic character has chi(2) = -1, so the eigenvalues are purely imaginary;
  * Brandt matrices: rows of B(l) sum to l + 1, B(l) diag(|Aut(E_j)|) is symmetric,
    B(l) commutes with B(2), and the number of j-invariants is the Eichler mass count;
  * level one: the exact characteristic polynomial agrees with PARI's mfheckemat for
    k <= 120, and T_2 on S_12 has the eigenvalue tau(2) = -24.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent.parent / "assets" / "data" / "limits"

GAMMA1_N = 101
QUADRATIC_NS = [101, 149, 173, 197, 269, 317, 389, 461, 557, 653, 773, 941, 1109,
                1301, 1549, 1861, 2213, 2621, 3109, 3701]
# Beyond 251: primes N for which (N - 1)/2 has no large prime factor, so that the fields
# Q(chi) stay small enough for PARI. N = 503 ((N - 1)/2 = 251 is prime) does not finish;
# N = 769 (Q(chi) of degree 128, dim S_2(769, chi) = 64) needs more than 6.5 GB. N = 661
# needs about 6 GB: run with PARI_STACK=6000000000.
PLANAR_NS = [29, 37, 53, 61, 73, 89, 101, 113, 131, 151, 173, 197, 223, 251,
             337, 421, 433, 577, 661]
DIAGONAL_PS = None       # filled in by diagonal_primes()
WEIGHT_KS = list(range(24, 961, 24))
WEIGHT_ELLS = (2, 3, 5)
BINS = 40                # histogram bins on [-2, 2]

pari = None


def start_pari():
    global pari
    import cypari2
    pari = cypari2.Pari()
    # PARI's stack: 3 GB by default; the largest planar levels need more, e.g.
    # PARI_STACK=6000000000 python3 _scripts/limit_laws.py planar
    pari.allocatemem(int(os.environ.get("PARI_STACK", 3 * 10 ** 9)))
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".gp", delete=False) as f:
        f.write(BRANDT_GP + LEVEL1_GP)
    pari(f'read("{f.name}")')
    Path(f.name).unlink()


# ---- laws ------------------------------------------------------------------------

def km(ell, x):
    """Density of the normalised Kesten-McKay law mu_l (of lambda/sqrt l) on [-2, 2]."""
    x = np.asarray(x, dtype=float)
    st = np.sqrt(np.maximum(4 - x * x, 0)) / (2 * np.pi)
    if ell == math.inf:
        return st
    return st * (ell + 1) * ell / ((ell + 1) ** 2 - ell * x * x)


_GRID = np.linspace(-2, 2, 40001)


def cdf(ell):
    f = km(ell, _GRID)
    c = np.concatenate([[0.0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(_GRID))])
    return c / c[-1]


def ks(values, ell):
    """Kolmogorov distance between the empirical law of values (in [-2, 2]) and mu_l."""
    v = np.sort(np.asarray(values, dtype=float))
    m = len(v)
    F = np.interp(v, _GRID, cdf(ell))
    i = np.arange(m)
    return float(max(np.max(np.abs((i + 1) / m - F)), np.max(np.abs(i / m - F))))


def hist(values):
    counts, _ = np.histogram(np.clip(values, -2, 2), bins=BINS, range=(-2, 2))
    return [int(c) for c in counts]


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)


# ---- spaces with character ------------------------------------------------------

def even_labels(N):
    return [a for a in range(1, N) if math.gcd(a, N) == 1
            and not int(pari(f"zncharisodd(znstar({N},1),{a})"))]


def chi_value(N, a, n):
    """chi_a(n) as a complex number, chi_a the character with Conrey label a mod N."""
    G = pari(f"znstar({N},1)")
    o = int(pari.charorder(G, pari(f"znconreylog(znstar({N},1),{a})")))
    v = pari(f"chareval(znstar({N},1), znconreylog(znstar({N},1),{a}), {n}, [exp(2*Pi*I/{o}), {o}])")
    return complex(v)


def eigen_chi(N, chi, ell=2):
    """Eigenvalues of T_l on S_2(N, chi) as complex numbers; chi is a PARI character."""
    mf = pari.mfinit([N, 2, chi], 1)
    d = int(pari.mfdim(mf))
    if d == 0:
        return np.zeros(0, dtype=complex)
    E = pari.mfembed(mf, pari.mfheckemat(mf, ell))
    A = np.array([[complex(E[j][i]) for j in range(d)] for i in range(d)])
    return np.linalg.eigvals(A)


def gamma1_spectrum(N, ell=2):
    """[(a, chi(l), eigenvalues)] over the even characters chi_a mod N."""
    out = []
    for a in even_labels(N):
        lam = eigen_chi(N, pari(f"Mod({a},{N})"), ell)
        c = chi_value(N, a, ell)
        if len(lam):
            check(np.max(np.abs(lam)) <= 2 * ell ** 0.5 + 1e-6, f"|lambda| > 2 sqrt l at N={N}, a={a}")
            check(np.max(np.abs(lam - c * np.conj(lam))) < 1e-6, f"lambda != chi(l) conj(lambda) at N={N}, a={a}")
        out.append((a, c, lam))
    total = sum(len(l) for _, _, l in out)
    check(total == (N - 5) * (N - 7) // 24, f"dim S_2(Gamma_1({N})) = {total}")
    # conj(chi) has Conrey label a^-1 mod N; its eigenvalues are the conjugates.
    by = {a: lam for a, _, lam in out}
    for a, _, lam in out:
        b = pow(a, -1, N)
        check(same_multiset(np.conj(lam), by[b]), f"eigenvalues for conj(chi) at N={N}, a={a}")
    return out


def same_multiset(u, v, tol=1e-6):
    """Whether two lists of complex numbers agree up to order and rounding."""
    if len(u) != len(v):
        return False
    left = list(v)
    for z in u:
        k = min(range(len(left)), key=lambda i: abs(left[i] - z))
        if abs(left[k] - z) > tol:
            return False
        left.pop(k)
    return True


def make_gamma1():
    N = GAMMA1_N
    rows = []
    for a, c, lam in gamma1_spectrum(N):
        rows.append({"a": a, "chi": [round(c.real, 6), round(c.imag, 6)],
                     "eig": [[round(z.real, 4), round(z.imag, 4)] for z in lam]})
    order = int(pari(f"znorder(Mod(2,{N}))"))
    lines = order if order % 2 else order // 2
    return {"N": N, "ell": 2, "lines": lines, "chars": rows}


def quadratic_row(N, y):
    """The row of quadratic.json for N from y = lambda/(i sqrt 2) over the eigenvalues
    lambda of T_2 on S_2(N, (./N)). Also used by magma_import.py."""
    y = np.asarray(y, dtype=float)
    check(np.allclose(np.sort(y), np.sort(-y), atol=1e-6), f"spectrum not symmetric at N = {N}")
    return {"N": N, "dim": len(y), "hist": hist(y), "ks": round(ks(y, 2), 4)}


def make_quadratic():
    rows = []
    for N in QUADRATIC_NS:
        check(N % 8 == 5 and int(pari(f"isprime({N})")), f"{N} is not a prime = 5 mod 8")
        chi = pari(f"[znstar({N},1), [({N}-1)/2]]")
        two = complex(pari(f"chareval(znstar({N},1), [({N}-1)/2], 2, [-1, 2])"))
        check(abs(two + 1) < 1e-12, f"chi(2) != -1 for N = {N}")
        lam = eigen_chi(N, chi)
        check(np.max(np.abs(lam.real)) < 1e-6, f"a real part at N = {N}")
        check(np.max(np.abs(lam)) <= 2 * 2 ** 0.5 + 1e-6, f"|lambda| > 2 sqrt 2 at N = {N}")
        rows.append(quadratic_row(N, lam.imag / 2 ** 0.5))
        print(f"  quadratic N={N}: dim {rows[-1]['dim']}, KS {rows[-1]['ks']}", flush=True)
    return {"ell": 2, "bins": BINS, "rows": rows}


def gamma1_by_orbits(N, ell=2):
    """The eigenvalues of T_l on S_2(Gamma_1(N)), N prime, one Galois orbit of characters
    at a time: [(o, c, eigenvalues)] over the even characters chi, o the order of chi and
    c = chi(l).

    The even characters of order o (o | (N-1)/2) form one Galois orbit. PARI computes T_l
    on S_2(N, chi) once, for one chi of order o, as a matrix over Q(chi) = Q[t]/Phi(t);
    substituting each complex root of Phi for t gives T_l on S_2(N, chi') for every chi'
    in the orbit. This is much faster than one mfinit per character (N = 251: seconds
    instead of minutes). chi(l) is read off the spectrum: lambda = chi(l) conj(lambda)
    for every eigenvalue lambda, and it must be an o-th root of unity.
    """
    h = (N - 1) // 2
    out = []
    for o in [d for d in range(1, h + 1) if h % d == 0]:
        # znstar(N, 1) is cyclic of order N - 1; exponent (N - 1)/o gives order o, and the
        # character is even because o divides h.
        mf = pari.mfinit([N, 2, pari(f"[znstar({N},1), [{(N - 1) // o}]]")], 1)
        d = int(pari.mfdim(mf))
        P = pari.mfparams(mf)[4]
        deg = int(pari.poldegree(P))
        roots = [complex(r) for r in pari.polroots(P)] if deg > 1 else [1.0]
        check(len(roots) == int(pari.eulerphi(o)), f"orbit of order {o} at N = {N}")
        if d == 0:
            continue
        M = pari.liftall(pari.mfheckemat(mf, ell))
        C = np.zeros((max(deg, 1), d, d))
        for i in range(d):
            for j in range(d):
                e = M[j][i]                      # cypari2: M[j] is column j
                if deg <= 1:
                    C[0, i, j] = float(e)
                else:
                    for k, c in enumerate(pari.Vecrev(e)):
                        C[k, i, j] = float(c)
        for z in roots:
            lam = np.linalg.eigvals(np.tensordot(z ** np.arange(C.shape[0]), C, axes=1))
            check(np.max(np.abs(lam)) <= 2 * ell ** 0.5 + 1e-6, f"|lambda| > 2 sqrt l at N={N}, o={o}")
            # The matrices have large rational entries, so the eigenvalues are only good to
            # about 1e-6 here (they are stored to 1e-3); hence the tolerance 1e-4.
            big = lam[np.argmax(np.abs(lam))]
            c = big / np.conj(big)
            check(np.max(np.abs(lam - c * np.conj(lam))) < 1e-4, f"lambda != c conj(lambda) at N={N}, o={o}")
            check(abs(c ** o - 1) < 1e-4, f"chi(l) is not an o-th root of unity at N={N}, o={o}")
            out.append((o, c, lam))
    total = sum(len(l) for _, _, l in out)
    check(total == (N - 5) * (N - 7) // 24, f"dim S_2(Gamma_1({N})) = {total}")
    return out


PLANAR_SAMPLE = 3000     # points drawn per N; the statistics use all of them


def make_planar():
    # The orbit method agrees with the one-character-at-a-time method of gamma1_spectrum.
    a = np.concatenate([lam for _, _, lam in gamma1_spectrum(GAMMA1_N)])
    b = np.concatenate([lam for _, _, lam in gamma1_by_orbits(GAMMA1_N)])
    check(same_multiset(a, b, 1e-8), f"the two methods disagree at N = {GAMMA1_N}")
    rows = []
    for N in PLANAR_NS:
        t = time.time()
        z = np.concatenate([lam for _, _, lam in gamma1_by_orbits(N)]) / 2 ** 0.5
        rows.append(planar_row(N, z))
        print(f"  planar N={N}: {rows[-1]['dim']} points, radial KS {rows[-1]['ksr']:.4f}, "
              f"|m20| {rows[-1]['mom']['20']}, m11 {rows[-1]['mom']['11']}, {time.time() - t:.0f}s", flush=True)
    return {"ell": 2, "sample": PLANAR_SAMPLE, "rows": rows}


def planar_row(N, z):
    """The row of planar.json for N from z = lambda/sqrt 2 over all the eigenvalues lambda
    of T_2 on S_2(Gamma_1(N)). Also used by magma_import.py."""
    z = np.asarray(z, dtype=complex)
    r = np.abs(z)
    # the limit: |z| distributed as |X|, X ~ mu_2; the CDF of |X| at r is 2 F(r) - 1
    rs = np.sort(r)
    m = len(rs)
    F = 2 * np.interp(rs, _GRID, cdf(2)) - 1
    i = np.arange(m)
    ks_r = float(max(np.max(np.abs((i + 1) / m - F)), np.max(np.abs(i / m - F))))
    mom = {f"{p}{q}": round(abs(complex(np.mean(z ** p * np.conj(z) ** q))), 4)
           for p, q in ((2, 0), (4, 0), (3, 1), (6, 0), (1, 1), (2, 2))}
    rhist, _ = np.histogram(np.minimum(r, 2), bins=20, range=(0, 2))
    order = next(k for k in range(1, N) if pow(2, k, N) == 1)
    shown = z
    if m > PLANAR_SAMPLE:   # a fixed random sample, the same on every run
        shown = z[np.sort(np.random.default_rng(N).choice(m, PLANAR_SAMPLE, replace=False))]
    return {"N": N, "dim": m, "lines": order if order % 2 else order // 2,
            "ksr": round(ks_r, 4), "mom": mom, "rhist": [int(c) for c in rhist],
            "pts": [v for w in shown for v in (round(w.real, 3), round(w.imag, 3))]}


# ---- Brandt matrices ---------------------------------------------------------------

BRANDT_GP = r"""
ssj(p) = {
  my(g = ffgen(ffinit(p,2),'a), D = [-3,-4,-7,-8,-11,-19,-43,-67,-163],
     J = [0,1728,-3375,8000,-32768,-884736,-884736000,-147197952000,-262537412640768000],
     j0, P2, seen, todo, out);
  for(i=1,#D, if(kronecker(D[i],p) != 1, j0 = J[i]*g^0; break));
  P2 = polmodular(2)*Mod(1,p);
  out = List([j0]); seen = Map(); mapput(seen, j0, 1); todo = List([j0]);
  while(#todo, my(j = todo[#todo]); listpop(todo);
    my(F = factor(subst(P2,'x,j)));
    for(k=1,#F~, my(r = -polcoef(F[k,1],0)/polcoef(F[k,1],1));
      if(!mapisdefined(seen,r), mapput(seen,r,#out+1); listput(out,r); listput(todo,r))));
  [Vec(out), seen];
}
brandt(S, l, p) = {
  my(js = S[1], seen = S[2], n = #js, Pl = polmodular(l)*Mod(1,p), B = matrix(n,n));
  for(i=1,n, my(F = factor(subst(Pl,'x,js[i])));
    for(k=1,#F~, if(poldegree(F[k,1]) != 1, error("Phi_l(j, Y) does not split"));
      my(r = -polcoef(F[k,1],0)/polcoef(F[k,1],1), c);
      if(!mapisdefined(seen,r,&c), error("not supersingular"));
      B[i,c] += F[k,2]));
  B;
}
autw(S) = apply(j -> if(j == 0, 3, if(j == 1728, 2, 1)), S[1]);
"""


def n_supersingular(p):
    return p // 12 + {1: 0, 5: 1, 7: 1, 11: 2}[p % 12]


def brandt_eigenvalues(p, ell):
    """Eigenvalues of T_l on S_2(Gamma_0(p)): those of B(l) without l + 1."""
    pari(f"S = ssj({p})")
    n = int(pari("#S[1]"))
    check(n == n_supersingular(p), f"number of supersingular j at p = {p}")
    w = np.array([int(x) for x in pari("autw(S)")], dtype=float)
    mats = {}
    for l in sorted({2, ell}):
        B = pari(f"brandt(S, {l}, {p})")
        A = np.array([[int(B[j][i]) for j in range(n)] for i in range(n)], dtype=float)
        check(np.all(A.sum(axis=1) == l + 1), f"row sums of B({l}) at p = {p}")
        check(np.array_equal(A * w[None, :], (A * w[None, :]).T), f"B({l}) D not symmetric at p = {p}")
        mats[l] = A
    return brandt_spectrum(p, ell, mats, w)


def brandt_spectrum(p, ell, mats, w):
    """The eigenvalues of T_l on S_2(Gamma_0(p)) from the Brandt matrices mats = {2: B(2),
    l: B(l)} (rows sum to l + 1) and the weights w_i = |Aut(E_i)|/2, after the checks.
    Also used by magma_import.py."""
    n = len(w)
    check(n == n_supersingular(p), f"number of supersingular j at p = {p}")
    for l, A in mats.items():
        check(np.all(A.sum(axis=1) == l + 1), f"row sums of B({l}) at p = {p}")
        check(np.array_equal(A * w[None, :], (A * w[None, :]).T), f"B({l}) D not symmetric at p = {p}")
    check(np.array_equal(mats[2] @ mats[ell], mats[ell] @ mats[2]), f"B(2), B({ell}) do not commute at p = {p}")
    s = np.sqrt(w)
    S = mats[ell] * w[None, :] / s[:, None] / s[None, :]
    lam = np.sort(np.linalg.eigvalsh(S))
    check(abs(lam[-1] - (ell + 1)) < 1e-7, f"top eigenvalue at p = {p}")
    lam = lam[:-1]
    check(np.max(np.abs(lam)) <= 2 * ell ** 0.5 + 1e-7, f"Ramanujan bound at p = {p}, l = {ell}")
    return lam


def diagonal_primes(count=28, lo=101, hi=20011):
    out = []
    for t in np.linspace(math.log(lo), math.log(hi), count):
        q = int(pari(f"nextprime({int(round(math.exp(t)))})"))
        if not out or q > out[-1]:
            out.append(q)
    return out


def ell_for(p):
    """The largest prime <= (log p)^2."""
    n = int(math.log(p) ** 2)
    while n > 2 and not all(n % d for d in range(2, int(n ** 0.5) + 1)):
        n -= 1
    return n


def diagonal_row(p, ell, lam):
    """The row of diagonal.json from the eigenvalues lam of T_l on S_2(Gamma_0(p)).
    Also used by magma_import.py."""
    check(ell == ell_for(p), f"l = {ell} is not the largest prime <= (log {p})^2")
    x = np.asarray(lam, dtype=float) / ell ** 0.5
    return {"p": p, "ell": ell, "dim": len(x), "hist": hist(x),
            "ks": round(ks(x, math.inf), 4), "ksl": round(ks(x, ell), 4)}


def make_diagonal():
    rows = []
    for p in DIAGONAL_PS or diagonal_primes():
        t = time.time()
        ell = ell_for(p)
        rows.append(diagonal_row(p, ell, brandt_eigenvalues(p, ell)))
        print(f"  diagonal p={p} l={ell}: dim {rows[-1]['dim']}, KS to ST {rows[-1]['ks']}, "
              f"to mu_l {rows[-1]['ksl']}, {time.time() - t:.0f}s", flush=True)
    return {"bins": BINS, "rows": rows}


# ---- level one ----------------------------------------------------------------------

LEVEL1_GP = r"""
level1dim(k) = if(k%2 || k < 12, 0, k\12 - (k%12 == 2));
level1hecke(k, l) = {
  my(d = level1dim(k), n = l*(d+1) + 2, q = 'q,
     E4 = 1 + 240*sum(m=1,n,sigma(m,3)*q^m) + O(q^(n+1)),
     E6 = 1 - 504*sum(m=1,n,sigma(m,5)*q^m) + O(q^(n+1)),
     D = q*eta(q + O(q^(n+1)))^24, B = vector(d), M = matrix(d,d), C = matrix(d,d));
  for(i=1,d, my(r = k - 12*i, a = -1, b);
    forstep(bb = r\6, 0, -1, if((r - 6*bb) % 4 == 0, a = (r - 6*bb)/4; b = bb; break));
    if(a < 0, error("no E_4^a E_6^b of weight ", r));
    B[i] = D^i * E4^a * E6^b);
  for(i=1,d, for(m=1,d,
    M[m,i] = polcoef(B[i], l*m) + if(m % l == 0, l^(k-1)*polcoef(B[i], m/l), 0);
    C[m,i] = polcoef(B[i], m)));
  matsolve(C, M);
}
level1eig(k, l) = {
  my(P = charpoly(level1hecke(k, l)), R);
  localprec(max(100, ceil(k*log(l)/log(10)) + 60));
  R = polroots(P);
  if(vecmax(apply(z -> abs(imag(z)), R)) > 1e-20 * l^((k-1)/2), error("non-real root"));
  vector(#R, i, real(R[i]) / l^((k-1)/2.));
}
"""


def level1_dim(k):
    return 0 if k % 2 or k < 12 else k // 12 - (k % 12 == 2)


def weight_row(k, l, x):
    """The row of weight.json from x = lambda / l^((k-1)/2) over the eigenvalues lambda of
    T_l on S_k(SL_2(Z)). Also used by magma_import.py."""
    x = np.asarray(x, dtype=float)
    check(len(x) == level1_dim(k), f"dim S_{k}(SL_2(Z)) = {len(x)}")
    check(np.max(np.abs(x)) <= 2 + 1e-9, f"Deligne bound at k = {k}, l = {l}")
    return {"k": k, "dim": len(x), "hist": hist(x), "ks": round(ks(x, l), 4)}


def make_weight():
    check(int(pari("level1hecke(12,2)[1,1]")) == -24, "T_2 on S_12 is not tau(2) = -24")
    for k in range(12, 121, 2):
        if int(pari(f"level1dim({k})")):
            for l in WEIGHT_ELLS:
                check(bool(pari(f"charpoly(level1hecke({k},{l})) == charpoly(mfheckemat(mfinit([1,{k}],1),{l}))")),
                      f"level one T_{l} at k = {k}")
    out = {"bins": BINS, "ells": {}}
    for l in WEIGHT_ELLS:
        rows = []
        for k in WEIGHT_KS:
            rows.append(weight_row(k, l, [float(v) for v in pari(f"level1eig({k},{l})")]))
        print(f"  weight l={l}: k up to {WEIGHT_KS[-1]}, dim {rows[-1]['dim']}, KS {rows[-1]['ks']}", flush=True)
        out["ells"][str(l)] = rows
    return out


MAKERS = {"gamma1": make_gamma1, "quadratic": make_quadratic, "planar": make_planar,
          "diagonal": make_diagonal, "weight": make_weight}


def main(argv):
    start_pari()
    names = argv or list(MAKERS)
    OUT.mkdir(parents=True, exist_ok=True)
    for name in names:
        t = time.time()
        print(f"{name} ...", flush=True)
        data = MAKERS[name]()
        (OUT / f"{name}.json").write_text(json.dumps(data, separators=(",", ":")) + "\n")
        print(f"{name}: written in {time.time() - t:.0f}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
