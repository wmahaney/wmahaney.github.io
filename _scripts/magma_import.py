#!/usr/bin/env python3
"""Check the output of _scripts/magma/harvest.m and merge it into the site's data stores.

    python3 _scripts/magma_import.py harvest.txt [more.txt ...]      # check, then write
    python3 _scripts/magma_import.py --dry-run harvest.txt           # check only
    python3 _scripts/magma_import.py --selftest                      # see below

Needs Python 3 and numpy (and cypari2 for --selftest only). The rows are built by the
same functions that build them from PARI (quadratic_row, planar_row, diagonal_row,
weight_row and brandt_spectrum in limit_laws.py; spectra_row in isogeny_graphs.py),
so the histograms, distances and checks are identical whichever program computed the
eigenvalues.

Stores and records (the record format is described at the top of harvest.m):

  limits/planar.json      planar records, grouped by N
  limits/quadratic.json   quadratic records
  limits/weight.json      weight records
  limits/diagonal.json    brandt records for l = 2 and l = the largest prime <= (log p)^2
  isogeny/spectra.json    brandt records for l = 2, 3 and 5

Checks, on top of those in the row functions (the import stops at the first failure,
before anything is written):
  * planar: |lambda| <= 2 sqrt 2; lambda = c conj(lambda) with c = chi(2), |c| = 1 and
    c^o = 1; each orbit (o) appears with phi(o) distinct conjugates k, all of the same
    dimension; the dimensions add up to the genus (N - 5)(N - 7)/24 of X_1(N);
  * brandt: all n rows are present once; the matrix, or its transpose, has rows summing
    to l + 1; the weights w (B(l)_ij w_j = B(l)_ji w_i) are found from B(2), lie in
    {1, 2, 3}, and make every B(l) D symmetric; the B(l) commute;
  * a record for a key already in a store (same N, p or k) must agree with it: the same
    dimension, Kolmogorov distances within 2e-4, and histograms differing in at most 2
    counts per bin (eigenvalues on a bin edge). The new row then replaces the old one.
    This compares Magma's numbers with PARI's wherever the two overlap.

--selftest writes harvest.m-style records for small cases from PARI (via limit_laws.py),
with one Brandt matrix transposed, imports them into a copy of the stores in a temporary
directory, and checks that every row matches the stored one. It tests this importer, not
harvest.m.
"""
from __future__ import annotations

import json
import math
import shutil
import sys
import tempfile
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import limit_laws as LL          # noqa: E402
import isogeny_graphs as IG      # noqa: E402

DATA = HERE.parent / "assets" / "data"


def fail(msg):
    raise SystemExit("IMPORT FAILED: " + msg)


def phi(n):
    out, m, q = n, n, 2
    while q * q <= m:
        if m % q == 0:
            while m % q == 0:
                m //= q
            out -= out // q
        q += 1
    return out - out // m if m > 1 else out


# ---- parsing ---------------------------------------------------------------------------

def parse(paths):
    planar = defaultdict(list)      # N -> [(o, k, d, c, lam)]
    quadratic, weight = {}, {}      # N -> y; (l, k) -> x
    brandt = defaultdict(dict)      # p -> {l: (n, {i: [(j, v)]})}
    for path in paths:
        for lineno, line in enumerate(Path(path).read_text().splitlines(), 1):
            f = line.split()
            if not f or f[0].startswith("#"):
                continue
            where = f"{path}:{lineno}"
            try:
                if f[0] == "planar":
                    N, o, k, d = map(int, f[1:5])
                    c = complex(float(f[5]), float(f[6]))
                    v = np.array(f[7:], dtype=float)
                    if len(v) != 2 * d:
                        fail(f"{where}: {len(v) // 2} eigenvalues, expected {d}")
                    planar[N].append((o, k, d, c, v[0::2] + 1j * v[1::2]))
                elif f[0] == "quadratic":
                    N, d = int(f[1]), int(f[2])
                    y = np.array(f[3:], dtype=float)
                    if len(y) != d:
                        fail(f"{where}: {len(y)} values, expected {d}")
                    quadratic[N] = y
                elif f[0] == "weight":
                    l, k, d = map(int, f[1:4])
                    x = np.array(f[4:], dtype=float)
                    if len(x) != d:
                        fail(f"{where}: {len(x)} values, expected {d}")
                    weight[(l, k)] = x
                elif f[0] == "brandt":
                    p, l, n, i = map(int, f[1:5])
                    rows = brandt[p].setdefault(l, (n, {}))[1]
                    if brandt[p][l][0] != n or not 1 <= i <= n or i in rows:
                        fail(f"{where}: row {i} of B({l}) at p = {p} is repeated or out of range")
                    rows[i] = [tuple(map(int, e.split(","))) for e in f[5:]]
                else:
                    fail(f"{where}: unknown record {f[0]!r}")
            except ValueError as e:
                fail(f"{where}: {e}")
    return planar, quadratic, weight, brandt


# ---- checks that need the raw data --------------------------------------------------------

def planar_points(N, recs):
    h = (N - 1) // 2
    seen = set()
    dims = defaultdict(set)
    total = 0
    for o, k, d, c, lam in recs:
        if h % o:
            fail(f"planar N={N}: order {o} does not divide (N-1)/2")
        if math.gcd(k, o) != 1 or (o, k % o) in seen:
            fail(f"planar N={N}: conjugate k={k} of order {o} is repeated or not prime to {o}")
        seen.add((o, k % o))
        dims[o].add(d)
        if abs(abs(c) - 1) > 1e-8 or abs(c ** o - 1) > 1e-7:
            fail(f"planar N={N}, o={o}, k={k}: chi(2) = {c} is not an o-th root of unity")
        if np.max(np.abs(lam)) > 2 * 2 ** 0.5 + 1e-8:
            fail(f"planar N={N}, o={o}, k={k}: |lambda| > 2 sqrt 2")
        if np.max(np.abs(lam - c * np.conj(lam))) > 1e-7:
            fail(f"planar N={N}, o={o}, k={k}: lambda != chi(2) conj(lambda)")
        total += d
    for o, ds in dims.items():
        if len(ds) != 1:
            fail(f"planar N={N}: conjugates of order {o} have different dimensions {ds}")
        if sum(1 for (oo, _) in seen if oo == o) != phi(o):
            fail(f"planar N={N}: order {o} has not all of its {phi(o)} conjugates")
    if total != (N - 5) * (N - 7) // 24:
        fail(f"planar N={N}: {total} eigenvalues, but the genus of X_1(N) is {(N - 5) * (N - 7) // 24}")
    return np.concatenate([lam for *_, lam in recs]) / 2 ** 0.5


def brandt_matrices(p, recs):
    """{l: B(l)} with rows summing to l + 1, and the weights w, from the brandt records."""
    if 2 not in recs:
        fail(f"brandt p={p}: B(2) is missing (it fixes the weights)")
    ns = {n for n, _ in recs.values()}
    if len(ns) != 1:
        fail(f"brandt p={p}: different dimensions {ns}")
    n = ns.pop()
    mats = {}
    for l, (_, rows) in recs.items():
        if len(rows) != n:
            fail(f"brandt p={p}: B({l}) has {len(rows)} of its {n} rows")
        A = np.zeros((n, n), dtype=np.int64)
        for i, entries in rows.items():
            for j, v in entries:
                A[i - 1, j - 1] += v
        if not np.all(A.sum(axis=1) == l + 1):
            A = A.T.copy()
        if not np.all(A.sum(axis=1) == l + 1):
            fail(f"brandt p={p}: neither the rows nor the columns of B({l}) sum to {l + 1}")
        mats[l] = A
    # weights: B_ij w_j = B_ji w_i along the edges of the (connected) 2-isogeny graph
    A = mats[2]
    w = [None] * n
    w[0] = Fraction(1)
    todo = [0]
    while todo:
        i = todo.pop()
        for j in np.nonzero(A[i])[0]:
            if w[j] is None:
                if A[j, i] == 0:
                    fail(f"brandt p={p}: B(2) has an edge {i}->{j} without a reverse edge")
                w[j] = w[i] * int(A[j, i]) / int(A[i, j])
                todo.append(int(j))
    if any(x is None for x in w):
        fail(f"brandt p={p}: the 2-isogeny graph is not connected")
    base = min(w)
    w = [x / base for x in w]
    if any(x not in (1, 2, 3) for x in w):
        fail(f"brandt p={p}: weights {sorted(set(w))} are not in {{1, 2, 3}}")
    return {l: M.astype(float) for l, M in mats.items()}, np.array([float(x) for x in w])


# ---- merging ---------------------------------------------------------------------------------

def agree(old, new, keys, what):
    """Fail unless the stored row old and the new row agree (see the docstring)."""
    if old["dim"] != new["dim"]:
        fail(f"{what}: dimension {new['dim']} here, {old['dim']} in the store")
    for k in keys:
        a, b = old[k], new[k]
        pairs = list(a.items()) if isinstance(a, dict) else [(None, a)]
        for key, va in pairs:
            vb = b[key] if key is not None else b
            if isinstance(va, list):
                if len(va) != len(vb) or max(abs(x - y) for x, y in zip(va, vb)) > 2:
                    fail(f"{what}: histogram {key or k} differs from the store")
            elif abs(va - vb) > 2e-4:
                fail(f"{what}: {key or k} = {vb} here, {va} in the store")


def merge(rows, new, key, keys, what):
    out = {r[key]: r for r in rows}
    for r in new:
        if r[key] in out:
            agree(out[r[key]], r, keys, f"{what} {key}={r[key]}")
        out[r[key]] = r
    return [out[k] for k in sorted(out)]


def run(paths, data=DATA, dry=False):
    planar, quadratic, weight, brandt = parse(paths)
    stores = {}
    report = []

    if planar:
        path = data / "limits" / "planar.json"
        store = json.loads(path.read_text())
        new = [LL.planar_row(N, planar_points(N, recs)) for N, recs in sorted(planar.items())]
        # pts are a random sample: compare the statistics, not the points
        store["rows"] = merge(store["rows"], new, "N", ["ksr", "rhist"], "planar")
        stores[path] = store
        report.append(f"planar: N = {', '.join(str(r['N']) for r in new)}")

    if quadratic:
        path = data / "limits" / "quadratic.json"
        store = json.loads(path.read_text())
        new = []
        for N, y in sorted(quadratic.items()):
            if N % 8 != 5:
                fail(f"quadratic N={N}: N is not 5 mod 8")
            new.append(LL.quadratic_row(N, y))
        store["rows"] = merge(store["rows"], new, "N", ["ks", "hist"], "quadratic")
        stores[path] = store
        report.append(f"quadratic: N = {', '.join(str(r['N']) for r in new)}")

    if weight:
        path = data / "limits" / "weight.json"
        store = json.loads(path.read_text())
        by_l = defaultdict(list)
        for (l, k), x in sorted(weight.items()):
            by_l[l].append(LL.weight_row(k, l, x))
        for l, new in by_l.items():
            store["ells"][str(l)] = merge(store["ells"].get(str(l), []), new, "k", ["ks", "hist"], f"weight l={l}")
            report.append(f"weight l={l}: k = {new[0]['k']}..{new[-1]['k']} ({len(new)})")
        stores[path] = store

    diag, spec = [], []
    for p, recs in sorted(brandt.items()):
        mats, w = brandt_matrices(p, recs)
        ell = LL.ell_for(p)
        if ell in mats:
            diag.append(LL.diagonal_row(p, ell, LL.brandt_spectrum(p, ell, {2: mats[2], ell: mats[ell]}, w)))
        if {2, 3, 5} <= set(mats):
            for l in (3, 5):
                if not np.array_equal(mats[2] @ mats[l], mats[l] @ mats[2]):
                    fail(f"brandt p={p}: B(2) and B({l}) do not commute")
            spec.append(IG.spectra_row(p, list(2 * w), {l: mats[l] for l in (2, 3, 5)}))
        if ell not in mats and not {2, 3, 5} <= set(mats):
            fail(f"brandt p={p}: neither l = {ell} (diagonal) nor l = 2, 3, 5 (spectra) is complete")
    if diag:
        path = data / "limits" / "diagonal.json"
        store = json.loads(path.read_text())
        store["rows"] = merge(store["rows"], diag, "p", ["ks", "ksl", "hist"], "diagonal")
        stores[path] = store
        report.append(f"diagonal: p = {', '.join(str(r['p']) for r in diag)}")
    if spec:
        path = data / "isogeny" / "spectra.json"
        store = json.loads(path.read_text())
        store["primes"] = merge(store["primes"], spec, "p", ["ks", "hist"], "spectra")
        stores[path] = store
        report.append(f"spectra: p = {', '.join(str(r['p']) for r in spec)}")

    for line in report:
        print(line)
    if dry:
        print("dry run: nothing written")
        return
    for path, store in stores.items():
        path.write_text(json.dumps(store, separators=(",", ":")) + "\n")
        print(f"wrote {path.relative_to(data.parent.parent)}")


# ---- selftest ------------------------------------------------------------------------------

def num(x):
    return f"{x:.12g}"


def selftest():
    LL.start_pari()
    pari = LL.pari
    lines = ["# selftest records, written from PARI"]
    for N in (29, 101):                                    # planar
        ks = defaultdict(int)
        for o, c, lam in LL.gamma1_by_orbits(N):
            ks[o] += 1
            k = [k for k in range(1, max(o, 1) + 1) if math.gcd(k, o) == 1][ks[o] - 1]
            vals = [num(t) for z in lam for t in (z.real, z.imag)]
            lines.append(" ".join(["planar", str(N), str(o), str(k), str(len(lam)),
                                   num(c.real), num(c.imag)] + vals))
    for N in (101, 149):                                   # quadratic
        lam = LL.eigen_chi(N, pari(f"[znstar({N},1), [({N}-1)/2]]"))
        y = lam.imag / 2 ** 0.5
        lines.append(" ".join(["quadratic", str(N), str(len(y))] + [num(t) for t in y]))
    for k in (24, 48, 120):                                # weight
        x = [float(v) for v in pari(f"level1eig({k},2)")]
        lines.append(" ".join(["weight", "2", str(k), str(len(x))] + [num(t) for t in x]))
    for p, ells, transpose in ((101, (2, 3, 5, 19), False), (1061, (2, 47), True)):   # brandt
        pari(f"S = ssj({p})")
        n = int(pari("#S[1]"))
        for l in ells:
            B = pari(f"brandt(S, {l}, {p})")
            A = np.array([[int(B[j][i]) for j in range(n)] for i in range(n)])
            if transpose:          # as if the matrix acted on the other side
                A = A.T
            for i in range(n):
                ent = [f"{j + 1},{A[i, j]}" for j in range(n) if A[i, j]]
                lines.append(" ".join(["brandt", str(p), str(l), str(n), str(i + 1)] + ent))
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        rec = tmp / "records.txt"
        rec.write_text("\n".join(lines) + "\n")
        copy = tmp / "data"
        shutil.copytree(DATA, copy)
        run([rec], data=copy)                             # fails if a row disagrees with the store
        # and the merged stores still hold exactly the keys they had
        for rel, key, field in (("limits/planar.json", "N", "rows"), ("limits/quadratic.json", "N", "rows"),
                                ("limits/diagonal.json", "p", "rows"), ("isogeny/spectra.json", "p", "primes")):
            a = [r[key] for r in json.loads((DATA / rel).read_text())[field]]
            b = [r[key] for r in json.loads((copy / rel).read_text())[field]]
            assert a == b, f"{rel}: keys changed"
    print("selftest passed: PARI's records import to the rows already stored")


def main(argv):
    if "--selftest" in argv:
        selftest()
        return
    paths = [a for a in argv if not a.startswith("--")]
    if not paths:
        raise SystemExit(__doc__)
    run(paths, dry="--dry-run" in argv)


if __name__ == "__main__":
    main(sys.argv[1:])
