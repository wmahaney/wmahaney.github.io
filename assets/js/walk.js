// Random walks on the isogeny graphs: the transition matrix, built from the adjacency
// matrix. Pure functions: nothing here draws or touches the page, so they are tested on
// their own (see _tests/walk.test.mjs; run "node --test _tests/*.test.mjs").
//
// Notation. A is the adjacency matrix of a directed multigraph: A[v][w] is the number of
// edges from v to w. For the supersingular l-isogeny graph, A[v][w] counts the subgroups
// C of E_v of order l with E_v / C isomorphic to E_w, so every row sums to l + 1.
// Matrices are stored by rows, each row a Map from a column index to a nonzero entry.

// ---- adjacency -----------------------------------------------------------------

// The two counts of a link in the data files (assets/data/isogeny/). A link is
//   [s, t, ell, m]        m edges from s to t and, unless the graph is directed, m back;
//   [s, t, ell, m, back]  m edges from s to t and "back" edges from t to s (they differ
//                         only next to a vertex with automorphisms other than +-1).
// A loop [s, s, ell, m] is m edges from s to itself.
export function linkCounts(link, directed = false) {
  const [s, t, ell, m, back] = link;
  if (s === t) return { s, t, ell, out: m, back: 0 };
  return { s, t, ell, out: m, back: back !== undefined ? back : directed ? 0 : m };
}

// The adjacency matrix of n vertices from a list of links.
export function adjacency(n, links, directed = false) {
  const A = Array.from({ length: n }, () => new Map());
  const add = (v, w, k) => { if (k) A[v].set(w, (A[v].get(w) || 0) + k); };
  for (const link of links) {
    const { s, t, out, back } = linkCounts(link, directed);
    add(s, t, out);
    add(t, s, back);
  }
  return A;
}

export const rowSum = (row) => { let x = 0; for (const v of row.values()) x += v; return x; };

// ---- walks ---------------------------------------------------------------------

// The simple random walk P0 = D^{-1} A, D the diagonal matrix of row sums: from v, each
// outgoing edge is equally likely. For the l-isogeny graph this is A / (l + 1). A vertex
// with no outgoing edge gets an empty row.
export function simpleWalk(A) {
  return A.map((row) => {
    const d = rowSum(row), out = new Map();
    for (const [w, a] of row) out.set(w, a / d);
    return out;
  });
}

// The closed 1-balls around the vertices in "centres": each centre, and every vertex
// joined to a centre by an edge in either direction.
export function closedBalls(A, centres) {
  const ball = new Set(centres), isCentre = new Set(centres);
  A.forEach((row, v) => {
    for (const w of row.keys()) {
      if (isCentre.has(v)) ball.add(w);
      if (isCentre.has(w)) ball.add(v);
    }
  });
  return ball;
}

// THE PETERSSON WALK. On the supersingular curves E_1, ..., E_n in characteristic p, the
// pairing <E_i, E_j> = w_i if i = j and 0 otherwise, with w_i = |Aut(E_i)| / 2, makes the
// Hecke operator T (T E_i = the sum of the targets of the isogenies out of E_i, one per
// kernel) self-adjoint: <T E_i, E_j> = A[i][j] w_j = A[j][i] w_i (Kohel's thesis,
// chapter 7, after Mestre and Oesterle). The Petersson walk takes its
// transition probabilities from this pairing:
//
//     P[i][j] = <T E_i, E_j> / sum over k of <T E_i, E_k> = A[i][j] w_j / s_i,
//     s_i = sum over k of A[i][k] w_k.
//
// A[i][j] w_j is the number of isogenies from E_i to E_j, counted up to sign instead of
// up to the automorphisms of E_j, so the walk picks an isogeny out of E_i uniformly
// where the simple random walk picks a kernel. The two agree at every vertex with no
// neighbour of extra automorphisms, that is outside the closed 1-balls of j = 0 and
// j = 1728. The matrix A[i][j] w_j is symmetric, so the walk is reversible and its
// stationary vector is proportional to s_i (see isogenyCounts).
//
//   weight(v)  w_v, or any positive multiple of it (only ratios matter).
export function peterssonWalk(A, weight) {
  return A.map((row) => {
    let s = 0;
    for (const [w, a] of row) s += a * weight(w);
    const out = new Map();
    for (const [w, a] of row) out.set(w, (a * weight(w)) / s);
    return out;
  });
}

// s_i = sum over k of A[i][k] w_k: the isogenies out of E_i, counted up to sign. The
// stationary vector of the Petersson walk is proportional to it.
export function isogenyCounts(A, weight) {
  return A.map((row) => { let s = 0; for (const [w, a] of row) s += a * weight(w); return s; });
}

// w_v from the "aut" field of a graph in the data files: {vertex: |Aut|} for the
// vertices where |Aut| is not 2.
export const weightFrom = (aut) => (v) => ((aut && aut[v]) || 2) / 2;

// The walk the animation runs: the Petersson walk where the automorphism groups are
// known (aut is given: the full graphs without level structure), and the simple random
// walk on every other graph shown.
export function walkFor(A, aut = null) {
  return aut ? peterssonWalk(A, weightFrom(aut)) : simpleWalk(A);
}
export const walkName = (aut = null) => (aut ? "petersson" : "simple");

// ---- stationary vector ---------------------------------------------------------

// A solution of pi P = pi with entries summing to 1, by iterating pi -> pi (P + I) / 2
// (the lazy walk has the same stationary vectors and cannot oscillate). Returns null if
// some vertex has no outgoing edge, or if the iteration has not converged: the
// stationary vector is then not unique or not meaningful (a disconnected graph).
export function stationary(P, { tol = 1e-12, maxIter = 20000 } = {}) {
  const n = P.length;
  if (!n || P.some((row) => row.size === 0)) return null;
  let pi = new Float64Array(n).fill(1 / n);
  for (let it = 0; it < maxIter; it++) {
    const next = new Float64Array(n);
    for (let v = 0; v < n; v++) {
      const x = pi[v] / 2;
      next[v] += x;
      for (const [w, p] of P[v]) next[w] += x * p;
    }
    let diff = 0;
    for (let v = 0; v < n; v++) diff += Math.abs(next[v] - pi[v]);
    pi = next;
    if (diff < tol) return pi;
  }
  return null;
}

// How far pi is from stationary: the largest entry of |pi P - pi|.
export function stationaryResidual(P, pi) {
  const out = new Float64Array(P.length);
  P.forEach((row, v) => { for (const [w, p] of row) out[w] += pi[v] * p; });
  let worst = 0;
  for (let v = 0; v < P.length; v++) worst = Math.max(worst, Math.abs(out[v] - pi[v]));
  return worst;
}

// pi rescaled so that its most common value is 1 (the value at a vertex with no extra
// automorphisms, for the walks here), with the vertices where it differs from 1.
export function relativeStationary(pi, digits = 6) {
  const groups = new Map();              // values equal to "digits" digits: [count, sum]
  for (const x of pi) {
    const key = x.toPrecision(digits), g = groups.get(key) || [0, 0];
    g[0] += 1; g[1] += x;
    groups.set(key, g);
  }
  let unit = pi[0], best = 0;
  for (const [c, sum] of groups.values()) if (c > best) { best = c; unit = sum / c; }
  const values = Array.from(pi, (x) => x / unit);
  const special = [];
  values.forEach((x, v) => { if (Math.abs(x - 1) > 1e-4) special.push(v); });
  return { values, special };
}

// ---- sampling ------------------------------------------------------------------

// One step from a row of P: the vertex chosen by a uniform random number u in [0, 1).
export function sampleStep(row, u) {
  let acc = 0, last = null;
  for (const [w, p] of row) {
    acc += p;
    last = w;
    if (u < acc) return w;
  }
  return last;                 // rounding: u fell in the last sliver below 1
}
