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

// A walk that agrees with the simple random walk except on the closed 1-balls around
// the vertices in "centres", where row(v, context) gives the transition probabilities
// out of v as a Map from vertices to probabilities.
//   context = { A, P0, centres: Set, ball: Set, aut }   (aut is passed through as given)
// With no centres, or no rule, this is the simple random walk. Every row supplied by the
// rule must be a probability vector supported on the vertices of the graph.
export function buildWalk(A, { centres = [], row = null, aut = null } = {}) {
  const P0 = simpleWalk(A);
  if (!row || !centres.length) return P0;
  const ball = closedBalls(A, centres);
  const context = { A, P0, centres: new Set(centres), ball, aut };
  return P0.map((r, v) => {
    if (!ball.has(v)) return r;
    const out = row(v, context);
    if (!(out instanceof Map)) throw new Error(`walk: no row for vertex ${v}`);
    let total = 0;
    for (const [w, x] of out) {
      if (!(Number.isInteger(w) && w >= 0 && w < A.length) || !(x >= 0)) {
        throw new Error(`walk: bad entry (${w}, ${x}) in the row of vertex ${v}`);
      }
      total += x;
    }
    if (Math.abs(total - 1) > 1e-9) throw new Error(`walk: the row of vertex ${v} sums to ${total}`);
    return out;
  });
}

// THE PETERSSON WALK. It agrees with the simple random walk at distance at least 2 from
// j = 0 and j = 1728, and differs from it on the closed 1-balls around those two
// vertices. Its rows there have NOT been confirmed by the owner of this site, so no rule
// is written here and the animation runs the simple random walk (see WALK below).
//
// To switch the animation to the Petersson walk, replace "null" by a function
//     (v, { A, P0, centres, ball, aut }) => Map(vertex -> probability)
// giving the row of each vertex v of the two closed 1-balls, add the expected
// stationary values to _tests/walk.test.mjs (the test marked "todo"), and change the
// "walk" entry of the widget's info panel in _partials/isogeny-widget.qmd.
export const PETERSSON_ROW = null;

// The walk the animation runs.
//   centres  the vertices with automorphisms other than +-1 (j = 0, 1728), or [] for
//            graphs where the walk is the simple random walk by definition.
export function walkFor(A, centres = [], aut = null) {
  return buildWalk(A, { centres, row: PETERSSON_ROW, aut });
}
export const WALK = PETERSSON_ROW ? "petersson" : "simple";

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
