// Tests for assets/js/walk.js, on the graph data the site ships.
//     node --test _tests/*.test.mjs
// They need nothing but Node (18 or newer).
import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync, readdirSync } from "node:fs";
import { fileURLToPath } from "node:url";
import path from "node:path";
import {
  adjacency, rowSum, simpleWalk, closedBalls, peterssonWalk, isogenyCounts, weightFrom,
  walkFor, walkName, stationary, stationaryResidual, relativeStationary, sampleStep, linkCounts,
} from "../assets/js/walk.js";

const DATA = path.join(path.dirname(fileURLToPath(import.meta.url)), "..", "assets", "data", "isogeny");
const load = (name) => JSON.parse(readFileSync(path.join(DATA, name), "utf8")).graph;
const DEGREE = { 2: 3, 23: 7, 235: 13 };                 // sum of l + 1 over the l used
const BASE = readdirSync(DATA).filter((f) => /^p\d+-(2|23|235)\.json$/.test(f));
// Small enough to test every vertex quickly.
const SMALL = BASE.filter((f) => Number(/^p(\d+)/.exec(f)[1]) <= 2000);
const centresOf = (g) => Object.keys(g.aut || {}).map(Number);
const close = (a, b, eps = 1e-9) => Math.abs(a - b) <= eps;

test("the data: every row of A sums to the number of kernels", () => {
  assert.ok(BASE.length > 0, "no graph data found");
  for (const f of BASE) {
    const g = load(f), kind = /-(\d+)\.json$/.exec(f)[1];
    const A = adjacency(g.nodes.length, g.links);
    A.forEach((row, v) => assert.equal(rowSum(row), DEGREE[kind], `${f}, vertex ${v}`));
  }
});

test("the data: A[v][w] |Aut w| = A[w][v] |Aut v|, and A is symmetric away from j = 0, 1728", () => {
  for (const f of SMALL) {
    const g = load(f), n = g.nodes.length, A = adjacency(n, g.links);
    const aut = (v) => (g.aut && g.aut[v]) || 2;
    for (let v = 0; v < n; v++) {
      for (const [w, a] of A[v]) {
        assert.equal(a * aut(w), (A[w].get(v) || 0) * aut(v), `${f}: ${v} -> ${w}`);
      }
    }
  }
});

test("the data: the vertices with extra automorphisms are j = 0 and j = 1728", () => {
  for (const f of SMALL) {
    const g = load(f), p = Number(/^p(\d+)/.exec(f)[1]);
    for (const [v, order] of Object.entries(g.aut || {})) {
      const j = Number(g.nodes[v]);
      assert.ok((order === 6 && j === 0) || (order === 4 && j === 1728 % p), `${f}: vertex ${v}`);
    }
    assert.equal(g.nodes.includes("0"), p % 3 === 2, `${f}: j = 0 is supersingular iff p = 2 mod 3`);
    assert.equal(g.nodes.includes(String(1728 % p)), p % 4 === 3, `${f}: j = 1728 iff p = 3 mod 4`);
  }
});

test("every row of P sums to 1", () => {
  for (const f of SMALL) {
    const g = load(f), A = adjacency(g.nodes.length, g.links);
    const P = walkFor(A, g.aut);
    P.forEach((row, v) => assert.ok(close(rowSum(row), 1), `${f}, vertex ${v}`));
  }
});

test("P agrees with the simple random walk outside the closed 1-balls of j = 0 and j = 1728", () => {
  for (const f of SMALL) {
    const g = load(f), A = adjacency(g.nodes.length, g.links);
    const centres = centresOf(g), ball = closedBalls(A, centres);
    const P = walkFor(A, g.aut), P0 = simpleWalk(A);
    P.forEach((row, v) => {
      if (ball.has(v)) return;
      assert.equal(row.size, P0[v].size, `${f}, vertex ${v}`);
      for (const [w, x] of row) assert.ok(close(x, P0[v].get(w)), `${f}: ${v} -> ${w}`);
    });
  }
});

test("Petersson walk: P[i][j] = A[i][j] w_j / s_i, and the matrix A[i][j] w_j is symmetric", () => {
  for (const f of SMALL) {
    const g = load(f), n = g.nodes.length, A = adjacency(n, g.links);
    const w = weightFrom(g.aut), s = isogenyCounts(A, w), P = walkFor(A, g.aut);
    assert.equal(walkName(g.aut), "petersson");
    for (let v = 0; v < n; v++) {
      for (const [u, a] of A[v]) {
        assert.equal(a * w(u), (A[u].get(v) || 0) * w(v), `${f}: ${v} -> ${u}`);
        assert.ok(close(P[v].get(u), (a * w(u)) / s[v]), `${f}: ${v} -> ${u}`);
      }
    }
  }
});

test("Petersson walk: the rows near j = 0 and j = 1728 in G_431(2)", () => {
  // Phi_2(0, Y) = (Y - 54000)^3 and Phi_2(1728, Y) = (Y - 1728)(Y - 287496)^2;
  // modulo 431: 54000 = 125, 1728 = 4, 287496 = 19.
  const g = load("p431-2.json"), A = adjacency(g.nodes.length, g.links);
  const at = (j) => g.nodes.indexOf(String(j));
  const [v0, v125, v1728, v19] = [0, 125, 4, 19].map(at);
  const P = walkFor(A, g.aut);
  assert.deepEqual([...A[v0]], [[v125, 3]]);
  assert.deepEqual([...P[v0]], [[v125, 1]]);                       // one neighbour
  assert.ok(close(P[v125].get(v0), 3 / 5));                        // 3 of the 5 isogenies
  for (const [u, x] of P[v125]) if (u !== v0) assert.ok(close(x, A[v125].get(u) / 5));
  assert.ok(close(P[v1728].get(v1728), 1 / 2));                    // the loop 1 + i, w = 2
  assert.ok(close(P[v1728].get(v19), 1 / 2));
  assert.ok(close(P[v19].get(v1728), 1 / 2));                      // 2 of the 4 isogenies
  const { values, special } = relativeStationary(stationary(P));
  assert.deepEqual(special, [v0, v125, v1728, v19].filter((v) => !close(values[v], 1, 1e-7)).sort((a, b) => a - b));
  assert.ok(close(values[v0], 1, 1e-7));
  assert.ok(close(values[v125], 5 / 3, 1e-7));
  assert.ok(close(values[v1728], 4 / 3, 1e-7));
  assert.ok(close(values[v19], 4 / 3, 1e-7));
});

test("Petersson walk: with no extra automorphisms (p = 1 mod 12) it is the simple random walk", () => {
  const files = SMALL.filter((f) => Number(/^p(\d+)/.exec(f)[1]) % 12 === 1);
  assert.ok(files.length > 0);
  for (const f of files) {
    const g = load(f), A = adjacency(g.nodes.length, g.links);
    assert.deepEqual(centresOf(g), []);
    assert.deepEqual(walkFor(A, g.aut).map((r) => [...r]), simpleWalk(A).map((r) => [...r]));
  }
});

test("without automorphism data the walk is the simple random walk", () => {
  const g = load("p431-2.json"), A = adjacency(g.nodes.length, g.links);
  assert.equal(walkName(null), "simple");
  assert.deepEqual(walkFor(A).map((r) => [...r]), simpleWalk(A).map((r) => [...r]));
});

test("the stationary vector solves pi P = pi", () => {
  for (const f of SMALL) {
    const g = load(f), A = adjacency(g.nodes.length, g.links);
    const P = walkFor(A, g.aut), pi = stationary(P);
    assert.ok(pi, `${f}: no stationary vector`);
    assert.ok(close(pi.reduce((a, b) => a + b, 0), 1));
    assert.ok(stationaryResidual(P, pi) < 1e-10, f);
  }
});

test("simple random walk: pi is proportional to 1 / |Aut| (1 generic, 1/2 at 1728, 1/3 at 0)", () => {
  for (const f of SMALL) {
    const g = load(f), A = adjacency(g.nodes.length, g.links);
    const { values, special } = relativeStationary(stationary(simpleWalk(A)));
    values.forEach((x, v) => assert.ok(close(x, 2 / ((g.aut && g.aut[v]) || 2), 1e-7), `${f}, vertex ${v}`));
    assert.deepEqual(special, centresOf(g).sort((a, b) => a - b), f);
  }
});

test("Petersson walk: pi is proportional to s_i, the isogenies out of E_i counted up to sign", () => {
  for (const f of SMALL) {
    const g = load(f), kind = /-(\d+)\.json$/.exec(f)[1], A = adjacency(g.nodes.length, g.links);
    const w = weightFrom(g.aut), s = isogenyCounts(A, w), total = s.reduce((a, b) => a + b, 0);
    const P = peterssonWalk(A, w), pi = stationary(P);
    pi.forEach((x, v) => assert.ok(close(x, s[v] / total, 1e-9), `${f}, vertex ${v}`));
    // s_i = (number of kernels) + 2 A[i][j=0] + A[i][j=1728]: 1 away from the two 1-balls
    const ball = closedBalls(A, centresOf(g));
    s.forEach((x, v) => {
      let extra = 0;
      for (const [u, a] of A[v]) extra += a * (w(u) - 1);
      assert.equal(x, DEGREE[kind] + extra);
      if (!ball.has(v)) assert.equal(x, DEGREE[kind], `${f}, vertex ${v}`);
    });
  }
});

// ---- level structure ----------------------------------------------------------

// The level-structure graphs and the number of level structures on one curve: the
// mass sum over v of 1/|Aut(v)| is (p - 1)/24 times that number.
const LEVEL_KINDS = { "2b5": 6, "23b5": 6, "2g15": 12, "3g2": 6, "3g4": 96, "5g4": 96 };
const LEVEL = readdirSync(DATA).filter((f) => {
  const m = /^p(\d+)-(\w+)\.json$/.exec(f);
  return m && LEVEL_KINDS[m[2]] && Number(m[1]) <= 1100;
});
const kindOf = (f) => /-(\w+)\.json$/.exec(f)[1];
const primeOf = (f) => Number(/^p(\d+)/.exec(f)[1]);
const orderOf = (g) => (v) => (g.aut[v] !== undefined ? g.aut[v] : g.autDefault || 2);

test("level data: every graph has its automorphism orders, and they satisfy the mass formula", () => {
  assert.ok(LEVEL.length > 0, "no level-structure data found");
  for (const f of LEVEL) {
    const g = load(f), n = g.nodes.length, aut = orderOf(g);
    assert.ok(g.aut, `${f}: no aut field`);
    let mass = 0;
    for (let v = 0; v < n; v++) mass += 1 / aut(v);
    assert.ok(close(mass, ((primeOf(f) - 1) / 24) * LEVEL_KINDS[kindOf(f)], 1e-9), `${f}: mass ${mass}`);
  }
});

test("level data: the Brandt relation, twisted by <2> on Gamma_1(5)", () => {
  for (const f of LEVEL) {
    const g = load(f), n = g.nodes.length, A = adjacency(n, g.links, !!g.directed), w = weightFrom(g.aut, g.autDefault);
    if (kindOf(f) === "2g15") {
      // The column sums of A w: sum over v of A[v][u] w_u = s_u, the identity behind pi ~ s.
      const s = isogenyCounts(A, w), col = new Array(n).fill(0);
      A.forEach((row) => { for (const [u, a] of row) col[u] += a * w(u); });
      col.forEach((x, u) => assert.ok(close(x, s[u]), `${f}: column ${u}`));
      continue;
    }
    for (let v = 0; v < n; v++) {
      for (const [u, a] of A[v]) assert.equal(a * w(u), (A[u].get(v) || 0) * w(v), `${f}: ${v} -> ${u}`);
    }
  }
});

test("level Petersson walk: rows sum to 1, and pi is proportional to s on the connected graphs", () => {
  for (const f of LEVEL) {
    const g = load(f), n = g.nodes.length, A = adjacency(n, g.links, !!g.directed);
    const w = weightFrom(g.aut, g.autDefault), P = walkFor(A, g.aut, g.autDefault);
    assert.equal(walkName(g.aut), "petersson");
    P.forEach((row, v) => assert.ok(close(rowSum(row), 1), `${f}, vertex ${v}`));
    if ((g.components || 1) > 1) continue;
    const s = isogenyCounts(A, w), total = s.reduce((a, b) => a + b, 0), pi = stationary(P);
    assert.ok(pi, `${f}: no stationary vector`);
    pi.forEach((x, v) => assert.ok(close(x, s[v] / total, 1e-8), `${f}, vertex ${v}`));
  }
});

test("level Petersson walk: equal to the simple walk where every vertex has the generic order", () => {
  for (const f of LEVEL) {
    const g = load(f), A = adjacency(g.nodes.length, g.links, !!g.directed);
    if (Object.keys(g.aut).length) continue;
    assert.deepEqual(walkFor(A, g.aut, g.autDefault).map((r) => [...r]), simpleWalk(A).map((r) => [...r]), f);
  }
});

test("sampling follows the row", () => {
  const row = new Map([[5, 0.25], [7, 0.5], [9, 0.25]]);
  assert.equal(sampleStep(row, 0), 5);
  assert.equal(sampleStep(row, 0.2499), 5);
  assert.equal(sampleStep(row, 0.25), 7);
  assert.equal(sampleStep(row, 0.7499), 7);
  assert.equal(sampleStep(row, 0.75), 9);
  assert.equal(sampleStep(row, 0.999999), 9);
});

test("links: both counts, loops, and directed graphs", () => {
  assert.deepEqual(linkCounts([0, 7, 2, 3, 1]), { s: 0, t: 7, ell: 2, out: 3, back: 1 });
  assert.deepEqual(linkCounts([1, 2, 2, 2]), { s: 1, t: 2, ell: 2, out: 2, back: 2 });
  assert.deepEqual(linkCounts([1, 1, 2, 1]), { s: 1, t: 1, ell: 2, out: 1, back: 0 });
  assert.deepEqual(linkCounts([1, 2, 2, 1], true), { s: 1, t: 2, ell: 2, out: 1, back: 0 });
  const A = adjacency(3, [[0, 1, 2, 3, 1], [1, 1, 2, 1], [1, 2, 2, 1]]);
  assert.deepEqual([...A[0]], [[1, 3]]);
  assert.deepEqual([...A[1]].sort(), [[0, 1], [1, 1], [2, 1]]);
});

test("a disconnected graph has no single stationary vector to report", () => {
  const files = readdirSync(DATA).filter((f) => /^p\d+-5g4\.json$/.test(f));
  if (!files.length) return;
  const g = load(files[0]), A = adjacency(g.nodes.length, g.links);
  const pi = stationary(simpleWalk(A));
  // It converges (to the uniform vector, by symmetry), so the page decides by the
  // number of components, not by this function.
  assert.ok(pi === null || stationaryResidual(simpleWalk(A), pi) < 1e-10);
});
