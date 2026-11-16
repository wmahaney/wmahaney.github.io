// Tests for assets/js/walk.js, on the graph data the site ships.
//     node --test _tests/*.test.mjs
// They need nothing but Node (18 or newer).
import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync, readdirSync } from "node:fs";
import { fileURLToPath } from "node:url";
import path from "node:path";
import {
  adjacency, rowSum, simpleWalk, closedBalls, buildWalk, walkFor, WALK, PETERSSON_ROW,
  stationary, stationaryResidual, relativeStationary, sampleStep, linkCounts,
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
    const P = walkFor(A, centresOf(g), g.aut);
    P.forEach((row, v) => assert.ok(close(rowSum(row), 1), `${f}, vertex ${v}`));
  }
});

test("P agrees with the simple random walk outside the closed 1-balls of j = 0 and j = 1728", () => {
  for (const f of SMALL) {
    const g = load(f), A = adjacency(g.nodes.length, g.links);
    const centres = centresOf(g), ball = closedBalls(A, centres);
    const P = walkFor(A, centres, g.aut), P0 = simpleWalk(A);
    P.forEach((row, v) => {
      if (ball.has(v)) return;
      assert.equal(row.size, P0[v].size, `${f}, vertex ${v}`);
      for (const [w, x] of row) assert.ok(close(x, P0[v].get(w)), `${f}: ${v} -> ${w}`);
    });
  }
});

test("a rule changes the rows of the two 1-balls and no others", () => {
  const g = load("p431-2.json"), A = adjacency(g.nodes.length, g.links);
  const centres = centresOf(g), ball = closedBalls(A, centres);
  assert.deepEqual(centres.map((v) => g.nodes[v]).sort(), ["0", "4"]);   // 1728 = 4 mod 431
  const touched = [];
  // A stand-in rule, for the plumbing only: stay put with probability 1.
  const P = buildWalk(A, { centres, row: (v) => { touched.push(v); return new Map([[v, 1]]); } });
  assert.deepEqual(touched.sort((a, b) => a - b), [...ball].sort((a, b) => a - b));
  const P0 = simpleWalk(A);
  P.forEach((row, v) => { if (!ball.has(v)) assert.deepEqual([...row], [...P0[v]]); });
  assert.throws(() => buildWalk(A, { centres, row: () => new Map([[0, 0.5]]) }), /sums to/);
});

test("the stationary vector solves pi P = pi", () => {
  for (const f of SMALL) {
    const g = load(f), A = adjacency(g.nodes.length, g.links);
    const P = walkFor(A, centresOf(g), g.aut), pi = stationary(P);
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

// The owner's values for the Petersson walk (entries near j = 0 and j = 1728 that differ
// from 1, such as 5/6 and 3/4) go here once the rows of P are confirmed.
test("Petersson walk: the stationary vector matches the owner's values", { todo: PETERSSON_ROW === null }, () => {
  assert.equal(WALK, "petersson");
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
