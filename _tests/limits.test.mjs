// Tests for assets/js/limits.js (the Limit laws page) and its data, assets/data/limits/.
// Run with:  node --test _tests/*.test.mjs
import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import {
  stDensity, kmDensity, kmstDistances, TV_CONST, KS_CONST, primesUpTo,
} from "../assets/js/limits.js";

const data = (name) => JSON.parse(readFileSync(new URL(`../assets/data/limits/${name}.json`, import.meta.url)));

// Integral of x^m f(x) over [-2, 2] by the midpoint rule.
function moment(f, m, n = 200000) {
  const h = 4 / n;
  let s = 0;
  for (let i = 0; i < n; i++) { const x = -2 + (i + 0.5) * h; s += x ** m * f(x) * h; }
  return s;
}

test("the laws are probability measures with the right second moments", () => {
  assert.ok(Math.abs(moment(stDensity, 0) - 1) < 1e-6);
  assert.ok(Math.abs(moment(stDensity, 2) - 1) < 1e-6);   // Sato-Tate: E[x^2] = 1
  for (const l of [2, 3, 5, 97, 1009]) {
    const f = (x) => kmDensity(l, x);
    assert.ok(Math.abs(moment(f, 0) - 1) < 1e-6, `mass of mu_${l}`);
    // closed walks of length 2 (and 4) at a vertex of the (l+1)-regular tree, over l (l^2)
    assert.ok(Math.abs(moment(f, 2) - (l + 1) / l) < 1e-6, `second moment of mu_${l}`);
    assert.ok(Math.abs(moment(f, 4) - (l + 1) * (2 * l + 1) / l ** 2) < 1e-6, `fourth moment of mu_${l}`);
    assert.ok(Math.abs(moment(f, 3)) < 1e-9, `mu_${l} is symmetric`);
  }
});

test("the distances to Sato-Tate decay like 3 sqrt 3 / (4 pi l) and 3 sqrt 3 / (8 pi l)", () => {
  const ells = primesUpTo(1009);
  assert.equal(ells.length, 169);
  assert.equal(ells.at(-1), 1009);
  let prev = Infinity;
  for (const l of ells) {
    const { tv } = kmstDistances(l);
    assert.ok(tv < prev, `total variation decreases at l = ${l}`);
    prev = tv;
  }
  for (const l of [1009, 997]) {
    const { tv, ks } = kmstDistances(l, 20000);
    assert.ok(Math.abs(tv * l / TV_CONST - 1) < 0.005, `tv at l = ${l}: ${tv * l} vs ${TV_CONST}`);
    assert.ok(Math.abs(ks * l / KS_CONST - 1) < 0.005, `ks at l = ${l}: ${ks * l} vs ${KS_CONST}`);
  }
});

test("gamma1.json: S_2(Gamma_1(101)) by character", () => {
  const d = data("gamma1");
  const N = d.N;
  assert.equal(d.chars.length, (N - 1) / 2);
  const all = d.chars.flatMap((c) => c.eig);
  assert.equal(all.length, (N - 5) * (N - 7) / 24);
  assert.equal(d.lines, 50);
  for (const c of d.chars) {
    const [cr, ci] = c.chi;
    assert.ok(Math.abs(Math.hypot(cr, ci) - 1) < 1e-5);
    for (const [re, im] of c.eig) {
      assert.ok(Math.hypot(re, im) <= 2 * Math.SQRT2 + 1e-3);
      // lambda = chi(2) conj(lambda)
      const pr = cr * re + ci * im, pi = ci * re - cr * im;
      assert.ok(Math.hypot(pr - re, pi - im) < 2e-3, `lambda = chi(2) conj(lambda) for a = ${c.a}`);
    }
  }
});

const sum = (a) => a.reduce((s, x) => s + x, 0);

test("quadratic.json, diagonal.json, weight.json: histograms count every eigenvalue", () => {
  for (const r of data("quadratic").rows) {
    assert.equal(r.N % 8, 5);
    assert.equal(sum(r.hist), r.dim);
  }
  const nss = (p) => Math.floor(p / 12) + { 1: 0, 5: 1, 7: 1, 11: 2 }[p % 12];
  for (const r of data("diagonal").rows) {
    assert.equal(r.dim, nss(r.p) - 1, `dim S_2(Gamma_0(${r.p}))`);
    assert.equal(sum(r.hist), r.dim);
    assert.ok(r.ell ** 1 <= Math.log(r.p) ** 2);
  }
  const dim1 = (k) => Math.floor(k / 12) - (k % 12 === 2 ? 1 : 0);
  for (const rows of Object.values(data("weight").ells)) {
    for (const r of rows) {
      assert.equal(r.dim, dim1(r.k), `dim S_${r.k}`);
      assert.equal(sum(r.hist), r.dim);
    }
  }
});

test("planar.json: the eigenvalues of T_2 on S_2(Gamma_1(N)), in the disc of radius 2", () => {
  const d = data("planar");
  let prev = 0;
  for (const r of d.rows) {
    assert.ok(r.N > prev, "N increases");
    prev = r.N;
    assert.equal(r.dim, (r.N - 5) * (r.N - 7) / 24, `genus of X_1(${r.N})`);
    // all of them, or a sample of d.sample when there are more
    assert.equal(r.pts.length / 2, Math.min(r.dim, d.sample), `points drawn at N = ${r.N}`);
    assert.equal(sum(r.rhist), r.dim, `radial histogram at N = ${r.N}`);
    for (let i = 0; i < r.pts.length; i += 2) assert.ok(Math.hypot(r.pts[i], r.pts[i + 1]) <= 2.001);
  }
  // the moments head for those of the limit: mean |z|^2 -> 3/2, mean z^2 -> 0
  const last = d.rows.at(-1);
  assert.ok(Math.abs(last.mom["11"] - 1.5) < 0.02);
  assert.ok(last.mom["20"] < 0.02);
});
