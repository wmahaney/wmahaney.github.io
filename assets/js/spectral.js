// The spectral-law animation: eigenvalues of the Hecke operator T_l on S_2(Gamma_0(p))
// for a growing list of primes p, as a histogram against the Kesten-McKay law.
// Data: assets/data/isogeny/spectra.json, made by _scripts/isogeny_graphs.py: for about
// 100 primes, a 60-bin histogram of the eigenvalues over [-2 sqrt l, 2 sqrt l] (the Brandt
// matrix of G_p(l) without the Eisenstein eigenvalue l + 1) and the Kolmogorov distance
// to the Kesten-McKay law. Between two primes the bars glide from one histogram to the next.
//
// Serre (1997): for fixed l, as p -> infinity the eigenvalues of T_l on S_2(Gamma_0(p))
// become equidistributed in [-2 sqrt l, 2 sqrt l] for the Kesten-McKay law of the
// (l + 1)-regular tree,
//     f(x) = (l + 1) / (2 pi) * sqrt(4 l - x^2) / ((l + 1)^2 - x^2).
//
// Sums. The T_l commute, so each eigenform has a tuple (lambda_2, lambda_3, lambda_5),
// and the adjacency matrix of G_p({2,3}) has the eigenvalues lambda_2 + lambda_3. Serre's
// theorem holds jointly: the tuples become equidistributed for the product measure
// mu_2 x mu_3 x mu_5. Pushing that forward under the sum gives the convolution
// mu_2 * mu_3 (* mu_5), which is the limit law drawn for "2, 3" and "2, 3, 5". For
// "2, 3" a small grid also shows the joint histogram of (lambda_2, lambda_3).

//
// "ℓ → ∞" (value "inf"): no eigenvalues, but the Kesten-McKay laws themselves, rescaled
// to [-2, 2], tending to the Sato-Tate law as l runs through the primes up to 1009. That
// picture is drawn by createKmSt in limits.js, which the Limit laws page uses too.
import { createKmSt } from "./limits.js";

const SVG_NS = "http://www.w3.org/2000/svg";
const el = (name, attrs = {}, parent) => {
  const node = document.createElementNS(SVG_NS, name);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
  if (parent) parent.appendChild(node);
  return node;
};
const text = (parent, x, y, str, attrs = {}) => {
  const t = el("text", { x, y, ...attrs }, parent);
  t.textContent = str;
  return t;
};
const SUB = { 2: "₂", 3: "₃", 5: "₅" };
const DISPLAY_BINS = 30;          // the 60 stored bins, merged in pairs

export function kmDensity(q, x) {
  const r = 4 * q - x * x;
  return r <= 0 ? 0 : ((q + 1) / (2 * Math.PI)) * Math.sqrt(r) / ((q + 1) ** 2 - x * x);
}

// The choices of the menu: one Hecke operator, or a sum of them.
export const LAWS = { 2: [2], 3: [3], 5: [5], 23: [2, 3], 235: [2, 3, 5] };

// The limit law for a list of primes: { a, density(x) } on [-a, a], a = sum of 2 sqrt(l).
// One prime: the Kesten-McKay density. Several: the convolution of the densities,
// by discrete convolution of midpoint samples on a grid of step h.
export function limitLaw(ells, h = 0.01) {
  if (ells.length === 1) {
    const q = ells[0];
    return { a: 2 * Math.sqrt(q), density: (x) => kmDensity(q, x) };
  }
  let mass = null, half = 0;
  for (const q of ells) {
    const m = Math.ceil((2 * Math.sqrt(q)) / h), w = new Float64Array(2 * m);
    let sum = 0;
    for (let i = 0; i < 2 * m; i++) { w[i] = kmDensity(q, (i - m + 0.5) * h); sum += w[i]; }
    for (let i = 0; i < 2 * m; i++) w[i] /= sum;
    if (!mass) mass = w;
    else {
      const out = new Float64Array(mass.length + w.length - 1);
      for (let i = 0; i < mass.length; i++) {
        const mi = mass[i];
        if (mi) for (let j = 0; j < w.length; j++) out[i + j] += mi * w[j];
      }
      mass = out;
    }
    half += m * h;
  }
  // mass[i] belongs to the cell centred at -half + (i + k / 2) h for k factors
  const shift = -half + (ells.length / 2) * h;
  const density = (x) => {
    const t = (x - shift) / h, i = Math.floor(t), f = t - i;
    const at = (k) => (k >= 0 && k < mass.length ? mass[k] / h : 0);
    return at(i) * (1 - f) + at(i + 1) * f;
  };
  return { a: ells.reduce((acc, q) => acc + 2 * Math.sqrt(q), 0), density };
}
const opName = (ells) => ells.map((q) => `T${SUB[q]}`).join(" + ");
const lawName = (ells) => (ells.length === 1 ? "the Kesten–McKay law"
  : `the convolution ${ells.map((q) => `μ${SUB[q]}`).join(" ∗ ")}`);
const edgeName = (ells) => ells.map((q) => `2√${q}`).join("+");

// still(): true while the animation is paused (a still picture is wanted).
export function createSpectral({ svg, base, W, H, onStatus, still }) {
  const g = el("g", { class: "wm-spec" }, svg);
  g.style.display = "none";
  let data = null, ell = 2, idx = 0, t0 = 0, holdUntil = 0, active = false;
  let bars = [], shown = [], from = [], to = [], tweenStart = 0, tweenMs = 1;
  let label = {}, inset = null, law = null, ells = [2], grid = null;
  // l -> infinity: the Kesten-McKay laws against Sato-Tate (createKmSt), with its own clock
  let kmst = null, kIdx = 0, kLast = 0, kHold = 0;
  const isInf = () => ell === "inf";
  function kmstGo(i, now, ms) {
    kIdx = i; kLast = now;
    kmst.go(i, now, ms, still());
  }
  function kmstStart(now) {
    kmst = createKmSt(g, { onStatus });
    kHold = 0;
    kmstGo(still() ? kmst.count - 1 : 0, now, 600);
  }

  async function ensure() {
    if (!data) {
      const r = await fetch(`${base}/spectra.json`);
      if (!r.ok) throw new Error(r.status);
      data = await r.json();
    }
  }

  // plot geometry
  const L = 56, R = W - 22, T = 128, B = H - 54;
  let X, Y, a, w;

  function heights(k) {
    const row = data.primes[k], counts = row.hist[ell], n = row.dim, out = [];
    const per = counts.length / DISPLAY_BINS;
    for (let b = 0; b < DISPLAY_BINS; b++) {
      let c = 0;
      for (let i = 0; i < per; i++) c += counts[b * per + i];
      out.push(c / (n * w));
    }
    return out;
  }

  // Axes, curve and empty bars: once per choice of l.
  function scaffold() {
    const q = ell;
    ells = LAWS[q] || [2];
    law = limitLaw(ells);
    a = law.a;
    w = (2 * a) / DISPLAY_BINS;
    const xmax = a * 1.1;
    X = (x) => L + ((x + xmax) / (2 * xmax)) * (R - L);
    let peak = 0;
    for (let i = 0; i <= 400; i++) peak = Math.max(peak, law.density(-a + (2 * a * i) / 400));
    const ymax = peak * 1.75;
    Y = (y) => B - (Math.min(y, ymax) / ymax) * (B - T);
    g.replaceChildren();
    bars = [];
    for (let k = 0; k < DISPLAY_BINS; k++) {
      const x0 = X(-a + k * w), x1 = X(-a + (k + 1) * w);
      bars.push(el("rect", { class: "wm-spec-bar", x: (x0 + 1).toFixed(1), y: B,
        width: Math.max(0.5, x1 - x0 - 2).toFixed(1), height: 0, rx: 2 }, g));
    }
    let d = "";
    for (let i = 0; i <= 300; i++) {
      const x = -a + (2 * a * i) / 300;
      d += `${i ? "L" : "M"}${X(x).toFixed(1)},${Y(law.density(x)).toFixed(1)}`;
    }
    el("path", { class: "wm-spec-curve", d }, g);
    el("line", { class: "wm-spec-axis", x1: L, x2: R, y1: B, y2: B }, g);
    const tickStep = a > 8 ? 2 : 1, tmax = Math.floor(xmax / tickStep) * tickStep;
    for (let t = -tmax; t <= tmax; t += tickStep) {
      el("line", { class: "wm-spec-tick", x1: X(t), x2: X(t), y1: B, y2: B + 5 }, g);
      text(g, X(t), B + 18, t < 0 ? `−${-t}` : String(t), { class: "wm-spec-label", "text-anchor": "middle" });
    }
    for (const s of [-1, 1]) {
      el("line", { class: "wm-spec-edge", x1: X(s * a), x2: X(s * a), y1: T - 6, y2: B }, g);
      text(g, X(s * a), T - 10, `${s < 0 ? "−" : ""}${ells.length > 1 ? "(" : ""}${edgeName(ells)}${ells.length > 1 ? ")" : ""}`,
        { class: "wm-spec-label", "text-anchor": "middle" });
    }
    text(g, (L + R) / 2, H - 14, `eigenvalue λ of ${opName(ells)}`, { class: "wm-spec-label", "text-anchor": "middle" });
    label = {
      p: text(g, L, 34, "", { class: "wm-spec-big" }),
      dim: text(g, L, 56, "", { class: "wm-spec-label" }),
      ks: text(g, L, 76, "", { class: "wm-spec-label" }),
      joint: text(g, L, 96, "", { class: "wm-spec-label" }),
    };
    // "2, 3": the joint histogram of (λ₂, λ₃), in the empty corner past the left edge.
    grid = null;
    const nb = data.jointBins;
    if (q === 23 && nb && data.primes[0].joint23) {
      const side = 84, cell = side / nb, gx = X(-a) + 12, gy = T + 20;
      text(g, gx, gy - 6, "joint law of (λ₂, λ₃)", { class: "wm-spec-small" });
      el("rect", { class: "wm-spec-inset", x: gx - 1, y: gy - 1, width: side + 2, height: side + 2 }, g);
      grid = [];
      for (let i = 0; i < nb; i++) {          // i: bin of λ₂ (across), j: bin of λ₃ (upwards)
        for (let j = 0; j < nb; j++) {
          grid.push(el("rect", { class: "wm-spec-cell", x: (gx + i * cell).toFixed(1),
            y: (gy + (nb - 1 - j) * cell).toFixed(1), width: cell.toFixed(1), height: cell.toFixed(1), opacity: 0 }, g));
        }
      }
      text(g, gx + side / 2, gy + side + 12, "λ₂ →", { class: "wm-spec-small", "text-anchor": "middle" });
      text(g, gx + side + 5, gy + side / 2, "λ₃ ↑", { class: "wm-spec-small" });
      label.corr = text(g, gx, gy + side + 26, "", { class: "wm-spec-small" });
    }
    // inset: the distance against p (log scale), for the primes shown so far
    const ix0 = W - 250, ix1 = W - 76, iy0 = 22, iy1 = 76;   // clear of the eye icon
    const lp = Math.log, pmin = lp(data.primes[0].p), pmax = lp(data.primes[data.primes.length - 1].p);
    const dmax = Math.max(...data.primes.map((row) => row.ks[q])) * 1.1;
    inset = {
      IX: (p) => ix0 + ((lp(p) - pmin) / (pmax - pmin)) * (ix1 - ix0),
      IY: (v) => iy1 - (v / dmax) * (iy1 - iy0),
    };
    el("rect", { class: "wm-spec-inset", x: ix0 - 34, y: iy0 - 12, width: ix1 - ix0 + 44, height: iy1 - iy0 + 34, rx: 8 }, g);
    text(g, ix0 - 26, iy0 - 1, "distance", { class: "wm-spec-small" });
    text(g, ix1, iy1 + 16, "p →", { class: "wm-spec-small", "text-anchor": "end" });
    text(g, ix0 - 26, iy1 + 16, `${data.primes[0].p}`, { class: "wm-spec-small" });
    inset.trace = el("path", { class: "wm-spec-trace" }, g);
    inset.dot = el("circle", { class: "wm-spec-dot wm-spec-now", r: 3.6 }, g);
    shown = new Array(DISPLAY_BINS).fill(0);
  }

  function paintBars(hs) {
    hs.forEach((h, k) => {
      const y = Y(h);
      bars[k].setAttribute("y", y.toFixed(1));
      bars[k].setAttribute("height", (B - y).toFixed(1));
    });
  }

  // Move to prime number k: labels and inset now, bars by a short glide.
  function goTo(k, now, ms) {
    idx = k;
    const row = data.primes[k], D = row.ks[ell];
    label.p.textContent = `p = ${row.p}`;
    label.dim.textContent = `${row.dim} eigenvalues = dim S₂(Γ₀(p))`;
    label.ks.textContent = `Kolmogorov distance ${D.toFixed(3)}`;
    if (ells.length > 1 && row.m2) {
      const limit = ells.reduce((acc, q) => acc + q + 1, 0);
      label.joint.textContent = `mean of λ² = ${row.m2[ell].toFixed(2)} (limit ${limit})`;
    } else label.joint.textContent = "";
    if (label.corr) {
      label.corr.textContent = row.corr && row.corr["23"] != null ? `correlation ${row.corr["23"].toFixed(3)}` : "";
    }
    if (grid && row.joint23) {
      const top = Math.max(1, ...row.joint23);
      row.joint23.forEach((c, k) => grid[k].setAttribute("opacity", (0.9 * c / top).toFixed(3)));
    }
    let path = "";
    for (let i = 0; i <= k; i++) {
      const r = data.primes[i];
      path += `${i ? "L" : "M"}${inset.IX(r.p).toFixed(1)},${inset.IY(r.ks[ell]).toFixed(1)}`;
    }
    inset.trace.setAttribute("d", path);
    inset.dot.setAttribute("cx", inset.IX(row.p).toFixed(1));
    inset.dot.setAttribute("cy", inset.IY(D).toFixed(1));
    from = shown.slice();
    to = heights(k);
    tweenStart = now;
    tweenMs = Math.max(1, ms * 0.85);
    if (still()) { shown = to.slice(); paintBars(shown); }
    const last = k === data.primes.length - 1;
    onStatus(`${opName(ells)} on S₂(Γ₀(${row.p})): ${row.dim} eigenvalues. ` +
      `Kolmogorov distance to ${lawName(ells)}: ${D.toFixed(3)}` +
      (last ? " (largest p shown; the animation starts again)." : "."));
  }

  return {
    async show(l) {
      ell = l === "inf" ? "inf" : Number(l);
      active = true;
      if (isInf()) {
        g.style.display = "";
        kmstStart(performance.now());
        return;
      }
      try { await ensure(); } catch {
        onStatus("Could not load the eigenvalues. Please try again.");
        return;
      }
      if (!active) return;
      g.style.display = "";
      scaffold();
      t0 = performance.now(); holdUntil = 0;
      // Paused: the largest prime, where the histogram is closest to the law.
      goTo(still() ? data.primes.length - 1 : 0, t0, 600);
    },
    // Play pressed: run through the primes from the first one.
    restart() {
      if (active && isInf() && kmst) { kmst.reset(); kHold = 0; kmstGo(0, performance.now(), 600); return; }
      if (!active || !data) return;
      shown = new Array(DISPLAY_BINS).fill(0);
      t0 = performance.now(); holdUntil = 0;
      goTo(0, t0, 600);
    },
    hide() { active = false; g.style.display = "none"; },
    setEll(l) {
      ell = l === "inf" ? "inf" : Number(l);
      if (active && isInf()) { kmstStart(performance.now()); return; }
      if (active && !data) { this.show(l); return; }
      if (active && data) {
        scaffold();
        goTo(idx, performance.now(), 600);
        if (still()) paintBars(shown);
      }
    },
    step(now, ms) {
      if (active && isInf() && kmst && !still()) {
        kmst.glide(now);
        if (kHold) {
          if (now >= kHold) { kHold = 0; kmst.reset(); kmstGo(0, now, ms); }
        } else if (now - kLast >= ms / 3) {   // 169 primes: three per period of the histograms
          if (kIdx < kmst.count - 1) {
            kmstGo(kIdx + 1, now, ms / 3);
            if (kIdx === kmst.count - 1) kHold = now + 8 * ms;
          }
        }
        return;
      }
      if (!active || !data || still() || isInf()) return;
      // the glide between two histograms
      const t = Math.min(1, (now - tweenStart) / tweenMs), e = t * t * (3 - 2 * t);
      shown = from.map((f, k) => f + (to[k] - f) * e);
      paintBars(shown);
      if (holdUntil) {
        if (now < holdUntil) return;
        holdUntil = 0; t0 = now;
        goTo(0, now, ms);
        return;
      }
      if (now - t0 < ms) return;
      t0 = now;
      if (idx < data.primes.length - 1) {
        goTo(idx + 1, now, ms);
        if (idx === data.primes.length - 1) holdUntil = now + 8 * ms;
      }
    },
  };
}
