// The spectral-law animation: eigenvalues of the Hecke operator T_l on S_2(Gamma_0(p))
// for a growing list of primes p, as a histogram against the Kesten-McKay law.
// Data: assets/data/isogeny/spectra.json, made by _scripts/isogeny_graphs.py
// (the Brandt matrices of G_p(l) with the Eisenstein eigenvalue l + 1 removed).
//
// Serre (1997): for fixed l, as p -> infinity the eigenvalues of T_l on S_2(Gamma_0(p))
// become equidistributed in [-2 sqrt l, 2 sqrt l] for the Kesten-McKay law of the
// (l + 1)-regular tree,
//     f(x) = (l + 1) / (2 pi) * sqrt(4 l - x^2) / ((l + 1)^2 - x^2).

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

export function kmDensity(q, x) {
  const r = 4 * q - x * x;
  return r <= 0 ? 0 : ((q + 1) / (2 * Math.PI)) * Math.sqrt(r) / ((q + 1) ** 2 - x * x);
}

// The Kesten-McKay distribution function on a fine grid (trapezoid rule, normalised).
function kmCdf(q) {
  const a = 2 * Math.sqrt(q), m = 4000, xs = [], cs = [0];
  for (let i = 0; i <= m; i++) xs.push(-a + (2 * a * i) / m);
  for (let i = 1; i <= m; i++) {
    cs.push(cs[i - 1] + ((kmDensity(q, xs[i - 1]) + kmDensity(q, xs[i])) / 2) * (xs[i] - xs[i - 1]));
  }
  const total = cs[m];
  return (x) => {
    if (x <= -a) return 0;
    if (x >= a) return 1;
    const t = ((x + a) / (2 * a)) * m, i = Math.min(m - 1, Math.floor(t));
    return (cs[i] + (cs[i + 1] - cs[i]) * (t - i)) / total;
  };
}

// sup_x |F_empirical(x) - F_KM(x)|
function kolmogorov(evs, F) {
  const s = [...evs].sort((u, v) => u - v), n = s.length;
  let d = 0;
  s.forEach((x, i) => { d = Math.max(d, Math.abs((i + 1) / n - F(x)), Math.abs(i / n - F(x))); });
  return d;
}

export function createSpectral({ svg, base, W, H, onStatus, reduced }) {
  const g = el("g", { class: "wm-spec" }, svg);
  g.style.display = "none";
  let data = null, ell = 2, idx = 0, t0 = 0, holdUntil = 0, active = false;
  const cdfs = {}, dists = {};

  async function ensure() {
    if (!data) {
      const r = await fetch(`${base}/spectra.json`);
      if (!r.ok) throw new Error(r.status);
      data = await r.json();
    }
  }
  function distances(l) {
    if (!dists[l]) {
      cdfs[l] = cdfs[l] || kmCdf(l);
      dists[l] = data.map((row) => kolmogorov(row.eigenvalues[l], cdfs[l]));
    }
    return dists[l];
  }

  // plot geometry
  const L = 56, R = W - 22, T = 128, B = H - 54;
  function draw() {
    const q = ell, a = 2 * Math.sqrt(q), row = data[idx], evs = row.eigenvalues[q];
    const n = evs.length, D = distances(q);
    const xmax = a * 1.1, X = (x) => L + ((x + xmax) / (2 * xmax)) * (R - L);
    let peak = 0;
    for (let i = 0; i <= 400; i++) peak = Math.max(peak, kmDensity(q, -a + (2 * a * i) / 400));
    const ymax = peak * 1.75, Y = (y) => B - (Math.min(y, ymax) / ymax) * (B - T);
    g.replaceChildren();

    // histogram of the eigenvalues, bins over [-2 sqrt l, 2 sqrt l]
    const nb = Math.max(10, Math.min(36, Math.round(1.7 * Math.sqrt(n)))), w = (2 * a) / nb;
    const counts = new Array(nb).fill(0);
    for (const x of evs) counts[Math.min(nb - 1, Math.max(0, Math.floor((x + a) / w)))] += 1;
    counts.forEach((c, k) => {
      if (!c) return;
      const h = c / (n * w), x0 = X(-a + k * w), x1 = X(-a + (k + 1) * w);
      el("rect", { class: "wm-spec-bar", x: (x0 + 1).toFixed(1), y: Y(h).toFixed(1),
        width: Math.max(0.5, x1 - x0 - 2).toFixed(1), height: (B - Y(h)).toFixed(1), rx: 2 }, g);
    });
    // the Kesten-McKay density
    let d = "";
    for (let i = 0; i <= 300; i++) {
      const x = -a + (2 * a * i) / 300;
      d += `${i ? "L" : "M"}${X(x).toFixed(1)},${Y(kmDensity(q, x)).toFixed(1)}`;
    }
    el("path", { class: "wm-spec-curve", d }, g);
    // axes and ticks
    el("line", { class: "wm-spec-axis", x1: L, x2: R, y1: B, y2: B }, g);
    for (let t = -Math.floor(xmax); t <= Math.floor(xmax); t++) {
      el("line", { class: "wm-spec-tick", x1: X(t), x2: X(t), y1: B, y2: B + 5 }, g);
      text(g, X(t), B + 18, t < 0 ? `−${-t}` : String(t), { class: "wm-spec-label", "text-anchor": "middle" });
    }
    for (const s of [-1, 1]) {
      el("line", { class: "wm-spec-edge", x1: X(s * a), x2: X(s * a), y1: T - 6, y2: B }, g);
      text(g, X(s * a), T - 10, `${s < 0 ? "−" : ""}2√${q}`, { class: "wm-spec-label", "text-anchor": "middle" });
    }
    text(g, (L + R) / 2, H - 14, `eigenvalue λ of T${"₂₃₅"[[2, 3, 5].indexOf(q)]}`, {
      class: "wm-spec-label", "text-anchor": "middle" });

    // header: p and the size of the space
    text(g, L, 34, `p = ${row.p}`, { class: "wm-spec-big" });
    text(g, L, 56, `${n} eigenvalues = dim S₂(Γ₀(p))`, { class: "wm-spec-label" });
    text(g, L, 76, `Kolmogorov distance ${D[idx].toFixed(3)}`, { class: "wm-spec-label" });

    // inset: the distance against p, for the primes shown so far
    const ix0 = W - 250, ix1 = W - 76, iy0 = 22, iy1 = 76;   // clear of the eye icon
    const lp = (p) => Math.log(p), pmin = lp(data[0].p), pmax = lp(data[data.length - 1].p);
    const dmax = Math.max(...D) * 1.1;
    const IX = (p) => ix0 + ((lp(p) - pmin) / (pmax - pmin)) * (ix1 - ix0);
    const IY = (v) => iy1 - (v / dmax) * (iy1 - iy0);
    el("rect", { class: "wm-spec-inset", x: ix0 - 34, y: iy0 - 12, width: ix1 - ix0 + 44, height: iy1 - iy0 + 34, rx: 8 }, g);
    text(g, ix0 - 26, iy0 - 1, "distance", { class: "wm-spec-small" });
    text(g, ix1, iy1 + 16, "p →", { class: "wm-spec-small", "text-anchor": "end" });
    let path = "";
    for (let k = 0; k <= idx; k++) path += `${k ? "L" : "M"}${IX(data[k].p).toFixed(1)},${IY(D[k]).toFixed(1)}`;
    el("path", { class: "wm-spec-trace", d: path }, g);
    for (let k = 0; k <= idx; k++) {
      el("circle", { class: k === idx ? "wm-spec-dot wm-spec-now" : "wm-spec-dot",
        cx: IX(data[k].p).toFixed(1), cy: IY(D[k]).toFixed(1), r: k === idx ? 3.6 : 2.2 }, g);
    }
    text(g, ix0 - 26, iy1 + 16, `${data[0].p}`, { class: "wm-spec-small" });

    onStatus(`T${"₂₃₅"[[2, 3, 5].indexOf(q)]} on S₂(Γ₀(${row.p})): ${n} eigenvalues. ` +
      `Kolmogorov distance to the Kesten–McKay law: ${D[idx].toFixed(3)}` +
      (idx === data.length - 1 ? " (largest p shown; the animation starts again)." : "."));
  }

  return {
    async show(l) {
      ell = Number(l);
      active = true;
      try { await ensure(); } catch {
        onStatus("Could not load the eigenvalues. Please try again.");
        return;
      }
      if (!active) return;
      g.style.display = "";
      idx = reduced ? data.length - 1 : 0;
      t0 = performance.now(); holdUntil = 0;
      draw();
    },
    hide() { active = false; g.style.display = "none"; },
    setEll(l) { ell = Number(l); if (active && data) draw(); },
    step(now, ms) {
      if (!active || !data || reduced) return;
      if (holdUntil) {
        if (now < holdUntil) return;
        holdUntil = 0; idx = 0; t0 = now; draw();
        return;
      }
      if (now - t0 < ms) return;
      t0 = now;
      if (idx < data.length - 1) {
        idx += 1;
        draw();
        if (idx === data.length - 1) holdUntil = now + 4 * ms;
      }
    },
  };
}
