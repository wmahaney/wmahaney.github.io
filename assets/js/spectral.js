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

export function createSpectral({ svg, base, W, H, onStatus, reduced }) {
  const g = el("g", { class: "wm-spec" }, svg);
  g.style.display = "none";
  let data = null, ell = 2, idx = 0, t0 = 0, holdUntil = 0, active = false;
  let bars = [], shown = [], from = [], to = [], tweenStart = 0, tweenMs = 1;
  let label = {}, inset = null;

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
    a = 2 * Math.sqrt(q);
    w = (2 * a) / DISPLAY_BINS;
    const xmax = a * 1.1;
    X = (x) => L + ((x + xmax) / (2 * xmax)) * (R - L);
    let peak = 0;
    for (let i = 0; i <= 400; i++) peak = Math.max(peak, kmDensity(q, -a + (2 * a * i) / 400));
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
      d += `${i ? "L" : "M"}${X(x).toFixed(1)},${Y(kmDensity(q, x)).toFixed(1)}`;
    }
    el("path", { class: "wm-spec-curve", d }, g);
    el("line", { class: "wm-spec-axis", x1: L, x2: R, y1: B, y2: B }, g);
    for (let t = -Math.floor(xmax); t <= Math.floor(xmax); t++) {
      el("line", { class: "wm-spec-tick", x1: X(t), x2: X(t), y1: B, y2: B + 5 }, g);
      text(g, X(t), B + 18, t < 0 ? `−${-t}` : String(t), { class: "wm-spec-label", "text-anchor": "middle" });
    }
    for (const s of [-1, 1]) {
      el("line", { class: "wm-spec-edge", x1: X(s * a), x2: X(s * a), y1: T - 6, y2: B }, g);
      text(g, X(s * a), T - 10, `${s < 0 ? "−" : ""}2√${q}`, { class: "wm-spec-label", "text-anchor": "middle" });
    }
    text(g, (L + R) / 2, H - 14, `eigenvalue λ of T${SUB[q]}`, { class: "wm-spec-label", "text-anchor": "middle" });
    label = {
      p: text(g, L, 34, "", { class: "wm-spec-big" }),
      dim: text(g, L, 56, "", { class: "wm-spec-label" }),
      ks: text(g, L, 76, "", { class: "wm-spec-label" }),
    };
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
    if (reduced) { shown = to.slice(); paintBars(shown); }
    const last = k === data.primes.length - 1;
    onStatus(`T${SUB[ell]} on S₂(Γ₀(${row.p})): ${row.dim} eigenvalues. ` +
      `Kolmogorov distance to the Kesten–McKay law: ${D.toFixed(3)}` +
      (last ? " (largest p shown; the animation starts again)." : "."));
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
      scaffold();
      t0 = performance.now(); holdUntil = 0;
      goTo(reduced ? data.primes.length - 1 : 0, t0, 600);
    },
    hide() { active = false; g.style.display = "none"; },
    setEll(l) {
      ell = Number(l);
      if (active && data) {
        scaffold();
        goTo(idx, performance.now(), 600);
        if (reduced) paintBars(shown);
      }
    },
    step(now, ms) {
      if (!active || !data || reduced) return;
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
