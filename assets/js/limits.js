// The Limit laws page (animations/limits.qmd): six pictures of weak convergence of
// spectra of Hecke operators. Statements and proofs: animations/about.qmd#limits.
//
//   kmst       the Kesten-McKay laws mu_l (of lambda / sqrt l) tend to the Sato-Tate law
//              as l -> infinity; computed here, no data file
//   star       T_2 on S_2(Gamma_1(101)), split by character: the eigenvalues lie on lines
//              through 0, and rotating by chi(2)^(-1/2) puts them all on the real line
//   quadratic  T_2 on S_2(N, (./N)), N = 5 (mod 8): purely imaginary eigenvalues, with
//              Im(lambda)/sqrt 2 tending to mu_2
//   planar     T_2 on all of S_2(Gamma_1(N)), N -> infinity: lambda / sqrt 2 in the plane
//   diagonal   T_l on S_2(Gamma_0(p)) with l and p growing together: Sato-Tate
//   weight     T_l on S_k(SL_2(Z)), k -> infinity: mu_l
//
// Data: assets/data/limits/*.json, made by _scripts/limit_laws.py.
// Colours come from the --wm-* variables, so light and dark mode both work.

const SVG_NS = "http://www.w3.org/2000/svg";
export const W = 600, H = 500;
const el = (name, attrs = {}, parent) => {
  const node = document.createElementNS(SVG_NS, name);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
  if (parent) parent.appendChild(node);
  return node;
};
// Text with subscripts written as T_ℓ or S_{k}: each becomes a lowered, smaller tspan.
function setText(t, str) {
  t.replaceChildren();
  const parts = String(str).split(/_(\{[^}]*\}|.)/u);
  parts.forEach((part, i) => {
    if (!part) return;
    if (i % 2 === 0) { t.appendChild(document.createTextNode(part)); return; }
    const s = el("tspan", { "baseline-shift": "sub", "font-size": "75%" }, t);
    s.textContent = part.startsWith("{") ? part.slice(1, -1) : part;
  });
}
const text = (parent, x, y, str, attrs = {}) => {
  const t = el("text", { x, y, ...attrs }, parent);
  setText(t, str);
  return t;
};
const SUBS = "₀₁₂₃₄₅₆₇₈₉";
export const sub = (n) => String(n).replace(/\d/g, (d) => SUBS[d]);
const fmt = (v, d = 3) => (v < 0 ? "−" : "") + Math.abs(v).toFixed(d);
const ease = (t) => t * t * (3 - 2 * t);

// ---- the laws, all normalised to [-2, 2] ----------------------------------------

// Sato-Tate: the semicircle law, the spectral measure of the infinite path's limit.
export function stDensity(x) {
  const r = 4 - x * x;
  return r <= 0 ? 0 : Math.sqrt(r) / (2 * Math.PI);
}
// Kesten-McKay mu_l for lambda / sqrt l: the semicircle times l(l+1) / ((l+1)^2 - l x^2).
// l may be any real number >= 1 (used for the glide between two primes).
export function kmDensity(l, x) {
  return stDensity(x) * l * (l + 1) / ((l + 1) ** 2 - l * x * x);
}
// sup over x of the density, for the vertical scale
function peak(f) {
  let m = 0;
  for (let i = 0; i <= 800; i++) m = Math.max(m, f(-2 + i / 200));
  return m;
}
// Total variation (half the L1 distance) and Kolmogorov distance between mu_l and the
// Sato-Tate law, by the trapezoid rule on a grid of 4001 points.
export function kmstDistances(l, n = 4000) {
  const h = 4 / n;
  let tv = 0, F = 0, S = 0, ks = 0, prev = null;
  for (let i = 0; i <= n; i++) {
    const x = -2 + i * h, a = kmDensity(l, x), b = stDensity(x);
    if (prev) {
      tv += (Math.abs(a - b) + Math.abs(prev[0] - prev[1])) * h / 2;
      F += (a + prev[0]) * h / 2; S += (b + prev[1]) * h / 2;
      ks = Math.max(ks, Math.abs(F - S));
    }
    prev = [a, b];
  }
  return { tv: tv / 2, ks };
}
// The first-order terms: tv ~ 3 sqrt 3 / (4 pi l), ks ~ 3 sqrt 3 / (8 pi l).
export const TV_CONST = (3 * Math.sqrt(3)) / (4 * Math.PI);
export const KS_CONST = (3 * Math.sqrt(3)) / (8 * Math.PI);

export function primesUpTo(n) {
  const out = [];
  for (let k = 2; k <= n; k++) {
    let prime = true;
    for (let d = 2; d * d <= k; d++) if (k % d === 0) { prime = false; break; }
    if (prime) out.push(k);
  }
  return out;
}

// ---- shared pieces -----------------------------------------------------------------

const L = 56, R = W - 22, T = 128, B = H - 54;
const X = (x) => L + ((x + 2.2) / 4.4) * (R - L);

function curvePath(f, Y, x0 = -2, x1 = 2, n = 300, Xf = X) {
  let d = "";
  for (let i = 0; i <= n; i++) {
    const x = x0 + ((x1 - x0) * i) / n;
    d += `${i ? "L" : "M"}${Xf(x).toFixed(1)},${Y(f(x)).toFixed(1)}`;
  }
  return d;
}
function axisX(g, Y0, xlabel) {
  el("line", { class: "wm-spec-axis", x1: L, x2: R, y1: Y0, y2: Y0 }, g);
  for (let t = -2; t <= 2; t++) {
    el("line", { class: "wm-spec-tick", x1: X(t), x2: X(t), y1: Y0, y2: Y0 + 5 }, g);
    text(g, X(t), Y0 + 18, t < 0 ? `−${-t}` : String(t), { class: "wm-spec-label", "text-anchor": "middle" });
  }
  if (xlabel) text(g, (L + R) / 2, H - 14, xlabel, { class: "wm-spec-label", "text-anchor": "middle" });
}
function headLabels(g) {
  return {
    big: text(g, L, 34, "", { class: "wm-spec-big" }),
    a: text(g, L, 56, "", { class: "wm-spec-label" }),
    b: text(g, L, 76, "", { class: "wm-spec-label" }),
    c: text(g, L, 96, "", { class: "wm-spec-label" }),
  };
}
// The inset in the top right corner: a distance against the parameter, log scale
// across, so far.
function makeInset(g, keys, vals, { title = "distance", xname = "p", logY = false, guide = null } = {}) {
  const ix0 = W - 206, ix1 = W - 76, iy0 = 22, iy1 = 76;   // clear of the eye icon
  const lk = keys.map(Math.log), kmin = Math.min(...lk), kmax = Math.max(...lk);
  const ys = vals.map((v) => (logY ? Math.log(v) : v));
  let vmax = Math.max(...ys), vmin = logY ? Math.min(...ys) : 0;
  if (!logY) vmax *= 1.1;
  const IX = (k) => ix0 + ((Math.log(k) - kmin) / (kmax - kmin || 1)) * (ix1 - ix0);
  const IY = (v) => iy1 - (((logY ? Math.log(v) : v) - vmin) / (vmax - vmin || 1)) * (iy1 - iy0);
  el("rect", { class: "wm-spec-inset", x: ix0 - 34, y: iy0 - 12, width: ix1 - ix0 + 44, height: iy1 - iy0 + 34, rx: 8 }, g);
  text(g, ix0 - 26, iy0 - 1, title, { class: "wm-spec-small" });
  text(g, ix1, iy1 + 16, `${xname} →`, { class: "wm-spec-small", "text-anchor": "end" });
  text(g, ix0 - 26, iy1 + 16, String(keys[0]), { class: "wm-spec-small" });
  if (guide) {
    let d = "";
    keys.forEach((k, i) => { d += `${i ? "L" : "M"}${IX(k).toFixed(1)},${IY(guide(k)).toFixed(1)}`; });
    el("path", { class: "wm-lim-guide", d }, g);
  }
  const trace = el("path", { class: "wm-spec-trace" }, g);
  const dot = el("circle", { class: "wm-spec-dot wm-spec-now", r: 3.6 }, g);
  return {
    upTo(i) {
      let d = "";
      for (let j = 0; j <= i; j++) d += `${j ? "L" : "M"}${IX(keys[j]).toFixed(1)},${IY(vals[j]).toFixed(1)}`;
      trace.setAttribute("d", d);
      dot.setAttribute("cx", IX(keys[i]).toFixed(1));
      dot.setAttribute("cy", IY(vals[i]).toFixed(1));
    },
  };
}

// A sequence of histograms on [-2, 2] gliding one into the next, against a limit law:
// used for the quadratic family, the diagonal limit and the weight aspect.
function histogramPlayer(g, { rows, law, extra, ymax, xlabel, inset, label, status, merge = 2 }) {
  const nb = rows[0].hist.length / merge, w = 4 / nb;   // stored bins, merged in groups
  const Y = (y) => B - (Math.min(y, ymax) / ymax) * (B - T);
  const bars = [];
  for (let k = 0; k < nb; k++) {
    const x0 = X(-2 + k * w), x1 = X(-2 + (k + 1) * w);
    bars.push(el("rect", { class: "wm-spec-bar", x: (x0 + 1).toFixed(1), y: B,
      width: Math.max(0.5, x1 - x0 - 2).toFixed(1), height: 0, rx: 2 }, g));
  }
  const extraPath = extra ? el("path", { class: "wm-lim-second" }, g) : null;
  el("path", { class: "wm-spec-curve", d: curvePath(law, Y) }, g);
  axisX(g, B, xlabel);
  for (const s of [-1, 1]) el("line", { class: "wm-spec-edge", x1: X(2 * s), x2: X(2 * s), y1: T - 6, y2: B }, g);
  const lab = headLabels(g);
  const ins = makeInset(g, rows.map(inset.key), rows.map(inset.val), inset);
  let shown = new Array(nb).fill(0), from = shown, to = shown, t0 = 0, tw = 1;
  const heights = (row) => {
    const n = row.hist.reduce((a, b) => a + b, 0), out = [];
    for (let b = 0; b < nb; b++) {
      let c = 0;
      for (let i = 0; i < merge; i++) c += row.hist[b * merge + i];
      out.push(c / (n * w));
    }
    return out;
  };
  const paint = (hs) => hs.forEach((h, k) => {
    const y = Y(h);
    bars[k].setAttribute("y", y.toFixed(1));
    bars[k].setAttribute("height", (B - y).toFixed(1));
  });
  return {
    count: rows.length,
    go(i, now, ms, still) {
      const row = rows[i];
      label(lab, row);
      ins.upTo(i);
      if (extraPath) extraPath.setAttribute("d", curvePath((x) => extra(row, x), Y));
      from = shown.slice(); to = heights(row); t0 = now; tw = Math.max(1, ms * 0.85);
      if (still) { shown = to.slice(); paint(shown); }
      status(row, i === rows.length - 1);
    },
    glide(now) {
      const t = Math.min(1, (now - t0) / tw), e = ease(t);
      shown = from.map((f, k) => f + (to[k] - f) * e);
      paint(shown);
    },
    reset() { shown = new Array(nb).fill(0); },
  };
}

// ---- kmst: Kesten-McKay to Sato-Tate -------------------------------------------------

// Draws mu_l for l running through the primes up to 1009 against the semicircle, into
// the group g. Used on the Limit laws page and by the "ℓ → ∞" entry of the spectral law
// (assets/js/spectral.js).
export function createKmSt(g, { onStatus = () => {} } = {}) {
  const ELLS = primesUpTo(1009);
  const dist = ELLS.map((l) => kmstDistances(l));
  const ymax = peak((x) => kmDensity(2, x)) * 1.12;
  const Y = (y) => B - (Math.min(y, ymax) / ymax) * (B - T);
  g.replaceChildren();
  const area = el("path", { class: "wm-lim-area" }, g);
  const line = el("path", { class: "wm-lim-line" }, g);
  el("path", { class: "wm-spec-curve", d: curvePath(stDensity, Y) }, g);
  axisX(g, B, "x = λ/√ℓ, λ in the spectrum of the (ℓ + 1)-regular tree");
  const lab = headLabels(g);
  const ins = makeInset(g, ELLS, dist.map((d) => d.tv), {
    title: "total variation (log)", xname: "ℓ", logY: true, guide: (l) => TV_CONST / l,
  });
  let idx = 0, q = 2, qFrom = 2, qTo = 2, t0 = 0, tw = 1;
  function draw(l) {
    const f = (x) => kmDensity(l, x);
    const top = curvePath(f, Y);
    line.setAttribute("d", top);
    area.setAttribute("d", `${top}L${X(2).toFixed(1)},${B}L${X(-2).toFixed(1)},${B}Z`);
  }
  function go(i, now, ms, still) {
    idx = i;
    const l = ELLS[i], d = dist[i];
    lab.big.textContent = `ℓ = ${l}`;
    lab.a.textContent = `total variation distance ${d.tv.toFixed(4)}`;
    lab.b.textContent = `Kolmogorov distance ${d.ks.toFixed(4)}`;
    lab.c.textContent = `(3√3/(4πℓ) = ${(TV_CONST / l).toFixed(4)}, 3√3/(8πℓ) = ${(KS_CONST / l).toFixed(4)})`;
    ins.upTo(i);
    qFrom = q; qTo = l; t0 = now; tw = Math.max(1, ms * 0.85);
    if (still) { q = l; draw(q); }
    onStatus(`μ${sub(l)} against the Sato–Tate law: total variation distance ${d.tv.toFixed(4)}, ` +
      `Kolmogorov distance ${d.ks.toFixed(4)}` + (i === ELLS.length - 1 ? " (ℓ = 1009, the last; it starts again)." : "."));
  }
  return {
    count: ELLS.length, ells: ELLS, go,
    glide(now) {
      const t = Math.min(1, (now - t0) / tw), e = ease(t);
      q = Math.exp(Math.log(qFrom) + (Math.log(qTo) - Math.log(qFrom)) * e);
      draw(q);
    },
    reset() { q = 2; },
    index: () => idx,
  };
}

// ---- star: S_2(Gamma_1(101)) by character ---------------------------------------------

function createStar(g, data, { onStatus }) {
  g.replaceChildren();
  const N = data.N;
  // Characters sorted by the angle of their line, arg(chi(2))/2 in [0, pi).
  const chars = data.chars.map((c) => {
    const arg = Math.atan2(c.chi[1], c.chi[0]);           // in (-pi, pi]
    const half = ((arg / 2) % Math.PI + Math.PI) % Math.PI; // the line's angle in [0, pi)
    return { ...c, half, w: arg / 2 };                     // w: a square root of chi(2) is e^{i w}
  }).sort((a, b) => a.half - b.half);
  const cx = W / 2, cy = 300, rad = 140;                   // |z| = 2 at radius rad
  const rowTop = 140, rowBot = 360, xs0 = 70, xs1 = W - 70;
  const RX = (x) => xs0 + ((x + 2) / 4) * (xs1 - xs0);
  el("circle", { class: "wm-lim-disk", cx, cy, r: rad }, g);
  const lines = el("g", { class: "wm-lim-lines" }, g);
  for (const c of chars) {
    const dx = Math.cos(c.half) * rad, dy = -Math.sin(c.half) * rad;
    el("line", { x1: (cx - dx).toFixed(1), y1: (cy - dy).toFixed(1), x2: (cx + dx).toFixed(1), y2: (cy + dy).toFixed(1) }, lines);
  }
  // The histogram of the rotated values, shown with the rows.
  const hb = 20, hw = 4 / hb, hTop = 380, hBot = 452;
  const all = [];
  const pts = [];
  chars.forEach((c, row) => {
    const yRow = rowTop + ((row + 0.5) / chars.length) * (rowBot - rowTop);
    for (const [re, im] of c.eig) {
      // z = lambda / sqrt 2; rotated: z e^{-i w}, real up to rounding
      const zr = re / Math.SQRT2, zi = im / Math.SQRT2;
      const x = zr * Math.cos(c.w) + zi * Math.sin(c.w);
      all.push(x);
      const dot = el("circle", { class: c.a === 1 ? "wm-lim-pt wm-lim-pt-triv" : "wm-lim-pt", r: 2.4 }, g);
      pts.push({ dot, sx: cx + (zr / 2) * rad, sy: cy - (zi / 2) * rad, rx: RX(x), ry: yRow });
    }
  });
  const hist = new Array(hb).fill(0);
  for (const x of all) hist[Math.min(hb - 1, Math.max(0, Math.floor((x + 2) / hw)))]++;
  const hmax = peak((x) => kmDensity(2, x)) * 1.15;
  const HY = (y) => hBot - (Math.min(y, hmax) / hmax) * (hBot - hTop);
  const hg = el("g", { class: "wm-lim-rowhist" }, g);
  hist.forEach((c, k) => {
    const y = HY(c / (all.length * hw));
    el("rect", { class: "wm-spec-bar", x: (RX(-2 + k * hw) + 1).toFixed(1), y: y.toFixed(1),
      width: (RX(-2 + (k + 1) * hw) - RX(-2 + k * hw) - 2).toFixed(1), height: (hBot - y).toFixed(1), rx: 2 }, hg);
  });
  el("path", { class: "wm-spec-curve", d: curvePath((x) => kmDensity(2, x), HY, -2, 2, 200, RX) }, hg);
  el("line", { class: "wm-spec-axis", x1: RX(-2), x2: RX(2), y1: hBot, y2: hBot }, hg);
  for (let t = -2; t <= 2; t++) text(hg, RX(t), hBot + 15, t < 0 ? `−${-t}` : String(t), { class: "wm-spec-small", "text-anchor": "middle" });
  const lab = headLabels(g);
  lab.big.textContent = `S${sub(2)}(Γ${sub(1)}(${N})) = ⊕ S${sub(2)}(${N}, χ)`;
  lab.a.textContent = `${all.length} eigenvalues λ of T₂, drawn as λ/√2; ${chars.length} even characters χ`;
  // t = 0: the plane; t = 1: rotated by chi(2)^(-1/2), one row per character.
  let t = 0, from = 0, to = 0, t0 = 0, tw = 1;
  function paint(s) {
    const e = ease(s);
    for (const p of pts) {
      p.dot.setAttribute("cx", (p.sx + (p.rx - p.sx) * e).toFixed(1));
      p.dot.setAttribute("cy", (p.sy + (p.ry - p.sy) * e).toFixed(1));
    }
    lines.style.opacity = String(1 - e);
    hg.style.opacity = String(e);
    g.querySelector(".wm-lim-disk").style.opacity = String(1 - e);
  }
  function label(rot) {
    lab.b.textContent = rot
      ? "rotated by χ(2)^(−1/2): every value is real; one row per character"
      : `on ${data.lines} lines through 0, one for each character (red: χ trivial)`;
    lab.c.textContent = rot ? "below: all of them together, against μ₂" : "";
  }
  return {
    count: 2,
    go(i, now, ms, still) {
      from = t; to = i; t0 = now; tw = Math.max(1, ms * 1.6);
      label(i === 1);
      if (still) { t = to; paint(t); }
      onStatus(i === 1
        ? `Rotated: each eigenvalue λ of T₂ on S₂(${N}, χ) times χ(2)^(−1/2)/√2 is real, and together they follow μ₂.`
        : `The eigenvalues of T₂ on S₂(Γ₁(${N})) in the plane: those on S₂(${N}, χ) lie on the line through 0 at angle arg χ(2)/2.`);
    },
    glide(now) { t = from + (to - from) * Math.min(1, (now - t0) / tw); paint(t); },
    reset() {},
    hold: 3,
  };
}

// ---- planar: S_2(Gamma_1(N)) as N grows -----------------------------------------------

function createPlanar(g, data, { onStatus }) {
  g.replaceChildren();
  const cx = 170, cy = 300, rad = 135;
  el("circle", { class: "wm-lim-disk", cx, cy, r: rad }, g);
  el("line", { class: "wm-spec-tick", x1: cx - rad, x2: cx + rad, y1: cy, y2: cy }, g);
  el("line", { class: "wm-spec-tick", x1: cx, x2: cx, y1: cy - rad, y2: cy + rad }, g);
  text(g, cx + rad + 4, cy + 4, "2", { class: "wm-spec-small" });
  text(g, cx, cy + rad + 16, "λ/√2 in ℂ", { class: "wm-spec-label", "text-anchor": "middle" });
  const dots = el("g", { class: "wm-lim-dots" }, g);
  // right: the radial law, |z| against the density 2 f(r) on [0, 2], f that of mu_2
  const px0 = 345, px1 = W - 24, pTop = 190, pBot = 430, nb = 20, bw = 2 / nb;
  const PX = (r) => px0 + (r / 2) * (px1 - px0);
  const radial = (r) => 2 * kmDensity(2, r);
  const ymax = peak((x) => (x >= 0 ? radial(x) : 0)) * 1.2;
  const PY = (y) => pBot - (Math.min(y, ymax) / ymax) * (pBot - pTop);
  const bars = [];
  for (let k = 0; k < nb; k++) {
    bars.push(el("rect", { class: "wm-spec-bar", x: (PX(k * bw) + 1).toFixed(1), y: pBot,
      width: (PX((k + 1) * bw) - PX(k * bw) - 2).toFixed(1), height: 0, rx: 2 }, g));
  }
  el("path", { class: "wm-spec-curve", d: curvePath(radial, PY, 0, 2, 200, PX) }, g);
  el("line", { class: "wm-spec-axis", x1: px0, x2: px1, y1: pBot, y2: pBot }, g);
  for (let r = 0; r <= 2; r += 0.5) {
    text(g, PX(r), pBot + 15, String(r), { class: "wm-spec-small", "text-anchor": "middle" });
  }
  text(g, (px0 + px1) / 2, pBot + 34, "|λ|/√2", { class: "wm-spec-label", "text-anchor": "middle" });
  text(g, px0, pTop - 12, "radius, against |X| for X ~ μ₂", { class: "wm-spec-small" });
  const lab = headLabels(g);
  const rows = data.rows;
  const ins = makeInset(g, rows.map((r) => r.N), rows.map((r) => r.ksr), { title: "radial distance", xname: "N" });
  let shown = new Array(nb).fill(0), from = shown, to = shown, t0 = 0, tw = 1;
  const paintBars = (hs) => hs.forEach((h, k) => {
    const y = PY(h);
    bars[k].setAttribute("y", y.toFixed(1));
    bars[k].setAttribute("height", (pBot - y).toFixed(1));
  });
  return {
    count: rows.length,
    go(i, now, ms, still) {
      const row = rows[i], P = row.pts, n = P.length / 2;
      const frag = document.createDocumentFragment();
      const counts = new Array(nb).fill(0);
      const rr = n > 1500 ? 1.3 : n > 600 ? 1.7 : 2.2;
      for (let j = 0; j < n; j++) {
        const zr = P[2 * j], zi = P[2 * j + 1];
        frag.appendChild(el("circle", { class: Math.abs(zi) < 1e-9 ? "wm-lim-pt wm-lim-pt-triv" : "wm-lim-pt",
          cx: (cx + (zr / 2) * rad).toFixed(1), cy: (cy - (zi / 2) * rad).toFixed(1), r: rr }));
        counts[Math.min(nb - 1, Math.floor(Math.hypot(zr, zi) / bw))]++;
      }
      dots.replaceChildren(frag);
      lab.big.textContent = `N = ${row.N}`;
      lab.a.textContent = `${n} eigenvalues of T₂ on S₂(Γ₁(${row.N}))`;
      lab.b.textContent = `on ${row.lines} lines through 0;  |mean z²| = ${row.mom["20"].toFixed(3)}`;
      lab.c.textContent = `mean |z|² = ${row.mom["11"].toFixed(3)}  (limit 3/2)`;
      ins.upTo(i);
      from = shown.slice(); to = counts.map((c) => c / (n * bw)); t0 = now; tw = Math.max(1, ms * 0.85);
      if (still) { shown = to.slice(); paintBars(shown); }
      onStatus(`T₂ on S₂(Γ₁(${row.N})): ${n} eigenvalues z = λ/√2 on ${row.lines} lines through 0. ` +
        `Radial Kolmogorov distance ${row.ksr.toFixed(3)}` + (i === rows.length - 1 ? " (the largest N shown; it starts again)." : "."));
    },
    glide(now) {
      const e = ease(Math.min(1, (now - t0) / tw));
      shown = from.map((f, k) => f + (to[k] - f) * e);
      paintBars(shown);
    },
    reset() { shown = new Array(nb).fill(0); },
  };
}

// ---- the three histogram sequences --------------------------------------------------

function createQuadratic(g, data, { onStatus }) {
  g.replaceChildren();
  return histogramPlayer(g, {
    rows: data.rows,
    law: (x) => kmDensity(2, x),
    ymax: peak((x) => kmDensity(2, x)) * 1.6,
    xlabel: "y = λ/(i√2), λ an eigenvalue of T₂ on S₂(N, (·/N))",
    inset: { key: (r) => r.N, val: (r) => r.ks, title: "distance", xname: "N" },
    label(lab, r) {
      lab.big.textContent = `N = ${r.N}`;
      lab.a.textContent = `${r.dim} eigenvalues = dim S₂(N, (·/N))`;
      lab.b.textContent = "(2/N) = −1, so all of them lie in iℝ";
      lab.c.textContent = `Kolmogorov distance to μ₂ ${r.ks.toFixed(3)}`;
    },
    status: (r, last) => onStatus(`T₂ on S₂(${r.N}, (·/N)): ${r.dim} purely imaginary eigenvalues. ` +
      `Kolmogorov distance of λ/(i√2) to μ₂: ${r.ks.toFixed(3)}` + (last ? " (the largest N shown; it starts again)." : ".")),
  });
}

function createDiagonal(g, data, { onStatus }) {
  g.replaceChildren();
  return histogramPlayer(g, {
    rows: data.rows,
    law: stDensity,
    extra: (r, x) => kmDensity(r.ell, x),
    ymax: (1 / Math.PI) * 1.75,
    xlabel: "x = λ/√ℓ, λ an eigenvalue of T_ℓ on S₂(Γ₀(p)), ℓ the largest prime ≤ (log p)²",
    inset: { key: (r) => r.p, val: (r) => r.ks, title: "distance to Sato–Tate", xname: "p" },
    label(lab, r) {
      lab.big.textContent = `p = ${r.p},  ℓ = ${r.ell}`;
      lab.a.textContent = `${r.dim} eigenvalues of T${sub(r.ell)} on S₂(Γ₀(${r.p}))`;
      lab.b.textContent = `Kolmogorov distance to Sato–Tate ${r.ks.toFixed(3)}`;
      lab.c.textContent = `to μ${sub(r.ell)} (dashed): ${r.ksl.toFixed(3)}`;
    },
    status: (r, last) => onStatus(`T${sub(r.ell)} on S₂(Γ₀(${r.p})): ${r.dim} eigenvalues. Kolmogorov distance to the ` +
      `Sato–Tate law ${r.ks.toFixed(3)}` + (last ? " (the largest p shown; it starts again)." : ".")),
  });
}

function createWeight(g, data, ell, { onStatus }) {
  g.replaceChildren();
  const rows = data.ells[String(ell)];
  return histogramPlayer(g, {
    rows,
    law: (x) => kmDensity(ell, x),
    ymax: peak((x) => kmDensity(ell, x)) * 1.6,
    xlabel: `x = λ/${ell}^((k−1)/2), λ an eigenvalue of T${sub(ell)} on S_k(SL₂(ℤ))`,
    inset: { key: (r) => r.k, val: (r) => r.ks, title: "distance", xname: "k" },
    label(lab, r) {
      lab.big.textContent = `k = ${r.k}`;
      setText(lab.a, `${r.dim} eigenvalue${r.dim === 1 ? "" : "s"} = dim S${sub(r.k)}(SL₂(ℤ))`);
      lab.b.textContent = `Kolmogorov distance to μ${sub(ell)} ${r.ks.toFixed(3)}`;
      lab.c.textContent = "";
    },
    status: (r, last) => onStatus(`T${sub(ell)} on S${sub(r.k)}(SL₂(ℤ)): ${r.dim} eigenvalues. Kolmogorov distance to ` +
      `μ${sub(ell)} ${r.ks.toFixed(3)}` + (last ? " (the largest weight shown; it starts again)." : ".")),
  });
}

// ---- the page widget -----------------------------------------------------------------

const MS = { calm: 1100, fast: 520, frantic: 220 };
const FILES = { star: "gamma1", quadratic: "quadratic", planar: "planar", diagonal: "diagonal", weight: "weight" };

function init(root, base) {
  const svg = root.querySelector("svg");
  const stats = root.querySelector(".wm-graph-stats");
  const modeSel = root.querySelector(".wm-lim-mode");
  const ellSel = root.querySelector(".wm-lim-ell select");
  const speedSel = root.querySelector(".wm-graph-speed");
  const playBtn = root.querySelector(".wm-graph-play");
  const slider = root.querySelector(".wm-lim-slider input");
  const eye = root.querySelector(".wm-graph-eye");
  const infoBlocks = [...root.querySelectorAll(".wm-anim-info [data-anim]")];
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  const g = el("g", {}, svg);
  const reduced = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  let playing = !reduced, visible = true, raf = 0;
  let player = null, idx = 0, last = 0, holdUntil = 0, buildId = 0;
  const cache = new Map();
  const load = (name) => {
    if (!cache.has(name)) {
      cache.set(name, fetch(`${base}/${name}.json`).then((r) => {
        if (!r.ok) throw new Error(r.status);
        return r.json();
      }).catch((e) => { cache.delete(name); throw e; }));
    }
    return cache.get(name);
  };
  const host = root.closest("[data-limit]");
  if (host && [...modeSel.options].some((o) => o.value === host.dataset.limit)) modeSel.value = host.dataset.limit;
  const onStatus = (msg) => { stats.textContent = msg; };
  const ms = () => MS[speedSel.value];

  async function build() {
    const id = ++buildId, mode = modeSel.value;
    root.dataset.mode = mode;
    for (const b of infoBlocks) b.hidden = b.dataset.anim !== mode;
    if (eye) eye.href = `${eye.dataset.base}#${mode === "kmst" ? "kmst" : mode}`;
    let made;
    try {
      if (mode === "kmst") made = createKmSt(g, { onStatus });
      else {
        const data = await load(FILES[mode]);
        if (id !== buildId) return;
        if (mode === "star") made = createStar(g, data, { onStatus });
        else if (mode === "planar") made = createPlanar(g, data, { onStatus });
        else if (mode === "quadratic") made = createQuadratic(g, data, { onStatus });
        else if (mode === "diagonal") made = createDiagonal(g, data, { onStatus });
        else made = createWeight(g, data, Number(ellSel.value), { onStatus });
      }
    } catch {
      onStatus("Could not load the eigenvalues. Please try again.");
      return;
    }
    player = made;
    if (slider) { slider.max = String(player.count - 1); slider.value = "0"; }
    // Paused: the last frame, where the picture is closest to the limit.
    goTo(playing ? 0 : player.count - 1);
  }
  function goTo(i) {
    idx = i;
    const now = performance.now();
    last = now; holdUntil = 0;
    player.go(i, now, ms(), !playing);
    if (slider) slider.value = String(i);
  }
  function frame(now) {
    raf = 0;
    if (!playing || !visible || document.hidden || !player) return;
    player.glide(now);
    if (holdUntil) {
      if (now >= holdUntil) { holdUntil = 0; player.reset(); goTo(0); }
    } else if (now - last >= ms() * (player.hold || 1)) {
      if (idx < player.count - 1) {
        goTo(idx + 1);
        if (idx === player.count - 1 && player.count > 2) holdUntil = now + 8 * ms();
      } else if (player.count === 2) goTo(0);
    }
    raf = requestAnimationFrame(frame);
  }
  function startLoop() { if (playing && !raf) raf = requestAnimationFrame(frame); }
  function setPlaying(on) {
    playing = on;
    playBtn.textContent = on ? "Pause" : "Play";
    playBtn.setAttribute("aria-label", on ? "Pause the animation" : "Play the animation");
    if (!on && raf) { cancelAnimationFrame(raf); raf = 0; }
    if (on) { last = performance.now(); startLoop(); }
  }
  modeSel.addEventListener("change", () => { build().then(startLoop); });
  ellSel.addEventListener("change", () => { if (modeSel.value === "weight") build().then(startLoop); });
  playBtn.addEventListener("click", () => setPlaying(!playing));
  if (slider) {
    slider.addEventListener("input", () => {
      if (!player) return;
      setPlaying(false);
      goTo(Number(slider.value));
      player.glide(performance.now() + 1e6);
    });
  }
  if ("IntersectionObserver" in window) {
    new IntersectionObserver(([e]) => { visible = e.isIntersecting; if (visible) startLoop(); }, { threshold: 0.05 }).observe(root);
  }
  document.addEventListener("visibilitychange", () => { if (!document.hidden) startLoop(); });
  setPlaying(playing);
  build().then(() => { startLoop(); root.classList.add("wm-graph-ready"); });
}

const root = typeof document !== "undefined" ? document.getElementById("wm-limits") : null;
if (root) init(root, new URL("../data/limits", import.meta.url).href);
