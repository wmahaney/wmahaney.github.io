// Interactive supersingular isogeny graphs for the home page.
// Data: assets/data/isogeny/ (index.json plus one p<prime>-<graph>.json per prime and
// graph, loaded when chosen), produced by _scripts/isogeny_graphs.py.
// Colours come from the --wm-* variables, so light and dark mode both work.
import {
  forceSimulation, forceLink, forceManyBody, forceCollide, forceX, forceY,
} from "./vendor/d3-force.min.js";
import { createSpectral } from "./spectral.js";
import {
  adjacency, linkCounts, walkFor, stationary, relativeStationary, sampleStep,
} from "./walk.js";

const W = 600, H = 500, PAD = 40;
const SVG_NS = "http://www.w3.org/2000/svg";
const KINDS = {
  "2": {
    defaultP: 431, ells: [2], bound: "2√2", boundVal: 2 * Math.SQRT2,
    caption: (p, n) => `The supersingular 2-isogeny graph in characteristic ${p}: ` +
      `${n} curves, each with three outgoing 2-isogenies.`,
  },
  "23": {
    defaultP: 211, ells: [2, 3], bound: "2√2 + 2√3",
    boundVal: 2 * Math.SQRT2 + 2 * Math.sqrt(3),
    caption: (p, n) => `The (2,3)-isogeny graph in characteristic ${p}: the same ${n} curves, ` +
      `joined by 2-isogenies and 3-isogenies.`,
  },
  "235": {
    defaultP: 157, ells: [2, 3, 5], bound: "2√2 + 2√3 + 2√5",
    boundVal: 2 * Math.SQRT2 + 2 * Math.sqrt(3) + 2 * Math.sqrt(5),
    caption: (p, n) => `The (2,3,5)-isogeny graph in characteristic ${p}: the same ${n} curves, ` +
      `joined by 2-, 3- and 5-isogenies.`,
  },
  "2b5": {
    defaultP: 101, level: true, ells: [2], bound: "2√2", boundVal: 2 * Math.SQRT2,
    caption: (p, n) => `The 2-isogeny graph with Borel level-5 structure in characteristic ${p}: ` +
      `${n} pairs (E, C), where C is a subgroup of E of order 5.`,
  },
  "23b5": {
    defaultP: 101, level: true, ells: [2, 3], bound: "2√2 + 2√3",
    boundVal: 2 * Math.SQRT2 + 2 * Math.sqrt(3),
    caption: (p, n) => `The (2,3)-isogeny graph with Borel level-5 structure in characteristic ${p}: ` +
      `${n} pairs (E, C), where C is a subgroup of E of order 5.`,
  },
  "2g15": {
    defaultP: 101, level: true, ells: [2], bound: "2√2", boundVal: 2 * Math.SQRT2,
    caption: (p, n) => `The 2-isogeny graph with Γ₁(5) structure in characteristic ${p}: ` +
      `${n} pairs (E, ±P), where P is a point of order 5. It is directed: the dual of an edge ` +
      "out of (E, ±P) leads back to (E, ±2P), so the walk follows the edges' directions.",
    tip: (l) => `(j(E), j(E/⟨P⟩)) = (${l[0]}, ${l[1]})`,
  },
  "3g2": {
    defaultP: 101, level: true, ells: [3], bound: "2√3", boundVal: 2 * Math.sqrt(3),
    caption: (p, n) => `The 3-isogeny graph with full level-2 structure in characteristic ${p}: ` +
      `${n} pairs (E, (P, Q)), (P, Q) an ordered basis of E[2]; these are the supersingular ` +
      "Legendre invariants λ, and there are (p − 1)/2 of them.",
    tip: (l) => `j(E) = ${l[0]}, λ = ${l[1]}`,
  },
  "3g4": {
    defaultP: 43, level: true, ells: [3], bound: null,
    caption: (p, n) => `The 3-isogeny graph with full level-4 structure in characteristic ${p}: ` +
      `${n} = 4(p − 1) pairs (E, (P, Q)), (P, Q) an ordered basis of E[4]. Colours: the two ` +
      "values ±i of the Weil pairing e₄(P, Q).",
    spectrum: () => "Connected and bipartite: a 3-isogeny sends e₄(P, Q) to its cube, so every edge " +
      "joins the two Weil-pairing classes, and −4 is an eigenvalue.",
    tip: (l) => `j(E) = ${l[0]}, Weil-pairing class ${l[1]}`,
  },
  "5g4": {
    defaultP: 43, level: true, ells: [5], bound: null,
    caption: (p, n) => `The 5-isogeny graph with full level-4 structure in characteristic ${p}: ` +
      `${n} = 4(p − 1) pairs (E, (P, Q)), (P, Q) an ordered basis of E[4]. Colours: the two ` +
      "values ±i of the Weil pairing e₄(P, Q).",
    spectrum: () => "Disconnected: a 5-isogeny sends e₄(P, Q) to its fifth power, which is itself, " +
      "so the two Weil-pairing classes are the two components, and the eigenvalue 6 appears twice.",
    tip: (l) => `j(E) = ${l[0]}, Weil-pairing class ${l[1]}`,
  },
};
// Why the F_p-subgraph (the vertices fixed by Frobenius) is empty, when the reason is a
// theorem. Write pi for the Frobenius of a model of E over F_p (or of a twist): it has
// trace 0, so pi^2 = -p. Proofs on the Read more page.
const EMPTY_WHY = {
  "2b5": (p) => (p % 5 === 2 || p % 5 === 3) ?
    `A fixed (E, C) needs Frobenius to preserve C, so to act on it by some λ with λ² ≡ −p (mod 5); ` +
    `there is no such λ when p ≡ ${p % 5} (mod 5).` : null,
  "2g15": (p) => (p % 5 !== 4) ?
    "A fixed (E, ±P) needs Frobenius to send P to ±P, so (±1)² ≡ −p (mod 5), i.e. p ≡ 4 (mod 5); " +
    `here p ≡ ${p % 5} (mod 5).` : null,
  "3g2": (p) => (p % 4 === 1) ?
    "A fixed (E, (P, Q)) needs all of E[2] to be fixed by Frobenius π, so (π − 1)/2 to be an " +
    "endomorphism; since π² = −p that needs p ≡ 3 (mod 4), and here p ≡ 1 (mod 4)." : null,
};
EMPTY_WHY["23b5"] = EMPTY_WHY["2b5"];
EMPTY_WHY["3g4"] = EMPTY_WHY["5g4"] = () =>
  "Frobenius sends e₄(P, Q) to e₄(P, Q)^p = e₄(P, Q)⁻¹ (here p ≡ 3 mod 4), so it swaps the two " +
  "Weil-pairing classes and fixes no vertex: for Γ(4) this subgraph is always empty.";

// The view menu's names for the three views.
const VIEW_NAME = { full: "full graph", spine: "F_p-subgraph", folded: "quotient by Frobenius",
                    over: "part over the F_p-subgraph of G_p" };

// The three ways to look at each graph (the data come from _scripts/isogeny_graphs.py).
//  full    every vertex and edge.
//  spine   the F_p-subgraph (also called the spine): the subgraph induced on the vertices
//          fixed by Frobenius, i.e. those vertices and every edge between them. Frobenius
//          acts on the level structure too, so with level structure this is smaller than
//          the set of vertices with j(E) in F_p (and it can be empty).
//  folded  the quotient by Frobenius (E, C) -> (E^(p), C^(p)): vertices are orbits of
//          vertices, edges are orbits of edges. The edge orbits need Frobenius on edges,
//          so they are computed by the script ("fold"), not here.
function derive(g, view) {
  const n = g.nodes.length, sigma = g.sigma, fixed = g.fp;
  if (view === "spine") {
    const keep = [], at = new Map();
    for (let i = 0; i < n; i++) if (fixed[i]) { at.set(i, keep.length); keep.push(i); }
    return {
      nodes: keep.map((i) => ({ label: g.nodes[i], members: [i], fp: true })),
      links: g.links.filter(([s, t]) => at.has(s) && at.has(t))
        .map(([s, t, ...rest]) => [at.get(s), at.get(t), ...rest]),
    };
  }
  if (view === "over") {
    // Level structure only: the vertices (E, .) whose curve E has j(E) in F_p, i.e. the
    // preimage of the F_p-subgraph of the graph without level structure, and the edges
    // between them. The vertices fixed by Frobenius (the F_p-subgraph proper) are gold,
    // and Frobenius itself is drawn as dashed links v - sigma(v).
    const jfp = g.jfp, keep = [], at = new Map();
    for (let i = 0; i < n; i++) if (jfp[i]) { at.set(i, keep.length); keep.push(i); }
    const frob = [];
    for (const i of keep) {
      const k = sigma[i];
      if (k !== i && at.has(k) && !(sigma[k] === i && k < i)) frob.push([at.get(i), at.get(k)]);
    }
    return {
      nodes: keep.map((i) => ({ label: g.nodes[i], members: [i], fp: fixed[i] })),
      links: g.links.filter(([s, t]) => at.has(s) && at.has(t))
        .map(([s, t, ...rest]) => [at.get(s), at.get(t), ...rest]),
      frob,
    };
  }
  if (view === "folded") {
    // Frobenius has order 2 (order 4 with Gamma_1(5) structure when p = +-2 mod 5):
    // its orbits have 1, 2 or 4 elements, each named by its smallest member.
    const at = new Map(), nodes = [];
    for (const cyc of cycles(sigma)) {
      for (const m of cyc) at.set(m, nodes.length);
      nodes.push({ label: g.nodes[cyc[0]], members: cyc, fp: cyc.length === 1, pair: cyc.length > 1 });
    }
    return { nodes, links: g.fold.map(([s, t, ...rest]) => [at.get(s), at.get(t), ...rest]) };
  }
  return { nodes: g.nodes.map((label, i) => ({ label, members: [i], fp: fixed[i] })), links: g.links };
}

// The cycles of a permutation, each starting at its smallest element.
function cycles(perm) {
  const seen = new Set(), out = [];
  for (let v = 0; v < perm.length; v++) {
    if (seen.has(v)) continue;
    const cyc = [];
    for (let w = v; !seen.has(w); w = perm[w]) { seen.add(w); cyc.push(w); }
    out.push(cyc);
  }
  return out;
}

// The j-invariant(s) to write next to a vertex. Labels are strings "a", "a + bi", "bi"
// (or pairs [j(E), j(E/C)] with level structure, where j(E) is used). A conjugate pair
// {a + bi, a - bi} is written a ± bi.
const jOf = (label) => (Array.isArray(label) ? label[0] : label);
function shortJ(nodes, nd) {
  const js = [...new Set(nd.members.map((m) => jOf(nodes[m])))];
  if (js.length === 1) return js[0];
  if (js.length > 2) return js.join(" / ");
  const [j1, j2] = js;
  const parse = (j) => {
    const m = /^(?:(\d+) \+ )?(\d*)i$/.exec(j);
    return m ? [m[1] || "0", Number(m[2] || 1)] : null;
  };
  const u = parse(j1), v = parse(j2);
  if (!u || !v || u[0] !== v[0]) return `${j1} / ${j2}`;
  const b = Math.min(u[1], v[1]), bi = b === 1 ? "i" : `${b}i`;
  return u[0] === "0" ? `±${bi}` : `${u[0]} ± ${bi}`;
}
// Labels are drawn only for graphs up to this many vertices; beyond that, hover instead.
const LABEL_MAX = 60;

function components(n, links) {
  const parent = [...Array(n).keys()];
  const find = (x) => { while (parent[x] !== x) x = parent[x] = parent[parent[x]]; return x; };
  for (const [s, t] of links) parent[find(s)] = find(t);
  return new Set(parent.map((_, i) => find(i))).size;
}

const el = (name, attrs = {}, parent) => {
  const node = document.createElementNS(SVG_NS, name);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
  if (parent) parent.appendChild(node);
  return node;
};

const html = (tag) => (text) => {
  const node = document.createElement(tag);
  node.textContent = text;
  return node;
};
const sub = html("sub"), sup = html("sup");
const fx = (v) => v.toFixed(1);

// A weight of the stationary vector, as a fraction when it is one with a small denominator.
function weight(x) {
  for (let q = 1; q <= 24; q++) {
    const pnum = Math.round(x * q);
    if (Math.abs(x - pnum / q) < 1e-6) return q === 1 ? String(pnum) : `${pnum}/${q}`;
  }
  return x.toFixed(3);
}

// How long one step of the random walk takes, in milliseconds.
const SPEEDS = { calm: 700, fast: 220, frantic: 90 };
// How long the spectral-law animation shows each prime.
const SPEC_MS = { calm: 900, fast: 450, frantic: 200 };
// How long the growing ball waits between layers, and how long it rests once it has
// found a conjugate pair before starting again somewhere else.
const BALL_MS = { calm: 1400, fast: 750, frantic: 330 };
const BALL_HOLD = 3200;

function init(root, index, base) {
  const svg = root.querySelector("svg");
  const tip = root.querySelector(".wm-graph-tip");
  const caption = root.querySelector(".wm-graph-caption");
  const stats = root.querySelector(".wm-graph-stats");
  const select = root.querySelector(".wm-graph-prime select");
  const trailSelect = root.querySelector(".wm-graph-trail-len");
  const speedSelect = root.querySelector(".wm-graph-speed");
  const viewSelect = root.querySelector(".wm-graph-view select");
  const animSelect = root.querySelector(".wm-graph-anim");
  const ballStatus = root.querySelector(".wm-graph-ball-status");
  const labelBox = root.querySelector(".wm-graph-labels input");
  const ellSelect = root.querySelector(".wm-graph-ell select");
  const legend5 = root.querySelector(".wm-graph-legend-5");
  const legendFrob = root.querySelector(".wm-graph-legend-frob");
  const legendCls = root.querySelector(".wm-graph-legend-cls");
  const eye = root.querySelector(".wm-graph-eye");
  const legendPair = root.querySelector(".wm-graph-legend-pair");
  const trailLabel = root.querySelector(".wm-trail-n");
  const buttons = [...root.querySelectorAll("[data-kind]")];
  const legend3 = root.querySelector(".wm-graph-legend-3");
  const legendFp = root.querySelector(".wm-graph-legend-fp");
  const fieldNote = root.querySelector(".wm-graph-field");
  const legendArrow = root.querySelector(".wm-graph-legend-arrow");
  const playBtn = root.querySelector(".wm-graph-play");
  const piBox = root.querySelector(".wm-anim-pi");
  const reduced = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  // Whether the animation is running. With "reduce motion" set, the page starts with a
  // still picture and the Play button starts the animation; otherwise it starts playing
  // and the button pauses it.
  let playing = !reduced;

  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  // A glow for the laser: blur the shape and lay the sharp shape on top. The filter
  // region is in user space, so straight (zero-height) segments still glow.
  const defs = el("defs", {}, svg);
  const glow = el("filter", { id: "wm-glow", filterUnits: "userSpaceOnUse", x: -40, y: -40,
    width: W + 80, height: H + 80 }, defs);
  el("feGaussianBlur", { in: "SourceGraphic", stdDeviation: 2.6, result: "b" }, glow);
  const merge = el("feMerge", {}, glow);
  for (const src of ["b", "b", "SourceGraphic"]) el("feMergeNode", { in: src }, merge);
  const gFrob = el("g", { class: "wm-g-frobs" }, svg);
  const gEdges = el("g", { class: "wm-g-edges" }, svg);
  const gNodes = el("g", { class: "wm-g-nodes" }, svg);
  const gLabels = el("g", { class: "wm-g-labels", "aria-hidden": "true" }, svg);
  // The laser: the beam from the vertex just left to the walker, a glow and a bright core.
  const beamGlow = el("path", { class: "wm-g-beam-glow" }, svg);
  const beamCore = el("path", { class: "wm-g-beam-core" }, svg);
  const walkerDot = el("circle", { class: "wm-g-walker", r: 5 }, svg);
  const spectral = createSpectral({
    svg, base, W, H, still: () => !playing,
    onStatus: (msg) => { if (mode === "spectral") stats.textContent = msg; },
  });

  const state = {
    kind: "2", p: {}, sim: null, nodes: [], strands: [], adj: [], pinned: null,
    view: null, dragging: false,
  };
  let visible = true, raf = 0, walker = null;
  // The walk's recent history: the strands of the last TRAIL steps, oldest first.
  // The newest is drawn boldest; older steps fade out.
  let TRAIL = Number(trailSelect.value), hopMs = SPEEDS[speedSelect.value];
  let trail = [], painted = [];
  // "walk" (a random walk with a fading trail) or "ball" (a ball grown layer by layer
  // from a random vertex until it contains a vertex and its Frobenius conjugate).
  const host = root.closest("[data-mode]");
  if (host && [...animSelect.options].some((o) => o.value === host.dataset.mode)) {
    animSelect.value = host.dataset.mode;
  }
  let mode = animSelect.value, ball = null;

  // Each prime's data is a separate small file, fetched the first time it is chosen.
  const cache = new Map();
  function load(p, kind) {
    const key = `p${p}-${kind}`;
    if (!cache.has(key)) {
      cache.set(key, fetch(`${base}/${key}.json`)
        .then((r) => { if (!r.ok) throw new Error(r.status); return r.json(); })
        .catch((e) => { cache.delete(key); throw e; }));
    }
    return cache.get(key);
  }

  // Stop whatever is running: the frame loop, the layout simulation, and every mark an
  // animation leaves on the picture. Every rebuild starts with this, so switching
  // animation (or graph, prime or view) never leaves the previous loop or its state.
  function stopAnimation() {
    if (raf) { cancelAnimationFrame(raf); raf = 0; }
    if (state.sim) state.sim.stop();
    spectral.hide();
    ball = null; walker = null;
    trail = []; paintTrail();
    beamGlow.removeAttribute("d"); beamCore.removeAttribute("d");
    beamGlow.style.display = beamCore.style.display = walkerDot.style.display = "none";
    if (piBox) piBox.textContent = "";
  }

  // The info panel. Each animation has its own block in the widget's markup
  // (.wm-anim-info [data-anim]); only the block of the running animation is shown.
  // On an Animations page the page title changes with it, and comes back when the
  // page's own animation is chosen again.
  const infoBlocks = [...root.querySelectorAll(".wm-anim-info [data-anim]")];
  const pageTitle = host ? document.querySelector("#title-block-header h1.title") : null;
  const ownTitle = pageTitle ? pageTitle.textContent : "", ownDocTitle = document.title;
  const ownMode = mode;
  function showInfo() {
    let block = null;
    for (const b of infoBlocks) { b.hidden = b.dataset.anim !== mode; if (!b.hidden) block = b; }
    if (!pageTitle || !block) return;
    const title = mode === ownMode ? ownTitle : block.dataset.title;
    pageTitle.textContent = title;
    document.title = mode === ownMode ? ownDocTitle : ownDocTitle.replace(ownTitle, title);
  }

  let buildId = 0;
  async function build() {
    const id = ++buildId;
    stopAnimation();
    root.dataset.mode = mode;
    showInfo();
    svg.classList.toggle("wm-spec-mode", mode === "spectral");
    if (mode === "spectral") {        // the spectral law: no graph is drawn
      hideTip(); showEmpty(null);
      ballStatus.hidden = true;
      caption.textContent = "Eigenvalues of the Hecke operator T_ℓ on weight-2 cusp forms " +
        "S₂(Γ₀(p)), for growing p, against the Kesten–McKay law (red).";
      stats.textContent = "";
      if (eye) {
        eye.href = `${eye.dataset.base}#spectral`;
        eye.dataset.tip = "Learn more about the spectral law";
        eye.setAttribute("aria-label", eye.dataset.tip);
      }
      await spectral.show(ellSelect.value);
      return;
    }
    const kind = state.kind, cfg = KINDS[kind];
    const cap = (index.caps && index.caps[kind]) || Infinity;
    const offered = (index.kind_primes && index.kind_primes[kind]) || index.primes;
    const primes = offered.filter((q) => q <= cap);
    if (!primes.includes(state.p[kind])) {
      state.p[kind] = primes.includes(cfg.defaultP) ? cfg.defaultP : primes[0];
    }
    // The growing ball runs on the full graph or on the quotient by Frobenius.
    const spineOpt = viewSelect.querySelector('option[value="spine"]');
    const overOpt = viewSelect.querySelector('option[value="over"]');
    spineOpt.disabled = mode === "ball";
    overOpt.hidden = overOpt.disabled = !cfg.level || mode === "ball";
    if (mode === "ball" && viewSelect.value === "spine") viewSelect.value = "full";
    if (overOpt.disabled && viewSelect.value === "over") viewSelect.value = cfg.level ? "full" : "spine";
    const view = viewSelect.value;
    // Some spines are empty (no vertex fixed by Frobenius): for example the level-3 graph
    // when p = 1 (mod 3). In the spine view those primes are marked in the menu, and
    // choosing one shows a short explanation with buttons to move on (see showEmpty).
    const empty = (q) => view === "spine" && index.spine && index.spine[q] && index.spine[q][kind] === 0;
    const p = state.p[kind];
    select.replaceChildren(...primes.map((q) => {
      const o = document.createElement("option");
      o.value = q; o.selected = q === p;
      o.textContent = empty(q) ? `${q} (empty F_p-subgraph)` : q;
      return o;
    }));
    const nearestNonempty = () => {
      const ok = primes.filter((q) => !empty(q));
      return ok.length ? ok.reduce((a, b) => (Math.abs(b - p) < Math.abs(a - p) ? b : a)) : null;
    };
    buttons.forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.kind === kind)));

    root.classList.add("wm-graph-loading");
    let entry;
    try { entry = await load(p, kind); } catch {
      if (id === buildId) {
        root.classList.remove("wm-graph-loading");
        caption.textContent = "Could not load this graph. Please try again.";
      }
      return;
    }
    if (id !== buildId) return;        // a newer choice has been made meanwhile
    root.classList.remove("wm-graph-loading");
    const g = entry.graph;
    const d = derive(g, view);
    const n = d.nodes.length;
    const nodes = d.nodes.map((nd, i) => ({
      ...nd, i,
      x: W / 2 + (W / 3) * Math.cos((2 * Math.PI * i) / n),
      y: H / 2 + (H / 3) * Math.sin((2 * Math.PI * i) / n),
    }));

    // One drawn strand per edge of the multigraph; strands between the same
    // pair of vertices fan out as curves.
    // A strand is [s, t, ell, both]. Where the numbers of edges s -> t and t -> s
    // agree, each pair of opposite edges is one plain strand ("both"), which can be
    // walked either way. Where they differ (next to a vertex with automorphisms other
    // than +-1, and everywhere in the directed Gamma_1(5) graph), every edge is its own
    // strand s -> t, drawn with an arrowhead and walked only forwards.
    const directed = !!g.directed;
    const groups = new Map();
    for (const link of d.links) {
      const { s, t, ell, out, back } = linkCounts(link, directed);
      const key = s <= t ? `${s}-${t}` : `${t}-${s}`;
      if (!groups.has(key)) groups.set(key, []);
      const list = groups.get(key);
      if (s !== t && out === back) {
        for (let k = 0; k < out; k++) list.push([s, t, ell, true]);
      } else {
        for (let k = 0; k < out; k++) list.push([s, t, ell, false]);
        for (let k = 0; k < back; k++) list.push([t, s, ell, false]);
      }
    }
    const strands = [], simLinks = [], adj = nodes.map(() => []);
    // Frobenius as dashed links (the "over" view); they also pull orbits together.
    gFrob.replaceChildren();
    const frobs = (d.frob || []).map(([a, b]) => {
      simLinks.push({ source: a, target: b });
      return { a, b, line: el("line", { class: "wm-g-frob" }, gFrob) };
    });
    for (const [key, list] of groups) {
      const [u, v] = key.split("-").map(Number);
      if (u !== v) simLinks.push({ source: u, target: v });
      list.forEach(([s, t, ell, both], k) => {
        // u, v: the endpoints in a fixed order, for drawing; s -> t: the direction;
        // arrow: drawn with an arrowhead at t (a loop has none)
        const strand = { s, t, u, v, ell, k, total: list.length, arrow: !both && s !== t };
        strands.push(strand);
        adj[s].push({ to: t, strand });
        if (both) adj[t].push({ to: s, strand });
      });
    }

    const r = n <= 20 ? 8 : n <= 45 ? 6.5 : n <= 110 ? 5 : n <= 300 ? 4 : n <= 700 ? 3 : 2.2;
    const rOf = (nd) => (nd.pair ? r * 1.35 : r);
    const arrowSize = n <= 45 ? 6 : n <= 110 ? 5 : n <= 300 ? 4.2 : 3;
    // The walk: its transition matrix, from the adjacency matrix of the graph shown.
    // "centres" are the vertices with automorphisms other than +-1, j = 0 and j = 1728,
    // known for the full graphs without level structure (see walk.js).
    const centres = view === "full" && g.aut ? Object.keys(g.aut).map(Number) : [];
    const P = walkFor(adjacency(n, d.links, directed), centres, g.aut);
    gEdges.replaceChildren();
    gNodes.replaceChildren();
    gLabels.replaceChildren();
    const labelled = n <= LABEL_MAX;
    labelBox.disabled = !labelled;
    labelBox.parentNode.title = labelled ? "" : `Labels are shown for graphs with at most ${LABEL_MAX} vertices; hover a vertex instead.`;
    for (const st of strands) {
      st.path = el("path", { class: `wm-g-edge wm-g-edge-${st.ell}` }, gEdges);
    }
    for (const nd of nodes) {
      nd.circle = el("circle", {
        class: `wm-g-node${nd.fp ? " wm-g-node-fp" : ""}${nd.pair ? " wm-g-node-pair" : ""}` +
          `${g.cls && view !== "folded" && g.cls[nd.members[0]] ? " wm-g-node-c1" : ""}`,
        r: rOf(nd), tabindex: "-1",
      }, gNodes);
      attachNodeEvents(nd);
      if (labelled) {
        nd.text = el("text", { class: "wm-g-label" }, gLabels);
        nd.text.textContent = shortJ(g.nodes, nd);
      }
    }
    gLabels.style.display = labelled && labelBox.checked ? "" : "none";
    Object.assign(state, { frobs, nodes, strands, adj, pinned: null, r, rOf, arrowSize, P,
                           graph: g, sigma: g.sigma, shownView: view, ballGraph: null });
    painted = [];
    hideTip();

    if (state.sim) state.sim.stop();
    const dist = Math.max(22, Math.min(90, 520 / Math.sqrt(n)));
    state.sim = forceSimulation(nodes)
      .force("link", forceLink(simLinks).distance(dist).strength(0.5))
      .force("charge", forceManyBody().strength(-Math.max(60, 5200 / n)).distanceMax(320))
      .force("collide", forceCollide((nd) => rOf(nd) + (n > 300 ? 1.5 : 3)))
      .force("x", forceX(W / 2).strength(0.06))
      .force("y", forceY(H / 2).strength(0.075))
      .alphaDecay(n > 700 ? 0.06 : n > 300 ? 0.03 : 0.035)
      .on("tick", draw)
      .stop();
    // Large graphs do most of their layout up front, so they settle in a moment.
    if (n > 700 && !reduced) state.sim.alpha(0.4);
    state.sim.tick(reduced ? 300 : n > 700 ? 160 : n > 300 ? 120 : 40);
    svg.classList.toggle("wm-dense", n > 300);
    state.view = null;          // snap the view to the new graph
    draw();
    if (!reduced && visible) state.sim.restart();

    const nFull = g.nodes.length;
    const nComp = components(n, d.links);
    // The random walk needs a connected graph: on a disconnected F_p-subgraph it is off.
    state.walkOff = (view === "spine" || view === "over") && (n === 0 || nComp > 1);
    if (view === "folded") {
      caption.textContent = cfg.caption(p, nFull).replace(/ in characteristic/, ", quotient by Frobenius, in characteristic")
        .replace(/: .*$/, `: ${n} vertices, the orbits of Frobenius on ${nFull} vertices.`);
      const bySize = {};
      for (const nd of d.nodes) bySize[nd.members.length] = (bySize[nd.members.length] || 0) + 1;
      stats.textContent = `${bySize[1] || 0} vertices fixed by Frobenius, ${bySize[2] || 0} orbits of size 2` +
        (bySize[4] ? `, ${bySize[4]} of size 4 (here Frobenius² is the diamond operator ⟨p⟩ ≠ 1)` : "") +
        " (drawn larger). Edges are Frobenius-orbits of edges.";
    } else if (view === "over") {
      const nFix = d.nodes.filter((nd) => nd.fp).length;
      caption.textContent = cfg.caption(p, nFull).replace(/ in characteristic/, ", over the F_p-subgraph, in characteristic")
        .replace(/: .*$/, `: the subgraph induced on the ${n} of ${nFull} vertices (E, ·) with j(E) in F_p. ` +
          `Gold: the ${nFix} fixed by Frobenius (the F_p-subgraph itself); dashed: Frobenius.`);
      stats.textContent = (nComp === 1 ? "This subgraph is connected." :
        `This subgraph has ${nComp} connected components.`) +
        ` Frobenius moves ${n - nFix} of its vertices.` +
        (state.walkOff && mode === "walk" ? " The random walk is off, since it is not connected." : "");
    } else if (view === "spine") {
      caption.textContent = cfg.caption(p, nFull).replace(/ in characteristic/, ", F_p-subgraph, in characteristic")
        .replace(/: .*$/, `: the subgraph induced on the ${n} of ${nFull} vertices fixed by Frobenius.`);
      stats.textContent = (n === 0 ? "The F_p-subgraph is empty." :
        nComp === 1 ? "The F_p-subgraph is connected." :
        nComp === n ? `The F_p-subgraph is totally disconnected: ${n} vertices and no edges between them.` :
        `The F_p-subgraph has ${nComp} connected components.`) +
        (state.walkOff && mode === "walk" && n > 0 ? " The random walk is off, since the graph is not connected." : "");
    } else {
      caption.textContent = cfg.caption(p, n);
      const lam = g.lambda2;
      stats.textContent = (cfg.spectrum ? cfg.spectrum() + " " : lam == null ? "" :
        `Second-largest eigenvalue modulus ${lam.toFixed(3)} ≤ ${cfg.bound} ≈ ${cfg.boundVal.toFixed(3)}. `) +
        `${g.fp.filter(Boolean).length} of the ${n} vertices are fixed by Frobenius (the F_p-subgraph).`;
    }
    legend3.hidden = !cfg.ells.includes(3);
    if (legend5) legend5.hidden = !cfg.ells.includes(5);
    if (legendFp) legendFp.hidden = !d.nodes.some((nd) => nd.fp);
    if (legendPair) legendPair.hidden = view !== "folded";
    if (legendFrob) legendFrob.hidden = view !== "over";
    if (legendCls) legendCls.hidden = !g.cls || view === "folded";
    if (legendArrow) legendArrow.hidden = !strands.some((st) => st.arrow);
    state.nComp = nComp;
    const nEdges = d.links.filter(([a, b]) => a !== b).length;
    showEmpty(view === "spine" && (n === 0 || nEdges === 0)
      ? { kind, p, n, alternative: nearestNonempty() } : null);
    const nr = entry.nonresidue;
    fieldNote.replaceChildren("Labels are j-invariants in F", sub("p"), "(i), i",
      sup("2"), ` = ${nr < 0 ? "−" + (-nr) : nr}.`);
    root.dataset.kind = kind;
    if (eye) {
      eye.href = `${eye.dataset.base}#graph-${kind}`;
      eye.dataset.tip = view === "full" ? "Learn more about this graph" :
        `Learn more about this graph and its ${VIEW_NAME[view]}`;
      eye.setAttribute("aria-label", eye.dataset.tip);
    }

    const movable = nodes.filter((nd) => adj[nd.i].length);
    const bloom = mode === "ball" || mode === "cycle";
    svg.classList.toggle("wm-ball-mode", bloom);
    ball = null;
    ballStatus.hidden = !bloom;
    if (bloom) {
      if (!playing) {             // no animation: show one finished ball
        startBloom(0);
        while (!ball.done) stepBloom(0);
      } else ballStatus.textContent = "";
    }
    // A fresh, uniformly random starting vertex every time a graph is drawn.
    const seed = new Uint32Array(1);
    crypto.getRandomValues(seed);
    walker = { at: movable.length ? movable[seed[0] % movable.length].i : 0,
               next: null, strand: null, t0: 0 };
    const walkOn = strands.length && mode === "walk" && !state.walkOff;
    walkerDot.style.display = walkOn ? "" : "none";
    beamGlow.style.display = beamCore.style.display = walkOn ? "" : "none";
    if (walkOn) { parkWalker(); showStationary(); }
  }

  // The walker at rest on its vertex: the still picture before Play is pressed.
  function parkWalker() {
    const nd = walker && walker.next == null ? state.nodes[walker.at] : null;
    if (!nd || nd.sx === undefined) return;
    walkerDot.setAttribute("cx", fx(nd.sx));
    walkerDot.setAttribute("cy", fx(nd.sy));
  }

  // The stationary vector pi of the walk (pi P = pi), computed numerically and shown
  // in the info panel scaled so that most vertices have weight 1.
  function showStationary() {
    if (!piBox) return;
    piBox.textContent = "";
    if (state.nComp !== 1 || state.nodes.length > 6000) return;
    const pi = stationary(state.P, { maxIter: 4000 });
    if (!pi) return;
    const { values, special } = relativeStationary(pi);
    const name = (v) => {
      const nd = state.nodes[v], j = shortJ(state.graph.nodes, nd);
      const order = state.shownView === "full" && state.graph.aut ? state.graph.aut[v] : 0;
      return `j = ${j}${order === 4 && j !== "1728" ? " (that is, 1728)" : ""}`;
    };
    if (!special.length) {
      piBox.textContent = "Stationary distribution: uniform. In the long run every vertex is visited equally often.";
    } else if (special.length <= 6) {
      piBox.textContent = "Stationary distribution, with the weight of an ordinary vertex set to 1: " +
        special.map((v) => `${weight(values[v])} at ${name(v)}`).join(", ") + ".";
    } else {
      let lo = Infinity, hi = -Infinity;
      for (const x of values) { lo = Math.min(lo, x); hi = Math.max(hi, x); }
      piBox.textContent = `Stationary distribution, with the most common weight set to 1: ` +
        `${special.length} of the ${values.length} vertices have another weight, between ` +
        `${weight(lo)} and ${weight(hi)}.`;
    }
  }

  // The simulation runs in its own coordinates; the view scales and centres it
  // so that the graph always fills the frame. (Isogeny graphs are expanders, so
  // their natural layout is a compact ball.)
  function fitView() {
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (const nd of state.nodes) {
      x0 = Math.min(x0, nd.x); x1 = Math.max(x1, nd.x);
      y0 = Math.min(y0, nd.y); y1 = Math.max(y1, nd.y);
    }
    const k = Math.min((W - 2 * PAD) / Math.max(x1 - x0, 1), (H - 2 * PAD) / Math.max(y1 - y0, 1), 2.4);
    const target = { k, tx: W / 2 - k * (x0 + x1) / 2, ty: H / 2 - k * (y0 + y1) / 2 };
    if (!state.view || reduced) state.view = target;
    else if (!state.dragging) {
      for (const key of ["k", "tx", "ty"]) state.view[key] += 0.18 * (target[key] - state.view[key]);
    }
    const v = state.view;
    for (const nd of state.nodes) { nd.sx = v.k * nd.x + v.tx; nd.sy = v.k * nd.y + v.ty; }
  }

  // Control point (or loop geometry) for a strand, shared by drawing and the walker.
  function geometry(st) {
    const a = state.nodes[st.u], b = state.nodes[st.v];
    if (st.u === st.v) {
      let ux = a.sx - W / 2, uy = a.sy - H / 2;
      const len = Math.hypot(ux, uy) || 1;
      ux /= len; uy /= len;
      const size = (20 + 9 * st.k) * 1.6, spread = 0.6;
      const cs = Math.cos(spread), sn = Math.sin(spread);
      return {
        loop: true, a,
        c1: [a.sx + size * (ux * cs - uy * sn), a.sy + size * (uy * cs + ux * sn)],
        c2: [a.sx + size * (ux * cs + uy * sn), a.sy + size * (uy * cs - ux * sn)],
      };
    }
    const off = (st.k - (st.total - 1) / 2) * 13;
    const dx = b.sx - a.sx, dy = b.sy - a.sy, len = Math.hypot(dx, dy) || 1;
    return { a, b, c: [(a.sx + b.sx) / 2 - (dy / len) * off, (a.sy + b.sy) / 2 + (dx / len) * off] };
  }

  function pointOn(geo, t, reverse) {
    if (reverse) t = 1 - t;
    const u = 1 - t;
    if (geo.loop) {
      const { a, c1, c2 } = geo;
      return [u * u * u * a.sx + 3 * u * u * t * c1[0] + 3 * u * t * t * c2[0] + t * t * t * a.sx,
              u * u * u * a.sy + 3 * u * u * t * c1[1] + 3 * u * t * t * c2[1] + t * t * t * a.sy];
    }
    const { a, b, c } = geo;
    return [u * u * a.sx + 2 * u * t * c[0] + t * t * b.sx, u * u * a.sy + 2 * u * t * c[1] + t * t * b.sy];
  }

  function draw() {
    fitView();
    const f = (v) => v.toFixed(1);
    for (const nd of state.nodes) {
      nd.circle.setAttribute("cx", f(nd.sx));
      nd.circle.setAttribute("cy", f(nd.sy));
    }
    if (gLabels.style.display !== "none") {
      const off = state.r + 2.5;
      for (const nd of state.nodes) {
        if (!nd.text) continue;
        nd.text.setAttribute("x", f(nd.sx + off));
        nd.text.setAttribute("y", f(nd.sy - off * 0.6));
      }
    }
    for (const fl of state.frobs || []) {
      const a = state.nodes[fl.a], b = state.nodes[fl.b];
      fl.line.setAttribute("x1", f(a.sx)); fl.line.setAttribute("y1", f(a.sy));
      fl.line.setAttribute("x2", f(b.sx)); fl.line.setAttribute("y2", f(b.sy));
    }
    for (const st of state.strands) {
      const g = geometry(st);
      st.path.setAttribute("d", g.loop
        ? `M${f(g.a.sx)},${f(g.a.sy)}C${f(g.c1[0])},${f(g.c1[1])} ${f(g.c2[0])},${f(g.c2[1])} ${f(g.a.sx)},${f(g.a.sy)}`
        : `M${f(g.a.sx)},${f(g.a.sy)}Q${f(g.c[0])},${f(g.c[1])} ${f(g.b.sx)},${f(g.b.sy)}` +
          (st.arrow ? arrowHead(st, g) : ""));
    }
    parkWalker();
  }

  // An open arrowhead at the target end of a directed strand, just outside the target's
  // circle. It is part of the strand's own path, so it is coloured, dimmed and
  // highlighted together with it.
  function arrowHead(st, g) {
    const { a, b, c } = g, toB = st.t === st.v;       // the head is at b, or at a
    const len = Math.hypot(b.sx - a.sx, b.sy - a.sy) || 1;
    const gap = state.rOf(toB ? b : a) + 1.5;
    const t = toB ? 1 - Math.min(0.45, gap / len) : Math.min(0.45, gap / len);
    const [x, y] = pointOn(g, t);
    let dx = 2 * (1 - t) * (c[0] - a.sx) + 2 * t * (b.sx - c[0]);   // direction of travel
    let dy = 2 * (1 - t) * (c[1] - a.sy) + 2 * t * (b.sy - c[1]);
    if (!toB) { dx = -dx; dy = -dy; }
    const norm = Math.hypot(dx, dy) || 1, h = state.arrowSize, w = 0.5 * h;
    dx /= norm; dy /= norm;
    const bx = x - dx * h, by = y - dy * h;
    return `M${fx(bx - dy * w)},${fx(by + dx * w)}L${fx(x)},${fx(y)}L${fx(bx + dy * w)},${fx(by - dx * w)}`;
  }

  // ---- random walk -------------------------------------------------------
  function positionWalker(now) {
    if (!walker || !state.nodes.length || state.walkOff) return;
    if (!playing) return;
    if (walker.next == null) {
      // The next vertex is drawn from the walker's row of the transition matrix P, and
      // from nothing else. Among parallel strands to that vertex, which one lights up
      // is a matter of drawing only.
      const row = state.P[walker.at];
      if (!row || !row.size) return;        // a vertex with no edge out (a tiny spine)
      const to = sampleStep(row, Math.random());
      const options = state.adj[walker.at].filter((o) => o.to === to);
      if (!options.length) return;
      const pick = options[Math.floor(Math.random() * options.length)];
      walker.next = pick.to; walker.strand = pick.strand; walker.t0 = now;
      walker.reverse = pick.strand.u !== walker.at;
      trail.push(pick.strand);
      if (trail.length > TRAIL) trail.shift();
      paintTrail();
    }
    let t = (now - walker.t0) / hopMs;
    if (t >= 1) {
      walker.at = walker.next; walker.next = null;
      t = 1;
    }
    const ease = 0.5 * t + 0.5 * t * t * (3 - 2 * t);   // quick, with a slight pause at each curve
    const geo = geometry(walker.strand);
    const [x, y] = pointOn(geo, ease, walker.reverse);
    walkerDot.setAttribute("cx", x.toFixed(1));
    walkerDot.setAttribute("cy", y.toFixed(1));
    // the beam: the part of the edge already travelled, shrinking from its tail
    const tail = Math.max(0, ease - 0.55);
    let d = "";
    for (let k = 0; k <= 14; k++) {
      const [bx, by] = pointOn(geo, tail + ((ease - tail) * k) / 14, walker.reverse);
      d += `${k ? "L" : "M"}${bx.toFixed(1)},${by.toFixed(1)}`;
    }
    beamGlow.setAttribute("d", d);
    beamCore.setAttribute("d", d);
  }
  // Strength 1 for the step being taken now, fading to 0.15 for the oldest in the window.
  // A strand walked more than once shows its most recent step.
  function paintTrail() {
    for (const st of painted) {
      st.path.classList.remove("wm-g-trail");
      st.path.style.removeProperty("--trail");
    }
    for (const st of painted) st.path.classList.remove("wm-g-hot");
    trail.forEach((st, k) => {
      const lin = (k + 1 + TRAIL - trail.length) / TRAIL;
      st.path.classList.add("wm-g-trail");
      st.path.classList.toggle("wm-g-hot", k >= trail.length - 3);
      st.path.style.setProperty("--trail", (0.15 + 0.85 * lin).toFixed(2));
    });
    painted = trail.slice();
  }
  // ---- "?" buttons: a small box with the relevant entries of the glossary ------
  const gloss = root.querySelector(".wm-graph-glossary");
  const pop = root.querySelector(".wm-graph-pop");
  let popFor = null;
  function closePop() {
    if (!popFor) return;
    popFor.setAttribute("aria-expanded", "false");
    pop.hidden = true;
    popFor = null;
  }
  function openPop(btn) {
    if (popFor === btn) { closePop(); return; }
    closePop();
    const dl = document.createElement("dl");
    for (const key of btn.dataset.terms.split(" ")) {
      const dt = gloss.querySelector(`dt[data-term="${key}"]`);
      if (dt) dl.append(dt.cloneNode(true), dt.nextElementSibling.cloneNode(true));
    }
    const close = button("×", closePop);
    close.className = "wm-graph-pop-close";
    close.setAttribute("aria-label", "Close");
    const all = button("All terms ↓", () => {
      closePop();
      gloss.open = true;
      gloss.scrollIntoView({ block: "nearest", behavior: reduced ? "auto" : "smooth" });
    });
    all.className = "wm-graph-pop-all";
    pop.replaceChildren(close, dl, all);
    pop.hidden = false;
    popFor = btn;
    btn.setAttribute("aria-expanded", "true");
    // Just under the button, kept inside the figure.
    const b = btn.getBoundingClientRect(), h = root.getBoundingClientRect();
    pop.style.top = `${b.bottom - h.top + 8}px`;
    pop.style.left = `${Math.max(0, Math.min(h.width - pop.offsetWidth, b.left - h.left - 24))}px`;
    close.focus({ preventScroll: true });
  }
  root.querySelectorAll(".wm-q").forEach((btn) => btn.addEventListener("click", (ev) => {
    ev.stopPropagation();
    openPop(btn);
  }));
  document.addEventListener("click", (ev) => { if (popFor && !pop.contains(ev.target)) closePop(); });
  document.addEventListener("keydown", (ev) => {
    if (ev.key === "Escape" && popFor) { const b = popFor; closePop(); b.focus(); }
  });

  // ---- empty or totally disconnected spines --------------------------------
  // With the spine induced on the vertices with j in F_p it is never empty (every p
  // has supersingular j-invariants in F_p), but the message is kept as a safeguard.
  const emptyBox = document.createElement("div");
  emptyBox.className = "wm-graph-empty";
  emptyBox.hidden = true;
  svg.parentNode.appendChild(emptyBox);
  function button(text, act) {
    const b = document.createElement("button");
    b.type = "button"; b.textContent = text;
    b.addEventListener("click", act);
    return b;
  }
  function showEmpty(info) {
    emptyBox.hidden = !info;
    emptyBox.classList.toggle("wm-graph-empty-top", !!info && info.n > 0);
    if (!info) return;
    const { kind, p, n, alternative } = info;
    const head = document.createElement("strong");
    const text = document.createElement("div");
    text.className = "wm-graph-empty-text";
    if (n === 0) {
      head.textContent = "Empty F_p-subgraph";
      text.append(`No vertex of this graph is fixed by Frobenius for p = ${p}, so the F`, sub("p"),
        "-subgraph is empty (not just disconnected).");
      const reason = EMPTY_WHY[kind] && EMPTY_WHY[kind](p);
      if (reason) {
        // The reason, folded away on small screens so the buttons stay in view.
        const why = document.createElement("details");
        why.open = window.matchMedia("(min-width: 576px)").matches;
        const summary = document.createElement("summary");
        summary.textContent = "Why?";
        why.append(summary, reason);
        text.append(why);
      }
    } else {
      head.textContent = "Totally disconnected F_p-subgraph";
      text.append(`For p = ${p} the F`, sub("p"), `-subgraph has ${n} vertices but no edges between them.`);
    }
    const go = document.createElement("div");
    go.className = "wm-graph-empty-go";
    if (alternative) {
      go.appendChild(button(`Show p = ${alternative}`, () => {
        state.p[state.kind] = alternative; build().then(startLoop);
      }));
    }
    go.appendChild(button("Full graph", () => { viewSelect.value = "full"; build().then(startLoop); }));
    go.appendChild(button("Folded graph", () => { viewSelect.value = "folded"; build().then(startLoop); }));
    emptyBox.replaceChildren(head, text, go);
  }

  // ---- growing ball -------------------------------------------------------
  // Breadth-first search, always run in the full graph ("upstairs") from a centre v, one
  // layer at a time, and drawn in the view shown.
  //  * Full graph: stop when the ball contains a vertex w and its Frobenius conjugate
  //    w^(p) != w. Frobenius preserves distances, so this happens exactly when B_r(v) and
  //    B_r(v^(p)) share a vertex off the F_p-subgraph.
  //  * Quotient by Frobenius: the ball downstairs is the image of the ball upstairs, and
  //    it remembers which lift it reached each vertex by. It stops when it reaches a
  //    vertex a second time by the other lift (a self-collision: w and w^(p) are both in
  //    the ball upstairs) or when it reaches the F_p-subgraph.
  const shortLabel = (i) => {
    const l = state.graph.nodes[i];
    return Array.isArray(l) ? `(${l[0]}, ${l[1]})` : `j = ${l}`;
  };
  function ballGraph() {
    if (state.ballGraph) return state.ballGraph;
    const g = state.graph, n = g.nodes.length, up = Array.from({ length: n }, () => new Set());
    for (const [s, t] of g.links) {
      if (s === t) continue;
      up[s].add(t);
      if (!g.directed) up[t].add(s);
    }
    const members = new Array(n);           // the Frobenius orbit of each vertex
    for (const cyc of cycles(g.sigma)) for (const m of cyc) members[m] = cyc;
    const orbitOf = new Array(n);
    for (const nd of state.nodes) for (const m of nd.members) orbitOf[m] = nd.i;
    const strandOf = new Map();
    for (const st of state.strands) {
      if (st.u === st.v) continue;
      const key = `${st.u}-${st.v}`;
      if (!strandOf.has(key)) strandOf.set(key, st);
    }
    state.ballGraph = { up: up.map((x) => [...x].sort((u, v) => u - v)), orbitOf, strandOf, members };
    return state.ballGraph;
  }
  function startBall(now, at) {
    const g = state.graph, quotient = state.shownView === "folded", fixed = g.fp;
    let center;
    if (at != null) center = state.nodes[at].members[0];
    else {
      const all = [...g.nodes.keys()], off = all.filter((i) => !fixed[i]);
      const pool = quotient && off.length ? off : all;
      const seed = new Uint32Array(1);
      crypto.getRandomValues(seed);
      center = pool[seed[0] % pool.length];
    }
    ball = { center, quotient, dist: new Map([[center, 0]]), parent: new Map(), layers: [[center]],
             r: 0, t0: now, done: null, holdUntil: 0 };
    if (quotient && fixed[center]) { ball.done = { fp: center }; ball.holdUntil = now + BALL_HOLD; }
    paintBall();
  }
  function stepBall(now) {
    const b = ball, g = state.graph, fixed = g.fp, next = [];
    for (const u of b.layers[b.r]) {
      for (const w of ballGraph().up[u]) {
        if (b.dist.has(w)) continue;
        b.dist.set(w, b.r + 1);
        b.parent.set(w, u);
        next.push(w);
      }
    }
    if (!next.length) {
      b.done = { none: true };
    } else {
      b.layers.push(next);
      b.r += 1;
      // Any new pair has a member in the new layer. Prefer the one closest to the centre.
      let pair = null, fp = null;
      for (const w of next) {
        // another member of w's Frobenius orbit already in the ball?
        for (const c of ballGraph().members[w]) {
          if (c === w || !b.dist.has(c)) continue;
          const score = b.dist.get(w) + b.dist.get(c);
          if (!pair || score < pair.score) pair = { w, c, score };
        }
        if (b.quotient && fp == null && fixed[w]) fp = w;
      }
      if (pair || fp != null) b.done = { pair, fp };
    }
    if (b.done) b.holdUntil = now + BALL_HOLD;
    paintBall();
  }
  // The two kinds of blooming share the clock: "ball" = blooming for conjugates,
  // "cycle" = blooming for collision.
  const startBloom = (now, at) => (mode === "cycle" ? startCycle(now, at) : startBall(now, at));
  const stepBloom = (now) => (mode === "cycle" ? stepCycle(now) : stepBall(now));
  function growBall(now) {
    if (!state.nodes.length) return;
    if (!ball) { startBloom(now); return; }
    if (ball.done) { if (now >= ball.holdUntil) startBloom(now); return; }
    if (now - ball.t0 >= BALL_MS[speedSelect.value]) { ball.t0 = now; stepBloom(now); }
  }

  // ---- blooming for collision ----------------------------------------------
  // Breadth-first search in the graph shown, one layer at a time, until an edge reaches
  // a vertex already in the ball other than along the edge it came in by: the ball
  // has closed a cycle (a loop is a cycle of length 1, two parallel edges one of
  // length 2). In the directed Gamma_1(5) graph the search follows arcs forwards, and
  // stops when two shortest paths from the centre meet.
  const nodeLabel = (i) => {
    const nd = state.nodes[i], l = shortJ(state.graph.nodes, nd);
    return Array.isArray(nd.label) && !nd.pair ? `(${nd.label[0]}, ${nd.label[1]})` : `j = ${l}`;
  };
  function startCycle(now, at) {
    const n = state.nodes.length;
    let center = at;
    if (center == null) {
      const seed = new Uint32Array(1);
      crypto.getRandomValues(seed);
      center = seed[0] % n;
    }
    ball = { kind: "cycle", center, dist: new Map([[center, 0]]), parent: new Map(),
             pstrand: new Map(), layers: [[center]], r: 0, t0: now, done: null, holdUntil: 0 };
    const loop = state.adj[center].find(({ to }) => to === center);
    if (loop) {
      ball.done = { u: center, w: center, strand: loop.strand, len: 1 };
      ball.holdUntil = now + BALL_HOLD;
    }
    paintCycle();
  }
  function cycleOf(b, u, w) {
    // the tree paths from u and w back to their last common ancestor
    const up = (x) => { const out = [x]; while (b.parent.has(x)) { x = b.parent.get(x); out.push(x); } return out; };
    const pu = up(u), pw = up(w), onW = new Set(pw);
    const lca = pu.find((x) => onW.has(x));
    const vs = [...pu.slice(0, pu.indexOf(lca) + 1), ...pw.slice(0, pw.indexOf(lca))];
    const ss = [];
    for (const x of [...pu.slice(0, pu.indexOf(lca)), ...pw.slice(0, pw.indexOf(lca))]) ss.push(b.pstrand.get(x));
    return { vs, ss, len: pu.indexOf(lca) + pw.indexOf(lca) + 1 };
  }
  function stepCycle(now) {
    const b = ball, next = [];
    let best = null;
    for (const u of b.layers[b.r]) {
      for (const { to: w, strand } of state.adj[u]) {
        if (strand === b.pstrand.get(u)) continue;          // the edge it came in by
        if (!b.dist.has(w)) {
          b.dist.set(w, b.r + 1);
          b.parent.set(w, u);
          b.pstrand.set(w, strand);
          next.push(w);
          continue;
        }
        if (strand === b.pstrand.get(w)) continue;          // the tree edge, seen from above
        const c = cycleOf(b, u, w);
        if (!best || c.len < best.len) best = { u, w, strand, len: c.len };
      }
    }
    if (next.length) { b.layers.push(next); b.r += 1; }
    if (best) b.done = best;
    else if (!next.length) b.done = { none: true };
    if (b.done) b.holdUntil = now + BALL_HOLD;
    paintCycle();
  }
  function paintCycle() {
    const b = ball, n = state.nodes.length;
    let cyc = { vs: [], ss: [] };
    if (b.done && !b.done.none) {
      cyc = b.done.len === 1 ? { vs: [b.center], ss: [] } : cycleOf(b, b.done.u, b.done.w);
      cyc.ss.push(b.done.strand);
    }
    const hit = new Set(b.done && !b.done.none ? [b.done.u, b.done.w] : []);
    const geo = new Set(cyc.ss);
    for (const nd of state.nodes) {
      const d = b.dist.get(nd.i), c = nd.circle;
      c.classList.toggle("wm-g-out", d === undefined);
      if (nd.text) nd.text.classList.toggle("wm-g-out", d === undefined);
      c.classList.toggle("wm-g-front", d !== undefined && d === b.r && b.r > 0 && !b.done);
      c.classList.toggle("wm-g-center", nd.i === b.center);
      c.classList.toggle("wm-g-hit", hit.has(nd.i));
      const r0 = nd.pair ? state.r * 1.35 : state.r;
      c.setAttribute("r", hit.has(nd.i) || nd.i === b.center ? r0 * 1.6 : r0);
    }
    for (const st of state.strands) {
      st.path.classList.toggle("wm-g-inball", b.dist.has(st.s) && b.dist.has(st.t));
      st.path.classList.toggle("wm-g-geo", geo.has(st));
    }
    const size = `${b.dist.size} of ${n} vertices`, start = nodeLabel(b.center);
    if (b.done && b.done.none) {
      ballStatus.textContent = `The ball around ${start} filled its component (${size}) without closing a cycle.`;
    } else if (b.done) {
      const what = b.done.len === 1 ? "a loop (a cycle of length 1)" :
        b.done.len === 2 ? "two parallel edges (a cycle of length 2)" : `a cycle of length ${b.done.len}`;
      ballStatus.textContent = state.graph.directed
        ? `Radius ${b.r}: two shortest paths from ${start} meet at ${nodeLabel(b.done.w)} (${size}), ` +
          `closing ${what} if directions are ignored. Bold: the two paths.`
        : `Radius ${b.r}: the ball around ${start} (${size}) collides with itself and closes ${what}. ` +
          "Red: where it closed; bold: the cycle.";
    } else {
      ballStatus.textContent = `Blooming around ${start}: radius ${b.r}, ${size}.`;
    }
  }
  function paintBall() {
    const b = ball, G = ballGraph(), o = G.orbitOf, n = state.nodes.length;
    const shown = new Map();                       // displayed vertex -> distance
    for (const [u, d] of b.dist) if (!shown.has(o[u]) || d < shown.get(o[u])) shown.set(o[u], d);
    const ends = [];                               // upstairs vertices where it stopped
    if (b.done && b.done.pair) ends.push(b.done.pair.w, b.done.pair.c);
    if (b.done && b.done.fp != null) ends.push(b.done.fp);
    const hit = new Set(ends.map((u) => o[u]));
    const geo = new Set();                         // the strands of the shortest paths
    for (const x of ends) {
      for (let y = x; b.parent.has(y); y = b.parent.get(y)) {
        const a = o[y], c = o[b.parent.get(y)];
        const st = G.strandOf.get(a < c ? `${a}-${c}` : `${c}-${a}`);
        if (st) geo.add(st);
      }
    }
    const centre = o[b.center];
    for (const nd of state.nodes) {
      const d = shown.get(nd.i), c = nd.circle;
      c.classList.toggle("wm-g-out", d === undefined);
      if (nd.text) nd.text.classList.toggle("wm-g-out", d === undefined);
      c.classList.toggle("wm-g-front", d !== undefined && d === b.r && b.r > 0 && !b.done);
      c.classList.toggle("wm-g-center", nd.i === centre);
      c.classList.toggle("wm-g-hit", hit.has(nd.i));
      const r0 = nd.pair ? state.r * 1.35 : state.r;
      c.setAttribute("r", hit.has(nd.i) || nd.i === centre ? r0 * 1.6 : r0);
    }
    for (const st of state.strands) {
      st.path.classList.toggle("wm-g-inball", shown.has(st.s) && shown.has(st.t));
      st.path.classList.toggle("wm-g-geo", geo.has(st));
    }
    const where = b.quotient ? " in the quotient by Frobenius" : "";
    const size = `${shown.size} of ${n} vertices${where}; √n ≈ ${Math.sqrt(n).toFixed(1)}`;
    const start = shortLabel(b.center);
    let msg;
    if (b.done && b.done.none) {
      msg = `The ball around ${start} filled its connected component (${size}) without stopping.`;
    } else if (b.done && b.r === 0) {
      msg = `The centre ${start} is already in the F_p-subgraph. Tap another vertex.`;
    } else if (b.done) {
      const parts = [];
      if (b.done.pair) {
        const { w, c } = b.done.pair;
        parts.push(b.quotient
          ? `meets itself: it reached the vertex {${shortLabel(w)}, ${shortLabel(c)}} by both lifts`
          : `contains the conjugate pair ${shortLabel(w)} and ${shortLabel(c)}`);
      }
      if (b.done.fp != null) parts.push(`reaches the F_p-subgraph at ${shortLabel(b.done.fp)}`);
      msg = `Radius ${b.r}: the ball around ${start} (${size}) ${parts.join(", and ")}. ` +
        "Red: where it stopped; bold: shortest paths from the centre.";
    } else {
      msg = `Growing a ball around ${start}${where}: radius ${b.r}, ${size}.`;
    }
    ballStatus.textContent = msg;
  }

  function frame(now) {
    raf = 0;
    if (!playing || !visible || document.hidden) return;
    if (mode === "ball" || mode === "cycle") growBall(now);
    else if (mode === "spectral") spectral.step(now, SPEC_MS[speedSelect.value]);
    else if (mode === "walk") positionWalker(now);
    raf = requestAnimationFrame(frame);
  }
  // The frame loop runs only while something moves: not in the plain graph viewer, not
  // when paused, and (see frame) not off screen or in a hidden tab.
  function startLoop() {
    if (playing && mode !== "view" && !raf) raf = requestAnimationFrame(frame);
  }
  function setPlaying(on) {
    playing = on;
    if (playBtn) {
      playBtn.textContent = on ? "Pause" : "Play";
      playBtn.setAttribute("aria-label", on ? "Pause the animation" : "Play the animation");
    }
    if (!on) { if (raf) { cancelAnimationFrame(raf); raf = 0; } return; }
    // Starting from a still picture: the blooming and the spectral law begin again.
    if (mode === "ball" || mode === "cycle") { if (ball && ball.done) ball = null; }
    else if (mode === "spectral") spectral.restart();
    else if (walker && walker.next != null) walker.t0 = performance.now();
    startLoop();
  }

  // ---- hover, tap and drag -------------------------------------------------
  function labelOf(nd) {
    const tip = KINDS[state.kind].tip;
    const one = (label) => (Array.isArray(label)
      ? (tip ? tip(label) : `(j(E), j(E/C)) = (${label[0]}, ${label[1]})`) : `j = ${label}`);
    if (!nd.pair) return one(nd.label);
    return `Frobenius orbit: ${nd.members.map((m) => one(state.graph.nodes[m])).join("  ~  ")}`;
  }
  function showTip(nd) {
    const nbrs = new Set(state.adj[nd.i].map((e) => e.to));
    nbrs.add(nd.i);
    const hot = new Set(state.adj[nd.i].map((e) => e.strand));
    state.nodes.forEach((m) => {
      m.circle.classList.toggle("wm-g-dim", !nbrs.has(m.i));
      if (m.text) m.text.classList.toggle("wm-g-dim", !nbrs.has(m.i));
    });
    state.strands.forEach((s) => s.path.classList.toggle("wm-g-dim", !hot.has(s)));
    nd.circle.classList.add("wm-g-focus");
    tip.replaceChildren(labelOf(nd), ...(nd.fp ? [" · fixed by Frobenius"] : []));
    tip.hidden = false;
    const box = svg.getBoundingClientRect(), host = root.getBoundingClientRect();
    const px = box.left - host.left + (nd.sx / W) * box.width;
    const py = box.top - host.top + (nd.sy / H) * box.height;
    tip.style.left = `${Math.max(8, Math.min(host.width - tip.offsetWidth - 8, px - tip.offsetWidth / 2))}px`;
    tip.style.top = `${Math.max(0, py - tip.offsetHeight - 14)}px`;
  }
  function hideTip() {
    tip.hidden = true;
    state.nodes.forEach((m) => {
      m.circle.classList.remove("wm-g-dim", "wm-g-focus");
      if (m.text) m.text.classList.remove("wm-g-dim");
    });
    state.strands.forEach((s) => s.path.classList.remove("wm-g-dim"));
  }
  function attachNodeEvents(nd) {
    const c = nd.circle;
    let dragging = false, moved = false;
    const toLocal = (ev) => {       // pointer position in simulation coordinates
      const box = svg.getBoundingClientRect(), v = state.view;
      return [(((ev.clientX - box.left) / box.width) * W - v.tx) / v.k,
              (((ev.clientY - box.top) / box.height) * H - v.ty) / v.k];
    };
    c.addEventListener("pointerenter", () => { if (state.pinned == null) { hideTip(); showTip(nd); } });
    c.addEventListener("pointerleave", () => { if (state.pinned == null && !dragging) hideTip(); });
    c.addEventListener("pointerdown", (ev) => {
      dragging = state.dragging = true; moved = false;
      c.setPointerCapture(ev.pointerId);
      [nd.fx, nd.fy] = [nd.x, nd.y];
      ev.preventDefault();
    });
    c.addEventListener("pointermove", (ev) => {
      if (!dragging) return;
      const [x, y] = toLocal(ev);
      if (Math.hypot(x - nd.x, y - nd.y) * state.view.k > 3) moved = true;
      [nd.fx, nd.fy] = [x, y];
      if (reduced) { nd.x = x; nd.y = y; draw(); } else state.sim.alphaTarget(0.25).restart();
      if (!tip.hidden) showTip(nd);
    });
    const end = () => {
      if (!dragging) return;
      dragging = state.dragging = false;
      nd.fx = nd.fy = null;
      if (!reduced) state.sim.alphaTarget(0);
      if (!moved && (mode === "ball" || mode === "cycle")) {   // a tap grows a ball from that vertex
        startBloom(performance.now(), nd.i);
        if (!playing) while (!ball.done) stepBloom(0);
        return;
      }
      if (!moved) {            // a tap or click pins the label
        hideTip();
        if (state.pinned === nd.i) state.pinned = null;
        else { state.pinned = nd.i; showTip(nd); }
      }
    };
    c.addEventListener("pointerup", end);
    c.addEventListener("pointercancel", end);
  }
  svg.addEventListener("pointerdown", (ev) => {
    if (ev.target === svg && state.pinned != null) { state.pinned = null; hideTip(); }
  });

  // ---- controls ---------------------------------------------------------------
  buttons.forEach((b) => b.addEventListener("click", () => {
    state.kind = b.dataset.kind; build().then(startLoop);
  }));
  select.addEventListener("change", () => {
    state.p[state.kind] = Number(select.value); build().then(startLoop);
  });
  trailSelect.addEventListener("change", () => {
    TRAIL = Number(trailSelect.value);
    trailLabel.textContent = TRAIL;
    while (trail.length > TRAIL) trail.shift();
    paintTrail();
  });
  viewSelect.addEventListener("change", () => { build().then(startLoop); });
  ellSelect.addEventListener("change", () => spectral.setEll(ellSelect.value));
  labelBox.addEventListener("change", () => {
    gLabels.style.display = labelBox.checked && !labelBox.disabled ? "" : "none";
    draw();
  });
  animSelect.addEventListener("change", () => { mode = animSelect.value; build().then(startLoop); });
  speedSelect.addEventListener("change", () => { hopMs = SPEEDS[speedSelect.value]; });
  if (playBtn) {
    playBtn.addEventListener("click", () => setPlaying(!playing));
    setPlaying(playing);
  }

  // Only animate while the figure is on screen and the tab is visible.
  if ("IntersectionObserver" in window) {
    new IntersectionObserver(([entry]) => {
      visible = entry.isIntersecting;
      if (!state.sim) return;
      if (visible) {
        if (!reduced && mode !== "spectral" && state.sim.alpha() > 0.02) state.sim.restart();
        startLoop();
      }
      else state.sim.stop();
    }, { threshold: 0.05 }).observe(root);
  }
  document.addEventListener("visibilitychange", () => { if (!document.hidden) startLoop(); });

  trailLabel.textContent = TRAIL;
  build().then(() => { startLoop(); root.classList.add("wm-graph-ready"); });
}

// The data sits at assets/data/isogeny/, next to this script's folder, so the same
// widget works on any page.
const root = document.getElementById("wm-graph");
if (root) {
  const base = new URL("../data/isogeny", import.meta.url).href;
  fetch(`${base}/index.json`)
    .then((r) => { if (!r.ok) throw new Error(r.status); return r.json(); })
    .then((index) => init(root, index, base))
    .catch(() => root.classList.add("wm-graph-failed"));
}
