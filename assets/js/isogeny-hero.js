// Interactive supersingular isogeny graphs for the home page.
// Data: assets/data/isogeny/ (index.json plus one p<prime>.json per prime, loaded when
// chosen), produced by _scripts/isogeny_graphs.py.
// Colours come from the --wm-* variables, so light and dark mode both work.
import {
  forceSimulation, forceLink, forceManyBody, forceCollide, forceX, forceY,
} from "./vendor/d3-force.min.js";

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
  "2b3": {
    defaultP: 157, levelCap: true, ells: [2], bound: "2√2", boundVal: 2 * Math.SQRT2,
    caption: (p, n) => `The 2-isogeny graph with Borel level-3 structure in characteristic ${p}: ` +
      `${n} pairs (E, C), where C is a subgroup of E of order 3.`,
  },
  "23b5": {
    defaultP: 101, cap5: true, ells: [2, 3], bound: "2√2 + 2√3",
    boundVal: 2 * Math.SQRT2 + 2 * Math.sqrt(3),
    caption: (p, n) => `The (2,3)-isogeny graph with Borel level-5 structure in characteristic ${p}: ` +
      `${n} pairs (E, C), where C is a subgroup of E of order 5.`,
  },
};

// The three ways to look at each graph.
//  full    every vertex.
//  folded  the quotient by Frobenius (E, C) -> (E^(p), C^(p)): conjugate pairs become one vertex.
//  spine   only the vertices fixed by Frobenius (the F_p-rational ones), and the edges among them.
function derive(g, view) {
  const n = g.nodes.length, sigma = g.sigma;
  const base = (i) => ({ label: g.nodes[i], members: [i] });
  if (view === "spine") {
    const keep = [], at = new Map();
    for (let i = 0; i < n; i++) if (g.fp[i]) { at.set(i, keep.length); keep.push(i); }
    return {
      nodes: keep.map((i) => ({ ...base(i), fp: true })),
      links: g.links.filter(([s, t]) => at.has(s) && at.has(t))
        .map(([s, t, ell, m]) => [at.get(s), at.get(t), ell, m]),
    };
  }
  if (view === "folded") {
    const at = new Map(), nodes = [];
    for (let i = 0; i < n; i++) {
      if (sigma[i] < i) continue;                    // its partner comes first
      const pair = sigma[i] !== i;
      at.set(i, nodes.length);
      if (pair) at.set(sigma[i], nodes.length);
      nodes.push({ label: g.nodes[i], members: pair ? [i, sigma[i]] : [i], fp: !pair, pair });
    }
    // An edge of the quotient is a Frobenius-orbit of edges of the original graph.
    const seen = new Set(), links = [];
    const norm = (a, b) => (a <= b ? [a, b] : [b, a]);
    for (const [s, t, ell, m] of g.links) {
      const [a, b] = norm(s, t), [c, d] = norm(sigma[s], sigma[t]);
      const key = `${ell}|${a < c || (a === c && b <= d) ? `${a}-${b}` : `${c}-${d}`}`;
      if (seen.has(key)) continue;
      seen.add(key);
      links.push([...norm(at.get(s), at.get(t)), ell, m]);
    }
    return { nodes, links };
  }
  return { nodes: g.nodes.map((_, i) => ({ ...base(i), fp: g.fp[i] })), links: g.links };
}

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

// How long one step of the random walk takes, in milliseconds.
const SPEEDS = { calm: 700, fast: 220, frantic: 90 };
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
  const legendPair = root.querySelector(".wm-graph-legend-pair");
  const trailLabel = root.querySelector(".wm-trail-n");
  const buttons = [...root.querySelectorAll("[data-kind]")];
  const legend3 = root.querySelector(".wm-graph-legend-3");
  const legendFp = root.querySelector(".wm-graph-legend-fp");
  const fieldNote = root.querySelector(".wm-graph-field");
  const reduced = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  const gEdges = el("g", { class: "wm-g-edges" }, svg);
  const gNodes = el("g", { class: "wm-g-nodes" }, svg);
  const walkerDot = el("circle", { class: "wm-g-walker", r: 5 }, svg);

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
  function load(p) {
    if (!cache.has(p)) {
      cache.set(p, fetch(`${base}/p${p}.json`)
        .then((r) => { if (!r.ok) throw new Error(r.status); return r.json(); })
        .catch((e) => { cache.delete(p); throw e; }));
    }
    return cache.get(p);
  }

  let buildId = 0;
  async function build() {
    const id = ++buildId;
    const kind = state.kind, cfg = KINDS[kind];
    const cap = cfg.levelCap ? index.level_max_p : cfg.cap5 ? index.level5_max_p : Infinity;
    const primes = index.primes.filter((q) => q <= cap);
    if (!primes.includes(state.p[kind])) {
      state.p[kind] = primes.includes(cfg.defaultP) ? cfg.defaultP : primes[0];
    }
    // The growing ball needs the whole graph; the view menu is set aside meanwhile.
    const view = mode === "ball" ? "full" : viewSelect.value;
    viewSelect.disabled = mode === "ball";
    // Some spines are empty (no vertex fixed by Frobenius): for example the level-3 graph
    // when p = 1 (mod 3). In the spine view those primes are marked in the menu, and
    // choosing one shows a short explanation with buttons to move on (see showEmpty).
    const empty = (q) => view === "spine" && index.spine && index.spine[q] && index.spine[q][kind] === 0;
    const p = state.p[kind];
    select.replaceChildren(...primes.map((q) => {
      const o = document.createElement("option");
      o.value = q; o.selected = q === p;
      o.textContent = empty(q) ? `${q} (empty spine)` : q;
      return o;
    }));
    const nearestNonempty = () => {
      const ok = primes.filter((q) => !empty(q));
      return ok.length ? ok.reduce((a, b) => (Math.abs(b - p) < Math.abs(a - p) ? b : a)) : null;
    };
    buttons.forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.kind === kind)));

    root.classList.add("wm-graph-loading");
    let entry;
    try { entry = await load(p); } catch {
      if (id === buildId) {
        root.classList.remove("wm-graph-loading");
        caption.textContent = "Could not load this graph. Please try again.";
      }
      return;
    }
    if (id !== buildId) return;        // a newer choice has been made meanwhile
    root.classList.remove("wm-graph-loading");
    const g = entry.graphs[kind];
    const d = derive(g, view);
    const n = d.nodes.length;
    const nodes = d.nodes.map((nd, i) => ({
      ...nd, i,
      x: W / 2 + (W / 3) * Math.cos((2 * Math.PI * i) / n),
      y: H / 2 + (H / 3) * Math.sin((2 * Math.PI * i) / n),
    }));

    // One drawn strand per edge of the multigraph; strands between the same
    // pair of vertices fan out as curves.
    const groups = new Map();
    for (const [s, t, ell, m] of d.links) {
      const key = `${s}-${t}`;
      if (!groups.has(key)) groups.set(key, []);
      for (let k = 0; k < m; k++) groups.get(key).push(ell);
    }
    const strands = [], simLinks = [], adj = nodes.map(() => []);
    for (const [key, ells] of groups) {
      const [s, t] = key.split("-").map(Number);
      if (s !== t) simLinks.push({ source: s, target: t });
      ells.forEach((ell, k) => {
        const strand = { s, t, ell, k, total: ells.length };
        strands.push(strand);
        adj[s].push({ to: t, strand });
        if (s !== t) adj[t].push({ to: s, strand });
      });
    }

    const r = n <= 20 ? 8 : n <= 45 ? 6.5 : n <= 110 ? 5 : n <= 300 ? 4 : n <= 700 ? 3 : 2.2;
    const rOf = (nd) => (nd.pair ? r * 1.35 : r);
    gEdges.replaceChildren();
    gNodes.replaceChildren();
    for (const st of strands) {
      st.path = el("path", { class: `wm-g-edge wm-g-edge-${st.ell}` }, gEdges);
    }
    for (const nd of nodes) {
      nd.circle = el("circle", {
        class: `wm-g-node${nd.fp ? " wm-g-node-fp" : ""}${nd.pair ? " wm-g-node-pair" : ""}`,
        r: rOf(nd), tabindex: "-1",
      }, gNodes);
      attachNodeEvents(nd);
    }
    Object.assign(state, { nodes, strands, adj, pinned: null, r, graph: g, sigma: g.sigma });
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

    const nFixed = g.fp.filter(Boolean).length, nFull = g.nodes.length;
    const nComp = components(n, d.links);
    if (view === "folded") {
      caption.textContent = cfg.caption(p, nFull).replace(/ in characteristic/, " folded by Frobenius, in characteristic")
        .replace(/: .*$/, `: ${n} vertices, the orbits of Frobenius on ${nFull} vertices.`);
      stats.textContent = `${nFixed} fixed by Frobenius (the spine) and ${(nFull - nFixed) / 2} conjugate pairs. ` +
        `Edges are Frobenius-orbits of edges.`;
    } else if (view === "spine") {
      caption.textContent = cfg.caption(p, nFull).replace(/ in characteristic/, " restricted to its spine, in characteristic")
        .replace(/: .*$/, `: the ${n} vertices defined over F_p, of ${nFull}, and the edges between them.`);
      stats.textContent = n === 0 ? "The spine is empty." :
        nComp === 1 ? "The spine is connected." :
        nComp === n ? `The spine is totally disconnected: ${n} vertices and no edges between them.` :
        `The spine has ${nComp} connected components.`;
    } else {
      caption.textContent = cfg.caption(p, n);
      const lam = g.lambda2;
      stats.textContent = (lam == null ? "" :
        `Second-largest eigenvalue modulus ${lam.toFixed(3)} ≤ ${cfg.bound} ≈ ${cfg.boundVal.toFixed(3)}. `) +
        `${nFixed} of the ${n} vertices are defined over F_p (the spine).`;
    }
    legend3.hidden = !cfg.ells.includes(3);
    if (legendFp) legendFp.hidden = !d.nodes.some((nd) => nd.fp);
    if (legendPair) legendPair.hidden = view !== "folded";
    const nEdges = d.links.filter(([a, b]) => a !== b).length;
    showEmpty(view === "spine" && (n === 0 || nEdges === 0)
      ? { kind, p, n, alternative: nearestNonempty() } : null);
    const nr = entry.nonresidue;
    fieldNote.replaceChildren("Labels are j-invariants in F", sub("p"), "(i), i",
      sup("2"), ` = ${nr < 0 ? "−" + (-nr) : nr}.`);
    root.dataset.kind = kind;

    const movable = nodes.filter((nd) => adj[nd.i].length);
    root.dataset.mode = mode;
    svg.classList.toggle("wm-ball-mode", mode === "ball");
    ball = null;
    ballStatus.hidden = mode !== "ball";
    if (mode === "ball") {
      trail = []; paintTrail();
      if (reduced) {              // no animation: show one finished ball
        startBall(0);
        while (!ball.done) stepBall(0);
      } else ballStatus.textContent = "";
    }
    walker = { at: movable.length ? movable[Math.floor(Math.random() * movable.length)].i : 0,
               next: null, strand: null, t0: 0 };
    trail = [];
    walkerDot.style.display = reduced || !strands.length || mode === "ball" ? "none" : "";
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
    const a = state.nodes[st.s], b = state.nodes[st.t];
    if (st.s === st.t) {
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
    for (const st of state.strands) {
      const g = geometry(st);
      st.path.setAttribute("d", g.loop
        ? `M${f(g.a.sx)},${f(g.a.sy)}C${f(g.c1[0])},${f(g.c1[1])} ${f(g.c2[0])},${f(g.c2[1])} ${f(g.a.sx)},${f(g.a.sy)}`
        : `M${f(g.a.sx)},${f(g.a.sy)}Q${f(g.c[0])},${f(g.c[1])} ${f(g.b.sx)},${f(g.b.sy)}`);
    }
  }

  // ---- random walk -------------------------------------------------------
  function positionWalker(now) {
    if (!walker || !state.nodes.length) return;
    if (reduced) return;
    if (walker.next == null) {
      const options = state.adj[walker.at];
      if (!options.length) return;          // a vertex with no edges (a tiny spine)
      const pick = options[Math.floor(Math.random() * options.length)];
      walker.next = pick.to; walker.strand = pick.strand; walker.t0 = now;
      walker.reverse = pick.strand.s !== walker.at;
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
    const [x, y] = pointOn(geometry(walker.strand), ease, walker.reverse);
    walkerDot.setAttribute("cx", x.toFixed(1));
    walkerDot.setAttribute("cy", y.toFixed(1));
  }
  // Strength 1 for the step being taken now, fading to 0.15 for the oldest in the window.
  // A strand walked more than once shows its most recent step.
  function paintTrail() {
    for (const st of painted) {
      st.path.classList.remove("wm-g-trail");
      st.path.style.removeProperty("--trail");
    }
    trail.forEach((st, k) => {
      const lin = (k + 1 + TRAIL - trail.length) / TRAIL;
      st.path.classList.add("wm-g-trail");
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
  // Why a spine can be empty: a supersingular curve over F_p (p > 3) has trace 0, so
  // its Frobenius pi satisfies pi^2 = -p. A subgroup of order N (N = 3 or 5, prime)
  // that Frobenius preserves is an eigenline of pi on E[N], whose eigenvalue would be a
  // square root of -p mod N. There is none for p = 1 (mod 3), or p = 2, 3 (mod 5).
  const LEVEL = { "2b3": 3, "23b5": 5 };
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
    const N = LEVEL[kind];
    const head = document.createElement("strong");
    const text = document.createElement("div");
    text.className = "wm-graph-empty-text";
    if (n === 0) {
      head.textContent = "Empty spine";
      text.append(`No vertex is fixed by Frobenius for p = ${p}, so the spine is empty, not just disconnected.`);
      if (N) {
        // The reason, folded away on small screens so the buttons stay in view.
        const why = document.createElement("details");
        why.open = window.matchMedia("(min-width: 576px)").matches;
        const summary = document.createElement("summary");
        summary.textContent = "Why?";
        why.append(summary, "On a curve over F", sub("p"), ` Frobenius satisfies π² = −p, and −p has no ` +
          `square root mod ${N} when p ≡ ${p % N} (mod ${N}), so no subgroup of order ${N} is ` +
          "stable under Frobenius.");
        text.append(why);
      }
    } else {
      head.textContent = "Totally disconnected spine";
      text.append(`For p = ${p} the spine has ${n} vertices but no edges between them.`);
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
  // Breadth-first search from a random centre v, one layer at a time. Frobenius is
  // an automorphism of the graph, so the ball B_r(v) contains a pair {w, w^(p)}
  // exactly when B_r(v) and B_r(v^(p)) share a vertex w not fixed by Frobenius.
  const shortLabel = (i) => {
    const l = state.graph.nodes[i];
    return Array.isArray(l) ? `(${l[0]}, ${l[1]})` : `j = ${l}`;
  };
  function startBall(now, at) {
    const n = state.nodes.length;
    const center = at ?? Math.floor(Math.random() * n);
    ball = { center, dist: new Map([[center, 0]]), parent: new Map(), layers: [[center]],
             r: 0, t0: now, done: null, holdUntil: 0 };
    paintBall();
  }
  function stepBall(now) {
    const b = ball, sigma = state.sigma, next = [];
    for (const u of b.layers[b.r]) {
      for (const { to, strand } of state.adj[u]) {
        if (b.dist.has(to)) continue;
        b.dist.set(to, b.r + 1);
        b.parent.set(to, { from: u, strand });
        next.push(to);
      }
    }
    if (!next.length) {
      b.done = { none: true };
    } else {
      b.layers.push(next);
      b.r += 1;
      // Any new pair has a member in the new layer. Prefer the one closest to the centre.
      let best = null;
      for (const w of next) {
        const c = sigma[w];
        if (c === w || !b.dist.has(c)) continue;
        const score = b.dist.get(w) + b.dist.get(c);
        if (!best || score < best.score) best = { w, c, score };
      }
      b.done = best;
    }
    if (b.done) b.holdUntil = now + BALL_HOLD;
    paintBall();
  }
  function growBall(now) {
    if (!state.nodes.length) return;
    if (!ball) { startBall(now); return; }
    if (ball.done) { if (now >= ball.holdUntil) startBall(now); return; }
    if (now - ball.t0 >= BALL_MS[speedSelect.value]) { ball.t0 = now; stepBall(now); }
  }
  function paintBall() {
    const b = ball, n = state.nodes.length;
    const hit = b.done && !b.done.none ? [b.done.w, b.done.c] : [];
    const geo = new Set();               // the strands of the two shortest paths
    for (const x of hit) {
      for (let y = x; b.parent.has(y); y = b.parent.get(y).from) geo.add(b.parent.get(y).strand);
    }
    for (const nd of state.nodes) {
      const d = b.dist.get(nd.i), c = nd.circle;
      c.classList.toggle("wm-g-out", d === undefined);
      c.classList.toggle("wm-g-front", d !== undefined && d === b.r && b.r > 0 && !b.done);
      c.classList.toggle("wm-g-center", nd.i === b.center);
      c.classList.toggle("wm-g-hit", hit.includes(nd.i));
      const r0 = nd.pair ? state.r * 1.35 : state.r;
      c.setAttribute("r", hit.includes(nd.i) || nd.i === b.center ? r0 * 1.6 : r0);
    }
    for (const st of state.strands) {
      st.path.classList.toggle("wm-g-inball", b.dist.has(st.s) && b.dist.has(st.t));
      st.path.classList.toggle("wm-g-geo", geo.has(st));
    }
    const size = b.dist.size, root2 = Math.sqrt(n).toFixed(1);
    if (b.done && b.done.none) {
      ballStatus.textContent = `The ball around ${shortLabel(b.center)} filled the graph without ` +
        "meeting a conjugate pair (every vertex it reached is defined over F_p).";
    } else if (b.done) {
      ballStatus.textContent = `Radius ${b.r}: the ball around ${shortLabel(b.center)} ` +
        `(${size} of ${n} vertices; √n ≈ ${root2}) now contains the conjugate pair ` +
        `${shortLabel(b.done.w)} and ${shortLabel(b.done.c)}. The pair is in red; bold: shortest paths from the centre to each.`;
    } else {
      ballStatus.textContent = `Growing a ball around ${shortLabel(b.center)}: radius ${b.r}, ` +
        `${size} of ${n} vertices.`;
    }
  }

  function frame(now) {
    raf = 0;
    if (!visible || document.hidden) return;
    if (mode === "ball") growBall(now); else positionWalker(now);
    raf = requestAnimationFrame(frame);
  }
  function startLoop() {
    if (!reduced && !raf) raf = requestAnimationFrame(frame);
  }

  // ---- hover, tap and drag -------------------------------------------------
  function labelOf(nd) {
    const one = (label) => (Array.isArray(label)
      ? `(j(E), j(E/C)) = (${label[0]}, ${label[1]})` : `j = ${label}`);
    if (!nd.pair) return one(nd.label);
    const other = state.graph.nodes[nd.members[1]];
    return `{E, E⁽ᵖ⁾}: ${one(nd.label)}  ~  ${one(other)}`;
  }
  function showTip(nd) {
    const nbrs = new Set(state.adj[nd.i].map((e) => e.to));
    nbrs.add(nd.i);
    const hot = new Set(state.adj[nd.i].map((e) => e.strand));
    state.nodes.forEach((m) => m.circle.classList.toggle("wm-g-dim", !nbrs.has(m.i)));
    state.strands.forEach((s) => s.path.classList.toggle("wm-g-dim", !hot.has(s)));
    nd.circle.classList.add("wm-g-focus");
    tip.replaceChildren(labelOf(nd), ...(nd.fp ? [" · defined over F", sub("p")] : []));
    tip.hidden = false;
    const box = svg.getBoundingClientRect(), host = root.getBoundingClientRect();
    const px = box.left - host.left + (nd.sx / W) * box.width;
    const py = box.top - host.top + (nd.sy / H) * box.height;
    tip.style.left = `${Math.max(8, Math.min(host.width - tip.offsetWidth - 8, px - tip.offsetWidth / 2))}px`;
    tip.style.top = `${Math.max(0, py - tip.offsetHeight - 14)}px`;
  }
  function hideTip() {
    tip.hidden = true;
    state.nodes.forEach((m) => m.circle.classList.remove("wm-g-dim", "wm-g-focus"));
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
      if (!moved && mode === "ball") {   // in ball mode a tap grows a ball from that vertex
        startBall(performance.now(), nd.i);
        if (reduced) while (!ball.done) stepBall(0);
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
  animSelect.addEventListener("change", () => { mode = animSelect.value; build().then(startLoop); });
  speedSelect.addEventListener("change", () => { hopMs = SPEEDS[speedSelect.value]; });

  // Only animate while the figure is on screen and the tab is visible.
  if ("IntersectionObserver" in window) {
    new IntersectionObserver(([entry]) => {
      visible = entry.isIntersecting;
      if (!state.sim) return;
      if (visible) { if (!reduced && state.sim.alpha() > 0.02) state.sim.restart(); startLoop(); }
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
