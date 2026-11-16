// Interactive supersingular isogeny graphs for the home page.
// Data: assets/data/isogeny-graphs.json, produced by _scripts/isogeny_graphs.py.
// Colours come from the --wm-* variables, so light and dark mode both work.
import {
  forceSimulation, forceLink, forceManyBody, forceCollide, forceX, forceY,
} from "./vendor/d3-force.min.js";

const W = 600, H = 500, PAD = 40;
const SVG_NS = "http://www.w3.org/2000/svg";
const KINDS = {
  "2": {
    defaultP: 431, maxP: Infinity, ells: [2], bound: "2√2", boundVal: 2 * Math.SQRT2,
    caption: (p, n) => `The supersingular 2-isogeny graph in characteristic ${p}: ` +
      `${n} curves, each with three outgoing 2-isogenies.`,
  },
  "23": {
    defaultP: 211, maxP: Infinity, ells: [2, 3], bound: "2√2 + 2√3",
    boundVal: 2 * Math.SQRT2 + 2 * Math.sqrt(3),
    caption: (p, n) => `The (2,3)-isogeny graph in characteristic ${p}: the same ${n} curves, ` +
      `joined by 2-isogenies and 3-isogenies.`,
  },
  "2b3": {
    defaultP: 157, maxP: 431, ells: [2], bound: "2√2", boundVal: 2 * Math.SQRT2,
    caption: (p, n) => `The 2-isogeny graph with Borel level-3 structure in characteristic ${p}: ` +
      `${n} pairs (E, C), where C is a subgroup of E of order 3.`,
  },
};

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

function init(root, data) {
  const svg = root.querySelector("svg");
  const tip = root.querySelector(".wm-graph-tip");
  const caption = root.querySelector(".wm-graph-caption");
  const stats = root.querySelector(".wm-graph-stats");
  const select = root.querySelector("select");
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
  const TRAIL = 10;
  let trail = [];

  function build() {
    const kind = state.kind, cfg = KINDS[kind];
    const primes = data.primes.filter((q) => q <= cfg.maxP);
    if (!primes.includes(state.p[kind])) {
      state.p[kind] = primes.includes(cfg.defaultP) ? cfg.defaultP : primes[0];
    }
    const p = state.p[kind];
    select.replaceChildren(...primes.map((q) => {
      const o = document.createElement("option");
      o.value = q; o.textContent = q; o.selected = q === p;
      return o;
    }));
    buttons.forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.kind === kind)));

    const entry = data.by_prime[p], g = entry.graphs[kind];
    const n = g.nodes.length;
    const nodes = g.nodes.map((label, i) => ({
      i, label, fp: g.fp[i],
      x: W / 2 + (W / 3) * Math.cos((2 * Math.PI * i) / n),
      y: H / 2 + (H / 3) * Math.sin((2 * Math.PI * i) / n),
    }));

    // One drawn strand per edge of the multigraph; strands between the same
    // pair of vertices fan out as curves.
    const groups = new Map();
    for (const [s, t, ell, m] of g.links) {
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

    const r = n <= 20 ? 8 : n <= 45 ? 6.5 : n <= 110 ? 5 : 4;
    gEdges.replaceChildren();
    gNodes.replaceChildren();
    for (const st of strands) {
      st.path = el("path", { class: `wm-g-edge wm-g-edge-${st.ell}` }, gEdges);
    }
    for (const nd of nodes) {
      nd.circle = el("circle", {
        class: `wm-g-node${nd.fp ? " wm-g-node-fp" : ""}`, r, tabindex: "-1",
      }, gNodes);
      attachNodeEvents(nd);
    }
    Object.assign(state, { nodes, strands, adj, pinned: null, r });
    hideTip();

    if (state.sim) state.sim.stop();
    const dist = Math.max(22, Math.min(90, 520 / Math.sqrt(n)));
    state.sim = forceSimulation(nodes)
      .force("link", forceLink(simLinks).distance(dist).strength(0.5))
      .force("charge", forceManyBody().strength(-Math.max(60, 5200 / n)).distanceMax(320))
      .force("collide", forceCollide(r + 3))
      .force("x", forceX(W / 2).strength(0.06))
      .force("y", forceY(H / 2).strength(0.075))
      .alphaDecay(0.035)
      .on("tick", draw)
      .stop();
    state.sim.tick(reduced ? 300 : 40);
    state.view = null;          // snap the view to the new graph
    draw();
    if (!reduced && visible) state.sim.restart();

    caption.textContent = cfg.caption(p, n);
    const lam = g.lambda2;
    stats.textContent = lam == null ? "" :
      `Second-largest eigenvalue modulus ${lam.toFixed(3)} ≤ ${cfg.bound} ≈ ${cfg.boundVal.toFixed(3)}.`;
    legend3.hidden = !cfg.ells.includes(3);
    if (legendFp) legendFp.hidden = !g.fp.some(Boolean);
    const nr = entry.nonresidue;
    fieldNote.replaceChildren("Labels are j-invariants in F", sub("p"), "(i), i",
      sup("2"), ` = ${nr < 0 ? "−" + (-nr) : nr}.`);
    root.dataset.kind = kind;

    walker = { at: Math.floor(Math.random() * n), next: null, strand: null, t0: 0 };
    trail = [];
    walkerDot.style.display = reduced ? "none" : "";
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
  const HOP_MS = 850;
  function positionWalker(now) {
    if (!walker || !state.nodes.length) return;
    if (reduced) return;
    if (walker.next == null) {
      const options = state.adj[walker.at];
      const pick = options[Math.floor(Math.random() * options.length)];
      walker.next = pick.to; walker.strand = pick.strand; walker.t0 = now;
      walker.reverse = pick.strand.s !== walker.at;
      trail.push(pick.strand);
      if (trail.length > TRAIL) trail.shift();
      paintTrail();
    }
    let t = (now - walker.t0) / HOP_MS;
    if (t >= 1) {
      walker.at = walker.next; walker.next = null;
      t = 1;
    }
    const ease = t * t * (3 - 2 * t);
    const [x, y] = pointOn(geometry(walker.strand), ease, walker.reverse);
    walkerDot.setAttribute("cx", x.toFixed(1));
    walkerDot.setAttribute("cy", y.toFixed(1));
  }
  // Strength 1 for the step being taken now, down to 1/TRAIL for the oldest.
  // A strand walked more than once shows its most recent step.
  function paintTrail() {
    for (const st of state.strands) {
      st.path.classList.remove("wm-g-trail");
      st.path.style.removeProperty("--trail");
    }
    trail.forEach((st, k) => {
      st.path.classList.add("wm-g-trail");
      st.path.style.setProperty("--trail", ((k + 1 + TRAIL - trail.length) / TRAIL).toFixed(2));
    });
  }
  function frame(now) {
    raf = 0;
    if (!visible || document.hidden) return;
    positionWalker(now);
    raf = requestAnimationFrame(frame);
  }
  function startLoop() {
    if (!reduced && !raf) raf = requestAnimationFrame(frame);
  }

  // ---- hover, tap and drag -------------------------------------------------
  function labelOf(nd) {
    return Array.isArray(nd.label)
      ? `(j(E), j(E/C)) = (${nd.label[0]}, ${nd.label[1]})`
      : `j = ${nd.label}`;
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
    state.kind = b.dataset.kind; build(); startLoop();
  }));
  select.addEventListener("change", () => {
    state.p[state.kind] = Number(select.value); build(); startLoop();
  });

  // Only animate while the figure is on screen and the tab is visible.
  if ("IntersectionObserver" in window) {
    new IntersectionObserver(([entry]) => {
      visible = entry.isIntersecting;
      if (visible) { if (!reduced && state.sim.alpha() > 0.02) state.sim.restart(); startLoop(); }
      else state.sim.stop();
    }, { threshold: 0.05 }).observe(root);
  }
  document.addEventListener("visibilitychange", () => { if (!document.hidden) startLoop(); });

  build();
  startLoop();
  root.classList.add("wm-graph-ready");
}

const root = document.getElementById("wm-graph");
if (root) {
  fetch(root.dataset.src)
    .then((r) => { if (!r.ok) throw new Error(r.status); return r.json(); })
    .then((data) => init(root, data))
    .catch(() => root.classList.add("wm-graph-failed"));
}
