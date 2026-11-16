#!/usr/bin/env python3
"""Figures for animations/balls.qmd: balls grown until they contain a Frobenius
conjugate pair.

Reads the graph data in assets/data/isogeny/ (made by isogeny_graphs.py) and writes

  assets/img/ball/ball-examples.svg   four worked examples, drawn layer by layer
  animations/_ball-chart.qmd          a chart (inline SVG, so it has hover labels and
                                      follows the page) and a table of the same numbers

    python3 _scripts/ball_diagrams.py

The ball: breadth-first search from a vertex v, one layer at a time, stopping at the
first radius r where B_r(v) contains some w and its Frobenius conjugate sigma(w) != w.
Pure Python, no dependencies.
"""
from __future__ import annotations

import json
import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "assets" / "data" / "isogeny"
INK, MUTED, BLUE, RED, GOLD, PAPER, LINE, FAINT = (
    "#1E2430", "#5A6070", "#2B4FB8", "#B93A2B", "#E5A11C", "#FFFDF7", "#E3D9C3", "#F3EBD8")
NAMES = {"2": "G_p(2)", "23": "G_p({2,3})", "2b3": "G_p(2), level 3", "23b5": "G_p({2,3}), level 5"}
LONG = {"2": "2-isogeny graph", "23": "(2,3)-isogeny graph",
        "2b3": "2-isogeny graph, level 3", "23b5": "(2,3)-isogeny graph, level 5"}


def load(p, kind):
    g = json.loads((DATA / f"p{p}.json").read_text())["graphs"].get(kind)
    if g is None:
        return None
    n = len(g["nodes"])
    adj = [dict() for _ in range(n)]          # neighbour -> list of ells
    for s, t, ell, m in g["links"]:
        if s == t:
            continue
        adj[s].setdefault(t, []).append(ell)
        adj[t].setdefault(s, []).append(ell)
    return g, adj


def grow(adj, sigma, v):
    """Layers, parents, and the pair found (w, sigma(w)), or None."""
    dist, parent, layers = {v: 0}, {}, [[v]]
    while True:
        nxt = []
        for u in layers[-1]:
            for w in sorted(adj[u]):
                if w not in dist:
                    dist[w] = len(layers)
                    parent[w] = u
                    nxt.append(w)
        if not nxt:
            return layers, parent, dist, None
        layers.append(nxt)
        best = None
        for w in nxt:
            c = sigma[w]
            if c != w and c in dist:
                score = dist[w] + dist[c]
                if best is None or score < best[0]:
                    best = (score, w, c)
        if best:
            return layers, parent, dist, best[1:]


def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;")


# --------------------------------------------------------------------------
# worked examples
# --------------------------------------------------------------------------
PW, PH = 400, 430


def example_panel(p, kind, seed):
    g, adj = load(p, kind)
    sigma, n = g["sigma"], len(g["nodes"])
    rng = random.Random(seed)
    runs = {v: grow(adj, sigma, v) for v in range(n)}
    good = [v for v, r in runs.items() if r[3] and len(r[0]) - 1 >= 2 and len(r[2]) <= 45]
    if not good:
        good = [v for v, r in runs.items() if r[3]]
    v = rng.choice(sorted(good))
    layers, parent, dist, (w, c) = runs[v]
    R = len(layers) - 1

    # angular sectors from the BFS tree, proportional to the number of leaves below
    kids = {u: [] for u in dist}
    for x, par in parent.items():
        kids[par].append(x)
    leaves = {}

    def count(u):
        leaves[u] = sum(count(k) for k in kids[u]) or 1
        return leaves[u]

    count(v)
    ang = {v: 0.0}

    def place(u, a0, a1):
        ang[u] = (a0 + a1) / 2
        a = a0
        for k in sorted(kids[u], key=lambda k: -leaves[k]):
            span = (a1 - a0) * leaves[k] / leaves[u]
            place(k, a, a + span)
            a += span

    place(v, -math.pi / 2, 3 * math.pi / 2)
    cx, cy, rmax = PW / 2, 78 + (PH - 120) / 2, min(PW, PH - 120) / 2 - 22
    step = rmax / max(R, 1)
    pos = {u: (cx + dist[u] * step * math.cos(ang[u]), cy + dist[u] * step * math.sin(ang[u]))
           for u in dist}

    geo = set()
    for x in (w, c):
        while x in parent:
            geo.add(frozenset((x, parent[x])))
            x = parent[x]

    lab = g["nodes"]
    short = (lambda i: lab[i][0]) if isinstance(lab[0], list) else (lambda i: lab[i])
    out = [f'<rect x="0.75" y="0.75" width="{PW - 1.5}" height="{PH - 1.5}" rx="16" '
           f'fill="{PAPER}" stroke="{LINE}" stroke-width="1.5"/>',
           f'<text x="18" y="27" font-size="14" font-weight="700" fill="{INK}">{esc(LONG[kind])}, '
           f'p = {p}</text>',
           f'<text x="18" y="46" font-size="12" fill="{MUTED}">{n} vertices; √n ≈ {math.sqrt(n):.1f}</text>',
           f'<text x="18" y="63" font-size="12" fill="{MUTED}">pair found at radius {R}, '
           f'ball of {len(dist)} vertices</text>']
    for k in range(1, R + 1):
        out.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{k * step:.1f}" fill="none" '
                   f'stroke="{LINE}" stroke-width="1.2" stroke-dasharray="3 4"/>')
        out.append(f'<text x="{cx + 4:.1f}" y="{cy - k * step - 4:.1f}" font-size="10.5" '
                   f'fill="{MUTED}">r = {k}</text>')
    drawn = set()
    for u in dist:
        for x, ells in adj[u].items():
            if x not in dist or frozenset((u, x)) in drawn:
                continue
            drawn.add(frozenset((u, x)))
            if frozenset((u, x)) in geo:
                continue
            col = BLUE if 2 in ells else RED
            (x1, y1), (x2, y2) = pos[u], pos[x]
            out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
                       f'stroke="{col}" stroke-width="1.3" opacity="0.55"/>')
    for e in geo:
        a, b = tuple(e)
        (x1, y1), (x2, y2) = pos[a], pos[b]
        out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
                   f'stroke="{INK}" stroke-width="3.4" stroke-linecap="round"/>')
    for u, (x, y) in pos.items():
        fill = RED if u in (w, c) else GOLD if g["jfp"][u] else BLUE
        big = u in (w, c) or u == v
        r = 7.5 if big else 4.6
        ring = INK if big else PAPER
        out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}" stroke="{ring}" '
                   f'stroke-width="{2.5 if big else 1.3}"><title>{esc(short(u))} (distance {dist[u]})</title></circle>')
    for u, name in ((w, "j"), (c, "j^p")):
        x, y = pos[u]
        dx = 12 if x >= cx else -12
        anchor = "start" if dx > 0 else "end"
        out.append(f'<text x="{x + dx:.1f}" y="{y + 4:.1f}" font-size="11.5" text-anchor="{anchor}" '
                   f'fill="{INK}" stroke="{PAPER}" stroke-width="3" paint-order="stroke">'
                   f'{esc(short(u))}</text>')
    out.append(f'<text x="{cx:.1f}" y="{cy + 22:.1f}" font-size="11" text-anchor="middle" fill="{INK}" '
               f'stroke="{PAPER}" stroke-width="3" paint-order="stroke">start</text>')
    return "\n".join(out)


def examples():
    cases = [("2", 1009, 3), ("23", 431, 1), ("2b3", 431, 2), ("23b5", 211, 5)]
    cols, gap = 2, 16
    W = cols * PW + gap
    H = 2 * PH + gap + 36
    body = []
    for k, (kind, p, seed) in enumerate(cases):
        x, y = (k % cols) * (PW + gap), (k // cols) * (PH + gap)
        body.append(f'<g transform="translate({x},{y})">{example_panel(p, kind, seed)}</g>')
    y = H - 14
    items = [("line", BLUE, "2-isogeny"), ("line", RED, "3-isogeny"), ("dot", BLUE, "curve"),
             ("dot", GOLD, "j-invariant in F_p"), ("dot", RED, "the conjugate pair"),
             ("thick", INK, "shortest paths")]
    x, leg = 6, []
    for shape, col, text in items:
        if shape == "dot":
            leg.append(f'<circle cx="{x + 5}" cy="{y}" r="5" fill="{col}"/>')
            x += 15
        else:
            leg.append(f'<line x1="{x}" y1="{y}" x2="{x + 18}" y2="{y}" stroke="{col}" '
                       f'stroke-width="{3.4 if shape == "thick" else 2}"/>')
            x += 24
        leg.append(f'<text x="{x}" y="{y + 4}" font-size="12" fill="{MUTED}">{esc(text)}</text>')
        x += 7 * len(text) + 16
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" '
           f'font-family="Inter, system-ui, -apple-system, Segoe UI, sans-serif" role="img">\n'
           + "\n".join(body) + "\n" + "\n".join(leg) + "\n</svg>\n")
    out = ROOT / "assets" / "img" / "ball"
    out.mkdir(parents=True, exist_ok=True)
    (out / "ball-examples.svg").write_text(svg)
    print("ball-examples.svg:", ", ".join(f"{k} p={p}" for k, p, _ in cases))


# --------------------------------------------------------------------------
# the chart: mean ball size over all starting vertices, against n
# --------------------------------------------------------------------------
SERIES = [  # kind, colour, marker (secondary encoding: shape as well as colour)
    ("2", "#2B4FB8", "circle"), ("23", "#B93A2B", "square"),
    ("2b3", "#C88A0E", "triangle"), ("23b5", "#2E8B57", "diamond")]


def stats():
    primes = json.loads((DATA / "index.json").read_text())["primes"]
    rows = {}
    for kind, *_ in SERIES:
        rows[kind] = []
        for p in primes:
            got = load(p, kind)
            if not got:
                continue
            g, adj = got
            n = len(g["nodes"])
            sizes, radii = [], []
            for v in range(n):
                layers, _, dist, pair = grow(adj, g["sigma"], v)
                if pair:
                    sizes.append(len(dist))
                    radii.append(len(layers) - 1)
            rows[kind].append((p, n, sum(sizes) / len(sizes), sum(radii) / len(radii)))
    return rows


def marker(shape, x, y, col):
    s = 4.6
    if shape == "circle":
        return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{s}" fill="{col}" stroke="{PAPER}" stroke-width="2"/>'
    if shape == "square":
        return (f'<rect x="{x - s:.1f}" y="{y - s:.1f}" width="{2 * s:.1f}" height="{2 * s:.1f}" rx="1.5" '
                f'fill="{col}" stroke="{PAPER}" stroke-width="2"/>')
    if shape == "triangle":
        t = s * 1.25
        return (f'<path d="M{x:.1f},{y - t:.1f}L{x + t:.1f},{y + t * 0.8:.1f}L{x - t:.1f},{y + t * 0.8:.1f}Z" '
                f'fill="{col}" stroke="{PAPER}" stroke-width="2"/>')
    t = s * 1.3
    return (f'<path d="M{x:.1f},{y - t:.1f}L{x + t:.1f},{y:.1f}L{x:.1f},{y + t:.1f}L{x - t:.1f},{y:.1f}Z" '
            f'fill="{col}" stroke="{PAPER}" stroke-width="2"/>')


def chart(rows):
    W, H = 760, 440
    L, R, T, B = 64, 150, 58, 52
    x0, x1, y0, y1 = math.log10(6), math.log10(3000), math.log10(2.5), math.log10(150)
    X = lambda n: L + (math.log10(n) - x0) / (x1 - x0) * (W - L - R)
    Y = lambda v: H - B - (math.log10(v) - y0) / (y1 - y0) * (H - T - B)
    out = [f'<rect x="0.75" y="0.75" width="{W - 1.5}" height="{H - 1.5}" rx="16" fill="{PAPER}" '
           f'stroke="{LINE}" stroke-width="1.5"/>',
           f'<text x="20" y="28" font-size="15" font-weight="700" fill="{INK}">Size of the ball when it first '
           f'contains a conjugate pair</text>',
           f'<text x="20" y="46" font-size="12" fill="{MUTED}">mean over every starting vertex, '
           f'against the number of vertices n (both axes logarithmic)</text>']
    for v in (3, 10, 30, 100):
        out.append(f'<line x1="{L}" x2="{W - R}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" stroke="{FAINT}" stroke-width="1"/>')
        out.append(f'<text x="{L - 8}" y="{Y(v) + 4:.1f}" font-size="11" text-anchor="end" fill="{MUTED}">{v}</text>')
    for n in (10, 30, 100, 300, 1000, 3000):
        if math.log10(n) > x1 + 1e-9:
            continue
        out.append(f'<line x1="{X(n):.1f}" x2="{X(n):.1f}" y1="{T}" y2="{H - B}" stroke="{FAINT}" stroke-width="1"/>')
        out.append(f'<text x="{X(n):.1f}" y="{H - B + 17}" font-size="11" text-anchor="middle" fill="{MUTED}">{n}</text>')
    out.append(f'<line x1="{L}" x2="{W - R}" y1="{H - B}" y2="{H - B}" stroke="{LINE}" stroke-width="1.2"/>')
    out.append(f'<text x="{(L + W - R) / 2:.1f}" y="{H - 12}" font-size="12" text-anchor="middle" '
               f'fill="{MUTED}">n, the number of vertices</text>')
    out.append(f'<text x="18" y="{(T + H - B) / 2:.1f}" font-size="12" text-anchor="middle" fill="{MUTED}" '
               f'transform="rotate(-90 18 {(T + H - B) / 2:.1f})">vertices in the ball</text>')
    # reference lines c * sqrt(n)
    ends = []
    for c, name in ((1, "√n"), (2, "2√n")):
        na, nb = 6, 3000
        out.append(f'<line x1="{X(na):.1f}" y1="{Y(c * math.sqrt(na)):.1f}" x2="{X(nb):.1f}" '
                   f'y2="{Y(c * math.sqrt(nb)):.1f}" stroke="{MUTED}" stroke-width="1.2" stroke-dasharray="5 5"/>')
        ends.append((Y(c * math.sqrt(nb)), name, MUTED, None))
    for kind, col, shape in SERIES:
        pts = sorted((n, m, p, rad) for p, n, m, rad in rows[kind])
        d = "M" + " L".join(f"{X(n):.1f},{Y(m):.1f}" for n, m, *_ in pts)
        out.append(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="2" stroke-linejoin="round"/>')
        for n, m, p, rad in pts:
            out.append(f'<g>{marker(shape, X(n), Y(m), col)}<circle cx="{X(n):.1f}" cy="{Y(m):.1f}" r="11" '
                       f'fill="transparent"><title>{esc(LONG[kind])}, p = {p}: n = {n}, mean ball '
                       f'{m:.1f} vertices (√n = {math.sqrt(n):.1f}), mean radius {rad:.2f}</title></circle></g>')
        n, m, *_ = pts[-1]
        ends.append((Y(m), NAMES[kind], INK, (shape, col)))
    # direct labels at the right, pushed apart so they do not overlap
    ends.sort()
    ys = [e[0] for e in ends]
    for k in range(1, len(ys)):
        ys[k] = max(ys[k], ys[k - 1] + 15)
    for (y, name, col, mk), yy in zip(ends, ys):
        x = W - R + 12
        if mk:
            out.append(marker(mk[0], x + 4, yy - 4, mk[1]))
            out.append(f'<text x="{x + 14}" y="{yy:.1f}" font-size="11.5" fill="{col}">{esc(name)}</text>')
        else:
            out.append(f'<text x="{x}" y="{yy:.1f}" font-size="11.5" fill="{col}">{esc(name)}</text>')
    # legend row, under the subtitle (also says which shape is which)
    x = L
    leg = []
    for kind, col, shape in SERIES:
        leg.append(marker(shape, x + 5, T - 2, col))
        leg.append(f'<text x="{x + 15}" y="{T + 2}" font-size="11.5" fill="{MUTED}">{esc(NAMES[kind])}</text>')
        x += 15 + 7 * len(NAMES[kind]) + 20
    out += leg
    svg = (f'<svg class="wm-ball-chart" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
           f'font-family="Inter, system-ui, -apple-system, Segoe UI, sans-serif" role="img" '
           f'aria-label="Mean ball size against n on log-log axes; every graph follows a line of slope one half, '
           f'between about one and two and a half times the square root of n.">\n' + "\n".join(out) + "\n</svg>")
    table = ["| graph | p | n | mean radius | mean ball size | ball / √n |", "|---|--:|--:|--:|--:|--:|"]
    for kind, *_ in SERIES:
        for p, n, m, rad in rows[kind]:
            table.append(f"| {LONG[kind]} | {p} | {n} | {rad:.2f} | {m:.1f} | {m / math.sqrt(n):.2f} |")
    qmd = ("<!-- Made by _scripts/ball_diagrams.py; do not edit by hand. -->\n\n"
           "```{=html}\n<div class=\"wm-figure-scroll\">\n" + svg + "\n</div>\n```\n\n"
           "<details>\n<summary>The numbers in the chart</summary>\n\n" + "\n".join(table) +
           "\n\n</details>\n")
    (ROOT / "animations" / "_ball-chart.qmd").write_text(qmd)
    print("_ball-chart.qmd written")


if __name__ == "__main__":
    examples()
    chart(stats())
