#!/usr/bin/env python3
"""Draw the F_p-subgraphs (spines) of the isogeny graphs, for animations/about.qmd.

The spine (F_p-subgraph) of a graph is the subgraph induced on the vertices fixed by
Frobenius (which acts on the level structure too): those vertices and every edge
between them. This script computes each graph with
_scripts/isogeny_graphs.py (all of its consistency checks run), restricts to the
spine, and writes static SVG files to assets/img/spine/.

    python3 _scripts/spine_diagrams.py

Needs networkx (for layouts only): pip install networkx
"""
from __future__ import annotations

import math
import random
import sys
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent))
import isogeny_graphs as ig  # noqa: E402

OUT = Path(__file__).resolve().parent.parent / "assets" / "img" / "spine"
INK, MUTED, BLUE, RED, GOLD, PAPER, LINE = (
    "#1E2430", "#5A6070", "#2B4FB8", "#B93A2B", "#E5A11C", "#FFFDF7", "#E3D9C3")
GREEN, PURPLE = "#2E8B57", "#7A3E9D"
PW, PH = 380, 300          # panel size


# --------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------
_cache: dict = {}


def graphs(p: int, level5: bool) -> dict:
    if (p, level5) not in _cache:
        random.seed(0)
        ig.LEVEL5_MAX_P = max(ig.LEVEL5_MAX_P, p)
        _cache[(p, level5)] = ig.compute(p, level5=level5)["graphs"]
    return _cache[(p, level5)]


def spine(p: int, kind: str):
    """(multigraph on the fixed vertices, set of special vertices) for G_p(kind)."""
    d = graphs(p, kind.endswith("b5"))[kind]
    # The spine: the subgraph induced on the vertices fixed by Frobenius.
    fixed = [i for i, f in enumerate(d["fp"]) if f]
    keep = set(fixed)
    G = nx.MultiGraph()
    G.add_nodes_from(fixed)
    for s, t, ell, m in d["links"]:
        if s in keep and t in keep:
            for _ in range(m):
                G.add_edge(s, t, ell=ell)
    special = {}
    for i in fixed:
        lab = d["nodes"][i]
        j = lab[0] if isinstance(lab, list) else lab
        if j == "0":
            special[i] = "0"
        elif j == str(1728 % p):        # labels store j mod p
            special[i] = "1728"
    return G, special


def pick(kind: str, cond, lo=7, hi=46, start=300):
    """Smallest prime >= start with p satisfying cond and a spine of lo..hi vertices."""
    p = start
    while True:
        p += 1
        if not ig.is_prime(p) or not cond(p):
            continue
        n = spine(p, kind)[0].number_of_nodes()
        if lo <= n <= hi:
            return p


# --------------------------------------------------------------------------
# layout and drawing
# --------------------------------------------------------------------------
def layout_components(G):
    """Place each connected component, then shelf-pack them into a PW x PH panel."""
    simple = nx.Graph(G)
    simple.remove_edges_from(nx.selfloop_edges(simple))
    comps = sorted(nx.connected_components(simple), key=lambda c: (-len(c), min(c)))
    boxes = []
    for c in comps:
        S = simple.subgraph(c)
        n = len(c)
        if n == 1:
            pos, w, h = {next(iter(c)): (0.0, 0.0)}, 0.0, 0.0
        elif n == 2:
            a, b = sorted(c)
            pos, w, h = {a: (0.0, 0.0), b: (1.0, 0.0)}, 64.0, 0.0
        else:
            raw = nx.kamada_kawai_layout(S, scale=1.0)
            xs = [v[0] for v in raw.values()]
            ys = [v[1] for v in raw.values()]
            span = max(max(xs) - min(xs), max(ys) - min(ys)) or 1
            pos = {v: ((raw[v][0] - min(xs)) / span, (raw[v][1] - min(ys)) / span) for v in c}
            w = h = 26 + 30 * math.sqrt(n)
            w *= (max(xs) - min(xs)) / span or 1
            h *= (max(ys) - min(ys)) / span or 1
            pos = {v: ((raw[v][0] - min(xs)) / span * (26 + 30 * math.sqrt(n)) / w if w else 0,
                       (raw[v][1] - min(ys)) / span * (26 + 30 * math.sqrt(n)) / h if h else 0)
                   for v in c}
        boxes.append((c, pos, w, h))
    # shelf packing into rows of roughly square overall shape
    pad, gap, margin = 15, 24, 22
    width = PW - 2 * margin - 2 * pad
    area = sum((w + gap) * (h + gap) for _, _, w, h in boxes)
    row_w = min(width, max(max(w for _, _, w, _ in boxes), math.sqrt(area * 1.6)))
    placed, rects, x, y, rowh = {}, [], 0.0, 0.0, 0.0
    for c, pos, w, h in boxes:
        if x + w > row_w and x > 0:
            x, y, rowh = 0.0, y + rowh + gap, 0.0
        for v in c:
            placed[v] = (x + pos[v][0] * w, y + pos[v][1] * h)
        rects.append((x, y, w, h))
        x += w + gap
        rowh = max(rowh, h)
    total_h = y + rowh
    used_w = max(rx + w for rx, _, w, _ in rects)
    avail_h = PH - 58 - margin - 2 * pad
    k = min(avail_h / total_h if total_h else 9, width / used_w if used_w else 9, 2.3)
    ox = margin + pad + (width - used_w * k) / 2
    oy = 58 + pad + (avail_h - total_h * k) / 2
    boxes_out = [(ox + rx * k, oy + ry * k, w * k, h * k) for rx, ry, w, h in rects]
    return ({v: (ox + px * k, oy + py * k) for v, (px, py) in placed.items()},
            len(comps), k, boxes_out)


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;")


def panel(p: int, kind: str, note: str) -> str:
    G, special = spine(p, kind)
    n = G.number_of_nodes()
    out = [f'<rect x="0.75" y="0.75" width="{PW - 1.5}" height="{PH - 1.5}" rx="16" '
           f'fill="{PAPER}" stroke="{LINE}" stroke-width="1.5"/>']
    pos, ncomp, k, boxes = layout_components(G) if n else ({}, 0, 1, [])
    comp_word = "component" if ncomp == 1 else "components"
    out.append(f'<text x="18" y="26" font-size="14" font-weight="700" fill="{INK}">p = {p}'
               f'<tspan font-weight="400" fill="{MUTED}">  ·  {esc(note)}</tspan></text>')
    out.append(f'<text x="18" y="44" font-size="12" fill="{MUTED}">'
               f'{n} vertices, {G.number_of_edges()} edges, {ncomp} {comp_word}</text>')
    if not n:
        out.append(f'<text x="{PW / 2}" y="{PH / 2 + 16}" font-size="14" text-anchor="middle" '
                   f'fill="{MUTED}">empty</text>')
        return "\n".join(out)
    r = max(3.4, min(8.5, 6.2 * math.sqrt(k) + 1.4 - 0.06 * n))
    if ncomp > 1:       # a faint card behind each component
        pad = 15
        for bx, by, bw, bh in boxes:
            out.append(f'<rect x="{bx - pad:.1f}" y="{by - pad:.1f}" width="{bw + 2 * pad:.1f}" '
                       f'height="{bh + 2 * pad:.1f}" rx="12" fill="#F6EFDD" stroke="none"/>')
    seen = {}
    for s, t, data in G.edges(data=True):
        key = (min(s, t), max(s, t))
        seen.setdefault(key, []).append(data["ell"])
    for (s, t), ells in sorted(seen.items()):
        for idx, ell in enumerate(sorted(ells)):
            col = BLUE if ell == 2 else RED if ell == 3 else GREEN
            (x1, y1), (x2, y2) = pos[s], pos[t]
            if s == t:
                ux, uy = x1 - PW / 2, y1 - PH / 2
                ln = math.hypot(ux, uy) or 1
                ux, uy = ux / ln, uy / ln
                size = 17 + 7 * idx
                c1 = (x1 + size * (ux * 0.8 - uy * 0.6), y1 + size * (uy * 0.8 + ux * 0.6))
                c2 = (x1 + size * (ux * 0.8 + uy * 0.6), y1 + size * (uy * 0.8 - ux * 0.6))
                out.append(f'<path d="M{x1:.1f},{y1:.1f}C{c1[0]:.1f},{c1[1]:.1f} {c2[0]:.1f},'
                           f'{c2[1]:.1f} {x1:.1f},{y1:.1f}" fill="none" stroke="{col}" '
                           f'stroke-width="1.6" opacity="0.85"/>')
                continue
            off = (idx - (len(ells) - 1) / 2) * 10
            dx, dy = x2 - x1, y2 - y1
            ln = math.hypot(dx, dy) or 1
            cx, cy = (x1 + x2) / 2 - dy / ln * off, (y1 + y2) / 2 + dx / ln * off
            out.append(f'<path d="M{x1:.1f},{y1:.1f}Q{cx:.1f},{cy:.1f} {x2:.1f},{y2:.1f}" '
                       f'fill="none" stroke="{col}" stroke-width="1.8" opacity="0.85"/>')
    for v, (x, y) in sorted(pos.items()):
        tag = special.get(v)
        if tag == "1728":
            out.append(f'<rect x="{x - r:.1f}" y="{y - r:.1f}" width="{2 * r:.1f}" height="{2 * r:.1f}" '
                       f'fill="{PURPLE}" stroke="{PAPER}" stroke-width="1.5"/>')
        else:
            fill = INK if tag == "0" else GOLD
            out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" fill="{fill}" '
                       f'stroke="{PAPER}" stroke-width="1.5"/>')
    return "\n".join(out)


def legend(ells, has_special: bool, y: int) -> str:
    items = [(BLUE, "line", "2-isogeny")]
    if 3 in ells:
        items.append((RED, "line", "3-isogeny"))
    if 5 in ells:
        items.append((GREEN, "line", "5-isogeny"))
    items.append((GOLD, "dot", "fixed by Frobenius"))
    if has_special:
        items += [(INK, "dot", "j = 0"), (PURPLE, "sq", "j = 1728")]
    out, x = [], 8
    for col, shape, text in items:
        if shape == "line":
            out.append(f'<line x1="{x}" y1="{y}" x2="{x + 18}" y2="{y}" stroke="{col}" stroke-width="2.4"/>')
            x += 24
        elif shape == "dot":
            out.append(f'<circle cx="{x + 5}" cy="{y}" r="5" fill="{col}"/>')
            x += 16
        else:
            out.append(f'<rect x="{x}" y="{y - 5}" width="10" height="10" fill="{col}"/>')
            x += 16
        out.append(f'<text x="{x}" y="{y + 4}" font-size="12" fill="{MUTED}">{esc(text)}</text>')
        x += 7 * len(text) + 18
    return "\n".join(out)


def figure(name: str, kind: str, cases, cols=2, ells=(2,), has_special=True):
    """cases: list of (prime, note). Writes assets/img/spine/<name>.svg."""
    rows = math.ceil(len(cases) / cols)
    gap = 16
    W = cols * PW + (cols - 1) * gap
    H = rows * PH + (rows - 1) * gap + 34
    body = []
    for k, (p, note) in enumerate(cases):
        x, y = (k % cols) * (PW + gap), (k // cols) * (PH + gap)
        body.append(f'<g transform="translate({x},{y})">{panel(p, kind, note)}</g>')
    body.append(f'<g transform="translate(0,{H - 14})">{legend(ells, has_special, 0)}</g>')
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" '
           f'font-family="Inter, system-ui, -apple-system, Segoe UI, sans-serif" role="img">\n'
           + "\n".join(body) + "\n</svg>\n")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.svg").write_text(svg)
    print(f"{name}.svg: " + ", ".join(f"p={p}" for p, _ in cases))


def main():
    for old in OUT.glob("spine-*.svg"):
        old.unlink()
    # the 2-isogeny graph, by p mod 8
    cases = [(pick("2", lambda p, r=r: p % 8 == r, lo=6, hi=30), f"p ≡ {r} (mod 8)")
             for r in (1, 5, 3, 7)]
    figure("spine-2", "2", cases)
    classes = (1, 13, 5, 17, 7, 19, 11, 23)
    # the (2,3)- and (2,3,5)-isogeny graphs, by p mod 24
    for kind, ells in (("23", (2, 3)), ("235", (2, 3, 5))):
        cases = [(pick(kind, lambda p, r=r: p % 24 == r, lo=6, hi=30), f"p ≡ {r} (mod 24)")
                 for r in classes]
        figure(f"spine-{kind}", kind, cases, ells=ells)
    # level 5 (six times as many vertices, so smaller primes)
    for kind, ells in (("2b5", (2,)), ("23b5", (2, 3))):
        cases = [(pick(kind, lambda p, r=r: p % 24 == r, lo=8, hi=60, start=60), f"p ≡ {r} (mod 24)")
                 for r in classes]
        figure(f"spine-{kind}", kind, cases, ells=ells)


if __name__ == "__main__":
    main()
