"""Generate the overview illustration as smooth vector SVG, in the manner of a journal's
in-house scientific art: gradient shading with one light source, soft inner shadows at the
edges instead of ink outlines, saturated colour only on the elements in focus (the
leiomyoma, the edited cells), neutral stone tones for context, pale blue glass.

    python3 illustration/vector.py

Writes the <svg> into docs/index.html between <!-- design:start --> and <!-- design:end -->.
Geometry is in the figure's viewBox units (900 x 360).
"""
import math
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PAGE = HERE.parent / "docs" / "index.html"

# palette: stone tones for tissue (context), the site's vermilion for the focus, pale blue glass
STONE_HI, STONE, STONE_LO, STONE_LINE = "#ece7df", "#d2cabd", "#b3a897", "#9e9382"
PINK_HI, PINK, PINK_LO = "#f7e9e4", "#edd3cb", "#d4aba3"
VERM_HI, VERM, VERM_LO, VERM_LINE = "#f3a872", "#d55e00", "#a84600", "#7e3400"
GLASS, GLASS_LINE, MEDIUM_HI, MEDIUM = "#eef3f7", "#b5c4cf", "#e3eef5", "#c7dbe8"
INK, ARROW = "#1c1e26", "#4a4f5a"


# ----------------------------------------------------------------------------- geometry
def flatten(d, n=28):
    """An SVG path (absolute M, L, C, Q, Z) as a list of point lists."""
    tokens = re.findall(r"[MLCQZ]|-?\d*\.?\d+", d)
    polys, cur, pos, i = [], [], None, 0
    while i < len(tokens):
        t = tokens[i]
        if t == "M":
            if cur: polys.append(cur)
            pos = (float(tokens[i + 1]), float(tokens[i + 2])); cur = [pos]; i += 3
        elif t == "L":
            pos = (float(tokens[i + 1]), float(tokens[i + 2])); cur.append(pos); i += 3
        elif t == "C":
            p0 = np.array(pos); p1, p2, p3 = (np.array([float(tokens[i + k]), float(tokens[i + k + 1])]) for k in (1, 3, 5))
            for s in np.linspace(0, 1, n)[1:]:
                cur.append(tuple((1 - s) ** 3 * p0 + 3 * (1 - s) ** 2 * s * p1 + 3 * (1 - s) * s ** 2 * p2 + s ** 3 * p3))
            pos = tuple(p3); i += 7
        elif t == "Q":
            p0 = np.array(pos); p1, p2 = (np.array([float(tokens[i + k]), float(tokens[i + k + 1])]) for k in (1, 3))
            for s in np.linspace(0, 1, n)[1:]:
                cur.append(tuple((1 - s) ** 2 * p0 + 2 * (1 - s) * s * p1 + s ** 2 * p2))
            pos = tuple(p2); i += 5
        elif t == "Z":
            if cur: polys.append(cur); cur = []
            i += 1
        else:
            i += 1
    if cur: polys.append(cur)
    return polys


def stroke_poly(points, widths):
    """The outline of a tapered stroke along points, as a closed point list."""
    p = np.asarray(points, float)
    t = np.gradient(p, axis=0); t /= np.maximum(np.hypot(*t.T), 1e-9)[:, None]
    nrm = np.column_stack([-t[:, 1], t[:, 0]])
    w = np.asarray(widths, float)[:, None] / 2
    return [tuple(q) for q in np.vstack([p + nrm * w, (p - nrm * w)[::-1]])]


def spindle(cx, cy, length, width, angle, n=64):
    """A smooth muscle cell: a lens tapering to both ends, as a point list."""
    pts = []
    for a in np.linspace(0, 2 * math.pi, n, endpoint=False):
        x = math.cos(a) * length / 2
        y = math.sin(a) * width / 2 * (1 - 0.5 * abs(math.cos(a)) ** 2)
        c, s = math.cos(angle), math.sin(angle)
        pts.append((cx + x * c - y * s, cy + x * s + y * c))
    return pts


def smooth_path(points, closed=True):
    """A closed point list as an SVG path of Catmull-Rom-derived cubic Béziers."""
    p = np.asarray(points, float)
    n = len(p)
    out = [f"M{p[0][0]:.1f} {p[0][1]:.1f}"]
    rng_ = range(n) if closed else range(n - 1)
    for i in rng_:
        p0, p1, p2, p3 = p[(i - 1) % n], p[i], p[(i + 1) % n], p[(i + 2) % n]
        c1, c2 = p1 + (p2 - p0) / 6, p2 - (p3 - p1) / 6
        out.append(f"C{c1[0]:.1f} {c1[1]:.1f} {c2[0]:.1f} {c2[1]:.1f} {p2[0]:.1f} {p2[1]:.1f}")
    return " ".join(out) + (" Z" if closed else "")


def spiral(cx, cy, r0, r1, turns, n=120):
    pts = []
    for k in range(n):
        s = k / (n - 1)
        a, r = s * turns * 2 * math.pi, r0 + (r1 - r0) * s
        pts.append((cx + r * math.cos(a), cy + r * math.sin(a)))
    return smooth_path(pts, closed=False)


# ----------------------------------------------------------------------------- drawing
def defs():
    return f"""
  <defs>
    <marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">
      <path d="M0 0.8 L10 5 L0 9.2 Z" fill="{ARROW}"/>
    </marker>
    <radialGradient id="gStone" cx="38%" cy="28%" r="78%">
      <stop offset="0" stop-color="{STONE_HI}"/><stop offset="0.55" stop-color="{STONE}"/><stop offset="1" stop-color="{STONE_LO}"/>
    </radialGradient>
    <linearGradient id="gTube" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0" stop-color="{STONE_HI}"/><stop offset="0.6" stop-color="{STONE}"/><stop offset="1" stop-color="{STONE_LO}"/>
    </linearGradient>
    <radialGradient id="gOvary" cx="35%" cy="30%" r="80%">
      <stop offset="0" stop-color="{STONE_HI}"/><stop offset="0.7" stop-color="{STONE}"/><stop offset="1" stop-color="{STONE_LO}"/>
    </radialGradient>
    <linearGradient id="gCavity" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0" stop-color="{PINK_HI}"/><stop offset="1" stop-color="{PINK}"/>
    </linearGradient>
    <radialGradient id="gFibroid" cx="36%" cy="30%" r="75%">
      <stop offset="0" stop-color="{VERM_HI}"/><stop offset="0.5" stop-color="{VERM}"/><stop offset="1" stop-color="{VERM_LO}"/>
    </radialGradient>
    <radialGradient id="gSphere" cx="36%" cy="32%" r="72%">
      <stop offset="0" stop-color="{VERM_HI}"/><stop offset="0.55" stop-color="{VERM}"/><stop offset="1" stop-color="{VERM_LO}"/>
    </radialGradient>
    <linearGradient id="gGlass" x1="0" y1="0" x2="1" y2="0">
      <stop offset="0" stop-color="#ffffff" stop-opacity="0.9"/><stop offset="0.5" stop-color="{GLASS}"/><stop offset="1" stop-color="{GLASS}"/>
    </linearGradient>
    <linearGradient id="gMedium" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0" stop-color="{MEDIUM_HI}"/><stop offset="1" stop-color="{MEDIUM}"/>
    </linearGradient>
    <radialGradient id="gSheen" cx="50%" cy="50%" r="50%">
      <stop offset="0" stop-color="#ffffff" stop-opacity="0.55"/><stop offset="1" stop-color="#ffffff" stop-opacity="0"/>
    </radialGradient>
    <!-- a soft shadow just inside an edge: the shading that stands in for an outline -->
    <filter id="inTissue" x="-20%" y="-20%" width="140%" height="140%">
      <feComponentTransfer in="SourceAlpha" result="inv"><feFuncA type="table" tableValues="1 0"/></feComponentTransfer>
      <feGaussianBlur in="inv" stdDeviation="5" result="blur"/>
      <feOffset in="blur" dx="-2" dy="-3" result="off"/>
      <feComposite in="off" in2="SourceAlpha" operator="in" result="shade"/>
      <feFlood flood-color="#6e6355" flood-opacity="0.45" result="col"/>
      <feComposite in="col" in2="shade" operator="in" result="tinted"/>
      <feMerge><feMergeNode in="SourceGraphic"/><feMergeNode in="tinted"/></feMerge>
    </filter>
    <filter id="inCavity" x="-20%" y="-20%" width="140%" height="140%">
      <feComponentTransfer in="SourceAlpha" result="inv"><feFuncA type="table" tableValues="1 0"/></feComponentTransfer>
      <feGaussianBlur in="inv" stdDeviation="4" result="blur"/>
      <feOffset in="blur" dx="1" dy="2" result="off"/>
      <feComposite in="off" in2="SourceAlpha" operator="in" result="shade"/>
      <feFlood flood-color="#a8766e" flood-opacity="0.55" result="col"/>
      <feComposite in="col" in2="shade" operator="in" result="tinted"/>
      <feMerge><feMergeNode in="SourceGraphic"/><feMergeNode in="tinted"/></feMerge>
    </filter>
    <filter id="inGlass" x="-20%" y="-20%" width="140%" height="140%">
      <feComponentTransfer in="SourceAlpha" result="inv"><feFuncA type="table" tableValues="1 0"/></feComponentTransfer>
      <feGaussianBlur in="inv" stdDeviation="3" result="blur"/>
      <feComposite in="blur" in2="SourceAlpha" operator="in" result="shade"/>
      <feFlood flood-color="#7f93a3" flood-opacity="0.35" result="col"/>
      <feComposite in="col" in2="shade" operator="in" result="tinted"/>
      <feMerge><feMergeNode in="SourceGraphic"/><feMergeNode in="tinted"/></feMerge>
    </filter>
    <filter id="cast" x="-30%" y="-30%" width="160%" height="160%">
      <feDropShadow dx="2" dy="3" stdDeviation="2.5" flood-color="#4a3326" flood-opacity="0.3"/>
    </filter>
    <filter id="castSoft" x="-30%" y="-30%" width="160%" height="160%">
      <feDropShadow dx="1" dy="2" stdDeviation="1.6" flood-color="#4a3326" flood-opacity="0.25"/>
    </filter>
    <filter id="blur6" x="-50%" y="-50%" width="200%" height="200%"><feGaussianBlur stdDeviation="6"/></filter>
  </defs>"""


def uterus():
    body = ("M220 72 C 262 72, 300 90, 306 120 C 312 150, 296 196, 266 232 C 258 242, 254 250, 254 258 "
            "C 256 288, 250 304, 236 309 C 228 312, 212 312, 204 309 C 190 304, 184 288, 186 258 "
            "C 186 250, 182 242, 174 232 C 144 196, 128 150, 134 120 C 140 90, 178 72, 220 72 Z")
    cavity = ("M178 124 C 196 116, 244 116, 262 124 C 250 158, 232 200, 224 234 Q 220 242, 216 234 "
              "C 208 200, 190 158, 178 124 Z")
    parts = []
    # tubes and ovaries behind the body
    for sign in (1, -1):
        cx = lambda x: 220 + sign * (x - 220)
        centre = np.array([(cx(x), y) for x, y in flatten("M300 106 C 326 92, 348 92, 364 106 C 376 116, 382 130, 394 140")[0]])
        parts.append(f'<path d="{smooth_path(stroke_poly(centre, np.linspace(8, 12, len(centre))))}" fill="url(#gTube)" stroke="{STONE_LINE}" stroke-width="0.8" filter="url(#inTissue)"/>')
        end, heading = centre[-1], centre[-1] - centre[-5]
        base = math.atan2(heading[1], heading[0])
        fingers = []
        for k, spread in enumerate((-1.0, -0.5, 0.0, 0.5, 1.0)):
            a, L = base + spread, (13, 17, 19, 17, 13)[k]
            tip = end + np.array([math.cos(a), math.sin(a)]) * L
            mid = (end + tip) / 2 + np.array([-math.sin(a), math.cos(a)]) * 1.8 * (1 if spread >= 0 else -1)
            fingers.append(f'<path d="{smooth_path(stroke_poly(np.vstack([end, mid, tip]), [5.5, 3.6, 1.6]))}" fill="{STONE_HI}" stroke="{STONE_LINE}" stroke-width="0.6"/>')
        parts.append(f'<g>{"".join(fingers)}</g>')
        ox, oy, rot = cx(406), 165, -16 * sign
        parts.append(f'<g transform="rotate({rot} {ox} {oy})" filter="url(#inTissue)"><ellipse cx="{ox}" cy="{oy}" rx="11" ry="16" fill="url(#gOvary)" stroke="{STONE_LINE}" stroke-width="0.8"/></g>')
        for fx, fy, fr in ((-3.5, -6, 2.4), (3.5, -1, 2.0), (-1.5, 6, 2.6)):
            parts.append(f'<circle cx="{ox + fx}" cy="{oy + fy}" r="{fr}" fill="{STONE_HI}" stroke="{STONE_LO}" stroke-width="0.6"/>')
    # the body of the uterus, with its endometrial cavity and cervical canal
    parts.append(f'<path d="{body}" fill="url(#gStone)" stroke="{STONE_LINE}" stroke-width="0.9" filter="url(#inTissue)"/>')
    parts.append(f'<ellipse cx="196" cy="104" rx="44" ry="18" fill="url(#gSheen)" transform="rotate(-18 196 104)"/>')
    parts.append(f'<path d="{cavity}" fill="url(#gCavity)" stroke="{PINK_LO}" stroke-width="0.7" filter="url(#inCavity)"/>')
    parts.append(f'<rect x="217.6" y="242" width="4.8" height="60" rx="2.4" fill="{PINK}" stroke="{PINK_LO}" stroke-width="0.6"/>')
    # the leiomyoma: the focus, saturated, seated in the wall with a cast shadow and a whorled surface
    fx, fy, fr = 168, 190, 20
    parts.append(f'<g filter="url(#cast)"><circle cx="{fx}" cy="{fy}" r="{fr}" fill="url(#gFibroid)" stroke="{VERM_LO}" stroke-width="0.8"/></g>')
    parts.append(f'<path d="{spiral(fx - 1, fy + 1, 3, 15, 2.0)}" fill="none" stroke="{VERM_LINE}" stroke-opacity="0.22" stroke-width="1.1" stroke-linecap="round"/>')
    parts.append(f'<ellipse cx="{fx - 7}" cy="{fy - 8}" rx="6" ry="4" fill="#ffffff" fill-opacity="0.35" transform="rotate(-30 {fx - 7} {fy - 8})"/>')
    return "\n      ".join(parts)


def cell(cx, cy, L, Wd, ang, gid, small=False):
    """A spindle-shaped smooth muscle cell with a nucleus, lit from the upper left."""
    body = smooth_path(spindle(cx, cy, L, Wd, ang))
    nuc = smooth_path(spindle(cx, cy, L * 0.36, Wd * 0.6, ang))
    deg = math.degrees(ang)
    grad = (f'<linearGradient id="{gid}" gradientUnits="userSpaceOnUse" x1="{cx}" y1="{cy - Wd}" x2="{cx}" y2="{cy + Wd}" gradientTransform="rotate({deg:.1f} {cx} {cy})">'
            f'<stop offset="0" stop-color="{VERM_HI}"/><stop offset="0.5" stop-color="{VERM}"/><stop offset="1" stop-color="{VERM_LO}"/></linearGradient>')
    sw = 0.5 if small else 0.7
    shape = (f'<path d="{body}" fill="url(#{gid})" stroke="{VERM_LO}" stroke-width="{sw}"/>'
             f'<path d="{nuc}" fill="{VERM_LINE}" fill-opacity="0.85"/>')
    if not small:
        hx, hy = cx - L * 0.12 * math.cos(ang) - Wd * 0.18 * math.sin(ang), cy - L * 0.12 * math.sin(ang) + Wd * 0.18 * math.cos(ang)
        shape += f'<ellipse cx="{hx:.1f}" cy="{hy:.1f}" rx="{L * 0.14:.1f}" ry="{Wd * 0.14:.1f}" fill="#ffffff" fill-opacity="0.45" transform="rotate({deg:.1f} {hx:.1f} {hy:.1f})"/>'
    return grad, shape


def cultures():
    rng = np.random.default_rng(5)
    grads, parts = [], []
    # the derived, edited cells: a small cluster
    for i, (cx, cy, L, Wd, ang) in enumerate(((548, 204, 46, 13, -0.38), (570, 225, 40, 12, 0.55), (541, 233, 38, 11.5, 1.3))):
        g, sh = cell(cx, cy, L, Wd, ang, f"gCell{i}")
        grads.append(g); parts.append(f'<g filter="url(#castSoft)">{sh}</g>')
    # 2D culture: a dish in section with medium and a single layer of cells on its floor
    dish = "M650 258 L650 290 Q650 304 664 304 L788 304 Q802 304 802 290 L802 258 L796 258 L796 288 Q796 298 786 298 L666 298 Q656 298 656 288 L656 258 Z"
    parts.append(f'<path d="M656 282 L796 282 L796 288 Q796 298 786 298 L666 298 Q656 298 656 288 Z" fill="url(#gMedium)"/>')
    layer = []
    x = 664
    k = 0
    while x < 786:
        L = rng.uniform(15, 21)
        g, sh = cell(x + L / 2, 295.3 + rng.uniform(-0.4, 0.4), L, rng.uniform(4.2, 5.0), rng.uniform(-0.16, 0.16), f"gLay{k}", small=True)
        grads.append(g); layer.append(sh); x += L * rng.uniform(0.9, 1.02); k += 1
    parts.append("".join(layer))
    parts.append(f'<path d="{dish}" fill="url(#gGlass)" stroke="{GLASS_LINE}" stroke-width="0.8" filter="url(#inGlass)"/>')
    parts.append(f'<path d="M659 262 L659 286" stroke="#ffffff" stroke-opacity="0.9" stroke-width="1.6" stroke-linecap="round"/>')
    parts.append(f'<path d="M657 282 L795 282" stroke="#ffffff" stroke-opacity="0.7" stroke-width="1"/>')
    # 3D culture: a round-bottom well with medium and spheroids
    well = ("M658 56 L658 108 Q658 152 702 152 L758 152 Q802 152 802 108 L802 56 L796 56 L796 108 Q796 146 758 146 L702 146 "
            "Q664 146 664 108 L664 56 Z")
    parts.append(f'<path d="M664 92 L796 92 L796 108 Q796 146 758 146 L702 146 Q664 146 664 108 Z" fill="url(#gMedium)"/>')
    for i, (sx, sy, sr) in enumerate(((701, 121, 12.5), (748, 126, 11.5), (727, 96, 10.5))):
        parts.append(f'<ellipse cx="{sx + 1}" cy="{sy + sr + 2}" rx="{sr * 0.9}" ry="{sr * 0.22}" fill="#5f7384" fill-opacity="0.18" filter="url(#blur6)"/>')
        parts.append(f'<circle cx="{sx}" cy="{sy}" r="{sr}" fill="url(#gSphere)" stroke="{VERM_LO}" stroke-width="0.7"/>')
        # the aggregate's cells, hinted as a few darker arcs and a soft highlight
        for a, rr in ((0.9, 0.55), (2.4, 0.62), (4.0, 0.5), (5.3, 0.6)):
            px, py = sx + sr * rr * math.cos(a), sy + sr * rr * math.sin(a)
            parts.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{sr * 0.3:.1f}" fill="none" stroke="{VERM_LINE}" stroke-opacity="0.28" stroke-width="0.8"/>')
        parts.append(f'<ellipse cx="{sx - sr * 0.35}" cy="{sy - sr * 0.4}" rx="{sr * 0.32}" ry="{sr * 0.22}" fill="#ffffff" fill-opacity="0.45" transform="rotate(-35 {sx - sr * 0.35} {sy - sr * 0.4})"/>')
    parts.append(f'<path d="{well}" fill="url(#gGlass)" stroke="{GLASS_LINE}" stroke-width="0.8" filter="url(#inGlass)"/>')
    parts.append(f'<path d="M667 60 L667 104" stroke="#ffffff" stroke-opacity="0.9" stroke-width="1.6" stroke-linecap="round"/>')
    parts.append(f'<path d="M665 92 L795 92" stroke="#ffffff" stroke-opacity="0.7" stroke-width="1"/>')
    return "\n    ".join(grads), "\n      ".join(parts)


def labels():
    return f"""
      <g stroke="{ARROW}" stroke-width="1.8" fill="none" stroke-linecap="round">
        <path d="M282 238 C 340 262, 430 250, 512 222" marker-end="url(#arrow)"/>
        <path d="M580 238 C 602 264, 622 282, 642 290" marker-end="url(#arrow)"/>
        <path d="M576 200 C 600 170, 620 138, 644 112" marker-end="url(#arrow)"/>
      </g>
      <g stroke="{ARROW}" stroke-width="1.2" fill="none" stroke-linecap="round">
        <path d="M154 204 L 118 246"/>
        <path d="M186 268 L 138 306"/>
        <path d="M236 192 L 326 200"/>
      </g>
      <g fill="{INK}"><circle cx="154" cy="204" r="2.4"/><circle cx="186" cy="268" r="2.4"/><circle cx="236" cy="192" r="2.4"/></g>
      <g font-size="17" fill="{INK}">
        <text x="112" y="251" text-anchor="end">Leiomyoma</text>
        <text x="134" y="311" text-anchor="end">Myometrium</text>
        <text x="332" y="205">Endometrium</text>
        <text x="570" y="270" text-anchor="end">MED12-edited cells</text>
        <text x="726" y="334" text-anchor="middle">2D culture</text>
        <text x="730" y="44" text-anchor="middle">3D spheroids</text>
      </g>
      <g font-size="14" fill="#6b7080">
        <text x="570" y="292" text-anchor="end">derived from the myometrium</text>
      </g>"""


def main():
    cell_grads, culture_parts = cultures()
    svg = f"""<svg class="design-svg" viewBox="0 0 900 360" role="img" aria-labelledby="design-title" xmlns="http://www.w3.org/2000/svg">
      <title id="design-title">Illustration: a frontal section of the uterus with a leiomyoma in the myometrial wall; myometrial cells edited at MED12 codon 44, grown as a flat 2D culture and as 3D spheroids.</title>{defs()}
  <defs>
    {cell_grads}
  </defs>
      {uterus()}
      {culture_parts}{labels()}
    </svg>"""
    html = PAGE.read_text()
    start, end = html.index("<!-- design:start -->"), html.index("<!-- design:end -->")
    html = html[:start] + "<!-- design:start -->\n    " + svg + "\n    " + html[end:]
    PAGE.write_text(html)
    print(f"wrote the illustration into {PAGE.name}: {len(svg) // 1024} KB of SVG")


if __name__ == "__main__":
    main()
