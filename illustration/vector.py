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

# palette: stone tones for tissue (context), the site's orange for the focus, pale blue glass
STONE_HI, STONE, STONE_LO, STONE_LINE = "#ece7df", "#d2cabd", "#b3a897", "#9e9382"
PINK_HI, PINK, PINK_LO = "#f7e9e4", "#edd3cb", "#d4aba3"
ORANGE_HI, ORANGE, ORANGE_LO, ORANGE_LINE = "#f3a872", "#d55e00", "#a84600", "#7e3400"
PINK_C_HI, PINK_C, PINK_C_LO, PINK_C_LINE = "#e6b3cf", "#bf6a9b", "#8e4a73", "#5e2f4c"
BLUE_HI, BLUE, BLUE_LO, BLUE_LINE = "#7cc0ea", "#0072b2", "#004b78", "#003757"
KO_HI, KO, KO_LO, KO_LINE = "#7b7f8a", "#3a3d46", "#1c1e26", "#0b0c10"
CELL = {"stone": (STONE_HI, STONE, STONE_LO, STONE_LINE), "verm": (ORANGE_HI, ORANGE, ORANGE_LO, ORANGE_LINE),
        "pink": (PINK_C_HI, PINK_C, PINK_C_LO, PINK_C_LINE), "ko": (KO_HI, KO, KO_LO, KO_LINE)}
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
      <stop offset="0" stop-color="{ORANGE_HI}"/><stop offset="0.5" stop-color="{ORANGE}"/><stop offset="1" stop-color="{ORANGE_LO}"/>
    </radialGradient>
    <radialGradient id="gFibroidBlue" cx="36%" cy="30%" r="75%">
      <stop offset="0" stop-color="{BLUE_HI}"/><stop offset="0.5" stop-color="{BLUE}"/><stop offset="1" stop-color="{BLUE_LO}"/>
    </radialGradient>
    <radialGradient id="gSphereVerm" cx="36%" cy="32%" r="72%">
      <stop offset="0" stop-color="{ORANGE_HI}"/><stop offset="0.55" stop-color="{ORANGE}"/><stop offset="1" stop-color="{ORANGE_LO}"/>
    </radialGradient>
    <radialGradient id="gSpherePink" cx="36%" cy="32%" r="72%">
      <stop offset="0" stop-color="{PINK_C_HI}"/><stop offset="0.55" stop-color="{PINK_C}"/><stop offset="1" stop-color="{PINK_C_LO}"/>
    </radialGradient>
    <radialGradient id="gSphereStone" cx="36%" cy="32%" r="72%">
      <stop offset="0" stop-color="{STONE_HI}"/><stop offset="0.55" stop-color="{STONE}"/><stop offset="1" stop-color="{STONE_LO}"/>
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
    parts.append(f'<g filter="url(#cast)"><circle cx="{fx}" cy="{fy}" r="{fr}" fill="url(#gFibroid)" stroke="{ORANGE_LO}" stroke-width="0.8"/></g>')
    parts.append(f'<path d="{spiral(fx - 1, fy + 1, 3, 15, 2.0)}" fill="none" stroke="{ORANGE_LINE}" stroke-opacity="0.22" stroke-width="1.1" stroke-linecap="round"/>')
    parts.append(f'<ellipse cx="{fx - 7}" cy="{fy - 8}" rx="6" ry="4" fill="#ffffff" fill-opacity="0.35" transform="rotate(-30 {fx - 7} {fy - 8})"/>')
    # a second leiomyoma, MED12 wild-type, on the fundus
    bx, by, br = 206, 80, 13
    parts.append(f'<g filter="url(#cast)"><circle cx="{bx}" cy="{by}" r="{br}" fill="url(#gFibroidBlue)" stroke="{BLUE_LO}" stroke-width="0.8"/></g>')
    parts.append(f'<path d="{spiral(bx - 0.5, by + 0.5, 2, 9.5, 1.8)}" fill="none" stroke="{BLUE_LINE}" stroke-opacity="0.25" stroke-width="1" stroke-linecap="round"/>')
    parts.append(f'<ellipse cx="{bx - 4.5}" cy="{by - 5}" rx="4" ry="2.6" fill="#ffffff" fill-opacity="0.35" transform="rotate(-30 {bx - 4.5} {by - 5})"/>')
    return "\n      ".join(parts)


def cell(cx, cy, L, Wd, ang, gid, pal="verm", small=False):
    """A spindle-shaped smooth muscle cell with a nucleus, lit from the upper left."""
    hi, mid, lo, line = CELL[pal]
    body = smooth_path(spindle(cx, cy, L, Wd, ang))
    nuc = smooth_path(spindle(cx, cy, L * 0.36, Wd * 0.6, ang))
    deg = math.degrees(ang)
    grad = (f'<linearGradient id="{gid}" gradientUnits="userSpaceOnUse" x1="{cx}" y1="{cy - Wd}" x2="{cx}" y2="{cy + Wd}" gradientTransform="rotate({deg:.1f} {cx} {cy})">'
            f'<stop offset="0" stop-color="{hi}"/><stop offset="0.5" stop-color="{mid}"/><stop offset="1" stop-color="{lo}"/></linearGradient>')
    sw = 0.5 if small else 0.7
    shape = (f'<path d="{body}" fill="url(#{gid})" stroke="{lo}" stroke-width="{sw}"/>'
             f'<path d="{nuc}" fill="{line}" fill-opacity="0.85"/>')
    if not small:
        hx, hy = cx - L * 0.12 * math.cos(ang) - Wd * 0.18 * math.sin(ang), cy - L * 0.12 * math.sin(ang) + Wd * 0.18 * math.cos(ang)
        shape += f'<ellipse cx="{hx:.1f}" cy="{hy:.1f}" rx="{L * 0.14:.1f}" ry="{Wd * 0.14:.1f}" fill="#ffffff" fill-opacity="0.45" transform="rotate({deg:.1f} {hx:.1f} {hy:.1f})"/>'
    return grad, shape


LINES = (("WT", "stone"), ("KO", "ko"), ("Mut1", "verm"), ("Mut2", "pink"))
LINE_X = (610, 650, 690, 730)          # the four lines side by side, after the editing step


def sphere(sx, sy, sr, pal):
    hi, mid, lo, line = CELL[pal]
    gid = {"verm": "gSphereVerm", "pink": "gSpherePink", "stone": "gSphereStone"}[pal]
    parts = [f'<ellipse cx="{sx + 1}" cy="{sy + sr + 2}" rx="{sr * 0.9}" ry="{sr * 0.22}" fill="#5f7384" fill-opacity="0.18" filter="url(#blur6)"/>',
             f'<circle cx="{sx}" cy="{sy}" r="{sr}" fill="url(#{gid})" stroke="{lo}" stroke-width="0.7"/>']
    for a_, rr in ((0.9, 0.55), (2.4, 0.62), (4.0, 0.5), (5.3, 0.6)):
        px, py = sx + sr * rr * math.cos(a_), sy + sr * rr * math.sin(a_)
        parts.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{sr * 0.3:.1f}" fill="none" stroke="{line}" stroke-opacity="0.28" stroke-width="0.8"/>')
    parts.append(f'<ellipse cx="{sx - sr * 0.35}" cy="{sy - sr * 0.4}" rx="{sr * 0.32}" ry="{sr * 0.22}" fill="#ffffff" fill-opacity="0.45" transform="rotate(-35 {sx - sr * 0.35} {sy - sr * 0.4})"/>')
    return "".join(parts)


def cultures():
    rng = np.random.default_rng(5)
    grads, parts = [], []
    # cells taken from the myometrium: still the tissue's colour
    for i, (cx, cy, L, Wd, ang) in enumerate(((508, 204, 44, 12.5, -0.38), (530, 224, 38, 11.5, 0.55), (502, 232, 36, 11, 1.3))):
        g, sh = cell(cx, cy, L, Wd, ang, f"gDer{i}", "stone")
        grads.append(g); parts.append(f'<g filter="url(#castSoft)">{sh}</g>')
    # after CRISPR: the four lines, each in its colour
    for (name, pal), x in zip(LINES, LINE_X):
        g, sh = cell(x, 214, 36, 11.5, -0.3, f"gLine{name}", pal)
        grads.append(g); parts.append(f'<g filter="url(#castSoft)">{sh}</g>')
    # 2D culture: one dish stands for the four lines, each grown separately
    x0, x1 = 800, 962
    dish = (f"M{x0} 258 L{x0} 290 Q{x0} 304 {x0 + 14} 304 L{x1 - 14} 304 Q{x1} 304 {x1} 290 L{x1} 258 L{x1 - 6} 258 L{x1 - 6} 288 "
            f"Q{x1 - 6} 298 {x1 - 16} 298 L{x0 + 16} 298 Q{x0 + 6} 298 {x0 + 6} 288 L{x0 + 6} 258 Z")
    parts.append(f'<path d="M{x0 + 6} 282 L{x1 - 6} 282 L{x1 - 6} 288 Q{x1 - 6} 298 {x1 - 16} 298 L{x0 + 16} 298 Q{x0 + 6} 298 {x0 + 6} 288 Z" fill="url(#gMedium)"/>')
    layer, k, x = [], 0, x0 + 12
    seg = (x1 - x0 - 24) / 4
    for (name, pal), s0 in zip(LINES, range(4)):
        xx = x0 + 12 + s0 * seg + 2
        while xx < x0 + 12 + (s0 + 1) * seg - 4:
            L = rng.uniform(13, 17)
            g, sh = cell(xx + L / 2, 295.3 + rng.uniform(-0.4, 0.4), L, rng.uniform(4.0, 4.8), rng.uniform(-0.16, 0.16), f"gLay{k}", pal, small=True)
            grads.append(g); layer.append(sh); xx += L * rng.uniform(0.92, 1.02); k += 1
    parts.append("".join(layer))
    parts.append(f'<path d="{dish}" fill="url(#gGlass)" stroke="{GLASS_LINE}" stroke-width="0.8" filter="url(#inGlass)"/>')
    parts.append(f'<path d="M{x0 + 9} 262 L{x0 + 9} 286" stroke="#ffffff" stroke-opacity="0.9" stroke-width="1.6" stroke-linecap="round"/>')
    parts.append(f'<path d="M{x0 + 7} 282 L{x1 - 7} 282" stroke="#ffffff" stroke-opacity="0.7" stroke-width="1"/>')
    # 3D culture: a round-bottom well with the three lines grown as spheroids
    w0, w1 = 812, 950
    well = (f"M{w0} 56 L{w0} 108 Q{w0} 152 {w0 + 44} 152 L{w1 - 44} 152 Q{w1} 152 {w1} 108 L{w1} 56 L{w1 - 6} 56 L{w1 - 6} 108 "
            f"Q{w1 - 6} 146 {w1 - 44} 146 L{w0 + 44} 146 Q{w0 + 6} 146 {w0 + 6} 108 L{w0 + 6} 56 Z")
    parts.append(f'<path d="M{w0 + 6} 92 L{w1 - 6} 92 L{w1 - 6} 108 Q{w1 - 6} 146 {w1 - 44} 146 L{w0 + 44} 146 Q{w0 + 6} 146 {w0 + 6} 108 Z" fill="url(#gMedium)"/>')
    for (sx, sy, sr), pal in zip(((852, 121, 12.5), (902, 126, 11.5), (880, 97, 10.5)), ("stone", "pink", "verm")):
        parts.append(sphere(sx, sy, sr, pal))
    parts.append(f'<path d="{well}" fill="url(#gGlass)" stroke="{GLASS_LINE}" stroke-width="0.8" filter="url(#inGlass)"/>')
    parts.append(f'<path d="M{w0 + 9} 60 L{w0 + 9} 104" stroke="#ffffff" stroke-opacity="0.9" stroke-width="1.6" stroke-linecap="round"/>')
    parts.append(f'<path d="M{w0 + 7} 92 L{w1 - 7} 92" stroke="#ffffff" stroke-opacity="0.7" stroke-width="1"/>')
    return "\n    ".join(grads), "\n      ".join(parts)


def labels():
    line_labels = "".join(f'<text x="{x}" y="243" text-anchor="middle">{name}</text>' for (name, _), x in zip(LINES, LINE_X))
    return f"""
      <g stroke="{ARROW}" stroke-width="1.8" fill="none" stroke-linecap="round">
        <path d="M282 238 C 340 262, 410 250, 474 222" marker-end="url(#arrow)"/>
        <path d="M548 214 L 588 214" marker-end="url(#arrow)"/>
        <path d="M744 206 C 768 180, 786 140, 804 112" marker-end="url(#arrow)"/>
        <path d="M744 224 C 768 252, 782 272, 794 286" marker-end="url(#arrow)"/>
      </g>
      <g stroke="{ARROW}" stroke-width="1.2" fill="none" stroke-linecap="round">
        <path d="M154 204 L 118 246"/>
        <path d="M186 268 L 138 306"/>
        <path d="M236 192 L 326 200"/>
        <path d="M196 74 L 158 48"/>
      </g>
      <g fill="{INK}"><circle cx="154" cy="204" r="2.4"/><circle cx="186" cy="268" r="2.4"/><circle cx="236" cy="192" r="2.4"/><circle cx="196" cy="74" r="2.4"/></g>
      <g font-size="17" fill="{INK}">
        <text x="112" y="251" text-anchor="end">Leiomyoma</text>
        <text x="152" y="46" text-anchor="end">Leiomyoma</text>
        <text x="134" y="311" text-anchor="end">Myometrium</text>
        <text x="332" y="205">Endometrium</text>
        <text x="524" y="272" text-anchor="middle">Myometrial cells</text>
        <text x="880" y="334" text-anchor="middle">2D culture</text>
        <text x="880" y="44" text-anchor="middle">3D spheroids</text>
      </g>
      <g font-size="14" fill="#5c6170">
        <text x="112" y="269" text-anchor="end">MED12 G44D</text>
        <text x="152" y="64" text-anchor="end">MED12 wild-type</text>
        <text x="568" y="198" text-anchor="middle">CRISPR editing</text>
        {line_labels}
      </g>"""


def main():
    cell_grads, culture_parts = cultures()
    svg = f"""<svg class="design-svg" viewBox="0 0 1000 360" role="img" aria-labelledby="design-title" xmlns="http://www.w3.org/2000/svg">
      <title id="design-title">Illustration: a frontal section of the uterus with two leiomyomas in the myometrial wall, one carrying MED12 G44D and one MED12 wild-type; myometrial cells taken from it and edited with CRISPR into four lines (unedited WT, MED12 knock-out, and two codon-44 mutant clones Mut1 and Mut2), grown as 2D cultures and, for WT, Mut1 and Mut2, as 3D spheroids.</title>{defs()}
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
