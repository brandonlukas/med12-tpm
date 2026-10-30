"""Paint the overview illustration in watercolour: line and wash.

    uv run --with numpy --with opencv-python-headless --with pillow python illustration/paint.py

Geometry is written in the page's viewBox units (900 x 360) and painted at 2400 px wide.
Each shape becomes a soft mask; a wash is built from it with noise-displaced edges, pigment
granulation, pooling at the edges and wet blooms, then layered over a transparent canvas. A
loose pen line goes over the washes. Labels and arrows are typeset on the page as SVG, not
painted here. Writes docs/img/design.webp (transparent) and illustration/paint_check.png.
"""
import math
import re
from pathlib import Path

import cv2
import numpy as np
from PIL import Image, ImageDraw

HERE = Path(__file__).resolve().parent
VW, VH = 900, 360
S = 2400 / VW
W, H = int(VW * S), int(VH * S)
SS = 3                                   # supersampling for anti-aliased masks
rng = np.random.default_rng(11)

# the site's role colours as pigments, plus the glass of the dishes and the pen
TISSUE, TISSUE_DEEP, LINING = "#9E9398", "#6b6066", "#c4a2a8"
TUMOUR, TUMOUR_DEEP, TUMOUR_LINE = "#D55E00", "#9c3f00", "#7a3000"
GLASS, GLASS_DEEP = "#b9c6d1", "#8da0b0"
INK = "#2b2a30"


def hex_rgb(h):
    return np.array([int(h[i:i + 2], 16) for i in (1, 3, 5)], float) / 255


# ----------------------------------------------------------------------------- geometry
def flatten(d, n=28):
    """An SVG path (M, L, C, Q, Z; absolute) as a list of polygons in viewBox units."""
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
                p = (1 - s) ** 3 * p0 + 3 * (1 - s) ** 2 * s * p1 + 3 * (1 - s) * s ** 2 * p2 + s ** 3 * p3
                cur.append(tuple(p))
            pos = tuple(p3); i += 7
        elif t == "Q":
            p0 = np.array(pos); p1, p2 = (np.array([float(tokens[i + k]), float(tokens[i + k + 1])]) for k in (1, 3))
            for s in np.linspace(0, 1, n)[1:]:
                p = (1 - s) ** 2 * p0 + 2 * (1 - s) * s * p1 + s ** 2 * p2
                cur.append(tuple(p))
            pos = tuple(p2); i += 5
        elif t == "Z":
            if cur: polys.append(cur); cur = []
            i += 1
        else:
            i += 1
    if cur: polys.append(cur)
    return polys


def circle(cx, cy, r, n=72):
    return [[(cx + r * math.cos(a), cy + r * math.sin(a)) for a in np.linspace(0, 2 * math.pi, n, endpoint=False)]]


def spindle(cx, cy, length, width, angle, n=60):
    """A smooth muscle cell: a lens that tapers to both ends."""
    pts = []
    for a in np.linspace(0, 2 * math.pi, n, endpoint=False):
        x = math.cos(a) * length / 2
        y = math.sin(a) * width / 2 * (1 - 0.55 * abs(math.cos(a)) ** 2)
        c, s = math.cos(angle), math.sin(angle)
        pts.append((cx + x * c - y * s, cy + x * s + y * c))
    return [pts]


def stroke_polys(points, widths):
    """The outline of a stroke along points, with a width per point, as one polygon."""
    p = np.asarray(points, float)
    t = np.gradient(p, axis=0); t /= np.maximum(np.hypot(*t.T), 1e-9)[:, None]
    nrm = np.column_stack([-t[:, 1], t[:, 0]])
    w = np.asarray(widths, float)[:, None] / 2
    return [list(map(tuple, np.vstack([p + nrm * w, (p - nrm * w)[::-1]])))]


def densify(poly, step=1.0):
    """Points along a closed polygon every `step` viewBox units."""
    p = np.asarray(poly + [poly[0]], float)
    seg = np.hypot(*np.diff(p, axis=0).T)
    out = [p[0]]
    for a, b, L in zip(p[:-1], p[1:], seg):
        for s in np.linspace(0, 1, max(2, int(L / step)) + 1)[1:]:
            out.append(a + (b - a) * s)
    return np.array(out[:-1])


# ----------------------------------------------------------------------------- raster
def mask(polys, soften=0.0):
    """Anti-aliased float mask of filled polygons."""
    img = Image.new("L", (W * SS, H * SS), 0)
    dr = ImageDraw.Draw(img)
    for poly in polys:
        dr.polygon([(x * S * SS, y * S * SS) for x, y in poly], fill=255)
    m = np.asarray(img.resize((W, H), Image.LANCZOS), float) / 255
    if soften:
        m = cv2.GaussianBlur(m, (0, 0), soften)
    return np.clip(m, 0, 1)


def value_noise(cells, seed):
    r = np.random.default_rng(seed)
    g = r.random((max(2, int(H / W * cells)) + 1, cells + 1)).astype(np.float32)
    return cv2.resize(g, (W, H), interpolation=cv2.INTER_CUBIC)


def fractal(cells, octaves, seed, gain=0.55):
    n, amp, total = np.zeros((H, W), np.float32), 1.0, 0.0
    for k in range(octaves):
        n += amp * value_noise(cells * 2 ** k, seed + k); total += amp; amp *= gain
    n /= total
    return (n - n.min()) / (n.max() - n.min() + 1e-9)


def displace(m, amp, cells, seed):
    """Warp a mask with a smooth noise field: the ragged edge of a wash."""
    dx = (fractal(cells, 2, seed) - 0.5) * 2 * amp
    dy = (fractal(cells, 2, seed + 50) - 0.5) * 2 * amp
    ys, xs = np.mgrid[0:H, 0:W].astype(np.float32)
    return cv2.remap(m.astype(np.float32), xs + dx, ys + dy, cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT)


class Canvas:
    def __init__(self):
        self.rgb = np.zeros((H, W, 3), np.float32)
        self.a = np.zeros((H, W), np.float32)

    def over(self, colour, alpha):
        """Composite a flat colour with a per-pixel alpha onto the canvas (straight alpha)."""
        c = hex_rgb(colour).astype(np.float32)
        a = np.clip(alpha, 0, 1).astype(np.float32)
        out_a = a + self.a * (1 - a)
        num = c[None, None, :] * a[..., None] + self.rgb * (self.a * (1 - a))[..., None]
        self.rgb = np.where(out_a[..., None] > 1e-6, num / np.maximum(out_a, 1e-6)[..., None], self.rgb)
        self.a = out_a

    def save(self, path, preview):
        rgba = np.dstack([self.rgb, self.a])
        img = Image.fromarray((np.clip(rgba, 0, 1) * 255).astype(np.uint8), "RGBA")
        img.save(path, quality=86, method=6)
        white = Image.new("RGBA", img.size, (255, 255, 255, 255))
        Image.alpha_composite(white, img).convert("RGB").save(preview)


def wash(canvas, m, colour, alpha, seed, deep=None, pool=0.45, pool_width=16, grain=0.4,
         wobble=10, cells=10, blooms=2, blur=1.6):
    """One wash of pigment from a mask: ragged edge, granulation, pooling at the edge, blooms."""
    r = np.random.default_rng(seed)
    m1 = displace(m, wobble, cells, seed)
    m1 = cv2.GaussianBlur(m1, (0, 0), blur)
    # pigment settles unevenly: fine granulation and a few lighter blooms
    g = fractal(90, 3, seed + 7)
    density = 1 - grain + grain * g ** 1.4
    for _ in range(blooms):
        ys, xs = np.nonzero(m1 > 0.5)
        if not len(xs): break
        k = r.integers(len(xs)); cx, cy = xs[k], ys[k]
        rad = r.uniform(0.08, 0.2) * math.sqrt(len(xs))
        yy, xx = np.mgrid[0:H, 0:W]
        density *= 1 - 0.38 * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * rad ** 2)).astype(np.float32)
    canvas.over(colour, m1 * alpha * density)
    if deep and pool:
        hard = (m1 > 0.5).astype(np.uint8)
        dist = cv2.distanceTransform(hard, cv2.DIST_L2, 5)
        edge = np.exp(-dist / pool_width) * (hard > 0)
        edge = cv2.GaussianBlur(edge.astype(np.float32), (0, 0), 1.2)
        canvas.over(deep, edge * alpha * pool * (0.7 + 0.3 * g))


def pen(canvas, polys, width=3.2, alpha=0.88, wobble=1.6, closed=True, colour=INK, seed=3):
    """A loose pen line: slightly off the edge, pressure varying, drawn with round joins."""
    r = np.random.default_rng(seed)
    img = Image.new("L", (W * SS, H * SS), 0)
    dr = ImageDraw.Draw(img)
    for poly in polys:
        pts = densify(poly, 1.2) if closed else np.asarray(poly, float)
        n = len(pts)
        ph = r.uniform(0, 6.28, 3)
        k = np.arange(n) / n * 2 * math.pi
        off = wobble * (0.6 * np.sin(2 * k + ph[0]) + 0.3 * np.sin(5 * k + ph[1]) + 0.2 * np.sin(11 * k + ph[2]))
        t = np.gradient(pts, axis=0); t /= np.maximum(np.hypot(*t.T), 1e-9)[:, None]
        nrm = np.column_stack([-t[:, 1], t[:, 0]])
        q = pts + nrm * off[:, None]
        wv = width * (0.78 + 0.22 * np.sin(3 * k + ph[1]) + 0.12 * np.sin(7 * k + ph[2]))
        if closed:
            q = np.vstack([q, q[:1]]); wv = np.append(wv, wv[0])
        for (x0, y0), (x1, y1), w0 in zip(q[:-1], q[1:], wv):
            dr.line([(x0 * S * SS, y0 * S * SS), (x1 * S * SS, y1 * S * SS)], fill=255, width=max(1, int(round(w0 * SS))))
            dr.ellipse([(x1 * S * SS - w0 * SS / 2, y1 * S * SS - w0 * SS / 2), (x1 * S * SS + w0 * SS / 2, y1 * S * SS + w0 * SS / 2)], fill=255)
    m = np.asarray(img.resize((W, H), Image.LANCZOS), np.float32) / 255
    canvas.over(colour, m * alpha)


def fibres(canvas, m, colour, period=9.0, alpha=0.11, seed=5):
    """Muscle fibres: broken lines parallel to the shape's edge, inside it."""
    hard = (m > 0.5).astype(np.uint8)
    dist = cv2.distanceTransform(hard, cv2.DIST_L2, 5)
    lines = 0.5 + 0.5 * np.cos(2 * math.pi * dist / (period * S / 2.667))
    lines = np.clip((lines - 0.72) / 0.28, 0, 1)
    broken = (fractal(40, 3, seed) > 0.42).astype(np.float32)
    broken = cv2.GaussianBlur(broken, (0, 0), 2)
    a = lines * broken * cv2.GaussianBlur(m, (0, 0), 3) * (dist > 6)
    canvas.over(colour, cv2.GaussianBlur(a, (0, 0), 0.8) * alpha)


def whorl(canvas, m, cx, cy, colour, alpha=0.14, period=6.5, seed=9):
    """The whorled cut surface of a leiomyoma: a loose spiral inside the nodule."""
    yy, xx = np.mgrid[0:H, 0:W].astype(np.float32)
    rr = np.hypot(xx - cx * S, yy - cy * S); th = np.arctan2(yy - cy * S, xx - cx * S)
    sp = 0.5 + 0.5 * np.cos(2 * math.pi * (rr / (period * S) - 0.9 * th / (2 * math.pi)))
    sp = np.clip((sp - 0.7) / 0.3, 0, 1)
    broken = cv2.GaussianBlur((fractal(60, 2, seed) > 0.4).astype(np.float32), (0, 0), 1.5)
    inner = cv2.erode((m > 0.5).astype(np.uint8), np.ones((9, 9), np.uint8)).astype(np.float32)
    canvas.over(colour, cv2.GaussianBlur(sp * broken * inner, (0, 0), 0.8) * alpha)


# ----------------------------------------------------------------------------- the scene
def main():
    cv = Canvas()

    # the uterus in frontal section: a domed fundus, the body narrowing to the isthmus,
    # a cylindrical cervix with a rounded end
    body = flatten("M220 72 C 262 72, 300 90, 306 120 C 312 150, 296 196, 266 232 C 258 242, 254 250, 254 258 "
                   "C 256 288, 250 304, 236 309 C 228 312, 212 312, 204 309 C 190 304, 184 288, 186 258 "
                   "C 186 250, 182 242, 174 232 C 144 196, 128 150, 134 120 C 140 90, 178 72, 220 72 Z")
    cavity = flatten("M178 124 C 196 116, 244 116, 262 124 C 250 158, 232 200, 224 234 Q 220 242, 216 234 "
                     "C 208 200, 190 158, 178 124 Z")
    canal = stroke_polys(np.array([(220, 242), (220, 302)]), [4.2, 4.2])

    tubes, fimbriae, ovaries = [], [], []
    for sign in (1, -1):
        cx = lambda x: 220 + sign * (x - 220)
        pts = np.array([(cx(x), y) for x, y in flatten("M302 108 C 326 94, 348 94, 364 108 C 376 118, 382 132, 394 142")[0]])
        tubes += stroke_polys(pts, np.linspace(7, 11, len(pts)))
        end, heading = pts[-1], pts[-1] - pts[-6]
        base = math.atan2(heading[1], heading[0])
        for k, spread in enumerate((-0.9, -0.45, 0.0, 0.45, 0.9)):
            a = base + spread
            L = (13, 17, 19, 17, 13)[k]
            tip = end + np.array([math.cos(a), math.sin(a)]) * L
            mid = (end + tip) / 2 + np.array([-math.sin(a), math.cos(a)]) * 1.5 * np.sign(spread or 1)
            fimbriae += stroke_polys(np.vstack([end, mid, tip]), [4.8, 3.2, 1.4])
        ox0, oy0 = cx(404), 164
        rot = math.radians(-18 * sign)
        ovaries += [[(ox0 + (11 * math.cos(a)) * math.cos(rot) - (16 * math.sin(a)) * math.sin(rot),
                      oy0 + (11 * math.cos(a)) * math.sin(rot) + (16 * math.sin(a)) * math.cos(rot))
                     for a in np.linspace(0, 2 * math.pi, 72, endpoint=False)]]

    m_body, m_cav, m_canal = mask(body), mask(cavity), mask(canal)
    m_tubes, m_fimb, m_ov = mask(tubes), mask(fimbriae), mask(ovaries)
    wall = np.clip(m_body - cv2.GaussianBlur(m_cav, (0, 0), 1.0) - cv2.GaussianBlur(m_canal, (0, 0), 0.8), 0, 1)
    fib_c = (168, 190)
    m_fib = mask(circle(*fib_c, 20))

    # washes, back to front
    wash(cv, m_fimb, TISSUE, 0.4, 21, deep=TISSUE_DEEP, pool=0.3, pool_width=4, wobble=1.5, cells=30, blooms=0, blur=1.0)
    wash(cv, m_tubes, TISSUE, 0.5, 22, deep=TISSUE_DEEP, pool=0.45, pool_width=6, wobble=3, cells=18, blooms=1)
    wash(cv, m_ov, TISSUE, 0.5, 23, deep=TISSUE_DEEP, pool=0.55, pool_width=6, wobble=2.5, cells=24, blooms=1, grain=0.5)
    wash(cv, wall, TISSUE, 0.5, 31, deep=TISSUE_DEEP, pool=0.62, pool_width=16, wobble=9, cells=9, blooms=3)
    wash(cv, wall, TISSUE, 0.28, 32, deep=None, wobble=6, cells=12, blooms=2, grain=0.55)
    # pigments mingling wet-in-wet: a rose and a cool grey drift through the wall in patches
    drift = cv2.GaussianBlur(wall, (0, 0), 4)
    cv.over("#b48f96", drift * np.clip(fractal(3, 2, 33) - 0.45, 0, 1) * 0.7)
    cv.over("#8d93a3", drift * np.clip(fractal(3, 2, 34) - 0.5, 0, 1) * 0.6)
    fibres(cv, wall, TISSUE_DEEP, period=11.5, alpha=0.08)
    # the endometrium: a soft lining just inside the cavity, and along the canal
    hard = (m_cav > 0.5).astype(np.uint8)
    dist_in = cv2.distanceTransform(hard, cv2.DIST_L2, 5)
    lining = np.exp(-dist_in / 5.5) * (hard > 0)
    cv.over(LINING, cv2.GaussianBlur(lining.astype(np.float32), (0, 0), 1.5) * 0.5)
    cv.over(LINING, (cv2.GaussianBlur(m_canal, (0, 0), 2.2) - m_canal).clip(0, 1) * 0.6)
    # the leiomyoma in the wall, with its whorled cut surface
    wash(cv, m_fib, TUMOUR, 0.7, 51, deep=TUMOUR_DEEP, pool=0.6, pool_width=11, wobble=3.5, cells=20, blooms=2, grain=0.35)
    wash(cv, mask(circle(fib_c[0] - 2, fib_c[1] + 3, 15)), TUMOUR, 0.26, 52, deep=None, wobble=3, cells=20, blooms=1)
    whorl(cv, m_fib, *fib_c, TUMOUR_LINE, alpha=0.15, period=6)

    # the derived cells: smooth muscle cells, spindle-shaped, in a small cluster
    cells, nuclei = [], []
    for cx0, cy0, L, Wd, ang in ((548, 204, 44, 12, -0.4), (568, 224, 38, 11, 0.55), (540, 232, 36, 10.5, 1.3)):
        cells += spindle(cx0, cy0, L, Wd, ang)
        nuclei += spindle(cx0, cy0, L * 0.36, Wd * 0.58, ang)
    for i, (c, nu) in enumerate(zip(cells, nuclei)):
        wash(cv, mask([c]), TUMOUR, 0.62, 61 + i, deep=TUMOUR_DEEP, pool=0.5, pool_width=5, wobble=1.4, cells=30, blooms=1, grain=0.3)
        cv.over(TUMOUR_DEEP, cv2.GaussianBlur(mask([nu]), (0, 0), 1.6) * 0.42)

    # 2D culture: a dish in section, medium, and a single layer of cells lying on the floor
    dish = flatten("M652 262 L 652 290 Q 652 302 664 302 L 788 302 Q 800 302 800 290 L 800 262")
    wash(cv, mask([[(654, 283), (798, 283), (798, 300), (654, 300)]]), GLASS, 0.3, 71, deep=GLASS_DEEP, pool=0.35, pool_width=8, wobble=3, cells=16, blooms=1, grain=0.3)
    layer, layer_nuclei = [], []
    x = 660
    while x < 790:
        L = rng.uniform(16, 22)
        ang = rng.uniform(-0.18, 0.18)
        layer += spindle(x + L / 2, 296 + rng.uniform(-0.5, 0.5), L, rng.uniform(4.0, 5.0), ang)
        layer_nuclei += spindle(x + L / 2, 296, L * 0.34, 2.6, ang)
        x += L * rng.uniform(0.86, 1.0)
    for i, (c, nu) in enumerate(zip(layer, layer_nuclei)):
        wash(cv, mask([c]), TUMOUR, 0.64, 81 + i, deep=TUMOUR_DEEP, pool=0.5, pool_width=3, wobble=0.7, cells=40, blooms=0, grain=0.3)
        cv.over(TUMOUR_DEEP, cv2.GaussianBlur(mask([nu]), (0, 0), 1.2) * 0.4)

    # 3D culture: a round-bottom well, medium, three spheroids with bumpy, shaded surfaces
    well = flatten("M660 60 L 660 108 Q 660 150 702 150 L 758 150 Q 800 150 800 108 L 800 60")
    medium3 = flatten("M662 96 L 662 108 Q 662 148 702 148 L 758 148 Q 798 148 798 108 L 798 96 Z")
    wash(cv, mask(medium3), GLASS, 0.3, 91, deep=GLASS_DEEP, pool=0.35, pool_width=10, wobble=4, cells=14, blooms=1, grain=0.3)
    spheroids = ((702, 122, 12.5), (748, 128, 11.5), (728, 94, 10.5))
    yy, xx = np.mgrid[0:H, 0:W].astype(np.float32)
    for i, (cx0, cy0, r0) in enumerate(spheroids):
        m_s = displace(mask(circle(cx0, cy0, r0)), 2.6, 80, 101 + i)
        wash(cv, m_s, TUMOUR, 0.72, 111 + i, deep=TUMOUR_DEEP, pool=0.75, pool_width=6, wobble=0.8, cells=40, blooms=1, grain=0.55)
        # the shadow side, lower right, and a soft highlight upper left
        shade = np.clip(((xx - cx0 * S) + (yy - cy0 * S)) / (2.2 * r0 * S) + 0.35, 0, 1)
        cv.over(TUMOUR_DEEP, cv2.GaussianBlur(m_s * shade, (0, 0), 3) * 0.35)
        hl = np.exp(-(((xx - (cx0 - r0 * 0.38) * S) ** 2 + (yy - (cy0 - r0 * 0.42) * S) ** 2) / (2 * (r0 * S * 0.3) ** 2)))
        cv.over("#ffffff", cv2.GaussianBlur(hl * m_s, (0, 0), 3) * 0.26)

    # pen: a loose line round the silhouettes, lighter on the small parts
    pen(cv, body, width=3.0, wobble=1.6, alpha=0.82, seed=1)
    pen(cv, tubes, width=2.0, wobble=0.9, alpha=0.7, seed=2)
    pen(cv, ovaries, width=2.0, wobble=0.9, alpha=0.7, seed=4)
    pen(cv, cavity, width=1.5, wobble=0.8, alpha=0.45, seed=5)
    pen(cv, [circle(*fib_c, 20)[0]], width=2.4, wobble=1.2, colour=TUMOUR_LINE, alpha=0.75, seed=6)
    for i, c in enumerate(cells):
        pen(cv, [c], width=1.8, wobble=0.5, colour=TUMOUR_LINE, alpha=0.7, seed=20 + i)
    pen(cv, dish, width=2.8, wobble=0.8, closed=False, seed=7)
    pen(cv, well, width=2.8, wobble=0.8, closed=False, seed=8)
    for i, c in enumerate(layer):
        pen(cv, [c], width=1.3, wobble=0.25, colour=TUMOUR_LINE, alpha=0.55, seed=30 + i)
    for i, (cx0, cy0, r0) in enumerate(spheroids):
        pen(cv, circle(cx0, cy0, r0, 90), width=1.8, wobble=1.3, colour=TUMOUR_LINE, alpha=0.55, seed=40 + i)

    out = HERE.parent / "docs" / "img"
    out.mkdir(exist_ok=True)
    cv.save(out / "design.webp", HERE / "paint_check.png")
    print("wrote", out / "design.webp", f"{(out / 'design.webp').stat().st_size // 1024} KB")


if __name__ == "__main__":
    main()
