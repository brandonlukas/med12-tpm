# The sketched half: one scene of simple solids, drawn in either of two
# finishes. Save this file beside ruled.py. Tested with Matplotlib 3.11
# and NumPy 2.
#
#   from ruled import *
#   from sketched import *
#
#   scene = Scene()
#   scene.box(0, 0, 0, 90, 6, 60)                       # the stage
#   scene.box(18, 6, 19, 54, 1.5, 22)                   # a glass slide
#   scene.disc(45, 7.5, 30, 5, PERSIMMON)               # the sample
#   with scene.object():                                # the objective,
#       scene.cylinder(45, 26, 30, 4, 8, r_top=9)       # from its lens up
#       scene.cylinder(45, 34, 30, 9, 22)
#       scene.cylinder(45, 56, 30, 11, 5)
#   scene.flow([(64, 7.5, 52), (88, 7.5, 52)])
#   scene.label("Objective", (62, 62), (51, 48, 30))
#   scene.label("Sample", (-36, 14), (41, 7.5, 30), ha="right")
#   scene.label("Stage scans", (86, -30), (80, 7.5, 54))
#
#   fig = plt.figure(figsize=(89 * MM, 60 * MM))        # for a figure
#   ax = drawing(fig)
#   scene.draw(ax)
#   fit(ax, *columns(1, 6, 89), top=1, bottom=0, letter="a")
#   fig.savefig("schematic.pdf")
#
#   fig, _ = slide("We scanned each sample under a 40× objective.",
#                  chart=False)                         # for a slide,
#   scene.draw(fig)                                     # either finish:
#   fig.savefig("method_line.png")                      # the line drawing
#
#   fig, _ = slide("We scanned each sample under a 40× objective.",
#                  chart=False)
#   scene.block((6, -8), (52, 70))                      # or the sketch
#   scene.sketch(fig, "Mynerve-Regular.ttf")
#   fig.savefig("method_sketch.png")
import contextlib

import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import to_hex, to_rgb
from matplotlib.font_manager import FontProperties

from ruled import INK, MARKER_GREYS, MUSTARD, PENCIL, fit, plt

YAW, PITCH = np.radians(-30), np.radians(22)         # the house view
TONES = dict(zip(("light", "mid", "dark"), MARKER_GREYS))
MARKER = "#5C6A73"       # one grey marker: on white, 20, 40 and 66 % of it
COVER = {"light": 0.2, "mid": 0.4, "dark": 0.66}     # give the three greys


def _shade(colour):
    """A colour under one stroke of the grey marker: the shaded side of a
    part drawn in colour."""
    k = COVER["light"]
    return to_hex(np.multiply(to_rgb(colour), 1 - k)
                  + np.multiply(to_rgb(MARKER), k))


def view(x, y, z):
    """The house three-quarter view: where a point of a scene lands on the
    page. x runs right, y up and z towards the viewer, so a box shows its
    top, its front (on the left) and its right side."""
    across = x * np.cos(YAW) + z * np.sin(YAW)
    depth = -x * np.sin(YAW) + z * np.cos(YAW)
    return across, y * np.cos(PITCH) - depth * np.sin(PITCH)


def drawing(fig):
    """Axes for a line drawing or a sketch: equal scale both ways, no
    frame, and limits that hug what is drawn. Draw in any unit, then call
    fit(), which places the drawing by its top left corner."""
    ax = fig.add_subplot()
    ax.set_aspect("equal")
    ax.set_anchor("NW")
    ax.margins(0)
    ax.axis("off")
    return ax


def leader(ax, text, at, to, ha="left"):
    """Label a part of a drawing: the text at one point, and a 0.5 pt line
    from it that ends in a dot on the part. Returns the label and the dot."""
    label = ax.annotate(text, to, xytext=at, ha=ha, va="center", zorder=4,
                        annotation_clip=False,
                        arrowprops={"arrowstyle": "-", "lw": 0.5,
                                    "color": INK, "shrinkA": 1.5,
                                    "shrinkB": 0})
    dot, = ax.plot(*to, "o", ms=1.5, color=INK, mew=0, zorder=4,
                   clip_on=False)
    return label, dot


def paste(fig, image):
    """Put a raster picture on a slide pixel for pixel, its top left
    corner in the top left of fig.box: for a sketch you rendered yourself
    from a photo or a CAD view. Give it on PAPER as a uint8 RGB array no
    larger than 1728 x 660 px, and save the slide at its own dpi."""
    if image.shape[0] > 660 or image.shape[1] > 1728:
        raise ValueError("the picture is larger than 1728 x 660 px")
    left, _, top, _ = fig.box
    fig.figimage(image, xo=round(left * 1920),
                 yo=round(top * 1080) - image.shape[0], origin="upper")


def _spline(points, n=10):
    """A smooth curve through a few points, as a dense row of points."""
    p = np.asarray(points, float)
    q = np.vstack([p[0], p, p[-1]])
    t = np.linspace(0, 1, n, endpoint=False)[:, None]
    return np.vstack([0.5 * (2 * b + (c - a) * t
                             + (2 * a - 5 * b + 4 * c - d) * t ** 2
                             + (3 * b - a - 3 * c + d) * t ** 3)
                      for a, b, c, d in zip(q, q[1:], q[2:], q[3:])]
                     + [p[-1:]])


def _hand(row, rng, u, over=1.4):
    """A row of points as a hand draws it. A straight edge runs past both
    corners, bows a little and wavers. A closed curve is drawn in one go
    and overlaps its own start. u is one point in drawing units."""
    p = np.asarray(row, float)
    wave = lambda x, f: np.sin(f * x + rng.uniform(0, 6.3))
    if len(p) == 2:
        a, b = p
        length = max(np.hypot(*(b - a)), 1e-9)
        t = (b - a) / length
        past = min(over * u, 0.16 * length) * rng.uniform(0.35, 1.15, 2)
        s = np.linspace(-past[0], length + past[1],
                        max(6, int((length + past.sum()) / (1.7 * u)) + 1))
        x = s / length
        bow = 0.005 * length * rng.uniform(-1, 1)
        off = (bow * np.sin(np.pi * x.clip(0, 1))
               + 0.22 * u * rng.uniform(-1, 1) * (2 * x - 1)
               + 0.13 * u * (0.6 * wave(x, 5.1) + 0.4 * wave(x, 11.3)))
        return a + np.outer(s, t) + np.outer(off, [-t[1], t[0]])
    closed = np.allclose(p[0], p[-1])
    if closed:
        n = len(p) - 1
        i = np.arange(int(n * rng.uniform(1.05, 1.12)) + 1)
        p = p[(i + rng.integers(n)) % n]
    x = np.linspace(0, 1, len(p))
    c = p.mean(axis=0)
    k = 1 + 0.012 * wave(x, 6.3) + 0.006 * wave(x, 12.6)
    if closed:
        k += 0.014 * rng.choice([-1, 1]) * (x - 0.5)
    return c + (p - c) * k[:, None]


def _pen(row, width, rng, u, ends=(0.5, 0.18)):
    """The outline of a pen stroke along a row of points: the pressure
    builds quickly, wavers and falls away towards the end."""
    p = np.asarray(row, float)
    s = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(p, axis=0).T))])
    tang = np.gradient(p, axis=0)
    tang /= np.maximum(np.hypot(*tang.T), 1e-12)[:, None]
    ease = lambda v: v.clip(0, 1) * (2 - v.clip(0, 1))
    total = max(s[-1], 1e-9)
    rise, fall = min(2.2 * u, total * 0.25), min(7 * u, total * 0.35)
    press = ((ends[0] + (1 - ends[0]) * ease(s / rise))
             * (ends[1] + (1 - ends[1]) * ease((total - s) / fall))
             * (1 + 0.1 * np.sin(4.4 * s / total + rng.uniform(0, 6.3))))
    side = np.column_stack([-tang[:, 1], tang[:, 0]])
    side *= (width * u / 2 * press)[:, None]
    return np.vstack([p + side, (p - side)[::-1]])


def _marker(region, angle, rng, u, nib=5.8, stray=0.45):
    """Chisel strokes that fill a region: rows along angle (in degrees),
    each a strip whose ends follow the edge of the region and stray a
    little over it. With a negative stray every stroke runs past the edge,
    for a region that is clipped. Returns the strips, and the fine pale
    streaks a marker leaves inside them."""
    c, s = np.cos(np.radians(angle)), np.sin(np.radians(angle))
    turn = np.array([[c, -s], [s, c]])        # strokes run along x
    q = np.asarray(region, float).reshape(-1, 2) @ turn
    if len(q) < 3:
        return [], []
    (x0, y0), (x1, y1) = q.min(axis=0), q.max(axis=0)
    nib = min(nib * u, max((y1 - y0) / 2.4, 1.2 * u))
    a, b = q, np.roll(q, -1, axis=0)

    def spans(y):                             # the row at y, inside the region
        cut = (a[:, 1] > y) != (b[:, 1] > y)
        x = np.sort(a[cut, 0] + (y - a[cut, 1]) * (b[cut, 0] - a[cut, 0])
                    / (b[cut, 1] - a[cut, 1]))
        return list(zip(x[::2], x[1::2]))

    rows = max(1, int(np.ceil((y1 - y0 - 0.16 * nib) / (0.84 * nib) - 1e-6)))
    pitch = (y1 - y0 - nib) / max(1, rows - 1)
    strips, streaks = [], []
    for r in range(rows):
        h = min(nib / 2 * rng.uniform(0.9, 1.08), (y1 - y0) / 2 * 1.02)
        yc = y0 + nib / 2 + r * pitch + rng.uniform(-0.06, 0.06) * nib
        if rows == 1:
            yc = (y0 + y1) / 2
        beside = spans(yc - 0.9 * h), spans(yc + 0.9 * h)
        for lo, hi in spans(yc):
            ends = []                         # where the nib's two edges end
            for row in beside:
                hit = [sp for sp in row if sp[1] > lo and sp[0] < hi]
                e = hit[0] if len(hit) == 1 else (lo, hi)
                ends.append((lo + (e[0] - lo).clip(-2 * nib, 2 * nib) / 0.9,
                             hi + (e[1] - hi).clip(-2 * nib, 2 * nib) / 0.9))
            (at, bt), (ab, bb) = ends
            ja, jb = stray * u, -stray * u
            if stray > 0:
                ja, jb = (rng.uniform(-stray, stray, 2) + [0.2, -0.2]) * u
            if max(bt - at, bb - ab) + jb - ja < 0.35 * nib:
                continue
            k0, k1 = rng.uniform(-0.1, 0.1, 2) * nib
            tilt = rng.uniform(-0.03, 0.03) * nib
            n = max(2, int((bt - at) / (3 * u)))
            wave = lambda: (np.linspace(-tilt, tilt, n) + 0.1 * u * np.sin(
                np.linspace(0, rng.uniform(2, 7), n) + rng.uniform(0, 6.3)))
            strips.append(np.vstack([
                np.column_stack([np.linspace(at + ja + k0, bt + jb + k1, n),
                                 yc - h + wave()]),
                np.column_stack([np.linspace(bb + jb - k1, ab + ja - k0, n),
                                 yc + h - wave()])]) @ turn.T)
            for _ in range(rng.integers(1, 4)):
                x = np.sort(rng.uniform(lo, hi, 2))
                y = yc + rng.uniform(-0.8, 0.8) * h
                w = rng.uniform(0.015, 0.05) * nib
                streaks.append(np.array([
                    [x[0], y - w], [x[1], y - w * 0.6], [x[1], y + w * 0.6],
                    [x[0], y + w]]) @ turn.T)
    return strips, streaks


class Scene:
    """One scene of simple solids, drawn in either finish from the same
    description: draw() makes the line drawing, for a figure or a slide,
    and sketch() the marker sketch, for a slide. Neither needs the other.
    x runs right, y up and z towards the viewer, in any unit.
    Add the parts from back to front and from the bottom up: each one
    hides what is behind it. A part is shaded in marker grey, or in the
    colour you give it: PERSIMMON for the thing the figure is about, GREY
    for its control. Where a point is asked for, three numbers are a place
    in the scene and two are a place on the page, as view() returns it.
    Each part returns where it lies on the page, as (left, bottom, right,
    top), to help place the labels.

    Lines come in two weights. The heavy one goes round the silhouette of
    each object; every other edge takes the thin one, including the edges
    where two parts of one object meet (see object()). An object less than
    a twelfth the size of the whole scene, a tube in a rack of twenty, is
    drawn in the thin line all round."""

    SMALL = 1 / 12

    def __init__(self, seed=1):
        self.parts, self.flows, self.labels, self.blocks = [], [], [], []
        self.seed = seed              # another seed, another hand
        self._object, self._count = None, 0

    @contextlib.contextmanager
    def object(self):
        """Parts added inside `with scene.object():` make one object, such
        as a tube and its cap: the heavy line goes round all of them, and
        where they meet the line is thin."""
        self._count += 1
        self._object = ("object", self._count)
        try:
            yield self
        finally:
            self._object = None

    def _objects(self):
        """The parts in drawing order, in runs that make one object each,
        and with each run whether its silhouette takes the heavy line."""
        page = np.vstack([np.vstack(p["edges"]) for p in self.parts
                          if "edges" in p] or [np.zeros((1, 2))])
        scene = max(np.ptp(page, axis=0))
        runs = []
        for part in self.parts:
            if runs and part.get("obj") is not None and \
                    part.get("obj") == runs[-1][0].get("obj"):
                runs[-1].append(part)
            else:
                runs.append([part])
        for run in runs:
            solid = [p for p in run if "edges" in p and not p.get("thin")]
            size = max(np.ptp(np.vstack([np.vstack(p["edges"])
                                         for p in solid]), axis=0)) \
                if solid else 0
            yield run, size >= self.SMALL * scene

    @staticmethod
    def _flat(points):
        return np.array([view(*p) if len(p) == 3 else p for p in points],
                        float)

    def _add(self, **part):
        part["obj"] = self._object
        self.parts.append(part)
        page = np.vstack(part["edges"])
        return tuple(float(v) for v in (*page.min(axis=0), *page.max(axis=0)))

    def box(self, x, y, z, w, h, d, colour=None):
        """A block: w wide, h high and d deep, from its back lower left
        corner at x, y, z."""
        p = lambda i, j, k: np.array(view(x + i * w, y + j * h, z + k * d))
        rim = [p(0, 0, 1), p(1, 0, 1), p(1, 0, 0), p(1, 1, 0), p(0, 1, 0),
               p(0, 1, 1)]
        near, far = p(1, 1, 1), p(0, 0, 0)
        return self._add(
            edges=[[a, b] for a, b in zip(rim, rim[1:] + rim[:1])],
            tones=[([rim[0], rim[1], near, rim[5]], colour or "light", 90),
                   ([rim[1], rim[2], rim[3], near],
                    _shade(colour) if colour else "mid", 90)],
            inside=[[near, rim[5]], [near, rim[1]], [near, rim[3]]],
            hidden=[[far, rim[2]], [far, rim[4]], [far, rim[0]]])

    def cylinder(self, x, y, z, r, h, r_top=None, colour=None):
        """An upright cylinder, the centre of its base at x, y, z. Give
        r_top for a cone or a taper; either radius may be 0."""
        r_top = r if r_top is None else r_top
        lean = np.arcsin(np.clip((r_top - r) / max(h, 1e-9) * np.tan(PITCH),
                                 -1, 1))

        def ring(radius, height, a0, a1, n=40):
            a = YAW + np.linspace(a0, a1, n)
            return np.array([view(x + radius * np.cos(t), height,
                                  z + radius * np.sin(t)) for t in a])

        def band(f0, f1, tone):       # down the side, from f0 to f1 across
            a0, a1 = (np.pi - lean) - np.array([f0, f1]) * (np.pi - 2 * lean)
            am = (a0 + a1) / 2
            way = ring(r_top, y + h, am, am, 1) - ring(r, y, am, am, 1)
            return (np.vstack([ring(r, y, a0, a1, 12),
                               ring(r_top, y + h, a1, a0, 12)]),
                    tone, np.degrees(np.arctan2(way[0, 1], way[0, 0])))

        front = ring(r, y, np.pi - lean, lean)
        back = ring(r_top, y + h, lean, -np.pi - lean)
        edges = [front, [front[-1], back[0]], back, [back[-1], front[0]]]
        axis = [view(x, y - 0.12 * h, z), view(x, y + 1.15 * h, z)]
        base = [ring(r, y, 0, 2 * np.pi, 49)] if r else []
        return self._add(
            guides=[[view(x, y, z), view(x, y + h, z)]] + base,
            edges=[e for e, there in zip(edges, (r, 1, r_top, 1)) if there],
            tones=([band(0, 1, colour), band(0.67, 1, _shade(colour))]
                   if colour else
                   [band(0.4, 1, "light"), band(0.67, 1, "mid")]),
            inside=[ring(r_top, y + h, 0, 2 * np.pi, 49)] if r_top else [],
            hidden=[axis] + base, round=True)

    def ball(self, x, y, z, r, colour=None):
        """A sphere with its centre at x, y, z."""
        c = np.array(view(x, y, z))
        a = np.radians(135) + np.linspace(0, 2 * np.pi, 73)
        circle = np.column_stack([np.cos(a), np.sin(a)])
        rim, lit = c + r * circle, c + r * np.array([-0.36, 0.38])

        def shade(k, tone):           # the ball outside a circle of light
            arc = lit + k * r * circle
            return (np.vstack([rim[np.hypot(*(rim - lit).T) > k * r],
                               arc[np.hypot(*(arc - c).T) < r][::-1]]),
                    tone, -40)

        return self._add(
            edges=[rim], round=True,
            tones=([(rim, colour, -40), shade(1.1, _shade(colour))]
                   if colour else
                   [shade(0.72, "light"), shade(1.1, "mid")]))

    def disc(self, x, y, z, r, tone="dark"):
        """A level circle on a surface: a well, a hole, an opening, a drop.
        tone is light, mid, dark, a colour or None."""
        rim = np.array([view(x + r * np.cos(t), y, z + r * np.sin(t))
                        for t in np.linspace(0, 2 * np.pi, 49)])
        return self._add(edges=[rim], thin=True, round=True,
                         tones=[(rim, tone, 0)] if tone else [])

    def face(self, points, tone="dark", outline=False):
        """A flat shape through the given points: a slot, a window, a
        sheet, a section. Its edge is thin unless it is a part of its own
        (outline=True). More than eight points are taken as a curve, and
        its edge is drawn as one line without corners."""
        p = list(self._flat(points))
        edges = [[a, b] for a, b in zip(p, p[1:] + p[:1])]
        return self._add(
            edges=edges if len(p) < 9 else [np.array(p + p[:1])],
            tones=[(p, tone, 90)] if tone else [], thin=not outline)

    def flow(self, points, behind=False):
        """Something that moves: a pencil-blue line through the points
        that ends in an open arrowhead. It is drawn over the parts, or with
        behind=True in its place among them, so that the parts added after
        it hide it: for a sample that goes into a slot, say."""
        row = _spline(self._flat(points))
        if behind:
            self.parts.append(dict(flow=row))
        else:
            self.flows.append(row)

    def label(self, text, at, to=None, ha="left", note=None):
        """A short label at one point, with a leader line to another. A
        sketch writes it in capitals. Give note to word it yourself where
        the case matters, as in "5 mL" or "10 µm". Handwriting is wider
        than type: if notes collide, empty scene.labels after draw() and
        place them again."""
        to = to if to is None else self._flat([to])[0]
        self.labels.append((text, self._flat([at])[0], to, ha,
                            note or text.upper()))

    def block(self, at, to):
        """Sketches only: one mustard block behind the parts, between two
        opposite corners. Keep it under a third of the area the scene
        covers, and clear of the lettering of the notes."""
        (x0, y0), (x1, y1) = self._flat([at, to])
        self.blocks.append([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])

    def draw(self, ax, passes=5):
        """The line drawing: flat greys, 1 pt silhouettes, 0.5 pt edges,
        0.75 pt flows and typeset labels. Give the axes of a figure panel
        and call fit() afterwards. Or give a slide made with slide(): the
        drawing then fills fig.box from its top left corner. For the steps
        of a build, give passes=1 to 4: construct (guides alone), outline,
        tone, colour, labels (with the flows). They keep their place.
        Returns the axes."""
        if not hasattr(ax, "plot"):           # a slide
            fig, ax = ax, drawing(ax)
            self.draw(ax, passes)
            fit(ax, *fig.box)
            return ax
        paper, later = ax.figure.get_facecolor(), []

        def step(level, *artists):            # drawn now, seen from its pass
            if level > passes:
                later.extend(artists)

        def arrow(row, z):
            step(5, *ax.plot(*row.T, color=PENCIL, lw=0.75, zorder=z))
            step(5, ax.annotate("", row[-1], xytext=row[-3], zorder=z,
                                arrowprops={"arrowstyle": "->", "lw": 0.75,
                                            "color": PENCIL, "shrinkA": 0,
                                            "shrinkB": 0,
                                            "mutation_scale": 7}))

        for run, heavy in self._objects():    # back to front
            if "flow" in run[0]:
                arrow(run[0]["flow"], 1)
                continue
            if heavy:                         # the silhouette: drawn wide,
                for part in run:              # then covered by the fills
                    if not part.get("thin"):  # inside, so 1 pt shows
                        step(2, ax.add_patch(plt.Polygon(
                            np.vstack(part["edges"]), fill=False, ec=INK,
                            lw=1.5, joinstyle="round")))
            for part in run:
                ax.add_patch(plt.Polygon(np.vstack(part["edges"]), fc=paper,
                                         lw=0))
                for shape, tone, _ in part["tones"]:
                    step(3 if tone in TONES else 4, ax.add_patch(plt.Polygon(
                        shape, lw=0, fc=TONES.get(tone, tone))))
            for part in run:
                inside = part.get("inside", [])
                if passes == 1:               # guides go with the outline
                    for row in (part.get("guides", part.get("hidden", []))
                                + part["edges"] + inside):
                        ax.plot(*np.transpose(row), color=PENCIL, lw=0.35,
                                alpha=0.8, zorder=1)
                for row in inside:
                    step(2, *ax.plot(*np.transpose(row), color=INK, lw=0.5,
                                     solid_capstyle="round", zorder=1))
                step(2, ax.add_patch(plt.Polygon(
                    np.vstack(part["edges"]), fill=False, ec=INK, lw=0.5,
                    joinstyle="round")))
        for row in self.flows:
            arrow(row, 3)
        for text, at, to, ha, _ in self.labels:
            if to is None:
                step(5, ax.text(*at, text, ha=ha, va="center", zorder=4))
            else:
                step(5, *leader(ax, text, at, to, ha))
        for artist in later:
            artist.set_alpha(0)
            if getattr(artist, "arrow_patch", None) is not None:
                artist.arrow_patch.set_alpha(0)
        ax.autoscale_view()
        return ax

    def sketch(self, fig, font, passes=5, box=None):
        """The marker sketch, for a slide made with slide(): construct
        (pencil), outline (pen), tone (grey marker), colour, labels (notes
        and flows). Give passes=1 to 4 for the steps of a build, as for
        draw(). font is the path of a handwriting face such as
        Mynerve-Regular.ttf. The sketch fills fig.box, or the box you give,
        from its top left corner. Returns the axes."""
        box = box or fig.box
        hand = FontProperties(fname=font, size=9.75)         # 52 px
        rows = [e for p in self.parts for e in p.get("edges") or [p["flow"]]]
        u = np.ptp(np.vstack(rows + self.flows)[:, 0]) / (
            (box[1] - box[0]) * fig.get_figwidth() * 72)
        for last in (0, 0, 0, 0, 1):  # draw, fit, draw again at that size
            ax = drawing(fig)
            self._sketch(ax, u, hand, passes)
            fit(ax, *box)
            now = fig.dpi / 72 / abs(ax.transData.get_matrix()[0, 0])
            if last or abs(now / u - 1) < 0.02:
                return ax
            u = now
            ax.remove()

    def _sketch(self, ax, u, hand, passes):
        rng = np.random.default_rng(self.seed)
        paper = to_rgb(ax.figure.get_facecolor())

        def lay(shapes, colour, level, clip=None, **more):
            """One layer of strokes. A layer above the pass asked for is
            laid down unseen, so that a build keeps its place."""
            if len(shapes):
                layer = ax.add_collection(PolyCollection(
                    shapes, facecolors=colour, edgecolors="none", **more))
                if clip is not None:
                    layer.set_clip_path(plt.Polygon(
                        clip, transform=ax.transData))
                if level > passes:
                    layer.set_alpha(0)

        def tone(shape, name, angle, nib=5.8, clip=None):
            grey = name in COVER
            ink = to_rgb(MARKER if grey else name)
            cover = COVER.get(name, 0.86)
            strips, streaks = _marker(shape, angle, rng, u, nib,
                                      0.45 if clip is None else -1.5)
            lay(strips, [(*ink, cover * rng.uniform(0.86, 1)) for _ in strips],
                3 if grey else 4, clip)
            lay(streaks, [(*paper, rng.uniform(0.08, 0.24)) for _ in streaks],
                3 if grey else 4, shape)

        def arrow(row, z):
            way = np.arctan2(*(row[-1] - row[-3])[::-1])
            barbs = [[row[-1] - 2.4 * u * np.array([np.cos(a), np.sin(a)]),
                      row[-1]] for a in (way + 0.42, way - 0.42)]
            lay([_pen(_hand(row, rng, u), 0.62, rng, u, (0.7, 0.5))]
                + [_pen(_hand(barb, rng, u, 0.2), 0.62, rng, u, (0.7, 0.7))
                   for barb in barbs], PENCIL, 5, zorder=z)

        for shape in self.blocks:
            tone(shape, MUSTARD, 90, nib=11)
        for run, heavy in self._objects():
            if "flow" in run[0]:
                arrow(run[0]["flow"], 1)
                continue
            if heavy:                         # the silhouette, wide; the
                lay([_pen(_hand(row, rng, u), 1.3, rng, u)    # fills cover
                     for part in run if not part.get("thin")  # its inside
                     for row in part["edges"]], INK, 2)
            for part in run:
                inside = part.get("inside", [])
                ax.add_patch(plt.Polygon(np.vstack(part["edges"]), fc=paper,
                                         lw=0))
                for shape, name, angle in part["tones"]:
                    tone(shape, name, angle,
                         clip=shape if part.get("round") else None)
                lay([_pen(_hand(row, rng, u, 5), 0.34, rng, u, (0.8, 0.6))
                     for row in part.get("hidden", []) + part["edges"]
                     + inside], PENCIL, 1, alpha=0.8)
            for part in run:                  # the pen, gone over twice
                twice = heavy and not part.get("thin")
                lay([_pen(_hand(row, rng, u), w, rng, u)
                     for row in part.get("inside", []) + part["edges"]
                     for w in ((0.48, 0.36) if twice else (0.48,))], INK, 2)
        for row in self.flows:
            arrow(row, 3)
        seen = float(passes >= 5)
        for _, at, to, ha, note in self.labels:
            bend = rng.uniform(0.08, 0.16) * rng.choice([-1, 1])
            line = {"arrowstyle": "-", "lw": 0.4, "color": INK, "shrinkA": 2,
                    "shrinkB": 0, "alpha": seen,
                    "connectionstyle": f"arc3,rad={bend:.2f}"}
            ax.annotate(note, at if to is None else to, xytext=at, ha=ha,
                        va="center", zorder=4, fontproperties=hand,
                        alpha=seen, rotation=rng.uniform(-2, 2),
                        linespacing=1, annotation_clip=False,
                        arrowprops=None if to is None else line)
            if to is not None:
                ax.add_patch(plt.Circle(to, 0.55 * u, fc=INK, lw=0, zorder=4,
                                        alpha=seen))
