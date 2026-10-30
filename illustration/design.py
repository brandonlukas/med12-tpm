"""The study-design drawing for the overview page, in the house line-drawing
finish (ruled.py and sketched.py beside this file).

    uv run --with matplotlib --with numpy python illustration/design.py

Writes docs/img/design.svg (live text, with a font fallback stack for the web),
illustration/design.pdf for reuse in slides or a paper, and design_check.png.
"""
import re
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.font_manager as fm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
# Arial from the Windows side of WSL, so the typeface is the house one
for name in ("arial.ttf", "arialbd.ttf", "ariali.ttf"):
    path = Path("/mnt/c/Windows/Fonts") / name
    if path.exists():
        fm.fontManager.addfont(str(path))

from ruled import *      # noqa: E402
from sketched import *   # noqa: E402
from sketched import _spline, YAW   # noqa: E402
import numpy as np   # noqa: E402

# The site's role colours, so the drawing means what the plots mean: the
# MED12-mutant material takes the plots' vermilion, the reference their grey.
MUTANT, REFERENCE = "#D55E00", "#8A8F9A"

# A coronal section of the uterus is a flat shape, so it is drawn as faces in the one
# plane that the house view shows without skew (perpendicular to the line of sight):
# a point (u, w) of the section lands on the page at (u, w). The cells and the two
# culture formats are solids in front of it.
E1 = np.array([np.cos(YAW), 0, np.sin(YAW)])
E2 = np.array([-0.187, 0.927, -0.324])
at = lambda u, w: tuple(u * E1 + w * E2)
curve = lambda pts: [at(u, w) for u, w in _spline(pts, 8)]

scene = Scene()
# the uterus: myometrium (the reference tissue, grey), the cavity, a leiomyoma in the wall
wall = [(0, -24), (-6, -22), (-11, -12), (-16, 0), (-21, 12), (-18, 21), (-8, 25.5), (0, 26),
        (8, 25.5), (18, 21), (21, 12), (16, 0), (11, -12), (6, -22), (0, -24)]
cavity = [(0, -7), (-5, 3), (-11, 13), (-5, 15.5), (0, 15), (5, 15.5), (11, 13), (5, 3), (0, -7)]
scene.face(curve(wall), tone=REFERENCE, outline=True)
scene.face(curve(cavity), tone="#FFFFFF")
scene.ball(*at(-14.5, 8), 5.4, colour=MUTANT)

# cells derived from the myometrium and edited at MED12 codon 44
for u, w in ((38.5, 3.6), (43, 1.4), (39.8, -3)):
    scene.ball(*at(u, w), 2.3, colour=MUTANT)

# 2D culture: a flat dish, the cells spread on its floor
dx, dy, dz = at(62, -17)
with scene.object():
    scene.cylinder(dx, dy, dz, 9.5, 2.6)
    scene.disc(dx, dy + 2.6, dz, 8.6, tone="light")
for ou, ow in ((-4.6, -1.8), (-1.2, -4.2), (2.6, -3.1), (-2.4, 2.4), (1.4, 0.6), (4.8, 2.2), (-5.4, 3.4)):
    scene.ball(dx + ou, dy + 3.5, dz + ow, 1.0, colour=MUTANT)

# 3D culture: a well with spheroids of the same cells
wx, wy, wz = at(90, -15)
with scene.object():
    scene.cylinder(wx, wy, wz, 7.2, 5)
    scene.disc(wx, wy + 5, wz, 6.3, tone="dark")
for ou, ow in ((-2.8, 0.4), (2.0, -2.0), (1.8, 2.6)):
    scene.ball(wx + ou, wy + 6.1, wz + ow, 2.1, colour=MUTANT)

# what moves: cells taken from the myometrium, then grown two ways
scene.flow([(21, 1), (28.5, 2.8), (34.6, 2.4)])
scene.flow([(42.3, -4.2), (49, -9.5), (54, -13)])
scene.flow([(45.2, 0.4), (64, -1.5), (83, -8.5)])

# five labels
scene.label("Leiomyoma", (-31, 22), (-15.5, 11.6), ha="right")
scene.label("Myometrium", (-31, -13), (-12.2, -13), ha="right")
scene.label("MED12-edited myometrial cells", (27, 19), (38.6, 5.6))
scene.label("2D culture", (46, -31), (54.5, -20))
scene.label("3D spheroids", (102, 3), (91.6, -6.6))

fig = plt.figure(figsize=(89 * MM, 27 * MM))
ax = drawing(fig)
scene.draw(ax)
box = fit(ax, *columns(1, 6, 89), top=1, bottom=0)
print("typeface:", typeface())
print("content box (left, right, top, bottom):", [round(v, 3) for v in box])

out = HERE.parent / "docs" / "img"
out.mkdir(exist_ok=True)
fig.savefig(out / "design.svg", transparent=True)
fig.savefig(HERE / "design.pdf")
fig.savefig(HERE / "design_check.png", dpi=300)
# the web has no Arial on every machine: give the SVG's text a fallback stack
svg = (out / "design.svg").read_text()
svg = re.sub(r"font-family:\s*[^;\"]+", "font-family: Arial, Helvetica, 'Liberation Sans', sans-serif", svg)
(out / "design.svg").write_text(svg)
print("wrote", out / "design.svg", "and", HERE / "design.pdf")
