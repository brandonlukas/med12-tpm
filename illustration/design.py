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

# The site's role colours, so the drawing means what the plots mean: the
# MED12-mutant material takes the plots' vermilion, the reference their grey.
MUTANT, REFERENCE = "#D55E00", "#8A8F9A"

scene = Scene()
scene.box(0, 0, 0, 94, 2, 46)                                # the bench

# the sequencer, at the back, with a screen on its front
with scene.object():
    scene.box(48, 2, 4, 26, 14, 14)
    scene.face([(51, 9.5, 18.05), (65, 9.5, 18.05), (65, 13.5, 18.05), (51, 13.5, 18.05)], tone="dark")

# tissue from a patient: a block of myometrium with a leiomyoma on it
scene.box(6, 2, 26, 24, 7, 15, colour=REFERENCE)
scene.ball(22, 12.2, 35, 4.4, colour=MUTANT)

# 2D culture: a dish of edited cells
with scene.object():
    scene.cylinder(47, 2, 34, 7.5, 2.4)
    scene.disc(47, 4.4, 34, 6.6, tone="light")
for dx, dz in ((-2.6, -1.4), (1.9, -2.7), (-0.5, 2.3), (3.3, 1.3)):
    scene.ball(47 + dx, 5.3, 34 + dz, 0.95, colour=MUTANT)

# 3D culture: a well with spheroids of the same cells
with scene.object():
    scene.cylinder(69, 2, 34, 6, 4)
    scene.disc(69, 6, 34, 5.2, tone="dark")
for dx, dz in ((-2.3, 0.3), (1.6, -1.7), (1.5, 2.1)):
    scene.ball(69 + dx, 6.9, 34 + dz, 1.7, colour=MUTANT)

# what moves: every sample into the sequencer; cells on into spheroid culture
scene.flow([(29, 7.5, 27), (51.5, 5.5, 18.6)])
scene.flow([(50, 4.6, 26.6), (58, 5, 18.6)])
scene.flow([(69, 6.5, 27.8), (65.5, 6.5, 18.6)])
scene.flow([(55.6, 4.6, 34), (62.2, 4.6, 34)])

# five labels
scene.label("Leiomyoma", (3, 22, 36), (19.5, 15.2, 34))
scene.label("Myometrium", (4, -5, 54), (16, 5.5, 41))
scene.label("MED12-edited cells", (42, -5, 56), (45, 4.4, 40.5))
scene.label("Spheroids", (80, -3, 50), (71.2, 8.4, 36.5))
scene.label("RNA-seq", (66, 25, 8), (62, 16, 10))

fig = plt.figure(figsize=(89 * MM, 38.6 * MM))
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
