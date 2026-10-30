# The ruled half: settings and helpers for charts, the layout grid and
# slides. sketched.py, beside this file, draws the explanations. Tested
# with Matplotlib 3.11. A chart on the right-hand columns of a figure:
#
#   from ruled import *
#
#   fig = plt.figure(figsize=(183 * MM, 60 * MM))
#   ax = fig.add_subplot()
#   line, = ax.plot(x, mean, "o-", color=PERSIMMON)
#   replicates(ax, x, runs, line)
#   label_series(ax, line, "Mutant")
#   reference(ax, [41.2, 47.9], [50, 50], "6.7 °C", INK, ha="left")
#   ruled(ax, xticks, yticks, "Temperature (°C)", "Activity (%)")
#   fit(ax, *columns(7, 12), top=1, bottom=0, letter="b")
#   fig.savefig("figure.pdf")
import re
import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation

INK, GREY = "#1B1A17", "#55524A"
PERSIMMON, PETROL, MUSTARD = "#D5501E", "#0079A6", "#E0A526"
PENCIL, PAPER = "#6C9DB8", "#FBFAF6"
MARKER_GREYS = "#DEE1E3", "#BEC3C7", "#939DA3"      # light, mid, dark
MM = 1 / 25.4                      # inches per millimetre

def _truetype(families):
    """Leave out the installed faces with CFF outlines, as in most .otf
    files: a PDF embeds them as TrueType, and the file comes out flawed.
    Families that are not installed stay, so an SVG still asks for them."""
    keep = []
    for family in families:
        try:
            path = mpl.font_manager.findfont(
                mpl.font_manager.FontProperties(family=family),
                fallback_to_default=False)
        except ValueError:             # not installed here
            keep.append(family)
            continue
        with open(path, "rb") as fh:
            if fh.read(4) != b"OTTO":
                keep.append(family)
    return keep


mpl.rcParams.update({
    # one grotesque at 7 pt, kept as live text in the PDF and SVG
    "font.family": "sans-serif",
    "font.sans-serif": _truetype(["Helvetica", "Arial", "Nimbus Sans",
                                  "TeX Gyre Heros", "Liberation Sans",
                                  "Arimo", "DejaVu Sans"]),
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    # indices in units, as in m$^{-2}$, stay in the same face
    "mathtext.fontset": "custom", "mathtext.rm": "sans",
    "mathtext.it": "sans:italic", "mathtext.bf": "sans:bold",
    "mathtext.cal": "sans",
    # ruled: hairline axes, ticks pointing out, no box, no grid
    "axes.linewidth": 0.5, "axes.edgecolor": INK, "axes.labelcolor": INK,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.color": INK, "ytick.color": INK, "text.color": INK,
    # a fallback order only: give each series its colour by role
    "axes.prop_cycle": mpl.cycler(color=[GREY, PERSIMMON, PETROL]),
    "lines.linewidth": 1.0, "lines.markersize": 3.6,
    "lines.markeredgecolor": "white", "lines.markeredgewidth": 0.7,
    "legend.frameon": False,
})


def typeface():
    """The family Matplotlib measures with and embeds in a PDF or PNG. It
    falls back without a word, so report this when it is not Helvetica or
    Arial. The list skips faces that would embed badly in a PDF. An SVG
    embeds no font: each viewer sets it in the first family it has."""
    found = mpl.font_manager.findfont(mpl.font_manager.FontProperties())
    return mpl.font_manager.get_font(found).family_name


def ruled(ax, xticks, yticks, xlabel, ylabel, xnames=None, xpad=0):
    """Call after plotting. Trim both axes to their ticks, detach them and
    set both titles horizontal. For categories, pass their names as xnames
    and half a step as xpad, and None as xlabel if the names say it all."""
    ax.set_xlim(xticks[0] - xpad, xticks[-1] + xpad)
    ax.set_ylim(yticks[0], yticks[-1])
    ax.set_xticks(xticks, xnames)
    ax.set_yticks(yticks)
    out = 4 + getattr(ax, "fan", 0)    # the y axis clears replicate dots
    for side, ticks, off in (("bottom", xticks, 4), ("left", yticks, out)):
        ax.spines[side].set_bounds(ticks[0], ticks[-1])
        ax.spines[side].set_position(("outward", off))
    for mark in (*ax.lines, *ax.collections):
        mark.set_clip_on(False)        # points on the end ticks stay whole
    if xlabel:
        ax.set_xlabel(xlabel, loc="right")
    fig = ax.figure
    fig.draw_without_rendering()       # measure the tick labels, so that
    wide = max(t.get_window_extent().width for t in ax.get_yticklabels())
    left = wide * 72 / fig.dpi + out + 6     # the y title starts flush
    up = ScaledTranslation(-left / 72, 8 / 72, fig.dpi_scale_trans)
    ax.text(0, 1, ylabel, transform=ax.transAxes + up, ha="left", va="bottom")


def replicates(ax, x, runs, line):
    """Every replicate as a 1.6 pt dot in the series colour, set to the
    left and right of its mean so that the marker never hides one.
    runs holds one sequence per replicate, each as long as x."""
    half = (len(runs) + 1) // 2
    step = min(1.4, 4 / max(1, half - 1))
    ax.fan = 3.5 + (half - 1) * step + 0.8           # how far the dots reach
    for i, run in enumerate(runs):
        side, k = (-1, half - 1 - i) if i < half else (1, i - half)
        dx = side * (3.5 + k * step) / 72
        shift = ScaledTranslation(dx, 0, ax.figure.dpi_scale_trans)
        ax.plot(x, run, "o", ms=1.6, mec="none", alpha=0.8,
                color=line.get_color(), zorder=line.get_zorder() - 1,
                transform=ax.transData + shift)


def label_series(ax, line, text, at=-1, side="right", dy=0):
    """Name a series in ink beside one of its points: the last one unless
    names would collide there. dy lifts the name, in points."""
    x, y = line.get_xdata()[at], line.get_ydata()[at]
    gap = max(9, getattr(ax, "fan", 0) + 3)          # clear of replicate dots
    ax.annotate(text, (x, y), xytext=(gap if side == "right" else -gap, dy),
                textcoords="offset points", va="center",
                ha="left" if side == "right" else "right")


def reference(ax, xs, ys, text, colour=MARKER_GREYS[2], ha="center"):
    """A 0.5 pt reference line with a short label above it. On a level
    line ha puts the label over the "left" end, the "center" or the
    "right" end: choose the place that is clear of the data. On any other
    line the label sits above the last point and ha says which way it
    runs. The dark marker grey marks context, INK a value you worked out
    between two series. Call it before ruled()."""
    ax.plot(xs, ys, color=colour, lw=0.5, zorder=1)
    ends = {"left": min(xs), "center": (min(xs) + max(xs)) / 2,
            "right": max(xs)}
    at = (ends[ha], ys[0]) if ys[0] == ys[-1] else (xs[-1], ys[-1])
    ax.annotate(text, at, xytext=(0, 3), textcoords="offset points",
                va="bottom", ha=ha)


def columns(first, last, width=183):
    """Left and right edge of grid columns first to last, as fractions of
    the figure width. The grid has twelve columns with 3 mm gutters, or
    six in a figure under 100 mm wide."""
    n = 12 if width > 100 else 6
    col = (width - 3 * (n - 1)) / n
    return ((first - 1) * (col + 3) / width,
            ((last - 1) * (col + 3) + col) / width)


def fit(ax, left, right, top, bottom, letter=None):
    """Call last. Move the axes until everything drawn in them, titles and
    names included, fills the box given in figure fractions. letter puts a
    bold panel letter in the top left corner, above the y title. Returns
    the box the content fills, as left, right, top, bottom. A drawing
    keeps its proportions, so it fills either the width or the height of
    the box. If it ends short of the right, make the figure taller, by
    the width you want over the width you got. If it ends short of the
    bottom, make the figure shorter, or give the drawing more columns."""
    fig = ax.figure
    if not ax.axison:
        ax.autoscale_view()            # a drawing's limits hug its shapes,
        for mark in (*ax.patches, *ax.lines, *ax.collections):
            if mark.get_clip_path() is None:
                mark.set_clip_on(False)    # so its outlines must not be cut
            else:                          # to them; a clip of its own stays
                mark.set_clip_box(None)
    if letter:
        if not hasattr(ax, "letter"):
            ax.letter = fig.text(0, 0, letter, fontsize=8, va="top",
                                 fontweight="bold")
        ax.letter.set_position((left, top))
        top -= 11 / 72 / fig.get_figheight()
    # half the heaviest outline hangs outside the shapes that are measured
    edge = max([p.get_linewidth() for p in ax.patches], default=0) / 144
    ex, ey = edge / fig.get_figwidth(), edge / fig.get_figheight()

    def measure():
        box = ax.get_tightbbox().transformed(fig.transFigure.inverted())
        return box.x0 - ex, box.x1 + ex, box.y1 + ey, box.y0 - ey

    seen = None
    for _ in range(24):                # a drawing settles in a few rounds
        x0, x1, y1, y0 = box = measure()
        if seen and max(abs(u - v) for u, v in zip(box, seen)) < 1e-7:
            break
        seen, now = box, ax.get_position()
        ax.set_position([now.x0 + left - x0, now.y0 + bottom - y0,
                         now.width + (right - left) - (x1 - x0),
                         now.height + (top - bottom) - (y1 - y0)])
    return measure()


def slide(headline, chart=True):
    """A 1920 x 1080 px slide with 96 px margins. The canvas is 5 x 2.8125
    in at 384 dpi, so every size above grows 5.33 times and 7 pt labels
    are 37 px. The headline is wrapped to fit between the margins and may
    take two lines. A number stays with its unit; put a line break in the
    headline to choose the place yourself. fig.box is the area under it.
    With a chart, finish with fit(ax, *fig.box), or with a narrower box if
    the full width would stretch the chart: (0.05, 0.68, *fig.box[2:])
    gives about 1200 px. With chart=False the area is left free for a
    drawing, and ax is None."""
    fig = plt.figure(figsize=(5, 2.8125), dpi=384, facecolor=PAPER)
    text = fig.text(0.05, 1 - 96 / 1080, "", fontsize=13.5, va="top",
                    fontweight="bold", linespacing=1.15)     # 72 px
    headline = re.sub(r"(?<=\d) (?=[^\W\d]|°|%)", " ", headline)
    text.set_text(headline)
    if text.get_window_extent().width > 1728 and "\n" not in headline:
        for count in range(16, len(headline)):       # two lines, balanced
            text.set_text(textwrap.fill(headline, count))
            if text.get_text().count("\n") < 2:
                break
    if (text.get_window_extent().width > 1728
            or text.get_text().count("\n") > 1):
        raise ValueError("the headline does not fit on two lines")
    fig.box = (0.05, 0.95, 0.7, 96 / 1080)           # left, right, top, bottom
    ax = fig.add_axes([0.12, 0.25, 0.7, 0.4], facecolor="none")
    if not chart:
        ax.remove()
    return fig, ax if chart else None
