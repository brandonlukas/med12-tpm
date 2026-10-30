#!/usr/bin/env python3
"""Draw the per-dataset PCA (docs/data/pca.json from pca.R) as three small SVG panels and
write them into docs/methods.html between <!-- pca:start --> and <!-- pca:end -->.

    python3 pca_svg.py
"""
import json
import re
from pathlib import Path

here = Path(__file__).resolve().parent
pca = json.load(open(here / "docs/data/pca.json"))
COLORS = {"MM": "#8A8F9A", "WT": "#8A8F9A", "LM_wtMED12": "#0072B2", "LM_mutMED12": "#D55E00", "KO": "#1c1e26", "Mut1": "#D55E00", "Mut2": "#BF6A9B"}
LABEL = {"MM": "MM", "LM_wtMED12": "LM WT", "LM_mutMED12": "LM mut"}
TITLE = {"JCI_tissue": "Tissue", "NC_2D": "2D cells", "NC_3D_spheroid": "3D spheroids"}
ORDER = {"JCI_tissue": ["MM", "LM_wtMED12", "LM_mutMED12"], "NC_2D": ["WT", "KO", "Mut1", "Mut2"], "NC_3D_spheroid": ["WT", "Mut1", "Mut2"]}
INK, MUTED, RULE = "#1c1e26", "#5c6170", "#e4e6ec"
W, H, L, R, T, B = 320, 290, 46, 12, 12, 44


def num(v):
    """A coordinate with a true minus sign, as elsewhere on the site."""
    return f"{v:.1f}".replace("-", "\u2212")


def sample_label(s):
    m = re.match(r"tissue_(\d+)_(LM|MM)$", s)
    if m: return f"Patient {m.group(1)}, {m.group(2)}"
    m = re.match(r"cell[23]D_(\w+?)_rep(\d+)$", s)
    return f"{m.group(1)}, replicate {m.group(2)}" if m else s


def panel(ds):
    d = pca[ds]; pts = d["samples"]; ve = d["explained"]
    xs, ys = [p["pc1"] for p in pts], [p["pc2"] for p in pts]
    span = max(max(xs) - min(xs), max(ys) - min(ys)) * 1.18   # equal scaling on both axes
    cx, cy = (max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2
    pw, ph = W - L - R, H - T - B
    k = min(pw, ph) / span
    X = lambda v: L + pw / 2 + (v - cx) * k
    Y = lambda v: T + ph / 2 - (v - cy) * k
    out = [f'<rect x="{L}" y="{T}" width="{pw}" height="{ph}" fill="none" stroke="{RULE}"/>',
           f'<line x1="{X(0):.1f}" y1="{T}" x2="{X(0):.1f}" y2="{T + ph}" stroke="{RULE}" stroke-dasharray="3 3"/>',
           f'<line x1="{L}" y1="{Y(0):.1f}" x2="{L + pw}" y2="{Y(0):.1f}" stroke="{RULE}" stroke-dasharray="3 3"/>']
    patient = lambda p: (re.match(r"tissue_(\d+)_", p["sample"]) or [None, ""])[1]
    if ds == "JCI_tissue":   # join each patient's two samples
        by = {}
        for p in pts: by.setdefault(patient(p), []).append(p)
        for pid, (a, b) in ((k, v) for k, v in by.items() if len(v) == 2):
            out.append(f'<line class="pair" data-patient="{pid}" x1="{X(a["pc1"]):.1f}" y1="{Y(a["pc2"]):.1f}" x2="{X(b["pc1"]):.1f}" y2="{Y(b["pc2"]):.1f}" stroke="#cfd3dc" stroke-width="1"/>')
    for p in pts:
        out.append(f'<circle cx="{X(p["pc1"]):.1f}" cy="{Y(p["pc2"]):.1f}" r="5.5" fill="{COLORS[p["group"]]}" stroke="#fff" stroke-width="1.2" '
                   f'data-label="{sample_label(p["sample"])}" data-group="{LABEL.get(p["group"], p["group"])}" data-color="{COLORS[p["group"]]}" '
                   f'data-pc="PC1 {num(p["pc1"])} · PC2 {num(p["pc2"])}" data-patient="{patient(p)}" tabindex="0" '
                   f'aria-label="{sample_label(p["sample"])}, {LABEL.get(p["group"], p["group"])}, PC1 {num(p["pc1"])}, PC2 {num(p["pc2"])}"/>')
    out.append(f'<text x="{L + pw / 2:.0f}" y="{H - 14}" text-anchor="middle" font-size="12" fill="{MUTED}">PC1 · {ve[0]:.0f}% of variation</text>')
    out.append(f'<text transform="translate(14 {T + ph / 2:.0f}) rotate(-90)" text-anchor="middle" font-size="12" fill="{MUTED}">PC2 · {ve[1]:.0f}%</text>')
    legend = " ".join(f'<span data-group="{LABEL.get(g, g)}"><span class="swatch" style="background:{COLORS[g]}"></span>{LABEL.get(g, g)}</span>' for g in ORDER[ds])
    return (f'<figure class="pca-panel"><h3>{TITLE[ds]}</h3>'
            f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="{TITLE[ds]}: samples by their first two principal components">{"".join(out)}</svg>'
            f'<p class="legend">{legend}</p></figure>')


html = '<div class="pca">' + "".join(panel(ds) for ds in ("JCI_tissue", "NC_2D", "NC_3D_spheroid")) + "</div>"
page = here / "docs/methods.html"
s = page.read_text()
a, b = s.index("<!-- pca:start -->"), s.index("<!-- pca:end -->")
page.write_text(s[:a] + "<!-- pca:start -->\n" + html + "\n" + s[b:])
print("wrote the PCA panels into methods.html")
