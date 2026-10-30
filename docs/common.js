// Shared by the Genes and Comparisons pages. Data layout: see build_site_data.py.

// Group colors (Okabe-Ito). A color is shared across datasets only when the state is the
// same kind of thing: the reference group is gray everywhere, and groups carrying a MED12 G44
// mutation are warm (LM mut and Mut1 orange, Mut2 a deeper pink than Okabe-Ito's so that
// it clears 3:1 on white). States unique to one dataset get their own color: LM without the
// mutation blue, the MED12 knock-out black.
const COLORS = {
  MM: "#8A8F9A", WT: "#8A8F9A",
  LM_wtMED12: "#0072B2", LM_mutMED12: "#D55E00",
  KO: "#1c1e26", Mut1: "#D55E00", Mut2: "#BF6A9B",
};
// higher / lower in a comparison: the RdBu red and blue, so neither hue is a group color
const UP = "#B2182B", DOWN = "#2166AC";
const INK = "#1c1e26", MUTED = "#5c6170", RULE = "#e4e6ec", GRID = "#eff0f4";
const FONT = 'system-ui, -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif';
const PLOT_FONT = { family: FONT, size: 13, color: INK };
const PLOT_CONFIG = { displayModeBar: false, responsive: true, doubleClick: "reset" };
const BASE_LAYOUT = {
  font: PLOT_FONT, plot_bgcolor: "#fff", paper_bgcolor: "#fff",
  hoverlabel: { bgcolor: "#fff", bordercolor: RULE, font: { family: FONT, size: 13, color: INK } },
};
const AXIS = { gridcolor: GRID, zeroline: false, ticks: "outside", tickcolor: "#c3c7d3", linecolor: RULE, fixedrange: true };

const SHORT_DS = { JCI_tissue: "Tissue", NC_2D: "2D cells", NC_3D_spheroid: "3D spheroids" };
const DATASET_TEXT = {
  JCI_tissue: "Leiomyoma (LM) and matched myometrium (MM) from 10 patients; 6 LM are <i>MED12</i> mutant (G44D), 4 wild-type (WT).",
  NC_2D: "An immortalized myometrial cell line, <i>MED12</i>-edited: unedited (WT), knock-out (KO), two G44N clones (Mut1, Mut2); 3 replicates each.",
  NC_3D_spheroid: "The same cells grown as spheroids; WT has 2 replicates.",
};
const GROUP_LABEL = { MM: "MM", LM_wtMED12: "LM WT", LM_mutMED12: "LM mut" };
const GROUP_TICK = {};
const CONTRAST_LABEL = {
  LM_vs_MM: "LM vs MM, all 10 pairs",
  LMmut_vs_MM: "LM mut vs MM",
  LMwt_vs_MM: "LM WT vs MM",
  LMmut_vs_LMwt: "LM mut vs LM WT",
};
// [numerator, denominator] in the words used on the page; log2FC > 0 means higher in the numerator
const CONTRAST_SIDES = {
  LM_vs_MM: ["LM", "MM"], LMmut_vs_MM: ["LM mut", "MM"], LMwt_vs_MM: ["LM WT", "MM"], LMmut_vs_LMwt: ["LM mut", "LM WT"],
  KO_vs_WT: ["KO", "WT"], Mut1_vs_WT: ["Mut1", "WT"], Mut2_vs_WT: ["Mut2", "WT"],
};
// the shortest unambiguous name, for the multi-gene "changes" column
const CONTRAST_SHORT = {
  LM_vs_MM: "LM vs MM", LMmut_vs_MM: "LM mut vs MM", LMwt_vs_MM: "LM WT vs MM", LMmut_vs_LMwt: "LM mut vs LM WT",
  KO_vs_WT: "KO", Mut1_vs_WT: "Mut1", Mut2_vs_WT: "Mut2",
};

const $ = (id) => document.getElementById(id);
const MINUS = "−";
const fmt = (x, d = 2) => (x == null ? "—" : (Math.abs(x) < 0.5 * 10 ** -d ? 0 : x).toFixed(d).replace("-", MINUS));
const fmtP = (p) => (p == null ? "—" : p === 0 ? "< 1e-300" : p < 1e-3 ? p.toExponential(1) : p < 0.1 ? p.toFixed(3) : p.toFixed(2));
const fmtTPM = (v) => (v == null ? "—" : v >= 100 ? Math.round(v).toLocaleString() : v >= 10 ? v.toFixed(1) : v.toFixed(2));
const ordinal = (n) => n + (n % 100 >= 11 && n % 100 <= 13 ? "th" : ["th", "st", "nd", "rd"][n % 10] || "th");
const groupLabel = (g) => GROUP_LABEL[g] || g;
const groupTick = (g) => GROUP_TICK[g] || GROUP_LABEL[g] || g;
const contrastLabel = (c) => CONTRAST_LABEL[c] || c.replace("_vs_", " vs ");
const contrastSides = (c) => CONTRAST_SIDES[c] || c.split("_vs_");
const contrastShort = (c) => CONTRAST_SHORT[c] || c.replace("_vs_", " vs ");
const contrastOption = (ds, c) => `${SHORT_DS[ds] || ds} · ${contrastLabel(c)}`;
const esc = (s) => String(s).replace(/[&<>"]/g, (ch) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[ch]);
const datasetNote = (id) => DATASET_TEXT[id] || "";
const genesLink = (names) => `genes.html#genes=${names.map(encodeURIComponent).join(",")}`;
const shortId = (id) => id.split(".")[0];

function sampleLabel(s) {
  let m = s.match(/^tissue_(\d+)_(LM|MM)$/);
  if (m) return `Patient ${m[1]}, ${m[2]}`;
  m = s.match(/^cell[23]D_(\w+?)_rep(\d+)$/);
  if (m) return `${m[1]}, replicate ${m[2]}`;
  return s;
}

// A y axis in TPM. On the log scale the spacing is log2(TPM + 1) but the ticks read in TPM,
// so every number on the page is in one unit. `top` is in axis units.
function tpmAxis(log, top) {
  const ax = { ...AXIS, title: { text: log ? "TPM (log scale)" : "TPM", standoff: 8 }, range: [-0.025 * top, top] };
  if (!log) return ax;
  let vals = [0, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000].filter((v) => Math.log2(v + 1) <= top);
  if (vals.length > 7) vals = vals.filter((v) => /^(0|10*)$/.test(String(v)));
  return { ...ax, tickvals: vals.map((v) => Math.log2(v + 1)), ticktext: vals.map((v) => v.toLocaleString()) };
}

// genes: [[id, name, type], ...]; lookups by symbol (may repeat) and by unversioned id
const Data = {
  meta: null, genes: null, byName: new Map(), byId: new Map(), chunks: new Map(),
  async load() {
    [this.meta, this.genes] = await Promise.all([
      fetch("data/meta.json").then((r) => r.json()),
      fetch("data/genes.json").then((r) => r.json()),
    ]);
    this.genes.forEach(([id, name], i) => {
      this.byId.set(shortId(id).toUpperCase(), i);
      const k = name.toUpperCase();
      if (!this.byName.has(k)) this.byName.set(k, []);
      this.byName.get(k).push(i);
    });
  },
  lookup(q) {
    q = q.trim().toUpperCase();
    if (!q) return [];
    if (this.byName.has(q)) return this.byName.get(q);
    const i = this.byId.get(q.split(".")[0]);
    return i == null ? [] : [i];
  },
  // up to n genes whose symbol starts with q: exact match first, then shortest
  suggest(q, n = 10) {
    q = q.trim().toUpperCase();
    if (!q) return [];
    const hits = [];
    for (let i = 0; i < this.genes.length && hits.length < 300; i++) {
      if (this.genes[i][1].toUpperCase().startsWith(q)) hits.push(i);
    }
    if (!hits.length) return this.lookup(q);
    const rank = (i) => (this.genes[i][1].toUpperCase() === q ? -1 : this.genes[i][1].length);
    return hits.sort((a, b) => rank(a) - rank(b) || this.genes[a][1].localeCompare(this.genes[b][1])).slice(0, n);
  },
  async gene(i) {
    const c = Math.floor(i / this.meta.chunk);
    if (!this.chunks.has(c)) this.chunks.set(c, fetch(`data/chunks/${c}.json`).then((r) => r.json()));
    return (await this.chunks.get(c))[this.genes[i][0]];
  },
};

// Gene search box with keyboard-navigable suggestions. onPick(i) gets a gene index.
function attachSearch(input, list, onPick, onMiss) {
  let sel = -1, matches = [];
  const render = () => {
    list.innerHTML = matches.map((i, k) => {
      const [id, name, type] = Data.genes[i];
      return `<li role="option" data-i="${i}" aria-selected="${k === sel}"><span>${esc(name)} <span class="dim small">${esc(type.replace(/_/g, " "))}</span></span><span class="id">${shortId(id)}</span></li>`;
    }).join("");
    list.hidden = !matches.length;
  };
  input.addEventListener("input", () => { matches = Data.suggest(input.value); sel = matches.length ? 0 : -1; render(); });
  input.addEventListener("keydown", (e) => {
    if (e.key === "ArrowDown" || e.key === "ArrowUp") {
      e.preventDefault();
      if (!matches.length) return;
      sel = (sel + (e.key === "ArrowDown" ? 1 : -1) + matches.length) % matches.length;
      render();
      list.children[sel]?.scrollIntoView({ block: "nearest" });
    } else if (e.key === "Enter") {
      e.preventDefault();
      const q = input.value;
      const hit = sel >= 0 ? matches[sel] : Data.lookup(q)[0];
      list.hidden = true; matches = []; sel = -1;
      if (hit != null) onPick(hit); else if (q.trim()) onMiss?.(q.trim());
    } else if (e.key === "Escape") list.hidden = true;
  });
  list.addEventListener("mousedown", (e) => {
    const li = e.target.closest("li");
    if (li) { e.preventDefault(); list.hidden = true; matches = []; onPick(+li.dataset.i); }
  });
  input.addEventListener("blur", () => setTimeout(() => (list.hidden = true), 120));
}

// Show a Reset button while a plot is drag-zoomed away from the view it was drawn with
// (el._home, set by the draw function); the button redraws it.
function zoomReset(el, btn, redraw) {
  if (el._zoomBound) return;
  el._zoomBound = true;
  const near = (a, b) => Math.abs(a - b) < 1e-6 * (1 + Math.abs(b));
  el.on("plotly_relayout", () => {
    const fx = el._fullLayout?.xaxis.range, fy = el._fullLayout?.yaxis.range, h = el._home;
    btn.hidden = !h || !fx || (near(fx[0], h.x[0]) && near(fx[1], h.x[1]) && near(fy[0], h.y[0]) && near(fy[1], h.y[1]));
  });
  btn.addEventListener("click", () => { btn.hidden = true; redraw(); });
}

// On a plot whose points can be clicked, the cursor becomes a hand only while it is over a point;
// elsewhere it keeps the plot's own cursor (the crosshair that means "drag to zoom").
function pointCursor(el) {
  if (el._cursorBound) return;
  el._cursorBound = true;
  const set = (c) => el.querySelectorAll(".drag").forEach((d) => { d.style.cursor = c; });
  el.on("plotly_hover", () => set("pointer"));
  el.on("plotly_unhover", () => set(""));
}

// Plotly is a 1 MB file; on a poor connection it can fail while the page loads. Try once more,
// and if it is still missing the pages show their tables without the plots.
const NO_PLOT = `<p class="noplot">The chart could not be drawn because the plotting library did not load.
  The numbers are complete; reload the page to try again.</p>`;
function ensurePlotly(src) {
  if (typeof Plotly !== "undefined") return Promise.resolve(true);
  return new Promise((resolve) => {
    const el = document.createElement("script");
    el.src = `${src}?retry=${Date.now()}`;
    el.onload = () => resolve(typeof Plotly !== "undefined");
    el.onerror = () => resolve(false);
    document.head.appendChild(el);
  });
}

// Swap a button's label for a moment, e.g. "Copy link" -> "Copied".
function flash(btn, text, ms = 1600) {
  const old = btn.textContent;
  btn.textContent = text;
  btn.disabled = true;
  setTimeout(() => { btn.textContent = old; btn.disabled = false; }, ms);
}

async function copyText(text, btn) {
  try {
    await navigator.clipboard.writeText(text);
    flash(btn, "Copied");
  } catch {
    window.prompt("Copy this link:", text);
  }
}

function saveText(text, filename, type = "text/csv") {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([text], { type }));
  a.download = filename;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}
