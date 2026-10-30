// Shared by the explorer and the DE page. Data layout: see build_site_data.py.
const COLORS = {           // reference groups grey; mutants vermilion across datasets
  MM: "#6B7F99", WT: "#6B7F99",
  LM_wtMED12: "#0072B2", LM_mutMED12: "#D55E00",
  KO: "#0072B2", Mut1: "#D55E00", Mut2: "#CC79A7",
};
const GROUP_LABEL = { MM: "MM", LM_wtMED12: "LM<br>MED12 WT", LM_mutMED12: "LM<br>MED12 mut" };
const CONTRAST_LABEL = {
  LM_vs_MM: "LM vs MM (all 10 pairs)",
  LMmut_vs_MM: "LM MED12 mut vs MM",
  LMwt_vs_MM: "LM MED12 WT vs MM",
  LMmut_vs_LMwt: "LM mut vs LM WT (shift vs own MM)",
};
const SHORT_DS = { JCI_tissue: "Tissue", NC_2D: "2D cells", NC_3D_spheroid: "3D spheroids" };

const $ = (id) => document.getElementById(id);
const fmt = (x, d = 2) => (x == null ? "—" : (Math.abs(x) < 0.5 * 10 ** -d ? 0 : x).toFixed(d));
const fmtP = (p) => (p == null ? "—" : p < 1e-3 ? p.toExponential(1) : p.toFixed(3));
const fmtTPM = (v) => (v == null ? "—" : v >= 100 ? Math.round(v).toLocaleString() : v >= 10 ? v.toFixed(1) : v.toFixed(2));
const groupLabel = (g, br = true) => (GROUP_LABEL[g] || g).replace("<br>", br ? "<br>" : " ");
const contrastLabel = (c) => CONTRAST_LABEL[c] || c.replace("_vs_", " vs ");
const esc = (s) => String(s).replace(/[&<>"]/g, (ch) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[ch]);
const explorerLink = (names) => `index.html#genes=${names.map(encodeURIComponent).join(",")}`;

const PLOT_FONT = { family: "Arial, sans-serif", size: 13, color: "#1a1a1a" };
const PLOT_CONFIG = { displayModeBar: false, responsive: true };

// genes: [[id, name, type], ...]; lookups by symbol (may repeat) and by unversioned id
const Data = {
  meta: null, genes: null, byName: new Map(), byId: new Map(), chunks: new Map(),
  async load() {
    [this.meta, this.genes] = await Promise.all([
      fetch("data/meta.json").then((r) => r.json()),
      fetch("data/genes.json").then((r) => r.json()),
    ]);
    this.genes.forEach(([id, name], i) => {
      this.byId.set(id.split(".")[0].toUpperCase(), i);
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
  // up to n genes whose symbol starts with q, shortest first
  suggest(q, n = 12) {
    q = q.trim().toUpperCase();
    if (!q) return [];
    const hits = [];
    for (let i = 0; i < this.genes.length && hits.length < 200; i++) {
      if (this.genes[i][1].toUpperCase().startsWith(q)) hits.push(i);
    }
    if (!hits.length) return this.lookup(q);
    return hits.sort((a, b) => this.genes[a][1].length - this.genes[b][1].length).slice(0, n);
  },
  async gene(i) {
    const c = Math.floor(i / this.meta.chunk);
    if (!this.chunks.has(c)) this.chunks.set(c, fetch(`data/chunks/${c}.json`).then((r) => r.json()));
    return (await this.chunks.get(c))[this.genes[i][0]];
  },
};

// Gene search box with keyboard-navigable suggestions; onPick(i) gets a gene index.
function attachSearch(input, list, onPick, onMiss) {
  let sel = -1, matches = [];
  const render = () => {
    list.innerHTML = matches.map((i, k) =>
      `<li role="option" data-i="${i}" aria-selected="${k === sel}"><span>${esc(Data.genes[i][1])}</span><span class="id">${Data.genes[i][0]}</span></li>`).join("");
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
      const hit = sel >= 0 ? matches[sel] : Data.lookup(input.value)[0];
      list.hidden = true; matches = []; sel = -1;
      if (hit != null) onPick(hit); else onMiss?.(input.value);
    } else if (e.key === "Escape") list.hidden = true;
  });
  list.addEventListener("mousedown", (e) => {
    const li = e.target.closest("li");
    if (li) { e.preventDefault(); list.hidden = true; matches = []; onPick(+li.dataset.i); }
  });
  input.addEventListener("blur", () => setTimeout(() => (list.hidden = true), 100));
}
