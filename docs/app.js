// Gene explorer. Data layout: see build_site_data.py.
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

// groups each bracket joins; X_vs_Y names default to [X, Y]. LM_vs_MM pools both
// LM groups, so it has no single bracket and stays in the table only.
const CONTRAST_GROUPS = {
  LM_vs_MM: null,
  LMmut_vs_MM: ["LM_mutMED12", "MM"],
  LMwt_vs_MM: ["LM_wtMED12", "MM"],
  LMmut_vs_LMwt: ["LM_mutMED12", "LM_wtMED12"],
};
const stars = (p) => (p < 0.001 ? "***" : p < 0.01 ? "**" : p < 0.05 ? "*" : "ns");

const $ = (id) => document.getElementById(id);
let meta, genes, byName = new Map(), byId = new Map();
const chunkCache = new Map();
let current = null;

const fmt = (x, d = 2) => (x == null ? "—" : (Math.abs(x) < 0.5 * 10 ** -d ? 0 : x).toFixed(d));
const fmtP = (p) => (p == null ? "—" : p < 1e-3 ? p.toExponential(1) : p.toFixed(3));
const groupLabel = (g) => GROUP_LABEL[g] || g;

async function init() {
  [meta, genes] = await Promise.all([
    fetch("data/meta.json").then((r) => r.json()),
    fetch("data/genes.json").then((r) => r.json()),
  ]);
  genes.forEach(([id, name], i) => {
    byId.set(id.split(".")[0].toUpperCase(), i);
    const k = name.toUpperCase();
    if (!byName.has(k)) byName.set(k, []);
    byName.get(k).push(i);
  });
  $("panels").innerHTML = meta.datasets.map((d) => `
    <section class="panel">
      <h2>${d.label}</h2>
      <p class="sub">${d.samples.length} samples</p>
      <div class="plot" id="plot-${d.id}"></div>
      <table id="de-${d.id}"></table>
    </section>`).join("");
  const q = new URLSearchParams(location.hash.slice(1)).get("gene") || "MED12";
  show(lookup(q)[0] ?? 0);
}

function lookup(q) {
  q = q.trim().toUpperCase();
  if (!q) return [];
  if (byName.has(q)) return byName.get(q);
  const id = byId.get(q.split(".")[0]);
  return id == null ? [] : [id];
}

async function getGene(i) {
  const c = Math.floor(i / meta.chunk);
  if (!chunkCache.has(c)) chunkCache.set(c, fetch(`data/chunks/${c}.json`).then((r) => r.json()));
  return (await chunkCache.get(c))[genes[i][0]];
}

async function show(i) {
  current = i;
  const [id, name] = genes[i];
  $("gene").value = name;
  $("status").textContent = "Loading…";
  history.replaceState(null, "", `#gene=${encodeURIComponent(name)}`);
  const g = await getGene(i);
  if (current !== i) return;
  $("status").innerHTML = `<b>${name}</b> · ${id}`;
  for (const d of meta.datasets) {
    plot(d, g.tpm[d.id], g.de[d.id]);
    table(d, g.de[d.id]);
  }
}

// deterministic jitter so points do not move between genes
const jitter = (k) => ((Math.abs(Math.sin(k * 12.9898) * 43758.5453) % 1) - 0.5) * 0.3;

function plot(d, tpm, de) {
  const log = document.querySelector('input[name="scale"]:checked').value === "log";
  const y = tpm.map((v) => (v == null ? null : log ? Math.log2(v + 1) : v));
  const xi = Object.fromEntries(d.groups.map((g, k) => [g, k]));
  const xs = d.samples.map((s, k) => xi[s.group] + jitter(k));
  const traces = [];

  // tissue: one grey line per patient, MM -> LM
  const pts = {};
  d.samples.forEach((s, k) => { if (s.patient) (pts[s.patient] ||= []).push(k); });
  for (const ks of Object.values(pts)) {
    if (ks.length !== 2) continue;
    traces.push({ x: ks.map((k) => xs[k]), y: ks.map((k) => y[k]), mode: "lines",
      line: { color: "#c8c8c8", width: 1.5 }, hoverinfo: "skip", showlegend: false });
  }
  for (const grp of d.groups) {
    const ks = d.samples.map((s, k) => k).filter((k) => d.samples[k].group === grp);
    traces.push({
      x: ks.map((k) => xs[k]), y: ks.map((k) => y[k]), mode: "markers", name: groupLabel(grp).replace("<br>", " "),
      marker: { color: COLORS[grp], size: 10, line: { color: "#fff", width: 1.5 } },
      text: ks.map((k) => `${d.samples[k].sample}<br>TPM ${fmt(tpm[k], tpm[k] < 10 ? 2 : 1)}`),
      hovertemplate: "%{text}<extra></extra>", showlegend: false,
    });
  }
  const ymax = Math.max(0, ...y.filter((v) => v != null)) || 1;
  const { shapes, annotations, top } = brackets(d, de, xi, ymax);
  Plotly.react(`plot-${d.id}`, traces, {
    shapes, annotations,
    margin: { l: 52, r: 8, t: 8, b: 48 },
    font: { family: "Arial, sans-serif", size: 13, color: "#1a1a1a" },
    xaxis: { tickvals: d.groups.map((_, k) => k), ticktext: d.groups.map(groupLabel),
      range: [-0.5, d.groups.length - 0.5], tickangle: 0, showgrid: false, zeroline: false, fixedrange: true },
    yaxis: { title: { text: log ? "log₂(TPM + 1)" : "TPM" }, range: [0, top], gridcolor: "#eeeeee",
      zeroline: false, fixedrange: true, ticks: "outside", tickcolor: "#bbb" },
    plot_bgcolor: "#fff", paper_bgcolor: "#fff", hoverlabel: { bgcolor: "#fff", font: { color: "#1a1a1a" } },
  }, { displayModeBar: false, responsive: true });
}

// One bracket per tested comparison, stacked so no two overlap; label = DESeq2 padj
// as stars. Filtered genes (no test) get no bracket.
function brackets(d, de, xi, ymax) {
  const items = [];
  for (const c of d.contrasts) {
    const pair = c in CONTRAST_GROUPS ? CONTRAST_GROUPS[c] : c.split("_vs_");
    const padj = de[c]?.[2];
    if (!pair || padj == null || !(pair[0] in xi) || !(pair[1] in xi)) continue;
    const [lo, hi] = [xi[pair[0]], xi[pair[1]]].sort((a, b) => a - b);
    items.push({ c, lo, hi, padj });
  }
  items.sort((a, b) => a.hi - a.lo - (b.hi - b.lo) || a.lo - b.lo);
  const levels = [];   // per level, the [lo, hi] spans already placed
  for (const it of items) {
    let k = levels.findIndex((spans) => spans.every(([l, h]) => it.hi <= l || it.lo >= h));
    if (k < 0) { k = levels.length; levels.push([]); }
    levels[k].push([it.lo, it.hi]);
    it.level = k;
  }
  const step = 0.13 * ymax, tick = 0.03 * ymax, inset = 0.08;
  const shapes = [], annotations = [];
  for (const it of items) {
    const yb = ymax * 1.08 + it.level * step;
    const [x1, x2] = [it.lo + inset, it.hi - inset];
    shapes.push({ type: "path", xref: "x", yref: "y", path: `M${x1},${yb - tick} L${x1},${yb} L${x2},${yb} L${x2},${yb - tick}`,
      line: { color: "#555", width: 1 } });
    const sig = it.padj < 0.05;
    annotations.push({ x: (it.lo + it.hi) / 2, y: yb, yanchor: "bottom", showarrow: false, text: stars(it.padj),
      font: { size: sig ? 14 : 11, color: sig ? "#1a1a1a" : "#767676" },
      hovertext: `${CONTRAST_LABEL[it.c] || it.c.replace("_vs_", " vs ")}: p<sub>adj</sub> ${fmtP(it.padj)}` });
  }
  return { shapes, annotations, top: ymax * 1.08 + levels.length * step + 0.04 * ymax };
}

function table(d, de) {
  const rows = d.contrasts.map((c) => {
    const [base, lfc, padj] = de[c] || [];
    const label = CONTRAST_LABEL[c] || c.replace(/_vs_/, " vs ");
    return `<tr><td>${label}</td><td>${fmt(lfc)}</td><td class="${padj != null && padj < 0.05 ? "sig" : ""}">${fmtP(padj)}</td><td>${base == null ? "—" : Math.round(base)}</td></tr>`;
  });
  $(`de-${d.id}`).innerHTML = `<tr><th>DESeq2</th><th>log₂FC</th><th>p<sub>adj</sub></th><th>baseMean</th></tr>${rows.join("")}`;
}

// search box with keyboard-navigable suggestions
let sel = -1, matches = [];
function suggest() {
  const q = $("gene").value.trim().toUpperCase();
  const ul = $("suggest");
  matches = [];
  if (q.length >= 1) {
    for (let i = 0; i < genes.length && matches.length < 12; i++) {
      if (genes[i][1].toUpperCase().startsWith(q)) matches.push(i);
    }
    if (!matches.length) matches = lookup(q);
    matches.sort((a, b) => genes[a][1].length - genes[b][1].length);
  }
  sel = matches.length ? 0 : -1;
  ul.innerHTML = matches.map((i, k) =>
    `<li role="option" data-i="${i}" aria-selected="${k === sel}"><span>${genes[i][1]}</span><span class="id">${genes[i][0]}</span></li>`).join("");
  ul.hidden = !matches.length;
}
function pick(i) { $("suggest").hidden = true; show(i); }
$("gene").addEventListener("input", suggest);
$("gene").addEventListener("keydown", (e) => {
  const items = [...$("suggest").children];
  if (e.key === "ArrowDown" || e.key === "ArrowUp") {
    e.preventDefault();
    if (!items.length) return;
    sel = (sel + (e.key === "ArrowDown" ? 1 : -1) + items.length) % items.length;
    items.forEach((li, k) => li.setAttribute("aria-selected", k === sel));
    items[sel].scrollIntoView({ block: "nearest" });
  } else if (e.key === "Enter") {
    const hit = sel >= 0 ? matches[sel] : lookup($("gene").value)[0];
    if (hit != null) pick(hit); else $("status").textContent = "Gene not found.";
  } else if (e.key === "Escape") $("suggest").hidden = true;
});
$("suggest").addEventListener("mousedown", (e) => {
  const li = e.target.closest("li");
  if (li) { e.preventDefault(); pick(+li.dataset.i); }
});
$("gene").addEventListener("blur", () => setTimeout(() => ($("suggest").hidden = true), 100));
document.querySelectorAll('input[name="scale"]').forEach((r) => r.addEventListener("change", () => current != null && show(current)));

init().catch((e) => { $("status").textContent = "Could not load data: " + e.message; });
