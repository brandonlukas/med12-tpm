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
    plot(d, g.tpm[d.id]);
    table(d, g.de[d.id]);
  }
}

// deterministic jitter so points do not move between genes
const jitter = (k) => ((Math.abs(Math.sin(k * 12.9898) * 43758.5453) % 1) - 0.5) * 0.3;

function plot(d, tpm) {
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
  Plotly.react(`plot-${d.id}`, traces, {
    margin: { l: 52, r: 8, t: 8, b: 48 },
    font: { family: "Arial, sans-serif", size: 13, color: "#1a1a1a" },
    xaxis: { tickvals: d.groups.map((_, k) => k), ticktext: d.groups.map(groupLabel),
      range: [-0.5, d.groups.length - 0.5], tickangle: 0, showgrid: false, zeroline: false, fixedrange: true },
    yaxis: { title: { text: log ? "log₂(TPM + 1)" : "TPM" }, rangemode: "tozero", gridcolor: "#eeeeee",
      zeroline: false, fixedrange: true, ticks: "outside", tickcolor: "#bbb" },
    plot_bgcolor: "#fff", paper_bgcolor: "#fff", hoverlabel: { bgcolor: "#fff", font: { color: "#1a1a1a" } },
  }, { displayModeBar: false, responsive: true });
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
