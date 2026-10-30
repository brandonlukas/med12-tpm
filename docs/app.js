// Gene explorer: one gene -> per-dataset plots with DESeq2 brackets; several genes ->
// TPM by gene, so genes can be compared with each other (what TPM is for).
const MAX_GENES = 12;
const PRESETS = {
  "Mediator kinase module": ["MED12", "MED12L", "MED13", "MED13L", "CCNC", "CDK8", "CDK19"],
  "Leiomyoma drivers": ["MED12", "HMGA2", "HMGA1", "FH", "COL4A5", "COL4A6", "RAD51B", "PLAG1"],
  "Hormone receptors": ["ESR1", "ESR2", "PGR", "AR", "GPER1"],
  "Housekeeping": ["ACTB", "GAPDH", "B2M", "TBP", "HPRT1", "PPIA"],
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

let selected = [];   // gene indices, in the order added
let renderId = 0;

const isLog = () => document.querySelector('input[name="scale"]:checked').value === "log";
const toY = (v, log) => (v == null ? null : log ? Math.log2(v + 1) : v);
// deterministic jitter in [-0.5, 0.5] so points do not move between genes
const jitter = (k) => ((Math.abs(Math.sin(k * 12.9898) * 43758.5453) % 1) - 0.5);

async function init() {
  await Data.load();
  $("presets").innerHTML = Object.keys(PRESETS).map((p) => `<button type="button" data-preset="${p}">${p}</button>`).join("");
  $("presets").addEventListener("click", (e) => {
    const p = e.target.dataset.preset;
    if (p) setGenes(PRESETS[p].map((n) => Data.lookup(n)[0]).filter((i) => i != null));
  });
  attachSearch($("gene"), $("suggest"), (i) => { $("gene").value = ""; addGene(i); },
    (q) => ($("status").textContent = `No gene called “${q}”.`));
  $("gene").addEventListener("keydown", (e) => {
    if (e.key === "Backspace" && !$("gene").value && selected.length) setGenes(selected.slice(0, -1));
  });
  $("chips").addEventListener("click", (e) => {
    const i = e.target.closest("[data-remove]")?.dataset.remove;
    if (i != null) setGenes(selected.filter((g) => g !== +i));
  });
  document.querySelectorAll('input[name="scale"]').forEach((r) => r.addEventListener("change", render));
  window.addEventListener("hashchange", fromHash);
  fromHash();
}

function fromHash() {
  const h = new URLSearchParams(location.hash.slice(1));
  const names = (h.get("genes") || h.get("gene") || "MED12").split(",");
  const idx = names.map((n) => Data.lookup(n)[0]).filter((i) => i != null);
  setGenes(idx.length ? idx : [Data.lookup("MED12")[0]], false);
}

function addGene(i) {
  if (selected.includes(i)) return;
  if (selected.length >= MAX_GENES) { $("status").textContent = `Up to ${MAX_GENES} genes at a time.`; return; }
  setGenes([...selected, i]);
}

function setGenes(idx, updateHash = true) {
  selected = [...new Set(idx)].slice(0, MAX_GENES);
  $("chips").innerHTML = selected.map((i) =>
    `<span class="chip">${esc(Data.genes[i][1])}<button type="button" data-remove="${i}" aria-label="Remove ${esc(Data.genes[i][1])}">×</button></span>`).join("");
  $("gene").placeholder = selected.length ? "Add a gene" : "Symbol or Ensembl ID, e.g. MED12";
  if (updateHash) history.replaceState(null, "", "#genes=" + selected.map((i) => encodeURIComponent(Data.genes[i][1])).join(","));
  render();
}

async function render() {
  const id = ++renderId;
  if (!selected.length) { $("panels").innerHTML = ""; $("status").textContent = "Add a gene to start."; return; }
  $("status").textContent = "Loading…";
  const data = await Promise.all(selected.map((i) => Data.gene(i)));
  if (id !== renderId) return;
  const multi = selected.length > 1;
  const [gid, gname, gtype] = Data.genes[selected[0]];
  $("panels").className = multi ? "multi" : "";
  $("status").innerHTML = multi
    ? `${selected.length} genes. TPM compares genes within a dataset, not across datasets.`
    : `<b>${esc(gname)}</b> · ${gid} · ${gtype.replace(/_/g, " ")}`;
  $("panels").innerHTML = Data.meta.datasets.map((d) => `
    <section class="panel">
      <h2>${d.label}</h2>
      <p class="sub">${d.samples.length} samples</p>
      <div class="plot" id="plot-${d.id}"></div>
      <div id="tab-${d.id}"></div>
    </section>`).join("");
  for (const d of Data.meta.datasets) {
    if (multi) { plotMulti(d, data); tableMulti(d, data); }
    else { plotOne(d, data[0]); tableOne(d, data[0]); }
  }
}

function yaxis(log, extra) {
  return { title: { text: log ? "log₂(TPM + 1)" : "TPM" }, gridcolor: "#eeeeee", zeroline: false,
    fixedrange: true, ticks: "outside", tickcolor: "#bbb", ...extra };
}
const LAYOUT = { font: PLOT_FONT, plot_bgcolor: "#fff", paper_bgcolor: "#fff",
  hoverlabel: { bgcolor: "#fff", font: { color: "#1a1a1a" } } };

function plotOne(d, g) {
  const log = isLog(), tpm = g.tpm[d.id];
  const y = tpm.map((v) => toY(v, log));
  const xi = Object.fromEntries(d.groups.map((grp, k) => [grp, k]));
  const xs = d.samples.map((s, k) => xi[s.group] + jitter(k) * 0.3);
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
    const ks = d.samples.map((_, k) => k).filter((k) => d.samples[k].group === grp);
    traces.push({
      x: ks.map((k) => xs[k]), y: ks.map((k) => y[k]), mode: "markers", showlegend: false,
      marker: { color: COLORS[grp], size: 10, line: { color: "#fff", width: 1.5 } },
      text: ks.map((k) => `${d.samples[k].sample}<br>TPM ${fmtTPM(tpm[k])}`), hovertemplate: "%{text}<extra></extra>",
    });
  }
  const ymax = Math.max(0, ...y.filter((v) => v != null)) || 1;
  const { shapes, annotations, top } = brackets(d, g.de[d.id], xi, ymax);
  Plotly.react(`plot-${d.id}`, traces, {
    ...LAYOUT, shapes, annotations, margin: { l: 52, r: 8, t: 8, b: 48 },
    xaxis: { tickvals: d.groups.map((_, k) => k), ticktext: d.groups.map((grp) => groupLabel(grp)),
      range: [-0.5, d.groups.length - 0.5], tickangle: 0, showgrid: false, zeroline: false, fixedrange: true },
    yaxis: yaxis(log, { range: [0, top] }),
  }, PLOT_CONFIG);
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
      hovertext: `${contrastLabel(it.c)}: p<sub>adj</sub> ${fmtP(it.padj)}` });
  }
  return { shapes, annotations, top: ymax * 1.08 + levels.length * step + 0.04 * ymax };
}

function tableOne(d, g) {
  const de = g.de[d.id], name = encodeURIComponent(Data.genes[selected[0]][1]);
  const rows = d.contrasts.map((c) => {
    const [base, lfc, padj] = de[c] || [];
    return `<tr><td><a href="de.html#c=${d.id}__${c}&gene=${name}" title="See this comparison">${contrastLabel(c)}</a></td>
      <td>${fmt(lfc)}</td><td class="${padj != null && padj < 0.05 ? "sig" : ""}">${fmtP(padj)}</td><td>${base == null ? "—" : Math.round(base)}</td></tr>`;
  });
  $(`tab-${d.id}`).innerHTML = `
    <table><tr><th>Expression percentile</th>${d.groups.map((grp) => `<th>${groupLabel(grp)}</th>`).join("")}</tr>
      <tr><td>vs coding genes</td>${g.pct[d.id].map((p) => `<td>${p}</td>`).join("")}</tr></table>
    <table><tr><th>DESeq2</th><th>log₂FC</th><th>p<sub>adj</sub></th><th>baseMean</th></tr>${rows.join("")}</table>`;
}

// several genes: genes on x, groups side by side within each gene
function plotMulti(d, data) {
  const log = isLog(), G = d.groups.length, slot = 0.72 / G;
  const traces = d.groups.map((grp, gi) => {
    const x = [], y = [], text = [];
    data.forEach((g, n) => {
      d.samples.forEach((s, k) => {
        if (s.group !== grp) return;
        const v = g.tpm[d.id][k];
        x.push(n + (gi - (G - 1) / 2) * slot + jitter(k + 31 * n) * slot * 0.5);
        y.push(toY(v, log));
        text.push(`${esc(Data.genes[selected[n]][1])} · ${s.sample}<br>TPM ${fmtTPM(v)}`);
      });
    });
    return { x, y, text, mode: "markers", name: groupLabel(grp, false), hovertemplate: "%{text}<extra></extra>",
      marker: { color: COLORS[grp], size: 8, line: { color: "#fff", width: 1 } } };
  });
  const shapes = selected.slice(1).map((_, n) => ({ type: "line", xref: "x", yref: "paper", x0: n + 0.5, x1: n + 0.5,
    y0: 0, y1: 1, line: { color: "#eeeeee", width: 1 } }));
  Plotly.react(`plot-${d.id}`, traces, {
    ...LAYOUT, shapes, margin: { l: 52, r: 8, t: 30, b: 40 },
    legend: { orientation: "h", x: 0, y: 1.02, yanchor: "bottom" },
    xaxis: { tickvals: selected.map((_, n) => n), ticktext: selected.map((i) => esc(Data.genes[i][1])),
      range: [-0.5, selected.length - 0.5], tickangle: 0, showgrid: false, zeroline: false, fixedrange: true },
    yaxis: yaxis(log, { rangemode: "tozero" }),
  }, PLOT_CONFIG);
}

function tableMulti(d, data) {
  const head = `<tr><th>Mean TPM (percentile)</th>${d.groups.map((grp) => `<th>${groupLabel(grp, false)}</th>`).join("")}</tr>`;
  const rows = data.map((g, n) => {
    const cells = d.groups.map((grp, gi) => {
      const v = d.samples.map((s, k) => (s.group === grp ? g.tpm[d.id][k] : null)).filter((x) => x != null);
      return `<td>${fmtTPM(v.reduce((a, b) => a + b, 0) / v.length)} <span class="muted">(${g.pct[d.id][gi]})</span></td>`;
    });
    const name = Data.genes[selected[n]][1];
    return `<tr><td><a href="${explorerLink([name])}">${esc(name)}</a></td>${cells.join("")}</tr>`;
  });
  $(`tab-${d.id}`).innerHTML = `<table>${head}${rows.join("")}</table>`;
}

init().catch((e) => { $("status").textContent = "Could not load data: " + e.message; });
