// Genes page. One gene: a panel per dataset with every sample, DESeq2 brackets and the
// stats. Several genes: TPM by gene, so genes can be compared with each other.
const MAX_GENES = 12;
const GENE_SETS = {
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
// the y axis always reaches 10 TPM, so a gene that is barely expressed looks that way
const FLOOR_TPM = 10;

let selected = [];   // gene indices, in the order added
let scale = "log";
let renderId = 0;
let current = [];    // data for the selected genes, for the CSV download
let lastMode = null; // "one" | "multi": plots are rebuilt when the layout switches
let undoSel = null;  // the selection a gene set replaced

const isLog = () => scale === "log";
const toY = (v) => (v == null ? null : isLog() ? Math.log2(v + 1) : v);
const yTop = (ymax) => Math.max(ymax * 1.05, toY(FLOOR_TPM));
// deterministic jitter in [-0.5, 0.5] so points do not move between genes
const jitter = (k) => ((Math.abs(Math.sin(k * 12.9898) * 43758.5453) % 1) - 0.5);
const mean = (v) => v.reduce((a, b) => a + b, 0) / v.length;
const status = (msg, html = false) => { const el = $("status"); el[html ? "innerHTML" : "textContent"] = msg; el.hidden = !msg; };

async function init() {
  if (typeof Plotly === "undefined") { status("The plotting library did not load. Reload the page."); return; }
  try { if (localStorage.getItem("scale") === "lin") scale = "lin"; } catch { /* storage unavailable */ }
  await Data.load();

  $("genesets").insertAdjacentHTML("beforeend", Object.keys(GENE_SETS).map((p) =>
    `<button type="button" data-set="${esc(p)}">${esc(p)}</button>`).join(""));
  $("genesets").addEventListener("click", (e) => {
    const p = e.target.closest("[data-set]")?.dataset.set;
    if (p) applySet(p);
  });
  $("status").addEventListener("click", (e) => {
    if (e.target.closest("[data-undo]") && undoSel) { const u = undoSel; undoSel = null; setGenes(u); }
  });

  attachSearch($("gene"), $("suggest"), (i) => { $("gene").value = ""; addGene(i); },
    (q) => status(`No gene called “${q}”. Try its official symbol, for example HMGA2.`));
  $("gene").addEventListener("keydown", (e) => {
    if (e.key === "Backspace" && !$("gene").value && selected.length) setGenes(selected.slice(0, -1));
  });
  $("chipbox").addEventListener("click", (e) => {
    const i = e.target.closest("[data-remove]")?.dataset.remove;
    if (i != null) setGenes(selected.filter((g) => g !== +i));
    else $("gene").focus();
  });

  document.querySelectorAll("[data-scale]").forEach((b) => b.addEventListener("click", () => setScale(b.dataset.scale)));
  setScale(scale, false);
  $("copy").addEventListener("click", (e) => copyText(location.href, e.currentTarget));
  $("csv").addEventListener("click", downloadCSV);

  $("panels").innerHTML = Data.meta.datasets.map((d) => `
    <section class="panel" aria-labelledby="h-${d.id}">
      <header><h2 id="h-${d.id}">${SHORT_DS[d.id] || d.label}</h2><p class="note">${datasetNote(d.id)}</p></header>
      <div class="plot" id="plot-${d.id}"></div>
      <div class="tables" id="tab-${d.id}"></div>
    </section>`).join("");
  // Plotly sizes a plot when it is created; redraw when the panel's width changes
  new ResizeObserver(() => document.querySelectorAll(".plot").forEach((el) => el.data && Plotly.Plots.resize(el))).observe($("panels"));

  window.addEventListener("hashchange", fromHash);
  fromHash();
}

function setScale(s, rerender = true) {
  scale = s;
  document.querySelectorAll("[data-scale]").forEach((b) => b.setAttribute("aria-checked", b.dataset.scale === s));
  try { localStorage.setItem("scale", s); } catch { /* storage unavailable */ }
  if (rerender && selected.length) render();
}

function fromHash() {
  const h = new URLSearchParams(location.hash.slice(1));
  const names = (h.get("genes") || h.get("gene") || "MED12").split(",");
  const idx = names.map((n) => Data.lookup(n)[0]).filter((i) => i != null);
  setGenes(idx.length ? idx : [Data.lookup("MED12")[0]], false);
}

function addGene(i) {
  if (selected.includes(i)) return;
  if (selected.length >= MAX_GENES) { status(`Up to ${MAX_GENES} genes at a time.`); return; }
  setGenes([...selected, i]);
}

// A gene set replaces the selection; the status line offers to undo that.
async function applySet(name) {
  const idx = GENE_SETS[name].map((n) => Data.lookup(n)[0]).filter((i) => i != null);
  const before = selected.slice();
  const same = before.length === idx.length && before.every((v, k) => v === idx[k]);
  await setGenes(idx);
  if (before.length && !same) {
    undoSel = before;
    const what = before.length === 1 ? `your gene (${esc(Data.genes[before[0]][1])})` : `your ${before.length} genes`;
    status(`Showing the ${esc(name)} set instead of ${what}. <button type="button" class="linklike" data-undo>Undo</button>`, true);
  }
}

function setGenes(idx, updateHash = true) {
  selected = [...new Set(idx)].slice(0, MAX_GENES);
  $("chips").innerHTML = selected.map((i) =>
    `<span class="chip">${esc(Data.genes[i][1])}<button type="button" data-remove="${i}" aria-label="Remove ${esc(Data.genes[i][1])}">×</button></span>`).join("");
  $("gene").placeholder = selected.length ? "Add another gene" : "Type a gene symbol, e.g. MED12";
  if (updateHash) history.replaceState(null, "", "#genes=" + selected.map((i) => encodeURIComponent(Data.genes[i][1])).join(","));
  return render();
}

async function render() {
  const id = ++renderId;
  if (!selected.length) {
    $("gene-head").hidden = true;
    status("Add a gene to start.");
    for (const d of Data.meta.datasets) { Plotly.purge(`plot-${d.id}`); $(`tab-${d.id}`).innerHTML = ""; }
    lastMode = null;
    return;
  }
  status("Loading…");
  let data;
  try { data = await Promise.all(selected.map((i) => Data.gene(i))); }
  catch (e) { if (id === renderId) status("The gene data did not load. Check your connection and try again."); return; }
  if (id !== renderId) return;
  current = data;
  const multi = selected.length > 1;
  $("gene-head").hidden = false;
  $("panels").classList.toggle("multi", multi);
  const mode = multi ? "multi" : "one";
  if (mode !== lastMode) document.querySelectorAll(".plot").forEach((el) => el.data && Plotly.purge(el));
  lastMode = mode;
  if (multi) {
    $("gene-name").textContent = `${selected.length} genes`;
    $("gene-meta").textContent = "Compare genes within a panel; the three datasets are not on the same scale.";
  } else {
    const [gid, gname, gtype] = Data.genes[selected[0]];
    $("gene-name").textContent = gname;
    $("gene-meta").textContent = `${shortId(gid)} · ${gtype.replace(/_/g, " ")}`;
  }
  try {
    for (const d of Data.meta.datasets) {
      if (multi) { plotMulti(d, data); tableMulti(d, data); }
      else { plotOne(d, data[0]); tableOne(d, data[0]); }
    }
    status("");
  } catch (e) {
    console.error(e);
    status("Could not draw the plots. Reload the page; if it happens again, please report it.");
  }
}

function plotOne(d, g) {
  const tpm = g.tpm[d.id];
  const y = tpm.map(toY);
  const xi = Object.fromEntries(d.groups.map((grp, k) => [grp, k]));
  const xs = d.samples.map((s, k) => xi[s.group] + jitter(k) * 0.32);
  const traces = [];
  // tissue: one grey line per patient, MM -> LM
  const pts = {};
  d.samples.forEach((s, k) => { if (s.patient) (pts[s.patient] ||= []).push(k); });
  for (const ks of Object.values(pts)) {
    if (ks.length !== 2) continue;
    traces.push({ x: ks.map((k) => xs[k]), y: ks.map((k) => y[k]), mode: "lines",
      line: { color: "#dde0e7", width: 1 }, hoverinfo: "skip", showlegend: false });
  }
  for (const grp of d.groups) {
    const ks = d.samples.map((_, k) => k).filter((k) => d.samples[k].group === grp);
    traces.push({
      x: ks.map((k) => xs[k]), y: ks.map((k) => y[k]), mode: "markers", showlegend: false,
      marker: { color: COLORS[grp], size: 10, line: { color: "#fff", width: 1.5 } },
      text: ks.map((k) => `${sampleLabel(d.samples[k].sample)}<br>TPM ${fmtTPM(tpm[k])}`),
      hovertemplate: "%{text}<extra></extra>",
    });
  }
  const ymax = yTop(Math.max(0, ...y.filter((v) => v != null)));
  const { shapes, annotations, top } = brackets(d, g.de[d.id], xi, ymax);
  Plotly.react(`plot-${d.id}`, traces, {
    ...BASE_LAYOUT, shapes, annotations, margin: { l: 56, r: 10, t: 10, b: 50 },
    xaxis: { ...AXIS, tickvals: d.groups.map((_, k) => k), ticktext: d.groups.map(groupTick),
      range: [-0.5, d.groups.length - 0.5], tickangle: 0, showgrid: false, ticks: "" },
    yaxis: tpmAxis(isLog(), top),
  }, PLOT_CONFIG);
}

// One bracket per tested comparison, stacked so that none overlap. The label is the
// DESeq2 padj as stars; a gene DESeq2 did not test gets no bracket.
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
  const step = 0.12 * ymax, tick = 0.03 * ymax, inset = 0.1;
  const shapes = [], annotations = [];
  for (const it of items) {
    const yb = ymax * 1.02 + it.level * step;
    const [x1, x2] = [it.lo + inset, it.hi - inset];
    shapes.push({ type: "path", xref: "x", yref: "y", path: `M${x1},${yb - tick} L${x1},${yb} L${x2},${yb} L${x2},${yb - tick}`,
      line: { color: MUTED, width: 1 } });
    const sig = it.padj < 0.05;
    annotations.push({ x: (it.lo + it.hi) / 2, y: yb, yanchor: "bottom", showarrow: false, text: stars(it.padj),
      font: { size: sig ? 15 : 11, color: sig ? INK : MUTED }, yshift: sig ? -2 : 1,
      hovertext: `${contrastLabel(it.c)}<br>p<sub>adj</sub> ${fmtP(it.padj)}` });
  }
  return { shapes, annotations, top: ymax * 1.02 + levels.length * step + 0.05 * ymax };
}

function groupCell(grp) {
  return `<span class="swatch" style="background:${COLORS[grp]}"></span>${esc(groupLabel(grp))}`;
}

function tableOne(d, g) {
  const tpm = g.tpm[d.id], de = g.de[d.id], name = encodeURIComponent(Data.genes[selected[0]][1]);
  const expr = d.groups.map((grp, gi) => {
    const v = d.samples.map((s, k) => (s.group === grp ? tpm[k] : null)).filter((x) => x != null);
    return `<tr><td>${groupCell(grp)}</td><td>${fmtTPM(mean(v))}</td><td>${ordinal(g.pct[d.id][gi])}</td></tr>`;
  });
  const diffs = d.contrasts.map((c) => {
    const [, lfc, padj] = de[c] || [];
    const label = `<a href="de.html#c=${d.id}__${c}&gene=${name}" title="See the whole comparison">${esc(contrastLabel(c))}</a>`;
    if (padj == null) return `<tr><td>${label}</td><td colspan="2" class="dim">not tested</td></tr>`;
    // no direction colour here: on this page orange and blue already mean sample groups
    const cls = padj < 0.05 ? "sig" : "";
    return `<tr><td>${label}</td><td class="${cls}">${fmt(lfc)}</td><td class="${cls}">${fmtP(padj)}</td></tr>`;
  });
  $(`tab-${d.id}`).innerHTML = `
    <table><caption>Expression by group</caption>
      <thead><tr><th scope="col">Group</th><th scope="col">Mean TPM</th><th scope="col">Percentile</th></tr></thead>
      <tbody>${expr.join("")}</tbody></table>
    <table><caption>Differences (DESeq2)</caption>
      <thead><tr><th scope="col">Comparison</th><th scope="col">log₂FC</th><th scope="col">p<sub>adj</sub></th></tr></thead>
      <tbody>${diffs.join("")}</tbody></table>`;
}

// several genes: genes on x, groups side by side within each gene
function plotMulti(d, data) {
  const G = d.groups.length, slot = 0.74 / G;
  const traces = d.groups.map((grp, gi) => {
    const x = [], y = [], text = [];
    data.forEach((g, n) => {
      d.samples.forEach((s, k) => {
        if (s.group !== grp) return;
        const v = g.tpm[d.id][k];
        x.push(n + (gi - (G - 1) / 2) * slot + jitter(k + 31 * n) * slot * 0.5);
        y.push(toY(v));
        text.push(`${esc(Data.genes[selected[n]][1])} · ${sampleLabel(s.sample)}<br>TPM ${fmtTPM(v)}`);
      });
    });
    return { x, y, text, mode: "markers", name: groupLabel(grp), hovertemplate: "%{text}<extra></extra>",
      marker: { color: COLORS[grp], size: 8, line: { color: "#fff", width: 1 } } };
  });
  const shapes = selected.slice(1).map((_, n) => ({ type: "line", xref: "x", yref: "paper", x0: n + 0.5, x1: n + 0.5,
    y0: 0, y1: 1, line: { color: RULE, width: 1 } }));
  // one star above a group where it differs from the reference group (padj < 0.05); a single
  // glyph so that neighbouring groups never collide. Comparisons that are not against the
  // reference stay in the table below
  const ref = d.groups[0], annotations = [];
  data.forEach((g, n) => {
    for (const c of d.contrasts) {
      const pair = c in CONTRAST_GROUPS ? CONTRAST_GROUPS[c] : c.split("_vs_");
      const padj = g.de[d.id][c]?.[2];
      if (!pair || pair[1] !== ref || !(pair[0] in Object.fromEntries(d.groups.map((x) => [x, 1]))) || padj == null || padj >= 0.05) continue;
      const gi = d.groups.indexOf(pair[0]);
      const ys = d.samples.map((s, k) => (s.group === pair[0] ? toY(g.tpm[d.id][k]) : null)).filter((v) => v != null);
      annotations.push({ x: n + (gi - (G - 1) / 2) * slot, y: Math.max(...ys), yanchor: "bottom", yshift: 4, showarrow: false,
        text: "*", font: { size: 16, color: INK }, hovertext: `${contrastLabel(c)}<br>p<sub>adj</sub> ${fmtP(padj)}` });
    }
  });
  const top = yTop(Math.max(0, ...traces.flatMap((t) => t.y.filter((v) => v != null)))) * (annotations.length ? 1.08 : 1);
  Plotly.react(`plot-${d.id}`, traces, {
    ...BASE_LAYOUT, shapes, annotations, margin: { l: 56, r: 10, t: 34, b: 44 },
    legend: { orientation: "h", x: 0, y: 1.02, yanchor: "bottom", font: { size: 13 }, itemclick: false, itemdoubleclick: false },
    xaxis: { ...AXIS, tickvals: selected.map((_, n) => n), ticktext: selected.map((i) => esc(Data.genes[i][1])),
      range: [-0.5, selected.length - 0.5], tickangle: "auto", showgrid: false, ticks: "" },
    yaxis: tpmAxis(isLog(), top),
  }, PLOT_CONFIG);
}

// which comparisons call the gene significant, as short labels with the direction
function changesCell(d, de) {
  if (!de || !Object.keys(de).length) return `<span class="dim">not tested</span>`;
  const hits = d.contrasts.filter((c) => de[c] && de[c][2] < 0.05)
    .map((c) => `${esc(contrastShort(c))} ${de[c][1] > 0 ? "↑" : "↓"}`);
  return hits.length ? hits.join(", ") : `<span class="dim">none</span>`;
}

function tableMulti(d, data) {
  const vsWT = d.contrasts.every((c) => c.endsWith("_vs_WT"));
  const head = `<tr><th scope="col">Gene</th>${d.groups.map((grp) => `<th scope="col">${groupCell(grp)}</th>`).join("")}
    <th scope="col" class="left">Significant changes${vsWT ? " vs WT" : ""}</th></tr>`;
  const rows = data.map((g, n) => {
    const cells = d.groups.map((grp, gi) => {
      const v = d.samples.map((s, k) => (s.group === grp ? g.tpm[d.id][k] : null)).filter((x) => x != null);
      return `<td>${fmtTPM(mean(v))} <span class="dim small">${ordinal(g.pct[d.id][gi])}</span></td>`;
    });
    const name = Data.genes[selected[n]][1];
    return `<tr><td><a href="${genesLink([name])}" title="Show only this gene">${esc(name)}</a></td>${cells.join("")}<td class="left">${changesCell(d, g.de[d.id])}</td></tr>`;
  });
  $(`tab-${d.id}`).innerHTML = `<div class="scroll"><table><caption>Mean TPM, with its percentile among protein-coding genes; changes at p<sub>adj</sub> &lt; 0.05</caption>
    <thead>${head}</thead><tbody>${rows.join("")}</tbody></table></div>`;
}

// long format: one row per sample and gene, for the selected genes
function downloadCSV() {
  const lines = ["dataset,sample,group,patient,gene,gene_id,TPM"];
  Data.meta.datasets.forEach((d) => d.samples.forEach((s, k) => selected.forEach((gi, n) => {
    lines.push([d.id, s.sample, s.group, s.patient ?? "", Data.genes[gi][1], Data.genes[gi][0], current[n].tpm[d.id][k] ?? ""].join(","));
  })));
  const stem = selected.map((i) => Data.genes[i][1]).join("_").slice(0, 60);
  saveText(lines.join("\n") + "\n", `TPM_${stem}.csv`);
}

init().catch((e) => status("The data did not load: " + e.message));
