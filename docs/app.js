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
let noPlots = false; // the plotting library did not load; tables only

const isLog = () => scale === "log";
const toY = (v) => (v == null ? null : isLog() ? Math.log2(v + 1) : v);
const yTop = (ymax) => Math.max(ymax * 1.05, toY(FLOOR_TPM));
// deterministic jitter in [-0.5, 0.5] so points do not move between genes
const jitter = (k) => ((Math.abs(Math.sin(k * 12.9898) * 43758.5453) % 1) - 0.5);
const mean = (v) => v.reduce((a, b) => a + b, 0) / v.length;
const status = (msg, html = false) => { const el = $("status"); el[html ? "innerHTML" : "textContent"] = msg; el.hidden = !msg; };
// the groups a comparison is about; the pooled LM vs MM involves every group
const pairFor = (d, c) => (c in CONTRAST_GROUPS ? CONTRAST_GROUPS[c] : c.split("_vs_")) || d.groups;

// Linked highlighting: hovering a table row, or one "KO ↓" entry, lights up the points it
// refers to and dims the rest. hl[dataset] records what each trace (or each point) stands for.
const hl = {};
function highlight(dsId, groups, gene) {
  const el = $(`plot-${dsId}`), h = hl[dsId];
  if (!el?.data || !h) return;
  const key = `${groups ? groups.join(",") : ""}|${gene ?? ""}`;
  if (h.key === key) return;
  h.key = key;
  const clear = !groups && gene == null;
  const on = (g, n) => (!groups || groups.includes(g)) && (gene == null || n === gene);
  if (h.mode === "one") {
    Plotly.restyle(el, { opacity: clear ? 1 : h.traces.map((g) => (g === null ? 0.3 : on(g, null) ? 1 : 0.15)) });
    if (h.brackets.length) {
      const lay = {};
      h.brackets.forEach((pair, i) => {
        const o = clear || pair.every((g) => on(g, null)) ? 1 : 0.2;
        lay[`shapes[${i}].opacity`] = o; lay[`annotations[${i}].opacity`] = o;
      });
      Plotly.relayout(el, lay);
    }
  } else if (clear) {
    Plotly.restyle(el, { "marker.opacity": 1 });
  } else {
    Plotly.restyle(el, { "marker.opacity": h.points.map((pts) => pts.map((p) => (on(p.g, p.n) ? 1 : 0.15))) });
  }
}

async function init() {
  try { if (localStorage.getItem("scale") === "lin") scale = "lin"; } catch { /* storage unavailable */ }
  const [plotly] = await Promise.all([ensurePlotly("vendor/plotly-basic-2.35.2.min.js"), Data.load()]);
  noPlots = !plotly;

  $("genesets").insertAdjacentHTML("beforeend", Object.keys(GENE_SETS).map((p) =>
    `<button type="button" data-set="${esc(p)}">${esc(p)}</button>`).join(""));
  $("genesets").addEventListener("click", (e) => {
    const p = e.target.closest("[data-set]")?.dataset.set;
    if (p) applySet(p);
  });
  $("status").addEventListener("click", (e) => {
    if (e.target.closest("[data-undo]") && undoSel) { const u = undoSel; undoSel = null; setGenes(u); }
    const t = e.target.closest("[data-try]")?.dataset.try;
    if (t) setGenes([Data.lookup(t)[0]]);
  });

  attachSearch($("gene"), $("suggest"), (i) => { $("gene").value = ""; addGene(i); },
    (q) => status(`No gene called “${q}”. Try its official symbol, for example HMGA2.`));
  $("gene").addEventListener("keydown", (e) => {
    if (e.key === "Backspace" && !$("gene").value && selected.length) setGenes(selected.slice(0, -1));
  });
  $("chipbox").addEventListener("click", (e) => {
    const i = e.target.closest("[data-remove]")?.dataset.remove;
    if (i != null) setGenes(selected.filter((g) => g !== +i));
    else if (!e.target.closest("#clear")) $("gene").focus();
  });
  $("clear").addEventListener("click", clearGenes);
  wireDrag();

  document.querySelectorAll("[data-scale]").forEach((b) => b.addEventListener("click", () => setScale(b.dataset.scale)));
  setScale(scale, false);
  $("copy").addEventListener("click", (e) => copyText(location.href, e.currentTarget));
  $("csv").addEventListener("click", downloadCSV);

  $("panels").innerHTML = Data.meta.datasets.map((d) => `
    <section class="panel" aria-labeledby="h-${d.id}">
      <header><h2 id="h-${d.id}">${SHORT_DS[d.id] || d.label}</h2><p class="note">${datasetNote(d.id)}</p></header>
      <div class="plot" id="plot-${d.id}"></div>
      <div class="tables" id="tab-${d.id}"></div>
    </section>`).join("");
  // Plotly sizes a plot when it is created; redraw when the panel's width changes
  if (!noPlots) new ResizeObserver(() => document.querySelectorAll(".plot").forEach((el) => el.data && Plotly.Plots.resize(el))).observe($("panels"));
  for (const d of Data.meta.datasets) {
    const tab = $(`tab-${d.id}`);
    const uncolor = () => tab.querySelectorAll("col.hl-col").forEach((c) => c.classList.remove("hl-col"));
    tab.addEventListener("mouseover", (e) => {
      uncolor();
      const th = e.target.closest("th[data-hl-groups]");   // a group header also shades its column
      if (th) th.closest("table").querySelectorAll("col")[th.cellIndex]?.classList.add("hl-col");
      const t = e.target.closest("[data-hl-groups],[data-hl-gene]");
      if (!t) { highlight(d.id, null, null); return; }
      highlight(d.id, t.dataset.hlGroups ? t.dataset.hlGroups.split(",") : null, t.dataset.hlGene != null ? +t.dataset.hlGene : null);
    });
    tab.addEventListener("mouseleave", () => { uncolor(); highlight(d.id, null, null); });
  }

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

// Reorder the genes by dragging a chip. Pointer events rather than native drag-and-drop, so the
// cursor stays a hand, the chip follows the pointer as a crisp copy, the other chips part to make
// room, and it works on touch screens. Alt+arrows move a focused chip; Delete removes it.
function renderChips() {
  $("chips").innerHTML = selected.map((i) =>
    `<span class="chip" data-i="${i}" tabindex="0" title="Drag to reorder">${esc(Data.genes[i][1])}<button type="button" data-remove="${i}" tabindex="-1" aria-label="Remove ${esc(Data.genes[i][1])}">×</button></span>`).join("");
}

function wireDrag() {
  const box = $("chips");
  let drag = null;
  const end = () => {
    if (!drag) return;
    drag.ghost?.remove();
    document.body.classList.remove("dragging-chip");
    drag.chip.classList.remove("placeholder");
    drag = null;
  };
  box.addEventListener("pointerdown", (e) => {
    const chip = e.target.closest(".chip");
    if (!chip || e.target.closest("button") || (e.pointerType === "mouse" && e.button !== 0)) return;
    drag = { chip, id: e.pointerId, x0: e.clientX, y0: e.clientY, moved: false };
    chip.setPointerCapture(e.pointerId);
  });
  box.addEventListener("pointermove", (e) => {
    if (!drag || e.pointerId !== drag.id) return;
    if (!drag.moved) {
      if (Math.hypot(e.clientX - drag.x0, e.clientY - drag.y0) < 5) return;   // a click, not a drag
      drag.moved = true;
      const r = drag.chip.getBoundingClientRect();
      drag.dx = e.clientX - r.left; drag.dy = e.clientY - r.top;
      drag.ghost = drag.chip.cloneNode(true);
      drag.ghost.className = "chip ghost"; drag.ghost.style.width = `${r.width}px`;
      document.body.appendChild(drag.ghost);
      document.body.classList.add("dragging-chip");
      drag.chip.classList.add("placeholder");
    }
    drag.ghost.style.transform = `translate(${e.clientX - drag.dx}px, ${e.clientY - drag.dy}px)`;
    // the placeholder moves among the other chips: before the first chip in the pointer's row whose
    // midpoint is right of the pointer, else after the last chip in that row
    const others = [...box.querySelectorAll(".chip:not(.placeholder)")];
    const inRow = (c) => { const r = c.getBoundingClientRect(); return e.clientY >= r.top - 4 && e.clientY <= r.bottom + 4; };
    const before = others.find((c) => inRow(c) && e.clientX < c.getBoundingClientRect().left + c.getBoundingClientRect().width / 2);
    if (before) { if (before !== drag.chip.nextElementSibling) box.insertBefore(drag.chip, before); }
    else { const row = others.filter(inRow); if (row.length && row.at(-1) !== drag.chip.previousElementSibling) row.at(-1).after(drag.chip); }
  });
  const finish = (e) => {
    if (!drag || e.pointerId !== drag.id) return;
    const moved = drag.moved;
    end();
    if (!moved) return;
    const order = [...box.querySelectorAll(".chip")].map((c) => +c.dataset.i);
    if (order.some((v, k) => v !== selected[k])) setGenes(order); else renderChips();
  };
  box.addEventListener("pointerup", finish);
  box.addEventListener("pointercancel", (e) => { if (drag && e.pointerId === drag.id) { end(); renderChips(); } });
  // keyboard: Alt+arrows move the focused chip, Delete or Backspace removes it
  box.addEventListener("keydown", (e) => {
    const chip = e.target.closest(".chip");
    if (!chip) return;
    const i = +chip.dataset.i, k = selected.indexOf(i);
    if (e.altKey && (e.key === "ArrowLeft" || e.key === "ArrowRight")) {
      const j = e.key === "ArrowLeft" ? k - 1 : k + 1;
      if (j < 0 || j >= selected.length) return;
      e.preventDefault();
      const order = selected.slice(); [order[k], order[j]] = [order[j], order[k]];
      setGenes(order).then(() => box.querySelector(`.chip[data-i="${i}"]`)?.focus());
    } else if (e.key === "Delete" || e.key === "Backspace") {
      e.preventDefault();
      setGenes(selected.filter((g) => g !== i)).then(() => ($("chips").querySelector(".chip") || $("gene")).focus());
    }
  });
}

// Clear every gene at once; the status line offers to undo that.
async function clearGenes() {
  const before = selected.slice();
  await setGenes([]);
  undoSel = before;
  status(`Cleared ${before.length} genes. <button type="button" class="linklike" data-undo>Undo</button>, or try <button type="button" class="linklike" data-try="MED12">MED12</button>.`, true);
  $("gene").focus();
}

function setGenes(idx, updateHash = true) {
  selected = [...new Set(idx)].slice(0, MAX_GENES);
  $("clear").hidden = selected.length < 2;
  renderChips();
  $("gene").placeholder = selected.length ? "Add another gene" : "Type a gene symbol, e.g. MED12";
  if (updateHash) history.replaceState(null, "", "#genes=" + selected.map((i) => encodeURIComponent(Data.genes[i][1])).join(","));
  return render();
}

async function render() {
  const id = ++renderId;
  if (!selected.length) {
    $("gene-head").hidden = true;
    status(`Add a gene to start, or try <button type="button" class="linklike" data-try="MED12">MED12</button>.`, true);
    for (const d of Data.meta.datasets) { if (!noPlots) Plotly.purge(`plot-${d.id}`); $(`tab-${d.id}`).innerHTML = ""; }
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
  if (mode !== lastMode && !noPlots) document.querySelectorAll(".plot").forEach((el) => el.data && Plotly.purge(el));
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
  if (noPlots) { $(`plot-${d.id}`).innerHTML = NO_PLOT; return; }
  const tpm = g.tpm[d.id];
  const y = tpm.map(toY);
  const xi = Object.fromEntries(d.groups.map((grp, k) => [grp, k]));
  const xs = d.samples.map((s, k) => xi[s.group] + jitter(k) * 0.32);
  const traces = [], what = [];   // what[i]: the group a trace shows, or null for a pairing line
  // tissue: one gray line per patient, MM -> LM
  const pts = {};
  d.samples.forEach((s, k) => { if (s.patient) (pts[s.patient] ||= []).push(k); });
  for (const ks of Object.values(pts)) {
    if (ks.length !== 2) continue;
    what.push(null);
    traces.push({ x: ks.map((k) => xs[k]), y: ks.map((k) => y[k]), mode: "lines",
      line: { color: "#cfd3dc", width: 1 }, hoverinfo: "skip", showlegend: false });
  }
  for (const grp of d.groups) {
    const ks = d.samples.map((_, k) => k).filter((k) => d.samples[k].group === grp);
    what.push(grp);
    traces.push({
      x: ks.map((k) => xs[k]), y: ks.map((k) => y[k]), mode: "markers", showlegend: false,
      marker: { color: COLORS[grp], size: 10, line: { color: "#fff", width: 1.5 } },
      text: ks.map((k) => `${sampleLabel(d.samples[k].sample)}<br>TPM ${fmtTPM(tpm[k])}`),
      hovertemplate: "%{text}<extra></extra>",
    });
  }
  const ymax = yTop(Math.max(0, ...y.filter((v) => v != null)));
  const { shapes, annotations, items, top } = brackets(d, g.de[d.id], xi, ymax);
  hl[d.id] = { mode: "one", traces: what, brackets: items.map((it) => it.pair) };
  Plotly.react(`plot-${d.id}`, traces, {
    ...BASE_LAYOUT, dragmode: false, shapes, annotations, margin: { l: 56, r: 10, t: 10, b: 50 },
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
    items.push({ c, pair, lo, hi, padj });
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
  return { shapes, annotations, items, top: ymax * 1.02 + levels.length * step + 0.05 * ymax };
}

function groupCell(grp) {
  return `<span class="swatch" style="background:${COLORS[grp]}"></span>${esc(groupLabel(grp))}`;
}

function tableOne(d, g) {
  const tpm = g.tpm[d.id], de = g.de[d.id], name = encodeURIComponent(Data.genes[selected[0]][1]);
  const expr = d.groups.map((grp, gi) => {
    const v = d.samples.map((s, k) => (s.group === grp ? tpm[k] : null)).filter((x) => x != null);
    return `<tr data-hl-groups="${grp}"><td>${groupCell(grp)}</td><td>${fmtTPM(mean(v))}</td><td>${ordinal(g.pct[d.id][gi])}</td></tr>`;
  });
  const diffs = d.contrasts.map((c) => {
    const [, lfc, padj] = de[c] || [];
    const label = `<a href="de.html#c=${d.id}__${c}&gene=${name}" title="See the whole comparison">${esc(contrastLabel(c))}</a>`;
    const row = `<tr data-hl-groups="${pairFor(d, c).join(",")}">`;
    if (padj == null) return `${row}<td>${label}</td><td colspan="2" class="dim">not tested</td></tr>`;
    // no direction color here: on this page orange and blue already mean sample groups
    const cls = padj < 0.05 ? "sig" : "";
    return `${row}<td>${label}</td><td class="${cls}">${fmt(lfc)}</td><td class="${cls}">${fmtP(padj)}</td></tr>`;
  });
  $(`tab-${d.id}`).innerHTML = `
    <table><caption>Expression by group</caption>
      <thead><tr><th scope="col">Group</th><th scope="col">Mean TPM</th><th scope="col">Percentile</th></tr></thead>
      <tbody>${expr.join("")}</tbody></table>
    <table><caption>Differences</caption>
      <thead><tr><th scope="col">Comparison</th><th scope="col">log₂FC</th><th scope="col">p<sub>adj</sub></th></tr></thead>
      <tbody>${diffs.join("")}</tbody></table>`;
}

// several genes: genes on x, groups side by side within each gene
function plotMulti(d, data) {
  if (noPlots) { $(`plot-${d.id}`).innerHTML = NO_PLOT; return; }
  const G = d.groups.length, slot = 0.74 / G;
  const points = [];   // per trace, what each point stands for: { g: group, n: gene index }
  const traces = d.groups.map((grp, gi) => {
    const x = [], y = [], text = [], pts = [];
    data.forEach((g, n) => {
      d.samples.forEach((s, k) => {
        if (s.group !== grp) return;
        const v = g.tpm[d.id][k];
        x.push(n + (gi - (G - 1) / 2) * slot + jitter(k + 31 * n) * slot * 0.5);
        y.push(toY(v));
        text.push(`${esc(Data.genes[selected[n]][1])} · ${sampleLabel(s.sample)}<br>TPM ${fmtTPM(v)}`);
        pts.push({ g: grp, n });
      });
    });
    points.push(pts);
    return { x, y, text, mode: "markers", name: groupLabel(grp), hovertemplate: "%{text}<extra></extra>",
      marker: { color: COLORS[grp], size: 8, line: { color: "#fff", width: 1 } } };
  });
  const shapes = selected.slice(1).map((_, n) => ({ type: "line", xref: "x", yref: "paper", x0: n + 0.5, x1: n + 0.5,
    y0: 0, y1: 1, line: { color: RULE, width: 1 } }));
  // one star above a group where it differs from the reference group (padj < 0.05), drawn as a
  // stroked marker so that it reads as a mark, not a speck. Comparisons that are not against
  // the reference stay in the table below
  const base = yTop(Math.max(0, ...traces.flatMap((t) => t.y.filter((v) => v != null))));
  const ref = d.groups[0], star = { x: [], y: [], text: [], pts: [] };
  data.forEach((g, n) => {
    for (const c of d.contrasts) {
      const pair = pairFor(d, c), padj = g.de[d.id][c]?.[2];
      if (pair.length !== 2 || pair[1] !== ref || !d.groups.includes(pair[0]) || padj == null || padj >= 0.05) continue;
      const gi = d.groups.indexOf(pair[0]);
      const ys = d.samples.map((s, k) => (s.group === pair[0] ? toY(g.tpm[d.id][k]) : null)).filter((v) => v != null);
      star.x.push(n + (gi - (G - 1) / 2) * slot);
      star.y.push(Math.max(...ys) + 0.055 * base);
      star.text.push(`${esc(Data.genes[selected[n]][1])}, ${contrastLabel(c)}<br>p<sub>adj</sub> ${fmtP(padj)}`);
      star.pts.push({ g: pair[0], n });
    }
  });
  if (star.x.length) {
    points.push(star.pts);
    traces.push({ x: star.x, y: star.y, text: star.text, mode: "markers", showlegend: false, hovertemplate: "%{text}<extra></extra>",
      marker: { symbol: "asterisk", size: 11, color: INK, line: { color: INK, width: 1.6 } } });
  }
  hl[d.id] = { mode: "multi", points };
  const top = base * (star.x.length ? 1.08 : 1);
  const byLabel = Object.fromEntries(d.groups.map((grp) => [groupLabel(grp), grp]));
  Plotly.react(`plot-${d.id}`, traces, {
    ...BASE_LAYOUT, dragmode: false, shapes, margin: { l: 56, r: 10, t: 34, b: 44 },
    legend: { orientation: "h", x: 0, y: 1.02, yanchor: "bottom", font: { size: 13 }, itemclick: false, itemdoubleclick: false },
    xaxis: { ...AXIS, tickvals: selected.map((_, n) => n), ticktext: selected.map((i) => esc(Data.genes[i][1])),
      range: [-0.5, selected.length - 0.5], tickangle: "auto", showgrid: false, ticks: "" },
    yaxis: tpmAxis(isLog(), top),
  }, PLOT_CONFIG);
  legendHover($(`plot-${d.id}`), (name) => highlight(d.id, name ? [byLabel[name]] : null, null));
}

// which comparisons call the gene significant, as short labels with the direction
function changesCell(d, de, n) {
  if (!de || !Object.keys(de).length) return `<span class="dim">not tested</span>`;
  const hits = d.contrasts.filter((c) => de[c] && de[c][2] < 0.05)
    .map((c) => `<span class="hl" data-hl-gene="${n}" data-hl-groups="${pairFor(d, c).join(",")}" title="padj ${fmtP(de[c][2])}">${esc(contrastShort(c))} ${de[c][1] > 0 ? "↑" : "↓"}</span>`);
  return hits.length ? hits.join(", ") : `<span class="dim">none</span>`;
}

function tableMulti(d, data) {
  const vsWT = d.contrasts.every((c) => c.endsWith("_vs_WT"));
  const head = `<tr><th scope="col">Gene</th>${d.groups.map((grp) => `<th scope="col" data-hl-groups="${grp}">${groupCell(grp)}</th>`).join("")}
    <th scope="col" class="left">Significant changes${vsWT ? " vs WT" : ""}</th></tr>`;
  const rows = data.map((g, n) => {
    const cells = d.groups.map((grp, gi) => {
      const v = d.samples.map((s, k) => (s.group === grp ? g.tpm[d.id][k] : null)).filter((x) => x != null);
      return `<td>${fmtTPM(mean(v))} <span class="dim small">${ordinal(g.pct[d.id][gi])}</span></td>`;
    });
    const name = Data.genes[selected[n]][1];
    return `<tr data-hl-gene="${n}"><td><a href="${genesLink([name])}" title="Show only this gene">${esc(name)}</a></td>${cells.join("")}<td class="left">${changesCell(d, g.de[d.id], n)}</td></tr>`;
  });
  const cols = "<col>".repeat(d.groups.length + 2);   // gene, one per group, changes
  $(`tab-${d.id}`).innerHTML = `<div class="scroll"><table><caption>Mean TPM, with its percentile among protein-coding genes; changes at p<sub>adj</sub> &lt; 0.05</caption>
    <colgroup>${cols}</colgroup><thead>${head}</thead><tbody>${rows.join("")}</tbody></table></div>`;
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
