// DE results: volcano, sortable table and cross-comparison agreement for one
// DESeq2 comparison at a time. Data: data/de/<dataset>__<contrast>.json (columns).
const UP = "#D55E00", DOWN = "#0072B2", NS = "#c9c9c9";
const PAGE = 50;
const cache = new Map();
const state = { c: null, rows: [], sort: { key: "rank", dir: 1 }, page: 0, pick: null };

const alpha = () => +$("alpha").value || 0.05;
const lfcMin = () => +$("lfcmin").value || 0;
const isSig = (r) => r.padj != null && r.padj < alpha() && Math.abs(r.lfc) >= lfcMin();
const dsOf = (c) => c.split("__")[0];
const fullLabel = (c) => `${SHORT_DS[dsOf(c)]}: ${contrastLabel(c.split("__")[1])}`;

async function load(c) {
  if (!cache.has(c)) {
    cache.set(c, fetch(`data/de/${c}.json`).then((r) => r.json()).then((d) => {
      const rows = d.g.map((g, k) => ({ i: g, name: Data.genes[g][1], base: d.base[k], lfc: d.lfc[k], mle: d.mle[k], p: d.p[k], padj: d.padj[k] }));
      const ranked = rows.filter((r) => r.p != null).sort((a, b) => a.p - b.p);
      ranked.forEach((r, k) => (r.rank = k + 1));
      rows.nTested = ranked.length;
      rows.byGene = new Map(rows.map((r) => [r.i, r]));
      return rows;
    }));
  }
  return cache.get(c);
}

async function init() {
  await Data.load();
  const opts = Data.meta.datasets.map((d) => `<optgroup label="${d.label}">${d.contrasts.map((c) =>
    `<option value="${d.id}__${c}">${contrastLabel(c)}</option>`).join("")}</optgroup>`).join("");
  $("contrast").innerHTML = opts;
  $("ax").innerHTML = opts; $("ay").innerHTML = opts;
  $("ax").value = "JCI_tissue__LMmut_vs_MM"; $("ay").value = "NC_2D__Mut1_vs_WT";

  const h = new URLSearchParams(location.hash.slice(1));
  if (h.get("c") && [...$("contrast").options].some((o) => o.value === h.get("c"))) $("contrast").value = h.get("c");
  const g = h.get("gene") && Data.lookup(h.get("gene"))[0];
  if (g != null) state.pick = g;

  $("contrast").addEventListener("change", () => showContrast());
  for (const id of ["alpha", "lfcmin"]) $(id).addEventListener("input", () => { refresh(); drawAgree(); });
  $("filter").addEventListener("input", () => { state.page = 0; drawTable(); });
  $("sigonly").addEventListener("change", () => { state.page = 0; drawTable(); });
  $("prev").addEventListener("click", () => { state.page--; drawTable(); });
  $("next").addEventListener("click", () => { state.page++; drawTable(); });
  $("download").addEventListener("click", download);
  $("detable").addEventListener("click", (e) => {
    const th = e.target.closest("th[data-key]");
    if (th) {
      const k = th.dataset.key;
      state.sort = { key: k, dir: state.sort.key === k ? -state.sort.dir : (k === "rank" || k === "name" ? 1 : -1) };
      state.page = 0; drawTable(); return;
    }
    const tr = e.target.closest("tr[data-i]");
    if (tr && !e.target.closest("a")) pick(+tr.dataset.i, false);
  });
  attachSearch($("find"), $("suggest"), (i) => { $("find").value = ""; pick(i, true); },
    (q) => { $("picked").hidden = false; $("picked").textContent = `No gene called “${q}”.`; });
  for (const id of ["ax", "ay"]) $(id).addEventListener("change", drawAgree);
  await showContrast();
  drawAgree();
}

async function showContrast() {
  state.c = $("contrast").value;
  state.rows = await load(state.c);
  state.page = 0;
  history.replaceState(null, "", `#c=${state.c}` + (state.pick != null ? `&gene=${encodeURIComponent(Data.genes[state.pick][1])}` : ""));
  refresh();
  if (state.pick != null) pick(state.pick, true);
}

function refresh() {
  const rows = state.rows;
  const up = rows.filter((r) => isSig(r) && r.lfc > 0).length, down = rows.filter((r) => isSig(r) && r.lfc < 0).length;
  $("summary").innerHTML = `<b>${fullLabel(state.c)}</b> · ${rows.nTested.toLocaleString()} genes tested ·
    <span style="color:${UP}">▲</span> ${up.toLocaleString()} up · <span style="color:${DOWN}">▼</span> ${down.toLocaleString()} down`;
  drawVolcano();
  drawTable();
}

function drawVolcano() {
  const rows = state.rows.filter((r) => r.padj != null && r.lfc != null);
  const finite = rows.filter((r) => r.padj > 0).map((r) => -Math.log10(r.padj));
  const cap = Math.max(10, ...finite);
  const yOf = (r) => (r.padj > 0 ? -Math.log10(r.padj) : cap);
  const groups = [
    { name: "not significant", color: NS, f: (r) => !isSig(r) },
    { name: "down", color: DOWN, f: (r) => isSig(r) && r.lfc < 0 },
    { name: "up", color: UP, f: (r) => isSig(r) && r.lfc > 0 },
  ];
  const traces = groups.map((g) => {
    const rs = rows.filter(g.f);
    return { type: "scattergl", mode: "markers", name: g.name, x: rs.map((r) => r.lfc), y: rs.map(yOf),
      customdata: rs.map((r) => r.i), text: rs.map((r) => r.name),
      marker: { color: g.color, size: g.color === NS ? 4 : 5, symbol: rs.map((r) => (r.padj > 0 ? "circle" : "triangle-up")) },
      hovertemplate: "<b>%{text}</b><br>log₂FC %{x:.2f}<extra></extra>" };
  });
  // label the 10 most significant genes (none on narrow screens, where they pile up), plus the picked one
  const nLabels = $("volcano").clientWidth < 600 ? 0 : 10;
  const top = rows.filter(isSig).sort((a, b) => a.padj - b.padj || Math.abs(b.lfc) - Math.abs(a.lfc)).slice(0, nLabels);
  const annotations = top.map((r) => ({ x: r.lfc, y: yOf(r), text: esc(r.name), showarrow: true, arrowhead: 0, arrowwidth: 1,
    arrowcolor: "#999", ax: r.lfc > 0 ? 24 : -24, ay: -14, font: { size: 11, color: "#1a1a1a" } }));
  const p = state.pick != null && state.rows.byGene.get(state.pick);
  if (p && p.padj != null) {
    traces.push({ type: "scatter", mode: "markers", x: [p.lfc], y: [yOf(p)], hoverinfo: "skip", showlegend: false,
      marker: { size: 16, color: "rgba(0,0,0,0)", line: { color: "#1a1a1a", width: 2 } } });
    annotations.push({ x: p.lfc, y: yOf(p), text: `<b>${esc(p.name)}</b>`, showarrow: true, arrowhead: 0, ax: 0, ay: -30,
      bgcolor: "#fff", font: { size: 13 } });
  }
  const t = lfcMin(), shapes = [
    { type: "line", xref: "paper", x0: 0, x1: 1, y0: -Math.log10(alpha()), y1: -Math.log10(alpha()), line: { color: "#999", width: 1, dash: "dot" } },
    ...(t > 0 ? [-t, t].map((x) => ({ type: "line", yref: "paper", y0: 0, y1: 1, x0: x, x1: x, line: { color: "#999", width: 1, dash: "dot" } })) : []),
  ];
  Plotly.react("volcano", traces, {
    font: PLOT_FONT, margin: { l: 56, r: 12, t: 30, b: 48 }, shapes, annotations,
    legend: { orientation: "h", x: 0, y: 1.02, yanchor: "bottom", traceorder: "reversed" },
    xaxis: { title: { text: "log₂ fold change (shrunken)" }, zeroline: true, zerolinecolor: "#ddd", gridcolor: "#f0f0f0", ticks: "outside", tickcolor: "#bbb" },
    yaxis: { title: { text: "−log₁₀ p<sub>adj</sub>" }, rangemode: "tozero", gridcolor: "#f0f0f0", zeroline: false, ticks: "outside", tickcolor: "#bbb" },
    plot_bgcolor: "#fff", paper_bgcolor: "#fff", hovermode: "closest", hoverlabel: { bgcolor: "#fff", font: { color: "#1a1a1a" } },
  }, PLOT_CONFIG);
  const el = $("volcano");
  el.removeAllListeners?.("plotly_click");
  el.on("plotly_click", (ev) => { const i = ev.points[0]?.customdata; if (i != null) pick(i, false); });
}

const COLS = [
  { key: "rank", label: "Rank" }, { key: "name", label: "Gene" }, { key: "base", label: "baseMean" },
  { key: "lfc", label: "log₂FC" }, { key: "mle", label: "log₂FC (unshrunken)" }, { key: "padj", label: "p<sub>adj</sub>" },
];

function filtered() {
  const q = $("filter").value.trim().toUpperCase();
  let rows = state.rows;
  if (q) rows = rows.filter((r) => r.name.toUpperCase().includes(q));
  if ($("sigonly").checked) rows = rows.filter(isSig);
  const { key, dir } = state.sort;
  return [...rows].sort((a, b) => {
    const x = a[key], y = b[key];
    if (x == null) return 1; if (y == null) return -1;   // empties last either way
    return (typeof x === "string" ? x.localeCompare(y) : x - y) * dir;
  });
}

function drawTable() {
  const rows = filtered();
  const pages = Math.max(1, Math.ceil(rows.length / PAGE));
  state.page = Math.min(Math.max(0, state.page), pages - 1);
  const head = `<tr>${COLS.map((c) => `<th data-key="${c.key}" class="${state.sort.key === c.key ? (state.sort.dir > 0 ? "asc" : "desc") : ""}">${c.label}</th>`).join("")}</tr>`;
  const body = rows.slice(state.page * PAGE, (state.page + 1) * PAGE).map((r) => {
    const sig = isSig(r), color = sig ? (r.lfc > 0 ? UP : DOWN) : null;
    return `<tr data-i="${r.i}" class="${r.i === state.pick ? "picked-row" : ""}">
      <td>${r.rank ?? "—"}</td><td><a href="${explorerLink([r.name])}">${esc(r.name)}</a></td>
      <td>${Math.round(r.base).toLocaleString()}</td>
      <td${color ? ` style="color:${color};font-weight:700"` : ""}>${fmt(r.lfc)}</td><td>${fmt(r.mle)}</td>
      <td class="${sig ? "sig" : ""}">${fmtP(r.padj)}</td></tr>`;
  }).join("");
  $("detable").innerHTML = head + body;
  $("count").textContent = `${rows.length.toLocaleString()} genes`;
  $("page").textContent = `page ${state.page + 1} of ${pages}`;
  $("prev").disabled = state.page === 0;
  $("next").disabled = state.page >= pages - 1;
}

// select a gene: describe its rank, ring it on the volcano, mark its row
function pick(i, jumpTable) {
  state.pick = i;
  const name = Data.genes[i][1];
  history.replaceState(null, "", `#c=${state.c}&gene=${encodeURIComponent(name)}`);
  const r = state.rows.byGene.get(i);
  $("picked").hidden = false;
  if (!r) {
    $("picked").innerHTML = `<b>${esc(name)}</b> was not tested in this comparison (too few counts). <a href="${explorerLink([name])}">Open in explorer</a>`;
  } else {
    const pct = r.rank ? ` (top ${Math.max(1, Math.ceil((r.rank / state.rows.nTested) * 100))}%)` : "";
    $("picked").innerHTML = `<b>${esc(name)}</b>: rank ${r.rank ? r.rank.toLocaleString() : "—"} of ${state.rows.nTested.toLocaleString()}${pct} ·
      log₂FC ${fmt(r.lfc)} · p<sub>adj</sub> ${fmtP(r.padj)} · baseMean ${Math.round(r.base).toLocaleString()}
      ${isSig(r) ? `· <span style="color:${r.lfc > 0 ? UP : DOWN};font-weight:700">${r.lfc > 0 ? "up" : "down"}</span>` : "· not significant"}
      · <a href="${explorerLink([name])}">Open in explorer</a>`;
    if (jumpTable) {
      $("filter").value = ""; $("sigonly").checked = false;
      state.sort = { key: "rank", dir: 1 };
      const k = filtered().indexOf(r);
      state.page = Math.floor(k / PAGE);
    }
  }
  drawVolcano();
  drawTable();
}

function download() {
  const rows = filtered();
  const lines = ["gene_id,gene_name,gene_type,baseMean,log2FC_shrunk,log2FC,pvalue,padj,rank"];
  for (const r of rows) {
    const [id, name, type] = Data.genes[r.i];
    lines.push([id, name, type, r.base, r.lfc, r.mle, r.p, r.padj, r.rank].map((v) => v ?? "").join(","));
  }
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([lines.join("\n") + "\n"], { type: "text/csv" }));
  a.download = `${state.c}${$("sigonly").checked ? "_significant" : ""}.csv`;
  a.click();
  URL.revokeObjectURL(a.href);
}

// Spearman correlation: Pearson on average ranks
function spearman(x, y) {
  const rank = (v) => {
    const o = v.map((x, k) => [x, k]).sort((a, b) => a[0] - b[0]), r = new Array(v.length);
    for (let s = 0; s < o.length;) {
      let e = s; while (e + 1 < o.length && o[e + 1][0] === o[s][0]) e++;
      for (let k = s; k <= e; k++) r[o[k][1]] = (s + e) / 2;
      s = e + 1;
    }
    return r;
  };
  const a = rank(x), b = rank(y), n = a.length, ma = (n - 1) / 2;
  let sab = 0, saa = 0, sbb = 0;
  for (let k = 0; k < n; k++) { const da = a[k] - ma, db = b[k] - ma; sab += da * db; saa += da * da; sbb += db * db; }
  return sab / Math.sqrt(saa * sbb);
}

async function drawAgree() {
  const [cx, cy] = [$("ax").value, $("ay").value];
  const [X, Y] = await Promise.all([load(cx), load(cy)]);
  const pairs = [];
  for (const rx of X) {
    const ry = Y.byGene.get(rx.i);
    if (ry && rx.lfc != null && ry.lfc != null) pairs.push([rx, ry]);
  }
  const cats = [
    { name: "neither significant", color: "#dddddd", f: (a, b) => !isSig(a) && !isSig(b) },
    { name: "significant in one", color: "#9aa5b8", f: (a, b) => isSig(a) !== isSig(b) },
    { name: "both, opposite direction", color: DOWN, f: (a, b) => isSig(a) && isSig(b) && Math.sign(a.lfc) !== Math.sign(b.lfc) },
    { name: "both, same direction", color: UP, f: (a, b) => isSig(a) && isSig(b) && Math.sign(a.lfc) === Math.sign(b.lfc) },
  ];
  const counts = cats.map((c) => pairs.filter(([a, b]) => c.f(a, b)).length);
  const traces = cats.map((c) => {
    const ps = pairs.filter(([a, b]) => c.f(a, b));
    return { type: "scattergl", mode: "markers", name: `${c.name} (${ps.length.toLocaleString()})`,
      x: ps.map(([a]) => a.lfc), y: ps.map(([, b]) => b.lfc), text: ps.map(([a]) => a.name), customdata: ps.map(([a]) => a.i),
      marker: { color: c.color, size: c.color === "#dddddd" ? 4 : 5 },
      hovertemplate: "<b>%{text}</b><br>x %{x:.2f} · y %{y:.2f}<extra></extra>" };
  });
  const rho = spearman(pairs.map(([a]) => a.lfc), pairs.map(([, b]) => b.lfc));
  $("agree-summary").innerHTML = `${pairs.length.toLocaleString()} genes tested in both · Spearman ρ ${rho.toFixed(2)} ·
    significant in both: ${counts[3].toLocaleString()} same direction, ${counts[2].toLocaleString()} opposite`;
  const lim = Math.max(1, ...pairs.flatMap(([a, b]) => [Math.abs(a.lfc), Math.abs(b.lfc)]));
  Plotly.react("agree", traces, {
    font: PLOT_FONT, margin: { l: 56, r: 12, t: 30, b: 48 },
    legend: { orientation: "h", x: 0, y: 1.02, yanchor: "bottom", traceorder: "reversed" },
    shapes: [{ type: "line", x0: -lim, y0: -lim, x1: lim, y1: lim, line: { color: "#bbb", width: 1, dash: "dot" } }],
    xaxis: { title: { text: `${fullLabel(cx)}, log₂FC` }, zeroline: true, zerolinecolor: "#ddd", gridcolor: "#f0f0f0", ticks: "outside", tickcolor: "#bbb" },
    yaxis: { title: { text: `${fullLabel(cy)}, log₂FC` }, zeroline: true, zerolinecolor: "#ddd", gridcolor: "#f0f0f0", ticks: "outside", tickcolor: "#bbb" },
    plot_bgcolor: "#fff", paper_bgcolor: "#fff", hovermode: "closest", hoverlabel: { bgcolor: "#fff", font: { color: "#1a1a1a" } },
  }, PLOT_CONFIG);
  const el = $("agree");
  el.removeAllListeners?.("plotly_click");
  el.on("plotly_click", (ev) => { const i = ev.points[0]?.customdata; if (i != null) pick(i, true); });
}

init().catch((e) => { $("summary").textContent = "Could not load data: " + e.message; });
