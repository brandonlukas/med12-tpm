// Comparisons page: volcano, ranked table and agreement between two comparisons.
// Data: data/de/<dataset>__<contrast>.json, one comparison as columns.
const NS = "#cfd2db";
// agreement categories: neither a direction nor a group hue
const AGREE = { same: "#2b3a8c", opposite: "#009E73", one: "#a3abb9", neither: "#e1e4ea" };
const PAGE = 25;
const cache = new Map();
const state = { c: null, rows: [], sort: { key: "rank", dir: 1 }, page: 0, pick: null };
let noPlots = false; // the plotting library did not load; the table still works

const alpha = () => +$("alpha").value || 0.05;
const lfcMin = () => +$("lfcmin").value || 0;
const isSig = (r) => r.padj != null && r.padj < alpha() && Math.abs(r.lfc) >= lfcMin();
const dsOf = (c) => c.split("__")[0];
const conOf = (c) => c.split("__")[1];
const fullLabel = (c) => contrastOption(dsOf(c), conOf(c));
const sides = (c) => contrastSides(conOf(c));
const nlog = (p) => -Math.log10(p);

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
  const [plotly] = await Promise.all([ensurePlotly("vendor/plotly-gl2d-2.35.2.min.js"), Data.load()]);
  noPlots = !plotly;
  const opts = Data.meta.datasets.map((d) => `<optgroup label="${esc(d.label)}">${d.contrasts.map((c) =>
    `<option value="${d.id}__${c}">${esc(contrastOption(d.id, c))}</option>`).join("")}</optgroup>`).join("");
  $("contrast").innerHTML = opts;
  $("ax").innerHTML = opts; $("ay").innerHTML = opts;
  $("ax").value = "JCI_tissue__LMmut_vs_MM"; $("ay").value = "NC_2D__Mut1_vs_WT";

  const h = new URLSearchParams(location.hash.slice(1));
  if (h.get("c") && [...$("contrast").options].some((o) => o.value === h.get("c"))) $("contrast").value = h.get("c");
  const g = h.get("gene") && Data.lookup(h.get("gene"))[0];
  if (g != null) state.pick = g;

  $("contrast").addEventListener("change", showContrast);
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
    if (tr) pick(+tr.dataset.i, false);
  });
  $("detable").addEventListener("keydown", (e) => {
    const tr = e.target.closest("tr[data-i]");
    if (tr && (e.key === "Enter" || e.key === " ")) { e.preventDefault(); pick(+tr.dataset.i, false); }
  });
  $("picked").addEventListener("click", (e) => { if (e.target.closest("[data-clear]")) clearPick(); });
  attachSearch($("find"), $("suggest"), (i) => { $("find").value = ""; pick(i, true); }, (q) => {
    $("picked").hidden = false;
    $("picked").innerHTML = `<span>No gene called “${esc(q)}”. Try its official symbol, for example HMGA2.</span>
      <span class="end"><button type="button" class="btn quiet" data-clear>Clear</button></span>`;
  });
  for (const id of ["ax", "ay"]) $(id).addEventListener("change", drawAgree);
  window.addEventListener("resize", () => { if (state.rows.length) drawVolcano(); });
  await showContrast();
  drawAgree();
}

async function showContrast() {
  state.c = $("contrast").value;
  state.rows = await load(state.c);
  state.page = 0;
  setHash();
  refresh();
  if (state.pick != null) pick(state.pick, true);
}

function setHash() {
  const gene = state.pick != null ? `&gene=${encodeURIComponent(Data.genes[state.pick][1])}` : "";
  history.replaceState(null, "", `#c=${state.c}${gene}`);
}

function refresh() {
  const rows = state.rows, [num] = sides(state.c);
  const up = rows.filter((r) => isSig(r) && r.lfc > 0).length, down = rows.filter((r) => isSig(r) && r.lfc < 0).length;
  $("summary").innerHTML = `<b>${esc(fullLabel(state.c))}</b> <span>${rows.nTested.toLocaleString()} genes tested</span>
    <span class="up">▲ ${up.toLocaleString()} higher in ${esc(num)}</span> <span class="down">▼ ${down.toLocaleString()} lower in ${esc(num)}</span>`;
  drawVolcano();
  drawTable();
}

// Greedy label placement in pixel space: candidates in priority order, each label at the
// first of a few offsets where it overlaps neither a placed label (with padding) nor the frame.
function placeLabels(cands, xr, yr, W, H, m, max = 14) {
  const pw = W - m.l - m.r, ph = H - m.t - m.b, pad = 4;
  const px = (x) => m.l + ((x - xr[0]) / (xr[1] - xr[0])) * pw;
  const py = (y) => m.t + (1 - (y - yr[0]) / (yr[1] - yr[0])) * ph;
  const placed = [], out = [];
  const tries = [[14, -12], [-14, -12], [16, 2], [-16, 2], [14, 14], [-14, 14], [0, -20], [0, 20]];
  for (const c of cands) {
    const w = c.name.length * 6.8 + 6, h = 15, cx = px(c.x), cy = py(c.y);
    for (const [dx, dy] of tries) {
      const left = dx > 0 ? cx + dx : dx < 0 ? cx + dx - w : cx - w / 2;
      const top = cy + dy - h / 2;
      if (left < m.l || left + w > W - m.r || top < m.t - 6 || top + h > H - m.b) continue;
      if (placed.some((b) => !(left + w + pad < b.left || b.left + b.w + pad < left || top + h + pad < b.top || b.top + b.h + pad < top))) continue;
      placed.push({ left, top, w, h });
      out.push({ ...c, ax: dx, ay: dy, xanchor: dx > 0 ? "left" : dx < 0 ? "right" : "center" });
      break;
    }
    if (out.length >= max) break;
  }
  return out;
}

function drawVolcano() {
  const el = $("volcano");
  if (noPlots) { el.innerHTML = NO_PLOT; return; }
  const rows = state.rows.filter((r) => r.padj != null && r.lfc != null);
  const [num, den] = sides(state.c);
  // Cap the y axis at 2.5x the 95th percentile of the significant genes (padj < 0.05) so that a
  // few genes at p ~ 1e-300 do not squash the rest; genes above the cap are drawn as triangles on it.
  const sv = rows.filter((r) => r.padj > 0 && r.padj < 0.05).map((r) => nlog(r.padj)).sort((a, b) => a - b);
  const vmax = Math.max(10, ...rows.filter((r) => r.padj > 0).map((r) => nlog(r.padj)));
  const cap = sv.length ? Math.min(vmax, Math.max(10, 2.5 * sv[Math.floor(sv.length * 0.95)])) : vmax;
  const capped = (r) => r.padj === 0 || nlog(r.padj) > cap;
  const yOf = (r) => (capped(r) ? cap : nlog(r.padj));
  const xmax = Math.max(1, ...rows.map((r) => Math.abs(r.lfc))) * 1.06;
  const xr = [-xmax, xmax], yr = [0, cap * 1.08];
  const margin = { l: 58, r: 14, t: 36, b: 52 };
  const colorOf = (r) => (isSig(r) ? (r.lfc > 0 ? UP : DOWN) : NS);
  const trace = (rs, extra) => ({
    type: "scattergl", mode: "markers", x: rs.map((r) => r.lfc), y: rs.map(yOf), text: rs.map((r) => r.name),
    customdata: rs.map((r) => [r.i, fmtP(r.padj)]),
    hovertemplate: "<b>%{text}</b><br>log₂FC %{x:.2f} · p<sub>adj</sub> %{customdata[1]}<extra></extra>", ...extra,
  });
  const plain = rows.filter((r) => !capped(r)), above = rows.filter(capped);
  const traces = [
    trace(plain.filter((r) => !isSig(r)), { name: "not significant", marker: { color: NS, size: 4 } }),
    trace(plain.filter((r) => isSig(r) && r.lfc < 0), { name: `lower in ${num}`, marker: { color: DOWN, size: 5.5 } }),
    trace(plain.filter((r) => isSig(r) && r.lfc > 0), { name: `higher in ${num}`, marker: { color: UP, size: 5.5 } }),
  ];
  if (above.length) traces.push(trace(above, { showlegend: false, marker: { color: above.map(colorOf), size: 7, symbol: "triangle-up" } }));

  // labels: the extremes first (capped genes by fold change, then the largest fold changes on
  // each side), then the most significant; the picked gene gets its own label
  const sig = rows.filter((r) => isSig(r) && r.i !== state.pick);
  const seen = new Set(), cands = [];
  for (const r of [
    ...sig.filter(capped).sort((a, b) => Math.abs(b.lfc) - Math.abs(a.lfc)).slice(0, 4),
    ...sig.filter((r) => r.lfc > 0).sort((a, b) => b.lfc - a.lfc).slice(0, 3),
    ...sig.filter((r) => r.lfc < 0).sort((a, b) => a.lfc - b.lfc).slice(0, 3),
    ...sig.filter((r) => !capped(r)).sort((a, b) => a.padj - b.padj).slice(0, 24),
  ]) if (!seen.has(r.i)) { seen.add(r.i); cands.push({ name: r.name, x: r.lfc, y: yOf(r) }); }
  const labels = el.clientWidth >= 600 ? placeLabels(cands, xr, yr, el.clientWidth, el.clientHeight, margin) : [];
  const annotations = labels.map((l) => ({ x: l.x, y: l.y, ax: l.ax, ay: l.ay, xanchor: l.xanchor, text: esc(l.name),
    showarrow: true, arrowhead: 0, arrowwidth: 1, arrowcolor: "#b4b8c4", font: { size: 11.5, color: INK } }));

  const p = state.pick != null && state.rows.byGene.get(state.pick);
  if (p && p.padj != null) {
    traces.push(ring(p.lfc, yOf(p)));
    annotations.push(ringLabel(p.name, p.lfc, yOf(p)));
  }
  const t = lfcMin(), a = nlog(alpha());
  const guide = { color: "#b4b8c4", width: 1, dash: "dot" };
  const shapes = [
    { type: "line", xref: "paper", x0: 0, x1: 1, y0: a, y1: a, line: guide },
    ...(t > 0 ? [-t, t].map((x) => ({ type: "line", yref: "paper", y0: 0, y1: 1, x0: x, x1: x, line: guide })) : []),
  ];
  Plotly.react(el, traces, {
    ...BASE_LAYOUT, margin, shapes, annotations, hovermode: "closest",
    legend: { orientation: "h", x: 0, y: 1.02, yanchor: "bottom", traceorder: "reversed", itemclick: false, itemdoubleclick: false },
    xaxis: { ...AXIS, fixedrange: false, range: xr, title: { text: `log₂ fold change, ${esc(num)} vs ${esc(den)}`, standoff: 8 },
      zeroline: true, zerolinecolor: RULE },
    yaxis: { ...AXIS, fixedrange: false, range: yr, title: { text: "−log₁₀ p<sub>adj</sub>", standoff: 8 } },
  }, PLOT_CONFIG);
  el._home = { x: [...xr], y: [...yr] };   // copies: Plotly edits the range arrays it was given when the user zooms
  zoomReset(el, $("volcano-reset"), drawVolcano);
  el.removeAllListeners?.("plotly_click");
  el.on("plotly_click", (ev) => { const i = ev.points[0]?.customdata?.[0]; if (i != null) pick(i, false); });
}

const ring = (x, y) => ({ type: "scatter", mode: "markers", x: [x], y: [y], hoverinfo: "skip", showlegend: false,
  marker: { size: 17, color: "rgba(0,0,0,0)", line: { color: INK, width: 2 } } });
const ringLabel = (name, x, y) => ({ x, y, text: `<b>${esc(name)}</b>`, showarrow: true, arrowhead: 0, arrowcolor: INK,
  ax: 0, ay: -34, bgcolor: "#fff", bordercolor: INK, borderpad: 3, font: { size: 13, color: INK } });

const COLS = [
  { key: "rank", label: "Rank" }, { key: "name", label: "Gene" },
  { key: "base", label: "Mean count", title: "Average normalized count over all samples in the dataset" },
  { key: "lfc", label: "log₂FC" }, { key: "padj", label: "p<sub>adj</sub>" },
];

// the tested genes, filtered and sorted as the table shows them
function filtered() {
  const q = $("filter").value.trim().toUpperCase();
  let rows = state.rows.filter((r) => r.p != null);
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
  const head = `<thead><tr>${COLS.map((c) => `<th scope="col" data-key="${c.key}" ${c.title ? `title="${c.title}"` : ""}
    class="${state.sort.key === c.key ? (state.sort.dir > 0 ? "asc" : "desc") : ""}">${c.label}</th>`).join("")}</tr></thead>`;
  const body = rows.slice(state.page * PAGE, (state.page + 1) * PAGE).map((r) => {
    const sig = isSig(r), dir = sig ? (r.lfc > 0 ? "up" : "down") : "";
    return `<tr data-i="${r.i}" tabindex="0" class="${r.i === state.pick ? "picked-row" : ""}">
      <td>${r.rank.toLocaleString()}</td><td>${esc(r.name)}</td>
      <td>${Math.round(r.base).toLocaleString()}</td>
      <td class="${dir ? `sig ${dir}` : ""}">${fmt(r.lfc)}</td>
      <td class="${sig ? "sig" : ""}">${fmtP(r.padj)}</td></tr>`;
  }).join("");
  $("detable").innerHTML = head + `<tbody>${body}</tbody>`;
  $("count").textContent = `${rows.length.toLocaleString()} genes`;
  $("page").textContent = `Page ${state.page + 1} of ${pages}`;
  $("prev").disabled = state.page === 0;
  $("next").disabled = state.page >= pages - 1;
}

// select a gene: say where it ranks, ring it on the plots, mark its row
function pick(i, jumpTable) {
  state.pick = i;
  setHash();
  const name = Data.genes[i][1], r = state.rows.byGene.get(i), [num] = sides(state.c);
  const link = `<a href="${genesLink([name])}">Show expression</a>`;
  const clear = `<button type="button" class="btn quiet" data-clear>Clear</button>`;
  $("picked").hidden = false;
  if (!r || r.p == null) {
    $("picked").innerHTML = `<b>${esc(name)}</b> <span>not tested in this comparison: too few counts</span>
      <span class="end">${link} ${clear}</span>`;
  } else {
    const pct = Math.max(1, Math.ceil((r.rank / state.rows.nTested) * 100));
    const call = r.padj == null ? `<span class="dim">no p<sub>adj</sub></span>`
      : isSig(r) ? `<span class="${r.lfc > 0 ? "up" : "down"}"><b>${r.lfc > 0 ? "higher" : "lower"} in ${esc(num)}</b></span>`
      : `<span class="dim">not significant</span>`;
    $("picked").innerHTML = `<b>${esc(name)}</b>
      <span>ranked ${r.rank.toLocaleString()} of ${state.rows.nTested.toLocaleString()} by p-value, top ${pct}%</span>
      <span>log₂FC ${fmt(r.lfc)}</span> <span>p<sub>adj</sub> ${fmtP(r.padj)}</span> ${call}
      <span class="end">${link} ${clear}</span>`;
    if (jumpTable) {
      $("filter").value = ""; $("sigonly").checked = false;
      state.sort = { key: "rank", dir: 1 };
      state.page = Math.floor(filtered().indexOf(r) / PAGE);
    }
  }
  drawVolcano();
  drawTable();
  drawAgree();
}

function clearPick() {
  state.pick = null;
  $("picked").hidden = true;
  setHash();
  drawVolcano();
  drawTable();
  drawAgree();
}

function download() {
  const rows = filtered();
  const lines = ["gene_id,gene_name,gene_type,baseMean,log2FC_shrunk,log2FC,pvalue,padj,rank"];
  for (const r of rows) {
    const [id, name, type] = Data.genes[r.i];
    lines.push([id, name, type, r.base, r.lfc, r.mle, r.p, r.padj, r.rank].map((v) => v ?? "").join(","));
  }
  saveText(lines.join("\n") + "\n", `${state.c}${$("sigonly").checked ? "_significant" : ""}.csv`);
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
  if (noPlots) { $("agree").innerHTML = NO_PLOT; $("agree-summary").textContent = ""; return; }
  const [cx, cy] = [$("ax").value, $("ay").value];
  const [X, Y] = await Promise.all([load(cx), load(cy)]);
  if (cx !== $("ax").value || cy !== $("ay").value) return;
  const pairs = [];
  for (const rx of X) {
    const ry = Y.byGene.get(rx.i);
    if (ry && rx.lfc != null && ry.lfc != null) pairs.push([rx, ry]);
  }
  const cats = [
    { name: "neither significant", color: AGREE.neither, f: (a, b) => !isSig(a) && !isSig(b) },
    { name: "significant in one", color: AGREE.one, f: (a, b) => isSig(a) !== isSig(b) },
    { name: "both, opposite directions", color: AGREE.opposite, f: (a, b) => isSig(a) && isSig(b) && Math.sign(a.lfc) !== Math.sign(b.lfc) },
    { name: "both, same direction", color: AGREE.same, f: (a, b) => isSig(a) && isSig(b) && Math.sign(a.lfc) === Math.sign(b.lfc) },
  ];
  const sets = cats.map((c) => pairs.filter(([a, b]) => c.f(a, b)));
  const xname = SHORT_DS[dsOf(cx)], yname = SHORT_DS[dsOf(cy)];
  const traces = cats.map((c, k) => {
    const ps = sets[k];
    return { type: "scattergl", mode: "markers", name: `${c.name} (${ps.length.toLocaleString()})`,
      x: ps.map(([a]) => a.lfc), y: ps.map(([, b]) => b.lfc), text: ps.map(([a]) => a.name), customdata: ps.map(([a]) => a.i),
      marker: { color: c.color, size: k === 0 ? 4 : 5.5 },
      hovertemplate: `<b>%{text}</b><br>${esc(xname)} %{x:.2f} · ${esc(yname)} %{y:.2f}<extra></extra>` };
  });
  const rho = pairs.length > 2 ? spearman(pairs.map(([a]) => a.lfc), pairs.map(([, b]) => b.lfc)) : NaN;
  $("agree-summary").innerHTML = cx === cy
    ? `The same comparison on both axes. Pick two different ones to compare them.`
    : `${pairs.length.toLocaleString()} genes tested in both · rank correlation ${rho.toFixed(2)} · significant in both:
       <b>${sets[3].length.toLocaleString()}</b> in the same direction, <b>${sets[2].length.toLocaleString()}</b> in opposite directions`;
  const lim = Math.max(1, ...pairs.flatMap(([a, b]) => [Math.abs(a.lfc), Math.abs(b.lfc)])) * 1.05;
  const annotations = [];
  const pk = state.pick != null && pairs.find(([a]) => a.i === state.pick);
  if (pk) { traces.push(ring(pk[0].lfc, pk[1].lfc)); annotations.push(ringLabel(pk[0].name, pk[0].lfc, pk[1].lfc)); }
  const el = $("agree");
  const axis = (title) => ({ ...AXIS, fixedrange: false, range: [-lim, lim], zeroline: true, zerolinecolor: RULE, title: { text: title, standoff: 8 } });
  Plotly.react(el, traces, {
    ...BASE_LAYOUT, margin: { l: 60, r: 14, t: 56, b: 56 }, hovermode: "closest", annotations,
    legend: { orientation: "h", x: 0, y: 1.02, yanchor: "bottom", traceorder: "reversed", itemclick: false, itemdoubleclick: false },
    shapes: [{ type: "line", x0: -lim, y0: -lim, x1: lim, y1: lim, line: { color: "#b4b8c4", width: 1, dash: "dot" } }],
    xaxis: axis(`${esc(fullLabel(cx))}, log₂FC`),
    yaxis: { ...axis(`${esc(fullLabel(cy))}, log₂FC`), scaleanchor: "x", scaleratio: 1 },
  }, PLOT_CONFIG);
  el._home = { x: [-lim, lim], y: [-lim, lim] };
  zoomReset(el, $("agree-reset"), drawAgree);
  el.removeAllListeners?.("plotly_click");
  el.on("plotly_click", (ev) => { const i = ev.points[0]?.customdata; if (i != null) pick(i, true); });
}

init().catch((e) => { $("summary").textContent = "The data did not load: " + e.message; });
