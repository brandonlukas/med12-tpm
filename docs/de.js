// Comparisons page: volcano, ranked table and agreement between two comparisons, for genes or
// for transcription-factor (TF) activity.
// Data: data/de/<dataset>__<contrast>.json (genes) and data/tf/<dataset>__<contrast>.json (TFs),
// one comparison as columns; data/tf/network.json, each TF's targets.
const NS = "#cfd2db";
// agreement categories: neither a direction nor a group hue
const AGREE = { same: "#2b3a8c", opposite: "#009E73", one: "#a3abb9", neither: "#e1e4ea" };
const PAGE = 25, TOP_TARGETS = 10;
const cache = new Map();
const state = { kind: "gene", c: null, rows: [], sort: { key: "rank", dir: 1 }, page: 0, pick: null, allTargets: false };
let noPlots = false; // the plotting library did not load; the table still works
let network = null;

const DEFAULTS = { alpha: SIG.padj, lfcmin: SIG.lfc };
// an empty or unreadable field falls back to its default; any number, including 0, is taken as typed
const num = (id) => { const v = parseFloat($(id).value); return Number.isFinite(v) ? v : DEFAULTS[id]; };
const alpha = () => num("alpha");
const lfcMin = () => num("lfcmin");
function resetThresholds() {
  $("alpha").value = DEFAULTS.alpha; $("lfcmin").value = DEFAULTS.lfcmin;
  refresh(); drawAgree();
}
// r.x is the plotted effect: log2 fold change for a gene, activity for a TF (which has no fold-change cut-off)
const isSig = (r) => r.padj != null && r.padj < alpha() && (r.tf || Math.abs(r.x) >= lfcMin());
const dsOf = (c) => c.split("__")[0];
const conOf = (c) => c.split("__")[1];
const fullLabel = (c) => contrastOption(dsOf(c), conOf(c));
const sides = (c) => contrastSides(dsOf(c), conOf(c));
const nlog = (p) => -Math.log10(p);
// the gene with exactly this symbol (Data.lookup ignores case)
const geneByName = (name) => Data.lookup(name).find((i) => Data.genes[i][1] === name) ?? null;
const hidePick = () => { $("picked").hidden = true; $("targets").hidden = true; };

// Rows of one comparison. A gene row's key is its gene index, a TF row's its name.
async function load(c, kind = state.kind) {
  const id = `${kind}:${c}`;
  if (!cache.has(id)) {
    cache.set(id, fetch(`data/${kind === "tf" ? "tf" : "de"}/${c}.json`).then((r) => r.json()).then((d) => {
      const rows = kind === "tf"
        ? d.tf.map((t, k) => ({ tf: true, key: t, name: t, x: d.act[k], n: d.n[k], p: d.p[k], padj: d.padj[k] }))
        : d.g.map((g, k) => ({ key: g, name: Data.genes[g][1], base: d.base[k], x: d.lfc[k], mle: d.mle[k], z: d.z[k], p: d.p[k], padj: d.padj[k] }));
      const ranked = rows.filter((r) => r.p != null).sort((a, b) => a.p - b.p);
      ranked.forEach((r, k) => (r.rank = k + 1));
      rows.nTested = ranked.length;
      rows.byKey = new Map(rows.map((r) => [r.key, r]));
      return rows;
    }));
  }
  return cache.get(id);
}

// TFs whose name contains what was typed: an exact match first, then those starting with it, then shortest
const TF_SEARCH = {
  suggest(q) {
    q = q.trim().toUpperCase();
    if (!q) return [];
    const hits = state.rows.filter((r) => r.name.toUpperCase().includes(q));
    const rank = (r) => (r.name.toUpperCase() === q ? -2 : r.name.toUpperCase().startsWith(q) ? -1 : 0);
    return hits.sort((a, b) => rank(a) - rank(b) || a.name.length - b.name.length).slice(0, 10).map((r) => r.key);
  },
  lookup: (q) => state.rows.find((r) => r.name.toUpperCase() === q.trim().toUpperCase())?.key,
  html: (key) => `<span>${esc(key)} <span class="dim small">${state.rows.byKey.get(key).n} targets measured</span></span>`,
};

// What differs between the gene view and the TF view. A gene's key is its gene index, a TF's its name.
const KINDS = {
  gene: {
    noun: "gene", verb: "tested", effect: "log₂FC", example: "HMGA2", lfcCut: true, ymin: 10, size: 5.5,
    up: (g) => `higher in ${g}`, down: (g) => `lower in ${g}`, axis: (a, b) => `log₂ fold change, ${a} vs ${b}`,
    miss: "Try its official symbol, for example HMGA2.", untested: "not tested in this comparison: too few counts",
    nameOf: (i) => Data.genes[i][1], keyOf: geneByName, keyAttr: (s) => +s, search: GENE_SEARCH,
    cols: [
      { key: "rank", label: "Rank" }, { key: "name", label: "Gene" },
      { key: "base", label: "Mean count", title: "Average normalized count over all samples in the dataset" },
      { key: "x", label: "log₂FC" }, { key: "padj", label: "p<sub>adj</sub>" },
    ],
    count: (r) => Math.round(r.base),
    csv: { head: "gene_id,gene_name,gene_type,baseMean,log2FC_shrunk,log2FC,pvalue,padj,rank", suffix: "",
      row: (r) => [...Data.genes[r.key], r.base, r.x, r.mle, r.p, r.padj, r.rank] },
  },
  tf: {
    noun: "TF", verb: "scored", effect: "activity", example: "ESR1", lfcCut: false, ymin: 3, size: 7,
    up: (g) => `more active in ${g}`, down: (g) => `less active in ${g}`, axis: (a, b) => `TF activity, ${a} vs ${b}`,
    miss: "TFs are named by their gene symbol, for example ESR1.", untested: "not scored in this comparison: fewer than 5 targets measured",
    nameOf: (t) => t, keyOf: (name) => (state.rows.byKey.has(name) ? name : null), keyAttr: (s) => s, search: TF_SEARCH,
    cols: [
      { key: "rank", label: "Rank" }, { key: "name", label: "TF" },
      { key: "n", label: "Targets", title: "Targets of the TF measured in this comparison" },
      { key: "x", label: "Activity" }, { key: "padj", label: "p<sub>adj</sub>" },
    ],
    count: (r) => r.n,
    csv: { head: "tf,activity,pvalue,padj,n_targets,rank", suffix: "_TF_activity", row: (r) => [r.name, r.x, r.p, r.padj, r.n, r.rank] },
  },
};
const K = () => KINDS[state.kind];
const isTF = () => state.kind === "tf";
const SEARCH = { suggest: (q) => K().search.suggest(q), lookup: (q) => K().search.lookup(q), html: (m) => K().search.html(m) };

async function init() {
  const [plotly] = await Promise.all([ensurePlotly("vendor/plotly-gl2d-2.35.2.min.js"), Data.load()]);
  noPlots = !plotly;
  const opts = Data.meta.datasets.map((d) => `<optgroup label="${esc(d.name)}">${d.contrasts.map(({ id: c }) =>
    `<option value="${d.id}__${c}">${esc(contrastOption(d.id, c))}</option>`).join("")}</optgroup>`).join("");
  $("contrast").innerHTML = opts;
  $("ax").innerHTML = opts; $("ay").innerHTML = opts;
  $("ax").value = "JCI_tissue__LMmut_vs_MM"; $("ay").value = "NC_2D__Mut1_vs_WT";

  const h = new URLSearchParams(location.hash.slice(1));
  if (h.get("c") && [...$("contrast").options].some((o) => o.value === h.get("c"))) $("contrast").value = h.get("c");
  if (h.get("view") === "tf" || h.get("tf")) state.kind = "tf";
  if (h.get("tf")) state.pick = h.get("tf");
  else { const g = h.get("gene") && Data.lookup(h.get("gene"))[0]; if (g != null) state.pick = g; }

  $("contrast").addEventListener("change", showContrast);
  document.querySelectorAll("[data-kind]").forEach((b) => b.addEventListener("click", () => setKind(b.dataset.kind)));
  for (const id of ["alpha", "lfcmin"]) $(id).addEventListener("input", () => { refresh(); drawAgree(); });
  $("thr-reset").addEventListener("click", resetThresholds);
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
    const tr = e.target.closest("tr[data-k]");
    if (tr) pick(K().keyAttr(tr.dataset.k), false);
  });
  $("detable").addEventListener("keydown", (e) => {
    const tr = e.target.closest("tr[data-k]");
    if (tr && (e.key === "Enter" || e.key === " ")) { e.preventDefault(); pick(K().keyAttr(tr.dataset.k), false); }
  });
  $("picked").addEventListener("click", (e) => { if (e.target.closest("[data-clear]")) clearPick(); });
  $("targets").addEventListener("click", (e) => {
    if (e.target.closest("[data-all]")) { state.allTargets = !state.allTargets; drawTargets(); return; }
    const t = e.target.closest("[data-gene]");
    if (t) setKind("gene", +t.dataset.gene);
  });
  attachSearch($("find"), $("suggest"), (key) => { $("find").value = ""; pick(key, true); }, (q) => {
    $("picked").hidden = false; $("targets").hidden = true;
    $("picked").innerHTML = `<span>No ${K().noun} called “${esc(q)}”. ${K().miss}</span>
      <span class="end"><button type="button" class="btn quiet" data-clear>Clear</button></span>`;
  }, SEARCH);
  for (const id of ["ax", "ay"]) $(id).addEventListener("change", drawAgree);
  window.addEventListener("resize", () => { if (state.rows.length) drawVolcano(); });
  showKind();
  await showContrast();
  drawAgree();
}

// Switch between genes and TFs. A selection carries over when the other view has the same name
// (ESR1 the gene and ESR1 the TF); `then` selects a given key in the new view instead.
async function setKind(kind, then) {
  if (kind === state.kind && then == null) return;
  const name = state.pick != null ? K().nameOf(state.pick) : null;
  state.kind = kind;
  state.sort = { key: "rank", dir: 1 }; state.page = 0; state.allTargets = false;
  $("filter").value = "";
  showKind();
  state.rows = await load(state.c);
  state.pick = then ?? (name == null ? null : K().keyOf(name));
  setHash();
  redraw();
}

function showKind() {
  document.querySelectorAll("[data-kind]").forEach((b) => b.setAttribute("aria-checked", b.dataset.kind === state.kind));
  document.querySelectorAll(".noun").forEach((el) => (el.textContent = K().noun));
  $("lfc-part").hidden = !K().lfcCut;
  $("find").placeholder = `e.g. ${K().example}`;
}

async function showContrast() {
  state.c = $("contrast").value;
  state.rows = await load(state.c);
  state.page = 0;
  setHash();
  redraw(false);
}

// After the comparison or the view changes: pick() draws everything when there is a selection
function redraw(agree = true) {
  summary();
  if (state.pick != null) { pick(state.pick, true); return; }
  hidePick(); drawVolcano(); drawTable();
  if (agree) drawAgree();
}

function setHash() {
  const sel = state.pick != null ? `&${state.kind}=${encodeURIComponent(K().nameOf(state.pick))}` : isTF() ? "&view=tf" : "";
  history.replaceState(null, "", `#c=${state.c}${sel}`);
}

function refresh() {
  summary();
  drawVolcano();
  drawTable();
}

function summary() {
  $("thr-reset").hidden = alpha() === DEFAULTS.alpha && lfcMin() === DEFAULTS.lfcmin;
  const rows = state.rows, [num] = sides(state.c);
  const up = rows.filter((r) => isSig(r) && r.x > 0).length, down = rows.filter((r) => isSig(r) && r.x < 0).length;
  $("summary").innerHTML = `<b>${esc(fullLabel(state.c))}</b> <span>${rows.nTested.toLocaleString()} ${K().noun}s ${K().verb}</span>
    <span class="up">↑ ${up.toLocaleString()} ${esc(K().up(num))}</span> <span class="down">↓ ${down.toLocaleString()} ${esc(K().down(num))}</span>`;
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
  const rows = state.rows.filter((r) => r.padj != null && r.x != null);
  const [num, den] = sides(state.c);
  // Cap the y axis at 2.5x the 95th percentile of the significant points (padj < 0.05) so that a
  // few genes at p ~ 1e-300 do not squash the rest; points above the cap are drawn as triangles on it.
  const sv = rows.filter((r) => r.padj > 0 && r.padj < 0.05).map((r) => nlog(r.padj)).sort((a, b) => a - b);
  const vmax = Math.max(K().ymin, ...rows.filter((r) => r.padj > 0).map((r) => nlog(r.padj)));
  const cap = sv.length ? Math.min(vmax, Math.max(10, 2.5 * sv[Math.floor(sv.length * 0.95)])) : vmax;
  const capped = (r) => r.padj === 0 || nlog(r.padj) > cap;
  const yOf = (r) => (capped(r) ? cap : nlog(r.padj));
  const xmax = Math.max(1, ...rows.map((r) => Math.abs(r.x))) * 1.06;
  const xr = [-xmax, xmax], yr = [0, cap * 1.08];
  const margin = { l: 58, r: 14, t: 36, b: 52 };
  const colorOf = (r) => (isSig(r) ? (r.x > 0 ? UP : DOWN) : NS);
  const size = K().size;
  const trace = (rs, extra) => ({
    type: "scattergl", mode: "markers", x: rs.map((r) => r.x), y: rs.map(yOf), text: rs.map((r) => r.name),
    customdata: rs.map((r) => [r.key, fmtP(r.padj)]),
    hovertemplate: `<b>%{text}</b><br>${K().effect} %{x:.2f} · p<sub>adj</sub> %{customdata[1]}<extra></extra>`, ...extra,
  });
  const plain = rows.filter((r) => !capped(r)), above = rows.filter(capped);
  const traces = [
    trace(plain.filter((r) => !isSig(r)), { name: "not significant", marker: { color: NS, size: size - 1.5 } }),
    trace(plain.filter((r) => isSig(r) && r.x < 0), { name: K().down(num), marker: { color: DOWN, size } }),
    trace(plain.filter((r) => isSig(r) && r.x > 0), { name: K().up(num), marker: { color: UP, size } }),
  ];
  if (above.length) traces.push(trace(above, { showlegend: false, marker: { color: above.map(colorOf), size: 7, symbol: "triangle-up" } }));

  // labels: the extremes first (capped points by effect, then the largest effects on each side),
  // then the most significant; the picked point gets its own label
  const sig = rows.filter((r) => isSig(r) && r.key !== state.pick && !/^ENSG\d+/.test(r.name));
  const seen = new Set(), cands = [];
  for (const r of [
    ...sig.filter(capped).sort((a, b) => Math.abs(b.x) - Math.abs(a.x)).slice(0, 4),
    ...sig.filter((r) => r.x > 0).sort((a, b) => b.x - a.x).slice(0, 3),
    ...sig.filter((r) => r.x < 0).sort((a, b) => a.x - b.x).slice(0, 3),
    ...sig.filter((r) => !capped(r)).sort((a, b) => a.padj - b.padj).slice(0, 24),
  ]) if (!seen.has(r.key)) { seen.add(r.key); cands.push({ name: r.name, x: r.x, y: yOf(r) }); }
  const labels = el.clientWidth >= 600 ? placeLabels(cands, xr, yr, el.clientWidth, el.clientHeight, margin) : [];
  const annotations = labels.map((l) => ({ x: l.x, y: l.y, ax: l.ax, ay: l.ay, xanchor: l.xanchor, text: esc(l.name),
    showarrow: true, arrowhead: 0, arrowwidth: 1, arrowcolor: "#b4b8c4", font: { size: 11.5, color: INK } }));

  const p = state.pick != null && state.rows.byKey.get(state.pick);
  if (p && p.padj != null) {
    traces.push(ring(p.x, yOf(p)));
    annotations.push(ringLabel(p.name, p.x, yOf(p)));
  }
  if (above.length) annotations.push({ xref: "paper", yref: "paper", x: 0, y: 1, xanchor: "left", yanchor: "top", showarrow: false,
    text: "▲ beyond the axis", font: { size: 11.5, color: MUTED } });
  const t = K().lfcCut ? lfcMin() : 0, a = nlog(alpha());
  const guide = { color: "#b4b8c4", width: 1, dash: "dot" };
  const shapes = [
    { type: "line", xref: "paper", x0: 0, x1: 1, y0: a, y1: a, line: guide },
    ...(t > 0 ? [-t, t].map((x) => ({ type: "line", yref: "paper", y0: 0, y1: 1, x0: x, x1: x, line: guide })) : []),
  ];
  Plotly.react(el, traces, {
    ...BASE_LAYOUT, margin, shapes, annotations, hovermode: "closest",
    legend: { orientation: "h", x: 0, y: 1.02, yanchor: "bottom", traceorder: "reversed", itemclick: false, itemdoubleclick: false },
    xaxis: { ...AXIS, fixedrange: false, range: xr, title: { text: K().axis(esc(num), esc(den)), standoff: 8 }, zeroline: true, zerolinecolor: RULE },
    yaxis: { ...AXIS, fixedrange: false, range: yr, title: { text: "−log₁₀ p<sub>adj</sub>", standoff: 8 } },
  }, PLOT_CONFIG);
  el._home = { x: [...xr], y: [...yr] };   // copies: Plotly edits the range arrays it was given when the user zooms
  zoomReset(el, $("volcano-reset"), drawVolcano);
  pointCursor(el);
  el.removeAllListeners?.("plotly_click");
  el.on("plotly_click", (ev) => { const k = ev.points[0]?.customdata?.[0]; if (k != null) pick(k, false); });
}

const ring = (x, y) => ({ type: "scatter", mode: "markers", x: [x], y: [y], hoverinfo: "skip", showlegend: false,
  marker: { size: 17, color: "rgba(0,0,0,0)", line: { color: INK, width: 2 } } });
const ringLabel = (name, x, y) => ({ x, y, text: `<b>${esc(name)}</b>`, showarrow: true, arrowhead: 0, arrowcolor: INK,
  ax: 0, ay: -34, bgcolor: "#fff", bordercolor: INK, borderpad: 3, font: { size: 13, color: INK } });

// the tested rows, filtered and sorted as the table shows them
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
  const rows = filtered(), cols = K().cols;
  const pages = Math.max(1, Math.ceil(rows.length / PAGE));
  state.page = Math.min(Math.max(0, state.page), pages - 1);
  const head = `<thead><tr>${cols.map((c) => `<th scope="col" data-key="${c.key}" ${c.title ? `title="${c.title}"` : ""}
    class="${state.sort.key === c.key ? (state.sort.dir > 0 ? "asc" : "desc") : ""}">${c.label}</th>`).join("")}</tr></thead>`;
  const body = rows.slice(state.page * PAGE, (state.page + 1) * PAGE).map((r) => {
    const sig = isSig(r), dir = sig ? (r.x > 0 ? "up" : "down") : "";
    return `<tr data-k="${esc(String(r.key))}" tabindex="0" class="${r.key === state.pick ? "picked-row" : ""}">
      <td>${r.rank.toLocaleString()}</td><td>${esc(r.name)}</td>
      <td>${K().count(r).toLocaleString()}</td>
      <td class="${dir ? `sig ${dir}` : ""}">${fmt(r.x)}</td>
      <td class="${sig ? "sig" : ""}">${fmtP(r.padj)}</td></tr>`;
  }).join("");
  $("detable").innerHTML = head + `<tbody>${body}</tbody>`;
  const nSig = state.rows.filter(isSig).length;
  $("count").textContent = $("sigonly").checked
    ? `${rows.length.toLocaleString()} significant ${K().noun}s`
    : `${rows.length.toLocaleString()} ${K().noun}s · ${nSig.toLocaleString()} significant`;
  $("page").textContent = `Page ${state.page + 1} of ${pages}`;
  $("prev").disabled = state.page === 0;
  $("next").disabled = state.page >= pages - 1;
}

// select a gene or TF: say where it ranks, ring it on the plots, mark its row, and for a TF list its targets
function pick(key, jumpTable) {
  state.pick = key;
  setHash();
  const name = K().nameOf(key), r = state.rows.byKey.get(key), [num] = sides(state.c);
  const gene = geneByName(name);
  const link = gene != null ? `<a href="${genesLink([name])}">Show expression</a>` : "";
  const clear = `<button type="button" class="btn quiet" data-clear>Clear</button>`;
  $("picked").hidden = false;
  if (!r || r.p == null) {
    $("picked").innerHTML = `<b>${esc(name)}</b> <span>${K().untested}</span>
      <span class="end">${link} ${clear}</span>`;
    $("targets").hidden = true;
  } else {
    const pct = Math.max(1, Math.ceil((r.rank / state.rows.nTested) * 100));
    const call = r.padj == null ? `<span class="dim">no p<sub>adj</sub></span>`
      : isSig(r) ? `<span class="${r.x > 0 ? "up" : "down"}"><b>${esc(r.x > 0 ? K().up(num) : K().down(num))}</b></span>`
      : `<span class="dim">not significant</span>`;
    $("picked").innerHTML = `<b>${esc(name)}</b>
      <span>ranked ${r.rank.toLocaleString()} of ${state.rows.nTested.toLocaleString()} by p-value, top ${pct}%</span>
      <span>${K().effect} ${fmt(r.x)}</span> <span>p<sub>adj</sub> ${fmtP(r.padj)}</span> ${call}
      <span class="end">${link} ${clear}</span>`;
    if (jumpTable) {
      $("filter").value = ""; $("sigonly").checked = false;
      state.sort = { key: "rank", dir: 1 };
      state.page = Math.floor((r.rank - 1) / PAGE);   // the table is now every row, by rank
    }
    if (isTF()) drawTargets(); else $("targets").hidden = true;
  }
  drawVolcano();
  drawTable();
  drawAgree();
}

// The picked TF's targets in this comparison, those pushing its score hardest first: a target's
// push is its test statistic, sign-flipped for a repressed target, in the direction of the activity.
async function drawTargets() {
  const tf = state.pick, c = state.c;
  network ??= fetch("data/tf/network.json").then((r) => r.json());
  const [net, genes] = await Promise.all([network, load(c, "gene")]);
  if (state.pick !== tf || state.c !== c || !isTF()) return;
  // a name shared by several genes counts once, the most expressed, as in the TF scores
  if (!genes.byName) {
    genes.byName = new Map();
    for (const g of genes) if (g.z != null && !(genes.byName.get(g.name)?.base > g.base)) genes.byName.set(g.name, g);
  }
  const [act, rep] = net[tf] ?? [[], []];
  const sign = Math.sign(state.rows.byKey.get(tf).x) || 1;
  const all = [...act.map((t) => [t, 1]), ...rep.map((t) => [t, -1])];
  const list = all.map(([t, w]) => ({ g: genes.byName.get(t), w })).filter(({ g }) => g)
    .map(({ g, w }) => ({ g, w, push: w * g.z * sign })).sort((a, b) => b.push - a.push);
  const shown = state.allTargets ? list : list.slice(0, TOP_TARGETS);
  const rows = shown.map(({ g, w }) => {
    const sig = isSig(g), dir = sig ? (g.x > 0 ? "up" : "down") : "";
    return `<tr><td><button type="button" class="linklike" data-gene="${g.key}" title="Select ${esc(g.name)} in the gene view">${esc(g.name)}</button></td>
      <td>${w > 0 ? "activates" : "represses"}</td><td class="${dir ? `sig ${dir}` : ""}">${fmt(g.x)}</td><td class="${sig ? "sig" : ""}">${fmtP(g.padj)}</td></tr>`;
  }).join("");
  const more = list.length > TOP_TARGETS
    ? `<button type="button" class="btn quiet" data-all>${state.allTargets ? `Show the top ${TOP_TARGETS}` : `Show all ${list.length}`}</button>` : "";
  $("targets").innerHTML = `<p class="small dim">${list.length} of ${all.length} targets of ${esc(tf)} measured in this
      comparison; those that push its score most come first. Click one to see it in the gene view.</p>
    <div class="scroll"><table><thead><tr><th>Target</th><th>Effect</th><th>log₂FC</th><th>p<sub>adj</sub></th></tr></thead>
    <tbody>${rows}</tbody></table></div>${more}`;
  $("targets").hidden = false;
}

function clearPick() {
  state.pick = null;
  hidePick();
  setHash();
  drawVolcano();
  drawTable();
  drawAgree();
}

function download() {
  const { head, row, suffix } = K().csv;
  const lines = [head, ...filtered().map((r) => row(r).map((v) => v ?? "").join(","))];
  saveText(lines.join("\n") + "\n", `${state.c}${suffix}${$("sigonly").checked ? "_significant" : ""}.csv`);
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
  const [cx, cy, kind] = [$("ax").value, $("ay").value, state.kind];
  const [X, Y] = await Promise.all([load(cx), load(cy)]);
  if (cx !== $("ax").value || cy !== $("ay").value || kind !== state.kind) return;
  const pairs = [];
  for (const rx of X) {
    const ry = Y.byKey.get(rx.key);
    if (ry && rx.x != null && ry.x != null) pairs.push([rx, ry]);
  }
  const cats = [
    { name: "neither significant", color: AGREE.neither, f: (a, b) => !isSig(a) && !isSig(b) },
    { name: "significant in one", color: AGREE.one, f: (a, b) => isSig(a) !== isSig(b) },
    { name: "both, opposite directions", color: AGREE.opposite, f: (a, b) => isSig(a) && isSig(b) && Math.sign(a.x) !== Math.sign(b.x) },
    { name: "both, same direction", color: AGREE.same, f: (a, b) => isSig(a) && isSig(b) && Math.sign(a.x) === Math.sign(b.x) },
  ];
  const sets = cats.map((c) => pairs.filter(([a, b]) => c.f(a, b)));
  const xname = dataset(dsOf(cx)).name, yname = dataset(dsOf(cy)).name;
  const traces = cats.map((c, k) => {
    const ps = sets[k];
    return { type: "scattergl", mode: "markers", name: `${c.name} (${ps.length.toLocaleString()})`,
      x: ps.map(([a]) => a.x), y: ps.map(([, b]) => b.x), text: ps.map(([a]) => a.name), customdata: ps.map(([a]) => a.key),
      marker: { color: c.color, size: k === 0 ? K().size - 1.5 : K().size },
      hovertemplate: `<b>%{text}</b><br>${esc(xname)} %{x:.2f} · ${esc(yname)} %{y:.2f}<extra></extra>` };
  });
  const rho = pairs.length > 2 ? spearman(pairs.map(([a]) => a.x), pairs.map(([, b]) => b.x)) : NaN;
  $("agree-summary").innerHTML = cx === cy
    ? `The same comparison on both axes. Pick two different ones to compare them.`
    : `${pairs.length.toLocaleString()} ${K().noun}s ${K().verb} in both · rank correlation ${rho.toFixed(2)} · significant in both:
       <b>${sets[3].length.toLocaleString()}</b> in the same direction, <b>${sets[2].length.toLocaleString()}</b> in opposite directions`;
  const lim = Math.max(1, ...pairs.flatMap(([a, b]) => [Math.abs(a.x), Math.abs(b.x)])) * 1.05;
  const annotations = [];
  const px = X.byKey.get(state.pick), py = Y.byKey.get(state.pick);
  if (px?.x != null && py?.x != null) { traces.push(ring(px.x, py.x)); annotations.push(ringLabel(px.name, px.x, py.x)); }
  const el = $("agree");
  const axis = (title) => ({ ...AXIS, fixedrange: false, range: [-lim, lim], zeroline: true, zerolinecolor: RULE, title: { text: title, standoff: 8 } });
  Plotly.react(el, traces, {
    ...BASE_LAYOUT, margin: { l: 60, r: 14, t: 56, b: 56 }, hovermode: "closest", annotations,
    legend: { orientation: "h", x: 0, y: 1.02, yanchor: "bottom", traceorder: "reversed", itemclick: false, itemdoubleclick: false },
    shapes: [{ type: "line", x0: -lim, y0: -lim, x1: lim, y1: lim, line: { color: "#b4b8c4", width: 1, dash: "dot" } }],
    xaxis: axis(`${esc(fullLabel(cx))}, ${K().effect}`),
    yaxis: { ...axis(`${esc(fullLabel(cy))}, ${K().effect}`), scaleanchor: "x", scaleratio: 1 },
  }, PLOT_CONFIG);
  el._home = { x: [-lim, lim], y: [-lim, lim] };
  zoomReset(el, $("agree-reset"), drawAgree);
  pointCursor(el);
  el.removeAllListeners?.("plotly_click");
  el.on("plotly_click", (ev) => { const k = ev.points[0]?.customdata; if (k != null) pick(k, true); });
}

init().catch((e) => { $("summary").textContent = "The data did not load: " + e.message; });
