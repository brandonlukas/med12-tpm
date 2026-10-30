// Overview page: a gene search that opens the Genes page, a "try a gene" strip with the real
// data for a few example genes, and the sample figure drawn from the sample table
// (data/meta.json), so that every dot is one sample.
const EXAMPLES = ["MED12", "HMGA2", "PGR", "ESR1", "CDK8"];
const FULL = { MM: "MM, myometrium", LM_wtMED12: "LM, MED12 wild-type", LM_mutMED12: "LM, MED12 G44D", WT: "WT", KO: "MED12 knock-out", Mut1: "Mut1", Mut2: "Mut2" };
const TICK = { MM: "MM", LM_wtMED12: "LM WT", LM_mutMED12: "LM mut" };
const R = 6, DOT = 20;
const reduced = () => window.matchMedia?.("(prefers-reduced-motion: reduce)").matches;
// deterministic jitter in [-0.5, 0.5], the same for a sample whichever gene is shown
const jitter = (k) => ((Math.abs(Math.sin(k * 12.9898) * 43758.5453) % 1) - 0.5);

// ---- the sample figure
const dot = (x, y, color, title, group, i) =>
  `<circle cx="${x}" cy="${y}" r="${R}" fill="${color}" stroke="#fff" stroke-width="1.5" data-group="${group}" style="animation-delay:${i * 22}ms"><title>${esc(title)}</title></circle>`;
const label = (x, y, text, anchor = "middle") =>
  `<text x="${x}" y="${y}" text-anchor="${anchor}" font-size="12" fill="${MUTED}">${esc(text)}</text>`;

// tissue: one column per patient, MM above LM, joined by a line
function tissueSvg(d) {
  const patients = [...new Set(d.samples.map((s) => s.patient))];
  const x0 = 36, w = x0 + patients.length * DOT, parts = [label(0, 30, "MM", "start"), label(0, 78, "LM", "start")];
  patients.forEach((p, i) => {
    const x = x0 + i * DOT + R, lm = d.samples.find((s) => s.patient === p && s.group !== "MM");
    parts.push(`<line x1="${x}" y1="${26 + R}" x2="${x}" y2="${74 - R}" stroke="${RULE}" stroke-width="1.5"/>`);
    parts.push(dot(x, 26, COLORS.MM, `Patient ${p}, myometrium`, "MM", i));
    parts.push(dot(x, 74, COLORS[lm.group], `Patient ${p}, leiomyoma (${groupLabel(lm.group)})`, lm.group, i + patients.length));
  });
  return `<svg width="${w}" height="92" viewBox="0 0 ${w} 92" role="img" aria-label="${patients.length} patients, each with a myometrium and a leiomyoma sample">${parts.join("")}</svg>`;
}

// cells: one column per group, replicates stacked
function cellsSvg(d) {
  const col = 56, parts = [];
  let i = 0;
  d.groups.forEach((g, k) => {
    const x = 20 + k * col;
    d.samples.filter((s) => s.group === g).forEach((s, j) => parts.push(dot(x, 14 + j * DOT, COLORS[g], sampleLabel(s.sample), g, i++)));
    parts.push(label(x, 88, g));
  });
  const w = 20 + d.groups.length * col - 16;
  return `<svg width="${w}" height="92" viewBox="0 0 ${w} 92" role="img" aria-label="${d.samples.length} samples in ${d.groups.length} groups">${parts.join("")}</svg>`;
}

function sampleBlocks(meta) {
  return meta.datasets.map((d) => {
    const counts = d.groups.map((g) => `<span data-group="${g}"><span class="swatch" style="background:${COLORS[g]}"></span>${esc(FULL[g] || g)} (${d.samples.filter((s) => s.group === g).length})</span>`);
    return `<section class="sample-block" aria-label="${esc(SHORT_DS[d.id])}">
      <h3>${esc(SHORT_DS[d.id])} <span class="dim small">· ${d.samples.length} samples</span></h3>
      ${d.id === "JCI_tissue" ? tissueSvg(d) : cellsSvg(d)}
      <p class="legend">${counts.join("")}</p>
      <p class="about">${datasetNote(d.id)}</p>
    </section>`;
  }).join("");
}

// hovering a legend entry fades the other groups' dots in that block
function wireLegends() {
  document.querySelectorAll(".sample-block").forEach((block) => {
    const dots = block.querySelectorAll("circle");
    block.querySelectorAll(".legend [data-group]").forEach((el) => {
      el.addEventListener("mouseenter", () => dots.forEach((c) => c.classList.toggle("dim", c.dataset.group !== el.dataset.group)));
      el.addEventListener("mouseleave", () => dots.forEach((c) => c.classList.remove("dim")));
    });
  });
}

// ---- try a gene: small TPM plots, one per dataset, that glide between example genes
const M = { w: 300, h: 190, l: 36, r: 8, t: 12, b: 34 };
const ly = (v) => Math.log2(v + 1);
const FLOOR = ly(10);
let shown = null;

function miniShell(d) {
  const { w, h, l, r, t, b } = M, pw = w - l - r, G = d.groups.length;
  const gx = (k) => l + ((k + 0.5) * pw) / G;
  const xi = Object.fromEntries(d.groups.map((g, k) => [g, gx(k)]));
  const parts = [`<g class="yticks"></g>`, `<line x1="${l}" y1="${t}" x2="${l}" y2="${h - b}" stroke="${RULE}"/>`];
  d.groups.forEach((g, k) => parts.push(`<text x="${gx(k)}" y="${h - b + 17}" text-anchor="middle" font-size="11" fill="${MUTED}">${esc(TICK[g] || g)}</text>`));
  const pairs = {};
  d.samples.forEach((s, k) => { if (s.patient) (pairs[s.patient] ||= []).push(k); });
  const xs = d.samples.map((s, k) => xi[s.group] + jitter(k) * Math.min(16, (pw / G) * 0.3));
  Object.values(pairs).forEach((ks) => { if (ks.length === 2) parts.push(`<line class="pair" data-a="${ks[0]}" data-b="${ks[1]}" x1="${xs[ks[0]].toFixed(1)}" x2="${xs[ks[1]].toFixed(1)}" y1="${h - b}" y2="${h - b}" stroke="${RULE}" stroke-width="1"/>`); });
  d.samples.forEach((s, k) => parts.push(`<circle data-k="${k}" cx="${xs[k].toFixed(1)}" cy="${h - b}" r="5" fill="${COLORS[s.group]}" stroke="#fff" stroke-width="1.2"><title></title></circle>`));
  return `<section class="mini-block"><h3>${esc(SHORT_DS[d.id])}</h3>
    <svg class="mini" viewBox="0 0 ${w} ${h}" role="img" aria-label="${esc(SHORT_DS[d.id])}: TPM of the chosen gene in every sample" id="mini-${d.id}">${parts.join("")}</svg>
    <p class="cap" id="cap-${d.id}"></p></section>`;
}

// move circles and pair lines to new y positions over ~600 ms (at once if motion is reduced)
function tween(svg, targets) {
  const els = [...svg.querySelectorAll("circle")], from = els.map((c) => +c.getAttribute("cy"));
  const lines = [...svg.querySelectorAll("line.pair")];
  const place = (f) => {
    els.forEach((c, k) => c.setAttribute("cy", (from[k] + (targets[k] - from[k]) * f).toFixed(1)));
    lines.forEach((ln) => { ln.setAttribute("y1", els[+ln.dataset.a].getAttribute("cy")); ln.setAttribute("y2", els[+ln.dataset.b].getAttribute("cy")); });
  };
  if (reduced() || from.every((v, k) => Math.abs(v - targets[k]) < 0.5)) { place(1); return; }
  const t0 = performance.now(), dur = 600, ease = (x) => 1 - Math.pow(1 - x, 3);
  const step = (now) => { const f = Math.min(1, (now - t0) / dur); place(ease(f)); if (f < 1) requestAnimationFrame(step); };
  requestAnimationFrame(step);
}

function miniUpdate(d, g, name) {
  const { h, l, t, b } = M, ph = h - t - b, svg = $(`mini-${d.id}`), tpm = g.tpm[d.id];
  // axis top: the next decade above the highest value (at least 10 TPM), so the dots sit under a labelled tick
  const vmax = Math.max(0, ...tpm.filter((v) => v != null));
  const top = Math.max(FLOOR, ly(10 ** Math.ceil(Math.log10(Math.max(vmax, 1))))) * 1.04;
  const y = (v) => h - b - (ly(v ?? 0) / top) * ph;
  svg.querySelectorAll("circle").forEach((c) => { const k = +c.dataset.k; c.querySelector("title").textContent = `${sampleLabel(d.samples[k].sample)} · TPM ${fmtTPM(tpm[k])}`; });
  const vals = [0, 1, 10, 100, 1000, 10000, 100000].filter((v) => ly(v) <= top);
  svg.querySelector(".yticks").innerHTML = vals.map((v) => `<line x1="${l - 3}" x2="${l}" y1="${y(v).toFixed(1)}" y2="${y(v).toFixed(1)}" stroke="#c3c7d3"/>
    <text x="${l - 6}" y="${(y(v) + 3.5).toFixed(1)}" text-anchor="end" font-size="10.5" fill="${MUTED}">${v.toLocaleString()}</text>`).join("")
    + `<text x="9" y="${(t + ph / 2).toFixed(1)}" text-anchor="middle" font-size="10.5" fill="${MUTED}" transform="rotate(-90 9 ${(t + ph / 2).toFixed(1)})">TPM (log)</text>`;
  tween(svg, tpm.map(y));
  // the caption: the comparison with the smallest padj
  const de = g.de[d.id], tested = d.contrasts.filter((c) => de[c]);
  let cap;
  if (!tested.length) cap = `<span class="dim">${esc(name)} was not tested here: too few counts.</span>`;
  else {
    const c = tested.reduce((a, b) => (de[a][2] <= de[b][2] ? a : b)), [, lfc, padj] = de[c], [num] = contrastSides(c);
    cap = padj < 0.05
      ? `${esc(contrastLabel(c))}: <b>${lfc > 0 ? "higher" : "lower"} in ${esc(num)}</b>, log₂FC ${fmt(lfc)}, p<sub>adj</sub> ${fmtP(padj)}`
      : `<span class="dim">No significant difference in any comparison (smallest p<sub>adj</sub> ${fmtP(padj)}).</span>`;
  }
  $(`cap-${d.id}`).innerHTML = cap;
}

async function showExample(name) {
  const i = Data.lookup(name)[0];
  if (i == null) return;
  shown = name;
  document.querySelectorAll("#examples button").forEach((b) => b.setAttribute("aria-pressed", b.dataset.gene === name));
  const g = await Data.gene(i);
  if (shown !== name) return;
  for (const d of Data.meta.datasets) miniUpdate(d, g, name);
  $("try-link").innerHTML = `<a href="${genesLink([name])}">Open ${esc(name)} in Genes</a> for the full view, with every comparison.`;
}

async function init() {
  const form = $("lookup"), q = $("q"), msg = $("lookup-msg");
  const go = (i) => { location.href = genesLink([Data.genes[i][1]]); };
  const miss = (text) => { msg.hidden = false; msg.textContent = `No gene called “${text}”. Try its official symbol, for example HMGA2.`; };
  form.addEventListener("submit", (e) => {
    e.preventDefault();
    if (!Data.genes) { location.href = `genes.html#genes=${encodeURIComponent(q.value.trim())}`; return; }
    const hit = Data.lookup(q.value)[0] ?? Data.suggest(q.value, 1)[0];
    if (hit != null) go(hit); else if (q.value.trim()) miss(q.value.trim());
  });
  await Data.load();
  attachSearch(q, $("suggest"), go, miss);

  $("samples").innerHTML = sampleBlocks(Data.meta);
  wireLegends();
  const n = Data.meta.datasets.reduce((a, d) => a + d.samples.length, 0), c = Data.meta.datasets.reduce((a, d) => a + d.contrasts.length, 0);
  $("glance").textContent = `${n} samples · ${Data.meta.datasets.length} datasets · ${c} comparisons · 62,298 GENCODE v44 genes`;

  $("examples").insertAdjacentHTML("beforeend", EXAMPLES.map((g) => `<button type="button" data-gene="${g}" aria-pressed="false">${g}</button>`).join(""));
  $("examples").addEventListener("click", (e) => { const g = e.target.closest("[data-gene]")?.dataset.gene; if (g) showExample(g); });
  $("minis").innerHTML = Data.meta.datasets.map(miniShell).join("");
  await showExample(EXAMPLES[0]);
}

init().catch((e) => { $("samples").innerHTML = `<p class="status">The sample table did not load: ${esc(e.message)}</p>`; });
