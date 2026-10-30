// Overview page: a gene search that opens the Genes page, and the sample figure, drawn from
// the real sample table (data/meta.json) so that every dot is one sample.
const CITE = {
  JCI_tissue: "Zuberi et al., JCI Insight 2023",
  NC_2D: "Buyukcelebi et al., Nat Commun 2023",
  NC_3D_spheroid: "Buyukcelebi et al., Nat Commun 2023",
};
const ABOUT = {
  JCI_tissue: "Leiomyoma (LM) and matched myometrium (MM) from 10 patients. Six of the tumours carry the MED12 G44D mutation, four are MED12 wild-type.",
  NC_2D: "hTERT-immortalized myometrial smooth muscle cells in 2D culture: unedited (WT), MED12 knocked out (KO), and two clones carrying a MED12 codon-44 mutation (Mut1, Mut2). Three replicates each.",
  NC_3D_spheroid: "The same cells grown as 3D spheroids: WT with two replicates, Mut1 and Mut2 with three.",
};
const FULL = { MM: "MM, myometrium", LM_wtMED12: "LM, MED12 wild-type", LM_mutMED12: "LM, MED12 G44D", WT: "WT", KO: "MED12 knock-out", Mut1: "Mut1", Mut2: "Mut2" };
const R = 6, DOT = 20;

const dot = (x, y, color, title) =>
  `<circle cx="${x}" cy="${y}" r="${R}" fill="${color}" stroke="#fff" stroke-width="1.5"><title>${esc(title)}</title></circle>`;
const label = (x, y, text, anchor = "middle") =>
  `<text x="${x}" y="${y}" text-anchor="${anchor}" font-size="12" fill="${MUTED}">${esc(text)}</text>`;

// tissue: one column per patient, MM above LM, joined by a line
function tissueSvg(d) {
  const patients = [...new Set(d.samples.map((s) => s.patient))];
  const x0 = 36, w = x0 + patients.length * DOT, parts = [label(0, 30, "MM", "start"), label(0, 78, "LM", "start")];
  patients.forEach((p, i) => {
    const x = x0 + i * DOT + R, mm = d.samples.find((s) => s.patient === p && s.group === "MM"), lm = d.samples.find((s) => s.patient === p && s.group !== "MM");
    parts.push(`<line x1="${x}" y1="${26 + R}" x2="${x}" y2="${74 - R}" stroke="${RULE}" stroke-width="1.5"/>`);
    parts.push(dot(x, 26, COLORS.MM, `Patient ${p}, myometrium`));
    parts.push(dot(x, 74, COLORS[lm.group], `Patient ${p}, leiomyoma (${groupLabel(lm.group)})`));
  });
  return `<svg width="${w}" height="92" viewBox="0 0 ${w} 92" role="img" aria-label="${patients.length} patients, each with a myometrium and a leiomyoma sample">${parts.join("")}</svg>`;
}

// cells: one column per group, replicates stacked
function cellsSvg(d) {
  const col = 56, parts = [];
  d.groups.forEach((g, k) => {
    const ks = d.samples.filter((s) => s.group === g);
    const x = 20 + k * col;
    ks.forEach((s, j) => parts.push(dot(x, 14 + j * DOT, COLORS[g], sampleLabel(s.sample))));
    parts.push(label(x, 88, g));
  });
  const w = 20 + d.groups.length * col - 16;
  return `<svg width="${w}" height="92" viewBox="0 0 ${w} 92" role="img" aria-label="${d.samples.length} samples in ${d.groups.length} groups">${parts.join("")}</svg>`;
}

function sampleBlocks(meta) {
  return meta.datasets.map((d) => {
    const counts = d.groups.map((g) => `<span><span class="swatch" style="background:${COLORS[g]}"></span>${esc(FULL[g] || g)} (${d.samples.filter((s) => s.group === g).length})</span>`);
    return `<section class="sample-block" aria-label="${esc(SHORT_DS[d.id])}">
      <h3>${esc(SHORT_DS[d.id])} <span class="dim small">· ${d.samples.length} samples</span></h3>
      ${d.id === "JCI_tissue" ? tissueSvg(d) : cellsSvg(d)}
      <p class="legend">${counts.join("")}</p>
      <p class="about">${esc(ABOUT[d.id] || "")} <span class="dim">${esc(CITE[d.id] || "")}.</span></p>
    </section>`;
  }).join("");
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
  const n = Data.meta.datasets.reduce((a, d) => a + d.samples.length, 0), c = Data.meta.datasets.reduce((a, d) => a + d.contrasts.length, 0);
  $("glance").textContent = `${n} samples · ${Data.meta.datasets.length} datasets · ${c} comparisons · 62,298 GENCODE v44 genes`;
}

init().catch((e) => { $("samples").innerHTML = `<p class="status">The sample table did not load: ${esc(e.message)}</p>`; });
