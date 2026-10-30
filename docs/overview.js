// Overview page: a gene search that opens the Genes page.
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
}

init().catch(() => { /* the form still submits to the Genes page */ });
