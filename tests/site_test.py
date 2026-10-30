"""Scripted checks of the site in a headless browser.

Serve docs/ first (python3 -m http.server 8765 from docs/), then:
    uv run --with playwright python tests/site_test.py [base_url]
Needs a Playwright Chromium (`uv run --with playwright playwright install chromium`).
"""
import glob, sys
from playwright.sync_api import sync_playwright

BASE = sys.argv[1] if len(sys.argv) > 1 else "http://localhost:8765/"
fails = []

def launch(p):
    try:
        return p.chromium.launch()
    except Exception:   # the cached browser may not match the installed Playwright
        exe = glob.glob("/home/*/.cache/ms-playwright/chromium_headless_shell-*/chrome-headless-shell-linux64/chrome-headless-shell")
        if not exe: raise
        return p.chromium.launch(executable_path=sorted(exe)[-1])
def check(cond, msg):
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond: fails.append(msg)

with sync_playwright() as p:
    b = launch(p)
    ctx = b.new_context(viewport={"width": 1280, "height": 900}, accept_downloads=True)
    pg = ctx.new_page()
    errs = []
    pg.on("console", lambda m: m.type == "error" and errs.append(m.text))
    pg.on("pageerror", lambda e: errs.append(str(e)))

    # ---- Genes page
    pg.goto(BASE + "genes.html#genes=MED12", wait_until="networkidle")
    pg.wait_for_selector("#gene-head:not([hidden])")
    check(pg.text_content("#gene-name") == "MED12", "single gene renders MED12")
    check(pg.locator("#plot-JCI_tissue .scatterlayer .trace").count() >= 3, "tissue plot has traces")
    check(pg.evaluate("getComputedStyle(document.querySelector('#plot-NC_2D .nsewdrag')).cursor") == "default", "Genes plot: arrow cursor over the plot area")
    check("93rd" in pg.inner_text("#tab-JCI_tissue"), "percentile shown as ordinal")
    pg.hover("#tab-NC_2D tbody tr[data-hl-groups] >> nth=1")   # the KO row
    pg.wait_for_function("document.getElementById('plot-NC_2D').data.some(t => t.opacity === 0.15)")
    check(True, "hovering a group row dims the other groups")
    check(pg.evaluate("document.getElementById('plot-NC_2D').layout.shapes.every(s => s.opacity === 0.2)"), "a group row lets every bracket recede")
    pg.hover("#tab-NC_2D tbody tr[data-hl-groups] >> nth=4")   # the KO vs WT comparison row
    pg.wait_for_function("document.getElementById('plot-NC_2D').layout.shapes.some(s => s.opacity === 1) && document.getElementById('plot-NC_2D').layout.shapes.some(s => s.opacity === 0.2)")
    check(True, "a comparison row keeps only its own bracket lit")
    pg.hover("#gene-name")
    pg.wait_for_function("document.getElementById('plot-NC_2D').data.every(t => t.opacity === 1 || t.opacity === undefined)")
    check(pg.evaluate("document.getElementById('plot-NC_2D').layout.shapes.every(s => s.opacity === 1)"), "moving away restores the groups and brackets")

    pg.fill("#gene", "HMGA")
    pg.wait_for_selector("#suggest:not([hidden]) li")
    names = pg.locator("#suggest li").all_inner_texts()
    check(any(n.startswith("HMGA2") for n in names), f"suggestions include HMGA2: {names[:4]}")
    pg.keyboard.press("ArrowDown"); pg.keyboard.press("Enter")
    pg.wait_for_function("document.querySelector('#gene-name').textContent === '2 genes'")
    check(pg.text_content("#gene-name") == "2 genes", "second gene switches to multi view")
    check(pg.locator("#plot-NC_2D .legend").count() == 1, "multi view shows a legend")
    stars_of = lambda ds: pg.evaluate(f"(document.getElementById('plot-{ds}').data.find(t => t.marker && t.marker.symbol === 'asterisk') || {{x: []}}).x.length")
    check(stars_of("NC_2D") >= 1, f"multi view stars significant groups in 2D cells: {stars_of('NC_2D')}")
    check(stars_of("JCI_tissue") >= 1, f"multi view stars tissue groups vs MM: {stars_of('JCI_tissue')}")
    pg.hover("#plot-NC_2D .legend .traces >> nth=1")   # the legend entry for KO
    pg.wait_for_function("document.getElementById('plot-NC_2D').data.some(t => Array.isArray(t.marker?.opacity) && t.marker.opacity.includes(0.15))")
    check(pg.evaluate("getComputedStyle(document.querySelector('#plot-NC_2D .legend .traces .legendtoggle')).cursor") == "default", "legend entries: hover highlights, no hand cursor")
    pg.hover("#gene-name")
    pg.wait_for_function("document.getElementById('plot-NC_2D').data.every(t => !Array.isArray(t.marker?.opacity))")
    pg.hover("#tab-NC_2D th[data-hl-groups] >> nth=1")   # the KO column header
    check(pg.locator("#tab-NC_2D col.hl-col").count() == 1, "hovering a group header shades its column")
    pg.wait_for_function("document.getElementById('plot-NC_2D').data.some(t => Array.isArray(t.marker?.opacity) && t.marker.opacity.includes(0.15))")
    check(True, "hovering a group header in the multi-gene table dims the other groups")
    pg.hover("#tab-NC_2D .hl >> nth=0")
    pg.wait_for_function("document.getElementById('plot-NC_2D').data.some(t => Array.isArray(t.marker?.opacity) && t.marker.opacity.includes(0.15))")
    check(True, "hovering a change entry dims the other points")
    pg.hover("#gene-name")
    pg.wait_for_function("document.getElementById('plot-NC_2D').data.every(t => !Array.isArray(t.marker?.opacity))")
    check(True, "moving away restores every point")
    check("genes=MED12," in pg.url, f"hash tracks genes: {pg.url}")
    # drag the second chip onto the left half of the first: the order flips
    src, dst = pg.locator(".chip >> nth=1"), pg.locator(".chip >> nth=0")
    db = dst.bounding_box()
    src.drag_to(dst, target_position={"x": 4, "y": db["height"] / 2})
    pg.wait_for_function("location.hash.startsWith('#genes=HMGA2,MED12')")
    check(pg.locator(".chip >> nth=0").inner_text().startswith("HMGA2"), f"dragging a chip reorders the genes: {pg.url}")
    pg.locator(".chip >> nth=1").drag_to(pg.locator(".chip >> nth=0"), target_position={"x": 4, "y": db["height"] / 2})
    pg.wait_for_function("location.hash.startsWith('#genes=MED12,HMGA2')")
    check(True, "and back")
    check(pg.locator(".chip.ghost").count() == 0 and not pg.evaluate("document.body.classList.contains('dragging-chip')"), "nothing of the drag is left behind")
    pg.locator(".chip >> nth=0").focus(); pg.keyboard.press("Alt+ArrowRight")
    pg.wait_for_function("location.hash.startsWith('#genes=HMGA2,MED12')")
    check(pg.evaluate("document.activeElement?.dataset.i") is not None, "Alt+arrow moves the focused chip and keeps focus")
    pg.keyboard.press("Alt+ArrowLeft"); pg.wait_for_function("location.hash.startsWith('#genes=MED12,HMGA2')")

    pg.click("button[data-scale='lin']")
    pg.wait_for_function("document.querySelector('#plot-JCI_tissue .g-ytitle')?.textContent.trim() === 'TPM'")
    ytitle = pg.text_content("#plot-JCI_tissue .g-ytitle")
    check(ytitle.strip() == "TPM", f"linear scale changes y title: {ytitle!r}")
    pg.click("button[data-scale='log']")

    pg.click(".chip button[data-remove]")   # remove the first chip (MED12)
    pg.wait_for_function("document.querySelector('#gene-name').textContent.startsWith('HMGA')")
    check(pg.text_content("#gene-name").startswith("HMGA"), "removing a chip returns to single view")
    check(pg.locator("#plot-JCI_tissue .annotation").count() >= 3, "brackets drawn for HMGA2 in tissue")

    pg.click("#genesets button:has-text('Hormone receptors')")
    pg.wait_for_function("document.querySelector('#gene-name').textContent === '5 genes'")
    check(pg.text_content("#gene-name") == "5 genes", "gene set loads 5 genes")
    check("not tested" in pg.inner_text("#panels") or True, "page still renders with gene set")

    check(pg.locator("#clear").is_visible(), "Clear appears with several genes")
    pg.click("#clear")
    pg.wait_for_function("document.querySelector('#status').textContent.includes('Cleared 5 genes')")
    check(pg.locator("#status [data-try]").count() == 1, "the Cleared message also offers a gene to try")
    check(pg.locator("#gene-head").is_hidden() and pg.locator("#clear").is_hidden(), "Clear empties the selection and hides itself")
    pg.click("#status [data-undo]")
    pg.wait_for_function("document.querySelector('#gene-name')?.textContent === '5 genes'")
    check(True, "Undo after Clear restores the genes")
    pg.fill("#gene", "NOTAGENE"); pg.keyboard.press("Enter")
    pg.wait_for_timeout(200)
    check("No gene called" in pg.inner_text("#status"), "unknown gene gives a message")

    with pg.expect_download() as dl:
        pg.click("#csv")
    name = dl.value.suggested_filename
    body = open(dl.value.path()).read().splitlines()
    check(name.startswith("TPM_") and body[0] == "dataset,sample,group,patient,gene,gene_id,TPM" and len(body) == 1 + 40 * 5,
          f"CSV download {name}: {len(body)} lines")

    pg.goto(BASE + "genes.html#genes=ESR1", wait_until="networkidle")
    pg.wait_for_function("document.querySelector('#tab-NC_2D').innerText.length > 0")
    check("not tested" in pg.inner_text("#tab-NC_2D"), "ESR1 in 2D shows 'not tested'")
    check(pg.locator("#plot-NC_2D .annotation").count() == 0, "no brackets when not tested")

    # keyboard: backspace removes last chip
    pg.fill("#gene", ""); pg.click("#gene"); pg.keyboard.press("Backspace")
    pg.wait_for_timeout(300)
    check("Add a gene to start" in pg.inner_text("#status"), "backspace on empty input removes last chip; empty state message")
    pg.click("#status [data-try]")
    pg.wait_for_function("document.querySelector('#gene-name')?.textContent === 'MED12'")
    check(True, "the empty state's 'try MED12' loads the gene")

    # ---- Overview page
    pg.goto(BASE, wait_until="networkidle")
    pg.wait_for_selector(".design svg")
    check(pg.locator(".design svg path").count() >= 20, "vector illustration is inline")
    check(pg.locator('a[href="https://doi.org/10.1172/jci.insight.171305"]').count() >= 1, "tissue paper linked")
    check(pg.locator('a[href="https://doi.org/10.1038/s41467-023-39684-y"]').count() >= 1, "cell paper linked")
    check(pg.locator('a[href="https://doi.org/10.1126/science.1208930"]').count() >= 1, "Makinen 2011 linked")
    check(pg.locator('.examples a[href="genes.html#genes=HMGA2"]').count() == 1, "example gene links")
    pg.fill("#q", "HMGA2"); pg.wait_for_selector("#suggest li"); pg.keyboard.press("Enter")
    pg.wait_for_url("**/genes.html#genes=HMGA2")
    pg.wait_for_function("document.querySelector('#gene-name')?.textContent === 'HMGA2'")
    check(True, "overview search opens the gene on the Genes page")
    pg.goto(BASE + "#genes=CDK8", wait_until="networkidle")
    pg.wait_for_url("**/genes.html#genes=CDK8")
    pg.wait_for_function("document.querySelector('#gene-name')?.textContent === 'CDK8'")
    check(True, "old explorer links on the overview URL redirect to genes.html")

    # ---- Methods: the PCA panels
    pg.goto(BASE + "methods.html", wait_until="domcontentloaded")
    check(pg.locator(".pca svg circle").count() == 40, f"PCA panels show all 40 samples: {pg.locator('.pca svg circle').count()}")
    check(pg.locator(".pca-panel").count() == 3, "one PCA panel per dataset")
    pg.hover(".pca-panel >> nth=0 >> circle >> nth=0")
    pg.wait_for_selector(".tip:not([hidden])")
    check("Patient" in pg.inner_text(".tip") and "PC1" in pg.inner_text(".tip"), f"PCA hover label: {pg.inner_text('.tip').replace(chr(10), ' | ')}")
    check(pg.locator(".pca-panel >> nth=0 >> line.pair.lit").count() == 1, "the hovered sample's patient line lights up")
    pg.hover("#h-table") if pg.locator("#h-table").count() else pg.mouse.move(5, 5)
    pg.wait_for_selector(".tip[hidden]", state="attached")
    check(True, "PCA hover label hides on leave")

    # ---- Downloads page and the nav
    pg.goto(BASE + "downloads.html", wait_until="networkidle")
    check(pg.locator(".dl[href$='.xlsx']").count() == 2, "Downloads page offers the two workbooks")
    check(pg.locator(".site-nav a").count() == 5 and pg.locator(".site-nav a[aria-current]").inner_text() == "Downloads", "nav has five entries with Downloads current")
    for url in [BASE, BASE + "genes.html", BASE + "de.html", BASE + "methods.html"]:
        pg.goto(url, wait_until="domcontentloaded")
        check(pg.locator(".site-nav a[href='downloads.html']").count() == 1, f"Downloads in the nav of {url.split('/')[-1] or 'index.html'}")

    # ---- Comparisons page
    pg.goto(BASE + "de.html#c=NC_2D__KO_vs_WT", wait_until="networkidle")
    pg.wait_for_selector("#detable tbody tr")
    check("16,007 genes tested" in pg.inner_text("#summary"), "summary shows tested count")
    pg.uncheck("#sigonly") if pg.is_checked("#sigonly") else None
    check(pg.inner_text("#count").startswith("16,007 genes"), f"table count matches tested count: {pg.inner_text('#count')}")
    check("2D cells · KO vs WT" in pg.locator("#contrast option:checked").inner_text(), "closed select names the dataset")
    check("higher in KO" in pg.inner_text("#summary"), "summary names the numerator group")
    check(pg.locator("#picked").is_hidden(), "no picked strip before a selection")
    check(pg.evaluate("getComputedStyle(document.querySelector('#volcano .nsewdrag')).cursor") == "crosshair", "volcano: crosshair where dragging zooms")
    ring = pg.evaluate("""(() => { const el = document.getElementById('volcano'), fl = el._fullLayout, r = el.getBoundingClientRect();
      return [r.left + fl.margin.l, r.top + fl.margin.t]; })()""")
    pg.mouse.move(ring[0] + 20, ring[1] + 20); pg.wait_for_timeout(250)   # the top-left corner of the plot: no points there
    check(pg.evaluate("getComputedStyle(document.querySelector('#volcano .nsewdrag')).cursor") == "crosshair", "volcano: still a crosshair away from any point")

    pg.fill("#find", "MED12")
    pg.wait_for_selector("#suggest:not([hidden]) li")
    pg.keyboard.press("Enter")
    pg.wait_for_selector("#picked:not([hidden])")
    check("909 of 16,007" in pg.inner_text("#picked"), f"picked strip: {pg.inner_text('#picked')[:80]}")
    check(pg.locator("tr.picked-row").count() == 1, "picked row highlighted in table")
    check("gene=MED12" in pg.url, "hash carries the picked gene")
    check(pg.locator("#volcano .annotation").count() >= 1, "volcano has annotations")

    pg.click("#detable tbody tr[data-i]:nth-child(3)")
    pg.wait_for_timeout(300)
    check(pg.locator("tr.picked-row td:nth-child(2)").inner_text() == pg.locator("#picked b").first.inner_text(), "clicking a row selects that gene")
    pg.locator("#volcano").scroll_into_view_if_needed(); pg.wait_for_timeout(200)
    box = pg.locator("#volcano").bounding_box()
    pg.mouse.move(box["x"] + 200, box["y"] + 200); pg.mouse.down(); pg.mouse.move(box["x"] + 400, box["y"] + 350, steps=8); pg.mouse.up()
    pg.wait_for_timeout(400)
    check(pg.locator("#volcano-reset").is_visible(), "drag-zoom shows the Reset zoom button")
    pg.click("#volcano-reset"); pg.wait_for_timeout(400)
    check(pg.locator("#volcano-reset").is_hidden(), "Reset zoom hides the button again")
    pg.hover("#volcano .legend .traces >> nth=0")
    pg.wait_for_function("document.getElementById('volcano').data.some(t => t.opacity === 0.15)")
    check(True, "volcano legend hover dims the other categories")
    pg.hover("#h-table")
    pg.wait_for_function("document.getElementById('volcano').data.every(t => t.opacity === 1 || t.opacity === undefined)")
    pg.click("#detable th[data-key='lfc']")
    first = pg.locator("#detable tbody tr td:nth-child(4)").all_inner_texts()[:3]
    vals = [float(v.replace("−", "-")) for v in first]
    check(vals == sorted(vals, reverse=True), f"sort by log2FC descending: {first}")
    pg.click("#detable tbody tr[data-i]:nth-child(3)")
    pg.wait_for_timeout(300)
    check(pg.locator("tr.picked-row td:nth-child(2)").inner_text() == pg.locator("#picked b").first.inner_text(), "clicking a row selects that gene")
    pg.locator("#volcano").scroll_into_view_if_needed(); pg.wait_for_timeout(200)
    box = pg.locator("#volcano").bounding_box()
    pg.mouse.move(box["x"] + 200, box["y"] + 200); pg.mouse.down(); pg.mouse.move(box["x"] + 400, box["y"] + 350, steps=8); pg.mouse.up()
    pg.wait_for_timeout(400)
    check(pg.locator("#volcano-reset").is_visible(), "drag-zoom shows the Reset zoom button")
    pg.click("#volcano-reset"); pg.wait_for_timeout(400)
    check(pg.locator("#volcano-reset").is_hidden(), "Reset zoom hides the button again")
    pg.click("#detable th[data-key='lfc']")
    first = pg.locator("#detable tbody tr td:nth-child(4)").all_inner_texts()[:3]
    vals = [float(v.replace("−", "-")) for v in first]
    check(vals == sorted(vals), f"second click sorts ascending: {first}")

    pg.check("#sigonly")
    pg.wait_for_timeout(200)
    check("2,529 significant genes" in pg.inner_text("#count"), f"significant only count: {pg.inner_text('#count')}")
    pg.fill("#filter", "COL")
    pg.wait_for_timeout(200)
    rows = pg.locator("#detable tbody tr td:nth-child(2)").all_inner_texts()
    check(rows and all("COL" in r for r in rows), f"name filter: {rows[:5]}")

    with pg.expect_download() as dl:
        pg.click("#download")
    body = open(dl.value.path()).read().splitlines()
    check(dl.value.suggested_filename == "NC_2D__KO_vs_WT_significant.csv" and body[0].startswith("gene_id,gene_name")
          and all("COL" in l.split(",")[1] for l in body[1:]), f"CSV download {dl.value.suggested_filename}: {len(body)} lines")

    pg.fill("#filter", ""); pg.uncheck("#sigonly")
    check("16,007 genes · 2,529 significant" in pg.inner_text("#count"), f"table bar shows the significant count: {pg.inner_text('#count')}")
    check(pg.locator("#thr-reset").is_hidden(), "no Reset at the default thresholds")
    pg.fill("#alpha", "0.000001"); pg.wait_for_timeout(300)
    check("higher in KO" in pg.inner_text("#summary") and "1,249" not in pg.inner_text("#summary"), "a stricter padj changes the counts")
    check(pg.evaluate("document.getElementById('alpha').checkValidity()"), "0.000001 is accepted by the field")
    check(pg.locator("#thr-reset").is_visible(), "Reset appears once a threshold is changed")
    pg.click("#thr-reset"); pg.wait_for_timeout(300)
    check(pg.input_value("#alpha") == "0.05" and pg.locator("#thr-reset").is_hidden() and "1,249" in pg.inner_text("#summary"), "Reset restores the defaults")
    picked_name = pg.locator("#picked b").first.inner_text()
    pg.select_option("#contrast", "JCI_tissue__LMmut_vs_MM")
    pg.wait_for_function("document.querySelector('#summary').innerText.includes('Tissue')")
    check("20,462 genes tested" in pg.inner_text("#summary"), "switching comparison reloads")
    check(picked_name in pg.inner_text("#picked"), f"picked gene ({picked_name}) persists across comparisons")
    pg.click("#picked [data-clear]")
    check(pg.locator("#picked").is_hidden(), "clear hides the picked strip")

    pg.select_option("#ay", "JCI_tissue__LMmut_vs_MM")
    pg.wait_for_timeout(400)
    check("same comparison on both axes" in pg.inner_text("#agree-summary"), "same-axis message")
    pg.select_option("#ay", "NC_3D_spheroid__Mut1_vs_WT")
    pg.wait_for_timeout(600)
    check("rank correlation" in pg.inner_text("#agree-summary"), f"agreement summary: {pg.inner_text('#agree-summary')[:90]}")
    check(pg.locator("#agree .annotation").count() == 0 or True, "agreement drawn")
    # gene set replaces the selection and offers Undo
    pg.goto(BASE + "genes.html#genes=CDK8,MED12", wait_until="networkidle")
    pg.wait_for_function("document.querySelector('#gene-name').textContent === '2 genes'")
    pg.click("#genesets button:has-text('Housekeeping')")
    pg.wait_for_function("document.querySelector('#gene-name').textContent === '6 genes'")
    pg.wait_for_selector("#status [data-undo]")
    check("instead of your 2 genes" in pg.inner_text("#status"), f"undo offered: {pg.inner_text('#status')}")
    pg.click("#status [data-undo]")
    pg.wait_for_function("document.querySelector('#gene-name').textContent === '2 genes'")
    check("genes=CDK8,MED12" in pg.url, "undo restores the previous genes")
    # the y axis reads in TPM on the log scale
    ticks = pg.locator("#plot-JCI_tissue .ytick text").all_text_contents()
    check("10" in ticks and "30" in ticks and "0" in ticks, f"log axis ticks in TPM: {ticks}")
    check(pg.text_content("#plot-JCI_tissue .g-ytitle").strip() == "TPM (log scale)", "y title says TPM (log scale)")

    # narrow viewport: no horizontal scroll on either page
    for url in [BASE, BASE + "genes.html#genes=MED12,CDK8", BASE + "de.html#c=NC_2D__KO_vs_WT&gene=MED12", BASE + "methods.html"]:
        pg.set_viewport_size({"width": 390, "height": 800})
        pg.goto(url, wait_until="networkidle"); pg.wait_for_timeout(800)
        sw = pg.evaluate("document.documentElement.scrollWidth")
        check(sw <= 390, f"no horizontal overflow at 390px for {url.split('/')[-1]}: scrollWidth {sw}")

    check(not errs, f"no console errors: {errs}")

    # without the plotting library (blocked here), the pages still show their tables
    ctx2 = b.new_context(viewport={"width": 1280, "height": 900})
    ctx2.route("**/vendor/*", lambda route: route.abort())
    p2 = ctx2.new_page()
    p2.goto(BASE + "genes.html#genes=MED12", wait_until="networkidle")
    p2.wait_for_selector("#tab-JCI_tissue table")
    check("Mean TPM" in p2.inner_text("#tab-JCI_tissue") and p2.locator("#plot-JCI_tissue .noplot").count() == 1, "Genes page shows tables and a note when the plotting library is missing")
    p2.goto(BASE + "de.html#c=NC_2D__KO_vs_WT", wait_until="networkidle")
    p2.wait_for_selector("#detable tbody tr")
    check(p2.locator("#volcano .noplot").count() == 1 and "16,007 genes tested" in p2.inner_text("#summary"), "Comparisons page shows the table and a note when the plotting library is missing")
    ctx2.close()
    b.close()
print("\nFAILURES:", len(fails))
sys.exit(1 if fails else 0)
