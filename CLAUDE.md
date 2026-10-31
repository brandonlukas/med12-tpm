# CLAUDE.md

Public MED12 RNA-seq explorer: a pipeline (nf-core/rnaseq → DESeq2 → TF activity) and a static site in `docs/`
served by GitHub Pages from `main`. Read `README.md` for the files and the run order; this file holds the rules and
the traps.

## Rules

- **Old numbers never move.** A data release (`v1.x`) adds datasets or analyses without recomputing earlier ones:
  each batch is its own nf-core run, `merge_salmon.R` places them side by side. After any rerun, check that
  earlier outputs are byte-identical. A new reference or pipeline version means reprocessing everything: `v2.0`.
- **`datasets.json` is the one registry** (datasets, comparisons, labels, release, batch folders, TF network
  file). Scripts and the site read it; never hard-code a dataset, comparison or file name elsewhere.
- **Genotypes:** never relabel from GEO. A correction needs evidence (`checks/`), and `samples.tsv` keeps the
  original label in its own column (`genotype_as_listed_on_GEO_(incorrect)`, `genotype_as_listed_in_paper`).
- **The TF network is frozen** (`networks/`, named in `datasets.json`). Never re-fetch silently; a new copy is a
  new data version.
- **One significance rule** for genes: p<sub>adj</sub> < 0.05 and |log2FC| ≥ 1 (`SIG` in `docs/common.js`).
  TFs: p<sub>adj</sub> < 0.05.
- **Generated page sections** sit between `<!-- name:start -->` and `<!-- name:end -->` markers
  (`page_tables.py`, `pca_svg.py`, `illustration/vector.py`). Edit the generator, not the page.

## Site writing

Plain language for biologists: no package or test names outside Methods, no code words, US spelling, colors named
plainly ("orange"). Less is more: cut repetition, prefer simple interactions to ornate ones. Group order is
MM → LM WT → LM mut. One meaning per color across the site.

## Changing the site

1. Edit `docs/`, then `python3 stamp.py` (content hashes on CSS/JS links; Pages caches for 10 minutes).
2. Serve `docs/` with `python3 -m http.server 8765` and run
   `uv run --with playwright python tests/site_test.py`. All checks must pass. Look at screenshots at desktop
   and 390 px width; nothing may scroll sideways.
3. Commit; push to `main` only when the user says so. Pages rebuilds in about a minute.

## Releasing data (`vX.Y`)

1. Bump `release` in `datasets.json`; add a `CHANGELOG.md` section (its first paragraph appears on Downloads).
2. Run the affected scripts, then `python3 page_tables.py` and `python3 stamp.py`; run the tests.
3. Merge to `main` and push. Download links point at `releases/download/<release>/`, so **every release carries
   every file**: `deliverables/*.csv` gzipped, `deliverables/de/*.csv` and `deliverables/tf/*.csv` gzipped (TF
   tables renamed `<dataset>__<contrast>_TF_activity.csv.gz`), `sample_info.csv`, the three workbooks, and the
   network file. Check every link on the built Downloads page resolves to an asset before creating the release.
4. `gh release create`; tag the site `site-YYYY-MM-DD`. Do not push a branch with the same name as a tag.

## Traps

- `pkill -f`/`pgrep -f` with a pattern that also appears in your own command kills your own shell (exit 144).
  Kill by PID.
- In Playwright, `goto` to the same page with only a new `#hash` does not reload it; call `reload()`.
- Test runs against the live site are flaky on this machine's network; rerun a timed-out step on its own before
  suspecting the site.
- 42-bp reads (Moyo, Leistico) are too short for minimap2; the checks align with bwa mem.
