"""Screenshot pages with one browser.
    uv run --with playwright python tests/screenshots.py out_dir "url|name|width" ...
"""
import sys, glob
from playwright.sync_api import sync_playwright

out = sys.argv[1]
jobs = [a.split("|") for a in sys.argv[2:]]
with sync_playwright() as p:
    try: b = p.chromium.launch()
    except Exception:
        exe = sorted(glob.glob("/home/*/.cache/ms-playwright/chromium_headless_shell-*/chrome-headless-shell-linux64/chrome-headless-shell"))[-1]
        b = p.chromium.launch(executable_path=exe)
    for url, name, w in jobs:
        pg = b.new_page(viewport={"width": int(w), "height": 900}, device_scale_factor=1)
        errs = []
        pg.on("console", lambda m: m.type == "error" and errs.append(m.text))
        pg.on("pageerror", lambda e: errs.append(str(e)))
        pg.goto(url, wait_until="networkidle")
        pg.wait_for_timeout(1200)
        pg.screenshot(path=f"{out}/{name}.png", full_page=True)
        print(name, "errors:", errs or "none")
        pg.close()
    b.close()
