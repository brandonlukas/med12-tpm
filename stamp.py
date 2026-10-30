#!/usr/bin/env python3
"""Stamp the site's own stylesheet and script links with a hash of their contents, so a
browser fetches a changed file instead of serving a cached one (GitHub Pages lets
browsers keep files for 10 minutes). Run before committing: python3 stamp.py"""
import hashlib
import re
from pathlib import Path

docs = Path(__file__).resolve().parent / "docs"
assets = ["style.css", "common.js", "app.js", "de.js", "overview.js", "methods.js"]
stamp = hashlib.sha1(b"".join((docs / a).read_bytes() for a in assets)).hexdigest()[:8]
pattern = re.compile(r'((?:href|src)=")(' + "|".join(map(re.escape, assets)) + r')(?:\?v=[0-9a-f]+)?"')
for page in docs.glob("*.html"):
    text = page.read_text()
    new = pattern.sub(lambda m: f'{m.group(1)}{m.group(2)}?v={stamp}"', text)
    if new != text:
        page.write_text(new)
print("stamp", stamp)
