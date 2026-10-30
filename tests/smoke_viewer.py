"""Browser smoke test of the viewer (installed Google Chrome via Playwright; dev only).

  python -m http.server 8777 --bind 127.0.0.1      # from the repo root, in another shell
  python tests/smoke_viewer.py [--url http://127.0.0.1:8777/docs/atlas/index.html] [--out tests/shots]

Loads the viewer, fails on console errors, selects a gene and a region, switches views and modes and
saves screenshots.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from playwright.sync_api import sync_playwright


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default="http://127.0.0.1:8777/docs/atlas/index.html")
    ap.add_argument("--out", default=str(Path(__file__).parent / "shots"))
    ap.add_argument("--gene", default="Prox1")
    ap.add_argument("--brain", default="yao638850")
    ap.add_argument("--region", default="DG")
    ap.add_argument("--expr", default="", help="expression host to test, e.g. the Hugging Face resolve URL")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    errors = []
    with sync_playwright() as p:
        b = p.chromium.launch(channel="chrome", headless=True, args=["--use-angle=swiftshader", "--enable-unsafe-swiftshader"])
        page = b.new_page(viewport={"width": 1500, "height": 900})
        page.on("console", lambda m: (print(f"[console.{m.type}] {m.text}"), errors.append(m.text) if m.type == "error" else None))
        page.on("pageerror", lambda e: (print(f"[pageerror] {e}"), errors.append(str(e))))
        t0 = time.time()
        extra = f"&expr={args.expr}" if args.expr else ""
        page.goto(f"{args.url}?brain={args.brain}{extra}", wait_until="domcontentloaded")
        page.wait_for_function("document.getElementById('loading').style.opacity === '0'", timeout=180_000)
        print(f"loaded in {time.time() - t0:.1f}s; hud: {page.inner_text('#hud')}")
        page.wait_for_timeout(1500)
        page.screenshot(path=str(out / "01_3d_types.png"))

        page.fill("#gene", args.gene)
        page.keyboard.press("Enter")
        page.wait_for_function(f"document.getElementById('genestatus').textContent.startsWith('{args.gene}:')", timeout=240_000)
        print("gene:", page.inner_text("#genestatus"))
        page.wait_for_timeout(1500)
        page.screenshot(path=str(out / "02_3d_gene_overlay.png"))

        page.fill("#region", args.region)
        page.wait_for_timeout(300)
        page.keyboard.press("Enter")
        page.wait_for_timeout(3000)
        print("profile:", page.inner_text("#profile")[:600].replace("\n", " | "))
        page.screenshot(path=str(out / "03_region_profile.png"))

        page.click("#view button[data-v='coronal']")
        page.wait_for_timeout(1000)
        page.wait_for_function("document.getElementById('genestatus').textContent.indexOf('template') < 0", timeout=240_000)
        page.wait_for_timeout(2500)
        page.screenshot(path=str(out / "04_coronal_overlay.png"))
        box = page.locator("#slice").bounding_box()
        page.mouse.move(box["x"] + box["width"] * 0.55, box["y"] + box["height"] * 0.4)
        page.wait_for_timeout(500)
        print("hover:", page.inner_text("#tip").replace("\n", " | ")[:400])

        page.click("#mode button[data-v='split']")
        page.wait_for_timeout(2500)
        page.screenshot(path=str(out / "05_coronal_split.png"))
        page.click("#view button[data-v='3d']")
        page.wait_for_timeout(1500)
        page.screenshot(path=str(out / "06_3d_split.png"))
        print("url:", page.url)
        print("gene regions:", page.inner_text("#generegions")[:300].replace("\n", " | "))
        root = args.url.split("/docs/")[0] + "/docs/"
        for name in ("index.html", "qc/index.html"):
            pg = b.new_page(viewport={"width": 1300, "height": 900})
            pg.on("console", lambda m: (print(f"[{name} console.{m.type}] {m.text}"), errors.append(m.text) if m.type == "error" else None))
            pg.on("pageerror", lambda e: (print(f"[{name} pageerror] {e}"), errors.append(str(e))))
            pg.goto(root + name, wait_until="networkidle")
            pg.wait_for_timeout(1500)
            pg.screenshot(path=str(out / f"page_{name.replace('/', '_').replace('.html', '')}.png"), full_page=True)
            print(f"{name}: title={pg.title()!r}")
        b.close()
    if errors:
        print(f"FAILED: {len(errors)} console errors")
        sys.exit(1)
    print("OK")


if __name__ == "__main__":
    main()
