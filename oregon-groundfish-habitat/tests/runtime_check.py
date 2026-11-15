"""Headless smoke test: python tests/runtime_check.py index.html /tmp/test.html  (needs playwright + chromium)"""
import os, sys, asyncio
from playwright.async_api import async_playwright
SRC = sys.argv[1]; TMP = sys.argv[2]
html = open(SRC, encoding="utf-8").read()
open(TMP, "w", encoding="utf-8").write(html)
async def main():
    async with async_playwright() as p:
        b = await p.chromium.launch(**({"executable_path": os.environ["CHROME"]} if os.environ.get("CHROME") else {}))
        for name, vp in [("desktop", (1440, 900)), ("tablet", (1000, 800)), ("phone", (400, 860))]:
            ctx = await b.new_context(viewport={"width": vp[0], "height": vp[1]}); pg = await ctx.new_page()
            errs = []
            pg.on("pageerror", lambda e: errs.append(str(e)))
            pg.on("console", lambda m: errs.append(m.text) if m.type == "error" and "fonts" not in m.text and "TUNNEL" not in m.text else None)
            await pg.goto("file://" + TMP); await pg.wait_for_timeout(1200)
            for sp in ["canary", "sablefish", "lingcod"]:
                await pg.click(f'[data-sp="{sp}"]'); await pg.wait_for_timeout(150)
            for lay in ["spawn", "change", "fishing", "temp", "habitat"]:
                await pg.click(f'[data-layer="{lay}"]'); await pg.wait_for_timeout(80)
            await pg.click('[data-seg="band"] button[data-v="all"]'); await pg.wait_for_timeout(150)
            ow = await pg.evaluate("document.documentElement.scrollWidth - window.innerWidth")
            info = await pg.evaluate("(() => { const c=document.getElementById('map'); const s=document.getElementById('stage'); return [c.style.width, c.style.height, s.style.height]; })()")
            lede = await pg.evaluate("document.querySelector('.lede').innerText")
            print(name, "overflow", ow, "canvas/stage", info, "errors", errs[:3]); print("   ", lede)
            await ctx.close()
        await b.close()
asyncio.run(main())
