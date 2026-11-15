"""Inline model outputs, geometry, parameters and the engine into one self-contained page.

Usage: python -I pipeline/05_build_app.py <data dir> <app dir> [<standalone index.html>]
Writes <app dir>/index.html (body fragment, as published to claude.ai) and, if a third path is
given, a full standalone HTML document for GitHub Pages with links to the repo and portfolio.
"""
import sys
import json
import pandas as pd
from pyproj import Transformer

DATA, APP = sys.argv[1], sys.argv[2]
to_utm = Transformer.from_crs(4326, 32610, always_xy=True)

mo = json.load(open(f"{DATA}/model_outputs.json"))
geo = json.load(open(f"{DATA}/geometry.json"))
params = json.load(open(f"{DATA}/species_params.json"))
hauls = pd.read_csv(f"{DATA}/hauls_oregon.csv")

g = mo["grid"]
grid = {
    "n": len(g["x"]),
    "x": [round(v) for v in g["x"]], "y": [round(v) for v in g["y"]],
    "depth": [round(v) for v in g["depth_m"]], "area": [round(v) for v in g["area_ha"]],
    "tclim": [round(v, 2) for v in g["t_clim"]],
    "lon": [round(v, 4) for v in g["lon"]], "lat": [round(v, 4) for v in g["lat"]],
}
species = {}
for k, o in mo["species"].items():
    species[k] = {kk: o[kk] for kk in ["n_hauls", "n_pos", "metrics", "t_range", "eta_p", "eta_g", "fp_t", "fg_t",
                                       "fp_draws", "fg_draws", "curves", "obs_by_year", "sigma2_log"]}
tg = mo["t_grid"]
haul_out = {
    "x": hauls.x.round().astype(int).tolist(), "y": hauls.y.round().astype(int).tolist(),
    "year": hauls.year.astype(int).tolist(), "depth": hauls.depth_m.round().astype(int).tolist(),
}
for k in species:
    haul_out[k] = hauls[f"{k}_cpue"].round(2).tolist()

# frames: split the coast at Heceta Head (44.14 N), as coastal charts often do
_, y_split = to_utm.transform(-124.13, 44.14)
x0, x1 = 308000, 468000
_, y_top = to_utm.transform(-124.0, 46.30)
_, y_bot = to_utm.transform(-124.0, 41.97)
frames = [
    {"id": "north", "label": "North coast", "sub": "Columbia River to Heceta Head", "x0": x0, "x1": x1,
     "y0": round(y_split - 6000), "y1": round(y_top)},
    {"id": "south", "label": "South coast", "sub": "Heceta Head to California", "x0": x0, "x1": x1,
     "y0": round(y_bot), "y1": round(y_split + 6000)},
]
features = []
for name, lon, lat in [("Heceta Bank", -125.02, 43.93), ("Stonewall Bank", -124.47, 44.52),
                       ("P A C I F I C   O C E A N", -125.22, 45.62)]:
    x, y = to_utm.transform(lon, lat)
    features.append({"name": name, "x": round(x), "y": round(y)})

payload = {
    "grid": grid, "species": species, "tgrid": {"t0": tg[0], "dt": round(tg[1] - tg[0], 3), "n": len(tg)},
    "tempAnomByYear": mo["temp_anomaly_by_year"], "years": mo["years"], "hauls": haul_out,
    "geo": geo, "params": params, "frames": frames, "features": features,
    "meta": {"nHauls": int(len(hauls)), "yearMin": int(hauls.year.min()), "yearMax": int(hauls.year.max()),
             "nCells": grid["n"], "built": pd.Timestamp.now().strftime("%Y-%m-%d")},
}
data_js = "const DATA = " + json.dumps(payload, separators=(",", ":")) + ";"

tpl = open(f"{APP}/template.html").read()
engine = open(f"{APP}/engine.js").read()
app = open(f"{APP}/app.js").read()
html = tpl.replace("/*__ENGINE__*/", engine).replace("/*__DATA__*/", data_js).replace("/*__APP__*/", app)
with open(f"{APP}/index.html", "w") as fh:
    fh.write(html.replace("    <!--__LINKS__-->\n", ""))
print(f"index.html {len(html) / 1e6:.2f} MB (data {len(data_js) / 1e6:.2f} MB)")

if len(sys.argv) > 3:
    REPO = "https://github.com/charliewiebe/habitat-portfolio/tree/main/oregon-groundfish-habitat"
    PORTFOLIO = "../"
    links = (f'    <nav class="links" aria-label="Project links"><a href="{REPO}">Code</a>'
             f'<a href="{PORTFOLIO}">Portfolio</a></nav>\n')
    head_end = html.index("</style>") + len("</style>")
    head, body = html[:head_end], html[head_end:]
    head = head.replace('<meta name="claude:viewport-fit-sides" content="cover">\n', "")
    reset = ("<style>:root{padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}"
             "body{margin:0}img{max-width:100%}[hidden]{display:none!important}</style>\n")
    doc = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
           '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
           + reset + head + "\n</head>\n<body>" + body.replace("    <!--__LINKS__-->\n", links) + "\n</body>\n</html>\n")
    with open(sys.argv[3], "w") as fh:
        fh.write(doc)
    print(f"standalone {sys.argv[3]} {len(doc) / 1e6:.2f} MB")
