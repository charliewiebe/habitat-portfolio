"""Map geometry in UTM 10N (EPSG:32610): land, depth contours, regulatory polygons, ports, graticule.

Usage: python -I pipeline/03_geometry.py <data dir>
"""
import sys
import json
import numpy as np
import pandas as pd
from pyproj import Transformer
from global_land_mask import globe
from skimage import measure
from shapely.geometry import Polygon, LineString, Point
from shapely.ops import unary_union
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

DATA = sys.argv[1]
to_utm = Transformer.from_crs(4326, 32610, always_xy=True)
LON0, LON1, LAT0, LAT1 = -126.0, -123.2, 41.75, 46.45

def dm(deg, minutes):
    return deg + minutes / 60.0

def utm_ring(lonlat):
    x, y = to_utm.transform([p[0] for p in lonlat], [p[1] for p in lonlat])
    return [[round(a, 0), round(b, 0)] for a, b in zip(x, y)]

geom = {}

# ---- land (GLOBE 30" land mask from global-land-mask) -------------------
res = 1 / 120
lats = np.arange(LAT0, LAT1, res)
lons = np.arange(LON0, LON1, res)
LL, LA = np.meshgrid(lons, lats)
mask = globe.is_land(LA, LL).astype(float)
pad = np.pad(mask, 1, constant_values=0)   # zero border => every land contour is closed
polys = []
for c in measure.find_contours(pad, 0.5):
    r, cidx = c[:, 0] - 1, c[:, 1] - 1
    lat = LAT0 + r * res
    lon = LON0 + cidx * res
    if len(c) < 20:
        continue
    poly = Polygon(zip(lon, lat))
    if not poly.is_valid:
        poly = poly.buffer(0)
    polys.append(poly)
land = unary_union(polys).simplify(0.0012)
land_rings = []
for p in getattr(land, "geoms", [land]):
    if p.area < 0.002:
        continue
    land_rings.append(utm_ring(list(p.exterior.coords)))
geom["land"] = land_rings
print(f"land rings: {len(land_rings)}, vertices {sum(len(r) for r in land_rings)}")

# ---- depth contours from the survey grid's depth field ------------------
g = pd.read_csv(f"{DATA}/grid_oregon.csv")
tri = mtri.Triangulation(g.x.values, g.y.values)
tri.set_mask(np.array([np.max(np.hypot(*(np.diff(np.r_[g.x.values[t], g.x.values[t][:1]]),
                                          np.diff(np.r_[g.y.values[t], g.y.values[t][:1]])))) > 4000
                       for t in tri.triangles]))
cs = plt.tricontour(tri, g.depth_m.values, levels=[100, 200, 500, 1000])
contours = []
for level, segs in zip(cs.levels, cs.allsegs):
    for seg in segs:
        if len(seg) < 6:
            continue
        ls = LineString(seg).simplify(400)
        contours.append({"depth": int(level), "pts": [[round(a), round(b)] for a, b in ls.coords]})
geom["contours"] = contours
print(f"depth contour lines: {len(contours)}")

# ---- Yelloweye Rockfish Conservation Areas off Oregon (50 CFR 660.70(g)-(m)) ----
W = lambda d, m: -dm(d, m)
YRCA = {
    "Tillamook": [(45, 40.96, 124, 27.52), (45, 40.96, 124, 19.99), (45, 34.44, 124, 14.48),
                  (45, 31.93, 124, 14.05), (45, 31.84, 124, 22.04), (45, 36.95, 124, 24.45),
                  (45, 38.89, 124, 25.92)],
    "Newport": [(44, 46.00, 124, 32.57), (44, 46.00, 124, 32.00), (44, 42.00, 124, 30.00),
                (44, 39.00, 124, 30.00), (44, 39.00, 124, 34.00), (44, 43.16, 124, 34.00),
                (44, 44.54, 124, 33.58)],
    "Stonewall Bank": [(44, 37.46, 124, 24.92), (44, 37.46, 124, 23.63), (44, 28.71, 124, 21.80),
                       (44, 28.71, 124, 24.10), (44, 31.42, 124, 25.47)],
    "Stonewall Bank Exp. 1": [(44, 41.76, 124, 30.02), (44, 41.73, 124, 21.60), (44, 25.25, 124, 16.94),
                              (44, 25.29, 124, 30.14)],
    "Stonewall Bank Exp. 2": [(44, 38.54, 124, 27.41), (44, 38.54, 124, 23.86), (44, 27.13, 124, 21.50),
                              (44, 27.13, 124, 26.89), (44, 31.30, 124, 28.35)],
    "Florence": [(44, 30.04, 124, 42.31), (44, 30.19, 124, 40.46), (44, 25.00, 124, 37.00),
                 (44, 25.00, 124, 45.00), (44, 26.71, 124, 45.00)],
    "Heceta Bank": [(44, 16.28, 124, 47.86), (44, 15.38, 124, 49.86), (44, 14.49, 124, 51.82),
                    (44, 14.01, 124, 52.88), (44, 13.47, 124, 54.08), (44, 12.72, 124, 54.07),
                    (44, 11.53, 124, 54.06), (44, 8.72, 124, 54.02), (44, 6.68, 124, 54.00),
                    (44, 5.34, 124, 53.10), (44, 2.88, 124, 53.96), (44, 2.18, 124, 54.29),
                    (44, 0.14, 124, 55.25), (43, 58.36, 124, 55.42), (43, 57.68, 124, 55.48),
                    (43, 56.66, 124, 55.45), (43, 56.65, 124, 55.49), (43, 56.64, 124, 56.53),
                    (43, 56.74, 124, 56.74), (43, 59.18, 124, 56.94), (44, 0.45, 124, 56.35),
                    (44, 2.34, 124, 55.49), (44, 4.81, 124, 55.65), (44, 6.45, 124, 55.78),
                    (44, 8.47, 124, 55.93), (44, 9.85, 124, 56.04), (44, 11.34, 124, 56.16),
                    (44, 12.92, 124, 56.28), (44, 14.06, 124, 55.10), (44, 15.32, 124, 53.79),
                    (44, 16.90, 124, 52.16), (44, 16.96, 124, 52.11), (44, 16.96, 124, 51.95),
                    (44, 17.02, 124, 48.02), (44, 17.02, 124, 47.47)],
}
CFR_PARA = {"Tillamook": "(g)", "Newport": "(h)", "Stonewall Bank": "(i)", "Stonewall Bank Exp. 1": "(j)",
            "Stonewall Bank Exp. 2": "(k)", "Florence": "(l)", "Heceta Bank": "(m)"}

# ---- Oregon marine reserves (OAR 141-142) --------------------------------
RESERVES = {
    "Cape Falcon": [(45, 48.00, 123, 58.10), (45, 46.80, 123, 58.17), (45, 45.79, 123, 58.41),
                    (45, 45.44, 123, 58.03), (45, 44.62, 123, 57.82), (45, 44.00, 123, 56.77),
                    (45, 44.00, 124, 1.50), (45, 48.00, 124, 1.50)],
    "Cascade Head": [(45, 2.00, 124, 0.52), (45, 1.98, 124, 0.65), (45, 1.53, 124, 0.74),
                     (44, 59.25, 124, 0.89), (44, 59.25, 124, 4.50), (45, 2.00, 124, 4.50)],
    "Otter Rock": [(44, 45.175, 124, 4.53), (44, 45.175, 124, 3.8583), (44, 44.7933, 124, 3.7833),
                   (44, 44.7283, 124, 3.7833), (44, 43.315, 124, 3.6567), (44, 43.315, 124, 4.2)],
    "Cape Perpetua": [(44, 17.00, 124, 6.73), (44, 16.82, 124, 6.69), (44, 16.15, 124, 6.56),
                      (44, 15.81, 124, 6.62), (44, 14.67, 124, 6.84), (44, 13.50, 124, 6.74),
                      (44, 13.50, 124, 11.00), (44, 17.00, 124, 11.00)],
    "Redfish Rocks": [(42, 42.54, 124, 29.64), (42, 42.96, 124, 27.78), (42, 41.52, 124, 27.18),
                      (42, 41.16, 124, 28.86)],
}
OAR = {"Cape Falcon": "141-142-0085", "Cascade Head": "141-142-0065", "Otter Rock": "141-142-0030",
       "Cape Perpetua": "141-142-0045", "Redfish Rocks": "141-142-0035"}

def to_lonlat(pts):
    return [(W(lo_d, lo_m), dm(la_d, la_m)) for la_d, la_m, lo_d, lo_m in pts]

gpts = np.column_stack([g.x.values, g.y.values])

def cells_inside(lonlat):
    poly_utm = Polygon(utm_ring(lonlat))
    # a cell counts as inside if its 2-km footprint overlaps the polygon by >= 25%
    inside = []
    for i, (x, y) in enumerate(gpts):
        if abs(x - poly_utm.centroid.x) > 60000 or abs(y - poly_utm.centroid.y) > 60000:
            continue
        cell = Polygon([(x - 1000, y - 1000), (x + 1000, y - 1000), (x + 1000, y + 1000), (x - 1000, y + 1000)])
        if cell.intersection(poly_utm).area >= 0.25 * cell.area:
            inside.append(i)
    return inside, poly_utm.area / 1e6

closures = []
for name, pts in YRCA.items():
    ll = to_lonlat(pts)
    cells, area = cells_inside(ll)
    closures.append({"id": "yrca_" + name.lower().replace(" ", "_").replace(".", ""), "name": f"{name} YRCA",
                     "kind": "yrca", "cite": f"50 CFR 660.70{CFR_PARA[name]}", "area_km2": round(area, 1),
                     "ring": utm_ring(ll), "cells": cells})
for name, pts in RESERVES.items():
    ll = to_lonlat(pts)
    cells, area = cells_inside(ll)
    closures.append({"id": "mr_" + name.lower().replace(" ", "_"), "name": f"{name} Marine Reserve",
                     "kind": "reserve", "cite": f"OAR {OAR[name]}", "area_km2": round(area, 1),
                     "ring": utm_ring(ll), "cells": cells})
for c in closures:
    print(f"  {c['name']:32s} {c['area_km2']:7.1f} km2  survey-grid cells {len(c['cells'])}")
geom["closures"] = closures

# ---- ports and labels ----------------------------------------------------
PORTS = [("Astoria", 46.188, -123.831), ("Garibaldi", 45.558, -123.911), ("Pacific City", 45.202, -123.962),
         ("Depoe Bay", 44.810, -124.063), ("Newport", 44.637, -124.053), ("Florence", 43.983, -124.100),
         ("Winchester Bay", 43.680, -124.180), ("Charleston", 43.345, -124.330), ("Bandon", 43.119, -124.408),
         ("Port Orford", 42.745, -124.497), ("Gold Beach", 42.407, -124.421), ("Brookings", 42.053, -124.284)]
px, py = to_utm.transform([p[2] for p in PORTS], [p[1] for p in PORTS])
geom["ports"] = [{"name": n, "x": round(a), "y": round(b)} for (n, _, _), a, b in zip(PORTS, px, py)]

# ---- graticule -----------------------------------------------------------
grat = []
for la in range(42, 47):
    lo = np.linspace(-126.5, -123.0, 30); x, y = to_utm.transform(lo, np.full_like(lo, la))
    grat.append({"label": f"{la}°N", "kind": "lat", "pts": [[round(a), round(b)] for a, b in zip(x, y)]})
for lo_ in [-126, -125.5, -125, -124.5, -124]:
    la = np.linspace(41.5, 46.6, 30); x, y = to_utm.transform(np.full_like(la, lo_), la)
    lab = f"{abs(lo_):g}°W"
    grat.append({"label": lab, "kind": "lon", "pts": [[round(a), round(b)] for a, b in zip(x, y)]})
geom["graticule"] = grat

with open(f"{DATA}/geometry.json", "w") as fh:
    json.dump(geom, fh)
print("wrote geometry.json")
