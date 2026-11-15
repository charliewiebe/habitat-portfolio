"""Build the Oregon haul table and prediction grid from surveyjoin's NWFSC files.

Usage: python -I pipeline/01_prepare.py <surveyjoin-data dir> <surveyjoin/data dir> <out dir>
"""
import sys
import numpy as np
import pandas as pd
import pyreadr
from pyproj import Transformer

SRC, GRID_DIR, OUT = sys.argv[1], sys.argv[2], sys.argv[3]

LAT_MIN, LAT_MAX = 42.0, 46.25          # California border to Columbia River mouth
SPECIES = {                             # ITIS codes in surveyjoin's species dictionary
    "lingcod": 167116,
    "canary": 166734,
    "sablefish": 167123,
}
to_utm = Transformer.from_crs(4326, 32610, always_xy=True)

# ---- hauls -------------------------------------------------------------
haul = pyreadr.read_r(f"{SRC}/nwfsc-haul.rds")[None]
haul = haul[haul.survey_name == "NWFSC.Combo"].copy()
# tow midpoint; fall back to the start position when the end position is implausible
# (one 2023 haul has its end longitude recorded with the sign dropped)
bad_end = ((haul.lat_start - haul.lat_end).abs() > 0.2) | ((haul.lon_start - haul.lon_end).abs() > 0.2)
haul["lat"] = np.where(bad_end, haul.lat_start, haul[["lat_start", "lat_end"]].mean(axis=1))
haul["lon"] = np.where(bad_end, haul.lon_start, haul[["lon_start", "lon_end"]].mean(axis=1))
print(f"hauls with implausible end positions (start used): {int(bad_end.sum())}")
haul = haul[(haul.lat >= LAT_MIN) & (haul.lat <= LAT_MAX)]
haul = haul.dropna(subset=["depth_m", "effort", "lat", "lon"])
haul["event_id"] = haul.event_id.astype("int64")
haul["year"] = haul.year.astype(int)
haul["x"], haul["y"] = to_utm.transform(haul.lon.values, haul.lat.values)

catch = pyreadr.read_r(f"{SRC}/nwfsc-catch.rds")[None]
catch["event_id"] = catch.event_id.astype("int64")
catch = catch.dropna(subset=["itis"])
catch["itis"] = catch.itis.astype(float).astype(int)

cols = ["event_id", "year", "date", "lat", "lon", "x", "y", "depth_m", "bottom_temp_c", "effort"]
out = haul[cols].copy()
for name, itis in SPECIES.items():
    c = catch[catch.itis == itis][["event_id", "catch_weight"]]
    c = c.groupby("event_id", as_index=False).catch_weight.sum()
    m = out[["event_id"]].merge(c, on="event_id", how="left")
    # surveyjoin's catch table is zero-filled for these species; a missing row means not recorded
    out[f"{name}_kg"] = m.catch_weight.values
    out[f"{name}_cpue"] = out[f"{name}_kg"] / out.effort           # kg per hectare

n0 = len(out)
out = out.dropna(subset=[f"{s}_kg" for s in SPECIES])
print(f"hauls: {n0} in Oregon band, {len(out)} with catch records, years {out.year.min()}-{out.year.max()}")
for s in SPECIES:
    pos = (out[f"{s}_kg"] > 0).mean()
    print(f"  {s:10s} encounter rate {pos:.3f}  positive hauls {(out[f'{s}_kg'] > 0).sum()}")
out.to_csv(f"{OUT}/hauls_oregon.csv", index=False)

# ---- prediction grid ---------------------------------------------------
grid = list(pyreadr.read_r(f"{GRID_DIR}/nwfsc_grid.rda").values())[0]
g = grid[(grid.survey == "NWFSC.Combo") & (grid.split_state == "O") & (grid.area > 0)].copy()
g = g[(g.lat >= LAT_MIN) & (g.lat <= LAT_MAX)]
g["x"], g["y"] = to_utm.transform(g.lon.values, g.lat.values)
g = g.rename(columns={"area": "area_ha"})[["lon", "lat", "x", "y", "depth_m", "area_ha"]].reset_index(drop=True)
print(f"grid cells (survey domain, Oregon): {len(g)}  depth {g.depth_m.min():.0f}-{g.depth_m.max():.0f} m")
g.to_csv(f"{OUT}/grid_oregon.csv", index=False)
