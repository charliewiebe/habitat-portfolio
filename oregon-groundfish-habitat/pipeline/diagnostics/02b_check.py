"""Diagnostics: predicted vs observed density at haul locations, by depth bin and year."""
import sys, json, numpy as np, pandas as pd
from scipy.spatial import cKDTree
D = sys.argv[1]
mo = json.load(open(f"{D}/model_outputs.json"))
g = pd.DataFrame(mo["grid"]); h = pd.read_csv(f"{D}/hauls_oregon.csv").dropna(subset=["bottom_temp_c"])
tree = cKDTree(g[["x", "y"]].values); dist, idx = tree.query(h[["x", "y"]].values)
print("haul->cell distance (km) median/95%:", np.round(np.percentile(dist, [50, 95]) / 1000, 2))
T = np.array(mo["t_grid"])
for sp, o in mo["species"].items():
    ti = np.clip(np.round((g.t_clim.values - T[0]) / 0.1).astype(int), 0, len(T) - 1)
    dens = 1 / (1 + np.exp(-(np.array(o["eta_p"]) + np.array(o["fp_t"])[ti]))) * np.exp(np.array(o["eta_g"]) + np.array(o["fg_t"])[ti])
    h["pred"] = dens[idx]
    h["dbin"] = pd.cut(h.depth_m, [55, 100, 150, 200, 300, 500, 800, 1300])
    t = h.groupby("dbin", observed=True).agg(obs=(f"{sp}_cpue", "mean"), pred=("pred", "mean"), n=("pred", "size"))
    print(f"\n{sp}  sigma2_log={o['sigma2_log']}  obs mean {h[f'{sp}_cpue'].mean():.2f}  pred@hauls {h.pred.mean():.2f}  max obs {h[f'{sp}_cpue'].max():.1f}")
    print(t.round(2).T.to_string())
