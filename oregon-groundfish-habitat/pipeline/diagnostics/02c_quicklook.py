import sys, json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
mo = json.load(open(f"{sys.argv[1]}/model_outputs.json")); out = sys.argv[2]
g = mo["grid"]; x = np.array(g["x"]) / 1000; y = np.array(g["y"]) / 1000; T = np.array(mo["t_grid"])
ti = np.clip(np.round((np.array(g["t_clim"]) - T[0]) / 0.1).astype(int), 0, len(T) - 1)
fig, ax = plt.subplots(1, 5, figsize=(16, 9))
sc = ax[0].scatter(x, y, c=g["t_clim"], s=3, cmap="RdYlBu_r", marker="s"); ax[0].set_title("bottom T clim"); plt.colorbar(sc, ax=ax[0])
sc = ax[1].scatter(x, y, c=np.log10(g["depth_m"]), s=3, marker="s"); ax[1].set_title("log10 depth"); plt.colorbar(sc, ax=ax[1])
for i, sp in enumerate(mo["species"]):
    o = mo["species"][sp]
    d = 1/(1+np.exp(-(np.array(o["eta_p"]) + np.array(o["fp_t"])[ti]))) * np.exp(np.array(o["eta_g"]) + np.array(o["fg_t"])[ti])
    sc = ax[2+i].scatter(x, y, c=np.log10(d + 1e-3), s=3, cmap="viridis", marker="s"); ax[2+i].set_title(sp + " log10 dens"); plt.colorbar(sc, ax=ax[2+i])
for a in ax: a.set_aspect("equal")
plt.tight_layout(); plt.savefig(f"{out}/quicklook_maps.png", dpi=80)
fig, ax = plt.subplots(2, 3, figsize=(13, 7))
for i, sp in enumerate(mo["species"]):
    c = mo["species"][sp]["curves"]
    for j, (k, lab) in enumerate([("depth", "depth (m)"), ("temp_p", "bottom T (C)")]):
        a = ax[j, i]; cc = c[k]; a.plot(cc["x"], cc["y"]); a.fill_between(cc["x"], cc["lo"], cc["hi"], alpha=.3)
        if k == "temp_p": a.plot(c["temp_g"]["x"], c["temp_g"]["y"], "r"); a.fill_between(c["temp_g"]["x"], c["temp_g"]["lo"], c["temp_g"]["hi"], color="r", alpha=.2)
        if k == "depth": a.set_xscale("log")
        a.set_title(f"{sp} {lab}")
plt.tight_layout(); plt.savefig(f"{out}/quicklook_curves.png", dpi=80)
