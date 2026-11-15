import sys, json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D, out = sys.argv[1], sys.argv[2]
geo = json.load(open(f"{D}/geometry.json")); mo = json.load(open(f"{D}/model_outputs.json"))
fig, ax = plt.subplots(1, 2, figsize=(12, 12))
for a, (xl, yl) in zip(ax, [((280e3, 480e3), (4640e3, 5140e3)), ((360e3, 440e3), (4880e3, 4980e3))]):
    a.scatter(mo["grid"]["x"], mo["grid"]["y"], c=np.log10(mo["grid"]["depth_m"]), s=2 if a is ax[0] else 30, marker="s", cmap="Blues")
    for r in geo["land"]: r = np.array(r); a.fill(r[:, 0], r[:, 1], color="#ddd", ec="k", lw=.5)
    for c in geo["contours"]: p = np.array(c["pts"]); a.plot(p[:, 0], p[:, 1], lw=.6, color={100: "c", 200: "b", 500: "m", 1000: "k"}[c["depth"]])
    for c in geo["closures"]: p = np.array(c["ring"] + c["ring"][:1]); a.plot(p[:, 0], p[:, 1], "r" if c["kind"] == "yrca" else "g", lw=1)
    for pt in geo["ports"]: a.plot(pt["x"], pt["y"], "ko", ms=3); a.text(pt["x"] + 2000, pt["y"], pt["name"], fontsize=7)
    a.set_xlim(*xl); a.set_ylim(*yl); a.set_aspect("equal")
plt.tight_layout(); plt.savefig(f"{out}/quicklook_geo.png", dpi=70)
