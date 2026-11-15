"""Independent per-recruit implementation to cross-check app/engine.js."""
import json, math, sys
P = json.load(open(sys.argv[1]))
def L(g, a):
    if g.get("a2") is None: v = g["Linf"] + (g["L1"] - g["Linf"]) * math.exp(-g["k"] * (a - g["a1"]))
    else: v = g["L1"] + (g["L2"] - g["L1"]) * (1 - math.exp(-g["k"] * (a - g["a1"]))) / (1 - math.exp(-g["k"] * (g["a2"] - g["a1"])))
    return max(v, 1)
def R(mu, th, b):
    from scipy.stats import nbinom
    p = th / (th + mu)
    return sum(1 - nbinom.cdf(k - 1, th, p) for k in range(1, b + 1))
for k, sp in P.items():
    lh = sp["lh"]; n = lh["amax"] + 1
    def sb(Frec, Fcom, bag, size):
        r, c = sp["rec"], sp["com"]; Rb = R(r["mu"], r["theta"], bag); g = (Rb + r["dmort"] * (r["mu"] - Rb)) / r["mu"]
        N, s = 1.0, 0.0
        for a in range(n):
            La, Lm = L(lh["growth"], a), L(lh["growth"], a + .5)
            W = lh["wl"]["a"] * La ** lh["wl"]["b"]; x = a if lh["mat"]["type"] == "age" else La
            mat = 1 / (1 + math.exp(lh["mat"]["slope"] * (x - lh["mat"]["x50"])))
            sr = 1 / (1 + math.exp(-(Lm - r["sel50"]) / r["selw"])); sc = 1 / (1 + math.exp(-(Lm - c["sel50"]) / c["selw"]))
            ret = 1 / (1 + math.exp(-(Lm - size) / 0.75)) if size > 0 else 1
            Z = lh["M"] + Frec * sr * (ret * g + (1 - ret) * r["dmort"]) + Fcom * sc
            if a == n - 1: N = N / (1 - math.exp(-Z))
            s += N * mat * W; N *= math.exp(-Z)
        return s
    s0 = sb(0, 0, sp["regs"]["bag"], sp["regs"]["size_cm"])
    Rb0 = R(sp["rec"]["mu"], sp["rec"]["theta"], sp["regs"]["bag"]); g0 = (Rb0 + sp["rec"]["dmort"] * (sp["rec"]["mu"] - Rb0)) / sp["rec"]["mu"]
    lo, hi = 0, 5
    for _ in range(80):
        mid = (lo + hi) / 2; s_ = sp["rec"]["share"]
        if sb(mid * s_ / g0, mid * (1 - s_), sp["regs"]["bag"], sp["regs"]["size_cm"]) / s0 > sp["status"]["spr"]: lo = mid
        else: hi = mid
    th = (lo + hi) / 2; s_ = sp["rec"]["share"]
    fr, fc = th * s_ / g0, th * (1 - s_)
    b0 = sp["regs"]["bag"]
    print(k, "theta", round(th, 4), "SPR bag-1", round(sb(fr, fc, b0 - 1, sp["regs"]["size_cm"]) / s0, 4), "nosize", round(sb(fr, fc, b0, 0) / s0, 4))
