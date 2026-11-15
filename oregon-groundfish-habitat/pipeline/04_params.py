"""Life-history, stock status, 2026 Oregon sport rules, and per-recruit assumptions.

Every value carries its source. Values marked "assumption" are modelling choices exposed
in the app's Advanced panel, not published estimates.
Usage: python -I pipeline/04_params.py <data dir>
"""
import sys
import json

DATA = sys.argv[1]
IN = 2.54

LINGCOD_SRC = "Taylor et al. 2021, Status of lingcod along the northern U.S. west coast in 2021 (PFMC)"
CANARY_SRC = "Langseth et al. 2023, Status of canary rockfish along the U.S. West Coast in 2023 (PFMC)"
SABLE_SRC = "Johnson et al. 2023, Status of sablefish along the U.S. West Coast in 2023 (PFMC)"
ODFW_2026 = "ODFW 2026 sport groundfish regulations (myodfw.com/sport-bottomfish-seasons)"

params = {
    "lingcod": {
        "name": "Lingcod", "sci": "Ophiodon elongatus", "niche": "Shelf, rocky reef to ~300 m",
        "lh": {
            "M": 0.418, "amax": 25,
            "growth": {"a1": 0.5, "L1": 12.668, "a2": 14.0, "L2": 105.084, "k": 0.152},
            "wl": {"a": 2.802e-6, "b": 3.2766},
            "mat": {"type": "age", "x50": 3.23, "slope": -1.5},
        },
        "status": {"spr": round(1 - 0.262 * 0.55, 3), "year": 2020, "target": 0.45,
                   "frac_unfished": 0.642, "frac_year": 2021,
                   "basis": "Relative fishing intensity (1-SPR)/(1-SPR45%) = 0.262 in 2020, northern stock"},
        "regs": {"bag": 3, "bag_prior": 2, "bag_label": "Daily bag limit", "bag_max": 6,
                 "size_in": 22, "size_max_in": 30,
                 "note": "3 fish daily (2 in 2022-2025), 22-inch minimum; separate from the general marine bag."},
        "rec": {"share": 0.60, "sel50": 45.0, "selw": 5.0, "mu": 1.2, "theta": 1.0, "dmort": 0.10},
        "com": {"sel50": 60.0, "selw": 5.0},
        "sources": {
            "M, growth, weight-length, A50": LINGCOD_SRC + " (Table 15; sections 2.3.2, 2.3.4, 3.2.5.1)",
            "status": LINGCOD_SRC + " (Tables iv, v)",
            "regs": ODFW_2026,
        },
        "assumed": ["maturity slope (assessment uses age-based maturity, slope not reported)",
                    "recreational share of current mortality", "selectivity curves",
                    "legal encounters per angler-trip (mean, dispersion)", "release mortality"],
    },
    "canary": {
        "name": "Canary rockfish", "sci": "Sebastes pinniger", "niche": "Outer shelf, 80-250 m, schooling",
        "lh": {
            "M": 0.078, "amax": 80,
            "growth": {"a1": 1.0, "L1": 11.38, "a2": None, "Linf": 57.9, "k": 0.143},
            "wl": {"a": 1.19e-5, "b": 3.09},
            "mat": {"type": "age", "x50": 10.87, "slope": -0.69},
        },
        "status": {"spr": round(1 - 0.469, 3), "year": 2022, "target": 0.50,
                   "frac_unfished": 0.351, "frac_year": 2023,
                   "basis": "Fishing intensity 1-SPR = 0.469 in 2022, near the SPR50% target"},
        "regs": {"bag": 2, "bag_prior": 1, "bag_label": "Daily sub-bag (within 4-fish marine bag)", "bag_max": 4,
                 "size_in": 0, "size_max_in": 0,
                 "note": "2-fish sub-bag inside the 4-fish general marine bag (1 fish in 2025); no minimum size."},
        "rec": {"share": 0.35, "sel50": 30.0, "selw": 3.0, "mu": 0.8, "theta": 0.5, "dmort": 0.20},
        "com": {"sel50": 35.0, "selw": 3.0},
        "sources": {
            "M": CANARY_SRC + " (female M of 0.0779 quoted for the base model in the maturity-ramp comparison; male fixed at 0.0643)",
            "growth, weight-length, functional maturity": CANARY_SRC + " (sections 2.4.2, 2.4.5, 2.4.6)",
            "status": CANARY_SRC + " (Table iv; reference points)",
            "regs": ODFW_2026,
        },
        "assumed": ["recreational share of current mortality", "selectivity curves",
                    "legal encounters per angler-trip (mean, dispersion)", "release mortality"],
    },
    "sablefish": {
        "name": "Sablefish", "sci": "Anoplopoma fimbria", "niche": "Upper slope, 200-1,200 m",
        "lh": {
            "M": 0.071, "amax": 70,
            "growth": {"a1": 0.5, "L1": 25.262, "a2": 30.0, "L2": 61.130, "k": 0.367},
            "wl": {"a": 3.315e-6, "b": 3.273},
            "mat": {"type": "length", "x50": 55.19, "slope": -0.421},
        },
        "status": {"spr": round(1 - 0.531 * 0.55, 3), "year": 2022, "target": 0.45,
                   "frac_unfished": 0.630, "frac_year": 2023,
                   "basis": "Relative fishing intensity (1-SPR)/(1-SPR45%) = 0.531 in 2022"},
        "regs": {"bag": 10, "bag_prior": 10, "bag_label": "Daily bag limit", "bag_max": 15,
                 "size_in": 0, "size_max_in": 0,
                 "note": "10 fish daily; no minimum size. Mostly a commercial fishery."},
        "rec": {"share": 0.03, "sel50": 40.0, "selw": 4.0, "mu": 0.4, "theta": 0.5, "dmort": 0.20},
        "com": {"sel50": 45.0, "selw": 4.0},
        "sources": {
            "M, growth, maturity": SABLE_SRC + " (Table 2; female length at age 30 = 61.1 cm)",
            "weight-length": "Haltuch et al. 2019 benchmark values, carried into the 2023 model",
            "status": SABLE_SRC + " (Tables iv, v)",
            "regs": ODFW_2026,
        },
        "assumed": ["reference age for L1 (0.5 yr)", "recreational share of current mortality",
                    "selectivity curves", "legal encounters per angler-trip (mean, dispersion)", "release mortality"],
    },
}
for p in params.values():
    p["regs"]["size_cm"] = round(p["regs"]["size_in"] * IN, 2)

with open(f"{DATA}/species_params.json", "w") as fh:
    json.dump(params, fh, indent=1)
print("wrote species_params.json")
