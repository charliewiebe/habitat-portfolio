"""Bottom-temperature climatology + delta-GAM habitat models for each species.

Usage: python -I pipeline/02_models.py <data dir>
Writes <data dir>/model_outputs.json and <data dir>/model_report.txt
"""
import sys
import json
import warnings
import numpy as np
import pandas as pd
from pygam import LinearGAM, LogisticGAM, s, te, f
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
DATA = sys.argv[1]
SPECIES = ["lingcod", "canary", "sablefish"]
REF_YEARS = list(range(2003, 2025))      # baseline = average year effect over the survey record
DEPTH_LIM = (55.0, 1280.0)               # sampled depth range of the survey design
T_GRID = np.round(np.arange(2.0, 14.01, 0.1), 2)
N_DRAWS = 60
rng = np.random.default_rng(42)

hauls = pd.read_csv(f"{DATA}/hauls_oregon.csv")
grid = pd.read_csv(f"{DATA}/grid_oregon.csv")
years = sorted(hauls.year.unique())
year_code = {y: i for i, y in enumerate(years)}
hauls["yc"] = hauls.year.map(year_code)
report = []

def log(msg):
    print(msg)
    report.append(msg)

def lin_parts(gam, X):
    """Per-term contributions on the link scale; asserts they sum to the linear predictor."""
    parts = [gam.partial_dependence(term=i, X=X) for i in range(len(gam.terms) - 1)]
    icpt = gam.coef_[-1]
    eta = gam._linear_predictor(X)
    assert np.allclose(np.sum(parts, axis=0) + icpt, eta, atol=1e-6)
    return parts, icpt

# ---- 1. bottom temperature climatology ---------------------------------
th = hauls.dropna(subset=["bottom_temp_c"]).copy()
Xt = np.column_stack([np.log(th.depth_m), th.x / 1000, th.y / 1000, th.yc])
tgam = LinearGAM(s(0, n_splines=20) + te(1, 2, n_splines=[6, 14]) + f(3)).gridsearch(
    Xt, th.bottom_temp_c.values, lam=np.logspace(-2, 3, 8), progress=False)
log(f"Temperature GAM: n={len(th)}, R2={tgam.statistics_['pseudo_r2']['explained_deviance']:.3f}")

gdep = grid.depth_m.clip(*DEPTH_LIM).values
def grid_X(yc):
    return np.column_stack([np.log(gdep), grid.x / 1000, grid.y / 1000, np.full(len(grid), yc)])

tclim = np.mean([tgam.predict(grid_X(year_code[y])) for y in REF_YEARS if y in year_code], axis=0)
grid["t_clim"] = tclim
log(f"  grid bottom-temp climatology {tclim.min():.2f}-{tclim.max():.2f} C")

# temperature anomaly by year (domain mean of year effect), for context in the app
yr_eff = {int(y): float(tgam.partial_dependence(term=2, X=grid_X(year_code[y])[:1])[0]) for y in years}
mean_eff = np.mean(list(yr_eff.values()))
temp_anom = {y: round(v - mean_eff, 3) for y, v in yr_eff.items()}

# ---- 2. delta-GAM habitat models ---------------------------------------
def design(df):
    return np.column_stack([np.log(df.depth_m.clip(*DEPTH_LIM)), df.bottom_temp_c,
                            df.x / 1000, df.y / 1000, df.yc])

def terms(ns_depth, ns_t, ns_xy):
    return (s(0, n_splines=ns_depth) + s(1, n_splines=ns_t) +
            te(2, 3, n_splines=list(ns_xy)) + f(4))

LAMS = np.logspace(-1, 3, 6)
md = hauls.dropna(subset=["bottom_temp_c"]).copy()
X = design(md)
lat_band = pd.qcut(md.y, 5, labels=False).values
kfold = rng.integers(0, 5, len(md))
outputs = {"species": {}, "t_grid": T_GRID.tolist(), "temp_anomaly_by_year": temp_anom}

for sp in SPECIES:
    y_pos = (md[f"{sp}_kg"] > 0).astype(int).values
    cpue = md[f"{sp}_cpue"].values
    n_pos = int(y_pos.sum())
    small = n_pos < 800
    tp = terms(10, 8, (5, 10) if not small else (4, 8))
    tg = terms(8, 6, (5, 8) if not small else (4, 6))

    # delta-lognormal: Bernoulli GAM for encounter, Gaussian GAM on log(CPUE | positive)
    pgam = LogisticGAM(tp).gridsearch(X, y_pos, lam=LAMS, progress=False)
    Xg_pos = X[y_pos == 1]
    ggam = LinearGAM(tg).gridsearch(Xg_pos, np.log(cpue[y_pos == 1]), lam=LAMS, progress=False)
    sigma2 = float(ggam.statistics_["scale"])            # residual variance on log scale

    # ---- cross-validation (lambda fixed at the full-data choice) ----
    def cv(fold_ids):
        pred_p = np.zeros(len(md)); pred_d = np.zeros(len(md))
        for k in np.unique(fold_ids):
            tr, te_ = fold_ids != k, fold_ids == k
            pm = LogisticGAM(tp, lam=pgam.lam).fit(X[tr], y_pos[tr])
            trp = tr & (y_pos == 1)
            gm = LinearGAM(tg, lam=ggam.lam).fit(X[trp], np.log(cpue[trp]))
            s2 = gm.statistics_["scale"]
            pred_p[te_] = pm.predict_mu(X[te_])
            pred_d[te_] = pred_p[te_] * np.exp(np.clip(gm.predict(X[te_]), -20, 20) + s2 / 2)
        auc = roc_auc_score(y_pos, pred_p)
        rho = spearmanr(pred_d, cpue).correlation
        return auc, rho

    auc_r, rho_r = cv(kfold)
    auc_s, rho_s = cv(lat_band)
    dev_p = pgam.statistics_["pseudo_r2"]["explained_deviance"]
    dev_g = ggam.statistics_["pseudo_r2"]["explained_deviance"]
    log(f"{sp}: n={len(md)}, positives={n_pos}; dev. explained presence={dev_p:.3f}, positive={dev_g:.3f}")
    log(f"   random 5-fold   AUC={auc_r:.3f}  Spearman(density)={rho_r:.3f}")
    log(f"   latitude-block  AUC={auc_s:.3f}  Spearman(density)={rho_s:.3f}")

    # temperature response, clamped to the 1st-99th pct of temperatures sampled
    t_lo, t_hi = np.percentile(md.bottom_temp_c, [1, 99])
    def t_curve(gam):
        Tc = np.clip(T_GRID, t_lo, t_hi)
        Xc = np.column_stack([np.full(len(Tc), np.log(200.0)), Tc, np.full(len(Tc), grid.x.mean() / 1000),
                              np.full(len(Tc), grid.y.mean() / 1000), np.zeros(len(Tc))])
        return gam.partial_dependence(term=1, X=Xc)

    fp_t, fg_t = t_curve(pgam), t_curve(ggam)

    # posterior draws of the temperature smooth (MVN on coefficients) for uncertainty bands
    def t_draws(gam, n=N_DRAWS):
        Tc = np.clip(T_GRID, t_lo, t_hi)
        Xc = np.column_stack([np.full(len(Tc), np.log(200.0)), Tc, np.full(len(Tc), grid.x.mean() / 1000),
                              np.full(len(Tc), grid.y.mean() / 1000), np.zeros(len(Tc))])
        B = gam._modelmat(Xc, term=1).toarray()
        idx = gam.terms.get_coef_indices(1)
        cov = gam.statistics_["cov"][np.ix_(idx, idx)]
        beta = rng.multivariate_normal(gam.coef_[idx], cov, size=n, method="eigh")
        cur = beta @ B.T
        assert np.allclose(B @ gam.coef_[idx], gam.partial_dependence(term=1, X=Xc), atol=1e-6)
        return cur

    fp_draws, fg_draws = t_draws(pgam), t_draws(ggam)
    tc_idx = np.clip(np.round((grid.t_clim.values - T_GRID[0]) / 0.1).astype(int), 0, len(T_GRID) - 1)

    # ---- grid predictions: per-year linear predictors without the temperature term ----
    def eta_by_year(gam):
        out = []
        for y in REF_YEARS:
            if y not in year_code:
                continue
            Xg = np.column_stack([np.log(gdep), grid.t_clim, grid.x / 1000, grid.y / 1000,
                                  np.full(len(grid), year_code[y])])
            parts, icpt = lin_parts(gam, Xg)
            out.append(icpt + parts[0] + parts[2] + parts[3])   # depth + space + year (no temp)
        return np.array(out)

    # Average over years on the response scale so the baseline is an arithmetic mean year.
    # Positive part: exp() factorises, so log-mean-exp is exact for any temperature.
    # Encounter part: chosen so the year-averaged probability is exact at climatology.
    ep_y, eg_y = eta_by_year(pgam), eta_by_year(ggam)
    eta_g = np.log(np.mean(np.exp(eg_y), axis=0)) + sigma2 / 2      # lognormal mean correction
    p_bar = np.mean(1 / (1 + np.exp(-(ep_y + fp_t[tc_idx]))), axis=0)
    p_bar = np.clip(p_bar, 1e-6, 1 - 1e-6)
    eta_p = np.log(p_bar / (1 - p_bar)) - fp_t[tc_idx]

    # baseline density at climatology (what the app reproduces at a 0 C anomaly)
    dens = (1 / (1 + np.exp(-(eta_p + fp_t[tc_idx])))) * np.exp(eta_g + fg_t[tc_idx])

    # depth response curves (link scale, with 95% CI) for the methods panel
    dgrid = np.exp(np.linspace(np.log(DEPTH_LIM[0]), np.log(DEPTH_LIM[1]), 60))
    Xd = np.column_stack([np.log(dgrid), np.full(60, 6.5), np.full(60, grid.x.mean() / 1000),
                          np.full(60, grid.y.mean() / 1000), np.zeros(60)])
    pd_d, ci_d = pgam.partial_dependence(term=0, X=Xd, width=0.95)
    Tq = np.linspace(t_lo, t_hi, 50)
    Xq = np.column_stack([np.full(50, np.log(200.0)), Tq, np.full(50, grid.x.mean() / 1000),
                          np.full(50, grid.y.mean() / 1000), np.zeros(50)])
    pd_t, ci_t = pgam.partial_dependence(term=1, X=Xq, width=0.95)
    gd_t, gci_t = ggam.partial_dependence(term=1, X=Xq, width=0.95)

    obs_by_year = md.assign(pos=y_pos).groupby("year").agg(enc=("pos", "mean"),
                                                            cpue=(f"{sp}_cpue", "mean"))
    outputs["species"][sp] = {
        "n_hauls": int(len(md)), "n_pos": n_pos, "sigma2_log": round(sigma2, 4),
        "metrics": {"dev_presence": round(dev_p, 3), "dev_positive": round(dev_g, 3),
                    "auc_random": round(auc_r, 3), "rho_random": round(rho_r, 3),
                    "auc_block": round(auc_s, 3), "rho_block": round(rho_s, 3)},
        "t_range": [round(float(t_lo), 2), round(float(t_hi), 2)],
        "eta_p": np.round(eta_p, 3).tolist(),
        "eta_g": np.round(eta_g, 3).tolist(),
        "fp_t": np.round(fp_t, 4).tolist(),
        "fg_t": np.round(fg_t, 4).tolist(),
        "fp_draws": np.round(fp_draws, 3).tolist(),
        "fg_draws": np.round(fg_draws, 3).tolist(),
        "dens_check_mean": float(np.average(dens, weights=grid.area_ha)),
        "curves": {
            "depth": {"x": np.round(dgrid, 1).tolist(), "y": np.round(pd_d, 3).tolist(),
                      "lo": np.round(ci_d[:, 0], 3).tolist(), "hi": np.round(ci_d[:, 1], 3).tolist()},
            "temp_p": {"x": np.round(Tq, 2).tolist(), "y": np.round(pd_t, 3).tolist(),
                       "lo": np.round(ci_t[:, 0], 3).tolist(), "hi": np.round(ci_t[:, 1], 3).tolist()},
            "temp_g": {"x": np.round(Tq, 2).tolist(), "y": np.round(gd_t, 3).tolist(),
                       "lo": np.round(gci_t[:, 0], 3).tolist(), "hi": np.round(gci_t[:, 1], 3).tolist()},
        },
        "obs_by_year": {int(k): [round(float(v.enc), 3), round(float(v.cpue), 3)]
                        for k, v in obs_by_year.iterrows()},
    }
    log(f"   grid mean density {outputs['species'][sp]['dens_check_mean']:.3f} kg/ha; "
        f"observed mean CPUE {cpue.mean():.3f} kg/ha; temp range used {t_lo:.2f}-{t_hi:.2f} C")

grid["t_clim"] = np.round(grid.t_clim, 3)
outputs["grid"] = {c: np.round(grid[c].values, 1 if c in ("x", "y", "depth_m") else 4).tolist()
                   for c in ["x", "y", "lon", "lat", "depth_m", "area_ha", "t_clim"]}
outputs["years"] = [int(y) for y in years]
with open(f"{DATA}/model_outputs.json", "w") as fh:
    json.dump(outputs, fh)
with open(f"{DATA}/model_report.txt", "w") as fh:
    fh.write("\n".join(report) + "\n")
grid.to_csv(f"{DATA}/grid_oregon.csv", index=False)
