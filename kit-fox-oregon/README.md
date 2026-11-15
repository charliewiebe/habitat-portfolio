# Kit Fox Oregon Outlook

Is there a connected block of habitat in Oregon large enough to hold a viable kit fox (*Vulpes macrotis*) population, now and as the climate changes to 2100?

**Open the app:** https://charliewiebe.github.io/habitat-portfolio/kit-fox-oregon/

**Portfolio:** https://charliewiebe.github.io/habitat-portfolio/

Kit foxes are listed as threatened in Oregon. The state is the northern edge of their range, in the salt-desert scrub of Harney and Malheur counties. The habitat model is trained on research-grade iNaturalist records from across the range outside Oregon: 714 records after cleaning, pulled on 2026-11-11. It corrects for observer effort, because records pile up where people live. It is then projected onto Oregon for 1970–2000 and for four CMIP6 periods to 2100, under five GCMs and two emissions scenarios. Oregon's own records are held out as an independent check. Suitable habitat is then grouped into dispersal-linked networks. Each network's carrying capacity is compared with a minimum viable population. In the app you can change the training set, scenario, period, suitability threshold, fox density, minimum viable population, dispersal distance and land-use cutoff, and the verdict, habitat totals and outlook chart recompute in the browser.

## Files

```
index.html                    the app (built by pipeline/07_export_web.R; self-contained, no server)
run_all.R                     runs stages 01–07 in order
pipeline/00_config.R          every setting: data sources, grid, predictors, GCMs, training sets, viability defaults
pipeline/01_occurrences.R     iNaturalist API pull and cleaning log; Oregon records set aside
pipeline/02_env_layers.R      WorldClim 2.1 + CMIP6 bioclim, SoilGrids, terrain ruggedness, ESA WorldCover,
                              GPWv4 population, travel time, GADM boundaries; aligned to a 2.5′ grid
pipeline/03_model_data.R      presence cells, accessible area, background sample, predictor selection, CV folds
pipeline/04_fit_models.R      GAM + lasso point-process models per training set; spatial block CV; importance; response curves
pipeline/05_project.R         range-wide baseline; Oregon window × GCM × SSP × period; MESS; Oregon record check
pipeline/06_habitat_viability.R  suitable and core habitat, dispersal networks, carrying capacity vs MVP, time series
pipeline/07_export_web.R      quantises the rasters and inlines them into app/template.html
pipeline/functions.R          models, habitat index, evaluation, MESS, network linking, viability
app/template.html             the web app's markup, styles and script
shiny/                        R Shiny + leaflet version with the same analysis, for local use after a pipeline run
data/                         iNaturalist pull (raw and cleaned), cleaning log, cleaned records as CSV
outputs/                      CV, importance, response curves, Oregon check, viability time series (CSV)
tests/                        web app vs pipeline numbers; headless Shiny check
```

Environmental rasters, fitted models and raster outputs (about 90 MB, plus about 2 GB of downloads) are not committed. Stages 02–06 rebuild them.

## Rebuild

```
cd kit-fox-oregon
Rscript run_all.R                  # or each stage: Rscript pipeline/01_occurrences.R, ...
node tests/page_check.js index.html
Rscript tests/shiny_check.R
Rscript -e 'shiny::runApp("shiny")'
```

R packages: `terra sf geodata mgcv glmnet igraph data.table jsonlite shiny bslib leaflet plotly DT`. Run every script from this folder.

On a 2-core machine, stage 02 takes about 40 minutes (about 45 global rasters of 0.2–0.6 GB). Stages 04 and 06 take about 10 minutes each, and the rest a few minutes. Stage 01 queries the live API, so a later run picks up newer records.

## Model

**Records.** Research-grade, wild iNaturalist records go through four filters:
- Records with obscured coordinates are dropped. iNaturalist randomises those within a 0.2° cell, about five times the grid size. This removes 93% of San Joaquin Valley records, which are obscured because that subspecies is federally endangered.
- Records with stated accuracy worse than 2 km are dropped.
- Records from before 1990 are dropped.
- Exact duplicates are removed.

1,306 records become 715 (see `data/occ_cleaning_log.csv`), reduced to one presence per grid cell.

**Predictors.** Predictors are chosen in a fixed priority order, keeping each one only if |r| < 0.7 with those already kept:

| Predictor | Source |
|---|---|
| Min temperature of coldest month | WorldClim 2.1 (bio6) |
| Max temperature of warmest month | WorldClim 2.1 (bio5) |
| Annual precipitation | WorldClim 2.1 (bio12) |
| Precipitation seasonality | WorldClim 2.1 (bio15) |
| Mean daily temperature range | WorldClim 2.1 (bio2) |
| Soil sand content, 30–60 cm | SoilGrids |
| Terrain ruggedness index | from 30″ elevation |
| Bare ground, shrub and cropland cover | ESA WorldCover 2020 |

**Observer effort.** Two effort covariates are fitted alongside the environment: log population density (GPWv4, 2020) and log travel time to a town of 20,000 or more. When predicting, they are held at their background medians, so the habitat surface is conditional on equal effort everywhere (Warton, Renner & Ramp 2013). The "No effort correction" training set shows what the adjustment changes. Without it, baseline habitat in Oregon falls from 17,630 to 10,936 km², because rural southeast Oregon is under-observed.

**Models.** A GAM with shrinkage smooths (mgcv) and a lasso with linear, quadratic and hinge features (glmnet) are fitted as infinitely weighted logistic regression, the presence–background equivalent of a Poisson point process. The background is 30,000 area-weighted points within 200 km of presences. Each model's prediction becomes a **habitat index**: the share of training presences predicted at or below a cell. So 0.10 is the P10 threshold, and 0.50 means "as good as the median occupied site". The ensemble averages the two indices.

**Training sets.** The cleaned records are dominated by the hot Mojave and Colorado deserts, while southeast Oregon is cold Great Basin desert. Four training sets show how much that matters:
- All records.
- Region-balanced: presence weights ∝ n_state^−0.5.
- Exclude San Joaquin: drops cismontane California presences and background. This is the default, because it has the best and most stable cross-validation.
- No effort correction.

| Training set | AUC | Boyce | Held-out presences ≥ P10 |
|---|---|---|---|
| All records | 0.819 | 0.716 | 84% |
| Region-balanced | 0.813 | 0.775 | 87% |
| Exclude San Joaquin | 0.829 | 0.884 | 86% |
| No effort correction | 0.833 | 0.836 | 85% |

Ensemble scores from spatial block cross-validation (250 km blocks, 5 folds).

**Change through time.** Bioclim variables are swapped for CMIP6 downscaled values from ACCESS-CM2, IPSL-CM6A-LR, MIROC6, MPI-ESM1-2-HR and UKESM1-0-LL under SSP2-4.5 and SSP5-8.5. Soils, terrain and land cover stay at present values, so the time series is the climate signal alone. Cells where a GCM's climate falls outside the training range (MESS < 0) are excluded from habitat totals.

**Viability.** The steps run in order:
1. A cell is suitable if its habitat index is ≥ 0.10, it is not extrapolated, and cropland plus built cover is under 50%.
2. Suitable cells form 8-connected patches.
3. Patches with any cells within 15 km of each other, centre to centre, join one network.
4. Each network's carrying capacity is K = area × 0.08 foxes/km². That density is the midpoint of the 0.04–0.12 foxes/km² measured in the Great Basin at Dugway, Utah.
5. A network is viable if K ≥ 100. Genetic work on southeast Oregon foxes estimated N ≈ 28 and Ne ≈ 14.7, so Ne/N ≈ 0.5, and reaching Ne ≥ 50 needs about 100 animals.

The 15 km dispersal distance is an assumption. Every assumption can be changed in the app.

## Results (Exclude San Joaquin, default assumptions)

| Period | Scenario | Suitable habitat in Oregon, km² (GCM range) | Largest network, foxes (GCM range) | GCMs with a viable network |
|---|---|---|---|---|
| 1970–2000 | baseline | 17,630 (10,621 in Harney + Malheur) | 752 | 4 viable networks |
| 2041–2060 | SSP2-4.5 | 57,731 (35,585–73,393) | 4,601 (2,184–5,838) | 5 of 5 |
| 2041–2060 | SSP5-8.5 | 63,873 (40,546–83,932) | 5,090 (2,561–6,630) | 5 of 5 |
| 2081–2100 | SSP2-4.5 | 71,462 (42,885–83,519) | 5,695 (2,649–6,593) | 5 of 5 |
| 2081–2100 | SSP5-8.5 | 94,318 (73,879–104,314) | 7,331 (5,818–8,043) | 5 of 5 |

At the baseline, the verdict holds even with a stricter suitability threshold (0.20) and half the density (0.04/km²): the largest network still holds about 125 foxes. No Oregon cell reaches "core" habitat (index ≥ 0.5), which fits a range-edge population. Four of the five Oregon iNaturalist records fall in suitable habitat. The miss is an obscured record near Ontario.

## Limits

- Coyotes, prey (kangaroo rats), cheatgrass invasion and fire are not modelled. Coyote density was linked to kit fox declines in the Great Basin, and Oregon's population is about 28 animals despite modelled room for hundreds. Habitat area may not be the binding constraint.
- Projected gains are driven by warmer winters at a cold range edge. Some of the new habitat is in central Oregon and the Columbia Basin, outside the species' historical range. Land cover is held fixed, so future conversion and fire would make the outlook worse.
- The climate baseline (1970–2000) predates most records (2015–2026).
- Networks are counted inside Oregon only. Habitat that continues into Nevada and Idaho is drawn on the map but not counted.
- The effort term absorbs any real response of kit foxes to people, such as the urban Bakersfield population.

## Sources

- iNaturalist API v1, research-grade *Vulpes macrotis* observations (taxon 42070), pulled 2026-11-11
- Fick & Hijmans 2017, *WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas*; WorldClim 2.1 CMIP6 downscaled bioclim
- Poggio et al. 2021, *SoilGrids 2.0*; ESA WorldCover 2020; CIESIN GPWv4.11 population density; Nelson et al. 2019, travel time to cities; GADM; all via the `geodata` R package
- Warton, Renner & Ramp 2013, *Model-based control of observer bias for the analysis of presence-only data in ecology*, PLoS ONE
- Hirzel et al. 2006 (continuous Boyce index); Elith, Kearney & Phillips 2010 (MESS)
- *Changes in kit fox–coyote–prey relationships in the Great Basin Desert, Utah*, Western North American Naturalist 67:389 (2007)
- *Genetic characterization of kit foxes at their northern range extent and monitoring recommendations*, Wildlife Society Bulletin (2018)
- ODFW, Oregon Conservation Strategy: kit fox

Public data only.
