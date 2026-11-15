# Oregon Groundfish Habitat

Where lingcod, canary rockfish and sablefish live off Oregon, and how bottom temperature, low-oxygen water and fishing rules change their spawning potential.

**Open the app:** https://charliewiebe.github.io/habitat-portfolio/oregon-groundfish-habitat/

**Portfolio:** https://charliewiebe.github.io/habitat-portfolio/

Habitat is modelled from 4,002 tows with bottom-temperature records, out of 4,072 tows of NOAA's West Coast Groundfish Bottom Trawl Survey (2003–2024) between the California border and the Columbia River. Habitat is predicted onto the survey's own 2-km sampling grid (7,858 Oregon cells). A per-recruit model calibrated to the latest Pacific Fishery Management Council stock assessments adds the fishing layer. Users can change bottom temperature, a low-oxygen scenario, bag and size limits, sport and commercial effort, and closures: the Yelloweye Rockfish Conservation Areas in 50 CFR 660.70, Oregon's five marine reserves, and a fathom-band closure. They can then see where habitat and spawning potential gain or lose relative to 2026 rules and average conditions.

## Files

```
index.html                  the app (built by pipeline/05_build_app.py; self-contained, no server)
pipeline/00_get_data.sh     fetch NWFSC haul, catch and grid tables from the public surveyjoin repos
pipeline/01_prepare.py      Oregon haul table (zero-filled catch, CPUE in kg/ha) and prediction grid, UTM 10N
pipeline/02_models.py       bottom-temperature climatology GAM; delta-lognormal habitat GAMs; cross-validation
pipeline/03_geometry.py     coastline, depth contours, YRCA and marine-reserve polygons digitized from regulation text
pipeline/04_params.py       life history, stock status, 2026 ODFW rules, labelled modelling assumptions
pipeline/05_build_app.py    inlines data, engine and UI into index.html
pipeline/diagnostics/       observed-vs-predicted by depth band, quick-look maps and response curves
app/engine.js               scenario math: per-recruit model, bag-limit arithmetic, calibration, density
app/app.js, template.html   map, controls, charts
data/                       outputs of 01–04 used by the app (hauls, grid, model outputs, geometry, parameters)
tests/                      JS engine vs independent Python per-recruit implementation; headless smoke test
```

## Rebuild

```
cd oregon-groundfish-habitat
pip install -r requirements.txt
bash pipeline/00_get_data.sh
python -I pipeline/01_prepare.py raw/surveyjoin-data raw/surveyjoin/data data
python -I pipeline/02_models.py  data
python -I pipeline/03_geometry.py data
python -I pipeline/04_params.py   data
python -I pipeline/05_build_app.py data app index.html
node tests/engine_check.js && python -I tests/engine_check.py data/species_params.json
```

## Model

**Habitat.** A delta-lognormal GAM for each species. One part models whether a tow catches the species; the other models log catch per hectare when it does. Both use s(log depth) + s(bottom temperature) + te(easting, northing) + a year effect, with smoothness chosen by grid search. The baseline is the average survey year on the response scale. A separate GAM fitted to tow temperatures (91.5% deviance explained) gives each grid cell a bottom-temperature climatology. The app's anomaly slider shifts that climatology. The temperature effect is held flat outside the 3.2–9.0 °C range the survey sampled, and 95% ranges come from 60 posterior draws of the temperature smooth.

| Species | Positive tows | AUC, random 5-fold | AUC, latitude blocks | Spearman (density), latitude blocks |
|---|---|---|---|---|
| Lingcod | 1,519 | 0.880 | 0.850 | 0.592 |
| Canary rockfish | 444 | 0.867 | 0.836 | 0.363 |
| Sablefish | 3,072 | 0.886 | 0.857 | 0.556 |

**Fishing.** A female equilibrium per-recruit model. Natural mortality, growth, weight-length, maturity and status-quo SPR come from the 2021 northern lingcod assessment and the 2023 canary rockfish and sablefish assessments. Current fishing mortality is calibrated so open water under 2026 rules reproduces the assessment SPR. Bag limits act through a negative-binomial count of legal fish encountered per angler-trip. Fish released over the bag or under the size limit die at the release-mortality rate. Closures remove the chosen fleet's mortality in covered cells, and displaced effort can move to open water, capped at double. The sport share, encounters per trip, release mortality, selectivity and the lingcod maturity slope are modelling assumptions. The app labels them as assumptions and lets you change them.

## Limits

- The survey does not sample shallower than 55 m or untrawlable rock, so nearshore reef habitat is under-represented.
- Cold bottom water on the Oregon shelf is usually freshly upwelled, low-oxygen water, so the temperature effect partly carries an oxygen signal. Read warming results as an association.
- Oxygen is not in this dataset. The low-oxygen control is a scenario overlay, not a fitted effect.
- Assessment estimates are coastwide (lingcod: north of 40°10′ N) and are applied to Oregon waters.
- Each cell's spawning ratio is computed locally, which assumes adults stay put.

## Sources

- NOAA NWFSC West Coast Groundfish Bottom Trawl Survey, via [surveyjoin](https://github.com/DFO-NOAA-Pacific/surveyjoin) and [surveyjoin-data](https://github.com/DFO-NOAA-Pacific/surveyjoin-data) v0.0.4
- Taylor et al. 2021, *Status of lingcod along the northern U.S. west coast in 2021* (PFMC)
- Langseth et al. 2023, *Status of canary rockfish along the U.S. West Coast in 2023* (PFMC)
- Johnson et al. 2023, *Status of sablefish along the U.S. West Coast in 2023* (PFMC)
- ODFW 2026 sport groundfish regulations; ODFW Stonewall Bank YRCA closure
- 50 CFR 660.70(g)–(m) (YRCA coordinates); OAR 141-142 (marine reserve boundaries)
- Coastline: GLOBE 30″ land mask via `global-land-mask`

Public data only.
