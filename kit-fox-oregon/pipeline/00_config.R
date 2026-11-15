# Kit fox (Vulpes macrotis) habitat suitability and reintroduction outlook for Oregon
# Shared configuration. Every stage sources this file; edit here, then re-run from the
# earliest stage your change affects (see README). Run every stage from the project root.

CFG <- list(

  # ---- occurrences: iNaturalist API v1 --------------------------------------------------
  inat_taxon_id     = 42070L,       # Vulpes macrotis (descendant taxa included)
  inat_oregon_place = 10L,          # iNat place id for Oregon
  inat_quality      = "research",
  max_pos_acc_m     = 2000,         # drop records with stated accuracy worse than this (NA kept)
  keep_obscured     = FALSE,        # obscured coords are random within a 0.2 deg cell (~5x the grid)
  min_year          = 1990,

  # ---- grid and extent ---------------------------------------------------------------------
  res_min   = 2.5,                        # arc-minutes; 2.5 or 5 (WorldClim, CMIP6, GPW all offer both)
  extent    = c(-125, -96.5, 20, 46.5),   # xmin, xmax, ymin, ymax: whole range + all of Oregon
  window    = c(-125, -114, 40, 46.5),    # focal window for projections (OR + N NV + SW ID + NE CA)
  data_dir  = "data",
  raw_dir   = "data/raw",                 # geodata download cache (global rasters; large)
  out_dir   = "outputs",
  keep_raw  = FALSE,                      # TRUE keeps the global CMIP6 files after cropping

  # ---- environmental predictors ------------------------------------------------------------
  # Greedy selection walks this list in order and keeps a variable only if |r| < max_abs_cor
  # with everything already kept (r computed over the training background).
  predictor_priority = c("bio6", "bio5", "bio12", "sand", "tri", "bare", "shrubs", "bio15",
                         "bio4", "cropland", "bio2", "clay", "grassland", "trees",
                         "bio18", "bio19", "bio1", "bio7", "bio14", "bio17"),
  max_abs_cor        = 0.7,
  max_predictors     = 10,
  soil_depth_cm      = 60,                # SoilGrids 30-60 cm layer: den-relevant depth
  landcover_vars     = c("bare", "shrubs", "grassland", "trees", "cropland", "built"),

  # ---- observer-effort (sampling-bias) covariates ------------------------------------------
  # Fitted alongside environment, then held at a constant when predicting, so the projected
  # surface reflects environment only (Warton, Renner & Ramp 2013, PLoS ONE 8:e79168).
  bias_vars   = c("log_pop", "log_travel"),
  pop_year    = 2020,                     # GPWv4 density year, matches the iNat record era
  travel_size = 7,                        # travel time to settlements >= 20k people (up = TRUE)

  # ---- climate futures (WorldClim 2.1 CMIP6 downscaled bioclim) ----------------------------
  gcms    = c("ACCESS-CM2", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", "UKESM1-0-LL"),
  ssps    = c("245", "585"),
  periods = c("2021-2040", "2041-2060", "2061-2080", "2081-2100"),
  period_mid = c(baseline = 1985, "2021-2040" = 2030, "2041-2060" = 2050,
                 "2061-2080" = 2070, "2081-2100" = 2090),

  # ---- training design ---------------------------------------------------------------------
  accessible_buffer_km = 200,             # background drawn within this distance of presences
  n_background  = 30000,                  # area-proportional background sample
  bg_weight     = 100,                    # IWLR background weight (Poisson point-process limit)
  gam_k         = 5,
  n_knots       = 6,                      # hinge knots per variable for the lasso model
  # Model sets. "all_nobias" exists only to show what the effort correction changes.
  model_sets = c("all", "balanced", "no_sjv", "all_nobias"),
  # Cismontane California (San Joaquin Valley, Carrizo Plain, Salinas/Cuyama valleys).
  # The no_sjv set drops presences AND background inside it, so the model learns the
  # cold-desert niche relevant to SE Oregon rather than the Mediterranean-lowland one.
  sjv_polygon = rbind(c(-122.2, 38.2), c(-121.0, 38.4), c(-119.9, 37.6), c(-119.0, 36.4),
                      c(-118.6, 35.4), c(-118.75, 34.85), c(-119.6, 34.7), c(-120.7, 34.9),
                      c(-121.6, 36.0), c(-122.2, 37.2), c(-122.2, 38.2)),
  balance_power = 0.5,                    # balanced set: presence weight ∝ n_state^-0.5
  cv_block_km   = 250,
  cv_folds      = 5,
  seed          = 42,

  # ---- habitat and reintroduction viability ------------------------------------------------
  # Habitat index = share of training presences whose predicted suitability is <= the cell's,
  # so 0.10 is the P10 rule and 0.50 means "as good as the median occupied site".
  suitable_q            = 0.10,
  core_q                = 0.50,
  unavailable_frac      = 0.50,           # cropland + built fraction above which a cell is lost
  exclude_extrapolation = TRUE,           # drop cells with MESS < 0 from habitat totals
  # 0.04-0.12 foxes/km2: Great Basin, Dugway UT (WNAN 67:389, 2007)
  density_km2           = 0.08,
  # SE Oregon: N ~ 28, Ne ~ 14.7 (Wildlife Soc Bull, Dec 2018) -> Ne/N ~ 0.5;
  # Ne >= 50 needs N ~ 100.
  mvp_n                 = 100,
  dispersal_km          = 15              # patches closer than this form one network (assumption)
)

dir.create(CFG$data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(CFG$out_dir,  showWarnings = FALSE, recursive = TRUE)
