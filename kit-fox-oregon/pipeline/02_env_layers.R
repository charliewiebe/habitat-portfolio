# Stage 02: download and align environmental layers (run where geodata.ucdavis.edu,
# www.worldclim.org and files.isric.org are reachable)
source("pipeline/00_config.R")
suppressPackageStartupMessages({ library(terra); library(geodata) })
options(timeout = 3600)   # global rasters take longer than the 60 s default
terraOptions(progress = 0)

env_dir <- file.path(CFG$data_dir, "env")
dir.create(file.path(env_dir, "future"), recursive = TRUE, showWarnings = FALSE)
dir.create(CFG$raw_dir, recursive = TRUE, showWarnings = FALSE)
raw <- CFG$raw_dir
ext_train <- ext(CFG$extent)
need <- function(x, what) { if (is.null(x)) stop("Download failed: ", what); x }

# ---- baseline climate defines the template grid ---------------------------------------------
bio <- crop(need(worldclim_global("bio", res = CFG$res_min, path = raw), "WorldClim bioclim"), ext_train)
names(bio) <- paste0("bio", 1:19)
tmpl <- bio[[1]]
to_tmpl <- function(x) project(crop(x, extend(ext_train, 0.5)), tmpl, method = "average")

# ---- soils (den substrate) ------------------------------------------------------------------
soil <- c(need(soil_world("sand", depth = CFG$soil_depth_cm, path = raw), "SoilGrids sand"),
          need(soil_world("clay", depth = CFG$soil_depth_cm, path = raw), "SoilGrids clay"))
soil <- to_tmpl(soil); names(soil) <- c("sand", "clay")

# ---- terrain ruggedness from 30" elevation, averaged to the grid ----------------------------
elv <- crop(need(worldclim_global("elev", res = 0.5, path = raw), "30s elevation"), extend(ext_train, 0.5))
tri <- project(terrain(elv, "TRI"), tmpl, method = "average"); names(tri) <- "tri"

# ---- land cover fractions (ESA WorldCover 2020) -----------------------------------------------
lc <- rast(lapply(CFG$landcover_vars, function(v) need(landcover(v, path = raw), paste("WorldCover", v))))
lc <- to_tmpl(lc); names(lc) <- CFG$landcover_vars

# ---- observer-effort covariates ---------------------------------------------------------------
pop <- need(population(CFG$pop_year, res = CFG$res_min, path = raw), "GPWv4 population density")
log_pop <- log1p(to_tmpl(pop)); names(log_pop) <- "log_pop"
tt <- need(travel_time("city", size = CFG$travel_size, up = TRUE, path = raw), "travel time")
tt <- classify(tt, cbind(65534.5, Inf, NA))
log_travel <- log1p(to_tmpl(tt)); names(log_travel) <- "log_travel"

current <- c(bio, soil, tri, lc, log_pop, log_travel)
current <- mask(current, app(current, anyNA), maskvalues = 1)
writeRaster(current, file.path(env_dir, "current.tif"), overwrite = TRUE)

# ---- boundaries -------------------------------------------------------------------------------
adm <- rbind(need(gadm("USA", level = 1, path = raw, resolution = 2), "GADM USA"),
             need(gadm("MEX", level = 1, path = raw, resolution = 2), "GADM MEX"))
adm <- crop(adm[, c("GID_0", "NAME_1")], ext_train)
cty <- need(gadm("USA", level = 2, path = raw, resolution = 2), "GADM USA counties")
hist_range <- cty[cty$NAME_1 == "Oregon" & cty$NAME_2 %in% c("Harney", "Malheur"), c("NAME_1", "NAME_2")]
gpkg <- file.path(env_dir, "boundaries.gpkg")
writeVector(adm, gpkg, layer = "admin1", overwrite = TRUE)
writeVector(hist_range, gpkg, layer = "hist_range", insert = TRUE)

# ---- CMIP6 futures (bioclim only; soils, terrain and land cover held at present) -------------
for (g in CFG$gcms) for (s in CFG$ssps) for (p in CFG$periods) {
  f <- file.path(env_dir, "future", sprintf("bioc_%s_ssp%s_%s.tif", g, s, p))
  if (file.exists(f)) next
  x <- cmip6_world(g, s, p, var = "bioc", res = CFG$res_min, path = raw)
  if (is.null(x)) { warning("Missing CMIP6 layer: ", g, " ssp", s, " ", p); next }
  src <- sources(x)
  x <- mask(crop(x, ext_train), current[[1]]); names(x) <- paste0("bio", 1:19)
  writeRaster(x, f, overwrite = TRUE)
  if (!CFG$keep_raw) unlink(src)
  message("saved ", basename(f))
}
