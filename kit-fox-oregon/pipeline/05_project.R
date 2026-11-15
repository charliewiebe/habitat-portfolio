# Stage 05: project every model set range-wide (baseline) and across the Oregon window for the
# baseline and each GCM x SSP x period; MESS extrapolation; Oregon external check
source("pipeline/00_config.R"); source("pipeline/functions.R")
md   <- readRDS(file.path(CFG$data_dir, "model_data.rds")); vars <- md$vars
env  <- rast(file.path(CFG$data_dir, "env", "current.tif"))
gpkg <- file.path(CFG$data_dir, "env", "boundaries.gpkg")
adm  <- vect(gpkg, layer = "admin1"); hr <- vect(gpkg, layer = "hist_range")
bio_vars <- grep("^bio", vars, value = TRUE)

chunked <- function(df, f, n = 50000) {
  idx <- split(seq_len(nrow(df)), ceiling(seq_len(nrow(df)) / n))
  unlist(lapply(idx, function(i) f(df[i, , drop = FALSE])), use.names = FALSE)
}
fill <- function(tmpl, cells, vals) { r <- rast(tmpl); r[cells] <- vals; r }

# ---- Oregon window: static layers -----------------------------------------------------------
env_w  <- crop(env, ext(CFG$window)); tmpl_w <- env_w[[1]]
cells_w <- which(!is.na(values(tmpl_w, mat = FALSE)))
base_w <- as.data.frame(env_w[cells_w])
static <- c(env_w[["cropland"]] + env_w[["built"]], env_w[["log_pop"]],
            rasterize(adm[adm$NAME_1 == "Oregon"], tmpl_w, field = 1),
            rasterize(hr, tmpl_w, field = 1), cellSize(tmpl_w, unit = "km"))
names(static) <- c("unavail", "log_pop", "oregon", "hist_range", "area_km2")
writeRaster(static, file.path(CFG$out_dir, "window_static.tif"), overwrite = TRUE)

# ---- futures (bioclim only; everything else held at present) ---------------------------------
fut <- list()
for (s in CFG$ssps) for (p in CFG$periods) for (g in CFG$gcms) {
  f <- file.path(CFG$data_dir, "env", "future", sprintf("bioc_%s_ssp%s_%s.tif", g, s, p))
  if (!file.exists(f)) { warning("missing ", basename(f)); next }
  fut[[sprintf("ssp%s__%s__%s", s, p, g)]] <- as.data.frame(crop(rast(f), ext(CFG$window))[cells_w])
}
message(length(fut), " future layers")

# ---- range-wide baseline ----------------------------------------------------------------------
tmpl  <- env[[1]]
cells <- which(!is.na(values(tmpl, mat = FALSE)))
full  <- as.data.frame(env[cells])
rw <- list()
for (s in CFG$model_sets) {
  m <- readRDS(file.path("models", paste0(s, ".rds")))
  message("== ", s)
  rw[[s]] <- fill(tmpl, cells, chunked(full, function(x) predict_index(m, x)$ens))
  if (s == "all")
    eff <- fill(tmpl, cells, log10(chunked(full, function(x) effort_multiplier(m, x))))

  ref <- md$sets[[s]]$d[, vars, drop = FALSE]
  idx <- list(baseline = predict_index(m, base_w)$ens)
  ms0 <- mess(base_w, ref, vars)
  mss <- list(baseline = ms0$mess)
  for (L in names(fut)) {
    nd <- base_w; for (v in bio_vars) nd[[v]] <- fut[[L]][[v]]
    idx[[L]] <- predict_index(m, nd)$ens
    mss[[L]] <- mess(nd, ref, vars)$mess
  }
  ir <- rast(lapply(idx, function(v) fill(tmpl_w, cells_w, v))); names(ir) <- names(idx)
  mr <- rast(lapply(mss, function(v) fill(tmpl_w, cells_w, v))); names(mr) <- names(mss)
  writeRaster(ir, file.path(CFG$out_dir, sprintf("window_index_%s.tif", s)), overwrite = TRUE)
  writeRaster(mr, file.path(CFG$out_dir, sprintf("window_mess_%s.tif", s)), overwrite = TRUE)
  mod <- fill(tmpl_w, cells_w, match(ms0$mod, vars))
  writeRaster(mod, file.path(CFG$out_dir, sprintf("window_mod_%s.tif", s)), overwrite = TRUE)
}
rw <- rast(rw); names(rw) <- CFG$model_sets
writeRaster(rw, file.path(CFG$out_dir, "rangewide_index.tif"), overwrite = TRUE)
names(eff) <- "log10_effort_multiplier"
writeRaster(eff, file.path(CFG$out_dir, "rangewide_effort.tif"), overwrite = TRUE)

# ---- external check: Oregon iNat records (never used in training) ---------------------------
oc <- readRDS(file.path(CFG$data_dir, "occ_oregon_all.rds"))
oc <- oc[lon > -124.7 & lon < -116.4 & lat > 41.9 & lat < 46.3]
chk <- rbindlist(lapply(seq_len(nrow(oc)), function(i) {
  r <- oc[i]
  vals <- if (isTRUE(r$obscured)) {
    x0 <- floor(r$lon / 0.2) * 0.2; y0 <- floor(r$lat / 0.2) * 0.2   # iNat 0.2-degree cell
    v <- as.data.frame(crop(rw, ext(x0, x0 + 0.2, y0, y0 + 0.2)))
    data.table(set = names(rw), index_mean = colMeans(v, na.rm = TRUE),
               index_max = apply(v, 2, max, na.rm = TRUE))
  } else {
    v <- unlist(extract(rw, cbind(r$lon, r$lat))[1, ])
    data.table(set = names(rw), index_mean = v, index_max = v)
  }
  cbind(r[, .(id, observed_on, lon, lat, obscured, url)], vals)
}))
fwrite(chk, file.path(CFG$out_dir, "oregon_check.csv"))
print(dcast(chk, id + obscured ~ set, value.var = "index_mean"))
