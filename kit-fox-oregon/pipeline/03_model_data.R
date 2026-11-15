# Stage 03: presences, accessible area, background, predictor selection, CV folds
source("pipeline/00_config.R"); source("pipeline/functions.R")
set.seed(CFG$seed)

env  <- rast(file.path(CFG$data_dir, "env", "current.tif"))
gpkg <- file.path(CFG$data_dir, "env", "boundaries.gpkg")
adm  <- vect(gpkg, layer = "admin1")
occ  <- readRDS(file.path(CFG$data_dir, "occ_clean.rds"))

tmpl   <- env[[1]]
land   <- values(!is.na(tmpl), mat = FALSE) == 1
area_v <- values(cellSize(tmpl, unit = "km"), mat = FALSE)
or_r   <- rasterize(adm[adm$NAME_1 == "Oregon"], tmpl, field = 1)
sjv_r  <- rasterize(vect(CFG$sjv_polygon, type = "polygons", crs = "EPSG:4326"), tmpl, field = 1)
or_v   <- !is.na(values(or_r, mat = FALSE))
sjv_v  <- !is.na(values(sjv_r, mat = FALSE))

# ---- presences ------------------------------------------------------------------------------
occ[, cell := cellFromXY(tmpl, cbind(lon, lat))]
occ <- occ[!is.na(cell) & land[cell]]
st  <- extract(adm, vect(occ, geom = c("lon", "lat"), crs = "EPSG:4326"))
st  <- st[!duplicated(st$id.y), ]
occ[, state := fifelse(is.na(st$NAME_1), "Unassigned", st$NAME_1)]
occ[, in_or  := or_v[cell] | in_oregon_place]
occ[, in_sjv := sjv_v[cell] | taxon %chin% "Vulpes macrotis mutica"]
train <- occ[in_or == FALSE]

cells_of <- function(dt) dt[, .(state = state[1], n_records = .N, in_sjv = any(in_sjv)), by = cell]
pres <- list(
  all      = cells_of(train),
  balanced = cells_of(train),
  no_sjv   = cells_of(train[in_sjv == FALSE])
)
pres$all_nobias <- pres$all
for (s in names(pres)) pres[[s]][, w := 1]
n_state <- pres$balanced[, .N, by = state]
pres$balanced[n_state, w := i.N^(-CFG$balance_power), on = "state"]
pres$balanced[, w := w / mean(w)]

# ---- accessible area and background ---------------------------------------------------------
accessible <- function(cells, drop_sjv) {
  r <- rast(tmpl); r[cells] <- 1
  d <- values(distance(r), mat = FALSE)
  which(land & d <= CFG$accessible_buffer_km * 1000 & !or_v & !(drop_sjv & sjv_v))
}
bg_for <- function(mcells) {
  if (length(mcells) <= CFG$n_background)
    return(data.table(cell = mcells, w = CFG$bg_weight * area_v[mcells] / mean(area_v[mcells])))
  data.table(cell = sample(mcells, CFG$n_background, prob = area_v[mcells]), w = CFG$bg_weight)
}
M_cells <- list(all = accessible(pres$all$cell, FALSE), no_sjv = accessible(pres$no_sjv$cell, TRUE))
M_cells$balanced <- M_cells$all; M_cells$all_nobias <- M_cells$all
bg_all <- bg_for(M_cells$all)
bg <- list(all = bg_all, balanced = bg_all, all_nobias = bg_all, no_sjv = bg_for(M_cells$no_sjv))

# ---- predictor selection on the "all" background --------------------------------------------
env_at <- function(cells) as.data.table(env[cells])
X <- env_at(bg$all$cell)
cand <- intersect(CFG$predictor_priority, names(X))
cand <- cand[vapply(cand, function(v) sd(X[[v]], na.rm = TRUE) > 0, TRUE)]
C <- abs(cor(X[, ..cand]))
vars <- character(0)
for (v in cand) {
  if (length(vars) >= CFG$max_predictors) break
  if (!length(vars) || all(C[v, vars] < CFG$max_abs_cor)) vars <- c(vars, v)
}
message("Predictors: ", paste(vars, collapse = ", "))

# ---- assemble modelling tables --------------------------------------------------------------
sets <- lapply(setNames(CFG$model_sets, CFG$model_sets), function(s) {
  p <- pres[[s]]; b <- bg[[s]]
  d <- rbind(data.table(cell = p$cell, pa = 1, w = p$w), data.table(cell = b$cell, pa = 0, w = b$w))
  xy <- xyFromCell(tmpl, d$cell)
  d[, `:=`(lon = xy[, 1], lat = xy[, 2])]
  d <- cbind(d, env_at(d$cell)[, c(vars, CFG$bias_vars), with = FALSE])
  d[, fold := assign_folds(lon, lat, pa == 1, CFG$cv_block_km, CFG$cv_folds, CFG$seed)]
  list(d = as.data.frame(d), bias = if (s == "all_nobias") character(0) else CFG$bias_vars,
       n_pres = nrow(p), n_records = sum(p$n_records))
})

M_r <- rast(tmpl, nlyrs = 2); names(M_r) <- c("all", "no_sjv")
M_r[[1]][M_cells$all] <- 1; M_r[[2]][M_cells$no_sjv] <- 1
writeRaster(M_r, file.path(CFG$out_dir, "accessible_area.tif"), overwrite = TRUE)

occ_out <- occ[, .(id, observed_on, year, lon, lat, state, taxon, in_or, in_sjv, url)]
saveRDS(list(vars = vars, cor = C, sets = sets, occ = occ_out,
             n_state = n_state[order(-N)]),
        file.path(CFG$data_dir, "model_data.rds"))
fwrite(data.table(variable = vars), file.path(CFG$out_dir, "selected_predictors.csv"))
print(rbindlist(lapply(names(sets), function(s)
  data.table(set = s, presence_cells = sets[[s]]$n_pres, records = sets[[s]]$n_records,
             background = sum(sets[[s]]$d$pa == 0)))))
