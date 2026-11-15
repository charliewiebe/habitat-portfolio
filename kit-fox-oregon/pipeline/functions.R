# Shared functions: models, indices, evaluation, extrapolation, viability
suppressPackageStartupMessages({ library(terra); library(data.table); library(mgcv); library(glmnet) })
terraOptions(progress = 0)
VAR_LABELS <- c(
  bio1 = "Annual mean temperature (°C)", bio2 = "Mean daily temperature range (°C)",
  bio3 = "Isothermality (daily range / annual range, %)", bio4 = "Temperature seasonality (SD × 100)",
  bio5 = "Max temperature of warmest month (°C)", bio6 = "Min temperature of coldest month (°C)",
  bio7 = "Annual temperature range (°C)", bio8 = "Mean temperature of wettest quarter (°C)",
  bio9 = "Mean temperature of driest quarter (°C)", bio10 = "Mean temperature of warmest quarter (°C)",
  bio11 = "Mean temperature of coldest quarter (°C)", bio12 = "Annual precipitation (mm)",
  bio13 = "Precipitation of wettest month (mm)", bio14 = "Precipitation of driest month (mm)",
  bio15 = "Precipitation seasonality (CV of monthly totals, %)", bio16 = "Precipitation of wettest quarter (mm)",
  bio17 = "Precipitation of driest quarter (mm)", bio18 = "Precipitation of warmest quarter (mm)",
  bio19 = "Precipitation of coldest quarter (mm)", sand = "Soil sand content, 30–60 cm",
  clay = "Soil clay content, 30–60 cm", tri = "Terrain ruggedness index",
  bare = "Bare ground cover (fraction)", shrubs = "Shrub cover (fraction)",
  grassland = "Grassland cover (fraction)", trees = "Tree cover (fraction)",
  cropland = "Cropland cover (fraction)", built = "Built-up cover (fraction)",
  log_pop = "Population density (log people/km²)", log_travel = "Travel time to a town of 20,000+ (log minutes)")
var_label <- function(v) ifelse(v %in% names(VAR_LABELS), VAR_LABELS[v], v)
ALBERS <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +datum=WGS84 +units=m"

# ---- weighted ECDF used for the presence-referenced habitat index -------------------------
make_ref <- function(eta, w = NULL) {
  if (is.null(w)) w <- rep(1, length(eta))
  o <- order(eta)
  list(x = eta[o], cw = cumsum(w[o]) / sum(w))
}
pct <- function(x, ref) {
  i <- findInterval(x, ref$x)
  out <- numeric(length(x)); out[i > 0] <- ref$cw[i[i > 0]]
  out[is.na(x)] <- NA
  out
}

# ---- lasso features: linear, quadratic, forward/reverse hinges (knots from background) ----
make_feature_spec <- function(df, vars, n_knots) {
  lapply(setNames(vars, vars), function(v) {
    x <- df[[v]]; mu <- mean(x); s <- sd(x)
    if (!is.finite(s) || s == 0) s <- 1
    z <- (x - mu) / s
    k <- if (n_knots > 0) unique(quantile(z, seq(0.1, 0.9, length.out = n_knots), names = FALSE)) else numeric(0)
    list(mu = mu, s = s, knots = k)
  })
}
build_X <- function(df, spec, hinge = TRUE) {
  do.call(cbind, lapply(names(spec), function(v) {
    sp <- spec[[v]]; z <- (df[[v]] - sp$mu) / sp$s
    m <- cbind(z, z^2); colnames(m) <- paste0(v, c("_l", "_q"))
    if (hinge && length(sp$knots)) {
      hf <- vapply(sp$knots, function(k) pmax(0, z - k), numeric(length(z)))
      hr <- vapply(sp$knots, function(k) pmax(0, k - z), numeric(length(z)))
      hf <- matrix(hf, nrow = length(z)); hr <- matrix(hr, nrow = length(z))
      colnames(hf) <- paste0(v, "_hf", seq_along(sp$knots))
      colnames(hr) <- paste0(v, "_hr", seq_along(sp$knots))
      m <- cbind(m, hf, hr)
    }
    m
  }))
}

# ---- model fitting: infinitely weighted logistic regression == Poisson point process --------
# pa = 1 presence / 0 background, w = weights (background weights large).
fit_gam <- function(d, vars, bias, k) {
  terms <- c(sprintf("s(%s, k = %d, bs = 'ts')", vars, k), sprintf("s(%s, k = 4)", bias))
  f <- as.formula(paste("pa ~", paste(terms, collapse = " + ")))
  withCallingHandlers(   # balanced set uses non-integer presence weights; that is intended
    bam(f, data = d, family = binomial(), weights = d$w, method = "fREML",
        discrete = TRUE, nthreads = max(1, parallel::detectCores() - 1)),
    warning = function(w) if (grepl("non-integer", conditionMessage(w))) invokeRestart("muffleWarning"))
}
fit_lasso <- function(d, vars, bias, n_knots) {
  bg <- d[d$pa == 0, ]
  spec_e <- make_feature_spec(bg, vars, n_knots)
  X <- build_X(d, spec_e); pf <- rep(1, ncol(X)); spec_b <- NULL
  if (length(bias)) {                       # effort terms enter unpenalised so they absorb effort
    spec_b <- make_feature_spec(bg, bias, 0)
    Xb <- build_X(d, spec_b, hinge = FALSE); X <- cbind(X, Xb); pf <- c(pf, rep(0, ncol(Xb)))
  }
  cv <- cv.glmnet(X, d$pa, family = "binomial", weights = d$w, alpha = 1,
                  penalty.factor = pf, nfolds = 5, nlambda = 40, lambda.min.ratio = 1e-3,
                  type.measure = "deviance")
  list(cv = cv, lambda = cv$lambda.min, spec_e = spec_e, spec_b = spec_b)
}

predict_eta <- function(algo, fit, nd, fix_bias) {
  if (length(fix_bias)) for (b in names(fix_bias)) nd[[b]] <- fix_bias[[b]]
  if (algo == "gam") return(as.numeric(predict(fit, nd, type = "link")))
  X <- build_X(nd, fit$spec_e)
  if (!is.null(fit$spec_b)) X <- cbind(X, build_X(nd, fit$spec_b, hinge = FALSE))
  as.numeric(predict(fit$cv, X, s = fit$lambda, type = "link"))
}

# Fit both algorithms; build presence-referenced indices and the ensemble reference
fit_set <- function(d, vars, bias, cfg) {
  m <- list(vars = vars, bias = bias)
  m$fix_bias <- if (length(bias)) vapply(bias, function(b) median(d[[b]][d$pa == 0]), 1) else NULL
  m$gam   <- fit_gam(d, vars, bias, cfg$gam_k)
  m$lasso <- fit_lasso(d, vars, bias, cfg$n_knots)
  pres <- d[d$pa == 1, ]
  m$ref_gam   <- make_ref(predict_eta("gam",   m$gam,   pres, m$fix_bias), pres$w)
  m$ref_lasso <- make_ref(predict_eta("lasso", m$lasso, pres, m$fix_bias), pres$w)
  raw <- (pct(predict_eta("gam", m$gam, pres, m$fix_bias), m$ref_gam) +
          pct(predict_eta("lasso", m$lasso, pres, m$fix_bias), m$ref_lasso)) / 2
  m$ref_ens <- make_ref(raw, pres$w)
  m
}
predict_index <- function(m, nd) {
  ig <- pct(predict_eta("gam",   m$gam,   nd, m$fix_bias), m$ref_gam)
  il <- pct(predict_eta("lasso", m$lasso, nd, m$fix_bias), m$ref_lasso)
  data.table(gam = ig, lasso = il, ens = pct((ig + il) / 2, m$ref_ens))
}
# Observer-effort multiplier: how much more often the fitted model expects records here
# than at the reference effort level (only defined for bias-corrected sets)
effort_multiplier <- function(m, nd) {
  if (!length(m$bias)) return(rep(NA_real_, nrow(nd)))
  exp(predict_eta("gam", m$gam, nd, NULL) - predict_eta("gam", m$gam, nd, m$fix_bias))
}

# ---- evaluation ------------------------------------------------------------------------------
auc <- function(p, b) {
  r <- rank(c(p, b)); n1 <- length(p); n0 <- length(b)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
# Continuous Boyce index (Hirzel et al. 2006, Ecol Model 199:142)
boyce <- function(p, b, n_win = 101, width = 0.1) {
  rng <- range(c(p, b), na.rm = TRUE); w <- width * diff(rng)
  if (w <= 0) return(NA_real_)
  mids <- seq(rng[1] + w / 2, rng[2] - w / 2, length.out = n_win)
  pe <- vapply(mids, function(x) {
    e <- mean(b >= x - w / 2 & b <= x + w / 2)
    if (e == 0) NA_real_ else mean(p >= x - w / 2 & p <= x + w / 2) / e
  }, 1)
  ok <- !is.na(pe)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(mids[ok], pe[ok], method = "spearman"))
}
# Spatial blocks in an Albers equal-area projection, dealt to folds so presence counts stay balanced
assign_folds <- function(lon, lat, is_pres, block_km, k, seed) {
  xy <- terra::project(cbind(lon, lat), from = "EPSG:4326", to = ALBERS)
  bid <- paste(floor(xy[, 1] / (block_km * 1000)), floor(xy[, 2] / (block_km * 1000)))
  set.seed(seed)
  cnt <- tapply(is_pres, bid, sum)
  ord <- names(sort(cnt + runif(length(cnt), 0, 0.5), decreasing = TRUE))
  fold_of <- setNames(rep_len(c(seq_len(k), rev(seq_len(k))), length(ord)), ord)
  unname(fold_of[bid])
}

# ---- extrapolation: MESS (Elith et al. 2010, Methods Ecol Evol 1:330) -----------------------
mess <- function(nd, ref, vars) {
  s <- vapply(vars, function(v) {
    r <- ref[[v]]; x <- nd[[v]]; mn <- min(r); mx <- max(r); rr <- mx - mn
    f <- 100 * ecdf(r)(x)
    ifelse(f == 0, (x - mn) / rr * 100,
    ifelse(f <= 50, 2 * f,
    ifelse(f < 100, 2 * (100 - f), (mx - x) / rr * 100)))
  }, numeric(nrow(nd)))
  s <- matrix(s, nrow = nrow(nd))
  list(mess = apply(s, 1, min), mod = vars[apply(s, 1, which.min)])
}

# ---- reintroduction viability ---------------------------------------------------------------
# Networks: 8-connected patches, merged when any cell of one lies within d_km (centre to centre,
# equal-area distance) of a cell of the other. Returns a raster of network ids.
link_patches <- function(pat, d_km) {
  pv <- values(pat, mat = FALSE); cells <- which(!is.na(pv))
  if (d_km <= 0 || length(unique(pv[cells])) < 2) return(pat)
  xy  <- terra::project(xyFromCell(pat, cells), from = crs(pat), to = ALBERS)
  pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y"), crs = ALBERS)
  nb  <- sf::st_is_within_distance(pts, dist = d_km * 1000)
  pid <- pv[cells]
  e <- unique(cbind(pid[rep(seq_along(nb), lengths(nb))], pid[unlist(nb)]))
  e <- e[e[, 1] < e[, 2], , drop = FALSE]
  ids <- sort(unique(pid))
  g <- igraph::graph_from_edgelist(cbind(match(e[, 1], ids), match(e[, 2], ids)), directed = FALSE)
  g <- igraph::add_vertices(g, length(ids) - igraph::vcount(g))
  comp <- igraph::components(g)$membership
  out <- rast(pat); out[cells] <- comp[match(pid, ids)]
  out
}

# idx: habitat index raster; unavail: cropland+built fraction; extrap_ok: logical raster or NULL
# region / hist: 1 inside Oregon / inside the historical range (Harney + Malheur), NA elsewhere
viability <- function(idx, unavail, extrap_ok, region, hist, area, p) {
  potential <- idx >= p$suitable_q
  if (isTRUE(p$exclude_extrapolation) && !is.null(extrap_ok)) potential <- potential & extrap_ok
  avail <- mask(potential & (unavail < p$unavailable_frac), region)
  pot_v <- values(mask(potential, region), mat = FALSE)
  av_v  <- values(avail, mat = FALSE)
  a_v   <- values(area, mat = FALSE)
  i_v   <- values(idx, mat = FALSE)
  h_v   <- values(hist, mat = FALSE)
  is_av <- which(av_v == 1)
  base <- data.table(
    potential_km2 = sum(a_v[which(pot_v == 1)]),
    suitable_km2  = sum(a_v[is_av]),
    core_km2      = sum(a_v[is_av][i_v[is_av] >= p$core_q]),
    hist_range_km2 = sum(a_v[is_av][!is.na(h_v[is_av])]),
    n_patches = 0L, n_networks = 0L, largest_network_km2 = 0, largest_network_K = 0,
    n_viable = 0L, K_viable_total = 0)
  suit <- classify(avail, cbind(0, NA))
  if (!length(is_av)) return(list(metrics = base, net = suit, tab = data.table()))
  pat <- patches(suit, directions = 8)
  base$n_patches <- length(unique(na.omit(values(pat, mat = FALSE))))
  net <- link_patches(pat, p$dispersal_km)
  tab <- data.table(net = values(net, mat = FALSE), km2 = a_v)[!is.na(net), .(km2 = sum(km2)), by = net]
  tab[, `:=`(K = km2 * p$density_km2)][, viable := K >= p$mvp_n]
  base[, `:=`(n_networks = nrow(tab), largest_network_km2 = max(tab$km2),
              largest_network_K = max(tab$K), n_viable = sum(tab$viable),
              K_viable_total = sum(tab$K[tab$viable]))]
  list(metrics = base, net = net, tab = tab)
}

layer_meta <- function(L) {
  if (L == "baseline") return(data.table(layer = L, ssp = "baseline", period = "1970-2000", gcm = "baseline"))
  s <- strsplit(L, "__", fixed = TRUE)[[1]]
  data.table(layer = L, ssp = s[1], period = s[2], gcm = s[3])
}

outlook_series <- function(idx_stack, mess_stack, static, p) {
  rbindlist(lapply(names(idx_stack), function(L) {
    ok <- if (!is.null(mess_stack)) mess_stack[[L]] >= 0 else NULL
    v <- viability(idx_stack[[L]], static[["unavail"]], ok, static[["oregon"]],
                   static[["hist_range"]], static[["area_km2"]], p)
    cbind(layer_meta(L), v$metrics)
  }))
}
