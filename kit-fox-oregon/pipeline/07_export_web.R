# Stage 07: build the self-contained web app. Rasters are quantised to 8 bits and embedded.
# Writes app/index.html (body fragment, as published to claude.ai) and index.html (a full
# standalone document for GitHub Pages with links to the code and the portfolio).
source("pipeline/00_config.R"); source("pipeline/functions.R")
suppressPackageStartupMessages(library(jsonlite))

O <- function(f) file.path(CFG$out_dir, f)
gpkg <- file.path(CFG$data_dir, "env", "boundaries.gpkg")
# 0-1 values are stored as floor(v * 240) so that every slider threshold (steps of 0.05) and every
# share of 5 GCMs falls exactly on a step: v >= t in R is then q >= t * 240 in the browser. 255 = no data.
Q8 <- 240
q8 <- function(v) { x <- floor(pmin(pmax(v, 0), 1) * Q8 + 1e-9); x[is.na(x)] <- 255; base64_enc(as.raw(x)) }
b8 <- function(v) base64_enc(as.raw(ifelse(is.na(v), 0L, 1L)))

static <- rast(O("window_static.tif"))
e <- as.vector(ext(static))
area_m <- matrix(values(static[["area_km2"]], mat = FALSE), nrow = nrow(static), byrow = TRUE)

lines_of <- function(v, tol, box) {
  v <- simplifyGeom(crop(v, box), tolerance = tol, preserveTopology = TRUE)
  if (!nrow(v)) return(list())
  g <- geom(as.lines(v))
  unname(lapply(split(seq_len(nrow(g)), paste(g[, "geom"], g[, "part"])),
                function(i) round(g[i, c("x", "y"), drop = FALSE], 3)))
}
adm <- vect(gpkg, layer = "admin1"); hr <- vect(gpkg, layer = "hist_range")
wbox <- ext(e[1] - 0.5, e[2] + 0.5, e[3] - 0.5, e[4] + 0.5)

set_labels <- c(all = "All records", balanced = "Region-balanced", no_sjv = "Exclude San Joaquin",
                all_nobias = "No effort correction")
sets <- lapply(setNames(CFG$model_sets, CFG$model_sets), function(s) {
  idx <- rast(O(sprintf("window_index_%s.tif", s))); ms <- rast(O(sprintf("window_mess_%s.tif", s)))
  L <- list(baseline = list(idx = q8(values(idx[["baseline"]], mat = FALSE)),
                            ext = q8(values(ms[["baseline"]] < 0, mat = FALSE))))
  for (sp in paste0("ssp", CFG$ssps)) for (p in CFG$periods) {
    ln <- grep(paste0("^", sp, "__", p, "__"), names(idx), value = TRUE)
    if (!length(ln)) next
    L[[paste(sp, p, sep = "__")]] <- list(idx = q8(values(mean(idx[[ln]]), mat = FALSE)),
                                          ext = q8(values(mean(ms[[ln]] < 0), mat = FALSE)),
                                          n_gcm = length(ln))
  }
  list(label = set_labels[[s]], layers = L)
})

# GCM range bands at the pipeline's default assumptions
ts <- fread(O("timeseries.csv"))
mets <- c("largest_network_K", "suitable_km2", "core_km2", "hist_range_km2", "n_viable")
bands <- lapply(split(ts, ts$set), function(d) {
  lapply(split(d[ssp != "baseline"], d[ssp != "baseline"]$ssp), function(x)
    lapply(setNames(mets, mets), function(m) x[, .(year = year[1], lo = min(get(m)), med = median(get(m)), hi = max(get(m))), by = period][order(year)]))
})

md  <- readRDS(file.path(CFG$data_dir, "model_data.rds"))
occ <- md$occ
inwin <- occ[in_or == FALSE & lon > e[1] & lon < e[2] & lat > e[3] & lat < e[4], .(lon = round(lon, 3), lat = round(lat, 3))]
orc <- fread(O("oregon_check.csv"))
or_obs <- unique(orc[, .(id, date = observed_on, lon = round(lon, 3), lat = round(lat, 3), obscured)])

rw  <- aggregate(rast(O("rangewide_index.tif")), 2, mean, na.rm = TRUE)
eff <- aggregate(rast(O("rangewide_effort.tif")), 2, mean, na.rm = TRUE)
eff_v <- values(eff, mat = FALSE); eff_lim <- 2   # log10 multiplier clamped to +/- 2 (100x)

cv <- fread(O("cv_summary.csv")); imp <- fread(O("importance.csv"))

D <- list(
  meta = list(q8 = Q8, window = e, nrow = nrow(static), ncol = ncol(static), res = res(static)[1],
              ssps = paste0("ssp", CFG$ssps), periods = CFG$periods, period_mid = as.list(CFG$period_mid),
              gcms = CFG$gcms, vars = md$vars, var_labels = as.list(VAR_LABELS),
              defaults = CFG[c("suitable_q", "core_q", "density_km2", "mvp_n", "dispersal_km",
                               "unavailable_frac", "exclude_extrapolation")],
              n_presence_cells = lapply(md$sets, `[[`, "n_pres"), n_records = nrow(occ[in_or == FALSE]),
              generated = format(Sys.Date())),
  rowArea = round(area_m[, 1], 4),
  oregon = b8(values(static[["oregon"]], mat = FALSE)),
  hist = b8(values(static[["hist_range"]], mat = FALSE)),
  unavail = q8(values(static[["unavail"]], mat = FALSE)),
  sets = sets, bands = bands,
  lines = list(states = lines_of(adm, 0.01, wbox), oregon = lines_of(adm[adm$NAME_1 == "Oregon"], 0.01, wbox),
               hist = lines_of(hr, 0.005, wbox)),
  obs = list(train = as.matrix(inwin), oregon = or_obs),
  rangewide = list(window = as.vector(ext(rw)), nrow = nrow(rw), ncol = ncol(rw),
                   sets = lapply(setNames(names(rw), names(rw)), function(s) q8(values(rw[[s]], mat = FALSE))),
                   effort = q8((pmin(pmax(eff_v, -eff_lim), eff_lim) + eff_lim) / (2 * eff_lim)), effort_lim = eff_lim,
                   lines = lines_of(adm, 0.05, ext(rw)),
                   train = as.matrix(occ[in_or == FALSE, .(lon = round(lon, 2), lat = round(lat, 2), sjv = in_sjv)])),
  cv = cv[algo == "ens", .(set, auc = round(auc, 3), boyce = round(boyce, 3), p10 = round(test_pres_above_p10, 3))],
  imp = lapply(split(imp, imp$set), function(d) d[order(-share), .(variable, share = round(share, 3))]),
  orcheck = dcast(orc, id + observed_on + obscured ~ set, value.var = "index_mean")
)

html <- paste(readLines("app/template.html", warn = FALSE, encoding = "UTF-8"), collapse = "\n")
html <- sub("/*__DATA__*/null", toJSON(D, auto_unbox = TRUE, digits = NA, na = "null"), html, fixed = TRUE)
LINKS_SLOT <- "    <!--__LINKS__-->\n"
write_utf8 <- function(x, f) { con <- file(f, "wb"); writeBin(charToRaw(enc2utf8(x)), con); close(con) }
write_utf8(sub(LINKS_SLOT, "", html, fixed = TRUE), "app/index.html")

REPO <- "https://github.com/charliewiebe/habitat-portfolio/tree/main/kit-fox-oregon"
links <- sprintf('    <nav class="links" aria-label="Project links"><a href="%s">Code</a><a href="../">Portfolio</a></nav>\n', REPO)
cut  <- regexpr("</style>", html, fixed = TRUE) + nchar("</style>")
head <- sub('<meta name="claude:viewport-fit-sides" content="cover">\n', "", substr(html, 1, cut - 1), fixed = TRUE)
body <- sub(LINKS_SLOT, links, substr(html, cut, nchar(html)), fixed = TRUE)
reset <- paste0("<style>:root{padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}",
                "body{margin:0}img{max-width:100%}[hidden]{display:none!important}</style>\n")
doc <- paste0('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n',
              '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n',
              reset, head, "\n</head>\n<body>", body, "\n</body>\n</html>\n")
write_utf8(doc, "index.html")
message(sprintf("wrote index.html (%.1f MB) and app/index.html", file.size("index.html") / 1e6))
