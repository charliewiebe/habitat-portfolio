# Stage 01: pull and clean kit fox records from the iNaturalist API
source("pipeline/00_config.R")
suppressPackageStartupMessages({ library(jsonlite); library(data.table) })

`%||%` <- function(a, b) if (is.null(a)) b else a

inat_page <- function(id_above) {
  url <- sprintf(paste0("https://api.inaturalist.org/v1/observations?taxon_id=%d",
                        "&quality_grade=%s&captive=false&order_by=id&order=asc",
                        "&per_page=200&id_above=%d"),
                 CFG$inat_taxon_id, CFG$inat_quality, id_above)
  for (i in 1:5) {
    res <- tryCatch(fromJSON(url, simplifyVector = FALSE), error = function(e) NULL)
    if (!is.null(res)) return(res$results)
    Sys.sleep(2^i)
  }
  stop("iNat API failed at id_above = ", id_above)
}

parse_obs <- function(o) {
  xy <- o$geojson$coordinates
  if (is.null(xy)) return(NULL)
  data.table(
    id               = o$id,
    observed_on      = o$observed_on %||% NA_character_,
    lon              = as.numeric(xy[[1]]),
    lat              = as.numeric(xy[[2]]),
    pos_acc          = as.numeric(o$positional_accuracy %||% NA),
    public_acc       = as.numeric(o$public_positional_accuracy %||% NA),
    obscured         = isTRUE(o$obscured),
    geoprivacy       = o$geoprivacy %||% NA_character_,
    taxon_geoprivacy = o$taxon_geoprivacy %||% NA_character_,
    in_oregon_place  = CFG$inat_oregon_place %in% unlist(o$place_ids),
    taxon            = o$taxon$name %||% NA_character_,
    url              = o$uri %||% NA_character_
  )
}

recs <- list(); last_id <- 0L
repeat {
  page <- inat_page(last_id)
  if (!length(page)) break
  recs[[length(recs) + 1]] <- rbindlist(lapply(page, parse_obs))
  last_id <- page[[length(page)]]$id
  message(sprintf("%d records (last id %d)", sum(vapply(recs, nrow, 1L)), last_id))
  Sys.sleep(1)
}
raw <- rbindlist(recs)
raw[, year := as.integer(substr(observed_on, 1, 4))]
saveRDS(raw, file.path(CFG$data_dir, "occ_raw.rds"))

steps <- list()
log_step <- function(dt, label) steps[[length(steps) + 1]] <<- data.table(step = label, n = nrow(dt))

d <- copy(raw);                                  log_step(d, "Research grade, wild, with coordinates")
if (!CFG$keep_obscured) d <- d[obscured == FALSE]; log_step(d, "Coordinates not obscured")
d <- d[is.na(pos_acc) | pos_acc <= CFG$max_pos_acc_m]
log_step(d, sprintf("Positional accuracy <= %s m or unstated", CFG$max_pos_acc_m))
d <- d[!is.na(year) & year >= CFG$min_year];     log_step(d, sprintf("Observed %d or later", CFG$min_year))
d <- unique(d, by = c("lon", "lat", "observed_on")); log_step(d, "Exact duplicates removed")

saveRDS(d, file.path(CFG$data_dir, "occ_clean.rds"))
fwrite(d[, .(id, observed_on, lon, lat, pos_acc, taxon, in_oregon_place, url)],
       file.path(CFG$data_dir, "occurrences_clean.csv"))
# Oregon records are never used for training; kept (obscured included) as an external check
saveRDS(raw[in_oregon_place == TRUE], file.path(CFG$data_dir, "occ_oregon_all.rds"))
fwrite(rbindlist(steps), file.path(CFG$data_dir, "occ_cleaning_log.csv"))
print(rbindlist(steps))
