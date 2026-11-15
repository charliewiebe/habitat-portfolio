# Stage 06: habitat totals, dispersal networks, carrying capacity vs MVP, through time
source("pipeline/00_config.R"); source("pipeline/functions.R")
static <- rast(file.path(CFG$out_dir, "window_static.tif"))

ts <- rbindlist(lapply(CFG$model_sets, function(s) {
  message("== ", s)
  idx <- rast(file.path(CFG$out_dir, sprintf("window_index_%s.tif", s)))
  ms  <- rast(file.path(CFG$out_dir, sprintf("window_mess_%s.tif", s)))
  cbind(set = s, outlook_series(idx, ms, static, CFG))
}))
ts[, year := CFG$period_mid[ifelse(ssp == "baseline", "baseline", period)]]
fwrite(ts, file.path(CFG$out_dir, "timeseries.csv"))

summ <- ts[, .(gcms = .N,
               suitable_km2 = median(suitable_km2),
               largest_K_median = median(largest_network_K),
               largest_K_range = sprintf("%.0f-%.0f", min(largest_network_K), max(largest_network_K)),
               share_gcms_viable = mean(n_viable > 0)),
           by = .(set, ssp, period)]
fwrite(summ, file.path(CFG$out_dir, "outlook_summary.csv"))
print(summ[set == "no_sjv"])
