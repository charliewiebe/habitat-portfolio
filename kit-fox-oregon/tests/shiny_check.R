# Headless check of the Shiny app: drives every input with shiny::testServer and renders every output.
# Usage, from the project root after the pipeline has run: Rscript tests/shiny_check.R
suppressPackageStartupMessages(library(shiny))
md <- readRDS("data/model_data.rds")
testServer("shiny", {
  session$setInputs(set = "no_sjv", ssp = "ssp245", period = "baseline", layer = "networks", show_obs = TRUE,
                    suitable_q = 0.1, core_q = 0.5, density = 0.08, mvp = 100, dispersal = 15, unavail = 0.5,
                    excl_extrap = TRUE, metric = "largest_network_K", rw_set = "all", rw_layer = "index",
                    imp_set = "no_sjv", rc_set = "no_sjv", rc_var = md$vars[1], recompute = 0)
  cat("baseline:", output$vb_suit, "|", output$vb_k, "|", output$vb_viable, "|", output$vb_verdict, "\n")
  for (per in c("2021-2040", "2081-2100")) for (s in c("ssp245", "ssp585")) {
    session$setInputs(period = per, ssp = s)
    cat(s, per, ":", output$vb_suit, "|", output$vb_k, "|", output$vb_verdict, "\n")
  }
  for (lay in c("index", "agree", "extrap", "unavail", "pop", "networks")) { session$setInputs(layer = lay); session$flushReact() }
  for (st in c("all", "balanced", "all_nobias")) { session$setInputs(set = st); session$flushReact() }
  session$setInputs(rw_layer = "effort", rw_set = "no_sjv"); session$flushReact()
  session$setInputs(set = "no_sjv", mvp = 50, dispersal = 25, recompute = 1); session$flushReact()
  for (m in c("largest_network_K", "suitable_km2", "n_viable", "hist_range_km2")) { session$setInputs(metric = m); invisible(output$ts_plot) }
  invisible(list(output$cv_tbl, output$imp_plot, output$or_tbl, output$yr_plot, output$clog_tbl,
                 output$state_tbl, output$map_title, output$map, output$rw_map))
  session$setInputs(rc_var = "bio6")
  stopifnot(grepl("Min temperature of coldest month", output$imp_plot),
            grepl("Min temperature of coldest month", output$rc_plot))
  cat("PASS: all outputs rendered\n")
})
